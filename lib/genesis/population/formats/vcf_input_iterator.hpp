#ifndef GENESIS_POPULATION_FORMATS_VCF_INPUT_ITERATOR_H_
#define GENESIS_POPULATION_FORMATS_VCF_INPUT_ITERATOR_H_

/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2021 Lucas Czech

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact:
    Lucas Czech <lczech@carnegiescience.edu>
    Department of Plant Biology, Carnegie Institution For Science
    260 Panama Street, Stanford, CA 94305, USA
*/

/**
 * @brief
 *
 * @file
 * @ingroup population
 */

#ifdef GENESIS_HTSLIB

#include "genesis/population/formats/hts_file.hpp"
#include "genesis/population/formats/vcf_header.hpp"
#include "genesis/population/formats/vcf_record.hpp"

#include "genesis/utils/core/thread_pool.hpp"

#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace genesis {
namespace population {

// =================================================================================================
//     VCF/BCF Input Iterator
// =================================================================================================

/**
 * @brief Iterate an input source and parse it as a VCF/BCF file.
 *
 * This simple wrapper allows easy iteration through the records/lines of a VCF/BCF file, and
 * takes care of setting up the HtsFile, VcfHeader, and VcfRecord.
 *
 * Basic usage:
 *
 *     auto it = VcfInputIterator( infile );
 *     while( it ) {
 *         // work with it.record() or it->...
 *         ++it;
 *     }
 *
 * or
 *
 *     for( auto it = VcfInputIterator( infile );  it; ++it ) {
 *         // work with it
 *     }
 *
 * For details on working with the records/lines, see VcfRecord and VcfFormatIterator.
 *
 * Caveat: The iterator is an input iterator that traverses a single VCF file in one go.
 * We internally use a buffer to speed up the reading asynchronously in the background.
 * However, as we do not copy the buffer when copying this iterator (for speed reasons, because
 * in typical iterator use cases that is just not necessary), any copy of this iterator might point
 * to old or invalid data. Hence, once an iterator is incremented, all copies of it are invalidated.
 *
 * This also means that the iterator is not thread safe: Incrementing an iterator or any copy of it
 * from multiple tasks leads to undefined behaviour, and hence needs to be synchronized externally.
 *
 * All this is typically not an issue, as an iterator is traversed in a single loop over a file.
 * Any copies are usually mere coincidental while preparing the input etc, but normally they are
 * not concurrently accessed in the actual loop over VcfRecord%s.
 */
class VcfInputIterator
{
public:

    // -------------------------------------------------------------------------
    //     Member Types
    // -------------------------------------------------------------------------

    using self_type         = VcfInputIterator;
    using value_type        = VcfRecord;
    using pointer           = value_type const*;
    using reference         = value_type const&;
    using difference_type   = std::ptrdiff_t;
    using iterator_category = std::input_iterator_tag;

    // -------------------------------------------------------------------------
    //     Constructors and Rule of Five
    // -------------------------------------------------------------------------

    /**
     * @brief Create a default instance, with no input. This is also the past-the-end iterator.
     */
    VcfInputIterator() = default;

    /**
     * @brief Create an instance that reads from an input file name.
     *
     * The optional parameter @p block_size sets the number of VcfRecord%s that are read
     * asynchronously into a buffer for speed improvements. This is mostly interesting for
     * window- or region-based analyses, where a certain number of records are needed to fill the
     * window, on which afterwards some (potentially time-consuming) operations and computations
     * are performed. In that time, an asynchronous thread can already read the next block of
     * VCF records. In these scenarios, it is best to chose a @p block_size that is larger than
     * the typical number of records per window/region that is being processed. E.g., if most
     * windows contain between 1200 and 1500 VcfRecord%s, a good @p block_size is 3000 or 5000,
     * so that subsequent windows can be filled quickly without having to wait for the reading.
     */
    explicit VcfInputIterator(
        std::string const& filename,
        size_t block_size = 1024
    )
        // Call the other constuctor, to avoid code duplication.
        : VcfInputIterator( filename, std::vector<std::string>{}, false, block_size )
    {}

    /**
     * @brief Create an instance that reads from an input file name.
     *
     * Additionally, this constructor takes a list of @p sample_names which are used as filter so
     * that only those samples (columns of the VCF records) are evaluated and accessible - or,
     * if @p inverse_sample_names is set to `true`, instead all <i>but</i> those samples.
     *
     * @copydetails VcfInputIterator( std::string const&, size_t )
     */
    VcfInputIterator(
        std::string const& filename,
        std::vector<std::string> const& sample_names,
        bool inverse_sample_names = false,
        size_t block_size = 1024
    )
        : filename_( filename )
        , block_size_( block_size )
        , file_(          std::make_shared<HtsFile>( filename ))
        , header_(        std::make_shared<VcfHeader>( *file_ ))
        , current_block_( std::make_shared<std::vector<VcfRecord>>() )
        , buffer_block_(  std::make_shared<std::vector<VcfRecord>>() )
        , thread_pool_(   std::make_shared<utils::ThreadPool>( 1 ))
        , future_(        std::make_shared<std::future<size_t>>() )
    {
        if( block_size_ == 0 ) {
            throw std::invalid_argument( "Invalid block_size == 0 for VcfInputIterator" );
        }

        // Above, we initialized part of the htslib-related data (file_, header_) to point to their
        // objects, as well as empty vectors for current_block_ and end_pos_, which will be filled in the
        // init_() call below. Furthermore, we started a thread_pool_ with exactly 1 thread,
        // and prepared the future that stores how many records were read in the prefetching.

        // Filter sample columns by their name.
        if( ! sample_names.empty() ) {
            header_->set_samples( sample_names, inverse_sample_names );
        }

        // Initialize the current_block_ and buffer_block_, and read the first block(s) of the file.
        init_();
    }

    ~VcfInputIterator() = default;

    VcfInputIterator( self_type const& ) = default;
    VcfInputIterator( self_type&& )      = default;

    self_type& operator= ( self_type const& ) = default;
    self_type& operator= ( self_type&& )      = default;

    // -------------------------------------------------------------------------
    //     Comparators
    // -------------------------------------------------------------------------

    /**
    * @brief Return true iff dereferencing is valid, i.e., iff there is a VCF record available.
    */
    explicit operator bool() const
    {
        assert( current_pos_ <= end_pos_ );
        return current_pos_ < end_pos_;
    }

    bool good() const
    {
        assert( current_pos_ <= end_pos_ );
        return current_pos_ < end_pos_;
    }

    // -------------------------------------------------------------------------
    //     Accessors
    // -------------------------------------------------------------------------

    std::string const& filename() const
    {
        return filename_;
    }

    HtsFile const& hts_file() const
    {
        // Here and below we assert the existence of a pointed-to object in the shared pointer,
        // which does not hold true if the default constructor of this class was used, in which
        // case any of these dereferencing functions here are supposed to be invalid - so, by
        // using assertions here, we can fail a bit more gracefully in such cases, instead of
        // getting seg faults due to nullptrs.
        assert( file_ );
        return *file_;
    }

    HtsFile& hts_file()
    {
        assert( file_ );
        return *file_;
    }

    VcfHeader const& header() const
    {
        assert( header_ );
        return *header_;
    }

    VcfHeader& header()
    {
        assert( header_ );
        return *header_;
    }

    VcfRecord const& record() const
    {
        assert( current_block_ );
        assert( current_pos_ < end_pos_ );
        return (*current_block_)[current_pos_];
    }

    VcfRecord& record()
    {
        assert( current_block_ );
        assert( current_pos_ < end_pos_ );
        return (*current_block_)[current_pos_];
    }

    VcfRecord const* operator ->() const
    {
        assert( current_block_ );
        assert( current_pos_ < end_pos_ );
        return &(*current_block_)[current_pos_];
    }

    VcfRecord* operator ->()
    {
        assert( current_block_ );
        assert( current_pos_ < end_pos_ );
        return &(*current_block_)[current_pos_];
    }

    VcfRecord const& operator*() const
    {
        assert( current_block_ );
        assert( current_pos_ < end_pos_ );
        return (*current_block_)[current_pos_];
    }

    VcfRecord& operator*()
    {
        assert( current_block_ );
        assert( current_pos_ < end_pos_ );
        return (*current_block_)[current_pos_];
    }

    // -------------------------------------------------------------------------
    //     Iteration
    // -------------------------------------------------------------------------

    self_type& operator ++ ()
    {
        increment_();
        return *this;
    }

    self_type& operator ++(int)
    {
        increment_();
        return *this;
    }

    bool operator==( self_type const& it ) const
    {
        // We want equality between iterators that are copies of each other, and inequality
        // for non copies that were not default constructed. For the default constructed iterator,
        // which serves as the past-the-end marker, we need a special case.

        // Test if either of the two was default constructed. If so, we want a non-default
        // constructed iterator also compare equal to a default constructed one only if it is done
        // reading data, in which case good() == false (whose code we copy here)
        if( ! file_ || ! it.file_ ) {
            return ( current_pos_ < end_pos_ ) == ( it.current_pos_ < it.end_pos_ );
        }

        // In all other cases, we have two normal iterators that we want to compare.
        // We only need to compare one of the pointers to make sure that the iterators point to the
        // same hts data (that is, one was created as a copy of the other), and assert the others.
        // We assert the rest (if file_ equals, then the others have to equal as well).
        assert(
            file_ != it.file_ ||
            (
                header_        == it.header_ &&
                current_block_ == it.current_block_ &&
                buffer_block_  == it.buffer_block_ &&
                thread_pool_   == it.thread_pool_ &&
                future_        == it.future_
            )
        );
        return file_ == it.file_;
    }

    bool operator!=( self_type const& it ) const
    {
        return !(*this == it);
    }

    // -------------------------------------------------------------------------
    //     Private Members
    // -------------------------------------------------------------------------

private:

    void init_()
    {
        assert( file_ );
        assert( header_ );
        assert( current_block_ );
        assert( buffer_block_ );

        // Init the records and create empty VcfRecord to read into.
        current_block_->reserve( block_size_ );
        buffer_block_->reserve( block_size_ );
        for( size_t i = 0; i < block_size_; ++i ) {
            current_block_->emplace_back( *header_ );
            buffer_block_->emplace_back( *header_ );
        }
        assert( current_block_->size() == block_size_ );
        assert( buffer_block_->size() == block_size_ );

        // Read the first block synchronously, so that there is data initialized to be dereferenced.
        end_pos_ = VcfInputIterator::read_block_( file_, current_block_, block_size_ );
        assert( current_pos_ == 0 );

        // If there is less data than the block size, the file is already done.
        // No need to start the async buffering.
        if( end_pos_ < block_size_ ) {
            return;
        }

        // Now start the worker thread to fill the buffer.
        fill_buffer_block_();
    }

    void increment_()
    {
        // Finish the reading (potentially waiting if not yet finished in the worker thread).
        // The future returns how much data there was to be read, which we use as our status.
        // After that, swap the buffer and start a new reading operation in the worker thread.

        // Move to the next element in the vector. If we are at the end of the record vector,
        // and if that vector was full (complete block size), there is more data, so start reading.
        ++current_pos_;
        if( current_pos_ == end_pos_ && end_pos_ == block_size_ ) {
            assert( future_ );
            assert( future_->valid() );

            // Get how many records were read into the buffer, which also waits for the reading
            // if necessary. After that, we can swao the buffer, start reading again, and set
            // or internal current location to the first element of the vector again.
            end_pos_ = future_->get();
            std::swap( buffer_block_, current_block_ );
            fill_buffer_block_();
            current_pos_ = 0;
        }
    }

    void fill_buffer_block_()
    {;
        assert( thread_pool_ );
        assert( future_ );

        // This function is only every called after we finished any previous operations,
        // so let's assert that the thread pool and future are in the states that we expect.
        assert( thread_pool_->load() == 0 );
        assert( ! future_->valid() );

        // In order to use lambda captures by copy for class member variables in C++11, we first
        // have to make local copies, and then capture those. Capturing the class members direclty
        // was only introduced later. Bit cumbersome, but gets the work done.
        auto file         = file_;
        auto buffer_block = buffer_block_;
        auto block_size   = block_size_;

        // The lambda returns the result of read_block_ call, that is, the number of records
        // that have been read, and which we later (in the future_) use to see how much data we got.
        *future_ = thread_pool_->enqueue(
            [ file, buffer_block, block_size ](){
                return VcfInputIterator::read_block_( file, buffer_block, block_size );
            }
        );
    }

    static size_t read_block_(
        std::shared_ptr<HtsFile> file,
        std::shared_ptr<std::vector<VcfRecord>> target,
        size_t block_size
    ) {
        // This is a static function that does not depend on the class member data, so that
        // we can use it from the future lambda in the thread pool above without having to worry
        // about lambda captures of `this` going extinct... which was an absolutely nasty bug to
        // find! Such a rookie mistake! For that reason, we also take all arguments as shared
        // pointers, so that they are kept alive while the thread pool is working.
        // However, once its done with its work, the function (the one that we give to the thread
        // pool with a lambda) is popped from the thread queue, so that the shared pointer can
        // be freed again - that is, we do not need to worry about the lambda keeping the shared
        // pointer from freeing its memory indefinitely.

        assert( file );
        assert( target->size() == block_size );

        // Read as long as there is data. Return number of read records.
        size_t i = 0;
        for( ; i < block_size; ++i ) {
            if( ! (*target)[i].read_next( *file ) ) {
                break;
            }
        }
        return i;
    }

    // -------------------------------------------------------------------------
    //     Data Members
    // -------------------------------------------------------------------------

private:

    std::string filename_;

    // We buffer block_size_ many vcf records, and within each block, iterate via current_pos_
    // from 0 (first element of the block) to end_pos_ (past the end counter).
    size_t block_size_ = 1024;
    size_t current_pos_ = 0;
    size_t end_pos_ = 0;

    // htslib structs. We use shared pointers here to allow copying this iterator.
    // Also, use a buffer into which we read asynchronously, and then swap with the current_block_
    std::shared_ptr<HtsFile>   file_;
    std::shared_ptr<VcfHeader> header_;
    std::shared_ptr<std::vector<VcfRecord>> current_block_;
    std::shared_ptr<std::vector<VcfRecord>> buffer_block_;

    // Thread pool to run the buffering in the background.
    // Also, store the future_ used to keep track of the background task. It returns the number of
    // record lines that have been read into the buffer (block_size_, or less at the end of the file).
    std::shared_ptr<utils::ThreadPool> thread_pool_;
    std::shared_ptr<std::future<size_t>> future_;

};

} // namespace population
} // namespace genesis

#endif // htslib guard
#endif // include guard
