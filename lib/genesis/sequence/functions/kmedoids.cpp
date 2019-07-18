/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2019 Lucas Czech and Pierre Barbera

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
    Lucas Czech <lucas.czech@h-its.org>
    Exelixis Lab, Heidelberg Institute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

/**
 * @brief
 *
 * @file
 * @ingroup sequence
 */

#include "genesis/sequence/functions/kmedoids.hpp"

#include "genesis/sequence/functions/distances.hpp"
#include "genesis/sequence/functions/functions.hpp"
#include "genesis/utils/math/common.hpp"
#include "genesis/utils/core/algorithm.hpp"

#include <cassert>
#include <stdexcept>
#include <vector>
#include <unordered_set>

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

namespace genesis {
namespace sequence {

// =================================================================================================
//     Settings
// =================================================================================================


// =================================================================================================
//    Sequence Kmeans
// =================================================================================================

/*
    Overload of the initialize function that additionally precomputes the pairwise distance matrix
 */
std::vector<SequenceKmedoids::BasePoint> SequenceKmedoids::initialize( SequenceSet const& data, size_t const k )
{
    size_t const N = data.size();

    medoid_set_.reserve( k );

    // initialize the PWD
    pairwise_distances_ = utils::Matrix<double>(N, N, 0.0);

    #pragma omp parallel for
    for ( size_t i = 0; i < N; ++i ) {
        for ( size_t j = i+1; j < N; ++j ) {
            auto dist = hamming_distance( data[i], data[j] );
            pairwise_distances_( i, j ) = dist;
            pairwise_distances_( j, i ) = dist;
        }
    }

    // create the "fake" data that the underlying kmeans impl. will work with
    // i.e. the indices to the pairwise distance matrix
    std::vector<BasePoint> fake_data( N );
    #pragma omp parallel for
    for ( size_t i = 0; i < N; ++i ) {
        fake_data[ i ] = i;
    }

    return fake_data;
}

bool SequenceKmedoids::data_validation( SequenceSet const& data ) const
{
    return is_alignment( data );
}

size_t SequenceKmedoids::run( SequenceSet const& data, size_t const k )
{
    // Validate the data. The function might also throw on its own, in order
    // to provide a more helpful message about what is actually invalid about the data.
    if( ! data_validation( data ) ) {
        throw std::runtime_error( "Invalid data: Sequences need to be of equal length." );
    }

    // run the udnerlying kmeans
    auto iterations = Kmeans::run( initialize( data, k ), k );

    // set the result medoids
    for( auto const& c : centroids() ) {
        medoids_.add( data[ c ] );
    }

    return iterations;
}

/*
 * Evaluates the cost of a given configuration of datapoints and their assignments to medoids
 */
double SequenceKmedoids::cost(
    std::vector<BasePoint>  const& data,
    std::vector<size_t> const& assignments,
    std::vector<BasePoint>  const& medoids
) const {
    assert( data.size() == assignments.size() );

    double cost = 0;
    for ( size_t i = 0; i < data.size(); ++i ) {
        cost += distance( data[i], medoids[ assignments[ i ] ] );
    }

    return cost;
}

void SequenceKmedoids::update_centroids(
    std::vector<BasePoint>  const& data,
    std::vector<size_t> const& assignments,
    std::vector<BasePoint>&        medoids
) {
    // This function is only called from within the run() function, which already
    // checks this condition. So, simply assert it here, instead of throwing.
    assert( data.size() == assignments.size() );

    auto const k = medoids.size();

    // get the cost of the current configuration
    double previous_cost = cost( data, assignments, medoids );

    // make a local copy of the assignments
    std::vector<size_t> proposal_assignments = assignments;

    // get an ID shorthand for the current medoids
    if( medoid_set_.size() != medoids.size() ) {
        medoid_set_.clear();
        for ( auto const& medoid : medoids ) {
            medoid_set_.insert( medoid );
        }
    }

    // for each current medoid
    for( size_t cur_medoid_id = 0; cur_medoid_id < k; ++cur_medoid_id ) {
        auto& medoid = medoids[ cur_medoid_id ];
        auto original_medoid = medoid;
        bool better_found = false;
        // find a new medoid that would improve the overall clustering
        for( size_t cur_non_medoid = 0; cur_non_medoid < data.size(); ++cur_non_medoid ) {
            // only try non-medoids
            // auto cur_non_medoid = data[ i ];
            if ( utils::contains( medoid_set_, cur_non_medoid ) ) {
                continue;
            }

            // propose a new configuration
            medoid = data[ cur_non_medoid ];
            assign_to_centroids( data, medoids, proposal_assignments );
            auto new_cost = cost( data, proposal_assignments, medoids );
            if( new_cost < previous_cost ) {
                // better config found for this medoid, continue with the next
                previous_cost = new_cost;
                medoid_set_.erase( cur_medoid_id );
                medoid_set_.insert( cur_non_medoid );
                better_found = true;
                break;
            }
        }

        if( not better_found ) {
            // revert medoid back
            medoid = original_medoid;
        }
    }
}

double SequenceKmedoids::distance( BasePoint const& lhs, BasePoint const& rhs ) const
{
    return pairwise_distances_( lhs, rhs );
}

} // namespace sequence
} // namespace genesis
