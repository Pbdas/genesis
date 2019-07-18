#ifndef GENESIS_SEQUENCE_FUNCTIONS_KMEANS_H_
#define GENESIS_SEQUENCE_FUNCTIONS_KMEANS_H_

/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2019 Lucas Czech, Pierre Barbera and HITS gGmbH

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
 * @brief A specialization of the generic Kmeans-class for the purpose of
 *        clustering `Sequence`s
 *
 * @file
 * @ingroup sequence
 */

#include "genesis/utils/math/kmeans.hpp"
#include "genesis/utils/containers/matrix.hpp"
#include "genesis/sequence/sequence.hpp"

#include <unordered_map>

namespace genesis {
namespace sequence {

// =================================================================================================
//     Forward Declarations
// =================================================================================================


// =================================================================================================
//     Mass Tree Kmeans
// =================================================================================================

class SequenceKmedoids
    : public utils::Kmeans< Sequence >
{
public:

    // -------------------------------------------------------------------------
    //     Typedefs and Constants
    // -------------------------------------------------------------------------

    using Point = Sequence;

    // -------------------------------------------------------------------------
    //     Constructors and Rule of Five
    // -------------------------------------------------------------------------

    SequenceKmedoids() = default;
    SequenceKmedoids( bool use_pwd ) : use_pwd_( use_pwd ) {};
    virtual ~SequenceKmedoids() override = default;

    SequenceKmedoids( SequenceKmedoids const& ) = default;
    SequenceKmedoids( SequenceKmedoids&& )      = default;

    SequenceKmedoids& operator= ( SequenceKmedoids const& ) = default;
    SequenceKmedoids& operator= ( SequenceKmedoids&& )      = default;

    // -------------------------------------------------------------------------
    //     Settings
    // -------------------------------------------------------------------------

    // -------------------------------------------------------------------------
    //     Default K-Means Functions
    // -------------------------------------------------------------------------

protected:
    virtual void initialize( std::vector<Point> const& data, size_t const k ) override;

private:

    virtual bool data_validation( std::vector<Point> const& data ) const override;

    virtual void update_centroids(
        std::vector<Point>  const& data,
        std::vector<size_t> const& assignments,
        std::vector<Point>&        centroids
    ) override;

    virtual double distance( Point const& lhs, Point const& rhs ) const override;

    double cost(
        std::vector<Point>  const& data,
        std::vector<size_t> const& assignments,
        std::vector<Point>  const& centroids
    ) const;

    // -------------------------------------------------------------------------
    //     Data
    // -------------------------------------------------------------------------

    // Pairwise distance matrix
    bool const use_pwd_ = true;
    utils::Matrix<double> pairwise_distances_;
    std::unordered_map<std::string, size_t> label_to_id_;
    std::unordered_set<size_t> medoid_set_;
};

} // namespace sequence
} // namespace genesis

#endif // include guard
