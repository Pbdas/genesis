#ifndef GENESIS_SEQUENCE_FUNCTION_DISTANCES_H_
#define GENESIS_SEQUENCE_FUNCTION_DISTANCES_H_

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
 * @brief A collection of distance function bewteen Sequences
 *
 * @file
 * @ingroup sequence
 */
#include "genesis/sequence/sequence.hpp"

#include <vector>
#include <stdexcept>

namespace genesis {

// =================================================================================================
//     Forward Declarations
// =================================================================================================

namespace sequence {

// =================================================================================================
//     Sequence to Sequence Distances
// =================================================================================================

/**
 * @brief Calculate the hamming distance between two Sequences of equal length
 */
double hamming_distance( Sequence const& lhs, Sequence const& rhs )
{
    // if ( not ( lhs.size() == rhs.size() ) ) {
    //     throw std::runtime_error{
    //         "Hamming distance can only be computed for sequences of equal length."
    //     };
    // }

    size_t distance = 0;
    for (size_t i = 0; i < lhs.size(); ++i) {
        distance += (lhs[i] != rhs[i]);
    }

    return distance;
}

} // namespace sequence
} // namespace genesis

#endif // include guard
