#ifndef GENESIS_PLACEMENT_FUNCTIONS_DIVERSITY_H_
#define GENESIS_PLACEMENT_FUNCTIONS_DIVERSITY_H_

/*
    Genesis - A toolkit for working with phylogenetic data.
    Copyright (C) 2014-2022 Pierre Barbera, Lucas Czech

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
 * @ingroup placement
 */

namespace genesis {

// =================================================================================================
//     Forward Declarations
// =================================================================================================

namespace placement {

    class Sample;

}

namespace placement {

// =================================================================================================
//     Within-Sample Diversity (alpha)
// =================================================================================================

/**
 * @brief      Function calculating the Balance Weighted Phylogenetic Diversity
 * of a given sample (from https://dx.doi.org/10.7717%2Fpeerj.157).
 *
 * @param      sample  The sample
 * @param[in]  theta   Value specifying whether the PD should be calculated as
 * normal (0.0), fully balance-weighted (1.0), or anything in between.
 *
 * @return     The BWPD of a sample
 */
double BWPD( Sample const& sample, double const theta );

/**
 * @brief      Function calculating the Phylogenetic Diversity (PD) of a given sample.
 * 
 *             
 *
 * @param      sample  The sample
 *
 * @return     Faith's PD of the given sample
 */
double PD( Sample const& sample );


// =================================================================================================
//     Between-Sample Diversity (beta)
// =================================================================================================

} // namespace placement
} // namespace genesis

#endif // include guard
