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

#include "genesis/placement/function/diversity.hpp"

#include "genesis/placement/sample.hpp"
#include "genesis/placement/function/functions.hpp"
#include "genesis/placement/function/helper.hpp"
#include "genesis/tree/function/functions.hpp"
#include "genesis/tree/iterator/postorder.hpp"

#include <cmath>
#include <vector>

#ifdef GENESIS_OPENMP
#   include <omp.h>
#endif

namespace genesis {
namespace placement {

// =================================================================================================
//     Local Helper Functions
// =================================================================================================

static bool equals_approx( double const a, double const b, double const epsilon = 1e-10 )
{
  return std::abs( a - b ) < epsilon;
}

// ==== Functions that can be used with BWPD_traversal ====

// These first two are currently unused, they implement the phylogenetic {quadratic} entropy.
// see the BWPD paper, where they use it.
static double __attribute__((unused)) phylo_entropy_g( double const x )
{
  if( equals_approx( x, 0.0 ) ) {
    return 0.0;
  }
  return x * std::log( x );
}

static double __attribute__((unused)) phylo_quad_entropy_g( double const x )
{
  return x * ( 1.0 - x );
}

static double step_function_g( double const x, double const theta )
{
  // special case: only consider edges in the spanning tree of the sample
  // (0 or 1 fraction means non-shared edge, which means the edge is not in the spanning)
  if( equals_approx( x, 0.0 ) or equals_approx( x, 1.0 ) ) {
    return 0.0;
  }

  return std::pow( 2 * std::min( x, 1.0 - x ), theta );
}

static size_t num_queries( std::vector< Pquery const* > const& pqs )
{
  size_t sum = 0;
  for( auto const& pq_ptr : pqs ) {
    sum += pq_ptr->name_size();
  }
  return sum;
}

static size_t num_queries( Sample const& sample )
{
  return total_name_count( sample );
}

/**
 * @brief      Traversal function to calculate the BWPD, taking a functional to
 * scale the individual values.
 *
 * @param      sample    The sample
 * @param[in]  g_func    The functional
 * @param[in]  args      The arguments to the functional
 *
 * @return     The BWPD of a sample
 */
template<class g_func_t, class... Args>
double
BWPD_traversal(Sample const& sample, g_func_t g_func, Args... args)
{
  /*  The goal is to iterate over all edges and calculate D_s(i), which is
     the fraction of reads in sample s that are in leaves (or edges /
     placements in our case) on the distal side of edge i. Then, the BWPD is the
     sum of edge length l(i) multiplied by [2min(D(i),1−D(i))]^θ. For θ=0, and
     evaluating only those edges that are in a samples spanning tree, this
     results in the somewhat classical FaithPD, except adapted to phylogenetic
     placement.
    */

  assert( theta >= 0.0 and theta <= 1.0 );

  // we start by first calculating D(i) for every edge in the tree
  // which we do bottom-up, post-order, dragging with us previous distal-side
  // counts, dividing them by the total count

  auto const& place_tree = sample.tree();
  auto pqrs_per_edge = pqueries_per_edge(sample, true);

  // if we only take the best hits, num_placements = num_pqueries, making the
  // total size, including multiplicities, simply the name count:
  double const total_queries = num_queries(sample);

  // vector to track the number of queries on the distal side of the edge
  std::vector<size_t> per_edge_distal_count(place_tree.edge_count());
  // same idea, but for D(i)
  std::vector<double> per_edge_D(place_tree.edge_count());

  double result = 0.0;

  // traverse the tree
  for (auto const& it : postorder(place_tree)) {
    // ensure last edge isn't visited twice
    if (it.is_last_iteration()) {
      continue;
    }

    auto const& edge = it.edge();
    auto const& edge_index = edge.index();

    // if this is a leaf, there cannot be any mass on the distal side, so we
    // skip
    if (is_leaf(edge)) {
      per_edge_distal_count[edge_index] = 0;
      per_edge_D[edge_index] = 0.0;
      continue;
    }
    // interior edges:
    // calculate the new distal count as the sum of distal counts and query
    // counts of the child edges
    auto const& node = it.node();
    auto const& lhs_edge = node.link().next().edge();
    auto const& rhs_edge = node.link().next().next().edge();
    auto const lhs_edge_index = lhs_edge.index();
    auto const rhs_edge_index = rhs_edge.index();

    per_edge_distal_count[edge_index] =
      per_edge_distal_count[lhs_edge_index] +
      num_queries(pqrs_per_edge[lhs_edge_index]) +
      per_edge_distal_count[rhs_edge_index] +
      num_queries(pqrs_per_edge[rhs_edge_index]);

    // calculate D(i) for this edge
    per_edge_D[edge_index] = per_edge_distal_count[edge_index] / total_queries;

    auto const branch_length = edge.data<PlacementEdgeData>().branch_length;

    // update the BWPD sum with the result of this edge
    result += branch_length * g_func(per_edge_D[edge_index], args...);
  }

  return result;
}

// =================================================================================================
//     Alpha Diversity
// =================================================================================================
double BWPD(
    Sample const&   sample,
    double const    theta
) {
    return BWPD_traversal( sample, step_function_g, theta );
}

double PD( Sample const& sample )
{
    // BWPD with theta = 0.0 is simply Faith's PD (but here implemented for placements)
    return BWPD( sample, 0.0 );
}

// =================================================================================================
//     Beta diversity
// =================================================================================================

} // namespace placement
} // namespace genesis
