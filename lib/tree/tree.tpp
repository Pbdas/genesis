/**
 * @brief Implementation of basic tree functions.
 *
 * For reasons of readability, in this implementation file, the template data types
 * NodeDataType and EdgeDataType are abbreviated using NDT and EDT, respectively.
 *
 * @file
 * @ingroup tree
 */

#include <assert.h>
#include <sstream>

#include "tree/distances.hpp"
#include "utils/logging.hpp"
#include "utils/utils.hpp"

namespace genesis {

// =============================================================================
//     Construction and Destruction
// =============================================================================

/**
 * @brief Copy constructor.
 *
 * This function creates all links, nodes and edges new, and shapes them so that the final
 * Tree has the same topology as the input Tree.
 *
 * The data belonging to the edges and nodes is copied using the copy assignment of the respective
 * template parameter classes EdgeDataType and NodeDataType. As this data might contain pointers and
 * other structures that need a deep copy, it is the responsibility of those data classes to make
 * sure its own data is copied correctly.
 *
 * In case of pointers, the class that is using the tree might need to do further copying after
 * calling this copy constructor. See Placements for an example.
 *
 * TODO Idea for a nice feature (not yet implemented):
 * The advantage of copying the topology only is that we are able to make this function completely
 * independend of the data, hence the `other` Tree does not need to share the same data types.
 * Some potential function declarations can be found in the header file tree.hpp.
 */
template <class NDT, class EDT>
Tree<NDT, EDT>::Tree (const Tree<NDT, EDT>& other)
{
    // preparation.
    clear();
    links_.resize(other.links_.size());
    nodes_.resize(other.nodes_.size());
    edges_.resize(other.edges_.size());

    // create all objects. we need two loops per array, because the pointers have to exist
    // in order to be used with each other.
    for (size_t i = 0; i < links_.size(); ++i) {
        links_[i] = make_unique<LinkType>();
        // links_[i]->data = other.links_[i]->data;
    }
    for (size_t i = 0; i < nodes_.size(); ++i) {
        nodes_[i] = make_unique<NodeType>();
        nodes_[i]->data = other.nodes_[i]->data;
    }
    for (size_t i = 0; i < edges_.size(); ++i) {
        edges_[i] = make_unique<EdgeType>();
        edges_[i]->data = other.edges_[i]->data;
    }

    // set all pointers for the topology in a second round of loops.
    for (size_t i = 0; i < links_.size(); ++i) {
        const auto& olink = other.links_[i];
        assert(olink->index_ == i);

        links_[i]->index_  = i;
        links_[i]->next_   = links_[olink->next_->index_].get();
        links_[i]->outer_  = links_[olink->outer_->index_].get();
        links_[i]->node_   = nodes_[olink->node_->index_].get();
        links_[i]->edge_   = edges_[olink->edge_->index_].get();
    }
    for (size_t i = 0; i < nodes_.size(); ++i) {
        const auto& onode = other.nodes_[i];
        assert(onode->index_ == i);

        nodes_[i]->index_  = i;
        nodes_[i]->link_   = links_[onode->link_->index_].get();
    }
    for (size_t i = 0; i < edges_.size(); ++i) {
        const auto& oedge = other.edges_[i];
        assert(oedge->index_ == i);

        edges_[i]->index_  = i;
        edges_[i]->link_p_ = links_[oedge->link_p_->index_].get();
        edges_[i]->link_s_ = links_[oedge->link_s_->index_].get();
    }
}

/**
 * @brief Assignment operator. Copies the topology, but not the data of a tree.
 *
 * See Tree copy constructor for more information.
 */
template <class NDT, class EDT>
Tree<NDT, EDT>& Tree<NDT, EDT>::operator = (const Tree<NDT, EDT>& other)
{
    // check for self-assignment. we want this explicitly, in order to avoid unnecessary copies of
    // the tree, which would mean loosing the data in process.
    if (&other == this) {
        return *this;
    }

    // the Tree tmp is a copy of the right hand side object (automatically created using the
    // copy constructor). we can thus simply swap the arrays, and upon leaving the function,
    // tmp is automatically destroyed, so that its arrays are cleared and the data freed.
    Tree<NDT, EDT> tmp(other);
    std::swap(links_, tmp.links_);
    std::swap(nodes_, tmp.nodes_);
    std::swap(edges_, tmp.edges_);
    return *this;
}

/**
 * @brief Destructor. Calls clear() to free all memory used by the tree and its substructures.
 */
template <class NDT, class EDT>
Tree<NDT, EDT>::~Tree ()
{
    // TODO totally not necessary to have that here! remove, run mem test!
    clear();
}

/**
 * @brief Deletes all data of the tree, including all links, nodes and edges.
 */
template <class NDT, class EDT>
void Tree<NDT, EDT>::clear()
{
    std::vector<std::unique_ptr<LinkType>>().swap(links_);
    std::vector<std::unique_ptr<NodeType>>().swap(nodes_);
    std::vector<std::unique_ptr<EdgeType>>().swap(edges_);
}

/**
 * @brief Swap.
 */
template <class NDT, class EDT>
void Tree<NDT, EDT>::swap (Tree<NDT, EDT>& other)
{
    std::swap(links_, other.links_);
    std::swap(nodes_, other.nodes_);
    std::swap(edges_, other.edges_);
}

/**
 * @brief Imports all elements of a tree.
 *
 * This function overwrites the topology and data of this tree with a given set of links, nodes
 * and edges. Use with care! No checks are done concerning the validity of the passed input.
 *
 * Caveat: Only the pointers to the tree elements are copied, not the elements themselves. Thus,
 * this function is not intended for creating a deep copy. It merely is a fast way to pass pointers
 * to tree elements.
 *
 * Therefore, the main usage of this function is to get a tree from different TreeProcessor objects
 * for reading trees from files.
 */
template <class NDT, class EDT>
void Tree<NDT, EDT>::import_content (LinkArray& links, NodeArray& nodes, EdgeArray& edges)
{
    using std::swap;

    clear();
    swap(links_, links);
    swap(nodes_, nodes);
    swap(edges_, edges);
}

/**
 * @brief Exports all elements of a tree.
 *
 * Caveat: Only the pointers to the tree elements are copied, not the elements themselves. Thus,
 * this function is not intended for creating a deep copy. It merely is a fast way to pass pointers
 * to tree elements.
 */
template <class NDT, class EDT>
void Tree<NDT, EDT>::export_content (LinkArray& links, NodeArray& nodes, EdgeArray& edges)
{
    links = links_;
    nodes = nodes_;
    edges = edges_;
}

// =============================================================================
//     Member Functions
// =============================================================================

/**
 * @brief Find a Node, given its name.
 *
 * TODO assumes that tree node has a name. not good.
 */
template <class NDT, class EDT>
typename Tree<NDT, EDT>::NodeType* Tree<NDT, EDT>::find_node(std::string name) const
{
    // TODO check first whether replacing underscores is necessary!
    name = string_replace_all(name, "_", " ");
    for (const auto& n : nodes_) {
        if (n->data.name == name) {
            return n.get();
        }
    }
    return nullptr;
}

/**
 * @brief Returns the highest rank of the nodes of the Tree.
 */
template <class NDT, class EDT>
int Tree<NDT, EDT>::max_rank() const
{
    int max = -1;
    for (size_t i = 0; i < nodes_.size(); ++i) {
        int rank = nodes_[i]->rank();
        if (rank == 1) {
            LOG_WARN << "Node with rank 1 found. This is a node without furcation, and usually "
                     << "indicates an error.";
        }
        max = std::max(rank, max);
    }
    return max;
}

/**
 * @brief Returns whether the Tree is bifurcating.
 */
template <class NDT, class EDT>
bool Tree<NDT, EDT>::is_bifurcating() const
{
    return max_rank() == 2;
}

/**
 * @brief Count the number of leaf nodes.
 */
template <class NDT, class EDT>
size_t Tree<NDT, EDT>::leaf_count() const
{
    size_t sum = 0;
    for (const auto& n : nodes_) {
        if (n->is_leaf()) {
            ++sum;
        }
    }
    return sum;
}

/**
 * @brief Count the number of inner nodes.
 */
template <class NDT, class EDT>
size_t Tree<NDT, EDT>::inner_count() const
{
    return nodes_.size() - leaf_count();
}

// =============================================================================
//     Comparison
// =============================================================================

/**
 * @brief Compares two trees for equality given a binary comparator functional.
 *
 * This function does a preorder traversal of both trees in parallel and calls a comparator
 * functional for each position of the iterator. It returns true iff the comparator is true for
 * every position.
 *
 * The comparator functional can be either a function pointer, function object, or lambda
 * expression.
 *
 * Furthermore, the trees are checked for equal topology: their elements (links, nodes, edges)
 * have to be equal in size and the rank of each node during the traversal has to be identical in
 * both trees. Those assumptions are made because two trees that do not have identical topology
 * are never considered equal for the purposes of this framework.
 */
template <class NDT, class EDT>
bool Tree<NDT, EDT>::equal(
    const TreeType& lhs,
    const TreeType& rhs,
    const std::function<bool (TreeType::ConstIteratorPreorder&, TreeType::ConstIteratorPreorder&)>
        comparator
) {
    // check array sizes
    if (lhs.links_.size() != rhs.links_.size() ||
        lhs.nodes_.size() != rhs.nodes_.size() ||
        lhs.edges_.size() != rhs.edges_.size()
    ) {
        return false;
    }

    // do a preorder traversal on both trees in parallel
    TreeType::ConstIteratorPreorder it_l = lhs.begin_preorder();
    TreeType::ConstIteratorPreorder it_r = rhs.begin_preorder();
    for (
        ;
        it_l != lhs.end_preorder() && it_r != rhs.end_preorder();
        ++it_l, ++it_r
    ) {
        if (it_l.node()->rank() != it_r.node()->rank() || !comparator(it_l, it_r)) {
            return false;
        }
    }

    // check whether we are done with both trees
    if (it_l != lhs.end_preorder() || it_r != rhs.end_preorder()) {
        return false;
    }

    return true;
}

/**
 * @brief Compares the tree to another one given a binary comparator functional.
 *
 * See static equal() function for more information.
 */
template <class NDT, class EDT>
bool Tree<NDT, EDT>::equal (
    const TreeType& other,
    const std::function<bool (TreeType::ConstIteratorPreorder&, TreeType::ConstIteratorPreorder&)>
        comparator
) const {
    return equal(*this, other, comparator);
}

/**
 * @brief Compares two trees for equality using the respective comparision operators for their nodes
 * and edges.
 *
 * This method is mainly a shortcut for the static equal function, where the comparator functional
 * is instanciated using the default comparision operators of the tree's data.
 */
template <class NDT, class EDT>
bool Tree<NDT, EDT>::equal (const TreeType& lhs, const TreeType& rhs)
{
    auto comparator = [] (
        TreeType::ConstIteratorPreorder& it_l,
        TreeType::ConstIteratorPreorder& it_r
    ) {
        return it_l.node() == it_r.node() && it_l.edge() == it_r.edge();
    };

    return equal(lhs, rhs, comparator);
}

/**
 * @brief Compares the tree to another one using the respective comparision operators for their nodes
 * and edges.
 *
 * See static equal() function for more information.
 */
template <class NDT, class EDT>
bool Tree<NDT, EDT>::equal(const TreeType& other) const
{
    return equal(*this, other);
}

/**
 * @brief Returns true iff both trees have an identical topology.
 *
 * The topology is considered identical only if the order of edges is also the same in both trees.
 * This means, although two trees might have the same number of leaves and branches, they might
 * still be not identical (with respect to this function) when the branches appear in a different
 * order or when the root sits at a different node.
 *
 * Thus, this function is mainly intended to check whether two trees have been produced from the
 * same input, for example from the same Newick file.
 */
template <class NDT, class EDT>
bool Tree<NDT, EDT>::identical_topology(const TreeType& lhs, const TreeType& rhs)
{
    auto comparator = [] (TreeType::ConstIteratorPreorder&, TreeType::ConstIteratorPreorder&)
    {
        return true;
    };

    return equal(lhs, rhs, comparator);
}

/**
 * @brief Returns true iff both trees have an identical topology.
 *
 * See static identical_topology() method for more information.
 */
template <class NDT, class EDT>
bool Tree<NDT, EDT>::identical_topology(const TreeType& other) const
{
    return identical_topology(*this, other);
}

// =============================================================================
//     Dump and Debug Functions
// =============================================================================

// TODO do all node->link_ links point to the root? same for all edge->primary?

/**
 * @brief Validate that all pointers of the tree elements (links, nodes, edges) to each other
 * are correct and that some other invariants are met.
 *
 * This check is a bit pedantic, but better safe than sorry.
 */
template <class NDT, class EDT>
bool Tree<NDT, EDT>::validate() const
{
    // check that the member arrays are valid: if at least one of them is empty, the tree is not
    // fully initialized, so either it is a new tree without any data (all arrays empty, valid),
    // or some are empty, but others not (not valid).
    if (links_.empty() || nodes_.empty() || edges_.empty()) {
        bool emp = links_.empty() && nodes_.empty() && edges_.empty();
        if (emp) {
            LOG_INFO << "Tree is empty.";
        } else {
            LOG_INFO << "Tree is not empty, but one of its data members is.";
        }
        return emp;
    }

    if (links_.front()->node_ != nodes_.front().get()) {
        LOG_INFO << "The first link does not correspond to the first node.";
        return false;
    }

    if (links_.front()->index() != 0 || links_.front()->node()->index() != 0) {
        LOG_INFO << "Root does not have index 0.";
        return false;
    }

    // Check Links.
    std::vector<size_t> links_to_edges(edges_.size(), 0);
    std::vector<size_t> links_to_nodes(nodes_.size(), 0);
    for (size_t i = 0; i < links_.size(); ++i) {
        // Check indices.
        if (i != links_[i]->index_) {
            LOG_INFO << "Link at index " << i << " has wrong index (" << links_[i]->index_ << ").";
            return false;
        }

        // Check next cycle and node.
        auto nl = links_[i].get();
        do {
            if (nl->node() != links_[i]->node()) {
                LOG_INFO << "Link at index " << nl->index_ << " points to wrong node.";
                return false;
            }
            nl = nl->next();
        } while(nl != links_[i].get());
        ++links_to_nodes[links_[i]->node()->index()];

        // Check outer cycle.
        if (links_[i]->outer()->outer() != links_[i].get()) {
            LOG_INFO << "Link at index " << i << " has wrong outer link.";
            return false;
        }

        // Check edge.
        auto edge = links_[i]->edge();
        if (edge->primary_link() != links_[i].get() && edge->secondary_link() != links_[i].get()) {
            LOG_INFO << "Link at index " << i << " has wrong edge pointer.";
            return false;
        }
        ++links_to_edges[links_[i]->edge()->index()];
    }

    // Check if all edges have been hit twice.
    for (size_t i = 0; i < edges_.size(); ++i) {
        if (links_to_edges[i] != 2) {
            LOG_INFO << "Edge at index " << i << " is not visited twice but " << links_to_edges[i]
                     << " times when traversing the links.";
            return false;
        }
    }

    // Check if all nodes have been hit as many times as their rank is.
    for (size_t i = 0; i < nodes_.size(); ++i) {
        if (links_to_nodes[i] != nodes_[i]->rank() + 1) {
            LOG_INFO << "Node at index " << i << " is not visited its rank + 1 ("
                     << nodes_[i]->rank() << " + 1 = " << nodes_[i]->rank() + 1
                     << ") times, but " << links_to_nodes[i] << " times when "
                     << "traversing the links.";
            return false;
        }
    }

    // Check Nodes.
    for (size_t i = 0; i < nodes_.size(); ++i) {
        // Check indices.
        if (i != nodes_[i]->index_) {
            LOG_INFO << "Node at index " << i << " has wrong index (" << nodes_[i]->index_ << ").";
            return false;
        }

        // Check link.
        if (nodes_[i]->link()->node() != nodes_[i].get()) {
            LOG_INFO << "Node at index " << i << " has wrong link.";
            return false;
        }
    }

    // Check Edges.
    for (size_t i = 0; i < edges_.size(); ++i) {
        // Check indices.
        if (i != edges_[i]->index_) {
            LOG_INFO << "Edge at index " << i << " has wrong index (" << edges_[i]->index_ << ").";
            return false;
        }

        // Check links.
        if (edges_[i]->primary_link()->edge() != edges_[i].get()) {
            LOG_INFO << "Edge at index " << i << " has wrong primary link.";
            return false;
        }
        if (edges_[i]->secondary_link()->edge() != edges_[i].get()) {
            LOG_INFO << "Edge at index " << i << " has wrong secondary link.";
            return false;
        }
    }

    // If we are here, all three arrays (links, nodes, edges) contain data, so we can start a full
    // traversal along all links.

    // Count, how many times the elements are hit while traversing.
    std::vector<size_t> it_links(links_.size(), 0);
    std::vector<size_t> it_edges(edges_.size(), 0);
    std::vector<size_t> it_nodes(nodes_.size(), 0);

    // Do the traversal. We do not use the iterator here, to keep it simple when testing.
    // (We want to validate the tree here, not the iterator.)
    LinkType* link = links_.front().get();
    do {
        ++it_links[link->index()];
        ++it_edges[link->edge()->index()];
        ++it_nodes[link->node()->index()];
        link = link->next_->outer_;
    } while (link != links_.front().get());

    // Check if all links have been hit once.
    for (size_t i = 0; i < links_.size(); i++) {
        if (it_links[i] != 1) {
            LOG_INFO << "Link at index " << i << " is not visited 1 but " << it_links[i]
                     << " times when iterating the tree.";
            return false;
        }
    }

    // Check if all edges have been hit twice.
    for (size_t i = 0; i < edges_.size(); ++i) {
        if (it_edges[i] != 2) {
            LOG_INFO << "Edge at index " << i << " is not visited 2 but " << it_edges[i]
                     << " times when iterating the tree.";
            return false;
        }
    }

    // Check if all nodes have been hit as many times as their rank is.
    for (size_t i = 0; i < nodes_.size(); ++i) {
        if (it_nodes[i] != nodes_[i]->rank() + 1) {
            LOG_INFO << "Node at index " << i << " is not visited "
                     << nodes_[i]->rank() << " + 1 = " << (nodes_[i]->rank() + 1)
                     << " times, but " << it_nodes[i] << " times when iterating the "
                     << "tree.";
            return false;
        }
    }

    return true;
}

/**
 * @brief Returns a simple text representation of the Tree, showing all nodes, edges and links
 * with their indices.
 */
template <class NDT, class EDT>
std::string Tree<NDT, EDT>::dump() const
{
    std::vector<int>    depth = node_depth_vector(*this);
    std::vector<size_t> done;
    std::ostringstream  out;

    // prepare link so that we point to the root link. this will ensure that the order of nodes
    // displayed by this funtion is the one expected by the user. usually, we would go into
    // the first branch immediately, but then there would be no way of first nicely displaying
    // the information about the root node. so we need to do it a bit more complex than the
    // usual iteration...
    LinkType* l = root_link();
    while (l->next() != root_link()) {
        l = l->next();
    }

    // do an euler tour traversal over all links. (we cannot use the iterator here, as
    // we need each link on its own, and not each node as the iterator gives)
    do {
        NodeType* n = l->node();
        std::string indent = std::string(4 * depth[n->index()], ' ');
        if (!contains(done, n->index())) {
            out << indent << "\033[1;31mNode " << n->index() << ": \"" << n->data.name << "\"\033[0m\n";
        }
        done.push_back(n->index());

        // dont display the next link when we are at the first iteration.
        if (l->next() == root_link()) {
            l = l->next();
        } else {
            out << indent;
            out << "    \033[34mLink " << l->index() << "\033[0m";
            l = l->next();
            out << " \033[32m>\033[0m \033[34mLink " << l->index() << "\033[0m\n";
        }

        out << indent;
        out << " -- \033[34mLink " << l->index() << "\033[0m";
        out << " -- \033[36mEdge " << l->edge()->index() << "\033[0m";
        l = l->outer();
        out << " --> \033[34mLink " << l->index() << "\033[0m\n";
    } while (l->next() != root_link());

    // output the last next link back to the root, because we skipped this in the loop
    // (the one that was skipped in the beginning).
    out << "    \033[34mLink " << l->index() << "\033[0m";
    l = l->next();
    out << " \033[32m>\033[0m \033[34mLink " << l->index() << "\033[0m\n";

    return out.str();
}

/**
 * @brief Returns lists of all nodes, edges and links including their indices and connections
 * with each other.
 */
template <class NDT, class EDT>
std::string Tree<NDT, EDT>::dump_lists() const
{
    std::ostringstream out;

    // nodes
    for (size_t i = 0; i < nodes_.size(); ++i) {
        out << "Node " << i
            << " \t Main Link: " << nodes_[i]->link_->index_
            << " \t " << nodes_[i]->dump() << "\n";
    }
    out << "\n";

    // edges
    for (size_t i = 0; i < edges_.size(); ++i) {
        out << "Edge " << i
            << " \t Link P: " << edges_[i]->link_p_->index_
            << " \t Link S: " << edges_[i]->link_s_->index_
            << " \t " << edges_[i]->dump() << "\n";
    }
    out << "\n";

    // links
    for (size_t i = 0; i < links_.size(); ++i) {
        out << "Link " << i
            << "  \t Next: "  << links_[i]->next_->index_
            << " \t Outer: " << links_[i]->outer_->index_
            << " \t Node: "  << links_[i]->node_->index_
            << " \t Edge: "  << links_[i]->edge_->index_
            << " \t " << links_[i]->dump()
            << "\n";
    }

    return out.str();
}

} // namespace genesis
