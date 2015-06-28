/**
 * @brief
 *
 * @file
 * @ingroup python
 */

#include <python/src/common.hpp>

#include "lib/tree/phyloxml_processor.hpp"

#include "lib/tree/default_tree.hpp"
#include "lib/placement/placement_tree.hpp"

// -------------------------------------------------------------------
//     Class Phyloxml Processor
// -------------------------------------------------------------------

template<class TreeType>
void Overload_PhyloxmlProcessor(boost::python::class_<genesis::PhyloxmlProcessor>& phyloxml_processor)
{
    phyloxml_processor
        .def(
            "to_file",
            ( bool ( ::genesis::PhyloxmlProcessor::* )( const std::string, const TreeType & ))( &::genesis::PhyloxmlProcessor::to_file ),
            ( boost::python::arg("fn"), boost::python::arg("tree") ),
            get_docstring("bool ::genesis::PhyloxmlProcessor::to_file (const std::string fn, const Tree< NDT, EDT > & tree)")
        )
        .def(
            "to_string",
            ( std::string ( ::genesis::PhyloxmlProcessor::* )( const TreeType & ))( &::genesis::PhyloxmlProcessor::to_string ),
            ( boost::python::arg("tree") ),
            get_docstring("std::string ::genesis::PhyloxmlProcessor::to_string (const Tree< NDT, EDT > & tree)")
        )
    ;
}

void BoostPythonExport_PhyloxmlProcessor()
{
    auto phyloxml_processor_ = boost::python::class_< ::genesis::PhyloxmlProcessor > ( "PhyloxmlProcessor" )

        // Public Member Variables

    ;

    /*
        TODO Phyloxml tree and other tree related classes are initilized using those hard coded lists
        of classes. instead, a more flexible approach could be to use variadic templates to achieve this:
        the tree export function takes the class types of all trees to be exported, defines the python
        bindings for them, and then delegates to all the other tree related exporters using the same
        template paramters.
    */
    Overload_PhyloxmlProcessor<genesis::PlacementTree>(phyloxml_processor_);
    Overload_PhyloxmlProcessor<genesis::DefaultTree>  (phyloxml_processor_);
}
