// see finite_cell_application/LICENSE.txt
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: hbui $
//   Date:                $Date: 12 Aug 2017 $
//   Revision:            $Revision: 1.0 $
//
//



// Project includes
#include "includes/element.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/polyfem_utility.h"


namespace Kratos
{

namespace Python
{

using namespace boost::python;

std::size_t PolyFEMUtility_GetLastNodeId(PolyFEMUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastNodeId(r_model_part);
}

std::size_t PolyFEMUtility_GetLastElementId(PolyFEMUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastElementId(r_model_part);
}

std::size_t PolyFEMUtility_GetLastConditionId(PolyFEMUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastConditionId(r_model_part);
}

std::size_t PolyFEMUtility_GetLastPropertiesId(PolyFEMUtility& rDummy, ModelPart& r_model_part)
{
    return rDummy.GetLastPropertiesId(r_model_part);
}

void PolyFEMUtility_AddElement(PolyFEMUtility& rDummy, ModelPart::ElementsContainerType& rpElements,
        Element::Pointer pElement)
{
    rDummy.AddElement(rpElements, pElement);
}

void PolyFEMApplication_AddCustomUtilitiesToPython()
{
    Condition::Pointer(PolyFEMUtility::*pointer_to_PyCreateCondition)(ModelPart&, const std::string&,
            const std::size_t&, Properties::Pointer, boost::python::list&) const = &PolyFEMUtility::PyCreateCondition;

    ModelPart::ElementsContainerType(PolyFEMUtility::*pointer_to_PyGetElements)(ModelPart&,
            boost::python::list&) const = &PolyFEMUtility::PyGetElements;

    void(PolyFEMUtility::*pointer_to_PyGetElements2)(ModelPart::ElementsContainerType&, ModelPart&,
            boost::python::list&) const = &PolyFEMUtility::PyGetElements;

    void(PolyFEMUtility::*pointer_to_Clean)(ModelPart&,
            ModelPart::ConditionsContainerType&, const int&) const = &PolyFEMUtility::Clean;

    void(PolyFEMUtility::*pointer_to_PrintGeometry)(Element::GeometryType::Pointer) const = &PolyFEMUtility::Print;

    class_<PolyFEMUtility, PolyFEMUtility::Pointer, boost::noncopyable>
    ("PolyFEMUtility", init<>())
    .def("CreateCondition", pointer_to_PyCreateCondition)
    .def("GetElements", pointer_to_PyGetElements)
    .def("GetElements", pointer_to_PyGetElements2)
    .def("Clean", pointer_to_Clean)
    .def("GetLastNodeId", &PolyFEMUtility_GetLastNodeId)
    .def("GetLastElementId", &PolyFEMUtility_GetLastElementId)
    .def("GetLastConditionId", &PolyFEMUtility_GetLastConditionId)
    .def("GetLastPropertiesId", &PolyFEMUtility_GetLastPropertiesId)
    .def("AddElement", &PolyFEMUtility_AddElement)
    .def("Print", pointer_to_PrintGeometry)
    .def("TestPolygonShapeFunction", &PolyFEMUtility::TestPolygonShapeFunction)
    .def("TestPolygonQuadrature", &PolyFEMUtility::TestPolygonQuadrature)
    ;
}

}  // namespace Python.

}  // namespace Kratos.

