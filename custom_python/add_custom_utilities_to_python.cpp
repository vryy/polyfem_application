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
#include "custom_utilities/poly_half_edge.h"
#include "custom_utilities/polytree_2d.h"
#include "custom_utilities/polytree_utility.h"


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

//////////////////////////////////////

template<std::size_t TDim>
std::size_t PolyVertex_Id(PolyVertex<TDim>& dummy)
{
    return dummy.Id();
}

template<std::size_t TDim>
void PolyVertex_setitem(PolyVertex<TDim>& dummy, int index, double value)
{
    if (index >= 0 && index < dummy.size())
    {
        dummy[index] = value;
    }
    else
    {
        PyErr_SetString(PyExc_IndexError, "PolyVertex: index out of range");
        throw_error_already_set();
    }
}

template<std::size_t TDim>
double PolyVertex_getitem(PolyVertex<TDim>& dummy, int index)
{
    if (index >= 0 && index < dummy.size())
    {
        return dummy[index];
    }
    else
    {
        PyErr_SetString(PyExc_IndexError, "PolyVertex: index out of range");
        throw_error_already_set();
    }
}

template<std::size_t TDim>
std::size_t PolyHalfEdge_Id(PolyHalfEdge<TDim>& dummy)
{
    return dummy.Id();
}

template<std::size_t TDim>
std::size_t PolyHalfEdge_HashCode(PolyHalfEdge<TDim>& dummy)
{
    return PolyHash<TDim>{}(dummy);
}

template<std::size_t TDim>
std::size_t PolyFace_HashCode(PolyFace<TDim>& dummy)
{
    return PolyHash<TDim>{}(dummy);
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

    //////////////////////////////////////////////////////////////////////////

    class_<PolyVertex<2>, PolyVertex<2>::Pointer, boost::noncopyable>
    ("PolyVertex2D", init<const std::size_t&, const double&, const double&>())
    .def("Id", &PolyVertex_Id<2>)
    .def("__getitem__", &PolyVertex_getitem<2>)
     .def("__setitem__", &PolyVertex_setitem<2>)
    .def(self_ns::str(self))
    ;

    class_<PolyVertex<3>, PolyVertex<3>::Pointer, boost::noncopyable>
    ("PolyVertex3D", init<const std::size_t&, const double&, const double&, const double&>())
    .def("Id", &PolyVertex_Id<3>)
    .def("__getitem__", &PolyVertex_getitem<3>)
     .def("__setitem__", &PolyVertex_setitem<3>)
    .def(self_ns::str(self))
    ;

    class_<PolyHalfEdge<2>, PolyHalfEdge<2>::Pointer, boost::noncopyable>
    ("PolyHalfEdge2D", init<PolyHalfEdge<2>::NodeType::Pointer, PolyHalfEdge<2>::NodeType::Pointer>())
    .def("HashCode", &PolyHalfEdge_HashCode<2>)
    .def(self_ns::str(self))
    ;

    class_<PolyFace<2>, PolyFace<2>::Pointer, boost::noncopyable>
    ("PolyFace2D", init<const std::size_t&>())
    .def("HashCode", &PolyFace_HashCode<2>)
    .def(self_ns::str(self))
    ;

    class_<PolyTree2D, PolyTree2D::Pointer, boost::noncopyable>
    ("PolyTree2D", init<>())
    .def("Synchronize", &PolyTree2D::Synchronize)
    .def(self_ns::str(self))
    ;

    //////////////////////////////////////////////////////////////////////////

    class_<PolyTreeUtility, PolyTreeUtility::Pointer, boost::noncopyable>
    ("PolyTreeUtility", init<>())
    .def("TestComputeVoronoiTesselation", &PolyTreeUtility::TestComputeVoronoiTesselation)
    .def("TestComputeReflectionPoints", &PolyTreeUtility::TestComputeReflectionPoints)
    .def("TestPolygonDecomposition", &PolyTreeUtility::TestPolygonDecomposition)
    .def(self_ns::str(self))
    ;

}

}  // namespace Python.

}  // namespace Kratos.

