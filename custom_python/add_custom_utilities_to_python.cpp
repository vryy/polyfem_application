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
#include "containers/data_value_container.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/polyfem_utility.h"
#include "custom_utilities/poly_half_edge.h"
#include "custom_utilities/polytree_2d.h"
#include "custom_utilities/polytree_utility.h"
#include "custom_utilities/polytree_sync_utility.h"
#include "python/pointer_vector_set_python_interface.h"


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

void PolyTreeSyncUtility_InitializeHalfEdges(PolyTreeSyncUtility& rDummy, ModelPart& r_model_part, PolyTree2D& r_tree)
{
    std::cout << "Constructing half-edge data structure from ModelPart " << r_model_part.Name() << std::endl;
    rDummy.InitializeHalfEdges(r_model_part, r_tree.Vertices(), r_tree.Edges(), r_tree.Faces(), r_tree.LastVertexId(), r_tree.LastFaceId());
    r_tree.SetInitialized();
    std::cout << "Constructing half-edge data structure completed " << std::endl;
    std::cout << "Number of vertices: " << r_tree.Vertices().size() << std::endl;
    std::cout << "Number of half-edges: " << r_tree.Edges().size() << std::endl;
    std::cout << "Number of faces: " << r_tree.Faces().size() << std::endl;
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
std::size_t PolyFace_Id(PolyFace<TDim>& dummy)
{
    return dummy.Id();
}

template<std::size_t TDim>
std::size_t PolyFace_HashCode(PolyFace<TDim>& dummy)
{
    return PolyHash<TDim>{}(dummy);
}

//////////////////////////////////////////////////

void PolyTree2D_WriteMatlabToFile(PolyTree2D& dummy, std::string filename,
        bool write_vertex_number, bool write_face_number)
{
    std::ofstream outfile;
    outfile.open(filename.c_str());
    dummy.WriteMatlab(outfile, write_vertex_number, write_face_number);
    outfile.close();
}

std::size_t PolyTree2D_LastVertexId(PolyTree2D& dummy)
{
    std::size_t LastVertexId = dummy.LastVertexId();
    return LastVertexId;
}

std::size_t PolyTree2D_LastFaceId(PolyTree2D& dummy)
{
    std::size_t LastFaceId = dummy.LastFaceId();
    return LastFaceId;
}

void PolyTree2D_ListVertex(PolyTree2D& dummy, const std::size_t Id)
{
    dummy.ListVertex(std::cout, Id);
}

void PolyTree2D_ListFace(PolyTree2D& dummy, const std::size_t Id)
{
    dummy.ListFace(std::cout, Id);
}

void PolyTree2D_ListVertices(PolyTree2D& dummy)
{
    dummy.ListVertices(std::cout);
}

void PolyTree2D_ListEdges(PolyTree2D& dummy)
{
    dummy.ListEdges(std::cout);
}

void PolyTree2D_ListFaces(PolyTree2D& dummy)
{
    dummy.ListFaces(std::cout);
}

PolyTree2D::VertexContainerType::Pointer PolyTree2D_GetVertices(PolyTree2D& dummy)
{
    return dummy.pVertices();
}

void PolyTree2D_SetVertices(PolyTree2D& dummy, PolyTree2D::VertexContainerType::Pointer pOtherVertices)
{
    dummy.SetVertices(pOtherVertices);
}

PolyTree2D::EdgeContainerType::Pointer PolyTree2D_GetEdges(PolyTree2D& dummy)
{
    return dummy.pEdges();
}

void PolyTree2D_SetEdges(PolyTree2D& dummy, PolyTree2D::EdgeContainerType::Pointer pOtherEdges)
{
    dummy.SetEdges(pOtherEdges);
}

PolyTree2D::FaceContainerType::Pointer PolyTree2D_GetFaces(PolyTree2D& dummy)
{
    return dummy.pFaces();
}

void PolyTree2D_SetFaces(PolyTree2D& dummy, PolyTree2D::FaceContainerType::Pointer pOtherFaces)
{
    dummy.SetFaces(pOtherFaces);
}

////////////////////////////////////////////////////

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
    ("PolyHalfEdge2D", init<PolyHalfEdge<2>::VertexType::Pointer, PolyHalfEdge<2>::VertexType::Pointer>())
    .def("HashCode", &PolyHalfEdge_HashCode<2>)
    .def(self_ns::str(self))
    ;

    class_<PolyFace<2>, PolyFace<2>::Pointer, boost::noncopyable>
    ("PolyFace2D", init<const std::size_t&>())
    .def("Id", &PolyFace_Id<2>)
    .def("HashCode", &PolyFace_HashCode<2>)
    .def(self_ns::str(self))
    ;

    class_<PolyTree2D, PolyTree2D::Pointer, bases<DataValueContainer>, boost::noncopyable>
    ("PolyTree2D", init<>())
    .def("LastVertexId", &PolyTree2D_LastVertexId)
    .def("LastFaceId", &PolyTree2D_LastFaceId)
    .add_property("Vertices", PolyTree2D_GetVertices, PolyTree2D_SetVertices)
    .add_property("Edges", PolyTree2D_GetEdges, PolyTree2D_SetEdges)
    .add_property("Faces", PolyTree2D_GetFaces, PolyTree2D_SetFaces)
    .def("VerticesArray", &PolyTree2D::Vertices, return_internal_reference<>())
    .def("EdgesArray", &PolyTree2D::Edges, return_internal_reference<>())
    .def("FacesArray", &PolyTree2D::Faces, return_internal_reference<>())
    .def("CreateFace", &PolyTree2D::CreateFace)
    .def("MarkFaceRefine", &PolyTree2D::MarkFaceRefine)
    .def("MarkFaceCoarsen", &PolyTree2D::MarkFaceCoarsen)
    .def("BeginRefineCoarsen", &PolyTree2D::BeginRefineCoarsen)
    .def("EndRefineCoarsen", &PolyTree2D::EndRefineCoarsen)
    .def("Validate", &PolyTree2D::Validate)
    .def("WriteMatlab", &PolyTree2D_WriteMatlabToFile)
    .def("ListVertex", &PolyTree2D_ListVertex)
    .def("ListFace", &PolyTree2D_ListFace)
    .def("ListVertices", &PolyTree2D_ListVertices)
    .def("ListEdges", &PolyTree2D_ListEdges)
    .def("ListFaces", &PolyTree2D_ListFaces)
    .def(self_ns::str(self))
    ;

    PointerVectorSetPythonInterface<PolyTree2D::VertexContainerType>::CreateInterface("PolyTree2DVertexArray");
    PointerVectorSetPythonInterface<PolyTree2D::EdgeContainerType>::CreateInterface("PolyTree2DEdgeArray");
    PointerVectorSetPythonInterface<PolyTree2D::FaceContainerType>::CreateInterface("PolyTree2DFaceArray");

    //////////////////////////////////////////////////////////////////////////

    class_<PolyTreeUtility, PolyTreeUtility::Pointer, boost::noncopyable>
    ("PolyTreeUtility", init<>())
    .def("TestComputeVoronoiTesselation", &PolyTreeUtility::TestComputeVoronoiTesselation)
    .def("TestComputeReflectionPoints", &PolyTreeUtility::TestComputeReflectionPoints)
    .def("TestPolygonDecomposition", &PolyTreeUtility::TestPolygonDecomposition)
    .def("TestClusterPoints1", &PolyTreeUtility::TestClusterPoints1)
    .def("TestClusterPoints2", &PolyTreeUtility::TestClusterPoints2)
    .def("TestClusterLengths1", &PolyTreeUtility::TestClusterLengths1)
    .def(self_ns::str(self))
    ;

    //////////////////////////////////////////////////////////////////////////

    class_<PolyTreeSyncUtility, PolyTreeSyncUtility::Pointer, boost::noncopyable>
    ("PolyTreeSyncUtility", init<>())
    .def("InitializeHalfEdges", &PolyTreeSyncUtility_InitializeHalfEdges)
    .def(self_ns::str(self))
    ;

}

}  // namespace Python.

}  // namespace Kratos.

