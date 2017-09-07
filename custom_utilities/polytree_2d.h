//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         polyfem_application/LICENSE.txt
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Hoang-Giang Bui
//  Date:            21 Aug 2017
//


#if !defined(KRATOS_POLY_TREE_H_INCLUDED )
#define  KRATOS_POLY_TREE_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <fstream>
#include <list>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "utilities/indexed_object.h"
#include "utilities/pair_indexed_object.h"
#include "containers/pointer_vector_set.h"
#include "containers/data_value_container.h"
#include "includes/model_part.h"
#include "poly_half_edge.h"


namespace Kratos
{

/// Short class definition.
/**
 * Polytree implementation in 2D
 * REF: H. Nguyen-Xuan et al, A polytree-based adaptive approach to limit analysis of cracked structures
 *
 * The polytree maintains a list of vertices and list of faces. It is kept to synchronize with the ModelPart.
 *
 * TODO:
 * 
 */
class PolyTree2D : public DataValueContainer
{
public:

    enum HalfEdgeStates
    {
        UNINITIALIZED = 0,
        INITIALIZED = 1,
        CACHED = 2,
        FINALIZED_READY = 3,
        FINALIZED = 4
    };

    /// Pointer definition of PolyTree2D
    KRATOS_CLASS_POINTER_DEFINITION(PolyTree2D);

    typedef PolyHalfEdge<2> EdgeType;
    typedef EdgeType::VertexType VertexType;
    typedef PolyFace<2> FaceType;

    // typedef std::map<std::size_t, VertexType::Pointer> VertexContainerType;
    typedef PointerVectorSet<VertexType, IndexedObject> VertexContainerType;

    typedef PairIndexedObject EdgeGetKeyType;
    typedef PointerVectorSet<EdgeType, EdgeGetKeyType, PairIndexedObjectCompare, PairIndexedObjectEqual> EdgeContainerType;

    // typedef std::list<FaceType::Pointer> FaceContainerType;
    typedef PointerVectorSet<FaceType, IndexedObject> FaceContainerType;

    /// Default constructor.
    PolyTree2D() : mLastVertexId(0), mLastFaceId(0) {}

    /// Destructor.
    virtual ~PolyTree2D()
    {
        std::cout << "PolyTree2D is destroyed" << std::endl;
    }

    /**
     * clear the internal data
     */
    void Clear();

    /**
     * @return the last vertex id of the polygon tree
     */
    std::size_t LastVertexId() const {return mLastVertexId;}

    /**
     * @return the last face id of the polygon tree
     */
    std::size_t LastFaceId() const {return mLastFaceId;}

    /**
     * Create a new face with Id // Do not use this in simulation, this is only for debugging
     * @param Id the id
     */
    FaceType::Pointer CreateFace(const std::size_t& Id);

    /**
     * Synchronize the half-edge data structure with ModelPart
     * @param r_model_part the input ModelPart
     */
    void Synchronize(ModelPart& r_model_part, bool forward);

    /**
     * Mark a face to refine
     * @param  face_index index of face to be refined
     * @return -1 if the face is inactive
     * @return -2 if the face does not exist
     * @return 0 is the face is marked successfully
     */
    int MarkFaceRefine(const std::size_t& face_index);

    /**
     * Mark a face to coarsen
     * @param face_index index of face to be coarsen
     * @return  -1 if the face is active
     * @return  -2 if the face does not exist
     * @return  0 if the face is marked successfully
     */
    int MarkFaceCoarsen(const std::size_t& face_index);

    /**
     * Cache refine operations
     */
    void BeginRefineCoarsen();

    /**
     * Finalize refine operations
     */
    void EndRefineCoarsen();

    /**
     * Perform various checking on the validity of the tree
     */
    void Validate() const;

    /**
     * Export the tree for visualization in Matlab
     * @param rOStream output stream
     * @param write_vertex_number flag to write the number to the vertex
     * @param write_face_number   flag to write the number to the face
     */
    void WriteMatlab(std::ostream& rOStream, const bool& write_vertex_number, const bool& write_face_number) const;

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "PolyTree2D";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << " Vertices:" << std::endl;
        for (VertexContainerType::const_iterator it = mVertexList.begin(); it != mVertexList.end(); ++it)
        {
            rOStream << "  " << *it << std::endl;
        }

        rOStream << " Edges:" << std::endl;
        for (EdgeContainerType::const_iterator it = mEdgeList.begin(); it != mEdgeList.end(); ++it)
        {
            rOStream << "  " << *it << std::endl;
        }

        rOStream << " Faces:" << std::endl;
        for (FaceContainerType::const_iterator it = mFaceList.begin(); it != mFaceList.end(); ++it)
        {
            rOStream << "  " << *it << std::endl;
        }
    }

private:

    /// Internal variables
    VertexContainerType mVertexList;
    EdgeContainerType mEdgeList;
    EdgeContainerType mNewEdgeList; // container to contain newly created half-edges, will be deleted after EndRefineCoarsen
    EdgeContainerType mCompositeEdgeList; // container to contain newly created composite half-edges for merge operation, will be deleted after EndRefineCoarsen
    FaceContainerType mFaceList;

    HalfEdgeStates m_half_edge_state;

    std::size_t mLastVertexId;
    std::size_t mLastFaceId;

    /**
     * Initialize the half edge data structure from ModelPart
     * @param r_model_part the input ModelPart
     * @param rVertexList  the output list of vertices
     * @param rEdgeList    the output list of edges
     * @param rFaceList    the output list of faces
     */
    void InitializeHalfEdges(ModelPart& r_model_part,
            VertexContainerType& rVertexList,
            EdgeContainerType& rEdgeList,
            FaceContainerType& rFaceList,
            std::size_t& LastVertexId, std::size_t& LastFaceId) const;

    /**
     * Refine an edge by adding a vertex in the middle
     * @param  pEdge the edge need to refine
     * @return       the newly created vertex
     */
    VertexType::Pointer RefineEdge(EdgeType::Pointer pEdge) const;

    /**
     * Merge the closed points on the composite edge of the face
     * @param rFace the input face
     * @param rVertexList the list of all vertices in the tree, new vertices will be added
     * @param rNewVertexList the list of newly created vertices
     * @param rEdgeList the list of all half-edges
     * @param rCompositeEdgeList the list of all composite edges
     * @param LastVertexId the last vertex id
     * @param alpha the merging parameter
     */
    void MergeEdges(EdgeType::Pointer pEdge,
            EdgeType::Pointer pOppositeEdge,
            VertexContainerType& rVertexList,
            VertexContainerType& rNewVertexList,
            EdgeContainerType& rEdgeList,
            EdgeContainerType& rCompositeEdgeList,
            std::size_t& LastVertexId,
            const double& alpha) const;

    /**
     * Refine the face and add the sub-faces to the face list
     * @param rFace        the face to be refined. It will then be removed from the face list
     * @param rVertexList  the vertex list. It will be added with the new vertices.
     * @param rEdgeList    the edge list. It will be added with the new half-edges.
     * @param rCompositeEdgeList    the composite edge list. It will be added with the new composite half-edges.
     * @param rFaceList    the face list. It will be added with the new faces.
     * @param LastVertexId the last vertex id. It will be increased accordingly during refinement.
     * @param LastFaceId   the last face id. It will be increased accordingly during refinement.
     */
    void RefineFace(FaceType::Pointer pFace,
            PolyTree2D::VertexContainerType& rVertexList,
            PolyTree2D::EdgeContainerType& rEdgeList,
            PolyTree2D::EdgeContainerType& rNewEdgeList,
            PolyTree2D::EdgeContainerType& rCompositeEdgeList,
            FaceContainerType& rFaceList,
            std::size_t& LastVertexId, std::size_t& LastFaceId) const;

    /**
     * Coarsen a face
     * @param rFace     the face need to be coarsen and removed
     * @param rFaceList the list of face
     */
    void CoarsenFace(FaceType& rFace, FaceContainerType& rFaceList) const;

    /**
     * Reconnect the edges in the case the point are merged
     * @param rEdgeList [description]
     */
    void ReconnectEdges(EdgeContainerType& rEdgeList) const;

    /**
     * Remove the zero edges and duplicated edges in the edge list
     * @param rEdgeList the processing edge list
     */
    void RemoveZeroAndDuplidateEdges(EdgeContainerType& rEdgeList) const;

    /**
     * Remove the lone edges in the edge list. The lone edge is the edge that point to an inactive face.
     * @param rEdgeList the processing edge list
     */
    void RemoveLoneEdges(EdgeContainerType& rEdgeList) const;

    /// Assignment operator.
    PolyTree2D& operator=(PolyTree2D const& rOther);

    /// Copy constructor.
    PolyTree2D(PolyTree2D const& rOther);
};

/// input stream function
inline std::istream& operator >> (std::istream& rIStream, PolyTree2D& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const PolyTree2D& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos.

#endif // KRATOS_POLY_TREE_2D_H_INCLUDED  defined

