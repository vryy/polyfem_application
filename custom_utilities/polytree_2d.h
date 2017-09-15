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
        CACHED_REFINE = 2,
        CACHED_COARSEN = 3,
        FINALIZED_READY = 4,
        FINALIZED = 5
    };

    /// Pointer definition of PolyTree2D
    KRATOS_CLASS_POINTER_DEFINITION(PolyTree2D);

    typedef PolyHalfEdge<2> EdgeType;
    typedef EdgeType::VertexType VertexType;
    typedef PolyFace<2> FaceType;

    typedef PointerVectorSet<VertexType, IndexedObject> VertexContainerType;

    typedef PairIndexedObject EdgeGetKeyType;
    typedef PointerVectorSet<EdgeType, EdgeGetKeyType, PairIndexedObjectCompare, PairIndexedObjectEqual> EdgeContainerType;

    typedef PointerVectorSet<FaceType, IndexedObject> FaceContainerType;

    /// Default constructor.
    PolyTree2D() : mLastVertexId(0), mLastFaceId(0)
    , mpVertexList(new VertexContainerType())
    , mpEdgeList(new EdgeContainerType())
    , mpFaceList(new FaceContainerType())
    {}

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
    std::size_t& LastVertexId() {return mLastVertexId;}
    const std::size_t& LastVertexId() const {return mLastVertexId;}

    /**
     * @return the last face id of the polygon tree
     */
    std::size_t& LastFaceId() {return mLastFaceId;}
    const std::size_t& LastFaceId() const {return mLastFaceId;}

    /// Get & Set for vertex container
    VertexContainerType& Vertices() {return *mpVertexList;}
    VertexContainerType::Pointer pVertices() {return mpVertexList;}
    void SetVertices(typename VertexContainerType::Pointer pOtherVertices) {mpVertexList = pOtherVertices;}

    /// Get & Set for edge container
    EdgeContainerType& Edges() {return *mpEdgeList;}
    EdgeContainerType::Pointer pEdges() {return mpEdgeList;}
    void SetEdges(typename EdgeContainerType::Pointer pOtherEdges) {mpEdgeList = pOtherEdges;}

    /// Get & Set for face container
    FaceContainerType& Faces() {return *mpFaceList;}
    FaceContainerType::Pointer pFaces() {return mpFaceList;}
    void SetFaces(typename FaceContainerType::Pointer pOtherFaces) {mpFaceList = pOtherFaces;}

    /**
     * Create a new face with Id // Do not use this in simulation, this is only for debugging
     * @param Id the id
     */
    FaceType::Pointer CreateFace(const std::size_t& Id);

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
     * Set the initialized state for the polytree
     */
    void SetInitialized() {m_half_edge_state = INITIALIZED;}

    /**
     * @return true if the tree is done with all refine/coarsen operations and ready to provide relevant information
     */
    bool IsFinalized() const {return m_half_edge_state == FINALIZED;}

    /**
     * Re-set the Id of all the vertices and the faces
     * @param rMapVertices the map from old vertex Id to new vertex Id
     * @param rMapFaces    the map from old face Id to new face Id
     */
    void Renumber(std::map<std::size_t, std::size_t>& rMapVertices,
            std::map<std::size_t, std::size_t>& rMapFaces);

    /**
     * Perform various checking on the validity of the tree
     */
    void Validate() const;

    /**
     * List the face
     * @param rOStream output stream
     * @param Id id of the face
     */
    void ListVertex(std::ostream& rOStream, const std::size_t& Id) const;

    /**
     * List the face
     * @param rOStream output stream
     * @param Id id of the face
     */
    void ListFace(std::ostream& rOStream, const std::size_t& Id) const;

    /**
     * List all vertices of the polytree
     * @param rOStream output stream
     */
    void ListVertices(std::ostream& rOStream) const;

    /**
     * List all edges of the polytree
     * @param rOStream output stream
     */
    void ListEdges(std::ostream& rOStream) const;

    /**
     * List all faces of the polytree
     * @param rOStream output stream
     */
    void ListFaces(std::ostream& rOStream) const;

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
        for (VertexContainerType::const_iterator it = mpVertexList->begin(); it != mpVertexList->end(); ++it)
        {
            rOStream << "  " << *it << std::endl;
        }

        rOStream << " Edges:" << std::endl;
        for (EdgeContainerType::const_iterator it = mpEdgeList->begin(); it != mpEdgeList->end(); ++it)
        {
            rOStream << "  " << *it << std::endl;
        }

        rOStream << " Faces:" << std::endl;
        for (FaceContainerType::const_iterator it = mpFaceList->begin(); it != mpFaceList->end(); ++it)
        {
            rOStream << "  " << *it << std::endl;
        }
    }

private:

    /// Internal variables
    VertexContainerType::Pointer mpVertexList;
    EdgeContainerType::Pointer mpEdgeList;
    EdgeContainerType mNewEdgeList; // container to contain newly created half-edges, will be deleted after EndRefineCoarsen
    EdgeContainerType mCompositeEdgeList; // container to contain newly created composite half-edges for merge operation, will be deleted after EndRefineCoarsen
    FaceContainerType::Pointer mpFaceList;

    HalfEdgeStates m_half_edge_state;

    std::size_t mLastVertexId;
    std::size_t mLastFaceId;

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

    /**
     * Remove the lone vertices in the vertex list. The lone vertex is the vertex that has no edge.
     * @param rVertexList the processing vertex list
     */
    void RemoveLoneVertices(VertexContainerType& rVertexList) const;

    /**
     * Extract the edges of the face and all the sub-faces
     * @param pFace     the input face
     * @param rFaceList the list of face
     */
    void ExtractEdges(FaceType::Pointer pFace, EdgeContainerType& rEdgeList) const;

    /**
     * Extract the face and all the sub-faces
     * @param pFace     the input face
     * @param rFaceList the list of face
     */
    void ExtractFaces(FaceType::Pointer pFace, FaceContainerType& rFaceList) const;

    /**
     * Extract the list of edges around a face. This function shall only be called after EndRefineCoarsen
     * @param pFace            the input face
     * @param rOuterEdgeList   the list of all outer edges. This is organized to make a closed loop.
     * @param rInnerEdgeList   the list of all inner edges (unorganized)
     * @param rInnerVertexList the list of all inner vertices (unorganized)
     * @param rInnerFaceList   the list of all inner faces (unorganized)
     */
    void ExtractEdgesAndVertices(FaceType::Pointer pFace,
            EdgeContainerType& rOuterEdgeList,
            EdgeContainerType& rInnerEdgeList,
            VertexContainerType& rInnerVertexList,
            FaceContainerType& rInnerFaceList) const;

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

