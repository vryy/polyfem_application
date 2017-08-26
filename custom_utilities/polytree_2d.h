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
#include "utilities/indexed_object.h"
#include "containers/pointer_vector_set.h"
#include "includes/model_part.h"
#include "poly_half_edge.h"


namespace Kratos
{

/// Short class definition.
/** Polytree implementation in 2D
*/
class PolyTree2D
{
public:

    enum HalfEdgeStates
    {
        UNINITIALIZED = 0,
        INITIALIZED = 1,
        CACHED = 2,
        FINALIZED = 3
    };

    /// Pointer definition of PolyTree2D
    KRATOS_CLASS_POINTER_DEFINITION(PolyTree2D);

    typedef PolyVertex<2> VertexType;
    typedef PolyHalfEdge<2> EdgeType;
    typedef PolyFace<2> FaceType;

    // typedef std::map<std::size_t, VertexType::Pointer> VertexContainerType;
    typedef PointerVectorSet<VertexType, IndexedObject> VertexContainerType;

    // typedef std::list<FaceType::Pointer> FaceContainerType;
    typedef PointerVectorSet<FaceType, IndexedObject> FaceContainerType;

    /// Default constructor.
    PolyTree2D() {}

    /// Destructor.
    virtual ~PolyTree2D() {}

    /**
     * clear the internal data
     */
    void Clear()
    {
        // clear all the containers
        mVertexList.clear();
        mFaceList.clear();

        m_half_edge_state = UNINITIALIZED;
    }

    /**
     * Synchronize the half-edge data structure with ModelPart
     * @param r_model_part the input ModelPart
     */
    void Synchronize(ModelPart& r_model_part, bool forward)
    {
        if (forward == true)
        {
            // initialize half-edges structure
            std::cout << "Constructing half-edge data structure from ModelPart " << r_model_part.Name() << std::endl;
            InitializeHalfEdges(r_model_part, mVertexList, mFaceList);
            std::cout << "Constructing half-edge data structure completed " << std::endl;
            std::cout << "Number of vertices: " << mVertexList.size() << std::endl;
            std::cout << "Number of faces: " << mFaceList.size() << std::endl;

            m_half_edge_state = INITIALIZED;
        }
        else
        {
            if (m_half_edge_state == FINALIZED)
            {
                std::cout << "Synchronizing half-edge data structure to ModelPart " << r_model_part.Name() << std::endl;
                
                // first phase: synchronizing forward, check for each element in the ModelPart if it is changed, then add new nodes and replace by new element accordingly
                std::vector<std::size_t> ChangedElementIds;
                ModelPart::ElementsContainerType NewChangedElements;
                for (ModelPart::ElementsContainerType::ptr_iterator it = r_model_part.Elements().ptr_begin();
                        it != r_model_part.Elements().ptr_end(); ++it)
                {
                    FaceContainerType::iterator it_face = mFaceList.find((*it)->Id());
                    if (it_face != mFaceList.end())
                    {
                        if (it_face->IsChanged())
                        {
                            // TODO create new nodes and element
                            ChangedElementIds.push_back((*it)->Id());
                        }
                    }
                    else
                    {
                        // if the face id is not in the face list. It's probably removed from coarsen operation. We mark to remove the respective element from the ModelPart.
                    }
                }

                // second phase: synchronizing backward, check for new faces and add to the ModelPart, also creating new nodes
                for (FaceContainerType::iterator it = mFaceList.begin(); it != mFaceList.end(); ++it)
                {
                    ModelPart::ElementsContainerType::iterator it_elem = r_model_part.Elements().find(it->Id());
                    if (it_elem == r_model_part.Elements().end())
                    {
                        // TODO create new nodes and element
                    }
                    else
                    {
                        // the respective face is already handled in the first phase so we SHALL do nothing here
                    }
                }
            }
        }

    }

    /**
     * Mark a face to refine
     */
    void MarkFaceRefine(const std::size_t& face_index)
    {
        FaceContainerType::iterator it = mFaceList.find(face_index);
        if (it != mFaceList.end())
        {
            it->SetRefine(true);
            it->SetCoarsen(false);
        }
        m_half_edge_state = CACHED;
    }

    /**
     * Mark a face to coarsen
     */
    void MarkFaceCoarsen(const std::size_t& face_index)
    {
        FaceContainerType::iterator it = mFaceList.find(face_index);
        if (it != mFaceList.end())
        {
            it->SetRefine(false);
            it->SetCoarsen(true);
        }
        m_half_edge_state = CACHED;
    }

    /**
     * Cache refine operations
     */
    void BeginRefineCoarsen()
    {
        if (m_half_edge_state == CACHED)
        {
            for (FaceContainerType::iterator it = mFaceList.begin(); it != mFaceList.end(); ++it)
            {
                if (it->IsRefined())
                {
                    if (it->IsLeaf())
                    {
                        // TODO refine this face
                    }
                    else
                    {
                        // Refine all the sub-faces
                    }

                    // turn off the refine flag
                    it->SetRefine(false);
                }

                if (it->IsCoarsen())
                {
                    if (!it->IsLeaf())
                    {
                        // TODO remove all the sub-faces
                    }
                    else
                    {
                        // shall do nothing because this is a leaf face
                    }

                    // turn off the coarsen flag
                    it->SetCoarsen(false);
                }
            }
        }
    }

    /**
     * Finalize refine operations
     */
    void EndRefineCoarsen()
    {
        // merge the composite edges. This is to reduce number of closed points in the polygon mesh.
        for (FaceContainerType::iterator it = mFaceList.begin(); it != mFaceList.end(); ++it)
        {
            MergeCompositeEdges(*it);
        }

        // align the list of vertices from the list of faces
        mVertexList.clear();
        for (FaceContainerType::iterator it = mFaceList.begin(); it != mFaceList.end(); ++it)
        {
            it->AddVertices(mVertexList);
        }
        mVertexList.Unique();

        m_half_edge_state = FINALIZED;
    }

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
    }

private:

    /// Internal variables
    VertexContainerType mVertexList;
    FaceContainerType mFaceList;

    HalfEdgeStates m_half_edge_state;

    /**
     * Initialize the half edge data structure from ModelPart
     * @param r_model_part the input ModelPart
     * @param rVertexList  the output list of vertices
     * @param rEdgeList    the output list of edges
     * @param rFaceList    the output list of faces
     */
    void InitializeHalfEdges(ModelPart& r_model_part,
            VertexContainerType& rVertexList, FaceContainerType& rFaceList) const
    {
        typedef ModelPart::ElementType::GeometryType GeometryType;

        // create half-edge vertices
        for (ModelPart::NodeIterator it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd() ; ++it)
        {
            VertexType::Pointer pNewVertex = boost::make_shared<VertexType>(it->Id(), it->X(), it->Y());
            rVertexList.insert(rVertexList.begin(), pNewVertex);
        }
        rVertexList.Unique();
        std::cout << "Create vertices complete, " << rVertexList.size() << " vertices are created" << std::endl;

        // create half-edge edges and faces
        ModelPart::ElementsContainerType& pElements = r_model_part.Elements();

        typedef std::map<std::pair<std::size_t, std::size_t>, std::vector<EdgeType::Pointer> > EdgeMapType;
        EdgeMapType HalfEdgesMap;
        std::size_t NumberOfHalfEdges = 0;

        for (ModelPart::ElementsContainerType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
            GeometryType& r_geom = (*it)->GetGeometry();

            std::vector<EdgeType::Pointer> pLocalEdges(r_geom.size());

            for (std::size_t i = 0; i < r_geom.size(); ++i)
            {
                std::size_t this_node = i;
                std::size_t next_node = (i < r_geom.size()-1) ? i+1 : 0;
                // std::cout << "at node " << r_geom[this_node].Id() << std::endl;

                // create new half-edge
                VertexType::Pointer pVertex1 = rVertexList(r_geom[this_node].Id());
                VertexType::Pointer pVertex2 = rVertexList(r_geom[next_node].Id());
                EdgeType::Pointer pNewHalfEdge = boost::make_shared<EdgeType>(pVertex1, pVertex2);
                // std::cout << "half-edge btw vertex " << pVertex1->Id() << " and " << pVertex2->Id() << " is created" << std::endl;

                // add to edge map
                std::pair<std::size_t, std::size_t> key;
                key.first = std::min(r_geom[this_node].Id(), r_geom[next_node].Id());
                key.second = std::max(r_geom[this_node].Id(), r_geom[next_node].Id());
                HalfEdgesMap[key].push_back(pNewHalfEdge);

                // add to the list of local edges
                pLocalEdges[i] = pNewHalfEdge;

                // set the edge for the vertex
                pVertex1->pSetEdge(pNewHalfEdge);
                // std::cout << "half-edge for vertex " << pVertex1->Id() << " is set" << std::endl;
                // KRATOS_WATCH(*pVertex1->pEdge())

                // insert to list of half-edges
                ++NumberOfHalfEdges;
            }

            // create new face and insert to list
            FaceType::Pointer pNewFace = boost::make_shared<FaceType>((*it)->Id());
            pNewFace->pSetEdge(pLocalEdges[0]);
            rFaceList.insert(rFaceList.end(), pNewFace);
            // std::cout << "face for element " << (*it)->Id() << " is created" << std::endl;

            // KRATOS_WATCH(pLocalEdges.size())
            for (std::size_t i = 0; i < pLocalEdges.size(); ++i)
            {
                std::size_t prev_edge = (i != 0) ? i - 1 : pLocalEdges.size()-1;
                std::size_t next_edge = (i != pLocalEdges.size()-1) ? i + 1 : 0;
                // std::cout << "at local edge " << i << ", prev = " << prev_edge << ", next = " << next_edge << std::endl;
 
                pLocalEdges[i]->pSetPrevEdge(pLocalEdges[prev_edge]);
                pLocalEdges[i]->pSetNextEdge(pLocalEdges[next_edge]);
                // std::cout << "at local edge " << i << ", set prev/next is complete" << std::endl;

                pLocalEdges[i]->pSetFace(pNewFace);
                // std::cout << "at local edge " << i << ", set face is complete" << std::endl;
            }
            // std::cout << "prev/next for half-edge at face" << pNewFace->Id() << " is created" << std::endl;
        }
        std::cout << "Create half-edges complete, " << NumberOfHalfEdges << " half-edges are created" << std::endl;
        std::cout << "Create faces complete, " << rFaceList.size() << " faces are created" << std::endl;

        // assign opposite edge for each half-edge
        for (EdgeMapType::iterator it = HalfEdgesMap.begin(); it != HalfEdgesMap.end(); ++it)
        {
            if (it->second.size() == 2)
            {
                it->second[0]->pSetOppositeEdge(it->second[1]);
                it->second[1]->pSetOppositeEdge(it->second[0]);

                // verify if the half-edge opposite edge is valid
                if ( ( it->second[0]->OppositeEdge().Node1().Id() != it->second[1]->OppositeEdge().Node2().Id() )
                  || ( it->second[1]->OppositeEdge().Node1().Id() != it->second[0]->OppositeEdge().Node2().Id() ) )
                {
                    std::stringstream ss;
                    ss << "The half-edge opposite is invalid at " << *(it->second[0]) << " and "
                       << *(it->second[1]) << std::endl;
                    KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
                }
            }
            else if(it->second.size() > 2)
            {
                KRATOS_THROW_ERROR(std::logic_error, "Something wrong with the half edges", "")
            }
        }
    }

    /**
     * Refine an edge by adding a vertex in the middle
     * @param  pEdge the edge need to refine
     * @return       the newly created vertex
     */
    VertexType::Pointer RefineEdge(EdgeType::Pointer pEdge) const
    {
        // Create new vertex
        double x = 0.5 * (pEdge->Node1()[0] + pEdge->Node2()[0]);
        double y = 0.5 * (pEdge->Node1()[1] + pEdge->Node2()[1]);
        VertexType::Pointer pVertex = boost::make_shared<VertexType>(0, x, y);

        // Create new half-edge
        EdgeType::Pointer pNewEdge1 = boost::make_shared<EdgeType>(pEdge->pNode1(), pVertex);
        EdgeType::Pointer pNewEdge2 = boost::make_shared<EdgeType>(pVertex, pEdge->pNode2());
        EdgeType::Pointer pNewOppositeEdge1 = boost::make_shared<EdgeType>(pVertex, pEdge->pNode1());
        EdgeType::Pointer pNewOppositeEdge2 = boost::make_shared<EdgeType>(pEdge->pNode2(), pVertex);

        // assign edge for first side
        pEdge->PrevEdge().pSetNextEdge(pNewEdge1);

        pNewEdge1->pSetPrevEdge(pEdge->pPrevEdge());
        pNewEdge1->pSetNextEdge(pNewEdge2);
        pNewEdge1->pSetOppositeEdge(pNewOppositeEdge1);
        pNewEdge1->pSetFace(pEdge->pFace());

        pNewEdge2->pSetPrevEdge(pNewEdge1);
        pNewEdge2->pSetNextEdge(pEdge->pNextEdge());
        pNewEdge2->pSetOppositeEdge(pNewOppositeEdge2);
        pNewEdge2->pSetFace(pEdge->pFace());

        pEdge->NextEdge().pSetPrevEdge(pNewEdge2);

        // assign edge for the opposite side
        pEdge->OppositeEdge().NextEdge().pSetPrevEdge(pNewOppositeEdge1);

        pNewOppositeEdge1->pSetPrevEdge(pNewOppositeEdge2);
        pNewOppositeEdge1->pSetNextEdge(pEdge->OppositeEdge().pNextEdge());
        pNewOppositeEdge1->pSetOppositeEdge(pNewEdge1);
        pNewOppositeEdge1->pSetFace(pEdge->OppositeEdge().pFace());

        pNewOppositeEdge2->pSetPrevEdge(pEdge->OppositeEdge().pPrevEdge());
        pNewOppositeEdge2->pSetNextEdge(pNewOppositeEdge1);
        pNewOppositeEdge2->pSetOppositeEdge(pNewEdge2);
        pNewOppositeEdge2->pSetFace(pEdge->OppositeEdge().pFace());

        pEdge->OppositeEdge().PrevEdge().pSetNextEdge(pNewOppositeEdge2);

        return pVertex;
    }

    /**
     * Merge the closed points on the composite edge of the face
     * @param rFace the input face
     */
    void MergeCompositeEdges(FaceType& rFace) const
    {
        if (rFace.IsLeaf())
        {
            EdgeType::Pointer pFirstEdge = rFace.pEdge();
            EdgeType::Pointer pEdge = pFirstEdge;

            do
            {
                if (pEdge->IsComposite())
                {
                    // TODO merge composite edges

                }
                pEdge = pFirstEdge->pNextEdge();
            } while (pEdge != pEdge);
        }
        else
        {
            for (std::size_t i = 0; i < rFace.NumberOfSubFaces(); ++i)
            {
                MergeCompositeEdges(rFace.SubFace(i));
            }
        }
    }

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
    rOStream << " ";
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos.

#endif // KRATOS_POLY_TREE_2D_H_INCLUDED  defined

