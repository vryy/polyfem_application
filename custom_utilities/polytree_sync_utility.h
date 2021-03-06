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
//  Date:            16 Sep 2017
//


#if !defined(KRATOS_POLYTREE_SYNC_UTILITY_H_INCLUDED )
#define  KRATOS_POLYTREE_SYNC_UTILITY_H_INCLUDED



// System includes
#include <iostream>


// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

#define DEBUG_SYNCHRONIZE

namespace Kratos
{
///@addtogroup PolyFEMApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** class for auxiliary routines
*/
class PolyTreeSyncUtility
{
public:
    ///@name Type Definitions
    ///@{

    enum SynchronizeStates
    {
        UNSYNCHRONIZED = 0,
        SYNCHRONIZED = 1
    };

    /// Pointer definition of PolyTreeSyncUtility
    KRATOS_CLASS_POINTER_DEFINITION(PolyTreeSyncUtility);

    typedef Element::GeometryType GeometryType;
    typedef GeometryType::CoordinatesArrayType CoordinatesArrayType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PolyTreeSyncUtility()
    : m_synchronize_state(UNSYNCHRONIZED)
    {}

    /// Destructor.
    virtual ~PolyTreeSyncUtility()
    {}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Initialize the half edge data structure for the polytree_2d from ModelPart
     * It is noted that the vertex coordinated is initialized to the original coordinates of the node. To synchronize the coordinates, use SetDeformed
     * @param r_model_part the input ModelPart
     * @param rVertexList  the output list of vertices
     * @param rEdgeList    the output list of edges
     * @param rFaceList    the output list of faces
     */
    template<class TVertexContainerType, class TEdgeContainerType, class TFaceContainerType>
    void InitializeHalfEdges(ModelPart& r_model_part,
            TVertexContainerType& rVertexList,
            TEdgeContainerType& rEdgeList,
            TFaceContainerType& rFaceList,
            std::size_t& LastVertexId, std::size_t& LastFaceId)
    {
        typedef ModelPart::ElementType::GeometryType GeometryType;
        typedef typename TVertexContainerType::data_type VertexType;
        typedef typename TEdgeContainerType::data_type EdgeType;
        typedef typename TFaceContainerType::data_type FaceType;

        // create half-edge vertices
        LastVertexId = 0;
        for (ModelPart::NodeIterator it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd() ; ++it)
        {
            typename VertexType::Pointer pNewVertex = boost::make_shared<VertexType>(it->Id(), it->X0(), it->Y0());
            rVertexList.insert(rVertexList.begin(), pNewVertex);

            if (it->Id() > LastVertexId)
                LastVertexId = it->Id();
        }
        rVertexList.Unique();
        std::cout << "Create vertices complete, " << rVertexList.size() << " vertices are created" << std::endl;

        // create half-edge edges and faces
        ModelPart::ElementsContainerType& pElements = r_model_part.Elements();

        typedef std::map<std::pair<std::size_t, std::size_t>, std::vector<typename EdgeType::Pointer> > EdgeMapType;
        EdgeMapType HalfEdgesMap;
        std::size_t NumberOfHalfEdges = 0;

        LastFaceId = 0;
        for (ModelPart::ElementsContainerType::ptr_iterator it = pElements.ptr_begin(); it != pElements.ptr_end(); ++it)
        {
            GeometryType& r_geom = (*it)->GetGeometry();

            std::vector<typename EdgeType::Pointer> pLocalEdges(r_geom.size());

            for (std::size_t i = 0; i < r_geom.size(); ++i)
            {
                std::size_t this_node = i;
                std::size_t next_node = (i < r_geom.size()-1) ? i+1 : 0;
                // std::cout << "at node " << r_geom[this_node].Id() << std::endl;

                // create new half-edge
                typename VertexType::Pointer pVertex1 = rVertexList(r_geom[this_node].Id());
                typename VertexType::Pointer pVertex2 = rVertexList(r_geom[next_node].Id());
                typename EdgeType::Pointer pNewHalfEdge = boost::make_shared<EdgeType>(pVertex1, pVertex2);
                // std::cout << "half-edge btw vertex " << pVertex1->Id() << " and " << pVertex2->Id() << " is created" << std::endl;

                // add to edge map
                std::pair<std::size_t, std::size_t> key;
                key.first = std::min(r_geom[this_node].Id(), r_geom[next_node].Id());
                key.second = std::max(r_geom[this_node].Id(), r_geom[next_node].Id());
                HalfEdgesMap[key].push_back(pNewHalfEdge);

                // add to global container
                rEdgeList.insert(rEdgeList.begin(), pNewHalfEdge);

                // add to the list of local edges
                pLocalEdges[i] = pNewHalfEdge;

                // set the edge for the vertex
                pVertex1->pSetEdge(pNewHalfEdge);
        //         // std::cout << "half-edge for vertex " << pVertex1->Id() << " is set" << std::endl;
        //         // KRATOS_WATCH(*pVertex1->pEdge())

                // insert to list of half-edges
                ++NumberOfHalfEdges;
            }

            // unique operation (removing duplicated half-edges), here it should remove nothing
            rEdgeList.Unique();

            // create new face and insert to list
            typename FaceType::Pointer pNewFace = boost::make_shared<FaceType>((*it)->Id());
            pNewFace->pSetEdge(pLocalEdges[0]);
            pNewFace->SetState(FaceType::NORMAL);
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

            if ((*it)->Id() > LastFaceId)
                LastFaceId = (*it)->Id();
        }
        std::cout << "Create half-edges complete, " << NumberOfHalfEdges << " half-edges are created" << std::endl;
        std::cout << "Create faces complete, " << rFaceList.size() << " faces are created" << std::endl;

        // assign opposite edge for each half-edge
        for (typename EdgeMapType::iterator it = HalfEdgesMap.begin(); it != HalfEdgesMap.end(); ++it)
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

        // after the first initialization, the state is assumed to be SYNCHRONIZED
        m_synchronize_state = SYNCHRONIZED;
    }


    /**
     * Synchronize the vertex coordinates to the current coordinates of the nodes in model_part. It is useful for ALE simulation.
     */
    template<class TTreeType>
    void SetDeformed(ModelPart& r_model_part, TTreeType& r_tree) const
    {
        typedef typename TTreeType::VertexContainerType VertexContainerType;
        typedef typename TTreeType::EdgeContainerType EdgeContainerType;
        typedef typename TTreeType::FaceContainerType FaceContainerType;

        typedef typename TTreeType::VertexType VertexType;
        typedef typename TTreeType::EdgeType EdgeType;
        typedef typename TTreeType::FaceType FaceType;

        for (typename VertexContainerType::iterator it = r_tree.Vertices().begin(); it != r_tree.Vertices().end(); ++it)
        {
            ModelPart::NodesContainerType::const_iterator it_node = r_model_part.Nodes().find(it->Id());

            if (it_node == r_model_part.Nodes().end())
            {
                std::stringstream ss;
                ss << "Node " << it->Id() << " does not exist in the model_part";
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
            }
            else
            {
                it->SetCoordinate(0, it_node->X());
                it->SetCoordinate(1, it_node->Y());
            }
        }
    }


    /**
     * Begin synchronization the nodes & elements in the model_part and the vertices and faces of the polytree after coarsen/refinement
     * In the beginning phase, the new nodes and elements are created, but the old ones won't be deleted. They must exist to transfer the physical data.
     */
    template<class TTreeType>
    void BeginSynchronize(ModelPart& r_model_part, const TTreeType& r_tree)
    {
        std::cout << "<<<<<<<<<<<<<<<<" << __FUNCTION__ << " begin" << ">>>>>>>>>>>>>>>>" << std::endl;

        if (m_synchronize_state == UNSYNCHRONIZED)
            return;

        if (m_synchronize_state == SYNCHRONIZED)
            m_synchronize_state = UNSYNCHRONIZED;

        typedef typename TTreeType::VertexContainerType VertexContainerType;
        typedef typename TTreeType::EdgeContainerType EdgeContainerType;
        typedef typename TTreeType::FaceContainerType FaceContainerType;

        typedef typename TTreeType::VertexType VertexType;
        typedef typename TTreeType::EdgeType EdgeType;
        typedef typename TTreeType::FaceType FaceType;

        if (!r_tree.IsReadyForSynchronization())
        {
            KRATOS_THROW_ERROR(std::logic_error, __LINE__, ": The polytree is not ready for synchronization");
            return;
        }

        // collect all the new vertices of the tree that does not exist in the model_part. In addition, find out what element the new node belongs to, so that the correct interpolation factor is computed.
        VertexContainerType pNewVertices;
        for (typename FaceContainerType::const_iterator it = r_tree.Faces().begin(); it != r_tree.Faces().end(); ++it)
        {
            if (it->State() == FaceType::NEW_BORN)
            {
                typename EdgeType::Pointer pEdge = it->pEdge();
                if (pEdge == NULL)
                    KRATOS_THROW_ERROR(std::logic_error, "The face has no edge. There are something wrong", "")
                typename EdgeType::Pointer pFirstEdge = pEdge;

                do
                {
                    pNewVertices.push_back(pEdge->pNode1());
                    pEdge = pEdge->pNextEdge();
                } while(pEdge != pFirstEdge);
            }
        }
        pNewVertices.Unique();

        // for all the new vertices, find out the parent face
        std::map<std::size_t, std::size_t> MapVertexParentFace;
        std::vector<std::size_t> ListNewNodes;
        std::vector<std::size_t> ListMovedNodes;
        for (typename VertexContainerType::const_iterator it = pNewVertices.begin(); it != pNewVertices.end(); ++it)
        {
            // check if the node existed in the model_part
            ModelPart::NodesContainerType::const_iterator it_node = r_model_part.Nodes().find(it->Id());
            if (it_node == r_model_part.Nodes().end())
            {
                typename FaceType::Pointer pFace = it->pEdge()->pFace();
                if (pFace->State() == FaceType::NEW_BORN)
                {
                    MapVertexParentFace[it->Id()] = pFace->pParent()->Id();
                    ListNewNodes.push_back(it->Id());
                }
                else if (pFace->State() == FaceType::CHANGED)
                {
                    MapVertexParentFace[it->Id()] = pFace->Id();
                    ListMovedNodes.push_back(it->Id());
                }
            }
        }

        #ifdef DEBUG_SYNCHRONIZE
        std::cout << "MapVertexParentFace:";
        for (std::map<std::size_t, std::size_t>::iterator it = MapVertexParentFace.begin(); it != MapVertexParentFace.end(); ++it)
        {
            std::cout << " (" << it->first << " " << it->second << ")";
        }
        std::cout << std::endl;
        #endif

        // in case of new node: create new nodes and add to model_part
        for (std::vector<std::size_t>::iterator it = ListNewNodes.begin(); it != ListNewNodes.end(); ++it)
        {
            // create and insert new node
            typename VertexType::Pointer pVertex = pNewVertices(*it);
            ModelPart::NodeType::Pointer pNewNode = r_model_part.CreateNewNode(pVertex->Id(), (*pVertex)[0], (*pVertex)[1], 0.0);
            #ifdef DEBUG_SYNCHRONIZE
            std::cout << *pNewNode << " is created" << std::endl;
            #endif
            mpNewNodes.insert(mpNewNodes.begin(), pNewNode);

            /* transfer data from parent element to new node */
            // find the parent element
            ModelPart::ElementsContainerType::const_iterator it_elem = r_model_part.Elements().find(MapVertexParentFace[*it]);
            if (it_elem == r_model_part.Elements().end())
            {
                KRATOS_THROW_ERROR(std::logic_error, "The parent element does not exist for node", pNewNode->Id())
            }

            mNewNodeParentElement[pNewNode->Id()] = it_elem->Id();
        }

        for (std::vector<std::size_t>::iterator it = ListMovedNodes.begin(); it != ListMovedNodes.end(); ++it)
        {
            // firstly check if the node is new after merging or the existing one
            ModelPart::NodesContainerType::iterator it_node = r_model_part.Nodes().find(*it);

            if (it_node == r_model_part.Nodes().end())
            {
                // this node existed in model_part, no transfer is required
            }
            else
            {
                // this node is created after merging operation
                typename TTreeType::ComposingVerticesContainerType::const_iterator it_composing_vertices
                        = r_tree.ComposingVertices().find(*it);
                if (it_composing_vertices != r_tree.ComposingVertices().end())
                {
                    const std::vector<std::size_t>& composing_vertices = it_composing_vertices->second;

                    mComposingVerticesData[it_composing_vertices->first] = composing_vertices;
                }
                else
                {
                    KRATOS_THROW_ERROR(std::logic_error, "The vertex is not contained in ComposingVertices information. Something's wrong", "")
                }
            }
        }

        // create new element and add to model_part
        for (typename FaceContainerType::const_iterator it = r_tree.Faces().begin(); it != r_tree.Faces().end(); ++it)
        {
            if (it->State() == FaceType::NEW_BORN)
                // a new born element is the element created after refinement. It will inherit the stresses and internal variables from the parent element.
            {
                // get the parent element
                const std::size_t& parent_element_id = it->pParent()->Id();
                ModelPart::ElementsContainerType::const_iterator it_elem = r_model_part.Elements().find(parent_element_id);

                if (it_elem == r_model_part.Elements().end())
                {
                    std::stringstream ss;
                    ss << "The parent element " << parent_element_id << " is not found in the Kratos database. Here it is expected that new born element must have some parent";
                    KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
                }

                // extract the parent element name and properties
                Properties::Pointer pProperties = it_elem->pGetProperties();
                Element::Pointer pNewElement = pCreateElementFrom(*it, *it_elem, r_model_part.Nodes(), pProperties);
                #ifdef DEBUG_SYNCHRONIZE
                std::cout << "Element " << pNewElement->Id()
                          << " of type " << typeid(*pNewElement).name()
                          << ", geometry " << typeid(*(pNewElement->pGetGeometry())).name()
                          << ", geometry type = " << pNewElement->GetGeometry().GetGeometryType()
                          << ", geometry size = " << pNewElement->GetGeometry().size()
                          << ", connectivity:";
                for (std::size_t i = 0; i < pNewElement->GetGeometry().size(); ++i)
                    std::cout << " " << pNewElement->GetGeometry()[i].Id();
                std::cout << ", is created" << std::endl;
                #endif
                mpNewElements.push_back(pNewElement);
                mNewElementParentElement[pNewElement->Id()] = parent_element_id;

                // set the state of the face to normal
                it->SetState(FaceType::NORMAL);
            }
            else if (it->State() == FaceType::CHANGED)
                // a changed element is the element that has nodes changed. It is the result of the refinement where nodes are moved or collapsed.
            {
                // delete the old element and create the new element
                const std::size_t& element_id = it->Id();
                ModelPart::ElementsContainerType::const_iterator it_elem = r_model_part.Elements().find(element_id);

                if (it_elem == r_model_part.Elements().end())
                {
                    std::stringstream ss;
                    ss << "The parent element " << element_id << " is not found in the database. Here it is expected that a changed element must still exist in the model_part";
                    KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
                }

                Properties::Pointer pProperties = it_elem->pGetProperties();
                Element::Pointer pNewElement = pCreateElementFrom(*it, *it_elem, r_model_part.Nodes(), pProperties);

                #ifdef DEBUG_SYNCHRONIZE
                std::cout << "Element " << pNewElement->Id()
                          << " of type " << typeid(*pNewElement).name()
                          << ", geometry " << typeid(*(pNewElement->pGetGeometry())).name()
                          << ", geometry type = " << pNewElement->GetGeometry().GetGeometryType()
                          << ", geometry size = " << pNewElement->GetGeometry().size()
                          << ", connectivity:";
                for (std::size_t i = 0; i < pNewElement->GetGeometry().size(); ++i)
                    std::cout << " " << pNewElement->GetGeometry()[i].Id();
                std::cout << ", is created" << std::endl;
                #endif
                mpChangedElements.push_back(pNewElement);

                // set the state of the face to normal
                it->SetState(FaceType::NORMAL);
            }
            else if (it->State() == FaceType::RESURRECTED)
                // a resurrected element is the element that is created after coarsen. It is refined and put to inactive state before. Now it is activated again, and its children elements will be deleted.
            {
                // get the children element
                ModelPart::ElementsContainerType p_children;
                for (std::size_t i = 0; i < it->HistoricalInnerFaces().size(); ++i)
                {
                    std::size_t child_element_id = it->HistoricalInnerFaces()[i];
                    ModelPart::ElementsContainerType::const_iterator it_elem = r_model_part.Elements().find(child_element_id);

                    if (it_elem == r_model_part.Elements().end())
                    {
                        std::stringstream ss;
                        ss << "The child element " << child_element_id << " is not found in the Kratos database. Here it is expected that a resurrected element must have some children";
                        KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
                    }

                    p_children.push_back(r_model_part.Elements()(child_element_id));
                }

                if (p_children.size() == 0)
                    KRATOS_THROW_ERROR(std::logic_error, "There are no inner element. Something is wrong.", "")

                // extract one child element name and properties
                Properties::Pointer pProperties = p_children.begin()->pGetProperties();
                Element::Pointer pNewElement = pCreateElementFrom(*it, *p_children.begin(), r_model_part.Nodes(), pProperties);
                #ifdef DEBUG_SYNCHRONIZE
                std::cout << "Element " << pNewElement->Id()
                          << " of type " << typeid(*pNewElement).name()
                          << ", geometry " << typeid(*(pNewElement->pGetGeometry())).name()
                          << ", geometry type = " << pNewElement->GetGeometry().GetGeometryType()
                          << ", geometry size = " << pNewElement->GetGeometry().size()
                          << ", connectivity:";
                for (std::size_t i = 0; i < pNewElement->GetGeometry().size(); ++i)
                    std::cout << " " << pNewElement->GetGeometry()[i].Id();
                std::cout << ", is created" << std::endl;
                #endif
                mpResurrectedElements.push_back(pNewElement);
                mResurrectedElementData[pNewElement->Id()] = it->HistoricalInnerFaces();

                // set the state of the face to normal
                it->SetState(FaceType::NORMAL);

                // clear the historical inner faces after used
                it->ClearHistoricalInnerFaces();
            }
        }

        std::cout << "<<<<<<<<<<<<<<<<" << __FUNCTION__ << " end" << ">>>>>>>>>>>>>>>>" << std::endl;

        return;
    }

    /**
     * Transfer the nodal variable. It's important that the node must be in its original position (i.e MoveMeshFlag = False), so that the local coordinates is found correctly.
     */
    template<class TVariableType>
    void TransferNodalVariable(ModelPart& r_model_part, TVariableType& rThisVariable)
    {
        this->TransferNodalVariableImpl(r_model_part, rThisVariable);
    }

    /**
     * End the synchronization, the old nodes and elements will be erased, as well as the internal data of this class.
     */
    template<class TTreeType>
    void EndSynchronize(ModelPart& r_model_part, const TTreeType& r_tree)
    {
        std::cout << "<<<<<<<<<<<<<<<<" << __FUNCTION__ << " begin" << ">>>>>>>>>>>>>>>>" << std::endl;

        if (m_synchronize_state != UNSYNCHRONIZED)
        {
            std::cout << "<<<<<<<<<<<<<<<<Calling " << __FUNCTION__ << " is not valid" << ">>>>>>>>>>>>>>>>" << std::endl;
            return;
        }

        typedef typename TTreeType::VertexContainerType VertexContainerType;
        typedef typename TTreeType::EdgeContainerType EdgeContainerType;
        typedef typename TTreeType::FaceContainerType FaceContainerType;

        typedef typename TTreeType::VertexType VertexType;
        typedef typename TTreeType::EdgeType EdgeType;
        typedef typename TTreeType::FaceType FaceType;

        // remove the lone nodes
        std::vector<std::size_t> list_removed_nodes_lonely;
        for (typename VertexContainerType::const_iterator it = r_tree.Vertices().begin(); it != r_tree.Vertices().end(); ++it)
        {
            if (it->IsLonely())
                list_removed_nodes_lonely.push_back(it->Id());
        }

        // remove nodes that are in model_part but not in the tree. Those nodes are deleted from the tree from coarsen operation
        std::vector<std::size_t> list_removed_nodes_delete;
        for (ModelPart::NodesContainerType::const_iterator it = r_model_part.Nodes().begin(); it != r_model_part.Nodes().end(); ++it)
        {
            typename VertexContainerType::const_iterator it_vertex = r_tree.Vertices().find(it->Id());

            if (it_vertex == r_tree.Vertices().end())
                list_removed_nodes_delete.push_back(it->Id());
        }

        // remove the nodes eventually
        for (std::vector<std::size_t>::iterator it = list_removed_nodes_lonely.begin(); it != list_removed_nodes_lonely.end(); ++it)
        {
            r_model_part.RemoveNode(*it);
            #ifdef DEBUG_SYNCHRONIZE
            std::cout << "Node " << *it << " is removed from model_part because it's lonely" << std::endl;
            #endif
        }
        for (std::vector<std::size_t>::iterator it = list_removed_nodes_delete.begin(); it != list_removed_nodes_delete.end(); ++it)
        {
            r_model_part.RemoveNode(*it);
            #ifdef DEBUG_SYNCHRONIZE
            std::cout << "Node " << *it << " is removed from model_part because it's deleted from the tree" << std::endl;
            #endif
        }

        r_model_part.Nodes().Unique(); // sort the node container

        #ifdef DEBUG_SYNCHRONIZE
        std::cout << "Start to remove elements" << std::endl;
        KRATOS_WATCH(mpNewElements.size())
        KRATOS_WATCH(mpChangedElements.size())
        KRATOS_WATCH(mpResurrectedElements.size())
        #endif

        // for sub-elements: add new element to model part
        for (ModelPart::ElementsContainerType::ptr_iterator it = mpNewElements.ptr_begin(); it != mpNewElements.ptr_end(); ++it)
        {
            r_model_part.Elements().push_back(*it);
            std::cout << "Element " << (*it)->Id() << " is added to the model_part" << std::endl;
        }

        r_model_part.Elements().Unique(); // sort the element container

        // for changed elements: remove the old and replace by new one
        for (ModelPart::ElementsContainerType::ptr_iterator it = mpChangedElements.ptr_begin(); it != mpChangedElements.ptr_end(); ++it)
        {
            r_model_part.RemoveElement((*it)->Id());
            std::cout << "Removed element " << (*it)->Id() << std::endl;
            r_model_part.Elements().push_back(*it);
            std::cout << "Element " << (*it)->Id() << " is added" << std::endl;
        }

        r_model_part.Elements().Unique(); // sort the element container

        // for resurrected elements: remove the children and add the resurrected one
        for (ModelPart::ElementsContainerType::ptr_iterator it = mpResurrectedElements.ptr_begin(); it != mpResurrectedElements.ptr_end(); ++it)
        {
            // remove the children elements from model_part
            for (std::size_t i = 0; i < mResurrectedElementData[(*it)->Id()].size(); ++i)
                r_model_part.RemoveElement(mResurrectedElementData[(*it)->Id()][i]);

            r_model_part.Elements().push_back(*it);
        }

        r_model_part.Elements().Unique(); // sort the element container

        // remove the inactive elements
        for (typename FaceContainerType::const_iterator it = r_tree.Faces().begin(); it != r_tree.Faces().end(); ++it)
        {
            if (it->State() == FaceType::SLEEP)
            {
                const std::size_t& element_id = it->Id();
                ModelPart::ElementsContainerType::const_iterator it_elem = r_model_part.Elements().find(element_id);

                if (it_elem == r_model_part.Elements().end())
                {
                    // std::stringstream ss;
                    // ss << "The parent element " << element_id << " is not found in the database.";
                    // KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
                    continue;
                }

                r_model_part.RemoveElement(element_id);
                #ifdef DEBUG_SYNCHRONIZE
                std::cout << "The element " << element_id << " is removed from the model_part" << std::endl;
                #endif
            }
        }

        r_model_part.Elements().Unique(); // sort the element container

        m_synchronize_state = SYNCHRONIZED;

        // clean the internal data structure
        mpNewNodes.clear();
        std::cout << "New nodes container is cleared" << std::endl;
        mNewNodeParentElement.clear();
        mComposingVerticesData.clear();

        mpNewElements.clear();
        std::cout << "New elements container is cleared" << std::endl;
        mNewElementParentElement.clear();
        mpChangedElements.clear();
        std::cout << "Changed elements container is cleared" << std::endl;
        mpResurrectedElements.clear();
        mResurrectedElementData.clear();
        std::cout << "Resurrected elements container is cleared" << std::endl;

        std::cout << "<<<<<<<<<<<<<<<<" << __FUNCTION__ << " end" << ">>>>>>>>>>>>>>>>" << std::endl;

        return;
    }

    ///@}
    ///@name Access
    ///@{

    ModelPart::NodesContainerType GetNewNodes()
    {
        return mpNewNodes;
    }

    ModelPart::ElementsContainerType GetNewElements()
    {
        return mpNewElements;
    }

    ModelPart::ElementsContainerType GetChangedElements()
    {
        return mpChangedElements;
    }

    ModelPart::ElementsContainerType GetResurrectedElements()
    {
        return mpResurrectedElements;
    }

    /**
     * Probe the status of the polytree
     */
    template<class TTreeType>
    void QueryState(typename TTreeType::Pointer pTree) const
    {
        typedef typename TTreeType::VertexContainerType VertexContainerType;
        typedef typename TTreeType::EdgeContainerType EdgeContainerType;
        typedef typename TTreeType::FaceContainerType FaceContainerType;

        typedef typename TTreeType::VertexType VertexType;
        typedef typename TTreeType::EdgeType EdgeType;
        typedef typename TTreeType::FaceType FaceType;

        std::cout << "Face states:" << std::endl;
        for (typename FaceContainerType::const_iterator it = pTree->Faces().begin(); it != pTree->Faces().end(); ++it)
        {
            if (it->State() == FaceType::NEW_BORN)
            {
                std::cout << "  Face " << it->Id() << " is new" << std::endl;
            }
            else if (it->State() == FaceType::CHANGED)
            {
                std::cout << "  Face " << it->Id() << " is changed" << std::endl;
            }
            else if (it->State() == FaceType::NORMAL)
            {
                std::cout << "  Face " << it->Id() << " is normal" << std::endl;
            }
            else if (it->State() == FaceType::SLEEP)
            {
                std::cout << "  Face " << it->Id() << " is sleep" << std::endl;
            }
        }
    }


    void ListModelPart(ModelPart& r_model_part) const
    {
        std::cout << "Nodes:" << std::endl;
        for (ModelPart::NodesContainerType::const_iterator it = r_model_part.Nodes().begin(); it != r_model_part.Nodes().end(); ++it)
        {
            std::cout << "  " << it->Id() << ": (" << it->X0() << " " << it->Y0() << " " << it->Z0() << ")" << std::endl;
        }

        std::cout << "Elements:" << std::endl;
        for (ModelPart::ElementsContainerType::const_iterator it = r_model_part.Elements().begin(); it != r_model_part.Elements().end(); ++it)
        {
            std::cout << "  " << it->Id() << ":";
            for (std::size_t i = 0; i < it->GetGeometry().size(); ++i)
                std::cout << " " << it->GetGeometry()[i].Id();
            std::cout << std::endl;
        }
        // for ( ModelPart::MeshType::ElementIterator element_iterator = r_model_part.GetMesh().ElementsBegin();
        //             element_iterator != r_model_part.GetMesh().ElementsEnd(); ++element_iterator)
        // {
        //     std::cout << "  " << element_iterator->Id() << ":";
        //     for (std::size_t i = 0; i < element_iterator->GetGeometry().size(); ++i)
        //         std::cout << " " << element_iterator->GetGeometry()[i].Id();
        //     std::cout << std::endl;
        // }
    }

    ///@}
    ///@name Inquiry
    ///@{


    /**
     * Compare a model part and a polytree (only vertices and faces)
     * @param r_model_part input model_part
     * @param r_tree       input polytree
     */
    template<class TTreeType>
    int CompareForward(ModelPart& r_model_part, TTreeType& r_tree) const
    {
        typedef typename TTreeType::VertexContainerType VertexContainerType;
        typedef typename TTreeType::EdgeContainerType EdgeContainerType;
        typedef typename TTreeType::FaceContainerType FaceContainerType;

        typedef typename TTreeType::VertexType VertexType;
        typedef typename TTreeType::EdgeType EdgeType;
        typedef typename TTreeType::FaceType FaceType;

        int error = 0;

        // check the node in the model_part
        for (typename VertexContainerType::const_iterator it = r_tree.Vertices().begin(); it != r_tree.Vertices().end(); ++it)
        {
            if (it->pEdge() != NULL) // only check the non-lonely node
            {
                ModelPart::NodesContainerType::const_iterator it_node = r_model_part.Nodes().find(it->Id());

                if (it_node == r_model_part.Nodes().end())
                {
                    std::cout << "The vertex " << it->Id() << " does not have respective node in the model_part" << std::endl;
                    error |= 0x0001;
                    continue;
                }
            }
        }

        // check the element in the model_part
        for (typename FaceContainerType::const_iterator it = r_tree.Faces().begin(); it != r_tree.Faces().end(); ++it)
        {
            if (it->IsActive())
            {
                ModelPart::ElementsContainerType::const_iterator it_elem = r_model_part.Elements().find(it->Id());

                if (it_elem == r_model_part.Elements().end())
                {
                    std::cout << "The face " << it->Id() << " does not have respective element in the model_part" << std::endl;
                    error |= 0x0010;
                    continue;
                }

                // check the element connectivity
                std::vector<std::size_t> connectivity = it->GetConnectivity();
                if (connectivity.size() != it_elem->GetGeometry().size())
                {
                    std::cout << "The element " << it_elem->Id() << " and face " << it->Id() << " has different number of nodes" << std::endl;
                    error |= 0x0100;
                    continue;
                }
                for (std::size_t i = 0; i < connectivity.size(); ++i)
                {
                    if (it_elem->GetGeometry()[i].Id() != connectivity[i])
                    {
                        std::cout << "The connectivity of the element is different at " << i << std::endl;
                        error |= 0x0100;
                        continue;
                    }
                }
            }
        }

        std::cout << __FUNCTION__ << " is successfully.";
        if (error == 0)
        {
            std::cout << " The polytree is contained in the model_part." << std::endl;
        }
        else
        {
            std::cout << " The polytree is NOT contained in the model_part." << std::endl;
        }

        return error;
    }


    /**
     * Compare a model part and a polytree (only nodes and elements)
     * @param r_model_part input model_part
     * @param r_tree       input polytree
     */
    template<class TTreeType>
    int CompareBackward(ModelPart& r_model_part, TTreeType& r_tree) const
    {
        typedef typename TTreeType::VertexContainerType VertexContainerType;
        typedef typename TTreeType::EdgeContainerType EdgeContainerType;
        typedef typename TTreeType::FaceContainerType FaceContainerType;

        typedef typename TTreeType::VertexType VertexType;
        typedef typename TTreeType::EdgeType EdgeType;
        typedef typename TTreeType::FaceType FaceType;

        int error = 0;

        // check the node in the model_part
        for (ModelPart::NodesContainerType::const_iterator it = r_model_part.Nodes().begin(); it != r_model_part.Nodes().end(); ++it)
        {
            typename VertexContainerType::const_iterator it_vertex = r_tree.Vertices().find(it->Id());

            if (it_vertex == r_tree.Vertices().end())
            {
                std::cout << "The node " << it->Id() << " does not have respective vertex in the polytree" << std::endl;
                error |= 0x0001;
                continue;
            }
        }

        // check the element in the model_part
        for (ModelPart::ElementsContainerType::const_iterator it = r_model_part.Elements().begin(); it != r_model_part.Elements().end(); ++it)
        {
            typename FaceContainerType::const_iterator it_face = r_tree.Faces().find(it->Id());

            if (it_face == r_tree.Faces().end())
            {
                std::cout << "The element " << it->Id() << " does not have respective face in the polytree" << std::endl;
                error |= 0x0010;
                continue;
            }

            // check the element connectivity
            std::vector<std::size_t> connectivity = it_face->GetConnectivity();
            if (connectivity.size() != it->GetGeometry().size())
            {
                std::cout << "The element " << it->Id() << " and face " << it_face->Id() << " has different number of nodes" << std::endl;
                error |= 0x0100;
                continue;
            }
            for (std::size_t i = 0; i < connectivity.size(); ++i)
            {
                if (it->GetGeometry()[i].Id() != connectivity[i])
                {
                    std::cout << "The connectivity of the element is different at " << i << std::endl;
                    error |= 0x0100;
                    continue;
                }
            }
        }

        std::cout << __FUNCTION__ << " is successfully.";
        if (error == 0)
        {
            std::cout << " The model_part is contained in the polytree." << std::endl;
        }
        else
        {
            std::cout << " The model_part is NOT contained in the polytree." << std::endl;
        }

        return error;
    }


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "PolyTree Synchronization Utility";
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


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:

    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    SynchronizeStates m_synchronize_state;

    ModelPart::NodesContainerType mpNewNodes;
    std::map<std::size_t, std::size_t> mNewNodeParentElement;
    std::map<std::size_t, std::vector<std::size_t> > mComposingVerticesData;

    ModelPart::ElementsContainerType mpNewElements;
    std::map<std::size_t, std::size_t> mNewElementParentElement;
    ModelPart::ElementsContainerType mpChangedElements;
    ModelPart::ElementsContainerType mpResurrectedElements;
    std::map<std::size_t, std::vector<std::size_t> > mResurrectedElementData;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * Create element from a face
     */
    template<class FaceType, class NodesContainerType>
    Element::Pointer pCreateElementFrom(FaceType& r_face,
            Element const& r_clone_element,
            NodesContainerType& rThisNodes,
            Properties::Pointer pProperties) const
    {
        // get the connectivity from the face
        std::vector<std::size_t> connectivity = r_face.GetConnectivity();

        // get the array of nodes
        Element::NodesArrayType temp_element_nodes;
        for (std::size_t i = 0; i < connectivity.size(); ++i)
            temp_element_nodes.push_back( *(FindKey(rThisNodes, connectivity[i], "Node").base()) );

        // create new element
        return r_clone_element.Create(r_face.Id(), temp_element_nodes, pProperties);
    }

    template<class TContainerType, class TKeyType>
    typename TContainerType::iterator FindKey(TContainerType& ThisContainer , TKeyType ThisKey, std::string ComponentName) const
    {
        typename TContainerType::iterator i_result;
        if((i_result = ThisContainer.find(ThisKey)) == ThisContainer.end())
        {
            std::stringstream buffer;
            buffer << ComponentName << " #" << ThisKey << " is not found.";
            KRATOS_THROW_ERROR(std::invalid_argument, buffer.str(), "");
        }

        return i_result;
    }

    void TransferNodalVariableImpl(ModelPart& r_model_part, Variable<double>& rThisVariable)
    {
        // transfer data for new node
        for (std::map<std::size_t, std::size_t>::iterator it = mNewNodeParentElement.begin();
                it != mNewNodeParentElement.end(); ++it)
        {
            std::size_t node_id = it->first;
            std::size_t parent_element_id = it->second;

            ModelPart::NodeType::Pointer pNewNode = r_model_part.Nodes()(node_id);
            Element::Pointer pParentElement = r_model_part.Elements()(parent_element_id);

            // find the local coordinates of the node within the parent element
            CoordinatesArrayType local_coords;
            local_coords = pParentElement->GetGeometry().PointLocalCoordinates(local_coords, *pNewNode);
            #ifdef DEBUG_SYNCHRONIZE
            std::cout << "local coordinates of new node " << pNewNode->Id() << " within parent element " << pParentElement->Id() << ": " << local_coords << std::endl;
            std::cout << " node on element " << pParentElement->Id() << ":" << std::endl;
            for (unsigned int i = 0; i < pParentElement->GetGeometry().size(); ++i)
            {
                std::cout << "  " << pParentElement->GetGeometry()[i].Id() << ": " << pParentElement->GetGeometry()[i].X0() << " " << pParentElement->GetGeometry()[i].Y0() << std::endl;
            }
            CoordinatesArrayType global_coords;
            global_coords = pParentElement->GetGeometry().GlobalCoordinates(global_coords, local_coords);
            std::cout << " recheck the global coordinates: " << global_coords << std::endl;
            #endif

            // transfer the data to new node
            double new_value = 0.0;
            Vector ShapeFunctionValues;
            ShapeFunctionValues = pParentElement->GetGeometry().ShapeFunctionsValues(ShapeFunctionValues, local_coords);
            for (std::size_t j = 0; j < pParentElement->GetGeometry().size(); ++j)
            {
                new_value += ShapeFunctionValues[j] * pParentElement->GetGeometry()[j].GetSolutionStepValue(rThisVariable);
            }
            pNewNode->SetValue(rThisVariable, new_value);
            #ifdef DEBUG_SYNCHRONIZE
            std::cout << " Node " << pNewNode->Id() << " has new " << rThisVariable.Name() << " value: " << new_value << std::endl;
            #endif
        }

        // transfer data from composing Vertices
        for (std::map<std::size_t, std::vector<size_t> >::iterator it = mComposingVerticesData.begin();
                it != mComposingVerticesData.end(); ++it)
        {
            std::size_t node_id = it->first;
            std::vector<std::size_t> composing_vertices = it->second;

            ModelPart::NodeType::Pointer pNewNode = r_model_part.Nodes()(node_id);

            double new_value = 0.0;
            for (std::size_t j = 0; j < composing_vertices.size(); ++j)
            {
                new_value += r_model_part.Nodes()[composing_vertices[j]].GetSolutionStepValue(rThisVariable);
            }
            new_value /= composing_vertices.size();
            pNewNode->SetValue(rThisVariable, new_value);
            #ifdef DEBUG_SYNCHRONIZE
            std::cout << " Node " << pNewNode->Id() << " has new " << rThisVariable.Name() << " value: " << new_value << std::endl;
            #endif
        }
    }

    void TransferNodalVariableImpl(ModelPart& r_model_part, Variable<array_1d<double, 3> >& rThisVariable)
    {
        // transfer data for new node
        for (std::map<std::size_t, std::size_t>::iterator it = mNewNodeParentElement.begin();
                it != mNewNodeParentElement.end(); ++it)
        {
            std::size_t node_id = it->first;
            std::size_t parent_element_id = it->second;

            ModelPart::NodeType::Pointer pNewNode = r_model_part.Nodes()(node_id);
            Element::Pointer pParentElement = r_model_part.Elements()(parent_element_id);

            // find the local coordinates of the node within the parent element
            CoordinatesArrayType local_coords;
            local_coords = pParentElement->GetGeometry().PointLocalCoordinates(local_coords, *pNewNode);
            #ifdef DEBUG_SYNCHRONIZE
            std::cout << "local coordinates of new node " << pNewNode->Id() << " within parent element " << pParentElement->Id() << ": " << local_coords << std::endl;
            std::cout << " node on element " << pParentElement->Id() << ":" << std::endl;
            for (unsigned int i = 0; i < pParentElement->GetGeometry().size(); ++i)
            {
                std::cout << "  " << pParentElement->GetGeometry()[i].Id() << ": " << pParentElement->GetGeometry()[i].X0() << " " << pParentElement->GetGeometry()[i].Y0();
                std::cout << ", " << rThisVariable.Name() << ": " << pParentElement->GetGeometry()[i].GetSolutionStepValue(rThisVariable);
                std::cout << std::endl;
            }
            CoordinatesArrayType global_coords;
            global_coords = pParentElement->GetGeometry().GlobalCoordinates(global_coords, local_coords);
            std::cout << " recheck the global coordinates: " << global_coords << std::endl;
            #endif

            // transfer the data to new node
            array_1d<double, 3> new_value;
            noalias(new_value) = ZeroVector(3);
            Vector ShapeFunctionValues;
            ShapeFunctionValues = pParentElement->GetGeometry().ShapeFunctionsValues(ShapeFunctionValues, local_coords);
            for (std::size_t j = 0; j < pParentElement->GetGeometry().size(); ++j)
            {
                new_value += ShapeFunctionValues[j] * pParentElement->GetGeometry()[j].GetSolutionStepValue(rThisVariable);
            }
            pNewNode->SetValue(rThisVariable, new_value);
            #ifdef DEBUG_SYNCHRONIZE
            std::cout << " Node " << pNewNode->Id() << " has new " << rThisVariable.Name() << " value: " << new_value << std::endl;
            #endif
        }

        // transfer data from composing Vertices
        for (std::map<std::size_t, std::vector<size_t> >::iterator it = mComposingVerticesData.begin();
                it != mComposingVerticesData.end(); ++it)
        {
            std::size_t node_id = it->first;
            std::vector<std::size_t> composing_vertices = it->second;

            ModelPart::NodeType::Pointer pNewNode = r_model_part.Nodes()(node_id);

            array_1d<double, 3> new_value;
            noalias(new_value) = ZeroVector(3);
            for (std::size_t j = 0; j < composing_vertices.size(); ++j)
            {
                new_value += r_model_part.Nodes()[composing_vertices[j]].GetSolutionStepValue(rThisVariable);
            }
            new_value /= composing_vertices.size();
            pNewNode->SetValue(rThisVariable, new_value);
            #ifdef DEBUG_SYNCHRONIZE
            std::cout << " Node " << pNewNode->Id() << " has new " << rThisVariable.Name() << " value: " << new_value << std::endl;
            #endif
        }
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Unaccessible methods
    ///@{

    /// Assignment operator.
    PolyTreeSyncUtility& operator=(PolyTreeSyncUtility const& rOther);

    /// Copy constructor.
    PolyTreeSyncUtility(PolyTreeSyncUtility const& rOther);


    ///@}

}; // Class PolyTreeSyncUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, PolyTreeSyncUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const PolyTreeSyncUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#undef DEBUG_SYNCHRONIZE

#endif // KRATOS_POLYTREE_SYNC_UTILITY_H_INCLUDED  defined
