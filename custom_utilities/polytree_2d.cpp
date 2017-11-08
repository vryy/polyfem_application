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

// System includes


// External includes


// Project includes
#include "polytree_2d.h"
#include "polytree_utility.h"
#include "polyfem_application/polyfem_application.h"

// #define DEBUG_REFINE
// #define DEBUG_COARSEN

namespace Kratos
{

void PolyTree2D::Clear()
{
    // clear all the containers
    mpVertexList->clear();
    mpEdgeList->clear();
    mpFaceList->clear();
    mNewEdgeList.clear();
    mCompositeEdgeList.clear();
    mComposingVertices.clear();
    m_half_edge_state = UNINITIALIZED;
}

PolyTree2D::FaceType::Pointer PolyTree2D::CreateFace(const std::size_t& Id)
{
    FaceType::Pointer pFace = boost::make_shared<FaceType>(Id);
    mpFaceList->insert(mpFaceList->begin(), pFace);
    return pFace;
}

int PolyTree2D::MarkFaceRefine(const std::size_t& face_index)
{
    FaceContainerType::iterator it = mpFaceList->find(face_index);
    if (it != mpFaceList->end())
    {
        if (it->IsActive())
        {
            it->SetRefine(true);
            it->SetCoarsen(false);
            #ifdef DEBUG_REFINE
            std::cout << *it << " is marked to be refined" << std::endl;
            #endif // DEBUG_REFINE
            if (m_half_edge_state == CACHED_COARSEN)
            {
                KRATOS_THROW_ERROR(std::logic_error, "The refine and coarsen operations are not allowed at the same time", "")
            }
            else
                m_half_edge_state = CACHED_REFINE;
            return 0;
        }
        else
            return -1;
    }
    else
        return -2;
}

int PolyTree2D::MarkFaceCoarsen(const std::size_t& face_index)
{
    FaceContainerType::iterator it = mpFaceList->find(face_index);
    if (it != mpFaceList->end())
    {
        if (!it->IsActive())
        {
            it->SetRefine(false);
            it->SetCoarsen(true);
            #ifdef DEBUG_COARSEN
            std::cout << *it << " is marked to be coarsen" << std::endl;
            #endif // DEBUG_COARSEN
            if (m_half_edge_state == CACHED_REFINE)
            {
                KRATOS_THROW_ERROR(std::logic_error, "The refine and coarsen operations are not allowed at the same time", "")
            }
            else
                m_half_edge_state = CACHED_COARSEN;
            return 0;
        }
        else
            return -1;
    }
    else
        return -2;
}

void PolyTree2D::BeginRefineCoarsen()
{
    std::cout << __FUNCTION__ << " starts" << std::endl;

    if (m_half_edge_state == CACHED_REFINE || m_half_edge_state == CACHED_COARSEN)
    {
        // firstly get all the indices of the face
        std::set<std::size_t> all_face_indices;
        for (FaceContainerType::iterator it = mpFaceList->begin(); it != mpFaceList->end(); ++it)
        {
            all_face_indices.insert(it->Id());
        }

        // then refine those face in that indices
        for (std::set<std::size_t>::iterator is = all_face_indices.begin(); is != all_face_indices.end(); ++is)
        {
            FaceContainerType::iterator it = mpFaceList->find(*is);
            if (it == mpFaceList->end())
                continue;

            FaceType::Pointer pFace = (*mpFaceList)(*is);

            if (pFace->IsRefined())
            {
                if (pFace->IsLeaf())
                {
                    RefineFace(pFace, *mpVertexList, *mpEdgeList, mNewEdgeList, mCompositeEdgeList, *mpFaceList, mLastVertexId, mLastFaceId);
                    mpEdgeList->Sort(); // sort to ensure find working correctly
                }
                // here we shall not refine the sub-faces. That can create enormous refinement if user call MarkFaceRefine again after BeginRefineCoarsen
                // else
                // {
                //     for (std::size_t i = 0; i < it->NumberOfSubFaces(); ++i)
                //     {
                //         RefineFace(it->SubFace(i), *mpFaceList, mLastVertexId, mLastFaceId);
                //     }
                // }

                // turn off the refine flag
                pFace->SetRefine(false);

                // set this face as inactive
                pFace->SetActive(false);
                pFace->SetState(FaceType::SLEEP);
            }

            if (pFace->IsCoarsen())
            {
                if (!pFace->IsLeaf())
                {
                    // collect all the edges around the face
                    EdgeContainerType pOuterEdges;
                    EdgeContainerType pInnerEdges;
                    VertexContainerType pInnerVertices;
                    FaceContainerType pInnerFaces;

                    ExtractEdgesAndVertices(pFace, pOuterEdges, pInnerEdges, pInnerVertices, pInnerFaces);

                    #ifdef DEBUG_COARSEN
                    std::cout << "list of inner faces of " << *pFace << ":" << std::endl;
                    for (FaceContainerType::const_iterator it = pInnerFaces.begin(); it != pInnerFaces.end(); ++it)
                        std::cout << "  " << *it << std::endl;
                    std::cout << std::endl;
                    #endif // DEBUG_COARSEN

                    #ifdef DEBUG_COARSEN
                    std::cout << "list of outer edges of " << *pFace << ":" << std::endl;
                    for (EdgeContainerType::const_iterator it = pOuterEdges.begin(); it != pOuterEdges.end(); ++it)
                        std::cout << "  " << *it << std::endl;
                    std::cout << std::endl;
                    #endif // DEBUG_COARSEN

                    // remove all the sub-faces and edges, also the inner vertices
                    std::vector<EdgeGetKeyType::result_type> removed_edge_keys;
                    for (EdgeContainerType::const_iterator it = pInnerEdges.begin(); it != pInnerEdges.end(); ++it)
                        removed_edge_keys.push_back(EdgeGetKeyType()(*it));
                    for (std::size_t i = 0; i < removed_edge_keys.size(); ++i)
                    {
                        mpEdgeList->erase(removed_edge_keys[i]);
                        #ifdef DEBUG_COARSEN
                        std::cout << "  edge(" << removed_edge_keys[i].first << " " << removed_edge_keys[i].second << ") is removed" << std::endl;
                        #endif // DEBUG_COARSEN
                    }
                    mpEdgeList->Sort(); // sort after erasing

                    std::vector<std::size_t> removed_vertex_keys;
                    for (VertexContainerType::const_iterator it = pInnerVertices.begin(); it != pInnerVertices.end(); ++it)
                        removed_vertex_keys.push_back(it->Id());
                    for (std::size_t i = 0; i < removed_vertex_keys.size(); ++i)
                    {
                        mpVertexList->erase(removed_vertex_keys[i]);
                        #ifdef DEBUG_COARSEN
                        std::cout << "  vertex " << removed_vertex_keys[i] << " is removed" << std::endl;
                        #endif // DEBUG_COARSEN
                    }
                    mpVertexList->Sort();

                    std::vector<std::size_t> removed_face_keys;
                    removed_face_keys.reserve(pInnerFaces.size());
                    for (FaceContainerType::const_iterator it = pInnerFaces.begin(); it != pInnerFaces.end(); ++it)
                        removed_face_keys.push_back(it->Id());
                    for (std::size_t i = 0; i < removed_face_keys.size(); ++i)
                    {
                        mpFaceList->erase(removed_face_keys[i]);
                        #ifdef DEBUG_COARSEN
                        std::cout << "  face " << removed_face_keys[i] << " is removed" << std::endl;
                        #endif // DEBUG_COARSEN
                    }
                    mpFaceList->Sort();

                    // clear the sub-faces, and also remove it from the *mpFaceList
                    std::set<std::size_t> removed_sub_faces;
                    pFace->ClearSubFaces(removed_sub_faces);

                    for (std::set<std::size_t>::iterator it = removed_sub_faces.begin(); it != removed_sub_faces.end(); ++it)
                    {
                        mpFaceList->erase(*it);
                    }

                    mpFaceList->Sort();

                    // create a new face from the round edges, also re-set the next edge for the outer edges
                    pFace->pSetEdge(*(pOuterEdges.ptr_begin()));
                    pFace->SetActive(true);
                    pFace->SetState(FaceType::RESURRECTED);
                    for (EdgeContainerType::ptr_iterator it = pOuterEdges.ptr_begin(); it != pOuterEdges.ptr_end(); ++it)
                    {
                        EdgeContainerType::ptr_iterator it2 = it+1;
                        if (it2 == pOuterEdges.ptr_end()) it2 = pOuterEdges.ptr_begin();
                        (*it)->pSetNextEdge(*it2);
                        (*it2)->pSetPrevEdge(*it);
                        (*it)->pSetFace(pFace);
                    }

                    #ifdef DEBUG_COARSEN
                    std::cout << *pFace << " is modified after coarsening" << std::endl;
                    #endif // DEBUG_COARSEN

                    // re-set the edge of the vertex
                    // this will also re-set the edge for valid vertex, but it's OK. It will take more time and more code to select the valid/invalid here
                    EdgeType::Pointer pEdge = pFace->pEdge();
                    EdgeType::Pointer pFirstEdge = pEdge;
                    do
                    {
                        pEdge->pNode1()->pSetEdge(pEdge);
                        pEdge = pEdge->pNextEdge();
                    } while (pEdge != pFirstEdge);

                    // set historical inner faces
                    pFace->SetHistoricalInnerFaces(removed_face_keys);

                    #ifdef DEBUG_COARSEN
                    std::cout << *pFace << " vertices:" << std::endl;
                    pEdge = pFace->pEdge();
                    pFirstEdge = pEdge;
                    do
                    {
                        std::cout << "  " << *(pEdge->pNode1()) << std::endl;
                        pEdge = pEdge->pNextEdge();
                    } while (pEdge != pFirstEdge);

                    std::cout << *pFace << " edges:" << std::endl;
                    pEdge = pFace->pEdge();
                    pFirstEdge = pEdge;
                    do
                    {
                        std::cout << "  " << *pEdge << std::endl;
                        pEdge = pEdge->pNextEdge();
                    } while (pEdge != pFirstEdge);
                    #endif // DEBUG_COARSEN
                }
                else
                {
                    // shall do nothing because this is a leaf face
                }

                // turn off the coarsen flag
                pFace->SetCoarsen(false);
            }
        }

        mpVertexList->Unique();
        mpEdgeList->Unique();
        mNewEdgeList.Unique();
        mCompositeEdgeList.Unique();
        mpFaceList->Unique();

        m_half_edge_state = FINALIZED_READY;
    }
}

void PolyTree2D::EndRefineCoarsen()
{
    std::cout << __FUNCTION__ << " starts" << std::endl;
    if (m_half_edge_state == FINALIZED_READY)
    {
        #ifdef DEBUG_REFINE
        std::cout << "  Faces before merging:" << std::endl;
        for (FaceContainerType::const_iterator it = mpFaceList->begin(); it != mpFaceList->end(); ++it)
        {
            std::cout << "    Face " << it->Id() << " edges: " << it->EdgeInfo() << std::endl;
        }

        std::cout << "  Edges before merging:" << std::endl;
        for (EdgeContainerType::ptr_const_iterator it = mpEdgeList->ptr_begin(); it != mpEdgeList->ptr_end(); ++it)
        {
            std::cout << "    " << *it << ": " << *(*it) << std::endl;
        }

        std::cout << "  New edges before merging:" << std::endl;
        for (EdgeContainerType::ptr_const_iterator it = mNewEdgeList.ptr_begin(); it != mNewEdgeList.ptr_end(); ++it)
        {
            std::cout << "    " << *it << ": " << *(*it) << std::endl;
        }

        std::cout << "  Composite edges before merging:" << std::endl;
        for (EdgeContainerType::const_iterator it = mCompositeEdgeList.begin(); it != mCompositeEdgeList.end(); ++it)
        {
            std::cout << "    " << *it << std::endl;
        }

        std::cout << "  Vertices before merging:" << std::endl;
        for (VertexContainerType::const_iterator it = mpVertexList->begin(); it != mpVertexList->end(); ++it)
        {
            std::cout << "    " << *it << std::endl;
        }
        #endif // DEBUG_REFINE

        // merge the composite edges. This is to reduce number of near points in the polygon mesh.
        // here we loop through the list of composite edges and make the pair with its opposite
        std::vector<std::pair<EdgeType::Pointer, EdgeType::Pointer> > EdgePairs;

        std::vector<std::pair<EdgeGetKeyType::result_type, EdgeGetKeyType::result_type> > EdgePairKeys1;
        std::vector<std::pair<EdgeGetKeyType::result_type, EdgeGetKeyType::result_type> > EdgePairKeys2;
        for (EdgeContainerType::const_iterator it = mCompositeEdgeList.begin();
                it != mCompositeEdgeList.end(); ++it)
        {
            EdgeGetKeyType::result_type key = EdgeGetKeyType()(*it);
            EdgeGetKeyType::result_type reversed_key(key.second, key.first);

            EdgeContainerType::const_iterator it2 = mCompositeEdgeList.find(reversed_key);
            if (it2 != mCompositeEdgeList.end())
            {
                EdgePairKeys1.push_back(std::pair<EdgeGetKeyType::result_type, EdgeGetKeyType::result_type>(key, reversed_key));
            }
            else
            {
                EdgeContainerType::const_iterator it3 = mpEdgeList->find(reversed_key);
                if (it3 != mpEdgeList->end())
                {
                    EdgePairKeys2.push_back(std::pair<EdgeGetKeyType::result_type, EdgeGetKeyType::result_type>(key, reversed_key));
                }
            }
        }

        // for (std::vector<std::pair<EdgeGetKeyType::result_type, EdgeGetKeyType::result_type> >::iterator it = EdgePairKeys1.begin();
        //         it != EdgePairKeys1.end(); ++it)
        // {
        //     EdgePairs.push_back(std::pair<EdgeType::Pointer, EdgeType::Pointer>(mCompositeEdgeList(it->first), mCompositeEdgeList(it->second)));
        // }

        // for (std::vector<std::pair<EdgeGetKeyType::result_type, EdgeGetKeyType::result_type> >::iterator it = EdgePairKeys2.begin();
        //         it != EdgePairKeys2.end(); ++it)
        // {
        //     EdgePairs.push_back(std::pair<EdgeType::Pointer, EdgeType::Pointer>(mCompositeEdgeList(it->first), (*mpEdgeList)(it->second)));
        // }

        for (std::vector<std::pair<EdgeGetKeyType::result_type, EdgeGetKeyType::result_type> >::iterator it = EdgePairKeys1.begin();
                it != EdgePairKeys1.end(); ++it)
        {
            std::pair<EdgeType::Pointer, EdgeType::Pointer> NewPair = std::pair<EdgeType::Pointer, EdgeType::Pointer>(mCompositeEdgeList(it->first), mCompositeEdgeList(it->second));
            std::pair<EdgeType::Pointer, EdgeType::Pointer> ReversedPair = std::pair<EdgeType::Pointer, EdgeType::Pointer>(mCompositeEdgeList(it->second), mCompositeEdgeList(it->first));

            if (std::find(EdgePairs.begin(), EdgePairs.end(), NewPair) == EdgePairs.end()
                && std::find(EdgePairs.begin(), EdgePairs.end(), ReversedPair) == EdgePairs.end())
            {
                EdgePairs.push_back(NewPair);
            }
        }

        for (std::vector<std::pair<EdgeGetKeyType::result_type, EdgeGetKeyType::result_type> >::iterator it = EdgePairKeys2.begin();
                it != EdgePairKeys2.end(); ++it)
        {
            std::pair<EdgeType::Pointer, EdgeType::Pointer> NewPair = std::pair<EdgeType::Pointer, EdgeType::Pointer>(mCompositeEdgeList(it->first), (*mpEdgeList)(it->second));
            std::pair<EdgeType::Pointer, EdgeType::Pointer> ReversedPair = std::pair<EdgeType::Pointer, EdgeType::Pointer>((*mpEdgeList)(it->second), mCompositeEdgeList(it->first));

            if (std::find(EdgePairs.begin(), EdgePairs.end(), NewPair) == EdgePairs.end()
                && std::find(EdgePairs.begin(), EdgePairs.end(), ReversedPair) == EdgePairs.end())
            {
                EdgePairs.push_back(NewPair);
            }
        }

        #ifdef DEBUG_REFINE
        std::cout << "Edge pairs to be merged:" << std::endl;
        for (std::size_t i = 0; i < EdgePairs.size(); ++i)
        {
            std::cout << " " << *(EdgePairs[i].first) << " is to be merged with" << std::endl
                      << "  " << *(EdgePairs[i].second) << std::endl;
        }
        std::cout << " ----------------" << std::endl;
        #endif // DEBUG_REFINE

        double alpha = 0.1;
        if (this->Has(MERGE_PARAMETER))
            alpha = this->GetValue(MERGE_PARAMETER);
        #ifdef DEBUG_REFINE
        KRATOS_WATCH(alpha)
        #endif // DEBUG_REFINE

        VertexContainerType pNewVertices;
        EdgeContainerType pNewEdges;
        EdgeContainerType pNewCompositeEdges;
        for (std::size_t i = 0; i < EdgePairs.size(); ++i)
        {
            bool is_composite = EdgePairs[i].second->IsComposite();

            MergeEdges(EdgePairs[i].first, EdgePairs[i].second, *mpVertexList, pNewVertices, mComposingVertices, pNewEdges, pNewCompositeEdges, mLastVertexId, alpha);

            // if the opposite edge is turned composite, we have to move them to other list
            if (!is_composite && EdgePairs[i].second->IsComposite())
            {
                pNewCompositeEdges.insert(pNewCompositeEdges.begin(), EdgePairs[i].second);
                mpEdgeList->erase(EdgeGetKeyType()(*(EdgePairs[i].second)));
            }
        }

        #ifdef DEBUG_REFINE
        std::cout << "  Vertices after merging:" << std::endl;
        for (VertexContainerType::const_iterator it = mpVertexList->begin(); it != mpVertexList->end(); ++it)
        {
            std::cout << "    " << *it << std::endl;
        }
        #endif // DEBUG_REFINE

        // clear the temporary memory
        EdgePairs.clear();

        // add the new edges from merge operation
        mNewEdgeList.insert(pNewEdges.ptr_begin(), pNewEdges.ptr_end());

        // remove the new edges which turns composite
        for (EdgeContainerType::iterator it = pNewCompositeEdges.begin(); it != pNewCompositeEdges.end(); ++it)
        {
            mNewEdgeList.erase(EdgeGetKeyType()(*it));
        }

        // find the zero edge and try to re-set the prev/next edge
        ReconnectEdges(mNewEdgeList);

        // remove zero and duplicated
        RemoveZeroAndDuplidateEdges(mNewEdgeList);

        // find the zero edge and try to re-set the prev/next edge
        ReconnectEdges(*mpEdgeList);

        // remove zero and duplicated
        RemoveZeroAndDuplidateEdges(*mpEdgeList);

        #ifdef DEBUG_REFINE
        std::cout << "  mNewEdgeList after merging:" << std::endl;
        for (EdgeContainerType::ptr_iterator it = mNewEdgeList.ptr_begin(); it != mNewEdgeList.ptr_end(); ++it)
        {
            std::cout << "    " << *it << ": " << *(*it) << std::endl;
        }

        std::cout << "  mpEdgeList after merging:" << std::endl;
        for (EdgeContainerType::ptr_iterator it = mpEdgeList->ptr_begin(); it != mpEdgeList->ptr_end(); ++it)
        {
            EdgeGetKeyType::result_type key = EdgeGetKeyType()(*(*it));
            EdgeContainerType::const_iterator it2 = mpEdgeList->find(key);
            std::cout << "    " << *it << ": " << *(*it) << ", key: (" << key.first << " " << key.second << "), found = " << (it2!=mpEdgeList->end()) << std::endl;
        }
        #endif // DEBUG_REFINE

        // add the composite edges to list
        mCompositeEdgeList.insert(pNewCompositeEdges.ptr_begin(), pNewCompositeEdges.ptr_end());

        // add the newly created half-edges to list, ignoring the ones already existed in edge list.
        for (EdgeContainerType::ptr_iterator it = mNewEdgeList.ptr_begin(); it != mNewEdgeList.ptr_end(); ++it)
        {
            EdgeGetKeyType::result_type key = EdgeGetKeyType()(*(*it));
            EdgeContainerType::const_iterator it2 = mpEdgeList->find(key);

            if (it2 != mpEdgeList->end())
            {
                EdgeType::Pointer pEdge = (*mpEdgeList)(key);
                pEdge->pSetNextEdge((*it)->pNextEdge());
                pEdge->pSetPrevEdge((*it)->pPrevEdge());
                pEdge->pSetFace((*it)->pFace());

                (*it)->pPrevEdge()->pSetNextEdge(pEdge);
                (*it)->pNextEdge()->pSetPrevEdge(pEdge);
            }
            else
            {
                mpEdgeList->insert(mpEdgeList->begin(), *it);
            }
        }

        // finally remove the remaining half-edges that turns to composite
        for (EdgeContainerType::iterator it = pNewCompositeEdges.begin(); it != pNewCompositeEdges.end(); ++it)
        {
            mpEdgeList->erase(EdgeGetKeyType()(*it));
        }

        #ifdef DEBUG_REFINE
        std::cout << "  Composite edges after merging:" << std::endl;
        for (EdgeContainerType::iterator it = mCompositeEdgeList.begin(); it != mCompositeEdgeList.end(); ++it)
        {
            std::cout << "    " << *it << std::endl;
        }

        std::cout << "  Edges after merging (final):" << std::endl;
        for (EdgeContainerType::ptr_iterator it = mpEdgeList->ptr_begin(); it != mpEdgeList->ptr_end(); ++it)
        {
            std::cout << "    " << *it << ": " << *(*it) << std::endl;
        }
        #endif // DEBUG_REFINE

        // this operation is normally not required. We put it here to stop the code if something wrong happens
        #ifdef DEBUG_REFINE
        std::cout << "  Unique operations are called" << std::endl;
        #endif // DEBUG_REFINE
        mpEdgeList->Unique();
        mpVertexList->Unique();
        #ifdef DEBUG_REFINE
        std::cout << "  Unique operations completed" << std::endl;
        #endif // DEBUG_REFINE

        // for all node, if the half-edge at the node is composite, then it shall have new half-edge
        for (VertexContainerType::iterator it = mpVertexList->begin(); it != mpVertexList->end(); ++it)
        {
            if (it->pEdge() != NULL) // TODO shall we check that? If we don't check, we have to remove all residual vertices
            {
                if (it->pEdge()->IsComposite())
                {
                    EdgeType::Pointer pEdge = it->pEdge()->pSubEdge(0);
                    #ifdef DEBUG_REFINE
                    std::cout << *it << " is set with ";
                    #endif // DEBUG_REFINE
                    it->pSetEdge(pEdge);
                    #ifdef DEBUG_REFINE
                    std::cout << *pEdge << std::endl;
                    #endif // DEBUG_REFINE
                }
            }
        }

        // for all refined face, now it shall not have any half-edge
        // in addition, if the face has no edge, it must not be active
        for (FaceContainerType::iterator it = mpFaceList->begin(); it != mpFaceList->end(); ++it)
        {
            if (!it->IsActive())
            {
                it->pSetEdge(NULL);
            }
        }

        // set the edge for the new vertices, because edge of the newly created vertices is removed later on
        for (VertexContainerType::iterator it = pNewVertices.begin(); it != pNewVertices.end(); ++it)
        {
            #ifdef DEBUG_REFINE
            std::cout << "     setting edge for " << *it << std::endl;
            #endif // DEBUG_REFINE
            EdgeType::Pointer pEdge = it->pEdge();
            EdgeGetKeyType::result_type key = EdgeGetKeyType()(*pEdge);
            EdgeContainerType::const_iterator it2 = mpEdgeList->find(key);

            if (it2 != mpEdgeList->end())
            {
                it->pSetEdge((*mpEdgeList)(key));
            }
            else
            {
                std::stringstream ss;
                ss << "Can't find edge with key (" << key.first << " " << key.second << ")";
                KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
            }

            // loop around each new vertex and say that the face is changed
            EdgeType::Pointer pThisEdge = pEdge;
            do
            {
                if (pThisEdge->pFace()->State() != FaceType::NEW_BORN)
                    pThisEdge->pFace()->SetState(FaceType::CHANGED);

                if (pThisEdge->pPrevEdge()->pOppositeEdge() == NULL)
                    break;
                else
                    pThisEdge = pThisEdge->pPrevEdge()->pOppositeEdge();
            } while (pThisEdge != pEdge);

            if (pEdge->pOppositeEdge() == NULL)
                continue;

            pThisEdge = pEdge->pOppositeEdge();
            do
            {
                if (pThisEdge->pFace()->State() != FaceType::NEW_BORN)
                    pThisEdge->pFace()->SetState(FaceType::CHANGED);

                if (pThisEdge->pNextEdge()->pOppositeEdge() == NULL)
                    break;
                else
                    pThisEdge = pThisEdge->pNextEdge()->pOppositeEdge();
            } while (pThisEdge != pEdge->pOppositeEdge());
        }

        // remove the composite edges
        #ifdef DEBUG_REFINE
        std::cout << "  mCompositeEdgeList is going to be cleared" << std::endl;
        #endif // DEBUG_REFINE
        mCompositeEdgeList.clear();

        // remove the newly created edge list to save memory
        #ifdef DEBUG_REFINE
        std::cout << "  mNewEdgeList is going to be cleared" << std::endl;
        #endif // DEBUG_REFINE
        mNewEdgeList.clear();

        // mark all the edges that belongs to the inactive face and remove it
        RemoveLoneEdges(*mpEdgeList);

        #ifdef DEBUG_REFINE
        std::cout << "  Faces after merging:" << std::endl;
        for (FaceContainerType::const_iterator it = mpFaceList->begin(); it != mpFaceList->end(); ++it)
        {
            std::cout << "    Face " << it->Id() << " edges:" << std::endl;
            EdgeType::Pointer pFirstEdge = it->pEdge();
            if (pFirstEdge != NULL)
            {
                EdgeType::Pointer pEdge = pFirstEdge;
                do
                {
                    std::cout << "      " << pEdge << ": " << *pEdge << std::endl;
                    pEdge = pEdge->pNextEdge();
                } while (pEdge != pFirstEdge);
            }
            // std::cout << "     " << *it << std::endl;
        }
        #endif // DEBUG_REFINE

        /*************FOR COARSENING***********/
        // TODO adjust vertex location to have all convex polygon. It happens that after the coarsening some polygons turn to be concave.

        m_half_edge_state = FINALIZED;
    }
    std::cout << __FUNCTION__ << " completed" << std::endl;
}

PolyTree2D::VertexType::Pointer PolyTree2D::RefineEdge(PolyTree2D::EdgeType::Pointer pEdge) const
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

void PolyTree2D::MergeEdges(PolyTree2D::EdgeType::Pointer pEdge,
        PolyTree2D::EdgeType::Pointer pOppositeEdge,
        PolyTree2D::VertexContainerType& rVertexList,
        PolyTree2D::VertexContainerType& rNewVertexList,
        std::map<std::size_t, std::vector<std::size_t> >& rComposingVertices,
        PolyTree2D::EdgeContainerType& rEdgeList,
        PolyTree2D::EdgeContainerType& rCompositeEdgeList,
        std::size_t& LastVertexId,
        const double& alpha) const
{
    #ifdef DEBUG_REFINE
    std::cout << "  MergeEdges starts on " << *pEdge << std::endl
              << "   and " << *pOppositeEdge << std::endl;
    #endif // DEBUG_REFINE

    if (pEdge->IsComposite())
    {
        if (pOppositeEdge->IsComposite())
        {
            #ifdef DEBUG_REFINE
            std::cout << "    Merge " << *pEdge << " with its opposite" << std::endl;
            #endif // DEBUG_REFINE

            // collect all the points on two edges
            std::map<std::size_t, double> map_distance;
            std::map<std::size_t, VertexType::Pointer> map_vertex;
            VertexType::Pointer pFirstVertex = pEdge->pSubEdge(0)->pNode1();
            std::vector<std::size_t> edge_indices;
            std::vector<std::size_t> opposite_edge_indices;
            for (std::size_t i = 0; i < pEdge->NumberOfSubEdges(); ++i)
            {
                VertexType::Pointer pVertex1 = pEdge->pSubEdge(i)->pNode1();
                double dist = sqrt(pow((*pVertex1)[0] - (*pFirstVertex)[0], 2) + pow((*pVertex1)[1] - (*pFirstVertex)[1], 2));
                map_distance[pVertex1->Id()] = dist;
                map_vertex[pVertex1->Id()] = pVertex1;
                edge_indices.push_back(pVertex1->Id());
                if (i == pEdge->NumberOfSubEdges()-1)
                    edge_indices.push_back(pEdge->pSubEdge(i)->pNode2()->Id());
            }

            for (std::size_t i = 0; i < pOppositeEdge->NumberOfSubEdges(); ++i)
            {
                VertexType::Pointer pVertex1 = pOppositeEdge->pSubEdge(i)->pNode1();
                double dist = sqrt(pow((*pVertex1)[0] - (*pFirstVertex)[0], 2) + pow((*pVertex1)[1] - (*pFirstVertex)[1], 2));
                map_distance[pVertex1->Id()] = dist;
                map_vertex[pVertex1->Id()] = pVertex1;
                opposite_edge_indices.push_back(pVertex1->Id());
                if (i == pOppositeEdge->NumberOfSubEdges())
                    opposite_edge_indices.push_back(pOppositeEdge->pSubEdge(i)->pNode2()->Id());
            }

            // sort the distance
            std::vector<double> dist_list;
            std::vector<std::size_t> id_list;
            for (std::map<std::size_t, double>::iterator it = map_distance.begin(); it != map_distance.end(); ++it)
            {
                dist_list.push_back(it->second);
                id_list.push_back(it->first);
            }

            // find the cluster of the length
            std::vector<std::vector<std::size_t> > cluster;
            std::vector<std::size_t> sorted_indices;
            PolyTreeUtility::ClusterLengths(cluster, sorted_indices, dist_list, alpha);

            for (std::size_t i = 0; i < cluster.size(); ++i)
            {
                for (std::size_t j = 0; j < cluster[i].size(); ++j)
                    cluster[i][j] = id_list[cluster[i][j]];
            }

            #ifdef DEBUG_REFINE
            std::cout << "    sub-edges on forward edge:" << std::endl;
            for (std::size_t i = 0; i < pEdge->NumberOfSubEdges(); ++i)
                std::cout << "     " << *pEdge->pSubEdge(i) << std::endl;
            std::cout << "    sub-edges on reversed edge:" << std::endl;
            for (std::size_t i = 0; i < pOppositeEdge->NumberOfSubEdges(); ++i)
                std::cout << "     " << *pOppositeEdge->pSubEdge(i) << std::endl;
            #endif // DEBUG_REFINE

            #ifdef DEBUG_REFINE
            std::cout << "     cluster:" << std::endl;
            for (std::size_t i = 0; i < cluster.size(); ++i)
            {
                std::cout << "      " << (i+1) << ":";
                for (std::size_t j = 0; j < cluster[i].size(); ++j)
                    std::cout << " " << cluster[i][j];
                std::cout << std::endl;
            }
            #endif // DEBUG_REFINE

            // find the indices of the length
            for (std::size_t i = 0; i < sorted_indices.size(); ++i)
                sorted_indices[i] = id_list[sorted_indices[i]];

            #ifdef DEBUG_REFINE
            std::cout << "     sorted indices:";
            for (std::size_t i = 0; i < sorted_indices.size(); ++i)
                std::cout << " " << sorted_indices[i];
            std::cout << std::endl;
            #endif // DEBUG_REFINE

            // create a list of potential forward edges and reversed edges
            std::vector<std::pair<std::size_t, std::size_t> > forward_edges;
            std::vector<std::pair<std::size_t, std::size_t> > reversed_edges;

            for (std::size_t i = 0, j = sorted_indices.size(); i < sorted_indices.size()-1, j > 1; ++i, --j)
            {
                forward_edges.push_back(std::pair<std::size_t, std::size_t>(sorted_indices[i], sorted_indices[i+1]));
                reversed_edges.push_back(std::pair<std::size_t, std::size_t>(sorted_indices[j-1], sorted_indices[j-2]));
            }

            #ifdef DEBUG_REFINE
            std::cout << "     forward_edges:";
            for (std::size_t i = 0; i < forward_edges.size(); ++i)
                std::cout << " (" << forward_edges[i].first << " " << forward_edges[i].second << ")";
            std::cout << std::endl;
            std::cout << "     reversed_edges:";
            for (std::size_t i = 0; i < reversed_edges.size(); ++i)
                std::cout << " (" << reversed_edges[i].first << " " << reversed_edges[i].second << ")";
            std::cout << std::endl;
            #endif // DEBUG_REFINE

            // find out which edge belong to which original edge
            std::vector<std::size_t> which_edge;
            PolyTreeUtility::FindEdge(which_edge, sorted_indices, edge_indices);

            #ifdef DEBUG_REFINE
            std::cout << "     forward_edges original edge:";
            for (std::size_t i = 0; i < which_edge.size(); ++i)
                std::cout << " " << i << ":" << which_edge[i];
            std::cout << std::endl;
            #endif // DEBUG_REFINE

            std::vector<std::size_t> which_reversed_edge;
            std::vector<std::size_t> reversed_sorted_indices = sorted_indices;
            std::reverse(reversed_sorted_indices.begin(), reversed_sorted_indices.end());
            PolyTreeUtility::FindEdge(which_reversed_edge, reversed_sorted_indices, opposite_edge_indices);

            #ifdef DEBUG_REFINE
            std::cout << "     reversed_edges original edge:";
            for (std::size_t i = 0; i < which_reversed_edge.size(); ++i)
                std::cout << " " << i << ":" << which_reversed_edge[i];
            std::cout << std::endl;
            #endif // DEBUG_REFINE

            // for each point cluster, choose a point, or create the new point
            for (std::size_t i = 0; i < cluster.size(); ++i)
            {
                VertexType::Pointer pVertex;
                if (cluster[i].size()%2 == 0)
                {
                    // create a new point
                    double x = 0.0, y = 0.0;
                    for (std::size_t j = 0; j < cluster[i].size(); ++j)
                    {
                        VertexType::Pointer pTmpVertex = rVertexList(cluster[i][j]);
                        x += (*pTmpVertex)[0];
                        y += (*pTmpVertex)[1];
                    }
                    x /= cluster[i].size();
                    y /= cluster[i].size();
                    pVertex = boost::make_shared<VertexType>(++LastVertexId, x, y);
                    rComposingVertices[pVertex->Id()] = cluster[i];
                    #ifdef DEBUG_REFINE
                    std::cout << "     New " << *pVertex << " is created" << std::endl;
                    #endif // DEBUG_REFINE
                    rVertexList.insert(rVertexList.begin(), pVertex);
                }
                else
                {
                    // pick the middle point
                    std::size_t point_index = (cluster[i].size()-1)/2;
                    pVertex = rVertexList(cluster[i][point_index]);
                }

                for (std::size_t j = 0; j < cluster[i].size(); ++j)
                {
                    for (std::size_t k = 0; k < forward_edges.size(); ++k)
                    {
                        if (forward_edges[k].first == cluster[i][j])
                        {
                            forward_edges[k].first = pVertex->Id();
                        }
                        if (forward_edges[k].second == cluster[i][j])
                        {
                            forward_edges[k].second = pVertex->Id();
                        }
                    }
                    for (std::size_t k = 0; k < reversed_edges.size(); ++k)
                    {
                        if (reversed_edges[k].first == cluster[i][j])
                        {
                            reversed_edges[k].first = pVertex->Id();
                        }
                        if (reversed_edges[k].second == cluster[i][j])
                        {
                            reversed_edges[k].second = pVertex->Id();
                        }
                    }
                }

                // add to the list of new vertices
                rNewVertexList.insert(rNewVertexList.begin(), pVertex);

                // loop half-edges around the vertex and set the new vertex
                bool vertex_is_set = false;
                EdgeType::Pointer pVertexEdge;
                for (std::size_t j = 0; j < cluster[i].size(); ++j)
                {
                    VertexType::Pointer pTmpVertex = rVertexList(cluster[i][j]);
                    EdgeType::Pointer pFirstEdge = pTmpVertex->pEdge();
                    pTmpVertex->pSetEdge(NULL);

                    EdgeType::Pointer pThisEdge = pFirstEdge;
                    do
                    {
                        pThisEdge->pSetNode1(pVertex);
                        pThisEdge->pPrevEdge()->pSetNode2(pVertex);

                        // set the edge for the vertex to an edge that is not on this line
                        if (std::find(sorted_indices.begin(), sorted_indices.end(), pThisEdge->Node2().Id()) == sorted_indices.end()
                            && !vertex_is_set)
                        {
                            pVertexEdge = pThisEdge;
                            vertex_is_set = true;
                        }

                        if (pThisEdge->pPrevEdge()->pOppositeEdge() == NULL)
                            break;
                        else
                            pThisEdge = pThisEdge->pPrevEdge()->pOppositeEdge();
                    } while (pThisEdge != pFirstEdge);

                    if (pFirstEdge->pOppositeEdge() == NULL)
                        continue;

                    pThisEdge = pFirstEdge->pOppositeEdge();
                    do
                    {
                        pThisEdge->pSetNode2(pVertex);
                        pThisEdge->pNextEdge()->pSetNode1(pVertex);

                        // set the edge for the vertex to an edge that is not on this line
                        if (std::find(sorted_indices.begin(), sorted_indices.end(), pThisEdge->pNextEdge()->Node2().Id()) == sorted_indices.end()
                            && !vertex_is_set)
                        {
                            pVertexEdge = pThisEdge->pNextEdge();
                            vertex_is_set = true;
                        }

                        if (pThisEdge->pNextEdge()->pOppositeEdge() == NULL)
                            break;
                        else
                            pThisEdge = pThisEdge->pNextEdge()->pOppositeEdge();
                    } while (pThisEdge != pFirstEdge->pOppositeEdge());
                }

                if (vertex_is_set)
                {
                    #ifdef DEBUG_REFINE
                    std::cout << "     " << *pVertex << " is set with edge ";
                    #endif // DEBUG_REFINE
                    pVertex->pSetEdge(pVertexEdge);
                    #ifdef DEBUG_REFINE
                    std::cout << *pVertexEdge << std::endl;
                    #endif // DEBUG_REFINE
                }
            }

            #ifdef DEBUG_REFINE
            std::cout << "    sub-edges on forward edge after merging:" << std::endl;
            for (std::size_t i = 0; i < pEdge->NumberOfSubEdges(); ++i)
                std::cout << "     " << *pEdge->pSubEdge(i) << std::endl;
            std::cout << "    sub-edges on reversed edge after merging:" << std::endl;
            for (std::size_t i = 0; i < pOppositeEdge->NumberOfSubEdges(); ++i)
                std::cout << "     " << *pOppositeEdge->pSubEdge(i) << std::endl;
            #endif // DEBUG_REFINE

            #ifdef DEBUG_REFINE
            std::cout << "     forward_edges after merging:";
            for (std::size_t i = 0; i < forward_edges.size(); ++i)
                std::cout << " (" << forward_edges[i].first << " " << forward_edges[i].second << ")";
            std::cout << std::endl;
            std::cout << "     reversed_edges after merging:";
            for (std::size_t i = 0; i < reversed_edges.size(); ++i)
                std::cout << " (" << reversed_edges[i].first << " " << reversed_edges[i].second << ")";
            std::cout << std::endl;
            #endif // DEBUG_REFINE

            // reset the forward edges and reversed edges if needed
            std::vector<std::pair<std::size_t, std::size_t> > new_forward_edges;
            std::vector<std::size_t> new_which_edge;
            for (std::size_t i = 0; i < forward_edges.size(); ++i)
            {
                if (forward_edges[i].first != forward_edges[i].second)
                {
                    new_forward_edges.push_back(forward_edges[i]);
                    new_which_edge.push_back(which_edge[i]);
                }
            }

            std::vector<std::pair<std::size_t, std::size_t> > new_reversed_edges;
            std::vector<std::size_t> new_which_reversed_edge;
            for (std::size_t i = 0; i < reversed_edges.size(); ++i)
            {
                if (reversed_edges[i].first != reversed_edges[i].second)
                {
                    new_reversed_edges.push_back(reversed_edges[i]);
                    new_which_reversed_edge.push_back(which_reversed_edge[i]);
                }
            }

            #ifdef DEBUG_REFINE
            std::cout << "     new_forward_edges:";
            for (std::size_t i = 0; i < new_forward_edges.size(); ++i)
                std::cout << " (" << new_forward_edges[i].first << " " << new_forward_edges[i].second << ")(" << new_which_edge[i] << ")";
            std::cout << std::endl;
            std::cout << "     new_reversed_edges:";
            for (std::size_t i = 0; i < new_reversed_edges.size(); ++i)
                std::cout << " (" << new_reversed_edges[i].first << " " << new_reversed_edges[i].second << ")(" << new_which_reversed_edge[i] << ")";
            std::cout << std::endl;
            #endif // DEBUG_REFINE

            // turn sub-edges to composite edges if needed (on forward side)
            std::vector<std::size_t> number_of_sub_edges_on_sub_edge(pEdge->NumberOfSubEdges());
            std::fill(number_of_sub_edges_on_sub_edge.begin(), number_of_sub_edges_on_sub_edge.end(), 0);
            for (std::size_t i = 0; i < new_which_edge.size(); ++i)
            {
                ++number_of_sub_edges_on_sub_edge[new_which_edge[i]];
            }

            #ifdef DEBUG_REFINE
            std::cout << "     number_of_sub_edges_on_sub_edge:";
            for (std::size_t i = 0; i < number_of_sub_edges_on_sub_edge.size(); ++i)
                std::cout << " " << number_of_sub_edges_on_sub_edge[i];
            std::cout << std::endl;
            #endif // DEBUG_REFINE

            std::vector<EdgeType::Pointer> pForwardEdges;
            for (std::size_t i = 0; i < pEdge->NumberOfSubEdges(); ++i)
            {
                if (number_of_sub_edges_on_sub_edge[i] > 1)
                {
                    // create new edges
                    std::vector<EdgeType::Pointer> pNewSubEdges;
                    for (std::size_t j = 0; j < new_forward_edges.size(); ++j)
                    {
                        if (new_which_edge[j] == i)
                        {
                            EdgeType::Pointer pNewEdge = boost::make_shared<EdgeType>(rVertexList(new_forward_edges[j].first), rVertexList(new_forward_edges[j].second));
                            pNewEdge->pSetFace(pEdge->pSubEdge(i)->pFace());
                            pNewSubEdges.push_back(pNewEdge);
                            rEdgeList.insert(rEdgeList.begin(), pNewEdge);
                            pForwardEdges.push_back(pNewEdge);
                        }
                    }

                    // set prev/next for new sub-edges
                    for (std::size_t j = 0; j < pNewSubEdges.size(); ++j)
                    {
                        if (j < pNewSubEdges.size()-1)
                        {
                            pNewSubEdges[j]->pSetNextEdge(pNewSubEdges[j+1]);
                        }
                        else
                        {
                            pNewSubEdges[j]->pSetNextEdge(pEdge->pSubEdge(i)->pNextEdge());
                            pEdge->pSubEdge(i)->pNextEdge()->pSetPrevEdge(pNewSubEdges[j]);
                        }

                        if (j > 0)
                        {
                            pNewSubEdges[j]->pSetPrevEdge(pNewSubEdges[j-1]);
                        }
                        else
                        {
                            pNewSubEdges[j]->pSetPrevEdge(pEdge->pSubEdge(i)->pPrevEdge());
                            pEdge->pSubEdge(i)->pPrevEdge()->pSetNextEdge(pNewSubEdges[j]);
                        }
                    }

                    // change the edge of the face if the face has edge pointer to the sub-edge
                    if (pEdge->pSubEdge(i)->pFace()->pEdge() == pEdge->pSubEdge(i))
                    {
                        pEdge->pSubEdge(i)->pFace()->pSetEdge(pNewSubEdges[0]);
                    }

                    // change the sub-edge to composite edge, so they will be deleted later on
                    #ifdef DEBUG_REFINE
                    std::cout << "     " << *pEdge->pSubEdge(i) << " turns to composite:";
                    #endif // DEBUG_REFINE
                    pEdge->pSubEdge(i)->SetSubEdges(pNewSubEdges);
                    rCompositeEdgeList.insert(rCompositeEdgeList.begin(), pEdge->pSubEdge(i));
                    #ifdef DEBUG_REFINE
                    std::cout << " " << *pEdge->pSubEdge(i) << std::endl;
                    #endif // DEBUG_REFINE

                    #ifdef DEBUG_REFINE
                    for (std::size_t i = 0; i < pNewSubEdges.size(); ++i)
                        std::cout << "     " << *pNewSubEdges[i] << " is created" << std::endl;
                    #endif // DEBUG_REFINE
                }
                else
                {
                    pForwardEdges.push_back(pEdge->pSubEdge(i));
                }
            }

            // turn sub-edges to composite edges if needed (on reversed side)
            std::vector<std::size_t> number_of_sub_edges_on_sub_edge_reversed(pOppositeEdge->NumberOfSubEdges());
            std::fill(number_of_sub_edges_on_sub_edge_reversed.begin(), number_of_sub_edges_on_sub_edge_reversed.end(), 0);
            for (std::size_t i = 0; i < new_which_reversed_edge.size(); ++i)
            {
                ++number_of_sub_edges_on_sub_edge_reversed[new_which_reversed_edge[i]];
            }

            #ifdef DEBUG_REFINE
            std::cout << "     number_of_sub_edges_on_sub_edge_reversed:";
            for (std::size_t i = 0; i < number_of_sub_edges_on_sub_edge_reversed.size(); ++i)
                std::cout << " " << number_of_sub_edges_on_sub_edge_reversed[i];
            std::cout << std::endl;
            #endif // DEBUG_REFINE

            std::vector<EdgeType::Pointer> pReverseEdges;
            for (std::size_t i = 0; i < pOppositeEdge->NumberOfSubEdges(); ++i)
            {
                if (number_of_sub_edges_on_sub_edge_reversed[i] > 1)
                {
                    // create new edges
                    std::vector<EdgeType::Pointer> pNewSubEdges;
                    for (std::size_t j = 0; j < new_reversed_edges.size(); ++j)
                    {
                        if (new_which_reversed_edge[j] == i)
                        {
                            EdgeType::Pointer pNewEdge = boost::make_shared<EdgeType>(rVertexList(new_reversed_edges[j].first), rVertexList(new_reversed_edges[j].second));
                            pNewEdge->pSetFace(pOppositeEdge->pSubEdge(i)->pFace());
                            pNewSubEdges.push_back(pNewEdge);
                            rEdgeList.insert(rEdgeList.begin(), pNewEdge);
                            pReverseEdges.push_back(pNewEdge);
                        }
                    }

                    // set prev/next for new reversed sub-edges
                    for (std::size_t j = 0; j < pNewSubEdges.size(); ++j)
                    {
                        if (j < pNewSubEdges.size()-1)
                        {
                            pNewSubEdges[j]->pSetNextEdge(pNewSubEdges[j+1]);
                        }
                        else
                        {
                            pNewSubEdges[j]->pSetNextEdge(pOppositeEdge->pSubEdge(i)->pNextEdge());
                            pOppositeEdge->pSubEdge(i)->pNextEdge()->pSetPrevEdge(pNewSubEdges[j]);
                        }

                        if (j > 0)
                        {
                            pNewSubEdges[j]->pSetPrevEdge(pNewSubEdges[j-1]);
                        }
                        else
                        {
                            pNewSubEdges[j]->pSetPrevEdge(pOppositeEdge->pSubEdge(i)->pPrevEdge());
                            pOppositeEdge->pSubEdge(i)->pPrevEdge()->pSetNextEdge(pNewSubEdges[j]);
                        }
                    }

                    // change the edge of the face if the face has edge pointer to the sub-edge
                    if (pOppositeEdge->pSubEdge(i)->pFace()->pEdge() == pOppositeEdge->pSubEdge(i))
                    {
                        pOppositeEdge->pSubEdge(i)->pFace()->pSetEdge(pNewSubEdges[0]);
                    }

                    // change the sub-edge to composite edge, so they will be deleted later on
                    #ifdef DEBUG_REFINE
                    std::cout << "     " << *pOppositeEdge->pSubEdge(i) << " turns to composite:";
                    #endif // DEBUG_REFINE
                    pOppositeEdge->pSubEdge(i)->SetSubEdges(pNewSubEdges);
                    rCompositeEdgeList.insert(rCompositeEdgeList.begin(), pOppositeEdge->pSubEdge(i));
                    #ifdef DEBUG_REFINE
                    std::cout << " " << *pOppositeEdge->pSubEdge(i) << std::endl;
                    #endif // DEBUG_REFINE

                    #ifdef DEBUG_REFINE
                    for (std::size_t i = 0; i < pNewSubEdges.size(); ++i)
                        std::cout << "     " << *pNewSubEdges[i] << " is created" << std::endl;
                    #endif // DEBUG_REFINE
                }
                else
                {
                    pReverseEdges.push_back(pOppositeEdge->pSubEdge(i));
                }
            }

            #ifdef DEBUG_REFINE
            std::cout << "     list of forward edges before setting opposite (only half-edges, no composite):" << std::endl;
            for (std::size_t i = 0; i < pForwardEdges.size(); ++i)
                std::cout << "      " << *pForwardEdges[i] << std::endl;
            std::cout << "     list of reversed edges before setting opposite (only half-edges, no composite):" << std::endl;
            for (std::size_t i = 0; i < pReverseEdges.size(); ++i)
                std::cout << "      " << *pReverseEdges[i] << std::endl;
            #endif // DEBUG_REFINE

            // set opposite edge again
            // we also modify the forward edges a bit to use the existing edges if having
            for (std::size_t i = 0, j = pReverseEdges.size(); i < pForwardEdges.size()-1, j > 0; ++i, --j)
            {
                if (pForwardEdges[i]->pOppositeEdge() == NULL)
                {
                    pForwardEdges[i]->pSetOppositeEdge(pReverseEdges[j-1]);
                }
                else
                {
                    if (pReverseEdges[j-1]->pFace()->pEdge() == pReverseEdges[j-1])
                        pReverseEdges[j-1]->pFace()->pSetEdge(pForwardEdges[i]->pOppositeEdge());
                    pReverseEdges[j-1] = pForwardEdges[i]->pOppositeEdge();
                }

                if (pReverseEdges[j-1]->pOppositeEdge() == NULL)
                {
                    pReverseEdges[j-1]->pSetOppositeEdge(pForwardEdges[i]);
                }
                else
                {
                    if (pForwardEdges[i]->pFace()->pEdge() == pForwardEdges[i])
                        pForwardEdges[i]->pFace()->pSetEdge(pReverseEdges[j-1]->pOppositeEdge());
                    pForwardEdges[i] = pReverseEdges[j-1]->pOppositeEdge();
                }
            }

            #ifdef DEBUG_REFINE
            std::cout << "     list of forward edges (only half-edges, no composite):" << std::endl;
            for (std::size_t i = 0; i < pForwardEdges.size(); ++i)
                std::cout << "      " << *pForwardEdges[i] << std::endl;
            std::cout << "     list of reversed edges (only half-edges, no composite):" << std::endl;
            for (std::size_t i = 0; i < pReverseEdges.size(); ++i)
                std::cout << "      " << *pReverseEdges[i] << std::endl;
            #endif // DEBUG_REFINE

            #ifdef DEBUG_REFINE
            std::cout << "     The edge is " << *pEdge << std::endl;
            std::cout << "     The opposite edge is " << *pOppositeEdge << std::endl;
            #endif // DEBUG_REFINE

            // turn the original sub-edges to composite if needed
            for (std::size_t i = 0; i < pEdge->NumberOfSubEdges(); ++i)
            {
                if (pEdge->pSubEdge(i)->pOppositeEdge() != NULL)
                {
                    if (pEdge->pSubEdge(i)->IsComposite() && !pEdge->pSubEdge(i)->pOppositeEdge()->IsComposite())
                    {
                        std::vector<EdgeType::Pointer> pTmpEdges;
                        for (std::size_t j = pEdge->pSubEdge(i)->NumberOfSubEdges(); j > 0; --j)
                        {
                            pTmpEdges.push_back(pEdge->pSubEdge(i)->pSubEdge(j-1)->pOppositeEdge());
                        }
                        #ifdef DEBUG_REFINE
                        std::cout << "     " << *pEdge->pSubEdge(i)->pOppositeEdge() << " turns to composite:";
                        #endif // DEBUG_REFINE
                        pEdge->pSubEdge(i)->pOppositeEdge()->SetSubEdges(pTmpEdges);
                        #ifdef DEBUG_REFINE
                        std::cout << " " << *pEdge->pSubEdge(i)->pOppositeEdge() << std::endl;
                        #endif // DEBUG_REFINE
                    }
                }
            }

            for (std::size_t i = 0; i < pOppositeEdge->NumberOfSubEdges(); ++i)
            {
                if (pOppositeEdge->pSubEdge(i)->pOppositeEdge() != NULL)
                {
                    if (pOppositeEdge->pSubEdge(i)->IsComposite() && !pOppositeEdge->pSubEdge(i)->pOppositeEdge()->IsComposite())
                    {
                        std::vector<EdgeType::Pointer> pTmpEdges;
                        for (std::size_t j = pOppositeEdge->pSubEdge(i)->NumberOfSubEdges(); j > 0; --j)
                        {
                            pTmpEdges.push_back(pOppositeEdge->pSubEdge(i)->pSubEdge(j-1)->pOppositeEdge());
                        }
                        #ifdef DEBUG_REFINE
                        std::cout << "     " << *pOppositeEdge->pSubEdge(i)->pOppositeEdge() << " turns to composite:";
                        #endif // DEBUG_REFINE
                        pOppositeEdge->pSubEdge(i)->pOppositeEdge()->SetSubEdges(pTmpEdges);
                        #ifdef DEBUG_REFINE
                        std::cout << " " << *pOppositeEdge->pSubEdge(i)->pOppositeEdge() << std::endl;
                        #endif // DEBUG_REFINE
                    }
                }
            }

            // re-set the sub-edges for forward edge
            pEdge->SetSubEdges(pForwardEdges);

            // re-set the sub-edges for the reversed edge
            pOppositeEdge->SetSubEdges(pReverseEdges);

            #ifdef DEBUG_REFINE
            std::cout << "     The edge now is " << *pEdge << std::endl;
            std::cout << "     The opposite edge now is " << *pOppositeEdge << std::endl;
            #endif // DEBUG_REFINE
        }
        else
        {
            #ifdef DEBUG_REFINE
            std::cout << "    " << *pEdge << " generates half-edges for opposite" << std::endl;
            #endif // DEBUG_REFINE

            std::vector<EdgeType::Pointer> pReverseEdges;
            for (std::size_t i = pEdge->NumberOfSubEdges(); i > 0 ; --i)
            {
                EdgeType::Pointer pSubEdge = pEdge->pSubEdge(i-1);
                VertexType::Pointer pVertex1 = pSubEdge->pNode2();
                VertexType::Pointer pVertex2 = pSubEdge->pNode1();
                EdgeType::Pointer pNewHalfEdge = boost::make_shared<EdgeType>(pVertex1, pVertex2);
                pNewHalfEdge->pSetOppositeEdge(pSubEdge);
                pNewHalfEdge->pSetFace(pOppositeEdge->pFace()); // set the face for opposite half-edge
                pSubEdge->pSetOppositeEdge(pNewHalfEdge);
                pReverseEdges.push_back(pNewHalfEdge);
                rEdgeList.insert(rEdgeList.begin(), pNewHalfEdge);
            }

            // make the opposite edge also composite edge
            pOppositeEdge->SetSubEdges(pReverseEdges);

            for (std::size_t i = 0; i < pReverseEdges.size(); ++i)
            {
                if (i < pReverseEdges.size()-1)
                {
                    pReverseEdges[i]->pSetNextEdge(pReverseEdges[i+1]);
                }
                else
                {
                    pReverseEdges[i]->pSetNextEdge(pOppositeEdge->pNextEdge());
                    pOppositeEdge->pNextEdge()->pSetPrevEdge(pReverseEdges[i]);
                }

                if (i > 0)
                {
                    pReverseEdges[i]->pSetPrevEdge(pReverseEdges[i-1]);
                }
                else
                {
                    pReverseEdges[i]->pSetPrevEdge(pOppositeEdge->pPrevEdge());
                    pOppositeEdge->pPrevEdge()->pSetNextEdge(pReverseEdges[i]);
                }
            }

            // change the state of the opposite face
            if (pOppositeEdge->pFace()->State() != FaceType::NEW_BORN)
                pOppositeEdge->pFace()->SetState(FaceType::CHANGED);

            #ifdef DEBUG_REFINE
            for (std::size_t i = 0; i < pReverseEdges.size(); ++i)
            {
                std::cout << "     " << *pReverseEdges[i] << " on reversed side is created" << std::endl;
            }
            #endif // DEBUG_REFINE
        }
    }
    else
    {
        if (pOppositeEdge->IsComposite())
        {
            #ifdef DEBUG_REFINE
            std::cout << "    Opposite generates half-edges for " << *pEdge << std::endl;
            #endif // DEBUG_REFINE

            std::vector<EdgeType::Pointer> pForwardEdges;
            for (std::size_t i = pOppositeEdge->NumberOfSubEdges(); i > 0 ; --i)
            {
                EdgeType::Pointer pSubEdge = pOppositeEdge->pSubEdge(i-1);
                VertexType::Pointer pVertex1 = pSubEdge->pNode2();
                VertexType::Pointer pVertex2 = pSubEdge->pNode1();
                EdgeType::Pointer pNewHalfEdge = boost::make_shared<EdgeType>(pVertex1, pVertex2);
                pNewHalfEdge->pSetOppositeEdge(pSubEdge);
                pNewHalfEdge->pSetFace(pEdge->pFace()); // set the face for half-edge
                pSubEdge->pSetOppositeEdge(pNewHalfEdge);
                pForwardEdges.push_back(pNewHalfEdge);
                rEdgeList.insert(rEdgeList.begin(), pNewHalfEdge);
            }

            // make this edge also composite edge
            pEdge->SetSubEdges(pForwardEdges);

            for (std::size_t i = 0; i < pForwardEdges.size(); ++i)
            {
                if (i < pForwardEdges.size()-1)
                {
                    pForwardEdges[i]->pSetNextEdge(pForwardEdges[i+1]);
                }
                else
                {
                    pForwardEdges[i]->pSetNextEdge(pEdge->pNextEdge());
                    pEdge->pNextEdge()->pSetPrevEdge(pForwardEdges[i]);
                }

                if (i > 0)
                {
                    pForwardEdges[i]->pSetPrevEdge(pForwardEdges[i-1]);
                }
                else
                {
                    pForwardEdges[i]->pSetPrevEdge(pEdge->pPrevEdge());
                    pEdge->pPrevEdge()->pSetNextEdge(pForwardEdges[i]);
                }
            }

            // change the state of the face
            if (pEdge->pFace()->State() != FaceType::NEW_BORN)
                pEdge->pFace()->SetState(FaceType::CHANGED);

            #ifdef DEBUG_REFINE
            for (std::size_t i = 0; i < pForwardEdges.size(); ++i)
            {
                std::cout << "     " << *pForwardEdges[i] << " on forward side is created" << std::endl;
            }
            #endif // DEBUG_REFINE
        }
        else
        {
            // here do nothing
        }
    }

    if (pEdge->IsComposite())
    {
        // if this face has pointer to the edge, then it shall point to the sub-edge
        if (pEdge->pFace()->pEdge() == pEdge)
        {
            pEdge->pFace()->pSetEdge(pEdge->pSubEdge(0));
        }
    }

    if (pOppositeEdge->IsComposite())
    {
        // if the opposite face has pointer to the opposite edge, then it shall point to the sub-edge
        if (pOppositeEdge->pFace()->pEdge() == pOppositeEdge)
        {
            pOppositeEdge->pFace()->pSetEdge(pOppositeEdge->pSubEdge(0));
        }
    }

    #ifdef DEBUG_REFINE
    std::cout << "  MergeEdges completed" << std::endl;
    #endif // DEBUG_REFINE
}

void PolyTree2D::RefineFace(PolyTree2D::FaceType::Pointer pFace,
        PolyTree2D::VertexContainerType& rVertexList,
        PolyTree2D::EdgeContainerType& rEdgeList,
        PolyTree2D::EdgeContainerType& rNewEdgeList,
        PolyTree2D::EdgeContainerType& rCompositeEdgeList,
        PolyTree2D::FaceContainerType& rFaceList,
        std::size_t& LastVertexId, std::size_t& LastFaceId) const
{
    std::cout << __FUNCTION__ << " " << pFace->Id() << " starts" << std::endl;
    std::cout << *pFace << " will be refined" << std::endl;

    // create the polygon data
    std::vector<std::vector<double> > polygon; // coordinates of the polygon used for Voronoi tessellation
    std::vector<VertexType::Pointer> pVertices; // list of vertices of the polygon
    std::vector<EdgeType::Pointer> pEdges; // list of edges of the polygon (can be composite edge, but the half-edge relation is not required)

    // only extract the vertices of the polygon
    pFace->ExtractMinimalPolygon(polygon, pVertices, pEdges);

    // set the opposite edges for minimal polygon edges
    for (std::size_t i = 0; i < pEdges.size(); ++i)
    {
        // set the face
        pEdges[i]->pSetFace(pFace);

        // set the opposite edge
        if (!pEdges[i]->IsComposite())
        {
            // here we don't need to set opposite edge, because the single edge shall have opposite already
            // EdgeGetKeyType::result_type key(pEdges[i]->Id1(), pEdges[i]->Id2());
        }
        else
        {
            for (std::size_t j = 0; j < pEdges[i]->NumberOfSubEdges(); ++j)
            {
                EdgeGetKeyType::result_type key(pEdges[i]->pSubEdge(j)->Id1(), pEdges[i]->pSubEdge(j)->Id2());
                EdgeGetKeyType::result_type reversed_key(key.second, key.first);
                EdgeContainerType::const_iterator it = rEdgeList.find(reversed_key);
                if (it != rEdgeList.end())
                {
                    pEdges[i]->pSubEdge(j)->pSetOppositeEdge(rEdgeList(reversed_key));
                }
                else
                    KRATOS_THROW_ERROR(std::logic_error, "The opposite edge is not found. That's strange.", "")
            }
        }
    }

    #ifdef DEBUG_REFINE
    std::cout << "polygon:" << std::endl;
    for (std::size_t i = 0; i < polygon.size(); ++i)
    {
        std::cout << "  " << i+1 << ": " << polygon[i][0] << " " << polygon[i][1] << std::endl;
    }
    std::cout << "vertices of the polygon:" << std::endl;
    for (std::size_t i = 0; i < pVertices.size(); ++i)
    {
        std::cout << "  " << i+1 << ": " << *pVertices[i] << std::endl;
    }
    std::cout << "edges of the polygon:" << std::endl;
    for (std::size_t i = 0; i < pEdges.size(); ++i)
    {
        std::cout << "  " << i+1 << ": " << *pEdges[i] << std::endl;
    }
    #endif // DEBUG_REFINE

    // extract the reversed composite edges
    for (std::size_t i = 0; i < pEdges.size(); ++i)
    {
        if (pEdges[i]->IsComposite())
        {
            typename EdgeType::Pointer pNewEdge = boost::make_shared<EdgeType>(pEdges[i]->pSubEdge(pEdges[i]->NumberOfSubEdges()-1)->pNode2(), pEdges[i]->pSubEdge(0)->pNode1());
            std::vector<EdgeType::Pointer> pSubEdges;
            for (std::size_t j = pEdges[i]->NumberOfSubEdges(); j > 0; --j)
            {
                pSubEdges.push_back(pEdges[i]->pSubEdge(j-1)->pOppositeEdge());
            }
            pNewEdge->SetSubEdges(pSubEdges);
            // std::cout << *pNewEdge << " is created" << std::endl;
            pNewEdge->pSetFace(pEdges[i]->pSubEdge(0)->pOppositeEdge()->pFace());
            rCompositeEdgeList.insert(rCompositeEdgeList.begin(), pNewEdge);
        }
    }

    // call utility to decompose polygon
    std::vector<std::vector<double> > compatible_coordinates;
    std::vector<std::vector<std::size_t> > voronoi_connectivities;
    double pro_constant = 1.5;      // TODO parameterize this
    std::size_t max_iter = 100;     // TODO parameterize this
    double tol = 1.0e-4;            // TODO parameterize this
    #ifdef DEBUG_REFINE
    std::size_t debug_level = 3;    // TODO parameterize this
    #else
    std::size_t debug_level = 0;    // TODO parameterize this
    #endif // DEBUG_REFINE
    PolyTreeUtility::ComputePolygonDecomposition(compatible_coordinates,
            voronoi_connectivities, polygon, pro_constant, max_iter, tol, debug_level);

    #ifdef DEBUG_REFINE
    std::cout << "Voronoi coordinates:" << std::endl;
    for (std::size_t i = 0; i < compatible_coordinates.size(); ++i)
    {
        std::cout << "  " << i+1 << ": " << compatible_coordinates[i][0] << " " << compatible_coordinates[i][1] << std::endl;
    }

    std::cout << "voronoi_connectivities:" << std::endl;
    for (std::size_t i = 0; i < voronoi_connectivities.size(); ++i)
    {
        std::cout << " ";
        for (std::size_t j = 0; j < voronoi_connectivities[i].size(); ++j)
        {
            std::cout << " " << voronoi_connectivities[i][j];
        }
        std::cout << std::endl;
    }
    #endif // DEBUG_REFINE

    // encode the vertices in the appropriate way
    std::vector<int> vertex_stat;
    std::vector<int> vertex_info;
    std::vector<int> edge_info;
    double detect_tol = 1.0e-6; // TODO parameterize this
    PolyTreeUtility::EncodeVertices(vertex_stat, vertex_info, edge_info, compatible_coordinates, polygon, 1, detect_tol);

    #ifdef DEBUG_REFINE
    std::cout << "vertex_stat:" << std::endl;
    for (std::size_t i = 0; i < compatible_coordinates.size(); ++i)
    {
        std::cout << "  " << (i+1) << ": " << vertex_stat[i];

        if (vertex_stat[i] == 0)
            std::cout << ", is inner";
        else if (vertex_stat[i] == 1)
            std::cout << ", on vertex " << vertex_info[i];
        else if (vertex_stat[i] == 2)
            std::cout << ", on edge " << vertex_info[i] << ", index " << edge_info[i];

        std::cout << std::endl;
    }
    #endif // DEBUG_REFINE

    // create new vertices
    std::map<std::size_t, VertexType::Pointer> VertexMap; // this is the map from Voronoi vertices (compatible_coordinates) to the new VertexType in the tree
    for (std::size_t i = 0; i < compatible_coordinates.size(); ++i)
    {
        if (vertex_stat[i] == 0 || vertex_stat[i] == 2)
        {
            VertexType::Pointer pNewVertex = boost::make_shared<VertexType>(++LastVertexId, compatible_coordinates[i][0], compatible_coordinates[i][1]);
            #ifdef DEBUG_REFINE
            std::cout << "Vertex " << pNewVertex->Id() << " is created at " << (*pNewVertex)[0] << "," << (*pNewVertex)[1] << std::endl;
            #endif // DEBUG_REFINE
            VertexMap[i] = pNewVertex;
            rVertexList.insert(rVertexList.begin(), pNewVertex);
        }
        else if (vertex_stat[i] == 1)
        {
            VertexMap[i] = pVertices[vertex_info[i]];
        }
        else if (vertex_stat[i] == 3)
        {
            KRATOS_THROW_ERROR(std::logic_error, "The vertex lies on edge but outside. Something is wrong with the vertices", "")
        }
    }

    // create new half-edges for the internal sub-faces
    typedef std::pair<std::size_t, std::size_t> EdgeKey;
    typedef std::map<EdgeKey, EdgeType::Pointer> EdgeMapType;
    EdgeMapType EdgeMap; // this is the map for half-edge
    std::vector<FaceType::Pointer> pNewFaces;
    for (std::size_t i = 0; i < voronoi_connectivities.size(); ++i)
    {
        // create new face
        FaceType::Pointer pNewFace = boost::make_shared<FaceType>(++LastFaceId);
        rFaceList.insert(rFaceList.begin(), pNewFace);
        pNewFaces.push_back(pNewFace);
        pNewFace->pSetParent(pFace);

        // create new half-edges
        std::vector<EdgeType::Pointer> pNewEdges;
        #ifdef DEBUG_REFINE
        KRATOS_WATCH(voronoi_connectivities[i].size())
        #endif // DEBUG_REFINE
        for (std::size_t j = 0; j < voronoi_connectivities[i].size(); ++j)
        {
            std::size_t this_vertex = voronoi_connectivities[i][j]-1;
            std::size_t next_vertex = (j==voronoi_connectivities[i].size()-1) ? voronoi_connectivities[i][0]-1 : voronoi_connectivities[i][j+1]-1;
            EdgeType::Pointer pEdge = boost::make_shared<EdgeType>(VertexMap[this_vertex], VertexMap[next_vertex]);
            #ifdef DEBUG_REFINE
            std::cout << j << ": " << *pEdge << " is created" << std::endl;
            #endif // DEBUG_REFINE
            rNewEdgeList.insert(rNewEdgeList.begin(), pEdge);
            pEdge->pSetFace(pNewFace);
            if (j == 0) pNewFace->pSetEdge(pEdge);
            pNewEdges.push_back(pEdge);
            std::pair<std::size_t, std::size_t> key(VertexMap[this_vertex]->Id(), VertexMap[next_vertex]->Id());
            if (vertex_stat[this_vertex] == 0 || vertex_stat[this_vertex] == 2) VertexMap[this_vertex]->pSetEdge(pEdge);
            EdgeMap[key] = pEdge;
        }

        // set the previous and next edge
        for (std::size_t j = 0; j < pNewEdges.size(); ++j)
        {
            std::size_t prev_edge = (j==0) ? pNewEdges.size()-1 : j-1;
            std::size_t next_edge = (j==pNewEdges.size()-1) ? 0 : j+1;
            pNewEdges[j]->pSetNextEdge(pNewEdges[next_edge]);
            pNewEdges[j]->pSetPrevEdge(pNewEdges[prev_edge]);
        }
    }

    // set the sub-faces
    pFace->SetSubFaces(pNewFaces);

    #ifdef DEBUG_REFINE
    std::cout << "EdgeMap:";
    for (EdgeMapType::iterator it = EdgeMap.begin(); it != EdgeMap.end(); ++it)
    {
        std::cout << " (" << it->first.first << "-" << it->first.second << ")";
    }
    std::cout << std::endl;
    #endif // DEBUG_REFINE

    // set the opposite edges
    for (EdgeMapType::iterator it = EdgeMap.begin(); it != EdgeMap.end(); ++it)
    {
        std::size_t first_vertex = it->first.first;
        std::size_t second_vertex = it->first.second;

        EdgeKey opposite_key(second_vertex, first_vertex);
        EdgeMapType::iterator it2 = EdgeMap.find(opposite_key);
        if (it2 != EdgeMap.end())
        {
            it->second->pSetOppositeEdge(it2->second);
            #ifdef DEBUG_REFINE
            std::cout << "Half-edge (" << it->second->Id1() << "-" << it->second->Id2()
                      << ") is set with the opposite (" << it2->second->Id1() << "-" << it2->second->Id2() << ")" << std::endl;
            #endif // DEBUG_REFINE
        }
    }

    // create new composite half-edges
    std::vector<std::vector<std::size_t> > composite_edges(polygon.size());
    for (std::size_t i = 0; i < compatible_coordinates.size(); ++i)
    {
        if (vertex_stat[i] == 2)
        {
            std::size_t edge_index = static_cast<std::size_t>(vertex_info[i]);
            std::size_t edge_point_index = static_cast<std::size_t>(edge_info[i]);
            if (edge_point_index+2 > composite_edges[edge_index].size())
            {
                composite_edges[edge_index].resize(edge_point_index+2);
                composite_edges[edge_index][edge_point_index] = i;
            }
            else
                composite_edges[edge_index][edge_point_index] = i;
        }
    }
    for (std::size_t i = 0; i < compatible_coordinates.size(); ++i)
    {
        if (vertex_stat[i] == 1)
        {
            std::size_t vertex_index = static_cast<std::size_t>(vertex_info[i]);
            std::size_t this_edge = vertex_index;
            std::size_t prev_edge = (vertex_index==0) ? polygon.size()-1 : vertex_index-1;

            if (composite_edges[this_edge].size() == 0)
                composite_edges[this_edge].resize(2);
            composite_edges[this_edge][0] = i;

            if (composite_edges[prev_edge].size() == 0)
            {
                composite_edges[prev_edge].resize(2);
                composite_edges[prev_edge][1] = i;
            }
            else
            {
                composite_edges[prev_edge].back() = i;
            }
        }
    }

    if (debug_level > 3)
    {
        for (std::size_t i = 0; i < composite_edges.size(); ++i)
        {
            std::cout << "composite_edges[" << i << "]:";
            for (std::size_t j = 0; j < composite_edges[i].size(); ++j)
                std::cout << " " << VertexMap[composite_edges[i][j]]->Id();
            std::cout << std::endl;
        }
    }

    // when we refine the face, we make composite edges around that face
    for (std::size_t i = 0; i < composite_edges.size(); ++i)
    {
        if (composite_edges[i].size() > 2)
        {
            std::vector<EdgeType::Pointer> pSubEdges;
            for (std::size_t j = 0; j < composite_edges[i].size()-1; ++j)
            {
                EdgeKey key(VertexMap[composite_edges[i][j]]->Id(), VertexMap[composite_edges[i][j+1]]->Id());
                EdgeMapType::iterator it2 = EdgeMap.find(key);
                if (it2 != EdgeMap.end())
                {
                    pSubEdges.push_back(it2->second);
                }
                else
                {
                    std::stringstream ss;
                    ss << "The edge (" << key.first << "-" << key.second << ") does not exist";
                    KRATOS_THROW_ERROR(std::logic_error, ss.str(), "")
                }
            }
            // KRATOS_WATCH(pSubEdges.size())

            // turn the corresponding edge to composite
            pEdges[i]->SetSubEdges(pSubEdges);
            rCompositeEdgeList.insert(rCompositeEdgeList.begin(), pEdges[i]);
            rEdgeList.erase(EdgeGetKeyType()(*pEdges[i]));
        }
    }

    std::cout << __FUNCTION__ << " " << pFace->Id() << " completed" << std::endl;
}

void PolyTree2D::CoarsenFace(PolyTree2D::FaceType& rFace, PolyTree2D::FaceContainerType& rFaceList) const
{

}

void PolyTree2D::ReconnectEdges(PolyTree2D::EdgeContainerType& rEdgeList) const
{
    for (EdgeContainerType::ptr_const_iterator it = rEdgeList.ptr_begin(); it != rEdgeList.ptr_end(); ++it)
    {
        EdgeGetKeyType::result_type key = EdgeGetKeyType()(*(*it));
        if (key.first == key.second)
        {
            EdgeType::Pointer pPrevEdge = (*it)->pPrevEdge();
            EdgeType::Pointer pNextEdge = (*it)->pNextEdge();
            pPrevEdge->pSetNextEdge(pNextEdge);
            pNextEdge->pSetPrevEdge(pPrevEdge);
            if ((*it)->pFace()->pEdge() == *it) (*it)->pFace()->pSetEdge(pNextEdge);
        }
    }
}

void PolyTree2D::RemoveZeroAndDuplidateEdges(PolyTree2D::EdgeContainerType& rEdgeList) const
{
    rEdgeList.Unique(); // remove the duplicated edge, sort is required to ensure find working correctly.

    // remove the zero edges (the one with Id1 == Id2)
    std::set<EdgeGetKeyType::result_type> removed_keys;
    for (EdgeContainerType::const_iterator it = rEdgeList.begin(); it != rEdgeList.end(); ++it)
    {
        EdgeGetKeyType::result_type key = EdgeGetKeyType()(*it);
        if (key.first == key.second)
            removed_keys.insert(key);
    }

    for (std::set<EdgeGetKeyType::result_type>::iterator it = removed_keys.begin(); it != removed_keys.end(); ++it)
    {
        rEdgeList.erase(*it);
    }

    rEdgeList.Sort(); // sort after erasing
}

void PolyTree2D::RemoveLoneEdges(PolyTree2D::EdgeContainerType& rEdgeList) const
{
    std::set<EdgeGetKeyType::result_type> removed_keys;
    for (EdgeContainerType::const_iterator it = rEdgeList.begin(); it != rEdgeList.end(); ++it)
    {
        if (it->pFace() == NULL)
        {
            removed_keys.insert(EdgeGetKeyType()(*it));
        }
        else
        {
            if (!it->pFace()->IsActive())
                removed_keys.insert(EdgeGetKeyType()(*it));
        }
    }

    for (std::set<EdgeGetKeyType::result_type>::iterator it = removed_keys.begin(); it != removed_keys.end(); ++it)
    {
        rEdgeList.erase(*it);
    }

    rEdgeList.Sort(); // sort after erasing
}

void PolyTree2D::RemoveLoneVertices(PolyTree2D::VertexContainerType& rVertexList) const
{
    std::set<std::size_t> removed_keys;
    for (VertexContainerType::const_iterator it = rVertexList.begin(); it != rVertexList.end(); ++it)
    {
        if (it->pEdge() == NULL)
        {
            removed_keys.insert(it->Id());
        }
    }

    for (std::set<std::size_t>::iterator it = removed_keys.begin(); it != removed_keys.end(); ++it)
    {
        rVertexList.erase(*it);
    }

    rVertexList.Sort(); // sort after erasing
}

void PolyTree2D::ExtractEdges(PolyTree2D::FaceType::Pointer pFace,
        PolyTree2D::EdgeContainerType& rEdgeList) const
{
    if (pFace->IsLeaf())
    {
        EdgeType::Pointer pFirstEdge = pFace->pEdge();
        EdgeType::Pointer pEdge = pFirstEdge;

        do
        {
            rEdgeList.insert(rEdgeList.begin(), pEdge);
            pEdge = pEdge->pNextEdge();
        } while(pEdge != pFirstEdge);
    }
    else
    {
        for (std::size_t i = 0; i < pFace->NumberOfSubFaces(); ++i)
            ExtractEdges(pFace->pSubFace(i), rEdgeList);
    }
}

void PolyTree2D::ExtractFaces(PolyTree2D::FaceType::Pointer pFace,
        PolyTree2D::FaceContainerType& rFaceList) const
{
    if (pFace->IsLeaf())
    {
        rFaceList.insert(rFaceList.begin(), pFace);
    }
    else
    {
        for (std::size_t i = 0; i < pFace->NumberOfSubFaces(); ++i)
            ExtractFaces(pFace->pSubFace(i), rFaceList);
    }
}

void PolyTree2D::ExtractEdgesAndVertices(PolyTree2D::FaceType::Pointer pFace,
        PolyTree2D::EdgeContainerType& rOuterEdgeList,
        PolyTree2D::EdgeContainerType& rInnerEdgeList,
        PolyTree2D::VertexContainerType& rInnerVertexList,
        PolyTree2D::FaceContainerType& rInnerFaceList) const
{
    if (pFace->IsLeaf())
    {
        // if the face is the leaf node of the polytree, everything is simple. To get the edges, one just need to iterate around the face.
        rOuterEdgeList.clear();

        EdgeType::Pointer pFirstEdge = pFace->pEdge();
        EdgeType::Pointer pEdge = pFirstEdge;

        do
        {
            rOuterEdgeList.push_back(pEdge);
            pEdge = pEdge->pNextEdge();
        } while (pEdge != pFirstEdge);

        rInnerEdgeList.clear();
        rInnerVertexList.clear();
        rInnerFaceList.clear();
        rInnerFaceList.push_back(pFace);
    }
    else
    {
        // in the case that the face is a parent node, thing is more complicated. Firstly all the edges will be collected, and then the inner edges shall be removed.
        // firstly collect all the faces
        rInnerFaceList.clear();
        ExtractFaces(pFace, rInnerFaceList);

        // get all the edges
        EdgeContainerType pAllEdges;
        ExtractEdges(pFace, pAllEdges);
        pAllEdges.Unique();

        // extract the inner edges and outer edges
        EdgeContainerType pOuterEdges;
        rInnerEdgeList.clear();
        rInnerVertexList.clear();
        for (EdgeContainerType::ptr_iterator it = pAllEdges.ptr_begin(); it != pAllEdges.ptr_end(); ++it)
        {
            EdgeGetKeyType::result_type key = EdgeGetKeyType()(*(*it));
            EdgeGetKeyType::result_type reversed_key(key.second, key.first);

            EdgeContainerType::const_iterator it2 = pAllEdges.find(reversed_key);
            if (it2 == pAllEdges.end())
            {
                pOuterEdges.insert(pOuterEdges.begin(), *it);
            }
            else
            {
                rInnerEdgeList.insert(rInnerEdgeList.begin(), *it);
            }

            rInnerVertexList.insert(rInnerVertexList.begin(), (*it)->pNode1());
            rInnerVertexList.insert(rInnerVertexList.begin(), (*it)->pNode2());
        }

        rInnerVertexList.Unique();

        // extract the outer edge list, so that it makes a closed loop
        EdgeType::Pointer pFirstEdge = *(pOuterEdges.ptr_begin());
        EdgeType::Pointer pEdge = pFirstEdge;

        rOuterEdgeList.clear();
        do
        {
            rOuterEdgeList.push_back(pEdge);
            EdgeGetKeyType::result_type key = EdgeGetKeyType()(*pEdge);

            bool found = false;
            for (EdgeContainerType::ptr_iterator it = pOuterEdges.ptr_begin(); it != pOuterEdges.ptr_end(); ++it)
            {
                if (*it != pEdge)
                {
                    EdgeGetKeyType::result_type key2 = EdgeGetKeyType()(*(*it));
                    if (key2.first == key.second)
                    {
                        found = true;
                        pEdge = *it;
                        break;
                    }
                }
            }

            if (!found)
                KRATOS_THROW_ERROR(std::logic_error, "The next edge is not found. Something is wrong", "")
        } while(pEdge != pFirstEdge);

        // extract all the inner vertices
        for (EdgeContainerType::ptr_const_iterator it = rOuterEdgeList.ptr_begin(); it != rOuterEdgeList.ptr_end(); ++it)
        {
            rInnerVertexList.erase((*it)->pNode1()->Id());
            rInnerVertexList.erase((*it)->pNode2()->Id());
        }
        rInnerVertexList.Sort();
    }
}

void PolyTree2D::Renumber(std::map<std::size_t, std::size_t>& rMapVertices,
            std::map<std::size_t, std::size_t>& rMapFaces)
{
    mLastVertexId = 0;
    for (VertexContainerType::iterator it = mpVertexList->begin(); it != mpVertexList->end(); ++it)
    {
        ++mLastVertexId;
        rMapVertices[it->Id()] = mLastVertexId;
        it->SetId(mLastVertexId);
    }

    // TODO shall we re-order the vertices to make it more direct-solver-efficiency?

    mLastFaceId = 0;
    for (FaceContainerType::iterator it = mpFaceList->begin(); it != mpFaceList->end(); ++it)
    {
        ++mLastFaceId;
        rMapFaces[it->Id()] = mLastFaceId;
        it->SetId(mLastFaceId);
    }
}

void PolyTree2D::Validate() const
{
    // validate the vertex
    const std::size_t max_iters = 10000;
    std::size_t iter;
    #ifdef DEBUG_REFINE
    std::cout << __FUNCTION__ << " initiated" << std::endl;
    #endif // DEBUG_REFINE

    for (VertexContainerType::ptr_const_iterator it = mpVertexList->ptr_begin(); it != mpVertexList->ptr_end(); ++it)
    {
        VertexType::Pointer pVertex = *it;
        EdgeType::Pointer pFirstEdge = pVertex->pEdge();
        if (pFirstEdge == NULL)
        {
            // if a vertex has no edge, then all the edge has to be searched if something contains this vertex
            for (EdgeContainerType::ptr_const_iterator it2 = mpEdgeList->ptr_begin(); it2 != mpEdgeList->ptr_end(); ++it2)
            {
                EdgeType::Pointer pEdge = *it2;
                if (pEdge->Node1().Id() == pVertex->Id() || pEdge->Node2().Id() == pVertex->Id())
                {
                    std::stringstream ss;
                    ss << *pEdge << " contains vertex " << *pVertex << ", but the vertex has no edge. It must be wrong during edge assignment";
                    KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, ss.str())
                }
            }

            // no edge contains this vertex; this vertex must be a lone one, skip it.
            continue;
        }

        // the vertex is valid if the edge of the vertex starting from the vertex
        if (pFirstEdge->Node1().Id() != pVertex->Id())
            KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "the edge and the vertex is incompatible")

        // iterate forwardly
        // the vertex is valid if the iteration makes a closed forward loop
        EdgeType::Pointer pThisEdge = pFirstEdge;
        iter = 0;
        do
        {
            pThisEdge = pThisEdge->pPrevEdge()->pOppositeEdge();
            if (pThisEdge == NULL)
                break;
            ++iter;
        } while (pThisEdge != pFirstEdge && iter < max_iters);

        if (iter >= max_iters)
            KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, ": The max iteration is exceeded")

        // iterate reversely
        // the vertex is valid if the iteration makes a closed reverse loop
        pThisEdge = pFirstEdge->pOppositeEdge();
        if (pThisEdge == NULL) // this happens when this edge is a boundary edge
            continue;
        iter = 0;
        do
        {
            pThisEdge = pThisEdge->pNextEdge()->pOppositeEdge();
            if (pThisEdge == NULL)
                break;
            ++iter;
        } while (pThisEdge != pFirstEdge->pOppositeEdge());

        if (iter >= max_iters)
            KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, ": The max iteration is exceeded")
    }

    #ifdef DEBUG_REFINE
    std::cout << "Validate vertices successfully" << std::endl;
    #endif // DEBUG_REFINE

    // validate the half-edge
    for (EdgeContainerType::ptr_const_iterator it = mpEdgeList->ptr_begin(); it != mpEdgeList->ptr_end(); ++it)
    {
        if ((*it)->IsComposite())
            KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "The composite edge is not allowed in mpEdgeList")

        EdgeType::Pointer pFirstEdge = *it;

        // the half-edge is valid if the iteration around a face makes a closed loop
        EdgeType::Pointer pThisEdge = pFirstEdge;
        iter = 0;
        do
        {
            if (pThisEdge->pNextEdge() == NULL)
            {
                std::stringstream ss;
                ss << *pThisEdge << " has no next edge. That must not happen";
                KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, ss.str())
            }
            if (pThisEdge->Node2().Id() != pThisEdge->pNextEdge()->Node1().Id())
                KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "The node id of next edge is not compatible")
            pThisEdge = pThisEdge->pNextEdge();
            ++iter;
        } while (pThisEdge != pFirstEdge && iter < max_iters);

        if (iter >= max_iters)
            KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, ": The max iteration is exceeded")

        // the half-edge is valid if it's not a composite edge
        if ((*it)->IsComposite())
        {
            KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, ": The half-edge is composite")
        }

        // the half-edge is valid if the opposite edge has the reversed vertex ids
        if ((*it)->pOppositeEdge() != NULL)
        {
            if (!((*it)->Node1().Id() == (*it)->pOppositeEdge()->Node2().Id() && (*it)->Node2().Id() == (*it)->pOppositeEdge()->Node1().Id()))
            {
                KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, ": The half-edge opposite has incompatible vertices")
            }
        }
    }

    #ifdef DEBUG_REFINE
    std::cout << "Validate half-edges successfully" << std::endl;
    #endif // DEBUG_REFINE

    // validate the face
    for (FaceContainerType::ptr_const_iterator it = mpFaceList->ptr_begin(); it != mpFaceList->ptr_end(); ++it)
    {
        EdgeType::Pointer pFirstEdge = (*it)->pEdge();
        if ((*it)->IsActive())
        {
            if (pFirstEdge == NULL)
            {
                std::stringstream ss;
                ss << *(*it) << " is active but its edge is null";
                KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, ss.str())
            }
        }
        else
            continue;

        // the face is valid if the iteration of its half-edge around a face makes a closed loop
        EdgeType::Pointer pThisEdge = pFirstEdge;
        iter = 0;
        do
        {
            if (pThisEdge->pFace() != *it)
                KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, ": The face has an incompatible edge")

            pThisEdge = pThisEdge->pNextEdge();
            ++iter;
        } while (pThisEdge != pFirstEdge && iter < max_iters);

        if (iter >= max_iters)
            KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, ": The max iteration is exceeded")

        for (std::size_t i = 0; i < (*it)->NumberOfSubFaces(); ++i)
        {
            if ((*it)->pSubFace(i)->pParent() != (*it))
            {
                KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, ": The sub-face does not point to the correct parent face")
            }
        }
    }

    #ifdef DEBUG_REFINE
    std::cout << "Validate faces successfully" << std::endl;
    #endif // DEBUG_REFINE

    std::cout << "The polytree is validated successfully" << std::endl;
}

void PolyTree2D::ListVertex(std::ostream& rOStream, const std::size_t& Id) const
{
    VertexContainerType::const_iterator it = mpVertexList->find(Id);
    if (it == mpVertexList->end())
    {
        rOStream << "The vertex " << Id << " does not exist" << std::endl;
    }
    else
    {
        rOStream << *it << std::endl;
    }
}

void PolyTree2D::ListFace(std::ostream& rOStream, const std::size_t& Id) const
{
    FaceContainerType::const_iterator it = mpFaceList->find(Id);
    if (it == mpFaceList->end())
    {
        rOStream << "The face " << Id << " does not exist" << std::endl;
    }
    else
    {
        rOStream << *it << std::endl;
    }
}

void PolyTree2D::ListVertices(std::ostream& rOStream) const
{
    for (VertexContainerType::const_iterator it = mpVertexList->begin(); it != mpVertexList->end(); ++it)
        rOStream << *it << std::endl;
}

void PolyTree2D::ListEdges(std::ostream& rOStream) const
{
    for (EdgeContainerType::const_iterator it = mpEdgeList->begin(); it != mpEdgeList->end(); ++it)
        rOStream << *it << std::endl;
}

void PolyTree2D::ListFaces(std::ostream& rOStream) const
{
    for (FaceContainerType::const_iterator it = mpFaceList->begin(); it != mpFaceList->end(); ++it)
        rOStream << *it << std::endl;
}

void PolyTree2D::WriteMatlab(std::ostream& rOStream,
        const bool& write_vertex_number, const bool& write_face_number, const bool& include_lone_node) const
{
    std::size_t max_number_of_nodes = 0;
    for (FaceContainerType::const_iterator it = mpFaceList->begin(); it != mpFaceList->end(); ++it)
    {
        if (it->IsActive())
            if (it->NumberOfNodes() > max_number_of_nodes)
                max_number_of_nodes = it->NumberOfNodes();
    }

    rOStream << "Vertices = [";
    std::vector<std::size_t> vertex_ids;
    std::map<std::size_t, std::size_t> map_vertex; // map from vertex id to matlab row
    std::size_t row = 0;
    for (VertexContainerType::const_iterator it = mpVertexList->begin(); it != mpVertexList->end(); ++it)
    {
        if (!include_lone_node && it->IsLonely()) continue;
        rOStream << (*it)[0] << " " << (*it)[1] << ";" << std::endl;
        map_vertex[it->Id()] = ++row;
        if (write_vertex_number) vertex_ids.push_back(it->Id());
    }
    rOStream << "];" << std::endl;

    rOStream << "Faces = [";
    std::vector<std::size_t> face_ids;
    for (FaceContainerType::const_iterator it = mpFaceList->begin(); it != mpFaceList->end(); ++it)
    {
        if (it->IsActive())
        {
            // rOStream << it->VertexInfo();
            typename EdgeType::Pointer pFirstEdge = it->pEdge();
            typename EdgeType::Pointer pEdge = pFirstEdge;

            do
            {
                rOStream << map_vertex[pEdge->pNode1()->Id()] << " ";
                pEdge = pEdge->pNextEdge();
            } while (pEdge != pFirstEdge);

            // appending NaN information
            for (std::size_t i = 0; i < max_number_of_nodes - it->NumberOfNodes(); ++i)
                rOStream << " NaN";
            rOStream << ";" << std::endl;
            if (write_face_number) face_ids.push_back(it->Id());
        }
    }
    rOStream << "];" << std::endl;

    rOStream << "patch('Faces',Faces,'Vertices',Vertices,'FaceColor','w');" << std::endl;

    if (write_vertex_number)
    {
        rOStream << "vertex_ids = [";
        for (std::size_t i = 0; i < vertex_ids.size(); ++i)
            rOStream << " " << vertex_ids[i];
        rOStream << "];" << std::endl;

        rOStream << "for i = 1:size(Vertices,1)\n";
        rOStream << "    text(Vertices(i,1),Vertices(i,2),num2str(vertex_ids(i)),'color','b');\n";
        rOStream << "end\n";
    }

    if (write_face_number)
    {
        rOStream << "face_ids = [";
        for (std::size_t i = 0; i < face_ids.size(); ++i)
            rOStream << " " << face_ids[i];
        rOStream << "];" << std::endl;

        rOStream << "for i = 1:size(Faces,1)\n";
        rOStream << "    x = 0.0;\n";
        rOStream << "    y = 0.0;\n";
        rOStream << "    n = 0;\n";
        rOStream << "    for j = 1:size(Faces,2)\n";
        rOStream << "        if ~isnan(Faces(i,j))\n";
        rOStream << "            x = x + Vertices(Faces(i,j),1);\n";
        rOStream << "            y = y + Vertices(Faces(i,j),2);\n";
        rOStream << "            n = n+1;\n";
        rOStream << "        end\n";
        rOStream << "    end\n";
        rOStream << "    x = x/n;\n";
        rOStream << "    y = y/n;\n";
        rOStream << "    text(x,y,num2str(face_ids(i)),'color','r')\n";
        rOStream << "end\n";
    }
}

void PolyTree2D::WriteMdpa(std::ostream& rOStream,
        const std::string& ele_prefix, const std::size_t& prop_id) const
{
    rOStream << "//mdpa for KRATOS structural application using polygon elements" << std::endl;
    rOStream << "//(c) 2017 Hoang Giang Bui, Ruhr-University Bochum" << std::endl << std::endl;

    rOStream << "Begin ModelPartData" << std::endl;
    rOStream << "End ModelPartData" << std::endl << std::endl;

    rOStream << "Begin Properties " << prop_id << std::endl;
    rOStream << "End Properties" << std::endl << std::endl;

    rOStream << "Begin Nodes" << std::endl;
    for (VertexContainerType::const_iterator it = mpVertexList->begin(); it != mpVertexList->end(); ++it)
    {
        rOStream << it->Id() << "\t" << (*it)[0] << "\t" << (*it)[1] << "\t0.0" << std::endl;
    }
    rOStream << "End Nodes" << std::endl << std::endl;

    typedef std::map<std::size_t, std::vector<FaceContainerType::const_iterator> > ElementBinType;
    ElementBinType ElementBin;
    for (FaceContainerType::const_iterator it = mpFaceList->begin(); it != mpFaceList->end(); ++it)
    {
        if (it->IsActive())
        {
            ElementBin[it->NumberOfNodes()].push_back(it);
        }
    }

    for (ElementBinType::iterator it = ElementBin.begin(); it != ElementBin.end(); ++it)
    {
        rOStream << "Begin Elements " << ele_prefix << "Polygon" << it->first << "N" << std::endl;
        for (std::size_t i = 0; i < it->second.size(); ++i)
        {
            FaceContainerType::const_iterator it2 = it->second[i];
            rOStream << "\t" << it2->Id() << "\t" << prop_id;

            typename EdgeType::Pointer pFirstEdge = it2->pEdge();
            typename EdgeType::Pointer pEdge = pFirstEdge;

            do
            {
                rOStream << "\t" << pEdge->pNode1()->Id();
                pEdge = pEdge->pNextEdge();
            } while (pEdge != pFirstEdge);

            rOStream << std::endl;
        }
        rOStream << "End Elements" << std::endl << std::endl;
    }
}

}  // namespace Kratos.

#undef DEBUG_REFINE
#undef DEBUG_COARSEN
