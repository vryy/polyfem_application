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


namespace Kratos
{
///@addtogroup FiniteCellApplication
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

    /// Pointer definition of PolyTreeSyncUtility
    KRATOS_CLASS_POINTER_DEFINITION(PolyTreeSyncUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PolyTreeSyncUtility() {}

    /// Destructor.
    virtual ~PolyTreeSyncUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Initialize the half edge data structure for the polytree_2d from ModelPart
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
            std::size_t& LastVertexId, std::size_t& LastFaceId) const
    {
        typedef ModelPart::ElementType::GeometryType GeometryType;
        typedef typename TVertexContainerType::data_type VertexType;
        typedef typename TEdgeContainerType::data_type EdgeType;
        typedef typename TFaceContainerType::data_type FaceType;

        // create half-edge vertices
        LastVertexId = 0;
        for (ModelPart::NodeIterator it = r_model_part.NodesBegin(); it != r_model_part.NodesEnd() ; ++it)
        {
            typename VertexType::Pointer pNewVertex = boost::make_shared<VertexType>(it->Id(), it->X(), it->Y());
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
    }


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
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


#endif // KRATOS_POLYTREE_SYNC_UTILITY_H_INCLUDED  defined
