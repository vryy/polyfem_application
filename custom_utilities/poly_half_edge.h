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


#if !defined(KRATOS_POLY_HALF_EDGE_H_INCLUDED )
#define  KRATOS_POLY_HALF_EDGE_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <fstream>
#include <array>
#include <cassert>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"

// #define TRACE_DESTRUCTION

namespace Kratos
{

template<std::size_t TDim> class PolyVertex;
template<std::size_t TDim> class PolyHalfEdge;
// template<std::size_t TDim> class PolyHalfEdgeComposite;
template<std::size_t TDim> class PolyFace;

template<std::size_t TDim>
struct PolyHash
{
    std::size_t operator() (PolyHalfEdge<TDim> const& rThis) const 
    {
        std::size_t h;
        std::size_t h1 = std::hash<std::size_t>{} (rThis.Node1().Id());
        std::size_t h2 = std::hash<std::size_t>{} (rThis.Node2().Id());
        h = h1 ^ (h2 << 1);
        if (rThis.pNextEdge() != NULL)
            h = h ^ (PolyHash<TDim>{} (rThis.NextEdge()) << 1);
        if (rThis.pPrevEdge() != NULL)
            h = h ^ (PolyHash<TDim>{} (rThis.PrevEdge()) << 1);
        if (rThis.pOppositeEdge() != NULL)
            h = h ^ (PolyHash<TDim>{} (rThis.OppositeEdge()) << 1);
        return h; 
    }

    std::size_t operator()(PolyFace<TDim> const& rThis) const 
    {
        std::size_t h;
        std::size_t h1 = std::hash<std::size_t>{}(rThis.Id());
        std::size_t h2 = PolyHash<TDim>{}(rThis.Edge());
        h = h1 ^ (h2 << 1);
        return h; 
    }
}; // struct PolyHash


/** 
 * A Polytree vertex in n-dimensional space
 */
template<std::size_t TDim>
class PolyVertex : public std::array<double, TDim>
{
public:

    /// Pointer definition of PolyVertex
    KRATOS_CLASS_POINTER_DEFINITION(PolyVertex);

    typedef PolyHalfEdge<TDim> EdgeType;

    /// Empty constructor for serializer
    PolyVertex() : mId(0) {}

    /// Constructor for PointerVectorSet
    PolyVertex(const std::size_t& Id) : mId(Id) {}

    /// Default constructor.
    PolyVertex(const std::size_t& Id, const double& X, const double& Y)
    : mId(Id)
    {
        if(TDim > 0)
            (*this)[0] = X;

        if(TDim > 1)
            (*this)[1] = Y;
    }

    PolyVertex(const std::size_t& Id, const double& X, const double& Y, const double& Z)
    : mId(Id)
    {
        if(TDim > 0)
            (*this)[0] = X;

        if(TDim > 1)
            (*this)[1] = Y;

        if(TDim > 2)
            (*this)[2] = Z;
    }

    /// Destructor.
    virtual ~PolyVertex()
    {
        #ifdef TRACE_DESTRUCTION
        std::cout << *this << " is destroyed" << std::endl;
        #endif
    }

    /// Returns the id of the associated node
    const std::size_t& Id() const {return mId;}

    /// Set the coordinate in direction i
    void SetCoordinate(const std::size_t& i, const double& x)
    {
        #ifdef NDEBUG
        assert(i < TDim);
        #endif
        (*this)[i] = x;
    }

    /// Get/Set for the underlying half-edge
    EdgeType& Edge() {return *pEdge();}
    const EdgeType& Edge() const {return *pEdge();}
    const typename EdgeType::Pointer pEdge() const {return mpEdge.lock();} // use this for get
    typename EdgeType::Pointer pEdge() {return mpEdge.lock();} // use this for get
    void pSetEdge(typename EdgeType::Pointer pEdge) {mpEdge = pEdge;} // use this for set

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "PolyVertex";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        // rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << Id() << " (";
        for(std::size_t i = 0; i < TDim; ++i)
            rOStream << " " << (*this)[i];
        rOStream << ")";
        rOStream << ", edge: ";
        if(!mpEdge.expired())
            rOStream << *mpEdge.lock();
        else
            rOStream << "null";
    }

private:

    std::size_t mId;

    typename EdgeType::WeakPointer mpEdge; // a half-edge containing this node

    /// Assignment operator.
    PolyVertex& operator=(PolyVertex const& rOther);

    /// Copy constructor.
    PolyVertex(PolyVertex const& rOther);

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
    }

    virtual void load(Serializer& rSerializer)
    {
    }

}; // Class PolyVertex

/**
 * Class representing an edge in the polytree data structure. This class uses half-edge data structure for underlying edge definition.
 * An edge containing references to two nodes and the next and previous half edge inside the face and the reference to the opposite edge.
 */
template<std::size_t TDim>
class PolyHalfEdge
{
public:
    /// Pointer definition of PolyHalfEdge
    KRATOS_CLASS_POINTER_DEFINITION(PolyHalfEdge);

    typedef PolyVertex<TDim> VertexType;
    // typedef PolySuperVertex<TDim> VertexType;
    typedef PolyHalfEdge<TDim> EdgeType;
    typedef PolyFace<TDim> FaceType;

    /// Empty constructor for serializer
    PolyHalfEdge() {}

    /// Constructor for pointer vector set
    PolyHalfEdge(const std::pair<std::size_t, std::size_t>& key) {}

    /// Default constructor.
    /// Note here: Node1->Node2 is understood as edge direction
    PolyHalfEdge(typename VertexType::Pointer pNode1, typename VertexType::Pointer pNode2)
    : mpNode1(pNode1), mpNode2(pNode2)
    {}

    /// Destructor.
    virtual ~PolyHalfEdge()
    {
        #ifdef TRACE_DESTRUCTION
        std::cout << *this << " is destroyed" << std::endl;
        #endif
    }

    std::size_t HashCode() const
    {
//        std::hash<EdgeType> hasher;
//        return hasher(*this);
        return this->operator()(*this);
    }

    void SetSubEdges(std::vector<EdgeType::Pointer>& pEdges)
    {
        if (pEdges.size() > 0)
        {
            mpSubEdges.resize(pEdges.size());
            for (std::size_t i = 0; i < pEdges.size(); ++i)
                mpSubEdges[i] = pEdges[i];
            mpNode1 = pEdges.front()->pNode1();
            mpNode2 = pEdges.back()->pNode2();
        }
        else
            KRATOS_THROW_ERROR(std::logic_error, "Number of given sub-edges is zero", "")
    }

    const std::size_t& Id1() const {return pNode1()->Id();}
    const std::size_t& Id2() const {return pNode2()->Id();}

    virtual VertexType& Node1() {return *pNode1();}
    virtual const VertexType& Node1() const {return *pNode1();}
    virtual typename VertexType::Pointer pNode1() {return mpNode1.lock();}
    virtual const typename VertexType::Pointer pNode1() const {return mpNode1.lock();}
    void pSetNode1(typename VertexType::Pointer pNode) {mpNode1 = pNode;}

    virtual VertexType& Node2() {return *pNode2();}
    virtual const VertexType& Node2() const {return *pNode2();}
    virtual typename VertexType::Pointer pNode2() {return mpNode2.lock();}
    virtual const typename VertexType::Pointer pNode2() const {return mpNode2.lock();}
    void pSetNode2(typename VertexType::Pointer pNode) {mpNode2 = pNode;}

    EdgeType& NextEdge() {return *pNextEdge();}
    const EdgeType& NextEdge() const {return *pNextEdge();}
    const EdgeType::Pointer pNextEdge() const {return mpNextEdge.lock();} // use this for get
    EdgeType::Pointer pNextEdge() {return mpNextEdge.lock();} // use this for get
    void pSetNextEdge(EdgeType::Pointer pEdge) {mpNextEdge = pEdge;} // use this for set

    EdgeType& PrevEdge() {return *pPrevEdge();}
    const EdgeType& PrevEdge() const {return *pPrevEdge();}
    const EdgeType::Pointer pPrevEdge() const {return mpPrevEdge.lock();} // use this for get
    EdgeType::Pointer pPrevEdge() {return mpPrevEdge.lock();} // use this for get
    void pSetPrevEdge(EdgeType::Pointer pEdge) {mpPrevEdge = pEdge;} // use this for set

    EdgeType& OppositeEdge() {return *pOppositeEdge();}
    const EdgeType& OppositeEdge() const {return *pOppositeEdge();}
    const EdgeType::Pointer pOppositeEdge() const {return mpOppositeEdge.lock();} // use this for get
    EdgeType::Pointer pOppositeEdge() {return mpOppositeEdge.lock();} // use this for get
    void pSetOppositeEdge(EdgeType::Pointer pEdge) {mpOppositeEdge = pEdge;} // use this for set

    FaceType& Face() {return *pFace();}
    const FaceType& Face() const {return *pFace();}
    const typename FaceType::Pointer pFace() const {return mpFace.lock();} // use this for get
    typename FaceType::Pointer pFace() {return mpFace.lock();} // use this for get
    void pSetFace(typename FaceType::Pointer pFace) {mpFace = pFace;} // use this for set

    /// Check if this edge is a composite edge
    virtual bool IsComposite() const {return mpSubEdges.size() != 0;}

    std::size_t NumberOfSubEdges() const {return mpSubEdges.size();}
    EdgeType& SubEdge(const std::size_t& i) {return *pSubEdge(i);}
    const EdgeType& SubEdge(const std::size_t& i) const {return *pSubEdge(i);}
    const typename EdgeType::Pointer pSubEdge(const std::size_t& i) const {return mpSubEdges[i].lock();} // use this for get
    typename EdgeType::Pointer pSubEdge(const std::size_t& i) {return mpSubEdges[i].lock();} // use this for get

    /// Get the number of nodes
    std::size_t NumberOfNodes() const
    {
        if (IsComposite())
        {
            return NumberOfSubEdges() + 1;
        }
        else
        {
            return 2;
        }
    }

    /// Check if this edge is collinear with another edge
    bool IsCollinear(typename EdgeType::Pointer pEdge) const
    {
        if (IsComposite() || pEdge->IsComposite())
        {
            KRATOS_THROW_ERROR(std::logic_error, "One of the edge is composite. It's not allowed to check collinearity of the composite edge with other.", "")
            return false;
        }
        else
        {
            double a1 = (Node2()[1] - Node1()[1]) / (Node2()[0] - Node1()[0]);
            double a2 = (pEdge->Node2()[1] - pEdge->Node1()[1]) / (pEdge->Node2()[0] - pEdge->Node1()[0]);
            return fabs(fabs(a1) - fabs(a2)) < 1.0e-6; // TODO shall we parameterize this
        }
    }

    /// Get the string representing the vertices on the edge
    std::string VertexInfo() const
    {
        std::stringstream ss;

        if (IsComposite())
        {
            ss << "(";
            for (std::size_t i = 0; i < NumberOfSubEdges(); ++i)
            {
                if (pSubEdge(i) != NULL)
                    ss << pSubEdge(i)->pNode1()->Id() << " ";
                else
                    ss << "null ";
            }
            if (mpSubEdges.back().lock() != NULL)
                ss << mpSubEdges.back().lock()->pNode2()->Id() << ")";
            else
                ss << "null)";
        }
        else
        {
            ss << "(" << pNode1()->Id() << " " << pNode2()->Id() << ")";
        }

        return ss.str();
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        if (!IsComposite())
            return "PolyHalfEdge";
        else
            return "PolyHalfEdgeComposite";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << VertexInfo();

        rOStream << ", next edge: ";
        if(pNextEdge() != NULL)
            rOStream << pNextEdge()->VertexInfo();
        else
            rOStream << "null";

        rOStream << ", prev edge: ";
        if(pPrevEdge() != NULL)
            rOStream << pPrevEdge()->VertexInfo();
        else
            rOStream << "null";

        rOStream << ", opposite edge: ";
        if(pOppositeEdge() != NULL)
            rOStream << pOppositeEdge()->VertexInfo();
        else
            rOStream << "null";

        rOStream << ", face: ";
        if(pFace() != NULL)
            rOStream << pFace()->Id();
        else
            rOStream << "null";

        if (IsComposite())
        {
            rOStream << ", sub-edges:";
            for (std::size_t i = 0; i < NumberOfSubEdges(); ++i)
            {
                if (pSubEdge(i) != NULL)
                    rOStream << " " << pSubEdge(i)->VertexInfo();
                else
                    rOStream << " null";
            }
        }
    }

private:

    typename EdgeType::WeakPointer mpPrevEdge;
    typename EdgeType::WeakPointer mpNextEdge;
    typename EdgeType::WeakPointer mpOppositeEdge;

    typename VertexType::WeakPointer mpNode1;
    typename VertexType::WeakPointer mpNode2;

    typename FaceType::WeakPointer mpFace;

    typename std::vector<EdgeType::WeakPointer> mpSubEdges;

    /// Assignment operator.
    PolyHalfEdge& operator=(PolyHalfEdge const& rOther)
    {
        mpPrevEdge = rOther.mpPrevEdge;
        mpNextEdge = rOther.mpNextEdge;
        mpOppositeEdge = rOther.mpOppositeEdge;
        mpNode1 = rOther.mpNode1;
        mpNode2 = rOther.mpNode2;
        mpSubEdges = rOther.mpSubEdges;
    }

    /// Copy constructor.
    PolyHalfEdge(PolyHalfEdge const& rOther)
    {
        mpPrevEdge = rOther.mpPrevEdge;
        mpNextEdge = rOther.mpNextEdge;
        mpOppositeEdge = rOther.mpOppositeEdge;
        mpNode1 = rOther.mpNode1;
        mpNode2 = rOther.mpNode2;
        mpSubEdges = rOther.mpSubEdges;
    }

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
    }

    virtual void load(Serializer& rSerializer)
    {
    }
}; // Class PolyHalfEdge

/**
 * A PolyTree Face
 * The polytree face contains a list of sub-faces, which results from refinement.
 */
template<std::size_t TDim>
class PolyFace
{
public:

    /// Pointer definition of PolyFace
    KRATOS_CLASS_POINTER_DEFINITION(PolyFace);

    typedef PolyHalfEdge<TDim> EdgeType;
    typedef typename EdgeType::VertexType VertexType;
    typedef PolyFace<TDim> FaceType;

    /// Empty constructor for serializer.
    PolyFace() : mId(0), m_is_active(true), m_is_refined(false), m_is_coarsen(false), m_is_changed(false) {}

    /// Default constructor.
    PolyFace(const std::size_t& Id)
    : mId(Id), m_is_active(true), m_is_refined(false), m_is_coarsen(false), m_is_changed(false)
    {}

    /// Destructor.
    virtual ~PolyFace()
    {
        #ifdef TRACE_DESTRUCTION
        std::cout << *this << " is destroyed" << std::endl;
        #endif
    }

    /// Returns the id of the associated face
    const std::size_t& Id() const {return mId;}

    /// Get/Set for the underlying half-edge
    EdgeType& Edge() {return *pEdge();}
    const EdgeType& Edge() const {return *pEdge();}
    const typename EdgeType::Pointer pEdge() const {return mpEdge.lock();} // use this for get
    typename EdgeType::Pointer pEdge() {return mpEdge.lock();} // use this for get
    void pSetEdge(typename EdgeType::Pointer pEdge) {mpEdge = pEdge;} // use this for set

    /// Get/Set for parent face
    FaceType& Parent() {return *pParent();}
    const FaceType& Parent() const {return *pParent();}
    const typename FaceType::Pointer pParent() const {return mpParent.lock();} // use this for get
    typename FaceType::Pointer pParent() {return mpParent.lock();} // use this for get
    void pSetParent(typename FaceType::Pointer pParent) {mpParent = pParent;} // use this for set

    /// Set the sub-faces
    void SetSubFaces(std::vector<typename FaceType::Pointer> pSubFaces)
    {
        mpSubFaces.resize(pSubFaces.size());
        for (std::size_t i = 0; i < pSubFaces.size(); ++i)
            mpSubFaces[i] = pSubFaces[i];
    }

    /// Check if this tree is a leaf
    bool IsLeaf() const {return mpSubFaces.size() == 0;}

    std::size_t NumberOfSubFaces() const {return mpSubFaces.size();}
    FaceType& SubFace(const std::size_t& i) {return *pSubFace(i);}
    const FaceType& SubFace(const std::size_t& i) const {return *pSubFace(i);}
    const FaceType::Pointer pSubFace(const std::size_t& i) const {return mpSubFaces[i].lock();} // use this for get
    FaceType::Pointer pSubFace(const std::size_t& i) {return mpSubFaces[i].lock();} // use this for get

    bool IsChanged() const {return m_is_changed;}
    void SetChange(const bool& is_changed) {m_is_changed = is_changed;}

    bool IsRefined() const {return m_is_refined;}
    void SetRefine(const bool& is_refined) {m_is_refined = is_refined;}

    bool IsCoarsen() const {return m_is_coarsen;}
    void SetCoarsen(const bool& is_coarsen) {m_is_coarsen = is_coarsen;}

    bool IsActive() const {return m_is_active;}
    void SetActive(const bool& is_active) {m_is_active = is_active;}

    /**
     * Extract and add all the vertices of the face to the container
     * @param rVertexList the output vertex container
     */
    template<class VertexContainerType>
    void AddVerticesTo(VertexContainerType& rVertexList) const
    {
        if (IsLeaf())
        {
            typename EdgeType::Pointer pFirstEdge = pEdge();
            typename EdgeType::Pointer pEdge = pFirstEdge;

            do
            {
                if (pEdge->IsComposite())
                {
                    for (std::size_t i = 0; i < pEdge->NumberOfSubEdges(); ++i)
                    {
                        typename VertexType::Pointer pNode = pEdge->pSubEdge(i)->pNode1();
                        #ifdef NDEBUG
                        assert(pNode != NULL);
                        #endif
                        rVertexList.insert(rVertexList.begin(), pNode);
                    }
                }
                else
                {
                    typename VertexType::Pointer pNode = pEdge->pNode1();
                    #ifdef NDEBUG
                    assert(pNode != NULL);
                    #endif
                    rVertexList.insert(rVertexList.begin(), pNode);
                }
                pEdge = pEdge->pNextEdge();
            } while (pEdge != pFirstEdge);
        }
        else
        {
            for (std::size_t i = 0; i < NumberOfSubFaces(); ++i)
            {
                SubFace(i).AddVerticesTo(rVertexList);
            }
        }
    }

    std::size_t NumberOfNodes() const
    {
        std::size_t number_of_nodes = 0;
        std::size_t number_of_edges = 0;

        typename EdgeType::Pointer pFirstEdge = pEdge();
        typename EdgeType::Pointer pEdge = pFirstEdge;

        do
        {
            ++number_of_edges;
            number_of_nodes += pEdge->NumberOfNodes();
            pEdge = pEdge->pNextEdge();
        } while (pEdge != pFirstEdge);

        return number_of_nodes - number_of_edges;
    }

    /// Get the number of edges of the polygon. If the face is leaf, it is also the number of nodes.
    std::size_t NumberOfEdges() const
    {
        std::size_t number_of_edges = 0;

        typename EdgeType::Pointer pFirstEdge = pEdge();
        typename EdgeType::Pointer pEdge = pFirstEdge;

        do
        {
            ++number_of_edges;
            pEdge = pEdge->pNextEdge();
        } while (pEdge != pFirstEdge);

        return number_of_edges;
    }

    
    /**
     * Extract the minimal polygon from this face. This polygon does not contain middle nodes on edge.
     * Note that this function requires the face to be complete, i.e. the half-edges make a closed loop. Do not call this function after BeginRefineCoarsen() and before EndRefineCoarsen().
     * @param polygon   the output polygon coordinates
     * @param pVertices the vertices of the polygon
     * @param pEdges    the edges of the polygon. In case the edge contains middle node, the composite edge will be created.
     */
    void ExtractMinimalPolygon(std::vector<std::vector<double> >& polygon,
            std::vector<typename VertexType::Pointer>& pVertices,
            std::vector<typename EdgeType::Pointer>& pEdges) const
    {
        typename EdgeType::Pointer pFirstEdge = pEdge();
        typename EdgeType::Pointer pEdge = pFirstEdge;
        typename EdgeType::Pointer pNextEdge = pEdge->pNextEdge();
        std::vector<typename EdgeType::Pointer> pSubEdges;

        pSubEdges.push_back(pEdge);
        do
        {
            if (pEdge->IsCollinear(pNextEdge))
            {
                pSubEdges.push_back(pNextEdge);
                pNextEdge = pNextEdge->pNextEdge();
            }
            else
            {
                std::vector<double> point(2);
                point[0] = pEdge->Node1()[0];
                point[1] = pEdge->Node1()[1];
                polygon.push_back(point);
                pVertices.push_back(pEdge->pNode1());
                pEdge = pNextEdge;
                pNextEdge = pNextEdge->pNextEdge();

                if (pSubEdges.size() > 1)
                {
                    typename EdgeType::Pointer pNewEdge = boost::make_shared<EdgeType>(pSubEdges.front()->pNode1(), pSubEdges.back()->pNode2());
                    pNewEdge->SetSubEdges(pSubEdges);
                    pEdges.push_back(pNewEdge);
                }
                else
                    pEdges.push_back(pSubEdges[0]);

                pSubEdges.clear();
                pSubEdges.push_back(pEdge);
            }
        } while (pNextEdge != pFirstEdge->pNextEdge());
    }

    /// Get the string representing the vertices on the polygon (not counting the one on the edge)
    std::string VertexInfo() const
    {
        std::stringstream ss;

        typename EdgeType::Pointer pFirstEdge = pEdge();
        if (pFirstEdge == NULL)
        {
            ss << "null";
            return ss.str();
        }
        typename EdgeType::Pointer pEdge = pFirstEdge;

        do
        {
            ss << pEdge->pNode1()->Id() << " ";
            pEdge = pEdge->pNextEdge();
        } while (pEdge != pFirstEdge);

        return ss.str();
    }

    /// Get the string representing the vertices on the edge
    std::string EdgeInfo() const
    {
        std::stringstream ss;

        typename EdgeType::Pointer pFirstEdge = pEdge();
        if (pFirstEdge == NULL)
        {
            ss << "null";
            return ss.str();
        }
        typename EdgeType::Pointer pEdge = pFirstEdge;

        ss << "( ";
        do
        {
            ss << pEdge->VertexInfo() << " ";
            pEdge = pEdge->pNextEdge();
        } while (pEdge != pFirstEdge);
        ss << ")";

        return ss.str();
    }

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "PolyFace";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << Id();
// KRATOS_WATCH(__LINE__)
        if (IsActive())
            rOStream << ", active,";
        else
            rOStream << ", inactive,";
// KRATOS_WATCH(__LINE__)
        rOStream << " Edge " << EdgeInfo();
// KRATOS_WATCH(__LINE__)
        rOStream << ", half-edge ";
        if(pEdge() != NULL)
            rOStream << pEdge()->VertexInfo();
        else
            rOStream << "null";
// KRATOS_WATCH(__LINE__)
        rOStream << ", parent:";
        if(pParent() != NULL)
            rOStream << pParent()->Id();
        else
            rOStream << "null";
// KRATOS_WATCH(__LINE__)
        rOStream << ", sub-faces:";
        for(std::size_t i = 0; i < NumberOfSubFaces(); ++i)
        {
            rOStream << " " << SubFace(i).Id();
        }
    }

    /// Assignment operator.
    PolyFace& operator=(PolyFace const& rOther)
    {
        mId = rOther.mId;
        m_is_changed = rOther.m_is_changed;
        m_is_refined = rOther.m_is_refined;
        mpEdge = rOther.mpEdge;
    }

    /// Copy constructor.
    PolyFace(PolyFace const& rOther)
    {
        mId = rOther.mId;
        m_is_changed = rOther.m_is_changed;
        m_is_refined = rOther.m_is_refined;
        mpEdge = rOther.mpEdge;
    }

    /// Compare operator
    friend bool operator==(const PolyFace& Left, const PolyFace& Right)
    {
        return (Left.mId == Right.mId);
    }

    friend bool operator<(const PolyFace& Left, const PolyFace& Right)
    {
        return (Left.mId < Right.mId);
    }

private:

    std::size_t mId;

    bool m_is_changed; // flag to mark if this face is changed

    bool m_is_refined; // flag to mark if this face is about to be refined
    
    bool m_is_coarsen; // flag to mark if this face is about to be coarsen

    bool m_is_active;  // flag to indicate if this face is active/inactive. For example, if the face is refined, it will become inactive

    // bool m_is_removed; // TODO consider this flag

    typename EdgeType::WeakPointer mpEdge; // a half-edge inside this face

    typename FaceType::WeakPointer mpParent; // pointer to the parent face

    std::vector<typename FaceType::WeakPointer> mpSubFaces; // list of sub-faces (in polytree)
    // if mpSubFaces.size() is nonzero then the mpEdge must be NULL, because at this time the face is not at the bottom of the tree anymore.

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
    }

    virtual void load(Serializer& rSerializer)
    {
    }
}; // Class PolyFace



/// input stream function
template<std::size_t TDim>
inline std::istream& operator >> (std::istream& rIStream, PolyVertex<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template<std::size_t TDim>
inline std::ostream& operator << (std::ostream& rOStream, const PolyVertex<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " ";
    rThis.PrintData(rOStream);

    return rOStream;
}

// /// input stream function
// template<std::size_t TDim>
// inline std::istream& operator >> (std::istream& rIStream, PolySuperVertex<TDim>& rThis)
// {
//     return rIStream;
// }

// /// output stream function
// template<std::size_t TDim>
// inline std::ostream& operator << (std::ostream& rOStream, const PolySuperVertex<TDim>& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << " ";
//     rThis.PrintData(rOStream);

//     return rOStream;
// }

/// input stream function
template<std::size_t TDim>
inline std::istream& operator >> (std::istream& rIStream, PolyHalfEdge<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template<std::size_t TDim>
inline std::ostream& operator << (std::ostream& rOStream, const PolyHalfEdge<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " ";
    rThis.PrintData(rOStream);

    return rOStream;
}

/// input stream function
template<std::size_t TDim>
inline std::istream& operator >> (std::istream& rIStream, PolyFace<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template<std::size_t TDim>
inline std::ostream& operator << (std::ostream& rOStream, const PolyFace<TDim>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " ";
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos.

#undef TRACE_DESTRUCTION

#endif // KRATOS_POLY_HALF_EDGE_H_INCLUDED  defined
