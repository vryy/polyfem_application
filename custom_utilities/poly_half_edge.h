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


namespace Kratos
{

template<std::size_t TDim> class PolyVertex;
template<std::size_t TDim> class PolyHalfEdge;
template<std::size_t TDim> class PolyHalfEdgeComposite;
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

    std::size_t operator() (PolyHalfEdgeComposite<TDim> const& rThis) const 
    {
        std::size_t h;
        for (std::size_t i = 0; i < rThis.NumberOfSubEdges(); ++i)
        {
            if (i == 0)
            {
                h = PolyHash<TDim>{} (rThis.SubEdge(i));
            }
            else
            {
                h = h ^ (PolyHash<TDim>{} (rThis.SubEdge(i)) << 1);
            }
        }
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


/// Short class definition.
/** A Polytree node in n-dimensional space
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
    virtual ~PolyVertex() {}

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
    EdgeType& Edge() {return *mpEdge;}
    const EdgeType& Edge() const {return *mpEdge;}
    const typename EdgeType::Pointer pEdge() const {return mpEdge;} // use this for get
    typename EdgeType::Pointer pEdge() {return mpEdge;} // use this for get
    void pSetEdge(typename EdgeType::Pointer pEdge) {mpEdge = pEdge;} // use this for set

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "PolyVertex";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << Id() << " (";
        for(std::size_t i = 0; i < TDim; ++i)
            rOStream << " " << (*this)[i];
        rOStream << ")";
        rOStream << ", edge: ";
        if(pEdge() != NULL)
            rOStream << PolyHash<TDim>{}(Edge());
        else
            rOStream << "null";
    }

private:

    std::size_t mId;

    typename EdgeType::Pointer mpEdge; // a half-edge containing this node

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

/// Short class definition.
/** class representing an edge in the polytree data structure. This class uses half-edge data structure for underlying edge definition.
 * An edge containing references to two nodes and the next and previous half edge inside the face and the reference to the opposite edge.
 */
template<std::size_t TDim>
class PolyHalfEdge
{
public:
    /// Pointer definition of PolyHalfEdge
    KRATOS_CLASS_POINTER_DEFINITION(PolyHalfEdge);

    typedef PolyVertex<TDim> NodeType;
    typedef PolyHalfEdge<TDim> EdgeType;
    typedef PolyFace<TDim> FaceType;

    /// Default constructor.
    /// Note here: Node1->Node2 is understood as edge direction
    PolyHalfEdge(typename NodeType::Pointer pNode1, typename NodeType::Pointer pNode2)
    : mpNode1(pNode1), mpNode2(pNode2)
    {}

    /// Destructor.
    virtual ~PolyHalfEdge() {}

    std::size_t HashCode() const
    {
//        std::hash<EdgeType> hasher;
//        return hasher(*this);
        return this->operator()(*this);
    }

    virtual NodeType& Node1() {return *mpNode1;}
    virtual const NodeType& Node1() const {return *mpNode1;}
    virtual typename NodeType::Pointer pNode1() {return mpNode1;}
    virtual const typename NodeType::Pointer pNode1() const {return mpNode1;}

    virtual NodeType& Node2() {return *mpNode2;}
    virtual const NodeType& Node2() const {return *mpNode2;}
    virtual typename NodeType::Pointer pNode2() {return mpNode2;}
    virtual const typename NodeType::Pointer pNode2() const {return mpNode2;}

    EdgeType& NextEdge() {return *mpNextEdge;}
    const EdgeType& NextEdge() const {return *mpNextEdge;}
    const EdgeType::Pointer pNextEdge() const {return mpNextEdge;} // use this for get
    EdgeType::Pointer pNextEdge() {return mpNextEdge;} // use this for get
    void pSetNextEdge(EdgeType::Pointer pEdge) {mpNextEdge = pEdge;} // use this for set

    EdgeType& PrevEdge() {return *mpPrevEdge;}
    const EdgeType& PrevEdge() const {return *mpPrevEdge;}
    const EdgeType::Pointer pPrevEdge() const {return mpPrevEdge;} // use this for get
    EdgeType::Pointer pPrevEdge() {return mpPrevEdge;} // use this for get
    void pSetPrevEdge(EdgeType::Pointer pEdge) {mpPrevEdge = pEdge;} // use this for set

    EdgeType& OppositeEdge() {return *mpOppositeEdge;}
    const EdgeType& OppositeEdge() const {return *mpOppositeEdge;}
    const EdgeType::Pointer pOppositeEdge() const {return mpOppositeEdge;} // use this for get
    EdgeType::Pointer pOppositeEdge() {return mpOppositeEdge;} // use this for get
    void pSetOppositeEdge(EdgeType::Pointer pEdge) {mpOppositeEdge = pEdge;} // use this for set

    FaceType& Face() {return *mpFace;}
    const FaceType& Face() const {return *mpFace;}
    const typename FaceType::Pointer pFace() const {return mpFace;} // use this for get
    typename FaceType::Pointer pFace() {return mpFace;} // use this for get
    void pSetFace(typename FaceType::Pointer pFace) {mpFace = pFace;} // use this for set

    /// Check if this edge is a composite edge
    virtual bool IsComposite() const {return false;}

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "PolyHalfEdge";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "Node: (" << mpNode1->Id() << ", " << mpNode2->Id() << ")";

        rOStream << ", next edge: ";
        if(pNextEdge() != NULL)
            rOStream << "(" << NextEdge().Node1().Id() << ", " << NextEdge().Node2().Id() << ")";
        else
            rOStream << "null";

        rOStream << ", prev edge: ";
        if(pPrevEdge() != NULL)
            rOStream << "(" << PrevEdge().Node1().Id() << ", " << PrevEdge().Node2().Id() << ")";
        else
            rOStream << "null";

        rOStream << ", opposite edge: ";
        if(pOppositeEdge() != NULL)
            rOStream << "(" << OppositeEdge().Node1().Id() << ", " << OppositeEdge().Node2().Id() << ")";
        else
            rOStream << "null";

        rOStream << ", face: ";
        if(pFace() != NULL)
            rOStream << Face().Id();
        else
            rOStream << "null";
    }

private:

    typename EdgeType::Pointer mpPrevEdge;
    typename EdgeType::Pointer mpNextEdge;
    typename EdgeType::Pointer mpOppositeEdge;

    typename NodeType::Pointer mpNode1;
    typename NodeType::Pointer mpNode2;

    typename FaceType::Pointer mpFace;

    /// Assignment operator.
    PolyHalfEdge& operator=(PolyHalfEdge const& rOther)
    {
        mpPrevEdge = rOther.mpPrevEdge;
        mpNextEdge = rOther.mpNextEdge;
        mpOppositeEdge = rOther.mpOppositeEdge;
        mpNode1 = rOther.mpNode1;
        mpNode2 = rOther.mpNode2;
    }

    /// Copy constructor.
    PolyHalfEdge(PolyHalfEdge const& rOther)
    {
        mpPrevEdge = rOther.mpPrevEdge;
        mpNextEdge = rOther.mpNextEdge;
        mpOppositeEdge = rOther.mpOppositeEdge;
        mpNode1 = rOther.mpNode1;
        mpNode2 = rOther.mpNode2;
    }

}; // Class PolyHalfEdge

/// Short class definition.
/** class representing an edge comprised of consecutive half-edge segment.
 */
template<std::size_t TDim>
class PolyHalfEdgeComposite : public PolyHalfEdge<TDim>
{
public:
    /// Pointer definition of PolyHalfEdge
    KRATOS_CLASS_POINTER_DEFINITION(PolyHalfEdgeComposite);

    typedef PolyHalfEdge<TDim> BaseType;
    typedef typename BaseType::NodeType NodeType;
    typedef typename BaseType::EdgeType EdgeType;
    typedef typename BaseType::FaceType FaceType;

    /// Default constructor.
    PolyHalfEdgeComposite(std::vector<typename EdgeType::Pointer> pEdges)
    {
        if(pEdges.size() > 0)
            mpSubEdges = pEdges;
        else
            KRATOS_THROW_ERROR(std::logic_error, "Empty edge list given to the constructor", "")
    }

    PolyHalfEdgeComposite(typename EdgeType::Pointer pEdge1, typename EdgeType::Pointer pEdge2)
    {
        mpSubEdges.push_back(pEdge1);
        mpSubEdges.push_back(pEdge2);
    }

    /// Destructor.
    virtual ~PolyHalfEdgeComposite() {}

    std::size_t HashCode() const
    {
        return this->operator()(*this);
    }

    virtual NodeType& Node1() {return mpSubEdges.front()->Node1();}
    virtual const NodeType& Node1() const {return mpSubEdges.front()->Node1();}
    virtual typename NodeType::Pointer pNode1() {return mpSubEdges.front()->pNode1();}
    virtual const typename NodeType::Pointer pNode1() const {return mpSubEdges.front()->pNode1();}

    virtual NodeType& Node2() {return mpSubEdges.back()->Node2();}
    virtual const NodeType& Node2() const {return mpSubEdges.back()->Node2();}
    virtual typename NodeType::Pointer pNode2() {return mpSubEdges.back()->pNode2();}
    virtual const typename NodeType::Pointer pNode2() const {return mpSubEdges.back()->pNode2();}

    std::size_t NumberOfSubEdges() const {return mpSubEdges.size();}
    EdgeType& SubEdge(const std::size_t& i) {return *mpSubEdges[i];}
    const EdgeType& SubEdge(const std::size_t& i) const {return *mpSubEdges[i];}
    const typename EdgeType::Pointer pSubEdge(const std::size_t& i) const {return mpSubEdges[i];} // use this for get
    typename EdgeType::Pointer pSubEdge(const std::size_t& i) {return mpSubEdges[i];} // use this for get

    virtual bool IsComposite() const {return true;}

    /// Turn back information as a string.
    virtual std::string Info() const
    {
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
        rOStream << "sub-edges:";
        for(std::size_t i = 0; i < NumberOfSubEdges(); ++i)
        {
            rOStream << " (" << SubEdge(i).Node1().Id() << ", " << SubEdge(i).Node2().Id() << ")";
        }
    }

private:

    std::vector<typename EdgeType::Pointer> mpSubEdges;

    /// Assignment operator.
    PolyHalfEdgeComposite& operator=(PolyHalfEdgeComposite const& rOther)
    {
        mpSubEdges = rOther.mpSubEdges;
    }

    /// Copy constructor.
    PolyHalfEdgeComposite(PolyHalfEdgeComposite const& rOther)
    {
        mpSubEdges = rOther.mpSubEdges;
    }

}; // Class PolyHalfEdgeComposite

/// Short class definition.
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

    typedef PolyVertex<TDim> NodeType;
    typedef PolyHalfEdge<TDim> EdgeType;
    typedef PolyFace<TDim> FaceType;

    /// Empty constructor for serializer.
    PolyFace() : mId(0) {}

    /// Default constructor.
    PolyFace(const std::size_t& Id)
    : mId(Id)
    {}

    /// Destructor.
    virtual ~PolyFace() {}

    /// Returns the id of the associated face
    const std::size_t& Id() const {return mId;}

    /// Get/Set for the underlying half-edge
    EdgeType& Edge() {return *mpEdge;}
    const EdgeType& Edge() const {return *mpEdge;}
    const typename EdgeType::Pointer pEdge() const {return mpEdge;} // use this for get
    typename EdgeType::Pointer pEdge() {return mpEdge;} // use this for get
    void pSetEdge(typename EdgeType::Pointer pEdge) {mpEdge = pEdge;} // use this for set

    /// Add a face to list
    void AddFace(typename FaceType::Pointer pFace)
    {
        mpSubFaces.push_back(pFace);
        mpEdge = NULL;
    }

    /// Check if this tree is a leaf
    bool IsLeaf() const {return mpSubFaces.size() == 0;}

    std::size_t NumberOfSubFaces() const {return mpSubFaces.size();}
    FaceType& SubFace(const std::size_t& i) {return *mpSubFaces[i];}
    const FaceType& SubFace(const std::size_t& i) const {return *mpSubFaces[i];}
    const FaceType::Pointer pSubFace(const std::size_t& i) const {return mpSubFaces[i];} // use this for get
    FaceType::Pointer pSubFace(const std::size_t& i) {return mpSubFaces[i];} // use this for get

    bool IsChanged() const {return m_is_changed;}
    void SetChange(const bool& is_changed) {m_is_changed = is_changed;}

    bool IsRefined() const {return m_is_refined;}
    void SetRefine(const bool& is_refined) {m_is_refined = is_refined;}

    bool IsCoarsen() const {return m_is_coarsen;}
    void SetCoarsen(const bool& is_coarsen) {m_is_coarsen = is_coarsen;}

    /**
     * Extract and add all the vertices of the face to the container
     * @param rVertexList the output vertex container
     */
    template<class VertexContainerType>
    void AddVertices(VertexContainerType& rVertexList) const
    {
        if (IsLeaf())
        {
            typename EdgeType::Pointer pFirstEdge = pEdge();
            typename EdgeType::Pointer pEdge = pFirstEdge;

            do
            {
                if (pEdge->IsComposite())
                {
                    typename PolyHalfEdgeComposite<TDim>::Pointer pEdgeComposite =
                            boost::dynamic_pointer_cast<PolyHalfEdgeComposite<TDim> >(pEdge);
                    for (std::size_t i = 0; i < pEdgeComposite->NumberOfSubEdges(); ++i)
                    {
                        typename NodeType::Pointer pNode = pEdgeComposite->pSubEdge(i)->pNode1();
                        #ifdef NDEBUG
                        assert(pNode != NULL);
                        #endif
                        rVertexList.insert(rVertexList.begin(), pNode);
                    }
                }
                else
                {
                    typename NodeType::Pointer pNode = pEdge->pNode1();
                    #ifdef NDEBUG
                    assert(pNode != NULL);
                    #endif
                    rVertexList.insert(rVertexList.begin(), pNode);
                }
                pEdge = pFirstEdge->pNextEdge();
            } while (pEdge != pEdge);
        }
        else
        {
            for (std::size_t i = 0; i < NumberOfSubFaces(); ++i)
            {
                SubFace(i).AddVertices(rVertexList);
            }
        }
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
        rOStream << ", edge: ";
        if(pEdge() != NULL)
            rOStream << PolyHash<TDim>{}(Edge());
        else
            rOStream << "null";
        rOStream << std::endl;

        rOStream << ", sub-faces:";
        for(std::size_t i = 0; i < NumberOfSubFaces(); ++i)
        {
            rOStream << " " << SubFace(i).Id();
        }
        rOStream << std::endl;
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

    typename EdgeType::Pointer mpEdge; // a half-edge inside this face

    std::vector<typename FaceType::Pointer> mpSubFaces; // list of sub-faces (in polytree)
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
inline std::istream& operator >> (std::istream& rIStream, PolyHalfEdgeComposite<TDim>& rThis)
{
    return rIStream;
}

/// output stream function
template<std::size_t TDim>
inline std::ostream& operator << (std::ostream& rOStream, const PolyHalfEdgeComposite<TDim>& rThis)
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


#endif // KRATOS_POLY_HALF_EDGE_H_INCLUDED  defined
