//  License: see polyfem_application/LICENSE.txt
//
//  Main authors:    Hoang Giang Bui
//  Date:   11 Aug 2017
//


#if !defined(KRATOS_POLYGON_H_INCLUDED )
#define  KRATOS_POLYGON_H_INCLUDED

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "geometries/line_2d_2.h"
#include "utilities/math_utils.h"
#include "custom_quadrature/polygon_gauss_legendre_integration_points.h"


namespace Kratos
{
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

/**
 * Polygon geometry in 2D using Wachspress shape function.
 REF: Paulino's PolyTop
 */

template<class TPointType, std::size_t TnVertices>
class Polygon : public Geometry<TPointType>
{
public:
    ///@}
    ///@name Type Definitions
    ///@{

    /**
     * Geometry as base class.
     */
    typedef Geometry<TPointType> BaseType;

    /**
     * Type of edge geometry
     */
    typedef Line2D2<TPointType> EdgeType;

    /**
     * Pointer definition of Polygon
     */
    KRATOS_CLASS_POINTER_DEFINITION( Polygon );

    /**
     * Integration methods implemented in geometry.
     */
    typedef GeometryData::IntegrationMethod IntegrationMethod;

    /**
     * A Vector of counted pointers to Geometries.
     * Used for returning edges of the geometry.
     */
    typedef typename BaseType::GeometriesArrayType GeometriesArrayType;

    /**
     * Redefinition of template parameter TPointType.
     */
    typedef TPointType PointType;

    /**
     * Type used for indexing in geometry class.
     * std::size_t used for indexing
     * point or integration point access methods and also all other
     * methods which need point or integration point index.
     */
    typedef typename BaseType::IndexType IndexType;

    /**
     * This type is used to return size or dimension in
     * geometry. Dimension, WorkingDimension, PointsNumber and
     * ... return this type as their results.
     */
    typedef typename BaseType::SizeType SizeType;

    /**
     * Array of counted pointers to point.
     * This type used to hold geometry's points.
     */
    typedef  typename BaseType::PointsArrayType PointsArrayType;

    /**
     * Array of coordinates. Can be Nodes, Points or IntegrationPoints
     */
    typedef typename BaseType::CoordinatesArrayType CoordinatesArrayType;

    /**
     * This type used for representing an integration point in geometry.
     * This integration point is a point with an additional weight component.
     */
    typedef typename BaseType::IntegrationPointType IntegrationPointType;

    /**
     * A Vector of IntegrationPointType which used to hold
     * integration points related to an integration
     * method.
     * IntegrationPoints functions used this type to return
     * their results.
     */
    typedef typename BaseType::IntegrationPointsArrayType IntegrationPointsArrayType;

    /**
     * A Vector of IntegrationPointsArrayType which used to hold
     * integration points related to different integration method
     * implemented in geometry.
     */
    typedef typename BaseType::IntegrationPointsContainerType IntegrationPointsContainerType;

    /**
     * A third order tensor used as shape functions' values
     * container.
     */
    typedef typename BaseType::ShapeFunctionsValuesContainerType ShapeFunctionsValuesContainerType;

    /**
     * A fourth order tensor used as shape functions' local
     * gradients container in geometry.
     */
    typedef typename BaseType::ShapeFunctionsLocalGradientsContainerType ShapeFunctionsLocalGradientsContainerType;

    /**
     * A third order tensor to hold jacobian matrices evaluated at
     * integration points. Jacobian and InverseOfJacobian functions
     * return this type as their result.
     */
    typedef typename BaseType::JacobiansType JacobiansType;

    /**
     * A third order tensor to hold shape functions' local
     * gradients. ShapefunctionsLocalGradients function return this
     * type as its result.
     */
    typedef typename BaseType::ShapeFunctionsGradientsType ShapeFunctionsGradientsType;

    /**
     * A third order tensor to hold shape functions' local second derivatives.
     * ShapefunctionsLocalGradients function return this
     * type as its result.
     */
    typedef typename BaseType::ShapeFunctionsSecondDerivativesType
    ShapeFunctionsSecondDerivativesType;

    /**
    * A third order tensor to hold shape functions' local second derivatives.
    * ShapefunctionsLocalGradients function return this
    * type as its result.
    */
    typedef typename BaseType::ShapeFunctionsThirdDerivativesType
    ShapeFunctionsThirdDerivativesType;

    /**
     * Type of the normal vector used for normal to edges in geomety.
     */
    typedef typename BaseType::NormalType NormalType;

    ///@}
    ///@name Life Cycle
    ///@{

    Polygon( const PointsArrayType& ThisPoints )
        : BaseType( ThisPoints, &msGeometryData )
    {
        if ( this->PointsNumber() != TnVertices )
        {
            std::stringstream ss;
            ss << "Invalid points number. Expected " << TnVertices << ", given " << this->PointsNumber();
            KRATOS_THROW_ERROR( std::invalid_argument, ss.str(), "" );
        }
    }

    /**
     * Copy constructor.
     * Construct this geometry as a copy of given geometry.
     *
     * @note This copy constructor does not copy the points and new
     * geometry shares points with given source geometry. It is
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    Polygon( Polygon const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Copy constructor from a geometry with other point type.
     * Construct this geometry as a copy of given geometry which
     * has different type of points. The given goemetry's
     * TOtherPointType* must be implicity convertible to this
     * geometry PointType.
     *
     * @note This copy constructor does not copy the points and new
     * geometry shares points with given source geometry. It is
     * obvious that any change to this new geometry's point affect
     * source geometry's points too.
     */
    template<class TOtherPointType> Polygon( Polygon<TOtherPointType, TnVertices> const& rOther )
        : BaseType( rOther )
    {
    }

    /**
     * Destructor. Does nothing!!!
     */
    virtual ~Polygon() {}

    GeometryData::KratosGeometryFamily GetGeometryFamily()
    {
        return GeometryData::Kratos_Polygon;
    }

    GeometryData::KratosGeometryType GetGeometryType()
    {
        // REMARKS: naming convention for polygon: https://www.mathsisfun.com/geometry/polygons.html
        if(TnVertices == 3)
            return GeometryData::Kratos_Tritagon;
        else if(TnVertices == 4)
            return GeometryData::Kratos_Tetragon;
        else if(TnVertices == 5)
            return GeometryData::Kratos_Pentagon;
        else if(TnVertices == 6)
            return GeometryData::Kratos_Hexagon;
        else if(TnVertices == 7)
            return GeometryData::Kratos_Heptagon;
        else if(TnVertices == 8)
            return GeometryData::Kratos_Octagon;
        else if(TnVertices == 9)
            return GeometryData::Kratos_Nonagon;
        else if(TnVertices == 10)
            return GeometryData::Kratos_Decagon;
        else if(TnVertices == 11)
            return GeometryData::Kratos_Hendecagon;
        else if(TnVertices == 12)
            return GeometryData::Kratos_Dodecagon;
        else if(TnVertices == 13)
            return GeometryData::Kratos_Triskaidecagon;
        else if(TnVertices == 14)
            return GeometryData::Kratos_Tetrakaidecagon;
        else if(TnVertices == 15)
            return GeometryData::Kratos_Pentadecagon;
        else if(TnVertices == 16)
            return GeometryData::Kratos_Hexakaidecagon;
        else if(TnVertices == 17)
            return GeometryData::Kratos_Heptadecagon;
        else if(TnVertices == 18)
            return GeometryData::Kratos_Octakaidecagon;
        else if(TnVertices == 19)
            return GeometryData::Kratos_Enneadecagon;
    }

    ///@}
    ///@name Operators
    ///@{

    /**
     * Assignment operator.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    Polygon& operator=( const Polygon& rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    /**
     * Assignment operator for geometries with different point type.
     *
     * @note This operator don't copy the points and this
     * geometry shares points with given source geometry. It's
     * obvious that any change to this geometry's point affect
     * source geometry's points too.
     *
     * @see Clone
     * @see ClonePoints
     */
    template<class TOtherPointType>
    Polygon& operator=( Polygon<TOtherPointType, TnVertices> const & rOther )
    {
        BaseType::operator=( rOther );
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    typename BaseType::Pointer Create( PointsArrayType const& ThisPoints ) const
    {
        switch ( ThisPoints.size() )
        {
            case 3:
                return typename BaseType::Pointer( new Polygon<TPointType, 3>( ThisPoints ) );
            case 4:
                return typename BaseType::Pointer( new Polygon<TPointType, 4>( ThisPoints ) );
            case 5:
                return typename BaseType::Pointer( new Polygon<TPointType, 5>( ThisPoints ) );
            case 6:
                return typename BaseType::Pointer( new Polygon<TPointType, 6>( ThisPoints ) );
            case 7:
                return typename BaseType::Pointer( new Polygon<TPointType, 7>( ThisPoints ) );
            case 8:
                return typename BaseType::Pointer( new Polygon<TPointType, 8>( ThisPoints ) );
            case 9:
                return typename BaseType::Pointer( new Polygon<TPointType, 9>( ThisPoints ) );
            case 10:
                return typename BaseType::Pointer( new Polygon<TPointType, 10>( ThisPoints ) );
            case 11:
                return typename BaseType::Pointer( new Polygon<TPointType, 11>( ThisPoints ) );
            case 12:
                return typename BaseType::Pointer( new Polygon<TPointType, 12>( ThisPoints ) );
            case 13:
                return typename BaseType::Pointer( new Polygon<TPointType, 13>( ThisPoints ) );
            case 14:
                return typename BaseType::Pointer( new Polygon<TPointType, 14>( ThisPoints ) );
            case 15:
                return typename BaseType::Pointer( new Polygon<TPointType, 15>( ThisPoints ) );
            case 16:
                return typename BaseType::Pointer( new Polygon<TPointType, 16>( ThisPoints ) );
            case 17:
                return typename BaseType::Pointer( new Polygon<TPointType, 17>( ThisPoints ) );
            case 18:
                return typename BaseType::Pointer( new Polygon<TPointType, 18>( ThisPoints ) );
            case 19:
                return typename BaseType::Pointer( new Polygon<TPointType, 19>( ThisPoints ) );
            default:
                KRATOS_THROW_ERROR(std::logic_error, ThisPoints.size(), "-gon geometry can't be created")
//                return typename BaseType::Pointer( new Polygon( ThisPoints ) );
        }
    }

    virtual Geometry< Point<3> >::Pointer Clone() const
    {
        Geometry< Point<3> >::PointsArrayType NewPoints;

        //making a copy of the nodes TO POINTS (not Nodes!!!)

        for ( IndexType i = 0 ; i < this->Points().size() ; i++ )
            NewPoints.push_back( this->Points()[i] );

        //creating a geometry with the new points
        Geometry< Point<3> >::Pointer p_clone( new Polygon< Point<3>, TnVertices >( NewPoints ) );

        p_clone->ClonePoints();

        return p_clone;
    }

    /**
     * returns the local coordinates of all nodes of the current geometry
     * @param rResult a Matrix object that will be overwritten by the result
     * @return the local coordinates of all nodes
     */
    virtual Matrix& PointsLocalCoordinates( Matrix& rResult ) const
    {
        noalias( rResult ) = ZeroMatrix( TnVertices, 2 );
        // TODO
//        rResult( 0, 0 ) = -1.0;
//        rResult( 0, 1 ) = -1.0;
//        rResult( 1, 0 ) =  1.0;
//        rResult( 1, 1 ) = -1.0;
//        rResult( 2, 0 ) =  1.0;
//        rResult( 2, 1 ) =  1.0;
//        rResult( 3, 0 ) = -1.0;
//        rResult( 3, 1 ) =  1.0;
        return rResult;
    }

    //lumping factors for the calculation of the lumped mass matrix
    virtual Vector& LumpingFactors( Vector& rResult ) const
    {
        if(rResult.size() != TnVertices)
            rResult.resize( TnVertices, false );
        std::fill( rResult.begin(), rResult.end(), 1.00 / static_cast<double>(TnVertices) );
        return rResult;
    }

    ///@}
    ///@name Information
    ///@{

    /**
     * This method calculates and returns Length or charactereistic
     * length of this geometry depending on it's dimension.
     * For one dimensional geometry for example Line it returns
     * length of it and for the other geometries it gives Characteristic
     * length otherwise.
     * In the current geometry this function returns the determinant of
     * jacobian
     *
     * @return double value contains length or Characteristic
     * length
     * @see Area()
     * @see Volume()
     * @see DomainSize()
     */
    /**
     * :TODO: could be replaced by something more suitable
     * (comment by janosch)
     */
    virtual double Length() const
    {
        //return sqrt(fabs( DeterminantOfJacobian(PointType())));
        double length = 0.000;
        length = sqrt( fabs( Area() ) );
        return length;

    }

    /** This method calculates and returns area or surface area of
     * this geometry depending to it's dimension. For one dimensional
     * geometry it returns zero, for two dimensional it gives area
     * and for three dimensional geometries it gives surface area.
     *
     * @return double value contains area or surface
     * area.N
     * @see Length()
     * @see Volume()
     * @see DomainSize()
     */
    /**
     * :TODO: could be replaced by something more suitable
     * (comment by janosch)
     */
    virtual double Area() const
    {
        // TODO check
        Vector temp;
        this->DeterminantOfJacobian( temp, msGeometryData.DefaultIntegrationMethod() );
        const IntegrationPointsArrayType& integration_points = this->IntegrationPoints( msGeometryData.DefaultIntegrationMethod() );
        double Area = 0.00;

        for ( unsigned int i = 0; i < integration_points.size(); i++ )
        {
            Area += temp[i] * integration_points[i].Weight();
        }

        return Area;
    }


    /** This method calculates and returns length, area or volume of
     * this geometry depending to it's dimension. For one dimensional
     * geometry it returns its length, for two dimensional it gives area
     * and for three dimensional geometries it gives its volume.
     *
     * @return double value contains length, area or volume.
     * @see Length()
     * @see Area()
     * @see Volume()
     */
    /**
     * :TODO: could be replaced by something more suitable
     * (comment by janosch)
     */
    virtual double DomainSize() const
    {
        return fabs( this->DeterminantOfJacobian( PointType() ) ) * 0.5;
    }

    /**
     * Returns whether given arbitrary point is inside the Geometry
     */
    virtual bool IsInside( const CoordinatesArrayType& rPoint, CoordinatesArrayType& rResult )
    {
        // TODO check
        this->PointLocalCoordinates( rResult, rPoint );

//        if ( rResult[0] >= -1.0 && rResult[0] <= 1.0 )
//            if ( rResult[1] >= -1.0 && rResult[1] <= 1.0 )
//                return true;

        return false;
    }



    /** This method gives you number of all edges of this
    geometry. This method will gives you number of all the edges
    with one dimension less than this geometry. for example a
    triangle would return three or a tetrahedral would return
    four but won't return nine related to its six edge lines.

    @return SizeType containes number of this geometry edges.
    @see Edges()
    @see Edge()
    */
    virtual SizeType EdgesNumber() const
    {
        return TnVertices;
    }

    /** This method gives you all edges of this geometry. This
    method will gives you all the edges with one dimension less
    than this geometry. for example a triangle would return
    three lines as its edges or a tetrahedral would return four
    triangle as its edges but won't return its six edge
    lines by this method.

    @return GeometriesArrayType containes this geometry edges.
    @see EdgesNumber()
    @see Edge()
    */
    virtual GeometriesArrayType Edges( void )
    {
        GeometriesArrayType edges = GeometriesArrayType();
        for(unsigned int i = 0; i < TnVertices-1; ++i)
            edges.push_back( EdgeType( this->pGetPoint( i ), this->pGetPoint( i+1 ) ) );
        edges.push_back( EdgeType( this->pGetPoint( TnVertices-1 ), this->pGetPoint( 0 ) ) );
        return edges;
    }

    ///@}
    ///@name Shape Function
    ///@{

    /**
     * Calculates the value of a given shape function at a given point.
     *
     * @param ShapeFunctionIndex The number of the desired shape function
     * @param rPoint the given point in local coordinates at which the value of the shape
     * function is calculated
     *
     * @return the value of the shape function at the given point
     */
    virtual double ShapeFunctionValue( IndexType ShapeFunctionIndex,
                                       const CoordinatesArrayType& rPoint ) const
    {
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, "is not supported")
        return 0.0;
    }

    /** This method gives all non-zero shape functions values
    evaluated at the rCoordinates provided

    \note There is no control if the return vector is empty or not!

    @return Vector of values of shape functions \f$ F_{i} \f$
    where i is the shape function index (for NURBS it is the index
    of the local enumeration in the element).

    @see ShapeFunctionValue
    @see ShapeFunctionsLocalGradients
    @see ShapeFunctionLocalGradient
    */
    static Vector& ShapeFunctionsValuesImpl(Vector &rResults, const CoordinatesArrayType& rCoordinates)
    {
        if(rResults.size() != TnVertices)
            rResults.resize(TnVertices, false);

        Vector A(TnVertices+1);

        Matrix Tmp(3, 3);
        Tmp(0, 0) = rCoordinates[0];
        Tmp(0, 1) = rCoordinates[1];
        Tmp(0, 2) = 1.0;
        Tmp(1, 2) = 1.0;
        Tmp(2, 2) = 1.0;
        double x1, y1, x2, y2;

        for(unsigned int i = 0; i < TnVertices; ++i)
        {
            x1 = cos(2.0*M_PI*i/TnVertices);
            y1 = sin(2.0*M_PI*i/TnVertices);
            if(i != TnVertices-1)
            {
                x2 = cos(2.0*M_PI*(i+1)/TnVertices);
                y2 = sin(2.0*M_PI*(i+1)/TnVertices);
            }
            else
            {
                x2 = 1.0;
                y2 = 0.0;
            }
            Tmp(1, 0) = x1;
            Tmp(1, 1) = y1;
            Tmp(2, 0) = x2;
            Tmp(2, 1) = y2;

            A(i) = MathUtils<double>::Det(Tmp);
        }
        A(TnVertices) = A(0);

        Vector alpha(TnVertices);
        double sum_alpha = 0.0;
        for(unsigned int i = 0; i < TnVertices; ++i)
        {
            alpha(i) = 1.0 / (A(i)*A(i+1));
            sum_alpha += alpha(i);
        }

        for(unsigned int i = 0; i < TnVertices; ++i)
        {
            rResults(i) = alpha(i) / sum_alpha;
        }

        return rResults;
    }

    virtual Vector& ShapeFunctionsValues(Vector &rResults, const CoordinatesArrayType& rCoordinates) const
    {
        return ShapeFunctionsValuesImpl(rResults, rCoordinates);
    }

    /**
     * Calculates the gradients in terms of local coordinates
     * of all shape functions in a given point.
     *
     * @param rPoint the current point at which the gradients are calculated in local
     * coordinates
     * @return the gradients of all shape functions
     * \f$ \frac{\partial N^i}{\partial \xi_j} \f$
     */
    static Matrix& ShapeFunctionsLocalGradientsImpl( Matrix& rResult,
            const CoordinatesArrayType& rPoint )
    {
        rResult.resize( TnVertices, 2, false );
        //noalias( rResult ) = ZeroMatrix( TnVertices, 2 );

        Vector A(TnVertices+1);
        Matrix dA(TnVertices+1, 2);

        Matrix Tmp(3, 3);
        Tmp(0, 0) = rPoint[0];
        Tmp(0, 1) = rPoint[1];
        Tmp(0, 2) = 1.0;
        Tmp(1, 2) = 1.0;
        Tmp(2, 2) = 1.0;
        double x1, y1, x2, y2;

        for(unsigned int i = 0; i < TnVertices; ++i)
        {
            x1 = cos(2.0*M_PI*i/TnVertices);
            y1 = sin(2.0*M_PI*i/TnVertices);
            if(i != TnVertices-1)
            {
                x2 = cos(2.0*M_PI*(i+1)/TnVertices);
                y2 = sin(2.0*M_PI*(i+1)/TnVertices);
            }
            else
            {
                x2 = 1.0;
                y2 = 0.0;
            }
            Tmp(1, 0) = x1;
            Tmp(1, 1) = y1;
            Tmp(2, 0) = x2;
            Tmp(2, 1) = y2;

            A(i) = 0.5*MathUtils<double>::Det(Tmp);
            dA(i, 0) = 0.5*(y2 - y1);
            dA(i, 1) = 0.5*(x1 - x2);
        }
        A(TnVertices) = A(0);
        dA(TnVertices, 0) = dA(0, 0);
        dA(TnVertices, 1) = dA(0, 1);

        Vector alpha(TnVertices);
        Matrix dalpha(TnVertices, 2);
        double sum_alpha = 0.0;
        Vector sum_dalpha(2);
        sum_dalpha(0) = 0.0;
        sum_dalpha(1) = 0.0;
        for(unsigned int i = 0; i < TnVertices; ++i)
        {
            alpha(i) = 1.0 / (A(i)*A(i+1));
            dalpha(i, 0) = -alpha(i)*(dA(i, 0)/A(i) + dA(i+1, 0)/A(i+1));
            dalpha(i, 1) = -alpha(i)*(dA(i, 1)/A(i) + dA(i+1, 1)/A(i+1));
            sum_alpha += alpha(i);
            sum_dalpha(0) += dalpha(i, 0);
            sum_dalpha(1) += dalpha(i, 1);
        }

        double Ni;
        for(unsigned int i = 0; i < TnVertices; ++i)
        {
            Ni = alpha(i) / sum_alpha;
            rResult(i, 0) = (dalpha(i, 0) - Ni*sum_dalpha(0))/sum_alpha;
            rResult(i, 1) = (dalpha(i, 1) - Ni*sum_dalpha(1))/sum_alpha;
        }

        return rResult;
    }

    virtual Matrix& ShapeFunctionsLocalGradients( Matrix& rResult,
            const CoordinatesArrayType& rPoint ) const
    {
        return ShapeFunctionsLocalGradientsImpl(rResult, rPoint);
    }


    /**
     * returns the second order derivatives of all shape functions
     * in given arbitrary pointers
     * @param rResult a third order tensor which contains the second derivatives
     * @param rPoint the given point the second order derivatives are calculated in
     */
    virtual ShapeFunctionsSecondDerivativesType& ShapeFunctionsSecondDerivatives( ShapeFunctionsSecondDerivativesType& rResult, const CoordinatesArrayType& rPoint ) const
    {
        if ( rResult.size() != this->PointsNumber() )
        {
            // KLUDGE: While there is a bug in
            // ublas vector resize, I have to put this beside resizing!!
            ShapeFunctionsGradientsType temp( this->PointsNumber() );
            rResult.swap( temp );
        }

        for(unsigned int i = 0; i < this->PointsNumber(); ++i)
            rResult[i].resize( 2, 2 );

        // TODO
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, " is yet implemented");

        return rResult;
    }

    /**
     * returns the third order derivatives of all shape functions
     * in given arbitrary pointers
     * @param rResult a fourth order tensor which contains the third derivatives
     * @param rPoint the given point the third order derivatives are calculated in
     */
    virtual ShapeFunctionsThirdDerivativesType& ShapeFunctionsThirdDerivatives( ShapeFunctionsThirdDerivativesType& rResult, const CoordinatesArrayType& rPoint ) const
    {
        if ( rResult.size() != this->PointsNumber() )
        {
            // KLUDGE: While there is a bug in
            // ublas vector resize, I have to put this beside resizing!!
//                 ShapeFunctionsGradientsType
            ShapeFunctionsThirdDerivativesType temp( this->PointsNumber() );
            rResult.swap( temp );
        }

        for ( IndexType i = 0; i < rResult.size(); i++ )
        {
            boost::numeric::ublas::vector<Matrix> temp( this->PointsNumber() );
            rResult[i].swap( temp );
        }

        for ( unsigned int i = 0; i < this->PointsNumber(); i++ )
        {
            for ( unsigned int j = 0; j < 2; j++ )
            {
                rResult[i][j].resize( 2, 2 );
                noalias( rResult[i][j] ) = ZeroMatrix( 2, 2 );
            }
        }

        // TODO
        KRATOS_THROW_ERROR(std::logic_error, __FUNCTION__, " is yet implemented");

        return rResult;
    }

    ///@}
    ///@name Input and output
    ///@{

    /**
     * Turn back information as a string.
     *
     * @return String contains information about this geometry.
     * @see PrintData()
     * @see PrintInfo()
     */
    virtual std::string Info() const
    {
        std::stringstream ss;
        ss << "2 dimensional polygon with " << TnVertices << " nodes in 2D space";
        return ss.str();
    }

    /**
     * Print information about this object.
     * @param rOStream Stream to print into it.
     * @see PrintData()
     * @see Info()
     */
    virtual void PrintInfo( std::ostream& rOStream ) const
    {
        rOStream << "2 dimensional polygon with " << TnVertices << " nodes in 2D space";
    }

    /**
     * Print geometry's data into given stream.
     * Prints it's points
     * by the order they stored in the geometry and then center
     * point of geometry.
     *
     * @param rOStream Stream to print into it.
     * @see PrintInfo()
     * @see Info()
     */
    /**
     * :TODO: needs to be reviewed because it is not properly implemented yet
     * (comment by janosch)
     */
    virtual void PrintData( std::ostream& rOStream ) const
    {
        BaseType::PrintData( rOStream );
        std::cout << std::endl;
        Matrix jacobian;
        this->Jacobian( jacobian, PointType() );
        rOStream << "    Jacobian in the origin\t : " << jacobian;
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:
    /**
     * There are no protected members in class Polygon
     */

private:
    ///@name Static Member Variables
    ///@{
    static const GeometryData msGeometryData;

    ///@}
    ///@name Member Variables
    ///@{


    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save( Serializer& rSerializer ) const
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, PointsArrayType );
    }

    virtual void load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, PointsArrayType );
    }

    Polygon(): BaseType( PointsArrayType(), &msGeometryData ) {}


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * TODO: implemented but not yet tested
     */
    /**
     * Calculates the values of all shape function in all integration points.
     * Integration points are expected to be given in local coordinates
     * @param ThisMethod the current integration method
     * @return the matrix of values of every shape function in each integration point
     * :KLUDGE: number of points is hard-coded -> be careful if you want to copy and paste!
     */
    static Matrix CalculateShapeFunctionsIntegrationPointsValues(
        typename BaseType::IntegrationMethod ThisMethod )
    {
        IntegrationPointsContainerType all_integration_points =
            AllIntegrationPoints();
        IntegrationPointsArrayType integration_points =
            all_integration_points[ThisMethod];
        //number of integration points
        const int integration_points_number = integration_points.size();
        //number of nodes in current geometry
        const int points_number = TnVertices;
        //setting up return matrix
        Matrix shape_function_values( integration_points_number, points_number );
        //loop over all integration points

        Vector tmp_values;
        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            tmp_values = ShapeFunctionsValuesImpl( tmp_values, integration_points[pnt] );
            row( shape_function_values, pnt ) = tmp_values;
        }

        return shape_function_values;
    }

    /**
     * TODO: implemented but not yet tested
     */
    /**
     * Calculates the local gradients of all shape functions
     * in all integration points.
     * Integration points are expected to be given in local coordinates
     *
     * @param ThisMethod the current integration method
     * @return the vector of the gradients of all shape functions
     * in each integration point
     */
    static ShapeFunctionsGradientsType
    CalculateShapeFunctionsIntegrationPointsLocalGradients(
        typename BaseType::IntegrationMethod ThisMethod )
    {
        IntegrationPointsContainerType all_integration_points =
            AllIntegrationPoints();
        IntegrationPointsArrayType integration_points =
            all_integration_points[ThisMethod];
        //number of integration points
        const int integration_points_number = integration_points.size();
        ShapeFunctionsGradientsType d_shape_f_values( integration_points_number );
        //initialising container
        //std::fill(d_shape_f_values.begin(), d_shape_f_values.end(), Matrix(4,2));
        //loop over all integration points

        for ( int pnt = 0; pnt < integration_points_number; pnt++ )
        {
            Matrix result( TnVertices, 2 );
            noalias( result ) = ShapeFunctionsLocalGradientsImpl( result, integration_points[pnt] );
            d_shape_f_values[pnt] = result;
        }

        return d_shape_f_values;
    }

    /**
     * TODO: testing
     */
    static const IntegrationPointsContainerType AllIntegrationPoints()
    {
        IntegrationPointsContainerType integration_points =
        {
            {
                Quadrature<PolygonGaussLegendreIntegrationPoints<TnVertices, TriangleGaussLegendreIntegrationPoints1>, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<PolygonGaussLegendreIntegrationPoints<TnVertices, TriangleGaussLegendreIntegrationPoints2>, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<PolygonGaussLegendreIntegrationPoints<TnVertices, TriangleGaussLegendreIntegrationPoints3>, 2, IntegrationPoint<3> >::GenerateIntegrationPoints(),
                Quadrature<PolygonGaussLegendreIntegrationPoints<TnVertices, TriangleGaussLegendreIntegrationPoints4>, 2, IntegrationPoint<3> >::GenerateIntegrationPoints()
            }
        };
        return integration_points;
    }

    /**
     * TODO: testing
     */
    static const ShapeFunctionsValuesContainerType AllShapeFunctionsValues()
    {
        ShapeFunctionsValuesContainerType shape_functions_values =
        {
            {
                Polygon<TPointType, TnVertices>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_1 ),
                Polygon<TPointType, TnVertices>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_2 ),
                Polygon<TPointType, TnVertices>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_3 ),
                Polygon<TPointType, TnVertices>::CalculateShapeFunctionsIntegrationPointsValues( GeometryData::GI_GAUSS_4 )
            }
        };
        return shape_functions_values;
    }

    /**
     * TODO: testing
     */
    static const ShapeFunctionsLocalGradientsContainerType AllShapeFunctionsLocalGradients()
    {
        ShapeFunctionsLocalGradientsContainerType shape_functions_local_gradients =
        {
            {
                Polygon<TPointType, TnVertices>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_1 ),
                Polygon<TPointType, TnVertices>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_2 ),
                Polygon<TPointType, TnVertices>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_3 ),
                Polygon<TPointType, TnVertices>::CalculateShapeFunctionsIntegrationPointsLocalGradients( GeometryData::GI_GAUSS_4 )
            }
        };
        return shape_functions_local_gradients;
    }

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Private Friends
    ///@{

    template<class TOtherPointType, std::size_t T> friend class Polygon;

    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}
}; // Class Geometry

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{
/**
 * input stream functions
 */
template<class TPointType, std::size_t TnVertices> inline std::istream& operator >> (
    std::istream& rIStream,
    Polygon<TPointType, TnVertices>& rThis );
/**
 * output stream functions
 */
template<class TPointType, std::size_t TnVertices> inline std::ostream& operator << (
    std::ostream& rOStream,
    const Polygon<TPointType, TnVertices>& rThis )
{
    rThis.PrintInfo( rOStream );
    rOStream << std::endl;
    rThis.PrintData( rOStream );
    return rOStream;
}

///@}

template<class TPointType, std::size_t TnVertices> const
GeometryData Polygon<TPointType, TnVertices>::msGeometryData(
    2, 2, 2,
    GeometryData::GI_GAUSS_2,
    Polygon<TPointType, TnVertices>::AllIntegrationPoints(),
    Polygon<TPointType, TnVertices>::AllShapeFunctionsValues(),
    Polygon<TPointType, TnVertices>::AllShapeFunctionsLocalGradients()
);
}// namespace Kratos.

#endif // KRATOS_POLYGON_H_INCLUDED  defined 

