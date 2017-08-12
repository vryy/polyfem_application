//
//   Project Name:        Kratos
//   Last Modified by:    $Author:   hbui $
//   Date:                $Date:  11 Aug 2017 $
//   Revision:            $Revision:         1.0 $
//
//

#if !defined(KRATOS_POLYGON_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED )
#define  KRATOS_POLYGON_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED


// System includes
#include <cmath>

// External includes

// Project includes
#include "integration/quadrature.h"
#include "integration/triangle_gauss_legendre_integration_points.h"
#include "utilities/math_utils.h"


namespace Kratos
{

template<std::size_t TnVertices, class TriangleIntegrationPoints>
class KRATOS_API(KRATOS_CORE) PolygonGaussLegendreIntegrationPoints
{
public:
    KRATOS_CLASS_POINTER_DEFINITION(PolygonGaussLegendreIntegrationPoints);
    typedef std::size_t SizeType;

    static const unsigned int Dimension = 2;

    typedef IntegrationPoint<2> IntegrationPointType;

    typedef boost::array<IntegrationPointType, TnVertices * TriangleIntegrationPoints::IntegrationPointsNumber()> IntegrationPointsArrayType;

    typedef IntegrationPointType::PointType PointType;

    static SizeType IntegrationPointsNumber()
    {
        return TnVertices * TriangleIntegrationPoints::IntegrationPointsNumber();
    }

    static IntegrationPointsArrayType& IntegrationPoints()
    {
        msIntegrationPoints = CreateIntegrationPoints();
        return msIntegrationPoints;
    }

    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Polygon Gauss-Legendre quadrature for " << TnVertices << "-gon using " << TriangleIntegrationPoints::Info();
        return buffer.str();
    }

     static constexpr IntegrationPointsArrayType CreateIntegrationPoints()
     {
         IntegrationPointsArrayType ThisIntegrationPoints;

         typename TriangleIntegrationPoints::IntegrationPointsArrayType SampleIntegrationPoints =
             TriangleIntegrationPoints::IntegrationPoints();

         double x1, x2, y1, y2;
         double N1, N2, N3;
         double dN1dXi, dN1dEta, dN2dXi, dN2dEta, dN3dXi, dN3dEta;
         Matrix Jac(2, 2);

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

             for(unsigned int j = 0; j < SampleIntegrationPoints.size(); ++j)
             {
                 TriShape(N1, N2, N3, SampleIntegrationPoints[j]);
                 TriGradient(dN1dXi, dN1dEta, dN2dXi, dN2dEta, dN3dXi, dN3dEta, SampleIntegrationPoints[j]);

                 IntegrationPointType IntegrationPoint;

                 IntegrationPoint[0] = N1*0.0 + N2*x1 + N3*x2;
                 IntegrationPoint[1] = N1*0.0 + N2*y1 + N3*y2;
                 IntegrationPoint[2] = 0.0;

                 Jac(0, 0) = dN1dXi*0.0 + dN2dXi*x1 + dN3dXi*x2;
                 Jac(0, 1) = dN1dEta*0.0 + dN2dEta*x1 + dN3dEta*x2;
                 Jac(1, 0) = dN1dXi*0.0 + dN2dXi*y1 + dN3dXi*y2;
                 Jac(1, 1) = dN1dEta*0.0 + dN2dEta*y1 + dN3dEta*y2;

                 IntegrationPoint.Weight() = MathUtils<double>::Det(Jac) * SampleIntegrationPoints[j].Weight();

                 ThisIntegrationPoints[i*SampleIntegrationPoints.size() + j] = IntegrationPoint;
             }
         }

         return ThisIntegrationPoints;
     }

protected:

private:

    static IntegrationPointsArrayType msIntegrationPoints;    

    static void TriShape(double& N1, double& N2, double& N3, const IntegrationPointType& rPoint)
    {
        N1 = 1.0 - rPoint[0] - rPoint[1];
        N2 = rPoint[0];
        N3 = rPoint[1];
    }

    static void TriGradient(double& dN1dXi, double& dN1dEta,
                            double& dN2dXi, double& dN2dEta,
                            double& dN3dXi, double& dN3dEta,
                            const IntegrationPointType& rPoint)
    {
        dN1dXi = -1.0; dN1dEta = -1.0;
        dN2dXi = 1.0; dN2dEta = 0.0;
        dN3dXi = 0.0; dN3dEta = 1.0;
    }

}; // Class PolygonGaussLegendreIntegrationPoints

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}


}  // namespace Kratos.

#endif // KRATOS_POLYGON_GAUSS_LEGENDRE_INTEGRATION_POINTS_H_INCLUDED  defined 


