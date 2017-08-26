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
//  Date:            24 Aug 2017
//


#if !defined(KRATOS_POLYTREE_UTILITY_H_INCLUDED )
#define  KRATOS_POLYTREE_UTILITY_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>


// External includes
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/foreach.hpp>
#include <boost/progress.hpp>

#ifdef POLYFEM_USE_QHULL
#include "libqhullcpp/RboxPoints.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/Qhull.h"
#include "libqhullcpp/QhullQh.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullFacetSet.h"
#include "libqhullcpp/QhullLinkedList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullSet.h"
#include "libqhullcpp/QhullVertexSet.h"
#endif


// Project includes
#include "includes/define.h"


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
class PolyTreeUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of PolyTreeUtility
    KRATOS_CLASS_POINTER_DEFINITION(PolyTreeUtility);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    PolyTreeUtility() {}

    /// Destructor.
    virtual ~PolyTreeUtility() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /** 
     * Given a polygon, compute a decomposition of this polygon, using Voronoi tessellation and LLoyd algorithm
     * REF: H. Nguyen-Xuan et al, A polytree-based adaptive approach to limit analysis of cracked structures
     * @param compatible_coordinates        the Voronoi coordinates
     * @param voronoi_connectivities     connectivities of polygonal Voronoi cell (index-1 based)
     * @param pro_constant       proportional constant to compute reflection
     * @param max_iter           max number of iteration for LLoyd algorithm
     * @param tol                tolerance for LLoyd iteration
     * @param debug_level        debug level
     *                           debug_level = 0: no printing at all
     *                           debug_level = 1: print the error
     *                           debug_level = 2: print the error and number of points generated in the iteration
     *                           debug_level = 3: print as (2) and the matlab command to visualize the Voronoi cells
     *                           debug_level = 4: very detail printing
     */
    static void ComputePolygonDecomposition(std::vector<std::vector<double> >& compatible_coordinates,
            std::vector<std::vector<std::size_t> >& voronoi_connectivities,
            const std::vector<std::vector<double> >& polygon,
            const double& pro_constant,
            const std::size_t& max_iter,
            const double& tol,
            const int& debug_level)
    {
        // compute the center of the polygon
        std::vector<double> center(2);
        center[0] = 0.0;
        center[1] = 0.0;
        for (std::size_t i = 0; i < polygon.size(); ++i)
        {
            center[0] += polygon[i][0];
            center[1] += polygon[i][1];
        }
        center[0] /= polygon.size();
        center[1] /= polygon.size();

        if (debug_level > 3)
        {
            std::cout << "center point: " << center[0] << " " << center[1] << std::endl;
        }

        // compute the seed points
        std::vector<std::vector<double> > seed_points;
        seed_points.push_back(center);
        for (std::size_t i = 0; i < polygon.size(); ++i)
        {
            std::vector<double> point(2);
            point[0] = 0.5 * (center[0] + polygon[i][0]);
            point[1] = 0.5 * (center[1] + polygon[i][1]);
            seed_points.push_back(point);
        }

        if (debug_level > 3)
        {
            std::cout << "Initial seed points:" << std::endl;
            for (std::size_t i = 0; i < seed_points.size(); ++i)
            {
                std::cout << "  " << seed_points[i][0] << " " << seed_points[i][1] << std::endl;
            }
        }

        // LLoyd iteration
        std::size_t it = 0;
        bool converged = false;
        double error, polygon_area;
        std::vector<std::vector<double> > voronoi_coordinates;
        do
        {
            if (debug_level > 0)
                std::cout << "Iteration " << it+1 << std::endl;

            // compute the reflection points
            std::vector<std::vector<double> > reflection_points;
            polygon_area = ComputeReflectionPoints(reflection_points, seed_points, polygon, pro_constant);

            if (debug_level > 3)
            {
                std::cout << "  Reflection points:" << std::endl;
                for (std::size_t i = 0; i < reflection_points.size(); ++i)
                {
                    std::cout << "    " << reflection_points[i][0] << " " << reflection_points[i][1] << std::endl;
                }
            }

            // compute the Voronoi tessellation
            // we insert the inside seed point first in order to get first the compatible polygon
            std::vector<double> seed_coordinates(2*(seed_points.size() + reflection_points.size()));
            std::size_t cnt = 0;
            for (std::size_t i = 0; i < seed_points.size(); ++i)
            {
                seed_coordinates[2*cnt] = seed_points[i][0];
                seed_coordinates[2*cnt + 1] = seed_points[i][1];
                ++cnt;
            }
            for (std::size_t i = 0; i < reflection_points.size(); ++i)
            {
                seed_coordinates[2*cnt] = reflection_points[i][0];
                seed_coordinates[2*cnt + 1] = reflection_points[i][1];
                ++cnt;
            }

            if (debug_level > 1)
            {
                std::cout << "  Number of seed points for Voronoi tessellation: " << seed_coordinates.size()/2 << std::endl;
            }

            voronoi_coordinates.clear();
            voronoi_connectivities.clear();
            ComputeVoronoiTesselation(voronoi_coordinates, voronoi_connectivities, 2, seed_coordinates);

            if (debug_level > 3)
            {
                std::cout << "  Output of the Voronoi tessellation: " << voronoi_coordinates.size() << " vertices, "
                          << voronoi_connectivities.size() << " cells" << std::endl;

                std::cout << "  voronoi_coordinates:" << std::endl;
                for (std::size_t i = 0; i < voronoi_coordinates.size(); ++i)
                {
                    std::cout << "    " << voronoi_coordinates[i][0] << " " << voronoi_coordinates[i][1] << std::endl;
                }
                std::cout << "  voronoi cells:" << std::endl;
                for (std::size_t i = 0; i < voronoi_connectivities.size(); ++i)
                {
                    std::cout << "   ";
                    for (std::size_t j = 0; j < voronoi_connectivities[i].size(); ++j)
                    {
                        std::cout << " " << voronoi_connectivities[i][j];
                    }
                    std::cout << std::endl;
                }
            }

            // compute the center of each Voronoi cell and set them as new seed points
            // also compute the error
            voronoi_connectivities.resize(seed_points.size()); // only keep the first n+1 compatible polygons, according with the inside seed point
            error = 0.0;
            for (std::size_t i = 0; i < voronoi_connectivities.size(); ++i)
            {
                center[0] = 0.0;
                center[1] = 0.0;
                for (std::size_t j = 0; j < voronoi_connectivities[i].size(); ++j)
                {
                    center[0] += voronoi_coordinates[voronoi_connectivities[i][j]-1][0];
                    center[1] += voronoi_coordinates[voronoi_connectivities[i][j]-1][1];
                }
                center[0] /= voronoi_connectivities[i].size();
                center[1] /= voronoi_connectivities[i].size();
                error += pow(center[0] - seed_points[i][0], 2) + pow(center[1] - seed_points[i][1], 2);
                seed_points[i] = center;
            }
            error = sqrt(error / polygon_area);

            converged = (error < tol);

            if (debug_level > 0)
            {
                std::cout << "  Iteration " << it+1 << ", error = " << error << std::endl;
            }

            if (debug_level > 3)
            {
                std::cout << "  New seed points:" << std::endl;
                for (std::size_t i = 0; i < seed_points.size(); ++i)
                {
                    std::cout << "    " << seed_points[i][0] << " " << seed_points[i][1] << std::endl;
                }
            }

            if (debug_level > 1)
            {
                std::cout << "  Number of Voronoi points: " << voronoi_coordinates.size() << std::endl;
                std::cout << "  Number of Voronoi cells: " << voronoi_connectivities.size() << std::endl;
            }

            ++it;
        } while (it < max_iter && converged != true);

        if (debug_level > 3)
        {
            std::cout << "voronoi_connectivities:" << std::endl;
            for (std::size_t i = 0; i < voronoi_connectivities.size(); ++i)
            {
                std::cout << "  ";
                for (std::size_t j = 0; j < voronoi_connectivities[i].size(); ++j)
                    std::cout << " " << voronoi_connectivities[i][j];
                std::cout << std::endl;
            }
        }

        if (debug_level > 0)
        {
            if (!converged)
            {
                std::cout << "LLoyd iteration does not converge to tolerance " << tol
                          << " within " << it << " steps"
                          << ", current error = " << error
                          << std::endl;
            }
            else
            {
                std::cout << "LLoyd iteration converges after " << it << " steps"
                          << ", current error = " << error
                          << ", tolerance = " << tol
                          << std::endl;
            }
        }

        // reconstruct the Voronoi cells and vertices
        std::set<std::size_t> all_indices;
        for (std::size_t i = 0; i < voronoi_connectivities.size(); ++i)
        {
            all_indices.insert(voronoi_connectivities[i].begin(), voronoi_connectivities[i].end());
        }

        std::map<std::size_t, std::size_t> new_indices;
        std::size_t cnt = 0;
        compatible_coordinates.resize(all_indices.size());
        for (std::set<std::size_t>::iterator it = all_indices.begin(); it != all_indices.end(); ++it)
        {
            compatible_coordinates[cnt] = voronoi_coordinates[*it-1];
            new_indices[*it] = ++cnt;
        }

        std::size_t max_polygon_size = 0;
        for (std::size_t i = 0; i < voronoi_connectivities.size(); ++i)
        {
            if (voronoi_connectivities[i].size() > max_polygon_size)
                max_polygon_size = voronoi_connectivities[i].size();

            for (std::size_t j = 0; j < voronoi_connectivities[i].size(); ++j)
            {
                voronoi_connectivities[i][j] = new_indices[voronoi_connectivities[i][j]];
            }
        }

        if (debug_level > 2)
        {
            std::cout << "Matlab command:" << std::endl;
            std::cout << "Faces = [";
            for (std::size_t i = 0; i < voronoi_connectivities.size(); ++i)
            {
                for (std::size_t j = 0; j < voronoi_connectivities[i].size(); ++j)
                {
                    std::cout << " " << voronoi_connectivities[i][j];
                }
                for (std::size_t j = voronoi_connectivities[i].size(); j < max_polygon_size; ++j)
                {
                    std::cout << " NaN";
                }
                std::cout << ";" << std::endl;
            }
            std::cout << "];" << std::endl;

            std::cout << "Vertices = [";
            for (std::size_t i = 0; i < compatible_coordinates.size(); ++i)
            {
                std::cout << compatible_coordinates[i][0] << " " << compatible_coordinates[i][1] << ";" << std::endl;
            }
            std::cout << "];" << std::endl;

            std::cout << "patch('Faces',Faces,'Vertices',Vertices,'FaceColor','w');" << std::endl;
        }
    }


    /**
     * Compute the Voronoi tessellation of a list of points
     * @param voronoi_coordinates output coordinates of Voronoi vertices
     * @param connectivities      connectivities of Voronoi cell (including the Inf vertex) (index-1 based)
     *                            connectivities is arranged according to the input point
     * @param dim                 dimension of input data
     * @param coordinates         coordinates data of input points
     */
    static void ComputeVoronoiTesselation(std::vector<std::vector<double> >& voronoi_coordinates,
            std::vector<std::vector<std::size_t> >& connectivities,
            const std::size_t& dim, const std::vector<double>& coordinates)
    {
        #ifdef POLYFEM_USE_QHULL
        orgQhull::Qhull qhull;

        orgQhull::PointCoordinates points;
        points.setDimension(dim);
        points.append(coordinates);

        std::string qhullCommand = "v Qbb o"; // v: compute the Voronoi diagram
                                        // p: print the Voronoi vertices
                                        // o: print the Voronoi vertices and cells

        qhull.runQhull(points.comment().c_str(), points.dimension(), points.count(), &*points.coordinates(), qhullCommand.c_str());
        qhull.outputQhull();
        qhull.clearQhullMessage();

        orgQhull::QhullVertexList vertices = qhull.vertexList();
        orgQhull::QhullFacetList facets = qhull.facetList();

        std::map<std::size_t, std::size_t> map_facet;
        std::size_t cnt = 0;
        voronoi_coordinates.resize(facets.count());
        for (orgQhull::QhullFacetList::iterator it = facets.begin(); it != facets.end(); ++it)
        {
            if (!(*it).isGood()) continue;
            orgQhull::QhullFacet& f = *it;

            orgQhull::QhullPoint center = f.getCenter();
            double* coords = center.coordinates();
            // std::cout << "  " << center;
            voronoi_coordinates[cnt].resize(dim);
            for (std::size_t i = 0; i < dim; ++i)
                voronoi_coordinates[cnt][i] = coords[i];

            ++cnt;
            map_facet[f.id()] = cnt;
        }

        connectivities.resize(vertices.count());
        for (orgQhull::QhullVertexList::iterator it = vertices.begin(); it != vertices.end(); ++it)
        {
            if (!(*it).isValid()) continue;
            if (!(*it).neighborFacetsDefined()) continue;
            orgQhull::QhullVertex& v = *it;
            orgQhull::QhullFacetSet fs = v.neighborFacets();
            // std::cout << " " << fs.count() << ": ";
            for (orgQhull::QhullFacetSet::iterator it2 = fs.begin(); it2 != fs.end(); ++it2)
            {
    //            if (!(*it2).isGood()) continue;
                orgQhull::QhullFacet f = *it2;

                std::map<std::size_t, std::size_t>::iterator it3 = map_facet.find(f.id());

                if(it3 != map_facet.end())
                {
                    // std::cout << " " << it3->second;
                    connectivities[(*it).point().id()].push_back(it3->second);
                }
                else
                {
                    // std::cout << " 0";
                    connectivities[(*it).point().id()].push_back(0);
                }
            }
            // std::cout << std::endl;
        }
        #endif
    }

    /**
     * Given a polygon and seed points inside the polygon, the reflection points through the boundary of the polygon is computed.
     * Some constraints:
     *     + the seed points must be inside the polygon
     *     + the polygon must be convex. The points are organized counter-clockwise and are not closed.
     *     + the vertices are in 2D
     * REF: Talischi et al, PolyMesher: a general-purpose mesh generator for polygonal elements written in Matlab
     * @param reflection_points output reflection points
     * @param seed_points       seed points
     * @param polygon           polygon
     * @param pro_constant      constant c (Eq.(20))
     * @return                  area of the polygon
     */
    static double ComputeReflectionPoints(std::vector<std::vector<double> >& reflection_points,
            const std::vector<std::vector<double> >& seed_points,
            const std::vector<std::vector<double> >& polygon,
            const double& pro_constant)
    {
        reflection_points.clear();

        // compute the area of the polygon
        double area = ComputeArea(polygon);
        // KRATOS_WATCH(area)
        if (area < 0.0)
            KRATOS_THROW_ERROR(std::logic_error, "The polygon is clockwise", "")

        for (std::size_t k = 0; k < seed_points.size(); ++k)
        {
            // compute the distance from the seed point to each side of the polygon
            for (std::size_t i = 0; i < polygon.size(); ++i)
            {
                std::size_t j = (i == polygon.size()-1) ? 0 : i+1;

                // construct the line equation (ax + by + c = 0)
                // std::cout << "line (" << polygon[i][0] << ", " << polygon[i][1] << ") - ("
                //           << polygon[j][0] << ", " << polygon[j][1] << ")" << std::endl;
                double a = polygon[j][1] - polygon[i][1]; // y2 - y1
                double b = polygon[i][0] - polygon[j][0]; // x1 - x2
                double c = -a*polygon[i][0] - b*polygon[i][1]; // -ax1 - by1
                // std::cout << "line eq: " << a << "x + " << b << "y + " << c << " = 0" << std::endl;

                double aux1 = a*seed_points[k][0] + b*seed_points[k][1] + c;
                double aux2 = a*a + b*b;
                double d = fabs(aux1) / sqrt(aux2);
                // KRATOS_WATCH(d)
                if (d < pro_constant*sqrt(area))
                {
                    // create the reflection point
                    // REF: https://math.stackexchange.com/questions/1013230/how-to-find-coordinates-of-reflected-point
                    std::vector<double> point(2);
                    point[0] = seed_points[k][0] - 2.0*a*aux1/aux2;
                    point[1] = seed_points[k][1] - 2.0*b*aux1/aux2;
                    reflection_points.push_back(point);
                }
            }
        }

        return area;
    }

    /**
     * Test the Voronoi tessellation
     */
    void TestComputeVoronoiTesselation() const
    {
        std::cout << "<<<<<<<<<<TestComputeVoronoiTesselation" << std::endl;

        std::vector<double> data = {0.5, 0.0,
                                0.0, 0.5,
                                -0.5, -0.5,
                                -0.2, -0.1,
                                -0.1, 0.1,
                                0.1, -0.1,
                                0.1, 0.1};

        std::vector<std::vector<double> > coordinates;
        std::vector<std::vector<std::size_t> > connectivities;

        ComputeVoronoiTesselation(coordinates, connectivities, 2, data);

        std::cout << "coordinates:" << std::endl;
        std::cout << "  0: Inf Inf" << std::endl;
        for (std::size_t i = 0; i < coordinates.size(); ++i)
        {
            std::cout << "  " << i+1 << ": " << coordinates[i][0] << " " << coordinates[i][1] << std::endl;
        }

        std::cout << "connectivities:" << std::endl;
        for (std::size_t i = 0; i < connectivities.size(); ++i)
        {
            std::cout << " ";
            for (std::size_t j = 0; j < connectivities[i].size(); ++j)
            {
                std::cout << " " << connectivities[i][j];
            }
            std::cout << std::endl;
        }
    }

    /**
     * Test compute the reflection points
     */
    void TestComputeReflectionPoints() const
    {
        std::cout << "<<<<<<<<<<TestComputeReflectionPoints" << std::endl;

        std::vector<std::vector<double> > polygon;
        polygon.push_back(std::vector<double>{0.0, 0.0});
        polygon.push_back(std::vector<double>{1.0, 0.0});
        polygon.push_back(std::vector<double>{2.0, 2.0});
        polygon.push_back(std::vector<double>{1.0, 2.0});
        polygon.push_back(std::vector<double>{0.0, 1.0});

        std::vector<std::vector<double> > seed_points;
        seed_points.push_back(std::vector<double>{1.8369e-01, 7.0443e-01});

        std::vector<std::vector<double> > reflection_points;
        ComputeReflectionPoints(reflection_points, seed_points, polygon, 1.5);

        std::cout << "reflection_points:" << std::endl;
        for (std::size_t i = 0; i < reflection_points.size(); ++i)
        {
            std::cout << "  " << reflection_points[i][0] << " " << reflection_points[i][1] << std::endl;
        }
    }

    /**
     * Test compute the polygon decomposition
     */
    void TestPolygonDecomposition(const std::size_t& max_iter, const double& tol, const int& debug_level) const
    {
        std::cout << "<<<<<<<<<<TestPolygonDecomposition" << std::endl;

        std::vector<std::vector<double> > polygon;
        polygon.push_back(std::vector<double>{0.0, 0.0});
        polygon.push_back(std::vector<double>{1.0, 0.0});
        polygon.push_back(std::vector<double>{2.0, 2.0});
        polygon.push_back(std::vector<double>{1.0, 2.0});
        polygon.push_back(std::vector<double>{0.0, 1.0});

        std::vector<std::vector<double> > voronoi_coordinates;
        std::vector<std::vector<std::size_t> > compatible_connectivities;
        ComputePolygonDecomposition(voronoi_coordinates, compatible_connectivities, polygon, 1.5, max_iter, tol, debug_level);

        std::cout << "voronoi coordinates:" << std::endl;
        for (std::size_t i = 0; i < voronoi_coordinates.size(); ++i)
        {
            std::cout << "  " << i+1 << ": " << voronoi_coordinates[i][0] << " " << voronoi_coordinates[i][1] << std::endl;
        }

        std::cout << "compatible_connectivities:" << std::endl;
        for (std::size_t i = 0; i < compatible_connectivities.size(); ++i)
        {
            std::cout << " ";
            for (std::size_t j = 0; j < compatible_connectivities[i].size(); ++j)
            {
                std::cout << " " << compatible_connectivities[i][j];
            }
            std::cout << std::endl;
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
        return "PolyTree Utility";
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

    /**
     * Check is a point is inside a polygon
     * The polygon must be convex
     * REF: https://stackoverflow.com/questions/217578/how-can-i-determine-whether-a-2d-point-is-within-a-polygon
     * @param  point   the input point
     * @param  polygon the input polygon
     * @return         true if inside
     */
    static bool IsInside(const std::vector<double>& point, const std::vector<std::vector<double> >& polygon)
    {
        std::size_t i, j;
        bool c = false;
        for (i = 0, j = polygon.size()-1; i < polygon.size(); j = i++)
        {
            if ( ((polygon[i][1]>point[1]) != (polygon[j][1]>point[1])) && (point[0] < (polygon[j][0]-polygon[i][0]) * (point[1]-polygon[i][1]) / (polygon[j][1]-polygon[i][1]) + polygon[i][0]) )
               c = !c;
        }
        return c;
    }

    /**
     * Compute the (signed) area of the polygon
     * // REF: http://geomalgorithms.com/a01-_area.html#2D%20Polygons
     * @param  polygon the input polygon
     * @return         the signed area
     */
    static double ComputeArea(const std::vector<std::vector<double> >& polygon)
    {
        double area = 0.0;
        for (std::size_t i = 0; i < polygon.size(); ++i)
        {
            std::size_t j = (i == polygon.size()-1) ? 0 : i+1;
            area += 0.5*(polygon[i][0]*polygon[j][1] - polygon[j][0]*polygon[i][1]);
        }
        return area;
    }

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
    PolyTreeUtility& operator=(PolyTreeUtility const& rOther);

    /// Copy constructor.
    PolyTreeUtility(PolyTreeUtility const& rOther);


    ///@}

}; // Class PolyTreeUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, PolyTreeUtility& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const PolyTreeUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.


#endif // KRATOS_POLYTREE_UTILITY_H_INCLUDED  defined
