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

    struct Comparator
    {
        /**
         * Comparator for pair comparison
         * @param  l first pair
         * @param  r second pair
         * @return   true if l.first < r.first, otherwise false
         */
        inline bool operator() (const std::pair<double, std::size_t>& l, const std::pair<double, std::size_t>& r)
        {
            return l.first < r.first;
        }
    };

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
     * @param compatible_coordinates    the Voronoi coordinates
     * @param voronoi_connectivities    connectivities of polygonal Voronoi cell (index-1 based)
     * @param polygon                   the input polygon
     * @param pro_constant      proportional constant to compute reflection
     * @param max_iter          max number of iteration for LLoyd algorithm
     * @param tol               tolerance for LLoyd iteration
     * @param debug_level       debug level
     *                          debug_level = 0: no printing at all
     *                          debug_level = 1: print the error
     *                          debug_level = 2: print the error and number of points generated in the iteration
     *                          debug_level = 3: print as (2) and the matlab command to visualize the Voronoi cells
     *                          debug_level = 4: very detail printing
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

        if (debug_level > 1)
        {
            std::cout << "  Number of vertices on the polygon: " << polygon.size() << std::endl;
        }

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

            if (debug_level > 1)
            {
                std::cout << "  Number of reflection points: " << reflection_points.size() << std::endl;
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
            if (debug_level > 0)
                std::cout << "   ComputeVoronoiTesselation completed" << std::endl;

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

        // compute the signed area to rotate if needed
        std::vector<double> area;
        ComputeArea(area, voronoi_coordinates, voronoi_connectivities);
        if (debug_level > 3)
        {
            std::cout << "Area:" << std::endl;
            for (std::size_t i = 0; i < area.size(); ++i)
                std::cout << (i+1) << ": " << area[i] << std::endl;
        }

        for (std::size_t i = 0; i < area.size(); ++i)
        {
            if (area[i] < 0.0)
            {
                std::reverse(voronoi_connectivities[i].begin(), voronoi_connectivities[i].end());
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
     * Encode the vertex information on the polygon
     * @param stat        -1: outer;
     *                    0: inner;
     *                    1: on vertex;
     *                    2: on edge (inner)
     *                    3: on edge (outer)
     * @param vertex_info if stat=1, the vertex index;
     *                    stat=2, the edge index;
     *                    otherwise -1
     * @param edge_info   if stat=2, the index of the point on edge (start with 1);
     *                    otherwise -1
     * @param coordinates the input coordinates of the vertex need to encode
     * @param polygon     the input polygon
     * @param check_flag  check inside flag;
     *                    if 0, check as normal;
     *                    if 1, all the exceptional points (not on vertex and edge) are understood as inside;
     *                    if 2, outside.
     * @param tol         the tolerance for detection
     */
    static void EncodeVertices(std::vector<int>& stat,
            std::vector<int>& vertex_info,
            std::vector<int>& edge_info,
            const std::vector<std::vector<double> >& coordinates,
            const std::vector<std::vector<double> >& polygon,
            const int& check_flag,
            const double& tol)
    {
        stat.resize(coordinates.size());
        vertex_info.resize(coordinates.size());
        edge_info.resize(coordinates.size());

        std::fill(stat.begin(), stat.end(), -2);
        std::fill(vertex_info.begin(), vertex_info.end(), -1);
        std::fill(edge_info.begin(), edge_info.end(), -1);

        // check for point on vertex of the polygon
        for (std::size_t i = 0; i < coordinates.size(); ++i)
        {
            for (std::size_t j = 0; j < polygon.size(); ++j)
            {
                double dist = sqrt(pow(coordinates[i][0] - polygon[j][0], 2) + pow(coordinates[i][1] - polygon[j][1], 2));
                if (dist < tol)
                {
                    stat[i] = 1;
                    vertex_info[i] = j;
                    edge_info[i] = -1;
                }
            }
            // std::cout << "point " << (i+1) << ": " << stat[i] << std::endl;
        }

        // check for point on edge of the polygon
        std::vector<std::vector<std::size_t> > edge_with_points(polygon.size());
        for (std::size_t i = 0; i < coordinates.size(); ++i)
        {
            if (stat[i] == -2)
            {
                for (std::size_t j = 0; j < polygon.size(); ++j)
                {
                    std::size_t k = (j == polygon.size()-1) ? 0 : j+1;

                    // construct the line equation (ax + by + c = 0)
                    double a = polygon[k][1] - polygon[j][1]; // y2 - y1
                    double b = polygon[j][0] - polygon[k][0]; // x1 - x2
                    double c = -a*polygon[j][0] - b*polygon[j][1]; // -ax1 - by1
                    double length = sqrt(a*a + b*b);

                    double dist = fabs(a*coordinates[i][0] + b*coordinates[i][1] + c) / length;
                    // std::cout << "dist to edge " << j << ": " << dist << std::endl;
                    if (dist < tol)
                    {
                        if (sqrt(pow(coordinates[i][0] - polygon[j][0], 2) + pow(coordinates[i][1] - polygon[j][1], 2)) < length)
                        {
                            stat[i] = 2;
                            vertex_info[i] = j;
                            edge_with_points[j].push_back(i);
                        }
                        else
                        {
                            stat[i] = 3;
                            vertex_info[i] = j;
                        }
                    }
                }
            }
            // std::cout << "point " << (i+1) << ": " << stat[i] << std::endl;
            // std::cout << "point " << (i+1) << " vertex_info: " << vertex_info[i] << std::endl;
        }

        // check for index of the point on edge
        for (std::size_t i = 0; i < polygon.size(); ++i)
        {
            // std::cout << "edge_with_points[" << i << "].size(): " << edge_with_points[i].size() << std::endl;
            if (edge_with_points[i].size() != 0)
            {
                std::vector<std::pair<double, std::size_t> > dist_list(edge_with_points[i].size());
                for (std::size_t j = 0; j < edge_with_points[i].size(); ++j)
                {
                    double dist = sqrt(pow(coordinates[edge_with_points[i][j]][0] - polygon[i][0], 2) + pow(coordinates[edge_with_points[i][j]][1] - polygon[i][1], 2));
                    dist_list[j] = std::pair<double, std::size_t>(dist, edge_with_points[i][j]);
                }

                std::sort(dist_list.begin(), dist_list.end(), Comparator());

                for (std::size_t j = 0; j < dist_list.size(); ++j)
                {
                    edge_info[dist_list[j].second] = j+1;
                }
            }
        }

        // check for inside/outside of the polygon
        if (check_flag == 0)
        {
            for (std::size_t i = 0; i < coordinates.size(); ++i)
            {
                if (IsInside(coordinates[i], polygon))
                    stat[i] = 0;
                else
                    stat[i] = -1;
            }
        }
        else if (check_flag == 1)
        {
            for (std::size_t i = 0; i < coordinates.size(); ++i)
            {
                if (stat[i] == -2)
                    stat[i] = 0;
            }
        }
        else if (check_flag == 2)
        {
            for (std::size_t i = 0; i < coordinates.size(); ++i)
            {
                if (stat[i] == -2)
                    stat[i] = -1;
            }
        }
    }

    /**
     * Given a list of arranged points on line, find the clustering of the points
     * @param cluster the output point cluster
     * @param points  the input points (organized, e.g. on a line)
     * @param alpha   the cluster deciding parameter
     */
    static void ClusterPoints(std::vector<std::vector<std::size_t> >& cluster,
            const std::vector<std::vector<double> >& points, const double& alpha)
    {
        // compute the largest distance
        double length = sqrt(pow(points.back()[0] - points.front()[0], 2) + pow(points.back()[1] - points.front()[1], 2));

        // start to search the cluster from left to right
        cluster.clear();
        std::size_t this_point = 0;
        do
        {
            if (this_point >= points.size() - 1)
                break;

            std::vector<size_t> group;
            group.push_back(this_point);
            std::size_t next_point = this_point + 1;
            do
            {
                if (this_point >= points.size() - 1)
                    break;

                double dist = sqrt(pow(points[next_point][0] - points[this_point][0], 2) + pow(points[next_point][1] - points[this_point][1], 2));
                if (dist <= alpha*length)
                {
                    group.push_back(next_point);
                    ++this_point;
                    ++next_point;
                }
                else
                    break;
            } while(true);

            if (group.size() > 1)
                cluster.push_back(group);

            ++this_point;
        } while (true);
    }

    /**
     * Given a list of lengths (not needed to be organized), compute the clustering of the length according to the criteria
     * @param cluster output cluster
     * @param sorted_indices sorted indices based on length length
     * @param lengths list of lengths
     * @param alpha   the clustering parameter
     */
    static void ClusterLengths(std::vector<std::vector<std::size_t> >& cluster,
            std::vector<std::size_t>& sorted_indices,
            const std::vector<double>& lengths,
            const double& alpha)
    {
        // firstly sort the distance
        typedef std::vector<std::pair<double, std::size_t> > container_t;
        container_t sorted_lengths;
        for (std::size_t i = 0; i < lengths.size(); ++i)
            sorted_lengths.push_back(std::pair<double, std::size_t>{lengths[i], i});
        std::sort(sorted_lengths.begin(), sorted_lengths.end(), Comparator());

        sorted_indices.resize(lengths.size());
        for (std::size_t i = 0; i < lengths.size(); ++i)
            sorted_indices[i] = sorted_lengths[i].second;

        std::size_t first_index = sorted_indices.front();
        std::size_t last_index = sorted_indices.back();

        // compute the largest distance
        double length = sorted_lengths.back().first - sorted_lengths.front().first;

        cluster.clear();
        container_t::iterator this_point = sorted_lengths.begin();
        do
        {
            if (this_point == sorted_lengths.end())
                break;

            std::vector<size_t> group;

            // group shall not contain element from two ends
            if (this_point->second != first_index && this_point->second != last_index)
                group.push_back(this_point->second);
            container_t::iterator next_point = this_point+1;
            do
            {
                if (next_point == sorted_lengths.end())
                    break;

                double dist = next_point->first - this_point->first;
                if (dist <= alpha*length)
                {
                    // group shall not contain element from two ends
                    if (next_point->second != first_index && next_point->second != last_index)
                        group.push_back(next_point->second);
                    ++this_point;
                    ++next_point;
                }
                else
                    break;
            } while(true);

            if (group.size() > 1)
                cluster.push_back(group);

            ++this_point;
        } while (true);
    }

    /**
     * Given the indices and sub-indices, find out which edge in indices the edge in sub-indices belong to
     * Example:
     *   Let say indices = (19 53 20) and sub-indices = (19 53 36 20)
     *   then within sub-indices:
     *       + edge 0 = (19 53) belongs to edge 0 = (19 53) of indices
     *       + edge 1 = (53 36) belongs to edge 1 = (53 20) of indices
     *       + edge 1 = (36 20) belongs to edge 1 = (53 20) of indices
     * @param which_edge  the edge index, which_edge[new_edge] = original_edge
     * @param sub_indices the indices of sub-edges (incrementally)
     * @param indices     the indices of original edges
     */
    static void FindEdge(std::vector<std::size_t>& which_edge,
            const std::vector<std::size_t>& sub_indices,
            const std::vector<std::size_t>& indices)
    {
        std::map<std::size_t, int> index_code;

        for (std::size_t i = 0; i < sub_indices.size(); ++i)
            index_code[sub_indices[i]] = -1;

        for (std::size_t i = 0; i < indices.size(); ++i)
        {
            index_code[indices[i]] = i;
        }

        for (std::size_t i = 0; i < sub_indices.size(); ++i)
        {
            if (index_code[sub_indices[i]] == -1)
                index_code[sub_indices[i]] = index_code[sub_indices[i-1]];
        }

        which_edge.resize(sub_indices.size()-1);
        for (std::size_t i = 0; i < sub_indices.size()-1; ++i)
        {
            int code1 = index_code[sub_indices[i]];
            // int code2 = index_code[sub-indices[i+1]];

            which_edge[i] = code1;
        }
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

        std::cout << "Voronoi coordinates:" << std::endl;
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

        // encode the vertices in the appropriate way
        std::vector<int> stat;
        std::vector<int> vertex_info;
        std::vector<int> edge_info;
        double detect_tol = 1.0e-6; // TODO parameterize this
        PolyTreeUtility::EncodeVertices(stat, vertex_info, edge_info, voronoi_coordinates, polygon, 1, detect_tol);

        std::cout << "stat:" << std::endl;
        for (std::size_t i = 0; i < voronoi_coordinates.size(); ++i)
        {
            std::cout << "  " << (i+1) << ": " << stat[i];

            if (stat[i] == 0)
                std::cout << ", is inner";
            else if (stat[i] == 1)
                std::cout << ", on vertex " << vertex_info[i];
            else if (stat[i] == 2)
                std::cout << ", on edge " << vertex_info[i] << ", index " << edge_info[i];

            std::cout << std::endl;
        }
    }

    void TestClusterPoints1()
    {
        std::vector<std::vector<double> > points;
        points.push_back(std::vector<double>{0.0, 0.0});
        points.push_back(std::vector<double>{1.0, 0.0});
        points.push_back(std::vector<double>{2.0, 0.0});
        points.push_back(std::vector<double>{2.1, 0.0});
        points.push_back(std::vector<double>{3.0, 0.0});
        points.push_back(std::vector<double>{4.0, 0.0});

        double alpha = 0.1;

        std::vector<std::vector<std::size_t> > cluster;
        ClusterPoints(cluster, points, alpha);

        std::cout << "cluster:" << std::endl;
        for (std::size_t i = 0; i < cluster.size(); ++i)
        {
            std::cout << "  " << (i+1) << ":";
            for (std::size_t j = 0; j < cluster[i].size(); ++j)
                std::cout << " " << cluster[i][j];
            std::cout << std::endl;
        }
    }

    void TestClusterPoints2()
    {
        std::vector<std::vector<double> > points;
        points.push_back(std::vector<double>{0.0, 0.0}); //0
        points.push_back(std::vector<double>{1.0, 0.0}); //1
        points.push_back(std::vector<double>{2.0, 0.0}); //2
        points.push_back(std::vector<double>{2.1, 0.0}); //3
        points.push_back(std::vector<double>{2.2, 0.0}); //4
        points.push_back(std::vector<double>{3.0, 0.0}); //5
        points.push_back(std::vector<double>{3.1, 0.0}); //6
        points.push_back(std::vector<double>{4.0, 0.0}); //7

        double alpha = 0.1;

        std::vector<std::vector<std::size_t> > cluster;
        ClusterPoints(cluster, points, alpha);

        std::cout << "cluster:" << std::endl;
        for (std::size_t i = 0; i < cluster.size(); ++i)
        {
            std::cout << "  " << (i+1) << ":";
            for (std::size_t j = 0; j < cluster[i].size(); ++j)
                std::cout << " " << cluster[i][j];
            std::cout << std::endl;
        }
    }

    void TestClusterLengths1()
    {
        std::vector<double> lengths;
        lengths.push_back(0.0); //0
        lengths.push_back(1.1); //1
        lengths.push_back(2.0); //2
        lengths.push_back(1.0); //3
        lengths.push_back(1.2); //4
        lengths.push_back(3.0); //5
        lengths.push_back(4.0); //6
        lengths.push_back(4.1); //7

        double alpha = 0.1;

        std::vector<std::vector<std::size_t> > cluster;
        std::vector<std::size_t> sorted_indices;
        ClusterLengths(cluster, sorted_indices, lengths, alpha);

        std::cout << "cluster:" << std::endl;
        for (std::size_t i = 0; i < cluster.size(); ++i)
        {
            std::cout << "  " << (i+1) << ":";
            for (std::size_t j = 0; j < cluster[i].size(); ++j)
                std::cout << " " << cluster[i][j];
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

    /**
     * Compute the (signed) area of the list of polygons
     * // REF: http://geomalgorithms.com/a01-_area.html#2D%20Polygons
     * @param  polygon the input polygon
     * @return         the signed area
     */
    static void ComputeArea(std::vector<double>& area,
            const std::vector<std::vector<double> >& coordinates,
            const std::vector<std::vector<std::size_t> >& connectivities)
    {
        area.resize(connectivities.size());
        for (std::size_t i = 0; i < connectivities.size(); ++i)
        {
            area[i] = 0.0;
            for (std::size_t j = 0; j < connectivities[i].size(); ++j)
            {
                std::size_t k = (j == connectivities[i].size()-1) ? 0 : j+1;
                area[i] += 0.5*(coordinates[connectivities[i][j]-1][0]*coordinates[connectivities[i][k]-1][1]
                              - coordinates[connectivities[i][k]-1][0]*coordinates[connectivities[i][j]-1][1]);
            }
        }
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
