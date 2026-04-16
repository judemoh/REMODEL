#include <boost/graph/adjacency_list.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/graph/graphviz.hpp>

#include <fstream>
#include <array>
#include <iostream>

namespace bg = boost::geometry;

using Point3 = bg::model::point<double, 3, bg::cs::cartesian>;

using Graph = boost::adjacency_list<
    boost::vecS,
    boost::vecS,
    boost::undirectedS,
    Point3
>;

// -------- VTK writer function (OUTSIDE main) --------
void write_vtk_lines(const Graph& g, const std::string& path)
{
    std::ofstream out(path);

    out << "# vtk DataFile Version 3.0\n";
    out << "Cube graph\n";
    out << "ASCII\n";
    out << "DATASET POLYDATA\n";

    // Write points
    out << "POINTS " << boost::num_vertices(g) << " float\n";
    for (auto v : boost::make_iterator_range(boost::vertices(g))) {
        const auto& p = g[v];
        out << bg::get<0>(p) << " "
            << bg::get<1>(p) << " "
            << bg::get<2>(p) << "\n";
    }

    // Write edges as lines
    std::size_t E = boost::num_edges(g);
    out << "LINES " << E << " " << E * 3 << "\n";

    for (auto e : boost::make_iterator_range(boost::edges(g))) {
        auto u = boost::source(e, g);
        auto v = boost::target(e, g);
        out << "2 " << u << " " << v << "\n";
    }
}

int main()
{
    Graph g;

    // Build cube vertices
    std::array<Point3, 8> cube_pts;
    for (int z = 0; z <= 1; ++z)
        for (int y = 0; y <= 1; ++y)
            for (int x = 0; x <= 1; ++x) {
                int idx = x + 2*y + 4*z;
                cube_pts[idx] = Point3(x,y,z);
            }

    for (int i = 0; i < 8; ++i)
        boost::add_vertex(cube_pts[i], g);

    // Add cube edges
    for (int i = 0; i < 8; ++i) {
        for (int bit = 0; bit < 3; ++bit) {
            int j = i ^ (1 << bit);
            if (i < j)
                boost::add_edge(i, j, g);
        }
    }

    // -------- Call function INSIDE main --------
    write_vtk_lines(g, "cube.vtk");

    std::cout << "cube.vtk written\n";

    //----------plot graph--------------------//
    std::ofstream dotfile("cube.dot");
    boost::write_graphviz(dotfile, g);

}
