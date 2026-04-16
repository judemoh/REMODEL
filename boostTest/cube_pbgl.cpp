// Enable PBGL interfaces (must come before some BGL headers)
#include <boost/graph/use_mpi.hpp>

#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#include <boost/graph/distributed/adjacency_list.hpp>
#include <boost/graph/distributed/mpi_process_group.hpp>
#include <boost/graph/distributed/breadth_first_search.hpp>

#include <array>
#include <iostream>

namespace bgd = boost::graph::distributed;

// A simple vertex property (coordinates)
struct Point3d { double x, y, z; };

// Distributed undirected graph
using Graph = boost::adjacency_list<
    boost::vecS,
    boost::distributedS<bgd::mpi_process_group, boost::vecS>,
    boost::undirectedS,
    Point3d
>;

int main(int argc, char** argv) {
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;

    // Create 8 global vertices collectively (distributed across ranks)
    Graph g(8); // distributed adjacency_list supports this pattern :contentReference[oaicite:2]{index=2}

    // Set vertex coordinates on the owning rank
    for (int z = 0; z <= 1; ++z)
        for (int y = 0; y <= 1; ++y)
            for (int x = 0; x <= 1; ++x) {
                int idx = x + 2*y + 4*z;
                auto v = vertex(idx, g);

                // Only the owner should write the local property
                if (get(boost::vertex_owner, g, v) == process_id(g.process_group())) {
                    g[v] = Point3d{double(x), double(y), double(z)};
                }
            }

    // Easiest construction: add edges on rank 0, then synchronize.
    // (For big graphs, you’d partition edge insertion too.)
    if (process_id(g.process_group()) == 0) {
        for (int i = 0; i < 8; ++i) {
            for (int bit = 0; bit < 3; ++bit) {
                int j = i ^ (1 << bit); // flip x/y/z bit
                if (i < j) {
                    add_edge(vertex(i, g), vertex(j, g), g);
                }
            }
        }
    }

    // Finalize distributed structural updates
    synchronize(g.process_group()); // standard PBGL step after modifications :contentReference[oaicite:3]{index=3}

    // Run a distributed algorithm (BFS) as a sanity check
    boost::breadth_first_search(g, vertex(0, g), boost::visitor(boost::null_visitor()));

    // Print what each rank owns (local vertices iterator is per-rank in distributed graphs)
    world.barrier();
    for (int r = 0; r < world.size(); ++r) {
        world.barrier();
        if (world.rank() == r) {
            std::cout << "\n=== Rank " << r << " local vertices ===\n";
            for (auto v : boost::make_iterator_range(vertices(g))) {
                auto gi = get(boost::vertex_index, g, v); // global index
                const auto& p = g[v];
                std::cout << "global " << gi << " at (" << p.x << "," << p.y << "," << p.z << ")\n";
            }
            std::cout << std::flush;
        }
    }

    return 0;
}
