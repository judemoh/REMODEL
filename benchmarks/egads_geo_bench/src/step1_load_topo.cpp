// MPI gives us rank/size and is what we will use later for parallel runs.
// Even in "serial" tests, keeping MPI_Init/MPI_Finalize makes it easy to
// scale the same code to multiple ranks.
#include <mpi.h>

// Standard C headers for printing and exiting
#include <cstdio>
#include <cstdlib>

// C++ convenience type for filenames
#include <string>

// EGADS / EGADSlite C API header.
// This defines types like `ego` and functions like EG_open, EG_loadModel, etc.
#include "egads.h"


// ------------------------------
// Small helper: check EGADS return codes
// ------------------------------
// Almost every EGADS call returns an integer status code.
// EGADS_SUCCESS means OK; any other value is an error.
// We use this helper so code stays readable.
static void eg_check(int stat, const char* where) {
  if (stat != EGADS_SUCCESS) {
    std::fprintf(stderr, "EGADS error %d at %s\n", stat, where);
    std::exit(stat);
  }
}


// ------------------------------
// main: program entry point
// ------------------------------
int main(int argc, char** argv) {

  // 1) Start MPI
  // MPI_Init must come before MPI calls like MPI_Comm_rank.
  MPI_Init(&argc, &argv);

  // Each process has an integer "rank" (0..size-1).
  // Rank 0 is usually used for printing "global" output.
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // 2) Parse command-line args
  // We expect: ./step1 cube.egads
  if (argc < 2) {
    if (rank == 0) {
      std::fprintf(stderr, "Usage: %s cube.egads\n", argv[0]);
    }
    MPI_Finalize();
    return 1;
  }

  // Filename of the serialized EGADS model
  const char* file = argv[1];


  // 3) Create an EGADS/EGADSlite "context"
  //
  // In EGADS, a Context is like a handle to the geometry engine instance.
  // Most objects you create/load belong to a specific context.
  // You must call EG_open first.
  ego ctx = nullptr;

  // The "model" is the top-level container we load from disk.
  ego model = nullptr;

  // Time the load. MPI_Wtime gives wall-clock seconds (double).
  double t0 = MPI_Wtime();

  // Open the context
  eg_check(EG_open(&ctx), "EG_open");

  // Load the model from file (here: cube.egads)
  // bFlag=0 means default behavior.
  eg_check(EG_loadModel(ctx, 0, file, &model), "EG_loadModel");

  double tLoad = MPI_Wtime() - t0;


  // 4) Extract the model topology
  //
  // EG_getTopology returns:
  // - geom: geometry reference (often NULL/unused at this level)
  // - oclass/mtype: what kind of object this is
  // - limits: bounding box / range info depending on object class
  // - children: array of child objects
  // For a MODEL, the "children" are usually BODIES.
  ego geom = nullptr;
  ego* bodies = nullptr;   // will be allocated by EGADS
  int oclass = 0, mtype = 0, nchild = 0;
  int* senses = nullptr;   // senses correspond to children relationship
  double limits[6] = {0,0,0,0,0,0};

  eg_check(
    EG_getTopology(model, &geom, &oclass, &mtype,
                   limits, &nchild, &bodies, &senses),
    "EG_getTopology(model)"
  );

  // We expect at least one body in a valid model
  if (nchild < 1) {
    if (rank == 0) std::fprintf(stderr, "Model has no bodies!\n");
    EG_free(bodies);
    EG_close(ctx);
    MPI_Finalize();
    return 2;
  }

  // For this simple benchmark, we just take the first body.
  ego body = bodies[0];


  // 5) For the body, ask EGADS for its sub-topology:
  // faces, edges, and nodes (vertices).
  //
  // EG_getBodyTopos(body, src, oclass, &n, &objs)
  // - src = NULL means "all of them"
  // - oclass chooses FACE/EDGE/NODE
  int nfaces = 0, nedges = 0, nnodes = 0;
  ego *faces = nullptr, *edges = nullptr, *nodes = nullptr;

  eg_check(EG_getBodyTopos(body, nullptr, FACE, &nfaces, &faces), "EG_getBodyTopos(FACE)");
  eg_check(EG_getBodyTopos(body, nullptr, EDGE, &nedges, &edges), "EG_getBodyTopos(EDGE)");
  eg_check(EG_getBodyTopos(body, nullptr, NODE, &nnodes, &nodes), "EG_getBodyTopos(NODE)");


  // 6) Print results
  // We print only from rank 0 to avoid duplicate output under MPI.
  if (rank == 0) {
    std::printf("[LOAD] file=%s time_s=%.6f\n", file, tLoad);
    std::printf("[TOPO] bodies=%d faces=%d edges=%d nodes=%d\n",
                nchild, nfaces, nedges, nnodes);
  }


  // 7) Cleanup
  //
  // EGADS allocates memory for arrays like faces/edges/nodes/bodies.
  // We must free them with EG_free.
  EG_free(faces);
  EG_free(edges);
  EG_free(nodes);
  EG_free(bodies);

  // Close the context (this releases the model and all associated objects)
  EG_close(ctx);

  // Shut down MPI
  MPI_Finalize();
  return 0;
}
