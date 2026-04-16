#include <cstdio>
#include <cstdlib>
#include <vector>
#include <random>
#include <chrono>
#include "/Users/jude/Documents/REMODEL/EGADSlite/EGADS/include/egads.h"

/*
Question: How expensive is geometry interrogation when simulations repeatedly interact with CAD geometry?
The following questions are answered 1E6 times:
1) “What surface is closest to this point?”
2) “Where on the CAD surface does this mesh node lie?”
3) “What are the local surface coordinates here?”

Main questions to answer include:
1) Which operations are most expensive? (e.g. closest point search, inverse evaluation, etc.)
2) How does the cost depend on the geometry? 

Query that happens often: 
Inverse evaluation: Given a 3D point (x,y,z), find the corresponding parametric coordinates (u,v) on the surface.

We meausre:
1) cost of geometry interrogration (eueries per second)
2) Batching

A cylinder introduces:
1) Non-zero curvature 
2) Non-linear inverse mapping & realistic CAD behavior



What I did: 
1)Load a b-rep via Egadslite stream
2) Run representative geometry interrogration - Inverse evaluation: Given a 3D point (x,y,z), find the corresponding parametric coordinates (u,v) on the surface.
3) Measuring throughput (queries per second) and how it depends on batch size (number of queries processed together) & sensitivity to geometric complexity
4) Evidence to motivate geometry handling& locality-aware access

*/ 

// void produces a result that cannot be used elsewhere
static void eg_check(int stat, const char* where) {
  if (stat != EGADS_SUCCESS) {
    std::fprintf(stderr, "EGADS error %d at %s\n", stat, where);
    std::exit(stat);
  }
}

static std::vector<char> read_all_bytes(const char* path) {
  std::FILE* fp = std::fopen(path, "rb");
  if (!fp) { std::perror("fopen"); std::exit(2); }
// seek to end to get file size, then back to start for reading
  std::fseek(fp, 0, SEEK_END);
// compute file size
  long n = std::ftell(fp);
  std::fseek(fp, 0, SEEK_SET);
  if (n <= 0) { std::fprintf(stderr, "Empty file: %s\n", path); std::exit(3); }
  //  std::vector<char> reads entire file into a bite array
  std::vector<char> buf((size_t)n);
  size_t nr = std::fread(buf.data(), 1, buf.size(), fp);
  std::fclose(fp);
  if (nr != buf.size()) { std::fprintf(stderr, "Short read on %s\n", path); std::exit(4); }
  return buf;
}

// converts a time point to seconds since that time point
static double seconds_since(const std::chrono::high_resolution_clock::time_point& t0) {
  return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count();
}

// Pick the body that has the most faces (avoids picking wire bodies)
static ego pick_body_with_most_faces(ego* bodies, int nbodies) {
  ego best = bodies[0];
  int bestFaces = -1;

  for (int i = 0; i < nbodies; ++i) {
    // stores count of mumber of faces 
    int nfaces = 0; 
    // populate pointer array 
    ego* faces = nullptr;

    (void) EG_getBodyTopos(bodies[i], nullptr, FACE, &nfaces, &faces);
    if (nfaces > bestFaces) {
      bestFaces = nfaces;
      best = bodies[i];
    }
  }
  return best;
}

int main(int argc, char** argv) {
  //retrieve the stream path, tenary operator checks if user provided an argument, otherwise uses default path
  const char* streamFile = (argc >= 2) ? argv[1] : "../data/cylinder.egads.stream";
//  N =  number of query points to generate, default 500k, can be set as second argument to executable
  const long long N = (argc >= 3) ? std::atoll(argv[2]) : 500000;

  // Import model, evaluate model size
  std::vector<char> bytes = read_all_bytes(streamFile);
  const size_t nbyte = bytes.size();
  // pointer to byte data, used for EG_importModel (needs to know buffer length to parse stream correctly)
  const char* streamPtr = bytes.data();

  //FYI egi creates blind pointers and initializes both to nullptr - pass them to EGADS instead of working with underlying data structure
  ego ctx=nullptr, model=nullptr;
  // loads cad model from byte stream prepared earlier
  eg_check(EG_open(&ctx), "EG_open");
  eg_check(EG_importModel(ctx, nbyte, streamPtr, &model), "EG_importModel");

  // Model -> bodies
  ego geom=nullptr, *bodies=nullptr;
  //oclass = type of geometry (e.g. solid, surface, curve), mytype - definition depends on toclass,  senses = orientation of each body
  int oclass=0, mtype=0, nbodies=0, *senses=nullptr;
  double limits[6] = {0,0,0,0,0,0};
  eg_check(EG_getTopology(model, &geom, &oclass, &mtype, limits, &nbodies, &bodies, &senses),
           "EG_getTopology(model)");

  // Pick solid body
  ego body = pick_body_with_most_faces(bodies, nbodies);

  // Get faces
  int nfaces=0; ego* faces=nullptr;
  eg_check(EG_getBodyTopos(body, nullptr, FACE, &nfaces, &faces), "EG_getBodyTopos(FACE)");
  if (nfaces < 1) {
    std::fprintf(stderr, "No faces found on chosen body (nfaces=%d)\n", nfaces);
    EG_close(ctx);
    return 2;
  }

  // Use first face (fine for baseline; later we can choose curved face specifically)
  ego face0 = faces[0];

  // Build a base point p0 from mid-UV on this face
  double frange[4]; int periodic=0;
  eg_check(EG_getRange(face0, frange, &periodic), "EG_getRange(face)");

  double uv_mid[2] = {
    // 
    frange[0] + 0.5*(frange[1]-frange[0]),
    frange[2] + 0.5*(frange[3]-frange[2])
  };

  double eval[18];
  //eg_evaluate returns coordinates and derivatives on a curve, always returns 18 values 
  // 18-element eval array with the 3D position (indices 0–2), first partial derivatives with respect to U and V (indices 3–8), and second partial derivatives (indices 9–17). 
  eg_check(EG_evaluate(face0, uv_mid, eval), "EG_evaluate(face)");
  const double p0[3] = {eval[0], eval[1], eval[2]};

  // Generate query points near p0, normally distributed perturbations around reference point p0, with standard deviation eps.
  // These will be the input to EG_invEvaluate, which should return UV parameters near uv_mid.
  std::mt19937 rng(123u);
  std::normal_distribution<double> normal(0.0, 1.0);
  //eps is a scaling factor that controls the spread of the query points around p0. 
  //smaller eps means points will be closer to p0, while a larger eps means points will be more spread out. 
  // eps can affect the difficulty of the inverse evaluation problem, as points farther from p0 --> more iterations for convergence 
  // may  fail to converge if points too far from the surface.
  const double eps = 1e-4;
// tight cluster of query points around p0, with some random variation to avoid identical queries
// tight clustering stress-tests the convergence of the inverse evaluation, as many points will be close together and may require more iterations to resolve accurately


//allocate storage for N 3D query points (N= no queries performed)
//
  std::vector<double> Q(3*(size_t)N);
  for (long long i=0; i<N; ++i) {
    Q[3*i+0] = p0[0] + eps*normal(rng);
    Q[3*i+1] = p0[1] + eps*normal(rng);
    Q[3*i+2] = p0[2] + eps*normal(rng);
  }

  std::printf("[INVEVAL_SETUP] file=%s bytes=%zu bodies=%d chosen_faces=%d N=%lld\n",
              streamFile, nbyte, nbodies, nfaces, N);

// chunking a list of queries into batches of size B (variable) --> Effect of B on performance 

  std::vector<int> batches = {1, 8, 32, 128, 512, 2048};

  for (int B : batches) {
    double outUV[2];
    double outXYZ[3];

    auto t0 = std::chrono::high_resolution_clock::now();

    for (long long i=0; i<N; i += B) {
      long long iEnd = (i + B < N) ? (i + B) : N;
      for (long long j=i; j<iEnd; ++j) {
        // EG_invEvaluate signature here requires non-const double*
        double xyz[3] = { Q[3*j+0], Q[3*j+1], Q[3*j+2] };
        (void) EG_invEvaluate(face0, xyz, outUV, outXYZ);
      }
    }

    double dt = seconds_since(t0);
    double qps = (dt > 0.0) ? (double)N / dt : 0.0;
    std::printf("[INVEVAL] batch=%d time_s=%.6f qps=%.3e\n", B, dt, qps);
  }

  EG_close(ctx);
  return 0;
}
