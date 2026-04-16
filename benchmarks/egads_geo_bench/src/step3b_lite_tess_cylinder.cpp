#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include "egads.h"

static void eg_check(int stat, const char* where) {
  if (stat != EGADS_SUCCESS) {
    std::fprintf(stderr, "EGADS error %d at %s\n", stat, where);
    std::exit(stat);
  }
}

static std::vector<char> read_all_bytes(const char* path) {
  std::FILE* fp = std::fopen(path, "rb");
  if (!fp) { std::perror("fopen"); std::exit(2); }
  std::fseek(fp, 0, SEEK_END);
  long n = std::ftell(fp);
  std::fseek(fp, 0, SEEK_SET);
  if (n <= 0) { std::fprintf(stderr, "Empty file: %s\n", path); std::exit(3); }
  std::vector<char> buf((size_t)n);
  size_t nr = std::fread(buf.data(), 1, buf.size(), fp);
  std::fclose(fp);
  if (nr != buf.size()) { std::fprintf(stderr, "Short read on %s\n", path); std::exit(4); }
  return buf;
}

// Choose the body that has the most faces (avoids picking wire bodies)
static ego pick_body_with_most_faces(ego* bodies, int nbodies) {
  ego best = bodies[0];
  int bestFaces = -1;

  for (int i = 0; i < nbodies; ++i) {
    int nfaces = 0;
    ego* faces = nullptr;
    // do not EG_free(faces) to avoid double-free issue in this build
    (void) EG_getBodyTopos(bodies[i], nullptr, FACE, &nfaces, &faces);
    if (nfaces > bestFaces) {
      bestFaces = nfaces;
      best = bodies[i];
    }
  }
  return best;
}

int main(int argc, char** argv) {
  const char* streamFile = (argc >= 2) ? argv[1] : "../data/cylinder.egads.stream";

  // Import model
  std::vector<char> bytes = read_all_bytes(streamFile);
  const size_t nbyte = bytes.size();
  const char* streamPtr = bytes.data();

  ego ctx=nullptr, model=nullptr;
  eg_check(EG_open(&ctx), "EG_open");

  auto t0 = std::chrono::high_resolution_clock::now();
  eg_check(EG_importModel(ctx, nbyte, streamPtr, &model), "EG_importModel");
  auto t1 = std::chrono::high_resolution_clock::now();
  double tImport = std::chrono::duration<double>(t1 - t0).count();

  // Model -> bodies
  ego geom=nullptr, *bodies=nullptr;
  int oclass=0, mtype=0, nbodies=0, *senses=nullptr;
  double limits[6] = {0,0,0,0,0,0};
  eg_check(EG_getTopology(model, &geom, &oclass, &mtype, limits, &nbodies, &bodies, &senses),
           "EG_getTopology(model)");

  // Pick the solid body (most faces)
  ego body = pick_body_with_most_faces(bodies, nbodies);

  int nfaces=0; ego* faces=nullptr;
  eg_check(EG_getBodyTopos(body, nullptr, FACE, &nfaces, &faces), "EG_getBodyTopos(FACE)");

  std::printf("[IMPORT] file=%s bytes=%zu time_s=%.6f\n", streamFile, nbyte, tImport);
  std::printf("[TOPO] bodies=%d chosen_body_faces=%d\n", nbodies, nfaces);

  if (nfaces < 1) {
    std::fprintf(stderr, "No faces found on chosen body — cannot tessellate.\n");
    EG_close(ctx);
    return 2;
  }

  // Tessellation sweep
  std::vector<double> defl = {1e-2, 5e-3, 1e-3, 5e-4};

  for (double h : defl) {
    double params[3] = {0.0, h, 20.0};
    ego tess = nullptr;

    auto tt0 = std::chrono::high_resolution_clock::now();
    eg_check(EG_makeTessBody(body, params, &tess), "EG_makeTessBody");
    auto tt1 = std::chrono::high_resolution_clock::now();
    double tTess = std::chrono::duration<double>(tt1 - tt0).count();

    long long tri_total = 0;
    for (int fIndex = 1; fIndex <= nfaces; ++fIndex) {
      int npnt=0, ntri=0;
      const double *pxyz=nullptr, *puv=nullptr;
      const int *ptype=nullptr, *pindx=nullptr, *tris=nullptr, *tric=nullptr;

      eg_check(EG_getTessFace(tess, fIndex, &npnt, &pxyz, &puv, &ptype, &pindx,
                              &ntri, &tris, &tric),
               "EG_getTessFace");
      tri_total += ntri;
    }

    std::printf("[TESS] deflection=%g time_s=%.6f triangles=%lld\n", h, tTess, tri_total);
    EG_deleteObject(tess);
  }

  EG_close(ctx);
  return 0;
}
