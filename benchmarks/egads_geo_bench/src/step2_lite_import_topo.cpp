#include <cstdio>
#include <cstdlib>
#include <vector>
#include <chrono>
#include "egads.h"

// Helper - tells me if goemetry call fails
static void eg_check(int stat, const char* where) {
  // Normal completion without errors.
  if (stat != EGADS_SUCCESS) {
    std::fprintf(stderr, "EGADS error %d at %s\n", stat, where);
    std::exit(stat);
  }
}
// reading stream file into ememory
static std::vector<char> read_all_bytes(const char* path) {
  std::FILE* fp = std::fopen(path, "rb");
  if (!fp) { std::perror("fopen"); std::exit(2); }

  std::fseek(fp, 0, SEEK_END);
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

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "Usage: %s cube.egads.stream\n", argv[0]);
    return 1;
  }
  const char* streamFile = argv[1];

  // 1) Read stream bytes
  std::vector<char> bytes = read_all_bytes(streamFile);
  const size_t nbyte = bytes.size();
  const char* streamPtr = bytes.data();

  // 2) Import into EGADSlite
  ego ctx = nullptr, model = nullptr;
  eg_check(EG_open(&ctx), "EG_open");

  auto t0 = std::chrono::high_resolution_clock::now();
  eg_check(EG_importModel(ctx, nbyte, streamPtr, &model), "EG_importModel");
  auto t1 = std::chrono::high_resolution_clock::now();
  double tImport = std::chrono::duration<double>(t1 - t0).count();

  // 3) Topology query
  ego geom = nullptr;
  ego* bodies = nullptr;
  int oclass = 0, mtype = 0, nchild = 0;
  int* senses = nullptr;
  double limits[6] = {0,0,0,0,0,0};

  eg_check(EG_getTopology(model, &geom, &oclass, &mtype,
                          limits, &nchild, &bodies, &senses),
           "EG_getTopology(model)");

  ego body = bodies[0];

  int nfaces = 0, nedges = 0, nnodes = 0;
  ego *faces = nullptr, *edges = nullptr, *nodes = nullptr;
  eg_check(EG_getBodyTopos(body, nullptr, FACE, &nfaces, &faces), "EG_getBodyTopos(FACE)");
  eg_check(EG_getBodyTopos(body, nullptr, EDGE, &nedges, &edges), "EG_getBodyTopos(EDGE)");
  eg_check(EG_getBodyTopos(body, nullptr, NODE, &nnodes, &nodes), "EG_getBodyTopos(NODE)");

  std::printf("[IMPORT] file=%s bytes=%zu time_s=%.6f\n", streamFile, nbyte, tImport);
  std::printf("[TOPO] bodies=%d faces=%d edges=%d nodes=%d\n", nchild, nfaces, nedges, nnodes);
  // NOTE: This ESP/EGADSlite build is double-freeing something during cleanup.
  // For short benchmark runs, skip freeing returned arrays and just close.
  // Cleanup
  EG_close(ctx);
  return 0;
}
