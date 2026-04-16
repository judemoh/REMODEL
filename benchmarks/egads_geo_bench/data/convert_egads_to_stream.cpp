#include <cstdio>
#include <cstdlib>
#include "egads.h"

static void eg_check(int stat, const char* where) {
  if (stat != EGADS_SUCCESS) {
    std::fprintf(stderr, "EGADS error %d at %s\n", stat, where);
    std::exit(stat);
  }
}

int main(int argc, char** argv) {
  // Usage: ./convert_to_stream cube.egads cube.egads.stream
  if (argc < 3) {
    std::fprintf(stderr, "Usage: %s input.egads output.stream\n", argv[0]);
    return 1;
  }

  const char* inFile  = argv[1];
  const char* outFile = argv[2];

  ego ctx = nullptr, model = nullptr;

  eg_check(EG_open(&ctx), "EG_open");
  eg_check(EG_loadModel(ctx, 0, inFile, &model), "EG_loadModel");

  // Export model to stream bytes
  size_t nbyte = 0;
  char* stream = nullptr;
  eg_check(EG_exportModel(model, &nbyte, &stream), "EG_exportModel");

  std::FILE* fp = std::fopen(outFile, "wb");
  if (!fp) { std::perror("fopen"); return 2; }

  size_t nw = std::fwrite(stream, 1, nbyte, fp);
  std::fclose(fp);

  if (nw != nbyte) {
    std::fprintf(stderr, "Short write: %zu/%zu bytes\n", nw, nbyte);
    return 3;
  }

  EG_free(stream);
  EG_close(ctx);

  std::printf("[OK] Exported %zu bytes to %s\n", nbyte, outFile);
  return 0;
}
