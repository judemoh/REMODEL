#include <cstdio>
#include <cstdlib>
#include <vector>
#include <random>
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

static double seconds_since(const std::chrono::high_resolution_clock::time_point& t0) {
  return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - t0).count();
}

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "Usage: %s cube.egads.stream [Nqueries]\n", argv[0]);
    return 1;
  }
  const char* streamFile = argv[1];
  const long long N = (argc >= 3) ? std::atoll(argv[2]) : 500000; // default half-million

  // ---- Import model ----
  std::vector<char> bytes = read_all_bytes(streamFile);
  const size_t nbyte = bytes.size();
  const char* streamPtr = bytes.data();

  ego ctx=nullptr, model=nullptr;
  eg_check(EG_open(&ctx), "EG_open");
  eg_check(EG_importModel(ctx, nbyte, streamPtr, &model), "EG_importModel");

  // model -> body -> faces
  ego geom=nullptr, *bodies=nullptr;
  int oclass=0, mtype=0, nchild=0, *senses=nullptr;
  double limits[6] = {0,0,0,0,0,0};
  eg_check(EG_getTopology(model, &geom, &oclass, &mtype, limits, &nchild, &bodies, &senses),
           "EG_getTopology(model)");
  ego body = bodies[0];

  int nfaces=0; ego* faces=nullptr;
  eg_check(EG_getBodyTopos(body, nullptr, FACE, &nfaces, &faces), "EG_getBodyTopos(FACE)");
  ego face0 = faces[0];

  // ---- Choose a base point on the face using mid-UV ----
  double frange[4]; int periodic=0;
  eg_check(EG_getRange(face0, frange, &periodic), "EG_getRange(face)");

  double uv_mid[2] = {
    frange[0] + 0.5*(frange[1]-frange[0]),
    frange[2] + 0.5*(frange[3]-frange[2])
  };

  // Evaluate surface at uv_mid -> xyz in eval[0:3]
  double eval[18];
  eg_check(EG_evaluate(face0, uv_mid, eval), "EG_evaluate(face)");
  const double p0[3] = {eval[0], eval[1], eval[2]};

  // ---- Generate query points near p0 ----
  std::mt19937 rng(123u);
  std::normal_distribution<double> normal(0.0, 1.0);
  const double eps = 1e-4;

  std::vector<double> Q(3*(size_t)N);
  for (long long i=0; i<N; ++i) {
    Q[3*i+0] = p0[0] + eps*normal(rng);
    Q[3*i+1] = p0[1] + eps*normal(rng);
    Q[3*i+2] = p0[2] + eps*normal(rng);
  }

  std::printf("[INVEVAL_SETUP] file=%s N=%lld faces=%d\n", streamFile, N, nfaces);

  // ---- Benchmark: same work, different loop batching ----
  std::vector<int> batches = {1, 8, 32, 128, 512, 2048};

for (int B : batches) {
  double outUV[2];
  double outXYZ[3];

  auto t0 = std::chrono::high_resolution_clock::now();

  for (long long i=0; i<N; i += B) {
    long long iEnd = (i + B < N) ? (i + B) : N;
    for (long long j=i; j<iEnd; ++j) {
      // EG_invEvaluate wants a non-const double* in this ESP build,
      // so make a tiny mutable copy.
      double xyz[3] = { Q[3*j+0], Q[3*j+1], Q[3*j+2] };
      (void) EG_invEvaluate(face0, xyz, outUV, outXYZ);
    }
  }

  double dt = seconds_since(t0);
  double qps = (dt > 0.0) ? (double)N / dt : 0.0;
  std::printf("[INVEVAL] batch=%d time_s=%.6f qps=%.3e\n", B, dt, qps);
}



  // IMPORTANT: skip EG_free(...) to avoid the double-free issue you saw.
  EG_close(ctx);
  return 0;
}
