#  if __has_include (<benchmark/benchmark.h>)
#    include <benchmark/benchmark.h>
#    define MSYS_WITH_BENCHMARK
#  endif

#include "io.hxx"
#include "clone.hxx"
#include "dms/dms.hxx"

using namespace desres::msys;

#if defined MSYS_WITH_BENCHMARK

static void BM_SystemCreation(benchmark::State& state) {
    for (auto _ : state) {
        System::create();
    }
}

static void BM_Load_jnk1(benchmark::State& state) {
    for (auto _ : state) {
        Load("tests/files/jnk1.dms");
    }
}

static void BM_Load_jnk1_structure(benchmark::State& state) {
    for (auto _ : state) {
        Load("tests/files/jnk1.dms", true);
    }
}

static void BM_Clone_jnk1(benchmark::State& state) {
    auto mol = Load("tests/files/jnk1.dms");
    for (auto _ : state) {
        Clone(mol, mol->atoms());
    }
}
static void BM_Clone_jnk1_structure(benchmark::State& state) {
    auto mol = Load("tests/files/jnk1.dms", true);
    for (auto _ : state) {
        Clone(mol, mol->atoms());
    }
}

static void BM_dms_jnk1_read(benchmark::State& state) {
    for (auto _ : state) {
        Sqlite::read("tests/files/jnk1.dms");
    }
}

static void BM_dms_jnk1_particle(benchmark::State& state) {
    auto dms = Sqlite::read("tests/files/jnk1.dms");
    for (auto _ : state) {
        auto r = dms.fetch("particle", false);
        for (; r; r.next());
    }
}

static void BM_dms_jnk1_exclusion(benchmark::State& state) {
    auto dms = Sqlite::read("tests/files/jnk1.dms");
    for (auto _ : state) {
        auto r = dms.fetch("exclusion", true);
        for (; r; r.next());
    }
}

static void BM_dms_jnk1_stretch_term(benchmark::State& state) {
    auto dms = Sqlite::read("tests/files/jnk1.dms");
    for (auto _ : state) {
        auto r = dms.fetch("stretch_harm_term", true);
        for (; r; r.next());
    }
}

BENCHMARK(BM_SystemCreation);
BENCHMARK(BM_dms_jnk1_read)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_dms_jnk1_particle)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_dms_jnk1_exclusion)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_dms_jnk1_stretch_term)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_Load_jnk1)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_Clone_jnk1)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_Load_jnk1_structure)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_Clone_jnk1_structure)->Unit(benchmark::kMillisecond);

int main(int argc, char** argv) {
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
  return 0;
}

#endif
