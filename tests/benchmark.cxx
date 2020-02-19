#  if __has_include (<benchmark/benchmark.h>)
#    include <benchmark/benchmark.h>
#    define MSYS_WITH_BENCHMARK
#  endif

#include "io.hxx"
#include "clone.hxx"
#include "dms/dms.hxx"
#include "MsysThreeRoe.hpp"

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

static void BM_dms_jnk1_all(benchmark::State& state) {
    auto dms = Sqlite::read("tests/files/jnk1.dms");
    std::vector<std::string> tables;
    auto r = dms.fetch("sqlite_master");
    int col = r.column("name");
    for (; r; r.next()) {
        tables.push_back(r.get_str(col));
    }
    for (auto _ : state) {
        for (auto& t : tables) {
            auto r = dms.fetch(t);
            for (; r; r.next());
        }
    }
}

static void BM_dms_jnk1_structure(benchmark::State& state) {
    auto dms = Sqlite::read("tests/files/jnk1.dms");
    std::vector<std::string> tables {"particle", "bond", "cell", "provenance" };
    for (auto _ : state) {
        for (auto& t : tables) {
            auto r = dms.fetch(t);
            for (; r; r.next());
        }
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

static void BM_dms_water_name_text(benchmark::State& state) {
    auto dms = Sqlite::read("tests/files/water.db");
    for (auto _ : state) {
        auto r = dms.fetch("particle_text", true);
        int col = r.column("name");
        for (; r; r.next()) {
            const char* s = r.get_str(col);
            ThreeRoe(s, strlen(s)).Final();
        }
    }
}
static void BM_dms_water_name_ints(benchmark::State& state) {
    auto dms = Sqlite::read("tests/files/water.db");
    for (auto _ : state) {
        auto r = dms.fetch("particle_ints", true);
        int col = r.column("name");
        for (; r; r.next()) {
            r.get_int(col);
        }
    }
}



BENCHMARK(BM_SystemCreation);
BENCHMARK(BM_dms_jnk1_all)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_dms_jnk1_structure)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_dms_jnk1_particle)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_dms_jnk1_exclusion)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_dms_jnk1_stretch_term)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_Load_jnk1)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_Clone_jnk1)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_Load_jnk1_structure)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_Clone_jnk1_structure)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_dms_water_name_text)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_dms_water_name_ints)->Unit(benchmark::kMillisecond);

int main(int argc, char** argv) {
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
  return 0;
}

#endif
