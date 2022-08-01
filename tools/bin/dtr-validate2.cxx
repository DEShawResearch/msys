#include <sys/stat.h>
#include <signal.h>
#include <thread>
#include <msys/molfile/dtrplugin.hxx>
#include <msys/molfile/vmddir.h>
#include <msys/types.hxx>

namespace mf = desres::molfile;

static const char USAGE[] = R"USAGE(
usage: dtr-validate [-h] [--progress] [--checkpoint] [--energy]
                    [--parseable-output] [--parallelism PARALLELISM]
                    input_path [input_path [input_path ...]]

dtr-validate input.{dtr,atr,etr,stk} Check that each frame in each input file
is readable and provides position and box data.

positional arguments:
  input_path            dtr/atr/etr/stk to validate

optional arguments:
  -h, --help            show this help message and exit
  --progress            show progress bar
  --checkpoint          turn off delta_t checks
  --energy              turn off box, pos and atom based checks
  --parseable-output    show parseable output compatible with d_validate
  --parallelism PARALLELISM
                        Number of frames to process in parallel

)USAGE";

static bool isdir(std::string path) {
    struct stat stbuf[1];
    if (stat(path.data(), stbuf)) {
        MSYS_FAIL(strerror(errno));
    }
    return S_ISDIR(stbuf->st_mode);
}

static int64_t getsize(std::string path) {
    struct stat stbuf[1];
    if (stat(path.data(), stbuf)) {
        MSYS_FAIL(strerror(errno));
    }
    return stbuf->st_size;
}

static bool is_framefile(const char* basename) {
    bool ret = !strncmp(basename, "frame", 5);
    if (ret) {
        for (int i=0; i<9; i++) ret &= isdigit(basename[5+i]);
    }
    return ret;
}

static bool time_to_do(double now, double first,double interval) {
    int64_t F = round(first);
    int64_t I = round(interval);
    int64_t T = round(now);
    return ((T-F)>=0) && (I==0 || ((T-F) % I)==0);
}

int main(int argc, char *argv[]) {
    signal(SIGUSR1, SIG_IGN);

    if (argc < 2) {
        fprintf(stderr, "%s", USAGE);
        return 1;
    }
    int parallelism = 1;
    bool progress = false;
    bool energy = false;
    bool checkpoint = false;
    bool parseable_output = false;

    std::vector<std::string> paths;

    for (int i=1; i<argc; i++) {
        if (argv[i][0]!= '-') paths.push_back(argv[i]);
        if (!strcmp(argv[i], "--progress")) progress = true;
        if (!strcmp(argv[i], "--energy")) energy = true;
        if (!strcmp(argv[i], "--checkpoint")) checkpoint = true;
        if (!strcmp(argv[i], "--parseable-output")) parseable_output = true;
        if (!strcmp(argv[i], "--parallelism")) parallelism = atoi(argv[++i]);
    }
    if (parallelism < 1) {
        MSYS_FAIL("Invalid --parallelism argument parsed as " << parallelism);
    }

    for (auto& path : paths) {

        mf::DtrReader r(path, mf::DtrReader::RandomAccess);
        r.init();

        if (isdir(path)) {
            int64_t total_framebytes = 0;
            auto directory = vmd_opendir(path.data());
            while (auto entry=vmd_readdir(directory)) {
                if (is_framefile(entry)) {
                    total_framebytes += getsize(path + "/" + entry);
                }
            }
            vmd_closedir(directory);
            if (total_framebytes != r.total_bytes()) {
                MSYS_FAIL("total frame bytes " << total_framebytes <<
                              " vs expected from timekeys file " << r.total_bytes());
            }
        }

        int natoms = r.natoms();
        int nframes = r.size();
        double first_frame_time = 999;
        double last_frame_time = -99;
        double last_delta_t = -99;
        if (progress) {
            printf("path: %s frames: %-9d atoms: %-9d\n", path.data(), nframes, natoms);
        }

        std::vector<molfile_timestep_t> frames(parallelism);
        for (auto& f : frames) {
            memset(&f, 0, sizeof(f));
            if (!energy) {
                f.coords = new float[3*natoms];
            }
        }

        for (int start=0; start<nframes; start += parallelism) {
            std::vector<std::thread> threads;
            for (int j=0; j<parallelism && start+j < nframes; j++) {
                auto& f = frames[j];
                // write invalid data into box so that we can tell if something was written
                memset(f.unit_cell, 7, sizeof(f.unit_cell));
                threads.emplace_back(&mf::DtrReader::frame, &r, start+j, &f, nullptr);
            }
            for (auto& t : threads) {
                t.join();
            }
            for (int j=0; j<parallelism && start+j < nframes; j++) {
                auto i = start + j;
                if (progress) {
                    if (time_to_do(i, 0, nframes/50)) fprintf(stderr, ".");
                    if (time_to_do(i, 0, (nframes)/5)) fprintf(stderr, "%d%%", 100*(i+1)/nframes);
                }
                auto& f = frames[j];
                auto frame_time = f.physical_time;
                if (!energy) {
                    // The low level frame parsing code already checks that pos is read;
                    // however, we have to check for box ourselves for some reason.
                    if (std::vector<double>(f.unit_cell, f.unit_cell+9) == std::vector<double>(9, 7)) {
                        MSYS_FAIL(path << " is missing box data at frame " << i);
                    }
                }

                if (i==0) {
                    first_frame_time = frame_time;
                } else {
                    auto delta_t = frame_time - last_frame_time;
                    if (i>1 && !checkpoint) {
                        auto change_delta_t = std::abs(last_delta_t - delta_t);
                        if (change_delta_t >= 0.02) {
                            MSYS_FAIL("delta_t changed from " << last_delta_t << " to " << delta_t);
                        }
                    }
                    last_delta_t = delta_t;
                }
                last_frame_time = frame_time;
            }
        }
        if (progress) fprintf(stderr, "\n");

        for (auto& f : frames) {
            delete[] f.coords;
        }

        if (parseable_output) {
            printf("first_frame_time = %f\n", first_frame_time);
            printf("last_frame_time = %f\n", last_frame_time);
            if (nframes > 1) {
                printf("delta_t = %f\n", (last_frame_time - first_frame_time) / (nframes - 1));
                printf("num_frames = %d\n", nframes);
            }
        }
    }
    return 0;
}


