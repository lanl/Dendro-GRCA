/**
 * @file bssn_scaling_bench.cpp
 * @brief Strong/weak scaling benchmark for the BSSN hybrid path. Builds a fixed
 * grid (bbh puncture or uniform), runs N RK steps with the full per-step compute
 * but NO remesh and NO IO, and emits the solver's per-step profiling JSONL.
 * Built through the same flag system as bssnSolver (mirrors its config).
 *
 * Usage: bssnScalingBench <paramFile> [--mode strong|weak] [--grid bbh|uniform]
 *        [--steps N] [--warmup K] [--lev L] [--id-type N] [--prefix P]
 */

#include <cmath>
#include <cstdlib>
#include <ctime>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include "TreeNode.h"
#include "bssnCtx.h"
#include "gr.h"
#include "grUtils.h"
#include "mathUtils.h"
#include "mesh.h"
#include "meshUtils.h"
#include "mpi.h"
#include "octUtils.h"
#include "parameters.h"

namespace {

struct BenchOpts {
    std::string param_file;
    std::string mode   = "strong";  // strong | weak
    std::string grid   = "bbh";     // bbh | uniform
    unsigned int steps = 10;
    unsigned int warmup = 2;
    unsigned int lev   = 3;  // uniform-grid refinement level (8^lev elements)
    // init data; -1 = use the param file's BSSN_ID_TYPE (puncture). Note: a
    // uniform grid + puncture is only marginally stable -- prefer --grid bbh.
    int id_type        = -1;
    // weak-scaling per-rank grain (BSSN_DENDRO_GRAIN_SZ); -1 = use param file.
    int grain          = -1;
    std::string prefix = "bssn_bench";
    // print mesh/state fingerprints; digests must not vary with threads
    bool fingerprint   = false;
};

void usage(const char* a0) {
    std::cout
        << "Usage: " << a0 << " <paramFile> [options]\n"
        << "  --mode strong|weak   strong = fixed grid (vary ranks); weak = fix "
           "work/rank (default strong)\n"
        << "  --grid bbh|uniform   bbh = representative puncture grid; uniform = "
           "clean kernel scaling (default bbh)\n"
        << "  --steps N            timed RK steps (default 10)\n"
        << "  --warmup K           untimed warm-up RK steps (default 2)\n"
        << "  --lev L              uniform-grid depth, 8^L elements (default 3)\n"
        << "  --id-type N          override BSSN_ID_TYPE init data (default: "
           "use param file; 0/1 = puncture)\n"
        << "  --grain N            per-rank grain (BSSN_DENDRO_GRAIN_SZ); weak "
           "mode keeps ~N elems/rank\n"
        << "  --prefix P           output prefix for <P>_steps.jsonl (default "
           "bssn_bench)\n"
        << "  --fingerprint        print mesh + state digests and exit after "
           "steps; digests MUST be identical across OMP_NUM_THREADS\n";
}

BenchOpts parse(int argc, char** argv) {
    BenchOpts o;
    if (argc >= 2) o.param_file = argv[1];
    for (int i = 2; i < argc; i++) {
        const std::string a = argv[i];
        auto need           = [&](const char* name) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << name << " requires a value\n";
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
            return argv[++i];
        };
        if (a == "--mode")
            o.mode = need("--mode");
        else if (a == "--grid")
            o.grid = need("--grid");
        else if (a == "--steps")
            o.steps = (unsigned int)std::stoul(need("--steps"));
        else if (a == "--warmup")
            o.warmup = (unsigned int)std::stoul(need("--warmup"));
        else if (a == "--lev")
            o.lev = (unsigned int)std::stoul(need("--lev"));
        else if (a == "--id-type")
            o.id_type = std::stoi(need("--id-type"));
        else if (a == "--grain")
            o.grain = std::stoi(need("--grain"));
        else if (a == "--prefix")
            o.prefix = need("--prefix");
        else if (a == "--fingerprint")
            o.fingerprint = true;
        else
            std::cerr << "[bench] ignoring unknown arg: " << a << "\n";
    }
    return o;
}

}  // namespace

int main(int argc, char** argv) {
    // Hybrid path runs OpenMP inside compute regions, not around MPI calls.
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    const BenchOpts opt = parse(argc, argv);
    if (opt.param_file.empty()) {
        if (!rank) usage(argv[0]);
        MPI_Finalize();
        return 0;
    }

    // ---- parameters (same loader as the solver) ------------------------
    bssn::readParamFile(opt.param_file.c_str(), comm);
    if (opt.id_type >= 0) bssn::BSSN_ID_TYPE = (unsigned int)opt.id_type;
    if (opt.grain > 0) bssn::BSSN_DENDRO_GRAIN_SZ = (unsigned int)opt.grain;
    // Resolve the RHS schedule (and set the block-major first-touch flag for the
    // "balanced" mode) BEFORE the ctx allocates the unzip buffers below, so the
    // NUMA first-touch placement is established on the first allocation. Mirrors
    // the solver's post-readParamFile call. No-op unless DENDRO_HYBRID_OMP.
    bssn::set_rhs_omp_schedule();
    _InitializeHcurve(bssn::BSSN_DIM);
    m_uiMaxDepth                   = bssn::BSSN_MAXDEPTH;
    bssn::BSSN_PROFILE_FILE_PREFIX = opt.prefix;

    if (!rank)
        std::cout << "[bench] mode=" << opt.mode << " grid=" << opt.grid
                  << " steps=" << opt.steps << " warmup=" << opt.warmup
                  << " lev=" << opt.lev << " npes=" << npes
                  << " rk_type=" << bssn::BSSN_RK_TYPE << " (no remesh, no IO)"
                  << std::endl;

    // ---- build the grid ------------------------------------------------
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double, double, double, double*)> f_init =
        [](double x, double y, double z, double* var) {
            bssn::punctureData(x, y, z, var);
        };
    const unsigned int interpVars = bssn::BSSN_NUM_VARS;
    unsigned int varIndex[interpVars];
    for (unsigned int i = 0; i < interpVars; i++) varIndex[i] = i;

    if (opt.grid == "uniform") {
        // Uniform octree (8^lev elements); weak scaling grows lev with ranks.
        unsigned int lev_use = opt.lev;
        if (opt.mode == "weak")
            lev_use += (unsigned int)std::lround(std::log((double)npes) /
                                                 std::log(8.0));
        createRegularOctree(tmpNodes, lev_use, bssn::BSSN_DIM, m_uiMaxDepth,
                            comm);
        if (!rank)
            std::cout << "[bench] uniform grid: lev=" << lev_use << std::endl;
    } else {
        // Representative BBH grid: refine on puncture data like the solver.
        const unsigned int f2olmin =
            std::min(bssn::BSSN_BH1_MAX_LEV, bssn::BSSN_BH2_MAX_LEV);
        if (f2olmin < MAXDEAPTH_LEVEL_DIFF + 2) {
            if (!rank)
                std::cout << "BH min level must be > " << (MAXDEAPTH_LEVEL_DIFF + 2)
                          << std::endl;
            MPI_Abort(comm, 0);
        }
        function2Octree(f_init, bssn::BSSN_NUM_VARS, varIndex, interpVars,
                        tmpNodes, (f2olmin - MAXDEAPTH_LEVEL_DIFF - 2),
                        bssn::BSSN_WAVELET_TOL, bssn::BSSN_ELE_ORDER, comm);
    }

    ot::Mesh* mesh = ot::createMesh(
        tmpNodes.data(), tmpNodes.size(), bssn::BSSN_ELE_ORDER, comm, 1,
        ot::SM_TYPE::FDM, bssn::BSSN_DENDRO_GRAIN_SZ, bssn::BSSN_LOAD_IMB_TOL,
        bssn::BSSN_SPLIT_FIX);
    mesh->setDomainBounds(
        Point(bssn::BSSN_GRID_MIN_X, bssn::BSSN_GRID_MIN_Y, bssn::BSSN_GRID_MIN_Z),
        Point(bssn::BSSN_GRID_MAX_X, bssn::BSSN_GRID_MAX_Y, bssn::BSSN_GRID_MAX_Z));
    unsigned int lmin, lmax;
    mesh->computeMinMaxLevel(lmin, lmax);
    tmpNodes.clear();

    // Fixed grid: disables AMR + the ctx's collective init grid-converge (which
    // hangs at multi-rank). Mirrors gr_scaling.cpp.
    bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY = 1;

    // dx / dt (same expressions as the solver)
    bssn::BSSN_CURRENT_MIN_DX =
        ((bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
         ((1u << (m_uiMaxDepth - lmax)) / ((double)bssn::BSSN_ELE_ORDER)) /
         ((double)(1u << (m_uiMaxDepth))));
    bssn::BSSN_RK45_TIME_STEP_SIZE =
        bssn::BSSN_CFL_FACTOR * bssn::BSSN_CURRENT_MIN_DX;

    // weak scaling on the BBH grid: split so each rank keeps ~grain-size elems.
    ot::Mesh* pMesh = mesh;
    if (opt.mode == "weak" && opt.grid == "bbh")
        pMesh = bssn::weakScalingReMesh(mesh, npes);

    // ---- ctx + stepper (mirror the solver, incl. MSRK/TSRK) ------------
    bssn::BSSNCtx* bssnCtx = new bssn::BSSNCtx(pMesh);
    const RKType rkType    = (RKType)bssn::BSSN_RK_TYPE;
    ts::ETS<DendroScalar, bssn::BSSNCtx>* ets = nullptr;
    if (rkType == RKType::RK4_MSRK2_1 || rkType == RKType::RK4_MSRK2_2 ||
        rkType == RKType::RK4_MSRK3) {
        ts::ETSType v = (rkType == RKType::RK4_MSRK2_1) ? ts::ETSType::RK4_MSRK2_1
                        : (rkType == RKType::RK4_MSRK2_2)
                            ? ts::ETSType::RK4_MSRK2_2
                            : ts::ETSType::RK4_MSRK3;
        ets = new ts::ETS_MSRK<DendroScalar, bssn::BSSNCtx>(bssnCtx, v);
    } else if (rkType == RKType::RK6_TSRK) {
        ets = new ts::ETS_TSRK<DendroScalar, bssn::BSSNCtx>(bssnCtx);
    } else {
        // RKType is index-matched to ts::ETSType: cast covers every
        // single-step method, falling back to RK4 on an unknown value.
        ets = new ts::ETS<DendroScalar, bssn::BSSNCtx>(bssnCtx);
        if (ets->set_ets_coefficients((ts::ETSType)rkType) != 0)
            ets->set_ets_coefficients(ts::ETSType::RK4);
    }
    ets->set_evolve_vars(bssnCtx->get_evolution_vars());
    ets->init();

    {
        // Report element + block counts (per-rank grain is the MPI split knob).
        const ot::Mesh* cmesh = ets->get_mesh();
        const bool act        = cmesh->isActive();
        DendroIntL le = act ? cmesh->getNumLocalMeshElements() : 0;
        DendroIntL lb = act ? (DendroIntL)cmesh->getLocalBlockList().size() : 0;
        DendroIntL ge = 0, gb = 0;
        par::Mpi_Reduce(&le, &ge, 1, MPI_SUM, 0, comm);
        par::Mpi_Reduce(&lb, &gb, 1, MPI_SUM, 0, comm);
        // Mesh digests BEFORE any evolution: a mesh-construction change shows up
        // here, isolated from anything the RHS/unzip threading might do.
        if (opt.fingerprint) bssn::meshFingerprint(cmesh, "mesh");
        if (!rank)
            std::cout << "[bench] elements=" << ge << " (" << ge / npes
                      << "/rank)  blocks=" << gb << " (" << gb / npes
                      << "/rank)  lmin=" << lmin << " lmax=" << lmax << std::endl;
    }

    // ---- warm-up + timed loop: NO remesh, NO IO ------------------------
    const DendroIntL start_step = ets->curr_step();
    const DendroIntL last_step  = start_step + opt.warmup + opt.steps;
    while (ets->curr_step() < last_step) {
        const DendroIntL step = ets->curr_step();
        const DendroIntL rel  = step - start_step;

        bssn::BSSN_CURRENT_RK_COORD_TIME = ets->curr_time();
        bssn::BSSN_CURRENT_RK_STEP       = step;

        // Zero every timer right before the first TIMED step (discard warm-up).
        if (rel == opt.warmup) {
#if defined(__PROFILE_ETS__) && defined(__PROFILE_CTX__)
            ets->reset_pt();
            bssnCtx->reset_pt();
#endif
            bssn::timer::resetSnapshot();
        }

        ets->evolve();

        // One JSONL record per timed step; reset after so each snap is one step.
        if (rel >= opt.warmup) {
#ifdef ENABLE_DENDRO_PROFILE_COUNTERS
#if defined(__PROFILE_ETS__) && defined(__PROFILE_CTX__)
            bssn::timer::profileInfoJSON(bssn::BSSN_PROFILE_FILE_PREFIX.c_str(),
                                         ets->get_mesh(), (unsigned int)step,
                                         &ets->m_uiCtxpt, &bssnCtx->m_uiCtxpt);
            ets->reset_pt();
#else
            bssn::timer::profileInfoJSON(bssn::BSSN_PROFILE_FILE_PREFIX.c_str(),
                                         ets->get_mesh(), (unsigned int)step);
#endif
            bssn::timer::resetSnapshot();
#endif
        }

        if (!rank)
            std::cout << "[bench] step " << step
                      << (rel < opt.warmup ? " (warmup)" : " (timed)")
                      << " t=" << ets->curr_time() << std::endl;
    }

    if (opt.fingerprint) {
        // State digest AFTER the timed steps: proves the answers are identical,
        // not just the mesh. Bit-exact by construction across threads, so any
        // difference is a real race or a reordered reduction.
        DendroScalar* eVar[bssn::BSSN_NUM_VARS];
        bssnCtx->get_evolution_vars().to_2d(eVar);
        bssn::stateFingerprint(ets->get_mesh(), eVar, bssn::BSSN_NUM_VARS,
                               "evolved");
    }

    if (!rank)
        std::cout << "[bench] done: " << opt.steps << " timed steps -> "
                  << opt.prefix << "_steps.jsonl" << std::endl;

    // ---- cleanup -------------------------------------------------------
    ot::Mesh* final_mesh = bssnCtx->get_mesh();
    delete bssnCtx;
    delete final_mesh;
    delete ets;

    MPI_Finalize();
    return 0;
}
