// rhs_experimental.cpp -- experimental BSSN RHS: prod path + all dispatch
// (cascade/naive/AVX*) behind compile flags, for anything being actively
// tested. Built instead of src/rhs.cpp when any such flag is set (CMake
// BSSN_RHS_SRC). Keep src/rhs.cpp pristine prod.

#include "rhs.h"

#include "gr.h"
#include "hadrhs.h"
#include "parameters.h"
#ifdef DENDRO_USE_NEW_DERIVS
#include "derivatives.h"  // full DendroDerivatives type for filter_cako()
#endif

#if defined(BSSN_USE_CASCADE_AVX) || defined(BSSN_USE_CASCADE_AVX_FUSED) || \
    defined(BSSN_USE_CASCADE_AVX512) ||                                      \
    defined(BSSN_USE_CASCADE_AVX512_FUSED) ||                               \
    defined(BSSN_USE_CASCADE_IR_AVX512) || defined(BSSN_USE_CASCADE_IR_AVX) || \
    defined(BSSN_USE_CASCADE_IR_AVX512_FUSED) ||                            \
    defined(BSSN_USE_CASCADE_IR_AVX_FUSED)
// Generated cascade bodies define their own VEC typedef (scoped to their
// body {}) and #undef/#define the macros before use, so AVX2 and AVX-512
// variants can coexist in the same TU.
#include <immintrin.h>

#include <cmath>
#endif

using namespace std;
using namespace bssn;

void bssnRHS(double **uzipVarsRHS, const double **uZipVars,
             const ot::Block *blkList, unsigned int numBlocks,
             const double curr_time, const double **uZipConstVars) {
    // First-call trace: prints once per process at first entry so the
    // colleague can confirm the cascade RHS is actually reached at scale.
    {
        static bool s_first_rhs_call = true;
        if (s_first_rhs_call) {
            int rank_global = 0;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank_global);
            if (!rank_global)
                std::cout << "[trace] bssnRHS first call (t=" << curr_time
                          << ", numBlocks=" << numBlocks << ")" << std::endl;
            s_first_rhs_call = false;
        }
    }

    unsigned int offset;
    double ptmin[3], ptmax[3];
    unsigned int sz[3];
    unsigned int bflag;
    double dx, dy, dz;
    const Point pt_min(bssn::BSSN_COMPD_MIN[0], bssn::BSSN_COMPD_MIN[1],
                       bssn::BSSN_COMPD_MIN[2]);
    const Point pt_max(bssn::BSSN_COMPD_MAX[0], bssn::BSSN_COMPD_MAX[1],
                       bssn::BSSN_COMPD_MAX[2]);
    const unsigned int PW = bssn::BSSN_PADDING_WIDTH;

#ifdef BSSN_ENABLE_CUDA
    cuda::BSSNComputeParams bssnParams;
    bssnParams.BSSN_LAMBDA[0]    = bssn::BSSN_LAMBDA[0];
    bssnParams.BSSN_LAMBDA[1]    = bssn::BSSN_LAMBDA[1];
    bssnParams.BSSN_LAMBDA[2]    = bssn::BSSN_LAMBDA[2];
    bssnParams.BSSN_LAMBDA[3]    = bssn::BSSN_LAMBDA[3];

    bssnParams.BSSN_LAMBDA_F[0]  = bssn::BSSN_LAMBDA_F[0];
    bssnParams.BSSN_LAMBDA_F[1]  = bssn::BSSN_LAMBDA_F[1];

    bssnParams.BSSN_ETA_POWER[0] = bssn::BSSN_ETA_POWER[0];
    bssnParams.BSSN_ETA_POWER[1] = bssn::BSSN_ETA_POWER[1];

    bssnParams.ETA_R0            = bssn::ETA_R0;
    bssnParams.ETA_CONST         = bssn::ETA_CONST;
    bssnParams.ETA_DAMPING       = bssn::ETA_DAMPING;
    bssnParams.ETA_DAMPING_EXP   = bssn::ETA_DAMPING_EXP;
    bssnParams.KO_DISS_SIGMA     = bssn::KO_DISS_SIGMA;

    dim3 threadBlock(16, 16, 1);
    cuda::computeRHS(uzipVarsRHS, (const double **)uZipVars, blkList, numBlocks,
                     (const cuda::BSSNComputeParams *)&bssnParams, threadBlock,
                     pt_min, pt_max, 1, curr_time,
                     (const double **)uZipConstVars);
#else

    // Hybrid path: blocks are independent. Each thread uses its own deriv
    // workspace slab and its own DendroDerivatives clone (via active_derivs()),
    // so the block loop is race-free under DENDRO_HYBRID_OMP.
#ifdef DENDRO_HYBRID_OMP
// schedule(runtime) lets OMP_SCHEDULE pick the policy at run time without a
// rebuild: "dynamic,1" = fine-grained load balance (default); "static" =
// thread-owned contiguous block partition, which keeps each thread's private L2
// warm across RK stages (MPI-subdomain-like locality) at the cost of some load
// imbalance. A/B these to see which wins for a given mesh.
#pragma omp parallel for schedule(runtime)  \
    num_threads(bssn::BSSN_HYBRID_NTHREADS) \
    private(offset, sz, bflag, dx, dy, dz, ptmin, ptmax)
#endif
    for (unsigned int blk = 0; blk < numBlocks; blk++) {
        offset   = blkList[blk].getOffset();
        sz[0]    = blkList[blk].getAllocationSzX();
        sz[1]    = blkList[blk].getAllocationSzY();
        sz[2]    = blkList[blk].getAllocationSzZ();

        bflag    = blkList[blk].getBlkNodeFlag();

        dx       = blkList[blk].computeDx(pt_min, pt_max);
        dy       = blkList[blk].computeDy(pt_min, pt_max);
        dz       = blkList[blk].computeDz(pt_min, pt_max);

        ptmin[0] = GRIDX_TO_X(blkList[blk].getBlockNode().minX()) - PW * dx;
        ptmin[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().minY()) - PW * dy;
        ptmin[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ()) - PW * dz;

        ptmax[0] = GRIDX_TO_X(blkList[blk].getBlockNode().maxX()) + PW * dx;
        ptmax[1] = GRIDY_TO_Y(blkList[blk].getBlockNode().maxY()) + PW * dy;
        ptmax[2] = GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ()) + PW * dz;

        bssnrhs(uzipVarsRHS, (const double **)uZipVars, offset, ptmin, ptmax,
                sz, bflag, curr_time, (const double **)uZipConstVars);
    }
#endif
}

/*----------------------------------------------------------------------;
 *
 * vector form of RHS
 *
 *----------------------------------------------------------------------*/
void bssnrhs(double **unzipVarsRHS, const double **uZipVars,
             const unsigned int &offset, const double *pmin, const double *pmax,
             const unsigned int *sz, const unsigned int &bflag, const double t,
             const double **uZipConstVars) {
    const double *const alpha = &uZipVars[VAR::U_ALPHA][offset];
    const double *const chi   = &uZipVars[VAR::U_CHI][offset];
    const double *const K     = &uZipVars[VAR::U_K][offset];
    const double *const gt0   = &uZipVars[VAR::U_SYMGT0][offset];
    const double *const gt1   = &uZipVars[VAR::U_SYMGT1][offset];
    const double *const gt2   = &uZipVars[VAR::U_SYMGT2][offset];
    const double *const gt3   = &uZipVars[VAR::U_SYMGT3][offset];
    const double *const gt4   = &uZipVars[VAR::U_SYMGT4][offset];
    const double *const gt5   = &uZipVars[VAR::U_SYMGT5][offset];
    const double *const beta0 = &uZipVars[VAR::U_BETA0][offset];
    const double *const beta1 = &uZipVars[VAR::U_BETA1][offset];
    const double *const beta2 = &uZipVars[VAR::U_BETA2][offset];
    const double *const At0   = &uZipVars[VAR::U_SYMAT0][offset];
    const double *const At1   = &uZipVars[VAR::U_SYMAT1][offset];
    const double *const At2   = &uZipVars[VAR::U_SYMAT2][offset];
    const double *const At3   = &uZipVars[VAR::U_SYMAT3][offset];
    const double *const At4   = &uZipVars[VAR::U_SYMAT4][offset];
    const double *const At5   = &uZipVars[VAR::U_SYMAT5][offset];
    const double *const Gt0   = &uZipVars[VAR::U_GT0][offset];
    const double *const Gt1   = &uZipVars[VAR::U_GT1][offset];
    const double *const Gt2   = &uZipVars[VAR::U_GT2][offset];
    const double *const B0    = &uZipVars[VAR::U_B0][offset];
    const double *const B1    = &uZipVars[VAR::U_B1][offset];
    const double *const B2    = &uZipVars[VAR::U_B2][offset];

    double *const a_rhs       = &unzipVarsRHS[VAR::U_ALPHA][offset];
    double *const chi_rhs     = &unzipVarsRHS[VAR::U_CHI][offset];
    double *const K_rhs       = &unzipVarsRHS[VAR::U_K][offset];
    double *const gt_rhs00    = &unzipVarsRHS[VAR::U_SYMGT0][offset];
    double *const gt_rhs01    = &unzipVarsRHS[VAR::U_SYMGT1][offset];
    double *const gt_rhs02    = &unzipVarsRHS[VAR::U_SYMGT2][offset];
    double *const gt_rhs11    = &unzipVarsRHS[VAR::U_SYMGT3][offset];
    double *const gt_rhs12    = &unzipVarsRHS[VAR::U_SYMGT4][offset];
    double *const gt_rhs22    = &unzipVarsRHS[VAR::U_SYMGT5][offset];
    double *const b_rhs0      = &unzipVarsRHS[VAR::U_BETA0][offset];
    double *const b_rhs1      = &unzipVarsRHS[VAR::U_BETA1][offset];
    double *const b_rhs2      = &unzipVarsRHS[VAR::U_BETA2][offset];
    double *const At_rhs00    = &unzipVarsRHS[VAR::U_SYMAT0][offset];
    double *const At_rhs01    = &unzipVarsRHS[VAR::U_SYMAT1][offset];
    double *const At_rhs02    = &unzipVarsRHS[VAR::U_SYMAT2][offset];
    double *const At_rhs11    = &unzipVarsRHS[VAR::U_SYMAT3][offset];
    double *const At_rhs12    = &unzipVarsRHS[VAR::U_SYMAT4][offset];
    double *const At_rhs22    = &unzipVarsRHS[VAR::U_SYMAT5][offset];
    double *const Gt_rhs0     = &unzipVarsRHS[VAR::U_GT0][offset];
    double *const Gt_rhs1     = &unzipVarsRHS[VAR::U_GT1][offset];
    double *const Gt_rhs2     = &unzipVarsRHS[VAR::U_GT2][offset];
    double *const B_rhs0      = &unzipVarsRHS[VAR::U_B0][offset];
    double *const B_rhs1      = &unzipVarsRHS[VAR::U_B1][offset];
    double *const B_rhs2      = &unzipVarsRHS[VAR::U_B2][offset];

    // then the constraints (should be optimized out if not called)
    const double *const ham   = &uZipConstVars[VAR_CONSTRAINT::C_HAM][offset];
    const double *const mom0  = &uZipConstVars[VAR_CONSTRAINT::C_MOM0][offset];
    const double *const mom1  = &uZipConstVars[VAR_CONSTRAINT::C_MOM1][offset];
    const double *const mom2  = &uZipConstVars[VAR_CONSTRAINT::C_MOM2][offset];
    const double *const psi4_real =
        &uZipConstVars[VAR_CONSTRAINT::C_PSI4_REAL][offset];
    const double *const psi4_img =
        &uZipConstVars[VAR_CONSTRAINT::C_PSI4_IMG][offset];

    mem::memory_pool<double> *__mem_pool = &BSSN_MEM_POOL;

    const unsigned int nx                = sz[0];
    const unsigned int ny                = sz[1];
    const unsigned int nz                = sz[2];

    const double hx                      = (pmax[0] - pmin[0]) / (nx - 1);
    const double hy                      = (pmax[1] - pmin[1]) / (ny - 1);
    const double hz                      = (pmax[2] - pmin[2]) / (nz - 1);

    const unsigned int lambda[4]         = {BSSN_LAMBDA[0], BSSN_LAMBDA[1],
                                            BSSN_LAMBDA[2], BSSN_LAMBDA[3]};
    const double lambda_f[2]             = {BSSN_LAMBDA_F[0], BSSN_LAMBDA_F[1]};

    // for CAHD we need also need dt, dx_i, and dx_min
    const double dt                      = bssn::BSSN_RK45_TIME_STEP_SIZE;
    // dx_i is the current dx (which is our hx)
    const double dx_i                    = hx;
    // dx_min is the minimum resolution of the entire grid
    const double dx_min                  = bssn::BSSN_CURRENT_MIN_DX;
    // and then the new parameter
    const double BSSN_CAHD_C             = bssn::BSSN_CAHD_C;

    // data needed for the black hole location
    const double bhMass1                 = bssn::BH1.getBHMass();
    const double bhMass2                 = bssn::BH2.getBHMass();
    const double bh1x                    = bssn::BSSN_BH_LOC[0].x();
    const double bh1y                    = bssn::BSSN_BH_LOC[0].y();
    const double bh1z                    = bssn::BSSN_BH_LOC[0].z();
    const double bh2x                    = bssn::BSSN_BH_LOC[1].x();
    const double bh2y                    = bssn::BSSN_BH_LOC[1].y();
    const double bh2z                    = bssn::BSSN_BH_LOC[1].z();

    const unsigned int PW                = bssn::BSSN_PADDING_WIDTH;
    const unsigned int n                 = sz[0] * sz[1] * sz[2];

#ifdef DENDRO_USE_NEW_DERIVS
    // Decide once per block: does a puncture lie within it (interior bounds
    // expanded by NBLOCKS-1 block-widths)? If so, route this block's derivs
    // through the explicit fallback below. NBLOCKS==0 disables.
    bool puncture_block = false;
    if (bssn::BSSN_DERIV_PUNCTURE_EXPLICIT_NBLOCKS > 0) {
        const unsigned int Nb = bssn::BSSN_DERIV_PUNCTURE_EXPLICIT_NBLOCKS;
        for (unsigned int b = 0; b < 2 && !puncture_block; ++b) {
            const double bc[3] = {bssn::BSSN_BH_LOC[b].x(),
                                  bssn::BSSN_BH_LOC[b].y(),
                                  bssn::BSSN_BH_LOC[b].z()};
            bool inside = true;
            for (unsigned int d = 0; d < 3 && inside; ++d) {
                const double h_d  = (pmax[d] - pmin[d]) / (sz[d] - 1);
                const double lo_i = pmin[d] + PW * h_d;       // interior min
                const double hi_i = pmax[d] - PW * h_d;       // interior max
                const double mrg  = (Nb - 1) * (hi_i - lo_i);  // block-width margin
                if (bc[d] < lo_i - mrg || bc[d] > hi_i + mrg) inside = false;
            }
            puncture_block = inside;
        }
    }
#endif

    bssn::timer::t_deriv.start();

    const unsigned int BLK_SZ = n;
    // each thread gets its own workspace slab so the block loop can be threaded
    double *const deriv_base  = bssn::BSSN_DERIV_WORKSPACE
#ifdef DENDRO_HYBRID_OMP
        + (size_t)omp_get_thread_num() * bssn::BSSN_DERIV_WORKSPACE_STRIDE
#endif
        ;

    // clang-format off
    #include "bssnrhs_evar_derivs.h"
#ifdef DENDRO_USE_NEW_DERIVS
    if (puncture_block) set_block_explicit_derivs(true);
#endif
#if defined(BSSN_USE_CASCADE_AVX512_FUSED) ||      \
    defined(BSSN_USE_CASCADE_IR_AVX512_FUSED) ||   \
    defined(BSSN_USE_CASCADE_IR_AVX_FUSED)
    // Fused kernels (hand or IR): for bflag==0 blocks (interior), use
    // mixed-only precompute — wide and narrow fused kernels both work from
    // mixed-only arrays. Only bflag!=0 (boundary) blocks need the full
    // 138-array workspace (their non-fused fallback reads it).
    if (bflag == 0) {
        #include "bssnrhs_derivs_mixed_only.h"
    } else {
        #include "bssnrhs_derivs.h"
        #include "bssnrhs_derivs_adv.h"
    }
#elif defined(BSSN_USE_CASCADE_AVX_FUSED)
    // AVX2 fused: bflag==0 → mixed-only precompute; else full precompute
    if (bflag == 0) {
        #include "bssnrhs_derivs_mixed_only.h"
    } else {
        #include "bssnrhs_derivs.h"
        #include "bssnrhs_derivs_adv.h"
    }
#else
#ifdef BSSN_USE_GRAD_SET
    if (puncture_block) {
        #include "bssnrhs_derivs.h"
    } else {
        #include "bssnrhs_derivs_gradset.h"
    }
#else
    #include "bssnrhs_derivs.h"
#endif
    #include "bssnrhs_derivs_adv.h"
#endif
#ifdef DENDRO_USE_NEW_DERIVS
    // restore configured operators before KO/algebra (KO uses BSSN_DERIVS directly)
    if (puncture_block) set_block_explicit_derivs(false);
#endif
    // clang-format on

    bssn::timer::t_deriv.stop();

    // loop dep. removed allowing compiler to optmize for vectorization.
    // if (bssn::RIT_ETA_FUNCTION == 0) {
    //     // HAD eta function
    //     eta=ETA_CONST;
    //     if (r_coord >= ETA_R0) {
    //         eta *= pow( (ETA_R0/r_coord), ETA_DAMPING_EXP);
    //     }
    // }
    // cout << "begin loop" << endl;

    // the SSL parameters
    const double h_ssl   = bssn::BSSN_SSL_H;
    const double sig_ssl = bssn::BSSN_SSL_SIGMA;

    bssn::timer::t_rhs.start();
#ifdef BSSN_USE_CASCADE_AVX512_FUSED
#pragma message("BSSN: using AVX-512 fused cascade (8-wide)")
    if (bflag == 0 && (nx - 2 * PW) >= 8) {
// Wide interior: 8-wide AVX-512 fused
#include "bssn_cascade_avx512_fused_interior.inc.cpp"
    } else if (bflag == 0) {
// Narrow interior (width < 8): 4-wide AVX2 fused (reads mixed-only)
#include "bssn_cascade_avx_fused_interior.inc.cpp"
    } else {
// Boundary block (bflag != 0): full precompute already done,
// use non-fused AVX2 (reads full 138-array workspace)
#include "bssn_cascade_avx_interior.inc.cpp"
    }
#elif defined(BSSN_USE_CASCADE_AVX_FUSED)
#pragma message( \
    "BSSN: using AVX2-batched cascade with inline deriv stencils (fused)")
    if (bflag == 0) {
#include "bssn_cascade_avx_fused_interior.inc.cpp"
    } else {
// Boundary block: fused centered stencils are wrong at the 3 outer
// points; fall back to non-fused AVX cascade (reads the full 138-
// array workspace that bflag-aware deriv code populated above).
#include "bssn_cascade_avx_interior.inc.cpp"
    }
#elif defined(BSSN_USE_CASCADE_AVX)
#pragma message("BSSN: using AVX2-batched cascade RHS")
#include "bssn_cascade_avx_interior.inc.cpp"
#elif defined(BSSN_USE_CASCADE_IR_AVX512_FUSED)
#pragma message( \
    "BSSN: using IR-generated polynomial-cascade RHS (AVX-512, 8-wide, fused stencils)")
    // Fused IR body: inline 6th-order stencils for 1st/pure-2nd derivs,
    // mixed 2nds from the mixed-only pre-pass above. Boundary blocks
    // (bflag != 0) got the full pre-pass and use the non-fused IR path.
    if (bflag == 0 && (nx - 2 * PW) >= 8) {
#include "bssn_cascade_ir_avx512_fused_interior.inc.cpp"
    } else if (bflag == 0) {
#include "bssn_cascade_ir_avx2_fused_interior.inc.cpp"
    } else if ((nx - 2 * PW) >= 8) {
#include "bssn_cascade_ir_avx512_interior.inc.cpp"
    } else {
#include "bssn_cascade_ir_avx2_interior.inc.cpp"
    }
#elif defined(BSSN_USE_CASCADE_IR_AVX_FUSED)
#pragma message( \
    "BSSN: using IR-generated polynomial-cascade RHS (AVX2, 4-wide, fused stencils)")
    if (bflag == 0) {
#include "bssn_cascade_ir_avx2_fused_interior.inc.cpp"
    } else {
#include "bssn_cascade_ir_avx2_interior.inc.cpp"
    }
#elif defined(BSSN_USE_CASCADE_IR_AVX512)
#pragma message("BSSN: using IR-generated polynomial-cascade RHS (AVX-512, 8-wide)")
    // IR-generated body (codegen/cascade_emit.py): non-fused, reads the full
    // deriv workspace; SSL/CAHD gauge terms included via runtime coefficients
    // (h_ssl = 0 / BSSN_CAHD_C = 0 disable them without rebuild). Wide blocks
    // use the 8-wide kernel; narrower blocks the 4-wide twin.
    if ((nx - 2 * PW) >= 8) {
#include "bssn_cascade_ir_avx512_interior.inc.cpp"
    } else {
#include "bssn_cascade_ir_avx2_interior.inc.cpp"
    }
#elif defined(BSSN_USE_CASCADE_IR_AVX)
#pragma message("BSSN: using IR-generated polynomial-cascade RHS (AVX2, 4-wide)")
    // 4-wide twin of the IR_AVX512 path for AVX2-only hardware (AMD Zen <= 3,
    // Intel parts with AVX-512 fused off). Same generated body.
#include "bssn_cascade_ir_avx2_interior.inc.cpp"
#elif defined(BSSN_USE_CASCADE_AVX512)
#pragma message("BSSN: using AVX-512-batched cascade RHS (8-wide)")
    // Non-fused AVX-512: reads the full 138-array deriv workspace populated by
    // the (#else) pre-pass above, so it is valid for both interior and
    // boundary blocks. Wide blocks (interior x-extent >= 8) use the 8-wide
    // kernel; narrower blocks fall back to the 4-wide AVX2 non-fused kernel.
    if ((nx - 2 * PW) >= 8) {
#include "bssn_cascade_avx512_interior.inc.cpp"
    } else {
#include "bssn_cascade_avx_interior.inc.cpp"
    }
#else
    for (unsigned int k = PW; k < nz - PW; k++) {
        for (unsigned int j = PW; j < ny - PW; j++) {
            // clang-format off
#ifdef BSSN_ENABLE_AVX
  #ifdef __INTEL_COMPILER
    #pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
    #pragma ivdep
  #endif
#endif
            // clang-format on

            for (unsigned int i = PW; i < nx - PW; i++) {
                // pull current coordinates
                const double x        = pmin[0] + i * hx;
                const double y        = pmin[1] + j * hy;
                const double z        = pmin[2] + k * hz;

                // index for variables at current location
                const unsigned int pp = i + nx * (j + ny * k);

                // set up some commonly used distances for the eta fxn
                const double r_coord  = sqrt(x * x + y * y + z * z);
                const double dr1 =
                    sqrt((x - bh1x) * (x - bh1x) + (y - bh1y) * (y - bh1y) +
                         (z - bh1z) * (z - bh1z));
                const double dr2 =
                    sqrt((x - bh2x) * (x - bh2x) + (y - bh2y) * (y - bh2y) +
                         (z - bh2z) * (z - bh2z));
                // eta formulation
                // clang-format off
                // const double eta = bssn::ETA_CONST; // constant damping
                #include "eta_RIT.inc.cpp" // RIT's prescription
                // #include "eta_linear_inverse.inc.cpp"
                // #include "eta_tophat.inc.cpp"
                // #include "eta_outerfloor.inc.cpp"
                // #include "eta_G.inc.cpp"
                // #include "eta_tophat_grow.inc.cpp"
                // #include "eta_outerfloor_inpand.inc.cpp"
                // #include "eta_single.inc.cpp"
                // #include "eta_tophat_wide.inc.cpp"
                // #include "eta_causal_grow.inc.cpp"
                // #include "eta_causal_fade.inc.cpp"
                // #include "eta_global_ramp.inc.cpp"
                // #include "eta_causal_fade_tuned.inc.cpp"
                // clang-format on

                // clang-format off
#ifdef BSSN_USE_CASCADE
  #pragma message("BSSN: using cascade RHS (experimental)")
  #include "bssneqs_cascade.cpp"
#elif defined(BSSN_USE_NAIVE)
  #pragma message("BSSN: using math-faithful naive RHS (no CSE)")
  #include "bssneqs_naive.cpp"
#elif defined(BSSN_ENABLE_SSL_HD)
  #pragma message("BSSN: enabling both SSL and CAHD")
  // #include "bssn_eqns_SSL_HD.cpp"
  // #include "bssn_eqns_SSL_HD_HAM_INCLUDED.inc.cpp"
  #include "bssneqs_SSL_HD_dxsq.cpp" // use dx^2/(1+10*dx^2) in H-damping
  // #include "test_bssn_shock_56.inc.cpp" // shock-avoiding lapse: use 5/6 post-merger? 
  // #include "test_bssn_etaG.inc.cpp" // use eta_G exactly as LH23
  // #include "test_bssn_etaG_SSL_CAHD.inc.cpp" // use eta_G, evolve w/ B
  // #include "test_bssn_etaG_LH23_SSL_CAHD.inc.cpp" // LH23 formulation
  // #include "test_bssn_eta_G_adv_new.inc.cpp" // eta_G w/ advect
#else
  #pragma message( \
    "BSSN: SSL and HD is **NOT** enabled! Using original formalism!")
  #ifdef USE_ROCHESTER_GAUGE
    #pragma message("BSSN: using rochester gauge")
    #ifdef USE_ETA_FUNC
      #pragma message("BSSN: using function eta damping")
      #include "bssneqs_eta_func_rochester_gauge.cpp"
    #else
      #pragma message("BSSN: using const eta damping")
      #include "bssneqs_eta_const_rochester_gauge.cpp"
    #endif
  // else for USE_ROCHESTER_GAUGE
  #else
    #pragma message("BSSN: using standard gauge")
    #ifdef USE_ETA_FUNC
      #pragma message("BSSN: using function eta damping")
      #include "bssneqs_eta_func_standard_gauge.cpp"
    #else
      #pragma message("BSSN: using const eta damping")
      #include "bssneqs_eta_const_standard_gauge.cpp"
    #endif
  #endif
#endif
                // clang-format on
            }
        }
    }
#endif  // BSSN_USE_CASCADE_AVX* branches
    bssn::timer::t_rhs.stop();

    if (bflag != 0) {
        bssn::timer::t_bdyc.start();

        bssn::timer::t_rhs_a.start();
        bssn_bcs(a_rhs, alpha, grad_0_alpha, grad_1_alpha, grad_2_alpha, pmin,
                 pmax, 1.0, 1.0, sz, bflag);
        bssn::timer::t_rhs_a.stop();

        bssn::timer::t_rhs_chi.start();
        bssn_bcs(chi_rhs, chi, grad_0_chi, grad_1_chi, grad_2_chi, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn::timer::t_rhs_chi.stop();

        bssn::timer::t_rhs_K.start();
        bssn_bcs(K_rhs, K, grad_0_K, grad_1_K, grad_2_K, pmin, pmax, 1.0, 0.0,
                 sz, bflag);
        bssn::timer::t_rhs_K.stop();

        bssn::timer::t_rhs_b.start();
        bssn_bcs(b_rhs0, beta0, grad_0_beta0, grad_1_beta0, grad_2_beta0, pmin,
                 pmax, 1.0, 0.0, sz, bflag);
        bssn_bcs(b_rhs1, beta1, grad_0_beta1, grad_1_beta1, grad_2_beta1, pmin,
                 pmax, 1.0, 0.0, sz, bflag);
        bssn_bcs(b_rhs2, beta2, grad_0_beta2, grad_1_beta2, grad_2_beta2, pmin,
                 pmax, 1.0, 0.0, sz, bflag);
        bssn::timer::t_rhs_b.stop();

        bssn::timer::t_rhs_Gt.start();
        bssn_bcs(Gt_rhs0, Gt0, grad_0_Gt0, grad_1_Gt0, grad_2_Gt0, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(Gt_rhs1, Gt1, grad_0_Gt1, grad_1_Gt1, grad_2_Gt1, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(Gt_rhs2, Gt2, grad_0_Gt2, grad_1_Gt2, grad_2_Gt2, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn::timer::t_rhs_Gt.stop();

        bssn::timer::t_rhs_B.start();
        bssn_bcs(B_rhs0, B0, grad_0_B0, grad_1_B0, grad_2_B0, pmin, pmax, 1.0,
                 0.0, sz, bflag);
        bssn_bcs(B_rhs1, B1, grad_0_B1, grad_1_B1, grad_2_B1, pmin, pmax, 1.0,
                 0.0, sz, bflag);
        bssn_bcs(B_rhs2, B2, grad_0_B2, grad_1_B2, grad_2_B2, pmin, pmax, 1.0,
                 0.0, sz, bflag);
        bssn::timer::t_rhs_B.stop();

        bssn::timer::t_rhs_At.start();
        bssn_bcs(At_rhs00, At0, grad_0_At0, grad_1_At0, grad_2_At0, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs01, At1, grad_0_At1, grad_1_At1, grad_2_At1, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs02, At2, grad_0_At2, grad_1_At2, grad_2_At2, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs11, At3, grad_0_At3, grad_1_At3, grad_2_At3, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs12, At4, grad_0_At4, grad_1_At4, grad_2_At4, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn_bcs(At_rhs22, At5, grad_0_At5, grad_1_At5, grad_2_At5, pmin, pmax,
                 2.0, 0.0, sz, bflag);
        bssn::timer::t_rhs_At.stop();

        bssn::timer::t_rhs_gt.start();
        bssn_bcs(gt_rhs00, gt0, grad_0_gt0, grad_1_gt0, grad_2_gt0, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(gt_rhs01, gt1, grad_0_gt1, grad_1_gt1, grad_2_gt1, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs02, gt2, grad_0_gt2, grad_1_gt2, grad_2_gt2, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs11, gt3, grad_0_gt3, grad_1_gt3, grad_2_gt3, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn_bcs(gt_rhs12, gt4, grad_0_gt4, grad_1_gt4, grad_2_gt4, pmin, pmax,
                 1.0, 0.0, sz, bflag);
        bssn_bcs(gt_rhs22, gt5, grad_0_gt5, grad_1_gt5, grad_2_gt5, pmin, pmax,
                 1.0, 1.0, sz, bflag);
        bssn::timer::t_rhs_gt.stop();

        bssn::timer::t_bdyc.stop();
    }

    // KO dissipation, timed as one region in BOTH paths.
    //
    // This used to charge the two paths differently. The old path precomputes
    // per-direction KO derivatives via bssnrhs_ko_derivs.h and that include sat
    // inside t_deriv, while under DENDRO_USE_NEW_DERIVS filter_cako computes the
    // same derivatives internally, inside t_rhs_ko. t_deriv therefore carried KO
    // derivative cost on the old path and none on the new one, so deriv_t0 was
    // not comparable between the two -- and neither was deriv_t0 + rhs_ko_t0,
    // which gave the old path no KO term at all while charging the new path for
    // both the KO derivatives and their application.
    //
    // Now t_deriv covers scheme derivatives only, and t_rhs_ko covers all KO
    // work (derivatives and application) on both paths. t_rhs still lumps
    // interior + KO, so interior-only is (t_rhs - t_rhs_ko) as before.
    bssn::timer::t_rhs.start();
    bssn::timer::t_rhs_ko.start();
#ifndef DENDRO_USE_NEW_DERIVS
    // Old path: precompute per-direction KO derivatives into the grad_* slots.
    // With DENDRO_USE_NEW_DERIVS, filter_cako computes them internally instead.
#include "bssnrhs_ko_derivs.h"
#endif

    double sigma = KO_DISS_SIGMA;

#ifdef DENDRO_USE_NEW_DERIVS
    // KO dissipation through dendrolib's DendroDerivatives::filter_cako: it
    // computes the KO derivative of each variable internally and accumulates
    // coeff_field * (KOx + KOy + KOz) into the corresponding *_rhs. CAKO is
    // reproduced by passing a per-point sqrt(chi)*epsilon field (gauge vs
    // other); otherwise a uniform KO_DISS_SIGMA field. Scratch buffers are the
    // (now-free) interior-RHS derivative slots — never aliased with an evolved
    // variable or its *_rhs.
    {
        const unsigned int npts = nx * ny * nz;
        double* kf = grad_0_K;    // per-point coefficient field
        double* wx = grad_1_K;    // KO workspaces
        double* wy = grad_2_K;
        double* wz = grad_0_Gt0;
        auto ko_apply = [&](const double* const in, double* const rhs) {
            bssn::active_derivs()->filter_cako(in, rhs, wx, wy, wz, hx, hy, hz, kf,
                                           sz, bflag);
        };

        if (bssn::BSSN_CAKO_ENABLED) {
            for (unsigned int p = 0; p < npts; p++)
                kf[p] = sqrt(chi[p]) * bssn::BSSN_EPSILON_CAKO_GAUGE;
            ko_apply(alpha, a_rhs);
            ko_apply(beta0, b_rhs0);
            ko_apply(beta1, b_rhs1);
            ko_apply(beta2, b_rhs2);
            ko_apply(B0, B_rhs0);
            ko_apply(B1, B_rhs1);
            ko_apply(B2, B_rhs2);

            for (unsigned int p = 0; p < npts; p++)
                kf[p] = sqrt(chi[p]) * bssn::BSSN_EPSILON_CAKO_OTHER;
            ko_apply(gt0, gt_rhs00);
            ko_apply(gt1, gt_rhs01);
            ko_apply(gt2, gt_rhs02);
            ko_apply(gt3, gt_rhs11);
            ko_apply(gt4, gt_rhs12);
            ko_apply(gt5, gt_rhs22);
            ko_apply(chi, chi_rhs);
            ko_apply(At0, At_rhs00);
            ko_apply(At1, At_rhs01);
            ko_apply(At2, At_rhs02);
            ko_apply(At3, At_rhs11);
            ko_apply(At4, At_rhs12);
            ko_apply(At5, At_rhs22);
            ko_apply(K, K_rhs);
            ko_apply(Gt0, Gt_rhs0);
            ko_apply(Gt1, Gt_rhs1);
            ko_apply(Gt2, Gt_rhs2);
        } else {
            for (unsigned int p = 0; p < npts; p++) kf[p] = sigma;
            ko_apply(alpha, a_rhs);
            ko_apply(beta0, b_rhs0);
            ko_apply(beta1, b_rhs1);
            ko_apply(beta2, b_rhs2);
            ko_apply(B0, B_rhs0);
            ko_apply(B1, B_rhs1);
            ko_apply(B2, B_rhs2);
            ko_apply(gt0, gt_rhs00);
            ko_apply(gt1, gt_rhs01);
            ko_apply(gt2, gt_rhs02);
            ko_apply(gt3, gt_rhs11);
            ko_apply(gt4, gt_rhs12);
            ko_apply(gt5, gt_rhs22);
            ko_apply(chi, chi_rhs);
            ko_apply(At0, At_rhs00);
            ko_apply(At1, At_rhs01);
            ko_apply(At2, At_rhs02);
            ko_apply(At3, At_rhs11);
            ko_apply(At4, At_rhs12);
            ko_apply(At5, At_rhs22);
            ko_apply(K, K_rhs);
            ko_apply(Gt0, Gt_rhs0);
            ko_apply(Gt1, Gt_rhs1);
            ko_apply(Gt2, Gt_rhs2);
        }
    }
#else
    for (unsigned int k = PW; k < nz - PW; k++) {
        for (unsigned int j = PW; j < ny - PW; j++) {
#ifdef BSSN_ENABLE_AVX
#ifdef __INTEL_COMPILER
#pragma vector vectorlength(__RHS_AVX_SIMD_LEN__) vecremainder
#pragma ivdep
#endif
#endif
            for (unsigned int i = PW; i < nx - PW; i++) {
                const unsigned int pp = i + nx * (j + ny * k);

                // as part of the improved techniques paper
                // (https://arxiv.org/pdf/2404.01137.pdf) they mentioned scaling
                // the KO dissipation of the RHS by sqrt(chi) with a strong
                // amount for the gauge variaables and a smaller amount for the
                // non-gauge variables. This is an option the user can use.
                if (bssn::BSSN_CAKO_ENABLED) {
                    sigma              = sqrt(chi[pp]);
                    double sigma_gauge = sigma * bssn::BSSN_EPSILON_CAKO_GAUGE;
                    double sigma_other = sigma * bssn::BSSN_EPSILON_CAKO_OTHER;
                    a_rhs[pp] +=
                        sigma_gauge * (grad_0_alpha[pp] + grad_1_alpha[pp] +
                                       grad_2_alpha[pp]);
                    b_rhs0[pp] +=
                        sigma_gauge * (grad_0_beta0[pp] + grad_1_beta0[pp] +
                                       grad_2_beta0[pp]);
                    b_rhs1[pp] +=
                        sigma_gauge * (grad_0_beta1[pp] + grad_1_beta1[pp] +
                                       grad_2_beta1[pp]);
                    b_rhs2[pp] +=
                        sigma_gauge * (grad_0_beta2[pp] + grad_1_beta2[pp] +
                                       grad_2_beta2[pp]);

                    gt_rhs00[pp] +=
                        sigma_other *
                        (grad_0_gt0[pp] + grad_1_gt0[pp] + grad_2_gt0[pp]);
                    gt_rhs01[pp] +=
                        sigma_other *
                        (grad_0_gt1[pp] + grad_1_gt1[pp] + grad_2_gt1[pp]);
                    gt_rhs02[pp] +=
                        sigma_other *
                        (grad_0_gt2[pp] + grad_1_gt2[pp] + grad_2_gt2[pp]);
                    gt_rhs11[pp] +=
                        sigma_other *
                        (grad_0_gt3[pp] + grad_1_gt3[pp] + grad_2_gt3[pp]);
                    gt_rhs12[pp] +=
                        sigma_other *
                        (grad_0_gt4[pp] + grad_1_gt4[pp] + grad_2_gt4[pp]);
                    gt_rhs22[pp] +=
                        sigma_other *
                        (grad_0_gt5[pp] + grad_1_gt5[pp] + grad_2_gt5[pp]);

                    chi_rhs[pp] +=
                        sigma_other *
                        (grad_0_chi[pp] + grad_1_chi[pp] + grad_2_chi[pp]);

                    At_rhs00[pp] +=
                        sigma_other *
                        (grad_0_At0[pp] + grad_1_At0[pp] + grad_2_At0[pp]);
                    At_rhs01[pp] +=
                        sigma_other *
                        (grad_0_At1[pp] + grad_1_At1[pp] + grad_2_At1[pp]);
                    At_rhs02[pp] +=
                        sigma_other *
                        (grad_0_At2[pp] + grad_1_At2[pp] + grad_2_At2[pp]);
                    At_rhs11[pp] +=
                        sigma_other *
                        (grad_0_At3[pp] + grad_1_At3[pp] + grad_2_At3[pp]);
                    At_rhs12[pp] +=
                        sigma_other *
                        (grad_0_At4[pp] + grad_1_At4[pp] + grad_2_At4[pp]);
                    At_rhs22[pp] +=
                        sigma_other *
                        (grad_0_At5[pp] + grad_1_At5[pp] + grad_2_At5[pp]);

                    K_rhs[pp] += sigma_other *
                                 (grad_0_K[pp] + grad_1_K[pp] + grad_2_K[pp]);

                    Gt_rhs0[pp] +=
                        sigma_other *
                        (grad_0_Gt0[pp] + grad_1_Gt0[pp] + grad_2_Gt0[pp]);
                    Gt_rhs1[pp] +=
                        sigma_other *
                        (grad_0_Gt1[pp] + grad_1_Gt1[pp] + grad_2_Gt1[pp]);
                    Gt_rhs2[pp] +=
                        sigma_other *
                        (grad_0_Gt2[pp] + grad_1_Gt2[pp] + grad_2_Gt2[pp]);

                    B_rhs0[pp] += sigma_gauge * (grad_0_B0[pp] + grad_1_B0[pp] +
                                                 grad_2_B0[pp]);
                    B_rhs1[pp] += sigma_gauge * (grad_0_B1[pp] + grad_1_B1[pp] +
                                                 grad_2_B1[pp]);
                    B_rhs2[pp] += sigma_gauge * (grad_0_B2[pp] + grad_1_B2[pp] +
                                                 grad_2_B2[pp]);
                } else {
                    a_rhs[pp] += sigma * (grad_0_alpha[pp] + grad_1_alpha[pp] +
                                          grad_2_alpha[pp]);
                    b_rhs0[pp] += sigma * (grad_0_beta0[pp] + grad_1_beta0[pp] +
                                           grad_2_beta0[pp]);
                    b_rhs1[pp] += sigma * (grad_0_beta1[pp] + grad_1_beta1[pp] +
                                           grad_2_beta1[pp]);
                    b_rhs2[pp] += sigma * (grad_0_beta2[pp] + grad_1_beta2[pp] +
                                           grad_2_beta2[pp]);

                    gt_rhs00[pp] += sigma * (grad_0_gt0[pp] + grad_1_gt0[pp] +
                                             grad_2_gt0[pp]);
                    gt_rhs01[pp] += sigma * (grad_0_gt1[pp] + grad_1_gt1[pp] +
                                             grad_2_gt1[pp]);
                    gt_rhs02[pp] += sigma * (grad_0_gt2[pp] + grad_1_gt2[pp] +
                                             grad_2_gt2[pp]);
                    gt_rhs11[pp] += sigma * (grad_0_gt3[pp] + grad_1_gt3[pp] +
                                             grad_2_gt3[pp]);
                    gt_rhs12[pp] += sigma * (grad_0_gt4[pp] + grad_1_gt4[pp] +
                                             grad_2_gt4[pp]);
                    gt_rhs22[pp] += sigma * (grad_0_gt5[pp] + grad_1_gt5[pp] +
                                             grad_2_gt5[pp]);

                    chi_rhs[pp] += sigma * (grad_0_chi[pp] + grad_1_chi[pp] +
                                            grad_2_chi[pp]);

                    At_rhs00[pp] += sigma * (grad_0_At0[pp] + grad_1_At0[pp] +
                                             grad_2_At0[pp]);
                    At_rhs01[pp] += sigma * (grad_0_At1[pp] + grad_1_At1[pp] +
                                             grad_2_At1[pp]);
                    At_rhs02[pp] += sigma * (grad_0_At2[pp] + grad_1_At2[pp] +
                                             grad_2_At2[pp]);
                    At_rhs11[pp] += sigma * (grad_0_At3[pp] + grad_1_At3[pp] +
                                             grad_2_At3[pp]);
                    At_rhs12[pp] += sigma * (grad_0_At4[pp] + grad_1_At4[pp] +
                                             grad_2_At4[pp]);
                    At_rhs22[pp] += sigma * (grad_0_At5[pp] + grad_1_At5[pp] +
                                             grad_2_At5[pp]);

                    K_rhs[pp] +=
                        sigma * (grad_0_K[pp] + grad_1_K[pp] + grad_2_K[pp]);

                    Gt_rhs0[pp] += sigma * (grad_0_Gt0[pp] + grad_1_Gt0[pp] +
                                            grad_2_Gt0[pp]);
                    Gt_rhs1[pp] += sigma * (grad_0_Gt1[pp] + grad_1_Gt1[pp] +
                                            grad_2_Gt1[pp]);
                    Gt_rhs2[pp] += sigma * (grad_0_Gt2[pp] + grad_1_Gt2[pp] +
                                            grad_2_Gt2[pp]);

                    B_rhs0[pp] +=
                        sigma * (grad_0_B0[pp] + grad_1_B0[pp] + grad_2_B0[pp]);
                    B_rhs1[pp] +=
                        sigma * (grad_0_B1[pp] + grad_1_B1[pp] + grad_2_B1[pp]);
                    B_rhs2[pp] +=
                        sigma * (grad_0_B2[pp] + grad_1_B2[pp] + grad_2_B2[pp]);
                }
            }
        }
    }
#endif  // DENDRO_USE_NEW_DERIVS (KO via filter_cako vs explicit loop)

    bssn::timer::t_rhs_ko.stop();
    bssn::timer::t_rhs.stop();

    bssn::timer::t_deriv.start();
    bssn::timer::t_deriv.stop();

#if 0
        for (unsigned int m = 0; m < 24; m++) {
            std::cout<<"  || dtu("<<m<<")|| = "<<normLInfty(unzipVarsRHS[m] + offset, n)<<std::endl;
        }
#endif
}

void bssn_bcs(double *f_rhs, const double *f, const double *dxf,
              const double *dyf, const double *dzf, const double *pmin,
              const double *pmax, const double f_falloff,
              const double f_asymptotic, const unsigned int *sz,
              const unsigned int &bflag) {
    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    const double hx       = (pmax[0] - pmin[0]) / (nx - 1);
    const double hy       = (pmax[1] - pmin[1]) / (ny - 1);
    const double hz       = (pmax[2] - pmin[2]) / (nz - 1);

    const unsigned int PW = bssn::BSSN_PADDING_WIDTH;

    unsigned int ib       = PW;
    unsigned int jb       = PW;
    unsigned int kb       = PW;
    unsigned int ie       = sz[0] - PW;
    unsigned int je       = sz[1] - PW;
    unsigned int ke       = sz[2] - PW;

    double x, y, z;
    unsigned int pp;
    double inv_r;

    if (bflag & (1u << OCT_DIR_LEFT)) {
        double x = pmin[0] + ib * hx;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k * hz;
            for (unsigned int j = jb; j < je; j++) {
                y         = pmin[1] + j * hy;
                pp        = IDX(ib, j, k);
                inv_r     = 1.0 / sqrt(x * x + y * y + z * z);
                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
            }
        }
    }
    if (bflag & (1u << OCT_DIR_RIGHT)) {
        x = pmin[0] + (ie - 1) * hx;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k * hz;
            for (unsigned int j = jb; j < je; j++) {
                y         = pmin[1] + j * hy;
                pp        = IDX((ie - 1), j, k);
                inv_r     = 1.0 / sqrt(x * x + y * y + z * z);
                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
            }
        }
    }
    if (bflag & (1u << OCT_DIR_DOWN)) {
        y = pmin[1] + jb * hy;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k * hz;
            for (unsigned int i = ib; i < ie; i++) {
                x         = pmin[0] + i * hx;
                inv_r     = 1.0 / sqrt(x * x + y * y + z * z);
                pp        = IDX(i, jb, k);
                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
            }
        }
    }
    if (bflag & (1u << OCT_DIR_UP)) {
        y = pmin[1] + (je - 1) * hy;
        for (unsigned int k = kb; k < ke; k++) {
            z = pmin[2] + k * hz;
            for (unsigned int i = ib; i < ie; i++) {
                x         = pmin[0] + i * hx;
                inv_r     = 1.0 / sqrt(x * x + y * y + z * z);
                pp        = IDX(i, (je - 1), k);
                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
            }
        }
    }
    if (bflag & (1u << OCT_DIR_BACK)) {
        z = pmin[2] + kb * hz;
        for (unsigned int j = jb; j < je; j++) {
            y = pmin[1] + j * hy;
            for (unsigned int i = ib; i < ie; i++) {
                x         = pmin[0] + i * hx;
                inv_r     = 1.0 / sqrt(x * x + y * y + z * z);
                pp        = IDX(i, j, kb);
                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
            }
        }
    }
    if (bflag & (1u << OCT_DIR_FRONT)) {
        z = pmin[2] + (ke - 1) * hz;
        for (unsigned int j = jb; j < je; j++) {
            y = pmin[1] + j * hy;
            for (unsigned int i = ib; i < ie; i++) {
                x         = pmin[0] + i * hx;
                inv_r     = 1.0 / sqrt(x * x + y * y + z * z);
                pp        = IDX(i, j, (ke - 1));
                f_rhs[pp] = -inv_r * (x * dxf[pp] + y * dyf[pp] + z * dzf[pp] +
                                      f_falloff * (f[pp] - f_asymptotic));
            }
        }
    }
}

void max_spacetime_speeds(double *const lambda1max, double *const lambda2max,
                          double *const lambda3max, const double *const alpha,
                          const double *const beta1, const double *const beta2,
                          const double *const beta3, const double *const gtd11,
                          const double *const gtd12, const double *const gtd13,
                          const double *const gtd22, const double *const gtd23,
                          const double *const gtd33, const double *const chi,
                          const unsigned int *sz) {
    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    const unsigned int PW = bssn::BSSN_PADDING_WIDTH;
    unsigned int ib       = PW;
    unsigned int jb       = PW;
    unsigned int kb       = PW;
    unsigned int ie       = sz[0] - PW;
    unsigned int je       = sz[1] - PW;
    unsigned int ke       = sz[2] - PW;

    for (unsigned int k = kb; k < ke; k++) {
        for (unsigned int j = jb; j < je; j++) {
            for (unsigned int i = ib; i < ie; i++) {
                unsigned int pp = IDX(i, j, k);
                /* note: gtu is the inverse tilde metric. It should have detgtd
                 * = 1. So, for the purposes of calculating wavespeeds, I simple
                 * set detgtd = 1. */
                double gtu11    = gtd22[pp] * gtd33[pp] - gtd23[pp] * gtd23[pp];
                double gtu22    = gtd11[pp] * gtd33[pp] - gtd13[pp] * gtd13[pp];
                double gtu33    = gtd11[pp] * gtd22[pp] - gtd12[pp] * gtd12[pp];
                if (gtu11 < 0.0 || gtu22 < 0.0 || gtu33 < 0.0) {
                    std::cout << "Problem computing spacetime characteristics"
                              << std::endl;
                    std::cout << "gtu11 = " << gtu11 << ", gtu22 = " << gtu22
                              << ", gtu33 = " << gtu33 << std::endl;
                    gtu11 = 1.0;
                    gtu22 = 1.0;
                    gtu33 = 1.0;
                }
                double t1 = alpha[pp] * sqrt(gtu11 * chi[pp]);
                double t2 = alpha[pp] * sqrt(gtu22 * chi[pp]);
                double t3 = alpha[pp] * sqrt(gtu33 * chi[pp]);
                lambda1max[pp] =
                    std::max(abs(-beta1[pp] + t1), abs(-beta1[pp] - t1));
                lambda2max[pp] =
                    std::max(abs(-beta2[pp] + t2), abs(-beta2[pp] - t2));
                lambda3max[pp] =
                    std::max(abs(-beta3[pp] + t3), abs(-beta3[pp] - t3));
            }
        }
    }
}

/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
void freeze_bcs(double *f_rhs, const unsigned int *sz,
                const unsigned int &bflag) {
    const unsigned int nx = sz[0];
    const unsigned int ny = sz[1];
    const unsigned int nz = sz[2];

    const unsigned int PW = bssn::BSSN_PADDING_WIDTH;
    unsigned int ib       = PW;
    unsigned int jb       = PW;
    unsigned int kb       = PW;
    unsigned int ie       = sz[0] - PW;
    unsigned int je       = sz[1] - PW;
    unsigned int ke       = sz[2] - PW;

    unsigned int pp;

    if (bflag & (1u << OCT_DIR_LEFT)) {
        for (unsigned int k = kb; k < ke; k++) {
            for (unsigned int j = jb; j < je; j++) {
                pp        = IDX(ib, j, k);
                f_rhs[pp] = 0.0;
            }
        }
    }

    if (bflag & (1u << OCT_DIR_RIGHT)) {
        for (unsigned int k = kb; k < ke; k++) {
            for (unsigned int j = jb; j < je; j++) {
                pp        = IDX((ie - 1), j, k);
                f_rhs[pp] = 0.0;
            }
        }
    }

    if (bflag & (1u << OCT_DIR_DOWN)) {
        for (unsigned int k = kb; k < ke; k++) {
            for (unsigned int i = ib; i < ie; i++) {
                pp        = IDX(i, jb, k);
                f_rhs[pp] = 0.0;
            }
        }
    }

    if (bflag & (1u << OCT_DIR_UP)) {
        for (unsigned int k = kb; k < ke; k++) {
            for (unsigned int i = ib; i < ie; i++) {
                pp        = IDX(i, (je - 1), k);
                f_rhs[pp] = 0.0;
            }
        }
    }

    if (bflag & (1u << OCT_DIR_BACK)) {
        for (unsigned int j = jb; j < je; j++) {
            for (unsigned int i = ib; i < ie; i++) {
                pp        = IDX(i, j, kb);
                f_rhs[pp] = 0.0;
            }
        }
    }

    if (bflag & (1u << OCT_DIR_FRONT)) {
        for (unsigned int j = jb; j < je; j++) {
            for (unsigned int i = ib; i < ie; i++) {
                pp        = IDX(i, j, (ke - 1));
                f_rhs[pp] = 0.0;
            }
        }
    }
}

/*----------------------------------------------------------------------;
 *
 * HAD RHS
 *
 *----------------------------------------------------------------------*/
void call_HAD_rhs() { had_bssn_rhs_(); }

#if 0
/*--------------------------------------------------------------
 * Kerr-Schild data
 *--------------------------------------------------------------*/

void ks_initial_data(double x, double y, double z, double *u)
{

    u[VAR::U_ALPHA] = 0.0;
    u[VAR::U_BETA0] = 0.0;
    u[VAR::U_BETA1] = pi/5.0*cos(y)*sin(z+x);
    u[VAR::U_BETA2] = 4.0/17.0*sin(f2*x)*sin(z);

    u[VAR::U_B0] = 0.0;
    u[VAR::U_B1] = 0.0;
    u[VAR::U_B2] = 0.0;

    u[VAR::U_GT0] = Gamt_1;
    u[VAR::U_GT1] = Gamt_2;
    u[VAR::U_GT2] = Gamt_3;

    u[VAR::U_CHI] = 1.0 + exp(-4.0*cos(x)*sin(y));

    u[VAR::U_SYMGT0] = 1.00+0.2*sin(x+z)*cos(y);
    u[VAR::U_SYMGT3] = 1.00+0.2*cos(y)*cos(z+ x);
    u[VAR::U_SYMGT5] = 1.00 / ( u[VAR::U_SYMGT0] + u[VAR::U_SYMGT3]);
    u[VAR::U_SYMGT1] = 0.7*cos(x*x + y*y);
    u[VAR::U_SYMGT2] = 0.3*sin(z)*cos(x);
    u[VAR::U_SYMGT4] = -0.5*sin(x*x)*cos(y)*cos(z);

    u[VAR::U_K] = 5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)
                  +5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)
                  +0.4*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))
                  *exp(-4.0*cos(x)*sin(y))*cos(z);

    u[VAR::U_SYMAT0] = exp(-4.0*cos(x)*sin(y))*(cos(x) -0.3333333333*exp(4.0*cos(x)*sin(y)) *(1.0+0.2*sin(x))*(5.0*exp(-4.0*cos(x)*sin(y)) /(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y)) /(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y) +5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
    u[VAR::U_SYMAT1] = 1.0 + x*z/(0.1 + x*x + y*y + z*z);
    u[VAR::U_SYMAT2] = 1.3 - x*y/(3.0 + x*x + 2.0*y*y + z*z)*(x*x+z*z);
    u[VAR::U_SYMAT3] = exp(-4.0*cos(x)*sin(y))*(cos(y)-0.33333333330*exp(4*cos(x)*sin(y))*(1+0.2*cos(y))*(5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
    u[VAR::U_SYMAT4] = -1.0 + y*z/(1.0 + 3.0*x*x + y*y + z*z);
    u[VAR::U_SYMAT5] = exp(-4.0*cos(x)*sin(y))*(cos(z)-0.3333333333*exp(4*cos(x)*sin(y))/(1+0.2*sin(x))/(1+0.2*cos(y))*(5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));

}
#endif
/*----------------------------------------------------------------------;
 *
 *
 *
 *----------------------------------------------------------------------*/
