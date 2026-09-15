//
// Created by milinda on 1/16/19.
//

#include "dataUtils.h"

#include <cmath>
#include <iomanip>
#include <limits>
#include <utility>

#include "bssnCtx.h"
#include "parameters.h"

namespace bssn {

void extractBHCoords(const ot::Mesh* pMesh, const DendroScalar* var,
                     double tolerance, const Point* ptIn, unsigned int numPt,
                     Point* ptOut) {
    if ((pMesh->isActive())) {
        MPI_Comm commActive     = pMesh->getMPICommunicator();
        unsigned int rankActive = pMesh->getMPIRank();
        unsigned int npesActive = pMesh->getMPICommSize();

        double v_min = vecMin((DendroScalar*)(var + pMesh->getNodeLocalBegin()),
                              pMesh->getNumLocalMeshNodes(), commActive);
        par::Mpi_Bcast(&v_min, 1, 0, commActive);

        assert(numPt == 2);

        const double extraction_tol     = tolerance;  // 10*v_min;

        const ot::TreeNode* allElements = &(*(pMesh->getAllElements().begin()));
        const unsigned int* e2n         = &(*(pMesh->getE2NMapping().begin()));
        const unsigned int* e2n_dg = &(*(pMesh->getE2NMapping_DG().begin()));
        const unsigned int* cgToDg = &(*(pMesh->getCG2DGMap().begin()));

        unsigned int lookup        = 0;
        unsigned int ownerID, ii_x, jj_y, kk_z;

        ot::TreeNode tmpOct;
        const unsigned int eleOrder = pMesh->getElementOrder();
        double hx, x, y, z;

        std::vector<Point> ptList;
        for (unsigned int node = pMesh->getNodeLocalBegin();
             node < pMesh->getNodeLocalEnd(); node++) {
            if (var[node] < extraction_tol) {
                lookup = cgToDg[node];
                pMesh->dg2eijk(lookup, ownerID, ii_x, jj_y, kk_z);
                tmpOct = allElements[ownerID];
                hx     = (tmpOct.maxX() - tmpOct.minX()) / ((double)eleOrder);
                x      = tmpOct.minX() + ii_x * (hx);
                y      = tmpOct.minY() + jj_y * (hx);
                z      = tmpOct.minZ() + kk_z * (hx);

                ptList.push_back(
                    Point(GRIDX_TO_X(x), GRIDY_TO_Y(y), GRIDZ_TO_Z(z)));

                // std::cout<<" x: "<<GRIDX_TO_X(x)<<" y: "<<GRIDY_TO_Y(y)<<" z:
                // "<<GRIDZ_TO_Z(z)<<std::endl;
            }
        }

        std::vector<Point>* ptCluster = new std::vector<Point>[numPt];

        double min0, min1;
        unsigned int cID;
        for (unsigned int pt = 0; pt < ptList.size(); pt++) {
            cID  = 0;
            min0 = (ptIn[0] - ptList[pt]).abs();
            min1 = (ptIn[1] - ptList[pt]).abs();

            if (fabs(min0 - min1) <
                1e-6) {  // implies the point is closer to the both clusters
                ptCluster[0].push_back(ptList[pt]);
                ptCluster[1].push_back(ptList[pt]);
            } else if (min0 < min1) {
                ptCluster[0].push_back(ptList[pt]);
            } else {
                ptCluster[1].push_back(ptList[pt]);
            }
        }

        ptList.clear();
        Point* ptMean          = new Point[numPt];
        DendroIntL* ptCounts   = new DendroIntL[numPt];
        DendroIntL* ptCounts_g = new DendroIntL[numPt];

        for (unsigned int c = 0; c < numPt; c++) {
            ptMean[c] = Point(0, 0, 0);
            ptOut[c]  = Point(0, 0, 0);
        }

        for (unsigned int c = 0; c < numPt; c++) {
            ptCounts[c] = ptCluster[c].size();
            for (unsigned int pt = 0; pt < ptCluster[c].size(); pt++)
                ptMean[c] += ptCluster[c][pt];
        }

        par::Mpi_Allreduce(ptCounts, ptCounts_g, numPt, MPI_SUM, commActive);
        par::Mpi_Allreduce(ptMean, ptOut, numPt,
                           par::Mpi_datatype<Point>::_SUM(), commActive);

        for (unsigned int c = 0; c < numPt; c++)
            ptOut[c] /= (double)ptCounts_g[c];

        // if(pMesh->getMPIRank()==0)std::cout<<"bh1 in : "<<ptIn[0].x()<<",
        // "<<ptIn[0].y()<<", "<<ptIn[0].z()<<std::endl;
        // if(pMesh->getMPIRank()==0)std::cout<<"bh2 in : "<<ptIn[1].x()<<",
        // "<<ptIn[1].y()<<", "<<ptIn[1].z()<<std::endl;

        // if(pMesh->getMPIRank()==0)std::cout<<"bh1 out: "<<ptOut[0].x()<<",
        // "<<ptOut[0].y()<<", "<<ptOut[0].z()<<std::endl;
        // if(pMesh->getMPIRank()==0)std::cout<<"bh2 out: "<<ptOut[1].x()<<",
        // "<<ptOut[1].y()<<", "<<ptOut[1].z()<<std::endl;

        delete[] ptCounts;
        delete[] ptCounts_g;
        delete[] ptMean;
        delete[] ptCluster;
    }
}

void writeBHCoordinates(const ot::Mesh* pMesh, const Point* ptLocs,
                        unsigned int numPt, unsigned int timestep,
                        double time) {
    unsigned int rankGlobal = pMesh->getMPIRankGlobal();
    if (!rankGlobal) {
        std::ofstream fileGW;
        char fName[256];
        sprintf(fName, "%s_BHLocations.dat",
                bssn::BSSN_PROFILE_FILE_PREFIX.c_str());
        fileGW.open(fName, std::ofstream::app);

        // writes the header
        if (timestep == 0)
            fileGW << "TimeStep\t" << "time\t" << "bh1_x\t" << "bh1_y\t"
                   << "bh1_z\t" << "bh2_x\t" << "bh2_y\t" << "bh2_z\t"
                   << std::endl;

        fileGW << std::setprecision(std::numeric_limits<double>::max_digits10)
               << timestep << "\t" << time << "\t" << ptLocs[0].x() << "\t"
               << ptLocs[0].y() << "\t" << ptLocs[0].z() << "\t"
               << ptLocs[1].x() << "\t" << ptLocs[1].y() << "\t"
               << ptLocs[1].z() << std::endl;
        fileGW.close();
        return;
    }
}

// The BH-derived kinematics that used to live here as free functions now live
// in dendro_bh::BHHistory (dendrolib, include/bh_history.h), which maintains
// them incrementally. isRemeshBH consumes that object directly.

// standard logistic transition from y_init to y_final over timescale sigma,
// centered about x0
static inline double logisticTransition(double x, double x0, double sigma,
                                        double y_init = 0.0,
                                        double y_final = 1.0) {
    const double s = 1.0 / (1.0 + std::exp(-(x - x0) / sigma));
    return y_init + (y_final - y_init) * s;
}

// logistic transition of the onion AMR ratio, tuned to the gauge wave leaving
// the BHs (retarded time tau = t - r / gauge_speed)
static inline double onionRatioGaugeLogistic(
    double t, double r, double y_final = bssn::BSSN_AMR_R_RATIO,
    double y_init = 2.0, double t_transition = 50.0, double sigma = 5.0,
    double gauge_speed = std::sqrt(2.0)) {
    const double tau = t - r / gauge_speed;
    return logisticTransition(tau, t_transition, sigma, y_init, y_final);
}

// ---------------------------------------------------------------------------
// Unified parallel refinement driver (hybrid-OpenMP, turn-off-able).
//
// refineFlagsPass runs `decide` for every LOCAL element (block-parallel, with a
// per-thread WaveletEl + scratch) and stores the returned octant flag in
// refine_flags. `decide(ctx, currentFlag)` receives element geometry + a
// thread-safe wavelet-error helper and the element's current flag, so passes
// compose (e.g. a wavelet decision followed by a BH-proximity override).
// Physicists edit only the per-element lambda; the threading lives here.
//
// With DENDRO_HYBRID_OMP OFF this is a plain serial loop with one shared scratch
// set -> bit-identical results, no OpenMP constructs compiled in.
// ---------------------------------------------------------------------------
struct RefineElemCtx {
    ot::Mesh* mesh;
    const ot::TreeNode* node;  // &pNodes[ele]
    unsigned int ele;
    unsigned int blk;
    unsigned int level;
    bool isBdyOct;
    Point domain_pt;  // domain coord of the element's min corner
    double hx[3];     // domain grid spacing along each axis
    // thread-local wavelet machinery (owned by the driver). im1/im2 are the
    // per-thread scratch the thread-safe compute_wavelets_3D needs so threads
    // don't race on RefElement's shared im_vec1/im_vec2.
    wavelet::WaveletEl* wrefEl;
    double* eVecTmp;
    std::vector<double>* wCout;
    const unsigned int* isz;
    double* im1;
    double* im2;

    // Relative L2 wavelet-coefficient error of one unzipped variable at this
    // element (matches the original WAMR normalization).
    double waveletErrorL2(const double* unzipVar, bool relative = false) const {
        mesh->getUnzipElementalNodalValues(unzipVar, blk, ele, eVecTmp, true);
        wrefEl->compute_wavelets_3D((double*)eVecTmp, isz, *wCout, isBdyOct, im1,
                                    im2);
        // relative WAMR normalizes by the element's L-inf; absolute uses 1.0.
        double Linf = 1.0;
        if (relative)
            Linf = std::max(normLInfty(eVecTmp, isz[0] * isz[1] * isz[2]), 1e-2);
        return (normL2(wCout->data(), wCout->size())) / sqrt(wCout->size()) /
               Linf;
    }
};

template <typename DecideFn>
static void refineFlagsPass(ot::Mesh* pMesh,
                            std::vector<unsigned int>& refine_flags,
                            DecideFn&& decide) {
    if (!pMesh->isActive()) return;
    const unsigned int eleLocalBegin      = pMesh->getElementLocalBegin();
    const ot::TreeNode* pNodes            = pMesh->getAllElements().data();
    const RefElement* refEl               = pMesh->getReferenceElement();
    const std::vector<ot::Block>& blkList = pMesh->getLocalBlockList();
    const unsigned int eOrder             = pMesh->getElementOrder();
    const unsigned int n1                 = 2 * eOrder + 1;
    const unsigned int szpd               = n1 * n1 * n1;
    const unsigned int nPe = (eOrder + 1) * (eOrder + 1) * (eOrder + 1);
    const unsigned int isz[3]             = {n1, n1, n1};

    auto run_block = [&](size_t b, wavelet::WaveletEl& wrefEl, double* eVecTmp,
                         std::vector<double>& wCout, double* im1, double* im2) {
        for (unsigned int ele = blkList[b].getLocalElementBegin();
             ele < blkList[b].getLocalElementEnd(); ele++) {
            const double oct_dx =
                (1u << (m_uiMaxDepth - pNodes[ele].getLevel())) /
                (double(eOrder));
            Point oct_pt1(pNodes[ele].minX(), pNodes[ele].minY(),
                          pNodes[ele].minZ());
            Point oct_pt2(pNodes[ele].minX() + oct_dx,
                          pNodes[ele].minY() + oct_dx,
                          pNodes[ele].minZ() + oct_dx);
            Point dpt1, dpt2;
            pMesh->octCoordToDomainCoord(oct_pt1, dpt1);
            pMesh->octCoordToDomainCoord(oct_pt2, dpt2);
            Point dxd = dpt2 - dpt1;

            RefineElemCtx ctx;
            ctx.mesh      = pMesh;
            ctx.node      = &pNodes[ele];
            ctx.ele       = ele;
            ctx.blk       = (unsigned int)b;
            ctx.level     = pNodes[ele].getLevel();
            ctx.isBdyOct  = pMesh->isBoundaryOctant(ele);
            ctx.domain_pt = dpt1;
            ctx.hx[0]     = dxd.x();
            ctx.hx[1]     = dxd.y();
            ctx.hx[2]     = dxd.z();
            ctx.wrefEl    = &wrefEl;
            ctx.eVecTmp   = eVecTmp;
            ctx.wCout     = &wCout;
            ctx.isz       = isz;
            ctx.im1       = im1;
            ctx.im2       = im2;
            refine_flags[ele - eleLocalBegin] =
                decide(ctx, refine_flags[ele - eleLocalBegin]);
        }
    };

#ifdef DENDRO_HYBRID_OMP
#pragma omp parallel
    {
        wavelet::WaveletEl wrefEl_t((RefElement*)refEl);
        std::vector<double> eVecTmp_t(szpd), wCout_t(szpd);
        std::vector<double> im1_t(nPe), im2_t(nPe);
#pragma omp for schedule(dynamic, 1)
        for (size_t b = 0; b < blkList.size(); b++)
            run_block(b, wrefEl_t, eVecTmp_t.data(), wCout_t, im1_t.data(),
                      im2_t.data());
    }
#else
    wavelet::WaveletEl wrefEl_t((RefElement*)refEl);
    std::vector<double> eVecTmp_t(szpd), wCout_t(szpd);
    std::vector<double> im1_t(nPe), im2_t(nPe);
    for (size_t b = 0; b < blkList.size(); b++)
        run_block(b, wrefEl_t, eVecTmp_t.data(), wCout_t, im1_t.data(),
                  im2_t.data());
#endif
}

// setMeshRefinementFlags on active ranks + the (collective) is-oct-change
// allreduce on all ranks. Pairs with refineFlagsPass.
static bool commitRefineFlags(ot::Mesh* pMesh,
                              std::vector<unsigned int>& refine_flags) {
    bool isOctChange = false;
    if (pMesh->isActive())
        isOctChange = pMesh->setMeshRefinementFlags(refine_flags);
    bool isOctChange_g = false;
    MPI_Allreduce(&isOctChange, &isOctChange_g, 1, MPI_CXX_BOOL, MPI_LOR,
                  pMesh->getMPIGlobalCommunicator());
    return isOctChange_g;
}

bool isRemeshBH(ot::Mesh* pMesh, const Point* bhLoc,
                const dendro_bh::BHHistory& bhHistory,
                const dendro_aeh::AEH_BHaHAHA* ahFinder) {
    ////////////////////////////////////////////////////////////////////
    // set up booleans to check whether the grid has changed
    bool isOctChange   = false;
    bool isOctChange_g = false;

    // wired through for future use (see HOOK below); unused by current policy
    (void)ahFinder;

    if (pMesh->isActive()) {
        ////////////////////////////////////////////////////////////////
        // Read in data from mesh
        const unsigned int eleLocalBegin = pMesh->getElementLocalBegin();
        const unsigned int eleLocalEnd   = pMesh->getElementLocalEnd();
        Point d1, d2, temp;

        // initialize refinement flags
        std::vector<unsigned int> refine_flags;
        // NOTE: these are the previous refine flags,
        // they're usually "no change" on remesh
        std::vector<unsigned int> prev_refine_flags =
            pMesh->getAllRefinementFlags();

        ////////////////////////////////////////////////////////////////
        // Set up onion parameters
        // read in AMR radii
        const double r_near[2] = {bssn::BSSN_BH1_AMR_R, bssn::BSSN_BH2_AMR_R};
        // level offset immediately about the BHs
        // (we don't actually refine to BSSN_BH?_MAX_LEV; two less)
        const unsigned int LVL_OFF     = MAXDEAPTH_LEVEL_DIFF + 1;

        // consider BHs merged if punctures are less than this value
        const double BH_MERGED_SEP_TOL = 0.1;
        // distance btw the black holes
        const double dBH               = (bhLoc[0] - bhLoc[1]).abs();

        ////////////////////////////////////////////////////////////////
        // Inert-BH robustness. A (near-)massless puncture contributes
        // nothing to the spacetime, so we must never build refinement
        // around it -- doing so wastes grid on a phantom center and can
        // seed instability (e.g. a single BH run with BH2 MASS=0 parked
        // on-grid). A BH is "active" iff its mass exceeds this tolerance.
        // For a real binary both are active, so every guard below is a
        // no-op and this path stays bit-identical to the two-BH code.
        constexpr double BH_ACTIVE_MASS_TOL = 1.0e-6;
        const bool bh_active[2] = {bssn::BSSN_BH1_MASS > BH_ACTIVE_MASS_TOL,
                                   bssn::BSSN_BH2_MASS > BH_ACTIVE_MASS_TOL};
        const int n_active_bh =
            (bh_active[0] ? 1 : 0) + (bh_active[1] ? 1 : 0);
        // Orbital scale is meaningless with fewer than two active BHs;
        // don't let an inert BH's (phantom) separation drive R_orbit.
        const double dBH_eff = (n_active_bh < 2) ? 0.0 : dBH;

        // lower of two max depths
        const unsigned int refLevMin =
            std::min(bssn::BSSN_BH1_MAX_LEV, bssn::BSSN_BH2_MAX_LEV);

        ////////////////////////////////////////////////////////////////
        // Set up a few ORBIT parameters
        // grid width
        const double Delta_x  = bssn::BSSN_GRID_MAX_X - bssn::BSSN_GRID_MIN_X;
        // smallest GW extraction radius
        const double R_GW_min = GW::BSSN_GW_RADAII[0];
        // largest GW extraction radius
        const double R_GW_max = GW::BSSN_GW_RADAII[GW::BSSN_GW_NUM_RADAII - 1];
        // element order
        const unsigned int n_order = bssn::BSSN_ELE_ORDER;
        // hardcode ending ORBIT some time past merger
        constexpr double t_ring    = 100;  // ringdown time
        // ORBIT disable time
        const double bh_merge_time = bssn::BSSN_BH_MERGE_TIME;
        const double t_disable     = bh_merge_time + R_GW_max + t_ring;

        ////////////////////////////////////////////////////////////////
        // if(!pMesh->getMPIRank())
        //     std::cout<<"bh distance: "<<dBH<<std::endl;

        const ot::TreeNode* pNodes = pMesh->getAllElements().data();
        refine_flags.resize(pMesh->getNumLocalMeshElements(), OCT_NO_CHANGE);

        // kinematics come from the incremental BHHistory; only the clean
        // orbital wavelength is consumed below

        // refine pass: per-element + write-disjoint, so thread directly (gated);
        // d1/d2/temp are per-thread. Geometry-only (no wavelet scratch needed),
        // so a direct pragma is cleaner than the driver here. Bit-identical to
        // the serial path.
#ifdef DENDRO_HYBRID_OMP
#pragma omp parallel for schedule(dynamic, 1) private(d1, d2, temp)
#endif
        for (unsigned int ele = eleLocalBegin; ele < eleLocalEnd; ele++) {
            // calculate which region of the grid we're in
            const unsigned int ln = 1u
                                    << (m_uiMaxDepth - pNodes[ele].getLevel());

            // Set obnoxiously large values for minimum radii
            // (to be overwritten) so as we go through box edges
            // we find minimum distances to BHs and to grid center
            double r1_min = 100000 * Delta_x;  // distance to BH1
            double r2_min = 100000 * Delta_x;  // distance to BH2
            double r_min  = 100000 * Delta_x;  // distance to grid center

            // measure minimum distances between both BHs
            for (unsigned int kk = 0; kk < 2; kk++)
                for (unsigned int jj = 0; jj < 2; jj++)
                    for (unsigned int ii = 0; ii < 2; ii++) {
                        const double x      = pNodes[ele].minX() + ii * ln;
                        const double y      = pNodes[ele].minY() + jj * ln;
                        const double z      = pNodes[ele].minZ() + kk * ln;
                        const Point oct_mid = Point(x, y, z);
                        pMesh->octCoordToDomainCoord(oct_mid, temp);

                        // vectors pointing toward each BH
                        d1     = temp - bhLoc[0];
                        d2     = temp - bhLoc[1];
                        // update minimum distance from each BH
                        r1_min = std::min(r1_min, d1.abs());
                        r2_min = std::min(r2_min, d2.abs());
                        // update minimum distance from grid center
                        r_min  = std::min(r_min, temp.abs());
                    }

            ////////////////////////////////////////////////////////////
            // wkb 5 Sept 2024: make this into nice functions

            // set default of coarsening; use BHLB before any others!
            refine_flags[ele - eleLocalBegin] = OCT_COARSE;

            // function which ensures we're at least at a given level
            auto setLevelFloor                = [&](int l_min) {
                int currentLevel = pNodes[ele].getLevel();

                if (currentLevel < l_min) {
                    // if below desired refinement level, split
                    refine_flags[ele - eleLocalBegin] = OCT_SPLIT;
                } else if (currentLevel == l_min &&
                           refine_flags[ele - eleLocalBegin] == OCT_COARSE) {
                    // if at desired level, prevent coarsening
                    refine_flags[ele - eleLocalBegin] = OCT_NO_CHANGE;
                }
                // else: sufficiently refined
            };

            ////////////////////////////////////////////////////////////
            // wkb 29 Aug 2024: add refinement based on radius
            // if radius r <= r*, keep level l >= l*
            // 9/9/24: change to depend on current BH separation dist.

            // set up orbital scale
            const double m1      = bssn::BSSN_BH1_MASS;  // mass of BH1
            const double m2      = bssn::BSSN_BH2_MASS;  // mass of BH2
            // calculate relative distance to each BH
            const double f1      = m2 / (m1 + m2);
            const double f2      = m1 / (m1 + m2);
            const double f       = std::max(f1, f2);
            const double R_orbit = f * dBH_eff + 8;  // M; resolve scale
            const int l_orbit    = 2;  // desired refinement level within
            if (r_min <= R_orbit) {
                // set up orbital radius scale
                setLevelFloor(l_orbit);
            }

            ////////////////////////////////////////////////////////////
            // wkb 2 Dec 2024: Onion refinement about the BHs
            auto onionLevel = [R_orbit, l_orbit](double radius, double r_AMR,
                                                 int maxLevel,
                                                 double ratio = 2.0) -> int {
                if (radius > R_orbit) {
                    // don't enforce onion outside orbital radius
                    return 0;
                } else {
                    // Start with the AMR radius and maximum level
                    double currentRadius = r_AMR;
                    int currentLevel     = maxLevel;
                    // Loop until we reach or go below the orbit level
                    while (currentLevel > l_orbit) {
                        // If input radius w/i current radius,
                        // return current level requirement
                        if (radius <= currentRadius) {
                            return currentLevel;
                        }
                        // Otherwise increase the radius limit
                        currentRadius *= ratio;
                        // and decrement the refinement level
                        currentLevel--;
                    }
                    // Shouldn't be here.
                    return -1;
                }
            };

            ////////////////////////////////////////////////////////////
            // Onion refinement immediately about the black holes
            if (dBH > BH_MERGED_SEP_TOL) {
                // if not merged yet, handle BHs separately to set up onion.
                // smoothly transition each BH's onion ratio toward the target
                // as the gauge wave passes (logistic in retarded time).
                const double t_current = bssn::BSSN_CURRENT_RK_COORD_TIME;
                const double ratio_1 = onionRatioGaugeLogistic(t_current, r1_min);
                const double ratio_2 = onionRatioGaugeLogistic(t_current, r2_min);
                // refinement levels near each BH
                const int l_goal_0 =
                    onionLevel(r1_min, r_near[0],
                               bssn::BSSN_BH1_MAX_LEV - LVL_OFF, ratio_1);
                const int l_goal_1 =
                    onionLevel(r2_min, r_near[1],
                               bssn::BSSN_BH2_MAX_LEV - LVL_OFF, ratio_2);
                // only refine around active BHs; an inert (massless)
                // puncture must not spawn its own onion (bit-exact for
                // real binaries, where both are active).
                if (bh_active[0]) setLevelFloor(l_goal_0);
                if (bh_active[1]) setLevelFloor(l_goal_1);
            } else {
                // if merged, handle BHs together
                // calculate minimum distance to either BH
                const double rBH_min = std::min(r1_min, r2_min);
                // calculate outer radius to which we should refine
                // ensure it captures both BHs and is >= than before
                //
                // HOOK for redesigning this target with the measured radius:
                //   if (ahFinder && ahFinder->common_horizon_found())
                //       r_meas = ahFinder->get_common_horizon_mean_radius();
                // left as the mass-based proxy for now (behavior unchanged).
                const double rBH_lim =
                    std::max(std::max(r_near[0], r_near[1]), 1.55 * (m1 + m2));
                // calculate level floor due to onion structure
                const int N_horizon =
                    50;  // goal number of points across the horizon
                const int l_post = std::ceil(
                    2 + std::log2(Delta_x * std::floor((N_horizon + 1) / 2.) /
                                  rBH_lim / n_order));
                // could here do l_post = std::min(l_post,bssn::BSSN_MAXDEPTH)
                // for low-res runs to keep l_post <= max depth
                const int l_goal =
                    onionLevel(rBH_min, rBH_lim, l_post - LVL_OFF);
                // const int l_goal = onionLevel(rBH_min,rBH_lim,refLevMin -
                // LVL_OFF); set level floor
                setLevelFloor(l_goal);
            }

#ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES
            ////////////////////////////////////////////////////////////
            // ORBIT refinement
            // @wkb 6 Feb 2025

            // calculate retarded time
            const double t_ret = bssn::BSSN_CURRENT_RK_COORD_TIME - r_min;
            // calculate retarded orbital wavelength
            const double lam_clean = bhHistory.clean_wavelength_at(t_ret);

            // calculate wavelength goal based on goal sph. har. order m
            auto get_ell = [lam_clean, Delta_x, n_order](
                               double t_ret, unsigned int m = 8) -> int {
                // Get refinement level necessary from wavelength
                return static_cast<int>(std::ceil(
                    std::log2(Delta_x * m / (n_order * lam_clean / 2.))));
            };

            // min refinement level required from GWs
            int ell_star;

            /*
            // enable this code here to use a different Nyquist criteria
            // beyond the GW extraction radii (to cut down costs)
            if (r_min <= R_GW_max) {
              // if w/i region of capturing GWs, resolve highly
              ell_star = get_ell(t_ret, bssn::BSSN_NYQUIST_M);
            } else {
              // otherwise, only prevent major backreflections
              ell_star = get_ell(t_ret, bssn::BSSN_NYQUIST_M - 2);
            }
            */

            // goal spherical harmonic order m to refine to
            const bool using_nyquist = bssn::BSSN_NYQUIST_M > 0;
            const double speed       = std::sqrt(2);
            const bool past_of_end =
                std::abs(r_min - R_GW_max) <
                speed * (t_disable - bssn::BSSN_CURRENT_RK_COORD_TIME);
            if (using_nyquist && past_of_end) {
                // same ell_star everywhere
                ell_star = get_ell(t_ret, bssn::BSSN_NYQUIST_M);
                setLevelFloor(ell_star);
            }
#endif

            ////////////////////////////////////////////////////////////
            // @wkb 12 June 2024: Don't allow min depth violation
            setLevelFloor(bssn::BSSN_MINDEPTH);
        }

        isOctChange = pMesh->setMeshRefinementFlags(refine_flags);
    }

    // communicate changes to the grid
    bool isOctChanged_g;
    MPI_Allreduce(&isOctChange, &isOctChanged_g, 1, MPI_CXX_BOOL, MPI_LOR,
                  pMesh->getMPIGlobalCommunicator());
    return isOctChanged_g;
}

bool isRemeshEH(ot::Mesh* pMesh, const double** unzipVec, unsigned int vIndex,
                double refine_th, double coarsen_th, bool isOverwrite) {
    const unsigned int eleLocalBegin = pMesh->getElementLocalBegin();
    const unsigned int eleLocalEnd   = pMesh->getElementLocalEnd();
    bool isOctChange                 = false;
    bool isOctChange_g               = false;
    const unsigned int eOrder        = pMesh->getElementOrder();

    std::vector<unsigned int> refine_flags;

    // NOTE: this is what the flags are set to, note that
    const std::vector<unsigned int> prev_refine_flags =
        pMesh->getAllRefinementFlags();

    if (pMesh->isActive()) {
        ot::TreeNode* pNodes =
            (ot::TreeNode*)&(*(pMesh->getAllElements().begin()));

        if (isOverwrite)
            for (unsigned int ele = eleLocalBegin; ele < eleLocalEnd; ele++)
                pNodes[ele].setFlag(((OCT_NO_CHANGE << NUM_LEVEL_BITS) |
                                     pNodes[ele].getLevel()));

        const std::vector<ot::Block>& blkList = pMesh->getLocalBlockList();
        unsigned int sz[3];
        unsigned int ei[3];

        // refine test: per-element, write-disjoint -> thread the block loop
        // (sz/ei per-thread). Gated; serial + bit-identical when off.
#ifdef DENDRO_HYBRID_OMP
#pragma omp parallel for schedule(dynamic, 1) private(sz, ei)
#endif
        for (unsigned int b = 0; b < blkList.size(); b++) {
            const ot::TreeNode blkNode = blkList[b].getBlockNode();

            sz[0]                      = blkList[b].getAllocationSzX();
            sz[1]                      = blkList[b].getAllocationSzY();
            sz[2]                      = blkList[b].getAllocationSzZ();

            const unsigned int bflag   = blkList[b].getBlkNodeFlag();
            const unsigned int offset  = blkList[b].getOffset();

            const unsigned int regLev  = blkList[b].getRegularGridLev();
            const unsigned int eleIndexMax =
                (1u << (regLev - blkNode.getLevel())) - 1;
            const unsigned int eleIndexMin = 0;

            for (unsigned int ele = blkList[b].getLocalElementBegin();
                 ele < blkList[b].getLocalElementEnd(); ele++) {
                ei[0] = (pNodes[ele].getX() - blkNode.getX()) >>
                        (m_uiMaxDepth - regLev);
                ei[1] = (pNodes[ele].getY() - blkNode.getY()) >>
                        (m_uiMaxDepth - regLev);
                ei[2] = (pNodes[ele].getZ() - blkNode.getZ()) >>
                        (m_uiMaxDepth - regLev);

                if ((bflag & (1u << OCT_DIR_LEFT)) && ei[0] == eleIndexMin)
                    continue;
                if ((bflag & (1u << OCT_DIR_DOWN)) && ei[1] == eleIndexMin)
                    continue;
                if ((bflag & (1u << OCT_DIR_BACK)) && ei[2] == eleIndexMin)
                    continue;

                if ((bflag & (1u << OCT_DIR_RIGHT)) && ei[0] == eleIndexMax)
                    continue;
                if ((bflag & (1u << OCT_DIR_UP)) && ei[1] == eleIndexMax)
                    continue;
                if ((bflag & (1u << OCT_DIR_FRONT)) && ei[2] == eleIndexMax)
                    continue;

                // refine test.
                for (unsigned int k = 3; k < eOrder + 1 + 3; k++)
                    for (unsigned int j = 3; j < eOrder + 1 + 3; j++)
                        for (unsigned int i = 3; i < eOrder + 1 + 3; i++) {
                            if (unzipVec[vIndex]
                                        [offset +
                                         (ei[2] * eOrder + k) * sz[0] * sz[1] +
                                         (ei[1] * eOrder + j) * sz[0] +
                                         (ei[0] * eOrder + i)] < refine_th) {
                                if ((pNodes[ele].getLevel() +
                                     MAXDEAPTH_LEVEL_DIFF + 1) < m_uiMaxDepth)
                                    pNodes[ele].setFlag(
                                        ((OCT_SPLIT << NUM_LEVEL_BITS) |
                                         pNodes[ele].getLevel()));
                            }
                        }
            }
        }

        // coarsen test. Sibling groups stay within a block, so blocks are
        // independent -> thread the block loop (sz/ei per-thread). The refine
        // pass above completed first (separate loop). Gated; serial when off.
#ifdef DENDRO_HYBRID_OMP
#pragma omp parallel for schedule(dynamic, 1) private(sz, ei)
#endif
        for (unsigned int b = 0; b < blkList.size(); b++) {
            const ot::TreeNode blkNode = blkList[b].getBlockNode();

            sz[0]                      = blkList[b].getAllocationSzX();
            sz[1]                      = blkList[b].getAllocationSzY();
            sz[2]                      = blkList[b].getAllocationSzZ();

            const unsigned int bflag   = blkList[b].getBlkNodeFlag();
            const unsigned int offset  = blkList[b].getOffset();

            const unsigned int regLev  = blkList[b].getRegularGridLev();
            const unsigned int eleIndexMax =
                (1u << (regLev - blkNode.getLevel())) - 1;
            const unsigned int eleIndexMin = 0;

            if ((eleIndexMax == 0) || (bflag != 0))
                continue;  // this implies the blocks with only 1 child and
                           // boundary blocks.

            for (unsigned int ele = blkList[b].getLocalElementBegin();
                 ele < blkList[b].getLocalElementEnd(); ele++) {
                assert(pNodes[ele].getParent() ==
                       pNodes[ele + NUM_CHILDREN - 1].getParent());
                bool isCoarsen = true;

                for (unsigned int child = 0; child < NUM_CHILDREN; child++) {
                    if ((pNodes[ele + child].getFlag() >> NUM_LEVEL_BITS) ==
                        OCT_SPLIT) {
                        isCoarsen = false;
                        break;
                    }
                }

                if (isCoarsen && pNodes[ele].getLevel() > 1) {
                    bool coarse = true;
                    for (unsigned int child = 0; child < NUM_CHILDREN;
                         child++) {
                        ei[0] = (pNodes[ele + child].getX() - blkNode.getX()) >>
                                (m_uiMaxDepth - regLev);
                        ei[1] = (pNodes[ele + child].getY() - blkNode.getY()) >>
                                (m_uiMaxDepth - regLev);
                        ei[2] = (pNodes[ele + child].getZ() - blkNode.getZ()) >>
                                (m_uiMaxDepth - regLev);

                        for (unsigned int k = 3; k < eOrder + 1 + 3; k++)
                            for (unsigned int j = 3; j < eOrder + 1 + 3; j++)
                                for (unsigned int i = 3; i < eOrder + +3; i++) {
                                    if (!((refine_th <
                                           unzipVec[vIndex]
                                                   [offset +
                                                    (ei[2] * eOrder + k) *
                                                        sz[0] * sz[1] +
                                                    (ei[1] * eOrder + j) *
                                                        sz[0] +
                                                    (ei[0] * eOrder + i)]) &&
                                          (unzipVec[vIndex]
                                                   [offset +
                                                    (ei[2] * eOrder + k) *
                                                        sz[0] * sz[1] +
                                                    (ei[1] * eOrder + j) *
                                                        sz[0] +
                                                    (ei[0] * eOrder + i)] <=
                                           coarsen_th)))
                                        coarse = false;
                                }
                    }

                    if (coarse)
                        for (unsigned int child = 0; child < NUM_CHILDREN;
                             child++)
                            pNodes[ele + child].setFlag(
                                ((OCT_COARSE << NUM_LEVEL_BITS) |
                                 pNodes[ele].getLevel()));
                }

                ele = ele + NUM_CHILDREN - 1;
            }
        }

        for (unsigned int ele = eleLocalBegin; ele < eleLocalEnd; ele++)
            if ((pNodes[ele].getFlag() >> NUM_LEVEL_BITS) ==
                OCT_SPLIT)  // trigger remesh only when some refinement occurs
                            // (laid back remesh :)  )
            {
                isOctChange = true;
                break;
            }
    }

    bool isOctChanged_g;
    MPI_Allreduce(&isOctChange, &isOctChanged_g, 1, MPI_CXX_BOOL, MPI_LOR,
                  pMesh->getMPIGlobalCommunicator());
    // if(!m_uiGlobalRank) std::cout<<"is oct changed:
    // "<<isOctChanged_g<<std::endl;
    return isOctChanged_g;
}

bool isReMeshWAMR(
    ot::Mesh* pMesh, const double** unzippedVec, const unsigned int* varIds,
    const unsigned int numVars,
    std::function<double(double, double, double, double*)> wavelet_tol,
    double amr_coarse_fac) {
    // if(!(pMesh->isReMeshUnzip((const double
    // **)unzippedVec,varIds,numVars,wavelet_tol,bssn::BSSN_DENDRO_AMR_FAC)))
    //     return false;

    // if(bssn::BSSN_CURRENT_RK_COORD_TIME > 0 &&
    // bssn::BSSN_CURRENT_RK_COORD_TIME < 80)
    //     return bssn::isReMeshBHRadial(pMesh);

    std::vector<unsigned int> refine_flags;
    const double r_near[2] = {bssn::BSSN_BH1_AMR_R, bssn::BSSN_BH2_AMR_R};

    const unsigned int eleLocalBegin = pMesh->getElementLocalBegin();
    const unsigned int eleLocalEnd   = pMesh->getElementLocalEnd();
    bool isOctChange                 = false;
    bool isOctChange_g               = false;
    Point d1, d2, temp;

    const unsigned int eOrder = pMesh->getElementOrder();
    const double dBH          = (BSSN_BH_LOC[0] - BSSN_BH_LOC[1]).abs();
    const unsigned int refLevMin =
        std::min(bssn::BSSN_BH1_MAX_LEV, bssn::BSSN_BH2_MAX_LEV);

    // BH considered merged if the distance between punctures are less than the
    // specified value.
    const double BH_MERGED_SEP_TOL = 0.1;

    if (pMesh->isActive()) {
        if (!pMesh->getMPIRank()) printf("BH coord sep: %.8E \n", dBH);
        // std::cout<<"BH coord sep: "<<dBH<<std::endl;

        refine_flags.resize(pMesh->getNumLocalMeshElements(), OCT_NO_CHANGE);

        // Phase 1: WAMR wavelet decision. Per-element + parallel via the unified
        // driver; the physicist edits only this lambda. (Math identical to the
        // former serial loop: max relative wavelet-coeff L2 over refine vars,
        // with the early-bail, vs the per-element tolerance.)
        refineFlagsPass(
            pMesh, refine_flags,
            [&](const RefineElemCtx& c, unsigned int) -> unsigned int {
                const double tol_ele =
                    wavelet_tol(c.domain_pt.x(), c.domain_pt.y(),
                                c.domain_pt.z(), const_cast<double*>(c.hx));
                double l_max = 0.0;
                for (unsigned int v = 0; v < numVars; v++) {
                    const double w = c.waveletErrorL2(unzippedVec[varIds[v]]);
                    l_max          = std::max(l_max, w);
                    if (w > tol_ele) break;  // early bail
                }
                if (l_max > tol_ele)
                    return (unsigned int)OCT_SPLIT;
                else if (l_max < amr_coarse_fac * tol_ele)
                    return (unsigned int)OCT_COARSE;
                return (unsigned int)OCT_NO_CHANGE;
            });
        // end of WAMR core calculation.

        ////////////////////////////////////////////////////////////////
        // Below code enforces a certain level of refinement at the BHs,
        // overriding what's currently set by the wavelets. Per-element and
        // write-disjoint, so thread it directly (gated); d1/d2/temp are
        // per-thread. Bit-identical to the serial path.
        const ot::TreeNode* pNodes = pMesh->getAllElements().data();
#ifdef DENDRO_HYBRID_OMP
#pragma omp parallel for schedule(dynamic, 1) private(d1, d2, temp)
#endif
        for (unsigned int ele = eleLocalBegin; ele < eleLocalEnd; ele++) {
            // refine_flags[ele-eleLocalBegin] =
            // (pNodes[ele].getFlag()>>NUM_LEVEL_BITS); std::cout<<"ref flag:
            // "<<(pNodes[ele].getFlag()>>NUM_LEVEL_BITS)<<std::endl;
            // if(refine_flags[ele-eleLocalBegin]==OCT_SPLIT)
            pMesh->octCoordToDomainCoord(
                Point((double)pNodes[ele].minX(), (double)pNodes[ele].minY(),
                      (double)pNodes[ele].minZ()),
                temp);
            d1 = temp - BSSN_BH_LOC[0];
            d2 = temp - BSSN_BH_LOC[1];

            //@milinda: 11/21/2020 : Don't allow to violate the min depth
            if (pNodes[ele].getLevel() < bssn::BSSN_MINDEPTH) {
                refine_flags[ele - eleLocalBegin] = OCT_SPLIT;
            } else if (pNodes[ele].getLevel() == bssn::BSSN_MINDEPTH &&
                       refine_flags[ele - eleLocalBegin] == OCT_COARSE) {
                refine_flags[ele - eleLocalBegin] = OCT_NO_CHANGE;
            }

            // don't overide things away from puntures let wavelets handle that.
            if (d1.abs() > 10 && d2.abs() > 10)
                continue;
            else {
                const unsigned int ln =
                    1u << (m_uiMaxDepth - pNodes[ele].getLevel());
                const double hx = ln / (double)(eOrder);
                for (unsigned int k = 0; k < (eOrder + 1); k++)
                    for (unsigned int j = 0; j < (eOrder + 1); j++)
                        for (unsigned int i = 0; i < (eOrder + 1); i++) {
                            const double x      = pNodes[ele].minX() + k * hx;
                            const double y      = pNodes[ele].minY() + j * hx;
                            const double z      = pNodes[ele].minZ() + i * hx;
                            const Point oct_mid = Point(x, y, z);

                            pMesh->octCoordToDomainCoord(oct_mid, temp);

                            d1                     = temp - BSSN_BH_LOC[0];
                            d2                     = temp - BSSN_BH_LOC[1];

                            // std::cout<<"d1: "<<d1 <<
                            // "BHLOC_0:"<<BSSN_BH_LOC[0]<<std::endl;
                            // std::cout<<"d2: "<<d2<<std::endl;

                            const double rd1       = d1.abs();
                            const double rd2       = d2.abs();

                            const bool isNearTobh1 = (rd1 <= r_near[0]);
                            const bool isNearTobh2 = (rd2 <= r_near[1]);

                            const bool isMidNearTobh1 =
                                (rd1 > r_near[0] && rd1 <= 10.0 * r_near[0]);
                            const bool isMidNearTobh2 =
                                (rd2 > r_near[1] && rd1 <= 10.0 * r_near[1]);

                            const bool isFarTobh1 = (rd1 > 2.0 * r_near[0]);
                            const bool isFarTobh2 = (rd2 > 2.0 * r_near[1]);

                            if (dBH < BH_MERGED_SEP_TOL) {
                                if (isNearTobh1 || isNearTobh2) {
                                    // std::cout<<"d1:
                                    // "<<d1.abs()<<"BHLOC_0:"<<BSSN_BH_LOC[0]<<std::endl;
                                    // std::cout<<"d2:
                                    // "<<d2.abs()<<"BHLOC_1:"<<BSSN_BH_LOC[1]<<std::endl;

                                    if ((pNodes[ele].getLevel() +
                                         MAXDEAPTH_LEVEL_DIFF + 1) < refLevMin)
                                        refine_flags[ele - eleLocalBegin] =
                                            OCT_SPLIT;
                                    else if ((pNodes[ele].getLevel() +
                                              MAXDEAPTH_LEVEL_DIFF + 1) >
                                             refLevMin)
                                        refine_flags[ele - eleLocalBegin] =
                                            OCT_COARSE;
                                    else
                                        refine_flags[ele - eleLocalBegin] =
                                            OCT_NO_CHANGE;

                                    // if( ( pNodes[ele].getLevel() +
                                    // MAXDEAPTH_LEVEL_DIFF +1)== refLevMin )
                                    //     refine_flags[ele-eleLocalBegin] =
                                    //     OCT_NO_CHANGE;
                                    // else if(( pNodes[ele].getLevel() +
                                    // MAXDEAPTH_LEVEL_DIFF +1)> refLevMin)
                                    //     refine_flags[ele-eleLocalBegin] =
                                    //     OCT_COARSE;
                                }

                            } else {
                                if (bssn::BSSN_BH1_MAX_LEV == refLevMin) {
                                    if (isNearTobh1) {
                                        // std::cout<<"d1:
                                        // "<<d1.abs()<<"BHLOC_0:"<<BSSN_BH_LOC[0]<<"
                                        // rnear: "<<r_near[0]<<std::endl;
                                        if ((pNodes[ele].getLevel() +
                                             MAXDEAPTH_LEVEL_DIFF + 1) <
                                            bssn::BSSN_BH1_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_SPLIT;
                                        else if ((pNodes[ele].getLevel() +
                                                  MAXDEAPTH_LEVEL_DIFF + 1) >
                                                 bssn::BSSN_BH1_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_COARSE;
                                        else
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_NO_CHANGE;

                                    } else {
                                        if (refine_flags[ele - eleLocalBegin] ==
                                                OCT_SPLIT &&
                                            (pNodes[ele].getLevel() +
                                             MAXDEAPTH_LEVEL_DIFF + 1) ==
                                                bssn::BSSN_BH1_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_NO_CHANGE;
                                        else if ((pNodes[ele].getLevel() +
                                                  MAXDEAPTH_LEVEL_DIFF + 1) >
                                                 bssn::BSSN_BH1_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_COARSE;
                                    }

                                    // changes in bh 1 will get overidden by lev
                                    // 2
                                    if (isNearTobh2) {
                                        if ((pNodes[ele].getLevel() +
                                             MAXDEAPTH_LEVEL_DIFF + 1) <
                                            bssn::BSSN_BH2_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_SPLIT;
                                        else if ((pNodes[ele].getLevel() +
                                                  MAXDEAPTH_LEVEL_DIFF + 1) >
                                                 bssn::BSSN_BH2_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_COARSE;
                                        else
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_NO_CHANGE;
                                    }

                                } else {
                                    assert(bssn::BSSN_BH2_MAX_LEV == refLevMin);
                                    if (isNearTobh2) {
                                        // std::cout<<"d1:
                                        // "<<d1.abs()<<"BHLOC_0:"<<BSSN_BH_LOC[0]<<"
                                        // rnear: "<<r_near[0]<<std::endl;
                                        if ((pNodes[ele].getLevel() +
                                             MAXDEAPTH_LEVEL_DIFF + 1) <
                                            bssn::BSSN_BH2_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_SPLIT;
                                        else if ((pNodes[ele].getLevel() +
                                                  MAXDEAPTH_LEVEL_DIFF + 1) >
                                                 bssn::BSSN_BH2_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_COARSE;
                                        else
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_NO_CHANGE;

                                    } else {
                                        if (refine_flags[ele - eleLocalBegin] ==
                                                OCT_SPLIT &&
                                            (pNodes[ele].getLevel() +
                                             MAXDEAPTH_LEVEL_DIFF + 1) ==
                                                bssn::BSSN_BH2_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_NO_CHANGE;
                                        else if ((pNodes[ele].getLevel() +
                                                  MAXDEAPTH_LEVEL_DIFF + 1) >
                                                 bssn::BSSN_BH2_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_COARSE;
                                    }

                                    // changes in bh 2 will get overidden by lev
                                    // 1 which is the higher level than bh2.
                                    if (isNearTobh1) {
                                        if ((pNodes[ele].getLevel() +
                                             MAXDEAPTH_LEVEL_DIFF + 1) <
                                            bssn::BSSN_BH1_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_SPLIT;
                                        else if ((pNodes[ele].getLevel() +
                                                  MAXDEAPTH_LEVEL_DIFF + 1) >
                                                 bssn::BSSN_BH1_MAX_LEV)
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_COARSE;
                                        else
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_NO_CHANGE;
                                    }
                                }
                            }
                        }
            }
        }

        isOctChange = pMesh->setMeshRefinementFlags(refine_flags);
    }

    // communicate refinement between cores.
    MPI_Allreduce(&isOctChange, &isOctChange_g, 1, MPI_CXX_BOOL, MPI_LOR,
                  pMesh->getMPIGlobalCommunicator());
    return isOctChange_g;
}

// Don't use this : 1) This is expensive
// 2). It has a coarsening bug (doesn't coarsen properly)
bool isReMeshBHRadial(ot::Mesh* pMesh) {
    std::vector<unsigned int> refine_flags;
    const double r_near[2] = {bssn::BSSN_BH1_AMR_R, bssn::BSSN_BH2_AMR_R};

    const unsigned int eleLocalBegin = pMesh->getElementLocalBegin();
    const unsigned int eleLocalEnd   = pMesh->getElementLocalEnd();
    bool isOctChange                 = false;
    bool isOctChange_g               = false;
    Point d1, d2, temp;

    const unsigned int eOrder = pMesh->getElementOrder();
    const double dBH          = (BSSN_BH_LOC[0] - BSSN_BH_LOC[1]).abs();
    const unsigned int refLevMin =
        std::min(bssn::BSSN_BH1_MAX_LEV, bssn::BSSN_BH2_MAX_LEV);

    // BH considered merged if the distance between punctures are less than the
    // specified value.
    const double BH_MERGED_SEP_TOL        = 0.1;
    const unsigned int NUM_REFINE_SPHERES = 5;

    std::vector<double> bh1_amr_r;
    std::vector<double> bh2_amr_r;

    bh1_amr_r.push_back(0.0);
    bh2_amr_r.push_back(0.0);

    for (unsigned int i = 0; i < NUM_REFINE_SPHERES; i++) {
        const unsigned int rfac = (1u << i);
        bh1_amr_r.push_back(rfac * bssn::BSSN_BH1_AMR_R);
        bh2_amr_r.push_back(rfac * bssn::BSSN_BH2_AMR_R);
    }

    if (pMesh->isActive()) {
        refine_flags.resize(pMesh->getNumLocalMeshElements(), OCT_NO_CHANGE);
        const ot::TreeNode* pNodes = pMesh->getAllElements().data();

        for (unsigned int ele = eleLocalBegin; ele < eleLocalEnd; ele++) {
            // refine_flags[ele-eleLocalBegin] =
            // (pNodes[ele].getFlag()>>NUM_LEVEL_BITS); std::cout<<"ref flag:
            // "<<(pNodes[ele].getFlag()>>NUM_LEVEL_BITS)<<std::endl;
            // if(refine_flags[ele-eleLocalBegin]==OCT_SPLIT)
            pMesh->octCoordToDomainCoord(
                Point((double)pNodes[ele].minX(), (double)pNodes[ele].minY(),
                      (double)pNodes[ele].minZ()),
                temp);
            d1 = temp - BSSN_BH_LOC[0];
            d2 = temp - BSSN_BH_LOC[1];

            //@milinda: 11/21/2020 : Don't allow to violate the min depth
            if (pNodes[ele].getLevel() < bssn::BSSN_MINDEPTH) {
                refine_flags[ele - eleLocalBegin] = OCT_SPLIT;
            } else if (pNodes[ele].getLevel() == bssn::BSSN_MINDEPTH &&
                       refine_flags[ele - eleLocalBegin] == OCT_COARSE) {
                refine_flags[ele - eleLocalBegin] = OCT_NO_CHANGE;
            }

            const unsigned int ln = 1u
                                    << (m_uiMaxDepth - pNodes[ele].getLevel());
            const double hx = ln / (double)(eOrder);

            for (unsigned int k = 0; k < (eOrder + 1); k++)
                for (unsigned int j = 0; j < (eOrder + 1); j++)
                    for (unsigned int i = 0; i < (eOrder + 1); i++) {
                        const double x      = pNodes[ele].minX() + i * hx;
                        const double y      = pNodes[ele].minY() + j * hx;
                        const double z      = pNodes[ele].minZ() + k * hx;
                        const Point oct_mid = Point(x, y, z);

                        pMesh->octCoordToDomainCoord(oct_mid, temp);

                        d1               = temp - BSSN_BH_LOC[0];
                        d2               = temp - BSSN_BH_LOC[1];

                        // std::cout<<"d1: "<<d1 <<
                        // "BHLOC_0:"<<BSSN_BH_LOC[0]<<std::endl;
                        // std::cout<<"d2: "<<d2<<std::endl;

                        const double rd1 = d1.abs();
                        const double rd2 = d2.abs();

                        if (dBH < BH_MERGED_SEP_TOL) {
                            for (unsigned int rs = 1;
                                 rs < NUM_REFINE_SPHERES + 1; rs++) {
                                if ((rd1 > bh1_amr_r[rs - 1]) &&
                                    (rd1 <= bh1_amr_r[rs])) {
                                    if ((pNodes[ele].getLevel() +
                                         MAXDEAPTH_LEVEL_DIFF + 1) <
                                        std::max(refLevMin - (rs - 1),
                                                 bssn::BSSN_MINDEPTH))
                                        refine_flags[ele - eleLocalBegin] =
                                            OCT_SPLIT;
                                    else if ((pNodes[ele].getLevel() +
                                              MAXDEAPTH_LEVEL_DIFF + 1) >
                                             std::max(refLevMin - (rs - 1),
                                                      bssn::BSSN_MINDEPTH))
                                        refine_flags[ele - eleLocalBegin] =
                                            OCT_COARSE;
                                    else
                                        refine_flags[ele - eleLocalBegin] =
                                            OCT_NO_CHANGE;
                                }
                                if (rd1 >
                                    bh1_amr_r.back())  // note : It is important
                                                       // to keep this inside
                                                       // the for loop to ensure
                                                       // proper refinement.
                                    refine_flags[ele - eleLocalBegin] =
                                        OCT_COARSE;
                            }

                        } else {
                            if (bssn::BSSN_BH1_MAX_LEV == refLevMin) {
                                for (unsigned int rs = 1;
                                     rs < NUM_REFINE_SPHERES + 1; rs++) {
                                    if ((rd1 > bh1_amr_r[rs - 1]) &&
                                        (rd1 <= bh1_amr_r[rs])) {
                                        if ((pNodes[ele].getLevel() +
                                             MAXDEAPTH_LEVEL_DIFF + 1) <
                                            std::max(bssn::BSSN_BH1_MAX_LEV -
                                                         (rs - 1),
                                                     bssn::BSSN_MINDEPTH))
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_SPLIT;
                                        else if ((pNodes[ele].getLevel() +
                                                  MAXDEAPTH_LEVEL_DIFF + 1) >
                                                 std::max(
                                                     bssn::BSSN_BH1_MAX_LEV -
                                                         (rs - 1),
                                                     bssn::BSSN_MINDEPTH))
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_COARSE;
                                        else
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_NO_CHANGE;

                                    } else if (rd1 > bh1_amr_r.back())
                                        refine_flags[ele - eleLocalBegin] =
                                            OCT_COARSE;

                                    // overide by smaller bh - BH2 is smaller
                                    // from the depth parameter.
                                    if ((rd2 > bh2_amr_r[rs - 1]) &&
                                        (rd2 <= bh2_amr_r[rs])) {
                                        if ((pNodes[ele].getLevel() +
                                             MAXDEAPTH_LEVEL_DIFF + 1) <
                                            std::max(bssn::BSSN_BH2_MAX_LEV -
                                                         (rs - 1),
                                                     bssn::BSSN_MINDEPTH))
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_SPLIT;
                                        else if ((pNodes[ele].getLevel() +
                                                  MAXDEAPTH_LEVEL_DIFF + 1) >
                                                 std::max(
                                                     bssn::BSSN_BH2_MAX_LEV -
                                                         (rs - 1),
                                                     bssn::BSSN_MINDEPTH))
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_COARSE;
                                        else
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_NO_CHANGE;
                                    } else if (rd2 > bh2_amr_r.back())
                                        refine_flags[ele - eleLocalBegin] =
                                            OCT_COARSE;
                                }

                            } else {
                                assert(bssn::BSSN_BH2_MAX_LEV == refLevMin);
                                for (unsigned int rs = 1;
                                     rs < NUM_REFINE_SPHERES + 1; rs++) {
                                    if ((rd2 > bh2_amr_r[rs - 1]) &&
                                        (rd2 <= bh2_amr_r[rs])) {
                                        if ((pNodes[ele].getLevel() +
                                             MAXDEAPTH_LEVEL_DIFF + 1) <
                                            std::max(bssn::BSSN_BH2_MAX_LEV -
                                                         (rs - 1),
                                                     bssn::BSSN_MINDEPTH))
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_SPLIT;
                                        else if ((pNodes[ele].getLevel() +
                                                  MAXDEAPTH_LEVEL_DIFF + 1) >
                                                 std::max(
                                                     bssn::BSSN_BH2_MAX_LEV -
                                                         (rs - 1),
                                                     bssn::BSSN_MINDEPTH))
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_COARSE;
                                        else
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_NO_CHANGE;
                                    } else if (rd2 > bh2_amr_r.back())
                                        refine_flags[ele - eleLocalBegin] =
                                            OCT_COARSE;

                                    // overide by smaller bh - BH1 is smaller
                                    // from the depth parameter.
                                    if ((rd1 > bh1_amr_r[rs - 1]) &&
                                        (rd1 <= bh1_amr_r[rs])) {
                                        if ((pNodes[ele].getLevel() +
                                             MAXDEAPTH_LEVEL_DIFF + 1) <
                                            std::max(bssn::BSSN_BH1_MAX_LEV -
                                                         (rs - 1),
                                                     bssn::BSSN_MINDEPTH))
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_SPLIT;
                                        else if ((pNodes[ele].getLevel() +
                                                  MAXDEAPTH_LEVEL_DIFF + 1) >
                                                 std::max(
                                                     bssn::BSSN_BH1_MAX_LEV -
                                                         (rs - 1),
                                                     bssn::BSSN_MINDEPTH))
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_COARSE;
                                        else
                                            refine_flags[ele - eleLocalBegin] =
                                                OCT_NO_CHANGE;
                                    } else if (rd1 > bh1_amr_r.back())
                                        refine_flags[ele - eleLocalBegin] =
                                            OCT_COARSE;
                                }
                            }
                        }
                    }
        }
        isOctChange = pMesh->setMeshRefinementFlags(refine_flags);
    }

    MPI_Allreduce(&isOctChange, &isOctChange_g, 1, MPI_CXX_BOOL, MPI_LOR,
                  pMesh->getMPIGlobalCommunicator());
    return isOctChange_g;
}

bool addRemeshWAMR(
    ot::Mesh* pMesh, const double** unzippedVec, const unsigned int* varIds,
    const unsigned int numVars,
    std::function<double(double, double, double, double*)> wavelet_tol,
    double amr_coarse_fac, bool relative_WAMR) {
    // WKB Aug 2024: add new refinement helper to add on refinement
    // in regions where wavelets say so, but don't use wavelets
    // directly for de-refinement. This adds up to this equalling
    // isReMeshWAMR, but without calling 'COARSEN' (nor trying to do
    // its own BHLB method

    // pull in current refinement flags
    std::vector<unsigned int> refine_flags = pMesh->getAllRefinementFlags();

    bool isOctChange   = false;
    bool isOctChange_g = false;
    const double dBH   = (BSSN_BH_LOC[0] - BSSN_BH_LOC[1]).abs();

    if (pMesh->isActive()) {
        if (!pMesh->getMPIRank()) printf("BH coord sep: %.8E \n", dBH);
        // std::cout<<"BH coord sep: "<<dBH<<std::endl;

        // WAMR add-only decision: per element + parallel via the unified
        // driver; physicist edits only this lambda.
        refineFlagsPass(
            pMesh, refine_flags,
            [&](const RefineElemCtx& c, unsigned int curFlag) -> unsigned int {
                double tol_ele =
                    wavelet_tol(c.domain_pt.x(), c.domain_pt.y(),
                                c.domain_pt.z(), const_cast<double*>(c.hx));
                // amplify the tolerance (effectively disable WAMR) within
                // either BH's AMR radius.
                const double r_BH1 =
                    (c.domain_pt - BSSN_BH_LOC[0]).abs() / bssn::BSSN_BH1_AMR_R;
                const double r_BH2 =
                    (c.domain_pt - BSSN_BH_LOC[1]).abs() / bssn::BSSN_BH2_AMR_R;
                if (std::min(r_BH1, r_BH2) <= 1) tol_ele *= 1e12;

                double l_max = 0.0;
                for (unsigned int v = 0; v < numVars; v++) {
                    const double w =
                        c.waveletErrorL2(unzippedVec[varIds[v]], relative_WAMR);
                    l_max = std::max(l_max, w);
                    if (w > tol_ele) break;  // early bail
                }
                // add refinement; protect needed refinement from coarsening;
                // otherwise leave the existing flag untouched.
                if (l_max > tol_ele)
                    return (unsigned int)OCT_SPLIT;
                else if (curFlag == OCT_COARSE &&
                         l_max >= amr_coarse_fac * tol_ele)
                    return (unsigned int)OCT_NO_CHANGE;
                return curFlag;
            });

        // check whether any changes have been made to the grid here
        isOctChange = pMesh->setMeshRefinementFlags(refine_flags);
    }

    // check whether changes have been made *anywhere* on the grid
    MPI_Allreduce(&isOctChange, &isOctChange_g, 1, MPI_CXX_BOOL, MPI_LOR,
                  pMesh->getMPIGlobalCommunicator());
    return isOctChange_g;
}

}  // end of namespace bssn
