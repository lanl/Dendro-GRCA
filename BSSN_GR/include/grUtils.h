//
// Created by milinda on 7/25/17.
/**
 *@author Milinda Fernando
 *School of Computing, University of Utah
 *@brief Contains utility functions for BSSN simulation.
 */

#ifndef SFCSORTBENCH_GRUTILS_H
#define SFCSORTBENCH_GRUTILS_H

#include <iomanip>
#include <limits>

#include "block.h"
#include "dendroProfileParams.h"
#include "grDef.h"
#include "json.hpp"
#include "lebedev.h"
#include "mesh.h"
#include "meshUtils.h"
#include "parUtils.h"
#include "parameters.h"
#include "point.h"
#include "profile_params.h"
#include "swsh.h"

#define RAISE_ERROR(msg)                                                       \
    std::cout << "[Error]: " << __FILE__ << ":" << __LINE__ << " at function " \
              << __FUNCTION__ << " (" << msg << ")" << std::endl

#define USE_NEW_COMMUNICATION_FOR_BH_COMMS

using json = nlohmann::json;
namespace bssn {

/**
 * @brief: Prints compile information, including Git Hash and Compile Date
 *
 * @param[in] rank: MPI Rank
 * @param[in] arg_s: Vectorized string arguments to program
 */
void printGitInformation(int rank, std::vector<std::string> arg_s);

/**
 * @brief: Read the parameter file and initialize the variables in parameters.h
 * file. TOML only -- JSON parameter files are deprecated and rejected with an
 * error.
 * @param[in] fName: file name
 * @param[in] comm: MPI communicator.
 * */
void readParamFile(const char* fName, MPI_Comm comm);

/**
 * @brief dump the read parameter files.
 *
 * @param sout
 * @param root
 * @param comm
 */
void dumpParamFile(std::ostream& sout, int root, MPI_Comm comm);

/**
 * @brief A general initial data function wrapper
 *
 * This function is designed to be called by the grid initialization. It can
 * only take an x, y, and z position on the grid and fill out the proper
 * variables array. If more variables are needed, special logic is required.
 * As a reminder, the grid initialization done by Dendro's function2Octree
 * routine cannot accept additional variables, so approximations are
 * required for initial grid construction and refinement.
 *
 * However, in the bssnCtx object
 * you can add additional logic in init_grid() to more accurately fill the
 * grid with precise values (see BSSN_ID_TYPE 0 and how it calls
 * punctureData() in this function and then calls TwoPunctures() in
 * init_grid())
 *
 * @param xx_grid : the x coord (in octree/grid coords)
 * @param yy_grid : the y coord (in octree/grid coords)
 * @param zz_grid : the z coord (in octree/grid coords)
 * @param var : initialized BSSN variables for the grid at the specified
 * points
 */
void initialDataFunctionWrapper(const double xx_grid, const double yy_grid,
                                const double zz_grid, double* var);

/**
 * @brief Two puncture intiial data from HAD code/
 *
 * @param xx1 : x coord (octree coord)
 * @param yy1 : y coord (octree coord)
 * @param zz1 : z coord (octree coord)
 * @param var : initialized bssn variables for the grid points
 */
void punctureData(const double xx1, const double yy1, const double zz1,
                  double* var);

/**
 * @brief  evaluate HAD initial data from the physical coordinates of the
 * domain. Uses domain coords no coord transformation inside.
 * @param xx1 x coord
 * @param yy1 y coord
 * @param zz1 z coord
 * @param var Initialized data eval at (x,y,z) point.
 */
void punctureDataPhysicalCoord(const double xx, const double yy,
                               const double zz, double* var);

/**
 * @brief compute the static Kerr-Schild BH data
 *
 * @param xx1 : x coord
 * @param yy1 : y coord
 * @param zz1 : z coord
 * @param var : initialized bssn variables for the grid points
 */
void KerrSchildData(const double xx1, const double yy1, const double zz1,
                    double* var);

/**
 * @brief add artificial noise to the initial data.
 * @param xx1 : x coord
 * @param yy1 : y coord
 * @param zz1 : z coord
 * @param var : initialized bssn variables for the grid points
 */
void noiseData(const double xx1, const double yy1, const double zz1,
               double* var);

/**
 * @brief Flat space + smooth Gaussian-enveloped DETERMINISTIC noise (ID 7).
 * Ported from CCZ4's ccz4FlatNoiseInit so BSSN and CCZ4 can be compared on the
 * same, reproducible initial data. Amplitude set by BSSN_NOISE_AMP. CCZ4's
 * U_THETA has no BSSN counterpart and is dropped.
 * @param xx1,yy1,zz1 : grid coords (converted to physical internally)
 * @param var : initialized bssn variables for the grid point
 */
void gaussianNoiseData(const double xx1, const double yy1, const double zz1,
                       double* var);

/**
 * @brief initial flat space data (Minkowski)
 */
void minkowskiInitialData(const double xx1, const double yy1, const double zz1,
                          double* var);

/**
 * @brief fake initial data.
 * @param xx1 : x coord
 * @param yy1 : y coord
 * @param zz1 : z coord
 * @param var : initialized bssn variables for the grid points
 */
void kerrData(double x, double y, double z, double* u);

/**
 * @brief: Generates block adaptive octree for the given binary blockhole
 * problem.
 * @param[out] tmpNodes: created octree tmpNodes
 * @param[in] pt_min: block min point
 * @param[in] pt_max: block max point
 * @param[in] regLev: regular grid level
 * @param[in] maxDepth: maximum refinement level.
 * @param[in] comm: MPI communicator.
 * */
void blockAdaptiveOctree(std::vector<ot::TreeNode>& tmpNodes,
                         const Point& pt_min, const Point& pt_max,
                         const unsigned int regLev, const unsigned int maxDepth,
                         MPI_Comm comm);

/**
 * @brief fake initial data.
 * @param xx1 : x coord
 * @param yy1 : y coord
 * @param zz1 : z coord
 * @param var : initialized bssn variables for the grid points
 */
void fake_initial_data(double x, double y, double z, double* u);

/**
 * @brief: Generates block adaptive octree for the given binary blockhole
 * problem.
 * @param[out] tmpNodes: created octree tmpNodes
 * @param[in] pt_min: block min point
 * @param[in] pt_max: block max point
 * @param[in] regLev: regular grid level
 * @param[in] maxDepth: maximum refinement level.
 * @param[in] comm: MPI communicator.
 * */
void blockAdaptiveOctree(std::vector<ot::TreeNode>& tmpNodes,
                         const Point& pt_min, const Point& pt_max,
                         const unsigned int regLev, const unsigned int maxDepth,
                         MPI_Comm comm);

/**
 * @brief Compute the wavelet tolerance as a function of space.
 *
 * @param x : x coord.
 * @param y : y coord
 * @param z : z coord
 * @param tol_min : min. tolerance value.
 * @return double
 */
double computeWTol(double x, double y, double z, double tol_min);

/**
 * @brief Compute the wavelet tolerance as a function of space (uses actual
 * domain coordinates not octree coordinates)
 *
 * @param x : x coord.
 * @param y : y coord
 * @param z : z coord
 * @param hx : resolution in x,y,z
 * @return double
 */
double computeWTolDCoords(double x, double y, double z, double* hx);

/**
 * @breif: Compute L2 constraint norms.
 */
template <typename T>
double computeConstraintL2Norm(const T* constraintVec, const T* maskVec,
                               unsigned int lbegin, unsigned int lend,
                               MPI_Comm comm);

/**
 * @brief Compute L2 constraint norms.
 */
template <typename T>
double computeConstraintL2Norm(const ot::Mesh* mesh, const T* constraintVec,
                               const T* maskVector, T maskthreshoold);

/**
 * @breif write constraints to a file.
 */
template <typename T>
double extractConstraints(const ot::Mesh* mesh, const T** constraintVar,
                          const T* maskVec, double maskthreshoold,
                          unsigned int timestep, double stime);

/**@brief : write a block to binary*/
void writeBLockToBinary(const double** unzipVarsRHS, unsigned int offset,
                        const double* pmin, const double* pmax, double* bxMin,
                        double* bxMax, const unsigned int* sz,
                        unsigned int blkSz, double dxFactor,
                        const char* fprefix);

/**@brief returns the octant weight for LTS timestepping. */
unsigned int getOctantWeight(const ot::TreeNode* pNode);

/**
 * @brief Compute the BH locations based on the integration of the shift
 * vectors.
 *
 * @param in current locations of the BHs
 * @param out : evolved locations of the BH
 * @param zipVars : zip representation of the current evolution variables.
 * (assumes ghost synced)
 * @param dt : time step size.
 */
void computeBHLocations(const ot::Mesh* pMesh, const Point* in, Point* out,
                        double** zipVars, double dt);

/**
 * @brief Creates mesh object to perform weak scaling, for specified grain size
 * (elements)
 *
 * @param pMesh : current mesh object ptr
 * @param target_npes : Elemental grain size dersired
 * @return ot::Mesh*
 */
ot::Mesh* weakScalingReMesh(ot::Mesh* pMesh, unsigned int target_npes);

/**
 * @brief allocated all the bssn deriv vars for the larges block size
 *
 * @param pMesh
 */
void allocate_bssn_deriv_workspace(const ot::Mesh* pMesh, unsigned int s_fac);

/**@brief deallocates the bssn derivs var space*/
void deallocate_bssn_deriv_workspace();

/**@brief Resolve the OpenMP run-sched ICV for the production RHS block loop
 * (rhs.cpp uses schedule(runtime)). Precedence: parfile BSSN_HYBRID_RHS_SCHEDULE
 * > OMP_SCHEDULE env > historical dynamic,1. Idempotent; no-op unless built with
 * DENDRO_HYBRID_OMP. Called from allocate_bssn_deriv_workspace() (the chokepoint
 * every RHS entry point hits) so the default is dynamic,1 even for binaries that
 * never set OMP_SCHEDULE (e.g. the scaling bench). */
void set_rhs_omp_schedule();

std::tuple<std::string, std::string, std::string> encode_bh_locs(
    const std::vector<std::pair<Point, Point>>& bh_history,
    const std::vector<double>& bh_times);

std::tuple<std::vector<std::pair<Point, Point>>, std::vector<double>>
decode_bh_locs(const std::string& bh1_str, const std::string& bh2_str,
               const std::string& time_str);

/**
 * @brief Order-sensitive fingerprint of the mesh structure, one hash per
 * construction stage, printed on rank 0 as `[fingerprint] <tag> <field> ...`.
 *
 * Guards mesh-construction changes (e.g. threading buildE2EMap / buildE2NMap /
 * performBlocksSetup): the mesh is a deterministic function of the octree, so
 * these hashes MUST NOT change with OMP_NUM_THREADS or with DENDRO_HYBRID_OMP.
 * Hashing is deliberately order-sensitive -- a reordered E2N map is exactly the
 * silent failure we are hunting, and it would change every numeric result.
 *
 * Per-field (not one lumped hash) so a divergence localizes to the stage that
 * caused it: `e2n` differing while `e2e` matches points at buildE2NMap.
 *
 * All ranks must call (gathers to rank 0). Rank-ordered, so it also detects
 * work migrating between ranks.
 */
void meshFingerprint(const ot::Mesh* pMesh, const char* tag);

/**
 * @brief Order-sensitive fingerprint of evolved state bits (all DOFs, all vars).
 *
 * The numerics half of the mesh guard: `meshFingerprint` proves the mesh is
 * identical, this proves the answers are. Hashes raw bit patterns, so it is a
 * bit-exactness check, not a tolerance check -- which is the standard
 * `DENDRO_HYBRID_OMP=OFF` must meet.
 */
void stateFingerprint(const ot::Mesh* pMesh, const DendroScalar* const* vars,
                      unsigned int nVars, const char* tag);

}  // end of namespace bssn

namespace bssn {

namespace timer {

/**@brief initialize all the flop counters. */
void initFlops();

/**@brief clears the snapshot counter for time profiler variables*/
void resetSnapshot();

/**@brief reduce min mean max.
 * @param [in] stat: local time
 * @param [out] stat_g 0-min, 1-mean 2-max
 * */
template <typename T>
void computeOverallStats(T* stat, T* stat_g, MPI_Comm comm) {
    int rank, npes;
    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &rank);

    par::Mpi_Reduce(stat, stat_g, 1, MPI_MIN, 0, comm);
    par::Mpi_Reduce(stat, stat_g + 1, 1, MPI_SUM, 0, comm);
    par::Mpi_Reduce(stat, stat_g + 2, 1, MPI_MAX, 0, comm);
    stat_g[1] /= (npes);
}

/** @breif : printout the profile parameters. */
void profileInfo(const char* filePrefix, const ot::Mesh* pMesh);

/** @breif : printout the profile parameters (intermediate profile information).
 */
void profileInfoIntermediate(const char* filePrefix, const ot::Mesh* pMesh,
                             const unsigned int currentStep);

/** @brief Append one JSONL record to <filePrefix>_steps.jsonl with
 *  min/mean/max-across-ranks timings. Optional ets_ctxpt / app_ctxpt source
 *  per-stage and unzip/zip/rhs from dendrolib's m_uiCtxpt; the legacy
 *  bssn::timer hooks are dead under the ETS evolve path. */
void profileInfoJSON(const char* filePrefix, const ot::Mesh* pMesh,
                     const unsigned int currentStep,
                     const std::vector<profiler_t>* ets_ctxpt = nullptr,
                     const std::vector<profiler_t>* app_ctxpt = nullptr);

}  // namespace timer

}  // namespace bssn

namespace GW {
/**
 * @brief : debug function to write psi4 to interpolated to spheres.
 */
void psi4ShpereDump(const ot::Mesh* mesh, DendroScalar** cVar,
                    unsigned int timestep, double time);

}  // end of namespace GW

#include "grUtils.tcc"

#endif  // SFCSORTBENCH_GRUTILS_H
