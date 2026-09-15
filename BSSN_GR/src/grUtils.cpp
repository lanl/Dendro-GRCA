//
// Created by milinda on 7/26/17.
/**
 *@author Milinda Fernando
 *School of Computing, University of Utah
 *@brief Contains utility functions for BSSN simulation.
 */
//

#include "grUtils.h"

#include <mpi.h>

#include <cstdint>
#include <cstdio>
#include <cstdlib>  // std::getenv for the RHS-schedule resolver
#include <string>
#include <tuple>

#include "base.h"
#include "git_version_and_date.h"
#include "parameters.h"

namespace bssn {

template <typename T>
inline void quickPrintVector(const std::vector<T>& vec, std::ostream& sout,
                             const std::string& prefix) {
    sout << prefix << " : ";
    for (const auto& val : vec) {
        sout << val << " ";
    }
    sout << NRM << std::endl;
}

void printGitInformation(int rank, std::vector<std::string> arg_s) {
    if (!rank) {
        std::cout << YLW << "  COMPILED ON  -  " << compile_info::compileDate
                  << NRM << std::endl;
        std::cout << YLW << "  LATEST GIT HASH - " << compile_info::currGitHash
                  << compile_info::dirtyStatus << NRM << std::endl;
    }

    for (size_t ii = 1; ii < arg_s.size(); ++ii) {
        if (arg_s[ii] == "--compile-info") {
            if (!rank)
                std::cout << "Compile info only flag found, exiting..." << NRM
                          << std::endl
                          << std::endl;
            MPI_Finalize();
            exit(0);
        }
    }
}

void readParamFile(const char* fName, MPI_Comm comm) {
    int rank;
    MPI_Comm_rank(comm, &rank);

    const std::string fNameStr(fName);
    const std::string tomlSuffix = ".toml";

    const bool isToml =
        fNameStr.size() >= tomlSuffix.size() &&
        fNameStr.compare(fNameStr.size() - tomlSuffix.size(),
                         tomlSuffix.size(), tomlSuffix) == 0;

    // JSON parameter files are gone: the hand-written reader duplicated every
    // key already declared in parameters.cpp, so each new parameter had to be
    // added twice. Every .par.json in BSSN_GR/pars/ has a .toml twin.
    if (!isToml) {
        if (!rank) {
            std::cerr << RED << "ERROR: " << NRM << "'" << fNameStr
                      << "' is not a TOML parameter file." << std::endl
                      << "  JSON parameter files are DEPRECATED and no longer "
                         "read."
                      << std::endl
                      << "  Pass a '.toml' parameter file instead (see "
                         "BSSN_GR/pars/*.par.toml)."
                      << std::endl;
        }
        MPI_Barrier(comm);
        MPI_Abort(comm, 1);
    }

    readParamTOMLFile(fName, comm);
}

void dumpParamFile(std::ostream& sout, int root, MPI_Comm comm) {
    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    if (rank == root) {
        sout << "parameters read: " << std::endl;
        sout << YLW << "\tnpes :" << npes << NRM << std::endl;
        sout << YLW << "\tBSSN_DIM :" << bssn::BSSN_DIM << NRM << std::endl;
        sout << YLW << "\tBSSN_ELE_ORDER :" << bssn::BSSN_ELE_ORDER << NRM
             << std::endl;
        sout << YLW << "\tBSSN_PADDING_WIDTH :" << bssn::BSSN_PADDING_WIDTH
             << NRM << std::endl;
        sout << YLW << "\tBSSN_CFL_FACTOR :" << bssn::BSSN_CFL_FACTOR << NRM
             << std::endl;
        sout << YLW << "\tBSSN_IO_OUTPUT_FREQ :" << bssn::BSSN_IO_OUTPUT_FREQ
             << NRM << std::endl;
        sout << YLW << "\tBSSN_GW_EXTRACT_FREQ :" << bssn::BSSN_GW_EXTRACT_FREQ
             << NRM << std::endl;
        sout << YLW
             << "\tBSSN_REMESH_TEST_FREQ :" << bssn::BSSN_REMESH_TEST_FREQ
             << NRM << std::endl;
        sout << YLW << "\tBSSN_CHECKPT_FREQ :" << bssn::BSSN_CHECKPT_FREQ << NRM
             << std::endl;
        sout << YLW << "\tBSSN_RESTORE_SOLVER :" << bssn::BSSN_RESTORE_SOLVER
             << NRM << std::endl;
        sout << YLW << "\tBSSN_ENABLE_BLOCK_ADAPTIVITY :"
             << bssn::BSSN_ENABLE_BLOCK_ADAPTIVITY << NRM << std::endl;
        sout << YLW << "\tBSSN_VTU_FILE_PREFIX :" << bssn::BSSN_VTU_FILE_PREFIX
             << NRM << std::endl;
        sout << YLW
             << "\tBSSN_CHKPT_FILE_PREFIX :" << bssn::BSSN_CHKPT_FILE_PREFIX
             << NRM << std::endl;
        sout << YLW
             << "\tBSSN_PROFILE_FILE_PREFIX :" << bssn::BSSN_PROFILE_FILE_PREFIX
             << NRM << std::endl;
        sout << YLW
             << "\tBSSN_VTU_Z_SLICE_ONLY :" << bssn::BSSN_VTU_Z_SLICE_ONLY
             << NRM << std::endl;
        sout << YLW << "\tBSSN_IO_OUTPUT_GAP :" << bssn::BSSN_IO_OUTPUT_GAP
             << NRM << std::endl;
        sout << YLW << "\tBSSN_DENDRO_GRAIN_SZ :" << bssn::BSSN_DENDRO_GRAIN_SZ
             << NRM << std::endl;
        sout << YLW << "\tBSSN_ASYNC_COMM_K :" << bssn::BSSN_ASYNC_COMM_K << NRM
             << std::endl;
        sout << YLW << "\tBSSN_DENDRO_AMR_FAC :" << bssn::BSSN_DENDRO_AMR_FAC
             << NRM << std::endl;

        sout << YLW << "\tBSSN_DENDRO_AMR_FAC_POST_MERGER: "
             << bssn::BSSN_DENDRO_AMR_FAC_POST_MERGER << NRM << std::endl;

        sout << YLW << "\tBSSN_USE_WAVELET_TOL_FUNCTION :"
             << bssn::BSSN_USE_WAVELET_TOL_FUNCTION << NRM << std::endl;
        sout << YLW << "\tBSSN_WAVELET_TOL :" << bssn::BSSN_WAVELET_TOL << NRM
             << std::endl;
        sout << YLW << "\tBSSN_GW_REFINE_WTOL:" << bssn::BSSN_GW_REFINE_WTOL
             << NRM << std::endl;
        sout << YLW << "\tBSSN_WAVELET_TOL_MAX:" << bssn::BSSN_WAVELET_TOL_MAX
             << NRM << std::endl;
        sout << YLW << "\tBSSN_WAVELET_TOL_FUNCTION_R0: "
             << bssn::BSSN_WAVELET_TOL_FUNCTION_R0 << NRM << std::endl;
        sout << YLW << "\tBSSN_WAVELET_TOL_FUNCTION_R1: "
             << bssn::BSSN_WAVELET_TOL_FUNCTION_R1 << NRM << std::endl;
        sout << YLW << "\tBSSN_LOAD_IMB_TOL :" << bssn::BSSN_LOAD_IMB_TOL << NRM
             << std::endl;
        sout << YLW << "\tBSSN_RK_TIME_BEGIN :" << bssn::BSSN_RK_TIME_BEGIN
             << NRM << std::endl;
        sout << YLW << "\tBSSN_RK_TIME_END :" << bssn::BSSN_RK_TIME_END << NRM
             << std::endl;
        sout << YLW << "\tBSSN_RK_TYPE :" << bssn::BSSN_RK_TYPE << " ("
             << rk_type_name(bssn::BSSN_RK_TYPE) << ")" << NRM << std::endl;
        sout << YLW
             << "\tBSSN_RK45_TIME_STEP_SIZE :" << bssn::BSSN_RK45_TIME_STEP_SIZE
             << NRM << std::endl;
        sout << YLW
             << "\tBSSN_RK45_DESIRED_TOL :" << bssn::BSSN_RK45_DESIRED_TOL
             << NRM << std::endl;
        sout << YLW << "\tBSSN_COMPD_MIN : ( " << bssn::BSSN_COMPD_MIN[0]
             << " ," << bssn::BSSN_COMPD_MIN[1] << ","
             << bssn::BSSN_COMPD_MIN[2] << " )" << NRM << std::endl;
        sout << YLW << "\tBSSN_COMPD_MAX : ( " << bssn::BSSN_COMPD_MAX[0]
             << " ," << bssn::BSSN_COMPD_MAX[1] << ","
             << bssn::BSSN_COMPD_MAX[2] << " )" << NRM << std::endl;
        sout << YLW << "\tBSSN_BLK_MIN : ( " << bssn::BSSN_BLK_MIN_X << " ,"
             << bssn::BSSN_BLK_MIN_Y << "," << bssn::BSSN_BLK_MIN_Z << " )"
             << NRM << std::endl;
        sout << YLW << "\tBSSN_BLK_MAX : ( " << bssn::BSSN_BLK_MAX_X << " ,"
             << bssn::BSSN_BLK_MAX_Y << "," << bssn::BSSN_BLK_MAX_Z << " )"
             << NRM << std::endl;
        sout << YLW << "\tBSSN_OCTREE_MIN : ( " << bssn::BSSN_OCTREE_MIN[0]
             << " ," << bssn::BSSN_OCTREE_MIN[1] << ","
             << bssn::BSSN_OCTREE_MIN[2] << " )" << NRM << std::endl;
        sout << YLW << "\tBSSN_OCTREE_MAX : ( " << bssn::BSSN_OCTREE_MAX[0]
             << " ," << bssn::BSSN_OCTREE_MAX[1] << ","
             << bssn::BSSN_OCTREE_MAX[2] << " )" << NRM << std::endl;
        sout << YLW << "\tETA_CONST :" << bssn::ETA_CONST << NRM << std::endl;
        sout << YLW << "\tETA_R0 :" << bssn::ETA_R0 << NRM << std::endl;
        sout << YLW << "\tETA_DAMPING :" << bssn::ETA_DAMPING << NRM
             << std::endl;
        sout << YLW << "\tETA_DAMPING_EXP :" << bssn::ETA_DAMPING_EXP << NRM
             << std::endl;
        sout << YLW << "\tBSSN_ETA_R0 :" << bssn::BSSN_ETA_R0 << NRM
             << std::endl;
        sout << YLW << "\tBSSN_ETA_POWER : (" << bssn::BSSN_ETA_POWER[0] << " ,"
             << bssn::BSSN_ETA_POWER[1] << " )" << NRM << std::endl;
        sout << YLW << "\tRIT_ETA_FUNCTION: " << bssn::RIT_ETA_FUNCTION << NRM
             << std::endl;
        sout << YLW << "\tRIT_ETA_CENTRAL: " << bssn::RIT_ETA_CENTRAL << NRM
             << std::endl;
        sout << YLW << "\tRIT_ETA_WIDTH: " << bssn::RIT_ETA_WIDTH << NRM
             << std::endl;
        sout << YLW << "\tRIT_ETA_OUTER: " << bssn::RIT_ETA_OUTER << NRM
             << std::endl;
        sout << YLW << "\tBSSN_LAMBDA : (" << bssn::BSSN_LAMBDA[0] << " ,"
             << bssn::BSSN_LAMBDA[1] << "," << bssn::BSSN_LAMBDA[2] << " ,"
             << bssn::BSSN_LAMBDA[3] << " )" << NRM << std::endl;
        sout << YLW << "\tBSSN_LAMBDA_F : (" << bssn::BSSN_LAMBDA_F[0] << " ,"
             << bssn::BSSN_LAMBDA_F[1] << " )" << NRM << std::endl;
        sout << YLW << "\tBSSN_XI : (" << bssn::BSSN_XI[0] << " ,"
             << bssn::BSSN_XI[1] << " ," << bssn::BSSN_XI[2] << " )" << NRM
             << std::endl;
        sout << YLW << "\tCHI_FLOOR :" << bssn::CHI_FLOOR << NRM << std::endl;
        sout << YLW << "\tBSSN_TRK0 :" << bssn::BSSN_TRK0 << NRM << std::endl;
        sout << YLW << "\tDISSIPATION_TYPE :" << bssn::DISSIPATION_TYPE << NRM
             << std::endl;
        sout << YLW << "\tKO_DISS_SIGMA :" << bssn::KO_DISS_SIGMA << NRM
             << std::endl;

        sout << YLW << "\tBH1 MASS :" << bssn::BH1.getBHMass() << NRM
             << std::endl;
        sout << YLW << "\tBH1 POSITION (x,y,z) : (" << bssn::BH1.getBHCoordX()
             << ", " << bssn::BH1.getBHCoordY() << ", "
             << bssn::BH1.getBHCoordZ() << " )" << NRM << std::endl;
        sout << YLW << "\tBH1 VELOCITY (x,y,z) : (" << bssn::BH1.getVx() << ", "
             << bssn::BH1.getVy() << ", " << bssn::BH1.getVz() << " )" << NRM
             << std::endl;
        sout << YLW << "\tBH1 SPIN (||,theta,phi): ( " << bssn::BH1.getBHSpin()
             << ", " << bssn::BH1.getBHSpinTheta() << ", "
             << bssn::BH1.getBHSpinPhi() << " )" << NRM << std::endl;

        sout << YLW << "\tBH2 MASS :" << bssn::BH2.getBHMass() << NRM
             << std::endl;
        sout << YLW << "\tBH2 POSITION (x,y,z) : (" << bssn::BH2.getBHCoordX()
             << ", " << bssn::BH2.getBHCoordY() << ", "
             << bssn::BH2.getBHCoordZ() << " )" << NRM << std::endl;
        sout << YLW << "\tBH2 VELOCITY (x,y,z) : (" << bssn::BH2.getVx() << ", "
             << bssn::BH2.getVy() << ", " << bssn::BH2.getVz() << " )" << NRM
             << std::endl;
        sout << YLW << "\tBH2 SPIN (||,theta,phi): ( " << bssn::BH2.getBHSpin()
             << ", " << bssn::BH2.getBHSpinTheta() << ", "
             << bssn::BH2.getBHSpinPhi() << " )" << NRM << std::endl;

        sout << YLW << "\tBSSN_DIM :" << bssn::BSSN_DIM << NRM << std::endl;
        sout << YLW << "\tBSSN_MAXDEPTH :" << bssn::BSSN_MAXDEPTH << NRM
             << std::endl;
        sout << YLW << "\tBSSN_MINDEPTH :" << bssn::BSSN_MINDEPTH << NRM
             << std::endl;

        sout << YLW << "\tBSSN_NUM_REFINE_VARS :" << bssn::BSSN_NUM_REFINE_VARS
             << NRM << std::endl;
        sout << YLW << "\tBSSN_REFINE_VARIABLE_INDICES :[";
        for (unsigned int i = 0; i < bssn::BSSN_NUM_REFINE_VARS - 1; i++)
            sout << bssn::BSSN_REFINE_VARIABLE_INDICES[i] << ", ";
        sout << bssn::BSSN_REFINE_VARIABLE_INDICES[bssn::BSSN_NUM_REFINE_VARS -
                                                   1]
             << "]" << NRM << std::endl;

        sout << YLW << "\tBSSN_REFINEMENT_MODE :" << bssn::BSSN_REFINEMENT_MODE
             << NRM << std::endl;
        sout << YLW << "\tBSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE :"
             << bssn::BSSN_USE_SET_REF_MODE_FOR_INITIAL_CONVERGE << NRM
             << std::endl;

        sout << YLW << "\tBSSN_BH1_AMR_R: " << bssn::BSSN_BH1_AMR_R << NRM
             << std::endl;
        sout << YLW << "\tBSSN_BH2_AMR_R: " << bssn::BSSN_BH2_AMR_R << NRM
             << std::endl;
        sout << YLW << "\tBSSN_AMR_R_RATIO: " << bssn::BSSN_AMR_R_RATIO << NRM
             << std::endl;

        sout << YLW << "\tBSSN_BH1_CONSTRAINT_R:" << bssn::BSSN_BH1_CONSTRAINT_R
             << NRM << std::endl;
        sout << YLW << "\tBSSN_BH2_CONSTRAINT_R:" << bssn::BSSN_BH2_CONSTRAINT_R
             << NRM << std::endl;

        sout << YLW << "\tBSSN_BH1_MAX_LEV:" << bssn::BSSN_BH1_MAX_LEV << NRM
             << std::endl;
        sout << YLW << "\tBSSN_BH2_MAX_LEV:" << bssn::BSSN_BH2_MAX_LEV << NRM
             << std::endl;
        sout << YLW << "\tBSSN_INIT_GRID_ITER:" << bssn::BSSN_INIT_GRID_ITER
             << NRM << std::endl;

#ifdef BSSN_REFINE_BASE_EH
        sout << YLW << "\tBSSN_EH_REFINE_VAL  : " << bssn::BSSN_EH_REFINE_VAL
             << NRM << std::endl;
        sout << YLW << "\tBSSN_EH_COARSEN_VAL : " << bssn::BSSN_EH_COARSEN_VAL
             << NRM << std::endl;
#endif

        sout << YLW << "\tBSSN_NUM_EVOL_VARS_VTU_OUTPUT :"
             << bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT << NRM << std::endl;
        sout << YLW << "\tBSSN_VTU_OUTPUT_EVOL_INDICES :[";
        for (unsigned int i = 0; i < bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT - 1;
             i++)
            sout << bssn::BSSN_VTU_OUTPUT_EVOL_INDICES[i] << ", ";
        sout << bssn::BSSN_VTU_OUTPUT_EVOL_INDICES
                    [bssn::BSSN_NUM_EVOL_VARS_VTU_OUTPUT - 1]
             << "]" << NRM << std::endl;

        sout << YLW << "\tBSSN_NUM_CONST_VARS_VTU_OUTPUT :"
             << bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT << NRM << std::endl;
        sout << YLW << "\tBSSN_VTU_OUTPUT_CONST_INDICES :[";
        for (unsigned int i = 0; i < bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT - 1;
             i++)
            sout << bssn::BSSN_VTU_OUTPUT_CONST_INDICES[i] << ", ";
        sout << bssn::BSSN_VTU_OUTPUT_CONST_INDICES
                    [bssn::BSSN_NUM_CONST_VARS_VTU_OUTPUT - 1]
             << "]" << NRM << std::endl;

        sout << YLW << "\tTPID_TARGET_M_PLUS :" << TPID::target_M_plus << NRM
             << std::endl;
        sout << YLW << "\tTPID_TARGET_M_MINUS :" << TPID::target_M_minus << NRM
             << std::endl;
        sout << YLW << "\tTPID_PAR_B :" << TPID::par_b << NRM << std::endl;
        sout << YLW << "\tTPID_PAR_P_PLUS : ( " << TPID::par_P_plus[0] << ", "
             << TPID::par_P_plus[1] << ", " << TPID::par_P_plus[2] << " )"
             << NRM << std::endl;
        sout << YLW << "\tTPID_PAR_P_MINUS : ( " << TPID::par_P_minus[0] << ", "
             << TPID::par_P_minus[1] << ", " << TPID::par_P_minus[2] << " )"
             << NRM << std::endl;
        sout << YLW << "\tTPID_PAR_S_PLUS : ( " << TPID::par_S_plus[0] << ", "
             << TPID::par_S_plus[1] << ", " << TPID::par_S_plus[2] << " )"
             << NRM << std::endl;
        sout << YLW << "\tTPID_PAR_S_MINUS : ( " << TPID::par_S_minus[0] << ", "
             << TPID::par_S_minus[1] << ", " << TPID::par_S_minus[2] << " )"
             << NRM << std::endl;
        sout << YLW << "\tTPID_CENTER_OFFSET : ( " << TPID::center_offset[0]
             << ", " << TPID::center_offset[1] << ", " << TPID::center_offset[2]
             << " )" << NRM << std::endl;

        sout << YLW << "\tTPID_INITIAL_LAPSE_PSI_EXPONENT :"
             << TPID::initial_lapse_psi_exponent << NRM << std::endl;
        sout << YLW << "\tTPID_NPOINTS_A :" << TPID::npoints_A << NRM
             << std::endl;
        sout << YLW << "\tTPID_NPOINTS_B :" << TPID::npoints_B << NRM
             << std::endl;
        sout << YLW << "\tTPID_NPOINTS_PHI :" << TPID::npoints_phi << NRM
             << std::endl;
        sout << YLW << "\tTPID_GIVE_BARE_MASS :" << TPID::give_bare_mass << NRM
             << std::endl;
        sout << YLW << "\tINITIAL_LAPSE :" << TPID::initial_lapse << NRM
             << std::endl;
        sout << YLW << "\tTPID_SOLVE_MOMENTUM_CONSTRAINT :"
             << TPID::solve_momentum_constraint << NRM << std::endl;
        sout << YLW << "\tTPID_GRID_SETUP_METHOD :" << TPID::grid_setup_method
             << NRM << std::endl;
        sout << YLW << "\tTPID_VERBOSE :" << TPID::verbose << NRM << std::endl;
        sout << YLW << "\tTPID_ADM_TOL :" << TPID::adm_tol << NRM << std::endl;
        sout << YLW << "\tTPID_NEWTON_TOL :" << TPID::Newton_tol << NRM
             << std::endl;

        sout << YLW << "\tEXTRACTION_VAR_ID :" << BHLOC::EXTRACTION_VAR_ID
             << NRM << std::endl;
        sout << YLW << "\tEXTRACTION_TOL :" << BHLOC::EXTRACTION_TOL << NRM
             << std::endl;

        sout << YLW << "\tBSSN_GW_NUM_RADAII: " << GW::BSSN_GW_NUM_RADAII << NRM
             << std::endl;
        sout << YLW << "\tBSSN_GW_NUM_LMODES: " << GW::BSSN_GW_NUM_LMODES << NRM
             << std::endl;

        sout << YLW << "\tBSSN_GW_RADAII: {";
        for (unsigned int i = 0; i < GW::BSSN_GW_NUM_RADAII; i++)
            sout << " ," << GW::BSSN_GW_RADAII[i];
        sout << "}" << NRM << std::endl;

        sout << YLW << "\tBSSN_GW_L_MODES: {";
        for (unsigned int i = 0; i < GW::BSSN_GW_NUM_LMODES; i++)
            sout << " ," << GW::BSSN_GW_L_MODES[i];
        sout << "}" << NRM << std::endl;

        sout << YLW << "\tBSSN_KO_SIGMA_SCALE_BY_CONFORMAL: "
             << (bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL ? "true" : "false")
             << NRM << std::endl;
        if (bssn::BSSN_KO_SIGMA_SCALE_BY_CONFORMAL) {
            sout << YLW << "\t\tBSSN_PSILON_CAKO_GAUGE: "
                 << bssn::BSSN_EPSILON_CAKO_GAUGE << NRM << std::endl;
            sout << YLW << "\t\tBSSN_PSILON_CAKO_OTHER: "
                 << bssn::BSSN_EPSILON_CAKO_OTHER << NRM << std::endl;
        }
        sout << YLW << "\tBSSN_NYQUIST_M: " << bssn::BSSN_NYQUIST_M << NRM
             << std::endl;

        sout << YLW << "\tBSSN_SSL_H: " << bssn::BSSN_SSL_H << NRM << std::endl;
        sout << YLW << "\tBSSN_SSL_SIGMA: " << bssn::BSSN_SSL_SIGMA << NRM
             << std::endl;

        sout << GRN << "\t----- AEH PARAMETERS -----" << NRM << std::endl;
        sout << YLW << "\t\tAEH_SOLVER_FREQ: " << AEH::AEH_SOLVER_FREQ << NRM
             << std::endl;
        sout << YLW << "\t\tN_HORIZONS: " << AEH::N_HORIZONS << NRM
             << std::endl;
        sout << YLW
             << "\t\tN_RESOLUTIONS_MULTIGRID: " << AEH::N_RESOLUTIONS_MULTIGRID
             << NRM << std::endl;
        quickPrintVector(AEH::INITIAL_X_CENTER, sout,
                         std::string(YLW) + "\t\tINITIAL_X_CENTER");
        quickPrintVector(AEH::INITIAL_Y_CENTER, sout,
                         std::string(YLW) + "\t\tINITIAL_Y_CENTER");
        quickPrintVector(AEH::INITIAL_Z_CENTER, sout,
                         std::string(YLW) + "\t\tINITIAL_Z_CENTER");
        quickPrintVector(AEH::M_SCALE, sout, std::string(YLW) + "\t\tM_SCALE");
        quickPrintVector(AEH::CFL_FACTOR, sout,
                         std::string(YLW) + "\t\tCFL_FACTOR");
        quickPrintVector(AEH::THETA_L2_M_TOL, sout,
                         std::string(YLW) + "\t\tTHETA_L2_M_TOL");
        quickPrintVector(AEH::THETA_LINF_M_TOL, sout,
                         std::string(YLW) + "\t\tTHETA_LINF_M_TOL");
        quickPrintVector(AEH::ETA_DAMP_M, sout,
                         std::string(YLW) + "\t\tETA_DAMP_M");
        quickPrintVector(AEH::KO_STRENGTH, sout,
                         std::string(YLW) + "\t\tKO_STRENGTH");
        quickPrintVector(AEH::MAX_SEARCH_RADIUS, sout,
                         std::string(YLW) + "\t\tMAX_SEARCH_RADIUS");
        quickPrintVector(AEH::NR_INTERP_MAX, sout,
                         std::string(YLW) + "\t\tNR_INTERP_MAX");
        sout << YLW << "\t\tNTHETA_MAX: " << AEH::NTHETA_MAX << NRM
             << std::endl;
        quickPrintVector(AEH::NTHETA_ARRAY, sout,
                         std::string(YLW) + "\t\tNTHETA_ARRAY");
        sout << YLW << "\t\tNPHI_MAX: " << AEH::NTHETA_MAX << NRM << std::endl;
        quickPrintVector(AEH::NPHI_ARRAY, sout,
                         std::string(YLW) + "\t\tNPHI_ARRAY");
        sout << YLW << "\t\tAEH_SAVE_DIR: " << AEH::AEH_SAVE_DIR << NRM
             << std::endl;
        sout << YLW << "\t\tNUM_RESOLUTIONS_AFTER_FIND: "
             << AEH::NUM_RESOLUTIONS_AFTER_FIND << NRM << std::endl;
        sout << YLW
             << "\t\tENABLE_ETA_VARYING_ALG: " << AEH::ENABLE_ETA_VARYING_ALG
             << NRM << std::endl;
        sout << YLW << "\t\tVERBOSITY_LEVEL: " << AEH::VERBOSITY_LEVEL << NRM
             << std::endl;
    }
}

double CalTolHelper(const double t, const double r, const double rad[],
                    const double eps[], const double toffset) {
    const double R0     = rad[0];
    const double R1     = rad[1];
    const double RGW    = rad[2];
    const double tol    = eps[0];
    const double tolGW  = eps[1];
    const double tolMax = eps[2];
    const double WRR    = std::min(tolGW, tolMax);
    if (r < R0) {
        return tol;
    }
    if (t > (R1 + toffset)) {
        const double RR         = std::min(t - toffset, RGW + 10.0);
        const double WTolExpFac = (RR - R0) / log10(WRR / tol);
        return std::min(tolMax, tol * pow(10.0, ((r - R0) / WTolExpFac)));
    } else {
        const double WTolExpFac = (R1 - R0) / log10(WRR / tol);
        return std::min(tolMax, tol * pow(10.0, ((r - R0) / WTolExpFac)));
    }
}

void initialDataFunctionWrapper(const double xx_grid, const double yy_grid,
                                const double zz_grid, double* var) {
    // code to convert grid functions to non-grid, if the function doesn't
    // already support it const double xx = GRIDX_TO_X(xx_grid); const double yy
    // = GRIDY_TO_Y(yy_grid); const double zz = GRIDZ_TO_Z(zz_grid);

    switch (bssn::BSSN_ID_TYPE) {
        case 0:
            // NOTE: this is the TwoPunctures code! For **pure**
            // initialization done when building up the grid before pure
            // population, this calls bssn::punctureData which is simply the
            // two puncture intiial data from HAD code. In bssnCtx.cpp's
            // init_grid the bssn::BSSN_ID_TYPE switch will call the true
            // TwoPunctures code which is the proper TPID initial data
            bssn::punctureData(xx_grid, yy_grid, zz_grid, var);

            break;

        case 1:
            // ID 1 is the puncture data function on its own
            bssn::punctureData(xx_grid, yy_grid, zz_grid, var);

            break;

        case 2:
            // ID 2 is the KerrSchild Data
            bssn::KerrSchildData(xx_grid, yy_grid, zz_grid, var);

            break;

        case 3:
            // ID 3 is a noise data
            bssn::noiseData(xx_grid, yy_grid, zz_grid, var);

            break;

        case 4:
            // ID 4 is a fake initial data
            bssn::fake_initial_data(xx_grid, yy_grid, zz_grid, var);

            break;

        case 5:
        // minkowski initial data is flat space!
        bssn:
            minkowskiInitialData(xx_grid, yy_grid, zz_grid, var);

            break;
        case 6:
            // minkowski initial data is flat space!
            bssn::kerrData(xx_grid, yy_grid, zz_grid, var);

            break;
            // MORE CAN BE ADDED HERE

        default:
            int rank;
            int npes;
            MPI_Comm comm = MPI_COMM_WORLD;
            MPI_Comm_rank(comm, &rank);
            MPI_Comm_size(comm, &npes);
            if (!rank) {
                std::cerr << RED << "ERROR::: Invalid initial data ID: "
                          << bssn::BSSN_ID_TYPE << NRM << std::endl;
            }

            MPI_Abort(comm, 0);

            break;
    }
}

void punctureDataPhysicalCoord(const double xx, const double yy,
                               const double zz, double* var) {
    /* Define the Levi-Civita pseudotensor and Kronecker delta */
    double epijk[3][3][3];
    int i, j, k;

    for (k = 0; k < 3; k++) {
        for (j = 0; j < 3; j++) {
            for (i = 0; i < 3; i++) {
                epijk[k][j][i] = 0.0;
            }
        }
    }

    epijk[0][1][2] = 1.0;
    epijk[1][2][0] = 1.0;
    epijk[2][0][1] = 1.0;
    epijk[0][2][1] = -1.0;
    epijk[2][1][0] = -1.0;
    epijk[1][0][2] = -1.0;

    double deltaij[3][3];

    for (j = 0; j < 3; j++) {
        for (i = 0; i < 3; i++) {
            deltaij[j][i] = 0.0;
        }
    }

    deltaij[0][0] = 1.0;
    deltaij[1][1] = 1.0;
    deltaij[2][2] = 1.0;

    double x1, y1, z1, rv1;
    double x2, y2, z2, rv2;
    double vn1[3], vn2[3];

    double vpsibl;
    double v_u_corr;
    double amp_capj, amp_capr, l_r;
    double u0_j, u2_j, mu_j, p2_mu_j, v_u_j1;
    double v1, v2, v3, v4, vt1, vt2;

    int i1, i2, i3, i4;

    double amp_capp;
    double u0_p, u2_p, mu_p, p2_mu_p;
    double v_u_p1, v_u_c1, v_u_j2, v_u_p2;
    double v_u_c2, vpsibl_u, vpsibl_u2;

    constexpr double PUNCTURE_EPS = 1.0e-14;

    // BH1
    const double mass1 = BH1.getBHMass();
    const double bh1x  = BH1.getBHCoordX();
    const double bh1y  = BH1.getBHCoordY();
    const double bh1z  = BH1.getBHCoordZ();

    double vp1[3];
    vp1[0] = BH1.getVx();
    vp1[1] = BH1.getVy();
    vp1[2] = BH1.getVz();

    const double vp1tot =
        sqrt(vp1[0] * vp1[0] +
             vp1[1] * vp1[1] +
             vp1[2] * vp1[2]);

    const double spin1     = BH1.getBHSpin();
    const double spin1_th  = BH1.getBHSpinTheta();
    const double spin1_phi = BH1.getBHSpinPhi();

    double vs1[3];
    vs1[0] = spin1 * sin(spin1_th) * cos(spin1_phi);
    vs1[1] = spin1 * sin(spin1_th) * sin(spin1_phi);
    vs1[2] = spin1 * cos(spin1_th);

    // BH2
    const double mass2 = BH2.getBHMass();
    const double bh2x  = BH2.getBHCoordX();
    const double bh2y  = BH2.getBHCoordY();
    const double bh2z  = BH2.getBHCoordZ();

    double vp2[3];
    vp2[0] = BH2.getVx();
    vp2[1] = BH2.getVy();
    vp2[2] = BH2.getVz();

    const double vp2tot =
        sqrt(vp2[0] * vp2[0] +
             vp2[1] * vp2[1] +
             vp2[2] * vp2[2]);

    const double spin2     = BH2.getBHSpin();
    const double spin2_th  = BH2.getBHSpinTheta();
    const double spin2_phi = BH2.getBHSpinPhi();

    double vs2[3];
    vs2[0] = spin2 * sin(spin2_th) * cos(spin2_phi);
    vs2[1] = spin2 * sin(spin2_th) * sin(spin2_phi);
    vs2[2] = spin2 * cos(spin2_th);

    // Coordinates relative to BH1
    x1  = xx - bh1x;
    y1  = yy - bh1y;
    z1  = zz - bh1z;
    rv1 = sqrt(x1 * x1 + y1 * y1 + z1 * z1);

    // Coordinates relative to BH2
    x2  = xx - bh2x;
    y2  = yy - bh2y;
    z2  = zz - bh2z;
    rv2 = sqrt(x2 * x2 + y2 * y2 + z2 * z2);

    const bool at_puncture1 = rv1 <= PUNCTURE_EPS;
    const bool at_puncture2 = rv2 <= PUNCTURE_EPS;

    /*
     * At either puncture, the radial unit vector is undefined and the
     * Brill-Lindquist conformal factor diverges. The evolved variables,
     * however, have well defined limiting values:
     *
     *     alpha -> 0,
     *     chi   -> 0,
     *     A_ij  -> 0.
     *
     * Apply those limits directly to avoid 0/0 and inf*0 operations.
     */
    if (at_puncture1 || at_puncture2) {
        var[VAR::U_ALPHA] = CHI_FLOOR;
        var[VAR::U_CHI]   = CHI_FLOOR;
        var[VAR::U_K]     = 0.0;

        var[VAR::U_BETA0] = 0.0;
        var[VAR::U_BETA1] = 0.0;
        var[VAR::U_BETA2] = 0.0;

        var[VAR::U_GT0] = 0.0;
        var[VAR::U_GT1] = 0.0;
        var[VAR::U_GT2] = 0.0;

        var[VAR::U_B0] = 0.0;
        var[VAR::U_B1] = 0.0;
        var[VAR::U_B2] = 0.0;

        var[VAR::U_SYMGT0] = 1.0;  // XX
        var[VAR::U_SYMGT1] = 0.0;  // XY
        var[VAR::U_SYMGT2] = 0.0;  // XZ
        var[VAR::U_SYMGT3] = 1.0;  // YY
        var[VAR::U_SYMGT4] = 0.0;  // YZ
        var[VAR::U_SYMGT5] = 1.0;  // ZZ

        var[VAR::U_SYMAT0] = 0.0;  // XX
        var[VAR::U_SYMAT1] = 0.0;  // XY
        var[VAR::U_SYMAT2] = 0.0;  // XZ
        var[VAR::U_SYMAT3] = 0.0;  // YY
        var[VAR::U_SYMAT4] = 0.0;  // YZ
        var[VAR::U_SYMAT5] = 0.0;  // ZZ

        return;
    }

    // Radial unit vectors are safe after the puncture check
    vn1[0] = x1 / rv1;
    vn1[1] = y1 / rv1;
    vn1[2] = z1 / rv1;

    vn2[0] = x2 / rv2;
    vn2[1] = y2 / rv2;
    vn2[2] = z2 / rv2;

    // Initial data is related to arXiv:0711.1165
    // Brill-Lindquist conformal factor
    vpsibl = 1.0 + mass1 / (2.0 * rv1);
    vpsibl = vpsibl + mass2 / (2.0 * rv2);

    v_u_corr = 0.0;

    // BH1 spinning puncture correction
    if (fabs(spin1) > 1.0e-6) {
        amp_capj = 4.0 * spin1 / (mass1 * mass1);
        amp_capr = 2.0 * rv1 / mass1;
        l_r      = 1.0 / (1.0 + amp_capr);

        u0_j =
            (l_r +
             l_r * l_r +
             l_r * l_r * l_r -
             4.0 * l_r * l_r * l_r * l_r +
             2.0 * l_r * l_r * l_r * l_r * l_r) /
            40.0;

        u2_j = -pow(l_r, 5) / 20.0;

        mu_j = vn1[0] * vs1[0];
        mu_j = mu_j + vn1[1] * vs1[1];
        mu_j = (mu_j + vn1[2] * vs1[2]) / fabs(spin1);

        p2_mu_j = (3.0 * mu_j * mu_j - 1.0) / 2.0;

        v_u_j1 =
            amp_capj * amp_capj *
            (u0_j + u2_j * amp_capr * amp_capr * p2_mu_j);

        v_u_corr = v_u_corr + v_u_j1;
    }

    // BH1 boosted puncture correction
    if (vp1tot > 1.0e-6) {
        amp_capp = 2.0 * vp1tot / mass1;
        amp_capr = 2.0 * rv1 / mass1;
        l_r      = 1.0 / (1.0 + amp_capr);

        u0_p = l_r -
               2.0 * l_r * l_r +
               2.0 * pow(l_r, 3);

        u0_p =
            (u0_p -
             pow(l_r, 4) +
             0.20 * pow(l_r, 5)) *
            (5.0 / 32.0);

        u2_p =
            15.0 * l_r +
            132.0 * l_r * l_r +
            53.0 * pow(l_r, 3);

        u2_p =
            u2_p +
            96.0 * pow(l_r, 4) +
            82.0 * pow(l_r, 5);

        u2_p =
            u2_p +
            (84.0 / amp_capr) *
                (pow(l_r, 5) + log(l_r) / amp_capr);

        u2_p = u2_p / (80.0 * amp_capr);

        mu_p = vn1[0] * vp1[0] / vp1tot;
        mu_p = mu_p + vn1[1] * vp1[1] / vp1tot;
        mu_p = mu_p + vn1[2] * vp1[2] / vp1tot;

        p2_mu_p = (3.0 * pow(mu_p, 2) - 1.0) / 2.0;

        v_u_p1 =
            pow(amp_capp, 2) *
            (u0_p + u2_p * p2_mu_p);

        v_u_corr = v_u_corr + v_u_p1;
    }

    // BH1 spinning and boosted puncture correction
    if (vp1tot > 1.0e-6 && fabs(spin1) > 1.0e-6) {
        v1 =
            (vp1[1] * vs1[2] - vp1[2] * vs1[1]) * vn1[0];

        v1 =
            v1 +
            (vp1[2] * vs1[0] - vp1[0] * vs1[2]) * vn1[1];

        v1 =
            v1 +
            (vp1[0] * vs1[1] - vp1[1] * vs1[0]) * vn1[2];

        v1 = v1 * (16.0 / pow(mass1, 4)) * rv1;

        amp_capr = 2.0 * rv1 / mass1;
        l_r      = 1.0 / (1.0 + amp_capr);

        v2 = 1.0 +
             5.0 * amp_capr +
             10.0 * pow(amp_capr, 2);

        v_u_c1 = (v1 * v2 * pow(l_r, 5)) / 80.0;

        v_u_corr = v_u_corr + v_u_c1;
    }

    // BH2 spinning puncture correction
    if (fabs(spin2) > 1.0e-6) {
        amp_capj = 4.0 * spin2 / (mass2 * mass2);
        amp_capr = 2.0 * rv2 / mass2;
        l_r      = 1.0 / (1.0 + amp_capr);

        u0_j =
            (l_r +
             l_r * l_r +
             l_r * l_r * l_r -
             4.0 * l_r * l_r * l_r * l_r +
             2.0 * l_r * l_r * l_r * l_r * l_r) /
            40.0;

        u2_j = -pow(l_r, 5) / 20.0;

        mu_j = vn2[0] * vs2[0];
        mu_j = mu_j + vn2[1] * vs2[1];
        mu_j = (mu_j + vn2[2] * vs2[2]) / fabs(spin2);

        p2_mu_j = (3.0 * mu_j * mu_j - 1.0) / 2.0;

        v_u_j2 =
            amp_capj * amp_capj *
            (u0_j + u2_j * amp_capr * amp_capr * p2_mu_j);

        v_u_corr = v_u_corr + v_u_j2;
    }

    // BH2 boosted puncture correction
    if (vp2tot > 1.0e-6) {
        amp_capp = 2.0 * vp2tot / mass2;
        amp_capr = 2.0 * rv2 / mass2;
        l_r      = 1.0 / (1.0 + amp_capr);

        u0_p = l_r -
               2.0 * l_r * l_r +
               2.0 * pow(l_r, 3);

        u0_p =
            (u0_p -
             pow(l_r, 4) +
             0.20 * pow(l_r, 5)) *
            (5.0 / 32.0);

        u2_p =
            15.0 * l_r +
            132.0 * l_r * l_r +
            53.0 * pow(l_r, 3);

        u2_p =
            u2_p +
            96.0 * pow(l_r, 4) +
            82.0 * pow(l_r, 5);

        u2_p =
            u2_p +
            (84.0 / amp_capr) *
                (pow(l_r, 5) + log(l_r) / amp_capr);

        u2_p = u2_p / (80.0 * amp_capr);

        mu_p = vn2[0] * vp2[0] / vp2tot;
        mu_p = mu_p + vn2[1] * vp2[1] / vp2tot;
        mu_p = mu_p + vn2[2] * vp2[2] / vp2tot;

        p2_mu_p = (3.0 * pow(mu_p, 2) - 1.0) / 2.0;

        v_u_p2 =
            pow(amp_capp, 2) *
            (u0_p + u2_p * p2_mu_p);

        v_u_corr = v_u_corr + v_u_p2;
    }

    // BH2 spinning and boosted puncture correction
    if (vp2tot > 1.0e-6 && fabs(spin2) > 1.0e-6) {
        v1 =
            (vp2[1] * vs2[2] - vp2[2] * vs2[1]) * vn2[0];

        v1 =
            v1 +
            (vp2[2] * vs2[0] - vp2[0] * vs2[2]) * vn2[1];

        v1 =
            v1 +
            (vp2[0] * vs2[1] - vp2[1] * vs2[0]) * vn2[2];

        v1 = v1 * (16.0 / pow(mass2, 4)) * rv2;

        amp_capr = 2.0 * rv2 / mass2;
        l_r      = 1.0 / (1.0 + amp_capr);

        v2 = 1.0 +
             5.0 * amp_capr +
             10.0 * pow(amp_capr, 2);

        v_u_c2 = (v1 * v2 * pow(l_r, 5)) / 80.0;

        v_u_corr = v_u_corr + v_u_c2;
    }

    // Corrected conformal factors
    vpsibl_u  = vpsibl + v_u_corr;
    vpsibl_u2 = vpsibl + v_u_corr;

    var[VAR::U_ALPHA] =
        1.0 / (vpsibl_u * vpsibl_u);

    var[VAR::U_ALPHA] =
        std::max(var[VAR::U_ALPHA], CHI_FLOOR);

    var[VAR::U_CHI] =
        1.0 / pow(vpsibl_u, 4);

    if (var[VAR::U_CHI] < CHI_FLOOR) {
        var[VAR::U_CHI] = CHI_FLOOR;
    }

    var[VAR::U_K] = 0.0;

    var[VAR::U_BETA0] = 0.0;
    var[VAR::U_BETA1] = 0.0;
    var[VAR::U_BETA2] = 0.0;

    var[VAR::U_GT0] = 0.0;
    var[VAR::U_GT1] = 0.0;
    var[VAR::U_GT2] = 0.0;

    var[VAR::U_B0] = 0.0;
    var[VAR::U_B1] = 0.0;
    var[VAR::U_B2] = 0.0;

    var[VAR::U_SYMGT0] = 1.0;  // XX
    var[VAR::U_SYMGT1] = 0.0;  // XY
    var[VAR::U_SYMGT2] = 0.0;  // XZ
    var[VAR::U_SYMGT3] = 1.0;  // YY
    var[VAR::U_SYMGT4] = 0.0;  // YZ
    var[VAR::U_SYMGT5] = 1.0;  // ZZ

    for (i1 = 0; i1 < 3; i1++) {
        for (i2 = 0; i2 < 3; i2++) {
            /*
             * BH1 contribution to the physical conformal
             * trace-free extrinsic curvature.
             */
            v2 = 0.0;

            for (i3 = 0; i3 < 3; i3++) {
                for (i4 = 0; i4 < 3; i4++) {
                    vt1 =
                        epijk[i1][i3][i4] *
                        vs1[i3] *
                        vn1[i4] *
                        vn1[i2];

                    vt2 =
                        epijk[i2][i3][i4] *
                        vs1[i3] *
                        vn1[i4] *
                        vn1[i1];

                    v2 = v2 + vt1 + vt2;
                }
            }

            v3 =
                vp1[i1] * vn1[i2] +
                vp1[i2] * vn1[i1];

            vt1 = 0.0;

            for (i3 = 0; i3 < 3; i3++) {
                vt1 = vt1 + vp1[i3] * vn1[i3];
            }

            vt1 =
                vt1 *
                (vn1[i1] * vn1[i2] -
                 deltaij[i1][i2]);

            v3 = v3 + vt1;

            v1 =
                3.0 /
                (pow(vpsibl_u2, 6) *
                 pow(rv1, 3));

            v4 =
                v1 *
                (v2 + (rv1 / 2.0) * v3);

            /*
             * BH2 contribution to the physical conformal
             * trace-free extrinsic curvature.
             */
            v2 = 0.0;

            for (i3 = 0; i3 < 3; i3++) {
                for (i4 = 0; i4 < 3; i4++) {
                    vt1 =
                        epijk[i1][i3][i4] *
                        vs2[i3] *
                        vn2[i4] *
                        vn2[i2];

                    vt2 =
                        epijk[i2][i3][i4] *
                        vs2[i3] *
                        vn2[i4] *
                        vn2[i1];

                    v2 = v2 + vt1 + vt2;
                }
            }

            v3 =
                vp2[i1] * vn2[i2] +
                vp2[i2] * vn2[i1];

            vt1 = 0.0;

            for (i3 = 0; i3 < 3; i3++) {
                vt1 = vt1 + vp2[i3] * vn2[i3];
            }

            vt1 =
                vt1 *
                (vn2[i1] * vn2[i2] -
                 deltaij[i1][i2]);

            v3 = v3 + vt1;

            v1 =
                3.0 /
                (pow(vpsibl_u2, 6) *
                 pow(rv2, 3));

            v4 =
                v4 +
                v1 *
                    (v2 + (rv2 / 2.0) * v3);

            if (i1 == 0 && i2 == 0) {
                var[VAR::U_SYMAT0] = v4;  // XX
            } else if (i1 == 0 && i2 == 1) {
                var[VAR::U_SYMAT1] = v4;  // XY
            } else if (i1 == 0 && i2 == 2) {
                var[VAR::U_SYMAT2] = v4;  // XZ
            } else if (i1 == 1 && i2 == 1) {
                var[VAR::U_SYMAT3] = v4;  // YY
            } else if (i1 == 1 && i2 == 2) {
                var[VAR::U_SYMAT4] = v4;  // YZ
            } else if (i1 == 2 && i2 == 2) {
                var[VAR::U_SYMAT5] = v4;  // ZZ
            }
        }
    }
}
void punctureData(const double xx1, const double yy1, const double zz1,
                  double* var) {
    const double xx = GRIDX_TO_X(xx1);
    const double yy = GRIDY_TO_Y(yy1);
    const double zz = GRIDZ_TO_Z(zz1);

    punctureDataPhysicalCoord(xx, yy, zz, var);
}
void kerrData(const double xx1, const double yy1, const double zz1,
              double* var) {
    const double xx = GRIDX_TO_X(xx1);
    const double yy = GRIDY_TO_Y(yy1);
    const double zz = GRIDZ_TO_Z(zz1);

    // parameters for the BH (mass, location, spin parameter)
    double M        = BH1.getBHMass();
    double bh1x     = BH1.getBHCoordX();
    double bh1y     = BH1.getBHCoordY();
    double bh1z     = BH1.getBHCoordZ();
    double spin1    = BH1.getBHSpin();

    // coordinates relative to the center of the BH
    double x        = xx - bh1x;
    double y        = yy - bh1y;
    double z        = zz - bh1z;

    // locating as a radial form
    double r        = sqrt(x * x + y * y + z * z);

    // HL : Angular momentum parameter will be added as param file after
    // testing
    double a        = spin1;

    double gtd[3][3], Atd[3][3];
    double alpha, Gamt[3];
    double Chi, TrK, Betau[3];

#include "Kerr.cpp"
#include "kerr_vars.cpp"
}

void KerrSchildData(const double xx1, const double yy1, const double zz1,
                    double* var) {
    const double xx = GRIDX_TO_X(xx1);
    const double yy = GRIDY_TO_Y(yy1);
    const double zz = GRIDZ_TO_Z(zz1);

    // parameters for the BH (mass, location, spin parameter)
    double M        = BH1.getBHMass();
    double bh1x     = BH1.getBHCoordX();
    double bh1y     = BH1.getBHCoordY();
    double bh1z     = BH1.getBHCoordZ();
    double spin1    = BH1.getBHSpin();

    // coordinates relative to the center of the BH
    double x        = xx - bh1x;
    double y        = yy - bh1y;
    double z        = zz - bh1z;

    // locating as a radial form
    double r        = sqrt(x * x + y * y + z * z);

    // HL : Angular momentum parameter will be added as param file after testing
    double a        = spin1;

    double gtd[3][3], Atd[3][3];
    double alpha, Gamt[3];
    double Chi, TrK, Betau[3];

#include "ks_vars.cpp"
#include "ksinit.cpp"

    var[VAR::U_ALPHA]  = alpha;
    var[VAR::U_CHI]    = Chi;
    var[VAR::U_K]      = TrK;

    var[VAR::U_BETA0]  = Betau[0];
    var[VAR::U_BETA1]  = Betau[1];
    var[VAR::U_BETA2]  = Betau[2];

    var[VAR::U_GT0]    = Gamt[0];
    var[VAR::U_GT1]    = Gamt[1];
    var[VAR::U_GT2]    = Gamt[2];

    var[VAR::U_B0]     = 0.0;
    var[VAR::U_B1]     = 0.0;
    var[VAR::U_B2]     = 0.0;

    var[VAR::U_SYMGT0] = gtd[0][0];
    var[VAR::U_SYMGT1] = gtd[0][1];
    var[VAR::U_SYMGT2] = gtd[0][2];
    var[VAR::U_SYMGT3] = gtd[1][1];
    var[VAR::U_SYMGT4] = gtd[1][2];
    var[VAR::U_SYMGT5] = gtd[2][2];

    var[VAR::U_SYMAT0] = Atd[0][0];
    var[VAR::U_SYMAT1] = Atd[0][1];
    var[VAR::U_SYMAT2] = Atd[0][2];
    var[VAR::U_SYMAT3] = Atd[1][1];
    var[VAR::U_SYMAT4] = Atd[1][2];
    var[VAR::U_SYMAT5] = Atd[2][2];

    // std::cout<<"KS init data: (x,y,z) = ( "<<x<<", "<<y<<", "<<z<<"), alpha =
    // "<<alpha<<std::endl;

#if 0
            //BSSN vars for Kerr-Schild
            var[VAR::U_ALPHA] = sqrt(rv1/(2.0*M+rv1));
            var[VAR::U_CHI] = 1.0/pow(1.0+2.0*M/rv1, 1.0/3.0);
            var[VAR::U_K] = 2.0*M*sqrt(rv1/(2.0*M+rv1))*(rv1+3.0*M)/(rv1*rv1*(2.0*M+rv1));

            var[VAR::U_BETA0] = 2.0*M*x1/(rv1*(2.0*M+rv1));
            var[VAR::U_BETA1] = 2.0*M*y1/(rv1*(2.0*M+rv1));
            var[VAR::U_BETA2] = 2.0*M*z1/(rv1*(2.0*M+rv1));

            var[VAR::U_GT0] = pow(2,8.0/3.0)*M*x1*(M+rv1*rv1+(3*M*M*rv1+x1*x1+y1*y1)/5.0)/(pow(5.0,1.0/3.0)*rv1*pow(M/5.0+rv1,2.0)*pow(M/(5.0*rv1)+1.0,2.0/3.0));
            var[VAR::U_GT1] = pow(2,8.0/3.0)*M*y1*(M+rv1*rv1+(3*M*M*rv1+x1*x1+y1*y1)/5.0)/(pow(5.0,1.0/3.0)*rv1*pow(M/5.0+rv1,2.0)*pow(M/(5.0*rv1)+1.0,2.0/3.0));
            var[VAR::U_GT2] = pow(2,8.0/3.0)*M*z1*(M+rv1*rv1+(3*M*M*rv1+x1*x1+y1*y1)/5.0)/(pow(5.0,1.0/3.0)*rv1*pow(M/5.0+rv1,2.0)*pow(M/(5.0*rv1)+1.0,2.0/3.0));

            var[VAR::U_B0] = M*x1/(500.0*rv1*(M/5.0+rv1));
            var[VAR::U_B1] = M*y1/(500.0*rv1*(M/5.0+rv1));
            var[VAR::U_B2] = M*z1/(500.0*rv1*(M/5.0+rv1));

            var[VAR::U_SYMGT0] = (M*x1*x1/2.0+rv1*rv1)/pow(10,5.0/3.0)*rv1*rv1*pow(M/(5.0*rv1)+1.0,1.0/3.0); //XX
            var[VAR::U_SYMGT1] = M*x1*y1/(50.0*pow(10.0,2.0/3.0)*rv1*rv1*pow(M/(5.0*rv1)+1.0,1.0/3.0)); //XY
            var[VAR::U_SYMGT2] = M*x1*z1/(50.0*pow(10.0,2.0/3.0)*rv1*rv1*pow(M/(5.0*rv1)+1.0,1.0/3.0)); //XZ
            var[VAR::U_SYMGT3] = (M*y1*y1/2.0+rv1*rv1)/pow(10,5.0/3.0)*rv1*rv1*pow(M/(5.0*rv1)+1.0,1.0/3.0); //YY
            var[VAR::U_SYMGT4] = M*y1*z1/(50.0*pow(10.0,2.0/3.0)*rv1*rv1*pow(M/(5.0*rv1)+1.0,1.0/3.0)); //YZ
            var[VAR::U_SYMGT5] = (M*z1*z1/2.0+rv1*rv1)/pow(10,5.0/3.0)*rv1*rv1*pow(M/(5.0*rv1)+1.0,1.0/3.0); //ZZ

            var[VAR::U_SYMAT0] = (sqrt(rv1/(M/5.0+rv1))*(rv1-(1.0/5.0+M/(10.0*rv1))*x1*x1)-M*sqrt(rv1/(M/5.0+rv1))*(3.0*M/10.0+rv1)*(1.0/10.0+M*x1*x1/(50.0*rv1*rv1))/(150.0*sqrt(10.0)*(M/5.0+rv1)*rv1))/pow(0.1+M/(50.0*rv1),1.0/3.0); //XX
            var[VAR::U_SYMAT1] = -(((1.0/5.0+M/(10.0*rv1))*sqrt(rv1/(M/5.0+rv1))*x1*y1)/(50.0*sqrt(10.0)*M*rv1)+(M*M*sqrt(rv1/(M/5.0+rv1))*(3.0*M/10.0*rv1)*x1*y1)/(750.0*sqrt(10.0)*rv1*rv1*(M/5.0*rv1)))/pow(0.1+M/(50.0*rv1),1.0/3.0); //XY
            var[VAR::U_SYMAT2] = -(((1.0/5.0+M/(10.0*rv1))*sqrt(rv1/(M/5.0+rv1))*x1*y1)/(50.0*sqrt(10.0)*M*rv1)+(M*M*sqrt(rv1/(M/5.0+rv1))*(3.0*M/10.0*rv1)*x1*z1)/(750.0*sqrt(10.0)*rv1*rv1*(M/5.0*rv1)))/pow(0.1+M/(50.0*rv1),1.0/3.0); //XZ
            var[VAR::U_SYMAT3] = (sqrt(rv1/(M/5.0+rv1))*(rv1-(1.0/5.0+M/(10.0*rv1))*y1*y1)-M*sqrt(rv1/(M/5.0+rv1))*(3.0*M/10.0+rv1)*(1.0/10.0+M*x1*x1/(50.0*rv1*rv1))/(150.0*sqrt(10.0)*(M/5.0+rv1)*rv1))/pow(0.1+M/(50.0*rv1),1.0/3.0); //YY
            var[VAR::U_SYMAT4] = -(((1.0/5.0+M/(10.0*rv1))*sqrt(rv1/(M/5.0+rv1))*x1*y1)/(50.0*sqrt(10.0)*M*rv1)+(M*M*sqrt(rv1/(M/5.0+rv1))*(3.0*M/10.0*rv1)*y1*z1)/(750.0*sqrt(10.0)*rv1*rv1*(M/5.0*rv1)))/pow(0.1+M/(50.0*rv1),1.0/3.0); //YZ
            var[VAR::U_SYMAT5] = (sqrt(rv1/(M/5.0+rv1))*(rv1-(1.0/5.0+M/(10.0*rv1))*z1*z1)-M*sqrt(rv1/(M/5.0+rv1))*(3.0*M/10.0+rv1)*(1.0/10.0+M*x1*x1/(50.0*rv1*rv1))/(150.0*sqrt(10.0)*(M/5.0+rv1)*rv1))/pow(0.1+M/(50.0*rv1),1.0/3.0); //ZZ
#endif
}

void noiseData(const double xx1, const double yy1, const double zz1,
               double* var) {
    // const double xx=GRIDX_TO_X(xx1);
    // const double yy=GRIDY_TO_Y(yy1);
    // const double zz=GRIDZ_TO_Z(zz1);

    // call random number generator between -1 and 1)
    double random_variable[30];
    int i;
    for (i = 0; i < 30; i++) {
        random_variable[i] = 2.0 * rand() / ((double)RAND_MAX) - 1.0;
    }

    // set a (uniform) amplitude for the noise
    double noise_amp   = bssn::BSSN_NOISE_AMP;

    var[VAR::U_ALPHA]  = 1.0 + noise_amp * random_variable[0];

    var[VAR::U_CHI]    = 1.0 + noise_amp * random_variable[1];

    var[VAR::U_K]      = noise_amp * random_variable[2];

    var[VAR::U_GT0]    = noise_amp * random_variable[3];
    var[VAR::U_GT1]    = noise_amp * random_variable[4];
    var[VAR::U_GT2]    = noise_amp * random_variable[5];

    var[VAR::U_BETA0]  = noise_amp * random_variable[6];
    var[VAR::U_BETA1]  = noise_amp * random_variable[7];
    var[VAR::U_BETA2]  = noise_amp * random_variable[8];

    var[VAR::U_B0]     = noise_amp * random_variable[9];
    var[VAR::U_B1]     = noise_amp * random_variable[10];
    var[VAR::U_B2]     = noise_amp * random_variable[11];

    var[VAR::U_SYMGT0] = 1.0 + noise_amp * random_variable[12];  // XX
    var[VAR::U_SYMGT1] = noise_amp * random_variable[13];        // XY
    var[VAR::U_SYMGT2] = noise_amp * random_variable[14];        // XZ
    var[VAR::U_SYMGT3] = 1.0 + noise_amp * random_variable[15];  // YY
    var[VAR::U_SYMGT4] = noise_amp * random_variable[16];        // YZ
    var[VAR::U_SYMGT5] = 1.0 + noise_amp * random_variable[17];  // ZZ

    var[VAR::U_SYMAT0] = noise_amp * random_variable[18];  // XX
    var[VAR::U_SYMAT1] = noise_amp * random_variable[19];  // XY
    var[VAR::U_SYMAT2] = noise_amp * random_variable[20];  // XZ
    var[VAR::U_SYMAT3] = noise_amp * random_variable[21];  // YY
    var[VAR::U_SYMAT4] = noise_amp * random_variable[22];  // YZ
    var[VAR::U_SYMAT5] = noise_amp * random_variable[23];  // ZZ
}

void fake_initial_data(double xx1, double yy1, double zz1, double* u) {
    /* const double x=GRIDX_TO_X(xx1);
     const double y=GRIDY_TO_Y(yy1);
     const double z=GRIDZ_TO_Z(zz1);


     const double pi = acos(-1.0);
     const double f1 = 31.0/17.0;
     const double f2 = 37.0/11.0;

     u[VAR::U_ALPHA] = 1.0 - 0.25*sin(f1*x);
     //u[F_ALPHA][pp] = 1.0;
     u[VAR::U_BETA0] = 4.0/17.0*sin(x)*cos(z);
     u[VAR::U_BETA1] = pi/5.0*cos(y)*sin(z+x);
     u[VAR::U_BETA2] = 4.0/17.0*sin(f2*x)*sin(z);

     u[VAR::U_B0] = 31.0*x*cos(f1*z+y);
     u[VAR::U_B1] = 7.0*y*sin(f1*x+y) + 3.0*cos(z);
     u[VAR::U_B2] = 5.0*z*cos(f1*x+y) + 7.0*sin(z+y+x) + 1.0;

     u[VAR::U_GT0] = 5.0*cos(x)/(10.0*sin(x+z)+26.0-1.0*cos(x*z)*cos(x));
     u[VAR::U_GT1] = -5.0*sin(y)/(25.0+10.0*cos(y+z)+cos(y)*cos(y*z));
     u[VAR::U_GT2] = -5.0*sin(z)/(25.0+10.0*cos(y+x)+cos(y*x)*cos(z));

     u[VAR::U_CHI] = 1.0 + exp(-4.0*cos(x)*sin(y));
     //u[F_CHI][pp] = 2.0;

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

     u[VAR::U_SYMAT0] = exp(-4.0*cos(x)*sin(y))*(cos(x)
     -0.3333333333*exp(4.0*cos(x)*sin(y))
     *(1.0+0.2*sin(x))*(5.0*exp(-4.0*cos(x)*sin(y))
     /(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))
     /(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)
     +5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
     u[VAR::U_SYMAT1] = 1.0 + x*z/(0.1 + x*x + y*y + z*z);
     u[VAR::U_SYMAT2] = 1.3 - x*y/(3.0 + x*x + 2.0*y*y + z*z)*(x*x+z*z);
     u[VAR::U_SYMAT3] =
     exp(-4.0*cos(x)*sin(y))*(cos(y)-0.33333333330*exp(4*cos(x)*sin(y))*(1+0.2*cos(y))*(5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
     u[VAR::U_SYMAT4] = -1.0 + y*z/(1.0 + 3.0*x*x + y*y + z*z);
     u[VAR::U_SYMAT5] =
     exp(-4.0*cos(x)*sin(y))*(cos(z)-0.3333333333*exp(4*cos(x)*sin(y))/(1+0.2*sin(x))/(1+0.2*cos(y))*(5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
     */

    const double x  = GRIDX_TO_X(xx1);
    const double y  = GRIDY_TO_Y(yy1);
    const double z  = GRIDZ_TO_Z(zz1);

    const double pi = acos(-1.0);
    const double f1 = 31.0 / 17.0;
    const double f2 = 37.0 / 11.0;

    u[VAR::U_ALPHA] = 1.0 - 0.25 * sin(f1 * x);
    // u[F_ALPHA][pp] = 1.0;
    u[VAR::U_BETA0] = 4.0 / 17.0 * sin(x) * cos(z);
    u[VAR::U_BETA1] = pi / 5.0 * cos(y) * sin(z + x);
    u[VAR::U_BETA2] = 4.0 / 17.0 * sin(f2 * x) * sin(z);

    u[VAR::U_B0]    = 31.0 * x * cos(f1 * z + y);
    u[VAR::U_B1]    = 7.0 * y * sin(f1 * x + y) + 3.0 * cos(z);
    u[VAR::U_B2]    = 5.0 * z * cos(f1 * x + y) + 7.0 * sin(z + y + x) + 1.0;

    u[VAR::U_GT0] =
        5.0 * cos(x) / (10.0 * sin(x + z) + 26.0 - 1.0 * cos(x * z) * cos(x));
    u[VAR::U_GT1] =
        -5.0 * sin(y) / (25.0 + 10.0 * cos(y + z) + cos(y) * cos(y * z));
    u[VAR::U_GT2] =
        -5.0 * sin(z) / (25.0 + 10.0 * cos(y + x) + cos(y * x) * cos(z));

    u[VAR::U_CHI]    = 1.0 + exp(-4.0 * cos(x) * sin(y));
    // u[F_CHI][pp] = 2.0;

    u[VAR::U_SYMGT0] = 1.00 + 0.2 * sin(x + z) * cos(y);
    u[VAR::U_SYMGT3] = 1.00 + 0.2 * cos(y) * cos(z + x);
    u[VAR::U_SYMGT5] = 1.00 / (u[VAR::U_SYMGT0] + u[VAR::U_SYMGT3]);
    u[VAR::U_SYMGT1] = 0.07 * (2.0 + cos(x * x + y * y));
    u[VAR::U_SYMGT2] = 0.1 * (3.0 + sin(z) * cos(x));
    u[VAR::U_SYMGT4] = 0.15 * (1.751 - sin(x * x) * cos(y) * cos(z));

    u[VAR::U_K] = 5.0 * exp(-4.0 * cos(x) * sin(y)) / (5.0 + sin(x)) * cos(x) +
                  5.0 * exp(-4.0 * cos(x) * sin(y)) / (5.0 + cos(y)) * cos(y) +
                  0.4 * (25.0 + 5.0 * cos(y) + 5.0 * sin(x) + sin(x) * cos(y)) *
                      exp(-4.0 * cos(x) * sin(y)) * cos(z);
    u[VAR::U_K] *= 0.01234;

    u[VAR::U_SYMAT0] =
        exp(-4.0 * cos(x) * sin(y)) *
        (cos(x) -
         0.3333333333 * exp(4.0 * cos(x) * sin(y)) * (1.0 + 0.2 * sin(x)) *
             (5.0 * exp(-4.0 * cos(x) * sin(y)) / (5.0 + sin(x)) * cos(x) +
              5.0 * exp(-4.0 * cos(x) * sin(y)) / (5.0 + cos(y)) * cos(y) +
              0.04 * (25.0 + 5.0 * cos(y) + 5.0 * sin(x) + sin(x) * cos(y)) *
                  exp(-4.0 * cos(x) * sin(y)) * cos(z)));
    u[VAR::U_SYMAT1] = 1.0 + x * z / (0.1 + x * x + y * y + z * z);
    u[VAR::U_SYMAT2] =
        1.3 - x * y / (3.0 + x * x + 2.0 * y * y + z * z) * (x * x + z * z);
    u[VAR::U_SYMAT3] =
        exp(-4.0 * cos(x) * sin(y)) *
        (cos(y) -
         0.33333333330 * exp(4 * cos(x) * sin(y)) * (1 + 0.2 * cos(y)) *
             (5.0 * exp(-4.0 * cos(x) * sin(y)) / (5.0 + sin(x)) * cos(x) +
              5.0 * exp(-4.0 * cos(x) * sin(y)) / (5.0 + cos(y)) * cos(y) +
              0.04 * (25.0 + 5.0 * cos(y) + 5.0 * sin(x) + sin(x) * cos(y)) *
                  exp(-4.0 * cos(x) * sin(y)) * cos(z)));
    u[VAR::U_SYMAT4] = -1.0 + y * z / (1.0 + 3.0 * x * x + y * y + z * z);
    u[VAR::U_SYMAT5] =
        exp(-4.0 * cos(x) * sin(y)) *
        (cos(z) -
         0.3333333333 * exp(4 * cos(x) * sin(y)) / (1 + 0.2 * sin(x)) /
             (1 + 0.2 * cos(y)) *
             (5.0 * exp(-4.0 * cos(x) * sin(y)) / (5.0 + sin(x)) * cos(x) +
              5.0 * exp(-4.0 * cos(x) * sin(y)) / (5.0 + cos(y)) * cos(y) +
              0.04 * (25.0 + 5.0 * cos(y) + 5.0 * sin(x) + sin(x) * cos(y)) *
                  exp(-4.0 * cos(x) * sin(y)) * cos(z)));

    /* Enforce BSSN constraints */
    double gtd[3][3], Atd[3][3];

    gtd[0][0]              = u[VAR::U_SYMGT0];
    gtd[0][1]              = u[VAR::U_SYMGT1];
    gtd[0][2]              = u[VAR::U_SYMGT2];
    gtd[1][0]              = gtd[0][1];
    gtd[1][1]              = u[VAR::U_SYMGT3];
    gtd[1][2]              = u[VAR::U_SYMGT4];
    gtd[2][0]              = gtd[0][2];
    gtd[2][1]              = gtd[1][2];
    gtd[2][2]              = u[VAR::U_SYMGT5];

    Atd[0][0]              = u[VAR::U_SYMAT0];
    Atd[0][1]              = u[VAR::U_SYMAT1];
    Atd[0][2]              = u[VAR::U_SYMAT2];
    Atd[1][0]              = Atd[0][1];
    Atd[1][1]              = u[VAR::U_SYMAT3];
    Atd[1][2]              = u[VAR::U_SYMAT4];
    Atd[2][0]              = Atd[0][2];
    Atd[2][1]              = Atd[1][2];
    Atd[2][2]              = u[VAR::U_SYMAT5];

    const double one_third = 1.0 / 3.0;
    double det_gtd =
        gtd[0][0] * (gtd[1][1] * gtd[2][2] - gtd[1][2] * gtd[1][2]) -
        gtd[0][1] * gtd[0][1] * gtd[2][2] +
        2.0 * gtd[0][1] * gtd[0][2] * gtd[1][2] -
        gtd[0][2] * gtd[0][2] * gtd[1][1];

    if (det_gtd < 0.0) {
        /* FIXME What to do here? The metric is not physical. Do we reset the
         * metric to be flat? */
        gtd[0][0] = 1.0;
        gtd[0][1] = 0.0;
        gtd[0][2] = 0.0;
        gtd[1][0] = 0.0;
        gtd[1][1] = 1.0;
        gtd[1][2] = 0.0;
        gtd[2][0] = 0.0;
        gtd[2][1] = 0.0;
        gtd[2][2] = 1.0;
        det_gtd   = 1.0;
    }
    double det_gtd_to_neg_third = 1.0 / pow(det_gtd, one_third);

    for (unsigned int j = 0; j < 3; j++) {
        for (unsigned int i = 0; i < 3; i++) {
            gtd[i][j] *= det_gtd_to_neg_third;
        }
    }

    det_gtd = gtd[0][0] * (gtd[1][1] * gtd[2][2] - gtd[1][2] * gtd[1][2]) -
              gtd[0][1] * gtd[0][1] * gtd[2][2] +
              2.0 * gtd[0][1] * gtd[0][2] * gtd[1][2] -
              gtd[0][2] * gtd[0][2] * gtd[1][1];

    double detgt_m1 = det_gtd - 1.0;

    if (fabs(detgt_m1) > 1.0e-6) {
        std::cout.precision(14);
        std::cout << "enforce_bssn_constraint: det(gtd) != 1. det="
                  << std::fixed << det_gtd << std::endl;
        std::cout << "      gtd(1,1)=" << gtd[0][0] << std::endl;
        std::cout << "      gtd(1,2)=" << gtd[0][1] << std::endl;
        std::cout << "      gtd(1,3)=" << gtd[0][2] << std::endl;
        std::cout << "      gtd(2,2)=" << gtd[1][1] << std::endl;
        std::cout << "      gtd(2,3)=" << gtd[1][2] << std::endl;
        std::cout << "      gtd(3,3)=" << gtd[2][2] << std::endl;
    }

    double gtu[3][3];
    double idet_gtd = 1.0 / det_gtd;
    gtu[0][0] = idet_gtd * (gtd[1][1] * gtd[2][2] - gtd[1][2] * gtd[1][2]);
    gtu[0][1] = idet_gtd * (-gtd[0][1] * gtd[2][2] + gtd[0][2] * gtd[1][2]);
    gtu[0][2] = idet_gtd * (gtd[0][1] * gtd[1][2] - gtd[0][2] * gtd[1][1]);
    gtu[1][0] = gtu[0][1];
    gtu[1][1] = idet_gtd * (gtd[0][0] * gtd[2][2] - gtd[0][2] * gtd[0][2]);
    gtu[1][2] = idet_gtd * (-gtd[0][0] * gtd[1][2] + gtd[0][1] * gtd[0][2]);
    gtu[2][0] = gtu[0][2];
    gtu[2][1] = gtu[1][2];
    gtu[2][2] = idet_gtd * (gtd[0][0] * gtd[1][1] - gtd[0][1] * gtd[0][1]);

    /* Require Atd to be traceless. */
    double one_third_trace_Atd =
        one_third *
        (Atd[0][0] * gtu[0][0] + Atd[1][1] * gtu[1][1] + Atd[2][2] * gtu[2][2] +
         2.0 * (Atd[0][1] * gtu[0][1] + Atd[0][2] * gtu[0][2] +
                Atd[1][2] * gtu[1][2]));

    Atd[0][0] -= one_third_trace_Atd * gtd[0][0];
    Atd[0][1] -= one_third_trace_Atd * gtd[0][1];
    Atd[0][2] -= one_third_trace_Atd * gtd[0][2];
    Atd[1][1] -= one_third_trace_Atd * gtd[1][1];
    Atd[1][2] -= one_third_trace_Atd * gtd[1][2];
    Atd[2][2] -= one_third_trace_Atd * gtd[2][2];

    double tr_A = Atd[0][0] * gtu[0][0] + Atd[1][1] * gtu[1][1] +
                  Atd[2][2] * gtu[2][2] +
                  2.0 * (Atd[0][1] * gtu[0][1] + Atd[0][2] * gtu[0][2] +
                         Atd[1][2] * gtu[1][2]);

    if (fabs(tr_A) > 1.0e-6) {
        std::cout << "enforce_bssn_constraint: tr_A != 0. tr_A=" << tr_A
                  << std::endl;
        std::cout << "      Atd(1,1)=" << Atd[0][0] << std::endl;
        std::cout << "      Atd(1,2)=" << Atd[0][1] << std::endl;
        std::cout << "      Atd(1,3)=" << Atd[0][2] << std::endl;
        std::cout << "      Atd(2,2)=" << Atd[1][1] << std::endl;
        std::cout << "      Atd(2,3)=" << Atd[1][2] << std::endl;
        std::cout << "      Atd(3,3)=" << Atd[2][2] << std::endl;
    }

    u[VAR::U_SYMAT0] = Atd[0][0];
    u[VAR::U_SYMAT1] = Atd[0][1];
    u[VAR::U_SYMAT2] = Atd[0][2];
    u[VAR::U_SYMAT3] = Atd[1][1];
    u[VAR::U_SYMAT4] = Atd[1][2];
    u[VAR::U_SYMAT5] = Atd[2][2];

    u[VAR::U_SYMGT0] = gtd[0][0];
    u[VAR::U_SYMGT1] = gtd[0][1];
    u[VAR::U_SYMGT2] = gtd[0][2];
    u[VAR::U_SYMGT3] = gtd[1][1];
    u[VAR::U_SYMGT4] = gtd[1][2];
    u[VAR::U_SYMGT5] = gtd[2][2];
}

void minkowskiInitialData(const double xx1, const double yy1, const double zz1,
                          double* var) {
    // Flat space initialization!
    var[VAR::U_ALPHA]  = 1;  // lapse
    var[VAR::U_CHI]    = 1;  // chi
    var[VAR::U_K]      = 0;  // trace K
    var[VAR::U_GT0]    = 0;  // Gt0
    var[VAR::U_GT1]    = 0;  // Gt1
    var[VAR::U_GT2]    = 0;  // Gt2
    var[VAR::U_BETA0]  = 0;  // shift 0
    var[VAR::U_BETA1]  = 0;  // shift 1
    var[VAR::U_BETA2]  = 0;  // shift 2
    var[VAR::U_B0]     = 0;  // gaugeB0
    var[VAR::U_B1]     = 0;  // gaugeB1
    var[VAR::U_B2]     = 0;  // gaugeB2
    var[VAR::U_SYMGT0] = 1;  // gt11
    var[VAR::U_SYMGT1] = 0;  // gt12
    var[VAR::U_SYMGT2] = 0;  // gt13
    var[VAR::U_SYMGT3] = 1;  // gt22
    var[VAR::U_SYMGT4] = 0;  // gt23
    var[VAR::U_SYMGT5] = 1;  // gt33
    var[VAR::U_SYMAT0] = 0;  // At11
    var[VAR::U_SYMAT1] = 0;  // At12
    var[VAR::U_SYMAT2] = 0;  // At13
    var[VAR::U_SYMAT3] = 0;  // At22
    var[VAR::U_SYMAT4] = 0;  // At23
    var[VAR::U_SYMAT5] = 0;  // At33
}

void blockAdaptiveOctree(std::vector<ot::TreeNode>& tmpNodes,
                         const Point& pt_min, const Point& pt_max,
                         const unsigned int regLev, const unsigned int maxDepth,
                         MPI_Comm comm) {
    int rank, npes;
    MPI_Comm_size(comm, &npes);
    MPI_Comm_rank(comm, &rank);

    double pt_g_min[3];
    double pt_g_max[3];

    pt_g_min[0] = X_TO_GRIDX(pt_min.x());
    pt_g_min[1] = Y_TO_GRIDY(pt_min.y());
    pt_g_min[2] = Z_TO_GRIDZ(pt_min.z());

    pt_g_max[0] = X_TO_GRIDX(pt_max.x());
    pt_g_max[1] = Y_TO_GRIDY(pt_max.y());
    pt_g_max[2] = Z_TO_GRIDZ(pt_max.z());

    assert(pt_g_min[0] >= 0 && pt_g_min[0] <= (1u << maxDepth));
    assert(pt_g_min[1] >= 0 && pt_g_min[1] <= (1u << maxDepth));
    assert(pt_g_min[2] >= 0 && pt_g_min[2] <= (1u << maxDepth));

    assert(pt_g_max[0] >= 0 && pt_g_max[0] <= (1u << maxDepth));
    assert(pt_g_max[1] >= 0 && pt_g_max[1] <= (1u << maxDepth));
    assert(pt_g_max[2] >= 0 && pt_g_max[2] <= (1u << maxDepth));

    unsigned int xRange_b, xRange_e;
    unsigned int yRange_b = pt_g_min[1], yRange_e = pt_g_max[1];
    unsigned int zRange_b = pt_g_min[2], zRange_e = pt_g_max[2];

    xRange_b =
        pt_g_min[0];  //(rank*(pt_g_max[0]-pt_g_min[0]))/npes + pt_g_min[0];
    xRange_e =
        pt_g_max[1];  //((rank+1)*(pt_g_max[0]-pt_g_min[0]))/npes + pt_g_min[0];

    unsigned int stepSz = 1u << (maxDepth - regLev);

    /* std::cout<<" x min: "<<xRange_b<<" x_max: "<<xRange_e<<std::endl;
     std::cout<<" y min: "<<yRange_b<<" y_max: "<<yRange_e<<std::endl;
     std::cout<<" z min: "<<zRange_b<<" z_max: "<<zRange_e<<std::endl;*/

    for (unsigned int x = xRange_b; x < xRange_e; x += stepSz)
        for (unsigned int y = yRange_b; y < yRange_e; y += stepSz)
            for (unsigned int z = zRange_b; z < zRange_e; z += stepSz) {
                if (x >= (1u << maxDepth)) x = x - 1;
                if (y >= (1u << maxDepth)) y = y - 1;
                if (z >= (1u << maxDepth)) z = z - 1;

                tmpNodes.push_back(
                    ot::TreeNode(x, y, z, regLev, m_uiDim, maxDepth));
            }

    return;
}

double computeWTol(double x, double y, double z, double tolMin) {
    double origin[3];
    origin[0] = (double)(1u << bssn::BSSN_MAXDEPTH - 1);
    origin[1] = (double)(1u << bssn::BSSN_MAXDEPTH - 1);
    origin[2] = (double)(1u << bssn::BSSN_MAXDEPTH - 1);

    double r =
        sqrt(GRIDX_TO_X(x) * GRIDX_TO_X(x) + GRIDY_TO_Y(y) * GRIDY_TO_Y(y) +
             GRIDZ_TO_Z(z) * GRIDZ_TO_Z(z));

    const double tolMax = bssn::BSSN_WAVELET_TOL_MAX;
    const double R0     = bssn::BSSN_WAVELET_TOL_FUNCTION_R0;
    const double R1     = bssn::BSSN_WAVELET_TOL_FUNCTION_R1;

    if (bssn::BSSN_USE_WAVELET_TOL_FUNCTION == 1) {
        return std::min(
            tolMax, std::max(tolMin, (tolMax - tolMin) / (R1 - R0) * (r - R0) +
                                         tolMin));
    } else {
        return tolMin;
    }
}

double computeWTolDCoords(double x, double y, double z, double* hx) {
    // set up a few useful values for computing the wavelet tolerances
    // element order: how many points we have in each grid
    const unsigned int eleOrder = bssn::BSSN_ELE_ORDER;
    // current simulation time
    const double T_CURRENT      = bssn::BSSN_CURRENT_RK_COORD_TIME;
    // radius (from the center of the grid)
    const double r              = sqrt(x * x + y * y + z * z);
    // distance between the BHs
    const double dbh = (bssn::BSSN_BH_LOC[0] - bssn::BSSN_BH_LOC[1]).abs();
    // set up grid point for relative distances to each BH
    Point grid_p(x, y, z);
    // distance from BH0
    const double dbh0 = (grid_p - bssn::BSSN_BH_LOC[0]).abs();
    // distance from BH1
    const double dbh1 = (grid_p - bssn::BSSN_BH_LOC[1]).abs();

    if (bssn::BSSN_USE_WAVELET_TOL_FUNCTION == 1) {
        const double tolMax = bssn::BSSN_WAVELET_TOL_MAX;
        const double tolMin = bssn::BSSN_WAVELET_TOL;

        const double R0     = bssn::BSSN_BH1_AMR_R;
        const double R1     = bssn::BSSN_BH2_AMR_R;

        // R_Max is defined based on the initial separation.
        const double R_MAX =
            (bssn::BH1.getBHCoord() - bssn::BH2.getBHCoord()).abs() + R0 + R1;

#ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES
        if (dbh < 0.1) {
            if ((dbh0 > R_MAX) && (dbh1 > R_MAX)) {
                if (r < (GW::BSSN_GW_RADAII[GW::BSSN_GW_NUM_RADAII - 1] + 10))
                    return BSSN_GW_REFINE_WTOL;
                else
                    return tolMax;
            } else {
                for (unsigned int k = 0; k < (eleOrder + 1); k++)
                    for (unsigned int j = 0; j < (eleOrder + 1); j++)
                        for (unsigned int i = 0; i < (eleOrder + 1); i++) {
                            const double xx = x + i * hx[0];
                            const double yy = y + j * hx[1];
                            const double zz = z + k * hx[2];

                            const Point grid_pp(xx, yy, zz);

                            const double dd0 =
                                (grid_pp - bssn::BSSN_BH_LOC[0]).abs();
                            const double dd1 =
                                (grid_pp - bssn::BSSN_BH_LOC[1]).abs();

                            // std::cout<<"x : "<<x<<" y: "<<y<<" z: "<<z<<"
                            // dd0: "<<dd0<<" dd1: "<<dd1<<" hx:
                            // "<<hx[0]<<std::endl;

                            if (dd0 < R0 || dd1 < R1) {
                                // std::cout<<"x : "<<x<<" y: "<<y<<" z: "<<z<<"
                                // dd0: "<<dd0<<" dd1: "<<dd1<<" hx:
                                // "<<hx[0]<<std::endl;
                                return tolMin;
                            }
                        }

                if (r < (GW::BSSN_GW_RADAII[GW::BSSN_GW_NUM_RADAII - 1] + 10))
                    return BSSN_GW_REFINE_WTOL;
                else
                    return tolMax;
            }

        } else {
            if ((dbh0 > R_MAX) &&
                (dbh1 >
                 R_MAX))  // no need to check individual points in the element
                return tolMax;
            else {
                for (unsigned int k = 0; k < (eleOrder + 1); k++)
                    for (unsigned int j = 0; j < (eleOrder + 1); j++)
                        for (unsigned int i = 0; i < (eleOrder + 1); i++) {
                            const double xx = x + i * hx[0];
                            const double yy = y + j * hx[1];
                            const double zz = z + k * hx[2];

                            const Point grid_pp(xx, yy, zz);

                            const double dd0 =
                                (grid_pp - bssn::BSSN_BH_LOC[0]).abs();
                            const double dd1 =
                                (grid_pp - bssn::BSSN_BH_LOC[1]).abs();

                            // std::cout<<"x : "<<x<<" y: "<<y<<" z: "<<z<<"
                            // dd0: "<<dd0<<" dd1: "<<dd1<<" hx:
                            // "<<hx[0]<<std::endl;

                            if (dd0 < R0 || dd1 < R1) {
                                // std::cout<<"x : "<<x<<" y: "<<y<<" z: "<<z<<"
                                // dd0: "<<dd0<<" dd1: "<<dd1<<" hx:
                                // "<<hx[0]<<std::endl;
                                return tolMin;
                            }
                        }

                //@milinda 21/11/2020 - smooth transition of the wtol.
                if (dbh0 < dbh1)
                    return std::min(
                        tolMax,
                        std::max(tolMin, (tolMax - tolMin) / (R_MAX - R0) *
                                                 (dbh0 - R0) +
                                             tolMin));
                else
                    return std::min(
                        tolMax,
                        std::max(tolMin, (tolMax - tolMin) / (R_MAX - R1) *
                                                 (dbh1 - R1) +
                                             tolMin));
            }
        }

#else
        if ((dbh0 > R_MAX) &&
            (dbh1 >
             R_MAX))  // no need to check individual points in the element
            return tolMax;
        else {
            for (unsigned int k = 0; k < (eleOrder + 1); k++)
                for (unsigned int j = 0; j < (eleOrder + 1); j++)
                    for (unsigned int i = 0; i < (eleOrder + 1); i++) {
                        const double xx = x + i * hx[0];
                        const double yy = y + j * hx[1];
                        const double zz = z + k * hx[2];

                        const Point grid_pp(xx, yy, zz);

                        const double dd0 =
                            (grid_pp - bssn::BSSN_BH_LOC[0]).abs();
                        const double dd1 =
                            (grid_pp - bssn::BSSN_BH_LOC[1]).abs();

                        // std::cout<<"x : "<<x<<" y: "<<y<<" z: "<<z<<" dd0:
                        // "<<dd0<<" dd1: "<<dd1<<" hx: "<<hx[0]<<std::endl;

                        if (dd0 < R0 || dd1 < R1) {
                            // std::cout<<"x : "<<x<<" y: "<<y<<" z: "<<z<<"
                            // dd0: "<<dd0<<" dd1: "<<dd1<<" hx:
                            // "<<hx[0]<<std::endl;
                            return tolMin;
                        }
                    }

            //@milinda 21/11/2020 - smooth transition of the wtol.
            if (dbh0 < dbh1)
                return std::min(
                    tolMax, std::max(tolMin, (tolMax - tolMin) / (R_MAX - R0) *
                                                     (dbh0 - R0) +
                                                 tolMin));
            else
                return std::min(
                    tolMax, std::max(tolMin, (tolMax - tolMin) / (R_MAX - R1) *
                                                     (dbh1 - R1) +
                                                 tolMin));
        }
#endif

    } else if (bssn::BSSN_USE_WAVELET_TOL_FUNCTION == 2) {
#ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES
        if (dbh < 0.1) {
            const double R0 =
                (GW::BSSN_GW_RADAII[GW::BSSN_GW_NUM_RADAII - 1] + 10);
            const double R1 =
                (GW::BSSN_GW_RADAII[GW::BSSN_GW_NUM_RADAII - 1] + 20);

            const double tolMin =
                std::min(bssn::BSSN_WAVELET_TOL, bssn::BSSN_GW_REFINE_WTOL);
            const double tolMax = bssn::BSSN_WAVELET_TOL_MAX;

            return std::min(
                tolMax,
                std::max(tolMin,
                         (tolMax - tolMin) / (R1 - R0) * (r - R0) + tolMin));

        } else {
            const double R0     = bssn::BSSN_WAVELET_TOL_FUNCTION_R0;
            const double R1     = bssn::BSSN_WAVELET_TOL_FUNCTION_R1;
            const double tolMax = bssn::BSSN_WAVELET_TOL_MAX;
            const double tolMin = bssn::BSSN_WAVELET_TOL;
            return std::min(
                tolMax,
                std::max(tolMin,
                         (tolMax - tolMin) / (R1 - R0) * (r - R0) + tolMin));
        }
#else
        const double R0     = bssn::BSSN_WAVELET_TOL_FUNCTION_R0;
        const double R1     = bssn::BSSN_WAVELET_TOL_FUNCTION_R1;
        const double tolMax = bssn::BSSN_WAVELET_TOL_MAX;
        const double tolMin = bssn::BSSN_WAVELET_TOL;
        return std::min(
            tolMax, std::max(tolMin, (tolMax - tolMin) / (R1 - R0) * (r - R0) +
                                         tolMin));
#endif

    } else if (bssn::BSSN_USE_WAVELET_TOL_FUNCTION == 3) {
        const double GW_R_SAFETY_FAC = 10.0;
        const double TIME_OFFSET_FAC = 5.0;

        if (T_CURRENT > bssn::BSSN_WAVELET_TOL_FUNCTION_R1 + TIME_OFFSET_FAC) {
            const double RR =
                std::min(T_CURRENT - TIME_OFFSET_FAC,
                         GW::BSSN_GW_RADAII[GW::BSSN_GW_NUM_RADAII - 1] + 10.0);
            const double GW_TOL = bssn::BSSN_GW_REFINE_WTOL;
            const double W_RR   = std::min(GW_TOL, bssn::BSSN_WAVELET_TOL_MAX);
            const double R0     = bssn::BSSN_WAVELET_TOL_FUNCTION_R0;
            const double WTOL_EXP = 10.0;
            // const double R1 = bssn::BSSN_WAVELET_TOL_FUNCTION_R1;
            const double WTOL_EXP_FAC =
                (RR - R0) / std::log10(W_RR / bssn::BSSN_WAVELET_TOL);
            if (r < R0)
                return bssn::BSSN_WAVELET_TOL;
            else
                return std::min(bssn::BSSN_WAVELET_TOL_MAX,
                                ((std::pow(WTOL_EXP, (r - R0) / WTOL_EXP_FAC)) *
                                 bssn::BSSN_WAVELET_TOL));

        } else {
            const double R0       = bssn::BSSN_WAVELET_TOL_FUNCTION_R0;
            const double R1       = bssn::BSSN_WAVELET_TOL_FUNCTION_R1;
            const double GW_TOL   = bssn::BSSN_GW_REFINE_WTOL;
            const double WTOL_EXP = 10.0;
            const double W_RR = std::min(GW_TOL, bssn::BSSN_WAVELET_TOL_MAX);
            const double WTOL_EXP_FAC =
                (R1 - R0) / std::log10(W_RR / bssn::BSSN_WAVELET_TOL);

            if (r < R0)
                return bssn::BSSN_WAVELET_TOL;
            else
                return std::min(bssn::BSSN_WAVELET_TOL_MAX,
                                ((std::pow(WTOL_EXP, (r - R0) / WTOL_EXP_FAC)) *
                                 bssn::BSSN_WAVELET_TOL));
        }

    } else if (bssn::BSSN_USE_WAVELET_TOL_FUNCTION == 4) {
        const double GW_R_SAFETY_FAC = 10.0;
        const double TIME_OFFSET_FAC = 20.0;

        if (T_CURRENT > bssn::BSSN_WAVELET_TOL_FUNCTION_R1 + TIME_OFFSET_FAC) {
            const double RR =
                std::min(T_CURRENT - TIME_OFFSET_FAC,
                         GW::BSSN_GW_RADAII[GW::BSSN_GW_NUM_RADAII - 1] + 10.0);
            const double GW_TOL = bssn::BSSN_GW_REFINE_WTOL;
            const double W_RR   = std::min(GW_TOL, bssn::BSSN_WAVELET_TOL_MAX);
            const double R0     = bssn::BSSN_WAVELET_TOL_FUNCTION_R0;
            const double WTOL_EXP = 10.0;
            // const double R1 = bssn::BSSN_WAVELET_TOL_FUNCTION_R1;
            const double WTOL_EXP_FAC =
                (RR - R0) / std::log10(W_RR / bssn::BSSN_WAVELET_TOL);
            if (r < R0)
                return bssn::BSSN_WAVELET_TOL;
            else
                return std::min(bssn::BSSN_WAVELET_TOL_MAX,
                                ((std::pow(WTOL_EXP, (r - R0) / WTOL_EXP_FAC)) *
                                 bssn::BSSN_WAVELET_TOL));

        } else {
            const double R01      = 3.0 * bssn::BSSN_BH1_MASS;
            const double R02      = 3.0 * bssn::BSSN_BH2_MASS;
            const double R11      = 4.0 * R01;
            const double R12      = 4.0 * R02;
            const double GW_TOL   = bssn::BSSN_GW_REFINE_WTOL;
            const double WTOL_EXP = 10.0;
            const double W_RR = std::min(GW_TOL, bssn::BSSN_WAVELET_TOL_MAX);
            const double W1 =
                (R11 - R01) / std::log10(W_RR / bssn::BSSN_WAVELET_TOL);
            const double W2 =
                (R12 - R02) / std::log10(W_RR / bssn::BSSN_WAVELET_TOL);

            if (dbh0 < R01 || dbh1 < R02)
                return bssn::BSSN_WAVELET_TOL;
            else {
                double minbheps =
                    std::min(((std::pow(WTOL_EXP, (dbh0 - R01) / W1)) *
                              bssn::BSSN_WAVELET_TOL),
                             ((std::pow(WTOL_EXP, (dbh1 - R02) / W2)) *
                              bssn::BSSN_WAVELET_TOL));
                return std::min(bssn::BSSN_WAVELET_TOL_MAX, minbheps);
            }
        }

    } else if (bssn::BSSN_USE_WAVELET_TOL_FUNCTION == 5) {
        Point grid_p(x, y, z);
        const double d1      = (grid_p - bssn::BSSN_BH_LOC[0]).abs();
        const double d2      = (grid_p - bssn::BSSN_BH_LOC[1]).abs();
        const double m1      = bssn::BSSN_BH1_MASS;
        const double m2      = bssn::BSSN_BH2_MASS;
        const double toffset = 16.0;

        const double eps[3]  = {bssn::BSSN_WAVELET_TOL,
                                bssn::BSSN_GW_REFINE_WTOL,
                                bssn::BSSN_WAVELET_TOL_MAX};
        double rad[3];
        rad[0]    = 3.0 * m1;
        rad[1]    = 4.0 * rad[0];
        rad[2]    = GW::BSSN_GW_RADAII[GW::BSSN_GW_NUM_RADAII - 1];

        double e1 = CalTolHelper(T_CURRENT, d1, rad, eps, toffset);

        rad[0]    = 3.0 * m2;
        rad[1]    = 4.0 * rad[0];
        rad[2]    = GW::BSSN_GW_RADAII[GW::BSSN_GW_NUM_RADAII - 1];
        double e2 = CalTolHelper(T_CURRENT, d2, rad, eps, toffset);

        return std::min(e1, e2);

    } else if (bssn::BSSN_USE_WAVELET_TOL_FUNCTION == 6) {
        // WKB Aug 2024
        // use different sensitivities for regions of spacetime which
        // are causally connected to the BHs as cf regions which are
        // spacelike, causally disconnected from the BHs.

        ////////////////////////////////////////////////////////////////
        // set up constants used in this function

        // (max) orbital radius; use strictest refinement here
        const double R_orbit = 8;
        // outermost GW extraction radius
        const double R_min   = GW::BSSN_GW_RADAII[0];
        const double R_max   = GW::BSSN_GW_RADAII[GW::BSSN_GW_NUM_RADAII - 1];

        // expected lapse wave tail length (M) + backreflections
        const double L       = 120;
        // calculate the time after which a given radius's relationship
        // with the grid center is both time-like & clean of lapse noise
        const double t_lim   = std::max(r, (r + L) / std::sqrt(2));

        // wavelet tolerance in acausal (or dirty or unneeded) regions.
        const double eps_disable = bssn::BSSN_WAVELET_TOL_MAX;
        // time to fade from eps_disable to eps_goal
        const double t_fade      = 100;

        ////////////////////////////////////////////////////////////////
        // set up goal resolution to hit in causal clean regions
        // linearly interpolate log tolerances vs log radii

        double eps_goal = eps_disable;  // default value in outer regions
        if (r <= R_orbit) {
            eps_goal = bssn::BSSN_WAVELET_TOL;
        } else if (r <= R_min) {  // log falloff
            // power we're raising the next expression to, scaling out radius
            const double pwr =
                std::log(r / R_orbit) / std::log(R_min / R_orbit);
            // goal wavelet tolerance at end times
            eps_goal =
                bssn::BSSN_WAVELET_TOL *
                std::pow(bssn::BSSN_GW_REFINE_WTOL / bssn::BSSN_WAVELET_TOL,
                         pwr);
        } else if (r <= R_max) {  // plateau
            eps_goal = bssn::BSSN_GW_REFINE_WTOL;
        }

        ////////////////////////////////////////////////////////////////
        // return time-delayed & smoothed wavelet tolerance

        if (T_CURRENT < t_lim) {  // in spacelike or dirty region
            // return max permissible wavelet tolerance
            // effectively disabling / kneecapping WAMR
            return eps_disable;
        } else if (T_CURRENT > t_lim + t_fade) {  // in clean timelike region
            // return standard wavelet tolerance
            return eps_goal;
        } else {  // in transition region
            // return linear transition between log tolerance values
            // slope of transition region
            const double slope = std::log10(eps_goal / eps_disable) / t_fade;
            const double lg_eps =
                std::log10(eps_disable) + slope * (T_CURRENT - t_lim);
            return std::pow(10.0, lg_eps);
        }
    } else {
        // return global wavelet tolerance, irrespective of position
        return bssn::BSSN_WAVELET_TOL;
    }
}

void writeBLockToBinary(const double** unzipVarsRHS, unsigned int offset,
                        const double* pmin, const double* pmax, double* bxMin,
                        double* bxMax, const unsigned int* sz,
                        unsigned int blkSz, double dxFactor,
                        const char* fprefix) {
    const unsigned int nx       = sz[0];
    const unsigned int ny       = sz[1];
    const unsigned int nz       = sz[2];

    const unsigned int ib       = 3;
    const unsigned int jb       = 3;
    const unsigned int kb       = 3;

    const unsigned int ie       = nx - 3;
    const unsigned int je       = ny - 3;
    const unsigned int ke       = nz - 3;

    const unsigned int blkInlSz = (nx - 3) * (ny - 3) * (nz - 3);

    double hx                   = (pmax[0] - pmin[0]) / (nx - 1);
    double hy                   = (pmax[1] - pmin[1]) / (ny - 1);
    double hz                   = (pmax[2] - pmin[2]) / (nz - 1);

    const double dx = (bssn::BSSN_COMPD_MAX[0] - bssn::BSSN_COMPD_MIN[0]) *
                      (1.0 / (double)(1u << bssn::BSSN_MAXDEPTH));
    unsigned int level = bssn::BSSN_MAXDEPTH - ((unsigned int)(hx / dx) - 1);

    MPI_Comm comm      = MPI_COMM_WORLD;

    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);
    // std::cout<<"ranl: "<<rank<<"npes: "<<npes<<std::endl;

    // std::cout<<"nx: "<<nx<<" level: "<<level<<" hx: "<<hx<<" dx:
    // "<<dx<<std::endl;

    if ((hx > (dxFactor * dx)) ||
        (pmin[0] < bxMin[0] || pmin[1] < bxMin[1] || pmin[2] < bxMin[2]) ||
        (pmax[0] > bxMax[0] || pmax[1] > bxMax[1] || pmax[2] > bxMax[2]))
        return;

    double* blkInternal = new double[blkInlSz];
    for (unsigned int var = 0; var < bssn::BSSN_NUM_VARS; var++) {
        char fName[256];
        sprintf(fName, "%s_%s_n_%d_r_%d_p_%d.bin", fprefix,
                bssn::BSSN_VAR_NAMES[var], nx, rank, npes);
        FILE* outfile = fopen(fName, "w");
        if (outfile == NULL) {
            std::cout << fName << " file open failed " << std::endl;
        }

        for (unsigned int k = kb; k < ke; k++)
            for (unsigned int j = jb; j < je; j++)
                for (unsigned int i = ib; i < ie; i++)
                    blkInternal[k * (ny - 3) * (nx - 3) + j * (nx - 3) + i] =
                        unzipVarsRHS[var]
                                    [offset + k * (ny * nx) + j * (ny) + i];

        fwrite(blkInternal, sizeof(double), blkInlSz,
               outfile);  // write out the number of elements.
        fclose(outfile);
    }

    delete[] blkInternal;
}

unsigned int getOctantWeight(const ot::TreeNode* pNode) {
    return (1u << (3 * pNode->getLevel())) * 1;
}

void computeBHLocations(const ot::Mesh* pMesh, const Point* in, Point* out,
                        double** zipVars, double dt) {
    // TODO: this should be easy to adjust based on how many bh's we actually
    // have at some point
    const unsigned int num_bhs              = 2;
    const unsigned int num_data_per_bh      = 3;
    const unsigned int total_points_to_comm = num_bhs * num_data_per_bh;

    dendro::logger::debug("[BH] Computing BH locations");
    MPI_Comm commActive = pMesh->getMPICommunicator();

    Point grid_limits[2];
    Point domain_limits[2];

    grid_limits[0]   = Point(bssn::BSSN_OCTREE_MIN[0], bssn::BSSN_OCTREE_MIN[1],
                             bssn::BSSN_OCTREE_MIN[2]);
    grid_limits[1]   = Point(bssn::BSSN_OCTREE_MAX[0], bssn::BSSN_OCTREE_MAX[1],
                             bssn::BSSN_OCTREE_MAX[2]);

    domain_limits[0] = Point(bssn::BSSN_COMPD_MIN[0], bssn::BSSN_COMPD_MIN[1],
                             bssn::BSSN_COMPD_MIN[2]);
    domain_limits[1] = Point(bssn::BSSN_COMPD_MAX[0], bssn::BSSN_COMPD_MAX[1],
                             bssn::BSSN_COMPD_MAX[2]);

    // gather the points for betas
    std::vector<double> beta0(num_bhs, 0.0);
    std::vector<double> beta1(num_bhs, 0.0);
    std::vector<double> beta2(num_bhs, 0.0);

    std::vector<unsigned int> validIndex_beta0;
    std::vector<unsigned int> validIndex_beta1;
    std::vector<unsigned int> validIndex_beta2;

    // IN is always the number of black holes
    std::vector<double> beta_interleaved(total_points_to_comm, 0.0);
    std::vector<double> bh_pts(total_points_to_comm, 0.0);

    for (unsigned int bhidx = 0; bhidx < num_bhs; ++bhidx) {
        // the points is the xyz direction, and fortunately beta also
        // cooresponds to xyz
        const unsigned int offset = bhidx * 3;
        bh_pts[offset]            = in[bhidx].x();
        bh_pts[offset + 1]        = in[bhidx].y();
        bh_pts[offset + 2]        = in[bhidx].z();
    }

    dendro::logger::debug(
        "[BH] Active processes will now call interpolateToCoords");

    if (pMesh->isActive()) {
        unsigned int activeRank = pMesh->getMPIRank();

        ot::da::interpolateToCoords(
            pMesh, zipVars[VAR::U_BETA0], bh_pts.data(), total_points_to_comm,
            grid_limits, domain_limits, beta0.data(), validIndex_beta0);
        ot::da::interpolateToCoords(
            pMesh, zipVars[VAR::U_BETA1], bh_pts.data(), total_points_to_comm,
            grid_limits, domain_limits, beta1.data(), validIndex_beta1);
        ot::da::interpolateToCoords(
            pMesh, zipVars[VAR::U_BETA2], bh_pts.data(), total_points_to_comm,
            grid_limits, domain_limits, beta2.data(), validIndex_beta2);

        assert(validIndex_beta0.size() == validIndex_beta1.size());
        assert(validIndex_beta1.size() == validIndex_beta2.size());

        for (unsigned int bhidx : validIndex_beta0) {
            // based on the bhidx we get the offset to fix it all up
            const unsigned int offset    = bhidx * 3;
            beta_interleaved[offset]     = beta0[bhidx];
            beta_interleaved[offset + 1] = beta1[bhidx];
            beta_interleaved[offset + 2] = beta2[bhidx];
        }
    }

    dendro::logger::debug(
        "[BH] interpolateToCoords finished, now prepping communication");

    std::vector<double> global_beta_interleaved(total_points_to_comm, 0.0);

#if 0
    dendro::logger::debug(
        "[BH] performing full sum for shift vector allreduce communication "
        "(MPI_Allreduce)...");
    MPI_Allreduce(beta_interleaved.data(), global_beta_interleaved.data(),
                  total_points_to_comm, MPI_DOUBLE, MPI_SUM,
                  pMesh->getMPIGlobalCommunicator());
    dendro::logger::debug("[BH] Allreduce sum of BH shift vectors complete.");
#else
    const int root_rank = 0;
    dendro::logger::debug(
        "[BH] performing full shift vector reduce communication "
        "(MPI_Reduce)...");
    //
    MPI_Reduce(beta_interleaved.data(), global_beta_interleaved.data(),
               total_points_to_comm, MPI_DOUBLE, MPI_SUM, root_rank,
               pMesh->getMPIGlobalCommunicator());

    dendro::logger::debug(
        "[BH] performing broadcast "
        "(MPI_Bcast)...");

    MPI_Bcast(global_beta_interleaved.data(), total_points_to_comm, MPI_DOUBLE,
              root_rank, pMesh->getMPIGlobalCommunicator());

    dendro::logger::debug("[BH] Reduce+Bcast of BH shift vectors complete.");
#endif

    // now final data is available to all processes
    for (unsigned int bh = 0; bh < num_bhs; bh++) {
        const unsigned int offset = bh * 3;
        const double shift_x      = global_beta_interleaved[offset + 0] * dt;
        const double shift_y      = global_beta_interleaved[offset + 1] * dt;
        const double shift_z      = global_beta_interleaved[offset + 2] * dt;

        out[bh] = Point(in[bh].x() - shift_x, in[bh].y() - shift_y,
                        in[bh].z() - shift_z);

        dendro::logger::info("[BH] Black Hole {} new position: [{}, {}, {}]",
                             bh, out[bh].x(), out[bh].y(), out[bh].z());
    }

    dendro::logger::debug("[BH] Finished computing BH locations!");

    return;
}

ot::Mesh* weakScalingReMesh(ot::Mesh* pMesh, unsigned int target_npes) {
    int rank, npes;
    MPI_Comm comm = pMesh->getMPIGlobalCommunicator();
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);

    if (target_npes > npes) {
        if (!rank)
            RAISE_ERROR("target npes "
                        << target_npes
                        << " is larger than global npes:" << npes);

        MPI_Abort(comm, 0);
    }

    const double R_RES_FAC         = 10;
    const double dr1               = bssn::BSSN_BH1_AMR_R / R_RES_FAC;
    const double dr2               = bssn::BSSN_BH2_AMR_R / R_RES_FAC;

    const DendroIntL ELE_SZ_REQ_GL = bssn::BSSN_DENDRO_GRAIN_SZ * target_npes;
    const unsigned int MAX_ITER    = 30;

    unsigned int iter_count        = 0;
    bool is_converged              = true;

    // Sequential split of each level, such that elemental grain size is
    // satiesfied.
    ot::Mesh* current_mesh         = NULL;

    do {
        if (current_mesh == NULL) current_mesh = pMesh;

        unsigned int LMIN, LMAX;
        current_mesh->computeMinMaxLevel(LMIN, LMAX);

        DendroIntL localSz  = current_mesh->getNumLocalMeshElements();
        DendroIntL globalSz = 0;

        par::Mpi_Allreduce(&localSz, &globalSz, 1, MPI_SUM, comm);

        DendroIntL req_l_splits =
            std::max(1ll, (ELE_SZ_REQ_GL - globalSz) / (14 * npes));

        std::vector<unsigned int> ref_flags;
        if (current_mesh->isActive()) {
            ref_flags.resize(current_mesh->getNumLocalMeshElements(),
                             OCT_NO_CHANGE);
            const unsigned int active_rank = current_mesh->getMPIRank();
            const unsigned int active_npes = current_mesh->getMPICommSize();
            const MPI_Comm active_comm     = current_mesh->getMPICommunicator();
            // const DendroIntL lsplit_b =
            // (active_rank)*req_l_splits/active_npes; const DendroIntL lsplit_e
            // = (active_rank+1)*req_l_splits/active_npes; const DendroIntL
            // num_split_rank = lsplit_e-lsplit_b;

            const ot::TreeNode* pNodes = current_mesh->getAllElements().data();

            // DendroIntL lcount =0;
            // for (unsigned int ele = current_mesh->getElementLocalBegin(); ele
            // < current_mesh->getElementLocalEnd(); ele++)
            // {
            //     if(pNodes[ele].getLevel() == LMAX-1)
            //         lcount++;
            // }

            // std::vector<DendroIntL> num_lev_counts;
            // std::vector<DendroIntL> num_lev_offsets;

            // num_lev_counts.resize(active_npes,0);
            // num_lev_offsets.resize(active_npes,0);

            // par::Mpi_Allgather(&lcount,num_lev_counts.data(),1,active_comm);
            // omp_par::scan(num_lev_counts.data(),num_lev_offsets.data(),active_npes);

            // auto it_upper =
            // std::upper_bound(num_lev_offsets.begin(),num_lev_offsets.end(),req_l_splits);
            // unsigned int valid_npes=active_npes;

            // if(it_upper != num_lev_offsets.end())
            //     valid_npes = std::distance(num_lev_offsets.begin(),it_upper)
            //     +1;

            // if (active_rank < valid_npes)
            // {
            //     unsigned int ele_offset =
            //     current_mesh->getElementLocalBegin(); lcount =0; for
            //     (unsigned int ele = current_mesh->getElementLocalBegin(); ele
            //     < current_mesh->getElementLocalEnd(); ele++)
            //     {
            //         if( lcount < req_l_splits  && pNodes[ele].getLevel() ==
            //         LMAX-1)
            //         {
            //             ref_flags[ele-ele_offset] = OCT_SPLIT;
            //             lcount++;
            //         }

            //     }

            // }
            unsigned int ele_offset    = current_mesh->getElementLocalBegin();
            DendroIntL lcount          = 0;
            for (unsigned int ele = current_mesh->getElementLocalBegin();
                 ele < current_mesh->getElementLocalEnd(); ele++) {
                if (lcount < req_l_splits &&
                    pNodes[ele].getLevel() == LMAX - 1) {
                    ref_flags[ele - ele_offset] = OCT_SPLIT;
                    lcount++;
                }
            }
        }

        bool is_refine   = current_mesh->setMeshRefinementFlags(ref_flags);
        bool is_refine_g = false;
        MPI_Allreduce(&is_refine, &is_refine_g, 1, MPI_CXX_BOOL, MPI_LOR, comm);
        if (is_refine_g) {
            ot::Mesh* new_mesh =
                current_mesh->ReMesh(bssn::BSSN_DENDRO_GRAIN_SZ);
            if (current_mesh == pMesh)
                current_mesh = new_mesh;
            else {
                std::swap(current_mesh, new_mesh);
                delete new_mesh;
            }

            localSz = current_mesh->getNumLocalMeshElements();
            par::Mpi_Allreduce(&localSz, &globalSz, 1, MPI_SUM, comm);

            if (!rank)
                std::cout << "weak scaling remesh iter : " << iter_count
                          << " num global elements : " << globalSz << std::endl;

            is_converged =
                (globalSz >=
                 ELE_SZ_REQ_GL);  //|| (fabs(globalSz/(double)npes -
                                  // bssn::BSSN_DENDRO_GRAIN_SZ)/(double)bssn::BSSN_DENDRO_GRAIN_SZ
                                  //< 0.1);//(globalSz >= ELE_SZ_REQ_GL);
            iter_count++;

        } else {
            // no remesh triggered hence, is_converged true;
            is_converged = true;
        }

    } while (!is_converged && iter_count < MAX_ITER);

    return current_mesh;
}

void allocate_bssn_deriv_workspace(const ot::Mesh* pMesh, unsigned int s_fac) {
    deallocate_bssn_deriv_workspace();

    if (!pMesh->isActive()) return;

    // gets the largest block size.
    const std::vector<ot::Block>& blkList = pMesh->getLocalBlockList();
    unsigned int max_blk_sz               = 0;
    for (unsigned int i = 0; i < blkList.size(); i++) {
        unsigned int blk_sz = blkList[i].getAllocationSzX() *
                              blkList[i].getAllocationSzY() *
                              blkList[i].getAllocationSzZ();
        if (blk_sz > max_blk_sz) max_blk_sz = blk_sz;
    }

    if (bssn::BSSN_DERIV_WORKSPACE != nullptr) {
        delete[] bssn::BSSN_DERIV_WORKSPACE;
        bssn::BSSN_DERIV_WORKSPACE = nullptr;
    }

    // One workspace slab per thread so the RHS block loop can be threaded under
    // DENDRO_HYBRID_OMP; each thread indexes WORKSPACE + tid * stride. Without
    // the flag this is a single slab (n_threads = 1) -> identical to before.
    bssn::BSSN_DERIV_WORKSPACE_STRIDE =
        (size_t)s_fac * max_blk_sz * bssn::BSSN_NUM_DERIVS;
#ifdef DENDRO_HYBRID_OMP
    const unsigned int n_threads = (unsigned int)omp_get_max_threads();
#else
    const unsigned int n_threads = 1;
#endif
    // Single source of truth for the threaded-region thread count: every
    // tid-indexed parallel region (RHS, constraints) pins num_threads() to this
    // so omp_get_thread_num() can never index past the slabs allocated here.
    bssn::BSSN_HYBRID_NTHREADS = n_threads;
    bssn::BSSN_DERIV_WORKSPACE =
        new double[(size_t)n_threads * bssn::BSSN_DERIV_WORKSPACE_STRIDE];
#ifdef DENDRO_HYBRID_OMP
    // NUMA first-touch: each thread zeroes its own slab so the pages bind to the
    // thread's local node (first-touch policy) instead of all landing on the
    // master's node -- matters on multi-socket machines.
    if (n_threads > 1) {
#pragma omp parallel num_threads(n_threads)
        {
            const unsigned int tid = (unsigned int)omp_get_thread_num();
            double* slab           = bssn::BSSN_DERIV_WORKSPACE +
                           (size_t)tid * bssn::BSSN_DERIV_WORKSPACE_STRIDE;
            for (size_t i = 0; i < bssn::BSSN_DERIV_WORKSPACE_STRIDE; i++)
                slab[i] = 0.0;
        }
    }

    // Resolve the RHS block-loop schedule (parfile > OMP_SCHEDULE > dynamic,1)
    // now that BSSN_HYBRID_NTHREADS is known. This is the chokepoint every RHS
    // entry point reaches, so the default stays dynamic,1 even for binaries that
    // never touch OMP_SCHEDULE (e.g. the scaling bench). Idempotent on remesh.
    set_rhs_omp_schedule();
#endif
}

void deallocate_bssn_deriv_workspace() {
    if (bssn::BSSN_DERIV_WORKSPACE != nullptr) {
        delete[] bssn::BSSN_DERIV_WORKSPACE;
        bssn::BSSN_DERIV_WORKSPACE = nullptr;
    }
}

void set_rhs_omp_schedule() {
#if defined(DENDRO_HYBRID_OMP) && defined(_OPENMP)
    // Precedence: explicit parfile BSSN_HYBRID_RHS_SCHEDULE ("kind[,chunk]") >
    // OMP_SCHEDULE env > historical default dynamic,1. Values "", "env",
    // "runtime" mean "do not override". A bare "static" (chunk 0) is the
    // NUMA-locality variant: each thread owns a contiguous block range that
    // matches the static first-touch of the unzip buffers (see the RHS NUMA Tax).
    const std::string& spec = bssn::BSSN_HYBRID_RHS_SCHEDULE;

    // "balanced" is the cost-balanced NUMA-aware path (rhs.cpp uses a manual
    // per-thread partition, not an OMP schedule kind) paired with a block-major
    // first-touch of the unzip buffers. Flag the first-touch here -- this runs
    // before the buffers are allocated (main after readParamFile; the scaling
    // bench likewise) so the placement is set on the FIRST allocation. Any other
    // value clears it, keeping the default flat first-touch.
    ot::g_padded_numa_first_touch = (spec == "balanced");
    if (spec == "balanced") {
        // consume uses the manual partition; leave the run-sched ICV at the
        // historical default so any other schedule(runtime) loops are unchanged.
        if (std::getenv("OMP_SCHEDULE") == nullptr)
            omp_set_schedule(omp_sched_dynamic, 1);
        return;
    }
    if (!spec.empty() && spec != "env" && spec != "runtime") {
        std::string kind_s  = spec;
        int chunk           = -1;  // <0 => use the kind's natural default below
        const size_t comma  = spec.find(',');
        if (comma != std::string::npos) {
            kind_s = spec.substr(0, comma);
            try {
                chunk = std::stoi(spec.substr(comma + 1));
            } catch (...) {
                chunk = -1;
            }
        }
        omp_sched_t kind = omp_sched_dynamic;
        bool ok          = true;
        if (kind_s == "static" || kind_s == "STATIC")
            kind = omp_sched_static;
        else if (kind_s == "dynamic" || kind_s == "DYNAMIC")
            kind = omp_sched_dynamic;
        else if (kind_s == "guided" || kind_s == "GUIDED")
            kind = omp_sched_guided;
        else if (kind_s == "auto" || kind_s == "AUTO")
            kind = omp_sched_auto;
        else
            ok = false;
        if (ok) {
            // chunk 0 => implementation default (even contiguous split for
            // static); keep the historical chunk 1 for an unspecified dynamic.
            if (chunk < 0) chunk = (kind == omp_sched_dynamic) ? 1 : 0;
            omp_set_schedule(kind, chunk);
        } else if (std::getenv("OMP_SCHEDULE") == nullptr) {
            omp_set_schedule(omp_sched_dynamic, 1);
        }
    } else if (std::getenv("OMP_SCHEDULE") == nullptr) {
        omp_set_schedule(omp_sched_dynamic, 1);
    }
#endif
}

std::tuple<std::string, std::string, std::string> encode_bh_locs(
    const std::vector<std::pair<Point, Point>>& bh_history,
    const std::vector<double>& bh_times) {
    std::vector<unsigned char> bh1_bytes;
    std::vector<unsigned char> bh2_bytes;

    for (const auto& pair : bh_history) {
        double coords1[3] = {pair.first.x(), pair.first.y(), pair.first.z()};
        double coords2[3] = {pair.second.x(), pair.second.y(), pair.second.z()};

        bh1_bytes.insert(
            bh1_bytes.end(), reinterpret_cast<unsigned char*>(coords1),
            reinterpret_cast<unsigned char*>(coords1) + sizeof(coords1));

        bh2_bytes.insert(
            bh2_bytes.end(), reinterpret_cast<unsigned char*>(coords2),
            reinterpret_cast<unsigned char*>(coords2) + sizeof(coords2));
    }

    // with bytes available, now we can encode with base91_encoding
    std::string bh1_str  = base<91>::encode(std::string(
        reinterpret_cast<const char*>(bh1_bytes.data()), bh1_bytes.size()));

    std::string bh2_str  = base<91>::encode(std::string(
        reinterpret_cast<const char*>(bh2_bytes.data()), bh2_bytes.size()));

    std::string time_str = base<91>::encode(
        std::string(reinterpret_cast<const char*>(bh_times.data()),
                    bh_times.size() * sizeof(double)));

    return std::make_tuple(bh1_str, bh2_str, time_str);
}

std::tuple<std::vector<std::pair<Point, Point>>, std::vector<double>>
decode_bh_locs(const std::string& bh1_str, const std::string& bh2_str,
               const std::string& time_str) {
    const size_t double_size = sizeof(double);
    // decode the strings
    std::string bh1_bytes    = base<91>::decode(bh1_str);
    std::string bh2_bytes    = base<91>::decode(bh2_str);
    std::string time_bytes   = base<91>::decode(time_str);

    const size_t num_entries = time_bytes.size() / double_size;

    // with the bytes back in place, we need to do a reinterpret cast for time
    std::vector<double> time_vector;
    std::vector<std::pair<Point, Point>> bh_locs;
    for (size_t i = 0; i < num_entries; ++i) {
        double value;
        std::memcpy(&value, time_bytes.data() + i * double_size, double_size);
        time_vector.push_back(value);

        double bh_temp[3];

        // create the point for b1
        std::memcpy(&bh_temp, bh1_bytes.data() + i * double_size * 3,
                    double_size * 3);
        Point bh1Pt = Point(bh_temp[0], bh_temp[1], bh_temp[2]);

        // then do it for bh2
        std::memcpy(&bh_temp, bh2_bytes.data() + i * double_size * 3,
                    double_size * 3);
        Point bh2Pt = Point(bh_temp[0], bh_temp[1], bh_temp[2]);

        bh_locs.push_back(std::make_pair(bh1Pt, bh2Pt));
    }

    return std::make_tuple(bh_locs, time_vector);
}

namespace {

// FNV-1a 64. Chosen for being order-sensitive and trivially reproducible across
// builds/compilers -- a reordered E2N map must produce a different digest.
// The offset basis was previously 1469598103934665603 -- 19 digits, one short of
// the real basis (a dropped trailing 7). It still hashed fine, but it meant the
// gate's "did we hash nothing?" guard, which greps for the standard basis
// 0xcbf29ce484222325, could never match. Digest VALUES change with this fix; that
// is harmless (they are only ever compared within one run set) but it voids any
// stored baseline. See the counts emitted by report() for the real non-vacuity check.
constexpr uint64_t FNV_OFFSET = 14695981039346656037ULL;  // 0xcbf29ce484222325
constexpr uint64_t FNV_PRIME  = 1099511628211ULL;         // 0x100000001b3

inline void hash_bytes(uint64_t& h, const void* p, size_t n) {
    const unsigned char* b = static_cast<const unsigned char*>(p);
    for (size_t i = 0; i < n; i++) {
        h ^= static_cast<uint64_t>(b[i]);
        h *= FNV_PRIME;
    }
}

// Fold every rank's local digest into one, IN RANK ORDER, so the result is
// deterministic and also detects data migrating between ranks. XOR/sum would be
// order-insensitive and would miss exactly that.
inline uint64_t combine_ranks(uint64_t local, MPI_Comm comm) {
    int rank, npes;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &npes);
    std::vector<uint64_t> all(rank == 0 ? npes : 0);
    MPI_Gather(&local, 1, MPI_UINT64_T, all.data(), 1, MPI_UINT64_T, 0, comm);
    if (rank) return 0;
    uint64_t h = FNV_OFFSET;
    for (int r = 0; r < npes; r++) hash_bytes(h, &all[r], sizeof(uint64_t));
    return h;
}

// NON-VACUITY, the honest way: report how many items were actually hashed,
// summed over ranks, alongside the digest.
//
// Why not sniff the un-hashed offset out of the printed value (what
// verify_bitexact.sh used to try)? Two reasons, both fatal: (1) the constant it
// grepped for never matched the one this file used (see FNV_OFFSET above), and
// (2) even with that fixed it could not work, because the printed value comes out
// of combine_ranks(), which re-hashes the gathered per-rank digests into a FRESH
// offset -- so the bare basis never appears in the output even when every rank
// hashed nothing. A digest of an empty vector is a perfectly ordinary-looking
// hash; the digest alone can never distinguish "identical" from "both empty".
// Counting the inputs can, so count them.
inline void report(const char* tag, const char* field, uint64_t local,
                   size_t local_n, MPI_Comm comm, int rank) {
    const uint64_t g = combine_ranks(local, comm);
    unsigned long long n_loc = (unsigned long long)local_n, n_glb = 0;
    MPI_Reduce(&n_loc, &n_glb, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, comm);
    if (!rank)
        printf("[fingerprint] %-10s %-8s %016lx %llu\n", tag, field,
               (unsigned long)g, n_glb);
}

}  // namespace

void meshFingerprint(const ot::Mesh* pMesh, const char* tag) {
    // Inactive ranks hold no mesh but MUST still reach the gathers below.
    MPI_Comm comm = pMesh->getMPIGlobalCommunicator();
    int rank;
    MPI_Comm_rank(comm, &rank);
    const bool act = pMesh->isActive();

    uint64_t h_ele = FNV_OFFSET, h_e2e = FNV_OFFSET, h_e2n = FNV_OFFSET,
             h_dg = FNV_OFFSET, h_blk = FNV_OFFSET, h_sm = FNV_OFFSET;
    size_t n_ele = 0, n_e2e = 0, n_e2n = 0, n_dg = 0, n_blk = 0, n_sm = 0;

    if (act) {
        const std::vector<ot::TreeNode>& ele = pMesh->getAllElements();
        n_ele = ele.size();
        for (size_t i = 0; i < ele.size(); i++) {
            const unsigned int c[4] = {ele[i].getX(), ele[i].getY(),
                                       ele[i].getZ(), ele[i].getLevel()};
            hash_bytes(h_ele, c, sizeof(c));
        }

        const std::vector<unsigned int>& e2e = pMesh->getE2EMapping();
        n_e2e = e2e.size();
        if (!e2e.empty())
            hash_bytes(h_e2e, e2e.data(), e2e.size() * sizeof(unsigned int));

        const std::vector<unsigned int>& e2n = pMesh->getE2NMapping();
        n_e2n = e2n.size();
        if (!e2n.empty())
            hash_bytes(h_e2n, e2n.data(), e2n.size() * sizeof(unsigned int));

        const std::vector<unsigned int>& dg = pMesh->getE2NMapping_DG();
        n_dg = dg.size();
        if (!dg.empty())
            hash_bytes(h_dg, dg.data(), dg.size() * sizeof(unsigned int));

        const std::vector<ot::Block>& blk = pMesh->getLocalBlockList();
        n_blk = blk.size();
        for (size_t i = 0; i < blk.size(); i++) {
            const ot::TreeNode bn = blk[i].getBlockNode();
            const unsigned int f[8] = {
                bn.getX(),       bn.getY(),
                bn.getZ(),       bn.getLevel(),
                blk[i].getAllocationSzX(), blk[i].getAllocationSzY(),
                blk[i].getAllocationSzZ(), blk[i].getRegularGridLev()};
            hash_bytes(h_blk, f, sizeof(f));
            const DendroIntL off = blk[i].getOffset();
            hash_bytes(h_blk, &off, sizeof(off));
        }

        // Scatter map: the ghost exchange's wiring. A reordered send/recv list
        // silently permutes ghost data.
        const std::vector<unsigned int>& sSM = pMesh->getSendNodeSM();
        const std::vector<unsigned int>& rSM = pMesh->getRecvNodeSM();
        n_sm = sSM.size() + rSM.size();
        if (!sSM.empty())
            hash_bytes(h_sm, sSM.data(), sSM.size() * sizeof(unsigned int));
        if (!rSM.empty())
            hash_bytes(h_sm, rSM.data(), rSM.size() * sizeof(unsigned int));
    }

    report(tag, "elements", h_ele, n_ele, comm, rank);
    report(tag, "e2e", h_e2e, n_e2e, comm, rank);
    report(tag, "e2n", h_e2n, n_e2n, comm, rank);
    report(tag, "e2n_dg", h_dg, n_dg, comm, rank);
    report(tag, "blocks", h_blk, n_blk, comm, rank);
    report(tag, "scattermap", h_sm, n_sm, comm, rank);

    // MESH THREADING CANARY -- deliberately NOT a [fingerprint] line.
    //
    // It must not join the compared digests: it is EXPECTED to vary with T (that
    // is the whole point), so folding it in would make the gate fail always.
    // It answers the question the digests structurally cannot: "did the Mesh ctor
    // actually thread?" Without it, a build with mesh threading disabled runs
    // identically single-threaded at every T, every digest matches trivially, and
    // the gate reports PASS -- the 9bb8f45 vacuous-pass mode, which was bit-exact
    // and silent for two days.
    //
    // bssn::BSSN_HYBRID_NTHREADS cannot serve here: it is BSSN's RHS/constraints
    // thread count, it is set AFTER a Mesh exists (it sizes off getLocalBlockList),
    // and it says nothing about dendrolib's ctor.
    if (!rank)
        printf("[meshcanary] %-10s omp_ctor_threads=%u\n", tag,
               ot::mesh_ctor_omp_threads);
}

void stateFingerprint(const ot::Mesh* pMesh, const DendroScalar* const* vars,
                      unsigned int nVars, const char* tag) {
    MPI_Comm comm = pMesh->getMPIGlobalCommunicator();
    int rank;
    MPI_Comm_rank(comm, &rank);

    uint64_t h  = FNV_OFFSET;
    size_t n_st = 0;
    if (pMesh->isActive() && vars) {
        // Local nodes only: ghost values are copies of some other rank's
        // local nodes, so hashing them would double-count and make the digest
        // depend on the partition rather than on the solution.
        const unsigned int nb = pMesh->getNodeLocalBegin();
        const unsigned int ne = pMesh->getNodeLocalEnd();
        for (unsigned int v = 0; v < nVars; v++)
            if (vars[v]) {
                hash_bytes(h, vars[v] + nb,
                           (size_t)(ne - nb) * sizeof(DendroScalar));
                n_st += (size_t)(ne - nb);
            }
    }
    report(tag, "state", h, n_st, comm, rank);
}

}  // end of namespace bssn

namespace bssn {

namespace timer {
void initFlops() {
    total_runtime.start();
    t_f2o.start();
    t_cons.start();
    t_cons_unzip.start();
    t_cons_kernel.start();
    t_cons_zipex.start();
    t_cons_deriv.start();
    t_cons_pts.start();
    t_bal.start();
    t_mesh.start();
    t_rkSolve.start();
    t_ghostEx_sync.start();
    t_unzip_sync.start();

    for (unsigned int i = 0; i < NUM_FACES; i++)
        dendro::timer::t_unzip_sync_face[i].start();

    dendro::timer::t_unzip_async_internal.start();
    dendro::timer::t_unzip_sync_edge.start();
    dendro::timer::t_unzip_sync_vtex.start();
    dendro::timer::t_unzip_p2c.start();
    dendro::timer::t_unzip_sync_nodalval.start();
    dendro::timer::t_unzip_sync_cpy.start();
    dendro::timer::t_unzip_sync_f_c1.start();
    dendro::timer::t_unzip_sync_f_c2.start();
    dendro::timer::t_unzip_sync_f_c3.start();

    t_unzip_async.start();
    dendro::timer::t_unzip_async_comm.start();

    dendro::timer::t_unzip_async_internal.start();
    dendro::timer::t_unzip_async_external.start();
    dendro::timer::t_unzip_async_comm.start();
    t_deriv.start();
    t_rhs.start();

    t_rhs_a.start();
    t_rhs_b.start();
    t_rhs_gt.start();
    t_rhs_chi.start();
    t_rhs_At.start();
    t_rhs_K.start();
    t_rhs_Gt.start();
    t_rhs_B.start();

    t_bdyc.start();

    t_zip.start();
    t_rkStep.start();
    t_isReMesh.start();
    t_gridTransfer.start();
    t_ioVtu.start();
    t_ioCheckPoint.start();
}

void resetSnapshot() {
    total_runtime.snapreset();
    t_f2o.snapreset();
    t_cons.snapreset();
    t_cons_unzip.snapreset();
    t_cons_kernel.snapreset();
    t_cons_zipex.snapreset();
    t_cons_deriv.snapreset();
    t_cons_pts.snapreset();
    t_bal.snapreset();
    t_mesh.snapreset();
    t_rkSolve.snapreset();
    t_ghostEx_sync.snapreset();
    t_unzip_sync.snapreset();

    for (unsigned int i = 0; i < NUM_FACES; i++)
        dendro::timer::t_unzip_sync_face[i].snapreset();

    dendro::timer::t_unzip_sync_internal.snapreset();
    dendro::timer::t_unzip_sync_edge.snapreset();
    dendro::timer::t_unzip_sync_vtex.snapreset();
    dendro::timer::t_unzip_p2c.snapreset();
    dendro::timer::t_unzip_sync_nodalval.snapreset();
    dendro::timer::t_unzip_sync_cpy.snapreset();

    dendro::timer::t_unzip_sync_f_c1.snapreset();
    dendro::timer::t_unzip_sync_f_c2.snapreset();
    dendro::timer::t_unzip_sync_f_c3.snapreset();

    t_unzip_async.snapreset();
    dendro::timer::t_unzip_async_internal.snapreset();
    dendro::timer::t_unzip_async_external.snapreset();
    dendro::timer::t_unzip_async_comm.snapreset();

    dendro::timer::t_ghost_pack.snapreset();
    dendro::timer::t_ghost_wait.snapreset();
    dendro::timer::t_ghost_unpack.snapreset();

    t_deriv.snapreset();
    t_rhs.snapreset();
    t_rhs_ko.snapreset();

    t_rhs_a.snapreset();
    t_rhs_b.snapreset();
    t_rhs_gt.snapreset();
    t_rhs_chi.snapreset();
    t_rhs_At.snapreset();
    t_rhs_K.snapreset();
    t_rhs_Gt.snapreset();
    t_rhs_B.snapreset();

    t_bdyc.snapreset();

    t_zip.snapreset();
    t_rkStep.snapreset();
    for (unsigned int s = 0; s < 6; s++) t_rkStage[s].snapreset();
    t_isReMesh.snapreset();
    t_gridTransfer.snapreset();
    t_ioVtu.snapreset();
    t_ioCheckPoint.snapreset();
}

void profileInfo(const char* filePrefix, const ot::Mesh* pMesh) {
    // commActive is only set on active ranks; the reductions below use it.
    if (!pMesh->isActive()) return;
    int activeRank, activeNpes, globalRank, globalNpes;

    MPI_Comm commActive;
    MPI_Comm commGlobal;

    if (pMesh->isActive()) {
        commActive = pMesh->getMPICommunicator();
        activeRank = pMesh->getMPIRank();
        activeNpes = pMesh->getMPICommSize();
    }

    globalRank = pMesh->getMPIRankGlobal();
    globalNpes = pMesh->getMPICommSizeGlobal();
    commGlobal = pMesh->getMPIGlobalCommunicator();

    double t_stat;
    double t_stat_g[3];

    const char separator = ' ';
    const int nameWidth  = 30;
    const int numWidth   = 10;

    char fName[256];
    std::ofstream outfile;

    DendroIntL localSz, globalSz;

    if (!activeRank) {
        sprintf(fName, "%s_final.prof", filePrefix);
        outfile.open(fName);
        if (outfile.fail()) {
            std::cout << fName << " file open failed " << std::endl;
            return;
        }

        outfile << "active npes : " << activeNpes << std::endl;
        outfile << "global npes : " << globalNpes << std::endl;
        outfile << "partition tol : " << bssn::BSSN_LOAD_IMB_TOL << std::endl;
        outfile << "wavelet tol : " << bssn::BSSN_WAVELET_TOL << std::endl;
        outfile << "maxdepth : " << bssn::BSSN_MAXDEPTH << std::endl;
    }

    MPI_Comm comm     = commActive;
    unsigned int rank = activeRank;

    localSz           = pMesh->getNumLocalMeshElements();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "Elements : " << globalSz << std::endl;

    localSz = pMesh->getNumLocalMeshNodes();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "DOG(zip) : " << globalSz << std::endl;

    localSz = pMesh->getDegOfFreedomUnZip();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "DOG(unzip) : " << globalSz << std::endl;

    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "step";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "min(s)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "mean(s)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "max(s)" << std::endl;

    t_stat = total_runtime.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "+runtime(s)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_f2o.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << " ++f2o";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_cons.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << " ++construction";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rkSolve.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << " ++rkSolve";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_bal.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --2:1 balance";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_mesh.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --mesh";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rkStep.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --rkstep";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ghostEx_sync.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --ghostExchge.";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_unzip_sync.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_unzip_async.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++unzip_async";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

#ifdef ENABLE_DENDRO_PROFILE_COUNTERS
    t_stat = dendro::timer::t_unzip_async_internal.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_internal";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_async_external.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_external";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_async_comm.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_comm (comm) ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;
#endif

    t_stat = t_deriv.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --deriv ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_bdyc.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --boundary con ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_zip.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --zip";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ioVtu.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --vtu";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ioCheckPoint.seconds;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --checkpoint";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    if (!rank) outfile.close();
}

void profileInfoIntermediate(const char* filePrefix, const ot::Mesh* pMesh,
                             const unsigned int currentStep) {
    // commActive is only set on active ranks; the reductions below use it.
    if (!pMesh->isActive()) return;
    int activeRank, activeNpes, globalRank, globalNpes;

    MPI_Comm commActive;
    MPI_Comm commGlobal;

    if (pMesh->isActive()) {
        commActive = pMesh->getMPICommunicator();
        activeRank = pMesh->getMPIRank();
        activeNpes = pMesh->getMPICommSize();
    }

    globalRank = pMesh->getMPIRankGlobal();
    globalNpes = pMesh->getMPICommSizeGlobal();
    commGlobal = pMesh->getMPIGlobalCommunicator();

    double t_stat;
    double t_stat_g[3];

    const char separator = ' ';
    const int nameWidth  = 30;
    const int numWidth   = 10;

    char fName[256];
    std::ofstream outfile;

    DendroIntL localSz, globalSz;

    DendroIntL ghostElements;
    DendroIntL localElements;

    DendroIntL ghostNodes;
    DendroIntL localNodes;

    DendroIntL totalSendNode;
    DendroIntL totalRecvNode;

    DendroIntL numCalls;

#ifdef BSSN_PROFILE_HUMAN_READABLE
    if (!activeRank) {
        sprintf(fName, "%s_im.prof", filePrefix);
        outfile.open(fName, std::fstream::app);
        if (outfile.fail()) {
            std::cout << fName << " file open failed " << std::endl;
            return;
        }

        outfile << "active npes : " << activeNpes << std::endl;
        outfile << "global npes : " << globalNpes << std::endl;
        outfile << "current step : " << currentStep << std::endl;
        outfile << "partition tol : " << bssn::BSSN_LOAD_IMB_TOL << std::endl;
        outfile << "wavelet tol : " << bssn::BSSN_WAVELET_TOL << std::endl;
        outfile << "maxdepth : " << bssn::BSSN_MAXDEPTH << std::endl;
    }

    MPI_Comm comm     = commActive;
    unsigned int rank = activeRank;

    localSz           = pMesh->getNumLocalMeshElements();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "Elements : " << globalSz << std::endl;

    localSz = pMesh->getNumLocalMeshNodes();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "DOG(zip) : " << globalSz << std::endl;

    localSz = pMesh->getDegOfFreedomUnZip();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << "DOG(unzip) : " << globalSz << std::endl;

    ghostElements =
        pMesh->getNumPreGhostElements() + pMesh->getNumPostGhostElements();
    localElements = pMesh->getNumLocalMeshElements();

    ghostNodes    = pMesh->getNumPreMeshNodes() + pMesh->getNumPostMeshNodes();
    localNodes    = pMesh->getNumLocalMeshNodes();

    if (!rank)
        outfile << "========================= MESH "
                   "==========================================================="
                   "============ "
                << std::endl;

    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "step";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "min(#)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "mean(#)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "max(#)" << std::endl;

    t_stat = ghostElements;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "ghost Elements";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = localElements;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "local Elements";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = ghostNodes;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "ghost Nodes";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = localNodes;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "local Nodes";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = pMesh->getGhostExcgTotalSendNodeCount();
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "send Nodes";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = pMesh->getGhostExcgTotalRecvNodeCount();
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "recv Nodes";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    if (!rank)
        outfile << "========================= RUNTIME "
                   "==========================================================="
                   "======== "
                << std::endl;
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "step";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "min(s)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "mean(s)";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "max(s)" << std::endl;

    /* t_stat=total_runtime.seconds;
    computeOverallStats(&t_stat,t_stat_g,comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<"+runtime(s)"; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


    t_stat=t_f2o.seconds;
    computeOverallStats(&t_stat,t_stat_g,comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<" ++f2o"; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


    t_stat=t_cons.seconds;
    computeOverallStats(&t_stat,t_stat_g,comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<" ++construction"; if(!rank)outfile << std::left
    << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


    t_stat=t_rkSolve.seconds;
    computeOverallStats(&t_stat,t_stat_g,comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<" ++rkSolve"; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;*/

    t_stat = t_bal.snap;
    // numCalls=t_bal.num_calls;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++2:1 balance";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_mesh.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++mesh";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rkStep.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++rkstep";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ghostEx_sync.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++ghostExchge.";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_unzip_sync.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++unzip_sync";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_unzip_async.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++unzip_async";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

#ifdef ENABLE_DENDRO_PROFILE_COUNTERS

    t_stat = dendro::timer::t_unzip_async_comm.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_comm_wait (comm) ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_nodalval.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_nodalVal";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_f_c1.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --t_unzip_sync_f_c1";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_f_c2.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --t_unzip_sync_f_c2";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_f_c3.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --t_unzip_sync_f_c3";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_cpy.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --t_unzip_sync_cpy";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_internal.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_internal";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[0].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_left";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[1].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_right";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[2].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_down";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[3].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_up";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[4].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_back";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_face[5].snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_face_front";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_edge.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_edge";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_sync_vtex.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_sync_vtex";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = dendro::timer::t_unzip_p2c.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --unzip_p2c";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;
#endif

    /*
    #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
    t_stat=t_unzip_async_internal.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<"  --unzip_internal"; if(!rank)outfile <<
    std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

    t_stat=t_unzip_async_external.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<"  --unzip_external"; if(!rank)outfile <<
    std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


    t_stat=t_unzip_async_comm.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator) <<"  --unzip_comm (comm) "; if(!rank)outfile <<
    std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
    if(!rank)outfile << std::left << std::setw(nameWidth) <<
    std::setfill(separator)<<t_stat_g[1]; if(!rank)outfile << std::left <<
    std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;
    #endif
    */
    t_stat = t_isReMesh.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++isReMesh";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_gridTransfer.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++gridTransfer";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_deriv.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++deriv ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++compute_rhs ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_a.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_a ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_b.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_b ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_gt.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_gt ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_chi.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_chi ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_At.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_At ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_K.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_K ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_Gt.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_Gt ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_rhs_B.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  --compute_rhs_B ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_bdyc.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++boundary con ";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_zip.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++zip";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ioVtu.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++vtu";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    t_stat = t_ioCheckPoint.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << "  ++checkpoint";
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[0];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[1];
    if (!rank)
        outfile << std::left << std::setw(nameWidth) << std::setfill(separator)
                << t_stat_g[2] << std::endl;

    if (!rank) outfile.close();
#else

    if (!activeRank) {
        sprintf(fName, "%s_im.prof", filePrefix);
        outfile.open(fName, std::fstream::app);
        if (outfile.fail()) {
            std::cout << fName << " file open failed " << std::endl;
            return;
        }

        // writes the header
        if (currentStep == 0)
            outfile << "step\t act_npes\t glb_npes\t part_tol\t wave_tol\t "
                       "maxdepth\t numOcts\t dof_zip\t dof_unzip\t"
                    << "element_ghost_min\t element_ghost_mean\t "
                       "element_ghost_max\t"
                    << "element_local_min\t element_local_mean\t "
                       "element_local_max\t"
                    << "nodes_local_min\t nodes_local_mean\t nodes_local|max\t"
                    << "send_nodes_min\t send_nodes_mean\t send_nodes_max\t"
                    << "recv_nodes_min\t recv_nodes_mean\t recv_nodes_max\t"
                    << "bal_min\t bal_mean\t bal_max\t"
                    << "mesh_min\t mesh_mean\t mesh_max\t"
                    << "rkstep_min\t rkstep_mean\t rkstep_max\t"
                    << "ghostEx_min\t ghostEx_mean\t ghostEx_max\t"
                    << "unzip_sync_min\t unzip_sync_mean\t unzip_sync_max\t"
                    << "unzip_async_min\t unzip_async_mean\t unzip_async_max\t"
                    << "unzip_async_wait_min\t unzip_async_wait_mean\t "
                       "unzip_async_wait_max\t"
                    << "isRemesh_min\t isRemesh_mean\t isRemesh_max\t"
                    << "GT_min\t GT_mean\t GT_max\t"
                    << "deriv_min\t deriv_mean\t deriv_max\t"
                    << "rhs_min\t rhs_mean\t rhs_max\t" << std::endl;
    }

    MPI_Comm comm     = commActive;
    unsigned int rank = activeRank;

    if (!rank) outfile << currentStep << "\t ";
    if (!rank) outfile << activeNpes << "\t ";
    if (!rank) outfile << globalNpes << "\t ";
    if (!rank) outfile << bssn::BSSN_LOAD_IMB_TOL << "\t ";
    if (!rank) outfile << bssn::BSSN_WAVELET_TOL << "\t ";
    if (!rank) outfile << bssn::BSSN_MAXDEPTH << "\t ";

    localSz = pMesh->getNumLocalMeshElements();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << globalSz << "\t ";

    localSz = pMesh->getNumLocalMeshNodes();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << globalSz << "\t ";

    localSz = pMesh->getDegOfFreedomUnZip();
    par::Mpi_Reduce(&localSz, &globalSz, 1, MPI_SUM, 0, comm);
    if (!rank) outfile << globalSz << "\t ";

    ghostElements =
        pMesh->getNumPreGhostElements() + pMesh->getNumPostGhostElements();
    localElements = pMesh->getNumLocalMeshElements();

    t_stat        = ghostElements;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = localElements;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    ghostNodes = pMesh->getNumPreMeshNodes() + pMesh->getNumPostMeshNodes();
    localNodes = pMesh->getNumLocalMeshNodes();

    /*t_stat=ghostNodes;
    computeOverallStats(&t_stat,t_stat_g,comm);
    if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t
    ";*/

    t_stat     = localNodes;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = pMesh->getGhostExcgTotalSendNodeCount();
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = pMesh->getGhostExcgTotalRecvNodeCount();
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_bal.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_mesh.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_rkStep.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_ghostEx_sync.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_unzip_sync.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_unzip_async.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = dendro::timer::t_unzip_async_comm.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_isReMesh.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_gridTransfer.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_deriv.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    t_stat = t_rhs.snap;
    computeOverallStats(&t_stat, t_stat_g, comm);
    if (!rank)
        outfile << t_stat_g[0] << "\t " << t_stat_g[1] << "\t " << t_stat_g[2]
                << "\t ";

    if (!rank) outfile << std::endl;
    if (!rank) outfile.close();
#endif
}

// JSONL profile emitter: one self-contained record per call to
// <prefix>_steps.jsonl. All ranks must call (MPI reduction inside); only
// active rank 0 writes. Consumer: plot_profile.py.
//
// ghost-comm is derivable as (unzip_wcomm - unzip) in the plotter.
#define BSSN_PROFILE_JSONL_SCHEMA_VERSION 2

namespace {
// All ranks must call; only rank 0 receives a non-null `os`.
inline void emit_timer_triplet(std::ostream* os, const char* key,
                               long double snap_seconds, MPI_Comm comm) {
    double t      = static_cast<double>(snap_seconds);
    double tg[3]  = {0, 0, 0};
    computeOverallStats(&t, tg, comm);
    if (os) {
        (*os) << "\"" << key << "\":{\"min\":" << tg[0]
              << ",\"mean\":" << tg[1] << ",\"max\":" << tg[2] << "}";
    }
}
}  // namespace

void profileInfoJSON(const char* filePrefix, const ot::Mesh* pMesh,
                     const unsigned int currentStep,
                     const std::vector<profiler_t>* ets_ctxpt,
                     const std::vector<profiler_t>* app_ctxpt) {
    if (!pMesh->isActive()) return;

    MPI_Comm comm     = pMesh->getMPICommunicator();
    int activeRank    = pMesh->getMPIRank();
    int activeNpes    = pMesh->getMPICommSize();
    int globalRank    = pMesh->getMPIRankGlobal();
    int globalNpes    = pMesh->getMPICommSizeGlobal();

    // File only opened on rank 0.
    std::ofstream out;
    std::ostream* os = nullptr;
    if (!activeRank) {
        char fName[256];
        std::snprintf(fName, sizeof(fName), "%s_steps.jsonl", filePrefix);
        out.open(fName, std::fstream::app);
        if (out.fail()) {
            std::cout << fName << " file open failed " << std::endl;
            // proceed without writing; we still need all ranks to do reductions
        } else {
            os = &out;
            out << std::scientific << std::setprecision(6);
        }
    }

    // ---- header / mesh ---------------------------------------------------
    DendroIntL localSz, globalSzElem = 0, globalSzZip = 0, globalSzUnzip = 0,
                        globalSzBlk = 0;
    localSz = pMesh->getNumLocalMeshElements();
    par::Mpi_Reduce(&localSz, &globalSzElem, 1, MPI_SUM, 0, comm);
    localSz = pMesh->getNumLocalMeshNodes();
    par::Mpi_Reduce(&localSz, &globalSzZip, 1, MPI_SUM, 0, comm);
    localSz = pMesh->getDegOfFreedomUnZip();
    par::Mpi_Reduce(&localSz, &globalSzUnzip, 1, MPI_SUM, 0, comm);
    localSz = (DendroIntL)pMesh->getLocalBlockList().size();
    par::Mpi_Reduce(&localSz, &globalSzBlk, 1, MPI_SUM, 0, comm);

    // min/max refinement level across ranks. Reduce over the ACTIVE comm only
    // (this routine already early-returns for inactive ranks) -- do NOT call
    // Mesh::computeMinMaxLevel here: it ends with a global Bcast that inactive
    // ranks would never reach, deadlocking the solver.
    unsigned int lmin_g = 0, lmax_g = 0;
    {
        const std::vector<ot::TreeNode>& allEle = pMesh->getAllElements();
        const unsigned int eb = pMesh->getElementLocalBegin();
        const unsigned int ee = pMesh->getElementLocalEnd();
        unsigned int lmin_l = (ee > eb) ? allEle[eb].getLevel() : 0xFFFFFFFFu;
        unsigned int lmax_l = (ee > eb) ? allEle[eb].getLevel() : 0u;
        for (unsigned int e = eb + 1; e < ee; e++) {
            const unsigned int lv = allEle[e].getLevel();
            if (lv < lmin_l) lmin_l = lv;
            if (lv > lmax_l) lmax_l = lv;
        }
        par::Mpi_Reduce(&lmin_l, &lmin_g, 1, MPI_MIN, 0, comm);
        par::Mpi_Reduce(&lmax_l, &lmax_g, 1, MPI_MAX, 0, comm);
    }

    if (os) {
        (*os) << "{"
              << "\"schema_version\":" << BSSN_PROFILE_JSONL_SCHEMA_VERSION
              << ",\"step\":" << currentStep
              << ",\"active_npes\":" << activeNpes
              << ",\"global_npes\":" << globalNpes
              << ",\"part_tol\":" << bssn::BSSN_LOAD_IMB_TOL
              << ",\"wavelet_tol\":" << bssn::BSSN_WAVELET_TOL
              << ",\"maxdepth\":" << bssn::BSSN_MAXDEPTH
              << ",\"num_elements\":" << globalSzElem
              << ",\"num_blocks\":" << globalSzBlk
              << ",\"num_zip_dof\":" << globalSzZip
              << ",\"num_unzip_dof\":" << globalSzUnzip
              << ",\"lmin\":" << lmin_g
              << ",\"lmax\":" << lmax_g
              << ",\"ele_order\":" << bssn::BSSN_ELE_ORDER
              << ",\"omp_threads\":" << bssn::BSSN_HYBRID_NTHREADS;
    }

    // Mirrors dendrolib_upstream/ODE/include/{ctx,ets}.h. If you reorder
    // either enum there, update here.
    enum CtxIdx {
        CTX_IS_REMESH    = 0,
        CTX_REMESH       = 1,
        CTX_GRID_TRASFER = 2,
        CTX_RHS          = 3,
        CTX_UNZIP_WCOMM  = 5,
        CTX_UNZIP        = 6,
        CTX_ZIP_WCOMM    = 7,
        CTX_ZIP          = 8,
    };
    enum EtsIdx {
        ETS_EVOLVE   = 0,
        ETS_STAGE_0  = 1,
    };

    // Legacy bssn::timer hooks don't tick under ETS; pick dendrolib values
    // when supplied.
    auto pick_app = [&](int ctx_idx, long double fallback) -> long double {
        return (app_ctxpt && ctx_idx < (int)app_ctxpt->size())
                   ? (*app_ctxpt)[ctx_idx].snap
                   : fallback;
    };
    auto pick_ets = [&](int ets_idx, long double fallback) -> long double {
        return (ets_ctxpt && ets_idx < (int)ets_ctxpt->size())
                   ? (*ets_ctxpt)[ets_idx].snap
                   : fallback;
    };

    // ---- phase block -----------------------------------------------------
    if (os) (*os) << ",\"phase\":{";
    bool first = true;
#define BSSN_JSONL_PHASE_RAW(name, snap_value) \
    do {                                                                  \
        if (os && !first) (*os) << ",";                                   \
        first = false;                                                    \
        emit_timer_triplet(os, name, snap_value, comm);                   \
    } while (0)

    BSSN_JSONL_PHASE_RAW("balance",       t_bal.snap);
    BSSN_JSONL_PHASE_RAW("mesh",          t_mesh.snap);
    BSSN_JSONL_PHASE_RAW("rk_step",
                         pick_ets(ETS_EVOLVE, t_rkStep.snap));
    // unzip = local interp (CTX::UNZIP); unzip_wcomm = the WHOLE Ctx::unzip,
    // i.e. unzip + pack + post + wait + unpack (ctx.h:664-796, UNZIP nested
    // inside). So unzip_wcomm - unzip is the ghost exchange, and the three
    // timers below break THAT down. Do not call unzip_wcomm a comm timer.
    BSSN_JSONL_PHASE_RAW("unzip",
                         pick_app(CTX_UNZIP, t_unzip_sync.snap));
    BSSN_JSONL_PHASE_RAW("unzip_wcomm",
                         pick_app(CTX_UNZIP_WCOMM, t_unzip_sync.snap));
    // Ghost exchange split (rank wall; started outside any parallel region).
    // ghost_pack/ghost_unpack are threaded over NEIGHBOUR RANKS, not nodes, so
    // their parallelism is capped by neighbour count -- raising T gives them
    // fewer ranks to spread over more threads while each rank packs MORE bytes
    // (per-rank ghost ~ (V/R)^(2/3) rises ~2.5x at T=4 even though the job's
    // aggregate ghost volume falls as R^(1/3)).
    // ghost_wait is an UPPER BOUND on wire time, never a measurement: it also
    // absorbs neighbour load imbalance.
    BSSN_JSONL_PHASE_RAW("ghost_pack",    dendro::timer::t_ghost_pack.snap);
    BSSN_JSONL_PHASE_RAW("ghost_wait",    dendro::timer::t_ghost_wait.snap);
    BSSN_JSONL_PHASE_RAW("ghost_unpack",  dendro::timer::t_ghost_unpack.snap);
    // rhs_wall: the whole RHS phase, rank wall time. CTX_RHS brackets bssnRHS()
    // from OUTSIDE the omp region (bssnCtx.cpp:200-218), so it covers every
    // block on every thread. Quote THIS for any cross-config comparison.
    //
    // The _t0 timers below are the sub-phase breakdown. They are started inside
    // the threaded block loop (rhs.cpp:53), and profiler_t::start/stop
    // early-return on worker threads (dendrolib/src/profiler.cpp:28,34) to avoid
    // racing the shared counter -- so each accumulates only THREAD 0's elapsed
    // time in that sub-phase. Because thread 0 runs concurrently with its peers,
    // that is ~the rank's wall time for the sub-phase WHEN THE THREADS ARE
    // BALANCED; under imbalance thread 0 need not be the critical path, so they
    // understate. The _t0 suffix marks that caveat -- it does not mean "divide
    // by T".
    //
    // Two traps, both measured (R=2, depth 9, gcc, 2026-07-16):
    //   1. t_rhs EXCLUDES t_deriv -- it starts AFTER the deriv calls
    //      (rhs.cpp:189-205, then :221). Hence rhs_eqn_t0, not rhs_t0. Treating
    //      the old "rhs" key as the whole RHS understates it by ~2x and makes
    //      any rk_step - rhs - uwcomm subtraction meaningless.
    //   2. t_rhs_ko is a SUBSET of t_rhs (profile_params.h:31), not a sibling.
    //      Adding it double-counts by ~12%.
    // The identity that does hold: deriv_t0 + rhs_eqn_t0 + bdyc_t0 == rhs_wall
    // (matched to <0.2% at T=1 and T=2). Use it to check these numbers.
    BSSN_JSONL_PHASE_RAW("rhs_wall",      pick_app(CTX_RHS, t_rhs.snap));
    BSSN_JSONL_PHASE_RAW("deriv_t0",      t_deriv.snap);
    BSSN_JSONL_PHASE_RAW("rhs_eqn_t0",    t_rhs.snap);
    BSSN_JSONL_PHASE_RAW("rhs_ko_t0",     t_rhs_ko.snap);
    BSSN_JSONL_PHASE_RAW("bdyc_t0",       t_bdyc.snap);
    // Constraint computation (unzip + block loop + zip); feeds GW extraction.
    BSSN_JSONL_PHASE_RAW("constraints",   t_cons.snap);
    BSSN_JSONL_PHASE_RAW("cons_unzip",    t_cons_unzip.snap);
    BSSN_JSONL_PHASE_RAW("cons_kernel",   t_cons_kernel.snap);
    BSSN_JSONL_PHASE_RAW("cons_zipex",    t_cons_zipex.snap);
    // Thread-0 samples; see profile_params.h.
    BSSN_JSONL_PHASE_RAW("cons_deriv_t0",  t_cons_deriv.snap);
    BSSN_JSONL_PHASE_RAW("cons_pts_t0",    t_cons_pts.snap);
    BSSN_JSONL_PHASE_RAW("zip",
                         pick_app(CTX_ZIP, t_zip.snap));
    BSSN_JSONL_PHASE_RAW("is_remesh",
                         pick_app(CTX_IS_REMESH, t_isReMesh.snap));
    BSSN_JSONL_PHASE_RAW("grid_transfer",
                         pick_app(CTX_GRID_TRASFER, t_gridTransfer.snap));
    BSSN_JSONL_PHASE_RAW("io_vtu",        t_ioVtu.snap);
    BSSN_JSONL_PHASE_RAW("io_checkpoint", t_ioCheckPoint.snap);
#undef BSSN_JSONL_PHASE_RAW
    if (os) (*os) << "}";

    // ---- unzip sub-phase breakdown --------------------------------------
    // dendrolib mesh.tcc ticks these inside Mesh::unzip(). The "sync"
    // naming is historical -- they tick under both ghost-exchange paths.
    {
        if (os) (*os) << ",\"unzip_breakdown\":{";
        bool first_ub = true;
#define BSSN_JSONL_UB(name, snap_value) \
    do {                                                                  \
        if (os && !first_ub) (*os) << ",";                                \
        first_ub = false;                                                 \
        emit_timer_triplet(os, name, snap_value, comm);                   \
    } while (0)

        long double face_sum = 0;
        for (unsigned int i = 0; i < NUM_FACES; i++)
            face_sum += dendro::timer::t_unzip_sync_face[i].snap;

        BSSN_JSONL_UB("internal", dendro::timer::t_unzip_sync_internal.snap);
        BSSN_JSONL_UB("p2c",      dendro::timer::t_unzip_p2c.snap);
        BSSN_JSONL_UB("face",     face_sum);
        BSSN_JSONL_UB("edge",     dendro::timer::t_unzip_sync_edge.snap);
        BSSN_JSONL_UB("vtex",     dendro::timer::t_unzip_sync_vtex.snap);
        BSSN_JSONL_UB("nodalval", dendro::timer::t_unzip_sync_nodalval.snap);
        BSSN_JSONL_UB("cpy",      dendro::timer::t_unzip_sync_cpy.snap);
        BSSN_JSONL_UB("async_internal",
                      dendro::timer::t_unzip_async_internal.snap);
        BSSN_JSONL_UB("async_external",
                      dendro::timer::t_unzip_async_external.snap);
        BSSN_JSONL_UB("async_comm",
                      dendro::timer::t_unzip_async_comm.snap);
#undef BSSN_JSONL_UB
        if (os) (*os) << "}";
    }

    // ---- per-RK-stage block ---------------------------------------------
    if (os) (*os) << ",\"rk_stage\":[";
    for (unsigned int s = 0; s < 6; s++) {
        if (os && s) (*os) << ",";
        const long double stage_snap =
            pick_ets(ETS_STAGE_0 + s, t_rkStage[s].snap);
        double t      = static_cast<double>(stage_snap);
        double tg[3]  = {0, 0, 0};
        computeOverallStats(&t, tg, comm);
        if (os) {
            (*os) << "{\"min\":" << tg[0] << ",\"mean\":" << tg[1]
                  << ",\"max\":" << tg[2] << "}";
        }
    }
    if (os) (*os) << "]";

    // ---- per-variable BC block (interior RHS is fused; see rhs.cpp) -----
    if (os) (*os) << ",\"rhs_var\":{";
    first = true;
#define BSSN_JSONL_VAR(name, timer_expr) \
    do {                                                                  \
        if (os && !first) (*os) << ",";                                   \
        first = false;                                                    \
        emit_timer_triplet(os, name, (timer_expr).snap, comm);             \
    } while (0)
    BSSN_JSONL_VAR("a",   t_rhs_a);
    BSSN_JSONL_VAR("b",   t_rhs_b);
    BSSN_JSONL_VAR("gt",  t_rhs_gt);
    BSSN_JSONL_VAR("chi", t_rhs_chi);
    BSSN_JSONL_VAR("At",  t_rhs_At);
    BSSN_JSONL_VAR("K",   t_rhs_K);
    BSSN_JSONL_VAR("Gt",  t_rhs_Gt);
    BSSN_JSONL_VAR("B",   t_rhs_B);
#undef BSSN_JSONL_VAR
    if (os) (*os) << "}";

    if (os) {
        (*os) << "}\n";
        out.close();
    }
}

}  // namespace timer

}  // namespace bssn

namespace GW {
void psi4ShpereDump(const ot::Mesh* mesh, DendroScalar** cVar,
                    unsigned int timestep, double time, double dtheta,
                    double dphi) {
    unsigned int rankGlobal   = mesh->getMPIRankGlobal();
    unsigned int npesGlobal   = mesh->getMPICommSizeGlobal();
    MPI_Comm commGlobal       = mesh->getMPIGlobalCommunicator();

    const unsigned int nTheta = (M_PI) / dtheta;
    const unsigned int nPhi   = (2 * M_PI) / dphi;
    const unsigned int numPts = nTheta * nPhi;

    unsigned int totalModes   = 0;
    for (unsigned int l = 0; l < BSSN_GW_NUM_LMODES; l++)
        totalModes += 2 * BSSN_GW_L_MODES[l] + 1;

    const unsigned int TOTAL_MODES = totalModes;

    DendroComplex* swsh_coeff =
        new DendroComplex[BSSN_GW_NUM_RADAII * TOTAL_MODES];
    DendroComplex* swsh_coeff_g =
        new DendroComplex[BSSN_GW_NUM_RADAII * TOTAL_MODES];

    std::vector<unsigned int> lmCounts;
    std::vector<unsigned int> lmOffset;

    lmCounts.resize(BSSN_GW_NUM_LMODES);
    lmOffset.resize(BSSN_GW_NUM_LMODES);

    for (unsigned int l = 0; l < BSSN_GW_NUM_LMODES; l++)
        lmCounts[l] = 2 * BSSN_GW_L_MODES[l] + 1;

    lmOffset[0] = 0;
    omp_par::scan(&(*(lmCounts.begin())), &(*(lmOffset.begin())),
                  BSSN_GW_NUM_LMODES);

    if (mesh->isActive()) {
        const unsigned int rankActive = mesh->getMPIRank();
        const unsigned int npesActive = mesh->getMPICommSize();

        std::vector<double> coords;
        coords.reserve(3 * numPts);

        std::vector<double> psi4_real;
        psi4_real.resize(numPts);

        std::vector<double> psi4_imag;
        psi4_imag.resize(numPts);

        Point grid_limits[2];
        Point domain_limits[2];

        grid_limits[0] =
            Point(bssn::BSSN_OCTREE_MIN[0], bssn::BSSN_OCTREE_MIN[1],
                  bssn::BSSN_OCTREE_MIN[2]);
        grid_limits[1] =
            Point(bssn::BSSN_OCTREE_MAX[0], bssn::BSSN_OCTREE_MAX[1],
                  bssn::BSSN_OCTREE_MAX[2]);

        domain_limits[0] =
            Point(bssn::BSSN_COMPD_MIN[0], bssn::BSSN_COMPD_MIN[1],
                  bssn::BSSN_COMPD_MIN[2]);
        domain_limits[1] =
            Point(bssn::BSSN_COMPD_MAX[0], bssn::BSSN_COMPD_MAX[1],
                  bssn::BSSN_COMPD_MAX[2]);

        std::vector<unsigned int> validIndex;

        for (unsigned int k = 0; k < BSSN_GW_NUM_RADAII; k++) {
            for (unsigned int i = 0; i < nTheta; i++)
                for (unsigned int j = 0; j < nPhi; j++) {
                    double x =
                        BSSN_GW_RADAII[k] * sin(j * dtheta) * cos(i * dphi);
                    double y =
                        BSSN_GW_RADAII[k] * sin(j * dtheta) * sin(i * dphi);
                    double z = BSSN_GW_RADAII[k] * cos(j * dtheta);

                    coords.push_back(x);
                    coords.push_back(y);
                    coords.push_back(z);
                }

            validIndex.clear();
            ot::da::interpolateToCoords(
                mesh, cVar[bssn::VAR_CONSTRAINT::C_PSI4_REAL],
                &(*(coords.begin())), coords.size(), grid_limits, domain_limits,
                &(*(psi4_real.begin())), validIndex);

            validIndex.clear();
            ot::da::interpolateToCoords(
                mesh, cVar[bssn::VAR_CONSTRAINT::C_PSI4_IMG],
                &(*(coords.begin())), coords.size(), grid_limits, domain_limits,
                &(*(psi4_imag.begin())), validIndex);
        }
    }
}

}  // namespace GW

