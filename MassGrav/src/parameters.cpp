//
// Created by milinda on 8/23/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief
*/
//

#include "parameters.h"

namespace massgrav
{

    unsigned int MASSGRAV_ELE_ORDER =4;
    unsigned int MASSGRAV_IO_OUTPUT_FREQ=10;
    unsigned int MASSGRAV_TIME_STEP_OUTPUT_FREQ=10;

    unsigned int MASSGRAV_GW_EXTRACT_FREQ=std::max(1u,MASSGRAV_IO_OUTPUT_FREQ>>1u);
    
    unsigned int MASSGRAV_REMESH_TEST_FREQ=10;
    double MASSGRAV_IO_OUTPUT_GAP=1.0;

    unsigned int MASSGRAV_USE_WAVELET_TOL_FUNCTION=0;
    double MASSGRAV_WAVELET_TOL=0.0001;
    double MASSGRAV_WAVELET_TOL_MAX=0.001;
    double MASSGRAV_WAVELET_TOL_FUNCTION_R0=10.0;
    double MASSGRAV_WAVELET_TOL_FUNCTION_R1=50.0;

    double MASSGRAV_CFL_FACTOR=0.1;

    double MASSGRAV_LOAD_IMB_TOL=0.1;
    unsigned int MASSGRAV_SPLIT_FIX=2;
    unsigned int MASSGRAV_ASYNC_COMM_K=4;
    double MASSGRAV_RK_TIME_BEGIN=0;
    double MASSGRAV_RK_TIME_END=10;

    double MASSGRAV_RK45_DESIRED_TOL=1e-6;

    unsigned int MASSGRAV_RK_TYPE;


    unsigned int MASSGRAV_CHECKPT_FREQ=10;
    unsigned int MASSGRAV_RESTORE_SOLVER=0;
    unsigned int MASSGRAV_ENABLE_BLOCK_ADAPTIVITY=0;

    double MASSGRAV_ETA_R0=1.31;
    double MASSGRAV_ETA_POWER[]={2.0,2.0};

    std::string MASSGRAV_VTU_FILE_PREFIX="massgrav_gr";
    std::string MASSGRAV_CHKPT_FILE_PREFIX="massgrav_cp";
    std::string MASSGRAV_PROFILE_FILE_PREFIX="massgrav_prof";


    BH BH1;
    BH BH2;
    unsigned int MASSGRAV_DIM=3;
    unsigned int MASSGRAV_MAXDEPTH=8;


    unsigned int MASSGRAV_ID_TYPE=0;

    double MASSGRAV_GRID_MIN_X=-50.0;
    double MASSGRAV_GRID_MAX_X=50.0;
    double MASSGRAV_GRID_MIN_Y=-50.0;
    double MASSGRAV_GRID_MAX_Y=50.0;
    double MASSGRAV_GRID_MIN_Z=-50.0;
    double MASSGRAV_GRID_MAX_Z=50.0;

    double MASSGRAV_BLK_MIN_X=-6.0;
    double MASSGRAV_BLK_MIN_Y=-6.0;
    double MASSGRAV_BLK_MIN_Z=-6.0;

    double MASSGRAV_BLK_MAX_X=6.0;
    double MASSGRAV_BLK_MAX_Y=6.0;
    double MASSGRAV_BLK_MAX_Z=6.0;

    double MASSGRAV_COMPD_MIN[3]={MASSGRAV_GRID_MIN_X,MASSGRAV_GRID_MIN_Y,MASSGRAV_GRID_MIN_Z};
    double MASSGRAV_COMPD_MAX[3]={MASSGRAV_GRID_MAX_X,MASSGRAV_GRID_MAX_Y,MASSGRAV_GRID_MAX_Z};

    double MASSGRAV_OCTREE_MIN[3]={0.0,0.0,0.0};
    double MASSGRAV_OCTREE_MAX[3]={(double)(1u<<MASSGRAV_MAXDEPTH),(double)(1u<<MASSGRAV_MAXDEPTH),(double)(1u<<MASSGRAV_MAXDEPTH)};

    //@note assumes the computational domain is a cube as well.
    double MASSGRAV_RK45_TIME_STEP_SIZE=MASSGRAV_CFL_FACTOR*(MASSGRAV_COMPD_MAX[0]-MASSGRAV_COMPD_MIN[0])*(1.0/(double)(1u<<MASSGRAV_MAXDEPTH));

    unsigned int MASSGRAV_LAMBDA[4]={1, 1, 1, 1};
    double MASSGRAV_LAMBDA_F[2]={1.0, 0.0};
    double MASSGRAV_TRK0=0.0;
    double ETA_CONST=2.0;
    double ETA_R0=50.0;
    double ETA_DAMPING=1.0;
    double ETA_DAMPING_EXP=1.0;
    double CHI_FLOOR=0.1;
    double KO_DISS_SIGMA=0.01;

    unsigned int DISSIPATION_TYPE=0;

    unsigned int MASSGRAV_DENDRO_GRAIN_SZ=1000;

    double MASSGRAV_DENDRO_AMR_FAC=0.1;

    unsigned int MASSGRAV_NUM_REFINE_VARS=1;
    unsigned int MASSGRAV_REFINE_VARIABLE_INDICES[MASSGRAV_NUM_VARS]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};

    unsigned int MASSGRAV_NUM_EVOL_VARS_VTU_OUTPUT=1;
    unsigned int MASSGRAV_NUM_CONST_VARS_VTU_OUTPUT=1;
    unsigned int MASSGRAV_VTU_OUTPUT_EVOL_INDICES[MASSGRAV_NUM_VARS]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
    unsigned int MASSGRAV_VTU_OUTPUT_CONST_INDICES[MASSGRAV_CONSTRAINT_NUM_VARS]={0,1,2,3,4,5};

    unsigned int MASSGRAV_XI[3]={0 , 0 , 0};

    unsigned int MASSGRAV_DISSIPATION_NC=0;

    unsigned int MASSGRAV_DISSIPATION_S=10;

    bool MASSGRAV_USE_FD_GRID_TRANSFER=false;

    double MASSGRAV_EH_REFINE_VAL = 0.3 ;
    double MASSGRAV_EH_COARSEN_VAL = 0.4;

    // by default use WAMR refinement. 
    RefinementMode MASSGRAV_REFINEMENT_MODE = RefinementMode::WAMR;

    bool MASSGRAV_VTU_Z_SLICE_ONLY = true;


}

namespace TPID {
  double target_M_plus=1.0;
  double target_M_minus=1.0;
  double par_m_plus=1.0;
  double par_m_minus=1.0;
  double par_b=4.0;
  double par_P_plus[3]={0.0,0.0,0.0};
  double par_P_minus[3]={0.0,0.0,0.0};
  double par_S_plus[3]={0.0,0.0,0.0};
  double par_S_minus[3]={0.0,0.0,0.0};
  double center_offset[3]={0.0,0.0,0.00014142135623730951};
  double initial_lapse_psi_exponent=-2;
  int npoints_A=30;
  int npoints_B=30;
  int npoints_phi=16;
  int give_bare_mass=0;
  int initial_lapse=2;
  int solve_momentum_constraint=1;
  int grid_setup_method=1;
  int verbose = 1;
  double adm_tol=1.0e-10;
  double Newton_tol=1.0e-10;
}


namespace BHLOC
{
    unsigned int EXTRACTION_VAR_ID=massgrav::VAR::U_ALPHA;
    double EXTRACTION_TOL=0.3;
}


namespace GW
{
    
    unsigned int MASSGRAV_GW_NUM_RADAII;
    
    unsigned int MASSGRAV_GW_NUM_LMODES;
    
    double MASSGRAV_GW_RADAII[MASSGRAV_GW_MAX_RADAII];
    
    unsigned int MASSGRAV_GW_SPIN=2;
    
    unsigned int MASSGRAV_GW_L_MODES[MASSGRAV_GW_MAX_LMODES];
}
