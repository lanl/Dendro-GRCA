//
// Created by milinda on 7/25/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief This file contains all the parameters related to QUADGRAV simulation.
*/
//

#ifndef SFCSORTBENCH_PARAMETERS_H
#define SFCSORTBENCH_PARAMETERS_H

#include "bh.h"
#include "grDef.h"
#include <string.h>
#include <iostream>

//Macro
#define QUADGRAV_EVOL

namespace quadgrav
{

    /**@brief element order*/
    extern unsigned int QUADGRAV_ELE_ORDER;

    /**@brief number of variables*/
    #ifdef QUADGRAV_EVOL

    static const unsigned int QUADGRAV_NUM_VARS=38;
    #else
    //Turn off QG terms
    static const unsigned int QUADGRAV_NUM_VARS=24;
    #endif
    /**@brief number of constraints variables*/
    static const unsigned int QUADGRAV_CONSTRAINT_NUM_VARS=6;

    /***@brief number of RK45 stages*/
    static const unsigned int QUADGRAV_RK45_STAGES=6;

    /***@brief number of RK4 stages*/
    static const unsigned int QUADGRAV_RK4_STAGES=4;

    /**@brief number of rk4 stages*/
    static const unsigned int QUADGRAV_RK3_STAGES=3;

    /**@brief: parameter used for adaptive time step update. */
    static const double QUADGRAV_SAFETY_FAC=0.8;

    /**@brief number of internal variables*/
    static const unsigned int QUADGRAV_NUM_VARS_INTENL=(QUADGRAV_RK45_STAGES+1)*QUADGRAV_NUM_VARS;

    /**@brief CFL stability number number (specifies how dt=QUADGRAV_CFL_FACTOR*dx)*/
    extern double QUADGRAV_CFL_FACTOR;

    /**@brief min bh domain add these to the parameter file.*/
    extern double QUADGRAV_COMPD_MIN[3];
    /**@brief min bh domain @todo add these to the parameter file. */
    extern double QUADGRAV_COMPD_MAX[3];

    /**@brief min coords of the OCTREE */
    extern double QUADGRAV_OCTREE_MIN[3];
    /**@brief max coords of the OCTREE */
    extern double QUADGRAV_OCTREE_MAX[3];

    /**@brief solution output frequency*/
    extern unsigned int QUADGRAV_IO_OUTPUT_FREQ;

    /**@brief Gravitational wave data extraction frequency*/
    extern unsigned int QUADGRAV_GW_EXTRACT_FREQ;

    /**@brief timestep norms out put freq.*/
    extern unsigned int QUADGRAV_TIME_STEP_OUTPUT_FREQ;

    /**@brief remesh test frequency*/
    extern unsigned int QUADGRAV_REMESH_TEST_FREQ;

    /**@brief checkpoint store frequency*/
    extern unsigned int QUADGRAV_CHECKPT_FREQ;

    /**@brief restore the solver from check point if set to 1. */
    extern unsigned int QUADGRAV_RESTORE_SOLVER;

    /**@brief use the block adaptivity and disable the AMR*/
    extern unsigned int QUADGRAV_ENABLE_BLOCK_ADAPTIVITY;

    /**@brief file prefix for VTU*/
    extern std::string QUADGRAV_VTU_FILE_PREFIX;

    /**@brief file prefix for write check point*/
    extern std::string QUADGRAV_CHKPT_FILE_PREFIX;

    /**@brief file prefix to write profile info.*/
    extern std::string QUADGRAV_PROFILE_FILE_PREFIX;

    /**@brief number of refine variables*/
    extern unsigned int QUADGRAV_NUM_REFINE_VARS;

    /**@brief indices of refine var ids*/
    extern unsigned int QUADGRAV_REFINE_VARIABLE_INDICES[QUADGRAV_NUM_VARS];

    /**@brief number of evolution variables written to vtu files*/
    extern unsigned int QUADGRAV_NUM_EVOL_VARS_VTU_OUTPUT;

    /**@brief number of constrint variables written to vtu files*/
    extern unsigned int QUADGRAV_NUM_CONST_VARS_VTU_OUTPUT;

    /**@brief evolution variable IDs written to vtu files*/
    extern unsigned int QUADGRAV_VTU_OUTPUT_EVOL_INDICES[QUADGRAV_NUM_VARS];

    /**@brief constraint variable IDs written to vtu files*/
    extern unsigned int QUADGRAV_VTU_OUTPUT_CONST_INDICES[QUADGRAV_CONSTRAINT_NUM_VARS];

    /**@brief solution output gap (instead of freq. we can use to output the solution if currentTime > lastIOOutputTime + QUADGRAV_IO_OUTPUT_GAP)*/
    extern  double QUADGRAV_IO_OUTPUT_GAP;

    /**@brief prefered grain sz to use when selecting active npes*/
    extern unsigned int QUADGRAV_DENDRO_GRAIN_SZ;

    /**@brief AMR coarsening factor (we coarsen if tol<QUADGRAV_DENDRO_AMR_FAC*QUADGRAV_WAVELET_TOL)*/
    extern double QUADGRAV_DENDRO_AMR_FAC;

    /**@brief wavelet tolerance value. */
    extern  double QUADGRAV_WAVELET_TOL;
    /**@brief load-imbalance tolerance value. */
    extern  double QUADGRAV_LOAD_IMB_TOL;
    /**@brief: Splitter fix value*/
    extern unsigned int QUADGRAV_SPLIT_FIX;

    /**@brief: async. communication at a time. (upper bound shoud be QUADGRAV_NUM_VARS) */
    extern unsigned int QUADGRAV_ASYNC_COMM_K;


    /**@brief simulation begin time. */
    extern double QUADGRAV_RK_TIME_BEGIN;
    /**@brief simulation end time*/
    extern double QUADGRAV_RK_TIME_END;
    /**@brief rk time step size. */
    extern double QUADGRAV_RK45_TIME_STEP_SIZE;

    /**@brief rk method type*/
    extern unsigned int QUADGRAV_RK_TYPE;

    /** desired tolerance value for the rk45 method (adaptive time stepping. )*/
    extern double QUADGRAV_RK45_DESIRED_TOL;

    /**@brief Black hole 1 */
    extern BH BH1;
    /**@brief Black hole 2 */
    extern BH BH2;

    /**@brief BBH initial data type */
    extern unsigned int QUADGRAV_ID_TYPE;

    /**@brief physical coordinates for grid, x_min */
    extern double QUADGRAV_GRID_MIN_X;

    /**@brief physical coordinates for grid, x_max */
    extern double QUADGRAV_GRID_MAX_X;

    /**@brief physical coordinates for grid, y_min */
    extern double QUADGRAV_GRID_MIN_Y;

    /**@brief physical coordinates for grid, y_max */
    extern double QUADGRAV_GRID_MAX_Y;

    /**@brief physical coordinates for grid, z_min */
    extern double QUADGRAV_GRID_MIN_Z;

    /**@brief physical coordinates for grid, z_max */
    extern double QUADGRAV_GRID_MAX_Z;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double QUADGRAV_BLK_MIN_X;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double QUADGRAV_BLK_MIN_Y;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double QUADGRAV_BLK_MIN_Z;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double QUADGRAV_BLK_MAX_X;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double QUADGRAV_BLK_MAX_Y;

    /**@brief physical coordinates for the blk adaptive x_min*/
    extern double QUADGRAV_BLK_MAX_Z;

    /**@brief: dimension of the grid*/
    extern unsigned int QUADGRAV_DIM;

    /**@brief: max refinement level*/
    extern unsigned int QUADGRAV_MAXDEPTH;

    /**@brief: lambda values for evolution */
    extern unsigned int QUADGRAV_LAMBDA[4];

    /**@brief: lambda values for evolution */
    extern double QUADGRAV_LAMBDA_F[2];

    /**@brief : parameters for eta_damping function */
    extern double QUADGRAV_ETA_R0;
    extern double QUADGRAV_ETA_POWER[2];

    /**@brief: lambda values for evolution */
    extern double QUADGRAV_TRK0;

    /**@brief: base value for eta in evolution */
    extern double ETA_CONST;

    /**@brief: eta_R0, radius where eta is damped for evolution */
    extern double ETA_R0;

    /**@brief: eta damping for evolution */
    extern double ETA_DAMPING;

    /**@brief: eta damping exponent for evolution */
    extern double ETA_DAMPING_EXP;

    /**@brief: chi floor value */
    extern double CHI_FLOOR;

    /**@brief: dissipation type */
    extern unsigned int DISSIPATION_TYPE;

    /**@brief: Kreiss-Oliger dissipation */
    extern double KO_DISS_SIGMA;

    /**@brief: QUADGRAV_USE_WAVELET_TOL_FUNCTION */
    extern unsigned int QUADGRAV_USE_WAVELET_TOL_FUNCTION;

    /**@brief: QUADGRAV_WAVELET_TOL_FUNCTION_R0 */
    extern double QUADGRAV_WAVELET_TOL_FUNCTION_R0;

    /**@brief: QUADGRAV_WAVELET_TOL_FUNCTION_R0 */
    extern double QUADGRAV_WAVELET_TOL_FUNCTION_R1;

    /**@brief: QUADGRAV_WAVELET_TOL_MAX */
    extern double QUADGRAV_WAVELET_TOL_MAX;

    /**@brief: eta function parameters*/
    extern double QUADGRAV_ETA_R0;

    /**@brief: eta function parameters (powers)*/
    extern double QUADGRAV_ETA_POWER[2];

    /**@brief: xi parameters to change between different gauge conditions*/
    extern unsigned int QUADGRAV_XI[3];

    /**@brief: @david can you please add some comeents for these parameters. */
    extern unsigned int QUADGRAV_DISSIPATION_NC;

    /**@brief: @david can you please add some comeents for these parameters. */
    extern unsigned int QUADGRAV_DISSIPATION_S;

    /**@brief: if true it will use finite differnce like grid transfer*/
    extern bool QUADGRAV_USE_FD_GRID_TRANSFER;

    /**@brief: tolerance for refinement based on EH */
    extern double QUADGRAV_EH_REFINE_VAL;
    
    /**@brief: tolerance for coarsen based on EH */
    extern double QUADGRAV_EH_COARSEN_VAL;

    /**@brief: refinement mode for the application*/
    extern RefinementMode QUADGRAV_REFINEMENT_MODE;
    
    /**@brief: if true output only the z slice*/
    extern bool QUADGRAV_VTU_Z_SLICE_ONLY;



}

namespace TPID {
  static const double TP_epsilon = 1.0e-6;
  static const int swap_xz = 0;
  static const int use_sources = 0;
  static const int rescale_sources = 0;
  static const int use_external_initial_guess = 0;
  static const int do_residuum_debug_output = 1;
  static const int do_initial_debug_output = 1;
  static const int multiply_old_lapse = 0;
  static const double TP_Tiny = 1.0e-15;
  static const double TP_Extend_Radius = 0.0;
  static const int Newton_maxit = 5;

  extern double target_M_plus;
  extern double target_M_minus;
  extern double par_m_plus;
  extern double par_m_minus;
  extern double par_b;
  extern double par_P_plus[3];
  extern double par_P_minus[3];
  extern double par_S_plus[3];
  extern double par_S_minus[3];
  extern double center_offset[3];
  extern double initial_lapse_psi_exponent;
  extern int npoints_A;
  extern int npoints_B;
  extern int npoints_phi;
  extern int give_bare_mass;
  extern int initial_lapse;
  extern int solve_momentum_constraint;
  extern int grid_setup_method;
  extern int verbose;
  extern double adm_tol;
  extern double Newton_tol;
}


/**@brief parameters related to BH location extraction*/
namespace BHLOC
{
    // Note: Current implementation of the BH location extraction is simple clustering property based on the initial location of the BH
    // Later we can add the code to compute the BH locations based on the event horizon extraction.
    
    
    /**@brief variable ID used for BH location extraction*/
    extern unsigned int EXTRACTION_VAR_ID;
    
    /**@brief tolerance for BH extraction*/
    extern double EXTRACTION_TOL;
    
    
}


/**@brief parameters related to GW extraction*/
namespace GW
{
    
    /**@brief max allowed radii values*/
    static const unsigned int QUADGRAV_GW_MAX_RADAII=20;
    
    /**@brief max allowed lmode values*/
    static const unsigned int QUADGRAV_GW_MAX_LMODES=8;

    static const unsigned int QUADGRAV_GW_LEBEDEV_PREC=25;
    
    /**@brief: number of different radius values for psi4 poly fit*/
    extern unsigned int QUADGRAV_GW_NUM_RADAII;
    
    /**@brief: number of L mode values*/
    extern unsigned int QUADGRAV_GW_NUM_LMODES;
    
    /**@brief: vallues of extraction radaii*/
    extern double QUADGRAV_GW_RADAII[QUADGRAV_GW_MAX_RADAII];
    
    /**@brief: value of the spin*/
    extern unsigned int QUADGRAV_GW_SPIN;
    
    /**@brief values for l modes in SWSH*/
    extern unsigned int QUADGRAV_GW_L_MODES[QUADGRAV_GW_MAX_LMODES];
    
    
}



#endif //SFCSORTBENCH_PARAMETERS_H
