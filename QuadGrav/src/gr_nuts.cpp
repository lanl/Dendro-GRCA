/**
 * @file gr_nuts.cpp
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief Solving QUADGRAV equations using non uniform time stepping. 
 * @version 0.1
 * @date 2020-03-04
 * 
 * @copyright Copyright (c) 2020, School of Computing, University of Utah. 
 * 
 */


#include "gr.h"
#include "grUtils.h"
#include "mpi.h"
#include "TreeNode.h"
#include "mesh.h"
#include <vector>
#include <iostream>
#include "rkQUADGRAV.h"
#include "octUtils.h"
#include "meshUtils.h"


int main (int argc, char** argv)
{
    // 0- NUTS 1-UTS
    unsigned int ts_mode=0;     
    
    if(argc<2)
    {
        std::cout<<"Usage: "<<argv[0]<<"paramFile TSMode(0){0-Spatially Adaptive Time Stepping(SATS, "<<GRN<<"default"<<NRM<<") , 1- Uniform Time Stepping.  }"<<std::endl;
        return 0;
    }
        
    if(argc>2)
        ts_mode = std::atoi(argv[2]);


    MPI_Init(&argc,&argv);
    MPI_Comm comm=MPI_COMM_WORLD;

    int rank,npes;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);


    // Print out CMAKE options
    if (!rank) {
        #ifdef QUADGRAV_COMPUTE_CONSTRAINTS
          std::cout<<GRN<<"  Compiled with QUADGRAV_COMPUTE_CONSTRAINTS"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without QUADGRAV_COMPUTE_CONSTRAINTS"<<NRM<<std::endl;
        #endif
        #ifdef QUADGRAV_ENABLE_VTU_CONSTRAINT_OUTPUT
          std::cout<<GRN<<"  Compiled with QUADGRAV_ENABLE_VTU_CONSTRAINT_OUTPUT"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without QUADGRAV_ENABLE_VTU_CONSTRAINT_OUTPUT"<<NRM<<std::endl;
        #endif
        #ifdef QUADGRAV_ENABLE_VTU_OUTPUT
          std::cout<<GRN<<"  Compiled with QUADGRAV_ENABLE_VTU_OUTPUT"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without QUADGRAV_ENABLE_VTU_OUTPUT"<<NRM<<std::endl;
        #endif
        #ifdef QUADGRAV_ETA_FUNCTION 
          std::cout<<GRN<<"  Compiled with  QUADGRAV_ETA_FUNCTION"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  QUADGRAV_ETA_FUNCTION"<<NRM<<std::endl;
        #endif
        #ifdef QUADGRAV_EXTRACT_BH_LOCATIONS 
          std::cout<<GRN<<"  Compiled with  QUADGRAV_EXTRACT_BH_LOCATIONS"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  QUADGRAV_EXTRACT_BH_LOCATIONS"<<NRM<<std::endl;
        #endif
        #ifdef QUADGRAV_EXTRACT_GRAVITATIONAL_WAVES 
          std::cout<<GRN<<"  Compiled with  QUADGRAV_EXTRACT_GRAVITATIONAL_WAVES"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  QUADGRAV_EXTRACT_GRAVITATIONAL_WAVES"<<NRM<<std::endl;
        #endif
        #ifdef QUADGRAV_EXTRACT_GRAVITATIONAL_WAVES 
          std::cout<<GRN<<"  Compiled with  QUADGRAV_EXTRACT_GRAVITATIONAL_WAVES"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  QUADGRAV_EXTRACT_GRAVITATIONAL_WAVES"<<NRM<<std::endl;
        #endif
        #ifdef QUADGRAV_GAUGE_ROCHESTER 
          std::cout<<GRN<<"  Compiled with  QUADGRAV_GAUGE_ROCHESTER"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  QUADGRAV_GAUGE_ROCHESTER"<<NRM<<std::endl;
        #endif
        #ifdef QUADGRAV_KERR_SCHILD_TEST 
          std::cout<<GRN<<"  Compiled with  QUADGRAV_KERR_SCHILD_TEST"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  QUADGRAV_KERR_SCHILD_TEST"<<NRM<<std::endl;
        #endif

        #ifdef QUADGRAV_REFINE_BASE_EH 
          std::cout<<GRN<<"  Compiled with  QUADGRAV_REFINE_BASE_EH"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  QUADGRAV_REFINE_BASE_EH"<<NRM<<std::endl;
        #endif

        #ifdef USE_FD_INTERP_FOR_UNZIP 
          std::cout<<GRN<<"  Compiled with  USE_FD_INTERP_FOR_UNZIP"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  USE_FD_INTERP_FOR_UNZIP"<<NRM<<std::endl;
        #endif

    }

    //1 . read the parameter file.
    if(!rank) std::cout<<" reading parameter file :"<<argv[1]<<std::endl;
    quadgrav::readParamFile(argv[1],comm);



    if(rank==1|| npes==1)
    {
        std::cout<<"parameters read: "<<std::endl;

        std::cout<<YLW<<"\tnpes :"<<npes<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_DIM :"<<quadgrav::QUADGRAV_DIM<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ELE_ORDER :"<<quadgrav::QUADGRAV_ELE_ORDER<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_CFL_FACTOR :"<<quadgrav::QUADGRAV_CFL_FACTOR<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_IO_OUTPUT_FREQ :"<<quadgrav::QUADGRAV_IO_OUTPUT_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_GW_EXTRACT_FREQ :"<<quadgrav::QUADGRAV_GW_EXTRACT_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_REMESH_TEST_FREQ :"<<quadgrav::QUADGRAV_REMESH_TEST_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_CHECKPT_FREQ :"<<quadgrav::QUADGRAV_CHECKPT_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_RESTORE_SOLVER :"<<quadgrav::QUADGRAV_RESTORE_SOLVER<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ENABLE_BLOCK_ADAPTIVITY :"<<quadgrav::QUADGRAV_ENABLE_BLOCK_ADAPTIVITY<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_VTU_FILE_PREFIX :"<<quadgrav::QUADGRAV_VTU_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_CHKPT_FILE_PREFIX :"<<quadgrav::QUADGRAV_CHKPT_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_PROFILE_FILE_PREFIX :"<<quadgrav::QUADGRAV_PROFILE_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_VTU_Z_SLICE_ONLY :"<<quadgrav::QUADGRAV_VTU_Z_SLICE_ONLY<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_IO_OUTPUT_GAP :"<<quadgrav::QUADGRAV_IO_OUTPUT_GAP<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_DENDRO_GRAIN_SZ :"<<quadgrav::QUADGRAV_DENDRO_GRAIN_SZ<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ASYNC_COMM_K :"<<quadgrav::QUADGRAV_ASYNC_COMM_K<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_DENDRO_AMR_FAC :"<<quadgrav::QUADGRAV_DENDRO_AMR_FAC<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_USE_WAVELET_TOL_FUNCTION :"<<quadgrav::QUADGRAV_USE_WAVELET_TOL_FUNCTION<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_USE_FD_GRID_TRANSFER :"<<quadgrav::QUADGRAV_USE_FD_GRID_TRANSFER<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_WAVELET_TOL :"<<quadgrav::QUADGRAV_WAVELET_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_WAVELET_TOL_MAX:"<<quadgrav::QUADGRAV_WAVELET_TOL_MAX<<NRM<<std::endl;
        std::cout<<YLW<<"\t:QUADGRAV_WAVELET_TOL_FUNCTION_R0: "<<quadgrav::QUADGRAV_WAVELET_TOL_FUNCTION_R0<<NRM<<std::endl;
        std::cout<<YLW<<"\t:QUADGRAV_WAVELET_TOL_FUNCTION_R1: "<<quadgrav::QUADGRAV_WAVELET_TOL_FUNCTION_R1<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_LOAD_IMB_TOL :"<<quadgrav::QUADGRAV_LOAD_IMB_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_RK_TIME_BEGIN :"<<quadgrav::QUADGRAV_RK_TIME_BEGIN<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_RK_TIME_END :"<<quadgrav::QUADGRAV_RK_TIME_END<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_RK_TYPE :"<<quadgrav::QUADGRAV_RK_TYPE<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_RK45_TIME_STEP_SIZE :"<<quadgrav::QUADGRAV_RK45_TIME_STEP_SIZE<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_RK45_DESIRED_TOL :"<<quadgrav::QUADGRAV_RK45_DESIRED_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_COMPD_MIN : ( :"<<quadgrav::QUADGRAV_COMPD_MIN[0]<<" ,"<<quadgrav::QUADGRAV_COMPD_MIN[1]<<","<<quadgrav::QUADGRAV_COMPD_MIN[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_COMPD_MAX : ( :"<<quadgrav::QUADGRAV_COMPD_MAX[0]<<" ,"<<quadgrav::QUADGRAV_COMPD_MAX[1]<<","<<quadgrav::QUADGRAV_COMPD_MAX[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_BLK_MIN : ( :"<<quadgrav::QUADGRAV_BLK_MIN_X<<" ,"<<quadgrav::QUADGRAV_BLK_MIN_Y<<","<<quadgrav::QUADGRAV_BLK_MIN_Z<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_BLK_MAX : ( :"<<quadgrav::QUADGRAV_BLK_MAX_X<<" ,"<<quadgrav::QUADGRAV_BLK_MAX_Y<<","<<quadgrav::QUADGRAV_BLK_MAX_Z<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_OCTREE_MIN : ( :"<<quadgrav::QUADGRAV_OCTREE_MIN[0]<<" ,"<<quadgrav::QUADGRAV_OCTREE_MIN[1]<<","<<quadgrav::QUADGRAV_OCTREE_MIN[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_OCTREE_MAX : ( :"<<quadgrav::QUADGRAV_OCTREE_MAX[0]<<" ,"<<quadgrav::QUADGRAV_OCTREE_MAX[1]<<","<<quadgrav::QUADGRAV_OCTREE_MAX[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tETA_CONST :"<<quadgrav::ETA_CONST<<NRM<<std::endl;
        std::cout<<YLW<<"\tETA_R0 :"<<quadgrav::ETA_R0<<NRM<<std::endl;
        std::cout<<YLW<<"\tETA_DAMPING :"<<quadgrav::ETA_DAMPING<<NRM<<std::endl;
        std::cout<<YLW<<"\tETA_DAMPING_EXP :"<<quadgrav::ETA_DAMPING_EXP<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ETA_R0 :"<<quadgrav::QUADGRAV_ETA_R0<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_ETA_POWER : ("<<quadgrav::QUADGRAV_ETA_POWER[0]<<" ,"<<quadgrav::QUADGRAV_ETA_POWER[1]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_LAMBDA : ("<<quadgrav::QUADGRAV_LAMBDA[0]<<" ,"<<quadgrav::QUADGRAV_LAMBDA[1]<<","<<quadgrav::QUADGRAV_LAMBDA[2]<<quadgrav::QUADGRAV_LAMBDA[3]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_LAMBDA_F : ("<<quadgrav::QUADGRAV_LAMBDA_F[0]<<" ,"<<quadgrav::QUADGRAV_LAMBDA_F[1]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_XI : ("<<quadgrav::QUADGRAV_XI[0]<<" ,"<<quadgrav::QUADGRAV_XI[1]<<" ,"<<quadgrav::QUADGRAV_XI[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tCHI_FLOOR :"<<quadgrav::CHI_FLOOR<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_TRK0 :"<<quadgrav::QUADGRAV_TRK0<<NRM<<std::endl;
        std::cout<<YLW<<"\tDISSIPATION_TYPE :"<<quadgrav::DISSIPATION_TYPE<<NRM<<std::endl;
        std::cout<<YLW<<"\tKO_DISS_SIGMA :"<<quadgrav::KO_DISS_SIGMA<<NRM<<std::endl;

        std::cout<<YLW<<"\tBH1 MASS :"<<quadgrav::BH1.getBHMass()<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH1 POSITION (x,y,z) : ("<<quadgrav::BH1.getBHCoordX()<<", "<<quadgrav::BH1.getBHCoordY()<<", "<<quadgrav::BH1.getBHCoordZ()<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH1 VELOCITY (x,y,z) : ("<<quadgrav::BH1.getVx()<<", "<<quadgrav::BH1.getVy()<<", "<<quadgrav::BH1.getVz()<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH1 SPIN (||,theta,phi): ( "<<quadgrav::BH1.getBHSpin()<<", "<<quadgrav::BH1.getBHSpinTheta()<<", "<<quadgrav::BH1.getBHSpinPhi()<<" )"<<NRM<<std::endl;

        std::cout<<YLW<<"\tBH2 MASS :"<<quadgrav::BH2.getBHMass()<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH2 POSITION (x,y,z) : ("<<quadgrav::BH2.getBHCoordX()<<", "<<quadgrav::BH2.getBHCoordY()<<", "<<quadgrav::BH2.getBHCoordZ()<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH2 VELOCITY (x,y,z) : ("<<quadgrav::BH2.getVx()<<", "<<quadgrav::BH2.getVy()<<", "<<quadgrav::BH2.getVz()<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH2 SPIN (||,theta,phi): ( "<<quadgrav::BH2.getBHSpin()<<", "<<quadgrav::BH2.getBHSpinTheta()<<", "<<quadgrav::BH2.getBHSpinPhi()<<" )"<<NRM<<std::endl;

        std::cout<<YLW<<"\tQUADGRAV_DIM :"<<quadgrav::QUADGRAV_DIM<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_MAXDEPTH :"<<quadgrav::QUADGRAV_MAXDEPTH<<NRM<<std::endl;

        std::cout<<YLW<<"\tQUADGRAV_NUM_REFINE_VARS :"<<quadgrav::QUADGRAV_NUM_REFINE_VARS<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_REFINE_VARIABLE_INDICES :[";
        for(unsigned int i=0;i<quadgrav::QUADGRAV_NUM_REFINE_VARS-1;i++)
            std::cout<<quadgrav::QUADGRAV_REFINE_VARIABLE_INDICES[i]<<", ";
        std::cout<<quadgrav::QUADGRAV_REFINE_VARIABLE_INDICES[quadgrav::QUADGRAV_NUM_REFINE_VARS-1]<<"]"<<NRM<<std::endl;

        std::cout<<YLW<<"\tQUADGRAV_REFINEMENT_MODE :"<<quadgrav::QUADGRAV_REFINEMENT_MODE<<NRM<<std::endl;

        #ifdef QUADGRAV_REFINE_BASE_EH
                std::cout<<YLW<<"\tQUADGRAV_EH_REFINE_VAL  : "<<quadgrav::QUADGRAV_EH_REFINE_VAL<<NRM<<std::endl;
                std::cout<<YLW<<"\tQUADGRAV_EH_COARSEN_VAL : "<<quadgrav::QUADGRAV_EH_COARSEN_VAL<<NRM<<std::endl;
        #endif 

        std::cout<<YLW<<"\tQUADGRAV_NUM_EVOL_VARS_VTU_OUTPUT :"<<quadgrav::QUADGRAV_NUM_EVOL_VARS_VTU_OUTPUT<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_VTU_OUTPUT_EVOL_INDICES :[";
        for(unsigned int i=0;i<quadgrav::QUADGRAV_NUM_EVOL_VARS_VTU_OUTPUT-1;i++)
            std::cout<<quadgrav::QUADGRAV_VTU_OUTPUT_EVOL_INDICES[i]<<", ";
        std::cout<<quadgrav::QUADGRAV_VTU_OUTPUT_EVOL_INDICES[quadgrav::QUADGRAV_NUM_EVOL_VARS_VTU_OUTPUT-1]<<"]"<<NRM<<std::endl;

        std::cout<<YLW<<"\tQUADGRAV_NUM_CONST_VARS_VTU_OUTPUT :"<<quadgrav::QUADGRAV_NUM_CONST_VARS_VTU_OUTPUT<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_VTU_OUTPUT_CONST_INDICES :[";
        for(unsigned int i=0;i<quadgrav::QUADGRAV_NUM_CONST_VARS_VTU_OUTPUT-1;i++)
            std::cout<<quadgrav::QUADGRAV_VTU_OUTPUT_CONST_INDICES[i]<<", ";
        std::cout<<quadgrav::QUADGRAV_VTU_OUTPUT_CONST_INDICES[quadgrav::QUADGRAV_NUM_CONST_VARS_VTU_OUTPUT-1]<<"]"<<NRM<<std::endl;


        std::cout<<YLW<<"\tTPID_TARGET_M_PLUS :"<<TPID::target_M_plus<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_TARGET_M_MINUS :"<<TPID::target_M_minus<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_PAR_B :"<<TPID::par_b<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_PAR_P_PLUS : ( "<<TPID::par_P_plus[0]<<", "<<TPID::par_P_plus[1]<<", "<<TPID::par_P_plus[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_PAR_P_MINUS : ( "<<TPID::par_P_minus[0]<<", "<<TPID::par_P_minus[1]<<", "<<TPID::par_P_minus[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_PAR_S_PLUS : ( "<<TPID::par_S_plus[0]<<", "<<TPID::par_S_plus[1]<<", "<<TPID::par_S_plus[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_PAR_S_MINUS : ( "<<TPID::par_S_minus[0]<<", "<<TPID::par_S_minus[1]<<", "<<TPID::par_S_minus[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_CENTER_OFFSET : ( "<<TPID::center_offset[0]<<", "<<TPID::center_offset[1]<<", "<<TPID::center_offset[2]<<" )"<<NRM<<std::endl;

        std::cout<<YLW<<"\tTPID_INITIAL_LAPSE_PSI_EXPONENT :"<<TPID::initial_lapse_psi_exponent<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_NPOINTS_A :"<<TPID::npoints_A<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_NPOINTS_B :"<<TPID::npoints_B<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_NPOINTS_PHI :"<<TPID::npoints_phi<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_GIVE_BARE_MASS :"<<TPID::give_bare_mass<<NRM<<std::endl;
        std::cout<<YLW<<"\tINITIAL_LAPSE :"<<TPID::initial_lapse<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_SOLVE_MOMENTUM_CONSTRAINT :"<<TPID::solve_momentum_constraint<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_GRID_SETUP_METHOD :"<<TPID::grid_setup_method<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_VERBOSE :"<<TPID::verbose<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_ADM_TOL :"<<TPID::adm_tol<<NRM<<std::endl;
        std::cout<<YLW<<"\tTPID_NEWTON_TOL :"<<TPID::Newton_tol<<NRM<<std::endl;
        
        
        std::cout<<YLW<<"\tEXTRACTION_VAR_ID :"<<BHLOC::EXTRACTION_VAR_ID<<NRM<<std::endl;
        std::cout<<YLW<<"\tEXTRACTION_TOL :"<<BHLOC::EXTRACTION_TOL<<NRM<<std::endl;


        std::cout<<YLW<<"\tQUADGRAV_GW_NUM_RADAII: "<<GW::QUADGRAV_GW_NUM_RADAII<<NRM<<std::endl;
        std::cout<<YLW<<"\tQUADGRAV_GW_NUM_LMODES: "<<GW::QUADGRAV_GW_NUM_LMODES<<NRM<<std::endl;

        std::cout<<YLW<<"\tQUADGRAV_GW_RADAII: {";
        for(unsigned int i=0;i<GW::QUADGRAV_GW_NUM_RADAII;i++)
            std::cout<<" ,"<<GW::QUADGRAV_GW_RADAII[i];
        std::cout<<"}"<<NRM<<std::endl;

        std::cout<<YLW<<"\tQUADGRAV_GW_L_MODES: {";
        for(unsigned int i=0;i<GW::QUADGRAV_GW_NUM_LMODES;i++)
            std::cout<<" ,"<<GW::QUADGRAV_GW_L_MODES[i];
        std::cout<<"}"<<NRM<<std::endl;

        

    }

    _InitializeHcurve(quadgrav::QUADGRAV_DIM);
    m_uiMaxDepth=quadgrav::QUADGRAV_MAXDEPTH;
    
    if(quadgrav::QUADGRAV_NUM_VARS%quadgrav::QUADGRAV_ASYNC_COMM_K!=0)
    {
        if(!rank) std::cout<<"[overlap communication error]: total QUADGRAV_NUM_VARS: "<<quadgrav::QUADGRAV_NUM_VARS<<" is not divisable by QUADGRAV_ASYNC_COMM_K: "<<quadgrav::QUADGRAV_ASYNC_COMM_K<<std::endl;
        MPI_Abort(comm,0);
    }

    if(quadgrav::QUADGRAV_GW_EXTRACT_FREQ> quadgrav::QUADGRAV_IO_OUTPUT_FREQ)
    {
      if(!rank) std::cout<<" QUADGRAV_GW_EXTRACT_FREQ  should be less QUADGRAV_IO_OUTPUT_FREQ "<<std::endl;
      MPI_Abort(comm,0);
    }


    //2. generate the initial grid.
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){quadgrav::punctureData(x,y,z,var);};
    std::function<double(double,double,double)> f_init_alpha=[](double x,double y,double z){ double var[24]; quadgrav::punctureData(x,y,z,var); return var[0];};
    //std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){quadgrav::KerrSchildData(x,y,z,var);};

    const unsigned int interpVars=quadgrav::QUADGRAV_NUM_VARS;
    unsigned int varIndex[interpVars];
    for(unsigned int i=0;i<quadgrav::QUADGRAV_NUM_VARS;i++)
        varIndex[i]=i;

    /*varIndex[0]=quadgrav::VAR::U_ALPHA;
    varIndex[1]=quadgrav::VAR::U_CHI;*/
    DendroIntL localSz,globalSz;
    double t_stat;
    double t_stat_g[3];

    quadgrav::timer::t_f2o.start();

    if(quadgrav::QUADGRAV_ENABLE_BLOCK_ADAPTIVITY)
    {
        if(!rank) std::cout<<YLW<<"Using block adaptive mesh. AMR disabled "<<NRM<<std::endl;
        const Point pt_min(quadgrav::QUADGRAV_BLK_MIN_X,quadgrav::QUADGRAV_BLK_MIN_Y,quadgrav::QUADGRAV_BLK_MIN_Z);
        const Point pt_max(quadgrav::QUADGRAV_BLK_MAX_X,quadgrav::QUADGRAV_BLK_MAX_Y,quadgrav::QUADGRAV_BLK_MAX_Z);

        quadgrav::blockAdaptiveOctree(tmpNodes,pt_min,pt_max,m_uiMaxDepth-2,m_uiMaxDepth,comm);
    }else
    {

        if(!rank) std::cout<<YLW<<"Using function2Octree. AMR enabled "<<NRM<<std::endl;
        function2Octree(f_init,quadgrav::QUADGRAV_NUM_VARS,varIndex,interpVars,tmpNodes,m_uiMaxDepth,quadgrav::QUADGRAV_WAVELET_TOL,quadgrav::QUADGRAV_ELE_ORDER,comm);
    }

    //std::vector<ot::TreeNode> f2Octants(tmpNodes);

    //ot::Mesh * mesh = ot::createMesh(tmpNodes.data(),tmpNodes.size(),quadgrav::QUADGRAV_ELE_ORDER,comm,1,ot::SM_TYPE::FDM,quadgrav::QUADGRAV_DENDRO_GRAIN_SZ,quadgrav::QUADGRAV_LOAD_IMB_TOL,quadgrav::QUADGRAV_SPLIT_FIX, quadgrav::getOctantWeight, (m_uiMaxDepth-MAXDEAPTH_LEVEL_DIFF-1));
    ot::Mesh * mesh = ot::createMesh(tmpNodes.data(),tmpNodes.size(),quadgrav::QUADGRAV_ELE_ORDER,comm,1,ot::SM_TYPE::FDM,quadgrav::QUADGRAV_DENDRO_GRAIN_SZ,quadgrav::QUADGRAV_LOAD_IMB_TOL,quadgrav::QUADGRAV_SPLIT_FIX, quadgrav::getOctantWeight, 0);
    DendroIntL lblocks = mesh->getLocalBlockList().size();
    DendroIntL gblocks =0; 
    par::Mpi_Reduce(&lblocks,&gblocks,1,MPI_SUM,0,comm);
    if(!rank)
      std::cout<<" number of blocks for coarset block level : "<<(m_uiMaxDepth-MAXDEAPTH_LEVEL_DIFF-1)<<" # blocks: "<<gblocks<<std::endl;
    unsigned int lmin, lmax;
    mesh->computeMinMaxLevel(lmin,lmax);    
    quadgrav::QUADGRAV_RK45_TIME_STEP_SIZE=quadgrav::QUADGRAV_CFL_FACTOR*((quadgrav::QUADGRAV_COMPD_MAX[0]-quadgrav::QUADGRAV_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) quadgrav::QUADGRAV_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))));
    par::Mpi_Bcast(&quadgrav::QUADGRAV_RK45_TIME_STEP_SIZE,1,0,comm);


    if(!rank)
      std::cout<<"ts_mode: "<<ts_mode<<std::endl;    

    if(ts_mode == 0)
    { 
        //NUTS
        assert(ot::test::isBlkFlagsValid(mesh));
        //std::cout<<"unzip : "<<mesh->getDegOfFreedomUnZip()<<" all : "<<mesh->getDegOfFreedomUnZip()*24<<" sz _per dof : "<<((mesh->getDegOfFreedomUnZip()*24)/24)<<std::endl;
        
        quadgrav::QUADGRAVCtx *  quadgravCtx = new quadgrav::QUADGRAVCtx(mesh); 
        ts::ExplicitNUTS<DendroScalar,quadgrav::QUADGRAVCtx>*  enuts = new ts::ExplicitNUTS<DendroScalar,quadgrav::QUADGRAVCtx>(quadgravCtx);

        //double * vec = mesh->createVector<double>(f_init_alpha);
        //bool state =ot::test::isSubScatterMapValid<double>(mesh,enuts->get_sub_scatter_maps(),vec);
        //std::cout<<" subSM valid : "<<state<<std::endl;
        //delete [] vec;
        
        std::vector<double> ld_stat_g;
        enuts->set_evolve_vars(quadgravCtx->get_evolution_vars());
        
        if((RKType)quadgrav::QUADGRAV_RK_TYPE == RKType::RK3)
            enuts->set_ets_coefficients(ts::ETSType::RK3);
        else if((RKType)quadgrav::QUADGRAV_RK_TYPE == RKType::RK4)
            enuts->set_ets_coefficients(ts::ETSType::RK4);
        else if((RKType)quadgrav::QUADGRAV_RK_TYPE == RKType::RK45)
            enuts->set_ets_coefficients(ts::ETSType::RK5);
        
        const unsigned int rank_global = enuts->get_global_rank();
        for(enuts->init(); enuts->curr_time() < quadgrav::QUADGRAV_RK_TIME_END ; enuts->evolve())
        {
            const DendroIntL step = enuts->curr_step();
            const DendroScalar time = enuts->curr_time();    

            const bool isActive = enuts->is_active();
            

            if(!rank_global)
                std::cout<<GRN<<"[Explicit NUTS]: Executing step :  "<<enuts->curr_step()<<std::setw(10)<<"\tcurrent time :"<<enuts->curr_time()<<std::setw(10)<<"\t dt(min):"<<enuts->get_dt_min()<<std::setw(10)<<"\t dt(max):"<<enuts->get_dt_max()<<std::setw(10)<<"\t"<<NRM<<std::endl;

            quadgravCtx->terminal_output();  

            bool isRemesh = false;    
            if( (step % quadgrav::QUADGRAV_REMESH_TEST_FREQ) == 0 )
                isRemesh = quadgravCtx->is_remesh();

            
            
            if(isRemesh)
            {
                if(!rank_global)
                    std::cout<<"[Explicit NUTS]: Remesh triggered "<<std::endl;;

                quadgravCtx->remesh(quadgrav::QUADGRAV_DENDRO_GRAIN_SZ, quadgrav::QUADGRAV_LOAD_IMB_TOL,quadgrav::QUADGRAV_SPLIT_FIX,true,false,false);
                quadgravCtx->terminal_output();

            }
            
        
            enuts->sync_with_mesh();

            if((step % quadgrav::QUADGRAV_IO_OUTPUT_FREQ) == 0 )
            quadgravCtx -> write_vtu();   

            if( (step % quadgrav::QUADGRAV_CHECKPT_FREQ) == 0 )
            quadgravCtx -> write_checkpt();
            
        }

        delete quadgravCtx->get_mesh();    
        delete quadgravCtx;

        delete enuts;

    }else if(ts_mode==1)
    { 
        //UTS
        quadgrav::QUADGRAVCtx *  quadgravCtx = new quadgrav::QUADGRAVCtx(mesh);
        ts::ETS<DendroScalar,quadgrav::QUADGRAVCtx>* ets = new ts::ETS<DendroScalar,quadgrav::QUADGRAVCtx>(quadgravCtx);
        ets->set_evolve_vars(quadgravCtx->get_evolution_vars());
        
        if((RKType)quadgrav::QUADGRAV_RK_TYPE == RKType::RK3)
            ets->set_ets_coefficients(ts::ETSType::RK3);
        else if((RKType)quadgrav::QUADGRAV_RK_TYPE == RKType::RK4)
            ets->set_ets_coefficients(ts::ETSType::RK4);
        else if((RKType)quadgrav::QUADGRAV_RK_TYPE == RKType::RK45)
            ets->set_ets_coefficients(ts::ETSType::RK5);


        for(ets->init(); ets->curr_time() < quadgrav::QUADGRAV_RK_TIME_END ; ets->evolve())
        {
            const DendroIntL   step = ets->curr_step();
            const DendroScalar time = ets->curr_time();    

            const bool isActive = ets->is_active();
            const unsigned int rank_global = ets->get_global_rank();

            if(!rank_global)
            std::cout<<"[ETS] : Executing step :  "<<ets->curr_step()<<"\tcurrent time :"<<ets->curr_time()<<"\t dt:"<<ets->ts_size()<<"\t"<<std::endl;

            quadgravCtx->terminal_output();  

            bool isRemesh = false;    
            if( (step % quadgrav::QUADGRAV_REMESH_TEST_FREQ) == 0 )
            isRemesh = quadgravCtx->is_remesh();

            if(isRemesh)
            {
                if(!rank_global)
                    std::cout<<"[ETS] : Remesh is triggered.  \n";

                quadgravCtx->remesh(quadgrav::QUADGRAV_DENDRO_GRAIN_SZ, quadgrav::QUADGRAV_LOAD_IMB_TOL,quadgrav::QUADGRAV_SPLIT_FIX,true,false,false);
                quadgravCtx->terminal_output();

            }
            
            ets->sync_with_mesh();

            if((step % quadgrav::QUADGRAV_IO_OUTPUT_FREQ) == 0 )
            quadgravCtx -> write_vtu();   

            if( (step % quadgrav::QUADGRAV_CHECKPT_FREQ) == 0 )
            quadgravCtx -> write_checkpt();
            
        }

        delete quadgravCtx->get_mesh();    
        delete quadgravCtx;
        delete ets;

    }else if(ts_mode ==2)
    {
        profiler_t t_rt;
        t_rt.clear();


        // perform a comparison test between ets and enuts. 
        quadgrav::QUADGRAVCtx *  quadgravCtx_enuts = new quadgrav::QUADGRAVCtx(mesh); 
        quadgrav::QUADGRAVCtx *  quadgravCtx_ets = new quadgrav::QUADGRAVCtx(mesh); 

        ts::ExplicitNUTS<DendroScalar,quadgrav::QUADGRAVCtx>*  enuts = new ts::ExplicitNUTS<DendroScalar,quadgrav::QUADGRAVCtx>(quadgravCtx_enuts);
        ts::ETS<DendroScalar,quadgrav::QUADGRAVCtx>*           ets   = new ts::ETS<DendroScalar,quadgrav::QUADGRAVCtx>(quadgravCtx_ets);


        ets   -> set_evolve_vars(quadgravCtx_ets->get_evolution_vars());
        enuts -> set_evolve_vars(quadgravCtx_enuts->get_evolution_vars());
        

        if((RKType)quadgrav::QUADGRAV_RK_TYPE == RKType::RK3)
        {
            ets   -> set_ets_coefficients(ts::ETSType::RK3);
            enuts -> set_ets_coefficients(ts::ETSType::RK3);

        }else if((RKType)quadgrav::QUADGRAV_RK_TYPE == RKType::RK4)
        {
            ets   -> set_ets_coefficients(ts::ETSType::RK4);
            enuts -> set_ets_coefficients(ts::ETSType::RK4);

        }else if((RKType)quadgrav::QUADGRAV_RK_TYPE == RKType::RK45)
        {
            ets   -> set_ets_coefficients(ts::ETSType::RK5);
            enuts -> set_ets_coefficients(ts::ETSType::RK5);
        }
        
        unsigned int num_steps = ( enuts->get_dt_max() / enuts->get_dt_min() );
        const unsigned int rank_global = ets->get_global_rank();

        if(!rank_global) 
          std::cout<<" num_steps: "<<num_steps<<std::endl;

        //enuts->dump_load_statistics(std::cout);

        
        t_rt.start();

        for(ets->init(); ets->curr_step() < num_steps ; ets->evolve())
        {
            if(!rank_global)
                std::cout<<"[ETS] : Executing step :  "<<ets->curr_step()<<"\tcurrent time :"<<ets->curr_time()<<"\t dt:"<<ets->ts_size()<<"\t"<<std::endl;

        }

        t_rt.stop();
        t_stat=t_rt.snap;
        quadgrav::timer::computeOverallStats(&t_stat, t_stat_g, comm);


        if(!rank_global)
                std::cout<<"[ETS] : Executing step :  "<<ets->curr_step()<<"\tcurrent time :"<<ets->curr_time()<<"\t dt:"<<ets->ts_size()<<"\t"<<std::endl;

        if(!rank_global)
          printf("[ETS] time (s): (min, mean, max): (%f, %f , %f)\n", t_stat_g[0], t_stat_g[1], t_stat_g[2]);

        
        DVec evar_ets = quadgravCtx_ets->get_evolution_vars();

        t_rt.snapreset();

        t_rt.start();
        for(enuts->init(); enuts->curr_step() < 1 ; enuts->evolve());
        t_rt.stop();

        t_stat=t_rt.snap;
        quadgrav::timer::computeOverallStats(&t_stat, t_stat_g, comm);

        if(!rank_global)
          std::cout<<GRN<<"[Explicit NUTS]: Executing step :  "<<enuts->curr_step()<<std::setw(10)<<"\tcurrent time :"<<enuts->curr_time()<<std::setw(10)<<"\t dt(min):"<<enuts->get_dt_min()<<std::setw(10)<<"\t dt(max):"<<enuts->get_dt_max()<<std::setw(10)<<"\t"<<NRM<<std::endl;

        if(!rank_global)
          printf("[ENUTS] time (s): (min, mean, max): (%f, %f , %f)\n", t_stat_g[0], t_stat_g[1], t_stat_g[2]);


        DVec evar_enuts = quadgravCtx_enuts->get_evolution_vars();

        DVec evar_diff;
        evar_diff.VecCopy(evar_ets,false);
               

        if(!rank)
          std::cout<<"\n\n\n"<<std::endl;

        for(unsigned int v=0; v < evar_diff.GetDof(); v++)
        {
          double min,max;
          evar_diff.VecMinMax(mesh, min, max, v);
          if(!rank)
            std::cout<<"[ETS] vec(min, max): ( "<<min<<" \t"<<", "<<max<<" ) "<<std::endl;
        }

        evar_diff.VecCopy(evar_enuts,true);

        if(!rank)
          std::cout<<"\n\n\n"<<std::endl;

        for(unsigned int v=0; v < evar_diff.GetDof(); v++)
        {
          double min,max;
          evar_diff.VecMinMax(mesh, min, max, v);
          if(!rank)
            std::cout<<"[ENUTS] vec(min, max): ( "<<min<<" \t"<<", "<<max<<" ) "<<std::endl;
        }


        evar_diff.VecFMA(mesh,evar_ets,1,-1,true);

        if(!rank)
          std::cout<<"\n\n\n"<<std::endl;


        for(unsigned int v=0; v < evar_diff.GetDof(); v++)
        {
          double min,max;
          evar_diff.VecMinMax(mesh, min, max, v);
          if(!rank)
            std::cout<<"[diff] vec(min, max): ( "<<min<<" \t"<<", "<<max<<" ) "<<std::endl;
        
        }


        evar_diff.VecDestroy();

        delete quadgravCtx_enuts->get_mesh();    
        delete quadgravCtx_enuts;
        delete quadgravCtx_ets;
        delete enuts;
        delete ets;

        


    }


    MPI_Finalize();
    return 0; 


}
