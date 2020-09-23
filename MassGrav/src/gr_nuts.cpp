/**
 * @file gr_nuts.cpp
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief Solving MASSGRAV equations using non uniform time stepping. 
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
#include "rkMASSGRAV.h"
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
        #ifdef MASSGRAV_COMPUTE_CONSTRAINTS
          std::cout<<GRN<<"  Compiled with MASSGRAV_COMPUTE_CONSTRAINTS"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without MASSGRAV_COMPUTE_CONSTRAINTS"<<NRM<<std::endl;
        #endif
        #ifdef MASSGRAV_ENABLE_VTU_CONSTRAINT_OUTPUT
          std::cout<<GRN<<"  Compiled with MASSGRAV_ENABLE_VTU_CONSTRAINT_OUTPUT"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without MASSGRAV_ENABLE_VTU_CONSTRAINT_OUTPUT"<<NRM<<std::endl;
        #endif
        #ifdef MASSGRAV_ENABLE_VTU_OUTPUT
          std::cout<<GRN<<"  Compiled with MASSGRAV_ENABLE_VTU_OUTPUT"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without MASSGRAV_ENABLE_VTU_OUTPUT"<<NRM<<std::endl;
        #endif
        #ifdef MASSGRAV_ETA_FUNCTION 
          std::cout<<GRN<<"  Compiled with  MASSGRAV_ETA_FUNCTION"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  MASSGRAV_ETA_FUNCTION"<<NRM<<std::endl;
        #endif
        #ifdef MASSGRAV_EXTRACT_BH_LOCATIONS 
          std::cout<<GRN<<"  Compiled with  MASSGRAV_EXTRACT_BH_LOCATIONS"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  MASSGRAV_EXTRACT_BH_LOCATIONS"<<NRM<<std::endl;
        #endif
        #ifdef MASSGRAV_EXTRACT_GRAVITATIONAL_WAVES 
          std::cout<<GRN<<"  Compiled with  MASSGRAV_EXTRACT_GRAVITATIONAL_WAVES"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  MASSGRAV_EXTRACT_GRAVITATIONAL_WAVES"<<NRM<<std::endl;
        #endif
        #ifdef MASSGRAV_EXTRACT_GRAVITATIONAL_WAVES 
          std::cout<<GRN<<"  Compiled with  MASSGRAV_EXTRACT_GRAVITATIONAL_WAVES"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  MASSGRAV_EXTRACT_GRAVITATIONAL_WAVES"<<NRM<<std::endl;
        #endif
        #ifdef MASSGRAV_GAUGE_ROCHESTER 
          std::cout<<GRN<<"  Compiled with  MASSGRAV_GAUGE_ROCHESTER"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  MASSGRAV_GAUGE_ROCHESTER"<<NRM<<std::endl;
        #endif
        #ifdef MASSGRAV_KERR_SCHILD_TEST 
          std::cout<<GRN<<"  Compiled with  MASSGRAV_KERR_SCHILD_TEST"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  MASSGRAV_KERR_SCHILD_TEST"<<NRM<<std::endl;
        #endif

        #ifdef MASSGRAV_REFINE_BASE_EH 
          std::cout<<GRN<<"  Compiled with  MASSGRAV_REFINE_BASE_EH"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  MASSGRAV_REFINE_BASE_EH"<<NRM<<std::endl;
        #endif

        #ifdef USE_FD_INTERP_FOR_UNZIP 
          std::cout<<GRN<<"  Compiled with  USE_FD_INTERP_FOR_UNZIP"<<NRM<<std::endl;
        #else
          std::cout<<RED<<"  Compiled without  USE_FD_INTERP_FOR_UNZIP"<<NRM<<std::endl;
        #endif

    }

    //1 . read the parameter file.
    if(!rank) std::cout<<" reading parameter file :"<<argv[1]<<std::endl;
    massgrav::readParamFile(argv[1],comm);



    if(rank==1|| npes==1)
    {
        std::cout<<"parameters read: "<<std::endl;

        std::cout<<YLW<<"\tnpes :"<<npes<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_DIM :"<<massgrav::MASSGRAV_DIM<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_ELE_ORDER :"<<massgrav::MASSGRAV_ELE_ORDER<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_CFL_FACTOR :"<<massgrav::MASSGRAV_CFL_FACTOR<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_IO_OUTPUT_FREQ :"<<massgrav::MASSGRAV_IO_OUTPUT_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_GW_EXTRACT_FREQ :"<<massgrav::MASSGRAV_GW_EXTRACT_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_REMESH_TEST_FREQ :"<<massgrav::MASSGRAV_REMESH_TEST_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_CHECKPT_FREQ :"<<massgrav::MASSGRAV_CHECKPT_FREQ<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_RESTORE_SOLVER :"<<massgrav::MASSGRAV_RESTORE_SOLVER<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_ENABLE_BLOCK_ADAPTIVITY :"<<massgrav::MASSGRAV_ENABLE_BLOCK_ADAPTIVITY<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_VTU_FILE_PREFIX :"<<massgrav::MASSGRAV_VTU_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_CHKPT_FILE_PREFIX :"<<massgrav::MASSGRAV_CHKPT_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_PROFILE_FILE_PREFIX :"<<massgrav::MASSGRAV_PROFILE_FILE_PREFIX<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_VTU_Z_SLICE_ONLY :"<<massgrav::MASSGRAV_VTU_Z_SLICE_ONLY<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_IO_OUTPUT_GAP :"<<massgrav::MASSGRAV_IO_OUTPUT_GAP<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_DENDRO_GRAIN_SZ :"<<massgrav::MASSGRAV_DENDRO_GRAIN_SZ<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_ASYNC_COMM_K :"<<massgrav::MASSGRAV_ASYNC_COMM_K<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_DENDRO_AMR_FAC :"<<massgrav::MASSGRAV_DENDRO_AMR_FAC<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_USE_WAVELET_TOL_FUNCTION :"<<massgrav::MASSGRAV_USE_WAVELET_TOL_FUNCTION<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_USE_FD_GRID_TRANSFER :"<<massgrav::MASSGRAV_USE_FD_GRID_TRANSFER<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_WAVELET_TOL :"<<massgrav::MASSGRAV_WAVELET_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_WAVELET_TOL_MAX:"<<massgrav::MASSGRAV_WAVELET_TOL_MAX<<NRM<<std::endl;
        std::cout<<YLW<<"\t:MASSGRAV_WAVELET_TOL_FUNCTION_R0: "<<massgrav::MASSGRAV_WAVELET_TOL_FUNCTION_R0<<NRM<<std::endl;
        std::cout<<YLW<<"\t:MASSGRAV_WAVELET_TOL_FUNCTION_R1: "<<massgrav::MASSGRAV_WAVELET_TOL_FUNCTION_R1<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_LOAD_IMB_TOL :"<<massgrav::MASSGRAV_LOAD_IMB_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_RK_TIME_BEGIN :"<<massgrav::MASSGRAV_RK_TIME_BEGIN<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_RK_TIME_END :"<<massgrav::MASSGRAV_RK_TIME_END<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_RK_TYPE :"<<massgrav::MASSGRAV_RK_TYPE<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_RK45_TIME_STEP_SIZE :"<<massgrav::MASSGRAV_RK45_TIME_STEP_SIZE<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_RK45_DESIRED_TOL :"<<massgrav::MASSGRAV_RK45_DESIRED_TOL<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_COMPD_MIN : ( :"<<massgrav::MASSGRAV_COMPD_MIN[0]<<" ,"<<massgrav::MASSGRAV_COMPD_MIN[1]<<","<<massgrav::MASSGRAV_COMPD_MIN[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_COMPD_MAX : ( :"<<massgrav::MASSGRAV_COMPD_MAX[0]<<" ,"<<massgrav::MASSGRAV_COMPD_MAX[1]<<","<<massgrav::MASSGRAV_COMPD_MAX[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_BLK_MIN : ( :"<<massgrav::MASSGRAV_BLK_MIN_X<<" ,"<<massgrav::MASSGRAV_BLK_MIN_Y<<","<<massgrav::MASSGRAV_BLK_MIN_Z<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_BLK_MAX : ( :"<<massgrav::MASSGRAV_BLK_MAX_X<<" ,"<<massgrav::MASSGRAV_BLK_MAX_Y<<","<<massgrav::MASSGRAV_BLK_MAX_Z<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_OCTREE_MIN : ( :"<<massgrav::MASSGRAV_OCTREE_MIN[0]<<" ,"<<massgrav::MASSGRAV_OCTREE_MIN[1]<<","<<massgrav::MASSGRAV_OCTREE_MIN[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_OCTREE_MAX : ( :"<<massgrav::MASSGRAV_OCTREE_MAX[0]<<" ,"<<massgrav::MASSGRAV_OCTREE_MAX[1]<<","<<massgrav::MASSGRAV_OCTREE_MAX[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tETA_CONST :"<<massgrav::ETA_CONST<<NRM<<std::endl;
        std::cout<<YLW<<"\tETA_R0 :"<<massgrav::ETA_R0<<NRM<<std::endl;
        std::cout<<YLW<<"\tETA_DAMPING :"<<massgrav::ETA_DAMPING<<NRM<<std::endl;
        std::cout<<YLW<<"\tETA_DAMPING_EXP :"<<massgrav::ETA_DAMPING_EXP<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_ETA_R0 :"<<massgrav::MASSGRAV_ETA_R0<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_ETA_POWER : ("<<massgrav::MASSGRAV_ETA_POWER[0]<<" ,"<<massgrav::MASSGRAV_ETA_POWER[1]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_LAMBDA : ("<<massgrav::MASSGRAV_LAMBDA[0]<<" ,"<<massgrav::MASSGRAV_LAMBDA[1]<<","<<massgrav::MASSGRAV_LAMBDA[2]<<massgrav::MASSGRAV_LAMBDA[3]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_LAMBDA_F : ("<<massgrav::MASSGRAV_LAMBDA_F[0]<<" ,"<<massgrav::MASSGRAV_LAMBDA_F[1]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_XI : ("<<massgrav::MASSGRAV_XI[0]<<" ,"<<massgrav::MASSGRAV_XI[1]<<" ,"<<massgrav::MASSGRAV_XI[2]<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tCHI_FLOOR :"<<massgrav::CHI_FLOOR<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_TRK0 :"<<massgrav::MASSGRAV_TRK0<<NRM<<std::endl;
        std::cout<<YLW<<"\tDISSIPATION_TYPE :"<<massgrav::DISSIPATION_TYPE<<NRM<<std::endl;
        std::cout<<YLW<<"\tKO_DISS_SIGMA :"<<massgrav::KO_DISS_SIGMA<<NRM<<std::endl;

        std::cout<<YLW<<"\tBH1 MASS :"<<massgrav::BH1.getBHMass()<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH1 POSITION (x,y,z) : ("<<massgrav::BH1.getBHCoordX()<<", "<<massgrav::BH1.getBHCoordY()<<", "<<massgrav::BH1.getBHCoordZ()<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH1 VELOCITY (x,y,z) : ("<<massgrav::BH1.getVx()<<", "<<massgrav::BH1.getVy()<<", "<<massgrav::BH1.getVz()<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH1 SPIN (||,theta,phi): ( "<<massgrav::BH1.getBHSpin()<<", "<<massgrav::BH1.getBHSpinTheta()<<", "<<massgrav::BH1.getBHSpinPhi()<<" )"<<NRM<<std::endl;

        std::cout<<YLW<<"\tBH2 MASS :"<<massgrav::BH2.getBHMass()<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH2 POSITION (x,y,z) : ("<<massgrav::BH2.getBHCoordX()<<", "<<massgrav::BH2.getBHCoordY()<<", "<<massgrav::BH2.getBHCoordZ()<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH2 VELOCITY (x,y,z) : ("<<massgrav::BH2.getVx()<<", "<<massgrav::BH2.getVy()<<", "<<massgrav::BH2.getVz()<<" )"<<NRM<<std::endl;
        std::cout<<YLW<<"\tBH2 SPIN (||,theta,phi): ( "<<massgrav::BH2.getBHSpin()<<", "<<massgrav::BH2.getBHSpinTheta()<<", "<<massgrav::BH2.getBHSpinPhi()<<" )"<<NRM<<std::endl;

        std::cout<<YLW<<"\tMASSGRAV_DIM :"<<massgrav::MASSGRAV_DIM<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_MAXDEPTH :"<<massgrav::MASSGRAV_MAXDEPTH<<NRM<<std::endl;

        std::cout<<YLW<<"\tMASSGRAV_NUM_REFINE_VARS :"<<massgrav::MASSGRAV_NUM_REFINE_VARS<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_REFINE_VARIABLE_INDICES :[";
        for(unsigned int i=0;i<massgrav::MASSGRAV_NUM_REFINE_VARS-1;i++)
            std::cout<<massgrav::MASSGRAV_REFINE_VARIABLE_INDICES[i]<<", ";
        std::cout<<massgrav::MASSGRAV_REFINE_VARIABLE_INDICES[massgrav::MASSGRAV_NUM_REFINE_VARS-1]<<"]"<<NRM<<std::endl;

        std::cout<<YLW<<"\tMASSGRAV_REFINEMENT_MODE :"<<massgrav::MASSGRAV_REFINEMENT_MODE<<NRM<<std::endl;

        #ifdef MASSGRAV_REFINE_BASE_EH
                std::cout<<YLW<<"\tMASSGRAV_EH_REFINE_VAL  : "<<massgrav::MASSGRAV_EH_REFINE_VAL<<NRM<<std::endl;
                std::cout<<YLW<<"\tMASSGRAV_EH_COARSEN_VAL : "<<massgrav::MASSGRAV_EH_COARSEN_VAL<<NRM<<std::endl;
        #endif 

        std::cout<<YLW<<"\tMASSGRAV_NUM_EVOL_VARS_VTU_OUTPUT :"<<massgrav::MASSGRAV_NUM_EVOL_VARS_VTU_OUTPUT<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_VTU_OUTPUT_EVOL_INDICES :[";
        for(unsigned int i=0;i<massgrav::MASSGRAV_NUM_EVOL_VARS_VTU_OUTPUT-1;i++)
            std::cout<<massgrav::MASSGRAV_VTU_OUTPUT_EVOL_INDICES[i]<<", ";
        std::cout<<massgrav::MASSGRAV_VTU_OUTPUT_EVOL_INDICES[massgrav::MASSGRAV_NUM_EVOL_VARS_VTU_OUTPUT-1]<<"]"<<NRM<<std::endl;

        std::cout<<YLW<<"\tMASSGRAV_NUM_CONST_VARS_VTU_OUTPUT :"<<massgrav::MASSGRAV_NUM_CONST_VARS_VTU_OUTPUT<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_VTU_OUTPUT_CONST_INDICES :[";
        for(unsigned int i=0;i<massgrav::MASSGRAV_NUM_CONST_VARS_VTU_OUTPUT-1;i++)
            std::cout<<massgrav::MASSGRAV_VTU_OUTPUT_CONST_INDICES[i]<<", ";
        std::cout<<massgrav::MASSGRAV_VTU_OUTPUT_CONST_INDICES[massgrav::MASSGRAV_NUM_CONST_VARS_VTU_OUTPUT-1]<<"]"<<NRM<<std::endl;


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


        std::cout<<YLW<<"\tMASSGRAV_GW_NUM_RADAII: "<<GW::MASSGRAV_GW_NUM_RADAII<<NRM<<std::endl;
        std::cout<<YLW<<"\tMASSGRAV_GW_NUM_LMODES: "<<GW::MASSGRAV_GW_NUM_LMODES<<NRM<<std::endl;

        std::cout<<YLW<<"\tMASSGRAV_GW_RADAII: {";
        for(unsigned int i=0;i<GW::MASSGRAV_GW_NUM_RADAII;i++)
            std::cout<<" ,"<<GW::MASSGRAV_GW_RADAII[i];
        std::cout<<"}"<<NRM<<std::endl;

        std::cout<<YLW<<"\tMASSGRAV_GW_L_MODES: {";
        for(unsigned int i=0;i<GW::MASSGRAV_GW_NUM_LMODES;i++)
            std::cout<<" ,"<<GW::MASSGRAV_GW_L_MODES[i];
        std::cout<<"}"<<NRM<<std::endl;

        

    }

    _InitializeHcurve(massgrav::MASSGRAV_DIM);
    m_uiMaxDepth=massgrav::MASSGRAV_MAXDEPTH;
    
    if(massgrav::MASSGRAV_NUM_VARS%massgrav::MASSGRAV_ASYNC_COMM_K!=0)
    {
        if(!rank) std::cout<<"[overlap communication error]: total MASSGRAV_NUM_VARS: "<<massgrav::MASSGRAV_NUM_VARS<<" is not divisable by MASSGRAV_ASYNC_COMM_K: "<<massgrav::MASSGRAV_ASYNC_COMM_K<<std::endl;
        MPI_Abort(comm,0);
    }

    if(massgrav::MASSGRAV_GW_EXTRACT_FREQ> massgrav::MASSGRAV_IO_OUTPUT_FREQ)
    {
      if(!rank) std::cout<<" MASSGRAV_GW_EXTRACT_FREQ  should be less MASSGRAV_IO_OUTPUT_FREQ "<<std::endl;
      MPI_Abort(comm,0);
    }


    //2. generate the initial grid.
    std::vector<ot::TreeNode> tmpNodes;
    std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){massgrav::punctureData(x,y,z,var);};
    std::function<double(double,double,double)> f_init_alpha=[](double x,double y,double z){ double var[24]; massgrav::punctureData(x,y,z,var); return var[0];};
    //std::function<void(double,double,double,double*)> f_init=[](double x,double y,double z,double*var){massgrav::KerrSchildData(x,y,z,var);};

    const unsigned int interpVars=massgrav::MASSGRAV_NUM_VARS;
    unsigned int varIndex[interpVars];
    for(unsigned int i=0;i<massgrav::MASSGRAV_NUM_VARS;i++)
        varIndex[i]=i;

    /*varIndex[0]=massgrav::VAR::U_ALPHA;
    varIndex[1]=massgrav::VAR::U_CHI;*/
    DendroIntL localSz,globalSz;
    double t_stat;
    double t_stat_g[3];

    massgrav::timer::t_f2o.start();

    if(massgrav::MASSGRAV_ENABLE_BLOCK_ADAPTIVITY)
    {
        if(!rank) std::cout<<YLW<<"Using block adaptive mesh. AMR disabled "<<NRM<<std::endl;
        const Point pt_min(massgrav::MASSGRAV_BLK_MIN_X,massgrav::MASSGRAV_BLK_MIN_Y,massgrav::MASSGRAV_BLK_MIN_Z);
        const Point pt_max(massgrav::MASSGRAV_BLK_MAX_X,massgrav::MASSGRAV_BLK_MAX_Y,massgrav::MASSGRAV_BLK_MAX_Z);

        massgrav::blockAdaptiveOctree(tmpNodes,pt_min,pt_max,m_uiMaxDepth-2,m_uiMaxDepth,comm);
    }else
    {

        if(!rank) std::cout<<YLW<<"Using function2Octree. AMR enabled "<<NRM<<std::endl;
        function2Octree(f_init,massgrav::MASSGRAV_NUM_VARS,varIndex,interpVars,tmpNodes,m_uiMaxDepth,massgrav::MASSGRAV_WAVELET_TOL,massgrav::MASSGRAV_ELE_ORDER,comm);
    }

    //std::vector<ot::TreeNode> f2Octants(tmpNodes);

    //ot::Mesh * mesh = ot::createMesh(tmpNodes.data(),tmpNodes.size(),massgrav::MASSGRAV_ELE_ORDER,comm,1,ot::SM_TYPE::FDM,massgrav::MASSGRAV_DENDRO_GRAIN_SZ,massgrav::MASSGRAV_LOAD_IMB_TOL,massgrav::MASSGRAV_SPLIT_FIX, massgrav::getOctantWeight, (m_uiMaxDepth-MAXDEAPTH_LEVEL_DIFF-1));
    ot::Mesh * mesh = ot::createMesh(tmpNodes.data(),tmpNodes.size(),massgrav::MASSGRAV_ELE_ORDER,comm,1,ot::SM_TYPE::FDM,massgrav::MASSGRAV_DENDRO_GRAIN_SZ,massgrav::MASSGRAV_LOAD_IMB_TOL,massgrav::MASSGRAV_SPLIT_FIX, massgrav::getOctantWeight, 0);
    DendroIntL lblocks = mesh->getLocalBlockList().size();
    DendroIntL gblocks =0; 
    par::Mpi_Reduce(&lblocks,&gblocks,1,MPI_SUM,0,comm);
    if(!rank)
      std::cout<<" number of blocks for coarset block level : "<<(m_uiMaxDepth-MAXDEAPTH_LEVEL_DIFF-1)<<" # blocks: "<<gblocks<<std::endl;
    unsigned int lmin, lmax;
    mesh->computeMinMaxLevel(lmin,lmax);    
    massgrav::MASSGRAV_RK45_TIME_STEP_SIZE=massgrav::MASSGRAV_CFL_FACTOR*((massgrav::MASSGRAV_COMPD_MAX[0]-massgrav::MASSGRAV_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) massgrav::MASSGRAV_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))));
    par::Mpi_Bcast(&massgrav::MASSGRAV_RK45_TIME_STEP_SIZE,1,0,comm);


    if(!rank)
      std::cout<<"ts_mode: "<<ts_mode<<std::endl;    

    if(ts_mode == 0)
    { 
        //NUTS
        assert(ot::test::isBlkFlagsValid(mesh));
        //std::cout<<"unzip : "<<mesh->getDegOfFreedomUnZip()<<" all : "<<mesh->getDegOfFreedomUnZip()*24<<" sz _per dof : "<<((mesh->getDegOfFreedomUnZip()*24)/24)<<std::endl;
        
        massgrav::MASSGRAVCtx *  massgravCtx = new massgrav::MASSGRAVCtx(mesh); 
        ts::ExplicitNUTS<DendroScalar,massgrav::MASSGRAVCtx>*  enuts = new ts::ExplicitNUTS<DendroScalar,massgrav::MASSGRAVCtx>(massgravCtx);

        //double * vec = mesh->createVector<double>(f_init_alpha);
        //bool state =ot::test::isSubScatterMapValid<double>(mesh,enuts->get_sub_scatter_maps(),vec);
        //std::cout<<" subSM valid : "<<state<<std::endl;
        //delete [] vec;
        
        std::vector<double> ld_stat_g;
        enuts->set_evolve_vars(massgravCtx->get_evolution_vars());
        
        if((RKType)massgrav::MASSGRAV_RK_TYPE == RKType::RK3)
            enuts->set_ets_coefficients(ts::ETSType::RK3);
        else if((RKType)massgrav::MASSGRAV_RK_TYPE == RKType::RK4)
            enuts->set_ets_coefficients(ts::ETSType::RK4);
        else if((RKType)massgrav::MASSGRAV_RK_TYPE == RKType::RK45)
            enuts->set_ets_coefficients(ts::ETSType::RK5);
        
        const unsigned int rank_global = enuts->get_global_rank();
        for(enuts->init(); enuts->curr_time() < massgrav::MASSGRAV_RK_TIME_END ; enuts->evolve())
        {
            const DendroIntL step = enuts->curr_step();
            const DendroScalar time = enuts->curr_time();    

            const bool isActive = enuts->is_active();
            

            if(!rank_global)
                std::cout<<GRN<<"[Explicit NUTS]: Executing step :  "<<enuts->curr_step()<<std::setw(10)<<"\tcurrent time :"<<enuts->curr_time()<<std::setw(10)<<"\t dt(min):"<<enuts->get_dt_min()<<std::setw(10)<<"\t dt(max):"<<enuts->get_dt_max()<<std::setw(10)<<"\t"<<NRM<<std::endl;

            massgravCtx->terminal_output();  

            bool isRemesh = false;    
            if( (step % massgrav::MASSGRAV_REMESH_TEST_FREQ) == 0 )
                isRemesh = massgravCtx->is_remesh();

            
            
            if(isRemesh)
            {
                if(!rank_global)
                    std::cout<<"[Explicit NUTS]: Remesh triggered "<<std::endl;;

                massgravCtx->remesh(massgrav::MASSGRAV_DENDRO_GRAIN_SZ, massgrav::MASSGRAV_LOAD_IMB_TOL,massgrav::MASSGRAV_SPLIT_FIX,true,false,false);
                massgravCtx->terminal_output();

            }
            
        
            enuts->sync_with_mesh();

            if((step % massgrav::MASSGRAV_IO_OUTPUT_FREQ) == 0 )
            massgravCtx -> write_vtu();   

            if( (step % massgrav::MASSGRAV_CHECKPT_FREQ) == 0 )
            massgravCtx -> write_checkpt();
            
        }

        delete massgravCtx->get_mesh();    
        delete massgravCtx;

        delete enuts;

    }else if(ts_mode==1)
    { 
        //UTS
        massgrav::MASSGRAVCtx *  massgravCtx = new massgrav::MASSGRAVCtx(mesh);
        ts::ETS<DendroScalar,massgrav::MASSGRAVCtx>* ets = new ts::ETS<DendroScalar,massgrav::MASSGRAVCtx>(massgravCtx);
        ets->set_evolve_vars(massgravCtx->get_evolution_vars());
        
        if((RKType)massgrav::MASSGRAV_RK_TYPE == RKType::RK3)
            ets->set_ets_coefficients(ts::ETSType::RK3);
        else if((RKType)massgrav::MASSGRAV_RK_TYPE == RKType::RK4)
            ets->set_ets_coefficients(ts::ETSType::RK4);
        else if((RKType)massgrav::MASSGRAV_RK_TYPE == RKType::RK45)
            ets->set_ets_coefficients(ts::ETSType::RK5);


        for(ets->init(); ets->curr_time() < massgrav::MASSGRAV_RK_TIME_END ; ets->evolve())
        {
            const DendroIntL   step = ets->curr_step();
            const DendroScalar time = ets->curr_time();    

            const bool isActive = ets->is_active();
            const unsigned int rank_global = ets->get_global_rank();

            if(!rank_global)
            std::cout<<"[ETS] : Executing step :  "<<ets->curr_step()<<"\tcurrent time :"<<ets->curr_time()<<"\t dt:"<<ets->ts_size()<<"\t"<<std::endl;

            massgravCtx->terminal_output();  

            bool isRemesh = false;    
            if( (step % massgrav::MASSGRAV_REMESH_TEST_FREQ) == 0 )
            isRemesh = massgravCtx->is_remesh();

            if(isRemesh)
            {
                if(!rank_global)
                    std::cout<<"[ETS] : Remesh is triggered.  \n";

                massgravCtx->remesh(massgrav::MASSGRAV_DENDRO_GRAIN_SZ, massgrav::MASSGRAV_LOAD_IMB_TOL,massgrav::MASSGRAV_SPLIT_FIX,true,false,false);
                massgravCtx->terminal_output();

            }
            
            ets->sync_with_mesh();

            if((step % massgrav::MASSGRAV_IO_OUTPUT_FREQ) == 0 )
            massgravCtx -> write_vtu();   

            if( (step % massgrav::MASSGRAV_CHECKPT_FREQ) == 0 )
            massgravCtx -> write_checkpt();
            
        }

        delete massgravCtx->get_mesh();    
        delete massgravCtx;
        delete ets;

    }else if(ts_mode ==2)
    {
        profiler_t t_rt;
        t_rt.clear();


        // perform a comparison test between ets and enuts. 
        massgrav::MASSGRAVCtx *  massgravCtx_enuts = new massgrav::MASSGRAVCtx(mesh); 
        massgrav::MASSGRAVCtx *  massgravCtx_ets = new massgrav::MASSGRAVCtx(mesh); 

        ts::ExplicitNUTS<DendroScalar,massgrav::MASSGRAVCtx>*  enuts = new ts::ExplicitNUTS<DendroScalar,massgrav::MASSGRAVCtx>(massgravCtx_enuts);
        ts::ETS<DendroScalar,massgrav::MASSGRAVCtx>*           ets   = new ts::ETS<DendroScalar,massgrav::MASSGRAVCtx>(massgravCtx_ets);


        ets   -> set_evolve_vars(massgravCtx_ets->get_evolution_vars());
        enuts -> set_evolve_vars(massgravCtx_enuts->get_evolution_vars());
        

        if((RKType)massgrav::MASSGRAV_RK_TYPE == RKType::RK3)
        {
            ets   -> set_ets_coefficients(ts::ETSType::RK3);
            enuts -> set_ets_coefficients(ts::ETSType::RK3);

        }else if((RKType)massgrav::MASSGRAV_RK_TYPE == RKType::RK4)
        {
            ets   -> set_ets_coefficients(ts::ETSType::RK4);
            enuts -> set_ets_coefficients(ts::ETSType::RK4);

        }else if((RKType)massgrav::MASSGRAV_RK_TYPE == RKType::RK45)
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
        massgrav::timer::computeOverallStats(&t_stat, t_stat_g, comm);


        if(!rank_global)
                std::cout<<"[ETS] : Executing step :  "<<ets->curr_step()<<"\tcurrent time :"<<ets->curr_time()<<"\t dt:"<<ets->ts_size()<<"\t"<<std::endl;

        if(!rank_global)
          printf("[ETS] time (s): (min, mean, max): (%f, %f , %f)\n", t_stat_g[0], t_stat_g[1], t_stat_g[2]);

        
        DVec evar_ets = massgravCtx_ets->get_evolution_vars();

        t_rt.snapreset();

        t_rt.start();
        for(enuts->init(); enuts->curr_step() < 1 ; enuts->evolve());
        t_rt.stop();

        t_stat=t_rt.snap;
        massgrav::timer::computeOverallStats(&t_stat, t_stat_g, comm);

        if(!rank_global)
          std::cout<<GRN<<"[Explicit NUTS]: Executing step :  "<<enuts->curr_step()<<std::setw(10)<<"\tcurrent time :"<<enuts->curr_time()<<std::setw(10)<<"\t dt(min):"<<enuts->get_dt_min()<<std::setw(10)<<"\t dt(max):"<<enuts->get_dt_max()<<std::setw(10)<<"\t"<<NRM<<std::endl;

        if(!rank_global)
          printf("[ENUTS] time (s): (min, mean, max): (%f, %f , %f)\n", t_stat_g[0], t_stat_g[1], t_stat_g[2]);


        DVec evar_enuts = massgravCtx_enuts->get_evolution_vars();

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

        delete massgravCtx_enuts->get_mesh();    
        delete massgravCtx_enuts;
        delete massgravCtx_ets;
        delete enuts;
        delete ets;

        


    }


    MPI_Finalize();
    return 0; 


}
