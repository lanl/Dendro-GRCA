//
// Created by milinda on 7/25/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Header file for the GR simulation.
*/
//

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


    if(argc<2)
        std::cout<<"Usage: "<<argv[0]<<" paramFile"<<std::endl;

    MPI_Init(&argc,&argv);
    MPI_Comm comm=MPI_COMM_WORLD;

    int rank,npes;
    MPI_Comm_rank(comm,&rank);
    MPI_Comm_size(comm,&npes);

    massgrav::timer::initFlops();

    massgrav::timer::total_runtime.start();

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

    ot::Mesh * mesh= ot::createMesh(tmpNodes.data(),tmpNodes.size(),massgrav::MASSGRAV_ELE_ORDER,comm,1,ot::SM_TYPE::FDM,massgrav::MASSGRAV_DENDRO_GRAIN_SZ,massgrav::MASSGRAV_LOAD_IMB_TOL,massgrav::MASSGRAV_SPLIT_FIX, massgrav::getOctantWeight);
    unsigned int lmin, lmax;
    mesh->computeMinMaxLevel(lmin,lmax);    
    massgrav::MASSGRAV_RK45_TIME_STEP_SIZE=massgrav::MASSGRAV_CFL_FACTOR*((massgrav::MASSGRAV_COMPD_MAX[0]-massgrav::MASSGRAV_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) massgrav::MASSGRAV_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth))));
    par::Mpi_Bcast(&massgrav::MASSGRAV_RK45_TIME_STEP_SIZE,1,0,comm);


    ode::solver::RK_MASSGRAV rk_massgrav(mesh,massgrav::MASSGRAV_RK_TIME_BEGIN,massgrav::MASSGRAV_RK_TIME_END,massgrav::MASSGRAV_RK45_TIME_STEP_SIZE,(RKType)massgrav::MASSGRAV_RK_TYPE);
    
    if(massgrav::MASSGRAV_RESTORE_SOLVER==1)
        rk_massgrav.restoreCheckPoint(massgrav::MASSGRAV_CHKPT_FILE_PREFIX.c_str(),comm);

    massgrav::timer::t_rkSolve.start();
    rk_massgrav.rkSolve();
    massgrav::timer::t_rkSolve.stop();

    massgrav::timer::total_runtime.stop();
    rk_massgrav.freeMesh();

    
    MPI_Finalize();
    return 0;

}
