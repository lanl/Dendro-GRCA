//
// Created by milinda on 7/26/17.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains utility functions for MASSGRAV simulation.
*/
//

#include "grUtils.h"

namespace massgrav
{

    void readParamFile(const char * fName,MPI_Comm comm)
    {


        json parFile;
        int rank,npes;
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        unsigned int vtu_len;
        unsigned int chp_len;
        unsigned int prf_len;

        if(!rank)
        {
            std::ifstream infile(fName);
            if(!infile) {std::cout<<fName<<" parameter file open failed "<<std::endl;}
            infile>>parFile;

            massgrav::MASSGRAV_IO_OUTPUT_FREQ=parFile["MASSGRAV_IO_OUTPUT_FREQ"];
            massgrav::MASSGRAV_REMESH_TEST_FREQ=parFile["MASSGRAV_REMESH_TEST_FREQ"];
            massgrav::MASSGRAV_CHECKPT_FREQ=parFile["MASSGRAV_CHECKPT_FREQ"];
            massgrav::MASSGRAV_IO_OUTPUT_GAP=parFile["MASSGRAV_IO_OUTPUT_GAP"];
            massgrav::MASSGRAV_VTU_FILE_PREFIX=parFile["MASSGRAV_VTU_FILE_PREFIX"].get<std::string>();
            massgrav::MASSGRAV_CHKPT_FILE_PREFIX=parFile["MASSGRAV_CHKPT_FILE_PREFIX"].get<std::string>();
            massgrav::MASSGRAV_PROFILE_FILE_PREFIX=parFile["MASSGRAV_PROFILE_FILE_PREFIX"].get<std::string>();
            massgrav::MASSGRAV_RESTORE_SOLVER=parFile["MASSGRAV_RESTORE_SOLVER"];
            massgrav::MASSGRAV_ENABLE_BLOCK_ADAPTIVITY=parFile["MASSGRAV_ENABLE_BLOCK_ADAPTIVITY"];
	        massgrav::MASSGRAV_ID_TYPE=parFile["MASSGRAV_ID_TYPE"];
            massgrav::MASSGRAV_BLK_MIN_X=parFile["MASSGRAV_BLK_MIN_X"];
            massgrav::MASSGRAV_BLK_MIN_Y=parFile["MASSGRAV_BLK_MIN_Y"];
            massgrav::MASSGRAV_BLK_MIN_Z=parFile["MASSGRAV_BLK_MIN_Z"];
            massgrav::MASSGRAV_BLK_MAX_X=parFile["MASSGRAV_BLK_MAX_X"];
            massgrav::MASSGRAV_BLK_MAX_Y=parFile["MASSGRAV_BLK_MAX_Y"];
            massgrav::MASSGRAV_BLK_MAX_Z=parFile["MASSGRAV_BLK_MAX_Z"];
            massgrav::MASSGRAV_DENDRO_GRAIN_SZ=parFile["MASSGRAV_DENDRO_GRAIN_SZ"];
            massgrav::MASSGRAV_ASYNC_COMM_K=parFile["MASSGRAV_ASYNC_COMM_K"];
            massgrav::MASSGRAV_DENDRO_AMR_FAC=parFile["MASSGRAV_DENDRO_AMR_FAC"];
            massgrav::MASSGRAV_LOAD_IMB_TOL=parFile["MASSGRAV_LOAD_IMB_TOL"];
            massgrav::MASSGRAV_RK_TIME_BEGIN=parFile["MASSGRAV_RK_TIME_BEGIN"];
            massgrav::MASSGRAV_RK_TIME_END=parFile["MASSGRAV_RK_TIME_END"];
            massgrav::MASSGRAV_RK_TYPE=parFile["MASSGRAV_RK_TYPE"];
            massgrav::MASSGRAV_RK45_TIME_STEP_SIZE=parFile["MASSGRAV_RK45_TIME_STEP_SIZE"];
            massgrav::MASSGRAV_RK45_DESIRED_TOL=parFile["MASSGRAV_RK45_DESIRED_TOL"];
            massgrav::MASSGRAV_DIM=parFile["MASSGRAV_DIM"];
            massgrav::MASSGRAV_MAXDEPTH=parFile["MASSGRAV_MAXDEPTH"];
            massgrav::BH1=BH((double)parFile["MASSGRAV_BH1"]["MASS"],(double)parFile["MASSGRAV_BH1"]["X"],(double)parFile["MASSGRAV_BH1"]["Y"],(double)parFile["MASSGRAV_BH1"]["Z"],(double)parFile["MASSGRAV_BH1"]["V_X"],(double)parFile["MASSGRAV_BH1"]["V_Y"],(double)parFile["MASSGRAV_BH1"]["V_Z"],(double)parFile["MASSGRAV_BH1"]["SPIN"],(double)parFile["MASSGRAV_BH1"]["SPIN_THETA"],(double)parFile["MASSGRAV_BH1"]["SPIN_PHI"]);
            massgrav::BH2=BH((double)parFile["MASSGRAV_BH2"]["MASS"],(double)parFile["MASSGRAV_BH2"]["X"],(double)parFile["MASSGRAV_BH2"]["Y"],(double)parFile["MASSGRAV_BH2"]["Z"],(double)parFile["MASSGRAV_BH2"]["V_X"],(double)parFile["MASSGRAV_BH2"]["V_Y"],(double)parFile["MASSGRAV_BH2"]["V_Z"],(double)parFile["MASSGRAV_BH2"]["SPIN"],(double)parFile["MASSGRAV_BH2"]["SPIN_THETA"],(double)parFile["MASSGRAV_BH2"]["SPIN_PHI"]);
            massgrav::MASSGRAV_GRID_MIN_X=parFile["MASSGRAV_GRID_MIN_X"];
            massgrav::MASSGRAV_GRID_MAX_X=parFile["MASSGRAV_GRID_MAX_X"];
            massgrav::MASSGRAV_GRID_MIN_Y=parFile["MASSGRAV_GRID_MIN_Y"];
            massgrav::MASSGRAV_GRID_MAX_Y=parFile["MASSGRAV_GRID_MAX_Y"];
            massgrav::MASSGRAV_GRID_MIN_Z=parFile["MASSGRAV_GRID_MIN_Z"];
            massgrav::MASSGRAV_GRID_MAX_Z=parFile["MASSGRAV_GRID_MAX_Z"];
            massgrav::ETA_CONST=parFile["ETA_CONST"];
            massgrav::ETA_R0=parFile["ETA_R0"];
            massgrav::ETA_DAMPING=parFile["ETA_DAMPING"];
            massgrav::ETA_DAMPING_EXP=parFile["ETA_DAMPING_EXP"];
            massgrav::MASSGRAV_LAMBDA[0]=(unsigned int) parFile["MASSGRAV_LAMBDA"]["MASSGRAV_LAMBDA_1"];
            massgrav::MASSGRAV_LAMBDA[1]=(unsigned int) parFile["MASSGRAV_LAMBDA"]["MASSGRAV_LAMBDA_2"];
            massgrav::MASSGRAV_LAMBDA[2]=(unsigned int) parFile["MASSGRAV_LAMBDA"]["MASSGRAV_LAMBDA_3"];
            massgrav::MASSGRAV_LAMBDA[3]=(unsigned int) parFile["MASSGRAV_LAMBDA"]["MASSGRAV_LAMBDA_4"];
            massgrav::MASSGRAV_LAMBDA_F[0]=parFile["MASSGRAV_LAMBDA_F"]["MASSGRAV_LAMBDA_F0"];
            massgrav::MASSGRAV_LAMBDA_F[1]=parFile["MASSGRAV_LAMBDA_F"]["MASSGRAV_LAMBDA_F1"];

            massgrav::MASSGRAV_XI[0] = (unsigned int ) parFile["MASSGRAV_XI"]["MASSGRAV_XI_0"];
            massgrav::MASSGRAV_XI[1] = (unsigned int ) parFile["MASSGRAV_XI"]["MASSGRAV_XI_1"];
            massgrav::MASSGRAV_XI[2] = (unsigned int ) parFile["MASSGRAV_XI"]["MASSGRAV_XI_2"];
            
            if(parFile.find("MASSGRAV_ELE_ORDER")!= parFile.end())
                massgrav::MASSGRAV_ELE_ORDER = parFile["MASSGRAV_ELE_ORDER"];
            
            massgrav::CHI_FLOOR=parFile["CHI_FLOOR"];
            massgrav::MASSGRAV_TRK0=parFile["MASSGRAV_TRK0"];
            if (parFile.find("DISSIPATION_TYPE") != parFile.end()) {
                massgrav::DISSIPATION_TYPE=parFile["DISSIPATION_TYPE"];
            }
            massgrav::KO_DISS_SIGMA=parFile["KO_DISS_SIGMA"];

  	        //Parameters for eta_damping function
	          massgrav::MASSGRAV_ETA_R0=parFile["MASSGRAV_ETA_R0"];
	          massgrav::MASSGRAV_ETA_POWER[0]=parFile["MASSGRAV_ETA_POWER"]["MASSGRAV_ETA_POWER_1"];
	          massgrav::MASSGRAV_ETA_POWER[1]=parFile["MASSGRAV_ETA_POWER"]["MASSGRAV_ETA_POWER_2"];

            massgrav::MASSGRAV_USE_WAVELET_TOL_FUNCTION=parFile["MASSGRAV_USE_WAVELET_TOL_FUNCTION"];
            massgrav::MASSGRAV_WAVELET_TOL=parFile["MASSGRAV_WAVELET_TOL"];
            massgrav::MASSGRAV_WAVELET_TOL_MAX=parFile["MASSGRAV_WAVELET_TOL_MAX"];
            massgrav::MASSGRAV_WAVELET_TOL_FUNCTION_R0=parFile["MASSGRAV_WAVELET_TOL_FUNCTION_R0"];
            massgrav::MASSGRAV_WAVELET_TOL_FUNCTION_R1=parFile["MASSGRAV_WAVELET_TOL_FUNCTION_R1"];

            massgrav::MASSGRAV_NUM_REFINE_VARS=parFile["MASSGRAV_NUM_REFINE_VARS"];
            for(unsigned int i=0;i<massgrav::MASSGRAV_NUM_REFINE_VARS;i++)
                massgrav::MASSGRAV_REFINE_VARIABLE_INDICES[i]=parFile["MASSGRAV_REFINE_VARIABLE_INDICES"][i];

            massgrav::MASSGRAV_NUM_EVOL_VARS_VTU_OUTPUT=parFile["MASSGRAV_NUM_EVOL_VARS_VTU_OUTPUT"];
            massgrav::MASSGRAV_NUM_CONST_VARS_VTU_OUTPUT=parFile["MASSGRAV_NUM_CONST_VARS_VTU_OUTPUT"];

            for(unsigned int i=0;i<massgrav::MASSGRAV_NUM_EVOL_VARS_VTU_OUTPUT;i++)
                massgrav::MASSGRAV_VTU_OUTPUT_EVOL_INDICES[i]=parFile["MASSGRAV_VTU_OUTPUT_EVOL_INDICES"][i];

            for(unsigned int i=0;i<massgrav::MASSGRAV_NUM_CONST_VARS_VTU_OUTPUT;i++)
                massgrav::MASSGRAV_VTU_OUTPUT_CONST_INDICES[i]=parFile["MASSGRAV_VTU_OUTPUT_CONST_INDICES"][i];

            if (parFile.find("MASSGRAV_CFL_FACTOR") != parFile.end()) {
                massgrav::MASSGRAV_CFL_FACTOR=parFile["MASSGRAV_CFL_FACTOR"];
            }

            if(parFile.find("MASSGRAV_VTU_Z_SLICE_ONLY") != parFile.end())
                massgrav::MASSGRAV_VTU_Z_SLICE_ONLY=parFile["MASSGRAV_VTU_Z_SLICE_ONLY"];

            if (parFile.find("MASSGRAV_GW_EXTRACT_FREQ") != parFile.end()) {
                massgrav::MASSGRAV_GW_EXTRACT_FREQ=parFile["MASSGRAV_GW_EXTRACT_FREQ"];
            }else
            {
                massgrav::MASSGRAV_GW_EXTRACT_FREQ=std::max(1u,massgrav::MASSGRAV_IO_OUTPUT_FREQ>>1u);
            }

           /* Parameters for TPID */
            TPID::target_M_plus=parFile["TPID_TARGET_M_PLUS"];
            TPID::target_M_minus=parFile["TPID_TARGET_M_MINUS"];
            TPID::par_m_plus=TPID::target_M_plus;
            TPID::par_m_minus=TPID::target_M_minus;
            TPID::par_b=parFile["TPID_PAR_B"];

            TPID::par_P_plus[0]=parFile["TPID_PAR_P_PLUS"]["X"];
            TPID::par_P_plus[1]=parFile["TPID_PAR_P_PLUS"]["Y"];
            TPID::par_P_plus[2]=parFile["TPID_PAR_P_PLUS"]["Z"];

            TPID::par_P_minus[0]=parFile["TPID_PAR_P_MINUS"]["X"];
            TPID::par_P_minus[1]=parFile["TPID_PAR_P_MINUS"]["Y"];
            TPID::par_P_minus[2]=parFile["TPID_PAR_P_MINUS"]["Z"];

            TPID::par_S_plus[0]=parFile["TPID_PAR_S_PLUS"]["X"];
            TPID::par_S_plus[1]=parFile["TPID_PAR_S_PLUS"]["Y"];
            TPID::par_S_plus[2]=parFile["TPID_PAR_S_PLUS"]["Z"];

            TPID::par_S_minus[0]=parFile["TPID_PAR_S_MINUS"]["X"];
            TPID::par_S_minus[1]=parFile["TPID_PAR_S_MINUS"]["Y"];
            TPID::par_S_minus[2]=parFile["TPID_PAR_S_MINUS"]["Z"];

            TPID::center_offset[0]=parFile["TPID_CENTER_OFFSET"]["X"];
            TPID::center_offset[1]=parFile["TPID_CENTER_OFFSET"]["Y"];
            TPID::center_offset[2]=parFile["TPID_CENTER_OFFSET"]["Z"];

            TPID::initial_lapse_psi_exponent=parFile["TPID_INITIAL_LAPSE_PSI_EXPONENT"];
            TPID::npoints_A=parFile["TPID_NPOINTS_A"];
            TPID::npoints_B=parFile["TPID_NPOINTS_B"];
            TPID::npoints_phi=parFile["TPID_NPOINTS_PHI"];

            TPID::give_bare_mass=parFile["TPID_GIVE_BARE_MASS"];
            TPID::initial_lapse=parFile["INITIAL_LAPSE"];
            TPID::solve_momentum_constraint=parFile["TPID_SOLVE_MOMENTUM_CONSTRAINT"];
            TPID::grid_setup_method=parFile["TPID_GRID_SETUP_METHOD"];
            TPID::verbose=parFile["TPID_VERBOSE"];
            TPID::adm_tol=parFile["TPID_ADM_TOL"];
            TPID::Newton_tol=parFile["TPID_NEWTON_TOL"];
            
            if (parFile.find("EXTRACTION_VAR_ID") != parFile.end()) {
                BHLOC::EXTRACTION_VAR_ID=parFile["EXTRACTION_VAR_ID"];
            }

            if (parFile.find("EXTRACTION_TOL") != parFile.end()) {
                BHLOC::EXTRACTION_TOL=parFile["EXTRACTION_TOL"];
            }


            vtu_len=MASSGRAV_VTU_FILE_PREFIX.size();
            chp_len=MASSGRAV_CHKPT_FILE_PREFIX.size();
            prf_len=MASSGRAV_PROFILE_FILE_PREFIX.size();


            GW::MASSGRAV_GW_NUM_RADAII=parFile["MASSGRAV_GW_NUM_RADAII"];
            GW::MASSGRAV_GW_NUM_LMODES=parFile["MASSGRAV_GW_NUM_LMODES"];

            for(unsigned int i=0;i<GW::MASSGRAV_GW_NUM_RADAII;i++)
                GW::MASSGRAV_GW_RADAII[i]=parFile["MASSGRAV_GW_RADAII"][i];

            for(unsigned int i=0;i<GW::MASSGRAV_GW_NUM_LMODES;i++)
                GW::MASSGRAV_GW_L_MODES[i]=parFile["MASSGRAV_GW_L_MODES"][i];


            if(parFile.find("MASSGRAV_USE_FD_GRID_TRANSFER")!=parFile.end())
            {
                massgrav::MASSGRAV_USE_FD_GRID_TRANSFER=parFile["MASSGRAV_USE_FD_GRID_TRANSFER"];
            }

            if(parFile.find("MASSGRAV_EH_COARSEN_VAL")!= parFile.end())
            {
                massgrav::MASSGRAV_EH_COARSEN_VAL = parFile["MASSGRAV_EH_COARSEN_VAL"];
            }

            if(parFile.find("MASSGRAV_EH_REFINE_VAL")!=parFile.end())
            {
                massgrav::MASSGRAV_EH_REFINE_VAL = parFile["MASSGRAV_EH_REFINE_VAL"];
            }

            if(parFile.find("MASSGRAV_REFINEMENT_MODE")!=parFile.end())
            {
                massgrav::MASSGRAV_REFINEMENT_MODE = static_cast<massgrav::RefinementMode>(parFile["MASSGRAV_REFINEMENT_MODE"]);
            }


        }

        par::Mpi_Bcast(&MASSGRAV_ELE_ORDER,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_IO_OUTPUT_FREQ,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_REMESH_TEST_FREQ,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_CHECKPT_FREQ,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_IO_OUTPUT_GAP,1,0,comm);

        par::Mpi_Bcast(&vtu_len,1,0,comm);
        par::Mpi_Bcast(&chp_len,1,0,comm);
        par::Mpi_Bcast(&prf_len,1,0,comm);

        par::Mpi_Bcast(&MASSGRAV_DENDRO_GRAIN_SZ,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_DENDRO_AMR_FAC,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_ASYNC_COMM_K,1,0,comm);
        par::Mpi_Bcast((int*)&MASSGRAV_REFINEMENT_MODE,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_GW_EXTRACT_FREQ,1,0,comm);

        char vtu_name[vtu_len+1];
        char chp_name[chp_len+1];
        char prf_name[prf_len+1];


        if(!rank)
        {
           for(unsigned int k=0;k<vtu_len;k++)
               vtu_name[k]=MASSGRAV_VTU_FILE_PREFIX[k];

            for(unsigned int k=0;k<chp_len;k++)
                chp_name[k]=MASSGRAV_CHKPT_FILE_PREFIX[k];

            for(unsigned int k=0;k<prf_len;k++)
                prf_name[k]=MASSGRAV_PROFILE_FILE_PREFIX[k];

            vtu_name[vtu_len]='\0';
            chp_name[chp_len]='\0';
            prf_name[prf_len]='\0';

        }


        MPI_Bcast(vtu_name,vtu_len+1,MPI_CHAR,0,comm);
        MPI_Bcast(chp_name,chp_len+1,MPI_CHAR,0,comm);
        MPI_Bcast(prf_name,prf_len+1,MPI_CHAR,0,comm);

        MASSGRAV_VTU_FILE_PREFIX=std::string(vtu_name);
        MASSGRAV_CHKPT_FILE_PREFIX=std::string(chp_name);
        MASSGRAV_PROFILE_FILE_PREFIX=std::string(prf_name);


        par::Mpi_Bcast(&MASSGRAV_RESTORE_SOLVER,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_ENABLE_BLOCK_ADAPTIVITY,1,0,comm);

        par::Mpi_Bcast(&MASSGRAV_USE_WAVELET_TOL_FUNCTION,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_WAVELET_TOL,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_WAVELET_TOL_MAX,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_WAVELET_TOL_FUNCTION_R0,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_WAVELET_TOL_FUNCTION_R1,1,0,comm);


        par::Mpi_Bcast(&MASSGRAV_CFL_FACTOR,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_LOAD_IMB_TOL,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_RK_TIME_BEGIN,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_RK_TIME_END,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_RK_TYPE,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_RK45_TIME_STEP_SIZE,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_RK45_DESIRED_TOL,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_DIM,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_MAXDEPTH,1,0,comm);

        MPI_Bcast(&(massgrav::BH1),sizeof(double)*10,MPI_BYTE,0,comm);
        MPI_Bcast(&(massgrav::BH2),sizeof(double)*10,MPI_BYTE,0,comm);

        par::Mpi_Bcast(&MASSGRAV_ID_TYPE,1,0,comm);

        par::Mpi_Bcast(&MASSGRAV_GRID_MIN_X,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_GRID_MAX_X,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_GRID_MIN_Y,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_GRID_MAX_Y,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_GRID_MIN_Z,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_GRID_MAX_Z,1,0,comm);

        par::Mpi_Bcast(&MASSGRAV_BLK_MIN_X,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_BLK_MIN_Y,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_BLK_MIN_Z,1,0,comm);

        par::Mpi_Bcast(&MASSGRAV_BLK_MAX_X,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_BLK_MAX_Y,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_BLK_MAX_Z,1,0,comm);

        MASSGRAV_OCTREE_MAX[0]=(double )(1u<<massgrav::MASSGRAV_MAXDEPTH);
        MASSGRAV_OCTREE_MAX[1]=(double )(1u<<massgrav::MASSGRAV_MAXDEPTH);
        MASSGRAV_OCTREE_MAX[2]=(double )(1u<<massgrav::MASSGRAV_MAXDEPTH);

        MASSGRAV_COMPD_MIN[0]=MASSGRAV_GRID_MIN_X;
        MASSGRAV_COMPD_MIN[1]=MASSGRAV_GRID_MIN_Y;
        MASSGRAV_COMPD_MIN[2]=MASSGRAV_GRID_MIN_Z;

        MASSGRAV_COMPD_MAX[0]=MASSGRAV_GRID_MAX_X;
        MASSGRAV_COMPD_MAX[1]=MASSGRAV_GRID_MAX_Y;
        MASSGRAV_COMPD_MAX[2]=MASSGRAV_GRID_MAX_Z;


        par::Mpi_Bcast(&ETA_CONST,1,0,comm);
        par::Mpi_Bcast(&ETA_R0,1,0,comm);
        par::Mpi_Bcast(&ETA_DAMPING,1,0,comm);
        par::Mpi_Bcast(&ETA_DAMPING_EXP,1,0,comm);

        par::Mpi_Bcast(&CHI_FLOOR,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_TRK0,1,0,comm);
        par::Mpi_Bcast(&DISSIPATION_TYPE,1,0,comm);
        par::Mpi_Bcast(&KO_DISS_SIGMA,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_ETA_R0,1,0,comm);

        MPI_Bcast(&(massgrav::MASSGRAV_LAMBDA),4,MPI_UNSIGNED,0,comm);
        MPI_Bcast(&(massgrav::MASSGRAV_LAMBDA_F),2,MPI_DOUBLE,0,comm);
        MPI_Bcast(&(massgrav::MASSGRAV_ETA_POWER),2,MPI_DOUBLE,0,comm);
        MPI_Bcast((massgrav::MASSGRAV_XI),3,MPI_INT,0,comm);


        par::Mpi_Bcast(&MASSGRAV_NUM_REFINE_VARS,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_NUM_EVOL_VARS_VTU_OUTPUT,1,0,comm);
        par::Mpi_Bcast(&MASSGRAV_NUM_CONST_VARS_VTU_OUTPUT,1,0,comm);


        if(MASSGRAV_NUM_REFINE_VARS>MASSGRAV_NUM_VARS){std::cout<<"Error[parameter file]: Number of refine variables should be less than number of MASSGRAV_NUM_VARS"<<std::endl; exit(0);}
        if(MASSGRAV_NUM_EVOL_VARS_VTU_OUTPUT>MASSGRAV_NUM_VARS){std::cout<<"Error[parameter file]: Number of evolution VTU variables should be less than number of MASSGRAV_NUM_VARS"<<std::endl; exit(0);}
        if(MASSGRAV_NUM_CONST_VARS_VTU_OUTPUT>MASSGRAV_CONSTRAINT_NUM_VARS){std::cout<<"Error[parameter file]: Number of constraint VTU variables should be less than number of MASSGRAV_CONSTRAINT_NUM_VARS"<<std::endl; exit(0);}

        par::Mpi_Bcast(MASSGRAV_REFINE_VARIABLE_INDICES,MASSGRAV_NUM_VARS,0,comm);
        par::Mpi_Bcast(MASSGRAV_VTU_OUTPUT_EVOL_INDICES,MASSGRAV_NUM_VARS,0,comm);
        par::Mpi_Bcast(MASSGRAV_VTU_OUTPUT_CONST_INDICES,MASSGRAV_CONSTRAINT_NUM_VARS,0,comm);

        par::Mpi_Bcast(&TPID::target_M_plus,1,0,comm);
        par::Mpi_Bcast(&TPID::target_M_minus,1,0,comm);
        par::Mpi_Bcast(&TPID::par_m_plus,1,0,comm);
        par::Mpi_Bcast(&TPID::par_m_minus,1,0,comm);
        par::Mpi_Bcast(&TPID::par_b,1,0,comm);

        MPI_Bcast(&(TPID::par_P_plus),3,MPI_DOUBLE,0,comm);
        MPI_Bcast(&(TPID::par_P_minus),3,MPI_DOUBLE,0,comm);
        MPI_Bcast(&(TPID::par_S_plus),3,MPI_DOUBLE,0,comm);
        MPI_Bcast(&(TPID::par_S_minus),3,MPI_DOUBLE,0,comm);

        MPI_Bcast(&(TPID::center_offset),3,MPI_DOUBLE,0,comm);

        par::Mpi_Bcast(&TPID::initial_lapse_psi_exponent,1,0,comm);
        par::Mpi_Bcast(&TPID::npoints_A,1,0,comm);
        par::Mpi_Bcast(&TPID::npoints_B,1,0,comm);
        par::Mpi_Bcast(&TPID::npoints_phi,1,0,comm);

        par::Mpi_Bcast(&TPID::give_bare_mass,1,0,comm);
        par::Mpi_Bcast(&TPID::initial_lapse,1,0,comm);
        par::Mpi_Bcast(&TPID::solve_momentum_constraint,1,0,comm);
        par::Mpi_Bcast(&TPID::grid_setup_method,1,0,comm);
        par::Mpi_Bcast(&TPID::verbose,1,0,comm);
        par::Mpi_Bcast(&TPID::adm_tol,1,0,comm);
        par::Mpi_Bcast(&TPID::Newton_tol,1,0,comm);
        
        
        par::Mpi_Bcast(&BHLOC::EXTRACTION_VAR_ID,1,0,comm);
        par::Mpi_Bcast(&BHLOC::EXTRACTION_TOL,1,0,comm);


        // gw extraction parameters.
        par::Mpi_Bcast(&GW::MASSGRAV_GW_NUM_RADAII,1,0,comm);
        par::Mpi_Bcast(&GW::MASSGRAV_GW_NUM_LMODES,1,0,comm);

        par::Mpi_Bcast(GW::MASSGRAV_GW_L_MODES,GW::MASSGRAV_GW_MAX_LMODES,0,comm);
        par::Mpi_Bcast(GW::MASSGRAV_GW_RADAII,GW::MASSGRAV_GW_MAX_RADAII,0,comm);
        par::Mpi_Bcast(&massgrav::MASSGRAV_USE_FD_GRID_TRANSFER,1,0,comm);
        par::Mpi_Bcast(&massgrav::MASSGRAV_VTU_Z_SLICE_ONLY,1,0,comm);
        par::Mpi_Bcast(&massgrav::MASSGRAV_EH_REFINE_VAL,1,0,comm);
        par::Mpi_Bcast(&massgrav::MASSGRAV_EH_COARSEN_VAL,1,0,comm);



    }



    void punctureData(const double xx1,const double yy1,const double zz1, double *var)
    {

        const double xx=GRIDX_TO_X(xx1);
        const double yy=GRIDY_TO_Y(yy1);
        const double zz=GRIDZ_TO_Z(zz1);

        /* Define the Levi-Cevita pseudo-tensor and Kroneckar delta */
        double epijk[3][3][3];
        int i,j,k;
        for (k=0;k<3;k++) {
            for (j=0;j<3;j++) {
                for (i=0;i<3;i++) {
                    epijk[k][j][i] = 0.0;
                }
            }
        }
        epijk[0][1][2] = 1.0;epijk[1][2][0] = 1.0;epijk[2][0][1] = 1.0;
        epijk[0][2][1] = -1.0;epijk[2][1][0] = -1.0;epijk[1][0][2] = -1.0;

        double deltaij[3][3];
        for (j=0;j<3;j++) {
            for (i=0;i<3;i++) {
                deltaij[j][i] = 0.0;
            }
        }

        deltaij[0][0] = 1.0;deltaij[1][1] = 1.0;deltaij[2][2] = 1.0;

        double x1,y1,z1,rv1;
        double x2,y2,z2,rv2;
        double vn1[3],vn2[3];

        double vpsibl;
        double v_u_corr,amp_capj,amp_capr,l_r,u0_j,u2_j,mu_j,p2_mu_j,v_u_j1;
        double v1,v2,v3,v4,vt1,vt2;

        int i1,i2,i3,i4;
        double amp_capp,u0_p,u2_p,mu_p,p2_mu_p;
        double v_u_p1,v_u_c1,v_u_j2,v_u_p2;
        double v_u_c2,vpsibl_u,vpsibl_u2;


        // bh 1
        double mass1 = BH1.getBHMass();
        double bh1x = BH1.getBHCoordX();
        double bh1y = BH1.getBHCoordY();
        double bh1z = BH1.getBHCoordZ();

        double vp1[3];
        vp1[0] = BH1.getVx();
        vp1[1] = BH1.getVy();
        vp1[2] = BH1.getVz();

        double vp1tot = sqrt( vp1[0]*vp1[0] + vp1[1]*vp1[1] + vp1[2]*vp1[2] );
        double spin1 = BH1.getBHSpin();
        double spin1_th = BH1.getBHSpinTheta();
        double spin1_phi = BH1.getBHSpinPhi();
        double vs1[3];

        vs1[0] = spin1*sin(spin1_th)*cos(spin1_phi);
        vs1[1] = spin1*sin(spin1_th)*sin(spin1_phi);
        vs1[2] = spin1*cos(spin1_th);

        // bh 2
        double mass2 = BH2.getBHMass();
        double bh2x =  BH2.getBHCoordX();
        double bh2y =  BH2.getBHCoordY();
        double bh2z =  BH2.getBHCoordZ();

        double vp2[3];
        vp2[0] = BH2.getVx();
        vp2[1] = BH2.getVy();
        vp2[2] = BH2.getVz();

        double vp2tot = sqrt( vp2[0]*vp2[0] + vp2[1]*vp2[1] + vp2[2]*vp2[2] );
        double spin2 = BH2.getBHSpin();
        double spin2_th = BH2.getBHSpinTheta();
        double spin2_phi = BH2.getBHSpinPhi();

        double vs2[3];
        vs2[0] = spin2*sin(spin2_th)*cos(spin2_phi);
        vs2[1] = spin2*sin(spin2_th)*sin(spin2_phi);
        vs2[2] = spin2*cos(spin2_th);


        // coordinates with respect to center of bh1
        x1 = xx - bh1x;
        y1 = yy - bh1y;
        z1 = zz - bh1z;

        //locating as a radial form
        rv1 = sqrt(x1*x1 + y1*y1 + z1*z1);
        vn1[0] = x1/rv1;
        vn1[1] = y1/rv1;
        vn1[2] = z1/rv1;

        //same as BH2
        x2 = xx - bh2x;
        y2 = yy - bh2y;
        z2 = zz - bh2z;

        rv2 = sqrt(x2*x2 + y2*y2 + z2*z2);
        vn2[0] = x2/rv2;
        vn2[1] = y2/rv2;
        vn2[2] = z2/rv2;

        //Initial data is related with the paper: http://arxiv.org/abs/0711.1165
        //Brill-Lindquist conformal factor
        vpsibl = 1.0 + mass1/(2.0*rv1);
        vpsibl = vpsibl + mass2/(2.0*rv2);

        v_u_corr = 0.0;
        // bh 1

        //For spinning puncture
        if ( fabs(spin1) > 1.e-6 ) {
            amp_capj = 4.0*spin1/(mass1*mass1);
            amp_capr = 2.0*rv1/mass1;
            l_r = 1.0/(1.0+amp_capr);
            u0_j = (l_r + l_r*l_r + l_r*l_r*l_r - 4.0*l_r*l_r*l_r*l_r + 2.0*l_r*l_r*l_r*l_r*l_r)/40.0;
            u2_j = -pow(l_r,5)/20.0;
            mu_j = vn1[0]*vs1[0];
            mu_j = mu_j + vn1[1]*vs1[1];
            mu_j = (mu_j + vn1[2]*vs1[2])/fabs(spin1);
            p2_mu_j = (3.0*mu_j*mu_j - 1.0)/2.0;
            v_u_j1 = amp_capj*amp_capj*(u0_j+u2_j*amp_capr*amp_capr*p2_mu_j);
            v_u_corr = v_u_corr + v_u_j1;
        }
        //For boosting puncture
        if (vp1tot > 1.e-6) {
            amp_capp = 2.0*vp1tot/mass1;
            amp_capr = 2.0*rv1/mass1;
            l_r = 1.0/(1.0 + amp_capr);
            u0_p = l_r - 2.0*l_r*l_r + 2.0*pow(l_r,3);
            u0_p = (u0_p - pow(l_r,4) + 0.20*pow(l_r,5))*(5.0/32.0);
            u2_p = 15.0*l_r + 132.0*l_r*l_r + 53.0*pow(l_r,3);
            u2_p = u2_p + 96.0*pow(l_r,4) + 82.0*pow(l_r,5);
            u2_p = u2_p + (84.0/amp_capr)*(pow(l_r,5)+log(l_r)/amp_capr);
            u2_p = (u2_p)/(80.0*amp_capr);
            mu_p =        vn1[0]*vp1[0]/vp1tot;
            mu_p = mu_p + vn1[1]*vp1[1]/vp1tot;
            mu_p = mu_p + vn1[2]*vp1[2]/vp1tot;
            p2_mu_p = (3.0*pow(mu_p,2) - 1.0)/2.0;
            v_u_p1 = pow(amp_capp,2)*(u0_p+u2_p*p2_mu_p);
            v_u_corr = v_u_corr + v_u_p1;
        }
        //For spinning boosted pucture
        if ( vp1tot > 1.e-6 && fabs(spin1) > 1.e-6 ) {
            v1 =      (vp1[1]*vs1[2]-vp1[2]*vs1[1])*vn1[0];
            v1 = v1 + (vp1[2]*vs1[0]-vp1[0]*vs1[2])*vn1[1];
            v1 = v1 + (vp1[0]*vs1[1]-vp1[1]*vs1[0])*vn1[2];
            v1 = v1*(16.0/pow(mass1,4))*rv1;

            amp_capr = 2.0*rv1/mass1;
            l_r = 1.0/(1.0 + amp_capr);

            v2 = 1.0 + 5.0*amp_capr + 10.0*pow(amp_capr,2);

            v_u_c1 = (v1*v2*pow(l_r,5))/80.0;
            v_u_corr = v_u_corr + v_u_c1;
        }
        // bh 2 same puncture as bh 1
        if ( fabs(spin2) > 1.e-6 ) {
            amp_capj = 4.0*spin2/(mass2*mass2);
            amp_capr = 2.0*rv2/mass2;
            l_r = 1.0/(1.0+amp_capr);
            u0_j = (l_r + l_r*l_r + l_r*l_r*l_r - 4.0*l_r*l_r*l_r*l_r + 2.0*l_r*l_r*l_r*l_r*l_r)/40.0;
            u2_j = -pow(l_r,5)/20.0;
            mu_j = vn2[0]*vs2[0];
            mu_j = mu_j + vn2[1]*vs2[1];
            mu_j = (mu_j + vn2[2]*vs2[2])/fabs(spin2);
            p2_mu_j = (3.0*mu_j*mu_j - 1.0)/2.0;
            v_u_j2 = amp_capj*amp_capj*(u0_j+u2_j*amp_capr*amp_capr*p2_mu_j);
            v_u_corr = v_u_corr + v_u_j2;
        }

        if (vp2tot > 1.e-6) {
            amp_capp = 2.0*vp2tot/mass2;
            amp_capr = 2.0*rv2/mass2;
            l_r = 1.0/(1.0 + amp_capr);
            u0_p = l_r - 2.0*l_r*l_r + 2.0*pow(l_r,3);
            u0_p = (u0_p - pow(l_r,4) + 0.20*pow(l_r,5))*(5.0/32.0);
            u2_p = 15.0*l_r + 132.0*l_r*l_r + 53.0*pow(l_r,3);
            u2_p = u2_p + 96.0*pow(l_r,4) + 82.0*pow(l_r,5);
            u2_p = u2_p + (84.0/amp_capr)*(pow(l_r,5)+log(l_r)/amp_capr);
            u2_p = (u2_p)/(80.0*amp_capr);
            mu_p =        vn2[0]*vp2[0]/vp2tot;
            mu_p = mu_p + vn2[1]*vp2[1]/vp2tot;
            mu_p = mu_p + vn2[2]*vp2[2]/vp2tot;
            p2_mu_p = (3.0*pow(mu_p,2) - 1.0)/2.0;
            v_u_p2 = pow(amp_capp,2)*(u0_p+u2_p*p2_mu_p);
            v_u_corr = v_u_corr + v_u_p2;
        }

        if ( vp2tot > 1.e-6 && fabs(spin2) > 1.e-6 ) {
            v1 =      (vp2[1]*vs2[2]-vp2[2]*vs2[1])*vn2[0];
            v1 = v1 + (vp2[2]*vs2[0]-vp2[0]*vs2[2])*vn2[1];
            v1 = v1 + (vp2[0]*vs2[1]-vp2[1]*vs2[0])*vn2[2];
            v1 = v1*(16.0/pow(mass2,4))*rv2;

            amp_capr = 2.0*rv2/mass2;
            l_r = 1.0/(1.0 + amp_capr);

            v2 = 1.0 + 5.0*amp_capr + 10.0*pow(amp_capr,2);

            v_u_c2 = (v1*v2*pow(l_r,5))/80.0;
            v_u_corr = v_u_corr + v_u_c2;
        }

        // vpsibl_u will be used for the conformal factor,
        vpsibl_u  = vpsibl + v_u_corr;
        // vpsibl_u2 is for the Aij terms...
        // ! since the corrections are first order...
        // ! adding half of the correction seems to give the best results...
        // ! update - do a fit for spin = 0.6...
        vpsibl_u2 = vpsibl + v_u_corr;

        var[VAR::U_ALPHA] = 1.0/(vpsibl_u*vpsibl_u);
        //std::cout<<"Alpha: "<<u[U_ALPHA]<<" vpsibl_u: "<< vpsibl_u<<std::endl;
        var[VAR::U_ALPHA] = std::max(var[VAR::U_ALPHA], CHI_FLOOR);

        v2 = 1.0/pow(vpsibl_u,4);
        var[VAR::U_CHI] = v2;

        if(var[VAR::U_CHI]<CHI_FLOOR)
            var[VAR::U_CHI]=CHI_FLOOR;

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

        var[VAR::U_SYMGT0] = 1.0; //XX
        var[VAR::U_SYMGT1] = 0.0; //XY
        var[VAR::U_SYMGT2] = 0.0; //XZ
        var[VAR::U_SYMGT3] = 1.0; //YY
        var[VAR::U_SYMGT4] = 0.0; //YZ
        var[VAR::U_SYMGT5] = 1.0; //ZZ

        for (i1=0;i1<3;i1++) {
            for (i2=0;i2<3;i2++) {
                // first BH
                v2 = 0.0;
                for (i3=0;i3<3;i3++) {
                    for (i4=0;i4<3;i4++) {
                        vt1 = epijk[i1][i3][i4]*vs1[i3]*vn1[i4]*vn1[i2];
                        vt2 = epijk[i2][i3][i4]*vs1[i3]*vn1[i4]*vn1[i1];
                        v2 = v2 + vt1 + vt2;
                    }
                }

                v3 = vp1[i1]*vn1[i2] + vp1[i2]*vn1[i1];
                vt1 = 0.0;
                for (i3=0;i3<3;i3++) {
                    vt1 = vt1 + vp1[i3]*vn1[i3];
                }
                vt1 = vt1*(vn1[i1]*vn1[i2] - deltaij[i1][i2]);
                v3 = v3 + vt1;

                v1 = 3.0/(pow(vpsibl_u2,6)*pow(rv1,3));
                v4 = v1*(v2+(rv1/2.0)*v3);

                // second BH
                v2 = 0.0;
                for (i3=0;i3<3;i3++) {
                    for (i4=0;i4<3;i4++) {
                        vt1 = epijk[i1][i3][i4]*vs2[i3]*vn2[i4]*vn2[i2];
                        vt2 = epijk[i2][i3][i4]*vs2[i3]*vn2[i4]*vn2[i1];
                        v2 = v2 + vt1 + vt2;
                    }
                }

                v3 = vp2[i1]*vn2[i2] + vp2[i2]*vn2[i1];
                vt1 = 0.0;
                for (i3=0;i3<3;i3++) {
                    vt1 = vt1 + vp2[i3]*vn2[i3];
                }
                vt1 = vt1*(vn2[i1]*vn2[i2] - deltaij[i1][i2]);
                v3 = v3 + vt1;

                v1 = 3.0/(pow(vpsibl_u2,6)*pow(rv2,3));
                v4 = v4 + v1*(v2+(rv2/2.0)*v3);

                if ( i1 == 0 && i2 == 0 ) {
                    var[VAR::U_SYMAT0] = v4; //XX
                } else if ( i1 == 0 && i2 == 1 ) {
                    var[VAR::U_SYMAT1] = v4; //XY
                } else if ( i1 == 0 && i2 == 2 ) {
                    var[VAR::U_SYMAT2] = v4; //XZ
                } else if ( i1 == 1 && i2 == 1 ) {
                    var[VAR::U_SYMAT3] = v4; //YY
                } else if ( i1 == 1 && i2 == 2 ) {
                    var[VAR::U_SYMAT4] = v4; //YZ
                } else if ( i1 == 2 && i2 == 2 ) {
                    var[VAR::U_SYMAT5] = v4; //ZZ
                }

            }
        }

    }
    
namespace trumpet_data { 

    const unsigned int KMAX = 5000;
    double dxsav;
    static double *xp;
    static double *yp[2];

    void trumpetData(const double xx1,const double yy1,const double zz1, double *var)
    {

        const double xx=GRIDX_TO_X(xx1);
        const double yy=GRIDY_TO_Y(yy1);
        const double zz=GRIDZ_TO_Z(zz1);

        double eps = 1.e-10; // tolerance on the integrator  
        double h1 = -1.e-5 ; // starting step off the critical point  
        double hmin = 1.e-25;// min allowed stepsize for adaptive integrator  
        
        double  bh_mass = 7.5 ; 
        double  alpha_c = 0.16227766016837933200 ; 
        double  C_sq = 1.5543095902183302040 ; 
        double  bigR_c = 1.5405694150420948330 * bh_mass ; 
        double  r_c = 0.30405997036 * bh_mass ; 

        static bool firstcall = true;

        const int neq = 2;
        static double *alpha0;
        static double *bigR0;
        static double *r_iso;

        static int np;

        const double third = 1.0 / 3.0 ;

        if (firstcall == true) {
            // solve ODE system 
            std::cout<<"trumpetData:  first call. Solve the ODE"<<std::endl;

            xp = new double [KMAX];
            yp[0] = new double [KMAX];
            yp[1] = new double [KMAX];
            alpha0 = new double [KMAX]; 
            bigR0 = new double [KMAX]; 
            r_iso = new double [KMAX]; 

            // integrate inward from r_C to r=0.
            double rmin = 0.0;   // min value for inward r integration
            double rmax = 1000.0;   // max value for outward r integration
            double rstart, ystart[2], dydr[2];
            bndcnd(h1, rstart, ystart, dydr);
            double rend = rmin + fabs(h1) ; 
            int nok, nbad;
            int kount;
            odeint(ystart,neq,rstart,rend,eps,h1,hmin,&nok,&nbad,derivs,rkqs,kount);

            int kountin = kount;
            for ( unsigned int i=0; i<kountin; i++ )
            {  
                r_iso [i] = xp[kountin - 1 - i] ; 
                alpha0[i] = yp[0][kountin - 1 - i] ; 
                bigR0 [i] = yp[1][kountin - 1 - i] ; 
                //std::cout<<"<in> xp = "<<xp[i]<<", i="<<i<<", alpha0 = "<<yp[0][i]<<", bigR0 = "<<yp[1][i]<<std::endl;
            }

            // integrate outwards from r_c to large r 
            h1 = fabs(h1) ; 
            rend = rmax ;
      
            std::cout<<"integrating outwards ... "<<std::endl;
            bndcnd(h1, rstart, ystart, dydr) ; 
            odeint(ystart,neq,rstart,rend,eps,h1,hmin,&nok,&nbad,derivs,rkqs,kount) ; 

            int kountout = kount ; 
            for ( unsigned int i=0; i<kountout; i++ )
            {  
              r_iso [i+kountin] = xp[i] ; 
              alpha0[i+kountin] = yp[0][i] ; 
              bigR0 [i+kountin] = yp[1][i] ; 
              //std::cout<<"<out> xp = "<<xp[i]<<", alpha0 = "<<yp[0][i]<<", bigR0 = "<<yp[1][i]<<std::endl;
            }
            np = kountin + kountout ;

            firstcall = false;
        }

        int nrby_indx = np / 2 ;
        double ax_eps = 1.0e-5 ;

        double tenh1 = 10.0 * fabs(h1) ;
        double rbar = sqrt( xx*xx + yy*yy + zz*zz ) ;
        if ( fabs(xx)<tenh1 && fabs(yy)<tenh1 && fabs(zz)<tenh1 ) { 
            rbar = tenh1 ; 
        }
        double alpha = interpolation4(r_iso,alpha0,np,rbar,&nrby_indx) ;  
        double bigR  = interpolation4(r_iso,bigR0,np,rbar,&nrby_indx) ;  

                if ( fabs(xx)<ax_eps && fabs(yy)<ax_eps && fabs(zz)<ax_eps)
                { rbar = sqrt( xx*xx + yy*yy + zz*zz + ax_eps*ax_eps) ; 
                } 

                double f0 = 1.0 - 2.0 * bh_mass / bigR ; 

                double Rsq_dalpha_dR = 4.0 * (   alpha*alpha - f0 
                                        - 0.5 * bh_mass / bigR ) 
                                    / ( alpha*alpha - 2.0*alpha - f0 ) 
                                    * bigR  ;

                double tmp_sqrt = sqrt( alpha*alpha - f0 ) ; 

                double tmp_chi = (rbar*rbar) / (bigR*bigR) ; 
        
                double tmp_trK =   2.0 * tmp_sqrt / bigR 
                          + ( alpha * Rsq_dalpha_dR - bh_mass ) 
                            / ( tmp_sqrt * bigR * bigR ) ; 

                double tmp_beta = tmp_sqrt / bigR ; 

                double tmp_Atilde = (   (   alpha * Rsq_dalpha_dR 
                                   - bh_mass 
                                 ) / tmp_sqrt 
                               - bigR * tmp_sqrt 
                             ) / bigR / bigR ;  






        var[VAR::U_ALPHA] = alpha;
        var[VAR::U_ALPHA] = std::max(var[VAR::U_ALPHA], CHI_FLOOR);

        var[VAR::U_CHI] = tmp_chi;

        if(var[VAR::U_CHI]<CHI_FLOOR)
            var[VAR::U_CHI]=CHI_FLOOR;

        var[VAR::U_K] = tmp_trK;

        var[VAR::U_BETA0] = xx * tmp_beta;
        var[VAR::U_BETA1] = yy * tmp_beta;
        var[VAR::U_BETA2] = zz * tmp_beta;

        var[VAR::U_GT0] = 0.0;
        var[VAR::U_GT1] = 0.0;
        var[VAR::U_GT2] = 0.0;

        var[VAR::U_B0] = 0.0;
        var[VAR::U_B1] = 0.0;
        var[VAR::U_B2] = 0.0;

        var[VAR::U_SYMGT0] = 1.0; //XX
        var[VAR::U_SYMGT1] = 0.0; //XY
        var[VAR::U_SYMGT2] = 0.0; //XZ
        var[VAR::U_SYMGT3] = 1.0; //YY
        var[VAR::U_SYMGT4] = 0.0; //YZ
        var[VAR::U_SYMGT5] = 1.0; //ZZ

        var[VAR::U_SYMAT0] = tmp_Atilde*( (xx/rbar)*(xx/rbar) - third);
        var[VAR::U_SYMAT1] = tmp_Atilde*xx * yy /(rbar * rbar);
        var[VAR::U_SYMAT2] = tmp_Atilde*xx * zz /(rbar * rbar);
        var[VAR::U_SYMAT3] = tmp_Atilde*( (yy/rbar)*(yy/rbar) - third);
        var[VAR::U_SYMAT4] = tmp_Atilde*yy * zz /(rbar * rbar);
        var[VAR::U_SYMAT5] = tmp_Atilde*( (zz/rbar)*(zz/rbar) - third);


    }

/*---------------------------------------------------------------------
 *  
 *  
 *  
 *---------------------------------------------------------------------*/ 
void bndcnd(double h, double &x, double y[], double dydx[])
{
    
    double bh_mass = 7.5 ; 
    double alpha_c = 0.16227766016837933200 ; 
    double bigR_c = 1.5405694150420948330 * bh_mass ; 
    double r_c = 0.30405997036 * bh_mass ; 

    double alpha = alpha_c ; 
    double bigR  = bigR_c ; 
    double r = r_c;

    double bigRsq = bigR * bigR;
    double alphasq = alpha * alpha;

    double dbigR_dr = alpha_c * bigR_c / r_c ; 

    double df0_dr = 2.0 * bh_mass / bigRsq * dbigR_dr;

    double tmp1 = 1.5 * bh_mass / bigRsq * dbigR_dr;

    double dalpha_dr =   0.25 * (r*df0_dr + 8.0*alphasq) / (alpha-1.0 )
                       - 0.25 * sqrt(   pow((r*df0_dr + 8.0*alphasq),2)
                                      + 32.0*alpha*(1.0-alpha)*r * tmp1
                                    ) / ( alpha - 1.0 );
    
    dalpha_dr = dalpha_dr/r;
    
    double ddbigR_drdr  =   ( dalpha_dr * bigR + alpha * dbigR_dr ) / r
                          - alpha * bigR / (r*r);
    
    double ddf0_drdr = - 4.0 * bh_mass * dbigR_dr * dbigR_dr / (bigR*bigRsq)
                       + 2.0 * bh_mass * ddbigR_drdr / bigRsq;
    
    double tmp2 = - 3.00 * bh_mass * dbigR_dr * dbigR_dr / (bigR*bigRsq)
                  + 1.50 * bh_mass * ddbigR_drdr / bigRsq;
    
    double ddalpha_drdr = - 2.0 * pow((r*dalpha_dr),3)
                          - 3.0 * r*dalpha_dr
                                * (   2.0 * r*dalpha_dr * (alpha - 1.0)
                                    - r*df0_dr )
                          + (r*dalpha_dr) * r*r * ddf0_drdr
                          +   (8.0 * r*dalpha_dr + 4.0*alpha)
                            * ( 2.0*alpha * r*dalpha_dr - r*tmp1)
                          + 4.0*alpha * (   2.0 *pow((r*dalpha_dr),2)
                                          - r*r * tmp2 );
        
    ddalpha_drdr =   ddalpha_drdr
                   / (   6.0 * pow(r,3) * (alpha-1.0) * dalpha_dr
                       - 2.0 * pow(r,3) * df0_dr
                       - 8.0 * r*r * alphasq );
    
    x = r + h;
    
    y[0] = alpha + h * dalpha_dr + 0.5 * h*h * ddalpha_drdr;
    y[1] =  bigR + h *  dbigR_dr + 0.5 * h*h *  ddbigR_drdr;
    
    dydx[0] = dalpha_dr + h * ddalpha_drdr;
    dydx[1] =  dbigR_dr + h *  ddbigR_drdr;
}
   


/*---------------------------------------------------------------------
 *  
 *  
 *  
 *---------------------------------------------------------------------*/ 
void derivs(double x, double y[], double dydx[])
{
    
    double alpha = y[0];
    double bigR = y[1];
    double r = x;
   
    double bh_mass = 7.5 ; 

    double alpha_sq = alpha * alpha;
    
    double alpha_sq_minus_1 = alpha_sq - 1.0;
    double M_over_R = bh_mass / bigR;
    
    double dalpha_dr = 4.0 * alpha * (alpha_sq_minus_1 + 1.5*M_over_R)
                           / (alpha_sq_minus_1 - 2.0*alpha + 2.0*M_over_R)
                           / r;
    
    double dbigR_dr  = alpha*bigR/r;
    
    dydx[0] = dalpha_dr;
    dydx[1] = dbigR_dr;
    
}
    
    
    
/*---------------------------------------------------------------------
 *  
 *  
 *  
 *---------------------------------------------------------------------*/ 
void hunt(double xx[], int n, double x, int *jlo)
{
    unsigned long jm,jhi,inc;
    int ascnd;
    
    ascnd=(xx[n-1] > xx[0]);
    //if (*jlo < 0 || *jlo > n) 
    if (*jlo < 0 || *jlo > n-1) {
        //*jlo=0;
        *jlo=-1;
        //jhi=n-1;
        jhi=n;
    } else {
        inc=1;
            if (x >= xx[*jlo] == ascnd) {
                if (*jlo == n-1) return;
                jhi=(*jlo)+1;
                while (x >= xx[jhi] == ascnd) {
                    *jlo=jhi;
                    inc += inc;
                    jhi=(*jlo)+inc;
                    if (jhi > n-1) {
                        jhi=n;
                    break;
                }
            }
        } else {
            if (*jlo == 0) {
                *jlo=-1;
                return;
            }
            jhi=(*jlo)--;
            while (x < xx[*jlo] == ascnd) {
                jhi=(*jlo);
                inc <<= 1;
                if (inc >= jhi) {
                    *jlo=-1;
                    break;
                }
                else *jlo=jhi-inc;
            }
        }
    }
    while (jhi-(*jlo) != 1) {
        jm=(jhi+(*jlo)) >> 1;
        if (x > xx[jm] == ascnd)
            *jlo=jm;
        else
            jhi=jm;
    }
} /* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */
   



/*---------------------------------------------------------------------
 *  
 *  
 *  
 *---------------------------------------------------------------------*/ 
void rkck(   double y[], double dydx[], int n, double x, double h, 
             double yout[], double yerr[], 
             void (*derivs)(double, double [], double [])   )
{
    int i;
    static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
        b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
        b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
        b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
        b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
        c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
        dc5 = -277.0/14336.0;
    double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
        dc4=c4-13525.0/55296.0,dc6=c6-0.25;

    double *ak2 = new double [n];
    double *ak3= new double [n];
    double *ak4= new double [n];
    double *ak5= new double [n];
    double *ak6= new double [n];
    double *ytemp= new double [n];
 
    for (i=0;i<n;i++)
        ytemp[i]=y[i]+b21*h*dydx[i];
    (*derivs)(x+a2*h,ytemp,ak2);
    for (i=0;i<n;i++)
        ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
    (*derivs)(x+a3*h,ytemp,ak3);
    for (i=0;i<n;i++)
        ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
    (*derivs)(x+a4*h,ytemp,ak4);
    for (i=0;i<n;i++)
        ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
    (*derivs)(x+a5*h,ytemp,ak5);
    for (i=0;i<n;i++)
        ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
    (*derivs)(x+a6*h,ytemp,ak6);
    for (i=0;i<n;i++)
        yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
    for (i=0;i<n;i++)
        yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
    delete [] ytemp;
    delete [] ak6;
    delete [] ak5;
    delete [] ak4;
    delete [] ak3;
    delete [] ak2;
} /* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */




/*---------------------------------------------------------------------
 *  
 *  
 *  
 *---------------------------------------------------------------------*/ 
void rkqs(   double y[], double dydx[], int n, double *x, double htry, 
             double eps, double yscal[], double *hdid, double *hnext,
             void (*derivs)(double, double [], double [])   )
{
    const double SAFETY = 0.9;
    const double PGROW = -0.2;
    const double PSHRNK = -0.25;
    const double ERRCON = 1.89e-4;
    
    void rkck( double y[], double dydx[], int n, double x, double h,
               double yout[], double yerr[], 
               void (*derivs)(double, double [], double []) );
    int i;
    double errmax, h, xnew;

    double *yerr = new double [n];
    double *ytemp = new double [n];
    h=htry;
    for (;;) {
        rkck(y,dydx,n,*x,h,ytemp,yerr,derivs);
        errmax=0.0;
        for (i=0;i<n;i++) errmax=std::max(errmax,fabs(yerr[i]/yscal[i]));
        errmax /= eps;
        if (errmax > 1.0) {
            h=SAFETY*h*pow(errmax,PSHRNK);
            if (h < 0.1*h) h *= 0.1;
            xnew=(*x)+h;
            if (xnew == *x) std::cerr<<"stepsize underflow in rkqs"<<std::endl;
            continue;
        } else {
            if (errmax > ERRCON) *hnext=SAFETY*h*pow(errmax,PGROW);
            else *hnext=5.0*h;
            *x += (*hdid=h);
            for (i=0;i<n;i++) y[i]=ytemp[i];
            break;
        }
    }
    delete [] ytemp;
    delete [] yerr;
} /* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */
    
    
/*---------------------------------------------------------------------
 *  
 *  
 *  
 *---------------------------------------------------------------------*/ 
void odeint(  double ystart[], int nvar, double x1, double x2, 
              double eps, double h1, double hmin, int *nok, int *nbad, 
              void (*derivs)( double, double [], double [] ),
              void  (*rkqs) ( double [], double [], int, double *, 
                              double, double, double [], double *, 
                              double *, 
                              void (*)( double, double [], double [] ) ), 
              int kount )
{
    int nstp,i;
    double xsav,x,hnext,hdid,h;
    const int MAXSTP = 10000;
    const double TINY = 1.0e-30;

    double *yscal = new double [nvar];
    double *y = new double [nvar];
    double *dydx = new double [nvar];
    x=x1;
    h=copysign(h1,x2-x1);
    *nok = (*nbad) = kount = 0;
    for (i=0;i<nvar;i++) y[i]=ystart[i];
    if (KMAX > 0) xsav=x-dxsav*2.0;
    for (nstp=0;nstp<MAXSTP;nstp++) {
        //std::cout<<"odeint: nstp="<<nstp<<", kount="<<kount<<std::endl;
        (*derivs)(x,y,dydx);
        for (i=0;i<nvar;i++)
            yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
        if (KMAX > 0 && kount < KMAX-1 && fabs(x-xsav) > fabs(dxsav)) {
            xp[kount]=x;
            for (i=0;i<nvar;i++) yp[i][kount]=y[i];
            xsav=x;
            kount++;
        }
        if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
        (*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
        if (hdid == h) ++(nok); else ++(nbad);
        if ((x-x2)*(x2-x1) >= 0.0) {
            for (i=0;i<nvar;i++) ystart[i]=y[i];
            if (KMAX) {
                xp[kount]=x;
                for (i=0;i<nvar;i++) yp[i][kount]=y[i];
                kount++;
            }
            delete [] dydx;
            delete [] y;
            delete [] yscal;
            return;
        }
        if (fabs(hnext) <= hmin) std::cerr<<"Step size too small in odeint"<<std::endl;
        h=hnext;
    }
    std::cerr<<"Too many steps in routine odeint"<<std::endl;
} 
/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */




/*---------------------------------------------------------------------
 *  
 *  
 *  
 *---------------------------------------------------------------------*/ 
double interpolation3( double xp[], double yp[], int np, double xb, 
                       int *n_nearest_pt )
{
    int  k ;      /* index of 1st point */
    int  m = 4;   /* degree of interpolation */
    double y ;    /* intermediate value */

    hunt ( xp, np, xb, n_nearest_pt );

    k = std::min  ( std::max( (*n_nearest_pt) - (m-1)/2 ,
                      1
                    ) ,
               np+1-m
             );
    
    double DBL_EPSILON = 1.e-12 ; 

    if ( xb==xp[k] || xb==xp[k+1] || xb==xp[k+2] || xb==xp[k+3] )
    { xb += 3.0 * DBL_EPSILON; }

    y=     (xb-xp[k+1])
         * (xb-xp[k+2])
         * (xb-xp[k+3])
         * yp[k  ] / (   (xp[k  ]-xp[k+1])
                       * (xp[k  ]-xp[k+2])
                       * (xp[k  ]-xp[k+3])  )

       +   (xb-xp[k  ])
         * (xb-xp[k+2])
         * (xb-xp[k+3])
         * yp[k+1] / (   (xp[k+1]-xp[k  ])
                       * (xp[k+1]-xp[k+2])
                       * (xp[k+1]-xp[k+3])  )

       +   (xb-xp[k  ])
         * (xb-xp[k+1])
         * (xb-xp[k+3])
         * yp[k+2] / (   (xp[k+2]-xp[k  ])
                       * (xp[k+2]-xp[k+1])
                       * (xp[k+2]-xp[k+3])  )

       +   (xb-xp[k  ])
         * (xb-xp[k+1])
         * (xb-xp[k+2])
         * yp[k+3] / (   (xp[k+3]-xp[k  ])
                       * (xp[k+3]-xp[k+1])
                       * (xp[k+3]-xp[k+2])  );

    return (y);
}



/*---------------------------------------------------------------------
 *  
 *  
 *  
 *---------------------------------------------------------------------*/ 
double interpolation4( double xx[], double yy[], int np, double xb, 
                       int *n_nearest_pt  )
{ 
    int  k    ;   /* index of 1st point */
    int  m = 5;   /* degree of interpolation */
    double y;     /* intermediate value */

    hunt ( xx, np, xb, n_nearest_pt );
    
    k = std::min ( std::max ( (*n_nearest_pt) - (m-1)/2 ,
                      1
                    ) , 
               np+1-m
             );

#if 0
    if ( fabs( xb - 1.0 ) < 1.e-13 )
    { // xb is 1.0 so just return with the corresponding y value -- no interp necessary 
        y = yy[np] ;
        return(y) ;
    }
    if ( fabs( xb ) < 1.e-13 )
    { // xb is zero so just return with the corresponding y value -- no interp necessary 
        y = yy[0] ;
        return(y) ;
    }
#endif 


    double DBL_EPSILON = 1.e-12 ; 
    
    if ( xb==xx[k] || xb==xx[k+1] || xb==xx[k+2] || xb==xx[k+3] || xb==xx[k+4] )
    { xb += 3.0 * DBL_EPSILON; }
    
    double xtmp0 = xb - xx[k  ] ;
    double xtmp1 = xb - xx[k+1] ;
    double xtmp2 = xb - xx[k+2] ;
    double xtmp3 = xb - xx[k+3] ;
    double xtmp4 = xb - xx[k+4] ;
    double xdiff01 = xx[k  ] - xx[k+1] ;
    double xdiff02 = xx[k  ] - xx[k+2] ;
    double xdiff03 = xx[k  ] - xx[k+3] ;
    double xdiff04 = xx[k  ] - xx[k+4] ;
    double xdiff12 = xx[k+1] - xx[k+2] ;
    double xdiff13 = xx[k+1] - xx[k+3] ;
    double xdiff14 = xx[k+1] - xx[k+4] ;
    double xdiff23 = xx[k+2] - xx[k+3] ;
    double xdiff24 = xx[k+2] - xx[k+4] ;
    double xdiff34 = xx[k+3] - xx[k+4] ;
  
    y=     xtmp1
         * xtmp2
         * xtmp3
         * xtmp4
         * yy[k  ] / (   xdiff01
                       * xdiff02
                       * xdiff03
                       * xdiff04  )
 
       -   xtmp0
         * xtmp2
         * xtmp3
         * xtmp4
         * yy[k+1] / (   xdiff01
                       * xdiff12
                       * xdiff13
                       * xdiff14  )

       +   xtmp0
         * xtmp1
         * xtmp3
         * xtmp4
         * yy[k+2] / (   xdiff02
                       * xdiff12
                       * xdiff23
                       * xdiff24  )

       -   xtmp0
         * xtmp1
         * xtmp2
         * xtmp4
         * yy[k+3] / (   xdiff03
                       * xdiff13
                       * xdiff23
                       * xdiff34  )
       +   xtmp0
         * xtmp1
         * xtmp2
         * xtmp3
         * yy[k+4] / (   xdiff04
                       * xdiff14
                       * xdiff24
                       * xdiff34  );
    
    return (y);
}



    }  // end of namespace trumpet_data
   

    void KerrSchildData(const double xx1,const double yy1,const double zz1, double *var)
    {

        const double xx=GRIDX_TO_X(xx1);
        const double yy=GRIDY_TO_Y(yy1);
        const double zz=GRIDZ_TO_Z(zz1);

        // parameters for the BH (mass, location, spin parameter)
        double M = BH1.getBHMass();
        double bh1x = BH1.getBHCoordX();
        double bh1y = BH1.getBHCoordY();
        double bh1z = BH1.getBHCoordZ();
        double spin1 = BH1.getBHSpin();

        // coordinates relative to the center of the BH
        double x = xx - bh1x;
        double y = yy - bh1y;
        double z = zz - bh1z;

        //locating as a radial form
        double r = sqrt(x*x + y*y + z*z);

	    //HL : Angular momentum parameter will be added as param file after testing
	    double a = spin1;

	    double gtd[3][3], Atd[3][3];
	    double alpha, Gamt[3];
	    double Chi, TrK, Betau[3];

        #include "ks_vars.cpp"
        #include "ksinit.cpp"
 	
        var[VAR::U_ALPHA] = alpha;
        var[VAR::U_CHI] = Chi;
        var[VAR::U_K] = TrK;

        var[VAR::U_BETA0] = Betau[0];
        var[VAR::U_BETA1] = Betau[1];
        var[VAR::U_BETA2] = Betau[2];

	    var[VAR::U_GT0] = Gamt[0];
	    var[VAR::U_GT1] = Gamt[1];
	    var[VAR::U_GT2] = Gamt[2];

        var[VAR::U_B0] = 0.0;
        var[VAR::U_B1] = 0.0;
        var[VAR::U_B2] = 0.0;
	
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

        //std::cout<<"KS init data: (x,y,z) = ( "<<x<<", "<<y<<", "<<z<<"), alpha = "<<alpha<<std::endl;

        #if 0
            //MASSGRAV vars for Kerr-Schild
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


    void FlatMinkowski(const double xx1,const double yy1,const double zz1, double *var)
    {

	    double gtd[3][3], Atd[3][3];
	    double alpha, Gamt[3];
	    double Chi, TrK, Betau[3];
 	
        var[VAR::U_ALPHA] = 1.0;
        var[VAR::U_CHI] = 0.0;
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
	
	    var[VAR::U_SYMGT0] = 1.0;
	    var[VAR::U_SYMGT1] = 0.0;
	    var[VAR::U_SYMGT2] = 0.0;
	    var[VAR::U_SYMGT3] = 1.0;
	    var[VAR::U_SYMGT4] = 0.0;
	    var[VAR::U_SYMGT5] = 1.0;

	    var[VAR::U_SYMAT0] = 0.0;
	    var[VAR::U_SYMAT1] = 0.0;
	    var[VAR::U_SYMAT2] = 0.0;
	    var[VAR::U_SYMAT3] = 0.0;
	    var[VAR::U_SYMAT4] = 0.0;
	    var[VAR::U_SYMAT5] = 0.0;
    }


    void noiseData(const double xx1,const double yy1,const double zz1, double *var)
    {

        //const double xx=GRIDX_TO_X(xx1);
        //const double yy=GRIDY_TO_Y(yy1);
        //const double zz=GRIDZ_TO_Z(zz1);

        // call random number generator between -1 and 1) 
        double random_variable[30] ; 
        int i ; 
        for ( i=0; i<30; i++ ) 
        {
          random_variable[i] = 2.0*rand()/((double)RAND_MAX) - 1.0; 
        }

        // set a (uniform) amplitude for the noise 
        double noise_amp = 1.0e-8 ; 

        var[VAR::U_ALPHA] = 1.0 + noise_amp * random_variable[0] ;  
        
        var[VAR::U_CHI] = 1.0 + noise_amp * random_variable[1] ;

        var[VAR::U_K] = noise_amp * random_variable[2] ;

        var[VAR::U_GT0] = noise_amp * random_variable[3] ;
        var[VAR::U_GT1] = noise_amp * random_variable[4] ;
        var[VAR::U_GT2] = noise_amp * random_variable[5] ;

        var[VAR::U_BETA0] = noise_amp * random_variable[6] ;
        var[VAR::U_BETA1] = noise_amp * random_variable[7] ;
        var[VAR::U_BETA2] = noise_amp * random_variable[8] ;

        var[VAR::U_B0] = noise_amp * random_variable[9] ;
        var[VAR::U_B1] = noise_amp * random_variable[10] ;
        var[VAR::U_B2] = noise_amp * random_variable[11] ;

        var[VAR::U_SYMGT0] = 1.0 + noise_amp * random_variable[12] ; //XX
        var[VAR::U_SYMGT1] =       noise_amp * random_variable[13] ; //XY
        var[VAR::U_SYMGT2] =       noise_amp * random_variable[14] ; //XZ
        var[VAR::U_SYMGT3] = 1.0 + noise_amp * random_variable[15] ; //YY
        var[VAR::U_SYMGT4] =       noise_amp * random_variable[16] ; //YZ
        var[VAR::U_SYMGT5] = 1.0 + noise_amp * random_variable[17] ; //ZZ

        var[VAR::U_SYMAT0] = noise_amp * random_variable[18] ; //XX
        var[VAR::U_SYMAT1] = noise_amp * random_variable[19] ; //XY
        var[VAR::U_SYMAT2] = noise_amp * random_variable[20] ; //XZ
        var[VAR::U_SYMAT3] = noise_amp * random_variable[21] ; //YY
        var[VAR::U_SYMAT4] = noise_amp * random_variable[22] ; //YZ
        var[VAR::U_SYMAT5] = noise_amp * random_variable[23] ; //ZZ

    }


    void fake_initial_data(double xx1, double yy1, double zz1, double *u)
    {

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

        u[VAR::U_SYMAT0] = exp(-4.0*cos(x)*sin(y))*(cos(x) -0.3333333333*exp(4.0*cos(x)*sin(y)) *(1.0+0.2*sin(x))*(5.0*exp(-4.0*cos(x)*sin(y)) /(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y)) /(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y) +5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
        u[VAR::U_SYMAT1] = 1.0 + x*z/(0.1 + x*x + y*y + z*z);
        u[VAR::U_SYMAT2] = 1.3 - x*y/(3.0 + x*x + 2.0*y*y + z*z)*(x*x+z*z);
        u[VAR::U_SYMAT3] = exp(-4.0*cos(x)*sin(y))*(cos(y)-0.33333333330*exp(4*cos(x)*sin(y))*(1+0.2*cos(y))*(5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
        u[VAR::U_SYMAT4] = -1.0 + y*z/(1.0 + 3.0*x*x + y*y + z*z);
        u[VAR::U_SYMAT5] = exp(-4.0*cos(x)*sin(y))*(cos(z)-0.3333333333*exp(4*cos(x)*sin(y))/(1+0.2*sin(x))/(1+0.2*cos(y))*(5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
        */

        const double x=GRIDX_TO_X(xx1);
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
        u[VAR::U_SYMGT1] = 0.07*(2.0 + cos(x*x + y*y));
        u[VAR::U_SYMGT2] = 0.1*(3.0 + sin(z)*cos(x));
        u[VAR::U_SYMGT4] = 0.15*(1.751 - sin(x*x)*cos(y)*cos(z));

        u[VAR::U_K] = 5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)
                      +5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)
                      +0.4*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))
                       *exp(-4.0*cos(x)*sin(y))*cos(z);
        u[VAR::U_K] *= 0.01234;

        u[VAR::U_SYMAT0] = exp(-4.0*cos(x)*sin(y))*(cos(x) -0.3333333333*exp(4.0*cos(x)*sin(y)) *(1.0+0.2*sin(x))*(5.0*exp(-4.0*cos(x)*sin(y)) /(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y)) /(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y) +5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
        u[VAR::U_SYMAT1] = 1.0 + x*z/(0.1 + x*x + y*y + z*z);
        u[VAR::U_SYMAT2] = 1.3 - x*y/(3.0 + x*x + 2.0*y*y + z*z)*(x*x+z*z);
        u[VAR::U_SYMAT3] = exp(-4.0*cos(x)*sin(y))*(cos(y)-0.33333333330*exp(4*cos(x)*sin(y))*(1+0.2*cos(y))*(5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));
        u[VAR::U_SYMAT4] = -1.0 + y*z/(1.0 + 3.0*x*x + y*y + z*z);
        u[VAR::U_SYMAT5] = exp(-4.0*cos(x)*sin(y))*(cos(z)-0.3333333333*exp(4*cos(x)*sin(y))/(1+0.2*sin(x))/(1+0.2*cos(y))*(5.0*exp(-4.0*cos(x)*sin(y))/(5.0+sin(x))*cos(x)+5.0*exp(-4.0*cos(x)*sin(y))/(5.0+cos(y))*cos(y)+0.04*(25.0+5.0*cos(y)+5.0*sin(x)+sin(x)*cos(y))*exp(-4.0*cos(x)*sin(y))*cos(z)));


        /* Enforce MASSGRAV constraints */
        double gtd[3][3], Atd[3][3];

        gtd[0][0] = u[VAR::U_SYMGT0];
        gtd[0][1] = u[VAR::U_SYMGT1];
        gtd[0][2] = u[VAR::U_SYMGT2];
        gtd[1][0] = gtd[0][1];
        gtd[1][1] = u[VAR::U_SYMGT3];
        gtd[1][2] = u[VAR::U_SYMGT4];
        gtd[2][0] = gtd[0][2];
        gtd[2][1] = gtd[1][2];
        gtd[2][2] = u[VAR::U_SYMGT5];

        Atd[0][0] = u[VAR::U_SYMAT0];
        Atd[0][1] = u[VAR::U_SYMAT1];
        Atd[0][2] = u[VAR::U_SYMAT2];
        Atd[1][0] = Atd[0][1];
        Atd[1][1] = u[VAR::U_SYMAT3];
        Atd[1][2] = u[VAR::U_SYMAT4];
        Atd[2][0] = Atd[0][2];
        Atd[2][1] = Atd[1][2];
        Atd[2][2] = u[VAR::U_SYMAT5];

        const double one_third = 1.0 / 3.0;
        double det_gtd  =  gtd[0][0]*( gtd[1][1]*gtd[2][2] - gtd[1][2]*gtd[1][2])
                           -  gtd[0][1]*gtd[0][1]*gtd[2][2]
                           +  2.0*gtd[0][1]*gtd[0][2]*gtd[1][2]
                           -  gtd[0][2]*gtd[0][2]*gtd[1][1];

        if (det_gtd < 0.0) {
            /* FIXME What to do here? The metric is not physical. Do we reset the metric to be flat? */
            gtd[0][0] = 1.0; gtd[0][1] = 0.0; gtd[0][2] = 0.0;
            gtd[1][0] = 0.0; gtd[1][1] = 1.0; gtd[1][2] = 0.0;
            gtd[2][0] = 0.0; gtd[2][1] = 0.0; gtd[2][2] = 1.0;
            det_gtd = 1.0;
        }
        double det_gtd_to_neg_third = 1.0 / pow(det_gtd, one_third);

        for (unsigned int j = 0; j < 3; j++)
        {
            for (unsigned int i = 0; i < 3; i++)
            {
                gtd[i][j] *= det_gtd_to_neg_third;
            }
        }

        det_gtd =   gtd[0][0]*( gtd[1][1]*gtd[2][2] - gtd[1][2]*gtd[1][2])
                    - gtd[0][1]*gtd[0][1]*gtd[2][2]
                    + 2.0*gtd[0][1]*gtd[0][2]*gtd[1][2]
                    - gtd[0][2]*gtd[0][2]*gtd[1][1];

        double detgt_m1 = det_gtd - 1.0;

        if (fabs(detgt_m1) > 1.0e-6) {
            std::cout.precision(14);
            std::cout<<"enforce_massgrav_constraint: det(gtd) != 1. det="<<std::fixed<<det_gtd<<std::endl;
            std::cout<<"      gtd(1,1)="<<gtd[0][0]<<std::endl;
            std::cout<<"      gtd(1,2)="<<gtd[0][1]<<std::endl;
            std::cout<<"      gtd(1,3)="<<gtd[0][2]<<std::endl;
            std::cout<<"      gtd(2,2)="<<gtd[1][1]<<std::endl;
            std::cout<<"      gtd(2,3)="<<gtd[1][2]<<std::endl;
            std::cout<<"      gtd(3,3)="<<gtd[2][2]<<std::endl;
        }

        double gtu[3][3];
        double idet_gtd = 1.0/det_gtd;
        gtu[0][0] = idet_gtd*(gtd[1][1]*gtd[2][2]-gtd[1][2]*gtd[1][2]);
        gtu[0][1] = idet_gtd*(-gtd[0][1]*gtd[2][2]+gtd[0][2]*gtd[1][2]);
        gtu[0][2] = idet_gtd*(gtd[0][1]*gtd[1][2]-gtd[0][2]*gtd[1][1]);
        gtu[1][0] = gtu[0][1];
        gtu[1][1] = idet_gtd*(gtd[0][0]*gtd[2][2]-gtd[0][2]*gtd[0][2]);
        gtu[1][2] = idet_gtd*(-gtd[0][0]*gtd[1][2]+gtd[0][1]*gtd[0][2]);
        gtu[2][0] = gtu[0][2];
        gtu[2][1] = gtu[1][2];
        gtu[2][2] = idet_gtd*(gtd[0][0]*gtd[1][1]-gtd[0][1]*gtd[0][1]);

        /* Require Atd to be traceless. */
        double one_third_trace_Atd =   one_third * (
                Atd[0][0]*gtu[0][0]
                + Atd[1][1]*gtu[1][1]
                + Atd[2][2]*gtu[2][2]
                + 2.0 * (   Atd[0][1]*gtu[0][1]
                            + Atd[0][2]*gtu[0][2]
                            + Atd[1][2]*gtu[1][2]  )
        );

        Atd[0][0] -= one_third_trace_Atd * gtd[0][0];
        Atd[0][1] -= one_third_trace_Atd * gtd[0][1];
        Atd[0][2] -= one_third_trace_Atd * gtd[0][2];
        Atd[1][1] -= one_third_trace_Atd * gtd[1][1];
        Atd[1][2] -= one_third_trace_Atd * gtd[1][2];
        Atd[2][2] -= one_third_trace_Atd * gtd[2][2];

        double tr_A =    Atd[0][0]*gtu[0][0]
                         + Atd[1][1]*gtu[1][1]
                         + Atd[2][2]*gtu[2][2]
                         + 2.0 * (   Atd[0][1]*gtu[0][1]
                                     + Atd[0][2]*gtu[0][2]
                                     + Atd[1][2]*gtu[1][2]  );

        if (fabs(tr_A) > 1.0e-6) {
            std::cout<<"enforce_massgrav_constraint: tr_A != 0. tr_A="<<tr_A<<std::endl;
            std::cout<<"      Atd(1,1)="<<Atd[0][0]<<std::endl;
            std::cout<<"      Atd(1,2)="<<Atd[0][1]<<std::endl;
            std::cout<<"      Atd(1,3)="<<Atd[0][2]<<std::endl;
            std::cout<<"      Atd(2,2)="<<Atd[1][1]<<std::endl;
            std::cout<<"      Atd(2,3)="<<Atd[1][2]<<std::endl;
            std::cout<<"      Atd(3,3)="<<Atd[2][2]<<std::endl;
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


    void blockAdaptiveOctree(std::vector<ot::TreeNode>& tmpNodes,const Point& pt_min,const Point & pt_max,const unsigned int regLev,const unsigned int maxDepth,MPI_Comm comm)
    {
        int rank,npes;
        MPI_Comm_size(comm,&npes);
        MPI_Comm_rank(comm,&rank);

        double pt_g_min[3];
        double pt_g_max[3];

        pt_g_min[0]=X_TO_GRIDX(pt_min.x());
        pt_g_min[1]=Y_TO_GRIDY(pt_min.y());
        pt_g_min[2]=Z_TO_GRIDZ(pt_min.z());

        pt_g_max[0]=X_TO_GRIDX(pt_max.x());
        pt_g_max[1]=Y_TO_GRIDY(pt_max.y());
        pt_g_max[2]=Z_TO_GRIDZ(pt_max.z());

        assert(pt_g_min[0]>=0 && pt_g_min[0]<=(1u<<maxDepth));
        assert(pt_g_min[1]>=0 && pt_g_min[1]<=(1u<<maxDepth));
        assert(pt_g_min[2]>=0 && pt_g_min[2]<=(1u<<maxDepth));

        assert(pt_g_max[0]>=0 && pt_g_max[0]<=(1u<<maxDepth));
        assert(pt_g_max[1]>=0 && pt_g_max[1]<=(1u<<maxDepth));
        assert(pt_g_max[2]>=0 && pt_g_max[2]<=(1u<<maxDepth));


        unsigned int xRange_b,xRange_e;
        unsigned int yRange_b=pt_g_min[1],yRange_e=pt_g_max[1];
        unsigned int zRange_b=pt_g_min[2],zRange_e=pt_g_max[2];

        xRange_b=pt_g_min[0];//(rank*(pt_g_max[0]-pt_g_min[0]))/npes + pt_g_min[0];
        xRange_e=pt_g_max[1];//((rank+1)*(pt_g_max[0]-pt_g_min[0]))/npes + pt_g_min[0];

        unsigned int stepSz=1u<<(maxDepth-regLev);

       /* std::cout<<" x min: "<<xRange_b<<" x_max: "<<xRange_e<<std::endl;
        std::cout<<" y min: "<<yRange_b<<" y_max: "<<yRange_e<<std::endl;
        std::cout<<" z min: "<<zRange_b<<" z_max: "<<zRange_e<<std::endl;*/


        for(unsigned int x=xRange_b;x<xRange_e;x+=stepSz)
            for(unsigned int y=yRange_b;y<yRange_e;y+=stepSz)
              for(unsigned int z=zRange_b;z<zRange_e;z+=stepSz)
	      {
		 if(x>=(1u<<maxDepth)) x=x-1;	
		 if(y>=(1u<<maxDepth)) y=y-1;	
		 if(z>=(1u<<maxDepth)) z=z-1;	
                 
		 tmpNodes.push_back(ot::TreeNode(x,y,z,regLev,m_uiDim,maxDepth));
	      }


        return ;


    }

    double computeWTol(double x,double y,double z,double tolMin)
    {
       double origin[3];
       origin[0]=(double)(1u<<massgrav::MASSGRAV_MAXDEPTH-1);
       origin[1]=(double)(1u<<massgrav::MASSGRAV_MAXDEPTH-1);
       origin[2]=(double)(1u<<massgrav::MASSGRAV_MAXDEPTH-1);

       double r = sqrt(  GRIDX_TO_X(x)*GRIDX_TO_X(x)
                        + GRIDY_TO_Y(y)*GRIDY_TO_Y(y)
                        + GRIDZ_TO_Z(z)*GRIDZ_TO_Z(z)
                        );

       const double tolMax = massgrav::MASSGRAV_WAVELET_TOL_MAX;
       const double R0 = massgrav::MASSGRAV_WAVELET_TOL_FUNCTION_R0;
       const double R1 = massgrav::MASSGRAV_WAVELET_TOL_FUNCTION_R1;

       if (massgrav::MASSGRAV_USE_WAVELET_TOL_FUNCTION == 1) {
         return std::min( tolMax, std::max( tolMin, (tolMax - tolMin)/(R1 - R0)*(r - R0) + tolMin));
       }
       else {
         return tolMin;
       }
    }



    void writeBLockToBinary(const double **unzipVarsRHS,unsigned int offset,const double *pmin, const double *pmax,double* bxMin,double * bxMax, const unsigned int *sz,unsigned int blkSz,double dxFactor,const char* fprefix)
    {


        const unsigned int nx = sz[0];
        const unsigned int ny = sz[1];
        const unsigned int nz = sz[2];

        const unsigned int ib=3;
        const unsigned int jb=3;
        const unsigned int kb=3;

        const unsigned int ie=nx-3;
        const unsigned int je=ny-3;
        const unsigned int ke=nz-3;

        const unsigned int blkInlSz=(nx-3)*(ny-3)*(nz-3);


        double hx = (pmax[0] - pmin[0]) / (nx - 1);
        double hy = (pmax[1] - pmin[1]) / (ny - 1);
        double hz = (pmax[2] - pmin[2]) / (nz - 1);

        const double dx=(massgrav::MASSGRAV_COMPD_MAX[0]-massgrav::MASSGRAV_COMPD_MIN[0])*(1.0/(double)(1u<<massgrav::MASSGRAV_MAXDEPTH));
        unsigned int level= massgrav::MASSGRAV_MAXDEPTH-((unsigned int)(hx/dx)-1);

        MPI_Comm comm=MPI_COMM_WORLD;

        int rank,npes;
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);
        //std::cout<<"ranl: "<<rank<<"npes: "<<npes<<std::endl;

        //std::cout<<"nx: "<<nx<<" level: "<<level<<" hx: "<<hx<<" dx: "<<dx<<std::endl;

        if((hx>(dxFactor*dx)) || (pmin[0]<bxMin[0] || pmin[1]<bxMin[1] || pmin[2]<bxMin[2] ) || (pmax[0]>bxMax[0] || pmax[1]>bxMax[1] || pmax[2]>bxMax[2]))
            return;

        double * blkInternal=new double[blkInlSz];
        for(unsigned int var=0;var<massgrav::MASSGRAV_NUM_VARS;var++)
        {
            char fName[256];
            sprintf(fName,"%s_%s_n_%d_r_%d_p_%d.bin",fprefix,massgrav::MASSGRAV_VAR_NAMES[var],nx,rank,npes);
            FILE* outfile = fopen(fName,"w");
            if(outfile==NULL) {std::cout<<fName<<" file open failed "<<std::endl;}

            for(unsigned int k=kb;k<ke;k++)
              for(unsigned int j=jb;j<je;j++)
                for(unsigned int i=ib;i<ie;i++)
                    blkInternal[k*(ny-3)*(nx-3)+j*(nx-3)+i]=unzipVarsRHS[var][offset+k*(ny*nx)+j*(ny)+i];


            fwrite(blkInternal,sizeof(double),blkInlSz,outfile); // write out the number of elements.
            fclose(outfile);

        }

        delete [] blkInternal;

    }


    unsigned int getOctantWeight(const ot::TreeNode* pNode)
    {
        return pNode->getLevel();
    }

}// end of namespace massgrav





namespace massgrav
{

    namespace timer
    {
        void initFlops()
        {
            total_runtime.start();
            t_f2o.start();
            t_cons.start();
            t_bal.start();
            t_mesh.start();
            t_rkSolve.start();
            t_ghostEx_sync.start();
            t_unzip_sync.start();

            for(unsigned int i=0;i<NUM_FACES;i++)
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

        void resetSnapshot()
        {

            total_runtime.snapreset();
            t_f2o.snapreset();
            t_cons.snapreset();
            t_bal.snapreset();
            t_mesh.snapreset();
            t_rkSolve.snapreset();
            t_ghostEx_sync.snapreset();
            t_unzip_sync.snapreset();

            for(unsigned int i=0;i<NUM_FACES;i++)
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

            t_deriv.snapreset();
            t_rhs.snapreset();

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
            t_isReMesh.snapreset();
            t_gridTransfer.snapreset();
            t_ioVtu.snapreset();
            t_ioCheckPoint.snapreset();

        }

        void profileInfo(const char* filePrefix,const ot::Mesh* pMesh)
        {


            int activeRank,activeNpes,globalRank,globalNpes;

            MPI_Comm commActive;
            MPI_Comm commGlobal;

            if(pMesh->isActive())
            {
                commActive=pMesh->getMPICommunicator();
                activeRank=pMesh->getMPIRank();
                activeNpes=pMesh->getMPICommSize();
            }

            globalRank=pMesh->getMPIRankGlobal();
            globalNpes=pMesh->getMPICommSizeGlobal();
            commGlobal=pMesh->getMPIGlobalCommunicator();

            double t_stat;
            double t_stat_g[3];

            const char separator    = ' ';
            const int nameWidth     = 30;
            const int numWidth      = 10;

            char fName[256];
            std::ofstream outfile;

            DendroIntL localSz,globalSz;

            if(!activeRank)
            {
                sprintf(fName,"%s_final.prof",filePrefix);
                outfile.open (fName);
                if(outfile.fail()) {std::cout<<fName<<" file open failed "<<std::endl; return ;}

                outfile<<"active npes : "<<activeNpes<<std::endl;
                outfile<<"global npes : "<<globalNpes<<std::endl;
                outfile<<"partition tol : "<<massgrav::MASSGRAV_LOAD_IMB_TOL<<std::endl;
                outfile<<"wavelet tol : "<<massgrav::MASSGRAV_WAVELET_TOL<<std::endl;
                outfile<<"maxdepth : "<<massgrav::MASSGRAV_MAXDEPTH<<std::endl;

            }

            MPI_Comm comm=commActive;
            unsigned int rank =activeRank;

            localSz=pMesh->getNumLocalMeshElements();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<"Elements : "<<globalSz<<std::endl;

            localSz=pMesh->getNumLocalMeshNodes();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<"DOG(zip) : "<<globalSz<<std::endl;

            localSz=pMesh->getDegOfFreedomUnZip();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<"DOG(unzip) : "<<globalSz<<std::endl;


            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "step";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "min(s)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "mean(s)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "max(s)"<<std::endl;


            t_stat=total_runtime.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"+runtime(s)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_f2o.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<" ++f2o";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_cons.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<" ++construction";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_rkSolve.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<" ++rkSolve";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_bal.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --2:1 balance";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_mesh.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --mesh";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_rkStep.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --rkstep";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_ghostEx_sync.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --ghostExchge.";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_unzip_sync.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_unzip_async.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++unzip_async";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

#ifdef ENABLE_DENDRO_PROFILE_COUNTERS
            t_stat=dendro::timer::t_unzip_async_internal.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_internal";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_async_external.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_external";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_async_comm.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_comm (comm) ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;
#endif

            t_stat=t_deriv.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --deriv ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_bdyc.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --boundary con ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;



            t_stat=t_zip.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --zip";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_ioVtu.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --vtu";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_ioCheckPoint.seconds;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --checkpoint";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            if(!rank) outfile.close();




        }


        void profileInfoIntermediate(const char* filePrefix,const ot::Mesh* pMesh,const unsigned int currentStep)
        {

            int activeRank,activeNpes,globalRank,globalNpes;

            MPI_Comm commActive;
            MPI_Comm commGlobal;

            if(pMesh->isActive())
            {
                commActive=pMesh->getMPICommunicator();
                activeRank=pMesh->getMPIRank();
                activeNpes=pMesh->getMPICommSize();
            }

            globalRank=pMesh->getMPIRankGlobal();
            globalNpes=pMesh->getMPICommSizeGlobal();
            commGlobal=pMesh->getMPIGlobalCommunicator();

            double t_stat;
            double t_stat_g[3];

            const char separator    = ' ';
            const int nameWidth     = 30;
            const int numWidth      = 10;

            char fName[256];
            std::ofstream outfile;

            DendroIntL localSz,globalSz;

            DendroIntL ghostElements;
            DendroIntL localElements;

            DendroIntL ghostNodes;
            DendroIntL localNodes;

            DendroIntL totalSendNode;
            DendroIntL totalRecvNode;

            DendroIntL numCalls;


#ifdef MASSGRAV_PROFILE_HUMAN_READABLE
            if(!activeRank)
            {
                sprintf(fName,"%s_im.prof",filePrefix);
                outfile.open (fName,std::fstream::app);
                if(outfile.fail()) {std::cout<<fName<<" file open failed "<<std::endl; return ;}

                outfile<<"active npes : "<<activeNpes<<std::endl;
                outfile<<"global npes : "<<globalNpes<<std::endl;
                outfile<<"current step : "<<currentStep<<std::endl;
                outfile<<"partition tol : "<<massgrav::MASSGRAV_LOAD_IMB_TOL<<std::endl;
                outfile<<"wavelet tol : "<<massgrav::MASSGRAV_WAVELET_TOL<<std::endl;
                outfile<<"maxdepth : "<<massgrav::MASSGRAV_MAXDEPTH<<std::endl;

            }

            MPI_Comm comm=commActive;
            unsigned int rank =activeRank;

            localSz=pMesh->getNumLocalMeshElements();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<"Elements : "<<globalSz<<std::endl;

            localSz=pMesh->getNumLocalMeshNodes();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<"DOG(zip) : "<<globalSz<<std::endl;

            localSz=pMesh->getDegOfFreedomUnZip();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<"DOG(unzip) : "<<globalSz<<std::endl;


            ghostElements=pMesh->getNumPreGhostElements()+pMesh->getNumPostGhostElements();
            localElements=pMesh->getNumLocalMeshElements();

            ghostNodes=pMesh->getNumPreMeshNodes()+pMesh->getNumPostMeshNodes();
            localNodes=pMesh->getNumLocalMeshNodes();

            if(!rank)outfile<<"========================= MESH ======================================================================= "<<std::endl;

            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "step";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "min(#)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "mean(#)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "max(#)"<<std::endl;

            t_stat=ghostElements;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"ghost Elements";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=localElements;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"local Elements";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=ghostNodes;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"ghost Nodes";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=localNodes;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"local Nodes";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=pMesh->getGhostExcgTotalSendNodeCount();
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"send Nodes";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=pMesh->getGhostExcgTotalRecvNodeCount();
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"recv Nodes";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            if(!rank)outfile<<"========================= RUNTIME =================================================================== "<<std::endl;
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "step";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "min(s)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "mean(s)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) << "max(s)"<<std::endl;




            /* t_stat=total_runtime.seconds;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"+runtime(s)";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_f2o.seconds;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<" ++f2o";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_cons.seconds;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<" ++construction";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_rkSolve.seconds;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<" ++rkSolve";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;*/


            t_stat=t_bal.snap;
            //numCalls=t_bal.num_calls;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++2:1 balance";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_mesh.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++mesh";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_rkStep.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++rkstep";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_ghostEx_sync.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++ghostExchge.";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_unzip_sync.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++unzip_sync";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_unzip_async.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++unzip_async";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;



            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS

            t_stat=dendro::timer::t_unzip_async_comm.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_comm_wait (comm) ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_nodalval.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_nodalVal";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_f_c1.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --t_unzip_sync_f_c1";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_f_c2.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --t_unzip_sync_f_c2";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_f_c3.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --t_unzip_sync_f_c3";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_cpy.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --t_unzip_sync_cpy";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=dendro::timer::t_unzip_sync_internal.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_internal";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_face[0].snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_face_left";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_face[1].snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_face_right";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_face[2].snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_face_down";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_face[3].snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_face_up";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_face[4].snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_face_back";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_face[5].snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_face_front";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=dendro::timer::t_unzip_sync_edge.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_edge";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_sync_vtex.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_sync_vtex";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=dendro::timer::t_unzip_p2c.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_p2c";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;
            #endif

            /*
            #ifdef ENABLE_DENDRO_PROFILE_COUNTERS
            t_stat=t_unzip_async_internal.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_internal";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_unzip_async_external.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_external";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_unzip_async_comm.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --unzip_comm (comm) ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;
            #endif
            */
            t_stat=t_isReMesh.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++isReMesh";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_gridTransfer.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++gridTransfer";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;



            t_stat=t_deriv.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++deriv ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++compute_rhs ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;



            t_stat=t_rhs_a.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_a ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs_b.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_b ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs_gt.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_gt ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_rhs_chi.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_chi ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs_At.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_At ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs_K.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_K ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs_Gt.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_Gt ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_rhs_B.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  --compute_rhs_B ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;



            t_stat=t_bdyc.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++boundary con ";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;



            t_stat=t_zip.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++zip";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            t_stat=t_ioVtu.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++vtu";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;

            t_stat=t_ioCheckPoint.snap;
            computeOverallStats(&t_stat, t_stat_g, comm);
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator) <<"  ++checkpoint";
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[0];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[1];
            if(!rank)outfile << std::left << std::setw(nameWidth) << std::setfill(separator)<<t_stat_g[2]<<std::endl;


            if(!rank) outfile.close();
#else

            if(!activeRank)
            {
                sprintf(fName,"%s_im.prof",filePrefix);
                outfile.open (fName,std::fstream::app);
                if(outfile.fail()) {std::cout<<fName<<" file open failed "<<std::endl; return ;}

                //writes the header
                if(currentStep==0)
                 outfile<<"step\t act_npes\t glb_npes\t part_tol\t wave_tol\t maxdepth\t numOcts\t dof_zip\t dof_unzip\t"<<\
                 "element_ghost_min\t element_ghost_mean\t element_ghost_max\t"<<\
                 "element_local_min\t element_local_mean\t element_local_max\t"<<\
                 "nodes_local_min\t nodes_local_mean\t nodes_local|max\t"<<\
                 "send_nodes_min\t send_nodes_mean\t send_nodes_max\t"<<\
                 "recv_nodes_min\t recv_nodes_mean\t recv_nodes_max\t"<<\
                 "bal_min\t bal_mean\t bal_max\t"<<\
                 "mesh_min\t mesh_mean\t mesh_max\t"<<\
                 "rkstep_min\t rkstep_mean\t rkstep_max\t"<<\
                 "ghostEx_min\t ghostEx_mean\t ghostEx_max\t"<<\
                 "unzip_sync_min\t unzip_sync_mean\t unzip_sync_max\t"<<\
                 "unzip_async_min\t unzip_async_mean\t unzip_async_max\t"<<\
                 "unzip_async_wait_min\t unzip_async_wait_mean\t unzip_async_wait_max\t"<<\
                 "isRemesh_min\t isRemesh_mean\t isRemesh_max\t"<<\
                 "GT_min\t GT_mean\t GT_max\t"<<\
                 "deriv_min\t deriv_mean\t deriv_max\t"<<\
                 "rhs_min\t rhs_mean\t rhs_max\t"<<std::endl;

            }

            MPI_Comm comm=commActive;
            unsigned int rank =activeRank;

            if(!rank) outfile<<currentStep<<"\t ";
            if(!rank) outfile<<activeNpes<<"\t ";
            if(!rank) outfile<<globalNpes<<"\t ";
            if(!rank) outfile<<massgrav::MASSGRAV_LOAD_IMB_TOL<<"\t ";
            if(!rank) outfile<<massgrav::MASSGRAV_WAVELET_TOL<<"\t ";
            if(!rank) outfile<<massgrav::MASSGRAV_MAXDEPTH<<"\t ";

            localSz=pMesh->getNumLocalMeshElements();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<globalSz<<"\t ";

            localSz=pMesh->getNumLocalMeshNodes();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<globalSz<<"\t ";

            localSz=pMesh->getDegOfFreedomUnZip();
            par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,comm);
            if(!rank)outfile<<globalSz<<"\t ";

            ghostElements=pMesh->getNumPreGhostElements()+pMesh->getNumPostGhostElements();
            localElements=pMesh->getNumLocalMeshElements();

            t_stat=ghostElements;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=localElements;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            ghostNodes=pMesh->getNumPreMeshNodes()+pMesh->getNumPostMeshNodes();
            localNodes=pMesh->getNumLocalMeshNodes();

            /*t_stat=ghostNodes;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";*/

            t_stat=localNodes;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=pMesh->getGhostExcgTotalSendNodeCount();
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=pMesh->getGhostExcgTotalRecvNodeCount();
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_bal.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_mesh.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_rkStep.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_ghostEx_sync.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_unzip_sync.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_unzip_async.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=dendro::timer::t_unzip_async_comm.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_isReMesh.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_gridTransfer.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_deriv.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            t_stat=t_rhs.snap;
            computeOverallStats(&t_stat,t_stat_g,comm);
            if(!rank) outfile<<t_stat_g[0]<<"\t "<<t_stat_g[1]<<"\t "<<t_stat_g[2]<<"\t ";

            if(!rank) outfile<<std::endl;
            if(!rank) outfile.close();
#endif




        }


    }

}



namespace GW
{
    void psi4ShpereDump(const ot::Mesh* mesh, DendroScalar ** cVar,unsigned int timestep,double time, double dtheta, double dphi)
    {
        
        unsigned int rankGlobal=mesh->getMPIRankGlobal();
        unsigned int npesGlobal=mesh->getMPICommSizeGlobal();
        MPI_Comm commGlobal=mesh->getMPIGlobalCommunicator();


        const unsigned int nTheta = (M_PI)/dtheta;
        const unsigned int nPhi = (2*M_PI)/dphi;
        const unsigned int numPts = nTheta * nPhi ;

        
        unsigned int totalModes=0;
        for(unsigned int l=0;l<MASSGRAV_GW_NUM_LMODES;l++)
            totalModes+=2*MASSGRAV_GW_L_MODES[l]+1;

        const unsigned int TOTAL_MODES=totalModes;

        DendroComplex * swsh_coeff = new DendroComplex[MASSGRAV_GW_NUM_RADAII*TOTAL_MODES];
        DendroComplex * swsh_coeff_g = new DendroComplex[MASSGRAV_GW_NUM_RADAII*TOTAL_MODES];

        std::vector<unsigned int> lmCounts;
        std::vector<unsigned int> lmOffset;

        lmCounts.resize(MASSGRAV_GW_NUM_LMODES);
        lmOffset.resize(MASSGRAV_GW_NUM_LMODES);

        for(unsigned int l=0;l<MASSGRAV_GW_NUM_LMODES;l++)
            lmCounts[l]=2*MASSGRAV_GW_L_MODES[l]+1;


        lmOffset[0]=0;
        omp_par::scan(&(*(lmCounts.begin())),&(*(lmOffset.begin())),MASSGRAV_GW_NUM_LMODES);

        if(mesh->isActive())
        {

            const unsigned int rankActive=mesh->getMPIRank();
            const unsigned int npesActive=mesh->getMPICommSize();

            std::vector<double> coords;
            coords.reserve(3*numPts);

            std::vector<double> psi4_real;
            psi4_real.resize(numPts);

            std::vector<double> psi4_imag;
            psi4_imag.resize(numPts);

            Point grid_limits[2];
            Point domain_limits[2];

            grid_limits[0] = Point(massgrav::MASSGRAV_OCTREE_MIN[0], massgrav::MASSGRAV_OCTREE_MIN[1], massgrav::MASSGRAV_OCTREE_MIN[2]);
            grid_limits[1] = Point(massgrav::MASSGRAV_OCTREE_MAX[0], massgrav::MASSGRAV_OCTREE_MAX[1], massgrav::MASSGRAV_OCTREE_MAX[2]);

            domain_limits[0] = Point(massgrav::MASSGRAV_COMPD_MIN[0], massgrav::MASSGRAV_COMPD_MIN[1], massgrav::MASSGRAV_COMPD_MIN[2]);
            domain_limits[1] = Point(massgrav::MASSGRAV_COMPD_MAX[0], massgrav::MASSGRAV_COMPD_MAX[1], massgrav::MASSGRAV_COMPD_MAX[2]);


            std::vector<unsigned int > validIndex;

            for(unsigned int k=0;k<MASSGRAV_GW_NUM_RADAII;k++)
            {
                for (unsigned int i=0; i< nTheta ; i++ )
                    for(unsigned int j=0; j< nPhi ; j++)
                    {
                        double x = MASSGRAV_GW_RADAII[k]*sin(j*dtheta) * cos(i*dphi) ;
                        double y = MASSGRAV_GW_RADAII[k]*sin(j*dtheta) * sin(i*dphi) ;
                        double z = MASSGRAV_GW_RADAII[k]*cos(j*dtheta);

                        coords.push_back(x);
                        coords.push_back(y);
                        coords.push_back(z);

                    }


                validIndex.clear();
                ot::da::interpolateToCoords(mesh,cVar[massgrav::VAR_CONSTRAINT::C_PSI4_REAL],&(*(coords.begin())),coords.size(), grid_limits, domain_limits , &(*(psi4_real.begin())),validIndex);

                validIndex.clear();
                ot::da::interpolateToCoords(mesh,cVar[massgrav::VAR_CONSTRAINT::C_PSI4_IMG],&(*(coords.begin())),coords.size(),  grid_limits, domain_limits ,&(*(psi4_imag.begin())),validIndex);

            }


        }





    }


}
