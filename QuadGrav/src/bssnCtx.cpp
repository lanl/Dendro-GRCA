/**
 * @file quadgravCtx.cpp
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief BSSN contex file. 
 * @version 0.1
 * @date 2019-12-20
 * 
 * School of Computing, University of Utah. 
 * @copyright Copyright (c) 2019
 * 
 */

#include "quadgravCtx.h"

namespace quadgrav
{
    BSSNCtx::BSSNCtx(ot::Mesh* pMesh) : Ctx()
    {   
        m_uiMesh = pMesh;
        // variable allocation for evolution variables
        m_uiEVar = this->create_vec(ts::CTXVType::EVOLUTION,true,false,false,BSSN_NUM_VARS);

        // variable allocation for constraint variables. 
        m_uiCVar = this->create_vec(ts::CTXVType::CONSTRAINT,true,false,false,BSSN_CONSTRAINT_NUM_VARS);

        // evars unzip in (0) evars unzip out (1) 
        m_uiEUnzip[0] = this->create_vec(ts::CTXVType::EVOLUTION,false,true,false,BSSN_NUM_VARS);
        m_uiEUnzip[1] = this->create_vec(ts::CTXVType::EVOLUTION,false,true,false,BSSN_NUM_VARS);

        // // constraint var unzip out (0)
        m_uiCUnzip[0] = this->create_vec(ts::CTXVType::CONSTRAINT,false,true,false,BSSN_CONSTRAINT_NUM_VARS); 

        m_uiTinfo._m_uiStep=0;
        m_uiTinfo._m_uiT = 0;
        m_uiTinfo._m_uiTb = quadgrav::BSSN_RK_TIME_BEGIN;
        m_uiTinfo._m_uiTe = quadgrav::BSSN_RK_TIME_END;
        m_uiTinfo._m_uiTh = quadgrav::BSSN_RK45_TIME_STEP_SIZE;     

        m_uiElementOrder = quadgrav::BSSN_ELE_ORDER;

        m_uiMinPt = Point(quadgrav::BSSN_GRID_MIN_X,quadgrav::BSSN_GRID_MIN_Y,quadgrav::BSSN_GRID_MIN_Z);
        m_uiMaxPt = Point(quadgrav::BSSN_GRID_MAX_X,quadgrav::BSSN_GRID_MAX_Y,quadgrav::BSSN_GRID_MAX_Z);
        
        return;

    }

    BSSNCtx::~BSSNCtx()
    {
        this->destroy_vec(m_uiEVar);
        this->destroy_vec(m_uiCVar);

        this->destroy_vec(m_uiEUnzip[0]);
        this->destroy_vec(m_uiEUnzip[1]);
        this->destroy_vec(m_uiCUnzip[0]);

    }

    int BSSNCtx::initialize()
    {
        if(quadgrav::BSSN_RESTORE_SOLVER)
        {
            this->restore_checkpt();
            return 0;
        }
        
        unsigned int nodeLookUp_CG;
        unsigned int nodeLookUp_DG;
        DendroScalar x,y,z,len;
        const ot::TreeNode * pNodes=&(*(m_uiMesh->getAllElements().begin()));
        unsigned int ownerID,ii_x,jj_y,kk_z;
        unsigned int eleOrder=m_uiMesh->getElementOrder();
        const unsigned int * e2n_cg=&(*(m_uiMesh->getE2NMapping().begin()));
        const unsigned int * e2n_dg=&(*(m_uiMesh->getE2NMapping_DG().begin()));
        const unsigned int nPe=m_uiMesh->getNumNodesPerElement();
        const unsigned int nodeLocalBegin=m_uiMesh->getNodeLocalBegin();
        const unsigned int nodeLocalEnd=m_uiMesh->getNodeLocalEnd();


        DendroScalar* var=new double[quadgrav::BSSN_NUM_VARS];
        DendroScalar** zipIn = new DendroScalar*[quadgrav::BSSN_NUM_VARS];
        m_uiEVar.Get2DArray(zipIn,true);

        DendroScalar mp, mm, mp_adm, mm_adm, E, J1, J2, J3;
        // set the TP communicator. 
        MPI_TP_COMM = m_uiMesh->getMPIGlobalCommunicator();

        for(unsigned int elem=m_uiMesh->getElementLocalBegin(); elem<m_uiMesh->getElementLocalEnd(); elem++)
        {
            for(unsigned int k=0; k<(eleOrder+1); k++)
                for(unsigned int j=0; j<(eleOrder+1); j++ )
                    for(unsigned int i=0; i<(eleOrder+1); i++)
                    {
                        nodeLookUp_CG=e2n_cg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                        if(nodeLookUp_CG>=nodeLocalBegin && nodeLookUp_CG<nodeLocalEnd)
                        {
                            nodeLookUp_DG=e2n_dg[elem*nPe+k*(eleOrder+1)*(eleOrder+1)+j*(eleOrder+1)+i];
                            m_uiMesh->dg2eijk(nodeLookUp_DG,ownerID,ii_x,jj_y,kk_z);
                            len=(double)(1u<<(m_uiMaxDepth-pNodes[ownerID].getLevel()));
                            x=pNodes[ownerID].getX()+ ii_x*(len/(eleOrder));
                            y=pNodes[ownerID].getY()+ jj_y*(len/(eleOrder));
                            z=pNodes[ownerID].getZ()+ kk_z*(len/(eleOrder));
                            
                            if (quadgrav::BSSN_ID_TYPE == 0) {
                                TwoPunctures((double)x,(double)y,(double)z,var,
                                            &mp, &mm, &mp_adm, &mm_adm, &E, &J1, &J2, &J3);
                            }
                            else if (quadgrav::BSSN_ID_TYPE == 1) {
                                quadgrav::punctureData((double)x,(double)y,(double)z,var);
                            }
                            else if (quadgrav::BSSN_ID_TYPE == 2) {
                                quadgrav::KerrSchildData((double)x,(double)y,(double)z,var);
                            }
                            else if (quadgrav::BSSN_ID_TYPE == 3) {
                                quadgrav::noiseData((double)x,(double)y,(double)z,var);
                            }
                            else if (quadgrav::BSSN_ID_TYPE == 4) {
                                quadgrav::fake_initial_data((double)x,(double)y,(double)z,var);
                            }
                            else {
                                std::cout<<"Unknown ID type"<<std::endl;
                            }
                            for(unsigned int v=0; v<quadgrav::BSSN_NUM_VARS; v++)
                                zipIn[v][nodeLookUp_CG]=var[v];
                        }

                    }
        }
        
        for(unsigned int node=m_uiMesh->getNodeLocalBegin(); node<m_uiMesh->getNodeLocalEnd(); node++)
            enforce_quadgrav_constraints(zipIn,node);


        delete [] var;
        delete [] zipIn;

        #ifdef BSSN_EXTRACT_BH_LOCATIONS
            m_uiBHLoc[0]=Point(quadgrav::BH1.getBHCoordX(),quadgrav::BH1.getBHCoordY(),quadgrav::BH1.getBHCoordZ());
            m_uiBHLoc[1]=Point(quadgrav::BH2.getBHCoordX(),quadgrav::BH2.getBHCoordY(),quadgrav::BH2.getBHCoordZ());
        #endif

        return 0;    
    }


    int BSSNCtx::finalize()
    {
        return 0;    
    }

    int BSSNCtx::rhs(DVec* in , DVec* out, unsigned int sz , DendroScalar time)
    {
        // all the variables should be packed together. 
        assert(sz==1);

        DendroScalar ** sVar;
        in[0].Get2DArray(sVar,false);

        for(unsigned int node=m_uiMesh->getNodeLocalBegin(); node< m_uiMesh->getNodeLocalEnd(); node++)
            enforce_quadgrav_constraints(sVar, node);
        
        this->unzip(in[0],m_uiEUnzip[0] , quadgrav::BSSN_ASYNC_COMM_K);
        
        DendroScalar ** unzipIn;
        DendroScalar **  unzipOut; 
        
        m_uiEUnzip[0].Get2DArray(unzipIn , false);
        m_uiEUnzip[1].Get2DArray(unzipOut ,false);

        const ot::Block* blkList = m_uiMesh->getLocalBlockList().data();
        const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();
        
        quadgravRHS(unzipOut,(const DendroScalar**)unzipIn,blkList,numBlocks);

        this->zip(m_uiEUnzip[1], out[0], quadgrav::BSSN_ASYNC_COMM_K);

        delete [] unzipIn;
        delete [] unzipOut;
        delete [] sVar;

        return 0;

    }


    int BSSNCtx::rhs_blkwise(DVec in, DVec out, const unsigned int* const blkIDs, unsigned int numIds, DendroScalar*  blk_time) const
    {

        DendroScalar ** unzipIn;
        DendroScalar **  unzipOut; 

        assert(in.GetDof() == out.GetDof());
        assert(in.IsUnzip() == out.IsUnzip());
        
        in.Get2DArray(unzipIn , false);
        out.Get2DArray(unzipOut ,false);

        const ot::Block* blkList = m_uiMesh->getLocalBlockList().data();
        const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();

        unsigned int offset;
        double ptmin[3], ptmax[3];
        unsigned int sz[3];
        unsigned int bflag;
        double dx,dy,dz;
        const Point pt_min(quadgrav::BSSN_COMPD_MIN[0],quadgrav::BSSN_COMPD_MIN[1],quadgrav::BSSN_COMPD_MIN[2]);
        const Point pt_max(quadgrav::BSSN_COMPD_MAX[0],quadgrav::BSSN_COMPD_MAX[1],quadgrav::BSSN_COMPD_MAX[2]);

        for(unsigned int i=0; i < numIds; i++)
        {
            const unsigned int blk = blkIDs[i];
            assert(blk<numBlocks);

            offset=blkList[blk].getOffset();
            sz[0]=blkList[blk].getAllocationSzX();
            sz[1]=blkList[blk].getAllocationSzY();
            sz[2]=blkList[blk].getAllocationSzZ();

            bflag=blkList[blk].getBlkNodeFlag();

            dx=blkList[blk].computeDx(pt_min,pt_max);
            dy=blkList[blk].computeDy(pt_min,pt_max);
            dz=blkList[blk].computeDz(pt_min,pt_max);

            ptmin[0]=GRIDX_TO_X(blkList[blk].getBlockNode().minX())-3*dx;
            ptmin[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().minY())-3*dy;
            ptmin[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ())-3*dz;

            ptmax[0]=GRIDX_TO_X(blkList[blk].getBlockNode().maxX())+3*dx;
            ptmax[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().maxY())+3*dy;
            ptmax[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ())+3*dz;

            #ifdef BSSN_RHS_STAGED_COMP
                quadgravrhs_sep(unzipOut, (const double **)unzipIn, offset, ptmin, ptmax, sz, bflag);
            #else
                quadgravrhs(unzipOut, (const double **)unzipIn, offset, ptmin, ptmax, sz, bflag);
            #endif

        }

        delete [] unzipIn;
        delete [] unzipOut;

        
    }

    int BSSNCtx::rhs_blk(const DendroScalar* in, DendroScalar* out, unsigned int dof, unsigned int local_blk_id, DendroScalar  blk_time) const 
    {
        //return 0;
        //std::cout<<"quadgrav_rhs"<<std::endl;
        DendroScalar **  unzipIn   = new DendroScalar*[dof];
        DendroScalar **  unzipOut  = new DendroScalar*[dof]; 

        const unsigned int blk = local_blk_id;

        const ot::Block* blkList = m_uiMesh->getLocalBlockList().data();
        const unsigned int numBlocks = m_uiMesh->getLocalBlockList().size();
        double ptmin[3], ptmax[3];
        unsigned int sz[3];
        unsigned int bflag;
        double dx,dy,dz;
        const Point pt_min(quadgrav::BSSN_COMPD_MIN[0],quadgrav::BSSN_COMPD_MIN[1],quadgrav::BSSN_COMPD_MIN[2]);
        const Point pt_max(quadgrav::BSSN_COMPD_MAX[0],quadgrav::BSSN_COMPD_MAX[1],quadgrav::BSSN_COMPD_MAX[2]);

        sz[0]=blkList[blk].getAllocationSzX();
        sz[1]=blkList[blk].getAllocationSzY();
        sz[2]=blkList[blk].getAllocationSzZ();

        const unsigned int NN  = sz[0]*sz[1]*sz[2];

        for(unsigned int v =0; v < dof; v++)
        {
            unzipIn[v] = (DendroScalar*) (in + v*NN);
            unzipOut[v] = (DendroScalar*) (out + v*NN);
        }


        bflag=blkList[blk].getBlkNodeFlag();

        dx=blkList[blk].computeDx(pt_min,pt_max);
        dy=blkList[blk].computeDy(pt_min,pt_max);
        dz=blkList[blk].computeDz(pt_min,pt_max);

        ptmin[0]=GRIDX_TO_X(blkList[blk].getBlockNode().minX())-3*dx;
        ptmin[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().minY())-3*dy;
        ptmin[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ())-3*dz;

        ptmax[0]=GRIDX_TO_X(blkList[blk].getBlockNode().maxX())+3*dx;
        ptmax[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().maxY())+3*dy;
        ptmax[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ())+3*dz;

        #ifdef BSSN_RHS_STAGED_COMP
            quadgravrhs_sep(unzipOut, (const DendroScalar **)unzipIn, 0, ptmin, ptmax, sz, bflag);
        #else
            quadgravrhs(unzipOut, (const DendroScalar **)unzipIn, 0, ptmin, ptmax, sz, bflag);
        #endif

        delete [] unzipIn;
        delete [] unzipOut;


    }


    int BSSNCtx::write_vtu()
    {

        unzip(m_uiEVar,m_uiEUnzip[0],BSSN_ASYNC_COMM_K);
        
        DendroScalar ** evolUnzipVar = NULL;
        DendroScalar ** consUnzipVar = NULL;
        DendroScalar ** consVar = NULL;
        DendroScalar ** evolVar = NULL;

        m_uiEUnzip[0].Get2DArray(evolUnzipVar,false);
        m_uiCUnzip[0].Get2DArray(consUnzipVar,false);
       
        m_uiEVar.Get2DArray(evolVar,false);
        m_uiCVar.Get2DArray(consVar,false);

        #ifdef BSSN_COMPUTE_CONSTRAINTS

            const std::vector<ot::Block> blkList=m_uiMesh->getLocalBlockList();

            unsigned int offset;
            double ptmin[3], ptmax[3];
            unsigned int sz[3];
            unsigned int bflag;
            double dx,dy,dz;
            const Point pt_min(quadgrav::BSSN_COMPD_MIN[0],quadgrav::BSSN_COMPD_MIN[1],quadgrav::BSSN_COMPD_MIN[2]);
            const Point pt_max(quadgrav::BSSN_COMPD_MAX[0],quadgrav::BSSN_COMPD_MAX[1],quadgrav::BSSN_COMPD_MAX[2]);

            for(unsigned int blk=0; blk<blkList.size(); blk++)
            {
                offset=blkList[blk].getOffset();
                sz[0]=blkList[blk].getAllocationSzX();
                sz[1]=blkList[blk].getAllocationSzY();
                sz[2]=blkList[blk].getAllocationSzZ();

                bflag=blkList[blk].getBlkNodeFlag();

                dx=blkList[blk].computeDx(pt_min,pt_max);
                dy=blkList[blk].computeDy(pt_min,pt_max);
                dz=blkList[blk].computeDz(pt_min,pt_max);

                ptmin[0]=GRIDX_TO_X(blkList[blk].getBlockNode().minX())-3*dx;
                ptmin[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().minY())-3*dy;
                ptmin[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().minZ())-3*dz;

                ptmax[0]=GRIDX_TO_X(blkList[blk].getBlockNode().maxX())+3*dx;
                ptmax[1]=GRIDY_TO_Y(blkList[blk].getBlockNode().maxY())+3*dy;
                ptmax[2]=GRIDZ_TO_Z(blkList[blk].getBlockNode().maxZ())+3*dz;

                physical_constraints(consUnzipVar, (const DendroScalar **) evolUnzipVar, offset, ptmin, ptmax, sz, bflag);
            }

            /*double consVecMin[quadgrav::BSSN_CONSTRAINT_NUM_VARS];
            double consVecMax[quadgrav::BSSN_CONSTRAINT_NUM_VARS];*/
            double constraintMaskedL2[quadgrav::BSSN_CONSTRAINT_NUM_VARS];
            
            this->zip(m_uiCUnzip[0],m_uiCVar,quadgrav::BSSN_ASYNC_COMM_K);
            


            quadgrav::extractConstraints(m_uiMesh,(const DendroScalar **)consVar,evolVar[BHLOC::EXTRACTION_VAR_ID],BHLOC::EXTRACTION_TOL,m_uiTinfo._m_uiStep);
            #ifndef BSSN_KERR_SCHILD_TEST
                #ifdef BSSN_EXTRACT_GRAVITATIONAL_WAVES
                GW::extractFarFieldPsi4(m_uiMesh,(const DendroScalar **)consVar,m_uiTinfo._m_uiStep,m_uiTinfo._m_uiT);
                #endif
            #endif

        #endif

        #ifdef BSSN_ENABLE_VTU_OUTPUT
            std::vector<std::string> pDataNames;
            const unsigned int numConstVars = quadgrav::BSSN_NUM_CONST_VARS_VTU_OUTPUT;
            const unsigned int numEvolVars  = quadgrav::BSSN_NUM_EVOL_VARS_VTU_OUTPUT;

            double *pData[(numConstVars+numEvolVars)];

            for(unsigned int i=0; i<numEvolVars; i++)
            {
                pDataNames.push_back(std::string(quadgrav::BSSN_VAR_NAMES[BSSN_VTU_OUTPUT_EVOL_INDICES[i]]));
                pData[i]=evolVar[BSSN_VTU_OUTPUT_EVOL_INDICES[i]];
            }


            for(unsigned int i=0; i<numConstVars; i++)
            {
                pDataNames.push_back(std::string(quadgrav::BSSN_CONSTRAINT_VAR_NAMES[BSSN_VTU_OUTPUT_CONST_INDICES[i]]));
                pData[numEvolVars+i]=consVar[BSSN_VTU_OUTPUT_CONST_INDICES[i]];
            }

            std::vector<char*> pDataNames_char;
            pDataNames_char.reserve(pDataNames.size());

            for(unsigned int  i = 0; i < pDataNames.size(); i++)
                pDataNames_char.push_back(const_cast<char*>(pDataNames[i].c_str()));

            const char * fDataNames[]= {"Time","Cycle"};
            const double fData[]= {m_uiTinfo._m_uiT,(double)m_uiTinfo._m_uiStep};

            char fPrefix[256];
            sprintf(fPrefix,"%s_%d",quadgrav::BSSN_VTU_FILE_PREFIX.c_str(),m_uiTinfo._m_uiStep);

            io::vtk::mesh2vtuFine(m_uiMesh,fPrefix,2,fDataNames,fData,(numEvolVars+numConstVars),(const char **)&pDataNames_char[0],(const double **)pData);
        #endif


        #ifdef BSSN_EXTRACT_BH_LOCATIONS
            quadgrav::writeBHCoordinates((const ot::Mesh *)m_uiMesh,(const Point *) m_uiBHLoc,2,m_uiTinfo._m_uiStep);
        #endif

        delete [] evolUnzipVar;
        delete [] consUnzipVar;
        delete [] evolVar;
        delete [] consVar;

        return 0;
    }

    int BSSNCtx::write_checkpt()
    {   
        if(m_uiMesh->isActive())
        {
            unsigned int cpIndex;
            (m_uiTinfo._m_uiStep %(2*quadgrav::BSSN_CHECKPT_FREQ)==0) ? cpIndex=0 : cpIndex=1; // to support alternate file writing.
            unsigned int rank=m_uiMesh->getMPIRank();
            unsigned int npes=m_uiMesh->getMPICommSize();

            DendroScalar** eVar = NULL;
            m_uiEVar.Get2DArray(eVar,false);
            

            char fName[256];
            const ot::TreeNode * pNodes=&(*(m_uiMesh->getAllElements().begin()+m_uiMesh->getElementLocalBegin()));
            sprintf(fName,"%s_octree_%d_%d.oct",quadgrav::BSSN_CHKPT_FILE_PREFIX.c_str(),cpIndex,rank);
            io::checkpoint::writeOctToFile(fName,pNodes,m_uiMesh->getNumLocalMeshElements());

            unsigned int numVars=quadgrav::BSSN_NUM_VARS;
            const char ** varNames=quadgrav::BSSN_VAR_NAMES;

            /*for(unsigned int i=0;i<numVars;i++)
            {
                sprintf(fName,"%s_%s_%d_%d.var",fNamePrefix,varNames[i],cpIndex,rank);
                io::checkpoint::writeVecToFile(fName,m_uiMesh,m_uiPrevVar[i]);
            }*/

            sprintf(fName,"%s_%d_%d.var",quadgrav::BSSN_CHKPT_FILE_PREFIX.c_str(),cpIndex,rank);
            io::checkpoint::writeVecToFile(fName,m_uiMesh,(const double **)eVar,quadgrav::BSSN_NUM_VARS);

            if(!rank)
            {
                sprintf(fName,"%s_step_%d.cp",quadgrav::BSSN_CHKPT_FILE_PREFIX.c_str(),cpIndex);
                std::cout<<"[BSSNCtx] \t writing checkpoint file : "<<fName<<std::endl;
                std::ofstream outfile(fName);
                if(!outfile) {
                    std::cout<<fName<<" file open failed "<<std::endl;
                    return 0;
                }

                json checkPoint;
                checkPoint["DENDRO_TS_TIME_BEGIN"]    = m_uiTinfo._m_uiTb;
                checkPoint["DENDRO_TS_TIME_END"]      = m_uiTinfo._m_uiTe;
                checkPoint["DENDRO_TS_ELEMENT_ORDER"] = m_uiElementOrder;

                checkPoint["DENDRO_TS_TIME_CURRENT"]    = m_uiTinfo._m_uiT;
                checkPoint["DENDRO_TS_STEP_CURRENT"]    = m_uiTinfo._m_uiStep;
                checkPoint["DENDRO_TS_TIME_STEP_SIZE"]  = m_uiTinfo._m_uiTh;
                checkPoint["DENDRO_TS_LAST_IO_TIME"]    = m_uiTinfo._m_uiT;

                checkPoint["DENDRO_TS_WAVELET_TOLERANCE"]=quadgrav::BSSN_WAVELET_TOL;
                checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"]=quadgrav::BSSN_LOAD_IMB_TOL;
                checkPoint["DENDRO_TS_NUM_VARS"]=numVars; // number of variables to restore.
                checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"]=m_uiMesh->getMPICommSize(); // (note that rank 0 is always active).
                
                checkPoint["DENDRO_BH1_X"]=m_uiBHLoc[0].x();
                checkPoint["DENDRO_BH1_Y"]=m_uiBHLoc[0].y();
                checkPoint["DENDRO_BH1_Z"]=m_uiBHLoc[0].z();
                
                
                checkPoint["DENDRO_BH2_X"]=m_uiBHLoc[1].x();
                checkPoint["DENDRO_BH2_Y"]=m_uiBHLoc[1].y();
                checkPoint["DENDRO_BH2_Z"]=m_uiBHLoc[1].z();
                

                outfile<<std::setw(4)<<checkPoint<<std::endl;
                outfile.close();

            }
            
            delete [] eVar;
        }
        return 0;
    }

    int BSSNCtx::restore_checkpt()
    {
        unsigned int numVars=0;
        std::vector<ot::TreeNode> octree;
        json checkPoint;

        int rank;
        int npes;
        MPI_Comm comm=m_uiMesh->getMPIGlobalCommunicator();
        MPI_Comm_rank(comm,&rank);
        MPI_Comm_size(comm,&npes);

        unsigned int activeCommSz;

        char fName[256];
        unsigned int restoreStatus=0;
        unsigned int restoreStatusGlobal=0; // 0 indicates successfully restorable.

        ot::Mesh* newMesh;
        unsigned int restoreStep[2];
        restoreStep[0]=0;
        restoreStep[1]=0;

        unsigned int restoreFileIndex=0;

        for(unsigned int cpIndex=0; cpIndex<2; cpIndex++) {

            restoreStatus=0;

            if(!rank)
            {
                sprintf(fName,"%s_step_%d.cp", quadgrav::BSSN_CHKPT_FILE_PREFIX.c_str() , cpIndex);
                std::ifstream infile(fName);
                if(!infile) {
                    std::cout<<fName<<" file open failed "<<std::endl;
                    restoreStatus=1;
                }


                if(restoreStatus==0)
                {
                    infile>>checkPoint;
                    m_uiTinfo._m_uiTb    = checkPoint["DENDRO_TS_TIME_BEGIN"];
                    m_uiTinfo._m_uiTe    = checkPoint["DENDRO_TS_TIME_END"];
                    m_uiTinfo._m_uiT     = checkPoint["DENDRO_TS_TIME_CURRENT"];
                    m_uiTinfo._m_uiStep  = checkPoint["DENDRO_TS_STEP_CURRENT"]; 
                    m_uiTinfo._m_uiTh    = checkPoint["DENDRO_TS_TIME_STEP_SIZE"];
                    m_uiElementOrder     = checkPoint["DENDRO_TS_ELEMENT_ORDER"];
                    
                    quadgrav::BSSN_WAVELET_TOL=checkPoint["DENDRO_TS_WAVELET_TOLERANCE"];
                    quadgrav::BSSN_LOAD_IMB_TOL=checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"];
                    
                    numVars=checkPoint["DENDRO_TS_NUM_VARS"];
                    activeCommSz=checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"];
                    
                    m_uiBHLoc[0]=Point((double)checkPoint["DENDRO_BH1_X"],(double)checkPoint["DENDRO_BH1_Y"],(double)checkPoint["DENDRO_BH1_Z"]);
                    m_uiBHLoc[1]=Point((double)checkPoint["DENDRO_BH2_X"],(double)checkPoint["DENDRO_BH2_Y"],(double)checkPoint["DENDRO_BH2_Z"]);
                    restoreStep[cpIndex]=m_uiTinfo._m_uiStep;

                    }
              }

            }

            if(!rank)
            {
                if(restoreStep[0]<restoreStep[1])
                    restoreFileIndex=1;
                else
                    restoreFileIndex=0;
            }

            par::Mpi_Bcast(&restoreFileIndex,1,0,comm);

            restoreStatus=0;
            octree.clear();
            if(!rank) std::cout<<"[BSSNCtx] :  Trying to restore from checkpoint index : "<<restoreFileIndex<<std::endl;
        
            if(!rank)
            {
                sprintf(fName,"%s_step_%d.cp", quadgrav::BSSN_CHKPT_FILE_PREFIX.c_str(), restoreFileIndex);
                std::ifstream infile(fName);
                if(!infile) {
                    std::cout<<fName<<" file open failed "<<std::endl;
                    restoreStatus=1;
                }


                if(restoreStatus==0)
                {
                    infile>>checkPoint;
                    m_uiTinfo._m_uiTb    = checkPoint["DENDRO_TS_TIME_BEGIN"];
                    m_uiTinfo._m_uiTe    = checkPoint["DENDRO_TS_TIME_END"];
                    m_uiTinfo._m_uiT     = checkPoint["DENDRO_TS_TIME_CURRENT"];
                    m_uiTinfo._m_uiStep  = checkPoint["DENDRO_TS_STEP_CURRENT"]; 
                    m_uiTinfo._m_uiTh    = checkPoint["DENDRO_TS_TIME_STEP_SIZE"];
                    m_uiElementOrder     = checkPoint["DENDRO_TS_ELEMENT_ORDER"];
                    
                    quadgrav::BSSN_WAVELET_TOL=checkPoint["DENDRO_TS_WAVELET_TOLERANCE"];
                    quadgrav::BSSN_LOAD_IMB_TOL=checkPoint["DENDRO_TS_LOAD_IMB_TOLERANCE"];
                    
                    numVars=checkPoint["DENDRO_TS_NUM_VARS"];
                    activeCommSz=checkPoint["DENDRO_TS_ACTIVE_COMM_SZ"];
                    
                    m_uiBHLoc[0]=Point((double)checkPoint["DENDRO_BH1_X"],(double)checkPoint["DENDRO_BH1_Y"],(double)checkPoint["DENDRO_BH1_Z"]);
                    m_uiBHLoc[1]=Point((double)checkPoint["DENDRO_BH2_X"],(double)checkPoint["DENDRO_BH2_Y"],(double)checkPoint["DENDRO_BH2_Z"]);
                    restoreStep[restoreFileIndex]=m_uiTinfo._m_uiStep;
                }


            }

            par::Mpi_Allreduce(&restoreStatus,&restoreStatusGlobal,1,MPI_MAX,comm);
            if(restoreStatusGlobal==1) 
            {
                if(!rank)
                    std::cout<<"[BSSNCtx] : Restore step failed, restore file corrupted. "<<std::endl;
                MPI_Abort(comm,0);
            }


            MPI_Bcast(&m_uiTinfo,sizeof(ts::TSInfo),MPI_BYTE,0,comm);
            par::Mpi_Bcast(&quadgrav::BSSN_WAVELET_TOL,1,0,comm);
            par::Mpi_Bcast(&quadgrav::BSSN_LOAD_IMB_TOL,1,0,comm);

            par::Mpi_Bcast(&numVars,1,0,comm);
            par::Mpi_Bcast(&m_uiElementOrder,1,0,comm);
            par::Mpi_Bcast(&activeCommSz,1,0,comm);
            
            par::Mpi_Bcast(m_uiBHLoc,2,0,comm);

                if(activeCommSz>npes)
                {
                    if(!rank)
                        std::cout<<" [BSSNCtx] : checkpoint file written from  a larger communicator than the current global comm. (i.e. communicator shrinking not allowed in the restore step. )"<<std::endl;
                    
                    MPI_Abort(comm,0);
                }



            bool isActive=(rank<activeCommSz);

            MPI_Comm newComm;
            par::splitComm2way(isActive,&newComm,comm);

            if(isActive) {

                int activeRank;
                int activeNpes;

                MPI_Comm_rank(newComm, &activeRank);
                MPI_Comm_size(newComm, &activeNpes);
                assert(activeNpes == activeCommSz);

                sprintf(fName, "%s_octree_%d_%d.oct", quadgrav::BSSN_CHKPT_FILE_PREFIX.c_str(),restoreFileIndex,activeRank);
                restoreStatus=io::checkpoint::readOctFromFile(fName, octree);
                assert(par::test::isUniqueAndSorted(octree, newComm));

            }

            par::Mpi_Allreduce(&restoreStatus,&restoreStatusGlobal,1,MPI_MAX,comm);
            if(restoreStatusGlobal==1) {

                if(!rank) std::cout<<"[BSSNCtx]: octree (*.oct) restore file is corrupted "<<std::endl;
                MPI_Abort(comm,0);
            }

            newMesh=new ot::Mesh(octree,1,m_uiElementOrder,activeCommSz,comm);
            // no need to transfer data only to resize the contex variables. 
            this->grid_transfer( newMesh , false, false, false);
            
            // only reads the evolution variables. 
            if(isActive) {

                int activeRank;
                int activeNpes;

                DendroScalar** inVec=NULL;
                m_uiEVar.Get2DArray(inVec,false);

                MPI_Comm_rank(newComm, &activeRank);
                MPI_Comm_size(newComm, &activeNpes);
                assert(activeNpes == activeCommSz);

                sprintf(fName,"%s_%d_%d.var",quadgrav::BSSN_CHKPT_FILE_PREFIX.c_str(),restoreFileIndex,activeRank);
                restoreStatus=io::checkpoint::readVecFromFile(fName,newMesh,inVec,quadgrav::BSSN_NUM_VARS);

                delete [] inVec;
            }

            MPI_Comm_free(&newComm);
            par::Mpi_Allreduce(&restoreStatus,&restoreStatusGlobal,1,MPI_MAX,comm);
            if(restoreStatusGlobal==1) {

                if(!rank) std::cout<<"[BSSNCtx]: varible (*.var) restore file currupted "<<std::endl;
                MPI_Abort(comm,0);
            }

            std::swap(m_uiMesh,newMesh);
            delete newMesh;
            
        unsigned int localSz=m_uiMesh->getNumLocalMeshElements();
        unsigned int totalElems=0;
        par::Mpi_Allreduce(&localSz, &totalElems ,1,MPI_SUM,comm);
        
        if(!rank) std::cout<<" checkpoint at step : "<<m_uiTinfo._m_uiStep<<"active Comm. sz: "<<activeCommSz<<" restore successful: "<<" restored mesh size: "<<totalElems<<std::endl;

        m_uiIsETSSynced = false;
        return 0;

    }

    int BSSNCtx::pre_timestep(DVec sIn)
    {

        #ifdef BSSN_EXTRACT_BH_LOCATIONS
            DendroScalar ** evar = new DendroScalar*[quadgrav::BSSN_NUM_VARS];
            Point bhLoc[2];
            sIn.Get2DArray(evar,true);
            quadgrav::extractBHCoords((const ot::Mesh *)m_uiMesh,(const DendroScalar*)evar[BHLOC::EXTRACTION_VAR_ID],BHLOC::EXTRACTION_TOL,(const Point *) m_uiBHLoc,2,(Point*)bhLoc);
            
            m_uiBHLoc[0] = bhLoc[0];
            m_uiBHLoc[1] = bhLoc[1];

            delete [] evar;
        #endif

        return 0;
    }

    int BSSNCtx::pre_stage(DVec  sIn)
    {
        
        return 0;
    }

    int BSSNCtx::post_stage(DVec sIn)
    { 
        return 0;
    }

    int BSSNCtx::post_timestep(DVec sIn)
    {
        return 0;
    }

    bool BSSNCtx::is_remesh()
    {
        bool isRefine = false;
        if(quadgrav::BSSN_ENABLE_BLOCK_ADAPTIVITY)
            return false;

        
        MPI_Comm comm = m_uiMesh->getMPIGlobalCommunicator();
        
        this->unzip(m_uiEVar,m_uiEUnzip[0],quadgrav::BSSN_ASYNC_COMM_K);

        DendroScalar** unzipVar;
        m_uiEUnzip[0].Get2DArray(unzipVar, false);

        unsigned int refineVarIds[quadgrav::BSSN_NUM_REFINE_VARS];
        for(unsigned int vIndex=0; vIndex<quadgrav::BSSN_NUM_REFINE_VARS; vIndex++)
            refineVarIds[vIndex]=quadgrav::BSSN_REFINE_VARIABLE_INDICES[vIndex];

        double wTol=quadgrav::BSSN_WAVELET_TOL;
            std::function<double(double,double,double)> waveletTolFunc =[wTol](double x,double y, double z) {
            return quadgrav::computeWTol(x,y,z,wTol);
        };
        

        if(quadgrav::BSSN_REFINEMENT_MODE == quadgrav::RefinementMode::WAMR) {
            isRefine=m_uiMesh->isReMeshUnzip((const double **)unzipVar,refineVarIds,quadgrav::BSSN_NUM_REFINE_VARS,waveletTolFunc,quadgrav::BSSN_DENDRO_AMR_FAC); 

        }else if(quadgrav::BSSN_REFINEMENT_MODE == quadgrav::RefinementMode::EH)
        {
            isRefine = quadgrav::isRemeshEH(m_uiMesh,(const double **)unzipVar,quadgrav::VAR::U_ALPHA,quadgrav::BSSN_EH_REFINE_VAL,quadgrav::BSSN_EH_COARSEN_VAL,true);

        }else if(quadgrav::BSSN_REFINEMENT_MODE == quadgrav::RefinementMode::EH_WAMR)
        {
            const bool isR1 = m_uiMesh->isReMeshUnzip((const double **)unzipVar,refineVarIds,quadgrav::BSSN_NUM_REFINE_VARS,waveletTolFunc,quadgrav::BSSN_DENDRO_AMR_FAC); 
            const bool isR2 = quadgrav::isRemeshEH(m_uiMesh,(const double **)unzipVar,quadgrav::VAR::U_ALPHA,quadgrav::BSSN_EH_REFINE_VAL,quadgrav::BSSN_EH_COARSEN_VAL,false);

            isRefine = (isR1 || isR2);
        }

        delete [] unzipVar;
        
        return isRefine;

    }

    int BSSNCtx::update_app_vars() 
    {
        m_uiEVar = m_uiEvolutionVar[0];
        m_uiCVar = m_uiConstrainedVar[0];

        m_uiEUnzip[0] = m_uiEvolutionUnzipVar[0];
        m_uiEUnzip[1] = m_uiEvolutionUnzipVar[1];

        m_uiCUnzip[0] = m_uiConstraintUnzipVar[0];
        
    }

    DVec BSSNCtx::get_evolution_vars()
    {
        return m_uiEVar;
    }
    
    DVec BSSNCtx::get_constraint_vars()
    {
        return m_uiCVar;
    }

    DVec BSSNCtx::get_primitive_vars()
    {
        return m_uiPVar;
    }

    int BSSNCtx::terminal_output()
    {
        DendroScalar min=0, max=0;
        m_uiEVar.VecMinMax(m_uiMesh,min,max,quadgrav::VAR::U_ALPHA);

        if(m_uiMesh->isActive())
        {
            if(!(m_uiMesh->getMPIRank()))
                std::cout<<"[BSSNCtx]:  "<<quadgrav::BSSN_VAR_NAMES[quadgrav::VAR::U_ALPHA]<<" (min,max) : \t ( "<<min<<", "<<max<<" ) "<<std::endl;
        }

        return 0; 
    }


}// end of namespace quadgrav. 
