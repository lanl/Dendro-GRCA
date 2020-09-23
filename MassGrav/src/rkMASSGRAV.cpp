//
// Created by milinda on 1/16/19.
//

/**
 * @brief contains RK time stepper for MASSGRAV equations.
 * @author Milinda Fernando
 * School of Computing, University of Utah
 *
 * */

#include "rkMASSGRAV.h"



namespace ode
{
namespace solver
{

RK_MASSGRAV::RK_MASSGRAV(ot::Mesh *pMesh, DendroScalar pTBegin, DendroScalar pTEnd,DendroScalar pTh,RKType rkType): RK(pMesh,pTBegin,pTEnd,pTh)
{

    m_uiRKType=rkType;

    m_uiBHLoc[0]=Point(massgrav::BH1.getBHCoordX(),massgrav::BH1.getBHCoordY(),massgrav::BH1.getBHCoordZ());
    m_uiBHLoc[1]=Point(massgrav::BH2.getBHCoordX(),massgrav::BH2.getBHCoordY(),massgrav::BH2.getBHCoordZ());


    // allocate memory for the variables.
    m_uiVar=new DendroScalar*[massgrav::MASSGRAV_NUM_VARS];
    for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
        m_uiVar[index]=m_uiMesh->createVector<DendroScalar>();

    m_uiPrevVar=new DendroScalar*[massgrav::MASSGRAV_NUM_VARS];
    for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
        m_uiPrevVar[index]=m_uiMesh->createVector<DendroScalar>();

    m_uiVarIm=new DendroScalar*[massgrav::MASSGRAV_NUM_VARS];
    for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
        m_uiVarIm[index]=m_uiMesh->createVector<DendroScalar>();

    if(m_uiRKType==RKType::RK3)
        m_uiNumRKStages=massgrav::MASSGRAV_RK3_STAGES;
    else if(m_uiRKType==RKType::RK4)
        m_uiNumRKStages=massgrav::MASSGRAV_RK4_STAGES;
    else if(m_uiRKType==RKType::RK45)
        m_uiNumRKStages=massgrav::MASSGRAV_RK45_STAGES;
    else
    {
        if(!(pMesh->getMPIRankGlobal()))
            std::cout<<"[RK Solver Error]: undefined rk solver type"<<std::endl;

        exit(0);
    }




    m_uiStage=new DendroScalar**[m_uiNumRKStages];
    for(unsigned int stage=0; stage<m_uiNumRKStages; stage++)
    {
        m_uiStage[stage]=new DendroScalar*[massgrav::MASSGRAV_NUM_VARS];
        for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
            m_uiStage[stage][index]=m_uiMesh->createVector<DendroScalar>();
    }

    m_uiUnzipVar=new DendroScalar*[massgrav::MASSGRAV_NUM_VARS];
    for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
        m_uiUnzipVar[index]=m_uiMesh->createUnZippedVector<DendroScalar>();


    m_uiUnzipVarRHS=new DendroScalar*[massgrav::MASSGRAV_NUM_VARS];
    for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
        m_uiUnzipVarRHS[index]=m_uiMesh->createUnZippedVector<DendroScalar>();



    // allocate memory for the constraint variables.
    m_uiConstraintVars=new DendroScalar*[massgrav::MASSGRAV_CONSTRAINT_NUM_VARS];
    for(unsigned int index=0; index<massgrav::MASSGRAV_CONSTRAINT_NUM_VARS; index++)
        m_uiConstraintVars[index]=m_uiMesh->createVector<DendroScalar>();


    m_uiUnzipConstraintVars=new DendroScalar*[massgrav::MASSGRAV_CONSTRAINT_NUM_VARS];
    for(unsigned int index=0; index<massgrav::MASSGRAV_CONSTRAINT_NUM_VARS; index++)
        m_uiUnzipConstraintVars[index]=m_uiMesh->createUnZippedVector<DendroScalar>();


    // mpi communication
    m_uiSendNodeBuf=new DendroScalar*[massgrav::MASSGRAV_ASYNC_COMM_K];
    m_uiRecvNodeBuf=new DendroScalar*[massgrav::MASSGRAV_ASYNC_COMM_K];

    m_uiSendReqs=new MPI_Request*[massgrav::MASSGRAV_ASYNC_COMM_K];
    m_uiRecvReqs=new MPI_Request*[massgrav::MASSGRAV_ASYNC_COMM_K];
    m_uiSendSts=new MPI_Status*[massgrav::MASSGRAV_ASYNC_COMM_K];
    m_uiRecvSts=new MPI_Status*[massgrav::MASSGRAV_ASYNC_COMM_K];

    for(unsigned int index=0; index<massgrav::MASSGRAV_ASYNC_COMM_K; index++)
    {
        m_uiSendNodeBuf[index]=NULL;
        m_uiRecvNodeBuf[index]=NULL;

        m_uiSendReqs[index]=NULL;
        m_uiRecvReqs[index]=NULL;
        m_uiSendSts[index]=NULL;
        m_uiRecvSts[index]=NULL;
    }

    if(m_uiMesh->isActive())
    {
        // allocate mpi comm. reqs and status
        for(unsigned int index=0; index<massgrav::MASSGRAV_ASYNC_COMM_K; index++)
        {
            if(m_uiMesh->getGhostExcgTotalSendNodeCount()!=0) m_uiSendNodeBuf[index]=new DendroScalar[m_uiMesh->getGhostExcgTotalSendNodeCount()];
            if(m_uiMesh->getGhostExcgTotalRecvNodeCount()!=0) m_uiRecvNodeBuf[index]=new DendroScalar[m_uiMesh->getGhostExcgTotalRecvNodeCount()];

            if(m_uiMesh->getSendProcListSize()!=0)
            {
                m_uiSendReqs[index]=new MPI_Request[m_uiMesh->getSendProcListSize()];
                m_uiSendSts[index]=new MPI_Status[m_uiMesh->getSendProcListSize()];
            }

            if(m_uiMesh->getRecvProcListSize()!=0)
            {
                m_uiRecvReqs[index]=new MPI_Request[m_uiMesh->getRecvProcListSize()];
                m_uiRecvSts[index]=new MPI_Status[m_uiMesh->getRecvProcListSize()];
            }


        }

    }


}

RK_MASSGRAV::~RK_MASSGRAV()
{
    for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
    {
        delete [] m_uiVar[index];
        delete [] m_uiPrevVar[index];
        delete [] m_uiVarIm[index];
        delete [] m_uiUnzipVar[index];
        delete [] m_uiUnzipVarRHS[index];


    }


    delete [] m_uiVar;
    delete [] m_uiPrevVar;
    delete [] m_uiVarIm;
    delete [] m_uiUnzipVar;
    delete [] m_uiUnzipVarRHS;

    for(unsigned int stage=0; stage<m_uiNumRKStages; stage++)
        for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
            delete [] m_uiStage[stage][index];

    for(unsigned int stage=0; stage<m_uiNumRKStages; stage++)
        delete [] m_uiStage[stage];

    delete [] m_uiStage;


    // deallocate memory for the constraint variables.
    for(unsigned int index=0; index<massgrav::MASSGRAV_CONSTRAINT_NUM_VARS; index++)
    {
        delete [] m_uiConstraintVars[index];
        delete [] m_uiUnzipConstraintVars[index];
    }


    delete [] m_uiConstraintVars;
    delete [] m_uiUnzipConstraintVars;


    // mpi communication
    for(unsigned int index=0; index<massgrav::MASSGRAV_ASYNC_COMM_K; index++)
    {
        delete [] m_uiSendNodeBuf[index];
        delete [] m_uiRecvNodeBuf[index];

        delete [] m_uiSendReqs[index];
        delete [] m_uiRecvReqs[index];

        delete [] m_uiSendSts[index];
        delete [] m_uiRecvSts[index];

    }

    delete [] m_uiSendNodeBuf;
    delete [] m_uiRecvNodeBuf;

    delete [] m_uiSendReqs;
    delete [] m_uiSendSts;
    delete [] m_uiRecvReqs;
    delete [] m_uiRecvSts;



}


void RK_MASSGRAV::applyInitialConditions(DendroScalar** zipIn)
{
    unsigned int nodeLookUp_CG;
    unsigned int nodeLookUp_DG;
    double x,y,z,len;
    const ot::TreeNode * pNodes=&(*(m_uiMesh->getAllElements().begin()));
    unsigned int ownerID,ii_x,jj_y,kk_z;
    unsigned int eleOrder=m_uiMesh->getElementOrder();
    const unsigned int * e2n_cg=&(*(m_uiMesh->getE2NMapping().begin()));
    const unsigned int * e2n_dg=&(*(m_uiMesh->getE2NMapping_DG().begin()));
    const unsigned int nPe=m_uiMesh->getNumNodesPerElement();
    const unsigned int nodeLocalBegin=m_uiMesh->getNodeLocalBegin();
    const unsigned int nodeLocalEnd=m_uiMesh->getNodeLocalEnd();


    double* var=new double[massgrav::MASSGRAV_NUM_VARS];

    double mp, mm, mp_adm, mm_adm, E, J1, J2, J3;
    // set the TP communicator. 
    MPI_TP_COMM = m_uiComm;

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
                        
                        if (massgrav::MASSGRAV_ID_TYPE == 0) {
                            TwoPunctures((double)x,(double)y,(double)z,var,
                                         &mp, &mm, &mp_adm, &mm_adm, &E, &J1, &J2, &J3);
                        }
                        else if (massgrav::MASSGRAV_ID_TYPE == 1) {
                            massgrav::punctureData((double)x,(double)y,(double)z,var);
                        }
                        else if (massgrav::MASSGRAV_ID_TYPE == 2) {
                            massgrav::KerrSchildData((double)x,(double)y,(double)z,var);
                        }
                        else if (massgrav::MASSGRAV_ID_TYPE == 3) {
                            massgrav::noiseData((double)x,(double)y,(double)z,var);
                        }
                        else if (massgrav::MASSGRAV_ID_TYPE == 4) {
                            massgrav::fake_initial_data((double)x,(double)y,(double)z,var);
                        }
                        else if (massgrav::MASSGRAV_ID_TYPE == 5) {
                            massgrav::FlatMinkowski((double)x,(double)y,(double)z,var);
                        }
                        else {
                            std::cout<<"Unknown ID type"<<std::endl;
                        }
                        for(unsigned int v=0; v<massgrav::MASSGRAV_NUM_VARS; v++)
                            zipIn[v][nodeLookUp_CG]=var[v];


                    }

                }

    }


    for(unsigned int node=m_uiMesh->getNodeLocalBegin(); node<m_uiMesh->getNodeLocalEnd(); node++)
        enforce_massgrav_constraints(zipIn,node);


    delete [] var;

}

void RK_MASSGRAV::initialGridConverge()
{

    applyInitialConditions(m_uiPrevVar);

    bool isRefine=false;
    unsigned int oldElements,oldElements_g;
    unsigned int newElements,newElements_g;

    // refine based on all the variables
    const unsigned int refineNumVars=massgrav::MASSGRAV_NUM_REFINE_VARS;
    unsigned int refineVarIds[refineNumVars];
    for(unsigned int vIndex=0; vIndex<refineNumVars; vIndex++)
        refineVarIds[vIndex]=massgrav::MASSGRAV_REFINE_VARIABLE_INDICES[vIndex];

    double wTol=massgrav::MASSGRAV_WAVELET_TOL;
    std::function<double(double,double,double)> waveletTolFunc =[wTol](double x,double y, double z) {
        return massgrav::computeWTol(x,y,z,wTol);
    };
    unsigned int iterCount=1;
    do
    {

        #ifdef RK_SOLVER_OVERLAP_COMM_AND_COMP
        unzipVars_async(m_uiPrevVar,m_uiUnzipVar);
        #else
        performGhostExchangeVars(m_uiPrevVar);
        unzipVars(m_uiPrevVar,m_uiUnzipVar);
        #endif


        if(massgrav::MASSGRAV_ENABLE_BLOCK_ADAPTIVITY)
            isRefine=false;
        else
            isRefine=massgrav::isRemeshEH(m_uiMesh,(const double **)m_uiUnzipVar,massgrav::VAR::U_ALPHA,massgrav::MASSGRAV_EH_REFINE_VAL,massgrav::MASSGRAV_EH_COARSEN_VAL);
            //m_uiMesh->isReMeshUnzip((const DendroScalar **)m_uiUnzipVar,refineVarIds,refineNumVars,waveletTolFunc,massgrav::MASSGRAV_DENDRO_AMR_FAC);

        if(isRefine)
        {
            ot::Mesh* newMesh=m_uiMesh->ReMesh(massgrav::MASSGRAV_DENDRO_GRAIN_SZ,massgrav::MASSGRAV_LOAD_IMB_TOL,massgrav::MASSGRAV_SPLIT_FIX);

            oldElements=m_uiMesh->getNumLocalMeshElements();
            newElements=newMesh->getNumLocalMeshElements();

            par::Mpi_Allreduce(&oldElements,&oldElements_g,1,MPI_SUM,m_uiMesh->getMPIGlobalCommunicator());
            par::Mpi_Allreduce(&newElements,&newElements_g,1,MPI_SUM,m_uiMesh->getMPIGlobalCommunicator());

            if(!(m_uiMesh->getMPIRankGlobal()))std::cout<<"initial grid iteration : "<<iterCount<<" old mesh: "<<oldElements_g<<" new mesh: "<<newElements_g<<std::endl;


            // performs the inter-grid transfer
            intergridTransferVars(m_uiPrevVar,newMesh,massgrav::MASSGRAV_USE_FD_GRID_TRANSFER);

            for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
            {
                delete [] m_uiVar[index];
                delete [] m_uiVarIm[index];
                delete [] m_uiUnzipVar[index];
                delete [] m_uiUnzipVarRHS[index];

                m_uiVar[index]=NULL;
                m_uiVarIm[index]=NULL;
                m_uiUnzipVar[index]=NULL;
                m_uiUnzipVarRHS[index]=NULL;

                m_uiVar[index]=newMesh->createVector<DendroScalar>();
                m_uiVarIm[index]=newMesh->createVector<DendroScalar>();
                m_uiUnzipVar[index]=newMesh->createUnZippedVector<DendroScalar>();
                m_uiUnzipVarRHS[index]=newMesh->createUnZippedVector<DendroScalar>();


            }

            for(unsigned int stage=0; stage<m_uiNumRKStages; stage++)
                for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
                {
                    delete [] m_uiStage[stage][index];
                    m_uiStage[stage][index]=NULL;
                    m_uiStage[stage][index]=newMesh->createVector<DendroScalar>();
                }

            // deallocate constraint vars allocate them for the new mesh.
            for(unsigned int index=0; index<massgrav::MASSGRAV_CONSTRAINT_NUM_VARS; index++)
            {
                delete [] m_uiConstraintVars[index];
                delete [] m_uiUnzipConstraintVars[index];

                m_uiConstraintVars[index]=newMesh->createVector<DendroScalar>();
                m_uiUnzipConstraintVars[index]=newMesh->createUnZippedVector<DendroScalar>();

            }


            std::swap(newMesh,m_uiMesh);
            delete newMesh;

            #ifdef RK_SOLVER_OVERLAP_COMM_AND_COMP
            // reallocates mpi resources for the the new mesh. (this will deallocate the old resources)
            reallocateMPIResources();
            #endif

            if(m_uiMesh->isActive())
            {
                DendroScalar l_min=vecMin(m_uiPrevVar[massgrav::VAR::U_ALPHA]+m_uiMesh->getNodeLocalBegin(),(m_uiMesh->getNumLocalMeshNodes()),m_uiMesh->getMPICommunicator());
                DendroScalar l_max=vecMax(m_uiPrevVar[massgrav::VAR::U_ALPHA]+m_uiMesh->getNodeLocalBegin(),(m_uiMesh->getNumLocalMeshNodes()),m_uiMesh->getMPICommunicator());
                if(!(m_uiMesh->getMPIRank())) {
                    std::cout << "transfer completed:    ||VAR::U_ALPHA|| (min, max) : ("<<l_min<<", "<<l_max<<" ) "<<std::endl;
                }

            }

            iterCount+=1;

        }

        applyInitialConditions(m_uiPrevVar);

    } while(isRefine && (newElements_g!=oldElements_g));


}

void RK_MASSGRAV::reallocateMPIResources()
{
    for(unsigned int index=0; index<massgrav::MASSGRAV_ASYNC_COMM_K; index++)
    {
        delete [] m_uiSendNodeBuf[index];
        delete [] m_uiRecvNodeBuf[index];

        delete [] m_uiSendReqs[index];
        delete [] m_uiRecvReqs[index];

        delete [] m_uiSendSts[index];
        delete [] m_uiRecvSts[index];

    }

    for(unsigned int index=0; index<massgrav::MASSGRAV_ASYNC_COMM_K; index++)
    {
        m_uiSendNodeBuf[index]=NULL;
        m_uiRecvNodeBuf[index]=NULL;

        m_uiSendReqs[index]=NULL;
        m_uiRecvReqs[index]=NULL;
        m_uiSendSts[index]=NULL;
        m_uiRecvSts[index]=NULL;
    }

    if(m_uiMesh->isActive())
    {
        // allocate mpi comm. reqs and status
        for(unsigned int index=0; index<massgrav::MASSGRAV_ASYNC_COMM_K; index++)
        {
            if(m_uiMesh->getGhostExcgTotalSendNodeCount()!=0) m_uiSendNodeBuf[index]=new DendroScalar[m_uiMesh->getGhostExcgTotalSendNodeCount()];
            if(m_uiMesh->getGhostExcgTotalRecvNodeCount()!=0) m_uiRecvNodeBuf[index]=new DendroScalar[m_uiMesh->getGhostExcgTotalRecvNodeCount()];

            if(m_uiMesh->getSendProcListSize()!=0)
            {
                m_uiSendReqs[index]=new MPI_Request[m_uiMesh->getSendProcListSize()];
                m_uiSendSts[index]=new MPI_Status[m_uiMesh->getSendProcListSize()];
            }

            if(m_uiMesh->getRecvProcListSize()!=0)
            {
                m_uiRecvReqs[index]=new MPI_Request[m_uiMesh->getRecvProcListSize()];
                m_uiRecvSts[index]=new MPI_Status[m_uiMesh->getRecvProcListSize()];
            }


        }

    }


}

void RK_MASSGRAV::writeToVTU(DendroScalar **evolZipVarIn, DendroScalar ** constrZipVarIn, unsigned int numEvolVars,unsigned int numConstVars,const unsigned int * evolVarIndices, const unsigned int * constVarIndices, bool zslice)
{
    massgrav::timer::t_ioVtu.start();

    std::vector<std::string> pDataNames;
    double *pData[(numConstVars+numEvolVars)];

    for(unsigned int i=0; i<numEvolVars; i++)
    {
        pDataNames.push_back(std::string(massgrav::MASSGRAV_VAR_NAMES[evolVarIndices[i]]));
        pData[i]=evolZipVarIn[evolVarIndices[i]];
    }


    for(unsigned int i=0; i<numConstVars; i++)
    {
        pDataNames.push_back(std::string(massgrav::MASSGRAV_CONSTRAINT_VAR_NAMES[constVarIndices[i]]));
        pData[numEvolVars+i]=constrZipVarIn[constVarIndices[i]];
    }

    std::vector<char*> pDataNames_char;
    pDataNames_char.reserve(pDataNames.size());

    for(unsigned int  i = 0; i < pDataNames.size(); i++)
        pDataNames_char.push_back(const_cast<char*>(pDataNames[i].c_str()));

    const char * fDataNames[]= {"Time","Cycle"};
    const double fData[]= {m_uiCurrentTime,(double)m_uiCurrentStep};

    char fPrefix[256];
    sprintf(fPrefix,"%s_%d",massgrav::MASSGRAV_VTU_FILE_PREFIX.c_str(),m_uiCurrentStep);

    if(zslice)
    {
        unsigned int s_val[3]= {1u<<(m_uiMaxDepth-1), 1u<<(m_uiMaxDepth-1), 1u<<(m_uiMaxDepth-1)};
        unsigned int s_norm[3] ={0,0,1};
        io::vtk::mesh2vtu_slice(m_uiMesh,s_val, s_norm, fPrefix,2,fDataNames,fData,(numEvolVars+numConstVars),(const char **)&pDataNames_char[0],(const double **)pData);

    }else{
        io::vtk::mesh2vtuFine(m_uiMesh,fPrefix,2,fDataNames,fData,(numEvolVars+numConstVars),(const char **)&pDataNames_char[0],(const double **)pData);
    }
    

    massgrav::timer::t_ioVtu.stop();

}

void RK_MASSGRAV::performGhostExchangeVars(DendroScalar** zipIn)
{

    massgrav::timer::t_ghostEx_sync.start();

    for(unsigned int v=0; v<massgrav::MASSGRAV_NUM_VARS; v++)
        m_uiMesh->performGhostExchange(zipIn[v]);

    massgrav::timer::t_ghostEx_sync.stop();

}

void RK_MASSGRAV::intergridTransferVars(DendroScalar**& zipIn, const ot::Mesh* pnewMesh,bool useUnzip)
{
    massgrav::timer::t_gridTransfer.start();

    if(!useUnzip)
    {
        for(unsigned int v=0; v<massgrav::MASSGRAV_NUM_VARS; v++)
            m_uiMesh->interGridTransfer(zipIn[v],pnewMesh);
    }else
    {
        unzipVars_async(zipIn,m_uiUnzipVar);
        for(unsigned int v=0; v<massgrav::MASSGRAV_NUM_VARS; v++)
        {
            m_uiMesh->interGridTransferUnzip(m_uiUnzipVar[v],zipIn[v],pnewMesh);
        }
        
    }



    

    massgrav::timer::t_gridTransfer.stop();

}

void RK_MASSGRAV::unzipVars(DendroScalar** zipIn, DendroScalar** uzipOut)
{
    massgrav::timer::t_unzip_sync.start();

    for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
        m_uiMesh->unzip(zipIn[index],uzipOut[index]);

    massgrav::timer::t_unzip_sync.stop();

}

void RK_MASSGRAV::unzipVars_async(DendroScalar ** zipIn, DendroScalar **uzipOut)
{

    massgrav::timer::t_unzip_async.start();

    for(unsigned int var=0; var<massgrav::MASSGRAV_NUM_VARS; var+=massgrav::MASSGRAV_ASYNC_COMM_K) {

        for(unsigned int i=0; (i<massgrav::MASSGRAV_ASYNC_COMM_K); i++)
            m_uiMesh->ghostExchangeStart(zipIn[var+i],m_uiSendNodeBuf[i],m_uiRecvNodeBuf[i],m_uiSendReqs[i],m_uiRecvReqs[i]);

        for(unsigned int i=0; (i<massgrav::MASSGRAV_ASYNC_COMM_K); i++)
        {
            m_uiMesh->ghostExchangeRecvSync(zipIn[var + i], m_uiRecvNodeBuf[i],m_uiRecvReqs[i], m_uiRecvSts[i]);
            m_uiMesh->unzip(zipIn[var+i],uzipOut[var+i]);
        }

        for(unsigned int i=0; (i<massgrav::MASSGRAV_ASYNC_COMM_K); i++)
            m_uiMesh->ghostExchangeSendSync(m_uiSendReqs[i], m_uiSendSts[i]);

    }

    massgrav::timer::t_unzip_async.stop();


}


void RK_MASSGRAV::zipVars(DendroScalar** uzipIn, DendroScalar** zipOut)
{
    massgrav::timer::t_zip.start();

    for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
        m_uiMesh->zip(uzipIn[index],zipOut[index]);

    massgrav::timer::t_zip.stop();

}


void RK_MASSGRAV::applyBoundaryConditions()
{

}

void RK_MASSGRAV::performSingleIteration()
{

    char frawName[256];


    if(m_uiMesh->isActive())
    {
        double current_t=m_uiCurrentTime;
        double current_t_adv=current_t;
        #if 0
            sprintf(frawName,"rkU_%d",3*m_uiCurrentStep);
            io::varToRawData((const ot::Mesh*)m_uiMesh,(const double **)m_uiPrevVar,massgrav::MASSGRAV_NUM_VARS,NULL,frawName);
        #endif

        #ifdef RK_SOLVER_OVERLAP_COMM_AND_COMP
        unzipVars_async(m_uiPrevVar,m_uiUnzipVar);
        #else
        //1. perform ghost exchange.
        performGhostExchangeVars(m_uiPrevVar);
        //2. unzip all the variables.
        unzipVars(m_uiPrevVar,m_uiUnzipVar);
        #endif


        int rank =m_uiMesh->getMPIRank();



        const unsigned int nodeLocalBegin=m_uiMesh->getNodeLocalBegin();
        const unsigned int nodeLocalEnd=m_uiMesh->getNodeLocalEnd();

        const std::vector<ot::Block> blkList=m_uiMesh->getLocalBlockList();
        unsigned int offset;
        double ptmin[3], ptmax[3];
        unsigned int sz[3];
        unsigned int bflag;
        double dx,dy,dz;
        const Point pt_min(massgrav::MASSGRAV_COMPD_MIN[0],massgrav::MASSGRAV_COMPD_MIN[1],massgrav::MASSGRAV_COMPD_MIN[2]);
        const Point pt_max(massgrav::MASSGRAV_COMPD_MAX[0],massgrav::MASSGRAV_COMPD_MAX[1],massgrav::MASSGRAV_COMPD_MAX[2]);

        if(m_uiRKType==RKType::RK3)
        {   
            // initial unzip and ghost exchange happens at the rkSolve class.  
            massgravRHS(m_uiUnzipVarRHS,(const DendroScalar **)m_uiUnzipVar,&(*(blkList.begin())),blkList.size());
            zipVars(m_uiUnzipVarRHS,m_uiStage[0]);
                        
            for(unsigned int node=nodeLocalBegin; node<nodeLocalEnd; node++)
            {

               for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
               {
                   m_uiStage[0][index][node]=m_uiPrevVar[index][node] + m_uiT_h * m_uiStage[0][index][node];
               }
               enforce_massgrav_constraints(m_uiStage[0], node);
            }
            
            #if 0            
            sprintf(frawName,"rkU_%d",3*m_uiCurrentStep + 1);
            io::varToRawData((const ot::Mesh*)m_uiMesh,(const double **)m_uiStage[0],massgrav::MASSGRAV_NUM_VARS,NULL,frawName);
            #endif            
            
            #ifdef RK_SOLVER_OVERLAP_COMM_AND_COMP
                unzipVars_async(m_uiStage[0],m_uiUnzipVar);
            #else
                performGhostExchangeVars(m_uiStage[0]);
                unzipVars(m_uiStage[0],m_uiUnzipVar);
            #endif
            
            massgravRHS(m_uiUnzipVarRHS,(const DendroScalar **)m_uiUnzipVar,&(*(blkList.begin())),blkList.size());
            zipVars(m_uiUnzipVarRHS,m_uiStage[1]);
            
            for(unsigned int node=nodeLocalBegin; node<nodeLocalEnd; node++)
            {

               for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
               {
                   m_uiStage[1][index][node]=(0.75)*m_uiPrevVar[index][node] + 
                   0.25*m_uiStage[0][index][node]+m_uiT_h * 0.25 *  m_uiStage[1][index][node];
               }
               enforce_massgrav_constraints(m_uiStage[1], node);
            }
            
            #if 0
            sprintf(frawName,"rkU_%d",3*m_uiCurrentStep + 2);
            io::varToRawData((const ot::Mesh*)m_uiMesh,(const double **)m_uiStage[1],massgrav::MASSGRAV_NUM_VARS,NULL,frawName);
            #endif

            #ifdef RK_SOLVER_OVERLAP_COMM_AND_COMP
                unzipVars_async(m_uiStage[1],m_uiUnzipVar);
            #else
                performGhostExchangeVars(m_uiStage[1]);
                unzipVars(m_uiStage[1],m_uiUnzipVar);
            #endif            
            
            
                
            massgravRHS(m_uiUnzipVarRHS,(const DendroScalar **)m_uiUnzipVar,&(*(blkList.begin())),blkList.size());
            zipVars(m_uiUnzipVarRHS,m_uiVar);
            
            
            for(unsigned int node=nodeLocalBegin; node<nodeLocalEnd; node++)
            {

               for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
               {
                   m_uiVar[index][node]=(1.0/3.0)*m_uiPrevVar[index][node] + 
                   (2.0/3.0)*m_uiStage[1][index][node]+m_uiT_h * (2.0/3.0) *  m_uiVar[index][node];
               }
               enforce_massgrav_constraints(m_uiVar, node);
            }
            
            /*
             * 
             *  !!!! this does not work for some reason . 
            // stage 0   f(u_k)
            massgravRHS(m_uiUnzipVarRHS,(const DendroScalar **)m_uiUnzipVar,&(*(blkList.begin())),blkList.size());    
            zipVars(m_uiUnzipVarRHS,m_uiStage[0]);
            #if 0
            sprintf(frawName,"rk3_step_%d_stage_%d",m_uiCurrentStep,0);
            io::varToRawData((const ot::Mesh*)m_uiMesh,(const double **)m_uiStage[0],massgrav::MASSGRAV_NUM_VARS,NULL,frawName);    
            #endif            
            // u1
            for(unsigned int node=nodeLocalBegin; node<nodeLocalEnd; node++)
            {

               for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
               {
                   m_uiVarIm[index][node]=m_uiPrevVar[index][node] + m_uiT_h * m_uiStage[0][index][node];
               }
               enforce_massgrav_constraints(m_uiVarIm, node);
            }
            
            #if 0
            sprintf(frawName,"rkU_%d",3*m_uiCurrentStep+(1));
            io::varToRawData((const ot::Mesh*)m_uiMesh,(const double **)m_uiVarIm,massgrav::MASSGRAV_NUM_VARS,NULL,frawName);
            #endif            
            current_t_adv=current_t+m_uiT_h;

            #ifdef RK_SOLVER_OVERLAP_COMM_AND_COMP
                unzipVars_async(m_uiVarIm,m_uiUnzipVar);
            #else
                performGhostExchangeVars(m_uiVarIm);
                unzipVars(m_uiVarIm,m_uiUnzipVar);
            #endif

                
            // stage 1
            massgravRHS(m_uiUnzipVarRHS,(const DendroScalar **)m_uiUnzipVar,&(*(blkList.begin())),blkList.size());    
            zipVars(m_uiUnzipVarRHS,m_uiStage[1]);
            #if 0      
            sprintf(frawName,"rk3_step_%d_stage_%d",m_uiCurrentStep,1);
            io::varToRawData((const ot::Mesh*)m_uiMesh,(const double **)m_uiStage[1],massgrav::MASSGRAV_NUM_VARS,NULL,frawName);
            #endif

            // u2
            for(unsigned int node=nodeLocalBegin; node<nodeLocalEnd; node++)
            {

               for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
               {
                   m_uiVarIm[index][node]=m_uiPrevVar[index][node] + m_uiT_h * 0.25 * (m_uiStage[0][index][node] + m_uiStage[1][index][node]);
               }
               enforce_massgrav_constraints(m_uiVarIm, node);
            }
            
            #if 0
            sprintf(frawName,"rkU_%d",3*m_uiCurrentStep+(2));
            io::varToRawData((const ot::Mesh*)m_uiMesh,(const double **)m_uiVarIm,massgrav::MASSGRAV_NUM_VARS,NULL,frawName);
            #endif            
            current_t_adv=current_t+m_uiT_h;

            #ifdef RK_SOLVER_OVERLAP_COMM_AND_COMP
                unzipVars_async(m_uiVarIm,m_uiUnzipVar);
            #else
                performGhostExchangeVars(m_uiVarIm);
                unzipVars(m_uiVarIm,m_uiUnzipVar);
            #endif
            
            // stage 2.     
            massgravRHS(m_uiUnzipVarRHS,(const DendroScalar **)m_uiUnzipVar,&(*(blkList.begin())),blkList.size());    
            zipVars(m_uiUnzipVarRHS,m_uiStage[2]);

            #if 0
            sprintf(frawName,"rk3_step_%d_stage_%d",m_uiCurrentStep,2);
            io::varToRawData((const ot::Mesh*)m_uiMesh,(const double **)m_uiStage[1],massgrav::MASSGRAV_NUM_VARS,NULL,frawName);
            #endif            
            
            // u_(k+1)
            for(unsigned int node=nodeLocalBegin; node<nodeLocalEnd; node++)
            {

               for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
               {
                   m_uiVar[index][node]=m_uiPrevVar[index][node] + m_uiT_h *( (1.0/6.0)*m_uiStage[0][index][node] + (1.0/6.0)*m_uiStage[1][index][node] + (2.0/3.0)*m_uiStage[2][index][node]);
               }
               enforce_massgrav_constraints(m_uiVar, node);
            }
         */   
            
        } else if (m_uiRKType==RKType::RK4)
        {   // rk4 solver
            //std::cout<<"rk4"<<std::endl;
            for(unsigned int stage=0; stage<(massgrav::MASSGRAV_RK4_STAGES-1); stage++)
            {


            #ifdef DEBUG_RK_SOLVER
                if(!rank)std::cout<<" stage: "<<stage<<" begin: "<<std::endl;
                for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
                    ot::test::isUnzipNaN(m_uiMesh,m_uiUnzipVar[index]);
            #endif

            #ifdef MASSGRAV_ENABLE_CUDA
                cuda::MASSGRAVComputeParams massgravParams;
                massgravParams.MASSGRAV_LAMBDA[0]=massgrav::MASSGRAV_LAMBDA[0];
                massgravParams.MASSGRAV_LAMBDA[1]=massgrav::MASSGRAV_LAMBDA[1];
                massgravParams.MASSGRAV_LAMBDA[2]=massgrav::MASSGRAV_LAMBDA[2];
                massgravParams.MASSGRAV_LAMBDA[3]=massgrav::MASSGRAV_LAMBDA[3];

                massgravParams.MASSGRAV_LAMBDA_F[0]=massgrav::MASSGRAV_LAMBDA_F[0];
                massgravParams.MASSGRAV_LAMBDA_F[1]=massgrav::MASSGRAV_LAMBDA_F[1];

                massgravParams.MASSGRAV_ETA_POWER[0]=massgrav::MASSGRAV_ETA_POWER[0];
                massgravParams.MASSGRAV_ETA_POWER[1]=massgrav::MASSGRAV_ETA_POWER[1];

                massgravParams.ETA_R0=massgrav::ETA_R0;
                massgravParams.ETA_CONST=massgrav::ETA_CONST;
                massgravParams.ETA_DAMPING=massgrav::ETA_DAMPING;
                massgravParams.ETA_DAMPING_EXP=massgrav::ETA_DAMPING_EXP;
                massgravParams.KO_DISS_SIGMA=massgrav::KO_DISS_SIGMA;

                dim3 threadBlock(16,16,1);
                cuda::computeRHS(m_uiUnzipVarRHS,(const double **)m_uiUnzipVar,&(*(blkList.begin())),blkList.size(),(const cuda::MASSGRAVComputeParams*) &massgravParams,threadBlock,pt_min,pt_max,1);
            #else

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


            #ifdef MASSGRAV_RHS_STAGED_COMP

                    massgravrhs_sep(m_uiUnzipVarRHS, (const double **)m_uiUnzipVar, offset, ptmin, ptmax, sz, bflag);
            #else
                    massgravrhs(m_uiUnzipVarRHS, (const double **)m_uiUnzipVar, offset, ptmin, ptmax, sz, bflag);
            #endif


                }
            #endif

            #ifdef DEBUG_RK_SOLVER
                if(!rank)std::cout<<" stage: "<<stage<<" af rhs UNZIP RHS TEST:"<<std::endl;
                for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
                    ot::test::isUnzipInternalNaN(m_uiMesh,m_uiUnzipVarRHS[index]);
            #endif



                zipVars(m_uiUnzipVarRHS,m_uiStage[stage]);


            #ifdef DEBUG_RK_SOLVER
                for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
                    if(seq::test::isNAN(m_uiStage[stage][index]+m_uiMesh->getNodeLocalBegin(),m_uiMesh->getNumLocalMeshNodes()))
                        std::cout<<" var: "<<index<<" contains nan af zip  stage: "<<stage<<std::endl;
            #endif

                /*for(unsigned int index=0;index<massgrav::MASSGRAV_NUM_VARS;index++)
                     for(unsigned int node=nodeLocalBegin;node<nodeLocalEnd;node++)
                             m_uiStage[stage][index][node]*=m_uiT_h;*/



                for(unsigned int node=nodeLocalBegin; node<nodeLocalEnd; node++)
                {
                    for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
                    {
                        m_uiVarIm[index][node]=m_uiPrevVar[index][node];
                        //if(!rank && index==0 && node==0) std::cout<<"rk stage: "<<stage<<" im coef: "<<stage<<" value: "<<RK4_U[stage+1]<<std::endl;
                        m_uiVarIm[index][node] += (RK4_U[stage + 1] * m_uiT_h * m_uiStage[stage][index][node]);

                    }
                    enforce_massgrav_constraints(m_uiVarIm, node);
                }


                current_t_adv=current_t+RK4_T[stage+1]*m_uiT_h;
            #ifdef RK_SOLVER_OVERLAP_COMM_AND_COMP
                unzipVars_async(m_uiVarIm,m_uiUnzipVar);
            #else
                performGhostExchangeVars(m_uiVarIm);
                unzipVars(m_uiVarIm,m_uiUnzipVar);
            #endif


            }

            current_t_adv=current_t+RK4_T[(massgrav::MASSGRAV_RK4_STAGES-1)]*m_uiT_h;


            #ifdef DEBUG_RK_SOLVER
            if(!rank)std::cout<<" stage: "<<(massgrav::MASSGRAV_RK4_STAGES-1)<<" begin: "<<std::endl;

            for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
                ot::test::isUnzipNaN(m_uiMesh,m_uiUnzipVar[index]);
            #endif

            #ifdef MASSGRAV_ENABLE_CUDA
            cuda::MASSGRAVComputeParams massgravParams;
            massgravParams.MASSGRAV_LAMBDA[0]=massgrav::MASSGRAV_LAMBDA[0];
            massgravParams.MASSGRAV_LAMBDA[1]=massgrav::MASSGRAV_LAMBDA[1];
            massgravParams.MASSGRAV_LAMBDA[2]=massgrav::MASSGRAV_LAMBDA[2];
            massgravParams.MASSGRAV_LAMBDA[3]=massgrav::MASSGRAV_LAMBDA[3];

            massgravParams.MASSGRAV_LAMBDA_F[0]=massgrav::MASSGRAV_LAMBDA_F[0];
            massgravParams.MASSGRAV_LAMBDA_F[1]=massgrav::MASSGRAV_LAMBDA_F[1];

            massgravParams.MASSGRAV_ETA_POWER[0]=massgrav::MASSGRAV_ETA_POWER[0];
            massgravParams.MASSGRAV_ETA_POWER[1]=massgrav::MASSGRAV_ETA_POWER[1];

            massgravParams.ETA_R0=massgrav::ETA_R0;
            massgravParams.ETA_CONST=massgrav::ETA_CONST;
            massgravParams.ETA_DAMPING=massgrav::ETA_DAMPING;
            massgravParams.ETA_DAMPING_EXP=massgrav::ETA_DAMPING_EXP;
            massgravParams.KO_DISS_SIGMA=massgrav::KO_DISS_SIGMA;

            dim3 threadBlock(16,16,1);
            cuda::computeRHS(m_uiUnzipVarRHS,(const double **)m_uiUnzipVar,&(*(blkList.begin())),blkList.size(),(const cuda::MASSGRAVComputeParams*) &massgravParams,threadBlock,pt_min,pt_max,1);
            #else


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


            #ifdef MASSGRAV_RHS_STAGED_COMP

                massgravrhs_sep(m_uiUnzipVarRHS, (const double **)m_uiUnzipVar, offset, ptmin, ptmax, sz, bflag);
            #else
                massgravrhs(m_uiUnzipVarRHS, (const double **)m_uiUnzipVar, offset, ptmin, ptmax, sz, bflag);
            #endif


            }

        #endif

        #ifdef DEBUG_RK_SOLVER
            if(!rank)std::cout<<" stage: "<<(massgrav::MASSGRAV_RK4_STAGES-1)<<" af rhs UNZIP RHS TEST:"<<std::endl;
            for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
                ot::test::isUnzipInternalNaN(m_uiMesh,m_uiUnzipVarRHS[index]);
        #endif

            zipVars(m_uiUnzipVarRHS,m_uiStage[(massgrav::MASSGRAV_RK4_STAGES-1)]);

            for(unsigned int node=nodeLocalBegin; node<nodeLocalEnd; node++)
            {

                for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
                {
                    m_uiVar[index][node]=m_uiPrevVar[index][node];
                    for(unsigned int s=0; s<(massgrav::MASSGRAV_RK4_STAGES); s++)
                    {
                        m_uiVar[index][node]+=(RK4_C[s]*m_uiT_h*m_uiStage[s][index][node]);
                    }


                }
                enforce_massgrav_constraints(m_uiVar, node);

            }

        } else if (m_uiRKType==RKType::RK45)
        {   // rk45 solver

            //std::cout<<"rk45"<<std::endl;

            bool repeatStep;
            double n_inf_max=0.0;
            double n_inf_max_g=0;

            do {

                repeatStep=false;
                n_inf_max=0;
                n_inf_max_g=0;

                if(m_uiMesh->isActive())
                {
                    double current_t=m_uiCurrentTime;
                    double current_t_adv=current_t;
                    #ifdef RK_SOLVER_OVERLAP_COMM_AND_COMP
                    unzipVars_async(m_uiPrevVar,m_uiUnzipVar);
                    #else
                    //1. perform ghost exchange.
                    performGhostExchangeVars(m_uiPrevVar);

                    //2. unzip all the variables.
                    unzipVars(m_uiPrevVar,m_uiUnzipVar);
                    #endif


                    int rank =m_uiMesh->getMPIRank();



                    const unsigned int nodeLocalBegin=m_uiMesh->getNodeLocalBegin();
                    const unsigned int nodeLocalEnd=m_uiMesh->getNodeLocalEnd();

                    const std::vector<ot::Block> blkList=m_uiMesh->getLocalBlockList();

                    unsigned int offset;
                    double ptmin[3], ptmax[3];
                    unsigned int sz[3];
                    unsigned int bflag;
                    double dx,dy,dz;
                    const Point pt_min(massgrav::MASSGRAV_COMPD_MIN[0],massgrav::MASSGRAV_COMPD_MIN[1],massgrav::MASSGRAV_COMPD_MIN[2]);
                    const Point pt_max(massgrav::MASSGRAV_COMPD_MAX[0],massgrav::MASSGRAV_COMPD_MAX[1],massgrav::MASSGRAV_COMPD_MAX[2]);


                    for(unsigned int stage=0; stage<(massgrav::MASSGRAV_RK45_STAGES-1); stage++)
                    {


                        #ifdef DEBUG_RK_SOLVER
                        if(!rank)std::cout<<" stage: "<<stage<<" begin: "<<std::endl;
                        for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
                            ot::test::isUnzipNaN(m_uiMesh,m_uiUnzipVar[index]);
                        #endif


                        #ifdef MASSGRAV_ENABLE_CUDA
                        cuda::MASSGRAVComputeParams massgravParams;
                        massgravParams.MASSGRAV_LAMBDA[0]=massgrav::MASSGRAV_LAMBDA[0];
                        massgravParams.MASSGRAV_LAMBDA[1]=massgrav::MASSGRAV_LAMBDA[1];
                        massgravParams.MASSGRAV_LAMBDA[2]=massgrav::MASSGRAV_LAMBDA[2];
                        massgravParams.MASSGRAV_LAMBDA[3]=massgrav::MASSGRAV_LAMBDA[3];

                        massgravParams.MASSGRAV_LAMBDA_F[0]=massgrav::MASSGRAV_LAMBDA_F[0];
                        massgravParams.MASSGRAV_LAMBDA_F[1]=massgrav::MASSGRAV_LAMBDA_F[1];

                        massgravParams.MASSGRAV_ETA_POWER[0]=massgrav::MASSGRAV_ETA_POWER[0];
                        massgravParams.MASSGRAV_ETA_POWER[1]=massgrav::MASSGRAV_ETA_POWER[1];

                        massgravParams.ETA_R0=massgrav::ETA_R0;
                        massgravParams.ETA_CONST=massgrav::ETA_CONST;
                        massgravParams.ETA_DAMPING=massgrav::ETA_DAMPING;
                        massgravParams.ETA_DAMPING_EXP=massgrav::ETA_DAMPING_EXP;
                        massgravParams.KO_DISS_SIGMA=massgrav::KO_DISS_SIGMA;

                        dim3 threadBlock(16,16,1);
                        cuda::computeRHS(m_uiUnzipVarRHS,(const double **)m_uiUnzipVar,&(*(blkList.begin())),blkList.size(),(const cuda::MASSGRAVComputeParams*) &massgravParams,threadBlock,pt_min,pt_max,1);
                        #else

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


                            #ifdef MASSGRAV_RHS_STAGED_COMP

                            massgravrhs_sep(m_uiUnzipVarRHS, (const double **)m_uiUnzipVar, offset, ptmin, ptmax, sz, bflag);
                            #else
                            massgravrhs(m_uiUnzipVarRHS, (const double **)m_uiUnzipVar, offset, ptmin, ptmax, sz, bflag);
                            #endif


                        }
                        #endif


                        #ifdef DEBUG_RK_SOLVER
                        if(!rank)std::cout<<" stage: "<<stage<<" af rhs UNZIP RHS TEST:"<<std::endl;
                        for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
                            ot::test::isUnzipInternalNaN(m_uiMesh,m_uiUnzipVarRHS[index]);
                        #endif



                        zipVars(m_uiUnzipVarRHS,m_uiStage[stage]);


                    #ifdef DEBUG_RK_SOLVER
                        for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
                            if(seq::test::isNAN(m_uiStage[stage][index]+m_uiMesh->getNodeLocalBegin(),m_uiMesh->getNumLocalMeshNodes()))
                                std::cout<<" var: "<<index<<" contains nan af zip  stage: "<<stage<<std::endl;
                    #endif

                        /*for(unsigned int index=0;index<massgrav::MASSGRAV_NUM_VARS;index++)
                            for(unsigned int node=nodeLocalBegin;node<nodeLocalEnd;node++)
                                    m_uiStage[stage][index][node]*=m_uiT_h;*/



                        for(unsigned int node=nodeLocalBegin; node<nodeLocalEnd; node++)
                        {
                            for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
                            {
                                m_uiVarIm[index][node]=m_uiPrevVar[index][node];
                                for(unsigned int s=0; s<(stage+1); s++) {
                                    //if(!rank && index==0 && node==0) std::cout<<"rk stage: "<<stage<<" im coef: "<<s<<" value: "<<RK_U[stage+1][s]<<std::endl;
                                    m_uiVarIm[index][node] += (RK_U[stage + 1][s] * m_uiT_h * m_uiStage[s][index][node]);
                                }
                            }
                            enforce_massgrav_constraints(m_uiVarIm, node);
                        }


                        current_t_adv=current_t+RK_T[stage+1]*m_uiT_h;
                    #ifdef RK_SOLVER_OVERLAP_COMM_AND_COMP
                        unzipVars_async(m_uiVarIm,m_uiUnzipVar);
                    #else
                        performGhostExchangeVars(m_uiVarIm);
                        unzipVars(m_uiVarIm,m_uiUnzipVar);
                    #endif


                    }

                    current_t_adv=current_t+RK_T[(massgrav::MASSGRAV_RK45_STAGES-1)];


                    #ifdef DEBUG_RK_SOLVER
                    if(!rank)std::cout<<" stage: "<<(massgrav::MASSGRAV_RK45_STAGES-1)<<" begin: "<<std::endl;

                    for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
                        ot::test::isUnzipNaN(m_uiMesh,m_uiUnzipVar[index]);
                    #endif



                    #ifdef MASSGRAV_ENABLE_CUDA
                    cuda::MASSGRAVComputeParams massgravParams;
                    massgravParams.MASSGRAV_LAMBDA[0]=massgrav::MASSGRAV_LAMBDA[0];
                    massgravParams.MASSGRAV_LAMBDA[1]=massgrav::MASSGRAV_LAMBDA[1];
                    massgravParams.MASSGRAV_LAMBDA[2]=massgrav::MASSGRAV_LAMBDA[2];
                    massgravParams.MASSGRAV_LAMBDA[3]=massgrav::MASSGRAV_LAMBDA[3];

                    massgravParams.MASSGRAV_LAMBDA_F[0]=massgrav::MASSGRAV_LAMBDA_F[0];
                    massgravParams.MASSGRAV_LAMBDA_F[1]=massgrav::MASSGRAV_LAMBDA_F[1];

                    massgravParams.MASSGRAV_ETA_POWER[0]=massgrav::MASSGRAV_ETA_POWER[0];
                    massgravParams.MASSGRAV_ETA_POWER[1]=massgrav::MASSGRAV_ETA_POWER[1];

                    massgravParams.ETA_R0=massgrav::ETA_R0;
                    massgravParams.ETA_CONST=massgrav::ETA_CONST;
                    massgravParams.ETA_DAMPING=massgrav::ETA_DAMPING;
                    massgravParams.ETA_DAMPING_EXP=massgrav::ETA_DAMPING_EXP;
                    massgravParams.KO_DISS_SIGMA=massgrav::KO_DISS_SIGMA;

                    dim3 threadBlock(16,16,1);
                    cuda::computeRHS(m_uiUnzipVarRHS,(const double **)m_uiUnzipVar,&(*(blkList.begin())),blkList.size(),(const cuda::MASSGRAVComputeParams*) &massgravParams,threadBlock,pt_min,pt_max,1);
                #else

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


                        #ifdef MASSGRAV_RHS_STAGED_COMP

                        massgravrhs_sep(m_uiUnzipVarRHS, (const double **)m_uiUnzipVar, offset, ptmin, ptmax, sz, bflag);
                        #else
                        massgravrhs(m_uiUnzipVarRHS, (const double **)m_uiUnzipVar, offset, ptmin, ptmax, sz, bflag);
                        #endif


                    }

                #endif

                #ifdef DEBUG_RK_SOLVER
                    if(!rank)std::cout<<" stage: "<<(massgrav::MASSGRAV_RK45_STAGES-1)<<" af rhs UNZIP RHS TEST:"<<std::endl;

                    for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
                        ot::test::isUnzipInternalNaN(m_uiMesh,m_uiUnzipVarRHS[index]);
                #endif

                    zipVars(m_uiUnzipVarRHS,m_uiStage[(massgrav::MASSGRAV_RK45_STAGES-1)]);

                    /*for(unsigned int index=0;index<massgrav::MASSGRAV_NUM_VARS;index++)
                        for(unsigned int node=nodeLocalBegin;node<nodeLocalEnd;node++)
                            m_uiStage[(massgrav::MASSGRAV_RK45_STAGES-1)][index][node]*=m_uiT_h;*/



                    for(unsigned int node=nodeLocalBegin; node<nodeLocalEnd; node++)
                    {

                        for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
                        {
                            m_uiVarIm[index][node]=m_uiPrevVar[index][node];
                            for(unsigned int s=0; s<(massgrav::MASSGRAV_RK45_STAGES-1); s++)
                            {
                                if(s==1) continue;
                                m_uiVarIm[index][node]+=(RK_4_C[s]*m_uiT_h*m_uiStage[s][index][node]);
                            }


                            m_uiVar[index][node]=m_uiPrevVar[index][node];
                            for(unsigned int s=0; s<massgrav::MASSGRAV_RK45_STAGES; s++)
                            {   if(s==1) continue; // because rk coef is zero.
                                m_uiVar[index][node]+=(RK_5_C[s]*m_uiT_h*m_uiStage[s][index][node]);
                            }


                        }

                        enforce_massgrav_constraints(m_uiVarIm, node);
                        enforce_massgrav_constraints(m_uiVar, node);

                    }


                    // update the m_uiTh bases on the normed diff between m_uiVarIm, m_uiVar.

                    double n_inf;
                    n_inf_max=0;
                    for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
                    {
                        n_inf=normLInfty(m_uiVarIm[index]+nodeLocalBegin,m_uiVar[index]+nodeLocalBegin,(nodeLocalEnd-nodeLocalBegin),m_uiMesh->getMPICommunicator());
                        if(n_inf>n_inf_max)
                            n_inf_max=n_inf;
                    }

                }

                // below all reduction is act as an barrier for inactive procs.
                par::Mpi_Allreduce(&n_inf_max,&n_inf_max_g,1,MPI_MAX,m_uiComm);
                n_inf_max=n_inf_max_g;

                if(n_inf_max>massgrav::MASSGRAV_RK45_DESIRED_TOL)
                {
                    repeatStep=true;
                    m_uiT_h=massgrav::MASSGRAV_SAFETY_FAC*m_uiT_h*(pow(fabs(massgrav::MASSGRAV_RK45_DESIRED_TOL/n_inf_max),0.25));
                    if(!m_uiMesh->getMPIRankGlobal()) std::cout<<" repeat : "<<m_uiCurrentStep<<" with : "<<m_uiT_h<<std::endl;

                } else
                {
                    repeatStep=false;
                    m_uiT_h=massgrav::MASSGRAV_SAFETY_FAC*m_uiT_h*(pow(fabs(massgrav::MASSGRAV_RK45_DESIRED_TOL/n_inf_max),0.20));
                }



            } while(repeatStep);




        }






    }

    m_uiMesh->waitAll();

    m_uiCurrentStep++;
    m_uiCurrentTime+=m_uiT_h;


}


void RK_MASSGRAV::rkSolve()
{

    if(m_uiCurrentStep==0)
    {
        //applyInitialConditions(m_uiPrevVar);
        initialGridConverge();
    }

    bool isRefine=true;
    unsigned int oldElements,oldElements_g;
    unsigned int newElements,newElements_g;

    // refine based on all the variables
    const unsigned int refineNumVars=massgrav::MASSGRAV_NUM_REFINE_VARS;
    unsigned int refineVarIds[refineNumVars];
    for(unsigned int vIndex=0; vIndex<refineNumVars; vIndex++)
        refineVarIds[vIndex]=massgrav::MASSGRAV_REFINE_VARIABLE_INDICES[vIndex];

    double wTol=massgrav::MASSGRAV_WAVELET_TOL;
    std::function<double(double,double,double)> waveletTolFunc =[wTol](double x,double y, double z) {
        return massgrav::computeWTol(x,y,z,wTol);
    };
    Point bhLoc[2];

    double l_min,l_max;
    for(double t=m_uiCurrentTime; t<m_uiTimeEnd; t=t+m_uiT_h)
    {

        // checkpoint the previous solution value before going to the next step.
        massgrav::timer::t_ioCheckPoint.start();
        if((m_uiMesh->isActive()) && (m_uiCurrentStep%massgrav::MASSGRAV_CHECKPT_FREQ)==0)
            storeCheckPoint(massgrav::MASSGRAV_CHKPT_FILE_PREFIX.c_str());
        massgrav::timer::t_ioCheckPoint.stop();

        // write sol to vtu.
        if((m_uiMesh->isActive()) && (m_uiCurrentStep%massgrav::MASSGRAV_TIME_STEP_OUTPUT_FREQ)==0)
        {

            l_min=vecMin(m_uiPrevVar[massgrav::VAR::U_ALPHA]+m_uiMesh->getNodeLocalBegin(),(m_uiMesh->getNumLocalMeshNodes()),m_uiMesh->getMPICommunicator());
            l_max=vecMax(m_uiPrevVar[massgrav::VAR::U_ALPHA]+m_uiMesh->getNodeLocalBegin(),(m_uiMesh->getNumLocalMeshNodes()),m_uiMesh->getMPICommunicator());
            if(!m_uiMesh->getMPIRank()) {
                std::cout << "executing step: " << m_uiCurrentStep << " dt: " << m_uiT_h << " rk_time : "<< m_uiCurrentTime << std::endl;
                std::cout << "\t ||VAR::U_ALPHA|| (min, max) : ("<<l_min<<", "<<l_max<<" ) "<<std::endl;
            }

            massgrav::timer::profileInfoIntermediate(massgrav::MASSGRAV_PROFILE_FILE_PREFIX.c_str(),m_uiMesh,m_uiCurrentStep);

        }

        if((m_uiCurrentStep%massgrav::MASSGRAV_TIME_STEP_OUTPUT_FREQ)==0)
            massgrav::timer::resetSnapshot();


        if((m_uiCurrentStep%massgrav::MASSGRAV_REMESH_TEST_FREQ)==0)
        {

            #ifdef RK_SOLVER_OVERLAP_COMM_AND_COMP
            unzipVars_async(m_uiPrevVar,m_uiUnzipVar);
            #else
            performGhostExchangeVars(m_uiPrevVar);
            //isRefine=m_uiMesh->isReMesh((const double **)m_uiPrevVar,refineVarIds,refineNumVars,massgrav::MASSGRAV_WAVELET_TOL);
            unzipVars(m_uiPrevVar,m_uiUnzipVar);
            #endif

            //                     char fPrefix[256];
            //                     sprintf(fPrefix,"%s_wavelets_%d",massgrav::MASSGRAV_VTU_FILE_PREFIX.c_str(),m_uiCurrentStep);
            //                     io::vtk::waveletsToVTU(m_uiMesh,(const double **)m_uiPrevVar,(const double**)m_uiUnzipVar,refineVarIds,refineNumVars,(const char*)fPrefix);

        #ifdef DEBUG_RK_SOLVER
            if(m_uiMesh->isActive())
            {
                if(!m_uiMesh->getMPIRank())std::cout<<" isRemesh Unzip : "<<std::endl;
                for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
                    ot::test::isUnzipNaN(m_uiMesh,m_uiUnzipVar[index]);

            }

        #endif
            massgrav::timer::t_isReMesh.start();
            //const double r[3] ={3.0,3.0,3.0};
            if(massgrav::MASSGRAV_ENABLE_BLOCK_ADAPTIVITY)
                isRefine=false;
            else
            {
                if(massgrav::MASSGRAV_REFINEMENT_MODE == massgrav::RefinementMode::WAMR)
                {
                    isRefine=m_uiMesh->isReMeshUnzip((const double **)m_uiUnzipVar,refineVarIds,refineNumVars,waveletTolFunc,massgrav::MASSGRAV_DENDRO_AMR_FAC); 

                }else if(massgrav::MASSGRAV_REFINEMENT_MODE == massgrav::RefinementMode::EH)
                {
                    isRefine = massgrav::isRemeshEH(m_uiMesh,(const double **)m_uiUnzipVar,massgrav::VAR::U_ALPHA,massgrav::MASSGRAV_EH_REFINE_VAL,massgrav::MASSGRAV_EH_COARSEN_VAL,true);

                }else if(massgrav::MASSGRAV_REFINEMENT_MODE == massgrav::RefinementMode::EH_WAMR)
                {
                    const bool isR1 = m_uiMesh->isReMeshUnzip((const double **)m_uiUnzipVar,refineVarIds,refineNumVars,waveletTolFunc,massgrav::MASSGRAV_DENDRO_AMR_FAC); 
                    const bool isR2 = massgrav::isRemeshEH(m_uiMesh,(const double **)m_uiUnzipVar,massgrav::VAR::U_ALPHA,massgrav::MASSGRAV_EH_REFINE_VAL,massgrav::MASSGRAV_EH_COARSEN_VAL,false);

                    isRefine = (isR1 || isR2);
                    
                }else
                {
                    std::cout<<" Error : "<<__func__<<" invalid refinement mode specified "<<std::endl;
                    MPI_Abort(m_uiComm,0);
                }
                
                
            }
                
                
            massgrav::timer::t_isReMesh.stop();

            if(isRefine)
            {


            #ifdef DEBUG_IS_REMESH
                unsigned int rank=m_uiMesh->getMPIRankGlobal();
                MPI_Comm globalComm=m_uiMesh->getMPIGlobalCommunicator();
                std::vector<ot::TreeNode> unChanged;
                std::vector<ot::TreeNode> refined;
                std::vector<ot::TreeNode> coarsened;
                std::vector<ot::TreeNode> localBlocks;

                const ot::Block* blkList=&(*(m_uiMesh->getLocalBlockList().begin()));
                for(unsigned int ele=0; ele<m_uiMesh->getLocalBlockList().size(); ele++)
                {
                    localBlocks.push_back(blkList[ele].getBlockNode());
                }


                const ot::TreeNode * pNodes=&(*(m_uiMesh->getAllElements().begin()));
                for(unsigned int ele=m_uiMesh->getElementLocalBegin(); ele<m_uiMesh->getElementLocalEnd(); ele++)
                {
                    if((pNodes[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_NO_CHANGE)
                    {
                        unChanged.push_back(pNodes[ele]);
                    } else if((pNodes[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT)
                    {
                        refined.push_back(pNodes[ele]);
                    } else
                    {
                        assert((pNodes[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_COARSE);
                        coarsened.push_back(pNodes[ele]);
                    }
                }

                char fN1[256];
                char fN2[256];
                char fN3[256];
                char fN4[256];

                sprintf(fN1,"unchanged_%d",m_uiCurrentStep);
                sprintf(fN2,"refined_%d",m_uiCurrentStep);
                sprintf(fN3,"coarsend_%d",m_uiCurrentStep);
                sprintf(fN4,"blocks_%d",m_uiCurrentStep);

                DendroIntL localSz=unChanged.size();
                DendroIntL globalSz;
                par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,globalComm);
                if(!rank) std::cout<<" total unchanged: "<<globalSz<<std::endl;

                localSz=refined.size();
                par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,globalComm);
                if(!rank) std::cout<<" total refined: "<<globalSz<<std::endl;


                localSz=coarsened.size();
                par::Mpi_Reduce(&localSz,&globalSz,1,MPI_SUM,0,globalComm);
                if(!rank) std::cout<<" total coarsend: "<<globalSz<<std::endl;


                io::vtk::oct2vtu(&(*(unChanged.begin())),unChanged.size(),fN1,globalComm);
                io::vtk::oct2vtu(&(*(refined.begin())),refined.size(),fN2,globalComm);
                io::vtk::oct2vtu(&(*(coarsened.begin())),coarsened.size(),fN3,globalComm);
                io::vtk::oct2vtu(&(*(localBlocks.begin())),localBlocks.size(),fN4,globalComm);

            #endif
                massgrav::timer::t_mesh.start();
                ot::Mesh* newMesh=m_uiMesh->ReMesh(massgrav::MASSGRAV_DENDRO_GRAIN_SZ,massgrav::MASSGRAV_LOAD_IMB_TOL,massgrav::MASSGRAV_SPLIT_FIX);
                massgrav::timer::t_mesh.stop();

                oldElements=m_uiMesh->getNumLocalMeshElements();
                newElements=newMesh->getNumLocalMeshElements();

                par::Mpi_Reduce(&oldElements,&oldElements_g,1,MPI_SUM,0,m_uiMesh->getMPIGlobalCommunicator());
                par::Mpi_Reduce(&newElements,&newElements_g,1,MPI_SUM,0,newMesh->getMPIGlobalCommunicator());

                if(!(m_uiMesh->getMPIRankGlobal()))std::cout<<"step : "<<m_uiCurrentStep<<" time : "<<m_uiCurrentTime<<" old mesh: "<<oldElements_g<<" new mesh: "<<newElements_g<<std::endl;


                // performs the inter-grid transfer
                intergridTransferVars(m_uiPrevVar,newMesh,massgrav::MASSGRAV_USE_FD_GRID_TRANSFER);

                for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
                {
                    delete [] m_uiVar[index];
                    delete [] m_uiVarIm[index];
                    delete [] m_uiUnzipVar[index];
                    delete [] m_uiUnzipVarRHS[index];

                    m_uiVar[index]=NULL;
                    m_uiVarIm[index]=NULL;
                    m_uiUnzipVar[index]=NULL;
                    m_uiUnzipVarRHS[index]=NULL;

                    m_uiVar[index]=newMesh->createVector<DendroScalar>();
                    m_uiVarIm[index]=newMesh->createVector<DendroScalar>();
                    m_uiUnzipVar[index]=newMesh->createUnZippedVector<DendroScalar>();
                    m_uiUnzipVarRHS[index]=newMesh->createUnZippedVector<DendroScalar>();


                }

                for(unsigned int stage=0; stage<m_uiNumRKStages; stage++)
                    for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
                    {
                        delete [] m_uiStage[stage][index];
                        m_uiStage[stage][index]=NULL;
                        m_uiStage[stage][index]=newMesh->createVector<DendroScalar>();
                    }

                // deallocate constraint vars allocate them for the new mesh.
                for(unsigned int index=0; index<massgrav::MASSGRAV_CONSTRAINT_NUM_VARS; index++)
                {
                    delete [] m_uiConstraintVars[index];
                    delete [] m_uiUnzipConstraintVars[index];

                    m_uiConstraintVars[index]=newMesh->createVector<DendroScalar>();
                    m_uiUnzipConstraintVars[index]=newMesh->createUnZippedVector<DendroScalar>();

                }




                std::swap(newMesh,m_uiMesh);
                delete newMesh;

                if(m_uiCurrentStep == 0)
                    applyInitialConditions(m_uiPrevVar);


                unsigned int lmin, lmax;
                m_uiMesh->computeMinMaxLevel(lmin,lmax);    
                massgrav::MASSGRAV_RK45_TIME_STEP_SIZE=massgrav::MASSGRAV_CFL_FACTOR*((massgrav::MASSGRAV_COMPD_MAX[0]-massgrav::MASSGRAV_COMPD_MIN[0])*((1u<<(m_uiMaxDepth-lmax))/((double) massgrav::MASSGRAV_ELE_ORDER))/((double)(1u<<(m_uiMaxDepth)))); 

                #ifdef RK_SOLVER_OVERLAP_COMM_AND_COMP
                    // reallocates mpi resources for the the new mesh. (this will deallocate the old resources)
                    reallocateMPIResources();
                #endif

                if(m_uiMesh->isActive())
                {
                    l_min=vecMin(m_uiPrevVar[massgrav::VAR::U_ALPHA]+m_uiMesh->getNodeLocalBegin(),(m_uiMesh->getNumLocalMeshNodes()),m_uiMesh->getMPICommunicator());
                    l_max=vecMax(m_uiPrevVar[massgrav::VAR::U_ALPHA]+m_uiMesh->getNodeLocalBegin(),(m_uiMesh->getNumLocalMeshNodes()),m_uiMesh->getMPICommunicator());
                    if(!(m_uiMesh->getMPIRank())) {
                        std::cout << "transfer completed:    ||VAR::U_ALPHA|| (min, max) : ("<<l_min<<", "<<l_max<<" ) "<<std::endl;
                    }

                }

                


            }

        }




        if((m_uiCurrentStep % massgrav::MASSGRAV_GW_EXTRACT_FREQ)==0)
        {

        #ifdef RK_SOLVER_OVERLAP_COMM_AND_COMP
            unzipVars_async(m_uiPrevVar,m_uiUnzipVar);
        #else
            performGhostExchangeVars(m_uiPrevVar);
            unzipVars(m_uiPrevVar,m_uiUnzipVar);
        #endif


        #ifdef MASSGRAV_COMPUTE_CONSTRAINTS

            const std::vector<ot::Block> blkList=m_uiMesh->getLocalBlockList();

            unsigned int offset;
            double ptmin[3], ptmax[3];
            unsigned int sz[3];
            unsigned int bflag;
            double dx,dy,dz;
            const Point pt_min(massgrav::MASSGRAV_COMPD_MIN[0],massgrav::MASSGRAV_COMPD_MIN[1],massgrav::MASSGRAV_COMPD_MIN[2]);
            const Point pt_max(massgrav::MASSGRAV_COMPD_MAX[0],massgrav::MASSGRAV_COMPD_MAX[1],massgrav::MASSGRAV_COMPD_MAX[2]);

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

                physical_constraints(m_uiUnzipConstraintVars, (const DendroScalar **) m_uiUnzipVar, offset, ptmin, ptmax, sz, bflag);
            }

            /*double consVecMin[massgrav::MASSGRAV_CONSTRAINT_NUM_VARS];
            double consVecMax[massgrav::MASSGRAV_CONSTRAINT_NUM_VARS];*/
            double constraintMaskedL2[massgrav::MASSGRAV_CONSTRAINT_NUM_VARS];
            for(unsigned int index=0; index<massgrav::MASSGRAV_CONSTRAINT_NUM_VARS; index++)
            {
                m_uiMesh->zip(m_uiUnzipConstraintVars[index],m_uiConstraintVars[index]);
                m_uiMesh->performGhostExchange(m_uiConstraintVars[index]);
            }

            massgrav::extractConstraints(m_uiMesh,(const DendroScalar **)m_uiConstraintVars,m_uiPrevVar[BHLOC::EXTRACTION_VAR_ID],BHLOC::EXTRACTION_TOL,m_uiCurrentStep);
        #ifndef MASSGRAV_KERR_SCHILD_TEST
            #ifdef MASSGRAV_EXTRACT_GRAVITATIONAL_WAVES
            GW::extractFarFieldPsi4(m_uiMesh,(const DendroScalar **)m_uiConstraintVars,m_uiCurrentStep,m_uiCurrentTime);
            #endif
        #endif

    #endif

        #ifdef MASSGRAV_ENABLE_VTU_OUTPUT
            if((m_uiCurrentStep % massgrav::MASSGRAV_IO_OUTPUT_FREQ) ==0)
                writeToVTU(m_uiPrevVar,m_uiConstraintVars,massgrav::MASSGRAV_NUM_EVOL_VARS_VTU_OUTPUT,massgrav::MASSGRAV_NUM_CONST_VARS_VTU_OUTPUT,massgrav::MASSGRAV_VTU_OUTPUT_EVOL_INDICES,massgrav::MASSGRAV_VTU_OUTPUT_CONST_INDICES,massgrav::MASSGRAV_VTU_Z_SLICE_ONLY);
        #endif


        #ifdef MASSGRAV_EXTRACT_BH_LOCATIONS
            massgrav::writeBHCoordinates((const ot::Mesh *)m_uiMesh,(const Point *) m_uiBHLoc,2,m_uiCurrentStep);
        #endif




        }



        massgrav::timer::t_rkStep.start();

        performSingleIteration();

        massgrav::timer::t_rkStep.stop();

        #ifdef MASSGRAV_EXTRACT_BH_LOCATIONS
            massgrav::extractBHCoords((const ot::Mesh *)m_uiMesh,(const DendroScalar*)m_uiVar[BHLOC::EXTRACTION_VAR_ID],BHLOC::EXTRACTION_TOL,(const Point *) m_uiBHLoc,2,(Point*)bhLoc);
            m_uiBHLoc[0]=bhLoc[0];
            m_uiBHLoc[1]=bhLoc[1];
        #endif


        std::swap(m_uiVar,m_uiPrevVar);
        //massgrav::artificial_dissipation(m_uiMesh,m_uiPrevVar,massgrav::MASSGRAV_NUM_VARS,massgrav::MASSGRAV_DISSIPATION_NC,massgrav::MASSGRAV_DISSIPATION_S,false);
        //if(m_uiCurrentStep==1) break;

    }


}


void RK_MASSGRAV::storeCheckPoint(const char * fNamePrefix)
{


    if(m_uiMesh->isActive())
    {
        unsigned int cpIndex;
        (m_uiCurrentStep%(2*massgrav::MASSGRAV_CHECKPT_FREQ)==0) ? cpIndex=0 : cpIndex=1; // to support alternate file writing.
        unsigned int rank=m_uiMesh->getMPIRank();
        unsigned int npes=m_uiMesh->getMPICommSize();

        char fName[256];
        const ot::TreeNode * pNodes=&(*(m_uiMesh->getAllElements().begin()+m_uiMesh->getElementLocalBegin()));
        sprintf(fName,"%s_octree_%d_%d.oct",fNamePrefix,cpIndex,rank);
        io::checkpoint::writeOctToFile(fName,pNodes,m_uiMesh->getNumLocalMeshElements());

        unsigned int numVars=massgrav::MASSGRAV_NUM_VARS;
        const char ** varNames=massgrav::MASSGRAV_VAR_NAMES;

        /*for(unsigned int i=0;i<numVars;i++)
        {
            sprintf(fName,"%s_%s_%d_%d.var",fNamePrefix,varNames[i],cpIndex,rank);
            io::checkpoint::writeVecToFile(fName,m_uiMesh,m_uiPrevVar[i]);
        }*/

        sprintf(fName,"%s_%d_%d.var",fNamePrefix,cpIndex,rank);
        io::checkpoint::writeVecToFile(fName,m_uiMesh,(const double **)m_uiPrevVar,massgrav::MASSGRAV_NUM_VARS);

        if(!rank)
        {
            sprintf(fName,"%s_step_%d.cp",fNamePrefix,cpIndex);
            std::cout<<"writing : "<<fName<<std::endl;
            std::ofstream outfile(fName);
            if(!outfile) {
                std::cout<<fName<<" file open failed "<<std::endl;
                return ;
            }

            json checkPoint;
            checkPoint["DENDRO_RK45_TIME_BEGIN"]=m_uiTimeBegin;
            checkPoint["DENDRO_RK45_TIME_END"]=m_uiTimeEnd;
            checkPoint["DENDRO_RK45_ELEMENT_ORDER"]=m_uiOrder;

            checkPoint["DENDRO_RK45_TIME_CURRENT"]=m_uiCurrentTime;
            checkPoint["DENDRO_RK45_STEP_CURRENT"]=m_uiCurrentStep;
            checkPoint["DENDRO_RK45_TIME_STEP_SIZE"]=m_uiT_h;
            checkPoint["DENDRO_RK45_LAST_IO_TIME"]=m_uiCurrentTime;

            checkPoint["DENDRO_RK45_WAVELET_TOLERANCE"]=massgrav::MASSGRAV_WAVELET_TOL;
            checkPoint["DENDRO_RK45_LOAD_IMB_TOLERANCE"]=massgrav::MASSGRAV_LOAD_IMB_TOL;
            checkPoint["DENDRO_RK45_NUM_VARS"]=numVars; // number of variables to restore.
            checkPoint["DENDRO_RK45_ACTIVE_COMM_SZ"]=m_uiMesh->getMPICommSize(); // (note that rank 0 is always active).
            
            checkPoint["DENDRO_BH1_X"]=m_uiBHLoc[0].x();
            checkPoint["DENDRO_BH1_Y"]=m_uiBHLoc[0].y();
            checkPoint["DENDRO_BH1_Z"]=m_uiBHLoc[0].z();
            
            
            checkPoint["DENDRO_BH2_X"]=m_uiBHLoc[1].x();
            checkPoint["DENDRO_BH2_Y"]=m_uiBHLoc[1].y();
            checkPoint["DENDRO_BH2_Z"]=m_uiBHLoc[1].z();
            

            outfile<<std::setw(4)<<checkPoint<<std::endl;
            outfile.close();

        }


    }



}


void RK_MASSGRAV::restoreCheckPoint(const char * fNamePrefix,MPI_Comm comm)
{
    unsigned int numVars=0;
    std::vector<ot::TreeNode> octree;
    json checkPoint;

    int rank;
    int npes;
    m_uiComm=comm;
    MPI_Comm_rank(m_uiComm,&rank);
    MPI_Comm_size(m_uiComm,&npes);

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
            sprintf(fName,"%s_step_%d.cp",fNamePrefix,cpIndex);
            std::ifstream infile(fName);
            if(!infile) {
                std::cout<<fName<<" file open failed "<<std::endl;
                restoreStatus=1;
            }


            if(restoreStatus==0)
            {
                infile>>checkPoint;
                m_uiTimeBegin=checkPoint["DENDRO_RK45_TIME_BEGIN"];
                m_uiTimeEnd=checkPoint["DENDRO_RK45_TIME_END"];

                m_uiOrder=checkPoint["DENDRO_RK45_ELEMENT_ORDER"];
                m_uiCurrentTime=checkPoint["DENDRO_RK45_TIME_CURRENT"];
                m_uiCurrentStep=checkPoint["DENDRO_RK45_STEP_CURRENT"];
                m_uiT_h=checkPoint["DENDRO_RK45_TIME_STEP_SIZE"];


                massgrav::MASSGRAV_WAVELET_TOL=checkPoint["DENDRO_RK45_WAVELET_TOLERANCE"];
                massgrav::MASSGRAV_LOAD_IMB_TOL=checkPoint["DENDRO_RK45_LOAD_IMB_TOLERANCE"];
                numVars=checkPoint["DENDRO_RK45_NUM_VARS"];
                activeCommSz=checkPoint["DENDRO_RK45_ACTIVE_COMM_SZ"];
                
                m_uiBHLoc[0]=Point((double)checkPoint["DENDRO_BH1_X"],(double)checkPoint["DENDRO_BH1_Y"],(double)checkPoint["DENDRO_BH1_Z"]);
                m_uiBHLoc[1]=Point((double)checkPoint["DENDRO_BH2_X"],(double)checkPoint["DENDRO_BH2_Y"],(double)checkPoint["DENDRO_BH2_Z"]);

                restoreStep[cpIndex]=m_uiCurrentStep;

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

    par::Mpi_Bcast(&restoreFileIndex,1,0,m_uiComm);

    for(unsigned int cpIndex=restoreFileIndex; cpIndex<2; cpIndex++)
    {
        restoreStatus=0;
        octree.clear();
        if(!rank) std::cout<<" Trying to restore from checkpoint index : "<<cpIndex<<std::endl;

        if(!rank)
        {
            sprintf(fName,"%s_step_%d.cp",fNamePrefix,cpIndex);
            std::ifstream infile(fName);
            if(!infile) {
                std::cout<<fName<<" file open failed "<<std::endl;
                restoreStatus=1;
            }


            if(restoreStatus==0)
            {
                infile>>checkPoint;
                m_uiTimeBegin=checkPoint["DENDRO_RK45_TIME_BEGIN"];
                m_uiTimeEnd=checkPoint["DENDRO_RK45_TIME_END"];

                m_uiOrder=checkPoint["DENDRO_RK45_ELEMENT_ORDER"];
                m_uiCurrentTime=checkPoint["DENDRO_RK45_TIME_CURRENT"];
                m_uiCurrentStep=checkPoint["DENDRO_RK45_STEP_CURRENT"];
                m_uiT_h=checkPoint["DENDRO_RK45_TIME_STEP_SIZE"];


                massgrav::MASSGRAV_WAVELET_TOL=checkPoint["DENDRO_RK45_WAVELET_TOLERANCE"];
                massgrav::MASSGRAV_LOAD_IMB_TOL=checkPoint["DENDRO_RK45_LOAD_IMB_TOLERANCE"];
                numVars=checkPoint["DENDRO_RK45_NUM_VARS"];
                activeCommSz=checkPoint["DENDRO_RK45_ACTIVE_COMM_SZ"];
                
                m_uiBHLoc[0]=Point((double)checkPoint["DENDRO_BH1_X"],(double)checkPoint["DENDRO_BH1_Y"],(double)checkPoint["DENDRO_BH1_Z"]);
                m_uiBHLoc[1]=Point((double)checkPoint["DENDRO_BH2_X"],(double)checkPoint["DENDRO_BH2_Y"],(double)checkPoint["DENDRO_BH2_Z"]);
            }


        }

        par::Mpi_Allreduce(&restoreStatus,&restoreStatusGlobal,1,MPI_MAX,m_uiComm);
        if(restoreStatusGlobal==1) continue;


        par::Mpi_Bcast(&m_uiTimeBegin,1,0,comm);
        par::Mpi_Bcast(&m_uiTimeEnd,1,0,comm);

        par::Mpi_Bcast(&m_uiCurrentTime,1,0,comm);
        par::Mpi_Bcast(&m_uiCurrentStep,1,0,comm);

        par::Mpi_Bcast(&m_uiT_h,1,0,comm);


        par::Mpi_Bcast(& massgrav::MASSGRAV_WAVELET_TOL,1,0,comm);
        par::Mpi_Bcast(&massgrav::MASSGRAV_LOAD_IMB_TOL,1,0,comm);

        par::Mpi_Bcast(&numVars,1,0,comm);
        par::Mpi_Bcast(&m_uiOrder,1,0,comm);
        par::Mpi_Bcast(&m_uiT_h,1,0,comm);

        par::Mpi_Bcast(&activeCommSz,1,0,comm);
        
        par::Mpi_Bcast(m_uiBHLoc,2,0,comm);

        if(activeCommSz>npes)
        {
            if(!rank)std::cout<<" [Error] : checkpoint file written from  a larger communicator than the current global comm. (i.e. communicator shrinking not allowed in the restore step. )"<<std::endl;
            exit(0);
        }



        bool isActive=(rank<activeCommSz);

        MPI_Comm newComm;
        par::splitComm2way(isActive,&newComm,m_uiComm);

        if(isActive) {

            int activeRank;
            int activeNpes;

            MPI_Comm_rank(newComm, &activeRank);
            MPI_Comm_size(newComm, &activeNpes);
            assert(activeNpes == activeCommSz);

            sprintf(fName, "%s_octree_%d_%d.oct", fNamePrefix,cpIndex,activeRank);
            restoreStatus=io::checkpoint::readOctFromFile(fName, octree);
            assert(par::test::isUniqueAndSorted(octree, newComm));

        }

        par::Mpi_Allreduce(&restoreStatus,&restoreStatusGlobal,1,MPI_MAX,m_uiComm);
        if(restoreStatusGlobal==1) {

            if(!rank) std::cout<<"[Error]: octree (*.oct) restore file currupted "<<std::endl;
            continue;
        }

        newMesh=new ot::Mesh(octree,1,m_uiOrder,activeCommSz,m_uiComm);

        for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
        {
            delete [] m_uiPrevVar[index];
            delete [] m_uiVar[index];
            delete [] m_uiVarIm[index];
            delete [] m_uiUnzipVar[index];
            delete [] m_uiUnzipVarRHS[index];

            m_uiPrevVar[index]=newMesh->createVector<DendroScalar>();
            m_uiVar[index]=newMesh->createVector<DendroScalar>();
            m_uiVarIm[index]=newMesh->createVector<DendroScalar>();
            m_uiUnzipVar[index]=newMesh->createUnZippedVector<DendroScalar>();
            m_uiUnzipVarRHS[index]=newMesh->createUnZippedVector<DendroScalar>();

        }

        for(unsigned int stage=0; stage<m_uiNumRKStages; stage++)
            for(unsigned int index=0; index<massgrav::MASSGRAV_NUM_VARS; index++)
            {
                delete [] m_uiStage[stage][index];
                m_uiStage[stage][index]=newMesh->createVector<DendroScalar>();
            }

        // deallocate constraint vars allocate them for the new mesh.
        for(unsigned int index=0; index<massgrav::MASSGRAV_CONSTRAINT_NUM_VARS; index++)
        {
            delete [] m_uiConstraintVars[index];
            delete [] m_uiUnzipConstraintVars[index];

            m_uiConstraintVars[index]=newMesh->createVector<DendroScalar>();
            m_uiUnzipConstraintVars[index]=newMesh->createUnZippedVector<DendroScalar>();

        }

        const char ** varNames=massgrav::MASSGRAV_VAR_NAMES;


        if(isActive) {

            int activeRank;
            int activeNpes;

            MPI_Comm_rank(newComm, &activeRank);
            MPI_Comm_size(newComm, &activeNpes);
            assert(activeNpes == activeCommSz);

            /*for(unsigned int i=0;i<numVars;i++)
            {
                sprintf(fName,"%s_%s_%d_%d.var",fNamePrefix,varNames[i],cpIndex,activeRank);
                restoreStatus=io::checkpoint::readVecFromFile(fName,newMesh,m_uiPrevVar[i]);
                if(restoreStatus==1) break;

            }*/

            sprintf(fName,"%s_%d_%d.var",fNamePrefix,cpIndex,activeRank);
            restoreStatus=io::checkpoint::readVecFromFile(fName,newMesh,m_uiPrevVar,massgrav::MASSGRAV_NUM_VARS);


        }
        MPI_Comm_free(&newComm);
        par::Mpi_Allreduce(&restoreStatus,&restoreStatusGlobal,1,MPI_MAX,m_uiComm);
        if(restoreStatusGlobal==1) {

            if(!rank) std::cout<<"[Error]: varible (*.var) restore file currupted "<<std::endl;
            continue;
        }

        std::swap(m_uiMesh,newMesh);
        delete newMesh;
        reallocateMPIResources();
        if(restoreStatusGlobal==0) break;

    }

    if(restoreStatusGlobal==1) {
        std::cout<<"rank: "<<rank<<"[Error]: rk solver restore error "<<std::endl;
        exit(0);
    }

    unsigned int localSz=m_uiMesh->getNumLocalMeshElements();
    unsigned int totalElems;

    par::Mpi_Allreduce(&localSz,&totalElems,1,MPI_SUM,m_uiComm);


    if(!rank) std::cout<<" checkpoint at step : "<<m_uiCurrentStep<<"active Comm. sz: "<<activeCommSz<<" restore successful: "<<" restored mesh size: "<<totalElems<<std::endl;
    return;


}


} // end of namespace solver

} // end of namespace ode
