//
// Created by milinda on 05/02/18.
/**
*@author Milinda Fernando
*School of Computing, University of Utah
*@brief Contains utility functions for MASSGRAV simulation.
*/

#include "mesh.h"
#include <math.h>
#include "daUtils.h"

namespace massgrav
{

    template <typename T>
    double computeConstraintL2Norm(const T* constraintVec, const T* maskVec, unsigned int lbegin, unsigned int lend,MPI_Comm comm)
    {

        double l2=0.0;
        double l2_g=0.0;
        const double MASK_THRESHOLD=1.0;
        for(unsigned int i=lbegin;i<lend;i++)
        {
            if(maskVec[i]<MASK_THRESHOLD)
                l2+=(constraintVec[i]*constraintVec[i]*maskVec[i]*maskVec[i]*maskVec[i]*maskVec[i]);
            else
                l2+=(constraintVec[i]*constraintVec[i]);
        }


        par::Mpi_Reduce(&l2,&l2_g,1,MPI_SUM,0,comm);

        return (sqrt(l2_g));


    }

    
    template <typename T>
    double computeConstraintL2Norm(const ot::Mesh* mesh, const T* constraintVec,const T* maskVector,T maskthreshoold)
    {
        double l2_g=0.0;
        if(mesh->isActive())
        {
            
            MPI_Comm comm=mesh->getMPICommunicator();
            const unsigned int eleLocalBegin=mesh->getElementLocalBegin();
            const unsigned int eleLocalEnd=mesh->getElementLocalEnd();
            
            //const ot::TreeNode* pNodes=&(*(mesh->getAllElements().begin()));
            //const unsigned int eleOrder=mesh->getElementOrder();
            //unsigned int ownerID,ii_x,jj_y,kk_z;
            //const unsigned int * e2n_cg=&(*(mesh->getE2NMapping().begin()));
            //const unsigned int * e2n_dg=&(*(mesh->getE2NMapping_DG().begin()));
            //const unsigned int nPe=mesh->getNumNodesPerElement();
            
            const unsigned int nodeLocalBegin=mesh->getNodeLocalBegin();
            const unsigned int nodeLocalEnd=mesh->getNodeLocalEnd();
                        
            double l2=0.0;
            
            DendroIntL localGridPts=0;
            DendroIntL globalGridPts=0;
            
            for(unsigned int node=nodeLocalBegin;node<nodeLocalEnd;node++)
            {
                if(maskVector[node]<maskthreshoold)
                    continue;
                
                 l2+=(constraintVec[node]*constraintVec[node]);
                 localGridPts++;
             
            }
            
            par::Mpi_Reduce(&l2,&l2_g,1,MPI_SUM,0,mesh->getMPICommunicator());
            par::Mpi_Reduce(&localGridPts,&globalGridPts,1,MPI_SUM,0,mesh->getMPICommunicator());


            if(!(mesh->getMPIRank()))
                l2_g=l2_g/(double)(globalGridPts);
            
        }
        
        return sqrt(l2_g);
        
        
    }

    template<typename T>
    double extractConstraints(const ot::Mesh* mesh,const T** constraintVar,const T* maskVec, double maskthreshoold,unsigned int timestep)
    {
        const unsigned int numConstraints=4;
        double constraintMaskedL2[numConstraints]; // remove the psi4

        unsigned int rankGlobal=mesh->getMPIRankGlobal();
        unsigned int npesGlobal=mesh->getMPICommSizeGlobal();
        MPI_Comm commGlobal=mesh->getMPIGlobalCommunicator();

        if(mesh->isActive())
        {

            unsigned int rankActive=mesh->getMPIRank();
            unsigned int npesActive=mesh->getMPICommSize();
            MPI_Comm commActive=mesh->getMPICommunicator();


            for(unsigned int index=0;index<numConstraints;index++)
            {
                constraintMaskedL2[index]=massgrav::computeConstraintL2Norm(mesh,constraintVar[index],maskVec,maskthreshoold);
                if(!rankActive)
                {
                    std::cout<<YLW<<"\tConstraint " <<massgrav::MASSGRAV_CONSTRAINT_VAR_NAMES[index]<< " L2 : ("<<constraintMaskedL2[index]<<" )"<<NRM<<std::endl;
                }

            }

            if(!rankActive)
            {

                std::ofstream fileGW;
                char fName[256];
                sprintf(fName,"%s_Constraints.dat",massgrav::MASSGRAV_PROFILE_FILE_PREFIX.c_str());
                fileGW.open (fName,std::ofstream::app);
                // writes the header
                if(timestep==0)
                    fileGW<<"TimeStep\t"<<" C_HAM\t"<<" C_MOM0\t"<<" C_MOM1\t"<<" C_MOM2\t"<<std::endl;

                fileGW<<timestep<<"\t"<<constraintMaskedL2[0]<<"\t"<<constraintMaskedL2[1]<<"\t"<<constraintMaskedL2[2]<<"\t"<<constraintMaskedL2[3]<<std::endl;
                fileGW.close();


            }



        }


        #if 0
        if(!rankGlobal)
        {
                for(unsigned int index=0;index<numConstraints;index++)
                {
                    if(constraintMaskedL2[index]>0.01)
                    {
                        if(massgrav::KO_DISS_SIGMA>0.06)
                            massgrav::KO_DISS_SIGMA=0.05;
                        else
                            massgrav::KO_DISS_SIGMA=0.10;

                        break;
                    }
                }

        }
        par::Mpi_Bcast(&massgrav::KO_DISS_SIGMA,1,0,commGlobal);
        #endif



    }


    template<typename T>
    void artificial_dissipation(ot::Mesh * pMesh , T** zipVars, unsigned int numVars, unsigned int nc, unsigned int s,bool isGhostEx)
    {

        if(!isGhostEx)
        {
            for(unsigned int var=0;var<numVars;var++)
                pMesh->performGhostExchange(zipVars[var]);
        }


        

        const unsigned int nPe = pMesh->getNumNodesPerElement();
        const unsigned int Nrp = pMesh->getElementOrder()+1;
        const unsigned int eOrder = pMesh -> getElementOrder();

        RefElement refEl(1,eOrder);
        refEl.computeFilterOp(nc,s);

        const unsigned int refElSz = refEl.getElementSz();
        const ot::TreeNode* allElements = &(*(pMesh->getAllElements().begin()));
        const unsigned int * e2n_cg = &(*(pMesh->getE2NMapping().begin()));
        const double * filterOp = refEl.getFr1D();

        T* nodalVec = new T[nPe*nPe];
        T* imV1 = new T[nPe*nPe];
        T* imV2 = new T[nPe*nPe];

        for(unsigned int var=0;var<numVars;var++)
        {
            for(unsigned int ele=pMesh->getElementLocalBegin();ele<pMesh->getElementLocalEnd();ele++)
            {
                pMesh->getElementNodalValues(zipVars[var],nodalVec,ele);

                DENDRO_TENSOR_IIAX_APPLY_ELEM(Nrp,filterOp,nodalVec,imV1);
                DENDRO_TENSOR_IAIX_APPLY_ELEM(Nrp,filterOp,imV1,imV2);
                DENDRO_TENSOR_AIIX_APPLY_ELEM(Nrp,filterOp,imV2,nodalVec);

                /*
                no need to scale by the jacobian
                const unsigned int szX = GRIDX_TO_X(allElements[ele].maxX()) - GRIDX_TO_X(allElements[ele].minX());
                const unsigned int szY = GRIDY_TO_Y(allElements[ele].maxY()) - GRIDY_TO_Y(allElements[ele].minY());
                const unsigned int szZ = GRIDZ_TO_Z(allElements[ele].maxZ()) - GRIDZ_TO_Z(allElements[ele].minZ());
                
                const double Jx = 1.0/(refElSz/(double (szX)));
                const double Jy = 1.0/(refElSz/(double (szY)));
                const double Jz = 1.0/(refElSz/(double (szZ)));
                
                for(unsigned int node=0;node<nPe;node++)
                    nodalVec[node]*= (Jx*Jy*Jz);*/

                
                for(unsigned int k=0; k < Nrp; k++)
                    for(unsigned int j=0; j < Nrp; j++)
                        for(unsigned int i=0; i < Nrp; i++)
                        {
                            if(pMesh->isNodeLocal(ele,i,j,k) && !(pMesh->isNodeHanging(ele,i,j,k)))
                                zipVars[var][e2n_cg[ele * nPe + k * Nrp * Nrp + j * Nrp + i]] = nodalVec[k * Nrp * Nrp + j * Nrp + i];
                        }

            }

        }


        delete [] nodalVec;
        delete [] imV1;
        delete [] imV2;


    }


} // end of namespace massgrav

