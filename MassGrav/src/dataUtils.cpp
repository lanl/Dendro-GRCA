//
// Created by milinda on 1/16/19.
//

#include "dataUtils.h"

namespace massgrav
{

    void extractBHCoords(const ot::Mesh* pMesh, const DendroScalar* var,double tolerance, const Point* ptIn, unsigned int numPt,Point* ptOut)
    {
        if((pMesh->isActive()))
        {

            unsigned int rankActive=pMesh->getMPIRank();
            unsigned int npesActive=pMesh->getMPICommSize();

            MPI_Comm commActive=pMesh->getMPICommunicator();

            const ot::TreeNode* allElements=&(*(pMesh->getAllElements().begin()));
            const unsigned int * e2n=&(*(pMesh->getE2NMapping().begin()));
            const unsigned int * e2n_dg=&(*(pMesh->getE2NMapping_DG().begin()));
            const unsigned int * cgToDg=&(*(pMesh->getCG2DGMap().begin()));

            unsigned int lookup=0;
            unsigned int ownerID, ii_x, jj_y, kk_z;

            ot::TreeNode tmpOct;
            const unsigned int eleOrder=pMesh->getElementOrder();
            double hx,x,y,z;


            std::vector<Point> ptList;
            for(unsigned int node=pMesh->getNodeLocalBegin();node<pMesh->getNodeLocalEnd();node++)
            {

                if(var[node]<tolerance)
                {
                    lookup=cgToDg[node];
                    pMesh->dg2eijk(lookup,ownerID,ii_x,jj_y,kk_z);
                    tmpOct=allElements[ownerID];
                    hx=(tmpOct.maxX()-tmpOct.minX())/((double) eleOrder);
                    x=tmpOct.minX() + ii_x*(hx);
                    y=tmpOct.minY() + jj_y*(hx);
                    z=tmpOct.minZ() + kk_z*(hx);

                    ptList.push_back(Point(GRIDX_TO_X(x),GRIDY_TO_Y(y),GRIDZ_TO_Z(z)));
                    
                    //std::cout<<" x: "<<GRIDX_TO_X(x)<<" y: "<<GRIDY_TO_Y(y)<<" z: "<<GRIDZ_TO_Z(z)<<std::endl;
                    

                }

            }

            
            
            

            std::vector<Point> * ptCluster=new std::vector<Point>[numPt];

            double min;
            unsigned int cID;
            for(unsigned int pt=0;pt<ptList.size();pt++)
            {
                cID=0;
                min=(ptIn[0]-ptList[pt]).abs();
                for(unsigned int c=1;c<numPt;c++)
                {
                    if(min>=(ptIn[c]-ptList[pt]).abs())
                    {
                        min=(ptIn[c]-ptList[pt]).abs();
                        cID=c;
                    }
                }

                ptCluster[cID].push_back(ptList[pt]);
            }


            ptList.clear();
            Point * ptMean=new Point[numPt];
            DendroIntL* ptCounts=new DendroIntL[numPt];
            DendroIntL* ptCounts_g=new DendroIntL[numPt];

            for(unsigned int c=0;c<numPt;c++)
            {
                ptMean[c]=Point(0,0,0);
                ptOut[c]=Point(0,0,0);

            }

            for(unsigned int c=0;c<numPt;c++)
            {
                ptCounts[c]=ptCluster[c].size();
                for(unsigned int pt=0;pt<ptCluster[c].size();pt++)
                    ptMean[c]+=ptCluster[c][pt];
            }

            
            
            par::Mpi_Allreduce(ptCounts,ptCounts_g,numPt,MPI_SUM,commActive);
            par::Mpi_Allreduce(ptMean,ptOut,numPt,par::Mpi_datatype<Point>::_SUM(),commActive);

            
            for(unsigned int c=0;c<numPt;c++)
               ptOut[c]/=(double)ptCounts_g[c];
            
            
            // if(pMesh->getMPIRank()==0)std::cout<<"bh1 in : "<<ptIn[0].x()<<", "<<ptIn[0].y()<<", "<<ptIn[0].z()<<std::endl;
            // if(pMesh->getMPIRank()==0)std::cout<<"bh2 in : "<<ptIn[1].x()<<", "<<ptIn[1].y()<<", "<<ptIn[1].z()<<std::endl;
            
            // if(pMesh->getMPIRank()==0)std::cout<<"bh1 out: "<<ptOut[0].x()<<", "<<ptOut[0].y()<<", "<<ptOut[0].z()<<std::endl;
            // if(pMesh->getMPIRank()==0)std::cout<<"bh2 out: "<<ptOut[1].x()<<", "<<ptOut[1].y()<<", "<<ptOut[1].z()<<std::endl;
            
            


            delete [] ptCounts;
            delete [] ptCounts_g;
            delete [] ptMean;
            delete [] ptCluster;

        }

    }

    void writeBHCoordinates(const ot::Mesh* pMesh,const Point* ptLocs, unsigned int numPt,unsigned int timestep)
    {
        
        unsigned int rankGlobal=pMesh->getMPIRankGlobal();
        if(!rankGlobal)
        {

            std::ofstream fileGW;
            char fName[256];
            sprintf(fName,"%s_BHLocations.dat",massgrav::MASSGRAV_PROFILE_FILE_PREFIX.c_str());
            fileGW.open (fName,std::ofstream::app);

            // writes the header
            if(timestep==0)
                fileGW<<"TimeStep\t"<<" bh1_x\t"<<" bh1_y\t"<<" bh1_z\t"<<" bh2_x\t"<<" bh2_y\t"<<" bh2_z\t"<<std::endl;

            fileGW<<timestep<<"\t"<<ptLocs[0].x()<<"\t"<<ptLocs[0].y()<<"\t"<<ptLocs[0].z()<<"\t"<<ptLocs[1].x()<<"\t"<<ptLocs[1].y()<<"\t"<<ptLocs[1].z()<<std::endl;
            fileGW.close();
            return;

        }
    }

    bool isRemeshBH(const ot::Mesh* pMesh, const Point* bhLoc, const double* r)
    {
        //Point bh_1 = Point(X_TO_GRIDX(bhLoc[0].x()),Y_TO_GRIDY(bhLoc[0].y()),Z_TO_GRIDZ(bhLoc[0].z()));
        //Point bh_2 = Point(X_TO_GRIDX(bhLoc[1].x()),Y_TO_GRIDY(bhLoc[1].y()),Z_TO_GRIDZ(bhLoc[1].z()));

        const double rx= r[0];
        const double ry= r[1];
        const double rz= r[2];

        const double rrx= 1.5 * r[0];
        const double rry= 1.5 * r[1];
        const double rrz= 1.5 * r[2];

        const unsigned int eleLocalBegin = pMesh->getElementLocalBegin();
        const unsigned int eleLocalEnd = pMesh->getElementLocalEnd();
        bool isOctChange=false;
        bool isOctChange_g =false;
        Point d1, d2;

        if(pMesh->isActive())
        {
            ot::TreeNode * pNodes = (ot::TreeNode*) &(*(pMesh->getAllElements().begin()));
            Point temp;

            for(unsigned int ele = eleLocalBegin; ele< eleLocalEnd; ele++)
                pNodes[ele].setFlag(((OCT_NO_CHANGE<<NUM_LEVEL_BITS)|pNodes[ele].getLevel()));

            // refine pass. 
            for(unsigned int ele = eleLocalBegin; ele< eleLocalEnd; ele++)
            {
                temp = Point( GRIDX_TO_X(pNodes[ele].minX()), GRIDY_TO_Y(pNodes[ele].minY()), GRIDZ_TO_Z(pNodes[ele].minZ()) ); 
                d1 = temp -bhLoc[0]; 
                d2 = temp -bhLoc[1];

                if( (fabs(d1.x()) <=rx && fabs(d1.y())<=ry && fabs(d1.z())<=rz ) || (fabs(d2.x()) <=rx && fabs(d2.y())<=ry && fabs(d2.z())<=rz) )
                {   
                    if((pNodes[ele].getLevel()+MAXDEAPTH_LEVEL_DIFF+1)<m_uiMaxDepth)
                        pNodes[ele].setFlag(((OCT_SPLIT<<NUM_LEVEL_BITS)|pNodes[ele].getLevel()));

                }
            
            }

            // coarsen pass. 
            for(unsigned int ele = eleLocalBegin; ele< eleLocalEnd; ele++)
            {
                if( ((ele + NUM_CHILDREN-1) < eleLocalEnd) && (pNodes[ele].getLevel()>1 && pNodes[ele].getParent() == pNodes[ele + NUM_CHILDREN-1].getParent() )) 
                {
                    bool isCoarsen=false;
                    for(unsigned int child=0; child < NUM_CHILDREN ; child ++)
                    {
                        temp = Point( GRIDX_TO_X(pNodes[ele + child].minX()), GRIDY_TO_Y(pNodes[ele + child].minY()), GRIDZ_TO_Z(pNodes[ele + child].minZ()) ); 
                        d1 = temp -bhLoc[0]; d2= temp -bhLoc[1];
                        if( ( (! ( (fabs(d1.x()) <=rx && fabs(d1.y())<=ry && fabs(d1.z())<=rz ) || (fabs(d2.x()) <=rx && fabs(d2.y())<=ry && fabs(d2.z())<=rz) ) )) && ( (fabs(d1.x()) <=rrx && fabs(d1.y())<=rry && fabs(d1.z())<=rrz ) || (fabs(d2.x()) <=rrx && fabs(d2.y())<=rry && fabs(d2.z())<=rrz ) ) && ((pNodes[ele+child].getFlag()>>NUM_LEVEL_BITS)!=OCT_SPLIT) )
                            isCoarsen=true;
                        else
                            isCoarsen=false;
                    }
                    
                    if(isCoarsen)
                    {
                        for(unsigned int child=0; child < NUM_CHILDREN ; child ++)
                            pNodes[ele + child].setFlag(((OCT_COARSE<<NUM_LEVEL_BITS)|pNodes[ele+child].getLevel()));

                        ele += NUM_CHILDREN-1;
                    }

                    
                        
                }
            }

            
            for(unsigned int ele=eleLocalBegin;ele<eleLocalEnd;ele++)
                if((pNodes[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT) // trigger remesh only when some refinement occurs
                { 
                    isOctChange=true;
                    break;
                }


        }

        bool isOctChanged_g;
        MPI_Allreduce(&isOctChange,&isOctChanged_g,1,MPI_CXX_BOOL,MPI_LOR,pMesh->getMPIGlobalCommunicator());
        //if(!m_uiGlobalRank) std::cout<<"is oct changed: "<<isOctChanged_g<<std::endl;
        return isOctChanged_g;

        

    }


    bool isRemeshEH(const ot::Mesh* pMesh, const double ** unzipVec, unsigned int vIndex, double refine_th, double coarsen_th, bool isOverwrite)
    {
        const unsigned int eleLocalBegin = pMesh->getElementLocalBegin();
        const unsigned int eleLocalEnd = pMesh->getElementLocalEnd();
        bool isOctChange=false;
        bool isOctChange_g =false;
        const unsigned int eOrder = pMesh->getElementOrder();

        if(pMesh->isActive())
        {
            ot::TreeNode * pNodes = (ot::TreeNode*) &(*(pMesh->getAllElements().begin()));
            
            if(isOverwrite)
            for(unsigned int ele = eleLocalBegin; ele< eleLocalEnd; ele++)
                pNodes[ele].setFlag(((OCT_NO_CHANGE<<NUM_LEVEL_BITS)|pNodes[ele].getLevel()));


            const std::vector<ot::Block> blkList = pMesh->getLocalBlockList();
            unsigned int sz[3];
            unsigned int ei[3];
            
            // refine test
            for(unsigned int b=0; b< blkList.size(); b++)
            {
                const ot::TreeNode blkNode = blkList[b].getBlockNode();

                sz[0]=blkList[b].getAllocationSzX();
                sz[1]=blkList[b].getAllocationSzY();
                sz[2]=blkList[b].getAllocationSzZ();

                const unsigned int bflag = blkList[b].getBlkNodeFlag();
                const unsigned int offset = blkList[b].getOffset();

                const unsigned int regLev=blkList[b].getRegularGridLev();
                const unsigned int eleIndexMax=(1u<<(regLev-blkNode.getLevel()))-1;
                const unsigned int eleIndexMin=0;

                for(unsigned int ele = blkList[b].getLocalElementBegin(); ele< blkList[b].getLocalElementEnd(); ele++)
                {
                    ei[0]=(pNodes[ele].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                    ei[1]=(pNodes[ele].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                    ei[2]=(pNodes[ele].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);

                    if((bflag &(1u<<OCT_DIR_LEFT)) && ei[0]==eleIndexMin)   continue;
                    if((bflag &(1u<<OCT_DIR_DOWN)) && ei[1]==eleIndexMin)   continue;
                    if((bflag &(1u<<OCT_DIR_BACK)) && ei[2]==eleIndexMin)   continue;

                    if((bflag &(1u<<OCT_DIR_RIGHT)) && ei[0]==eleIndexMax)  continue;
                    if((bflag &(1u<<OCT_DIR_UP)) && ei[1]==eleIndexMax)     continue;
                    if((bflag &(1u<<OCT_DIR_FRONT)) && ei[2]==eleIndexMax)  continue;

                    // refine test. 
                    for(unsigned int k=3; k< eOrder+1 +   3; k++)
                     for(unsigned int j=3; j< eOrder+1 +  3; j++)
                      for(unsigned int i=3; i< eOrder+1 + 3; i++)
                      {
                          if ( unzipVec[vIndex][offset + (ei[2]*eOrder + k)*sz[0]*sz[1] + (ei[1]*eOrder + j)*sz[0] + (ei[0]*eOrder + i)] < refine_th)
                          {
                            if( (pNodes[ele].getLevel() + MAXDEAPTH_LEVEL_DIFF +1) < m_uiMaxDepth  )
                                pNodes[ele].setFlag(((OCT_SPLIT<<NUM_LEVEL_BITS)|pNodes[ele].getLevel()));

                          }

                      }
                    
                    
                }

            }

            //coarsen test. 
            for(unsigned int b=0; b< blkList.size(); b++)
            {
                const ot::TreeNode blkNode = blkList[b].getBlockNode();

                sz[0]=blkList[b].getAllocationSzX();
                sz[1]=blkList[b].getAllocationSzY();
                sz[2]=blkList[b].getAllocationSzZ();

                const unsigned int bflag = blkList[b].getBlkNodeFlag();
                const unsigned int offset = blkList[b].getOffset();

                const unsigned int regLev=blkList[b].getRegularGridLev();
                const unsigned int eleIndexMax=(1u<<(regLev-blkNode.getLevel()))-1;
                const unsigned int eleIndexMin=0;

                if((eleIndexMax==0) || (bflag!=0)) continue; // this implies the blocks with only 1 child and boundary blocks.

                for(unsigned int ele = blkList[b].getLocalElementBegin(); ele< blkList[b].getLocalElementEnd(); ele++)
                {
                    assert(pNodes[ele].getParent()==pNodes[ele+NUM_CHILDREN-1].getParent());
                    bool isCoarsen =true;

                    for(unsigned int child=0;child<NUM_CHILDREN;child++)
                    {
                        if((pNodes[ele+child].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT)
                        {
                            isCoarsen=false;
                            break;
                        }

                    }

                    if(isCoarsen && pNodes[ele].getLevel()>1)
                    {
                        bool coarse = true;
                        for(unsigned int child=0;child<NUM_CHILDREN;child++)
                        {
                            ei[0]=(pNodes[ele + child].getX()-blkNode.getX())>>(m_uiMaxDepth-regLev);
                            ei[1]=(pNodes[ele + child].getY()-blkNode.getY())>>(m_uiMaxDepth-regLev);
                            ei[2]=(pNodes[ele + child].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLev);

                            for(unsigned int k=3; k< eOrder+1 + 3; k++)
                            for(unsigned int j=3; j< eOrder+1 +3; j++)
                             for(unsigned int i=3; i< eOrder+ +3; i++)
                             {
                                if ( !((refine_th  < unzipVec[vIndex][offset + (ei[2]*eOrder + k)*sz[0]*sz[1] + (ei[1]*eOrder + j)*sz[0] + (ei[0]*eOrder + i)]) &&  (unzipVec[vIndex][offset + (ei[2]*eOrder + k)*sz[0]*sz[1] + (ei[1]*eOrder + j)*sz[0] + (ei[0]*eOrder + i)] <=coarsen_th ))  )
                                    coarse = false;
                             }


                        }

                        if(coarse)
                            for(unsigned int child=0;child<NUM_CHILDREN;child++)
                                pNodes[ele+child].setFlag(((OCT_COARSE<<NUM_LEVEL_BITS)|pNodes[ele].getLevel()));



                        

                    }

                    ele = ele + NUM_CHILDREN-1;

                    
                    
                    
                }

            }



            for(unsigned int ele=eleLocalBegin;ele<eleLocalEnd;ele++)
                if((pNodes[ele].getFlag()>>NUM_LEVEL_BITS)==OCT_SPLIT) // trigger remesh only when some refinement occurs (laid back remesh :)  )
                { 
                    isOctChange=true;
                    break;
                }

            

        }

        bool isOctChanged_g;
        MPI_Allreduce(&isOctChange,&isOctChanged_g,1,MPI_CXX_BOOL,MPI_LOR,pMesh->getMPIGlobalCommunicator());
        //if(!m_uiGlobalRank) std::cout<<"is oct changed: "<<isOctChanged_g<<std::endl;
        return isOctChanged_g;




    }



}// end of namespace massgrav
