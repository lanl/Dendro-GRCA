/**
 * @file enuts.h
 * @author Milinda Fernando (milinda@cs.utah.edu)
 * @brief Spatially adaptive non uniform time stepping framework. 
 * @version 0.1
 * @date 2019-07-12
 * @copyright Copyright (c) 2019, School of Computing, University of Utah. 
 * 
 */
#pragma once
#include "mesh.h"
#include <functional>
#include "block.h"
#include "point.h"
#include <assert.h>
#include "ts.h"
#include "ets.h"
#include "meshUtils.h"
#include <bitset>
#include "subSM.h"
#include "oct2vtk.h"
#include <iostream>
#include "blkAsync.h"

namespace ts
{   


    /**
     * @brief : simple structure to support storing of a single time step (explicit methods)
     * note that the stages are numbered from 1 to m_uiNumStages. 
     * @tparam T : vector type
     */
    template <typename T>
    struct BlockTimeStep
    {

        public:
            /**@brief: ets stages
             * stage 0    : input vector
             * stage 1    : ets stage 1
             *            .
             *            .
             * stage p    : ets stage p
             * stage p+1  : output vector
            */
            std::vector< BlockAsyncVector<T> > _vec;

            /**@brief: rk stage*/
            unsigned int  _rks  = ETS_STAGE_DEFAULT;

            /**@brief: block time*/
            unsigned int  _time = LOOK_UP_TABLE_DEFAULT;

            /**
             * @brief allocate the Block time step vector. 
             * @param numVec : number of vectors per block. 
             * @param blkid  : block id
             * @param sz     : sizes of the each dimension
             * @param dof    : degrees of freedoms. 
             */
            void alloc_vec(unsigned int numVec, unsigned int blkid ,const unsigned int *sz , unsigned int dof=1)
            {
                _vec.resize(numVec);
                for(unsigned int k=0; k < _vec.size(); k++)
                    _vec[k].createVec(blkid,sz, false, BLK_ASYNC_VEC_MODE::BLK_UNZIP, dof);
            }

            /**
             * @brief deallocate the vectors. 
             */
            void dealloc_vec()
            {
                for(unsigned int k=0; k < _vec.size(); k++)
                    _vec[k].destroyVec();

                _vec.clear();
            }


    };


    
    /**
    * @brief Basic: class for performing non-uniform time stepping. (explicit time stepping)
    * In order to perform non uniform time stepping the octree to block decomposition
    * needed to be completed. We assume that the blocks are setup in the mesh class. 
    * 
    * Note: The stages are numbered from 1 to m_uiNumStages. 
    */
    template<typename T, typename Ctx>
    class ExplicitNUTS  : public ETS<T,Ctx>
    {
        using ETS<T,Ctx>::m_uiAppCtx;
        using ETS<T,Ctx>::m_uiNumStages;
        using ETS<T,Ctx>::m_uiEVar;
        using ETS<T,Ctx>::m_uiStVec;
        using ETS<T,Ctx>::m_uiEVecTmp;
        using ETS<T,Ctx>::m_uiTimeInfo;
        using ETS<T,Ctx>::m_uiAij;
        using ETS<T,Ctx>::m_uiBi;
        using ETS<T,Ctx>::m_uiCi;
        using ETS<T,Ctx>::m_uiType;
        
        protected:

            /**@brief: minimum level of the grid. */
            unsigned int m_uiLevMin;

            /**@brief: max level of the grid. */
            unsigned int m_uiLevMax;

            /**@brief: element to block map. */
            std::vector< unsigned int > m_uiE2B; 
            
            /**@brief: unzip vector for m_uiEVar*/
            DVec m_uiEvarUzip;

            /**@brief: DG vector for m_uiEVar*/
            DVec m_uiEvarDG;

            /**@brief: List of block async vector. */
            std::vector<ts::BlockTimeStep<T>> m_uiBVec;

            /**@brief : keep track of the element time. */
            std::vector<unsigned int> m_uiEleTime;

            /**@brief: explicit timer interpolation operators for ENUTS */
            ENUTSOp* m_uiECOp = NULL;


        private:

            /**@brief: allocates the data stuctures and initialize with the current mesh stores in the class. */
            void init_data_structures();

            /**@brief: freee the allocated data structures. */
            void free_data_structures();

            
        private:
            /**
             * @brief Allocates internal variables for the time stepper. 
             * @return int 
             */
            virtual int allocate_internal_vars();

            /**@brief: Deallocate internal variables. */
            virtual int deallocate_internal_vars(); 

            virtual void sync_blk_timestep(unsigned int blk, unsigned int rk_s);

            virtual void update_ele_timestep();


            

        public:

            /**@brief: constructor
             * Assumptions: Note that blocks result in from octree to block decomposition can be not 2:1 balanced.
             * But to perform NUTS, blocks should be 2:1 balanced. 
            */
            ExplicitNUTS(Ctx* ctx);

            /**@brief: default destructor */
            ~ExplicitNUTS();

            /**@brief: build all the required data structures for the non uniform TS*/
            void init();

            /**@brief : evolve the variables to the next coarsest time step (i.e. loop over the block sequence. ) Evolution in the sense of the  NUTS*/
            virtual void evolve();
            
            /**@brief: returns  a constant pointer to the sub scatter maps. */
            //inline ot::SubScatterMap* const get_sub_scatter_maps() const {return m_uiSubSM;}

            /**@brief: update all the internal data strucutures, with a new mesh pointer. */
            virtual int sync_with_mesh();

            /**@brief returns the dt min value */
            T get_dt_min() const { return m_uiTimeInfo._m_uiTh; } 
            
            /**@brief returns the dt max value */
            T get_dt_max() const { return m_uiTimeInfo._m_uiTh*(1u<<(m_uiLevMax-m_uiLevMin)); }

            /**@brief: prints the load balance statistis*/
            void  dump_load_statistics(std::ostream & sout) const ;
        
    };

    template<typename T, typename Ctx>
    ExplicitNUTS<T,Ctx>::ExplicitNUTS(Ctx* ctx) : ETS<T,Ctx>(ctx)
    {
        this->init_data_structures();
        
    }

    template<typename T, typename Ctx>
    ExplicitNUTS<T,Ctx>::~ExplicitNUTS()
    {
        this->deallocate_internal_vars();
        this->free_data_structures();

    }

    template<typename T, typename Ctx>
    void ExplicitNUTS<T,Ctx>::init_data_structures()
    {

         // identify dependent and independent blocks. 
        const ot::Mesh* pMesh = m_uiAppCtx->get_mesh();
        const bool isActive = pMesh->isActive();
        
        pMesh->computeMinMaxLevel(m_uiLevMin,m_uiLevMax);
        par::Mpi_Bcast(&m_uiLevMin,1,0,pMesh->getMPIGlobalCommunicator());
        par::Mpi_Bcast(&m_uiLevMax,1,0,pMesh->getMPIGlobalCommunicator());
        assert( (m_uiLevMin > 0 ) && (m_uiLevMax <= m_uiMaxDepth) );
        
        

    }

    template<typename T, typename Ctx>
    void ExplicitNUTS<T,Ctx>::free_data_structures()
    {
        
        return ;

    }
    
    template<typename T, typename Ctx>
    int ExplicitNUTS<T,Ctx>::allocate_internal_vars()
    {

        assert(m_uiNumStages>0);
        const unsigned int DOF = m_uiEVar.GetDof();

        m_uiStVec.resize(m_uiNumStages);
        
        // for(unsigned int i=0; i < m_uiNumStages; i++)
        //     m_uiStVec[i].VecCreate(m_uiAppCtx->get_mesh(), false , true, false , m_uiEVar.GetDof());

        for(unsigned int i=0; i < m_uiNumStages; i++)
            m_uiStVec[i].VecCreateDG(m_uiAppCtx->get_mesh(), true, m_uiEVar.GetDof());

        m_uiEVecTmp[0].VecCreate(m_uiAppCtx->get_mesh(), m_uiEVar.IsGhosted() , m_uiEVar.IsUnzip(), m_uiEVar.IsElemental() , m_uiEVar.GetDof());
        m_uiEVecTmp[1].VecCreate(m_uiAppCtx->get_mesh(), m_uiEVar.IsGhosted() , m_uiEVar.IsUnzip(), m_uiEVar.IsElemental() , m_uiEVar.GetDof());

        m_uiEvarUzip.VecCreate(m_uiAppCtx->get_mesh(), false , true, false ,m_uiEVar.GetDof());
        m_uiEvarDG.VecCreateDG(m_uiAppCtx->get_mesh(), true, m_uiEVar.GetDof());

        const ot::Mesh* pMesh = m_uiAppCtx->get_mesh();
        if(pMesh->isActive())
        {
            m_uiEleTime.resize(pMesh->getAllElements().size(),0);
            std::vector<ot::Block> blkList = pMesh->getLocalBlockList();
            m_uiBVec.resize(blkList.size());

            for(unsigned int blk =0; blk < blkList.size(); blk++)
            {
                const unsigned int sz[3] = { blkList[blk].getAllocationSzX(), blkList[blk].getAllocationSzY(), blkList[blk].getAllocationSzZ()};
                m_uiBVec[blk].alloc_vec(m_uiNumStages+2, blk, sz, DOF);
            }

        }

        m_uiECOp = new ENUTSOp(m_uiType);
        return 0;

        
    }

    template<typename T, typename Ctx>
    int ExplicitNUTS<T,Ctx>::deallocate_internal_vars()
    {
        for(unsigned int i=0; i < m_uiNumStages; i++)
        {
            this->m_uiStVec[i].VecDestroy();
        }

        this->m_uiStVec.clear();
        this->m_uiEVecTmp[0].VecDestroy();
        this->m_uiEVecTmp[1].VecDestroy();

        m_uiEvarUzip.VecDestroy();
        m_uiEvarDG.VecDestroy();


        for(unsigned int k=0; k < m_uiBVec.size(); k++)
        {
            for(unsigned int j=0; j < m_uiBVec[k]._vec.size(); j++)
                m_uiBVec[k]._vec[j].destroyVec();
            
            m_uiBVec[k]._vec.clear();
        }
        
        m_uiBVec.clear();

        delete m_uiECOp;
        return 0;

    }

    template<typename T, typename Ctx>
    int ExplicitNUTS<T,Ctx>::sync_with_mesh()
    {
        if(m_uiAppCtx -> is_ets_synced())
            return 0;

        this->deallocate_internal_vars();
        this->free_data_structures();

        this->init_data_structures();
        this->allocate_internal_vars();

        m_uiEVar = m_uiAppCtx->get_evolution_vars();
        m_uiAppCtx->set_ets_synced(true);

        return 0;


    }

    template<typename T, typename Ctx>
    void ExplicitNUTS<T,Ctx>::init()
    {

        m_uiAppCtx->initialize();
        m_uiTimeInfo = m_uiAppCtx->get_ts_info();
        allocate_internal_vars();

    }


    template<typename T, typename Ctx>
    void ExplicitNUTS<T,Ctx>::sync_blk_timestep(unsigned int blk, unsigned int rk_s)
    {

       
        ot::Mesh* pMesh = (ot::Mesh*)m_uiAppCtx->get_mesh();
        assert(rk_s>=1 && rk_s <= m_uiNumStages);
        
        if((!(pMesh->isActive())) || m_uiBVec[blk]._vec[rk_s].isSynced() )
            return;

        const unsigned int * etVec = m_uiEleTime.data();
        const ot::Block* blkList     =   pMesh->getLocalBlockList().data();
        const unsigned int regLevel  =   blkList[blk].getRegularGridLev();
        const ot::TreeNode* pNodes   =   pMesh->getAllElements().data();
        const ot::TreeNode blkNode   =   blkList[blk].getBlockNode();
        const unsigned int eOrder    =   pMesh->getElementOrder();
        const unsigned int nPe       =   pMesh->getNumNodesPerElement();

        MPI_Comm comm = pMesh->getMPICommunicator();
        const unsigned int paddWidth =   blkList[blk].get1DPadWidth();

        assert(paddWidth>0);
        assert(rk_s>=1);
        
        const unsigned int lx     =  blkList[blk].getAllocationSzX();
        const unsigned int ly     =  blkList[blk].getAllocationSzY();
        const unsigned int lz     =  blkList[blk].getAllocationSzZ();
        const unsigned int offset =  blkList[blk].getOffset(); 

        const unsigned int dgSz  =  pMesh->getAllElements().size() * pMesh->getNumNodesPerElement();
        const unsigned int cgSz  =  pMesh->getDegOfFreedom();

        const unsigned int* e2n  =  pMesh->getE2NMapping().data();
        const unsigned int* e2e  =  pMesh->getE2EMapping().data();


        const unsigned int bLev  =  pNodes[blkList[blk].getLocalElementBegin()].getLevel();
        const unsigned int bTime =  etVec[blkList[blk].getLocalElementBegin()];
        T dt;
        const T dt_f             =  (1u<<(m_uiLevMax-bLev-1))*(m_uiTimeInfo._m_uiTh);
        const T dt_c             =  (1u<<(m_uiLevMax-bLev+1))*(m_uiTimeInfo._m_uiTh);
        unsigned int tl=0;
        unsigned int lookUp;
        
        const unsigned int cSz[3] = { eOrder+1, eOrder+1, eOrder+1 };
        
        unsigned int fchild[4];
        unsigned int echild[2];

        

        const unsigned int dof = m_uiEVar.GetDof();

        std::vector<T> cVec;
        cVec.resize(nPe*rk_s);

        T* sV[rk_s];    // stage pointer for each variable v. 
        T* cVin[rk_s];  // correction vec. pointer in,
        T* cVout[rk_s]; // correction vec. pointer out,

        T * dgWVec = m_uiEvarDG.GetVecArray();
        T * cgWVec = m_uiEVecTmp[0].GetVecArray();
        T * uzWVec = m_uiEvarUzip.GetVecArray();
        T * dgStages[m_uiNumStages];

        

        
        const unsigned int m_uiBlkID = blk;
        ENUTSOp* Op            = m_uiECOp;
        
        for(unsigned int i=0; i < m_uiNumStages; i++)
            dgStages[i] = m_uiStVec[i].GetVecArray();

        T * m_uiVec = ( T *) m_uiBVec[blk]._vec[rk_s].data();
        
        
        for(unsigned int v = 0; v < dof; v++)
        {
            T * vVec = m_uiVec + v*lx*ly*lz;
            T * dgV  = dgWVec  + v*dgSz;  

            const unsigned int m_uiDof   = 1; //m_uiEVar.GetDof();
            
            for(unsigned int s=1; s <=rk_s; s++)
                sV[s-1] = dgStages[s-1] + v * dgSz;


            std::vector<unsigned int> unzipEids;    
            unzipEids.reserve(4*blkList[blk].getElemSz1D() * blkList[blk].getElemSz1D()* blkList[blk].getElemSz1D());

            for(unsigned int elem = blkList[m_uiBlkID].getLocalElementBegin(); elem < blkList[m_uiBlkID].getLocalElementEnd(); elem++)
            {
                const unsigned int ei=(pNodes[elem].getX()-blkNode.getX())>>(m_uiMaxDepth-regLevel);
                const unsigned int ej=(pNodes[elem].getY()-blkNode.getY())>>(m_uiMaxDepth-regLevel);
                const unsigned int ek=(pNodes[elem].getZ()-blkNode.getZ())>>(m_uiMaxDepth-regLevel);

                const unsigned int emin = 0;
                const unsigned int emax = (1u<<(regLevel-blkNode.getLevel()))-1;

                assert(etVec[elem] == bTime );

                for(unsigned int node =0; node < nPe; node++)
                    dgV[elem*nPe + node] = sV[rk_s-1][elem*nPe +  node];


                //std::cout<<"ele  internal coppied "<<std::endl;
                unzipEids.push_back(elem);


                // OCT_DIR_LEFT
                if(ei==emin)
                { 
                    const unsigned int dir = OCT_DIR_LEFT;
                    lookUp = e2e[elem*NUM_FACES + dir];
                    if(lookUp!=LOOK_UP_TABLE_DEFAULT)
                    {
                        
                        if(pNodes[lookUp].getLevel() == bLev )
                        {
                            tl = etVec[lookUp];
                            unzipEids.push_back(lookUp);

                            assert(bTime == tl);
                            // no corrections needs to be applied. 
                            for(unsigned int node =0; node < nPe; node++)
                                dgV[lookUp *nPe + node] = sV[rk_s-1][lookUp*nPe +  node];

                        }else if(pNodes[lookUp].getLevel() < bLev)
                        {

                            tl = etVec[lookUp];
                            unzipEids.push_back(lookUp);
                             // need to do coarser to finer correction. 
                            for(unsigned int s=1; s <= rk_s; s++ )
                            {
                                cVin[s-1]  = &sV[s-1][lookUp*nPe];
                                cVout[s-1] = cVec.data() + (s-1)*nPe; 
                            }
                            
                            if(tl==bTime)
                                dt = 0;
                            else if(tl == bTime + (1u<<(m_uiLevMax - bLev)))
                                dt = dt_c/2.0;
                            else
                            {
                                std::cout<<"ENUTS sync error at  "<<__func__<<" line : "<<__LINE__<<" finer and coarse time gap violated "<<std::endl;
                                assert(false); // finer block cannot exceed time of the coarset block. 
                                MPI_Abort(comm,0);
                            }
                                
                            Op->Ccf(cVout, (const T**)cVin, cSz, rk_s, dt_c, dt_f, dt, m_uiDof );
                            for(unsigned int node =0; node < nPe; node++)
                                dgV[lookUp *nPe + node] = cVout[rk_s-1][node];

                        }else
                        {
                            assert(pNodes[lookUp].getLevel() == bLev +1);
                            // finner to coaser correction. 

                            pMesh->getFinerFaceNeighbors(elem, dir, (unsigned int *)fchild);
                            assert(etVec[fchild[0]] == etVec[fchild[1]] && etVec[fchild[1]] == etVec[fchild[2]] && etVec[fchild[2]] == etVec[fchild[3]] );
                            unzipEids.push_back(fchild[0]); unzipEids.push_back(fchild[1]); unzipEids.push_back(fchild[2]); unzipEids.push_back(fchild[3]); 
                            tl = etVec[fchild[0]];
                            if(tl==bTime)
                                dt = 0;
                            else if(tl + (1u<<(m_uiLevMax-bLev-1)) == bTime)
                                dt = dt_c/2.0;
                            else
                            {
                                std::cout<<"ENUTS sync error at  "<<__func__<<" line : "<<__LINE__<<" finer and coarse time gap violated "<<std::endl;
                                assert(false); // finer block cannot exceed time of the coarset block. 
                                MPI_Abort(comm,0);
                            }

                            for(unsigned int m=0; m < (NUM_CHILDREN>>1u); m++)
                            {
                                for(unsigned int s=1; s <= rk_s; s++ )
                                {
                                    cVin[s-1]  = &sV[s-1][fchild[m]*nPe];
                                    cVout[s-1] = cVec.data() + (s-1)*nPe; 
                                }

                                Op->Cfc(cVout, (const T**)cVin, cSz, rk_s, dt_c, dt_f, dt, m_uiDof );
                                for(unsigned int node =0; node < nPe; node++)
                                    dgV[fchild[m] *nPe + node] = cVout[rk_s-1][node];
                        
                            }

                        }
                        
                    }
                        
                }

                //std::cout<<"ele  left "<<std::endl;

                // OCT_DIR_RIGHT
                if(ei==emax)
                {
                    const unsigned int dir = OCT_DIR_RIGHT;
                    lookUp = e2e[elem*NUM_FACES + dir];
                    if(lookUp!=LOOK_UP_TABLE_DEFAULT)
                    {
                        
                        if(pNodes[lookUp].getLevel() == bLev )
                        {
                            tl = etVec[lookUp];
                            unzipEids.push_back(lookUp);

                            assert(bTime == tl);

                            for(unsigned int node =0; node < nPe; node++)
                                dgV[lookUp *nPe + node] = sV[rk_s-1][lookUp*nPe +  node];

                        }else if(pNodes[lookUp].getLevel() < bLev)
                        {

                            tl = etVec[lookUp];
                            unzipEids.push_back(lookUp);
                            
                            // need to do coarser to finer correction. 
                            for(unsigned int s=1; s <= rk_s; s++ )
                            {
                                cVin[s-1]  = &sV[s-1][lookUp*nPe];
                                cVout[s-1] = cVec.data() + (s-1)*nPe; 
                            }
                            
                            if(tl==bTime)
                                dt = 0;
                            else if(tl == bTime + (1u<<(m_uiLevMax - bLev)))
                                dt = dt_c/2.0;
                            else
                            {
                                std::cout<<"ENUTS sync error at  "<<__func__<<" line : "<<__LINE__<<" finer and coarse time gap violated "<<std::endl;
                                assert(false); // finer block cannot exceed time of the coarset block. 
                                MPI_Abort(comm,0);
                            }

                            Op->Ccf(cVout, (const T**) cVin, cSz, rk_s, dt_c, dt_f, dt, m_uiDof );
                            for(unsigned int node =0; node < nPe; node++)
                                dgV[lookUp *nPe + node] = cVout[rk_s-1][node];


                        }else
                        {
                            assert(pNodes[lookUp].getLevel() == bLev +1);
                            // finner to coaser correction. 

                            pMesh->getFinerFaceNeighbors(elem, dir, (unsigned int *)fchild);
                            assert(etVec[fchild[0]] == etVec[fchild[1]] && etVec[fchild[1]] == etVec[fchild[2]] && etVec[fchild[2]] == etVec[fchild[3]] );
                            unzipEids.push_back(fchild[0]); unzipEids.push_back(fchild[1]); unzipEids.push_back(fchild[2]); unzipEids.push_back(fchild[3]); 

                            tl = etVec[fchild[0]];
                            if(tl==bTime)
                                dt = 0;
                            else if(tl + (1u<<(m_uiLevMax-bLev-1)) == bTime)
                                dt = dt_c/2.0;
                            else
                            {
                                std::cout<<"ENUTS sync error at  "<<__func__<<" line : "<<__LINE__<<" finer and coarse time gap violated "<<std::endl;
                                assert(false); // finer block cannot exceed time of the coarset block. 
                                MPI_Abort(comm,0);
                            }

                            for(unsigned int m=0; m < (NUM_CHILDREN>>1u); m++)
                            {
                                for(unsigned int s=1; s <= rk_s; s++ )
                                {
                                    cVin[s-1]  = &sV[s-1][fchild[m]*nPe];
                                    cVout[s-1] = cVec.data() + (s-1)*nPe; 
                                }

                                Op->Cfc(cVout, (const T**) cVin, cSz, rk_s, dt_c, dt_f, dt, m_uiDof );
                                for(unsigned int node =0; node < nPe; node++)
                                    dgV[fchild[m] *nPe + node] = cVout[rk_s-1][node];
                        
                            }

                        }
                        
                    }

                }

                //std::cout<<"ele  right "<<std::endl;

                if(ej==emin)
                {
                    const unsigned int dir = OCT_DIR_DOWN;
                    lookUp = e2e[elem*NUM_FACES + dir];
                    if(lookUp!=LOOK_UP_TABLE_DEFAULT)
                    {
                        
                        if(pNodes[lookUp].getLevel() == bLev )
                        {
                            tl = etVec[lookUp];
                            unzipEids.push_back(lookUp);
                            assert(bTime == tl);

                            for(unsigned int node =0; node < nPe; node++)
                                dgV[lookUp *nPe + node] = sV[rk_s-1][lookUp*nPe +  node];

                        }else if(pNodes[lookUp].getLevel() < bLev)
                        {

                            tl = etVec[lookUp];
                            unzipEids.push_back(lookUp);
                            
                            // need to do coarser to finer correction. 
                            for(unsigned int s=1; s <= rk_s; s++ )
                            {
                                cVin[s-1]  = &sV[s-1][lookUp*nPe];
                                cVout[s-1] = cVec.data() + (s-1)*nPe; 
                            }
                            
                            if(tl==bTime)
                                dt = 0;
                            else if(tl == bTime + (1u<<(m_uiLevMax - bLev)))
                                dt = dt_c/2.0;
                            else
                            {
                                std::cout<<"ENUTS sync error at  "<<__func__<<" line : "<<__LINE__<<" finer and coarse time gap violated "<<std::endl;
                                assert(false); // finer block cannot exceed time of the coarset block. 
                                MPI_Abort(comm,0);
                            }

                            Op->Ccf(cVout, (const T**) cVin, cSz, rk_s, dt_c, dt_f, dt, m_uiDof );
                            for(unsigned int node =0; node < nPe; node++)
                                dgV[lookUp *nPe + node] = cVout[rk_s-1][node];


                        }else
                        {
                            assert(pNodes[lookUp].getLevel() == bLev +1);
                            // finner to coaser correction. 

                            pMesh->getFinerFaceNeighbors(elem, dir, (unsigned int *)fchild);
                            assert(etVec[fchild[0]] == etVec[fchild[1]] && etVec[fchild[1]] == etVec[fchild[2]] && etVec[fchild[2]] == etVec[fchild[3]] );
                            unzipEids.push_back(fchild[0]); unzipEids.push_back(fchild[1]); unzipEids.push_back(fchild[2]); unzipEids.push_back(fchild[3]); 

                            tl = etVec[fchild[0]];
                            if(tl==bTime)
                                dt = 0;
                            else if(tl + (1u<<(m_uiLevMax-bLev-1)) == bTime)
                                dt = dt_c/2.0;
                            else
                            {
                                std::cout<<"ENUTS sync error at  "<<__func__<<" line : "<<__LINE__<<" finer and coarse time gap violated "<<std::endl;
                                assert(false); // finer block cannot exceed time of the coarset block. 
                                MPI_Abort(comm,0);
                            }

                            for(unsigned int m=0; m < (NUM_CHILDREN>>1u); m++)
                            {
                                for(unsigned int s=1; s <= rk_s; s++ )
                                {
                                    cVin[s-1]  = &sV[s-1][fchild[m]*nPe];
                                    cVout[s-1] = cVec.data() + (s-1)*nPe; 
                                }

                                Op->Cfc(cVout, (const T**) cVin, cSz, rk_s, dt_c, dt_f, dt, m_uiDof );
                                for(unsigned int node =0; node < nPe; node++)
                                    dgV[fchild[m] *nPe + node] = cVout[rk_s-1][node];
                        
                            }

                        }
                        
                    }

                }

                //std::cout<<"ele  down "<<std::endl;

                if(ej==emax)
                {
                    const unsigned int dir = OCT_DIR_UP;
                    lookUp = e2e[elem*NUM_FACES + dir];
                    if(lookUp!=LOOK_UP_TABLE_DEFAULT)
                    {
                        
                        if(pNodes[lookUp].getLevel() == bLev )
                        {
                            tl = etVec[lookUp];
                            unzipEids.push_back(lookUp);
                            assert(bTime == tl);

                            for(unsigned int node =0; node < nPe; node++)
                                dgV[lookUp *nPe + node] = sV[rk_s-1][lookUp*nPe +  node];

                        }else if(pNodes[lookUp].getLevel() < bLev)
                        {

                            tl = etVec[lookUp];
                            unzipEids.push_back(lookUp);
                            
                            // need to do coarser to finer correction. 
                            for(unsigned int s=1; s <= rk_s; s++ )
                            {
                                cVin[s-1]  = &sV[s-1][lookUp*nPe];
                                cVout[s-1] = cVec.data() + (s-1)*nPe; 
                            }
                            
                            if(tl==bTime)
                                dt = 0;
                            else if(tl == bTime + (1u<<(m_uiLevMax - bLev)))
                                dt = dt_c/2.0;
                            else
                            {
                                std::cout<<"ENUTS sync error at  "<<__func__<<" line : "<<__LINE__<<" finer and coarse time gap violated "<<std::endl;
                                assert(false); // finer block cannot exceed time of the coarset block. 
                                MPI_Abort(comm,0);
                            }

                            Op->Ccf(cVout, (const T**) cVin, cSz, rk_s, dt_c, dt_f, dt, m_uiDof );
                            for(unsigned int node =0; node < nPe; node++)
                                dgV[lookUp *nPe + node] = cVout[rk_s-1][node];


                        }else
                        {
                            assert(pNodes[lookUp].getLevel() == bLev +1);
                            // finner to coaser correction. 

                            pMesh->getFinerFaceNeighbors(elem, dir, (unsigned int *)fchild);
                            assert(etVec[fchild[0]] == etVec[fchild[1]] && etVec[fchild[1]] == etVec[fchild[2]] && etVec[fchild[2]] == etVec[fchild[3]] );
                            unzipEids.push_back(fchild[0]); unzipEids.push_back(fchild[1]); unzipEids.push_back(fchild[2]); unzipEids.push_back(fchild[3]); 

                            tl = etVec[fchild[0]];
                            if(tl==bTime)
                                dt = 0;
                            else if(tl + (1u<<(m_uiLevMax-bLev-1)) == bTime)
                                dt = dt_c/2.0;
                            else
                            {
                                std::cout<<"ENUTS sync error at  "<<__func__<<" line : "<<__LINE__<<" finer and coarse time gap violated "<<std::endl;
                                assert(false); // finer block cannot exceed time of the coarset block. 
                                MPI_Abort(comm,0);
                            }

                            for(unsigned int m=0; m < (NUM_CHILDREN>>1u); m++)
                            {
                                for(unsigned int s=1; s <= rk_s; s++ )
                                {
                                    cVin[s-1]  = &sV[s-1][fchild[m]*nPe];
                                    cVout[s-1] = cVec.data() + (s-1)*nPe; 
                                }

                                Op->Cfc(cVout, (const T**) cVin, cSz, rk_s, dt_c, dt_f, dt, m_uiDof );
                                for(unsigned int node =0; node < nPe; node++)
                                    dgV[fchild[m] *nPe + node] = cVout[rk_s-1][node];
                        
                            }

                        }
                        
                    }

                }

                //std::cout<<"ele  up "<<std::endl;

                if(ek==emin)
                {

                    const unsigned int dir = OCT_DIR_BACK;
                    lookUp = e2e[elem*NUM_FACES + dir];
                    if(lookUp!=LOOK_UP_TABLE_DEFAULT)
                    {
                        
                        if(pNodes[lookUp].getLevel() == bLev )
                        {
                            tl = etVec[lookUp];
                            unzipEids.push_back(lookUp);
                            assert(bTime == tl);

                            for(unsigned int node =0; node < nPe; node++)
                                dgV[lookUp *nPe + node] = sV[rk_s-1][lookUp*nPe +  node];

                        }else if(pNodes[lookUp].getLevel() < bLev)
                        {

                            tl = etVec[lookUp];
                            unzipEids.push_back(lookUp);
                            
                            // need to do coarser to finer correction. 
                            for(unsigned int s=1; s <= rk_s; s++ )
                            {
                                cVin[s-1]  = &sV[s-1][lookUp*nPe];
                                cVout[s-1] = cVec.data() + (s-1)*nPe; 
                            }
                            
                            if(tl==bTime)
                                dt = 0;
                            else if(tl == bTime + (1u<<(m_uiLevMax - bLev)))
                                dt = dt_c/2.0;
                            else
                            {
                                std::cout<<"ENUTS sync error at  "<<__func__<<" line : "<<__LINE__<<" finer and coarse time gap violated "<<std::endl;
                                assert(false); // finer block cannot exceed time of the coarset block. 
                                MPI_Abort(comm,0);
                            }

                            Op->Ccf(cVout, (const T**) cVin, cSz, rk_s, dt_c, dt_f, dt, m_uiDof );
                            for(unsigned int node =0; node < nPe; node++)
                                dgV[lookUp *nPe + node] = cVout[rk_s-1][node];


                        }else
                        {
                            assert(pNodes[lookUp].getLevel() == bLev +1);
                            // finner to coaser correction. 

                            pMesh->getFinerFaceNeighbors(elem, dir, (unsigned int *)fchild);
                            assert(etVec[fchild[0]] == etVec[fchild[1]] && etVec[fchild[1]] == etVec[fchild[2]] && etVec[fchild[2]] == etVec[fchild[3]] );
                            unzipEids.push_back(fchild[0]); unzipEids.push_back(fchild[1]); unzipEids.push_back(fchild[2]); unzipEids.push_back(fchild[3]); 

                            tl = etVec[fchild[0]];
                            if(tl==bTime)
                                dt = 0;
                            else if(tl + (1u<<(m_uiLevMax-bLev-1)) == bTime)
                                dt = dt_c/2.0;
                            else
                            {
                                std::cout<<"ENUTS sync error at  "<<__func__<<" line : "<<__LINE__<<" finer and coarse time gap violated "<<std::endl;
                                assert(false); // finer block cannot exceed time of the coarset block. 
                                MPI_Abort(comm,0);
                            }

                            for(unsigned int m=0; m < (NUM_CHILDREN>>1u); m++)
                            {
                                for(unsigned int s=1; s <= rk_s; s++ )
                                {
                                    cVin[s-1]  = &sV[s-1][fchild[m]*nPe];
                                    cVout[s-1] = cVec.data() + (s-1)*nPe; 
                                }

                                Op->Cfc(cVout, (const T**) cVin, cSz, rk_s, dt_c, dt_f, dt, m_uiDof );
                                for(unsigned int node =0; node < nPe; node++)
                                    dgV[fchild[m] *nPe + node] = cVout[rk_s-1][node];
                        
                            }

                        }
                        
                    }


                }

                //std::cout<<"ele back"<<std::endl;

                if(ek==emax)
                {

                    const unsigned int dir = OCT_DIR_FRONT;
                    lookUp = e2e[elem*NUM_FACES + dir];
                    if(lookUp!=LOOK_UP_TABLE_DEFAULT)
                    {
                        
                        if(pNodes[lookUp].getLevel() == bLev )
                        {   
                            tl = etVec[lookUp];
                            unzipEids.push_back(lookUp);
                            assert(bTime == tl);

                            for(unsigned int node =0; node < nPe; node++)
                                dgV[lookUp *nPe + node] = sV[rk_s-1][lookUp*nPe +  node];

                        }else if(pNodes[lookUp].getLevel() < bLev)
                        {

                            tl = etVec[lookUp];
                            unzipEids.push_back(lookUp);
                            
                            // need to do coarser to finer correction. 
                            for(unsigned int s=1; s <= rk_s; s++ )
                            {
                                cVin[s-1]  = &sV[s-1][lookUp*nPe];
                                cVout[s-1] = cVec.data() + (s-1)*nPe; 
                            }
                            
                            if(tl==bTime)
                                dt = 0;
                            else if(tl == bTime + (1u<<(m_uiLevMax - bLev)))
                                dt = dt_c/2.0;
                            else
                            {
                                std::cout<<"ENUTS sync error at  "<<__func__<<" line : "<<__LINE__<<" finer and coarse time gap violated "<<std::endl;
                                assert(false); // finer block cannot exceed time of the coarset block. 
                                MPI_Abort(comm,0);
                            }

                            Op->Ccf(cVout, (const T**) cVin, cSz, rk_s, dt_c, dt_f, dt, m_uiDof );
                            for(unsigned int node =0; node < nPe; node++)
                                dgV[lookUp *nPe + node] = cVout[rk_s-1][node];


                        }else
                        {
                            assert(pNodes[lookUp].getLevel() == bLev +1);
                            // finner to coaser correction. 

                            pMesh->getFinerFaceNeighbors(elem, dir, (unsigned int *)fchild);
                            assert(etVec[fchild[0]] == etVec[fchild[1]] && etVec[fchild[1]] == etVec[fchild[2]] && etVec[fchild[2]] == etVec[fchild[3]] );
                            unzipEids.push_back(fchild[0]); unzipEids.push_back(fchild[1]); unzipEids.push_back(fchild[2]); unzipEids.push_back(fchild[3]); 

                            tl = etVec[fchild[0]];
                            if(tl==bTime)
                                dt = 0;
                            else if(tl + (1u<<(m_uiLevMax-bLev-1)) == bTime)
                                dt = dt_c/2.0;
                            else
                            {
                                std::cout<<"ENUTS sync error at  "<<__func__<<" line : "<<__LINE__<<" finer and coarse time gap violated "<<std::endl;
                                assert(false); // finer block cannot exceed time of the coarset block. 
                                
                            }

                            for(unsigned int m=0; m < (NUM_CHILDREN>>1u); m++)
                            {
                                for(unsigned int s=1; s <= rk_s; s++ )
                                {
                                    cVin[s-1]  = &sV[s-1][fchild[m]*nPe];
                                    cVout[s-1] = cVec.data() + (s-1)*nPe; 
                                }

                                Op->Cfc(cVout, (const T**) cVin, cSz, rk_s, dt_c, dt_f, dt, m_uiDof );
                                for(unsigned int node =0; node < nPe; node++)
                                    dgV[fchild[m] *nPe + node] = cVout[rk_s-1][node];
                        
                            }

                        }
                        
                    }


                }

                //std::cout<<"ele front"<<std::endl;



                
            }

            const std::vector<unsigned int> blkEdge = blkList[m_uiBlkID].getBlk2DiagMap_vec();
            const std::vector<unsigned int> blkVert = blkList[m_uiBlkID].getBlk2VertexMap_vec();

            const unsigned int ele_1d = blkList[m_uiBlkID].getElemSz1D();

            // do for the edge corrections. 
            for(unsigned int dir =0; dir < NUM_EDGES; dir++)
            {
                for(unsigned int i=0; i < ele_1d; i++)
                {
                    if(blkEdge[dir*2*ele_1d + 2*i]!=LOOK_UP_TABLE_DEFAULT)
                    {
                        if(blkEdge[dir*2*ele_1d + 2*i] == blkEdge[dir*2*ele_1d + 2*i+1])
                        {

                            lookUp = blkEdge[dir*2*ele_1d + 2*i];
                            if(lookUp!=LOOK_UP_TABLE_DEFAULT)
                            {
                                
                                if(pNodes[lookUp].getLevel() == bLev )
                                {
                                    tl = etVec[lookUp];
                                    unzipEids.push_back(lookUp);
                                    assert(bTime == tl);

                                    for(unsigned int node =0; node < nPe; node++)
                                        dgV[lookUp *nPe + node] = sV[rk_s-1][lookUp*nPe +  node];

                                }else if(pNodes[lookUp].getLevel() < bLev)
                                {

                                    tl = etVec[lookUp];
                                    unzipEids.push_back(lookUp);
                                    
                                    // need to do coarser to finer correction. 
                                    for(unsigned int s=1; s <= rk_s; s++ )
                                    {
                                        cVin[s-1]  = &sV[s-1][lookUp*nPe];
                                        cVout[s-1] = cVec.data() + (s-1)*nPe; 
                                    }
                                    
                                    if(tl==bTime)
                                        dt = 0;
                                    else if(tl == bTime + (1u<<(m_uiLevMax - bLev)))
                                        dt = dt_c/2.0;
                                    else
                                        assert(false); // finer block cannot exceed time of the coarset block. 

                                    Op->Ccf(cVout, (const T**) cVin, cSz, rk_s, dt_c, dt_f, dt, m_uiDof );
                                    for(unsigned int node =0; node < nPe; node++)
                                        dgV[lookUp *nPe + node] = cVout[rk_s-1][node];


                                }else
                                    assert(false);
                                
                                
                            }


                        }else
                        {
                            lookUp = blkEdge[dir*2*ele_1d + 2*i];
                            echild[0] = blkEdge[dir*2*ele_1d + 2*i];
                            echild[1] = blkEdge[dir*2*ele_1d + 2*i+1];
                            assert(pNodes[lookUp].getLevel() == bLev +1);
                            assert(etVec[echild[0]] == etVec[echild[1]]);

                            unzipEids.push_back(echild[0]);
                            unzipEids.push_back(echild[1]);

                            tl = etVec[echild[0]];
                            if(tl==bTime)
                                dt = 0;
                            else if(tl + (1u<<(m_uiLevMax-bLev-1)) == bTime)
                                dt = dt_c/2.0;
                            else
                            {
                                std::cout<<"ENUTS sync error at  "<<__func__<<" line : "<<__LINE__<<" finer and coarse time gap violated "<<std::endl;
                                assert(false); // finer block cannot exceed time of the coarset block. 
                                MPI_Abort(comm,0);
                            }

                            for(unsigned int m=0; m < (NUM_CHILDREN>>2u); m++)
                            {
                                for(unsigned int s=1; s <= rk_s; s++ )
                                {
                                    cVin[s-1]  = &sV[s-1][echild[m]*nPe];
                                    cVout[s-1] = cVec.data() + (s-1)*nPe; 
                                }

                                Op->Cfc(cVout, (const T**) cVin, cSz, rk_s, dt_c, dt_f, dt, m_uiDof );
                                for(unsigned int node =0; node < nPe; node++)
                                    dgV[echild[m] *nPe + node] = cVout[rk_s-1][node];
                        
                            }

                        }
                    
                    }
                
                }

                //std::cout<<"ele diag"<<dir<<std::endl;
                
            }


            // do the vertex element corrections. 
            for(unsigned int dir =0; dir < NUM_CHILDREN; dir++)
            {
                lookUp = blkVert[dir];

                if(lookUp!=LOOK_UP_TABLE_DEFAULT)
                {
                    if(pNodes[lookUp].getLevel() == bLev )
                    {
                        tl = etVec[lookUp];
                        unzipEids.push_back(lookUp);
                        assert(bTime == tl);
                        

                        for(unsigned int node =0; node < nPe; node++)
                            dgV[lookUp *nPe + node] = sV[rk_s-1][lookUp*nPe +  node];

                    }else if(pNodes[lookUp].getLevel() < bLev)
                    {

                        tl = etVec[lookUp];
                        unzipEids.push_back(lookUp);
                        
                        // need to do coarser to finer correction. 
                        for(unsigned int s=1; s <= rk_s; s++ )
                        {
                            cVin[s-1]  = &sV[s-1][lookUp*nPe];
                            cVout[s-1] = cVec.data() + (s-1)*nPe; 
                        }
                        
                        if(tl==bTime)
                            dt = 0;
                        else if(tl == bTime + (1u<<(m_uiLevMax - bLev)))
                            dt = dt_c/2.0;
                        else
                        {
                            std::cout<<"ENUTS sync error at  "<<__func__<<" line : "<<__LINE__<<" finer and coarse time gap violated "<<std::endl;
                            assert(false); // finer block cannot exceed time of the coarset block. 
                            MPI_Abort(comm,0);
                        }

                        Op->Ccf(cVout, (const T**) cVin, cSz, rk_s, dt_c, dt_f, dt, m_uiDof );
                        for(unsigned int node =0; node < nPe; node++)
                            dgV[lookUp *nPe + node] = cVout[rk_s-1][node];


                    }else
                    {
                        assert(pNodes[lookUp].getLevel() == bLev +1);
                        // finner to coaser correction. 

                        tl = etVec[lookUp];
                        unzipEids.push_back(lookUp);

                        // need to do coarser to finer correction. 
                        for(unsigned int s=1; s <= rk_s; s++ )
                        {
                            cVin[s-1]  = &sV[s-1][lookUp*nPe];
                            cVout[s-1] = cVec.data() + (s-1)*nPe; 
                        }

                        if(tl==bTime)
                            dt = 0;
                        else if(tl + (1u<<(m_uiLevMax-bLev-1)) == bTime)
                            dt = dt_c/2.0;
                        else
                        {
                            std::cout<<"ENUTS sync error at  "<<__func__<<" line : "<<__LINE__<<" finer and coarse time gap violated "<<std::endl;
                            assert(false); // finer block cannot exceed time of the coarset block. 
                            MPI_Abort(comm,0);
                        }
                            

                        Op->Cfc(cVout, (const T**) cVin, cSz, rk_s, dt_c, dt_f, dt, m_uiDof );
                        for(unsigned int node =0; node < nPe; node++)
                            dgV[lookUp *nPe + node] = cVout[rk_s-1][node];

                        
                    }
                    
                }

                //std::cout<<"ele vert "<<dir<<std::endl;


            }

            // do unzip using the DG vector (Now the mesh class supports it. )
            pMesh->DG2CGVec(dgWVec,cgWVec,true,unzipEids.data(),unzipEids.size(),1);
            pMesh->unzip(cgWVec,uzWVec,&m_uiBlkID,1);
            
            // copy from uzWvec to m_uiVec. 
            std::memcpy(vVec,&uzWVec[offset],sizeof(T)*(lx*ly*lz));

        }

        m_uiBVec[blk]._vec[rk_s].mark_synced();
        return;


    }
    
    template<typename T, typename Ctx>
    void ExplicitNUTS<T,Ctx>::update_ele_timestep()
    {
        ot::Mesh* pMesh = m_uiAppCtx->get_mesh();

        if(!(pMesh->isActive()))
            return;

        std::vector<ot::Block> blkList = pMesh->getLocalBlockList();

        for(unsigned int ele = pMesh->getElementPreGhostBegin(); ele < pMesh->getElementPostGhostEnd(); ele ++)
            m_uiEleTime[ele]=0;

        for(unsigned int blk =0; blk < blkList.size(); blk++)
        {
            for(unsigned int ele = blkList[blk].getLocalElementBegin(); ele < blkList[blk].getLocalElementEnd(); ele ++)
                m_uiEleTime[ele] = m_uiBVec[blk]._time;
            
        }

    }
   
   
    template<typename T, typename Ctx>
    void ExplicitNUTS<T,Ctx>::evolve()
    {

        ot::Mesh* pMesh = m_uiAppCtx->get_mesh();
        m_uiTimeInfo = m_uiAppCtx->get_ts_info();
        const double current_t= m_uiTimeInfo._m_uiT;
        double current_t_adv=current_t;
        const double dt_finest = m_uiTimeInfo._m_uiTh;
        const double dt_coarset = (1u<<(m_uiLevMax-m_uiLevMin)) * m_uiTimeInfo._m_uiTh;
        // Assumption: m_uiEvar is time synced accross all the blocks. 

        if(pMesh->isActive())
        {

            const unsigned int rank = pMesh->getMPIRank();
            const unsigned int npes = pMesh->getMPICommSize();
            
            assert (m_uiEVar.GetDof()== m_uiEVecTmp[0].GetDof());
            assert (m_uiEVar.GetDof()== m_uiEVecTmp[1].GetDof());
            
            const unsigned int DOF = m_uiEVar.GetDof();
            MPI_Comm comm = pMesh->getMPICommunicator();

            const ot::TreeNode* const pNodes = pMesh->getAllElements().data();
            m_uiAppCtx->pre_timestep(m_uiEVar);

            const unsigned int  finest_t = 1u<<(m_uiLevMax - m_uiLevMax); // finest  time level.
            const unsigned int coarset_t = 1u<<(m_uiLevMax - m_uiLevMin); // coarset time level. 

            m_uiAppCtx->unzip(m_uiEVar,m_uiEvarUzip ,m_uiAppCtx->get_async_batch_sz());  // unzip m_uiEVec to m_uiEVarUnzip
            std::vector<ot::Block> blkList = pMesh->getLocalBlockList();

            // initialize the block async vectors
            for(unsigned int blk =0; blk < blkList.size(); blk++)
            {
                m_uiBVec[blk]._time = 0;
                m_uiBVec[blk]._rks  = 0;
                m_uiBVec[blk]._vec[0].copyFromUnzip(pMesh,m_uiEvarUzip.GetVecArray(), true, DOF);
                m_uiBVec[blk]._vec[0].mark_synced();

                for(unsigned int s=1;  s <= m_uiNumStages; s++)
                    m_uiBVec[blk]._vec[s].mark_unsynced();
            }

            this->update_ele_timestep();
            
            pMesh->readFromGhostBeginElementVec(m_uiEleTime.data(),1);
            pMesh->readFromGhostEndElementVec(m_uiEleTime.data(),1);

            for(unsigned int pt=0; pt<coarset_t; pt ++)
            {
                std::cout<<"[ENUTS] : pt: "<<pt<<" \n";
                for(unsigned int rk=1; rk <= m_uiNumStages; rk++ )
                {
                    const unsigned int BLK_S = rk;
                    for(unsigned int blk =0; blk < blkList.size(); blk++)
                    {
                        const unsigned int BLK_T = m_uiBVec[blk]._time;
                        const unsigned int NN   =  m_uiBVec[blk]._vec[0].getSz();
                        const unsigned int bLev =  pNodes[blkList[blk].getLocalElementBegin()].getLevel();
                        const unsigned int BLK_DT = 1u<<(m_uiLevMax - bLev);
                        if( pt% BLK_DT !=0 )
                            continue;

                        
                        //std::cout<<"[NUTS]: pt: "<<pt<<" blk: "<<blk<<" rk : "<<rk<<" step size: "<<BLK_DT<<std::endl;

                        assert(m_uiBVec[blk]._vec[rk-1].isSynced());

                        m_uiBVec[blk]._vec[m_uiNumStages+1].vecCopy(m_uiBVec[blk]._vec[0]); // copy the in vector. 
                        T* out_ptr = (T*)m_uiBVec[blk]._vec[m_uiNumStages+1].data();
                        
                        for(unsigned int s=1;  s < BLK_S; s++)
                        {
                            const T* stage_ptr = m_uiBVec[blk]._vec[s].data();
                            for(unsigned int n =0; n < DOF*NN; n++)
                                out_ptr[n] += m_uiAij[ (BLK_S-1) * m_uiNumStages + (s-1)]* stage_ptr[n];
                        }

                        m_uiBVec[blk]._vec[BLK_S].computeVec(m_uiAppCtx, m_uiBVec[blk]._vec[m_uiNumStages+1], current_t + dt_finest*BLK_T);
                        m_uiBVec[blk]._vec[BLK_S].mark_unsynced();
                        m_uiBVec[blk]._vec[BLK_S].zipDG(pMesh,m_uiStVec[BLK_S-1].GetVecArray(),DOF);
                        m_uiBVec[blk]._rks = rk;

                    }

                    // do the DG vec ghost sync. 
                    pMesh->readFromGhostBeginEleDGVec(m_uiStVec[BLK_S-1].GetVecArray(),DOF);
                    pMesh->readFromGhostEndEleDGVec(m_uiStVec[BLK_S-1].GetVecArray(),DOF);

                    for(unsigned int blk =0; blk < blkList.size(); blk++)
                    {
                        const unsigned int BLK_T = m_uiBVec[blk]._time;
                        const unsigned int NN   =  m_uiBVec[blk]._vec[0].getSz();
                        const unsigned int bLev =  pNodes[blkList[blk].getLocalElementBegin()].getLevel();
                        const unsigned int BLK_DT = 1u<<(m_uiLevMax - bLev);
                        if( pt% BLK_DT !=0 )
                            continue;

                        sync_blk_timestep(blk,BLK_S);
                        m_uiBVec[blk]._vec[BLK_S].mark_synced();
                        //std::cout<<"[NUTS]: pt :  "<<pt<<" blk : "<<blk<<" rk : "<<rk<<" step size: "<<BLK_DT<<" is synced: "<<m_uiBVec[blk]._vec[BLK_S].isSynced()<<std::endl;

                    }
                        
                    
                }

                // compute the time step vector and increment time. 
                for(unsigned int blk =0; blk < blkList.size(); blk++)
                {

                    const unsigned int BLK_T = m_uiBVec[blk]._time;
                    const unsigned int NN   =  m_uiBVec[blk]._vec[0].getSz();
                    const unsigned int bLev =  pNodes[blkList[blk].getLocalElementBegin()].getLevel();
                    const unsigned int BLK_DT = 1u<<(m_uiLevMax - bLev);
                    if( pt% BLK_DT !=0 )
                        continue;

                    T* out_ptr = (T*)m_uiBVec[blk]._vec[0].data();
                    for(unsigned int s=1;  s <= m_uiNumStages; s++)
                    {
                        assert(m_uiBVec[blk]._vec[s].isSynced());
                        const T* stage_ptr = m_uiBVec[blk]._vec[s].data();
                        for(unsigned int n =0; n < DOF*NN; n++)
                            out_ptr[n] += m_uiBi[(s-1)]*stage_ptr[n];
                    }

                    m_uiBVec[blk]._time += 1u<<(m_uiLevMax -bLev); 
                    m_uiBVec[blk]._rks=0;
                    m_uiBVec[blk]._vec[0].mark_synced();

                    for(unsigned int s=1;  s <= m_uiNumStages; s++)
                        m_uiBVec[blk]._vec[s].mark_unsynced();

                }


                update_ele_timestep();
                pMesh->readFromGhostBeginElementVec(m_uiEleTime.data(),1);
                pMesh->readFromGhostEndElementVec(m_uiEleTime.data(),1);

            
            }


        }

        pMesh->waitAll();
        m_uiAppCtx->increment_ts_info();
        m_uiTimeInfo = m_uiAppCtx->get_ts_info();
        pMesh->waitAll();


    }


}


