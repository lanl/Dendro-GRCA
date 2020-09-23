/**
 * @file ts.h
 * @author Milinda Fernando
 * @brief generic time integrator class for Dendro. 
 * @version 0.1
 * @date 2019-10-18
 * 
 * School of Computing, University of Utah
 * @copyright Copyright (c) 2019
 * 
 */

#pragma once
#include "dendro.h"
#include "mesh.h"
#include "ctx.h"
#include "dvec.h"
#include "ts.h"

namespace ts
{
    /**time stepper type
     * UTS uniform time stepper. 
     * UTS_ADAP: uniform over the grid but time step size changes over time. 
     * NUTS: spatially adaptive time stepping. 
     * NUTS_ADAP: NUTS where the smallest time step varies in time. 
     * 
    */
    enum TimeStepperType {UTS=0 , UTS_ADAP ,NUTS, NUTS_ADAP};

    /**
     * @brief ETS Flags (currently not used)
     */
    enum ETSFlags {FROM_T0 =0 , CHECKPT, CURR_STEP, CURR_TIME };

    /**
     * @brief General explicit time stepper class for Dendro-5.0
     * @tparam T 
     */
    template<typename T, typename Ctx>
    class ETS
    {
        protected: 

            /** @brief Application context. */
            Ctx* m_uiAppCtx;
            
            /**@brief: Time stepper type*/
            ETSType m_uiType;

            /**@brief: CFL factor*/
            DendroScalar m_uiCFL;

            /**@brief: time integrator coefficients for solution u*/
            DendroScalar* m_uiAij = NULL;

            /**@brief: time integrator coefficients for time */
            DendroScalar* m_uiBi = NULL;

            /**@brief: time integrator weights*/
            DendroScalar* m_uiCi = NULL;

            /**@brief: number of stages*/
            unsigned int m_uiNumStages;

            /**@brief: time step info*/
            TSInfo m_uiTimeInfo;

            /**@brief: evolution variables*/
            DVec m_uiEVar;
            
            /**@brief: stage vector*/
            std::vector<DVec> m_uiStVec;

            /**@brief: evolution temp vector*/
            DVec m_uiEVecTmp[2];
            

        private:
            /**
             * @brief Allocates internal variables for the time stepper. 
             * @return int 
             */
            virtual int allocate_internal_vars();

            /**@brief: Deallocate internal variables. */
            virtual int deallocate_internal_vars();
            

        public: 
            /**
             * @brief Construct a new ETS object
             * @param pMesh : underlying mesh data structure. 
             */
            ETS(Ctx* appCtx);
            
            /**
             * @brief Destroy the ETS object
             */
            ~ETS();

            /**
             * @brief Set the ets coefficients
             * $ k_i = f(u^{n} + \sum_{j=1}^{i-1} a_{i,j}*dt , t^{n} + b_i*th)$
             * $u^{n+1} = u^{n} + \sum_{m=1}^{num\_stages} k_m $
             * @param aij : time integrator coefficients for the f_rhs first term
             * @param bi : time integrator coefficients for the f_rhs second term
             * @param ci : time integrator weights
             * @param num_stages : number of stages 
             * @return int : is success return zero. 
             */
            int set_ets_coefficients(DendroScalar * aij , DendroScalar * bi ,DendroScalar * ci,unsigned int num_stages);

            /**
             * @brief : sets default ETS time integrator.  
             * @param [in] type: time integrator type. 
             */
            int set_ets_coefficients(ETSType type);

            /**
             * @brief Set the evolve vars for the ETS time stepper. 
             * 
             * @param eVar : evolution variables, multiple evolution variables should be added as one vector with multiple dof. 
             * @return int 
             */
            int set_evolve_vars(DVec eVar);
            
            /**@brief: initialize the ETS solver*/
            void init();

            /**@brief: returns the current time step*/
            inline DendroIntL curr_step() {return m_uiTimeInfo._m_uiStep;};

            /**@brief: returns the current time*/
            inline DendroScalar curr_time() {return m_uiTimeInfo._m_uiT;};

            /**@brief: */
            inline DendroScalar ts_size() const { return m_uiTimeInfo._m_uiTh;}

            /**@brief: returns the mesh is active. */
            inline bool is_active() const {return m_uiAppCtx->get_mesh()->isActive();}
            
            /**@brief: returns the active rank*/
            unsigned int get_active_rank() const ;
            
            /**@brief: returns the active npes*/
            unsigned int get_active_npes() const;
            
            /**@brief: returns the global rank*/
            unsigned int get_global_rank() const;

            /**@brief: return the global npes*/
            unsigned int get_global_npes() const;

            /**@brief: return the active communicator*/
            MPI_Comm get_active_comm () const;

            /**@brief: return the global communicator. */
            MPI_Comm get_global_comm() const;

            /**@brief: returns the underlying mesh data structure. */
            const ot::Mesh* get_mesh() const  { return m_uiAppCtx->get_mesh();} 

            /**@brief: returns the evolution variables. */ 
            inline DVec get_evolve_vars() const { return m_uiEVar; }

            /**@brief: returns the time step info*/
            inline TSInfo get_timestep_info() const {return m_uiTimeInfo;} 

            /**@brief: perform synchronizations with correct variable allocations for the new mesh: should be called after remeshing.  */
            virtual int sync_with_mesh();
            
            /**@brief: advance to next time step*/
            void evolve();
            
    };


    template<typename T, typename  Ctx>
    ETS<T , Ctx>::ETS(Ctx* appCtx)
    {   
        
        m_uiAppCtx = appCtx;
        m_uiAij = NULL;
        m_uiBi = NULL;
        m_uiCi = NULL;
        m_uiNumStages  = 0;
        m_uiTimeInfo = appCtx->get_ts_info();

        m_uiEVar = m_uiAppCtx->get_evolution_vars();

    }

    template<typename T, typename  Ctx>
    ETS<T , Ctx>::~ETS()
    {   
        return ;
    }

    template<typename T, typename  Ctx>
    int ETS<T , Ctx>::set_ets_coefficients(DendroScalar * aij , DendroScalar * bi ,DendroScalar * ci,unsigned int num_stages)
    {
        m_uiAij = aij;
        m_uiBi  = bi;
        m_uiCi  = ci;

        m_uiNumStages = num_stages;
        return 0;

    }

    template<typename T, typename Ctx>
    int ETS<T,Ctx> :: allocate_internal_vars()
    {
        m_uiStVec.resize(m_uiNumStages);
        for(unsigned int i=0; i < m_uiNumStages; i++)
            m_uiStVec[i].VecCreate(m_uiAppCtx->get_mesh(), m_uiEVar.IsGhosted() , m_uiEVar.IsUnzip(), m_uiEVar.IsElemental() , m_uiEVar.GetDof());
        
        m_uiEVecTmp[0].VecCreate(m_uiAppCtx->get_mesh(), m_uiEVar.IsGhosted() , m_uiEVar.IsUnzip(), m_uiEVar.IsElemental() , m_uiEVar.GetDof());

        m_uiEVecTmp[1].VecCreate(m_uiAppCtx->get_mesh(), m_uiEVar.IsGhosted() , m_uiEVar.IsUnzip(), m_uiEVar.IsElemental() , m_uiEVar.GetDof());
        return 0;
    }
    
    template<typename T, typename Ctx>
    int ETS<T,Ctx> :: deallocate_internal_vars()
    {
        for(unsigned int i=0; i < m_uiNumStages; i++)
            m_uiStVec[i].VecDestroy();

        m_uiStVec.clear();

        m_uiEVecTmp[0].VecDestroy();
        m_uiEVecTmp[1].VecDestroy();
        
        return 0; 
    }


    template<typename T, typename  Ctx>
    int ETS<T,Ctx>::set_ets_coefficients(ETSType type)
    {
        m_uiType = type ;

        if(type == ETSType::RK3)
        {
            m_uiNumStages = 3;

            static const DendroScalar ETS_C[] = {1.0/6.0,1.0/6.0,2.0/3.0};
            static const DendroScalar ETS_T[] = {0.0,1.0,1.0/2.0};
            static const DendroScalar ETS_U[] = { 0.0     ,0.0      , 0.0, 
                                                  1.0     ,0.0      , 0.0,
                                                  1.0/4.0 ,1.0/4.0  , 0.0 };


            m_uiCi  =  (DendroScalar*) ETS_T;
            m_uiBi  =  (DendroScalar*)ETS_C;
            m_uiAij =  (DendroScalar*)ETS_U;
        
        }else if( type == ETSType::RK4)
        {
            m_uiNumStages = 4; 
            
            static const DendroScalar ETS_C[] = {1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};
            static const DendroScalar ETS_T[] = {0,1.0/2.0,1.0/2.0,1.0};
            static const DendroScalar ETS_U[] = { 0.0   , 0.0       , 0.0   , 0.0,
                                                1.0/2.0 , 0.0       , 0.0   , 0.0,
                                                0.0     , 1.0/2.0   , 0.0   , 0.0,
                                                0.0     , 0.0       , 1.0   , 0.0 };

            m_uiCi  =  (DendroScalar*) ETS_T;
            m_uiBi  =  (DendroScalar*)ETS_C;
            m_uiAij =  (DendroScalar*)ETS_U;

        }else if ( type == ETSType::RK5)
        {

            return -1;
            // these coefficients looks wrong. - Milinda. (need to fix those and enable the rk5 method. )
            m_uiNumStages = 5;
            static const DendroScalar ETS_C[] = {7.0/90.0, 0, 32/90.0,12.0/90.0, 32.0/90.0, 7.0/90.0};
            static const DendroScalar ETS_T[] = {0,1.0/4.0, 1.0/4.0, 1.0/2.0, 3.0/4.0, 1.0};
            static const DendroScalar ETS_U[] = {   0.0       ,  0.0        , 0.0       , 0.0        , 0.0,
                                                    1.0/8.0   ,  1.0/8.0    , 0.0       , 0.0        , 0.0,
                                                    0.0       , -1.0/2.0    , 1.0       , 0.0        , 0.0,
                                                    3.0/16.0  ,  0.0        , 0.0       , 9.0/16.0   , 0.0,
                                                    -3.0/7.0  ,  2.0/7.0    , 12.0/7.0  , -12.0/7.0  , 8.0/7.0};
            
            m_uiCi  =  (DendroScalar*) ETS_T;
            m_uiBi  =  (DendroScalar*)ETS_C;
            m_uiAij =  (DendroScalar*)ETS_U;

        }else
            return -1;

        

        return 0;
    }

    template<typename T, typename  Ctx>
    void ETS<T,Ctx>::init()
    {
        m_uiAppCtx->initialize();
        m_uiTimeInfo = m_uiAppCtx->get_ts_info();
        allocate_internal_vars();

    }

    template<typename T, typename Ctx>
    int ETS<T,Ctx>::set_evolve_vars(DVec eVars)
    {
        m_uiEVar = eVars;
        return 0; 
    }

    template<typename T, typename  Ctx>
    void ETS<T, Ctx>::evolve()
    {
        const ot::Mesh* pMesh = m_uiAppCtx->get_mesh();
        m_uiTimeInfo = m_uiAppCtx->get_ts_info();
        const double current_t= m_uiTimeInfo._m_uiT;
        double current_t_adv=current_t;
        const double dt = m_uiTimeInfo._m_uiTh;

        if(pMesh->isActive())
        {
            
            int rank =pMesh->getMPIRank();

            const unsigned int nodeLocalBegin=pMesh->getNodeLocalBegin();
            const unsigned int nodeLocalEnd=pMesh->getNodeLocalEnd();

            const std::vector<ot::Block> blkList=pMesh->getLocalBlockList();
            unsigned int offset;
            double ptmin[3], ptmax[3];
            unsigned int sz[3];
            unsigned int bflag;
            double dx,dy,dz;
            
            m_uiAppCtx->pre_timestep(m_uiEVar);
            for(int stage=0; stage< m_uiNumStages ; stage++)
            {
                //m_uiEVecTmp[0].VecCopy(m_uiEVar,true);
                m_uiEVecTmp[0].VecFMA(m_uiAppCtx->get_mesh(),m_uiEVar,0,1.0,true); 
                
                for(int p = 0 ; p < stage;  p++ )
                    m_uiEVecTmp[0].VecFMA(m_uiAppCtx->get_mesh(), m_uiStVec[stage], 1, m_uiAij[(stage)*m_uiNumStages + p]*dt,true);

                current_t_adv=current_t+m_uiCi[stage] * dt;

                m_uiAppCtx -> pre_stage(m_uiStVec[stage]);
                m_uiAppCtx -> rhs(&m_uiEVecTmp[0], &m_uiStVec[stage], 1, current_t_adv);
                m_uiAppCtx -> post_stage(m_uiStVec[stage]);
                

            }

            for(unsigned int k=0; k< m_uiNumStages; k++)
                m_uiEVar.VecFMA(m_uiAppCtx->get_mesh(), m_uiStVec[k], 1, m_uiBi[k]*dt, true);
            
            m_uiAppCtx->post_timestep(m_uiEVar);
            

        }

        m_uiAppCtx->increment_ts_info();
        m_uiTimeInfo = m_uiAppCtx->get_ts_info();
        pMesh->waitAll();

        
    }
    
    template<typename T, typename Ctx>
    unsigned int ETS<T,Ctx>::get_active_rank() const {
        if(is_active())
            return m_uiAppCtx->get_mesh()->getMPIRank();
        else
            return get_global_rank();
    }

    template<typename T, typename Ctx>
    unsigned int ETS<T,Ctx>::get_active_npes() const {
        if(is_active())
            return m_uiAppCtx->get_mesh()->getMPICommSize();
        else
            return get_global_npes();
    }


    template<typename T, typename Ctx>
    unsigned int ETS<T,Ctx>::get_global_rank() const {
        return m_uiAppCtx->get_mesh()->getMPIRankGlobal();
    }

    template<typename T, typename Ctx>
    unsigned int ETS<T,Ctx>::get_global_npes() const {
        return m_uiAppCtx->get_mesh()->getMPICommSizeGlobal();
    }

    template<typename T, typename Ctx>
    MPI_Comm ETS<T,Ctx>:: get_active_comm() const 
    {
        if(is_active())
            return m_uiAppCtx->get_mesh()->getMPICommunicator();
        else
            return MPI_COMM_NULL;
    }

    template<typename T, typename Ctx>
    MPI_Comm ETS<T,Ctx>:: get_global_comm() const 
    {
        return m_uiAppCtx->get_mesh()->getMPIGlobalCommunicator();
    }

    template <typename T, typename Ctx >
    int ETS<T,Ctx>::sync_with_mesh()
    {
        if(m_uiAppCtx -> is_ets_synced())
            return 0;
        

        deallocate_internal_vars();
        allocate_internal_vars();
        m_uiEVar = m_uiAppCtx->get_evolution_vars();
        m_uiAppCtx->set_ets_synced(true);

        return 0;
    }


}// end of namespace ts
