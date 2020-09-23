/**
 * @file ctx.h
 * @author Milinda Fernando 
 * @brief : Application context which can be used for time stepping. 
 * @version 0.1
 * @date 2019-12-17
 * 
 * School of Computing, University of Utah. 
 * @copyright Copyright (c) 2019
 * 
 */
 #pragma once
 #include "dendro.h"
 #include "mesh.h"
 #include <vector>
 #include <iostream>
 #include "dvec.h"
 #include "ts.h"
 #include <functional>

namespace ts
{
    /**@brief: different variable types. */
    enum CTXVType {EVOLUTION=0, CONSTRAINT, PRIMITIVE};

    template<typename T, typename I>
    class Ctx{
        
        protected: 

            /**@brief : Mesh object  */
            ot::Mesh * m_uiMesh; 

            /**@brief: variable list related to the application*/
            std::vector< ot::DVector<T,I> > m_uiEvolutionVar;

            /**@brief: variable list related to the application*/
            std::vector< ot::DVector<T,I> > m_uiConstrainedVar;

            /**@brief: variable list related to the application*/
            std::vector< ot::DVector<T,I> > m_uiPrimitiveVar;

            /**@brief: Evolution unzip variable*/
            std::vector< ot::DVector<T,I> > m_uiEvolutionUnzipVar;
            
            /**@brief: Constraint unzip variables*/
            std::vector< ot::DVector<T,I> > m_uiConstraintUnzipVar;

            /**@brief: primitive unzip variables. */
            std::vector< ot::DVector<T,I> > m_uiPrimitiveUnzipVar;


            /**@brief: Total memory allocated for time stepper in B */
            DendroIntL m_uiMemAlloc;

            /**@brief: Total memory deallocated for time stepper B*/
            DendroIntL m_uiMemDeAlloc;

            /**@brief: time stepper info*/
            ts::TSInfo m_uiTinfo;

            /**@brief: element order for mesh generation*/
            unsigned int m_uiElementOrder;

            /** @brief: min point of the compute domain. */
            Point m_uiMinPt;
            
            /**@brief: max point of the compute domain */
            Point m_uiMaxPt;
            
            /**@brief: is true indicates that the ETS is synced with the current mesh. */
            bool m_uiIsETSSynced = true;

            
        public: 
        
            /**@brief: default constructor*/  
            Ctx(){};

            /**@brief: default destructor*/
            ~Ctx(){};

            /**@brief: returns the const ot::Mesh. */
            inline ot::Mesh* get_mesh() const {return m_uiMesh; } 
            
            /**@brief: returns the time stamp info, related to ets*/
            inline ts::TSInfo get_ts_info() const {return m_uiTinfo; }
            
            /**@breif: returns the ETS synced status*/
            inline bool is_ets_synced() const { return m_uiIsETSSynced;} 

            /**@brief:sets the ets sync status*/
            inline void set_ets_synced(bool s) {m_uiIsETSSynced = s; }

            /**@brief: initial solution*/
            virtual int initialize() {return 0;};
            
            /**@brief: right hand side computation v= F(u,t);
             * @param [in] in : input (u vector) for the evolution variables. 
             * @param [out] out: output (v vector) computed rhs. 
             * @param [in] sz: number of in and out vectors (dof)
             * @param [in] time : current time at the evolution. 
            */
            virtual int rhs(ot::DVector<T,I>* in , ot::DVector<T,I>* out, unsigned int sz , T time) {return 0;};

            /**
             * @brief block wise RHS. 
             * 
             * @param in : input vector (unzip version)
             * @param out : output vector (unzip version)
             * @param blkIDs : local block ids where the rhs is computed.  
             * @param sz : size of the block ids
             * @param blk_time : block time  corresponding to the block ids. 
             * @return int 
             */
            virtual int rhs_blkwise(ot::DVector<T,I> in , ot::DVector<T,I> out, const unsigned int* const blkIDs, unsigned int numIds, T*  blk_time) const {return 0;}

            /**
             * @brief compute the block for the rhs. 
             * 
             * @param in 
             * @param out 
             * @param local_blk_id 
             * @param blk_time 
             * @return int 
             */
            virtual int rhs_blk(const T* in, T* out, unsigned int dof, unsigned int local_blk_id, T  blk_time) const {return 0;};

            /**@brief: compute constraints. */
            virtual int compute_constraints() {return 0;}

            /**@brief: compute primitive variables. */
            virtual int compute_primitives() {return 0;}
            
            /**@brief: function execute before each stage
             * @param sIn: stage var in. 
            */
            virtual int pre_stage(ot::DVector<T,I> sIn)  {return 0;}; 

            /**@brief: function execute after each stage
             * @param sIn: stage var in. 
            */
            virtual int post_stage(ot::DVector<T,I> sIn)  {return 0;};

            /**@brief: function execute before each step*/
            virtual int pre_timestep(ot::DVector<T,I> sIn) {return 0;}; 

            /**@brief: function execute after each step*/
            virtual int post_timestep(ot::DVector<T,I> sIn) {return 0;};

            
            /**@brief: function execute after each step*/
            virtual bool is_remesh() {return false;};

            /**
             * @brief performs remesh if the remesh flags are true. 
             * 
             * @param grain_sz : grain size for the remesh 
             * @param ld_tol : load imbalance tolerance. 
             * @param sf_k : splitter fix value 
             * @param transferEvolution : if true transform the evolution variales. 
             * @param transferConstraint : if true transform the constraint variables. 
             * @param transferPrimitive : if true transform the primitive variables. 
             * @return int : return 0 if success. 
             */
            virtual int remesh(unsigned int grain_sz = DENDRO_DEFAULT_GRAIN_SZ, double ld_tol = DENDRO_DEFAULT_LB_TOL, unsigned int sf_k = DENDRO_DEFAULT_SF_K ,bool transferEvolution = true, bool transferConstraint = false, bool transferPrimitive = false ); 

            /**
             * @brief performs intergrid transfer for a given mesh variable. 
             * 
             * @param newMesh : Mesh object the variables needed to be transfered.
             * @param transferEvolution : if true transform the evolution variables.  
             * @param transferConstraint : if true transform the constraint variales. 
             * @param transferPrimitive : if true transform the primitive variables. 
             * @return int 
             */
            virtual int grid_transfer(ot::Mesh* newMesh, bool transferEvolution = true, bool transferConstraint = false, bool transferPrimitive = false);

            /**@brief: write vtu files and output related stuff. */
            virtual int write_vtu() {return 0;};

            /**@brief: writes checkpoint*/
            virtual int write_checkpt() {return 0;};

            /**@brief: restore from check point*/
            virtual int restore_checkpt() {return 0;};

            /**@brief: should be called for free up the contex memory. */
            virtual int finalize() {return 0;};

            /**@brief: Add variables to the time stepper*/
            virtual ot::DVector<T,I> create_vec(CTXVType type, bool isGhosted = false, bool isUnzip =false, bool isElemental = false, unsigned int dof=1);

            /**@brief: destroy a vector 
             * @param vec: DVector object 
             * @param type: DVector type. 
            */
            virtual void destroy_vec(ot::DVector<T,I>& vec, CTXVType type);

            /**@brief: destroy vector 
             * @param vec: DVector object
            */
            virtual void destroy_vec(ot::DVector<T,I>& vec);

            /**
             * @brief performs unzip operation, 
             * 
             * @param in : input zip vector. 
             * @param out : output unzip vector. 
             * @param async_k : async communicator. 
             */
            virtual void unzip(ot::DVector<T,I>& in , ot::DVector<T,I>& out, unsigned int async_k = 1);

            /**
             * @brief performs zip operation
             * 
             * @param in : unzip vector
             * @param out : zip vector. 
             */
            virtual void zip(ot::DVector<T,I>& in , ot::DVector<T,I>& out, unsigned int async_k = 1);

            /**@brief: pack and returns the evolution variables to one DVector*/
            virtual ot::DVector<T,I> get_evolution_vars() { ot::DVector<T,I> tmp; return tmp; };

            /**@brief: pack and returns the constraint variables to one DVector*/
            virtual ot::DVector<T,I> get_constraint_vars(){ ot::DVector<T,I> tmp; return tmp; };

            /**@brief: pack and returns the primitive variables to one DVector*/
            virtual ot::DVector<T,I> get_primitive_vars(){ ot::DVector<T,I> tmp; return tmp; };

            /**@brief: updates the time step information. */
            virtual int increment_ts_info();

            /**@brief: updates the application variables from the variable list. */
            virtual int update_app_vars() {return 0;};

            /**@brief: prints any messages to the terminal output. */
            virtual int terminal_output() {return 0;};

            /**@brief: returns the async communication batch size. */
            virtual unsigned int get_async_batch_sz() {return 1;};

            
    };

    template<typename T, typename I>
    ot::DVector<T,I> Ctx<T,I>::create_vec(CTXVType type, bool isGhosted, bool isUnzip, bool isElemental, unsigned int dof)
    {
        ot::DVector<T,I> tmp;
        tmp.VecCreate(m_uiMesh,isGhosted,isUnzip,isElemental,dof);

        if(type == CTXVType::EVOLUTION)
            (isUnzip) ? m_uiEvolutionUnzipVar.push_back(tmp) : m_uiEvolutionVar.push_back(tmp);
        else if(type == CTXVType::CONSTRAINT)
            (isUnzip) ? m_uiConstraintUnzipVar.push_back(tmp) : m_uiConstrainedVar.push_back(tmp);
        else if(type == CTXVType::PRIMITIVE)
            (isUnzip) ? m_uiPrimitiveUnzipVar.push_back(tmp) : m_uiPrimitiveVar.push_back(tmp);
        
        return tmp;
    }

    template<typename T,typename I>
    void Ctx<T,I>::destroy_vec(ot::DVector<T,I>& vec, CTXVType type)
    {
        int index=-1;
        if(type == CTXVType::EVOLUTION)
        {
            for(unsigned int i=0; i< m_uiEvolutionVar.size(); i++)
            {
                if(m_uiEvolutionVar[i] == vec)
                {
                    index =i;
                    break;
                }
            }

            if(index >=0)
                m_uiEvolutionVar.erase( m_uiEvolutionVar.begin() + index);

            
        }
        else if(type == CTXVType::CONSTRAINT)
        {
            for(unsigned int i=0; i< m_uiConstrainedVar.size(); i++)
            {
                if(m_uiConstrainedVar[i] == vec)
                {
                    index =i;
                    break;
                }
            }

            if(index >=0)
                m_uiConstrainedVar.erase(m_uiConstrainedVar.begin()+index);
        }
        else if(type == CTXVType::PRIMITIVE)
        {
            for(unsigned int i=0; i< m_uiPrimitiveVar.size(); i++)
            {
                if(m_uiPrimitiveVar[i] == vec)
                {
                    index =i;
                    break;
                }
            }

            if(index >=0)
                m_uiPrimitiveVar.erase(m_uiPrimitiveVar.begin()+index);

        }
        
    }


    template<typename T,typename I>
    void Ctx<T,I>::destroy_vec(ot::DVector<T,I>& vec)
    {
        int index=-1;
        
        for(unsigned int i=0; i< m_uiEvolutionVar.size(); i++)
        {
            if(m_uiEvolutionVar[i] == vec)
            {
                index =i;
                break;
            }
        }

        if(index >=0)
            m_uiEvolutionVar.erase(m_uiEvolutionVar.begin() +  index);

        index =-1;
        for(unsigned int i=0; i< m_uiConstrainedVar.size(); i++)
        {
            if(m_uiConstrainedVar[i] == vec)
            {
                index =i;
                break;
            }
        }

        if(index >=0)
            m_uiConstrainedVar.erase(m_uiConstrainedVar.begin() + index);

        index=-1;
        for(unsigned int i=0; i< m_uiPrimitiveVar.size(); i++)
        {
            if(m_uiPrimitiveVar[i] == vec)
            {
                index =i;
                break;
            }
        }

        if(index >=0)
            m_uiPrimitiveVar.erase( m_uiPrimitiveVar.begin() +  index);

        
    }
    
    template<typename T, typename I>
    void Ctx<T,I>::unzip(ot::DVector<T,I>& in , ot::DVector<T,I>& out, unsigned int async_k)
    {
        
        assert( (in.IsUnzip() == false) && (in.GetDof()== out.GetDof()) && (out.IsUnzip()==true) && (in.IsGhosted()==true) && async_k <= in.GetDof());
        
        const unsigned int dof = in.GetDof();
        T* in_ptr  = in.GetVecArray();
        T* out_ptr = out.GetVecArray();

        //std::cout<<" dof: "<<in.GetDof()<<std::endl;

        const unsigned int sz_per_dof_zip  = in.GetSizePerDof();
        const unsigned int sz_per_dof_uzip = out.GetSizePerDof();
        assert(sz_per_dof_uzip == m_uiMesh->getDegOfFreedomUnZip());
        assert(sz_per_dof_zip == m_uiMesh->getDegOfFreedom());

        for(unsigned int i=0 ; i < async_k; i++)
        {
            for(unsigned int j=((i*dof)/async_k); j < (((i+1)*dof)/async_k); j++)   
                m_uiMesh->readFromGhostBegin(in_ptr + j*sz_per_dof_zip, 1);

            
            if(i>0)
            {
                for(unsigned int j= (((i-1)*dof)/async_k); j < ((i*dof)/async_k); j++ )
                    m_uiMesh->unzip(in_ptr + j*sz_per_dof_zip, out_ptr + j*sz_per_dof_uzip);

            }

            for(unsigned int j=((i*dof)/async_k); j < (((i+1)*dof)/async_k); j++)
                m_uiMesh->readFromGhostEnd(in_ptr + j*sz_per_dof_zip, 1);

        }

        for(unsigned int j= (((async_k-1)*dof)/async_k); j < ((async_k*dof)/async_k); j++ )
            m_uiMesh->unzip(in_ptr + j*sz_per_dof_zip, out_ptr + j*sz_per_dof_uzip);

        
    }


    template<typename T, typename I>
    void Ctx<T,I>::zip(ot::DVector<T,I>& in , ot::DVector<T,I>& out, unsigned int async_k)
    {
        assert( (in.IsUnzip() == true) && (in.GetDof()== out.GetDof()) && (out.IsUnzip()==false) && (out.IsGhosted()==true));
        const unsigned int dof = in.GetDof();
        const unsigned int sz_per_dof_uzip = in.GetSizePerDof();
        const unsigned int sz_per_dof_zip = out.GetSizePerDof();

        T* in_ptr = in.GetVecArray();
        T* out_ptr = out.GetVecArray();

        assert(sz_per_dof_uzip == m_uiMesh->getDegOfFreedomUnZip());
        assert(sz_per_dof_zip == m_uiMesh->getDegOfFreedom());

        for(unsigned int i=0 ; i < async_k; i++)
        {
            
            for(unsigned int j=((i*dof)/async_k); j < (((i+1)*dof)/async_k); j++)
            {
                m_uiMesh->zip(in_ptr + j*sz_per_dof_uzip,out_ptr + j*sz_per_dof_zip);
            }
               
            if(i>0)
            {
                for(unsigned int j= (((i-1)*dof)/async_k); j < ((i*dof)/async_k); j++ )
                    m_uiMesh->readFromGhostBegin(out_ptr + j*sz_per_dof_zip,1);
            }
        }

        for(unsigned int j= (((async_k-1)*dof)/async_k); j < ((async_k*dof)/async_k); j++ )
            m_uiMesh->readFromGhostBegin(out_ptr + j*sz_per_dof_zip,1);

        for(unsigned int j=0; j < dof; j++)
            m_uiMesh->readFromGhostEnd(out_ptr + j*sz_per_dof_zip,1);
            
        
        return;

    }

    
    template<typename T,typename I>
    int Ctx<T,I>::remesh(unsigned int grain_sz, double ld_tol, unsigned int sf_k, bool transferEvolution, bool transferConstraint, bool transferPrimitive)
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
        
        ot::Mesh* newMesh=m_uiMesh->ReMesh(grain_sz,ld_tol,sf_k);
        
        DendroIntL oldElements=m_uiMesh->getNumLocalMeshElements();
        DendroIntL newElements=newMesh->getNumLocalMeshElements();

        DendroIntL oldElements_g, newElements_g;

        par::Mpi_Reduce(&oldElements,&oldElements_g,1,MPI_SUM,0,m_uiMesh->getMPIGlobalCommunicator());
        par::Mpi_Reduce(&newElements,&newElements_g,1,MPI_SUM,0,newMesh->getMPIGlobalCommunicator());

        if(!(m_uiMesh->getMPIRankGlobal()))
            std::cout<<"[Ctx]: step : "<<m_uiTinfo._m_uiStep<<"\ttime : "<<m_uiTinfo._m_uiT<<"\told mesh: "<<oldElements_g<<"\tnew mesh:"<<newElements_g<<std::endl;

        this-> grid_transfer(newMesh,transferEvolution, transferConstraint, transferPrimitive);
        this-> update_app_vars();


        std::swap(newMesh,m_uiMesh);
        delete newMesh;

        m_uiIsETSSynced = false;

        
        return 0; 

    }

    template<typename T, typename I>
    int Ctx<T,I>::grid_transfer(ot::Mesh* newMesh, bool transferEvolution, bool transferConstraint, bool transferPrimitive)
    {

        std::vector< ot::DVector<T,I> > eVars;
        std::vector< ot::DVector<T,I> > cVars;
        std::vector< ot::DVector<T,I> > pVars;

        std::vector< ot::DVector<T,I> > eVarsUnzip;
        std::vector< ot::DVector<T,I> > cVarsUnzip;
        std::vector< ot::DVector<T,I> > pVarsUnzip;

        const unsigned int numEVars = m_uiEvolutionVar.size();
        const unsigned int numCVars = m_uiConstrainedVar.size();
        const unsigned int numPVars = m_uiPrimitiveVar.size();

        for(unsigned int i=0; i< m_uiEvolutionVar.size(); i++)
        {
            ot::DVector < T,I > v1;
            v1.VecCreate(newMesh,m_uiEvolutionVar[i].IsGhosted(),m_uiEvolutionVar[i].IsUnzip(), m_uiEvolutionVar[i].IsElemental(),m_uiEvolutionVar[i].GetDof());
            eVars.push_back(v1);
        }

        for(unsigned int i=0; i< m_uiConstrainedVar.size(); i++)
        {
            ot::DVector < T,I > v1;
            v1.VecCreate(newMesh,m_uiConstrainedVar[i].IsGhosted(),m_uiConstrainedVar[i].IsUnzip(), m_uiConstrainedVar[i].IsElemental(),m_uiConstrainedVar[i].GetDof());
            cVars.push_back(v1);
        }

        for(unsigned int i=0; i< m_uiPrimitiveVar.size(); i++)
        {
            ot::DVector < T,I > v1;
            v1.VecCreate(newMesh,m_uiPrimitiveVar[i].IsGhosted(),m_uiPrimitiveVar[i].IsUnzip(), m_uiPrimitiveVar[i].IsElemental(),m_uiPrimitiveVar[i].GetDof());
            pVars.push_back(v1);
        }

        for(unsigned int i=0; i< m_uiEvolutionUnzipVar.size(); i++)
        {
            ot::DVector < T,I > v1;
            v1.VecCreate(newMesh,m_uiEvolutionUnzipVar[i].IsGhosted(),m_uiEvolutionUnzipVar[i].IsUnzip(), m_uiEvolutionUnzipVar[i].IsElemental(),m_uiEvolutionUnzipVar[i].GetDof());
            eVarsUnzip.push_back(v1);
        }

        for(unsigned int i=0; i< m_uiConstraintUnzipVar.size(); i++)
        {
            ot::DVector < T,I > v1;
            v1.VecCreate(newMesh,m_uiConstraintUnzipVar[i].IsGhosted(),m_uiConstraintUnzipVar[i].IsUnzip(), m_uiConstraintUnzipVar[i].IsElemental(),m_uiConstraintUnzipVar[i].GetDof());
            cVarsUnzip.push_back(v1);
        }

        for(unsigned int i=0; i< m_uiPrimitiveUnzipVar.size(); i++)
        {
            ot::DVector < T,I > v1;
            v1.VecCreate(newMesh,m_uiPrimitiveUnzipVar[i].IsGhosted(),m_uiPrimitiveUnzipVar[i].IsUnzip(), m_uiPrimitiveUnzipVar[i].IsElemental(),m_uiPrimitiveUnzipVar[i].GetDof());
            pVarsUnzip.push_back(v1);
        }


        if(transferEvolution)
        {
            for(unsigned int i=0; i< m_uiEvolutionVar.size(); i++)
            {
                assert(m_uiEvolutionVar[i].IsUnzip() == false && m_uiEvolutionVar[i].IsGhosted()==true);
                const unsigned int dof = m_uiEvolutionVar[i].GetDof();
                T* in  = m_uiEvolutionVar[i].GetVecArray();
                T* out = eVars[i].GetVecArray(); 
                
                const unsigned int nPDOF_old = m_uiEvolutionVar[i].GetSizePerDof();
                const unsigned int nPDOF_new = eVars[i].GetSizePerDof();

                for(unsigned int v=0; v < dof; v++)
                    m_uiMesh->interGridTransfer( (in + v*nPDOF_old) , (out + v*nPDOF_new) , newMesh, true);
            
            }
            
        }

        if(transferConstraint)
        {
            for(unsigned int i=0; i< m_uiConstrainedVar.size(); i++)
            {
                assert(m_uiConstrainedVar[i].IsUnzip() == false && m_uiConstrainedVar[i].IsGhosted()==true);
                const unsigned int dof = m_uiConstrainedVar[i].GetDof();
                T* in = m_uiConstrainedVar[i].GetVecArray();
                T* out = cVars[i].GetVecArray();
                
                const unsigned int nPDOF_old = m_uiConstrainedVar[i].GetSizePerDof();
                const unsigned int nPDOF_new = cVars[i].GetSizePerDof();

                for(unsigned int v=0; v < dof; v++)
                    m_uiMesh->interGridTransfer(in + v*nPDOF_old , out + v*nPDOF_new ,newMesh, true);
            
            }

            
        }

        if(transferPrimitive)
        {
            for(unsigned int i=0; i< m_uiPrimitiveVar.size(); i++)
            {
                assert(m_uiPrimitiveVar[i].IsUnzip() == false && m_uiPrimitiveVar[i].IsGhosted()==true);
                const unsigned int dof = m_uiPrimitiveVar[i].GetDof();
                T* in = m_uiPrimitiveVar[i].GetVecArray();
                T* out = pVars[i].GetVecArray();
                
                const unsigned int nPDOF_old = m_uiPrimitiveVar[i].GetSizePerDof();
                const unsigned int nPDOF_new = pVars[i].GetSizePerDof();

                for(unsigned int v=0; v < dof; v++)
                    m_uiMesh->interGridTransfer(in + v*nPDOF_old , out + v*nPDOF_new ,newMesh, true);
                
            }

        }

        std::swap(m_uiEvolutionVar,eVars);
        std::swap(m_uiConstrainedVar,cVars);
        std::swap(m_uiPrimitiveVar,pVars);
        std::swap(m_uiEvolutionUnzipVar,eVarsUnzip);
        std::swap(m_uiConstraintUnzipVar,cVarsUnzip);
        std::swap(m_uiPrimitiveUnzipVar,pVarsUnzip);
        
        eVars.clear();
        pVars.clear();
        cVars.clear();
        
        eVarsUnzip.clear();
        pVarsUnzip.clear();
        cVarsUnzip.clear();
       
    }


    template<typename T, typename I>
    int Ctx<T,I>::increment_ts_info()
    {
        m_uiTinfo._m_uiT += m_uiTinfo._m_uiTh;
        m_uiTinfo._m_uiStep +=1;
        
        return 0; 

    }


 } // end of namespace ts