/**
 * @file massgravCtx.h
 * @author Milinda Fernando
 * @brief Application context class for solving the Einstein equations in MASSGRAVKO formulation. 
 * @version 0.1
 * @date 2019-12-20
 * 
 * @copyright Copyright (c) 2019, University of Utah. 
 * 
 */

#pragma once
#include "ctx.h"
#include "parameters.h"
#include "grDef.h"
#include "grUtils.h"
#include "rhs.h"
#include "physcon.h"
#include "oct2vtk.h"
#include "checkPoint.h"
#include <iostream>
#include "parUtils.h"
#include "dataUtils.h"
//#include "gwExtract.h"
#include "physcon.h"
#include "massgrav_constraints.h"
#include "TwoPunctures.h"

namespace massgrav
{
    class MASSGRAVCtx : public ts::Ctx<DendroScalar, DendroIntL>
    {

        protected:
            /**@brief: evolution var (zip)*/
            DVec m_uiEVar;
            
            /**@brief: evolution var tmp (zip)*/
            DVec m_uiETmp;
            
            /**@brief: constraint var (zip)*/    
            DVec m_uiCVar;

            /**@brief: primitive var (zip)*/
            DVec m_uiPVar;

            /**@brief: Evolution var unzip 0 - unzip in , 1 - unzip out */
            DVec m_uiEUnzip[2];
            
            /**@brief: Constraint var unzip 0 - unzip in , 1 - unzip out */
            DVec m_uiCUnzip[2];

            /**@brief: Primitive var unzip 0 - unzip in , 1 - unzip out */
            DVec m_uiPUnzip[2];

            /**@brief: extracted BH locations*/
            Point m_uiBHLoc[2];

        public :

            /**@brief: default constructor*/
            MASSGRAVCtx(ot::Mesh* pMesh);

            /**@brief: default deconstructor*/
            ~MASSGRAVCtx();

            /**@brief: initial solution*/
            virtual int initialize();
            
            /**
             * @brief computes the MASSGRAV rhs 
             * 
             * @param in : zipped input
             * @param out : zipped output
             * @param sz  : number of variables. 
             * @param time : current time. 
             * @return int : status. (0) on success. 
             */
            virtual int rhs(DVec* in , DVec* out, unsigned int sz , DendroScalar time);

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
            virtual int rhs_blkwise(DVec in, DVec out, const unsigned int* const blkIDs, unsigned int numIds, DendroScalar*  blk_time) const;

            virtual int rhs_blk(const DendroScalar* in, DendroScalar* out, unsigned int dof ,unsigned int local_blk_id, DendroScalar  blk_time) const ;
            
            /**@brief: function execute before each stage
             * @param sIn: stage var in. 
            */
            virtual int pre_stage(DVec sIn); 

            /**@brief: function execute after each stage
             * @param sIn: stage var in. 
            */
            virtual int post_stage(DVec sIn);

            /**@brief: function execute before each step*/
            virtual int pre_timestep(DVec sIn); 

            /**@brief: function execute after each step*/
            virtual int post_timestep(DVec sIn);

            /**@brief: function execute after each step*/
            virtual bool is_remesh();

            /**@brief: write to vtu. */
            virtual int write_vtu();

            /**@brief: writes checkpoint*/
            virtual int write_checkpt();

            /**@brief: restore from check point*/
            virtual int restore_checkpt();

            /**@brief: should be called for free up the contex memory. */
            virtual int finalize();

            /**@brief: pack and returns the evolution variables to one DVector*/
            virtual DVec get_evolution_vars();

            /**@brief: pack and returns the constraint variables to one DVector*/
            virtual DVec get_constraint_vars();

            /**@brief: pack and returns the primitive variables to one DVector*/
            virtual DVec get_primitive_vars();

            /**@brief: updates the application variables from the variable list. */
            virtual int update_app_vars();

            /**@brief: prints any messages to the terminal output. */
            virtual int terminal_output();

            /**@brief: returns the async communication batch size. */
            virtual unsigned int get_async_batch_sz() {return massgrav::MASSGRAV_ASYNC_COMM_K;}




    };

}// end of namespace massgrav
