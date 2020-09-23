/**
 * @file aeh_solver.h
 * @author: Sharvaree Vadgama
 * @author: Milinda Fernando
 * School of Computing, University of Utah. 
 * @brief : To solve the Apperent event horizon solver for MASSGRAV formulation
 * @version 0.1
 * @date 2019-11-10
 * 
 * @copyright Copyright (c) 2019
 * 
 * 
 */
#pragma once

#include "mesh.h"
#include "mpi.h"
#include "grDef.h"

namespace massgrav
{

    enum AEHErrorType {SUCCESS, MAX_ITERATIONS_REACHED};

    template<typename T>
    void theta_rhs(T** s, T** uZipVars, unsigned int offset, const T *pmin, const T *pmax, const unsigned int *sz, const unsigned int& bflag)
    {
        const double *alpha = &uZipVars[VAR::U_ALPHA][offset];
        const double *chi = &uZipVars[VAR::U_CHI][offset];
        const double *K = &uZipVars[VAR::U_K][offset];
        const double *gt0 = &uZipVars[VAR::U_SYMGT0][offset];
        const double *gt1 = &uZipVars[VAR::U_SYMGT1][offset];
        const double *gt2 = &uZipVars[VAR::U_SYMGT2][offset];
        const double *gt3 = &uZipVars[VAR::U_SYMGT3][offset];
        const double *gt4 = &uZipVars[VAR::U_SYMGT4][offset];
        const double *gt5 = &uZipVars[VAR::U_SYMGT5][offset];
        const double *beta0 = &uZipVars[VAR::U_BETA0][offset];
        const double *beta1 = &uZipVars[VAR::U_BETA1][offset];
        const double *beta2 = &uZipVars[VAR::U_BETA2][offset];
        const double *At0 = &uZipVars[VAR::U_SYMAT0][offset];
        const double *At1 = &uZipVars[VAR::U_SYMAT1][offset];
        const double *At2 = &uZipVars[VAR::U_SYMAT2][offset];
        const double *At3 = &uZipVars[VAR::U_SYMAT3][offset];
        const double *At4 = &uZipVars[VAR::U_SYMAT4][offset];
        const double *At5 = &uZipVars[VAR::U_SYMAT5][offset];
        const double *Gt0 = &uZipVars[VAR::U_GT0][offset];
        const double *Gt1 = &uZipVars[VAR::U_GT1][offset];
        const double *Gt2 = &uZipVars[VAR::U_GT2][offset];
        const double *B0 = &uZipVars[VAR::U_B0][offset];
        const double *B1 = &uZipVars[VAR::U_B1][offset];
        const double *B2 = &uZipVars[VAR::U_B2][offset];

        // 1. Allocate memory for derivative variables. 
        

        // 2. Need to compute the derivatives for theta rhs. 


        // 3. Call the theta_rhs code generated from sympyGR


        // 4. Deallocate Derivatiove memory. 
        





    }



    /**
     * @brief :To solve the AEH using parabolic approach. 
     * @tparam T : type of the evolving var (double, float)
     * @param pMesh : underlying mesh from Dendro
     * @param s : surface normal to the Apparent Event Horizon(AEH) 3-vector (initial guess solution to $S_2$)
     * @param massgravVars : MASSGRAV variables at a given time
     * @param tol : tolerance for the time step iterations
     * @param max_iter : maximum number of iterations
     * @return int : error code (0 for success, )
     */
    template<typename T>
    int aeh_solver(const ot::Mesh* pMesh, T** s, T** massgravVars, T tol, unsigned int max_iter)
    {
        // Notes: We might need to change the mesh during the solver to adaptively capture the mesh. (Future work)

        // 1. compute TH based on the CFL of the grid. 
        double th;

        // 
        double T = th*max_iter;
        // 2. use ETS class to evolve s until converges. 
        for(double t=0; t < T; t+=th )
        {
            // use ETS to advance in time. 
            


            // check for the convergence criteria. 
            //if(converged ) break;


        }

        
    }

}

