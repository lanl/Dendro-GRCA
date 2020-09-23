/**
*@brief Contains utility functions for MASSGRAV simulation.
*/


#ifndef SFCSORTBENCH_GRUTILS_H
#define SFCSORTBENCH_GRUTILS_H

#include "point.h"
#include "parameters.h"
#include "mesh.h"
#include "block.h"
#include "parUtils.h"
#include "json.hpp"
#include "dendroProfileParams.h"
#include "profile_params.h"
#include "swsh.h"
#include "lebedev.h"
#include "grDef.h"


using json = nlohmann::json;
namespace massgrav
{


 /**
 * @brief internal variables needed for rk update.
 * */


 /**
  * @brief: Read the parameter file and initialize the variables in parameters.h file.
  * @param[in] fName: file name
  * @param[in] comm: MPI communicator.
  * */
  void readParamFile(const char * fName,MPI_Comm comm);


 /**
  * @brief Two puncture intiial data from HAD code/ 
  * 
  * @param xx1 : x coord
  * @param yy1 : y coord
  * @param zz1 : z coord
  * @param var : initialized massgrav variables for the grid points
  */
 void punctureData(const double xx1,const double yy1,const double zz1, double *var);
 
 /**
  * @brief compute the static Kerr-Schild BH data
  * 
  * @param xx1 : x coord
  * @param yy1 : y coord
  * @param zz1 : z coord
  * @param var : initialized massgrav variables for the grid points
  */
 void KerrSchildData(const double xx1,const double yy1,const double zz1, double *var);
 
 /**
  * @brief Compute Minkowski space in Carteisan coordinate
  * @param xx1 : x coord
  * @param yy1 : y coord
  * @param zz1 : z coord
  * @param var : initialized massgrav variables for the grid points
  */
 void FlatMinkowski(const double xx1,const double yy1,const double zz1, double *var);
 
 /**
  * @brief add artificial noise to the initial data.
  * @param xx1 : x coord
  * @param yy1 : y coord
  * @param zz1 : z coord
  * @param var : initialized massgrav variables for the grid points
  */
 void noiseData(const double xx1,const double yy1,const double zz1, double *var);

 /**
  * @brief fake initial data. 
  * @param xx1 : x coord
  * @param yy1 : y coord
  * @param zz1 : z coord
  * @param var : initialized massgrav variables for the grid points
  */
 void fake_initial_data(double x, double y, double z, double *u);

namespace trumpet_data 
{
 void trumpetData(const double xx1,const double yy1,const double zz1, double *var);
  void bndcnd(double h, double &x, double y[], double dydx[]) ;
void derivs(double x, double y[], double dydx[]) ;
void hunt(double xx[], int n, double x, int *jlo) ;
void rkck(   double y[], double dydx[], int n, double x, double h,
             double yout[], double yerr[],
             void (*derivs)(double, double [], double [])   ) ;
void rkqs(   double y[], double dydx[], int n, double *x, double htry,
             double eps, double yscal[], double *hdid, double *hnext,
             void (*derivs)(double, double [], double [])   ) ;
void odeint(  double ystart[], int nvar, double x1, double x2,
              double eps, double h1, double hmin, int *nok, int *nbad,
              void (*derivs)( double, double [], double [] ),
              void  (*rkqs) ( double [], double [], int, double *,
                              double, double, double [], double *,
                              double *,
                              void (*)( double, double [], double [] ) ),
              int kount ) ;
double interpolation3( double xp[], double yp[], int np, double xb,
                       int *n_nearest_pt ) ;
double interpolation4 ( double xp[], double yp[], int np, double xb,
                        int *n_nearest_pt  ) ;
}  // end of namespace trumpet_data


 /**
  * @brief: Generates block adaptive octree for the given binary blockhole problem.
  * @param[out] tmpNodes: created octree tmpNodes
  * @param[in] pt_min: block min point
  * @param[in] pt_max: block max point
  * @param[in] regLev: regular grid level
  * @param[in] maxDepth: maximum refinement level. 
  * @param[in] comm: MPI communicator. 
  * */
 void blockAdaptiveOctree(std::vector<ot::TreeNode>& tmpNodes,const Point& pt_min,const Point & pt_max,const unsigned int regLev,const unsigned int maxDepth,MPI_Comm comm);

 /**
  * @brief Compute the wavelet tolerance as a function of space. 
  * 
  * @param x : x coord. 
  * @param y : y coord
  * @param z : z coord
  * @param tol_min : min. tolerance value. 
  * @return double 
  */
 double computeWTol(double x,double y,double z,double tol_min);

   /**
   * @breif: Compute L2 constraint norms. 
   */
   template <typename T>
   double computeConstraintL2Norm(const T* constraintVec, const T* maskVec, unsigned int lbegin, unsigned int lend,MPI_Comm comm);

   /**
    * @breif: Compute L2 constraint norms. 
    */
   template <typename T>
   double computeConstraintL2Norm(const ot::Mesh* mesh, const T* constraintVec,const T* maskVector,T maskthreshoold);

   /**
    * @breif write constraints to a file. 
    */
   template<typename T>
   double extractConstraints(const ot::Mesh* mesh,const T** constraintVar,const T* maskVec, double maskthreshoold,unsigned int timestep);


   void writeBLockToBinary(const double **unzipVarsRHS,unsigned int offset,const double *pmin, const double *pmax,double* bxMin,double * bxMax, const unsigned int *sz,unsigned int blkSz,double dxFactor,const char* fprefix);


   /**
    * @brief Performs elemental based artificial dissipation based on DG book by Hesthavan
    * @param[in] mesh: const pointer to the mesh. 
    * @param[in/out] zipVars: gr variables that need to dissipation applied to.  
    * @param[in] numVars: number of variables
    * @param[in] nc: cutoff order should be 0,4 since we are using 4th order elements (during interpolations etc). 
    * @param[in] s: dissipation parameter s , 
    * @note nc, s increasing, will lower the dissipation. 
    */
   template<typename T>
   void artificial_dissipation(ot::Mesh * pMesh , T** zipVars, unsigned int numVars, unsigned int nc, unsigned int s,bool isGhostEx=false);


   unsigned int getOctantWeight(const ot::TreeNode* pNode);




}// end of namespace massgrav



namespace massgrav
{

    namespace timer
    {

        /**@brief initialize all the flop counters. */
        void initFlops();

        /**@brief clears the snapshot counter for time profiler variables*/
        void resetSnapshot();


        /**@brief reduce min mean max.
         * @param [in] stat: local time
         * @param [out] stat_g 0-min, 1-mean 2-max
        * */
       template<typename T>
       void computeOverallStats(T *stat, T *stat_g, MPI_Comm comm)
       {
           int rank,npes;
           MPI_Comm_size(comm,&npes);
           MPI_Comm_rank(comm,&rank);

           par::Mpi_Reduce(stat,stat_g,1,MPI_MIN,0,comm);
           par::Mpi_Reduce(stat,stat_g+1,1,MPI_SUM,0,comm);
           par::Mpi_Reduce(stat,stat_g+2,1,MPI_MAX,0,comm);
           stat_g[1]/=(npes);

       }


        /** @breif : printout the profile parameters. */
        void profileInfo(const char* filePrefix,const ot::Mesh* pMesh);

        /** @breif : printout the profile parameters (intermediate profile information). */
        void profileInfoIntermediate(const char* filePrefix,const ot::Mesh* pMesh,const unsigned int currentStep);


    }


}


namespace GW
{
    /**
    * @brief : debug function to write psi4 to interpolated to spheres. 
    */
    void psi4ShpereDump(const ot::Mesh* mesh, DendroScalar ** cVar,unsigned int timestep,double time );
    

}// end of namespace GW




#include "grUtils.tcc"


#endif //SFCSORTBENCH_GRUTILS_H
