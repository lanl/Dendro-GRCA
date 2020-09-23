//
// Created by milinda on 1/16/19.
//

/**
 * @author Milinda Fernando
 * School of Computing univerity of Utah
 * @brief Contain utility functions to perform processing simulation data and post processing.
 *
 */

#ifndef DENDRO_5_0_DATAUTILS_H
#define DENDRO_5_0_DATAUTILS_H

#include "mesh.h"
#include "TreeNode.h"
#include "parameters.h"
#include "point.h"
#include "grDef.h"


namespace quadgrav
{

    /**
     * @brief performs the BH coordinates extraction
     * @param[in] pMesh: current mesh data structure
     * @param[in] var: considering variable for BH extraction
     * @param[in] tolerance: tolerance value for coords extraction
     * @param[in] ptIn: previous known location of the BHs.
     * @param[in] numPt: number of points
     * @param[out] ptOut: new locations of the BHs. (only rank =0 of active comm will have the coords)
     * */
    void extractBHCoords(const ot::Mesh* pMesh, const DendroScalar* var,double tolerance, const Point* ptIn, unsigned int numPt,Point* ptOut);
    
    
    /**
     *@brief write the extracted BH coordinates to a file. 
     *@param[in] ptLocs: point locations
     *@param[in] numPt: number of points
     *@param[in] timestep: time step 
     */
    void writeBHCoordinates(const ot::Mesh* pMesh,const Point* ptLocs, unsigned int numPt,unsigned int timestep);


    /**
     * @brief refine based on the black hole location locations. 
     * @param[in] pMesh : pointer to the mesh. 
     * @param[in] bhLoc : BH location 
     * @param[in] r: radius to refine based on the bh location. (square block refinement. )
     */
    bool isRemeshBH(const ot::Mesh* pMesh, const Point* bhLoc, const double* r);

    /**
     * @brief refine only based on the alpha variable event horizon. 
     * 
     * @param pMesh : pointer to the mesh 
     * @param unzipVec : unzip vars. 
     * @param refine_th : refine tol for alpha
     * @param coarsen_th : coarsend threshold for alpha
     * @return true : is mesh need to be changed
     * @return false : otherwise. 
     */
    bool isRemeshEH(const ot::Mesh* pMesh, const double ** unzipVec, unsigned int vIndex, double refine_th, double coarsen_th, bool isOverwrite=true);




} // end of namespace quadgrav

#endif //DENDRO_5_0_DATAUTILS_H
