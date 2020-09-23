//
// Created by milinda on 10/6/16.
//

/**
 * @author: Milinda Fernando
 * School of Computing , University of Utah.
 *
 * @brief: Contains the functions related to visulizing octress, and mesh data in vtk format.
 *
 * */


#ifndef SFCSORTBENCH_VISUALIZE_H
#define SFCSORTBENCH_VISUALIZE_H

#include "bh.h"
#include "TreeNode.h"
#include <algorithm>

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#define VTK_HEXAHEDRON 12

/**
 * @author Milinda Fernando
 * @breif Writes the octree to the vtk file format with the function values. (currenlty all the 24 functions. )
 * @param[in] octree -> octree which needs to be written
 * @param[in] rank -> mpi rank
 * @param[in] function values corresponding to each octant.
 */

void octree2VTK(const std::vector<ot::TreeNode>& nodes, unsigned int rank, const double * u,std::string vtk_file_name);


#endif //SFCSORTBENCH_VISUALIZE_H
