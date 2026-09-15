// * Created by milinda on 8/10/18
/**
 * @author Milinda Fernando.
 * School of Computing, University of Utah.
 * @brief Computation of the rhs using cuda
 * Created by milinda on 8/10/18
 * */

#ifndef SFCSORTBENCH_CUDARHS_H
#define SFCSORTBENCH_CUDARHS_H

#include <device_launch_parameters.h>

#include "block_cu.h"
#include "bssn_rhs_deriv_mem_cuda.h"
#include "cudaUtils.cuh"
#include "cudaUtils.h"
#include "cuda_runtime.h"
#include "gpuTest.cuh"
#include "grDef.h"
#include "parameters.h"
#include "params_cu.h"
#include "profile_gpu.h"
#include "rhs.h"
#include "rhs_bssn.cuh"

namespace cuda {

/***
 * @brief performs kernel pre-launch tasks and launch the bssnrhs kernel
 *
 **/
void computeRHS(double** unzipVarsRHS, const double** uZipVars,
                const ot::Block* dendroBlockList, unsigned int numBlocks,
                const cuda::BSSNComputeParams* bssnPars, dim3 blockDim,
                const Point& pt_min, const Point& pt_max,
                unsigned int numStreams, const double curr_time,
                const double** uZipConstVars, unsigned int device = 0);

}  // namespace cuda

#endif  // SFCSORTBENCH_CUDARHS_H
