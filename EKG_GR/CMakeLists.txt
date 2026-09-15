# This is a SEPARATE recompilation of the application with 26 fields.
# Do not link the original 24-field bssn_common object library.
if(WITH_CUDA OR DENDRO_USE_NEW_DERIVS)
  message(FATAL_ERROR "EKS starter: use CPU and legacy explicit derivatives")
endif()
set(EKS_SRC
  src/parameters.cpp src/grUtils.cpp src/rhs.cpp src/sommerfeld_bcs.cpp
  src/derivs.cpp src/physcon.cpp src/profile_params.cpp src/bssn_constraints.cpp
  src/TPCoordTransf.cpp src/TPEquations.cpp src/TPFuncAndJacobian.cpp
  src/TPNewton.cpp src/TPUtilities.cpp src/TwoPunctures.cpp src/dataUtils.cpp
  src/hadrhs.F90 src/trumpet.cpp src/eks_config.cpp)
add_library(eks_common OBJECT ${EKS_SRC})
target_include_directories(eks_common PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/scripts ${CMAKE_CURRENT_SOURCE_DIR}/gencode)
if(DEFINED CHECK_GIT_INCLUDE_PATH)
  target_include_directories(eks_common PUBLIC ${CHECK_GIT_INCLUDE_PATH})
endif()
target_compile_definitions(eks_common PUBLIC
  BSSN_USE_6TH_ORDER_DERIVS BSSN_GW_USE_DENDROLIB
  BSSN_ENABLE_VTU_OUTPUT BSSN_COMPUTE_CONSTRAINTS
  BSSN_ENABLE_VTU_CONSTRAINT_OUTPUT)
target_compile_options(eks_common PUBLIC
  $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-fext-numeric-literals>)
target_link_libraries(eks_common PUBLIC
  dendro5 dendro_gr MPI::MPI_CXX GSL::gsl toml11::toml11
  dendro_git_version_and_date)
add_executable(eksSolver src/eks_main.cpp src/bssnAEH.cpp src/bssnCtx.cpp)
target_link_libraries(eksSolver PRIVATE eks_common m)
# Generated .cpp.inc files are textually included, NOT separate sources.
