# only QUADGRAV (GR) options
option(QUADGRAV_PROFILE_HUMAN_READABLE " enables human readable profile outputs" ON)
option(QUADGRAV_CONSEC_COMM_SELECT "rank selection performed using consecative ranks (works with arbitary number of procs.) if off ranks selected in binary tree fashion. (works only with powers of 2)" ON)
option(QUADGRAV_ENABLE_VTU_OUTPUT "Enables vtu output" ON)
option(QUADGRAV_COMPUTE_CONSTRAINTS "Enable constraint equation computations" OFF)
option(QUADGRAV_ENABLE_VTU_CONSTRAINT_OUTPUT "write constrints to vtu" ON)
option(QUADGRAV_ETA_FUNCTION "Use function as ETA damping" OFF)
option(QUADGRAV_KERR_SCHILD_TEST "compare with the analytical Kerr-Schild sol. for a single BH" OFF)
option(QUADGRAV_EXTRACT_BH_LOCATIONS "compute the BH locations as time evolves" ON)
option(QUADGRAV_REFINE_BASE_EH "enable refine based on BH event horizons" OFF)
option(QUADGRAV_EXTRACT_GRAVITATIONAL_WAVES "extract GW waves" OFF)


if(QUADGRAV_PROFILE_HUMAN_READABLE)
    add_definitions(-DQUADGRAV_PROFILE_HUMAN_READABLE)
endif()

if(QUADGRAV_CONSEC_COMM_SELECT)
    add_definitions(-DQUADGRAV_CONSEC_COMM_SELECT)
endif()

if(QUADGRAV_ENABLE_VTU_OUTPUT)
    add_definitions(-DQUADGRAV_ENABLE_VTU_OUTPUT)
endif()

if(QUADGRAV_COMPUTE_CONSTRAINTS)
    add_definitions(-DQUADGRAV_COMPUTE_CONSTRAINTS)
endif()

if(QUADGRAV_ENABLE_VTU_CONSTRAINT_OUTPUT)
    add_definitions(-DQUADGRAV_ENABLE_VTU_CONSTRAINT_OUTPUT)
endif()

if(QUADGRAV_ETA_FUNCTION)
    add_definitions(-DUSE_ETA_FUNC)
endif()

if(QUADGRAV_KERR_SCHILD_TEST)
    add_definitions(-DQUADGRAV_KERR_SCHILD_TEST)
endif()

if(QUADGRAV_EXTRACT_BH_LOCATIONS)
    add_definitions(-DQUADGRAV_EXTRACT_BH_LOCATIONS)
endif()


if(QUADGRAV_REFINE_BASE_EH)
    add_definitions(-DQUADGRAV_REFINE_BASE_EH)
endif()


#message(${KWAY})

# some debug flags. uncomment to print the debug information.
#add_definitions(-DDEBUG_RK_SOLVER)
#add_definitions(-DDEBUG_UPWIND_INTERP)
#add_definitions(-DDEBUG_DOWNWIND_INTERP)
#add_definitions(-DDEBUG_DERIVS_COMP)
#add_definitions(-DMATVEC_PROFILE)
#add_definitions(-DFEM_ACCUMILATE_ONES_TEST)
#add_definitions(-DNLSM_DIRICHLET_BDY)

set(QUADGRAV_INC include/parameters.h
        include/bh.h
        include/grUtils.h
        include/grUtils.tcc
        include/rhs.h
        include/derivs.h
        include/physcon.h
        include/profile_params.h
        include/quadgrav_constraints.h
        include/TPUtilities.h
        include/TwoPunctures.h
        include/dataUtils.h
        )

set(QUADGRAV_SRC src/rkQUADGRAV.cpp
           src/parameters.cpp
           src/grUtils.cpp
           src/rhs.cpp
           src/derivs.cpp
           src/physcon.cpp
           src/profile_params.cpp
           src/quadgrav_constraints.cpp
           src/TPCoordTransf.cpp
           src/TPEquations.cpp
           src/TPFuncAndJacobian.cpp
           src/TPNewton.cpp
           src/TPUtilities.cpp
           src/TwoPunctures.cpp
           src/dataUtils.cpp

        )


set(SOURCE_FILES src/gr.cpp include/gr.h ${QUADGRAV_INC} ${QUADGRAV_SRC})
add_executable(quadgravSolver ${SOURCE_FILES})
install(TARGETS quadgravSolver RUNTIME DESTINATION bin)

target_include_directories(quadgravSolver PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(quadgravSolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/include/test)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
target_include_directories(quadgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(quadgravSolver PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(quadgravSolver PRIVATE ${GSL_INCLUDE_DIRS})

target_link_libraries(quadgravSolver dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)

