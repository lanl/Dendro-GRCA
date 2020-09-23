# only MASSGRAV (GR) options
option(MASSGRAV_PROFILE_HUMAN_READABLE " enables human readable profile outputs" ON)
option(MASSGRAV_CONSEC_COMM_SELECT "rank selection performed using consecative ranks (works with arbitary number of procs.) if off ranks selected in binary tree fashion. (works only with powers of 2)" ON)
option(MASSGRAV_ENABLE_VTU_OUTPUT "Enables vtu output" ON)
option(MASSGRAV_COMPUTE_CONSTRAINTS "Enable constraint equation computations" OFF)
option(MASSGRAV_ENABLE_VTU_CONSTRAINT_OUTPUT "write constrints to vtu" ON)
option(MASSGRAV_ETA_FUNCTION "Use function as ETA damping" OFF)
option(MASSGRAV_KERR_SCHILD_TEST "compare with the analytical Kerr-Schild sol. for a single BH" OFF)
option(MASSGRAV_EXTRACT_BH_LOCATIONS "compute the BH locations as time evolves" ON)
option(MASSGRAV_REFINE_BASE_EH "enable refine based on BH event horizons" OFF)
option(MASSGRAV_EXTRACT_GRAVITATIONAL_WAVES "extract GW waves" OFF)


if(MASSGRAV_PROFILE_HUMAN_READABLE)
    add_definitions(-DMASSGRAV_PROFILE_HUMAN_READABLE)
endif()

if(MASSGRAV_CONSEC_COMM_SELECT)
    add_definitions(-DMASSGRAV_CONSEC_COMM_SELECT)
endif()

if(MASSGRAV_ENABLE_VTU_OUTPUT)
    add_definitions(-DMASSGRAV_ENABLE_VTU_OUTPUT)
endif()

if(MASSGRAV_COMPUTE_CONSTRAINTS)
    add_definitions(-DMASSGRAV_COMPUTE_CONSTRAINTS)
endif()

if(MASSGRAV_ENABLE_VTU_CONSTRAINT_OUTPUT)
    add_definitions(-DMASSGRAV_ENABLE_VTU_CONSTRAINT_OUTPUT)
endif()

if(MASSGRAV_ETA_FUNCTION)
    add_definitions(-DUSE_ETA_FUNC)
endif()

if(MASSGRAV_KERR_SCHILD_TEST)
    add_definitions(-DMASSGRAV_KERR_SCHILD_TEST)
endif()

if(MASSGRAV_EXTRACT_BH_LOCATIONS)
    add_definitions(-DMASSGRAV_EXTRACT_BH_LOCATIONS)
endif()


if(MASSGRAV_REFINE_BASE_EH)
    add_definitions(-DMASSGRAV_REFINE_BASE_EH)
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

set(MASSGRAV_INC include/parameters.h
        include/bh.h
        include/grUtils.h
        include/grUtils.tcc
        include/rhs.h
        include/derivs.h
        include/physcon.h
        include/profile_params.h
        include/massgrav_constraints.h
        include/TPUtilities.h
        include/TwoPunctures.h
        include/dataUtils.h
        )

set(MASSGRAV_SRC src/rkMASSGRAV.cpp
           src/parameters.cpp
           src/grUtils.cpp
           src/rhs.cpp
           src/derivs.cpp
           src/physcon.cpp
           src/profile_params.cpp
           src/massgrav_constraints.cpp
           src/TPCoordTransf.cpp
           src/TPEquations.cpp
           src/TPFuncAndJacobian.cpp
           src/TPNewton.cpp
           src/TPUtilities.cpp
           src/TwoPunctures.cpp
           src/dataUtils.cpp

        )


set(SOURCE_FILES src/gr.cpp include/gr.h ${MASSGRAV_INC} ${MASSGRAV_SRC})
add_executable(massgravSolver ${SOURCE_FILES})
install(TARGETS massgravSolver RUNTIME DESTINATION bin)

target_include_directories(massgravSolver PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(massgravSolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(massgravSolver PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
target_include_directories(massgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(massgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/include/test)
target_include_directories(massgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/examples/include)
target_include_directories(massgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/FEM/include)
target_include_directories(massgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/ODE/include)
target_include_directories(massgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/LinAlg/include)
target_include_directories(massgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/IO/vtk/include)
target_include_directories(massgravSolver PRIVATE ${CMAKE_SOURCE_DIR}/IO/zlib/inc)
target_include_directories(massgravSolver PRIVATE ${MPI_INCLUDE_PATH})
target_include_directories(massgravSolver PRIVATE ${GSL_INCLUDE_DIRS})

target_link_libraries(massgravSolver dendro5 ${LAPACK_LIBRARIES} ${MPI_LIBRARIES} ${GSL_LIBRARIES} m)

