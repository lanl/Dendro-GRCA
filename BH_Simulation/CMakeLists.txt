cmake_minimum_required(VERSION 2.8)
project(BH_Simulation)

find_package(OpenMP REQUIRED)
find_package(MPI REQUIRED)


if(OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


if(MPI_COMPILE_FLAGS)
    set(COMPILE_FLAGS "${COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set(LINK_FLAGS "${LINK_FLAGS} ${MPI_LINK_FLAGS}")
endif()


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

option(USE_64BIT_INDICES "Use 64-Bit indices. Reverts to 32-bit if turned off" ON)
option(HILBERT_ORDERING "Use Hilbert ordering instead of Morton" ON)
option(RUN_WEAK_SCALING "Run Entire Weak Scaling" ON)
option(ALLTOALLV_FIX "Use K-way all to all v" ON)
option(SPLITTER_SELECTION_FIX "Turn on Splitter Selection fix" ON)
option(DIM_2 "use the two dimentional sorting" OFF)
option(PROFILE_TREE_SORT "Profile the tree sort code, construction and balancing. " OFF)

#option(KWAY "K parameter for alltoallv_kway" 128)
set(KWAY 128 CACHE INT 128)
set(NUM_NPES_THRESHOLD 2 CACHE INT 2)
#message(${KWAY})

if(DIM_2)
    add_definitions(-DDIM_2)
endif()

if(PROFILE_TREE_SORT)
    add_definitions(-DPROFILE_TREE_SORT)
endif()


if(USE_64BIT_INDICES)
    add_definitions(-DUSE_64BIT_INDICES)
    #message('Configuring 64BIT indices')
endif()

if(HILBERT_ORDERING)
    add_definitions(-DHILBERT_ORDERING)
    #message (Configuring with HIlbert)
endif()

if(RUN_WEAK_SCALING)
    add_definitions(-DRUN_WEAK_SCALING)
    #message('Configuring with Weak Scaling')
endif()

if(ALLTOALLV_FIX)
    add_definitions(-DALLTOALLV_FIX)
    add_definitions(-DKWAY=${KWAY})
endif()

if(SPLITTER_SELECTION_FIX)
    add_definitions(-DSPLITTER_SELECTION_FIX)
    add_definitions(-DNUM_NPES_THRESHOLD=${NUM_NPES_THRESHOLD})
endif()


include_directories(${PROJECT_BINARY_DIR}
        ../include
        ../include/test
        ${MPI_INCLUDE_PATH}
        ../examples/include
        include
        )


set (COMMON_FILES ../include/seqUtils.h
        ../include/seqUtils.tcc
        ../include/ompUtils.h
        ../include/ompUtils.tcc
        ../include/dendro.h
        ../include/binUtils.h
        ../src/binUtils.cpp
        ../include/hcurvedata.h
        ../src/hcurvedata.cpp
        ../include/genPts_par.h
        ../src/genPts_par.cpp
        ../include/TreeNode.h
        ../src/TreeNode.cpp
        ../include/treenode2vtk.h
        ../src/treenode2vtk.cpp
        ../include/parUtils.h
        ../include/parUtils.tcc
        ../src/parUtils.cpp
        ../include/octUtils.h
        ../src/octUtils.cpp
        ../include/point.h
        ../src/point.cpp

        )

set(SOURCE_FILES include/bh.h
        include/interpolate.h
        include/visualize.h
        src/visualize.cpp
        src/bh.cpp
        src/bh_example.cpp
        ${COMMON_FILES}
        )

add_executable(BHOctreeTest ${SOURCE_FILES})
target_link_libraries(BHOctreeTest ${MPI_LIBRARIES} m)