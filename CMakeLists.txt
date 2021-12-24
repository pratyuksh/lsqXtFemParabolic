CMAKE_MINIMUM_REQUIRED( VERSION 3.15 )

message ("Recognized hostname:\t" $ENV{HOSTNAME})
message ("Recognized user:\t" $ENV{USER})

project(LsqXtFemParabolic)

# Manage dependencies
include ("cmake/GetFMT.cmake")
include ("cmake/GetJSON.cmake")
include("cmake/warning_flags.cmake")
include("cmake/clang_tidy.cmake")

if (($ENV{USER} MATCHES "pratyuksh"))
    
    message ("\n\n\t---------------------\t")
    message (    "\tBuild on local system\t")
    message (    "\t---------------------\t\n\n")
    
    option(USE_LIB_PARDISO "Use Pardiso library" ON)
    
    # Manage dependencies
    include ("cmake/GetGTest.cmake")
    
    # compile flags
    set(CMAKE_CXX_COMPILER "g++")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -O1 -fopenmp")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp -DNDEBUG -msse2 -march=native")
    message ("Compiling with:\t" ${CMAKE_CXX_COMPILER})
    #message ("CXX_DEBUG flags debug:\t" ${CMAKE_CXX_FLAGS_DEBUG})
    #message ("CXX_RELEASE flags release:\t" ${CMAKE_CXX_FLAGS_RELEASE})
    
    if (CMAKE_BUILD_TYPE STREQUAL "Debug") # Debug mode
        message ("CXX_DEBUG flags debug:\t" ${CMAKE_CXX_FLAGS_DEBUG})
        # add_compile_definitions(MYVERBOSE)
    endif (CMAKE_BUILD_TYPE STREQUAL "Debug")

    if (CMAKE_BUILD_TYPE STREQUAL "Release") # Release mode
        message ("CXX_RELEASE flags release:\t" ${CMAKE_CXX_FLAGS_RELEASE})
    endif (CMAKE_BUILD_TYPE STREQUAL "Release")
    
    # MPI libraries
    # set(MPI_INC_DIR "/usr/include/openmpi/")
    # set(MPI_C_LIB "/usr/lib/openmpi/lib/libmpi.so")
    # set(MPI_CXX_LIB "/usr/lib/openmpi/lib/libmpi_cxx.so")

    # external libraries
    set(PROJECT_LIB_DIR "$ENV{HOME}/Documents/softwares")
    
    # MFEM
    set(MFEM_DIR "${PROJECT_LIB_DIR}/mfem-4.1-custom/build-serial")
    set(MFEM_LIB_DIR "${PROJECT_LIB_DIR}/mfem-4.1-custom/build-serial")
    
    # Eigen
    set(EIGEN_INC_DIR "${PROJECT_LIB_DIR}/eigen-3.3.7")

    # Linear algebra
    message ("USE Pardiso library:\t" ${USE_LIB_PARDISO})
    add_library(linear_algebra INTERFACE)
    if(USE_LIB_PARDISO) # Pardiso library, else uses MKL Pardiso
        add_compile_definitions(LIB_PARDISO)
        include ("cmake/GetPardiso_local.cmake")
        target_link_libraries(linear_algebra INTERFACE ${PARDISO})
    endif(USE_LIB_PARDISO)
    unset(USE_LIB_PARDISO CACHE)
    target_link_libraries(linear_algebra INTERFACE -L${PROJECT_LIB_DIR}/intel_mkl/install/mkl/2021.3.0/lib/intel64
    -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core 
    -L${PROJECT_LIB_DIR}/intel_mkl/install/compiler/2021.3.0/linux/compiler/lib/intel64_lin -liomp5 -lpthread -lm -ldl)


elseif (($ENV{USER} MATCHES "prbansal"))
    
    MESSAGE("\n\n\t----------------------\t")
    MESSAGE(    "\tBuild on Euler cluster\t")
    MESSAGE(    "\t----------------------\t\n\n")
    
    option(USE_LIB_PARDISO "Use Pardiso library" ON)
    
    # Manage dependencies
    include ("cmake/GetGTest_euler.cmake")
    
    # compile flags
    #set(CMAKE_CXX_COMPILER "gcc")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -O1")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -msse2 -march=native")
    MESSAGE ("Compiling with:\t" ${CMAKE_CXX_COMPILER})
    #MESSAGE ("CXX_DEBUG flags debug:\t" ${CMAKE_CXX_FLAGS_DEBUG})
    #MESSAGE ("CXX_RELEASE flags release:\t" ${CMAKE_CXX_FLAGS_RELEASE})
    
    if (CMAKE_BUILD_TYPE STREQUAL "Debug") # Debug mode
        message ("CXX_DEBUG flags debug:\t" ${CMAKE_CXX_FLAGS_DEBUG})
        # add_compile_definitions(MYVERBOSE)
    endif (CMAKE_BUILD_TYPE STREQUAL "Debug")

    if (CMAKE_BUILD_TYPE STREQUAL "Release") # Release mode
        message ("CXX_RELEASE flags release:\t" ${CMAKE_CXX_FLAGS_RELEASE})
    endif (CMAKE_BUILD_TYPE STREQUAL "Release")
    
    # MPI libraries
    # set(MPI_INC_DIR "/cluster/apps/gcc-8.2.0/openmpi-4.0.2-vvr7fdofwljfy4qgkkhao36r5qx44hni/include")
    # set(MPI_C_LIB "/cluster/apps/gcc-8.2.0/openmpi-4.0.2-vvr7fdofwljfy4qgkkhao36r5qx44hni/lib/libmpi.so")
    # set(MPI_CXX_LIB "/cluster/apps/gcc-8.2.0/openmpi-4.0.2-vvr7fdofwljfy4qgkkhao36r5qx44hni/lib/libmpi_cxx.so")
    
    # external libraries
    set(PROJECT_LIB_DIR "$ENV{HOME}/Documents/libraries")
    
    # MFEM
    set(MFEM_DIR "${PROJECT_LIB_DIR}/mfem-4.1-custom/build-serial")
    set(MFEM_LIB_DIR "${PROJECT_LIB_DIR}/mfem-4.1-custom/build-serial")
    
    # Eigen
    set(EIGEN_INC_DIR "/cluster/apps/eigen/3.2.1/x86_64/gcc_4.8.2/serial/include/eigen3")
    
    #add_library(linear_algebra INTERFACE)
    #target_link_libraries(linear_algebra INTERFACE -L$ENV{OPENBLAS_ROOT}/lib -lopenblas)
    #target_link_libraries(linear_algebra INTERFACE -L$ENV{NETLIB_SCALAPACK_ROOT}/lib -lscalapack)
    
    # Linear algebra
    message ("USE Pardiso library:\t" ${USE_LIB_PARDISO})
    add_library(linear_algebra INTERFACE)
    if(USE_LIB_PARDISO) # Pardiso library, else uses MKL Pardiso
        add_compile_definitions(LIB_PARDISO)
        include ("cmake/GetPardiso_euler.cmake")
        target_link_libraries(linear_algebra INTERFACE ${PARDISO})
    endif(USE_LIB_PARDISO)
    unset(USE_LIB_PARDISO CACHE)
    target_link_libraries(linear_algebra INTERFACE
    -L/cluster/apps/intel/parallel_studio_xe_2018_r1/mkl/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
    -L/cluster/apps/intel/parallel_studio_xe_2018_r1/mkl/../compiler/lib/intel64 -liomp5 -lpthread -lm -ldl)

endif()

# message ("\nMPI Libraries:")
# message (${MPI_INC_DIR})
# message (${MPI_C_LIB})
# message (${MPI_CXX_LIB})

message ("\nInclude directories:")
message ("MFEM_DIR: " ${MFEM_DIR})

message ("\nLink libraries:")
message ("MFEM_LIB_DIR: " ${MFEM_LIB_DIR})

# include_directories(${MPI_INC_DIR}
#                     ${MFEM_DIR}
#                     ${EIGEN_INC_DIR})
include_directories(${MFEM_DIR}
                    ${EIGEN_INC_DIR})

link_directories(${MFEM_LIB_DIR})


##########
## Core ##
##########
add_library(Core OBJECT)
target_compile_features(Core PUBLIC cxx_std_17)
target_include_directories(Core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(Core PUBLIC FMT)
target_link_libraries(Core PUBLIC JSON)
target_link_libraries(Core PUBLIC stdc++fs)
target_link_libraries(Core PRIVATE warning_flags)
    
if(CLANG_TIDY_EXE)
    set_target_properties(Core PROPERTIES CXX_CLANG_TIDY "${DO_CLANG_TIDY}")
endif()


############
## MyMFEM ##
############

# MyMFEM; serial
add_library(MyMfem OBJECT)
target_compile_features(MyMfem PUBLIC cxx_std_17)
target_include_directories(MyMfem PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(MyMfem PUBLIC Core)
target_link_libraries(MyMfem PUBLIC linear_algebra)
target_link_libraries(MyMfem PRIVATE -fopenmp "${MFEM_LIB_DIR}/libmfem.a" -lmfem)

#############
## Pardiso ##
#############

add_library(Pardiso OBJECT)
target_compile_features(Pardiso PUBLIC cxx_std_17)
target_include_directories(Pardiso PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(Pardiso PUBLIC Core)
target_link_libraries(Pardiso PRIVATE linear_algebra)

##########
## Heat ##
##########

# Heat library
add_library(LibHeat OBJECT)
target_link_libraries(LibHeat PRIVATE Core)

# SparseHeat library
add_library(LibSparseHeat OBJECT)
target_link_libraries(LibSparseHeat PRIVATE Core)

# Heat, serial
add_executable(heat)
target_link_libraries(heat PRIVATE Core)
target_link_libraries(heat PRIVATE LibHeat)
target_link_libraries(heat PRIVATE MyMfem)
target_link_libraries(heat PRIVATE Pardiso)

# SparseHeat, serial
add_executable(sparse_heat)
target_link_libraries(sparse_heat PRIVATE Core)
target_link_libraries(sparse_heat PRIVATE LibHeat)
target_link_libraries(sparse_heat PRIVATE LibSparseHeat)
target_link_libraries(sparse_heat PRIVATE MyMfem)
target_link_libraries(sparse_heat PRIVATE Pardiso)


############
## MLMCMC ##
############

# add_library(LibMlmcmc OBJECT)
# target_link_libraries(LibMlmcmc PRIVATE Core)

# add_library(LibMlmcmcMpi OBJECT)
# target_link_libraries(LibMlmcmcMpi PRIVATE Core)

###############
## MCMC Heat ##
###############

# serial
# add_executable(mcmc_heat)
# target_link_libraries(mcmc_heat PRIVATE Core)
# target_link_libraries(mcmc_heat PRIVATE LibHeat)
# target_link_libraries(mcmc_heat PRIVATE MyMfem)
# target_link_libraries(mcmc_heat PRIVATE Pardiso)
# target_link_libraries(mcmc_heat PRIVATE LibMlmcmc)

# parallel
# add_executable(mcmc_heat_mpi)
# target_link_libraries(mcmc_heat_mpi PRIVATE Core)
# target_link_libraries(mcmc_heat_mpi PRIVATE LibHeat)
# target_link_libraries(mcmc_heat_mpi PRIVATE MyMfem)
# target_link_libraries(mcmc_heat_mpi PRIVATE Pardiso)
# target_link_libraries(mcmc_heat_mpi PRIVATE LibMlmcmc)
# target_link_libraries(mcmc_heat_mpi PRIVATE LibMlmcmcMpi)
# target_link_libraries(mcmc_heat_mpi PRIVATE "${MPI_C_LIB}" "${MPI_CXX_LIB}")


#################
## MLMCMC Heat ##
#################

# serial
# add_executable(mlmcmc_heat)
# target_link_libraries(mlmcmc_heat PRIVATE Core)
# target_link_libraries(mlmcmc_heat PRIVATE LibHeat)
# target_link_libraries(mlmcmc_heat PRIVATE MyMfem)
# target_link_libraries(mlmcmc_heat PRIVATE Pardiso)
# target_link_libraries(mlmcmc_heat PRIVATE LibMlmcmc)

# parallel
# add_executable(mlmcmc_heat_mpi)
# target_link_libraries(mlmcmc_heat_mpi PRIVATE Core)
# target_link_libraries(mlmcmc_heat_mpi PRIVATE LibHeat)
# target_link_libraries(mlmcmc_heat_mpi PRIVATE MyMfem)
# target_link_libraries(mlmcmc_heat_mpi PRIVATE Pardiso)
# target_link_libraries(mlmcmc_heat_mpi PRIVATE LibMlmcmc)
# target_link_libraries(mlmcmc_heat_mpi PRIVATE LibMlmcmcMpi)
# target_link_libraries(mlmcmc_heat_mpi PRIVATE "${MPI_C_LIB}" "${MPI_CXX_LIB}")

#############
## Testing ##
#############

# Testing, serial
add_executable(unit_tests)

target_link_libraries(unit_tests PRIVATE Core)
target_link_libraries(unit_tests PRIVATE LibHeat)
target_link_libraries(unit_tests PRIVATE LibSparseHeat)
target_link_libraries(unit_tests PRIVATE MyMfem)
target_link_libraries(unit_tests PRIVATE Pardiso)
target_link_libraries(unit_tests PRIVATE GTest)

#

# add sources
add_subdirectory(src)
add_subdirectory(tests)

#
