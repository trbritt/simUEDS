cmake_minimum_required(VERSION 3.5)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ./bin)

project(opensc)
enable_language(Fortran)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ./lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ./include)
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
 set(dialect "-Wuninitialized -fbacktrace -ffree-form -ffree-line-length-none -std=f2008 -fimplicit-none")
   #     set(dialect "-Wconversion -ffree-form -ffree-line-length-none   -std=f2008 -fimplicit-none")
 #set(dialect "-ffree-form -ffree-line-length-none   -std=f2008 -fimplicit-none")   
    set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone")
    set(bounds "-check bounds")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-C")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

enable_language(C)
#if(CMAKE_C_COMPILER_ID MATCHES "GNU")
# set(dialect "-Wall ")
#endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wuninitialized -pedantic -Wextra")


set(CMAKE_BUILD_TYPE Debug)

#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_Fortran_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_Fortran_FLAGS}")
#endif()


option(BUILD_KDTREE "Enable building kdtree2" ON)

if(BUILD_KDTREE)
    # add_subdirectory(kdtree/src-f90)
    add_custom_target(
        kdtree
        COMMAND make
        WORKING_DIRECTORY src/mpi/kdtree/src-f90
    )
endif()
#--------------------------
option(GetFFTW "GetFFTW" ON)
if(GetFFTW)

# FFTW2 
include(ExternalProject)
ExternalProject_Add(project_fftw
  URL http://www.fftw.org/fftw-2.1.5.tar.gz
  PREFIX ${CMAKE_CURRENT_BINARY_DIR}/include/fftw
  CONFIGURE_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/include/fftw/src/project_fftw/configure --enable-mpi --prefix=<INSTALL_DIR>
  BUILD_COMMAND make -j 12
  INSTALL_COMMAND make install
  PREFIX=.${CMAKE_CURRENT_BINARY_DIR}/fftw
)

ExternalProject_Get_Property(project_fftw install_dir)
message(STATUS "fftw install_dir: " ${install_dir})
add_library(fftw STATIC IMPORTED)
set_property(TARGET fftw PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libfftw.a)

add_library(fftw_mpi STATIC IMPORTED)
set_property(TARGET fftw_mpi PROPERTY IMPORTED_LOCATION ${install_dir}/lib/libfftw_mpi.a)

add_dependencies(fftw project_fftw)
include_directories(${install_dir}/include)


endif(GetFFTW)
#--------------------------


##find_library (FFTW_LIBRARIES NAMES libdfftw.a)
##message(STATUS "FFTW_LIBRARIES: " ${FFTW_LIBRARIES})

##find_path (FFTW_INCLUDES NAMES fftw)
##message(STATUS "FFTW_INCLUDES: " ${FFTW_INCLUDES})



find_package(MPI REQUIRED)
# set(HDF5_ROOT /home/trbritt/hdf5)
# find_package(HDF5 REQUIRED COMPONENTS Fortran )
# set(HDF5_Fortran_INCLUDE_DIRS /home/trbritt/hdf5-1.12.0/hdf5/include)
# set(HDF5_Fortran_LIBRARIES  -ldl -lm -Wl,-rpath -Wl,/home/trbritt/hdf5-1.12.0/hdf5/lib)
# set(HDF5_PREFER_PARALLEL true)

# message(STATUS "HDF5 Fortran Include: " ${HDF5_Fortran_INCLUDE_DIRS})
# message(STATUS "HDF5 Fortran Libraries: " ${HDF5_Fortran_LIBRARIES})
# message(STATUS "HDF5 Parallel? : " ${HDF5_IS_PARALLEL})
#--------------------------
#
# Preprocessor
#
option(MPIPARALLEL "use MPI" ON)
IF(MPIPARALLEL)
    ADD_DEFINITIONS(-DMPIPARALLEL)
ENDIF(MPIPARALLEL)

option(FFT_FFTW "use fftw" ON)
IF(FFT_FFTW)
    ADD_DEFINITIONS(-DFFT_FFTW)
ENDIF(FFT_FFTW)

#set (CMAKE_Fortran_FLAGS "${CMAKE_C_FLAGS} ${MPI_Fortran_FLAGS}")
set(MPI_COMPILER /home/trbritt/hdf5/bin/h5pfc)

message(STATUS "MPI_Fortran_INCLUDE_PATH: " ${MPI_Fortran_INCLUDE_PATH})
message(STATUS "MPI_Fortran_Libraries: " ${MPI_Fortran_LIBRARIES})

include_directories(SYSTEM ${MPI_Fortran_INCLUDE_PATH} ${HDF5_Fortran_INCLUDE_DIRS} src/mpi/code/ src/mpi/kdtree/src-f90)



#--------------------------
#
# Compile.
#
file(GLOB_RECURSE sources_sc src/mpi/code/fft/*.h  src/mpi/code/fft/*.c 
                             src/mpi/kdtree/src-f90/kdtree2.f90
                             src/mpi/code/data_movement_mod.f90
                             src/mpi/code/decomposition_mod.f90
                             src/mpi/code/deposition_mod.f90
                             src/mpi/code/fft_mod.f90
                             src/mpi/code/hdf5utils_mod.f90
                             src/mpi/code/fieldutils_mod.f90
                             src/mpi/code/numerical_distributions_mod.f90
                             src/mpi/code/open_spacecharge_mod.f90
                             src/mpi/code/opensc.f90
                             src/mpi/code/rf_mod.f90
                             src/mpi/code/propagate_mod.f90
                            )
set(CMAKE_Fortran_COMPILER /home/trbritt/hdf5-1.12.0/hdf5/bin/h5pfc)

add_executable(spacecharge ${sources_sc})
# add_library(spacecharge ${sources_sc})
target_link_libraries(spacecharge fftw fftw_mpi ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES} )

file(GLOB_RECURSE objects_sc build/*.mod)
file(COPY ${objects_sc} DESTINATION ../../../include)

# Link FFTW
## target_link_libraries(test_opensc ${FFTW_LIBRARIES} ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES} )

# enable_language(CXX)
# set(MPI_C_LIBRARIES ${MPI_C_LIBRARIES} PARENT_SCOPE)
# set(MPI_Fortran_LIBRARIES ${MPI_Fortran_LIBRARIES} PARENT_SCOPE)
# set(sources_sc ${sources_sc} PARENT_SCOPE)
# set(fftw ${fftw} PARENT_SCOPE)
# set(FFTW_LIBRARIES ${FFTW_LIBRARIES} PARENT_SCOPE)

#add_executable(test_linking test_linking.cpp ${sources_sc})
#target_link_libraries(test_linking fftw fftw_mpi ${MPI_C_LIBRARIES} ${MPI_Fortran_LIBRARIES})