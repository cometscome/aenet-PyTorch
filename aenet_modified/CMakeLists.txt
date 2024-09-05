cmake_minimum_required(VERSION 3.15)
project(aenet Fortran)

find_package(MPI    REQUIRED)

option(MKLUSE "Use MKL or not" on)
option(ACCUSE "Use MacOs Accelerate framework" OFF)

set(aenet_PATH "${CMAKE_CURRENT_SOURCE_DIR}/")

message(STATUS "extracting...")
execute_process(
      COMMAND ${CMAKE_COMMAND} -E tar xvf "${aenet_PATH}/lib/Lbfgsb.3.0.tar.gz"
      WORKING_DIRECTORY "${aenet_PATH}/lib/"
)
message(STATUS "extracted.")

set(MKL_INTERFACE lp64)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Fujitsu")
else()
  if(MKLUSE)
    find_package(MKL REQUIRED)
  else()
    find_package(LAPACK REQUIRED)
  endif()
endif()


include_directories(SYSTEM ${MPI_C_INCLUDE_PATH} ${MPI_CXX_INCLUDE_PATH})
add_definitions(-DPARALLEL)



add_library(lbfgsb STATIC )


execute_process(COMMAND ${CMAKE_Fortran_COMPILER} --version
                OUTPUT_VARIABLE GCC_VERSION_OUTPUT)
string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" GCC_VERSION ${GCC_VERSION_OUTPUT})

if(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    message(STATUS "Using Intel Fortran Compiler")
    
    add_compile_options(-some-intel-specific-flag)

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    message(STATUS "Using GCC Fortran Compiler")
    
    execute_process(COMMAND ${CMAKE_Fortran_COMPILER} --version
                    OUTPUT_VARIABLE GCC_VERSION_OUTPUT)
    string(REGEX MATCH "[0-9]+\\.[0-9]+\\.[0-9]+" GCC_VERSION ${GCC_VERSION_OUTPUT})

    if(GCC_VERSION VERSION_GREATER_EQUAL "10.0.0")
        message(STATUS "GCC version >= 10, adding -fallow-argument-mismatch flag")
        add_compile_options(-fallow-argument-mismatch)
    else()
        message(STATUS "GCC version < 10, not adding -fallow-argument-mismatch flag")
    endif()
endif()









add_custom_command(
  TARGET lbfgsb
  POST_BUILD
  COMMAND cp  $<TARGET_FILE:lbfgsb>  ${CMAKE_SOURCE_DIR}
)

set(LBFGSB ${CMAKE_CURRENT_SOURCE_DIR}/lib/Lbfgsb.3.0)

target_sources(lbfgsb
PRIVATE
${LBFGSB}/blas.f
${LBFGSB}/lbfgsb.f
${LBFGSB}/linpack.f
${LBFGSB}/timer.f
)



add_library(aenet STATIC)
target_link_libraries(aenet PRIVATE  lbfgsb)

set(aenet_dir ${aenet_PATH}/src_modified)

set(SOURCES
  ${aenet_dir}/ext/chebyshev.f90 
  ${aenet_dir}/ext/feedforward.f90  
  ${aenet_dir}/ext/io.f90  
  ${aenet_dir}/ext/lclist.f90  
  ${aenet_dir}/ext/sfbasis.f90  
  ${aenet_dir}/ext/sortlib.f90  
  ${aenet_dir}/ext/symmfunc.f90  
  ${aenet_dir}/ext/timing.f90
    ${aenet_dir}/ext/unittest.f90  
    ${aenet_dir}/ext/xsflib.f90
    ${aenet_dir}/aenet.f90 
  ${aenet_dir}/aeio.f90
  ${aenet_dir}/constants.f90
  ${aenet_dir}/geometry.f90
  ${aenet_dir}/input.f90
  ${aenet_dir}/optimize.f90
  ${aenet_dir}/parallel.F90
  ${aenet_dir}/potential.f90
  ${aenet_dir}/random.f90
  ${aenet_dir}/sfsetup.f90
  ${aenet_dir}/trainset.f90
  ${aenet_dir}/pytorchoutput.f90
  ${aenet_dir}/pytorchoutput_MPI.f90
  ${aenet_dir}/sub_generate.f90
  ${aenet_dir}/sub_trainbin2ASCII.f90
  ${aenet_dir}/sub_nnASCII2bin.f90
  ${aenet_dir}/trainset_MPI.f90
  ${aenet_dir}/sub_generate_MPI.f90
)


add_executable(aenet_generate.x ${aenet_dir}/generate.f90 ${SOURCES})
add_executable(aenet_generate_MPI.x ${aenet_dir}/generate_MPI.F90 ${SOURCES})


if(NOT MATHLIBDONE)
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "Fujitsu")
    target_link_libraries(aenet_generate.x PRIVATE lbfgsb ${LAPACK_LIBRARIES} )
    target_link_libraries(aenet_generate_MPI.x PRIVATE lbfgsb ${LAPACK_LIBRARIES} )

  else()
    if(MKLUSE)
        target_link_libraries(aenet_generate.x PRIVATE  lbfgsb MKL::MKL )
        target_link_libraries(aenet_generate_MPI.x PRIVATE  lbfgsb MKL::MKL )
    else()
    	target_link_libraries(aenet_generate.x PRIVATE lbfgsb ${LAPACK_LIBRARIES} )
      target_link_libraries(aenet_generate_MPI.x PRIVATE lbfgsb ${LAPACK_LIBRARIES} )
    endif()
      
  endif()
else()
  target_link_libraries(aenet_generate.x PRIVATE lbfgsb )
  target_link_libraries(aenet_generate_MPI.x PRIVATE lbfgsb )
endif() 

target_sources(aenet
PRIVATE
${SOURCES}
)



if(MPI_COMPILE_FLAGS)
  set_target_properties(aenet  PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")

endif()


include_directories(SYSTEM ${MPI_C_INCLUDE_PATH} ${MPI_CXX_INCLUDE_PATH})

target_link_libraries(aenet PRIVATE ${MPI_LIBRARIES})
target_include_directories(aenet PRIVATE ${MPI_INCLUDE_PATH})

target_link_libraries(aenet PRIVATE MPI::MPI_Fortran)

target_compile_definitions(aenet PUBLIC PARALLEL )



