cmake_minimum_required(VERSION 3.14.0)

################################################################################
#                                Compiler Flags                                #
################################################################################
# set compiler
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

# set build type to release if not specified (shouldn't be an issue)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# add compiler flags for cxx17, wAll and -O3 optimization
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

################################################################################
#                              Define ShorelineDDG                             #
################################################################################
# Set the version number
set(SHORELINE_MAJOR_VERSION 0)
set(SHORELINE_MINOR_VERSION 1)
set(SHORELINE_PATCH_VERSION 0)
set(PACKAGE_VERSION
  ${SHORELINE_MAJOR_VERSION}.${SHORELINE_MINOR_VERSION}.${SHORELINE_PATCH_VERSION}
)

# setup an optional install directory
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/install")

#define the shoreline project
project(ShorelineDDG)

# make compiling with VTK optional
option(BUILD_WITH_VTK "Build Shoreline with VTK functionality" ON)
option(BUILD_WITH_PARAVIEW "Build Shoreline with Paraview Catalyst features" ON)

# set directory information
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeFiles")
INCLUDE_DIRECTORIES(src)

# Set the header files for Shoreline
set(SHRLN_HEADERS
  src/MeshCourier.h
  src/MeshMapper.h
  src/MeshMill.h
  src/MeshMillDDG.h
  src/ShorelineMPI.h
  src/MeshCourierClient.h
  src/MeshCourierServer.h
)
# set the source files for Shoreline
set(SHRLN_SOURCE
  src/ShorelineMPI.cpp
  src/MeshMill.cpp
  src/MeshCourier.cpp
  src/MeshMapper.cpp
  src/MeshMill_SetupSolver.cpp
  src/MeshMill_ElemFormation.cpp
  src/MeshMill_ElemAssembly.cpp
  src/MeshMill_SolveMeshMotion.cpp
  src/MeshMill_ViewMatrix.cpp
  src/MeshMillDDG_ElemAssembly.cpp
  src/MeshMillDDG_ElemFormation.cpp
  src/MeshMillDDG_SetupSolver.cpp
  src/MeshMillDDG_SolveMeshMotion.cpp
  src/MeshCourierClient.cpp
  src/MeshCourierServer.cpp
)
# define the ShorelineDDG library
add_library(ShorelineDDG STATIC "${SHRLN_SOURCE}")

# define the runner executable
add_executable(runner.exe main.cpp)
add_executable(chefRunner.exe chefRunner.cpp)
# the SCOREC Numbering objects are a little weird, in that defining them
# alone does all the linking to the mesh that is required. As such, the compiler
# thinks the apf::Numbering variables are unused. So we suppress that warning
target_compile_options(runner.exe PRIVATE "-Wno-unused-variable")
target_compile_options(chefRunner.exe PRIVATE "-Wno-unused-variable")

################################################################################
#                            Add External Libraries                            #
################################################################################
# add PETSc
FIND_PACKAGE(PETSc REQUIRED)
INCLUDE_DIRECTORIES(${PETSC_INCLUDES})
ADD_DEFINITIONS(${PETSC_DEFINITIONS})

# add MKL
find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIRS})

# add SCOREC/core
find_package(SCOREC REQUIRED)

# create list of external libraries
if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(SHORELINE_LIBS ${MKL_LIBRARIES} ${PETSC_LIBRARIES}
    SCOREC::apf SCOREC::mds SCOREC::gmi SCOREC::core SCOREC::ph
  )
else()
set(SHORELINE_LIBS -Wl,--start-group ${MKL_LIBRARIES} -Wl,--end-group ${PETSC_LIBRARIES}
  SCOREC::apf SCOREC::mds SCOREC::gmi SCOREC::core SCOREC::ph
)
endif()

# add VTK
if(BUILD_WITH_VTK)
  message( "building with VTK")
  add_definitions(-DBUILD_WITH_VTK)
  find_package(VTK COMPONENTS CommonCore CommonDataModel IOXML)
  list(APPEND SHORELINE_LIBS ${VTK_LIBRARIES}) #VTK::PythonUsed 
endif()

# add ParaView Catalyst
if(BUILD_WITH_PARAVIEW)
  add_definitions(-DUSE_CATALYST)
  find_package(ParaView REQUIRED COMPONENTS Catalyst VTKExtensionsCore PythonCatalyst)
  list(APPEND SHORELINE_LIBS
    ParaView::Catalyst
    ParaView::VTKExtensionsCore
    ParaView::PythonCatalyst
  )
endif()

# link external libraries to ShorelineDDG
target_link_libraries(ShorelineDDG ${SHORELINE_LIBS})

################################################################################`
#                         Attach PHASTA... Maybe...                            #
################################################################################
option(BUILD_WITH_PHASTA off)

if(BUILD_WITH_PHASTA)
  enable_language(Fortran)
  set(PHASTA_SRC_DIR phasta CACHE FILEPATH "path to phasta source code")
  set(PHASTA_CHEF_ENABLED ON)
  add_subdirectory(${PHASTA_SRC_DIR} ${CMAKE_BINARY_DIR}/phasta)
  find_package(phastaIC PATHS ${CMAKE_BINARY_DIR})
  message(STATUS "PHASTAIC_LIBS ${PHASTAIC_LIBS}")
  add_executable(phastaChefRunner.exe phastaChefRunner.cpp)
  set_target_properties(phastaChefRunner.exe PROPERTIES LINKER_LANGUAGE Fortran)
  target_compile_options(phastaChefRunner.exe PRIVATE "-Wno-unused-variable")
  target_include_directories(phastaChefRunner.exe PRIVATE ${PHASTAIC_INCLUDE_DIRS})
  target_link_libraries(phastaChefRunner.exe ShorelineDDG "dl")
  target_link_libraries(phastaChefRunner.exe ${PHASTAIC_LIBS})
  target_link_libraries(phastaChefRunner.exe SCOREC::core)

  add_executable(shPhasta.exe shPhasta.cpp)
  set_target_properties(shPhasta.exe PROPERTIES LINKER_LANGUAGE Fortran)
  target_compile_options(shPhasta.exe PRIVATE "-Wno-unused-variable")
  target_include_directories(shPhasta.exe PRIVATE ${PHASTAIC_INCLUDE_DIRS})
  target_link_libraries(shPhasta.exe ShorelineDDG "dl")
  target_link_libraries(shPhasta.exe ${PHASTAIC_LIBS})
  target_link_libraries(shPhasta.exe SCOREC::core)
endif()

################################################################################
#                 Install Shoreline DDG to an install location                 #
################################################################################
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/install")
install(TARGETS ShorelineDDG
  EXPORT ShorelineTargets
  DESTINATION lib
)
install(FILES ${SHRLN_HEADERS} DESTINATION include)

install(EXPORT ShorelineTargets
  FILE ShorelineTargets.cmake
  DESTINATION lib/cmake/Shoreline
)

include(CMakePackageConfigHelpers)
# generate the config file that includes the exports
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/ShorelineConfig.cmake"
  INSTALL_DESTINATION "lib/cmake/example"
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO
)
# generate the version file for the config file
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/ShorelineConfigVersion.cmake"
  VERSION "${PACKAGE_VERSION}"
  COMPATIBILITY AnyNewerVersion
)

# install the configuration file
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/ShorelineConfig.cmake
  DESTINATION lib/cmake/Shoreline
)

export(EXPORT ShorelineTargets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/ShorelineTargets.cmake"
)
################################################################################
#                             activate VTK modules                             #
################################################################################
# vtk_module_autoinit is needed... apparently
# vtk_module_autoinit( TARGETS ShorelineDDG MODULES ${VTK_LIBRARIES})

################################################################################
#                          Link ShorelineDDG to Runner                         #
################################################################################
target_link_libraries( runner.exe ShorelineDDG "dl")
target_link_libraries( chefRunner.exe ShorelineDDG "dl" )
target_link_libraries(chefRunner.exe SCOREC::core)

add_executable(shClient.exe shClient.cpp)
add_executable(shServer.exe shServer.cpp)
target_link_libraries( shClient.exe ShorelineDDG "dl")
target_link_libraries( shServer.exe ShorelineDDG "dl")
target_compile_options(shServer.exe PRIVATE "-Wno-unused-variable")
