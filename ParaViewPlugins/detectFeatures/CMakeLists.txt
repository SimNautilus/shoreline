cmake_minimum_required(VERSION 3.16)
project(detectFeaturesFilter)

set(CMAKE_BUILD_TYPE "Release")

find_package(ParaView REQUIRED)

option(BUILD_SHARED_LIBS "Build shared libraries" ON)
# enable_testing()

# Required for python testing
# set(_vtk_build_TEST_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Testing/Temporary)

include(GNUInstallDirs)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

set("_paraview_plugin_default_${CMAKE_PROJECT_NAME}" ON)
paraview_plugin_scan(
  PLUGIN_FILES      "${CMAKE_CURRENT_SOURCE_DIR}/Plugin/paraview.plugin"
  PROVIDES_PLUGINS  plugins
  REQUIRES_MODULES  required_modules
  ENABLE_BY_DEFAULT ON)

message("the plugins are set to ${plugins}")

foreach (module IN LISTS required_modules)
  if (NOT TARGET "${module}")
    message(#FATAL_ERROR
      "Skipping example ${CMAKE_PROJECT_NAME}: Missing required module: "
      "${module}")
    return ()
  endif ()
endforeach ()

message("CMAKE_INSTALL_BINDIR = ${CMAKE_INSTALL_BINDIR}")
message("CMAKE_INSTALL_LIBDIR = ${CMAKE_INSTALL_LIBDIR}")
message("PARAVIEW_PLUGIN_SUBDIR = ${PARAVIEW_PLUGIN_SUBDIR}")
paraview_plugin_build(
  RUNTIME_DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY_DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  LIBRARY_SUBDIRECTORY "${PARAVIEW_PLUGIN_SUBDIR}"
  PLUGINS ${plugins})
