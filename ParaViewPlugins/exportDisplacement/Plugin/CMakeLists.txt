message("did you get here?")
paraview_add_plugin(exportDisplacementFilter
  VERSION "1.0"
  MODULES exportDisplacementFilters
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Filters/vtk.module")

message("Hey we ran paraview_add_plugin")

# option(BUILD_TESTING "Build Plugin Testing" OFF)
# if (BUILD_TESTING AND BUILD_SHARED_LIBS)
#   add_subdirectory(Testing)
# endif()
