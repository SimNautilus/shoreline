set(interfaces)
set(sources
  pqMyPolyLineWidget.cxx
  pqMyPolyLineWidget.h
)

paraview_plugin_add_property_widget(
  KIND GROUP_WIDGET
  TYPE "my_skeleton_widget"
  CLASS_NAME pqMyPolyLineWidget
  INTERFACES property_interfaces
  SOURCES property_sources)
list(APPEND interfaces
  ${property_interfaces})
list(APPEND sources
  ${property_sources})

paraview_add_plugin(HandleWidgets
  VERSION "1.0"
  UI_INTERFACES ${interfaces}
  UI_FILES pqMyPolyLineWidget.ui
  SOURCES ${sources}
  MODULES interactiveSkeletonFilters
  MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/Filters/vtk.module")
