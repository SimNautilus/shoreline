set(classes
vtkInteractiveSkeletonFilter
vtkMySkeletonWidget
vtkMySkeletonSource
vtkMySkeletonRepresentation)

vtk_module_add_module(interactiveSkeletonFilters
  FORCE_STATIC
  CLASSES ${classes})

#Find external libraries
find_package(LIBIGL REQUIRED igl::core igl::common)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

vtk_module_link(interactiveSkeletonFilters
  PUBLIC igl::core igl::common Eigen3::Eigen ${EMBREE_LIBRARIES}
)

paraview_add_server_manager_xmls(XMLS vtkInteractiveSkeletonFilter.xml)
