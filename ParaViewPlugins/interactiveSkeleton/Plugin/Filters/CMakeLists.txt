set(classes
vtkInteractiveSkeletonFilter
vtkMySkeletonWidget
vtkMySkeletonSource
vtkMySkeletonRepresentation)

vtk_module_add_module(interactiveSkeletonFilters
  FORCE_STATIC
  CLASSES ${classes})

#Find external libraries
# find_package(CGAL REQUIRED)

find_package(LIBIGL REQUIRED igl::core igl::common)
# find_package(Embree 3.12.1 REQUIRED)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
# add_definitions("-DCGAL_EIGEN3_ENABLED") # -DCGAL_HEADER_ONLY

vtk_module_link(interactiveSkeletonFilters
  PUBLIC igl::core igl::common Eigen3::Eigen ${EMBREE_LIBRARIES}
)
vtk_module_include(interactiveSkeletonFilters PUBLIC ${EMBREE_INCLUDES} "/Users/coreylee/Git/libigl_Test/libigl/external/embree/include" "/Users/coreylee/Git/libigl_Test/build/include")

paraview_add_server_manager_xmls(XMLS vtkInteractiveSkeletonFilter.xml)
