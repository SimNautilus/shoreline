set(classes
  vtkcreateSkeletonFilter
)

vtk_module_add_module(createSkeletonFilters
  # FORCE_STATIC # Using FORCE_STATIC build the vtk module statically into the plugin library, to avoid confusion when loading
  CLASSES ${classes})

#Find external package
find_package(LIBIGL REQUIRED)
find_package(CGAL REQUIRED)
find_package(Boost REQUIRED) 

vtk_module_include( createSkeletonFilters
  PRIVATE "/users/jeffhadley/Builds/src/libigl/include"
)
vtk_module_link( createSkeletonFilters PRIVATE CGAL::CGAL Eigen3::Eigen igl::core Boost::boost)

paraview_add_server_manager_xmls(
XMLS  createSkeletonFilter.xml)
