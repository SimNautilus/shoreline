set(classes
  vtkcreateSkeletonFilter
)

vtk_module_add_module(createSkeletonFilters
  # FORCE_STATIC # Using FORCE_STATIC build the vtk module statically into the plugin library, to avoid confusion when loading
  CLASSES ${classes})

#Find CGAL
find_package(LIBIGL REQUIRED)
find_package(CGAL REQUIRED)
find_package(Boost REQUIRED) 
# find_package (Eigen3 3.3 REQUIRED NO_MODULE)
# add_definitions("-DCGAL_EIGEN3_ENABLED")#-DCGAL_HEADER_ONLY ***Commented this out***

vtk_module_include( createSkeletonFilters
  PRIVATE "/users/jeffhadley/Builds/src/libigl/include"
)
vtk_module_link( createSkeletonFilters PRIVATE CGAL::CGAL Eigen3::Eigen igl::core Boost::boost)
	#PRIVATE CGAL::CGAL Eigen3::Eigen igl::core ${Boost_LIBRARIES})
paraview_add_server_manager_xmls(
XMLS  createSkeletonFilter.xml)
