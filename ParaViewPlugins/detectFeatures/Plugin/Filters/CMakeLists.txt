set(classes DetectFeaturesFilter)

vtk_module_add_module(detectFeaturesFilters
  # FORCE_STATIC # Using FORCE_STATIC build the vtk module statically into the plugin library, to avoid confusion when loading
  CLASSES ${classes})

find_package(CGAL REQUIRED)
vtk_module_link( detectFeaturesFilters PRIVATE  CGAL::CGAL )

paraview_add_server_manager_xmls(XMLS  DetectFeaturesFilter.xml)

message("hey we got here too")
