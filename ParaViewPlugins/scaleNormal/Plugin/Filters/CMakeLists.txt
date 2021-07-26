set(classes scaleNormalFilter)

vtk_module_add_module(scaleNormalFilters
  # FORCE_STATIC # Using FORCE_STATIC build the vtk module statically into the plugin library, to avoid confusion when loading
  CLASSES ${classes})

find_package(LIBIGL REQUIRED COMPONENTS core)
find_package (Eigen3 3.3 REQUIRED)
vtk_module_include(scaleNormalFilters PRIVATE "/users/jeffhadley/Builds/src/libigl/include/")
vtk_module_link(   scaleNormalFilters PRIVATE Eigen3::Eigen igl::core)

paraview_add_server_manager_xmls(XMLS  scaleNormalFilter.xml)
