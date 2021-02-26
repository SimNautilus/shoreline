set(classes scaleDirectionFilter)

vtk_module_add_module(scaleDirectionFilters
  # FORCE_STATIC # Using FORCE_STATIC build the vtk module statically into the plugin library, to avoid confusion when loading
  CLASSES ${classes})

find_package(LIBIGL REQUIRED COMPONENTS core)
find_package (Eigen3 3.3 REQUIRED)
vtk_module_include(scaleDirectionFilters PRIVATE "/Users/coreylee/Git/libigl/include/")
vtk_module_link(   scaleDirectionFilters PRIVATE Eigen3::Eigen igl::core)

paraview_add_server_manager_xmls(XMLS  scaleDirectionFilter.xml)
