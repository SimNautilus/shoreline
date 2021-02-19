set(classes
  vtkcreateSkeletonFilter)

vtk_module_add_module(createSkeletonFilters
  # FORCE_STATIC # Using FORCE_STATIC build the vtk module statically into the plugin library, to avoid confusion when loading
  CLASSES ${classes})

vtk_module_include(createSkeletonFilters PRIVATE "/Users/coreylee/Git/Pinocchio/Pinocchio")
vtk_module_link(createSkeletonFilters
  PRIVATE "/Users/coreylee/Git/Pinocchio/Pinocchio/libpinocchio.a")

paraview_add_server_manager_xmls(
XMLS  createSkeletonFilter.xml)
