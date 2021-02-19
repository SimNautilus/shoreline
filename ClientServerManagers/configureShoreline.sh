# this quick script is for building Shoreline on the viz nodes and
# probably other linux based systems that have PHASTA access
# the default options will build successfully on a mac with the requisite libs

cd build
cmake \
-DBUILD_WITH_PHASTA=ON \
-DBUILD_WITH_VTK=OFF \
-DBUILD_WITH_PARAVIEW=ON \
-DPHASTA_SRC_DIR=/users/nelson15/Git/covizPhasta/phasta-next \
-DPHASTA_INCOMPRESSIBLE=ON \
-DPHASTA_COMPRESSIBLE=OFF \
-DPHASTA_USE_LESLIB=OFF \
-DPHASTA_USE_SVLS=ON \
-DPHASTA_USE_PETSC=OFF \
-DParaView_DIR=/users/nelson15/LIBS/pvBuild/install/lib/cmake/paraview-5.8 \
-DVTK_DIR=/users/nelson15/LIBS/pvBuild/install/lib/cmake/paraview-5.8/vtk \
-DPHASTA_USE_CATALYST=ON \
-DUSE_CATALYST=ON \
-DPHASTA_ADAPTOR_LIB=/users/nelson15/Git/covizPhasta/build_adaptor/libPhastaAdaptor.a \
-DCMAKE_PREFIX_PATH=/usr/local/qt/5.11.0/lib/cmake/Qt5/ \
-DCMAKE_BUILD_TYPE=Release \
-DLESLIB=/users/matthb2/libles1.5/libles-debianjessie-gcc-ompi.a \
-DPHASTA_TESTING=OFF \
-DCASES=/users/nelson15/projects/ImmersiveSims/phastaChefTests \
..
