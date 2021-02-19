# this quick script is for building Shoreline on Summit and
# probably other linux based systems that have PHASTA access
# the default options will build successfully on a mac with the requisite libs


module purge
module load intel/17.4
module load impi/17.3
module load mkl/17.3
module load petsc/3.8.0
module load cmake/3.14.1

export PETSC_LIBRARIES=$CURC_PETSC_LIB
export PETSC_INCLUDES=$CURC_PETSC_INC

cd build
cmake \
-DBUILD_WITH_PHASTA=ON \
-DCMAKE_Fortran_COMPILER=mpiifort \
-DBUILD_WITH_VTK=OFF \
-DBUILD_WITH_PARAVIEW=OFF \
-DPHASTA_SRC_DIR=/projects/cone6342/Git/phasta-next \
-DPHASTA_INCOMPRESSIBLE=ON \
-DPHASTA_COMPRESSIBLE=OFF \
-DPHASTA_USE_LESLIB=OFF \
-DPHASTA_USE_SVLS=ON \
-DPHASTA_USE_PETSC=OFF \
-DPHASTA_TESTING=OFF \
-DPETSC_LIBRARIES=$CURC_PETSC_LIB \
-DPETSC_INCLUDES=$CURC_PETSC_INC \
-DCMAKE_Fortran_FLAGS="-O3 -xCORE-AVX2 -nofor_main" \
-DCMAKE_C_FLAGS="-O3 -xCORE-AVX2 -I$PETSC_INCLUDES -L $PETSC_LIBRARIES -lpetsc" \
-DCMAKE_CXX_FLAGS="-O3 -xCORE-AVX2 -I$PETSC_INCLUDES -L $PETSC_LIBRARIES -lpetsc" \
..


#-DLESLIB=/users/matthb2/libles1.5/libles-debianjessie-gcc-ompi.a \
#-DCASES=/users/nelson15/projects/ImmersiveSims/phastaChefTests \
