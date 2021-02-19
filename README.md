# shoreline
Tools to enable _in-situ_ interactive geometry modification of simulation domains

# Getting Setup

This system is architected to operate as a server and a client. Typically, the server will run on a remote computing resource while the client will run on the user's workstation or laptop. However, there is nothing preventing both the server and the client from running on the same hardware.

## ParaView Configuration
ParaView is the driving software for user interactivity and _in-situ_ visualization. It is important that the ParaView installed on the client workstation matches that of the ParaView install on the server compute resource.

ParaView can be cloned from https://gitlab.kitware.com/paraview/paraview. Version 5.8.1 is the only version tested, but later versions will likely work as well.
A build guide is available here https://gitlab.kitware.com/paraview/paraview/-/blob/master/Documentation/dev/build.md

Python is required for the catalyst components used in visualization. You will need CMake, Qt, Python and a set of compilers installed. Few optional components are required for this system. It is important to turn on the option PARAVIEW_BUILD_SHARED_LIBS.

ParaView must be configured to the same version on the client computer and the server computer to enable in-situ visualization via Catalyst.

## Building the ParaView Plugins
Currently, each plugin is built standalone. That is, each plugin will need to be setup with a build directory and then separately built. In the future, all of them could be built at once with the help of CMake. All plugins will require setting the environment variable ParaView_DIR, VTK_DIR, and Qt5_DIR.
Mileage will vary, but on Corey's macbook, these are set to:
Qt5_DIR=/Users/coreylee/Qt/5.13.1/clang_64/lib/cmake/Qt5
ParaView_DIR=/Users/coreylee/Git/paraview/install_5.6.0/lib/cmake/paraview-5.8
VTK_DIR=/Users/coreylee/Git/paraview/install_5.6.0/lib/cmake/paraview-5.8/vtk

A number of plugins require other external libraries for operation. The most common ones are [libIGL](https://libigl.github.io) and [CGAL](https://www.cgal.org). Each of these libraries have install instructions on their website. libIGL is quite straightforward to install, CGAL can be tricky. If things go well, the easiest course of installation for CGAL is through a package manager as recommended on their website (brew for macOS and apt-get for Linux).

To build a plugin, from a plugin's top level directory:
mkdir build
cd build
cmake ..
make -j

Note, of course, the build directory can be wherever you'd like it to be.

Plugins must then be turned on in ParaView using the plugin manager. The first time this is done, you will need to navigate to the shared library generated in the plugin build. You can then ask ParaView to auto-load the plugins and remember where they all are. If this is done correctly, the plugins should show up in the Filters/Nautilus dropdown menu.

## Building the Client-Server Managers
The client and server managers are built from the same project with options managing whether PHASTA is included (so that the client computer does not need access to the simulation code). To build this, the minimum requirements are the libraries [PETSc](https://www.mcs.anl.gov/petsc/), [SCOREC core](https://github.com/SCOREC/core), [MKL](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html), and VTK. VTK should be installed in the ParaView installation. However, the use of VTK in these apps is decoupled from ParaView proper, and can thus be from a different version of ParaView than will be used for catalyst. The VTK library in this context is just used by the Client app to shuttle data to ParaView through file I/O.

There are two options available that default to off. BUILD_WITH_PARAVIEW, BUILD_WITH_PHASTA, BUILD_WITH_VTK. The first two correspond to Catalyst and PHASTA capabilities as related to the Server app. These two are not required by the Client. BUILD_WITH_VTK corresponds to the above mentioned behavior of the Client app shuttling data to and from ParaView via file I/O. This setting is required to be on for the client and not required for the server. These different configurations should be considered based on where you are installing the client-server managers. As they share significant amounts of code, it makes sense to keep the together in one place, but the Client app has a very different setup to the Server app. With PHASTA enabled, all of the relevant PHASTA CMake options are required to be passed in as well. The above libraries and options do not interact with PHASTA. That is, you can turn off PETSc in PHASTA and it will still be required by the Server app (for volume mesh deformation).

On the server, from the top level of ClientServerManagers:
mkdir build
cd build
cmake -DBUILD_WITH_PARAVIEW=ON -DBUILD_WITH_PHASTA=ON ..
make -j

on the client, from the top level of ClientServerManagers:
mkdir build
cd build
cmake -DBUILD_WITH_VTK=ON ..
make -j

# Connecting Devices for Runtime
The system requires 2 ssh tunnels, one for geometry data, and one for Catalyst visualization. For example, on the client, boot up 3 terminals.
1. run ssh -L22222:viz003:22222 user@jumpgate-phasta.colorado.edu
   run mpirun -np X pvserver
   (This will connect to ParaView on the client for live visualization)
2. run ssh -L 8080:viz003:8080 user@jumpgate-phasta.colorado.edu (currently the 8080 port is hard coded for this)
   (this will be where the Server app is run. See section below)
3. Navigate to where you would like to run the Client application and wait until the Server app is running

# Running the System
A simulation directory should be setup as if it were prepared for Chef. That is, you need a mesh (partitioned .smb file(s)), a model (I stick with .dmg files to avoid needing Simmetrix in the build tree), an adapt.inp file, (i prefer to use a bcs.spj for BC specification), and a solver.inp. In the solver.inp, the total timesteps is NOT  the total time steps of the simulation, and has been hijacked to now represent the number of timesteps simulated between geometry deformations. from that directory, run the PHASTA server with mpi:

mpirun -np X /path/to/shPhasta.exe

This will then hang and wait for a client to connect.

On the client computer, run ./shClient.exe

When the Client is initialized, it will let you know if it successfully connected. Upon successful connection, the Server app will kick off the simulation and it will begin buzzing away.

The client ought to have at this point downloaded the surface mesh and wrote it to a .vtp file for manipulation in ParaView. Open a ParaView session and open the surface mesh. You can then modify the mesh with the installed plugins. When you are happy with the modifications, run the final plugin exportDisplacements and point it to the folder where the original vtp was saved by the Client (where you ran the app from). Then from terminnal 3, press the 's' key to tell the Client the mesh has been modified and the modification has been exported from ParaView. The client will read the displacement file and then hang as it tries to send the data on to the server. The server will eventually step out of the simulation routine and check if the client has data ready, and then initialize communication of the displacement field. It will then automatically deform the mesh.
