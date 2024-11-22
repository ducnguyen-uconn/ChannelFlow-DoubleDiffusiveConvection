## Installation
This note presents how to install packages required for Channelflow 2.0 on university HPC. First, we neet to install miniconda and create a new environment named `channelflow`. For each usage, we just do `conda activate channelflow` before doing anything regarding Channelflow. Second, Channelflow 2.0 requires some external libraries and we must install them within `channelflow` environment. Pls read this [instruction](https://kb.uconn.edu/space/SH/26079723879/Miniconda+Environment+Set+Up) for installing miniconda.
- MPI - a standardized and portable message-passing standard designed to function on parallel computing architectures. We can use OpenMPI or MPICH for this.
```bash
conda install conda-forge::openmpi
# conda install openmpi-mpicc
mpiexec --version
```
- Eigen3 - a library for linear algebra: matrices, vectors, numerical solvers, and related algorithms
```bash
conda install conda-forge::eigen
```
- FFTW3 - a library for computing the discrete Fourier transform (DFT) in one or more dimensions, of arbitrary input size, and of both real and complex data
```bash 
conda install conda-forge::fftw
```
- FFTW3-MPI - a sub-library of FFTW3 package for MPI
```bash 
conda install cryoem/label/archive::fftw-mpi
```
- NetCDF - a set of software libraries and self-describing, machine-independent data formats that support the creation, access, and sharing of array-oriented scientific data.
```bash 
conda install conda-forge::netcdf4
```
- Cmake
```bash
conda install conda-forge::cmake
cmake --version
```
- hdf5
```bash
conda install -c conda-forge hdf5
```
Maybe you also need these
```bash
conda install conda-forge::pthread-stubs
conda install -c conda-forge doxygen
```

Also, for each compiler, you may are required to change MPICXX compiler to specific compiler of the system, eg HPC. On local machine, common mpicxx compiler is set `/usr/bin/mpicxx`. On university HPC, you are not accessed to use mpicxx from miniconda when compiling global executable files. Thus, you can indicate `CMAKE_CXX_COMPILER` by specific compiler of system. For examples, if you're using UCONN HPC, you need to use compiler from `/gpfs/sharedfs1/admin/hpc2.0/apps/openmpi/5.0.5/bin/mpicxx`. You can find this via `which mpicxx`.

A sample installation, with all features enabled, might look like this:
```bash 
cmake ../channelflow -DCMAKE_CXX_COMPILER=/gpfs/sharedfs1/admin/hpc2.0/apps/openmpi/5.0.5/bin/mpicxx -DWITH_DDC=ON -DWITH_NSOLVER=ON -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/chflow -DCMAKE_CXX_FLAGS_RELEASE:STRING=" -fPIC -lfftw3 -lm -Wno-unused-variable -Wno-deprecated " -DWITH_SHARED=OFF -DWITH_HDF5CXX=OFF
make -j16
```