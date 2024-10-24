## Installation
This note present how to install packages which are required for Channelflow 2.0. First, we neet to install conda and create a new environment named `channelflow`. For each usages, we just do `conda activate channelflow` before do anything regarding Channelflow. Second, Channelflow 2.0 requires some external libraries and we must install them within `channelflow` environment.
- MPI - a standardized and portable message-passing standard designed to function on parallel computing architectures. We can use OpenMPI or MPICH for this.
```bash
conda install conda-forge::openmpi
conda install openmpi-mpicc
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
conda install conda-forge::pthread-stubs
conda install -c conda-forge doxygen