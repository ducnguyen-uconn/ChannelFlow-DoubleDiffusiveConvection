# Instructions for Channelflow 2.0
Channelflow 2.0 is version 2.0 of [Channelflow](https://www.channelflow.ch/) software based on language C++ and parallelized for numerical analysis of the incompressible flow in channel geometries. Channelflow 2.0 is developed by the Emergent Complexity in Physical Systems Laboratory ([ECPS](https://ecps.epfl.ch/)) 

## Installation
Channelflow 2.0 requires some external libraries:
- MPI - a standardized and portable message-passing standard designed to function on parallel computing architectures. We can use OpenMPI or MPICH for this.
```bash
sudo apt-get install libopenmpi-dev
```
- Eigen3 - a library for linear algebra: matrices, vectors, numerical solvers, and related algorithms
```bash
sudo apt-get install libeigen3-dev
```
- FFTW3 - a library for computing the discrete Fourier transform (DFT) in one or more dimensions, of arbitrary input size, and of both real and complex data
```bash 
sudo apt-get install libfftw3-dev 
```
- FFTW3-MPI - a sub-library of FFTW3 package for MPI
```bash 
sudo apt-get install libfftw3-mpi-dev
```
- NetCDF - a set of software libraries and self-describing, machine-independent data formats that support the creation, access, and sharing of array-oriented scientific data.
```bash 
sudo apt install netcdf-bin libnetcdff-dev
```

Now, let us build the Channelflow 2.0. We first need to clone official source code via [GitHub](https://github.com/epfl-ecps/channelflow).
```bash 
git clone https://github.com/epfl-ecps/channelflow.git
```
and build it by using syntax as follows
```bash 
mkdir build
cd build
cmake PATH_TO_SOURCE -DCMAKE_BUILD_TYPE=<debug/release> (configuration options)
make -j
make install
```

Sometimes, we can get a few errors. Here, I will list some common errors, which I got when installing, and how to fix it.

- If get errors regarding `undefined reference to "fftw_..."`, add below flag after `cmake`
```bash  
-DCMAKE_CXX_FLAGS_RELEASE:STRING=" -lfftw3 -lm "
```
- If get errors regarding `fPIC`, add `-fPIC` in above flag or rebuild FFTW with flag `--enable-shared`
```bash  
-DCMAKE_CXX_FLAGS_RELEASE:STRING=" -fPIC "
```
- If get errors regarding library finding (e.g. `libfftw3.a`, `libfftw3_mpi.a`) when compiling. We need to rebuild FFTW3 by hand with both flags `--enable-mpi` and `--enable-shared` after `./configure`. Pls read detailed instructions for installation [here](https://www.fftw.org/fftw3_doc/FFTW-MPI-Installation.html).

Channelflow supports, beneath other standard cmake flags, the following options


|Option                   | Values  | Default   | Description                                                       |
|:------------------------|:--------|:----------|:------------------------------------------------------------------|
|`-DCMAKE_INSTALL_PREFIX` | path    | usr/local | Installation path for make install                                |
|`-DUSE_MPI`              | ON/OFF  | ON        | Enable MPI                                                        |
|`-DWITH_SHARED`          | ON/OFF  | ON        | build shared channelflow and nsolver libraries                    |
|`-DWITH_STATIC`          | ON/OFF  | OFF       | build static libraries (also enables linking to static libraries) |
|`-DWITH_PYTHON`          | ON/OFF  | OFF       | build a python wrapper for flowfields, disabled by default because it requires boost-python |
|`-DWITH_HDF5CXX`         |  ON/OFF | OFF       | enable legacy .h5 file format (using HDF5 C++)                    |


A sample installation, with all features enabled, might look like this:
```bash 
cmake ../channelflow -DCMAKE_BUILD_TYPE=release -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx -DCMAKE_CXX_FLAGS_RELEASE:STRING=" -fPIC -lfftw3 -lm " -DCMAKE_INSTALL_PREFIX=/usr/local -DWITH_PYTHON=ON -DWITH_HDF5CXX=ON
make -j4
sudo make install
```