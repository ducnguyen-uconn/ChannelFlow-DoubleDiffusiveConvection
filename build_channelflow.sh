sudo apt-get update
sudo apt-get install libopenmpi-dev
sudo apt install libeigen3-dev
sudo apt-get install libfftw3-dev
sudo apt-get install libfftw3-mpi-dev
sudo apt install netcdf-bin libnetcdff-dev

git clone https://github.com/epfl-ecps/channelflow.git
mkdir build
cd build
cmake ../channelflow -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/user/local -DCMAKE_CXX_FLAGS_RELEASE:STRING=" -fPIC -lfftw3 -lm " -DWITH_HDF5CXX=On
make -j24
sudo make install
