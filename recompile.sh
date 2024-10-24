#### re-compile all code
cd build
cmake ../channelflow -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=~/channelflow-2.0 -DWITH_HDF5CXX=On
make -j4

#### copy exe file and run it
cd ../works
cp ../build/custom_duc/customtest ./
mpiexec -n 4 ./customtest