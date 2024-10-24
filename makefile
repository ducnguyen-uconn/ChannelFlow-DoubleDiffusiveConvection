CXX = /usr/bin/mpicxx
# CXX = /home/jms24002/miniforge3/envs/channelflow/bin/mpicxx
# CXX = /gpfs/sharedfs1/admin/hpc2.0/apps/openmpi/5.0.5/bin/mpicxx
TYPE = release
PREFIX = /user/local/
# PREFIX = ~/chflow/
CFLAGS = " -fPIC -lfftw3 -lm "
# CFLAGS = ""

DDC = ON
NSOLVER = ON
HDF5 = OFF

.PHONY: clean build
update:
	sudo apt-get update
	sudo apt-get install libopenmpi-dev
	sudo apt install libeigen3-dev
	sudo apt-get install libfftw3-dev
	sudo apt-get install libfftw3-mpi-dev
	sudo apt install netcdf-bin libnetcdff-dev

clone:
	rm -rf channelflow
	git clone https://github.com/epfl-ecps/channelflow.git

cloneilc:
	rm -rf channelflow
	git clone --single-branch --branch module/ilc https://github.com/epfl-ecps/channelflow.git

build:
	mkdir -p build
	cd build;\
	cmake ../channelflow -DCMAKE_CXX_COMPILER=$(CXX) -DWITH_DDC=$(DDC) -DWITH_NSOLVER=$(NSOLVER) -DCMAKE_BUILD_TYPE=$(TYPE) -DCMAKE_INSTALL_PREFIX=$(PREFIX) -DCMAKE_CXX_FLAGS_RELEASE:STRING=$(CFLAGS) -DWITH_SHARED=OFF -DWITH_HDF5CXX=$(HDF5);\
	make -j16
builduconn:
	mkdir -p build
	cd build;\
	cmake ../channelflow -DWITH_DDC=$(DDC) -DWITH_NSOLVER=$(NSOLVER) -DCMAKE_BUILD_TYPE=$(TYPE) -DCMAKE_INSTALL_PREFIX=$(PREFIX) -DCMAKE_CXX_FLAGS_RELEASE:STRING=$(CFLAGS) -DWITH_SHARED=OFF -DWITH_HDF5CXX=$(HDF5);\
	make -j16
buildpsc:
	mkdir -p build
	cd build;\
	cmake ../channelflow -DCMAKE_CXX_COMPILER=/jet/home/vnguyen9/miniconda3/envs/channelflow/bin/mpicxx -DWITH_NSOLVER=ON -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/chflow -DCMAKE_CXX_FLAGS_RELEASE:STRING=" -fPIC -lfftw3 -lm " -DWITH_SHARED=OFF -DWITH_HDF5CXX=OFF;\
	make -j16
# -DCMAKE_CXX_COMPILER=$(CXX)
run:
	cp ./CMakeLists.txt ./channelflow/CMakeLists.txt
	mkdir -p ./channelflow/modules/
	rm -rf ./channelflow/modules/ddc
	cp -r ./ddc ./channelflow/modules/ddc
	make build
	rm -rf data
# mpiexec -n 16 ./build/modules/ddc/validations/yang2021jfm_case3_2d
# mpiexec -n 16 ./build/modules/ddc/examples/2d_diffusive_convection
# mpiexec -n 16 ./build/modules/ddc/programs/ddc_simulateflow -Pr 10 -Ra 100000 -Le 100 -Rr 2 -dt 0.02 -dT 1 -T 100 -Nx 200 -Ny 81 -Nz 10 -Lx 2 -Lz 0.02 -nl "conv" -Ta 0 -Tb 1 -Sa 0 -Sb 1 
# mpiexec -n 16 ./build/modules/ddc/programs/ddc_findeigenvals -Pr 10 -Ra 100000 -Le 100 -Rr 2 -Nx 200 -Ny 81 -Nz 10 -Lx 2 -Lz 0.02 -Ta 0 -Tb 1 -Sa 0 -Sb 1 
# mpiexec -n 16 ./build/modules/ddc/programs/ddc_continuesoln -eqb

clean:
	rm -rf build
	rm -rf data

allclean:
	rm -rf channelflow
	rm -rf build
	rm -rf data


