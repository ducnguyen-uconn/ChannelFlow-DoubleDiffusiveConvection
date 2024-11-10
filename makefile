OUTPUT = 2d_diffusive_convection 2d_finger_convection 2d_finger_convection 3d_finger_convection 2d_finger data yang2021jfm_case3_2d binary_fluid_convection *flags.txt *args processinfo *.asc *.nc


.PHONY: clean update build
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
	cp ./CMakeLists.txt ./channelflow/CMakeLists.txt
	mkdir -p ./channelflow/modules/
	rm -rf ./channelflow/modules/ddc
	cp -r ./ddc ./channelflow/modules/ddc
	rm -rf build
	mkdir -p build
	cd build;\
	cmake ../channelflow -DCMAKE_CXX_COMPILER=/usr/bin/mpicxx -DWITH_ILC=ON -DWITH_DDC=ON -DWITH_NSOLVER=ON -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/user/local/ -DCMAKE_CXX_FLAGS_RELEASE:STRING=" -fPIC -lfftw3 -lm -Wno-unused-variable " -DWITH_SHARED=OFF -DWITH_HDF5CXX=OFF;\
	make -j16
builduconn:
	git clone https://github.com/epfl-ecps/channelflow.git
	cp ./CMakeLists.txt ./channelflow/CMakeLists.txt
	mkdir -p ./channelflow/modules/
	rm -rf ./channelflow/modules/ddc
	cp -r ./ddc ./channelflow/modules/ddc
	mkdir -p build
	cd build;\
	cmake ../channelflow -DCMAKE_CXX_COMPILER=/gpfs/sharedfs1/admin/hpc2.0/apps/openmpi/5.0.5/bin/mpicxx -DWITH_DDC=ON -DWITH_NSOLVER=ON -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/chflow -DCMAKE_CXX_FLAGS_RELEASE:STRING=" -fPIC -lfftw3 -lm -Wno-unused-variable -Wno-deprecated " -DWITH_SHARED=OFF -DWITH_HDF5CXX=OFF;\
	make -j16
buildpsc:
	git clone https://github.com/epfl-ecps/channelflow.git
	cp ./CMakeLists.txt ./channelflow/CMakeLists.txt
	mkdir -p ./channelflow/modules/
	rm -rf ./channelflow/modules/ddc
	cp -r ./ddc ./channelflow/modules/ddc
	mkdir -p build
	cd build;\
	cmake ../channelflow -DCMAKE_CXX_COMPILER=/jet/home/vnguyen9/miniconda3/envs/channelflow/bin/mpicxx -DWITH_NSOLVER=ON -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/chflow -DCMAKE_CXX_FLAGS_RELEASE:STRING=" -fPIC -lfftw3 -lm " -DWITH_SHARED=OFF -DWITH_HDF5CXX=OFF;\
	make -j16

run:
	mpiexec -n 16 ./build/modules/ddc/programs/ddc_findeigenvals -o "bfc/" -Pr 7 -Ra 1908 -Le 100 -Rs -0.1 -Nx 652 -Ny 33 -Nz 10 -Lx 14 -Lz 0.02 -Ta 1 -Tb 0 -Sa 0 -Sb 0 -dt 0.009 -dtmax 0.01 -T 100
runbfc: # binary fluid convection
	rm -rf binary_fluid_convection
	mpiexec -n 16 ./build/modules/ddc/programs/ddc_simulateflow -o "binary_fluid_convection/" -Pr 7 -Ra 1908 -Le 100 -Rs -0.1 -Nx 652 -Ny 33 -Nz 6 -Lx 14 -Lz 0.004 -Ta 1 -Tb 0 -Sa 0 -Sb 0 -dt 0.005 -dtmax 0.01 -dT 0.5 -T 100
example:
	rm -rf 2d_diffusive_convection
	mpiexec -n 6 ./build/modules/ddc/examples/2d_diffusive_convection
2d_finger_convection:
	rm -rf 2d_finger
	mpiexec -n 16 ./build/modules/ddc/programs/ddc_simulateflow -o "2d_finger/" -Pr 7 -Ra 10000 -Le 100 -Rr 2 -Nx 768 -Ny 385 -Nz 6 -Lx 2 -Lz 0.004 -Ta 0 -Tb 1 -Sa 0 -Sb 1 -dt 0.01 -dtmin 1e-9 -dtmax 0.02 -dT 0.02 -T 5
yang2021jfm_case3:
	rm -rf yang2021jfm_case3_2d
	mpiexec -n 16 ./build/modules/ddc/programs/ddc_simulateflow -m "yang2021jfm" -o "yang2021jfm_case3_2d/" -Pr 10 -Ra 1000000 -Le 100 -Rr 2 -Nx 384 -Ny 385 -Nz 6 -Lx 2 -Lz 0.004 -Ua -0.5 -Ub 0.5 -Ta 1 -Tb 0 -Sa 1 -Sb 0 -dt 1e-3 -dtmin 1e-9 -dtmax 1e-2 -dT 1e-0 -T 100
# mpiexec -n 16 ./build/modules/ddc/validations/yang2021jfm_case3_2d
eigen:
	mpiexec -n 16 ./build/programs/findeigenvals 
testilc:
	rm -rf data
	mpiexec -n 16 ./build/modules/ilc/programs/ilc_simulateflow -Pr 10 -Ra 1000000 -ilcUw 0.5 -dt 0.001 -dT 0.01 u t
	
# 
# mpiexec -n 16 ./build/modules/ddc/programs/ddc_simulateflow -Pr 10 -Ra 100000 -Le 100 -Rr 2 -dt 0.02 -dT 1 -T 100 -Nx 200 -Ny 81 -Nz 10 -Lx 2 -Lz 0.02 -nl "conv" -Ta 0 -Tb 1 -Sa 0 -Sb 1 
# mpiexec -n 16 ./build/modules/ddc/programs/ddc_findeigenvals -Pr 10 -Ra 100000 -Le 100 -Rr 2 -Nx 200 -Ny 81 -Nz 10 -Lx 2 -Lz 0.02 -Ta 0 -Tb 1 -Sa 0 -Sb 1 
# mpiexec -n 16 ./build/modules/ddc/programs/ddc_continuesoln -eqb

clean:
	rm -rf build data $(OUTPUT)

allclean:
	rm -rf channelflow build $(OUTPUT)


