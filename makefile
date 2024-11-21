OUTPUT = 2d_diffusive_convection 2d_finger_convection 2d_finger_convection \
		3d_finger_convection 2d_finger data yang2021jfm_case3_2d binary_fluid_convection \
		*flags.txt *args processinfo *.asc *.nc snapshots


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
	cmake ../channelflow -DCMAKE_CXX_COMPILER=/opt/packages/openmpi/gnu/5.0.3-gcc13.2.1-cpu/bin/mpicxx -DWITH_DDC=ON -DWITH_NSOLVER=ON -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX=/chflow -DCMAKE_CXX_FLAGS_RELEASE:STRING=" -fPIC -lfftw3 -lm -Wno-unused-variable -Wno-deprecated " -DWITH_SHARED=OFF -DWITH_HDF5CXX=OFF;\
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


yang2021jfm_case3_simulateflow:
	rm -rf yang2021jfm_case3_2d
	mpiexec -n 16 ./build/modules/ddc/programs/ddc_simulateflow -m "yang2021jfm" -o "yang2021jfm_case3_2d/" -Pr 10 -Ra 1000000 -Le 100 -Ri 1 -Rr 0.5 -Nx 192 -Ny 193 -Nz 6 -Lx 2 -ymin 0 -ymax 1 -Lz 0.004 -Ua -0.5 -Ub 0.5 -Ta 1 -Tb 0 -Sa 1 -Sb 0 -dt 5e-3 -dtmin 1e-9 -dtmax 5e-2 -dT 1e0 -T 100 -ys 1
yang2021jfm_case5_simulateflow:
	rm -rf yang2021jfm_case5_2d
	mpiexec -n 16 ./build/modules/ddc/programs/ddc_simulateflow -m "yang2021jfm_case5" -o "yang2021jfm_case5_2d/" -Pr 10 -Ra 1e7 -Le 100 -Ri 1 -Rr 0.5 -Nx 768 -Ny 385 -Nz 6 -Lx 4 -Lz 0.004 -Ua -0.5 -Ub 0.5 -Ta 1 -Tb 0 -Sa 1 -Sb 0 -dt 2.5e-3 -dtmin 1e-9 -dtmax 1e-2 -dT 10 -T 6000
yang2021jfm_case3_findsoln:
	mpiexec -n 16 ./build/modules/ddc/programs/ddc_findsoln -orb -Pr 10 -Ra 1000000 -Le 100 -Rr 2 -Ua -0.5 -Ub 0.5 -Ta 1 -Tb 0 -Sa 1 -Sb 0 -dt 1e-3 inputfields/u100 inputfields/t100 inputfields/s100

eigen:
	mpiexec -n 16 ./build/programs/findeigenvals 
testilc:
	rm -rf data
	mpiexec -n 16 ./build/modules/ilc/programs/ilc_simulateflow -Pr 10 -Ra 1000000 -ilcUw 0.5 -dt 0.001 -dT 0.01 u t
	
spcf: # stratified plane Couette flow
	rm -rf spcf
	mpiexec -n 16 ./build/modules/ddc/programs/ddc_simulateflow -o "spcf/" -Rey 1e3 -Pr 1 -Ri -1.22e-3 -ymin -1 -ymax 1 -Nx 64 -Ny 129 -Nz 64 -Lx 6.2831853 -Lz 6.2831853 -Ua -1 -Ub 1 -Ta 1 -Tb -1 -dt 1e-2 -dtmax 1e-1 -dT 1e1 -T 100 -ys 1
# mpiexec -n 16 ./build/modules/ddc/programs/ddc_findsoln -eqb -Rey 1e3 -Pr 1 -Ri -1.22e-3 -Ua -1 -Ub 1 -Ta 1 -Tb -1 -dt 1e-2 -ys 1 inputfields/u10 inputfields/t10 inputfields/t10
# mpiexec -n 16 ./build/modules/ddc/programs/ddc_edgetracking -Rey 1e3 -Pr 1 -Ri -1.22e-3 -Ua -1 -Ub 1 -Ta 1 -Tb -1 -dt 1e-2 -dtmax 1e-1 -ys 1 inputfields/u10 inputfields/t10 inputfields/t10

langham2020jfm_simulate:
	rm -rf langham2020jfm
	mpiexec -n 16 ./build/modules/ddc/programs/ddc_simulateflow -o "langham2020jfm/" -Rey 400 -Pr 1 -Ri 0.01 -ymin -1 -ymax 1 -Nx 64 -Ny 65 -Nz 32 -Lx 6.2831853 -Lz 3.14159 -Ua -1 -Ub 1 -Ta 1 -Tb -1 -dt 2.5e-2 -dtmax 1e-1 -dT 1e1 -T 1000 -ys 1
langham2020jfm_findsoln: 
# mpiexec -n 16 ./build/modules/ddc/programs/ddc_findsoln -orb -Rey 400 -Pr 1 -Ri 0.01 -Ua -1 -Ub 1 -Ta 1 -Tb -1 -dt 2.5e-2 -ys 1 langham2020jfm/u40 langham2020jfm/t40 langham2020jfm/t40
	mpiexec -n 16 ./build/modules/ddc/programs/ddc_findsoln -eqb -Rey 400 -Pr 1 -Ri 0.01 -Ua -1 -Ub 1 -Ta 1 -Tb -1 -dt 2.5e-2 -ys 1 langham2020jfm/u40 langham2020jfm/t40 langham2020jfm/t40

eaves2016jfm_simulate_p600:
	rm -rf eaves2016jfm_p600
	mpiexec -n 16 ./build/modules/ddc/programs/ddc_simulateflow -o "eaves2016jfm_p600/" -Rey 600 -Pr 300 -Ri 0.52 -ymin -1 -ymax 1 -Nx 256 -Ny 257 -Nz 6 -Lx 6.2831853 -Lz 0.004 -Ua -1 -Ub 1 -Ta 1 -Tb -1 -dt 2e-3 -dtmax 1e-1 -dT 1e1 -T 400 -ys 1

cf_ddc: # Couette flow [2π/1.14, 2, 4π/5]
	rm -rf couetteflow_ddc
	mpiexec -n 16 ./build/modules/ddc/programs/ddc_simulateflow -o "couetteflow_ddc/" -Rey 1000 -ymin -1 -ymax 1 -Nx 64 -Ny 65 -Nz 64 -Lx 5.5115660589 -Lz 2.51327412287 -Ua -1 -Ub 1 -dt 2e-2 -dtmax 5e-2 -dT 1e0 -T 10 -nl "conv" -ys 1

cf_chflow: # Couette flow [2π/1.14, 2, 4π/5]
	rm -rf couetteflow_chflow
	mpiexec -n 16 ./build/programs/simulateflow -o "couetteflow_chflow/" -R 400 -Uwall 1 -theta 0 -dt 2e-2 -dtmax 5e-2 -dT 1e1 -T 500 -nl "conv" u
cf_soln: # Couette flow [2π/1.14, 2, 4π/5]
	mpiexec -n 16 ./build/programs/findsoln -eqb -R 400 -Uwall 1 -theta 0 -dt 2e-2 couetteflow_chflow/u20


# mpiexec -n 16 ./build/modules/ddc/programs/ddc_simulateflow -m "yang2021jfm" -o "testing/" -Rey 1e3 -Pr 1 -Ra 1e5 -Le 100 -Rr 2 -Rs -0.1 -Ri 1e-3 -ymin -1 -ymax 1 -Nx 128 -Ny 65 -Nz 128 -Lx 6.2831853 -Lz 6.2831853 -Ua -1 -Ub 1 -Ta 1 -Tb -1 -Sa 0 -Sb 0 -dt 1e-2 -dtmin 1e-9 -dtmax 1e-1 -dT 1e1 -T 100 -ys 1
# mpiexec -n 16 ./build/modules/ddc/programs/ddc_simulateflow -m "yang2021jfm" -o "testing/" -Pr 10 -Ra 10000 -Le 100 -Rr 2 -Nx 128 -Ny 129 -Nz 6 -Lx 2 -Lz 0.004 -Ua -0.5 -Ub 0.5 -Ta 1 -Tb 0 -Sa 1 -Sb 0 -dt 1e-2 -dtmin 1e-9 -dtmax 1e-1 -dT 1e-0 -T 100 -ys 1
# mpiexec -n 16 ./build/modules/ddc/programs/ddc_findsoln -eqb -Pr 10 -Ra 1000000 -Le 100 -Rr 2 -Ua -0.5 -Ub 0.5 -Ta 1 -Tb 0 -Sa 1 -Sb 0 -dt 1e-3 inputfields/u6 inputfields/t6 inputfields/s6
# 
# mpiexec -n 16 ./build/modules/ddc/programs/ddc_simulateflow -Pr 10 -Ra 100000 -Le 100 -Rr 2 -dt 0.02 -dT 1 -T 100 -Nx 200 -Ny 81 -Nz 10 -Lx 2 -Lz 0.02 -nl "conv" -Ta 0 -Tb 1 -Sa 0 -Sb 1 
# mpiexec -n 16 ./build/modules/ddc/programs/ddc_findeigenvals -Pr 10 -Ra 100000 -Le 100 -Rr 2 -Nx 200 -Ny 81 -Nz 10 -Lx 2 -Lz 0.02 -Ta 0 -Tb 1 -Sa 0 -Sb 1 
# mpiexec -n 16 ./build/modules/ddc/programs/ddc_continuesoln -eqb

clean:
	rm -rf build data $(OUTPUT)

allclean:
	rm -rf channelflow build $(OUTPUT)


