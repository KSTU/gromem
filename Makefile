all: nvtmain

# Define Fortran compiler
#FC= /opt/intel/composerxe-2011.5.220/bin/ia32/ifort
#FC=/opt/pgi/linux86/10.6/bin/pgfortran  # gfortran
FC=gfortran
CC=gcc

nvtmain: tempprog.f90
	$(FC) -o tempprog tempprog.f90
	$(FC) -o memcv6 create6.f90
	$(FC) -o memcv5 create5.f90
	$(FC) -o dcvmd dcv.f90
	$(FC) -o denscalc denscalc.f90
	$(CC) -O2 -o velcalc velcalc.c -lm	#
	
#cuda_add.o: cuda_add.cu
#	/usr/local/cuda/bin/nvcc -c cuda_add.cu     #-deviceemu

clean: 
	rm tempprog tempprog.o
	rm memcv6 create6.o
	rm memcv5 create5.o
	rm dcvmd dcvmd.o
	rm denscalc denscalc.o
	rm velcalc velcalc.o
