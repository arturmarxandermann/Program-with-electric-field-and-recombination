#

.SUFFIXES: .f .f90

FC = ifort

#FFLAGS = -check all -g -traceback
FFLAGS = -O2 -c -parallel -xHost -ip -align
# O2 paraleliza corretamente

LIB_BLAS = -lmkl_blas95_lp64
LIB_LAPACK = -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
LIB_OMP = -liomp5 -lpthread

INCS_MKL = -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include/fftw

LIB  = $(LIB_BLAS) $(LIB_LAPACK) $(LIB_OMP) -lrt
INCS = $(INCS_MKL)

#-----------------------------------------------------------------------
# general rules
#-----------------------------------------------------------------------
SOURCE = types.o\
       constants.o\
     parameters.o\
     functions.o\
     overlap.o\
     system_hamiltonian.o\
     verlet.o\
     rdftensor.o\
          rkf.o\
     time_evolution.o\
       main.o\




a: $(SOURCE)
	-rm -f a
	$(FC) $(FC_ALL) $(INCS) -o a $(SOURCE) $(LIB)
	-rm -f *.log
.f.o:
	$(FC) -fpp -free $(FC_ALL) $(FFLAGS) $(INCS) -c $<
.f90.o:
	$(FC) -fpp -free $(FC_ALL) $(FFLAGS) $(INCS) -c $<
clean:
	-rm -f *.o *.mod; touch *.f *.f90;
#safe: FC_ALL += -g -traceback -check all -fstack-protector -assume protect_parens -implicitnone -warn all -warn -debug extended -check bounds -check uninit -ftrapuv -debug all -gen-interfaces -warn interfaces -fpe3 -fpe-all=3 -assume ieee_fpe_flags -ftz -mp -fp-model strict -fp-model precise -fp-speculation=off 
#safe: a 	
