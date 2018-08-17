################################################################
# $Id: Makefile,v 1.2 2007/09/19 23:14:02 frolov Exp frolov $
# DEFROST Makefile (Intel Fortran compiler)
################################################################

.SUFFIXES: .f .F .f90

ifndef MPIFC
MPIFC = mpif90
endif

PG_INFO=$(shell $(MPIFC) -V 2>&1 | grep pgf90)
ifneq ($(PG_INFO),)
USE_PG=1
endif

GCC_INFO=$(shell $(MPIFC) -v 2>&1 | grep "gcc version")
ifneq ($(GCC_INFO),)
USE_GNU=1
endif

FFLAGS_DEFS =
BASEINC = -I$(shell cd $(dir $(shell which fftw-wisdom))../include && pwd)

BASELIBDIR = $(shell cd $(dir $(shell which fftw-wisdom))../lib && pwd)
BASELIB = -L$(BASELIBDIR)
LDFLAGS = $(BASELIB)
export LD_RUN_PATH = $(BASELIBDIR)

ifdef USE_PG
FFLAGS = -O3 -fast -r8 # -Mconcur -mp
FFLAGS_COMMON = -Mpreprocess
FFLAGS_FREE = -Mfreeform -Mextend
FFLAGS_FIXED = -Mfixed -Mextend
LDFLAGS += 
else
ifdef USE_GNU
FFLAGS = -O3 -march=pentium4 -fdefault-real-8 -fdollar-ok # -fopenmp
FFLAGS_COMMON = -x f95-cpp-input 
FFLAGS_FREE = -ffree-line-length-none
FFLAGS_FIXED = -ffixed-form -ffixed-line-length-none
LDFLAGS += 
else
# Fortran compiler (adjust for your machine, -r8 and -fpp are mandatory)
FFLAGS = -O3 -ipo -xP -r8 -pc80 # -parallel -openmp
FFLAGS_COMMON = -fpp
FFLAGS_FREE = 
FFLAGS_FIXED = -132
LDFLAGS += -static-intel
endif
endif

FFTWINC = 
FFTWLIB = -lfftw3

NO_FFTW_THREADS=1
ifneq ($(shell ls $(BASELIBDIR) | grep libfftw3_threads),)
ifndef NO_FFTW_THREADS
FFTWLIB += -lfftw3_threads
FFLAGS_DEFS += -DUSE_FFTW_THREADS
endif
endif

FFT3DLIB = p3dfft/build/module.o p3dfft/build/fft_spec.o p3dfft/build/fft_init.o p3dfft/build/fft_exec.o p3dfft/build/wrap.o
FFLAGS_DEFS += -DFFTW -DDOUBLE_PREC

# SILO library (native VisIt data format)
ifdef USE_SILO
SILOINC = 
SILOLIB = -L/usr/local/lib64 -lsilo -lhdf5

FFLAGS_DEFS += -DSILO
endif


################################################################

L_VALUES = 400 # 100 200 400
N_VALUES = 1024 # 64 128 256 512 1024
LL_VALUES = 2.8125E-6 # 1.25E-6 2.8125E-6 5.0E-6 7.8125E-6 1.125E-5
GG_FACTS = 10.0 # 5.0 10.0
SEED_VALUES = 2 # 1 2 3 4 5
EXES = $(foreach L,$(L_VALUES),$(foreach N,$(N_VALUES),$(foreach LL,$(LL_VALUES),$(foreach GG,$(GG_FACTS),$(foreach SEED,$(SEED_VALUES),defrost_$L_$N_$(LL)_$(GG)_$(SEED))))))
OBJS = $(foreach E,$(EXES),$E.o)

all: $(EXES)

clean:
	rm -f *.bak gmon.out core *.o *.mod $(EXES)
	rm -f p3dfft/build/*.o

$(EXES): %: %.o $(FFT3DLIB)
	$(MPIFC) $(FFLAGS) $^ -o $@ $(LDFLAGS) $(FFTWLIB) $(SILOLIB)

################################################################

p3dfft/build/module.o: p3dfft/build/module.F p3dfft/build/fft_spec.o
p3dfft/build/fft_init.o: p3dfft/build/fft_init.F p3dfft/build/module.o p3dfft/build/fft_spec.o
p3dfft/build/fft_exec.o: p3dfft/build/fft_exec.F p3dfft/build/module.o p3dfft/build/fft_spec.o
p3dfft/build/wrap.o: p3dfft/build/wrap.f p3dfft/build/module.o
$(OBJS): p3dfft/build/module.o

$(OBJS): defrost.f90 parameters.inc model.inc
	$(MPIFC) $(BASEINC) -c $(FFLAGS_COMMON) $(FFLAGS_FREE) $(FFLAGS) $(FFLAGS_DEFS) -DL_VALUE=$(word 2,$(subst _, ,$@)).0 -DN_VALUE=$(word 3,$(subst _, ,$@)) -DLL_VALUE=$(word 4,$(subst _, ,$@)) -DGG_FACT=$(word 5,$(subst _, ,$@)) -DSEED=$(word 6,$(subst _, ,$(subst ., ,$@))) $(FFTWINC) $(SILOINC) $< -o $@

.F.o:
	$(MPIFC) $(BASEINC) -c $(FFLAGS_COMMON) $(FFLAGS_FIXED) $(FFLAGS) $(FFLAGS_DEFS) $(FFTWINC) $< -o $@

.f.o:
	$(MPIFC) $(BASEINC) -c $(FFLAGS_COMMON) $(FFLAGS_FIXED) $(FFLAGS) $(FFLAGS_DEFS) $(FFTWINC) $< -o $@

