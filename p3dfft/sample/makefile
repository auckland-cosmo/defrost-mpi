SHELL=/bin/csh

MACHINE = IBM_Power4
#MACHINE = IBM_BGL
#MACHINE = INTEL

P3DFFT_ROOT = ..
P3DFFT_LIB = p3dfft

########################################################################
ifeq ($(MACHINE),IBM_Power4)

CC = xlc_r 
FF = mpxlf95_r
FFLAGS = -O3 -qfixed=132 -qarch=pwr4 -qstrict -qtune=pwr4 -q64 -d 
LDFLAGS= $(FFLAGS) -lessl -L/usr/local/apps/mass -lmass 
# For FFTW use path to the installed FFTW library:
# -L/usr/local/apps/fftw301s/lib -lfftw3f 

endif

########################################################################
ifeq ($(MACHINE),IBM_BGL)

CC = xlc
FF = mpxlf90

FFLAGS = -O3 -qfixed=132 -qarch=440 -qstrict -qsave -qtune=440 -qcache=auto -d
AR = ar
ARFLAGS = -v -r -u
LDFLAGS= $(FFLAGS) -L/opt/ibmmath/essl/4.2/lib -lesslbg -L/usr/local/apps/lib -lmass 
# For FFTW use path to the installed FFTW library:
# -L/usr/local/apps/fftw312s/lib -lfftw3f 

endif

########################################################################
ifeq ($(MACHINE),INTEL)

FF = mpif90 

FFLAGS = -O3 -132 -fpp
LDFLAGS= $(FFLAGS) -L/usr/local/apps/fftw312s/lib -lfftw3f
AR = ar
ARFLAGS = -ruv

endif

##-----------------

all: driver1.o driver2.o 
	$(FF) -o test1.x driver1.o -L$(P3DFFT_ROOT)/lib -l$(P3DFFT_LIB) $(LDFLAGS)
	$(FF) -o test2.x driver2.o -L$(P3DFFT_ROOT)/lib -l$(P3DFFT_LIB) $(LDFLAGS)

driver1.o: driver1.f
driver2.o: driver2.f

.f.o: 
		$(FF) -c $(FFLAGS) -I$(P3DFFT_ROOT)/include $<

clean:
	/bin/rm driver*.o

