SHELL=/bin/csh

MACHINE = IBM_Power4
#MACHINE = IBM_BGL
#MACHINE = INTEL

P3DFFT_ROOT = ..
P3DFFT_LIB = libp3dfft.a

# ----------------

FFT3DLIB = module.o fft_spec.o fft_init.o fft_exec.o wrap.o

all: lib install
lib: $(FFT3DLIB)
	$(AR) $(ARFLAGS) $(P3DFFT_LIB) $(FFT3DLIB)	
install: 
	if(!(-e $(P3DFFT_ROOT))) mkdir $(P3DFFT_ROOT)
	if (!(-e $(P3DFFT_ROOT)/lib)) mkdir $(P3DFFT_ROOT)/lib	
	cp $(P3DFFT_LIB) $(P3DFFT_ROOT)/lib
	if (!(-e $(P3DFFT_ROOT)/include)) mkdir $(P3DFFT_ROOT)/include
	cp p3dfft.mod $(P3DFFT_ROOT)/include

module.o: module.F fft_spec.o
fft_init.o: fft_init.F module.o fft_spec.o
fft_exec.o: fft_exec.F module.o fft_spec.o
wrap.o: wrap.f module.o

.c.o: 
		$(CC) -c $(CFLAGS) *.c
.F.o:
		$(FF) $(DFLAGS) -c $(FFLAGS) $(INCL) $<
.f.o: 
#		$(CPP) $(CPPFLAGS) $< $(SRCDIR)/$<
#		$(FF) -c $(FFLAGS) $(SRCDIR)/$<
		$(FF) -c $(FFLAGS) $(INCL) $<
clean:
	/bin/rm $(FFT3DLIB) *.mod

########################################################################
ifeq ($(MACHINE),IBM_Power4) ##

CC = xlc_r 
FF = mpxlf95_r
FFLAGS = -O3 -qfixed=132 -qarch=pwr4 -qstrict -qtune=pwr4 -q64 -d 
AR = ar
ARFLAGS = -X64 -v -r -u
LDFLAGS= $(FFLAGS) -lessl -L/usr/local/apps/mass -lmass 
# For FFTW use path to the installed FFTW library:
# -L/usr/local/apps/fftw312s/lib -lfftw3f 

DFLAGS = -WF,"-P","-DESSL"
# For FFTW use "-DFFTW" insteadl of "-DESSL"
# For double precision add "-DDOUBLE_PREC"
# For 1D decomposition add "-DONED" above

INCL = -I/usr/local/apps/fftw312s/include 
# For FFTW add include file location, for example: 
# INCL = -I/usr/local/apps/fftw312s/include 

endif

########################################################################
ifeq ($(MACHINE),IBM_BGL) 

CC = xlc
FF = mpxlf90

FFLAGS = -O3 -qfixed=132 -qarch=440 -qstrict -qsave -qtune=440 -qcache=auto -d
AR = ar
ARFLAGS = -v -r -u
LDFLAGS= $(FFLAGS) -L/opt/ibmmath/essl/4.2/lib -lesslbg -L/usr/local/apps/lib -lmass 

DFLAGS = -WF,"-P","-DESSL"
# For FFTW use "-DFFTW" insteadl of "-DESSL"
# For double precision add "-DDOUBLE_PREC"
# For 1D decomposition add "-DONED" above

INCL = -I/usr/local/apps/fftw312s/include 
# For FFTW add include file location, for example: 
# INCL = -I/usr/local/apps/fftw312s/include 

endif

########################################################################
ifeq ($(MACHINE),INTEL)
FF = mpif90 

FFLAGS = -O3 -132 -fpp
LDFLAGS= $(FFLAGS) -L/usr/local/apps/fftw301s/lib -lfftw3f
AR = ar
ARFLAGS = -ruv

#For FFTW

INCL = -I/usr/local/apps/fftw301s/include 
DFLAGS = -Wp,"-P","-DFFTW"
# For double precision add "-DDOUBLE_PREC"
# For 1D decomposition add "-DONED" above

endif



