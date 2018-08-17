! This sample program illustrates the 
! use of P3DFFT library for highly scalable parallel 3D FFT. 
!
! This program initializes a 3D array with random numbers, then 
! performs forward transform, backward transform, and checks that 
! the results are correct, namely the same as in the start except 
! for a normalization factor. It can be used both as a correctness
! test and for timing the library functions. 
!
! Processor grid in this test is chosen to be square or close to square.
! For better performance, experiment with this setting, varying 
! iproc and jproc. In many cases, minimizing iproc gives best results. 
! Setting it to 1 corresponds to one-dimensional decomposition, which 
! is the best in case it's feasible. 
!
! To build this program, use one of the provided makefiles, modifying 
! it as needed for your compilers, flags and library locations. 
! Be sure to link the program with both the P3DFFT library and the underlying 
! 1D FFT library such as ESSL or FFTW. 
!
! If you have questions please contact Dmitry Pekurovsky, dmitry@sdsc.edu

      program fft3d

      use p3dfft
      implicit none
      include 'mpif.h'

      integer i,n,nx,ny,nz
      integer m,x,y,z
      integer fstatus
      logical flg_inplace

      real(mytype), dimension(:,:,:),  allocatable :: BEG,FIN
      complex(mytype), dimension(:,:,:),  allocatable :: AEND
      real(mytype) pi,twopi,sinyz,diff,cdiff,ccdiff,ans

      integer*8 Ntot
      real(mytype) factor
      real(mytype),dimension(:),allocatable:: sinx,siny,sinz
      real(8) rtime1,rtime2
      real(8) t1,t2,t3,t4,gt1,gt2,gt3,gt4,tp1,gtp1,gtcomm
      common /timers/ t1,t2,t3,t4,tp1
      integer ierr,nu,ndim,dims(2),nproc,proc_id
      integer istart(3),iend(3),isize(3)
      integer fstart(3),fend(3),fsize(3)
      integer iproc,jproc

      call MPI_INIT (ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD,proc_id,ierr)

      pi=atan(1.)*4.
      twopi=2.*pi

      t1 = 0.0
      t2 = 0.0
      t3 = 0.0
      t4 = 0.0
      tp1 = 0.0
      gt1=0.0
      gt2=0.0
      gt3=0.0
      gt4=0.0
      gtp1 = 0.0

      if (proc_id.eq.0) then 
         open (unit=3,file='stdin',status='old',
     &         access='sequential',form='formatted', iostat=fstatus)
         if (fstatus .eq. 0) then
            write(*, *) ' Reading from input file stdin'
         endif 
         ndim = 2

        read (3,*) nx, ny, nz, ndim,n
        write (*,*) "procs=",nproc," nx=",nx,
     &          " ny=", ny," nz=", nz,"ndim=",ndim," repeat=", n
        if(mytype .eq. 4) then
           print *,'Single precision version'
        else if(mytype .eq. 8) then
           print *,'Double precision version'
        endif
       endif

      call MPI_Bcast(nx,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(ny,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(nz,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(n,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(ndim,1, MPI_INTEGER,0,mpi_comm_world,ierr)

!    nproc is devided into a iproc x jproc stencle
!

      if(ndim .eq. 1) then
         dims(1) = 1
         dims(2) = nproc
      else if(ndim .eq. 2) then
         dims(1) = 0
         dims(2) = 0
         call MPI_Dims_create(nproc,2,dims,ierr)
         if(dims(1) .gt. dims(2)) then
            dims(1) = dims(2)
            dims(2) = nproc / dims(1)
         endif
      endif

      iproc = dims(1)
      jproc = dims(2)

      if(proc_id .eq. 0) then
         print *,'Using processor grid ',iproc,' x ',jproc
      endif

      call p3dfft_setup (dims,nx,ny,nz)
      call get_dims(istart,iend,isize,1)
      call get_dims(fstart,fend,fsize,2)

      allocate (sinx(nx))
      allocate (siny(ny))
      allocate (sinz(nz))
c
c initialize
c

c      print *,'Allocating BEG (',isize,istart,iend
      allocate (BEG(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)), stat=ierr)
      if(ierr .ne. 0) then
         print *,'Error ',ierr,' allocating array BEG'
      endif
c      print *,'Allocating AEND (',fsize,fstart,fend
      allocate (AEND(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)), stat=ierr)
      if(ierr .ne. 0) then
         print *,'Error ',ierr,' allocating array AEND'
      endif
      allocate (FIN(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)), stat=ierr)
      if(ierr .ne. 0) then
         print *,'Error ',ierr,' allocating array FIN'
      endif

c start with x-z slabs in physical space
c
      do z=istart(3),iend(3)
         do y=istart(2),iend(2)
            do x=istart(1),iend(1)
               call random_number(BEG(x,y,z))
            enddo
         enddo
      enddo

c
c transform from physical space to wavenumber space
c (XgYiZj to XiYjZg)
c then transform back to physical space
c (XiYjZg to XgYiZj)
c

      factor = 1.0d0/(nx*ny)
      factor = factor / nz
      Ntot = fsize(1)*fsize(2)*fsize(3)

      rtime1 = 0.0               
      do  m=1,n
         if(proc_id .eq. 0) then
            print *,'Iteration ',m
         endif

         rtime1 = rtime1 - MPI_wtime()
         call ftran_r2c (BEG,AEND)
         
         rtime1 = rtime1 + MPI_wtime()
         
c         print *,'Result of forward transform:'
c         call print_all(AEND,Ntot,proc_id)
         
         call mult_array(AEND, Ntot,factor)
         
         rtime1 = rtime1 - MPI_wtime()
         call p3dfft_btran_c2r (AEND,FIN)       
         rtime1 = rtime1 + MPI_wtime()
         
      end do

      call p3dfft_clean

      call MPI_Reduce(rtime1,rtime2,1,mpi_real8,MPI_MAX,0,
     &  MPI_COMM_WORLD,ierr)

      if (proc_id.eq.0) write(6,*)'proc_id, cpu time per loop',
     &   proc_id,rtime2/dble(n )

      call MPI_Reduce(t1,gt1,1,mpi_real8,MPI_MAX,0,
     &  MPI_COMM_WORLD,ierr)
      call MPI_Reduce(tp1,gtp1,1,mpi_real8,MPI_MAX,0,
     &  MPI_COMM_WORLD,ierr)
      call MPI_Reduce(t2,gt2,1,mpi_real8,MPI_MAX,0,
     &  MPI_COMM_WORLD,ierr)
      call MPI_Reduce(t3,gt3,1,mpi_real8,MPI_MAX,0,
     &  MPI_COMM_WORLD,ierr)
      call MPI_Reduce(t4,gt4,1,mpi_real8,MPI_MAX,0,
     &  MPI_COMM_WORLD,ierr)
      call MPI_Reduce(t1+t2+t3+t4,gtcomm,1,mpi_real8,MPI_MAX,0,
     &  MPI_COMM_WORLD,ierr)

      gt1 = gt1 / dble(n)
      gtp1 = gtp1 / dble(n)
      gt2 = gt2 / dble(n)
      gt3 = gt3 / dble(n)
      gt4 = gt4 / dble(n)
      gtcomm = gtcomm / dble(n)

      if(proc_id .eq. 0) then
         print *,'t1,tp1,t2,t3,t4: ',gt1,gtp1,gt2,gt3,gt4
         print *,'Total comm: ',gtcomm
      endif

      cdiff=0.0d0
      do 20 z=istart(3),iend(3)
         do 20 y=istart(2),iend(2)
            sinyz=siny(y)*sinz(z)
            do 20 x=istart(1),iend(1)
            ans=sinx(x)*sinyz
            if(cdiff .lt. abs(BEG(x,y,z)-FIN(x,y,z))) then
               cdiff = abs(BEG(x,y,z)-FIN(x,y,z))
c               print *,'x,y,z,cdiff=',x,y,z,cdiff
            endif
 20   continue
      call MPI_Reduce(cdiff,ccdiff,1,mpireal,MPI_MAX,0,
     &  MPI_COMM_WORLD,ierr)

      if (proc_id.eq.0) write (6,*) 'max diff at do 20=',ccdiff

      call MPI_FINALIZE (ierr)

      contains 
!=========================================================

      subroutine mult_array(X,nar,f)

      integer*8 nar,i
      complex(mytype) X(nar)
      real(mytype) f

      do i=1,nar
         X(i) = X(i) * f
      enddo

      return
      end subroutine
#line 291
      subroutine print_all(Ar,Nar,proc_id)

      use p3dfft

      integer x,y,z,proc_id
      integer(8) i,Nar
      complex(mytype) Ar(1,1,*)
      integer Fstart(3),Fend(3),Fsize(3)

      call get_dims(Fstart,Fend,Fsize,2)
      do i=1,Nar
         if(abs(Ar(1,1,i)) .gt. Nar *1.25e-4) then
            z = (i-1)/(Fsize(1)*Fsize(2))
            y = (i-1 - z * Fsize(1)*Fsize(2))/Fsize(1)
            x = i-1-z*Fsize(1)*Fsize(2) - y*Fsize(1)
            print *,'(',x+Fstart(1),y+Fstart(2),z+Fstart(3),') ',Ar(1,1,i)
         endif
      enddo

      return
      end subroutine

      end
