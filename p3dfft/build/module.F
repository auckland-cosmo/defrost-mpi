! This file is part of P3DFFT library
!
! Version 2.2 beta
!
! Copyright (C) 2006-2008 Dmitry Pekurovsky
!
!    P3DFFT is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation; either version 2 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program; if not, write to the Free Software
!    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

!----------------------------------------------------------------------------

! global variables
      module p3dfft

      include 'mpif.h'

      private

#ifdef DOUBLE_PREC
       integer, parameter,public :: mytype=8
       integer, parameter,public:: mpireal = MPI_DOUBLE_PRECISION
       integer,parameter,public:: mpicomplex = MPI_DOUBLE_COMPLEX
#else
       integer, parameter,public :: mytype=4
       integer, parameter,public:: mpireal = MPI_REAL
       integer,parameter,public:: mpicomplex = MPI_COMPLEX
#endif

      integer, public :: padd

      integer,save :: NX_fft,NY_fft,NZ_fft,numtasks,iproc,jproc
      integer,save :: ipid,jpid,taskid
      integer,save :: iistart,iiend,iisize,jjstart,jjsize,jjend
      integer,save ::jistart,kjstart,jisize,kjsize,jiend,kjend

      integer,save ::  nxh,nxhp

! mpi process info
!
      logical :: mpi_set=.false.
      integer, save :: mpi_comm_cart      
      integer, save :: mpi_comm_row, mpi_comm_col
      integer,save, dimension(:), allocatable :: iist,iien,iisz
      integer,save, dimension(:), allocatable :: jist,jien,jisz
      integer,save, dimension(:), allocatable :: jjst,jjen,jjsz
      integer,save, dimension(:), allocatable :: kjst,kjen,kjsz

! mpi derived data types for implementing alltoallv using send-recvs
      integer,save,dimension(:),allocatable:: IfSndCnts,IfSndStrt
      integer,save,dimension(:),allocatable:: IfRcvCnts,IfRcvStrt
      integer,save,dimension(:),allocatable:: KfSndCnts,KfSndStrt
      integer,save,dimension(:),allocatable:: KfRcvCnts,KfRcvStrt
      integer,save,dimension(:),allocatable:: JrSndCnts,JrSndStrt
      integer,save,dimension(:),allocatable:: JrRcvCnts,JrRcvStrt
      integer,save,dimension(:),allocatable:: KrSndCnts,KrSndStrt
      integer,save,dimension(:),allocatable:: KrRcvCnts,KrRcvStrt
      integer,save,dimension(:,:),allocatable:: status

      public :: get_dims,p3dfft_setup,ftran_r2c,btran_c2r,p3dfft_clean,print_buf

      contains

      subroutine get_dims(istart,iend,isize,conf)

      integer istart(3),iend(3),isize(3),conf
      
      if(.not. mpi_set) then
         print *,'P3DFFT error: call setup before other routines'
      else

      if(conf .eq. 1) then
         istart(1) = 1
         iend(1) = NX_fft
         isize(1) = NX_fft
         istart(2) = jistart
         iend(2) = jiend
         isize(2) = jisize
         istart(3) = kjstart
         iend(3) = kjend
         isize(3) = kjsize
      else if(conf .eq. 2) then
         istart(1) = iistart
         iend(1) = iiend
         isize(1) = iisize
         istart(2) = jjstart
         iend(2) = jjend
         isize(2) = jjsize
         istart(3) = 1
         iend(3) = NZ_fft
         isize(3) = NZ_fft
      endif

      endif
      end subroutine get_dims

c =========================================================
      subroutine p3dfft_setup(dims,nx,ny,nz)
c========================================================

      implicit none

      integer i,j,k,nx,ny,nz
      integer ierr, dims(2),  cartid(2)
      logical periodic(2),remain_dims(2)
      integer impid, ippid, jmpid, jppid,nm
      complex(mytype), allocatable :: A(:)
      real(mytype), allocatable :: B(:)

      mpi_set = .true.
      nx_fft = nx
      ny_fft = ny
      nz_fft = nz
      nxh=nx/2
      nxhp=nxh+1

      call MPI_COMM_SIZE (MPI_COMM_WORLD,numtasks,ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD,taskid,ierr)

      iproc = dims(1)
      jproc = dims(2)
      periodic(1) = .false.
      periodic(2) = .false.
! creating cartesian processor grid
      call MPI_Cart_create(MPI_COMM_WORLD,2,dims,periodic,
     &     .false.,mpi_comm_cart,ierr)
! Obtaining process ids with in the cartesian grid
      call MPI_Cart_coords(mpi_comm_cart,taskid,2,cartid,ierr)
! process with a linear id of 5 may have cartid of (3,1)
      ipid = cartid(1)
      jpid = cartid(2)
! here i is east-west j is north-south
! impid is west neighbour ippid is east neighbour and so on
      impid = ipid - 1
      ippid = ipid + 1
      jmpid = jpid - 1
      jppid = jpid + 1
!boundary processes
      if (ipid.eq.0) impid = MPI_PROC_NULL
      if (jpid.eq.0) jmpid = MPI_PROC_NULL
      if (ipid.eq.iproc-1) ippid = MPI_PROC_NULL
      if (jpid.eq.jproc-1) jppid = MPI_PROC_NULL
! using cart comworld create east-west(row) sub comworld
      remain_dims(1) = .true.
      remain_dims(2) = .false.
      call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_row,ierr)
! using cart comworld create north-south(column) sub comworld
      remain_dims(1) = .false.
      remain_dims(2) = .true.
      call MPI_Cart_sub(mpi_comm_cart,remain_dims,mpi_comm_col,ierr)

      allocate (iist(0:iproc-1))
      allocate (iisz(0:iproc-1))
      allocate (iien(0:iproc-1))
      allocate (jjst(0:jproc-1))
      allocate (jjsz(0:jproc-1))
      allocate (jjen(0:jproc-1))
      allocate (jist(0:iproc-1))
      allocate (jisz(0:iproc-1))
      allocate (jien(0:iproc-1))
      allocate (kjst(0:jproc-1))
      allocate (kjsz(0:jproc-1))
      allocate (kjen(0:jproc-1))
!
!Mapping 3-D data arrays onto 2-D process grid
! (nx+2,ny+2,nz) => (iproc,jproc)      
! 
      call MapDataToProc(nxhp,iproc,iist,iien,iisz)
c      call MapDataToProc(nxhp,iproc,xist,xien,xisz)
      call MapDataToProc(ny,iproc,jist,jien,jisz)
      call MapDataToProc(ny,jproc,jjst,jjen,jjsz)
      call MapDataToProc(nz,jproc,kjst,kjen,kjsz)


      iistart = iist(ipid)
      jjstart = jjst(jpid)
      jistart = jist(ipid)
      kjstart = kjst(jpid)
      iisize= iisz(ipid)
      jjsize= jjsz(jpid)
      jisize= jisz(ipid)
      kjsize= kjsz(jpid)
      iiend = iien(ipid)
      jjend = jjen(jpid)
      jiend = jien(ipid)
      kjend = kjen(jpid)

      padd = max(iisize*jjsize*nz_fft,iisize*ny_fft*kjsize) - nxhp*jisize*kjsize
      if(padd .le. 0) then 
         padd=0
      else
         if(mod(padd,nxhp*jisize) .eq. 0) then
            padd = padd / (nxhp*jisize)
         else
            padd = padd / (nxhp*jisize)+1
         endif

      endif

! Displacements and buffer counts for mpi_alltoallv

      allocate (IfSndStrt(0:iproc-1))
      allocate (IfSndCnts(0:iproc-1))     
      allocate (IfRcvStrt(0:iproc-1))
      allocate (IfRcvCnts(0:iproc-1))

      allocate (KfSndStrt(0:jproc-1))
      allocate (KfSndCnts(0:jproc-1))     
      allocate (KfRcvStrt(0:jproc-1))
      allocate (KfRcvCnts(0:jproc-1))

      allocate (JrSndStrt(0:jproc-1))
      allocate (JrSndCnts(0:jproc-1))     
      allocate (JrRcvStrt(0:jproc-1))
      allocate (JrRcvCnts(0:jproc-1))

      allocate (KrSndStrt(0:iproc-1))
      allocate (KrSndCnts(0:iproc-1))     
      allocate (KrRcvStrt(0:iproc-1))
      allocate (KrRcvCnts(0:iproc-1))

!   start pointers and types of send  for the 1st forward transpose
      do i=0,iproc-1
         IfSndStrt(i) = (iist(i) -1)* jisize*kjsize*mytype*2
         IfSndCnts(i) = iisz(i) * jisize*kjsize*mytype*2

!   start pointers and types of recv for the 1st forward transpose
         IfRcvStrt(i) = (jist(i) -1) * iisize*kjsize*mytype*2
         IfRcvCnts(i) = jisz(i) * iisize*kjsize*mytype*2
      end do

!   start pointers and types of send  for the 2nd forward transpose
      do i=0,jproc-1
         KfSndStrt(i) = (jjst(i) -1)*iisize*kjsize*mytype*2
         KfSndCnts(i) = iisize*kjsize*jjsz(i)*mytype*2

!   start pointers and types of recv for the 2nd forward transpose
         KfRcvStrt(i) = (kjst(i) -1) * iisize * jjsize*mytype*2
         KfRcvCnts(i) = iisize*jjsize*kjsz(i)*mytype*2
      end do

!   start pointers and types of send  for the 1st inverse transpose
      do i=0,jproc-1
         JrSndStrt(i) = (kjst(i) -1) * iisize * jjsize*mytype*2
         JrSndCnts(i) = iisize*jjsize*kjsz(i)*mytype*2

!   start pointers and types of recv for the 1st inverse transpose
         JrRcvStrt(i) = (jjst(i) -1)*iisize*kjsize*mytype*2
         JrRcvCnts(i) = jjsz(i) * iisize * kjsize*mytype*2
      end do

!   start pointers and types of send  for the 2nd inverse transpose
      do i=0,iproc-1
         KrSndStrt(i) = (jist(i) -1) * iisize*kjsize*mytype*2
         KrSndCnts(i) = jisz(i) * iisize*kjsize*mytype*2

!   start pointers and types of recv for the 2nd inverse transpose
         KrRcvStrt(i) = (iist(i) -1) * jisize*kjsize*mytype*2
         KrRcvCnts(i) = jisize*iisz(i)*kjsize*mytype*2
      enddo

      nm = iisize * ny_fft * (kjsize+padd)
      allocate(A(nm))
      allocate(B(nm*2))
      A = 0.0
      B = 0.0
      call init_fft(A,B,nm)
      deallocate(A)
      deallocate(B)

      end subroutine

!==================================================================       
      subroutine MapDataToProc (data,proc,st,en,sz)
c========================================================
!    
       implicit none
       integer data,proc,st(0:proc-1),en(0:proc-1),sz(0:proc-1)
       integer i,size,nl,nu

       size=data/proc
       nu = data - size * proc
       nl = proc - nu
       st(0) = 1
       sz(0) = size
       en(0) = size
       do i=1,nl-1
         st(i) = st(i-1) + size
         sz(i) = size
         en(i) = en(i-1) + size
      enddo
      size = size + 1
      do i=nl,proc-1
         st(i) = en(i-1) + 1
         sz(i) = size
         en(i) = en(i-1) + size
      enddo
      en(proc-1)= data 
      sz(proc-1)= data-st(proc-1)+1

      end subroutine

c========================================================

      subroutine init_fft(A,B,n1)

      use fft_spec
      implicit none

      complex(mytype) A(n1)
      real(mytype) B(n1*2)
      integer n1

      
      call init_work(nx_fft,ny_fft,nz_fft)
      call plan_f_r2c(B,nx_fft,A,nxhp,nx_fft,jisize*kjsize,.false.) 
      call plan_b_c2r(A,nxhp,B,nx_fft,nx_fft,jisize*kjsize,.false.) 
#ifdef STRIDE1
      call plan_f_c1(A,1,ny_fft,A,1,ny_fft,ny_fft,iisize*kjsize,.false.)
      call plan_b_c1(A,1,ny_fft,A,1,ny_fft,ny_fft,iisize*kjsize,.false.)
      call plan_f_c2(A,1,nz_fft, 
     &     A,1,nz_fft,nz_fft,iisize*jjsize,.false.)
      call plan_b_c2(A,1,nz_fft, 
     &     A,1,nz_fft,nz_fft,iisize*jjsize,.false.)
#else
      call plan_f_c1(A,iisize,1,A,iisize,1,ny_fft,iisize,.false.)
      call plan_b_c1(A,iisize,1,A,iisize,1,ny_fft,iisize,.false.)
      call plan_f_c2(A,iisize*jjsize, 1,
     &     A,iisize*jjsize, 1,nz_fft,iisize*jjsize,.false.)
      call plan_b_c2(A,iisize*jjsize, 1,
     &     A,iisize*jjsize, 1,nz_fft,iisize*jjsize,.false.)
#endif

      return
      end subroutine

c========================================================
c  3D FFT inverse transform with 2D domain decomposition
c
c  This version uses MPI_Alltoallv to exchange data between processors
c  In the second step, y planes are sent separately
c  The order of array elements in memory is the same in all stages: (x,y,z) 

c Input: XYZg - comlpex array, with y dimension contained entirely,
c               while x and y are block-distributed among processors in 2D grid
c Output: XgYZ - an array of real, x dimension is contained entirely within 
c processors memory  while z and y are block-distributed among processors 
c in 2D grid
c
c !!! CAUTION: In this version: all arrays may occupy the same memory space
c In this case their first elements should coincide. 
c Naturally, output overwrites input
c

      subroutine btran_c2r (XYZg,XgYZ)
c========================================================

      use fft_spec
      implicit none

      real(mytype) XgYZ(nx_fft,jistart:jiend,kjstart:kjend)
      complex(mytype) XYZg(iistart:iiend,jjstart:jjend,nz_fft)
      complex(mytype), allocatable :: XYgZ(:,:,:)

      integer x,y,z,i,nx,ny,nz
c      logical flg_inplace

      if(.not. mpi_set) then
         print *,'P3DFFT error: call setup before other routines'
         return
      endif
      nx = nx_fft
      ny = ny_fft
      nz = nz_fft

c      call init_work(nx,ny,nz)

      allocate(XYgZ(nxhp,jistart:jiend,kjstart:kjend+padd))

      if(jproc .gt. 1) then

         call init_b_c(XYZg, iisize*jjsize, 1,
     &        XYZg, iisize*jjsize, 1,nz,iisize*jjsize)
c Transform in z dimension for all x and y

         call exec_b_c2(XYZg, iisize*jjsize, 1,
     &        XYZg, iisize*jjsize, 1,nz,iisize*jjsize)
         
c Exhange z-y pencils for all x

         call bcomm1(XYZg,XYgZ)

      else
c In case of 1D decomposition do only FFT in Z without data exchange

         call init_b_c(XYZg,iisize*ny,1,XYgZ,
     &        iisize*ny, 1, nz,iisize*ny)
         call exec_b_c2(XYZg,iisize*ny,1,XYgZ,
     &        iisize*ny,1, nz,iisize*ny)
      endif

c      call clean_x3

      if(iproc .gt. 1) then
c
c Transform in y dimension for all x, one z-plane at a time
c
         call init_b_c(XYgZ,iisize,1,XYgZ,iisize,1,ny,iisize)
         
         do z=kjstart,kjend

               call btran_y_zplane(XYgZ,z-kjstart,iisize,kjsize,iisize,1,
     &            XYgZ,z-kjstart,iisize,kjsize,iisize,1,ny,iisize) 
            
         enddo

         call bcomm2(XYgZ,XYgZ)

      else

         call init_b_c(XYgZ,nxhp,1,XYgZ,nxhp,1,ny,nxhp)

c In case of 1D decomposition do only FFT in Y without data exchange
         do z=kjstart,kjend

            call btran_y_zplane(XYgZ,z-kjstart,nxhp,kjsize,nxhp,1,XYgZ,
     &           z-kjstart,nxhp,kjsize, nxhp,1,ny,nxhp)

         enddo

      endif

c      call clean_x2

c Perform Complex-to-real FFT in x dimension for all y and z

      call init_b_c2r(XYgZ,nxhp,XgYZ,nx,nx,jisize*kjsize)
      call exec_b_c2r(XYgZ,nxhp,XgYZ,nx,nx,jisize*kjsize)

c      call clean_x1

      deallocate(XYgZ)

c      call free_work

      return
      end subroutine



c========================================================
c 3D FFT Forward transform, 2D domain decomposition
c
c  This version uses MPI_Alltoallv to exchange data between processors
c  In the first step, y planes are sent separately
c  The order of array elements in memory is the same in all stages: (x,y,z) 
c
c Input: XgYZ - an array of real, x dimension is contained entirely within 
c  processors memory  while z and y are block-distributed among processors 
c  in 2D grid XgYZ 
c Output: XYZg - comlpex array, with y dimension contained entirely,
c               while x and y are block-distributed among processors in 2D grid
c
c !!! CAUTION: In this version: all arrays may occupy the same memory space
c In this case their first elements should coincide. 
c Naturally, output overwrites input

      subroutine ftran_r2c (XgYZ,XYZg)
c========================================================

      use fft_spec
      implicit none

      real(mytype) XgYZ(nx_fft,jistart:jiend,kjstart:kjend)
      complex(mytype) XYZg(iistart:iiend,jjstart:jjend,nz_fft)
      complex(mytype), allocatable :: XYgZ(:,:,:)
      integer x,y,z,i,err,nx,ny,nz
c      logical flg_inplace
      
      if(.not. mpi_set) then
         print *,'P3DFFT error: call setup before other routines'
         return
      endif

      nx = nx_fft
      ny = ny_fft
      nz = nz_fft
c
c FFT transform in x for all z and y
c

      allocate(XYgZ(nxhp,jisize,kjsize+padd),stat=err)
      if(err .ne. 0) then
         print *,'Error ',err,' allocating array XYgZ'
      endif

c      call init_work(nx,ny,nz)

      call init_f_r2c(XgYZ,nx,XYgZ,nxhp,nx,jisize*kjsize) 

      call exec_f_r2c(XgYZ,nx,XYgZ,nxhp,nx,jisize*kjsize)

      if(iproc .gt. 1) then
         call fcomm1(XYgZ,XYgZ)

         call init_f_c(XgYZ,iisize,1,XgYZ,iisize,1,ny,iisize)
         do z=1,kjsize
            call ftran_y_zplane(XYgZ,z-1,iisize,kjsize,iisize,1,
     &          XYgZ,z-1,iisize,kjsize,iisize,1,ny,iisize)
         enddo
      else
         call init_f_c(XYgZ,nxhp,1,XYZg,nxhp,1,ny,nxhp)
         do z=1,kjsize

c In case of 1D decomposition do only FFT in Y without data exchange
            call ftran_y_zplane(XYgZ,z-1,nxhp,kjsize,nxhp,1,
     &            XYgZ,z-1,nxhp,kjsize,nxhp,1,ny,nxhp)
         enddo

      endif

c      call clean_x1
c      call clean_x2

      if(jproc .gt. 1) then

         call fcomm2(XYgZ,XYZg)

c Transform in z for all x and y
         call init_f_c(XYZg,iisize*jjsize, 1,
     &   XYZg,iisize*jjsize, 1,nz,iisize*jjsize)

         call exec_f_c2(XYZg,iisize*jjsize, 1,
     &   XYZg,iisize*jjsize, 1,nz,iisize*jjsize)

      else
c In case of 1D decomposition do only FFT in Y without data exchange
         call init_f_c(XYZg,iisize*ny,1,XYZg,
     &    iisize*ny,1, nz,iisize*ny)
         call exec_f_c2(XYgZ,iisize*ny,1,XYZg,
     &    iisize*ny,1, nz,iisize*ny)
      endif
      
      deallocate(XYgZ)
c      call clean_x3
c      call free_work

      return
      end subroutine


c Communication Module
c 
c Contains 4 routines for forward and backward exchanges inrows and columns 
c of processors. Uses MPI_Alltoallv routine
c

c========================================================
c Transpose X and Y pencils
c ! Assume source and dest are different !

      subroutine fcomm1(source,dest)
c========================================================

      implicit none

      complex(mytype) source(nxhp,jisize,kjsize)
#ifdef STRIDE1
      complex(mytype) dest(ny_fft,iisize,kjsize)
#else
      complex(mytype) dest(iisize,ny_fft,kjsize)
#endif

      real(8) t1,t2,t3,t4,tp1
      common /timers/ t1,t2,t3,t4,tp1
      complex(mytype) buf1(jisize*nxhp*kjsize),buf2(iisize*ny_fft*kjsize)
      integer x,y,i,ierr,z,xs
      integer*8 position

   
c Pack the real part of nx/2 frequency into the imaginary part of 0 frequency
c (both have imaginary part zero due to conjugte symmetry of a real-to-complex 
c transform)

      t1 = t1 - MPI_Wtime() 

      tp1 = tp1 - MPI_Wtime() 

c Pack the send buffer for exchanging y and x (within a given z plane ) into sendbuf

      position = 1

#ifdef PACK
      do y=1,jisize
         source(1,y,z) = cmplx(real(source(1,y,z)),real(source(nxhp,y,z)),mytype)
      enddo
#endif
      do i=0,iproc-1
         do z=1,kjsize
            do y=1,jisize
               do x=iist(i),iien(i)
                  buf1(position) = source(x,y,z)
                  position = position +1
               enddo
            enddo
         enddo
      enddo
      tp1 = tp1 + MPI_Wtime() 
      
c Exchange the y-x buffers

      call mpi_alltoallv(buf1,IfSndCnts, IfSndStrt,mpi_byte, 
     &     buf2,IfRcvCnts, IfRcvStrt,mpi_byte,mpi_comm_row,ierr)

      position = 1
      tp1 = tp1 - MPI_Wtime() 
      do i=0,iproc-1
         do z=1,kjsize
            do y=jist(i),jien(i)
               do x=1,iisize
#ifdef STRIDE1
                  dest(y,x,z) = buf2(position)
#else
                  dest(x,y,z) = buf2(position)
#endif
                  position = position + 1
               enddo
            enddo
         enddo
      enddo
      tp1 = tp1 + MPI_Wtime() 
      t1 = MPI_Wtime() + t1

      return
      end subroutine

c========================================================
c Transpose Y and Z pencils

      subroutine fcomm2(source,dest)
c========================================================

      implicit none

#ifdef STRIDE1
      complex(mytype) source(ny_fft,iisize,kjsize)
      complex(mytype) dest(nz_fft,iisize,jjsize)
#else
      complex(mytype) source(iisize,ny_fft,kjsize)
      complex(mytype) dest(iisize,jjsize,nz_fft)
#endif
      real(8) t1,t2,t3,t4,tp1
      common /timers/ t1,t2,t3,t4,tp1
      complex(mytype) buf1(iisize*(kjsize+padd)*ny_fft)
#ifdef STRIDE1
      complex(mytype) buf2(iisize*(jjsize+padd)*nz_fft)
#endif
      integer x,z,y,i,ierr,xs,ys
      integer*8 position

c Pack send buffers for exchanging y and z for all x at once 
      position = 1
      t2 = t2 - MPI_Wtime()
 
      do i=0,jproc-1
#ifdef STRIDE1
         do y=jjst(i),jjen(i)
            do x=1,iisize
               do z=1,kjsize
                  buf1(position) = source(y,x,z)
                  position = position+1
               enddo
            enddo
         enddo
#else
         do z=1,kjsize
            do y=jjst(i),jjen(i)
               do x=1,iisize
                  buf1(position) = source(x,y,z)
                  position = position+1
               enddo
            enddo
         enddo
#endif
      enddo
      
c Exchange y-z buffers


#ifdef STRIDE1
      call mpi_alltoallv(buf1,KfSndCnts, KfSndStrt,mpi_byte,
     &     buf2,KfRcvCnts, KfRcvStrt,mpi_byte,mpi_comm_col,ierr)

      position = 1
      do i=0,jproc-1
         do y=1,jjsize
            do x=1,iisize
               do z=kjst(i),kjen(i)
                  dest(z,x,y) = buf2(position)
                  position = position +1
               enddo
            enddo
         enddo
      enddo

#else

      call mpi_alltoallv(buf1,KfSndCnts, KfSndStrt,mpi_byte,
     &     dest,KfRcvCnts, KfRcvStrt,mpi_byte,mpi_comm_col,ierr)
#endif
      t2 = MPI_Wtime() + t2
         
      return
      end subroutine


c========================================================
c Transpose back Z to Y pencils

      subroutine bcomm1 (source,dest)
c========================================================

      implicit none

#ifdef STRIDE1
      complex(mytype) source(nz_fft,iisize,jjsize)
      complex(mytype) dest(ny_fft,iisize,kjsize)
      complex(mytype) buf1(iisize*jjsize*nz_fft)
#else
      complex(mytype) source(iisize,jjsize,nz_fft)
      complex(mytype) dest(iisize,ny_fft,kjsize)
#endif
      real(8) t1,t2,t3,t4,tp1
      common /timers/ t1,t2,t3,t4,tp1
      complex(mytype) buf2(iisize*kjsize*ny_fft)
      integer x,y,z,i,ierr,xs,ys
      integer*8 position

      t3 = t3 - MPI_Wtime() 

#ifdef STRIDE1
      position=1
      do z=1,nz_fft
         do x=1,iisize
            do y=1,jjsize
               buf1(position) = source(z,x,y)
               position = position+1
            enddo
         enddo
      enddo

      call mpi_alltoallv(buf1,JrSndCnts, JrSndStrt,mpi_byte,
     &     buf2,JrRcvCnts, JrRcvStrt,mpi_byte,mpi_comm_col,ierr)
#else
      call mpi_alltoallv(source,JrSndCnts, JrSndStrt,mpi_byte,
     &     buf2,JrRcvCnts, JrRcvStrt,mpi_byte,mpi_comm_col,ierr)
#endif      
c Unpack receive buffers into XZgY

      position=1
      do i=0,jproc-1
         do z=1,kjsize
#ifdef STRIDE1
            do x=1,iisize
               do y=jjst(i),jjen(i)
                  dest(y,x,z) = buf2(position)
                  position = position+1
               enddo
            enddo
#else
            do y=jjst(i),jjen(i)
               do x=1,iisize
                  dest(x,y,z) = buf2(position)
                  position = position+1
               enddo
            enddo
#endif
         enddo
      enddo
      
      t3 = t3 + MPI_Wtime() 
      
      return
      end subroutine


c========================================================
c Transpose back Y to X pencils

      subroutine bcomm2(source,dest)
c========================================================

      implicit none

      complex(mytype) dest(nxhp,jisize,kjsize)
#ifdef STRIDE1
      complex(mytype) source(ny_fft,iisize,kjsize)
#else
      complex(mytype) source(iisize,ny_fft,kjsize)
#endif
      real(8) t1,t2,t3,t4,tp1
      common /timers/ t1,t2,t3,t4,tp1
      complex(mytype) buf1(iisize*ny_fft*kjsize)
      complex(mytype) buf2(jisize*nxhp*kjsize)
      integer x,y,z,i,ierr
      integer*8 position

c Exchange x-z buffers

      t4 = t4 - MPI_Wtime() 
      
      position = 1
      do i=0,iproc-1
         do z=1,kjsize
            do y=jist(i),jien(i)
               do x=1,iisize
#ifdef STRIDE1
                  buf1(position) = source(y,x,z)
#else
                  buf1(position) = source(x,y,z)
#endif
                  position = position+1
               enddo
            enddo
         enddo
      enddo         
           
      call mpi_alltoallv (buf1,KrSndCnts, KrSndStrt,
     &           mpi_byte, buf2,KrRcvCnts,KrRcvStrt,mpi_byte,   
     &           mpi_comm_row,ierr)

c Unpack receive buffers into XgYZ

      position=1
      do i=0,iproc-1
         do z=1,kjsize
            do y=1,jisize
               do x=iist(i),iien(i)
                  dest(x,y,z) = buf2(position)
                  position = position+1
               enddo
            enddo
         enddo
      enddo

c Move the real part of nx/2 frequency to its proper location to prepare
c for the inverse X transform

#ifdef PACK
         do y=1,jisize
            dest(nxhp,y,z)=cmplx(imag(dest(1,y,z)),0.0,mytype)
            dest(1,y,z) = cmplx(real(dest(1,y,z)),0.0,mytype)
         enddo
#endif
      
      t4 = t4 + MPI_Wtime() 

      return
      end subroutine

c========================================================
      subroutine p3dfft_clean

!!--------------------------------------------------------------
c Clean-up routines for FFTW

      use fft_spec

#ifdef FFTW
#ifdef DOUBLE_PREC
      call dfftw_destroy_plan(plan1_frc)      
      call dfftw_destroy_plan(plan1_bcr)      
      call dfftw_destroy_plan(plan1_fc)      
      call dfftw_destroy_plan(plan2_fc)      
      call dfftw_destroy_plan(plan1_bc)      
      call dfftw_destroy_plan(plan2_bc)      
#else
      call sfftw_destroy_plan(plan1_frc)      
      call sfftw_destroy_plan(plan1_bcr)      
      call sfftw_destroy_plan(plan1_fc)      
      call sfftw_destroy_plan(plan2_fc)      
      call sfftw_destroy_plan(plan1_bc)      
      call sfftw_destroy_plan(plan2_bc)      
#endif

#elif defined ESSL
      deallocate(caux1)      
      deallocate(caux2)
      deallocate(raux1)
      deallocate(raux2)
#endif

      return
      end subroutine

      end module







