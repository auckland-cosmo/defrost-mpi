! $Id: defrost.f90,v 2.0 2008/09/30 01:23:01 frolov Stab $
! [compile with: ifort -O3 -ipo -xT -r8 -pc80 -fpp defrost.f90 -lfftw3]

! Reheating code doing something...
! http://www.sfu.ca/physics/cosmology/defrost

! Copyright (C) 2007 Andrei Frolov <frolov@sfu.ca>
! Distributed under the terms of the GNU General Public License
! If you use this code for your research, please cite arXiv:0809.4904

! Ported to GNU compiler and MPI-enabled by Hal Finkel 2008-2010

program defrost

use mpi; use p3dfft;
#ifdef _OPENMP
use omp_lib
#endif

implicit none

include "fftw3.f"

#ifdef SILO
include "silo.inc"
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! some useful constants
real, parameter :: twopi = 6.2831853071795864769252867665590
real, parameter :: sqrt3 = 1.7320508075688772935274463415059

#include "parameters.inc"

! buffers holding variable names, statistics and spectra (on per-frame basis)
character(12) DVAR(fields+3); real CDF(n+1,fields+3), PSD(ns,fields+3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer l
integer :: mpierror, mpisize, mpirank, mpistatus(MPI_STATUS_SIZE)
integer :: dims(2)
integer :: istart(3), iend(3), isize(3)
integer :: sistart(3), siend(3), sisize(3)
integer :: fstart(3), fend(3), fsize(3)

integer :: zdivs, ydivs                ! Number of divisions in the z and y dimensions
integer :: lzdiv, lydiv                ! Local division number in z and y
integer :: leftzrank, leftyrank        ! The left rank numbers (cyclic)
integer :: rightzrank, rightyrank      ! The right rank numbers (cyclic)

! smp array samples all fields on a 3D grid for three subsequent time slices
real, allocatable :: smp(:,:,:,:,:), pp(:,:,:,:), tmp(:,:,:), tmp2(:,:,:), bov(:,:,:); complex, allocatable :: Fk(:,:,:)

! init MPI
call MPI_Init(mpierror)
call MPI_Comm_size(MPI_COMM_WORLD, mpisize, mpierror)
call MPI_Comm_rank(MPI_COMM_WORLD, mpirank, mpierror)

! total number of divisions, divide z up to 
ydivs = (mpisize/n) + 1; zdivs = (mpisize/ydivs)

! local division numbers
lydiv = mpirank/zdivs; lzdiv = mpirank - (lydiv*zdivs)

! right y rank
if (lydiv .eq. (ydivs - 1)) then
        rightyrank = lzdiv
else
        rightyrank = mpirank + zdivs
endif

! left y rank
if (lydiv .eq. 0) then
        leftyrank = (ydivs-1)*zdivs + lzdiv
else
        leftyrank = mpirank - zdivs
endif

! right z rank
if (lzdiv .eq. (zdivs - 1)) then
        rightzrank = mpirank - (zdivs - 1)
else
        rightzrank = mpirank + 1
endif

! left z rank
if (lzdiv .eq. 0) then
        leftzrank = mpirank + (zdivs - 1)
else
        leftzrank = mpirank - 1
endif

! setup the FFT dims array
dims(1) = ydivs; dims(2) = zdivs

! use threaded FFTW on SMP machines (link with -lfftw3_threads)
! note: on some machines it may be helpful to compile FFTW with --enable-openmp
#ifdef _OPENMP
#ifdef USE_FFTW_THREADS
call dfftw_init_threads
call dfftw_plan_with_nthreads(OMP_GET_MAX_THREADS())
#endif
#endif

! init FFT
call p3dfft_setup (dims, n, n, n)
call get_dims(istart, iend, isize, 1)
call get_dims(fstart, fend, fsize, 2)

! storage size and indicies
sistart(:) = istart(:) - 1
siend(:) = iend(:) + (p-n)
sisize(:) = siend(:) - sistart(:) + 1

allocate(smp(fields,sistart(1):siend(1),sistart(2):siend(2),sistart(3):siend(3),3))
allocate(pp(2,sistart(1):siend(1),sistart(2):siend(2),sistart(3):siend(3)))
allocate(tmp(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)))
allocate(tmp2(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)))
allocate(Fk(fstart(1):fend(1),fstart(2):fend(2),fstart(3):fend(3)))
!if (output$bov .and. nx /= n) allocate(bov(nx,nx,nx))
if (output$bov .and. ds /= 1) allocate(bov(isize(1)/ds,isize(2)/ds,isize(3)/ds))

#ifndef SEED
#define SEED 0
#endif

! initialize random number generator
call random_seed(SIZE=rssize)
allocate(rs(rssize))
rs=SEED+mpirank
call random_seed(PUT=rs)

! initialize and run simulation
call head(6, (/"t    ", "a    ", "H    ", "<rho>", "<P>  "/))
call init(smp(:,:,:,:,1), smp(:,:,:,:,2))

do l = 1,tt,3
        if (LA(2) .le. amax) then
                if (LA(2) .gt. aflip .and. fliplambda) then
                        lambdapf = -1.0
                end if
                if (LA(2) .gt. aflip .and. flipg) then
                        gpf = 0.0
                end if

                call step(l,   smp(:,:,:,:,1), smp(:,:,:,:,2), smp(:,:,:,:,3))
                call step(l+1, smp(:,:,:,:,2), smp(:,:,:,:,3), smp(:,:,:,:,1))
                call step(l+2, smp(:,:,:,:,3), smp(:,:,:,:,1), smp(:,:,:,:,2))
        end if
end do

call p3dfft_clean()
call MPI_Finalize(mpierror)

contains

#include "model.inc"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! periodic boundary conditions
subroutine wrap(up)
        real, dimension(fields,sistart(1):siend(1),sistart(2):siend(2),sistart(3):siend(3)) :: up

        ! the x dimension is non-distributed        
        up(:,0,:,:) = up(:,n,:,:); up(:,n+1,:,:) = up(:,1,:,:)

        ! the y dimension may be distributed
        if (ydivs .eq. 1) then
                up(:,:,0,:) = up(:,:,n,:); up(:,:,n+1,:) = up(:,:,1,:)
        else
                call MPI_Sendrecv(up(:,:,iend(2),:), fields*sisize(1)*sisize(3), MPI_DOUBLE_PRECISION, rightyrank, 0, &
                       up(:,:,istart(2)-1,:), fields*sisize(1)*sisize(3), MPI_DOUBLE_PRECISION, leftyrank, 0, &
                       MPI_COMM_WORLD, mpistatus, mpierror)
                call MPI_Sendrecv(up(:,:,istart(2),:), fields*sisize(1)*sisize(3), MPI_DOUBLE_PRECISION, leftyrank, 0, &
                       up(:,:,iend(2)+1,:), fields*sisize(1)*sisize(3), MPI_DOUBLE_PRECISION, rightyrank, 0, &
                       MPI_COMM_WORLD, mpistatus, mpierror)
        endif

        ! the z dimension is distributed
        call MPI_Sendrecv(up(:,:,:,iend(3)), fields*sisize(1)*sisize(2), MPI_DOUBLE_PRECISION, rightzrank, 0, &
               up(:,:,:,istart(3)-1), fields*sisize(1)*sisize(2), MPI_DOUBLE_PRECISION, leftzrank, 0, &
               MPI_COMM_WORLD, mpistatus, mpierror)
        call MPI_Sendrecv(up(:,:,:,istart(3)), fields*sisize(1)*sisize(2), MPI_DOUBLE_PRECISION, leftzrank, 0, &
               up(:,:,:,iend(3)+1), fields*sisize(1)*sisize(2), MPI_DOUBLE_PRECISION, rightzrank, 0, &
               MPI_COMM_WORLD, mpistatus, mpierror)
end subroutine wrap

! scalar field initial conditions
subroutine init(dn, hr)
        real, dimension(fields,sistart(1):siend(1),sistart(2):siend(2),sistart(3):siend(3)) :: dn, hr
        real, dimension(fields) :: m2eff

        m2eff = modelm2eff()
        m2eff = m2eff + m2phi - 2.25*H0**2
        
        call sample(tmp, -0.25, m2eff(phi)); tmp = initscale*tmp
        hr(phi,istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)) = tmp + phi0

        call sample(tmp, +0.25, m2eff(phi)); tmp = initscale*tmp
        dn(phi,istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)) = &
                hr(phi,istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)) - tmp*dt - dphi0*dt + ddphi0*dt**2/2.0
        
        call wrap(dn); call wrap(hr)
end subroutine init

! stencil operators (implemented as preprocessor macros)
#define HR(x,y,z) hr(:,i+(x),j+(y),k+(z))
#define GRAD2(x,y,z) sum((hr(:,i+(x),j+(y),k+(z))-hr(:,i,j,k))**2)

#define RANK0(O) (O(0,0,0))
#define RANK1(O) (O(-1,0,0) + O(1,0,0) + O(0,-1,0) + O(0,1,0) + O(0,0,-1) + O(0,0,1))
#define RANK2(O) (O(-1,-1,0) + O(1,-1,0) + O(-1,1,0) + O(1,1,0) + O(-1,0,-1) + O(1,0,-1) + O(-1,0,1) + O(1,0,1) + O(0,-1,-1) + O(0,1,-1) + O(0,-1,1) + O(0,1,1))
#define RANK3(O) (O(-1,-1,-1) + O(1,-1,-1) + O(-1,1,-1) + O(1,1,-1) + O(-1,-1,1) + O(1,-1,1) + O(-1,1,1) + O(1,1,1))

! pgf90 has a line length limit of 2048 which is too short for the expanded stencil
#ifdef __PGI
#define STENCIL(C,O) ((C/**/0) * RANK0(O) + (C/**/1) * RANK1(O) + (C/**/2) * RANK2(O) + (C/**/3) * RANK3(O))
#else
#ifdef __GFORTRAN__
#define STENCIL(C,O) ((C/**/0) * RANK0(O) + (C/**/1) * RANK1(O) + (C/**/2) * RANK2(O) + (C/**/3) * RANK3(O))
#else
#define STENCIL(C,O) ((C ## 0) * RANK0(O) + (C ## 1) * RANK1(O) + (C ## 2) * RANK2(O) + (C ## 3) * RANK3(O))
#endif
#endif

! scalar field evolution step
subroutine step(l, dn, hr, up)
        real, dimension(fields,sistart(1):siend(1),sistart(2):siend(2),sistart(3):siend(3)) :: dn, hr, up; integer i, j, k, l
        
        ! Laplacian operator stencils: traditional one and three isotropic variants
        ! stable for dx/dt > sqrt(3), sqrt(2), sqrt(21)/3, and 8/sqrt(30) respectively
        ! computational cost difference is insignificant for large grids
        
        !real, parameter :: c3 = 0.0, c2 = 0.0, c1 = 1.0, c0 = -6.0, cc = 1.0
        !real, parameter :: c3 = 0.0, c2 = 1.0, c1 = 2.0, c0 = -24.0, cc = 6.0
        !real, parameter :: c3 = 1.0, c2 = 0.0, c1 = 8.0, c0 = -56.0, cc = 12.0
        real, parameter :: c3 = 1.0, c2 = 3.0, c1 = 14.0, c0 = -128.0, cc = 30.0
        
        ! field energy (distributed for parallelization)
        real :: V, T, G
        real, dimension(n) :: PE, KE, GE
        
        ! various evolution operator coefficients
        real c, d, b0, b1, b2, b3, d1, d2, e1, e2, e3
        
        ! optional computations flags
        logical, parameter :: dumping = output .and. (output$bov .or. output$crv)
        logical, parameter :: needTii = dumping .and. (output$set .or. output$pot)
        logical checkpt; integer db, idx
        
        ! flat or expanding background
        !real, parameter :: a = 1.0, H = 0.0; real Q
        real a, H, Q, R; a = LA(2); H = 1.0/LH(2)
        
        ! all coefficients inside the loop are pre-calculated here
        d = 1.0 + 1.5*H*dt; c = cc * alpha**2 * a**2 * d
        b0 = 2.0/d + c0/c; b1 = c1/c; b2 = c2/c; b3 = c3/c
        d1 = -(1.0 - 1.5*H*dt)/d; d2 = -dt**2/d
        e1 = 1.0/(8.0*dt**2); e2 = 1.0/(4.0*(a*dx)**2*cc); e3 = e2/3.0
        
        ! initialize accumulators
        PE = 0.0; KE = 0.0; GE = 0.0
        checkpt = mod(l-1, nt) == 0
        db = 0; idx = 0
        
        ! discretized scalar field evolution step
        !$omp parallel do
        do k = istart(3),iend(3); do j = istart(2),iend(2); do i = istart(1),iend(1)
                ! scalar field potential derivatives are inlined here
                up(:,i,j,k) = STENCIL(b,HR) + d1 * dn(:,i,j,k) + &
                    d2 * ( modeldv(hr, i, j, k) ) * hr(:,i,j,k)
                
                ! scalar field potential multiplied by 2 is inlined here
                V = modelv(hr, i, j, k)
                T = sum((up(:,i,j,k)-dn(:,i,j,k))**2)

                PE(k) = PE(k) + V; KE(k) = KE(k) + T
                
                ! calculate density and pressure when needed
                if (checkpt) then
                        G = STENCIL(c,GRAD2); GE(k) = GE(k) + G
                        
                        if (needTii) then
                                pp(rho,i,j,k) = e1*T + e2*G + 0.5*V
                                pp(prs,i,j,k) = e1*T - e3*G - 0.5*V
                        end if
                end if
        end do; end do; end do

        call MPI_Allreduce(MPI_IN_PLACE, PE, n, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
        call MPI_Allreduce(MPI_IN_PLACE, KE, n, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
        call MPI_Allreduce(MPI_IN_PLACE, GE, n, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
        
        ! periodic boundary conditions
        call wrap(up)
        
        ! update expansion factors
        Q = sum(4.0*e1*KE - PE)/(6.0*n**3)
        R = LH(1) + (1.0 + Q * LH(2)**2) * dt
        LH = (/ LH(2), LH(1) + (1.0 + Q * R**2) * (2.0*dt) /)
        LA = (/ LA(2), LA(1) + (H*a) * (2.0*dt) /)
        
        ! dump simulation data
        if (checkpt) then
                if (mpirank .eq. 0) write (*,'(5g25.16e2)') (l-1)*dt, a, H, sum(e1*KE + e2*GE + 0.5*PE)/n**3, &
                        sum(e1*KE - e3*GE - 0.5*PE)/n**3
                
                if (dumping .and. output$any) db = fopen("frame", (l-1)/nt, (l-1)*dt)
                if (dumping .and. output$fld) then
                        Q = 1.0; if (oscale) Q = a**1.5
                        tmp = Q*hr(phi,istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
                        call dump(db, "phi", (l-1)/nt, (l-1)*dt, tmp, idx)
                end if
                if (dumping .and. output$set) then
                        Q = 1.0; if (oscale) Q = 1.0/(3.0*H**2)
                        tmp = Q*pp(rho,istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
                        call dump(db, "rho", (l-1)/nt, (l-1)*dt, tmp, idx)

                        tmp = Q*pp(prs,istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
                        call dump(db, "prs", (l-1)/nt, (l-1)*dt, tmp, idx)
                end if
                if (dumping .and. output$pot) then
                        Q = a**2/2.0; tmp = Q*pp(rho,istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
                        call laplace(tmp, tmp); call dump(db, "PSI", (l-1)/nt, (l-1)*dt, tmp, idx)
                end if
                if (idx > 0 .and. output$gnu) call fflush((l-1)*dt, idx)
                if (db /= 0) call fclose(db)
        end if
end subroutine step


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! sample Gaussian random field with power-law spectrum
subroutine sample(f, gamma, m2eff)
        real f(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)), gamma, m2eff; integer*8 plan

#ifndef DONT_FIX_ICS
#define KCUT_FAC 2.0
#else
#define KCUT_FAC 1.0
#endif
        
        integer, parameter :: os = 16, nos = n * os**2
        real, parameter :: dxos = dx/os, dkos = dk/(2*os), kcut = KCUT_FAC*nn*dk/2.0
        complex, parameter :: w = (0.0, twopi)
        
        real ker(nos), a(nn), p(nn)
        integer i, j, k, l; real kk, norm
        
        ! calculate (oversampled) radial profile of convolution kernel
        do k = 1,nos; kk = (k-0.5)*dkos
                ker(k) = kk*(kk**2 + m2eff)**gamma * exp(-(kk/kcut)**2)
        end do
        
        call dfftw_plan_r2r_1d(plan,nos,ker,ker,FFTW_RODFT10,FFTW_ESTIMATE)
        call dfftw_execute(plan); call dfftw_destroy_plan(plan)
        
        norm = 0.5/(n**3 * sqrt(twopi*dk**3) * mpl) * (dkos/dxos)
        do k = 1,nos; ker(k) = norm * ker(k)/k; end do
        
        ! initialize 3D convolution kernel (using linear interpolation of radial profile)
        !$omp parallel do
        do k = istart(3),iend(3); do j = istart(2),iend(2); do i = istart(1),iend(1)
                kk = sqrt(real(i-nn)**2 + real(j-nn)**2 + real(k-nn)**2) * os; l = floor(kk)
                
                if (l > 0) then
                        f(i,j,k) = ker(l) + (kk-l)*(ker(l+1)-ker(l))
                else
                        f(i,j,k) = (4.0*ker(1)-ker(2))/3.0
                end if
        end do; end do; end do
        
        ! convolve kernel with delta-correlated Gaussian noise
        call p3dfft_ftran_r2c(f, Fk)
        
        !$omp parallel do
        do k = fstart(3),fend(3); do j = fstart(2),fend(2)
                call random_number(a); call random_number(p)
                Fk(:,j,k) = sqrt(-2.0*log(a)) * exp(w*p) * Fk(:,j,k)
        end do; end do

        call p3dfft_btran_c2r(Fk, f)        
end subroutine sample

! solve Laplace equation $\Delta f = \rho$
subroutine laplace(f, rho)
        real f(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)), rho(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
        integer i, j, k; real :: ii, jj, kk
        
        real, parameter :: w = twopi/n
        
        ! Laplacian operator stencils: traditional one and three isotropic variants
        ! (corresponding to discretization used for field evolution equations above)
        
        !real, parameter :: c3 = 0.0, c2 = 0.0, c1 = 1.0, c0 = -3.0, cc = 0.5
        !real, parameter :: c3 = 0.0, c2 = 1.0, c1 = 1.0, c0 = -6.0, cc = 1.5
        !real, parameter :: c3 = 1.0, c2 = 0.0, c1 = 2.0, c0 = -7.0, cc = 1.5
        real, parameter :: c3 = 2.0, c2 = 3.0, c1 = 7.0, c0 = -32.0, cc = 7.5
        
        real, parameter :: c = cc * dx**2/real(n)**3

        call p3dfft_ftran_r2c(rho, Fk)
        
        !$omp parallel do
        do k = fstart(3),fend(3); kk = cos(w*(k-1))
        do j = fstart(2),fend(2); jj = cos(w*(j-1))
        do i = fstart(1),fend(1); ii = cos(w*(i-1))
                Fk(i,j,k) = c*Fk(i,j,k)/(c0 + c1*(ii+jj+kk) + c2*(ii*jj+ii*kk+jj*kk) + c3*ii*jj*kk)
        end do; end do; end do
        
        Fk(1,1,1) = 0.0

        call p3dfft_btran_c2r(Fk, f)        
end subroutine laplace

! one-sided power spectrum density estimator (output valid only on rank-0 node)
subroutine spectrum(f, S)
        real f(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3)), S(ns), W(ns)
        integer i, j, k, ii, jj, kk, l; real p, c(2)

        call p3dfft_ftran_r2c(f, Fk)
        
        S = 0.0; W = 0.0
        
        do k = fstart(3),fend(3); if (k <= nn) then; kk = k-1; else; kk = n+1-k; end if
        do j = fstart(2),fend(2); if (j <= nn) then; jj = j-1; else; jj = n+1-j; end if
        do i = 1,n;               if (i <= nn) then; ii = i-1; else; ii = n+1-i; end if
                if ((ii+1) .lt. fstart(1) .or. (ii+1) .gt. fend(1)) cycle

                p = sqrt(real(ii**2 + jj**2 + kk**2)); l = floor(p)
                
                c = (1.0 - (/l-p,l+1-p/)**2)**2
                
                S(l+1:l+2) = S(l+1:l+2) + c * Fk(ii+1,j,k)*conjg(Fk(ii+1,j,k))
                W(l+1:l+2) = W(l+1:l+2) + c
        end do; end do; end do

        call MPI_Allreduce(MPI_IN_PLACE, S, ns, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
        call MPI_Allreduce(MPI_IN_PLACE, W, ns, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
        
        where (W /= 0.0) S = S/W/real(n)**6
end subroutine spectrum


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! local merge sort
subroutine mergesort(f, ft)
        real, dimension(:), intent(inout) :: f
        real, dimension(1:size(f)), intent(inout) :: ft
        call mergesortinner(f, ft)
end subroutine mergesort

recursive subroutine mergesortinner(f, ft)
        real, dimension(:), intent(inout) :: f
        real, dimension(1:size(f)), intent(inout) :: ft

        integer fs, mid, cur, left, right 

        fs = size(f)
        if (fs .lt. 2) return

        mid = fs/2
        call mergesortinner(f(1:mid), ft(1:mid))
        call mergesortinner(f((mid+1):fs), ft((mid+1):fs))

        cur = 1; left = 1
        right = mid + 1

        do while (cur .le. fs)
                if (left .gt. mid) then
                        ft(cur) = f(right)
                        right = right + 1
                else if (right .gt. fs) then
                        ft(cur) = f(left)
                        left = left + 1
                else
                        if (f(left) .le. f(right)) then
                                ft(cur) = f(left)
                                left = left + 1
                        else
                                ft(cur) = f(right)
                                right = right + 1
                        end if
                end if

                cur = cur + 1
        end do

        f = ft
end subroutine mergesortinner

! distributed selection sort
! see: Distributed selectsort sorting algorithms on broadcast communication networks,
!      Jau-Hsiung Huang and Leonard Kleinrock, Parallel Computing 16 (1990) 183-190

subroutine distsortflat(f, ft, ft2, fs)
        integer fs
        real, dimension(1:fs), intent(inout) :: f
        real, dimension(1:fs) :: ft, ft2

        call distsort(f, ft, ft2)
end subroutine

subroutine distsort(f, ft, ft2)
        real, dimension(:), intent(inout) :: f
        real, dimension(1:size(f)) :: ft, ft2

        ! each node sorts its local list
        call mergesort(f, ft)

        ! if (mpisize .lt. 2) return

        call distsortinner(f, ft, ft2)
end subroutine distsort

subroutine distsortinner(f, ft, ft2)
        real, dimension(:), intent(inout) :: f
        real, dimension(1:size(f)) :: ft, ft2

        integer fs
        integer rrank, ri, cd, cs
        integer bsi(0:(mpisize-1)), bsl(0:(mpisize-1)), bsr(0:(mpisize-1))
        integer bsc(0:(mpisize-1)), gbsc(0:(mpisize-1),0:(mpisize-1))
        integer bsf(0:(mpisize-1)), gbsf(0:(mpisize-1)), bspn(0:(mpisize-1))
        integer bsrank(0:(mpisize-1))

        real bsv(0:(mpisize-1)), gbsv(0:(mpisize-1),0:(mpisize-1))

        fs = size(f)
        bspn = mpisize

        bsi = fs/2; bsl = 1; bsr = fs; bsf = 0; gbsf = 0

        ! each node is doing mpisize binary searches (for the starting value
        ! for each node)
        do while (sum(gbsf) .lt. mpisize)
                ! for each binary search, find the number of elements less than
                ! the current value

                ! note: the sorting needs to "break the tie" when duplicate
                ! values are in the list; as a result, it really sorts the
                ! tuples (value, rank, index), but only the values need be
                ! compared if they are different.

                ! get the values of the current index in each binary search
                do rrank = 0,(mpisize-1)
                        bsv(rrank) = f(bsi(rrank))
                end do

                ! send those values to all processors
                call MPI_Allgather(bsv, mpisize, MPI_DOUBLE_PRECISION, gbsv, mpisize*mpisize,&
                                   MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, mpierror)

                ! count the number of values less than those in gbsv
                gbsc = 0
                do cs = 1,fs
                        do rrank = 0,(mpisize-1); do ri = 0,(mpisize-1)
                                if (f(cs) .lt. gbsv(ri, rrank)) then
                                        gbsc(ri, rrank) = gbsc(ri, rrank) + 1
                                else if (f(cs) .eq. gbsv(ri, rrank) .and. mpirank .lt. rrank) then
                                        gbsc(ri, rrank) = gbsc(ri, rrank) + 1
                                else if (f(cs) .eq. gbsv(ri, rrank) .and. mpirank .eq. rrank .and. cs .lt. bsi(ri)) then
                                        gbsc(ri, rrank) = gbsc(ri, rrank) + 1
                                end if
                        end do; end do
                end do

                ! add the results, and send those counts back to each node
                call MPI_Reduce_scatter(gbsc, bsc, bspn, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierror)

                ! update each binary search
                do rrank = 0,(mpisize-1)
                        if (bsc(rrank) .eq. fs*rrank) then
                                ! this value is the target starting element, mark as found
                                bsf(rrank) = 1
                        else
                                if (bsc(rrank) .lt. fs*rrank) then
                                        ! search in the upper half of the current range
                                        bsl(rrank) = bsi(rrank)
                                else
                                        ! search in the lower half of the current range
                                        bsr(rrank) = bsi(rrank)
                                end if

                                bsi(rrank) = (bsl(rrank) + bsr(rrank))/2
                        end if
                end do

                ! count the number of complete searches
                call MPI_Allreduce(bsf, gbsf, mpisize, MPI_INTEGER, MPI_BOR, MPI_COMM_WORLD, mpierror)
        end do

        ! get the search-result values
        bsv = 0.0; bsrank = 0
        do rrank = 0,(mpisize-1)
                if (bsf(rrank) .eq. 1) then
                        bsv(rrank) = f(bsi(rrank))
                        bsrank(rrank) = rrank
                else
                        bsi(rrank) = 0
                end if
        end do

        ! send the values found to all nodes
        call MPI_Allreduce(MPI_IN_PLACE, bsv, mpisize, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
        call MPI_Allreduce(MPI_IN_PLACE, bsi, mpisize, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierror)
        call MPI_Allreduce(MPI_IN_PLACE, bsrank, mpisize, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, mpierror)

        ft2 = 0.0

        ! send all elements to other nodes, keep the ones which belong to this node
        cd = 1
        do rrank = 0,(mpisize-1)
                if (rrank .eq. mpirank) then
                        ft2 = f
                end if
                        
                call MPI_Bcast(ft2, fs, MPI_DOUBLE_PRECISION, rrank, MPI_COMM_WORLD, mpierror)

                ! copy values to keep into ft
                do cs = 1,fs
                        if (ft2(cs) .lt. bsv(mpirank)) cycle

                        if (ft2(cs) .eq. bsv(mpirank)) then
                                if (rrank .lt. bsrank(mpirank)) cycle

                                if (rrank .eq. bsrank(mpirank)) then
                                        if (cs .lt. bsi(mpirank)) cycle

                                        if (mpirank .lt. (mpisize-1)) then
                                                ! last node has no upper bound
                                                if (cs .ge. bsi(mpirank+1)) cycle
                                        end if
                                end if

                                if (mpirank .lt. (mpisize-1)) then
                                        ! last node has no upper bound
                                        if (rrank .ge. bsrank(mpirank+1)) cycle
                                end if
                        end if

                        if (mpirank .lt. (mpisize-1)) then
                                ! last node has no upper bound
                                if (ft2(cs) .ge. bsv(mpirank+1)) cycle
                        end if

                        if (cd .le. fs) then
                                ft(cd) = ft2(cs)
                                cd = cd + 1
                        else
                                call MPI_Abort(MPI_COMM_WORLD, 1, mpierror);
                                stop 'Too many elements for local list part in sort'
                        end if
                end do
        end do

        call mergesort(ft, ft2)
        f = ft
end subroutine distsortinner

! partially sort an array, finding (r:n:m)th smallest elements
! recursive subroutine sieve(f, n, m, r)
!         integer i, j, k, n, m, r; real p, f(n)
!         
!         if (r > n) return
!         
!         i = 1; k = n/2+1; j = n
!         
!         if (f(1) > f(n)) f((/1,n/)) = f((/n,1/))
!         if (f(1) > f(k)) f((/1,k/)) = f((/k,1/))
!         if (f(n) < f(k)) f((/n,k/)) = f((/k,n/))
!         
!         if (n < 4) return
!         
!         p = f(k)
!         
!         do while (i < j)
!                 i = i+1; do while (f(i) < p); i = i+1; end do
!                 j = j-1; do while (f(j) > p); j = j-1; end do
!                 if (i < j) f((/i,j/)) = f((/j,i/))
!         end do
!         
!         call sieve(f(1:i-1), i-1, m, r)
!         call sieve(f(j+1:n), n-j, m, modulo(r-j-1,m)+1)
! end subroutine sieve


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! identify yourself
subroutine head(fd, vars)
        integer(4) fd; character(*) :: vars(:)
        character(512) :: buffer; integer a, b, c, l
        character(*), parameter :: rev = "$Revision: 2.0 $"

        if (mpirank .ne. 0) return
        
        ! behold the horror that is Fortran string parsing
        a = index(rev, ": ") + 2
        b = index(rev, ".")
        c = index(rev, " $", .true.) - 1
        buffer = rev(a:b-1); read (buffer, '(i12)') a
        buffer = rev(b+1:c); read (buffer, '(i12)') b
        
        ! ID string
        write (fd,'(a,4(i0,a))') "# This is DEFROST revision ", a-1, ".", b, " (", fields, " fields, ", n, "^3 grid)"

#ifdef _OPENMP
        write (fd,'(a,i0)') "# Using OpenMP with max threads: ", OMP_GET_MAX_THREADS()
#endif

        write (fd,'(a,i0,a,i0)') "# MPI rank: ", mpirank, "; total nodes: ", mpisize
        write (fd,'(a,i0,a,i0,a)') "# Using ", ydivs, " y division(s) and ", zdivs, " z division(s)"
        write (fd,'(a,i0,a,i0)') "# Local y division: ", lydiv, "; local z division: ", lzdiv
        write (fd,'(a,i0,a,i0)') "# neighbor ranks: ", leftyrank, " <- y -> ", rightyrank
        write (fd,'(a,i0,a,i0)') "# neighbor ranks: ", leftzrank, " <- z -> ", rightzrank

        ! model summary
        call modelsummary(fd)
        
        ! variable list
        write (buffer,'(a,g12.12",",32(g24.12","))') "OUTPUT:", adjustr(vars); l = index(buffer, ',', .true.);
        write (fd,'(2a)') "# ", repeat('=', l)
        write (fd,'(2a)') "# ", buffer(1:l-1)//';'
        write (fd,'(2a)') "# ", repeat('=', l)
end subroutine head

! open frame
function fopen(file, frame, t)
        character(*) :: file; integer fopen, frame; real t

#ifdef SILO
        integer db, opts, e
        integer, parameter :: D = 3, nodes(D) = nx+1
        integer, parameter :: lut(nx) = (/0:nx-1/)*n/nx
        real, parameter :: mesh(nx+1) = (/lut*dx,n*dx/)
        
        character(256) :: buffer; write (buffer,'(a,a,i4.4,a)') file, '-', frame, '.silo'

        if (mpirank .ne. 0) return        
        if (.not. (output$bov .or. output$vis)) return
        
#define STR(string) string, len(string)
        e = dbcreate(buffer, len_trim(buffer), DB_CLOBBER, DB_LOCAL, STR("DEFROST frame"), DB_HDF5, db)
        
        e = dbmkoptlist(3, opts)
        e = dbaddiopt(opts, DBOPT_CYCLE, frame)
        e = dbadddopt(opts, DBOPT_DTIME, t)
        e = dbaddiopt(opts, DBOPT_COORDSYS, DB_CARTESIAN)
        
        e = dbputqm(db, STR("mesh"), STR("x"), STR("y"), STR("z"), mesh, mesh, mesh, nodes, D, DB_DOUBLE, DB_COLLINEAR, opts, e)
        e = dbfreeoptlist(opts)
        
        fopen = db
#else
        fopen = 0
#endif
end function fopen

! close frame
subroutine fclose(db)
        integer db, e
        
#ifdef SILO
        if (mpirank .ne. 0) return        
        e = dbclose(db)
#endif
end subroutine fclose

! output field configuration and its aggregates
subroutine dump(db, v, frame, t, f, idx)
        character(*) :: v; integer db, frame, k; real t, f(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
        integer, optional :: idx
        character(256) :: buffer; real avg, var, S(ns), P(n+1), X(n+1), PDF(2:n);
        integer :: i, j, rawfh, datasz, ii, jj, kk
        integer(kind=MPI_OFFSET_KIND) disp
        
        !integer, parameter :: D = 3, zones(D) = nx, lut(nx) = (/ (i, i=0, nx-1) /)*n/nx + 1; integer e
        
        ! output 3D box of values
        if (output$bov) then
#ifdef SILO
                if (nx /= n) then
                        bov = f(lut,lut,lut)
                        e = dbputqv1(db, STR(v), STR("mesh"), bov, zones, D, DB_F77NULL, 0, DB_DOUBLE, DB_ZONECENT, DB_F77NULL, e)
                else
                        e = dbputqv1(db, STR(v), STR("mesh"), f, zones, D, DB_F77NULL, 0, DB_DOUBLE, DB_ZONECENT, DB_F77NULL, e)
                end if
#else
                if (mpirank .eq. 0) then
                        write (buffer,'(a,a,i4.4,a)') v, '-', frame, '.bov'; open(10, file=buffer)
                end if

                write (buffer,'(a,a,i4.4,a)') v, '-', frame, '.raw'
                call MPI_File_open(MPI_COMM_WORLD, buffer, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, rawfh, mpierror)

                ! raw data
                datasz = (isize(1)*isize(2)*isize(3))/(ds**3)
                disp = mpirank * sizeof(t) * datasz
                call MPI_File_set_view(rawfh, disp, MPI_DOUBLE_PRECISION, MPI_DOUBLE_PRECISION, 'native', MPI_INFO_NULL, mpierror)

                if (ds .eq. 1) then
                        call MPI_File_write(rawfh, f, datasz, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, mpierror)
                else
                        do k = istart(3),iend(3); do j = istart(2),iend(2); do i = istart(1),iend(1)
                                ii = (i-istart(1))/ds + 1; jj = (j - istart(2))/ds + 1; kk = (k - istart(3))/ds + 1
                                bov(ii, jj, kk) = f(i,j,k)                                
                        end do; end do; end do

                        call MPI_File_write(rawfh, bov, datasz, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, mpierror)
                endif

                call MPI_File_close(rawfh, mpierror)

                if (mpirank .eq. 0) then
                        ! BOV header
                        write (10,'(a,g25.16e2)') "TIME:        ", t
                        write (10,'(a,3g25.16e2)') "BRICK_ORIGIN:", 0.0, 0.0, 0.0
                        write (10,'(a,3g25.16e2)') "BRICK_SIZE:  ", n*dx, n*dx, n*dx
                        write (10,'(2a)') "VARIABLE:        ", v
                        write (10,'(2a)') "DATA_FILE:       ", buffer
                        write (10,'(a,3i5)') "DATA_SIZE:     ", n/ds, n/ds, n/ds
                        write (10,'(2a)') "DATA_FORMAT:     ", "DOUBLE"
                        write (10,'(2a)') "DATA_ENDIAN:     ", "LITTLE"
                        write (10,'(2a)') "CENTERING:       ", "zonal"

                        close (10); 
                end if
#endif
        end if
        
        ! output spectra and statistics
        if (output$crv) then
#ifndef SILO
                ! in VisIt format, all curves go into a single file per field, per frame
                if (output$vis .and. mpirank .eq. 0) then
                        write (buffer,'(a,a,i4.4,a)') v, '-', frame, '.ult'; open(12, file=buffer)
                end if
#endif
                
                ! in gnuplot format, all fields and frames go into a single file per curve
                ! (to be fflush()ed every frame after all the fields are analyzed)
                if (present(idx)) then; idx = idx + 1; DVAR(idx) = v; end if
                
                ! output power spectrum
                if (output$psd) then
                        call spectrum(f, S); if (present(idx)) PSD(:,idx) = S
                        
                        if (output$vis .and. mpirank .eq. 0) then
#ifdef SILO
                                e = dbputcurve(db, STR("PSD_"//v), (/0:ns-1/)*dk, S, DB_DOUBLE, ns, DB_F77NULL, e)
                                e = dbputcurve(db, STR("log10_PSD_"//v), log10((/1:ns-1/)*dk), log10(S(2:ns)), DB_DOUBLE, ns-1, DB_F77NULL, e)
#else
                                write (12,'(a)') "# PSD"
                                do k = 1,ns; write (12,'(2g25.16e2)') (k-1)*dk, S(k); end do
                                write (12,'(2a)') "", ""
                                
                                write (12,'(a)') "# PSD [logarithmic]"
                                do k = 2,ns; write (12,'(2g25.16e2)') log10((k-1)*dk), log10(S(k)); end do
                                write (12,'(2a)') "", ""
#endif
                        end if
                end if
                
                ! output distribution of values
                if (output$cdf) then
                        call distsortflat(f, tmp, tmp2, size(f,1)*size(f,2)*size(f,3))

                        P = 0.0;
                        if (istart(1) .eq. 1 .and. istart(2) .eq. 1) then
                                P(istart(3):iend(3)) = f(1,1,istart(3):iend(3))
                        end if
                        
                        if (mpirank .eq. (mpisize-1)) then
                                P(n+1) = f(iend(1),iend(2),iend(3))
                        end if
                       
                        call MPI_Allreduce(MPI_IN_PLACE, P, n+1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierror)
 
                        if (present(idx)) CDF(:,idx) = P
                        avg = sum(P)/(n+1); var = sum((P-avg)**2)/n; X = (P-avg)/sqrt(2.0*var)
                        do k = 2,n; PDF(k) = (2.0/n)/(P(k+1)-P(k-1)); end do
                        
                        if (output$vis .and. mpirank .eq. 0) then
#ifdef SILO
                                e = dbputcurve(db, STR("CDF_"//v), P, (/0:n/)/real(n), DB_DOUBLE, n+1, DB_F77NULL, e)
                                e = dbputcurve(db, STR("PDF_"//v), P(2:n), PDF, DB_DOUBLE, n-1, DB_F77NULL, e)
                                e = dbputcurve(db, STR("gaussian_CDF_"//v), P, (1.0 + erf(X))/2.0, DB_DOUBLE, n+1, DB_F77NULL, e)
                                e = dbputcurve(db, STR("gaussian_PDF_"//v), P(2:n), exp(-X(2:n)**2)/sqrt(twopi*var), DB_DOUBLE, n-1, DB_F77NULL, e)
#else
                                write (12,'(a)') "# CDF"
                                do k = 1,n+1; write (12,'(2g25.16e2)') P(k), real(k-1)/n; end do
                                write (12,'(2a)') "", ""
                                
                                write (12,'(a)') "# CDF [gaussian]"
                                do k = 1,n+1; write (12,'(2g25.16e2)') P(k), (1.0 + erf(X(k)))/2.0; end do
                                write (12,'(2a)') "", ""
                                
                                write (12,'(a)') "# PDF"
                                do k = 2,n; write (12,'(2g25.16e2)') P(k), PDF(k); end do
                                write (12,'(2a)') "", ""
                                
                                write (12,'(a)') "# PDF [gaussian]"
                                do k = 2,n; write (12,'(2g25.16e2)') P(k), exp(-X(k)**2)/sqrt(twopi*var); end do
                                write (12,'(2a)') "", ""
#endif
                        end if
                end if
                
#ifndef SILO
                if (output$vis .and. mpirank .eq. 0) close (12)
#endif
        end if
end subroutine dump

! flush frame data into gnuplot-style curves
subroutine fflush(t, idx)
        real t; integer idx, k; logical o

        if (mpirank .ne. 0) return        

        ! output power spectrum
        if (output$psd) then
                inquire (31, opened=o)
                
                if (.not. o) then
                        open(31, file="PSD"); call head(31, (/"t           ", "k           ", DVAR(1:idx)/))
                end if
                
                do k = 1,ns; write (31,'(32g25.16e2)') t, (k-1)*dk, log10(PSD(k,1:idx)); end do
                write (31,'(2a)') "", ""; flush(31)
        end if
        
        ! output distribution of values
        if (output$cdf) then
                inquire (32, opened=o)
                
                if (.not. o) then
                        open(32, file="CDF"); call head(32, (/"t           ", "percentile  ", DVAR(1:idx)/))
                end if
                
                do k = 1,n+1; write (32,'(32g25.16e2)') t, real(k-1)/n, CDF(k,1:idx); end do
                write (32,'(2a)') "", ""; flush(32)
        end if
end subroutine fflush

end
