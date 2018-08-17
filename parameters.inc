#ifndef N_VALUE
#define N_VALUE 128
#endif

#ifndef L_VALUE
#define L_VALUE 100.0
#endif

! solver control parameters
integer, parameter :: n = N_VALUE                   ! sampled grid size (simulation cube is n^3 pts)
integer, parameter :: p = n+2                   ! padded grid size (>n, adjust for cache efficiency)
integer, parameter :: tt = 2**24                ! total number of time steps to take (i.e. runtime)
integer, parameter :: nn = n/2+1                ! Nyquist frequency (calculated, leave it alone)
integer, parameter :: ns = sqrt3*(n/2) + 2      ! highest wavenumber on 3D grid (leave it alone)

real, parameter :: alpha = 10.0                 ! dx/dt (be careful not to violate Courant condition)
real, parameter :: len = L_VALUE
real, parameter :: dx = len/n                  ! grid spacing   (physical grid size is n*dx)
real, parameter :: dt = dx/alpha                ! time step size (simulated timespan is tt*dt)
real, parameter :: dk = twopi/(n*dx)            ! frequency domain grid spacing (leave it alone)
real, parameter :: amax = 5.0			! the maximum value of a
real, parameter :: initscale = 1.0 ! 3e-7       ! the initial mode scaling

!logical, parameter :: leinit = .true.           ! use LE-style ICs
!real, parameter :: len0 = 100.0                 ! the base of the L ICs scaling

! largest distance is sqrt(3(n-1)^2) = sqrt(3)*(n-1)
!integer, parameter :: dmax = int(sqrt(3.0)*(n-1))
!logical, parameter :: direct2pt = .false.
!logical, parameter :: fft2pt = .true.

real, parameter :: aflip = 3.75
logical, parameter :: fliplambda = .false.
logical, parameter :: flipg = .false.

! output control parameters
!integer, parameter :: nx = n                  ! spatial grid is downsampled to nx^3 pts for output
integer, parameter :: ds = 2			! the down-sample factor
integer, parameter :: nt = 2**14                 ! simulation will be logged every nt time steps

logical, parameter :: output = .true.           ! set this to false to disable all file output at once
logical, parameter :: oscale = .true.           ! scale output variables to counter-act expansion

logical, parameter :: output$bov = .true.       ! output 3D data cube (storage-expensive)
logical, parameter :: output$psd = .true.       ! output power spectra (time-expensive)
logical, parameter :: output$cdf = .false.      ! output distributions (time-expensive)
logical, parameter :: output$fvar = .true.      ! output the field variances
!logical, parameter :: output$twoptcorr = .true. ! output the 2-pt. correlation function

logical, parameter :: output$fld = .true.       ! output scalar fields
logical, parameter :: output$set = .true.       ! output stress-energy tensor components
logical, parameter :: output$pot = .false.       ! output gravitatinal potential (expensive)
logical, parameter :: output$any = output$fld .or. output$set .or. output$pot

logical, parameter :: output$gnu = .true.       ! output curves in gnuplot format (sinle file)
logical, parameter :: output$vis = .true.       ! output curves in VisIt X-Y format (per frame)
logical, parameter :: output$crv = (output$psd .or. output$cdf) .and. (output$gnu .or. output$vis)

!logical, parameter :: output$upframe = .false.    ! output a raw frame of the next slice also
!logical, parameter :: output$lhla    = .false.    ! write LH and LA to a file

! WARNING: Make sure that dt is the same as in the run which produced these files!
!logical, parameter      :: input$fromfile  = .false.      ! read the initial conditions from a file
!character(*), parameter :: input$hrphifile = "phi_hr.raw" ! the 2nd slice of phi
!character(*), parameter :: input$dnphifile = "phi_dn.raw" ! the 1st slice of psi
!character(*), parameter :: input$hrpsifile = "psi_hr.raw" ! the 2nd slice of psi
!character(*), parameter :: input$dnpsifile = "psi_dn.raw" ! the 1st slice of psi
!character(*), parameter :: input$lhlafile  = "lhla.raw"   ! the file with LH and LA

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! fields are referred to by their symbolic aliases
integer, parameter :: fields = 1                ! total number of scalar fields being evolved
integer, parameter :: phi = 1, psi = 2          ! symbolic aliases for scalar fields
integer, parameter :: rho = 1, prs = 2          ! symbolic aliases for stress-energy

! potential and its derivatives are (separately) inlined in step()
! model summary in human-readable form is printed out in head()
! parameters of the scalar fields potential are defined here

#ifndef LL_VALUE
#define LL_VALUE 1.25E-6
#endif

#ifndef GG_FACT
#define GG_FACT 5.0
#endif

real, parameter :: m2phi = 1.0, m2psi = 0.0, g2 = 0.0, mpl = 2e5, lambda = -0.0, ll=LL_VALUE, gg=sqrt(GG_FACT)*ll
real, parameter :: KillField = 0.0, f=0.001

! Needs to be changed accordingly.

! initial conditions for homogeneous field component
real, parameter ::  phi0 =  SQRT((m2phi/mpl**2)*(3.0/5.0)*(ll/gg**2))
real, parameter :: dphi0 =  -SQRT((1./2.0)*phi0**2-(ll/4.0)*(mpl**2/m2phi)*phi0**4+(mpl**4/m2phi**2)*(gg**2/6.0)*phi0**6)
real, parameter ::    H0 =   -(1.0/SQRT(2.0))*dphi0 
real, parameter ::   dH0 =   -H0**2
real, parameter :: ddphi0 =  -3.*H0*dphi0-(phi0-ll*(mpl**2/m2phi)*phi0**3+(gg**2)*(mpl**4/m2phi**2)*phi0**5)
real, parameter ::  ddH0 = -(1./2.)*dphi0**2-H0**3

! scale factor and horizon size (sampled on two subsequent time slices)
real :: LA(2) = (/ 1.0 - H0*dt, 1.0 /)
real :: LH(2) = 1.0/(/ H0 - dH0*dt + ddH0*dt**2/2.0, H0 /)

real :: lambdapf = 1.0
real :: gpf = 1.0

integer :: rssize
integer, allocatable :: rs(:)

