!
!
! --> BETA code starting anew <--
!
! -- data storge also changed in this version to coincide with FFT requirement
! -- this version uses CVODE instead VODPK or DVODPK
! -- language is updated to F90
! -- check for SILO.INC file format also (might still be in F77 format)
!    -- corrected to include file SILO_F90.inc instead
! -- FFTW 2.1.5 is used here
! -- this is a three-field version, density and potential fluctuations and a mean field
! -- precision is set by compiler flags forcing all REAL and DOUBLE to be 8-byte wide
! -- indices of density or potential arrays reflect the indices of corresponding modes,
!    this eliminates the guess work and indexing routines
! -- tracers, flights, lypunov tracking are reinstated
! -- file outputs for postprocessing are binary files (little-endian according to gcc 4.6.2)

MODULE precmod
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  MODULE precmod : specifies the precision of REAL types
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE

      INTEGER, PARAMETER :: IDP = KIND(1.d0)

END MODULE precmod

MODULE fft_variables
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  MODULE fft_variables: accounts for variables required for FFT
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE
!C
!C   First, the size of the box in Fourier space:
!C
!C     kxmax    = highest mode number in the x direction
!C     kymax    = highest mode number in the y-direction
!C
!C     kxm1  = total number of modes in the x-direction
!C     kym1  = total number of modes in the y-direction
!C
!C     kxm2  = max x-mode number + 1
!C     kym2  = max y-mode number + 1
!C
!C     nn    = total number of modes to be solved (excluding the (0,0) mode)
!C
!C     writcount is a dummy variable
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCC                                                                  CCC
!CCC    !!!!!!! W A R N I N G !!!!!!!                                 CCC
!CCC                                                                  CCC
!CCC    KXMAX and KYMAX must be multiples of 2, 3, and 5              CCC
!CCC                                                                  CCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      INTEGER, PARAMETER :: &
!           kxmax = 10, kymax = 10, &
           kxmax = 40, kymax = 40, &
!           kxmax = 128, kymax = 128, &
!           kxmax = 220, kymax = 220, &
           kxm1 = 2*kxmax+1, kym1 = 2*kymax+1, &
           kxm2 = kxmax + 1, kym2 = kymax + 1, &
           nn = 6*(kym1*kxmax + kym2)                    ! two-fields (REAL and IMAG) + mean field modes (REAL and IMAG), a total of three fields
!           nn = 6*kym1*kxm2

      INTEGER, PARAMETER :: &                            ! parameters for FFT - accounting for aliasing (require 3*kxmax + 1 points)
!           nxreal = 128, nyreal = 128, &                 ! for a maximum of 20x20 modes
           nxreal = 256, nyreal = 256, &                 ! number of points in REAL-space, 256 for a maximum of 40x40 modes
!           nxreal = 343, nyreal = 343, &
!           nxreal = 1024, nyreal = 1024, &                 ! for a maximum of 170x170 modes
           nxcmplx = nxreal/2 + 1, nycmplx = nyreal      ! number of points in k-space (this is twice the one-field version)

!C constant definitions needed by FFTW ------------------------------------------------
      INTEGER, PARAMETER :: FFTW_FORWARD = 1,FFTW_BACKWARD = -1
      INTEGER, PARAMETER :: FFTW_REAL_TO_COMPLEX = -1,FFTW_COMPLEX_TO_REAL = 1
      INTEGER, PARAMETER :: FFTW_ESTIMATE = 0,FFTW_MEASURE = 1
      INTEGER, PARAMETER :: FFTW_OUT_OF_PLACE = 0,FFTW_IN_PLACE = 8,FFTW_USE_WISDOM = 16
      INTEGER, PARAMETER ::  FFTW_THREADSAFE = 128
      INTEGER, PARAMETER :: FFTW_TRANSPOSED_ORDER = 1,FFTW_NORMAL_ORDER = 0
!C ------------------------------------------------------------------------------------

      INTEGER(kind=8) :: plan1, plan2
      INTEGER(kind=8) :: plan1_1d, plan2_1d
      
END MODULE fft_variables


MODULE variables
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  MODULE variables: accounts for main variables in code
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE FFT_VARIABLES
      IMPLICIT NONE

!C
!C      Constants:
!C        -- values extracted from former two field code
!C       epsi    = fraction of trapped electrons
!C       rhos    = ion gyration radius
!C       emass   = electron mass
!C       nue     = effective collision frequency
!C       Lnc     = critical gradient

      REAL, PARAMETER :: pi = 4.d0*ATAN(1.d0) !3.14159265358979323846264338328d0
      REAL, PARAMETER :: twopi = 2.d0*pi
      COMPLEX, PARAMETER :: ic = CMPLX(0.d0,1.d0), czero = CMPLX(0.d0,0.d0)
!      REAL, PARAMETER:: epsi = 0.3 !0.001
!      REAL, PARAMETER:: rhos = 0.1 !0.01      ! this helps with time step convergence, controls scale length
!      REAL, PARAMETER:: emass = 1.0, nue = 5.0
      REAL :: epsi, rhos, emass, nue, rtepsi, rhos3

!      REAL, PARAMETER:: rtepsi = SQRT(epsi)
!      REAL, PARAMETER:: rhos3 = rhos*rhos*rhos
!      REAL, PARAMETER:: Lnc = 1.0E-1
 
!C	kxnorm	= normalization for the kx modes
!C	kynorm	= normalization for the ky modes
!C	ekkold	= an array to store old values of the amplitudes
!C	v1	= beta
!C	gama1	\                             -> low k
!C	gama2	 = k-dependent damping rates  -> mid k
!C	gama3	/                             -> high k
!C	amp	= used to define tola 
!C	ek	= the total (time integrated) energy
!C	en	= the total (time integrated) enstrphy
!C	ekold	= old value of the energy
!C	enold	= old value of the enstrophy
!C	t	= a time counter (REAL type)
!C	dt	= time step size
!C	tol	= tolerance measure
!C	gamk	= rate of growth of the total energy
!C	gamen	= rate of growth of the total enstrophy
!C	ratio	= how peaked an initial profile
!C	psi0	= for vortex and peaked profile initialization
!C	din     = \                            -> low k            (old variable is "d")
!C	dmid	=  = k-dependent driving terms -> mid k
!C	dout	= /                            -> high k
!C	flreal	= real part of flat profile initialization
!C	flima	= imag part of flat profile initialization
 
      REAL :: kxnorm,kynorm,v1,gama1,gama2,gama3,amp
!      REAL, DIMENSION(0:kxmax,-kymax:kymax):: ekkold
      REAL :: denen,dengamen,denekold,denenold,denek,dengamek                         ! density values
      REAL :: phien,phigamen,phiekold,phiek1old,phienold,phiek,phiek1,phigamek        ! phi values
      REAL :: denavgek, denavgekold                                                   ! mean field values
      REAL :: ek, ekold                                                               ! global energy values
      REAL :: t,dt,tol,phi0flow,ratio, psi0, din, dmid, dout, flREAL, flima           ! DIN=DMID=DOUT=LE if original DTEM equation is to be integrated.
      REAL :: dcoeff                                                                  ! variable to model old set of equations
      REAL :: dmeanf                                                                  ! to model mean field diffusion at high-k modes

!C	ncl	= time counter
!C      ncl_in  = previous time counter read-in from RESTART file

      INTEGER :: ncl, ncl_in
 
!C      ku,kl	= upper and lower cutoffs for flat band 
!C	kg1sq	= lower cut-off for k^2-dependent damping
!C   	kg3sq	= upper cut-off for k^2-dependent damping
!C	kd1sq	= lower cut-off for k^2-dependent drive
!C	kd3sq	= upper cut-off for k^2-dependent drive
!C	kinit	= lower cut-off for initialization
!C      kdmeanf = number of modes to apply mean field diffusion

      INTEGER :: ku,kl,kg1sq,kg3sq,kd1sq,kd3sq,kinit,kdmeanf
    
!C	den 	= the array which contains the stream function in k-space.
!C      phi     = array contains the potential function in k-space
!C      denavg  = density mean field (x = radial direction, y = poloidal direction), d<n>/dy = 0 only ky = 0 modes are valid
!C      phi0    = this field is for an external field

      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax) :: den, phi, denavg   ! these are the principal quantities!!!
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax) :: phi0

!C      phi0amp = amplitude of an external flow
      REAL :: phi0amp

!C      rLn       = length scale stored to write to file
!C      rdenavg1d = mean field feedback on the fluctuating fields
      
      REAL, DIMENSION(nxreal) :: rLn, rdenavg1d

!C      philinarr1  = array of terms multiplying phi in phi equation
!C      philinarr2  = array of terms multiplying den in phi equation
!C      philinarr3  = array of terms multiplying non-linear term in phi equation
!C      denlinarr   = array of terms miltiplying non-linear term in den equation
!C      densource   = external source to define external gradient
!C      rdensource  = source term in REAL-space

!      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax) :: philinarr1, philinarr2, philinarr3, philinarr4
!      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax) :: denlinarr, densource
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax) :: densource
      REAL, DIMENSION(nxreal,nyreal) :: rdensource

      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax) :: philinarr, damparr, dampdenavgarr, drivearr
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax) :: philengthlinarr
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax) :: linkdrive        ! stored array for output of the drive for fluctuation fields

      ! extras
!      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax) :: drivedenavgarr   ! temporary array to drive mean field at a certain mode
!      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax) :: linarr
!      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax) :: gradnarr
!C      gradnarr    = scale length of density gradient, using the average values to define scale length

!C	ntstep	= number of time steps
!C	irand	= number to feed random number generator

      INTEGER :: ntstep
      INTEGER(kind=8):: irand
 
!C       lp      = dummy constant (lp is a control on Pol nonlinear terms) (pol drift NL or NS)
!C       le      = dummy constant (le is a control on ExB nonlinear terms) (ExB NL)
!C       lm      = dummy constant (lm is a control on mean field convection)
      REAL :: lp,le,lm

!C    Theses terms pertain to the source and the mean field.
!C       Lnc     = critical gradient
!C       samp    = source amplitude
!C       sscale  = source horizontal scale
      REAL :: Lnc, samp, sscale

!C       tola the tolerance vector passed to VODPK matches element wise 
!C       with the Y(nn) vector

      REAL, DIMENSION(nn) :: tola

!C       nspecmodesx = number of modes in the kx direction
!C       nspecmodesy = number of modes in the ky direction
!C       specout     = storage array to output energy spectrum at specific modes
!C       specmodesx  = specific kx modes for output
!C       specmodesy  = specific ky modes for output
      INTEGER, PARAMETER :: nspecmodesx = 3, nspecmodesy = 3       ! these values are currently hard wired but can be modified through some sort of input
!      COMPLEX, DIMENSION(3,nspecmodesx,nspecmodesy) :: specout  ! currently set to 9 modes per field (9 modes * 3 fields = 27 modes)
      INTEGER, DIMENSION(nspecmodesx), PARAMETER :: specmodesx = (/ 1, INT(kxmax/2), kxmax /)
      INTEGER, DIMENSION(nspecmodesy), PARAMETER :: specmodesy = (/ 0, INT(kymax/2), kymax /)

END MODULE variables


MODULE prog_flags
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  MODULE prog_flags: accounts for all the code flags
!C     (e.g. fft routine, save big files, etc.)
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE

!C	lpout	= turn on(=1) or off(=0) nonlinearity at high k's
!C	lpmid	= turn on(=1) or off(=0) nonlinearity at mid k's
      INTEGER :: lpout,lpmid

!C	init	= flag for which initial conditions to use
      INTEGER :: init

!C	nwrt	= how often to write ENERGY files
!C	nwrs	= how often to write RESTART files for spectrum and particles
!C       nvisit  = how often to write to silo files - SD
!C       nfft    = which fft subroutine to use. fftw if nfft=1, matfft otherwise (matfft has been completely disabled!)
      INTEGER :: nwrt,nwrs,nvisit,nfft

!C	To save or not to save the BIG files: stream(3), vorticity(34), 
!C       ensemble(12), energ_vs_ak(14), energyspec(11)
      LOGICAL :: l_save

!C       l_firstrun = construct linear driving and damping terms (not needed)
!      LOGICAL :: l_firstrun

      CHARACTER(LEN=5):: BETA_version = 'x.x0'

      LOGICAL :: l_tracers, l_tracers_renew                              ! L_TRACERS = .T. uses tracers, L_TRACERS_RENEW = reinit tracers every NPT iter
      LOGICAL :: l_tracers_visit     ! output tracer trajectories to visit file 

!C       l_omp = enables OpenMP capability
      LOGICAL :: l_omp

!C       drivedamp = select model for drive and damp construction
      INTEGER :: drivedamp

!C       lcont = checks if the simulation is a continuation
      LOGICAL :: lcont

!C       l_specout = spectral output for specific modes
      LOGICAL :: l_specout

!C       l_extflow = flag for external flow
      LOGICAL :: l_extflow

!C       l_nltrans = flag for nonlinear transfer output of specific modes
      LOGICAL :: l_nltrans

END MODULE prog_flags


MODULE cvode_variables
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  MODULE cvode_variables: accounts for CVODE construction and solving
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE

      INTEGER(kind=8), DIMENSION(21):: iout
      REAL, DIMENSION(6):: rout
      INTEGER :: ipar
      REAL :: rpar

END MODULE cvode_variables


MODULE visit
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  MODULE visit: accounts for ViSit outputs
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      IMPLICIT NONE

      INCLUDE "silo_f90.inc"

      INTEGER :: nx_visit, ny_visit, nkx_visit, nky_visit, ivisit
      INTEGER :: ndims 
      INTEGER, DIMENSION(:), ALLOCATABLE:: xdims, kdims
      REAL, DIMENSION(:), ALLOCATABLE:: x_visit, y_visit, kx_visit, ky_visit


END MODULE visit


MODULE fileunits
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  MODULE fileunits: define external fileunits
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE

      INTEGER, PARAMETER :: &
           u_denpolavg = 585, u_phipolavg = 586, u_denavgpolavg = 587, &   ! files for poloidal averages
           u_Ln = 588, u_denavg1d = 589 ! inverse length scale, feedback field on fulctuating fields
      
      INTEGER, PARAMETER :: &
           u_energy = 10 ! file for energy

      INTEGER, PARAMETER :: &
           u_tracers = 645, &  ! binary file containing positions and velocities of all tracers
!           u_tposx = 645, u_tposy = 646, u_tvelx = 647, u_tvely = 648, &    ! files for tracers
           u_tlypmax = 649  ! file for maximum lypunov exponent

      INTEGER, PARAMETER :: &
           u_tflightsx = 731, u_tflightsy = 732  ! files for tracer flights

      INTEGER, PARAMETER :: &
           u_inputscopy = 8, & ! file for copy of input values
           u_specout = 432, &  ! file for spectrum output
           u_nltrans = 433     ! file for nonlinear transfer ouput

END MODULE fileunits



MODULE tracers
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  Particle initialization parameters
!C     mbox= number of marker particles
!C     nbox= number of pair particles
!C
!C     Total number of particles = mbox**2*nbox**2
!C
!C     nshift  = meridional shift of particle initial positions
!C     nvpart  = type of particle initialization
!C     vsize   = size of circle for nvpart.eq.3 initialization
!C     pspread = max size of delta-y_rms before reinitializing
!C     ncount  = reinitialization counter for tracers                  ! Moved from VARIABLES here.
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod

      IMPLICIT NONE

      INTEGER:: mbox, nbox, nvpart
      REAL:: vsize, pspread, nshift
      INTEGER:: nmax, npoint, npair, ixmax, jymax
      INTEGER :: lypcounter, itrack, iflightsx, iflightsy                         ! counters, ** "itrack" counts every time iteration minus the first iteration **
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: tposarr_old                         ! old tracer positions for VISIT output
      INTEGER, DIMENSION(:), ALLOCATABLE :: ixstart, jystart
      REAL, DIMENSION(:,:), ALLOCATABLE :: px0, py0, px, pvx0, pvy0, py, pvx, pvy ! position and velocities for pairs
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: px_tail, py_tail                     ! stores some histories of the positions
      REAL, DIMENSION(:), ALLOCATABLE :: lypmax, pdistinit                        ! initial tracers' distance
      REAL, DIMENSION(:), ALLOCATABLE :: deltax_old, deltay_old, &                ! variables for flights
           timex_old, timey_old, xpos_old, ypos_old
      REAL :: told                                                                ! stores time of previous step, used to advance tracers
      REAL :: smalldist                                                           ! this value is specified through the input file
      INTEGER :: lypmaxdtstep                                                     ! maximum lypunov exponent time step
      REAL, PARAMETER :: boxlen = 1                                               ! box length
!      REAL, DIMENSION(:,:), ALLOCATABLE:: vx, vy
!      COMPLEX, DIMENSION(:), ALLOCATABLE:: pt, vt, vtfix, pt0 
      INTEGER:: npos                                                    ! frequency with which tracer positions written on file = NPOS
      INTEGER:: ncount                                                  ! Moved here from VARIABLES. Set in INITALIZE_TRACERS.

END MODULE tracers

!!$MODULE convolve_module
!!$!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$!C
!!$!C  MODULE convolve_module: declares pointer values to use
!!$!C                          in CONVOLVE subroutine
!!$!C
!!$!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$      USE fft_variables
!!$      IMPLICIT NONE
!!$
!!$!C  To compensate for de-aliasing, dimensions are doubled for kx and ky dimensions.
!!$!C  although "2/3" rule applies to compensate for de-aliasing, FFTW transforms efficiently for multiples of "2";
!!$!C  hence, this is motivation for doubling the array size. Need to only transform half on modes, since the
!!$!C  others are just conjugates. The modes must be placed with zero modes on the boundaries.
!!$!C
!!$
!!$!C  Declare FFT data.
!!$ 
!!$      COMPLEX, DIMENSION(nxcmplx,nycmplx) :: &
!!$           kden, &               !stream function
!!$           kphi, &               !potential function
!!$           cnkx,cnky, &          !components of gradients of dn/dy
!!$           cphiwx, cphiwy, &     !components of the lapacian of potential
!!$           cphix, cphiy, &       !components of the gradient of potential
!!$           cPol,cExB             !Polarization and ExB non-linearity
!!$
!!$      REAL, DIMENSION(nxreal,nyreal) :: &
!!$           rden,rphi,nkx,nky,phiwx,phiwy,phix,phiy,Pol,ExB
!!$
!!$!C  The components of the k-vector.
!!$
!!$      REAL :: kx_conv(nxcmplx), ky_conv(nycmplx)
!!$
!!$END MODULE convolve_module


MODULE omp_variables
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  MODULE omp_variables: common variables for OpenMP
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE

      INTEGER :: nthreads,thread_id
!      INTEGER :: nthreads_req ! number of requested threads
!      INTEGER :: OMP_GET_NUM_PROCS,OMP_GET_THREAD_NUM,OMP_GET_MAX_THREADS
! ** cannot declare OpenMP functions in a module **

END MODULE omp_variables


!C BEGIN PROGRAM SECTION ********************************************************************************************************************************************************************
PROGRAM BETA
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  PROGRAM BETA: main program
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE prog_flags
      USE variables, ONLY:ntstep,ncl,ncl_in
      USE fileunits
      IMPLICIT NONE

      INTEGER :: i,j
      REAL :: start_time, end_time
      REAL :: time_init, time_fin

!C  Banner variables
      CHARACTER(LEN=50), PARAMETER:: banner = ' THIS IS BETA . VERSION: '
      CHARACTER(LEN=30):: computer
      INTEGER:: imon
      CHARACTER(LEN=10):: date0, time0, zone0
      CHARACTER(LEN=3), DIMENSION(12), SAVE:: months
      CHARACTER(LEN=40):: dateloc
      DATA months/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/



!C >>>>>> Disable this before submit to PACMAN! <<<<<<<<<<<<<<<
!      WRITE(*,*) "Hit enter to start BETA!"
!      READ(*,*)
!C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      computer = ' Darwin OSX'

!C  Print banner -------------------------------------------------------
      WRITE(6,*) ' ========================================================================================='
      CALL DATE_AND_TIME(date0,time0,zone0)
      READ(date0(5:6),'(i2)') imon
      WRITE(dateloc,100) months(imon),date0(7:8),date0(1:4),time0(1:2),time0(3:4),time0(5:6)
100   FORMAT(' DATE = ',a3,' ',a2,',',a4,' ',' TIME= ',2(a2,':'),a2)
      WRITE(6,'(1x,2a,/,1x,2a)') banner, BETA_version, computer, dateloc
      WRITE(6,*) " This is a two-field version."
      WRITE(6,*) ' ========================================================================================='
      WRITE(6,*)
!C  END printing banner ------------------------------------------------



!C  Initialize...
      CALL initial
      CALL INIT_VISIT

!      CALL energy        ! calculate initial energies
      IF (lcont .EQV. .FALSE.) THEN    ! if this is NOT a continuation run (or first run)
         WRITE(u_inputscopy,'(1x,2a,/,1x,2a)') banner, BETA_version, computer, dateloc
         CALL writerestart                                    ! write initial "restart" files
         CALL energy                                          ! calculate initial energies
         CALL writeit(u_energy)                               ! write initial energy to file
         CALL write_visit                                     ! write intial output
         CALL write_polavg(u_denpolavg,u_phipolavg,u_denavgpolavg,u_Ln,u_denavg1d)
         IF (l_specout .EQV. .TRUE.) CALL write_specout       ! write output for specific modes
         IF (l_nltrans .EQV. .TRUE.) CALL write_nltrans       ! write nonlinear transfer for specific modes
      END IF

!C  Iterate over time steps... also shift time step if this is a continuation run
      DO ncl = ncl_in+1, ntstep+ncl_in

         CALL timead

        ! write output to ViSit
         IF ( MOD(ncl,nvisit) == 0 ) THEN
            CALL write_visit
            CALL energy               ! calculates energy
            CALL writeit(u_energy)    ! write to energy file
            CALL writerestart         ! write "restart" file
            CALL write_polavg(u_denpolavg,u_phipolavg,u_denavgpolavg,u_Ln,u_denavg1d)
            IF (l_nltrans .EQV. .TRUE.) CALL write_nltrans    ! write nonlinear transfer for specific modes
         END IF

         IF ( (MOD(ncl,INT(nvisit/10)) == 0) .AND. (l_specout .EQV. .TRUE.) ) CALL write_specout    ! write output for specific modes at 1/10 the frequency of visit output

      END DO
      ncl = ncl - 1  ! counter is automatically advanced at end of DO loop

!C  Write final "restart" file
!      CALL energy
!      CALL writeit(u_energy)
      CALL writeinit
!      CALL write_polavg(u_denpolavg,u_phipolavg,u_denavgpolavg,u_Ln,u_denavg1d)
      CALL make_namelist

!C  Terminate all files
      CALL FINALIZE
      CALL que_hora_es(6)

END PROGRAM BETA
!C END PROGRAM SECTION ********************************************************************************************************************************************************************


SUBROUTINE initial
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                                                                      C
!C   subroutine initial                                                 C
!C                                                                      C
!C   Calls the initializing subroutines                                 C
!C   and writes the initail energy spectrum.                            C
!C                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE prog_flags
      USE variables
      USE fft_variables
      USE fileunits
      USE omp_variables
!      USE convolve_module
      IMPLICIT NONE

      INTEGER :: i, j, itime
      REAL :: start_time, end_time    !timing variables


!C  OpenMP functions
      INTEGER :: OMP_GET_NUM_PROCS,OMP_GET_THREAD_NUM,OMP_GET_MAX_THREADS

      WRITE(6,*) "Precision of REAL type is: ", PRECISION(pi)  ! print out machine precision in decimal places
      WRITE(6,*) "Range of REAL type is: ", RANGE(pi)
      WRITE(6,*) "Maximum exponent of REAL type is: ", MAXEXPONENT(pi)
      WRITE(6,*) "Minimum exponent of REAL type is: ", MINEXPONENT(pi)

      CALL check_cont(lcont)   !!! check for continuation

      ! set random variable
      CALL SYSTEM_CLOCK(itime)
      irand = -INT(itime)
      IF (lcont .EQV. .FALSE.) WRITE(u_inputscopy,*) 'IRAND=', irand

      CALL readit              !!! Read in parameters

      ! define constants from read-in parameters
      rtepsi = SQRT(epsi); rhos3 = rhos*rhos*rhos
!      dcoeff = 1./rtepsi ! defined a particular case

      CALL SetToleranceVector                         !!! Set the tolerance vector
      CALL makesource                                 !!! Precompute source terms
      CALL pert                                       !!! Set the inital perturbation
!      CALL energy                                    !!! Calculates the initial energies
      CALL solverinit                                 !!! Starts the solver routines
      CALL linterm_construct                          !!! Construct linear terms

      ! initialize optional outputs
      IF (l_tracers .EQV. .TRUE.) CALL tracers_init   !!! Initialize tracer parameters
      IF (l_specout .EQV. .TRUE.) CALL specoutinit    !!! Specify the modes for output
      IF (l_nltrans .EQV. .TRUE.) CALL nltransinit    !!! Opens file for nonlinear transfer data

!      denavg = denavg + densource   !sets initial vallue to source + initial flat shifted profile

!!$!C  Initialize the "k-vectors" -- these are used for the CONVOLVE routine
!!$
!!$      kx_conv(1:kxm2) = kxnorm*(/ (i,i=0,kxmax) /)
!!$      ky_conv(1:kym2) = kynorm*(/ (j,j=0,kymax) /)
!!$      ky_conv(nycmplx-kymax+1:nycmplx) = kynorm*(/ (j,j=-kymax,-1) /)


!C  Get OpenMP threads
!$      nthreads = OMP_GET_MAX_THREADS()
!$      IF (nthreads .LE. 0) nthreads = OMP_GET_NUM_PROCS()
!$      WRITE(6,*) "Number of threads to use is: ", nthreads
!$      WRITE(6,*) "Number of processors is: ", OMP_GET_NUM_PROCS()

!C   Starting FFTW initiallization -------------------------------------
      
      CALL CPU_TIME(start_time) !start timer
      CALL fftwinit
      CALL CPU_TIME(end_time)                                      !stop timer
      WRITE(6,*) "Time to create plans is:", end_time - start_time
!C   End FFTW initiallization ------------------------------------------


END SUBROUTINE initial


SUBROUTINE solverinit
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                                                                      C
!C   subroutine: solverinit                                             C
!C                                                                      C
!C   Initialize solver routines (CVODE)                                 C
!C                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE fft_variables
      USE variables
      USE cvode_variables

      IMPLICIT NONE

      INTEGER :: ierr
      REAL, DIMENSION(nn) :: y

!  ** WRITE statements give allocation errors in CVODE ** (issue solved, use INT8 instead!)
!  ** INT8 is important not to get memory errors. **

      ipar = 0; rpar = 0.0  ! solver variables

!C   Convert den to y
      CALL convert(den,phi,denavg,y)

!C  Begin CVODE initialization -----------------------------------------
!C  Initialize CVODE solver (serial version)
!C     call FNVINITS(KEY=1,NEQ=nn,IER=ierr)
      CALL FNVINITS(1,INT(nn,kind=8),ierr)
      IF (IERR .NE. 0) THEN
         WRITE(6,20) IERR
20       FORMAT(///' SUNDIALS_ERROR: FNVINITS returned IER = ', I5)
         CALL FINALIZE
         STOP
      ENDIF

!C   Allocate memory before solving
!C     CALL FCVMALLOC(T0=t, Y0=y, METH=2, ITMETH=1, IATOL=2, RTOL=tol,
!C    & ATOL=tola, IOUT=iout, ROUT=rout, IPAR=ipar, RPAR=rpar, IER=ierr)
      CALL FCVMALLOC(t, y, 2, 1, 2, tol, tola, iout, rout, ipar, rpar, ierr)
      IF (IERR .NE. 0) THEN
         WRITE(6,30) IERR
30       FORMAT(///' SUNDIALS_ERROR: FCVMALLOC returned IER = ', I5)
         CALL FINALIZE
         STOP
      ENDIF

!      WRITE(*,*) "MALLAOC Success!"

!   Use diagonal approximation for Jacobian
!      CALL FCVDIAG(ierr)

!   Trying another method since convergence was very pure for diagonal approximation
!   CALL FCVSPGMR(IPRETYPE, IGSTYPE, MAXL, DELT, IER)
!     IPRETYPE = preconditioner, 0 for no preconditioner
!     IGSTYPE  = Gram-Schimdt process, 1 for modified, 2 for classic
!     MAXL     = maximum Krylov subspace dimension
!     DELT     = linear convergence, 0.0 for default
      CALL FCVSPGMR(0,2,5,0.0,IERR)
      IF (IERR .NE. 0) THEN
         WRITE(6,40) IERR
40       FORMAT(///' SUNDIALS_ERROR: FCVSPGMR returned IER = ', I5)
         CALL FINALIZE
         STOP
      ENDIF

!   Set maximum internal steps
!      CALL FCVSETIIN("MAX_NSTEPS", 500, IERR)
!      IF (IERR .NE. 0) THEN
!         WRITE(6,40) IERR
!50       FORMAT(///' SUNDIALS_ERROR: FCVSETIIN returned IER = ', I5)
!         STOP
!      ENDIF

!C   End CVODE initialization ------------------------------------------
END SUBROUTINE solverinit


SUBROUTINE fftwinit
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                                                                      C
!C   subroutine: fftwinit                                               C
!C                                                                      C
!C   Initialize FFTW routines                                           C
!C                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE fft_variables
      USE prog_flags, ONLY: l_omp
      IMPLICIT NONE
      INTEGER :: ierr
      INTEGER(kind = 8) :: fftwomp_intflag

      WRITE(6,*) "creating plan..."

!C  Use FFTW_ESTIMATE for quick plan creation ----------------

      fftwomp_intflag = 0
!$      IF (l_omp .EQV. .TRUE.) THEN  ! switch FFTW routines
!$         fftwomp_intflag = FFTW_THREADSAFE
!$      WRITE(6,'(/a/)') "*** This message prints when OpenMP is enabled! ***"
!$      CALL fftw_f77_threads_init(ierr)
!$      ELSE
         WRITE(6,'(/a/)') "*** This is a serial run. ***"
!$      END IF
      CALL rfftw2d_f77_create_plan(plan2, nxreal, nyreal, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE + fftwomp_intflag)   ! used for two-fields
      CALL rfftw2d_f77_create_plan(plan1, nxreal, nyreal, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE + fftwomp_intflag)
      CALL rfftw2d_f77_create_plan(plan2_1d, nxreal, 1 , FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE + fftwomp_intflag)     ! using 2D form for 1D transforms, this is for the mean field
      CALL rfftw2d_f77_create_plan(plan1_1d, nxreal, 1, FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE + fftwomp_intflag)
!      CALL rfftw2d_f77_create_plan(plan2_1d,2*kxmax,1,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE)     ! using 2D form for 1D transforms
!      CALL rfftw2d_f77_create_plan(plan1_1d,2*kxmax,1,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE)

      WRITE(6,*) "plan1: ", plan1
      WRITE(6,*) "plan2: ", plan2

!C  ---------------------------------------------------------
END SUBROUTINE fftwinit

SUBROUTINE readit
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C     subroutine READIT
!C
!C     Reads in the parameter values.
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE variables
      USE tracers
      USE prog_flags
      USE fileunits
      IMPLICIT NONE
 
!C   Namelist input.
      NAMELIST /runp/ dt,ntstep,nwrt,nwrs,nvisit,tol,kxnorm,kynorm
      NAMELIST /physp/ v1,rhos,nue,emass,epsi,gama1,gama2,gama3,din,dout,dmid,dcoeff,dmeanf
      NAMELIST /physp2/ lp,le,lm,Lnc,lpmid,lpout
      NAMELIST /intilp/ amp,ku,kl,kg1sq,kg3sq,ratio
      NAMELIST /linzones/ kg1sq,kg3sq,kd1sq,kd3sq,kdmeanf,kinit
      NAMELIST /pertin/ init,psi0,flreal,flima, samp, sscale, phi0amp
      NAMELIST /parts/ mbox,nbox,nshift,nvpart,vsize,pspread,npos,smalldist,lypmaxdtstep
      NAMELIST /flags/ drivedamp,l_omp,l_save,l_tracers,l_tracers_renew,l_tracers_visit,l_specout, l_nltrans, l_extflow

      OPEN(4,file='inputs',status='old')

      READ (4,NML = runp)
      READ (4,NML = physp)
      READ (4,NML = physp2)
      READ (4,NML = intilp)
      READ (4,NML = linzones)
      READ (4,NML = pertin)
      READ (4,NML = parts)
      READ (4,NML = flags)

      CLOSE(4)

      ! write initial inputs to file for reference
      IF (lcont .EQV. .FALSE.) THEN
         WRITE (u_inputscopy,*)'xsize by ysize ',kxm1,' by ',kym1
         WRITE (u_inputscopy,runp)
         WRITE (u_inputscopy,physp)
         WRITE (u_inputscopy,physp2)
         WRITE (u_inputscopy,intilp)
         WRITE (u_inputscopy,linzones)
         WRITE (u_inputscopy,pertin)
         WRITE (u_inputscopy,parts)
         WRITE (u_inputscopy,flags)
      END IF

END SUBROUTINE readit



SUBROUTINE SetToleranceVector
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C     subroutine SetToleranceVector
!C  	  This routine will set the values of a tolerance vector 
!C     to be passed to VODPK. The layout follow the format of the
!C     Y(nn) vector, where each entry is inversly proportional
!C     to the wavenumber squared  (i.e. ATOL(i) = TOL/(kx^2+ky^2))
!C	  where k^2 = (kx^2+ky^2).  This yields tighter absolute tolerance
!C     in higher modes.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE variables, ONLY: tol,tola
      USE fft_variables, ONLY: kxmax,kymax,kym1,kym2,nn
      IMPLICIT NONE

      INTEGER :: i
      INTEGER :: kx, ky, kount, kountm, kountd, kountdm, kounts

      kount = kym2
!      kount = 0
!      kountm = kym1*kxmax + kym2
      kountm = nn/(2*3)
      kountd = 2*kountm
      kountdm = 3*kountm
      kounts = 4*kountm

!C  Do for kx = 0, -kymax <= ky <= 0
      tola(kym2) = tol ! (0,0) mode
      tola(1:kymax) = (/ (tol/REAL(ky*ky),ky=-kymax,-1) /)
      tola(kountm+1:kountm+kym2) = tola(1:kym2)

!C  Do for the rest
      DO ky = -kymax, kymax
         DO kx = 1, kxmax
            kount = kount + 1
            tola(kount) = tol/REAL(kx*kx+ky*ky)
!            tola(kountm+kount) = tola(kount)
         END DO
      END DO

!!$      DO ky = -kymax, 0
!!$         DO kx = 0, kxmax
!!$            kount = kount + 1
!!$            IF ((kx .EQ. 0) .AND. (ky .EQ. 0)) THEN
!!$               tola(kount) = tol
!!$            ELSE
!!$               tola(kount) = tol/REAL(kx*kx+ky*ky)
!!$            END IF
!!$         END DO
!!$      END DO
!!$
!!$      DO ky = 1, kymax
!!$         DO kx = 1, kxmax
!!$            kount = kount + 1
!!$            IF ((kx .EQ. 0) .AND. (ky .EQ. 0)) THEN
!!$               tola(kount) = tol
!!$            ELSE
!!$               tola(kount) = tol/REAL(kx*kx+ky*ky)
!!$            END IF
!!$         END DO
!!$      END DO

      tola(kountm+1:kountd) = tola(1:kountm)  !apply to imaginary part
      tola(kountd+1:kounts) = tola(1:kountd)  !apply the same to potential tolerances
      tola(kounts+1:nn)     = tola(1:kountd)  !apply the same to mean field tolerances

      OPEN(45,file="tola_dump",status="replace")
      WRITE(45,*) (tola(i),i=1,kountm)
      CLOSE(45)

END SUBROUTINE SetToleranceVector


SUBROUTINE pert
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  subroutine pert                                              12/16/94
!C
!C  Sets the initial perturbation.
!C 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE fft_variables
      USE variables
      USE prog_flags, ONLY: init, lcont
      USE fileunits
      IMPLICIT NONE
 
      REAL, EXTERNAL :: RAN2

      LOGICAL :: temp_logic, lexist
      COMPLEX :: temp_cmplx
      COMPLEX, DIMENSION(-kymax:kymax) :: temp_cmplx_vect
      INTEGER :: i, j
      REAL :: xk, yk, ak, a1, a2, a3
      REAL :: phase, small, zmodeavg
      CHARACTER(LEN=256) :: restartfile_str

!C   Initialization.

      small = TINY(0.0d0)
      den = CMPLX(0.0,0.0); phi = CMPLX(0.0,0.0); denavg = CMPLX(0.0,0.0)
      zmodeavg = 100.0d0      ! value for average of mean field

!C   Set the initial stream function perturbation according to INIT.
!C
!C   Value of INIT   Perturbation
!C   -------------   ------------
!C   0               Flat energy band initial profile (from kl to ku)
!C                   with random phases.
!C
!C   1               Completely flat amplitude profile (energy ~ k^2).
!C
!C   2               Peaked initial profile with random phases.
!C
!C   3               Peaked initial profile with random phases,
!C                   no energy in the k < kinit modes.
!C
!C   4               Read in amplitudes from file restart.bin file.
!C
!C   5               Peaked initial profile with random phases,
!C                   plus a jet or vortex read from file initspec.
!C
!C   6               Peaked initial profile with random phases,
!C                   no energy in the k < kinit modes, 
!C                   plus a jet or vortex read from file initspec.
!C
!C
!C   42              Zero perturbations
!C

      IF (lcont .EQV. .FALSE.) THEN
         
         WRITE(6,*) "init = ", init

         SELECT CASE (init)
         CASE (1)
            WRITE(6,*) 'Completely flat amplitude profile (energy ~ k^2)'
            WRITE(u_inputscopy,*)'Completely flat amplitude profile (energy ~ k^2)'
 
            DO i = 0, kxmax
               xk = kxnorm*i

               DO j = -kymax, kymax
                  yk = kynorm*j

                  ak = SQRT(xk*xk + yk*yk)
                  IF (ak .EQ. 0.0) ak = 1.0d0
                  IF ((ak .GE. kl) .AND. (ak .LT. ku)) THEN
                     den(i,j) = 2.*CMPLX(flreal,flima)
                     phi(i,j) = CMPLX(flreal,flima)
                  END IF

               END DO

            END DO

            !         den = CMPLX(0.0,0.0)
            !         phi = den
            denavg(0,0) = CMPLX(zmodeavg,0)

         CASE (3)
            WRITE(6,*) "Peaked initial profile with random phases."
            WRITE(u_inputscopy,*) "Peaked initial profile with random phases."

            DO j = -kymax, kymax
               yk = kynorm*j

               DO i = 0, kxmax
                  xk = kxnorm*i
                  ak = SQRT(xk*xk + yk*yk)
                  a1 = psi0/(1.d0 + ak**ratio)
                  a2 = a1*2.*SIN(2.d0*pi*ran2(irand))
                  a3 = SQRT(ABS(2.d0*a1*a1 - a2*a2))

                  den(i,j) = CMPLX(a3,a2)

               END DO

            END DO
         
            den(0,0) = CMPLX(0.0,0.0)   ! <-- no average initial flow
            den(0,1:kymax) = CONJG(den(0,-1:-kymax:-1))  ! enforce conjugacy

            phi = den
!            denavg = CMPLX(small,0)
!            denavg(0,0) = CMPLX(zmodeavg,0.d0)
!            denavg(1,0) = CMPLX(1.d0,0.d0)

            denavg = 5.0d0/samp * densource  ! gives an initial profile
            denavg(0,0) = CMPLX(zmodeavg,0.d0)
            !         denavg(0,0) = CMPLX(1000,0)
            denavg(0,1:kymax) = CONJG(denavg(0,-1:-kymax:-1))  ! enforce conjugacy

         CASE (4)
            WRITE(6,*) "Read in amplitudes from restart.bin file."
            WRITE(u_inputscopy,*) "Read in amplitudes from restart.bin file."

            ! a small error checking loop
            lexist = .FALSE.
            DO
               WRITE(6,'(/a)',advance='no') "Specify the restart.bin file (enter 'q' to quit): "
               READ(5,'(a)') restartfile_str
!               restartfile_str = "restart_betatf40.bin" ! for testing only!!
!               restartfile_str = "restart_betatf10.bin"
               IF (TRIM(ADJUSTL(restartfile_str)) == "q") THEN
                  WRITE(6,'(//a/)') "User quits from specifying restart file. BETA code aborted."
                  STOP
               END IF
               INQUIRE(file=TRIM(ADJUSTL(restartfile_str)), exist=lexist)
               IF (lexist .EQV. .TRUE.) EXIT ! break from loop 
            END DO
            
            WRITE(6,'(/a)', advance='no') "Specified restart.bin file: "
            WRITE(6,*) TRIM(ADJUSTL(restartfile_str))
            OPEN(unit=47, file=TRIM(ADJUSTL(restartfile_str)), status='old', form='unformatted')
            READ(47) den,phi,denavg
            CLOSE(47)

            den(0,1:kymax) = CONJG(den(0,-1:-kymax:-1))  ! enforce conjugacy
            phi(0,1:kymax) = CONJG(phi(0,-1:-kymax:-1))
            denavg(0,1:kymax) = CONJG(denavg(0,-1:-kymax:-1))

!            phi(0:kxmax,0) = czero ! special case only

!      OPEN(253,file="denavg_dump",status="replace")
!      WRITE(253,*) denavg(0:kxmax,0)
!      CLOSE(253)
!      STOP "Force STOP at pert"

         CASE (42)
            WRITE(6,*) "Zero perturbation case"
            WRITE(u_inputscopy,*) "Zero perturbation case"

            den = czero; phi = czero; denavg = czero
            denavg(0,0) = CMPLX(zmodeavg,0.d0)

         CASE DEFAULT

            CALL FINALIZE
            STOP 'WARNING: Must have INIT == 1'

         END SELECT

      ELSE
         WRITE(6,*) 'Read in amplitudes from file initspec'

         OPEN(unit=47,file='initspec.bin',status='old',form='unformatted')
         READ(47) denekold,denenold,phiekold,phiek1old,phienold
         READ(47) den,phi,denavg
         CLOSE(47)

!!$         OPEN(unit=47,file='initspec_den',status='old')
!!$         OPEN(unit=48,file='initspec_phi',status='old')
!!$         OPEN(unit=49,file='initspec_denavg',status='old')
!!$
!!$         DO j = -kymax, kymax
!!$            READ(47,*) (den(i,j),i=0,kxmax)
!!$            READ(48,*) (phi(i,j),i=0,kxmax)
!!$            READ(49,*) (denavg(i,j),i=0,kxmax)
!!$         END DO

!!$!C  disregards one half of the input
!!$         DO i=1,kxmax
!!$            READ(48,*) (temp_cmplx,j=-kymax,kymax)
!!$            READ(49,*) (temp_cmplx,j=-kymax,kymax)
!!$         END DO
!!$
!!$!C  retains one half of input
!!$         DO i = 0, kxmax
!!$            READ(48,*) (den(i,j),j=-kymax,kymax)
!!$            READ(49,*) (phi(i,j),j=-kymax,kymax)
!!$         END DO
      END IF

END SUBROUTINE pert



SUBROUTINE energy
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine energy                                           12/14/94
!C
!C   Calculates the total energy and enstrophy.
!C 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE variables
      IMPLICIT NONE
 
      INTEGER :: i, j
      REAL :: xk, yk, ak, aden, bden, es, aphi, bphi
      REAL :: denektemp, denentemp, phiektemp, phiek1temp, phientemp, denavgektemp
 
!C  Initialization
 
      denek = 0.0; denen = 0.0; es = TINY(0.0d0)
      phiek = 0.0; phiektemp = 0.0; phien = 0.0
      phiek1 = 0.0; phiek1temp = 0.0
      ek = 0.0
      denavgek = 0.0

!C  Reduced loops -----------------------------------------------
!C  -- first part calculates ky modes for kx = 0
!C  -- second part calculates the rest of modes

      xk = 0.0
      DO j = -kymax, -1
         yk = kynorm*j
         ak = xk*xk + yk*yk

         ! calculated different from one-field model
         denektemp = den(0,j)*CONJG(den(0,j))*2.0d0
         denentemp = den(0,j)*CONJG(den(0,j))*ak*2.0d0
         phiektemp = phi(0,j)*CONJG(phi(0,j))*2.0d0
         phiek1temp = phi(0,j)*CONJG(phi(0,j))*2.d0*(1.d0-rtepsi+ak*rhos*rhos)  ! <-- when need to calculate conserved quantity
         phientemp = phi(0,j)*CONJG(phi(0,j))*ak*2.0d0
         denavgektemp = denavg(0,j)*CONJG(denavg(0,j))*2.0d0

         denek = denek + denektemp; denen = denen + denentemp
         phiek = phiek + phiektemp; phien = phien + phientemp
         phiek1 = phiek1 + phiek1temp
         denavgek = denavgek + denavgektemp;
         ek = denek + phiek

      END DO

      DO j = -kymax, kymax
         yk = kynorm*j
         DO i = 1,kxmax
            xk = kxnorm*i
            ak = xk*xk + yk*yk

            denektemp = den(i,j)*CONJG(den(i,j))*2.0d0
            denentemp = den(i,j)*CONJG(den(i,j))*ak*2.0d0
            phiektemp = phi(i,j)*CONJG(phi(i,j))*2.0d0
            phiek1temp = phi(i,j)*CONJG(phi(i,j))*2.d0*(1.d0-rtepsi+ak*rhos*rhos)  ! <-- when need to calculate conserved quantity
            phientemp = phi(i,j)*CONJG(phi(i,j))*ak*2.0d0
            denavgektemp = denavg(i,j)*CONJG(denavg(i,j))*2.0d0

            denek = denek + denektemp; denen = denen + denentemp
            phiek = phiek + phiektemp; phien = phien + phientemp
            phiek1 = phiek1 + phiek1temp
            denavgek = denavgek + denavgektemp
            ek = denek + phiek

         END DO
      END DO
!C ---------------------------------------------------------------
 
!C  Calculate rates of energy and enstrophy production.
 
      aden = MAX((denen+denenold)/2.d0,es); bden = MAX((denek+denekold)/2.d0,es)
      aphi = MAX((phien+phienold)/2.d0,es); bphi = MAX((phiek+phiekold)/2.d0,es)

      dengamen = ((denen-denenold)/dt)/aden; dengamek = ((denek-denekold)/dt)/bden
      phigamen = ((phien-phienold)/dt)/aphi; phigamek = ((phiek-phiekold)/dt)/bphi
 
!C  Save old energies
 
      denekold = denek; denenold = denen
      phiekold = phiek; phiek1old = phiek1; phienold = phien
      denavgekold = denavgek
      ekold = ek
      

END SUBROUTINE energy



SUBROUTINE timead
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                                                                      C
!C  subroutine timead                                                   C
!C                                                                      C
!C      Calls the ODE solver vodpk.                                     C
!C                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE variables
      USE cvode_variables
      USE fft_variables
      USE tracers
      USE fileunits, ONLY: u_tracers
      USE prog_flags, ONLY: l_tracers,l_tracers_visit
      IMPLICIT NONE

      INTEGER :: ierr
      REAL :: tf
      REAL, DIMENSION(nn) :: y

      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax) :: dentempout,junkout,dendiff

      INTEGER :: i,j
!      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax) :: den_temp,phi_temp,denavg_temp

      tf = t + dt

!C   Convert den to y

!      WRITE(*,*) "denavg = ", denavg(:,0)

!      phi(0:kxmax,0) = czero ! special case only when zeroing out the ky = 0 modes

      CALL convert(den,phi,denavg,y)

      ! test convert routine
!      CALL invert(y,dentempout,junkout,junkout)
!      dendiff = den - dentempout
!      IF (SUM(dendiff) .NE. czero) THEN
!         STOP "Forced STOP in timead: CONVERT and INVERT function mismatch."
!      END IF

!      OPEN(253,file="den_invert_dump",status="replace")
!      DO j = -kymax, kymax
!         WRITE(253,*) dentempout(0:kxmax,j)
!      END DO
!      CLOSE(253)
!      STOP "Force STOP in timead"

!      denavg_temp = denavg
!      CALL invert(y,den_temp,phi_temp,denavg_temp)

!      OPEN(unit=64,file='denavg_dump',status='replace')
!      DO j = -kymax, kymax
!         WRITE(64,*) (ABS(denavg(i,j)-denavg_temp(i,j)),i=0,kxmax)
!      END DO

!      IF (SUM(SUM(ABS(den-den_temp),2),1) .NE. CMPLX(0.0,0.0)) WRITE(*,*) "den not equal"
!      IF (SUM(SUM(ABS(phi-phi_temp),2),1) .NE. CMPLX(0.0,0.0)) WRITE(*,*) "phi not equal"
!      IF (SUM(SUM(ABS(denavg-denavg_temp),2),1) .NE. 0.0) WRITE(*,*) "denavg not equal", SUM(ABS(denavg-denavg_temp),2)
!      STOP "Forced STOP!!!"

!C   Call CVODE to solve for a timestep
!C     CALL FCVODE(TOUT=tf, T=t, Y=y, ITASK=1, IER=ierr)
      CALL FCVODE(tf, t, y, 1, ierr)

!C   Check for problems (refer to reference for more outputs)
      WRITE(95,50) T, IOUT(3), IOUT(4), IOUT(9), ROUT(2)
50    FORMAT(/' t = ', E11.3, 3X, 'nsteps = ', I5,'  nfuncevals = ',I4,'  q = ', I2, '  h = ', E14.6)

      IF (IERR .NE. 0) THEN
         WRITE(6,60) IERR, IOUT(15)
60       FORMAT(///' SUNDIALS_ERROR: FCVODE returned IER = ', I5, /, &
              '                 Linear Solver returned IER = ', I5)
         CALL FCVFREE
         CALL FINALIZE
         STOP
      ENDIF

 
!C   Convert y back to den
 
      CALL invert(y,den,phi,denavg)
!      CALL check_meanfield(denavg,denavg)

!      OPEN(unit=567,file='dendiff',status='replace')
!      DO j = -kymax, kymax
!         WRITE(567,*) (ABS(denavg(i,j)-denavg_temp(i,j)),i=0,kxmax)
!      END DO
!      STOP "Forced stop!"

!C   Write tracers' data
!      CALL pwrite(tf,u_tposx,u_tposy,u_tvelx,u_tvely)
      itrack = itrack + 1      ! specify as one entry
!      CALL pwrite(itrack,tf,u_tracers)
      IF (l_tracers .EQV. .TRUE.) THEN 
         WRITE(u_tracers) itrack,tf,px,py,pvx,pvy
         IF (MOD(ncl,lypmaxdtstep) .EQ. 0) CALL plypmax   ! add iteration to the calculation of maximum lypunov exponent
      END IF

      ! update tracer positions for VISIT output (tracers are advanced with the fields coupled with the solver, not by the number of time iterations)
      IF (l_tracers_visit) CALL update_tracer_tail

END SUBROUTINE timead


SUBROUTINE convert(u,v,w,y)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  subroutine convert                                           10/19/11
!C
!C   Puts two complex two-dimensional arrays, u(kx,ky) and v(kx,ky)
!C   into a one-dimensional array, y.
!C 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE fft_variables
      IMPLICIT NONE
 
      INTEGER :: kxarr, kyarr, kount, kountm, kountd, kountdm, kounts, kountsm
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax), INTENT(in):: u, v, w
      REAL, DIMENSION(nn), INTENT(out):: y

!C  Initialization
      kount = kym2                          ! number of kx = 0 elements
!      kount = 0
!      kountm = kym1*kxmax + kym2            ! starting index for imaginary parts of density terms
      kountm = nn/(2*3)
      kountd = 2*kountm                     ! starting index for real parts of potential terms
      kountdm = 3*kountm                    ! starting index for imaginary parts of potential terms
      kounts = 4*kountm                     ! starting index for mean field
      kountsm = 5*kountm                    ! starting index for imaginary parts of mean filed terms

!C  Making densities occupy first portion
      y(1:kym2) = (/ (REAL( u(0,kyarr) ), kyarr = -kymax, 0) /)
      y(kountm+1:kountm+kym2) = (/ (AIMAG( u(0,kyarr) ), kyarr = -kymax, 0) /)
!C  Potentials occupy second portion
      y(kountd+1:kountd+kym2) = (/ (REAL( v(0,kyarr) ), kyarr = -kymax, 0) /)
      y(kountdm+1:kountdm+kym2) = (/ (AIMAG( v(0,kyarr) ), kyarr = -kymax, 0) /)
!C  Mean fields occupy third portion
      y(kounts+1:kounts+kym2) = (/ (REAL( w(0,kyarr) ), kyarr = -kymax, 0) /)
      y(kountsm+1:kountsm+kym2) = (/ (AIMAG( w(0,kyarr) ), kyarr = -kymax, 0) /)

!!$OMP PARALLEL DO PRIVATE(kount)
      DO kyarr = -kymax, kymax
         DO kxarr = 1, kxmax

            kount = kount + 1
!            kount = kym2 + (kxarr + kxmax*(kyarr+kymax))   ! this is done to use OpenMP directive
            y(kount) = REAL( u(kxarr,kyarr) )
            y(kountm+kount) = AIMAG( u(kxarr,kyarr) )
            y(kountd+kount) = REAL( v(kxarr,kyarr) )
            y(kountdm+kount) = AIMAG( v(kxarr,kyarr) )
            y(kounts+kount) = REAL( w(kxarr,kyarr) )
            y(kountsm+kount) = AIMAG( w(kxarr,kyarr) )

         END DO
      END DO
!!$OMP END PARALLEL DO

!C  Different arrangement
!!$      DO kyarr = -kymax, 0
!!$         DO kxarr = 0, kxmax
!!$
!!$            kount = kount + 1
!!$            y(kount) = REAL( u(kxarr,kyarr) )
!!$            y(kountm+kount) = AIMAG( u(kxarr,kyarr) )
!!$            y(kountd+kount) = REAL( v(kxarr,kyarr) )
!!$            y(kountdm+kount) = AIMAG( v(kxarr,kyarr) )
!!$            y(kounts+kount) = REAL( w(kxarr,kyarr) )
!!$            y(kountsm+kount) = AIMAG( w(kxarr,kyarr) )
!!$
!!$         END DO
!!$      END DO
!!$
!!$      DO kyarr = 1, kymax
!!$         DO kxarr = 1, kxmax
!!$
!!$            kount = kount + 1
!!$            y(kount) = REAL( u(kxarr,kyarr) )
!!$            y(kountm+kount) = AIMAG( u(kxarr,kyarr) )
!!$            y(kountd+kount) = REAL( v(kxarr,kyarr) )
!!$            y(kountdm+kount) = AIMAG( v(kxarr,kyarr) )
!!$            y(kounts+kount) = REAL( w(kxarr,kyarr) )
!!$            y(kountsm+kount) = AIMAG( w(kxarr,kyarr) )
!!$
!!$         END DO
!!$      END DO

      IF (kount .NE. kountm) THEN 
         CALL FINALIZE
         STOP "CONVERT routine: Number of data elements does not match."
      END IF

END SUBROUTINE convert



SUBROUTINE invert(y,u,v,w)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine invert                                          10/19/11
!C
!C   Puts a one-dimensional array into the 2-D complex arrays,
!C   u(kx,ky) and v(kx,ky).
!C 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE fft_variables
      USE variables, ONLY: czero
      IMPLICIT NONE

      INTEGER :: i, j
      INTEGER :: kxarr, kyarr, kount, kountm, kountd, kountdm, kounts, kountsm
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax), INTENT(out):: u, v, w
      REAL, DIMENSION(nn), INTENT(in):: y

!C  Initialization
      kount = kym2
!      kount = 0
      kountm = kym1*kxmax + kym2
      kountm = nn/(2*3)
      kountd = 2*kountm
      kountdm = 3*kountm
      kounts = 4*kountm
      kountsm = 5*kountm

!      u = czero; v = czero; w = czero

!C  Extract densities from first portion
      u(0,-kymax:0) = (/ (CMPLX( y(j), y(kountm+j) ),j=1,kym2) /)
      u(0,1:kymax) = CONJG( u(0,-1:-kymax:-1) )
!C  Extract potentials from second portion
      v(0,-kymax:0) = (/ (CMPLX( y(kountd+j), y(kountdm+j) ),j=1,kym2) /)
      v(0,1:kymax) = CONJG( v(0,-1:-kymax:-1) )
!C  Extract mean field from third portion
      w(0,-kymax:0) = (/ (CMPLX( y(kounts+j), y(kountsm+j) ),j=1,kym2) /)
      w(0,1:kymax) = CONJG( w(0,-1:-kymax:-1) )

!!$OMP PARALLEL DO PRIVATE(kount)
      DO kyarr = -kymax, kymax
         DO kxarr = 1, kxmax

            kount = kount + 1
!            kount = kym2 + (kxarr + kxmax*(kyarr+kymax))         ! for OpenMP directives
            u(kxarr,kyarr) = CMPLX( y(kount), y(kountm+kount) )
            v(kxarr,kyarr) = CMPLX( y(kountd+kount), y(kountdm+kount) )
            w(kxarr,kyarr) = CMPLX( y(kounts+kount), y(kountsm+kount) )

         END DO
      END DO
!!$OMP END PARALLEL DO

!C  Different arrangement
!!$      DO kyarr = -kymax, 0
!!$         DO kxarr = 0, kxmax
!!$
!!$            kount = kount + 1
!!$            u(kxarr,kyarr) = CMPLX( y(kount), y(kountm+kount) )
!!$            v(kxarr,kyarr) = CMPLX( y(kountd+kount), y(kountdm+kount) )
!!$            w(kxarr,kyarr) = CMPLX( y(kounts+kount), y(kountsm+kount) )
!!$
!!$         END DO
!!$      END DO
!!$
!!$      DO kyarr = 1, kymax
!!$         DO kxarr = 1, kxmax
!!$
!!$            kount = kount + 1
!!$            u(kxarr,kyarr) = CMPLX( y(kount), y(kountm+kount) )
!!$            v(kxarr,kyarr) = CMPLX( y(kountd+kount), y(kountdm+kount) )
!!$            w(kxarr,kyarr) = CMPLX( y(kounts+kount), y(kountsm+kount) )
!!$
!!$         END DO
!!$      END DO
!!$
!!$      u(0,1:kymax) = CONJG( u(0,-1:-kymax:-1) )  ! enforce conjugacy
!!$      v(0,1:kymax) = CONJG( v(0,-1:-kymax:-1) )
!!$      w(0,1:kymax) = CONJG( w(0,-1:-kymax:-1) )


END SUBROUTINE invert


SUBROUTINE FCVFUN(T, Y, YDOT, IPAR, RPAR, IER)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   SUBROUTINE FCVFUN: specifies right-hand side of diff.eq.
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE fft_variables
      IMPLICIT NONE

      INTEGER :: i, j
      REAL, INTENT(in):: t
      INTEGER, DIMENSION(:), INTENT(in):: ipar
      REAL, DIMENSION(:), INTENT(in):: rpar
      REAL, DIMENSION(nn), INTENT(in):: y
      REAL, DIMENSION(nn), INTENT(out):: ydot
      INTEGER, INTENT(out):: IER
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax):: u, v, denkp, phikp, w, denavgkp

!C  Convert y to den.
 
      CALL invert(y,u,v,w)
 
!C  Get y-prime.  (The right hand side)
 
 
!      OPEN(253,file="phi1_dump",status="replace")
!      WRITE(253,*) v
!      CLOSE(253)
!      STOP "Forced STOP in FCVFUN"

      CALL fden(t,u,v,w,denkp,phikp,denavgkp)
 
!C  Convert denkp back to yprime.

      CALL convert(denkp,phikp,denavgkp,ydot)

!C  Return success
      IER = 0


END SUBROUTINE FCVFUN


!!$SUBROUTINE FCVJTIMES (V, FJV, T, Y, FY, H, IPAR, RPAR, WORK, IER)
!!$!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$!C
!!$!C   SUBROUTINE FCVJTIMES: specifies an approximate Jacobian vector
!!$!C     -- currently, these are just linear terms
!!$!C
!!$!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!$      USE fft_variables, ONLY:nn,kxmax,kymax
!!$      USE variables, ONLY:philinarr1,philinarr2,denlinarr,nue,emass
!!$      IMPLICIT NONE
!!$  
!!$      REAL, DIMENSION(nn), INTENT(out):: FJV
!!$      REAL, DIMENSION(nn), INTENT(in):: V,T,Y,FY,WORK
!!$      REAL, DIMENSION(:), INTENT(in):: RPAR
!!$      REAL, INTENT(in):: H
!!$
!!$      INTEGER, DIMENSION(:), INTENT(in):: IPAR
!!$      INTEGER, INTENT(out):: IER
!!$
!!$      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax):: utemp, vtemp, denjac, phijac
!!$
!!$!C  transform into complex array
!!$      CALL invert(Y,utemp,vtemp)
!!$
!!$!C  apply all linear terms
!!$      phijac = philinarr1*vtemp + philinarr2*utemp
!!$      denjac = denlinarr*vtemp - nue*utemp/emass
!!$
!!$!C  transform back into vector form
!!$      CALL convert(denjac,phijac,FJV)
!!$
!!$!C  return success
!!$      IER = 0
!!$
!!$END SUBROUTINE FCVJTIMES



SUBROUTINE fden(tnow,uin,vin,win,denkp,phikp,denavgkp)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine fden                                             12/15/94
!C
!C   This calculates the derivative and also handles the particles.
!C 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE fft_variables
      USE variables
      USE prog_flags
      USE tracers, ONLY: told
      IMPLICIT NONE
 
      INTEGER :: i,j
      REAL, INTENT(in) :: tnow
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax), INTENT(IN):: uin, vin, win                   ! u = density, v = phi, w = mean field
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax) :: wtemp                                      ! temporary arrays, might need to be fixed later
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax), INTENT(OUT):: denkp, phikp, denavgkp
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax):: &
           dennl, phinl, denavgnl, philength, phisrc
!      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax):: vtemp

!!$      ! implementing a floor
!!$      CALL rfftwnd_f77_one_complex_to_real(plan2, kdenavg, rdenavg)
!!$      DO j = 1, nyreal
!!$         DO i = 1, nxreal
!!$            IF (rdenavg(i,j) .LT. 1.0E-10) THEN
!!$               rdenavg(i,j) = 1.0E-10
!!$            END IF
!!$         END DO
!!$      END DO
!!$      kdenavg = czero
!!$      CALL rfftwnd_f77_one_real_to_complex(plan1, rdenavg, kdenavg)
!!$      kdenavg = kdenavg/REAL(nxreal*nyreal)


      ! special case only ky = 0 modes are zeroed out
!      vtemp = v
!      vtemp(0:kxmax,0) = czero

!      OPEN(253,file="phi1_dump",status="replace")
!      WRITE(253,*) v
!      CLOSE(253)

!C   Check mean field for negative values first
!      w = win
!      wtemp = czero
      wtemp = win
!      IF (lm .NE. 0.0) CALL check_meanfield(win,wtemp)

!C   Calculate for non-linear terms.

!      CALL denscale(w,dengrad)                     ! uses current density values to get average local density scale length
      IF ((lp .NE. 0.0) .OR. (le .NE. 0.0) .OR. (lm .NE. 0.0)) THEN 
         CALL convolve(uin, vin+phi0, wtemp, dennl, phinl, denavgnl, philength, phisrc)    ! convolution routine for all non-linear terms
!      CALL convolve(u, vtemp, w, dennl, phinl, denavgnl, denlength, philength, phisrc)   ! special case only
      END IF
      IF (lm .EQ. 0.0) philength = philengthlinarr*vin
!      philength = philengthlinarr*vin   ! <-- over ride fixed gradient drive

!      WRITE(*,*) "dengrad = ", dengrad(0:kxmax,0)

!      phikp = 0.5*phinl*philinarr    ! <-- energy conserved, depends on how the energy summation is defined
!      denkp = 0.5*dennl*philinarr

      ! linear part test
!      phikp = (v1*denlength - v1*philength + nue*(u - v))*philinarr
!      denkp = -v1*denlength + nue*(v - u)
!      denavgkp = densource !- rhos3*denavgnl

!      phikp = ((v1*linarr*drivearr + nue)*u - (v1*linarr + nue + damparr)*v - 0.5*phinl)*philinarr    ! <-- this works for 40x40 modes
!      denkp = (nue)*v - (v1*linarr*drivearr + nue)*u + 0.5*dennl

      ! special case only when ky = 0 modes are zeroed out
!      phikp = (-1.d0*(1.d0 - rtepsi*dcoeff)*v1*philength - damparr*(vtemp + phi0) + nue*rtepsi*drivearr*u &   ! david's old equations
!           - nue*rtepsi*vtemp + 0.5d0*lp*dcoeff*phinl)*philinarr
!      denkp = -1.d0*v1*dcoeff*philength + nue*vtemp - nue*u - 0.5d0*le*dcoeff*dennl

      ! ** currently using these equations **
      phikp = (-1.d0*(1.d0 - rtepsi*dcoeff)*v1*philength - damparr*(vin + phi0)  &          ! david's old equations
            + nue*rtepsi*(drivearr*uin - vin) + 0.5d0*lp*dcoeff*phinl)*philinarr
      denkp = -1.d0*v1*dcoeff*philength + nue*(vin - uin) - 0.5d0*le*dcoeff*dennl

!      phikp = (-1.d0*(1.d0 - rtepsi)*v1*philength - damparr*(vin + phi0)  &          ! Variant #1 equations
!            + nue*(uin - rtepsi*vin) + 0.5d0*lp*phinl)*philinarr
!      denkp = -1.d0*rtepsi*v1*philength + nue*(rtepsi*vin - uin) - 0.5d0*le*dcoeff*dennl

!      phikp = (-1.d0*(1.d0 - rtepsi*dcoeff)*v1*philength - damparr*(vin + phi0)  &          ! Variant #2 equations
!            + nue*(uin - rtepsi*vin) + lp*phinl)*philinarr
!      denkp = -1.d0*v1*rtepsi*dcoeff*philength + nue*(rtepsi*vin - uin) - le*dcoeff*dennl

!      phikp = 0.001*v; denkp = 0.001*u   ! <-- linear test equations
!      phikp = rtepsi*u - v; denkp = v - u
!      phikp = -1.d0*damparr*v + nue*rtepsi*(drivearr*u - v)
!      denkp = nue*(v - u)

!      phikp = (-0.5d0*(1.d0 - rtepsi*dcoeff)*le*dennl - damparr*(v + phi0) + nue*rtepsi*drivearr*u &   ! <-- modified with source added
!           - nue*rtepsi*v + 0.5d0*lp*dcoeff*phinl)*philinarr
!      denkp = densource + nue*(v - u) - 0.5d0*le*dcoeff*dennl  ! <-- adding the source directly with the density perturbation

!      phikp(0:kxmax,0) = czero ! special case only when not evolving the ky = 0 modes

!      phikp = ( (rtepsi-1)*(densource + phisrc) - damparr*v + nue*rtepsi*drivearr*u + v1*rtepsi*dcoeff*philength &   ! modified equations with adhoc source
!           - nue*rtepsi*v - 0.5d0*lp*dcoeff*phinl)*philinarr
!      denkp = -1.d0*rtepsi*densource + nue*v - nue*u - 0.5d0*le*dcoeff*dennl

!      phikp = ( -1.0*(1.0 - rtepsi*dcoeff)*v1*philength + nue*(drivearr*u - rtepsi*v) - damparr*v - lp*phinl)*philinarr   ! <-- testing this, equations as derived
!      denkp = -1.0*v1*rtepsi*philength - nue*(u - rtepsi*v) - le*dennl

      IF (lm .NE. 0.0) denavgkp = densource - 0.5d0*lm*denavgnl - dampdenavgarr*wtemp  ! mean field needs to add a damping term
!      IF (lm .NE. 0.0) denavgkp = densource - lm*denavgnl - dampdenavgarr*win  ! Variant #2
!      IF (lm .NE. 0.0) denavgkp = czero  ! set source but not evolving mean field
!      IF (lm .NE. 0.0) denavgkp = -0.5d0*lm*denavgnl  ! mean field with no source
!      IF (lm .NE. 0.0) denavgkp = -0.5d0*lm*denavgnl - dampdenavgarr*w  ! mean field with damping only
!      IF (lm .NE. 0.0) denavgkp = densource  ! this is for testing the source

!      phikp = ( -1.0*(1.0 - rtepsi*dcoeff)*v1*philength + nue*(u - rtepsi*v) - damparr*v - lp*phinl)*philinarr   ! <-- Craddock's equations
!      denkp = -1.0*v1*rtepsi*dcoeff*philength - nue*(u - rtepsi*v) - le*dennl


!      IF (ncl .GT. 20) THEN  ! stop source after some iterations
!         denavgkp = -rhos3*denavgnl
!      ELSE
!         denavgkp = densource - rhos3*denavgnl
!      END IF

!C   Advance tracers with current field
      IF (l_tracers .EQV. .TRUE.) THEN
         CALL padvance(vin+phi0,tnow - told)
         told = tnow
      END IF

!!$      IF (ncl .EQ. 2) THEN
!!$         OPEN(253,file="phikp_dump",status="replace")
!!$         WRITE(253,*) phikp
!!$         CLOSE(253)
!!$         STOP "Forced STOP in fden"
!!$      END IF

END SUBROUTINE fden




SUBROUTINE linterm_construct
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine linterm_construct                                22/12/11
!C
!C   Construct linear portions of the propagator
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE variables
      USE prog_flags, ONLY: drivedamp, l_extflow
      IMPLICIT NONE

      LOGICAL :: temp_logic
      INTEGER :: i,j,kxarr,kyarr
      REAL :: xk, yk, ak, ak1, g1, d1, akmax, akmeanflim
      REAL :: kd1, kd3, kg1, kg3
!      REAL :: damp, drive
      
      WRITE(6,*) "Constructing array for linear terms..."

      philinarr = czero; damparr = czero; drivearr = czero; dampdenavgarr = czero;
      phi0 = czero

      ! extras
!      linarr = czero; linkdrive = czero; drivedenavgarr = czero

!C  Parameter checks
!      WRITE(6,'(/a)') "--> Parameter display... <--"
!      WRITE(6,*) "gama1 = ", gama1, " gama2 = ", gama2, " gama3 = ", gama3
!      WRITE(6,*) "din = ", din, " dmid = ", dmid, " dout = ", dout
!      WRITE(6,*) "rtepsi = ", rtepsi, " rhos = ", rhos, " nue = ", nue
!      WRITE(6,*) " v1 = ", v1
!      WRITE(6,*) "kxnorm = ", kxnorm, " kynorm = ", kynorm
!      WRITE(6,*) "kd1sq = ", kd1sq, "kd3sq = ", kd3sq
!      WRITE(6,'(a/)') "--> End parameter display <--"

      kg1 = SQRT(REAL(kg1sq)); kg3 = SQRT(REAL(kg3sq))  ! this is used for square drive and damp regions

      akmax = REAL(kxmax)*REAL(kxmax) + REAL(kymax)*REAL(kymax)
!      akmeanflim = REAL(kxmax-kdmeanf)*REAL(kxmax-kdmeanf) + REAL(kymax-kdmeanf)*REAL(kymax-kdmeanf)
      akmeanflim = kdmeanf
      
      ! external flow
      IF (l_extflow .EQV. .TRUE.) THEN            ! currently hardwired, but these quantities can be modified
         CALL make_extflow
      END IF

      DO kyarr = -kymax, kymax
         DO kxarr = 0, kxmax
            xk = kxnorm*kxarr
            yk = kynorm*kyarr
            ak = yk*yk + xk*xk
            ak1 = (1.0d0 - rtepsi + ak*rhos*rhos)      ! <-- new set of equations
!            ak1 = (1.0 - rtepsi)*(1 + ak*rhos*rhos)  ! <-- old set of equations
!            ak1 = ak
!            IF(ak1 .EQ. 0.) THEN   ! <-- this might not be necessary
!               ak1 = 1.0
!            ENDIF

!C  Phi terms

!C  Setup drive and damp terms
            SELECT CASE(drivedamp)

!C  Shell model for drive and damp (circular regions)
            CASE(1)
!               WRITE(6,'(a)') "Step functions for drive and damp."
               ! Specify damp regions
               IF (ak .LT. REAL(kg1sq)) THEN                        !kd1sq = lower cut-off for k^2-dependent drive (a circle)
                  damparr(kxarr,kyarr) = gama1
               ELSE IF (ak .GE. REAL(kg3sq)) THEN                   !kd3sq = upper cut-off for k^2-dependent drive
                  damparr(kxarr,kyarr) = gama3*ak*ak
               ELSE
                  damparr(kxarr,kyarr) = gama2*ak
               ENDIF

               ! Specify drive regions
               IF (ak .LT. REAL(kd1sq)) THEN
                  drivearr(kxarr,kyarr) = din
               ELSE IF (ak .GE. REAL(kd3sq)) THEN
                  drivearr(kxarr,kyarr) = dout
               ELSE
                  drivearr(kxarr,kyarr) = dmid
               END IF

               ! driving mean field for a certain region
!               IF (ak .GE. kxmax*kxmax) THEN
!                  drivedenavgarr(kxarr,kyarr) = gama3*ak*ak*(0.001d0 + ak*ak)
!               END IF

               ! no drive nor damp on (0,0) mode
               damparr(0,0) = czero; drivearr(0,0) = czero
               
!C  Construct the mean damping array
               IF (ak .GE. akmeanflim) THEN
                  dampdenavgarr(kxarr,kyarr) = dmeanf*ak  ! dissipation proportional to k^2
               END IF

!C  Shell model for drive and damp (square regions)
            CASE(2)
               ! Specify damp regions
               IF ((ABS(kyarr) .LE. kg1) .AND. (ABS(kxarr) .LE. kg1)) THEN ! (a square)
                  damparr(kxarr,kyarr) = gama1
               ELSE IF ((ABS(kyarr) .GT. kg3) .OR. (ABS(kxarr) .GT. kg3)) THEN
                  damparr(kxarr,kyarr) = gama3*ak*ak
               ELSE
                  damparr(kxarr,kyarr) = gama2*ak
               ENDIF

               ! Specify drive regions
               IF (ak .LT. REAL(kd1sq)) THEN
                  drivearr(kxarr,kyarr) = din
               ELSE IF (ak .GE. REAL(kd3sq)) THEN
                  drivearr(kxarr,kyarr) = dout
               ELSE
                  drivearr(kxarr,kyarr) = dmid
               END IF

               ! no drive nor damp on (0,0) mode
               damparr(0,0) = czero; drivearr(0,0) = czero

               ! damp the ky = 0 mode to counteract the self-consistent flows
!               IF (kyarr .EQ. 0) THEN
!                  damparr(kxarr,kyarr) = 0.5d0 !100.d0*gama1
!               END IF

!C  Construct the mean damping array
               IF (ak .GT. akmeanflim) THEN
!                  dampdenavgarr(kxarr,kyarr) = dmeanf*ak*(100.0d0 + 1.0d0*ak)   ! damping proportional to k^4
                  dampdenavgarr(kxarr,kyarr) = dmeanf*ak   ! damping proportional to k^2
               END IF

            CASE(3)
!C  Continuous variables
!               WRITE(6,'(a)') "Continuous functions for drive and damp."
               IF (ak .EQ. 0.0) THEN ! catching zeroth mode
                  damparr(0,0) = gama1
               ELSE
                  damparr(kxarr,kyarr) = gama1*ak + gama3*ak*ak  ! k^2 dominance in low-k, k^4 dominance in high-k
                  drivearr(kxarr,kyarr) = din/(1.0 + EXP(8.0*(4.0*SQRT(ak/akmax) - 1.0)))    ! driving is forced close to a step function
                  dampdenavgarr(kxarr,kyarr) = dmeanf/(1.0 + EXP(-8.0*(4.0*SQRT(ak/akmax) - 3.0)))  ! damping for the mean field
               END IF

            END SELECT

!            damparr(kxarr,kyarr) = (gama1 + gama2*ak + gama3*ak*ak)*ak

!            linarr(kxarr,kyarr) = - ic * yk
!            philinarr1(kxarr,kyarr) =  ic * rtepsi * yk !-1.0*damp*ak - ic*yk*v1 - nue*rtepsi <-- former phi term
!            philinarr2(kxarr,kyarr) = - ic * yk                 !drive*nue*rtepsi + ic*v1*rtepsi*din*yk <-- former den term
            philinarr(kxarr,kyarr) = 1.0/ak1

!C  y-derivative terms

            philengthlinarr(kxarr,kyarr) = - ic * yk

         END DO
      END DO


      WRITE(6,*) "Linear terms constructed."

!C  Linear arrays dump files
!!$      OPEN(96,file="philinarr1_dump",status="replace")
!!$      temp_logic = write_me(philinarr1,96)
!!$      CLOSE(96)
!!$      OPEN(96,file="philinarr2_dump",status="replace")
!!$      temp_logic = write_me(philinarr2,96)
!!$      CLOSE(96)
!!$      OPEN(96,file="philinarr3_dump",status="replace")
!!$      temp_logic = write_me(philinarr3,96)
!!$      CLOSE(96)
!!$      OPEN(96,file="denlinarr_dump",status="replace")
!!$      temp_logic = write_me(denlinarr,96)
!!$      CLOSE(96)

!!$      OPEN(96,file="damparr_dump",status="replace")
!!$      temp_logic = write_me(damparr,96)
!!$      CLOSE(96)
!!$      OPEN(96,file="drivearr_dump",status="replace")
!!$      temp_logic = write_me(drivearr,96)
!!$      CLOSE(96)
!!$      OPEN(96,file="philinarr_dump",status="replace")
!!$      temp_logic = write_me(philinarr,96)
!!$      CLOSE(96)
!!$      OPEN(96,file="philengthlinarr_dump",status="replace")
!!$      temp_logic = write_me(philengthlinarr,96)
!!$      CLOSE(96)

!      stop "Force STOP!"

CONTAINS
      LOGICAL FUNCTION write_me( some_arr,file_unit )
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C     LOGICAL FUNCTION write_me :
!C        -- writes linear array
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            USE fft_variables, ONLY: kxmax, kymax
            IMPLICIT NONE
            INTEGER :: i,j
            INTEGER, INTENT(in) :: file_unit
            COMPLEX, DIMENSION(0:kxmax,-kymax:kymax), INTENT(IN) :: some_arr
            
            DO j = -kymax, kymax
               WRITE(file_unit,*) (REAL(some_arr(i,j)*CONJG(some_arr(i,j))), i=0,kxmax)
!               WRITE(file_unit,*) (some_arr(i,j), i=0,kxmax)
            END DO

            write_me = .TRUE.

      END FUNCTION write_me

END SUBROUTINE linterm_construct



SUBROUTINE convolve(u, v, w, dennl, phinl, denavgnl, philength, phisrc)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine convolve                                         12/23/11
!C
!C     This subroutine handles the nonlinear convolutions
!C     for two non-linearities and one mean field.
!C     -- ExB => dennl, Pol => phinl
!C     -- average field => denavgnl
!C 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE variables
      USE prog_flags
      USE fft_variables
      USE omp_variables
!      USE convolve_module
      IMPLICIT NONE

      LOGICAL :: temp_logic

!C  Argument declarations.

      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax),INTENT(in):: u,v,w
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax),INTENT(out):: &
           dennl,phinl,denavgnl,philength,phisrc

!C  Declare FFT data.
 
      ! CMPLX arrays
      COMPLEX, DIMENSION(0:nxcmplx-1,0:nycmplx-1) :: &
           kden, &               ! stream function
           kphi, &               ! potential function
           kdenavg, &            ! mean field
           cnkx,cnky, &          ! components of gradients of dn/dy
           cphix, cphiy, &       ! components of the lapacian of potential
           cvex, cvey, &         ! components of ExB velocity
           cnavgkx, cnavgky, &   ! gradient of mean field dn/dx, dn/dy = 0
           cnavgkx1d, &          ! poloidal averages for mean field
           kdenavg1d, &
           cPol,cExB, &          ! Polarization and ExB non-linearity
           cConvt, &             ! Convective term for mean field
           kvdv                  ! feedback with 1/Ln terms

      COMPLEX, DIMENSION(0:nxcmplx-1,0:nycmplx-1) :: kdenavg_temp  ! temporary variable to compare how filtering is done

!      COMPLEX, DIMENSION(0:nxcmplx-1) ::  cnavgkx1d, kdenavg1d ! for gradient length scale calculation

      ! extras
!      COMPLEX, DIMENSION(0:nxcmplx-1) ::  cnavgky1d, cvey1d, cvex1d
      COMPLEX, DIMENSION(0:nxcmplx-1,0:nycmplx-1) :: kphicopy, ksource, cphisource

      ! REAL arrays
      REAL, DIMENSION(nxreal,nyreal) :: &
           nkx,nky,phix,phiy,vex,vey,navgkx,navgky,Pol,ExB,Convt,Ln,vdu,vdv,rdenavg, &
           denavg1d, navgkx1d

!      REAL, DIMENSION(nxreal) :: denavg1d, navgkx1d, navgky1d, dengradtemp

      ! extras
!      REAL, DIMENSION(nxreal) :: vey1d, vex1d, Convt1d
      REAL, DIMENSION(nxreal,nyreal) :: rphi, rsource, rphisource

!C  The components of the k-vector.

      REAL :: kx_conv(0:nxcmplx-1), ky_conv(0:nycmplx-1)

!C  Misc...
 
      INTEGER :: i, j
      REAL :: ksq
      CHARACTER(LEN=10) :: itstr

!C  Initialize the "k-vectors" -- these are used for the CONVOLVE routine

      kx_conv = 0.0; ky_conv = 0.0
      kx_conv(0:kxmax) = kxnorm*(/ (i,i=0,kxmax) /)
      ky_conv(0:kymax) = kynorm*(/ (j,j=0,kymax) /)
      ky_conv(nycmplx-kymax:nycmplx-1) = kynorm*(/ (j,j=-kymax,-1) /)

!C  Initialize arrays to zero

      kden = czero; kphi = czero
      cnkx = czero; cnky = czero
      cphix = czero; cphiy = czero
      cvex = czero; cvey = czero
      cPol = czero; cExB = czero
      kdenavg = czero; cnavgkx = czero; cnavgky = czero; cConvt = czero
      cnavgkx1d = czero; kdenavg1d = czero

      ! extras
!      cnavgky1d = czero; cvex1d = czero; cvey1d = czero
      ksource = czero

!C  Do the shuffle...

      kden(0:kxmax,0:kymax) = u(0:kxmax,0:kymax)
      kden(0:kxmax,nycmplx-kymax:nycmplx-1) = u(0:kxmax,-kymax:-1)

      kphi(0:kxmax,0:kymax) = v(0:kxmax,0:kymax)
      kphi(0:kxmax,nycmplx-kymax:nycmplx-1) = v(0:kxmax,-kymax:-1)

      kdenavg(0:kxmax,0:kymax) = w(0:kxmax,0:kymax)
      kdenavg(0:kxmax,nycmplx-kymax:nycmplx-1) = w(0:kxmax,-kymax:-1)

      ! extras
!      kphicopy = kphi
!      ksource(0:kxmax,0:kymax) = densource(0:kxmax,0:kymax)
!      ksource(0:kxmax,nycmplx-kymax:nycmplx-1) = densource(0:kxmax,-kymax:-1)

!      kdenavg1d(0:kxmax) = kdenavg(0:kxmax,0)
!      kdenavg1d(0:5) = kdenavg(0:5,0)  ! extracting only low k-modes

!      kdenavg(0,0) = denavg(0,0) ! fill in the (0,0) mode to preserve shift in the mean field profile

!C  Calculates stuff to be transformed to x-space . . .
!C   Order of i, j loop flipped - SD  
!C  The loops are redundant, they could be split to decrease computing time.

!!$OMP PARALLEL DO PRIVATE(ksq)
      DO j = 0, nycmplx-1
         DO i = 0, kxmax
            ksq = kx_conv(i)*kx_conv(i) + ky_conv(j)*ky_conv(j)

            cnkx(i,j) = ic * kx_conv(i) * kden(i,j)     ! fourier transform is d/dx -> -ikx
            cnky(i,j) = ic * ky_conv(j) * kden(i,j)

            cvex(i,j) = -1.d0 * ic * ky_conv(j) * kphi(i,j)     ! vx component carries an additional minus sign
            cvey(i,j) =         ic * kx_conv(i) * kphi(i,j)

            cphix(i,j) = -1.d0 * ic * kx_conv(i) * ksq * kphi(i,j)
            cphiy(i,j) = -1.d0 * ic * ky_conv(j) * ksq * kphi(i,j)

            cnavgkx(i,j) = ic * kx_conv(i) * kdenavg(i,j)
            cnavgky(i,j) = ic * ky_conv(j) * kdenavg(i,j)

            IF (j .EQ. 0) THEN
               cnavgkx1d(i,j) = ic * kx_conv(i) * kdenavg(i,j)
               kdenavg1d(i,j) = kdenavg(i,j)
            END IF

         END DO
      END DO
!!$OMP END PARALLEL DO

!!$      OPEN(253,file="kdenavg_dump",status="replace")
!!$      WRITE(253,*) kdenavg
!!$      CLOSE(253)
!!$
!!$      OPEN(253,file="cnavgkx1d_dump",status="replace")
!!$      WRITE(253,*) cnavgkx1d
!!$      CLOSE(253)
!!$
!!$      OPEN(253,file="kdenavg1d_dump",status="replace")
!!$      WRITE(253,*) kdenavg1d
!!$      CLOSE(253)
!!$
!!$      STOP "Forced STOP!"

      ! these are the terms for the gradient scale length -- filtering steps
!      cnavgkx1d(0:16) = cnavgkx(0:16,0)  ! d<n>/dx term, FFTW CMPLX -> REAL destroys CMPLX data
!      kdenavg1d(0:16) = kdenavg(0:16,0)  ! extracting only low k-modes

!      kdenavg1d(0:kxmax) = kdenavg(0:kxmax,0)  ! <-- working conditions
!      cnavgkx1d(0:kxmax) = cnavgkx(0:kxmax,0)
!      kdenavg1d(0:kxmax) = kden(0:kxmax,0)   ! <-- these are used for the modified equations
!      cnavgkx1d(0:kxmax) = cnkx(0:kxmax,0)

!      cnavgky1d(0:kxmax) = cnavgky(0:kxmax,0)
!      cvex1d(0:kxmax) = cvex(0:kxmax,0)
!      cvey1d(0:kxmax) = cvey(0:kxmax,0)

!      WRITE(*,*) cnavgkx

!C  Off to REAL space . . .

!$      IF (l_omp .EQV. .TRUE.) THEN ! switches to OpenMP commands
!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads, plan2, cnkx, nkx)
!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads, plan2, cnky, nky)
!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads, plan2, cphix, phix)
!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads, plan2, cphiy, phiy)
!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads, plan2, cvex, vex)
!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads, plan2, cvey, vey)
!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads, plan2, cnavgkx, navgkx)
!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads, plan2, cnavgky, navgky)
!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads, plan2, kdenavg, rdenavg)
!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads, plan2, kdenavg1d, denavg1d)
!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads, plan2, cnavgkx1d, navgkx1d)

      ! extra for source convolution
!!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads, plan2, kphicopy, rphi)
!!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads, plan2, ksource, rsource)

      ! extra arrrays for convolutions for kx ExB for mean field
!!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads, plan2_1d,cnavgky1d,navgky1d)
!!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads, plan2_1d,cvex1d,vex1d)
!!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads, plan2_1d,cvey1d,vey1d)


!$      ELSE ! uses serial commands

         CALL rfftwnd_f77_one_complex_to_real(plan2, cnkx, nkx)
         CALL rfftwnd_f77_one_complex_to_real(plan2, cnky, nky)
         
!!$      OPEN(253,file="phi_in_dump",status="replace")
!!$      WRITE(253,*) kphi
!!$      CLOSE(253)
!!$
!!$      OPEN(253,file="cphix_dump",status="replace")
!!$      WRITE(253,*) cphix
!!$      CLOSE(253)
!!$
!!$      OPEN(253,file="cphiy_dump",status="replace")
!!$      WRITE(253,*) cphiy
!!$      CLOSE(253)
!!$
!!$      OPEN(253,file="cvex_dump",status="replace")
!!$      WRITE(253,*) cvex
!!$      CLOSE(253)
!!$
!!$      OPEN(253,file="cvey_dump",status="replace")
!!$      WRITE(253,*) cvey
!!$      CLOSE(253)

         CALL rfftwnd_f77_one_complex_to_real(plan2, cphix, phix)
         CALL rfftwnd_f77_one_complex_to_real(plan2, cphiy, phiy)

         CALL rfftwnd_f77_one_complex_to_real(plan2, cvex, vex)
         CALL rfftwnd_f77_one_complex_to_real(plan2, cvey, vey)

         CALL rfftwnd_f77_one_complex_to_real(plan2, cnavgkx, navgkx)
         CALL rfftwnd_f77_one_complex_to_real(plan2, cnavgky, navgky)

         CALL rfftwnd_f77_one_complex_to_real(plan2, kdenavg, rdenavg)

         CALL rfftwnd_f77_one_complex_to_real(plan2,kdenavg1d,denavg1d)
         CALL rfftwnd_f77_one_complex_to_real(plan2,cnavgkx1d,navgkx1d)


!!$      OPEN(253,file="phix_dump",status="replace")
!!$      WRITE(253,*) phix
!!$      CLOSE(253)
!!$
!!$      OPEN(253,file="phiy_dump",status="replace")
!!$      WRITE(253,*) phiy
!!$      CLOSE(253)
!!$
!!$      OPEN(253,file="vex_dump",status="replace")
!!$      WRITE(253,*) vex
!!$      CLOSE(253)
!!$
!!$      OPEN(253,file="vey_dump",status="replace")
!!$      WRITE(253,*) vey
!!$      CLOSE(253)

         ! extras
!         CALL rfftwnd_f77_one_complex_to_real(plan2_1d,cnavgky1d,navgky1d)  ! these are used for 1D transforms for poloidal averages
!         CALL rfftwnd_f77_one_complex_to_real(plan2_1d,cvex1d,vex1d)
!         CALL rfftwnd_f77_one_complex_to_real(plan2_1d,cvey1d,vey1d)
!         CALL rfftwnd_f77_one_complex_to_real(plan2, kphicopy, rphi)
!         CALL rfftwnd_f77_one_complex_to_real(plan2, ksource, rsource)

!$      END IF

!         WRITE(*,*) "cnavgkx = ", cnavgkx(0:kxmax,0)
!         WRITE(*,*) "nxavg = ", nxavg

!C  Check for zero values in the mean field
      DO j = 1, nyreal
         DO i = 1, nxreal
            IF (rdenavg(i,j) .LT. 0.d0) THEN
               WRITE(*,'(/a//)') "Unphysical mean field."
               WRITE(*,*) "Value at: i = ", i, " j = ", j
               CALL FINALIZE
               STOP 'Process terminated due to negative value.'
            END IF
         END DO
      END DO

!C  Calculate the nonlinear terms
       
      Pol = vex * phix + vey * phiy
      ExB = vex * nkx + vey * nky
      Convt = vex * navgkx + vey * navgky

!      OPEN(253,file="rPol_dump",status="replace")
!      WRITE(253,*) Pol
!      CLOSE(253)

      ! extras
!      Convt = (rdensource + rhos3 * (vex * navgkx + vey * navgky)) / rdenavg ! <-- navgky is zero since the beginning
!      Convt1d = vex1d * navgkx1d + vey1d * navgky1d
!      rphisource = rphi * rsource

!C  Mean field flag
!      IF (lm .EQ. 0.0) THEN
!         dengradtemp = -1.0d0
!      ELSE

         ! catch zero Ln values
!         dengradtemp = denavg1d / navgkx1d          ! length values
!         dengradtemp = navgkx1d / denavg1d     ! inverse scale length 1/Ln
!         dengradtemp = denavg1d/nxavg + Lnc   ! adding a critical gradient
!         dengradtemp = EXP(-1.0*(denavg1d/nxavg - Lnc))
!         DO i = 1, nxreal
!            IF ( isnan(dengradtemp(i)) ) THEN   ! catch NaN values
!               dengradtemp(i) = HUGE(0.0d0)
!               dengradtemp(i) = TINY(0.0d0)  ! switching to calculate the drive differently
!               WRITE(*,*) "Setting dengradtemp(i) = 0.0 at i = ", i
!            dengradtemp(i) = ABS(nxavg(i))/denavg1d(i)  ! this quantity calculates 1/Ln
!            dengradtemp(i) = 1/(denavg1d(i)/nxavg(i) - Lnc)
!               dengradtemp(i) = 1.0d0/dengradtemp(i)
!            END IF
!         END DO
!      END IF
 
      rdenavg1d = denavg1d(1:nxreal,nyreal/2)  ! stores the feedback field
      rLn = navgkx1d(1:nxreal,nyreal/2)/denavg1d(1:nxreal,nyreal/2)     ! stores lengths to write to file

!      DO j = 1, nyreal
!         Ln(1:nxreal,j) = -1.0d0 * dengradtemp             ! replicate to form an array, this is supposed to be the average poloidal length scale, negative sign is for direction of gradient
!         Convt(1:nxreal,j) = Convt1d
!      END DO

!      vdv = -1.0d0 * vex * navgkx1d/denavg1d             ! divide the length scale, the negative sign is to add a negative sign from "vex" specified above
      vdv = -1.0d0 * vex * navgkx/rdenavg                ! try local gradient scale length instead
!      vdv = -1.0d0 * vex * navgkx/denavg1d 

!      IF ((MOD(ncl,500) .EQ. 0) .OR. (ncl .EQ. 1)) THEN
!         CALL get_iterstr(itstr,ncl,LEN(itstr))
!         OPEN(unit=897, file="dengrad_dump."//itstr, status="replace")
!         DO i = 1, nxreal
!            WRITE(897,*) Ln(i,1),denavg1d(i),nxavg(i)
!         END DO
!         CLOSE(897)
!      END IF
!      STOP


!C  Back to Fourier space . . .

!$      IF (l_omp .EQV. .TRUE.) THEN
!$         CALL rfftwnd_f77_threads_one_real_to_complex(nthreads, plan1, Pol, cPol)
!$         CALL rfftwnd_f77_threads_one_real_to_complex(nthreads, plan1, ExB, cExB)
!$         CALL rfftwnd_f77_threads_one_real_to_complex(nthreads, plan1, Convt, cConvt)
!$         CALL rfftwnd_f77_threads_one_real_to_complex(nthreads, plan1, vdv, kvdv)

      ! extras
!!$         CALL rfftwnd_f77_threads_one_real_to_complex(nthreads, plan1, rphisource, cphisource)

!$      ELSE

         CALL rfftwnd_f77_one_real_to_complex(plan1, Pol, cPol)
         CALL rfftwnd_f77_one_real_to_complex(plan1, ExB, cExB)
         CALL rfftwnd_f77_one_real_to_complex(plan1, Convt, cConvt)
         CALL rfftwnd_f77_one_real_to_complex(plan1, vdv, kvdv)

         ! extras
!         CALL rfftwnd_f77_one_real_to_complex(plan1, rphisource, cphisource)

!$      END IF


!C  Do the UN-shuffle...

!C   The kx >= 0, ky >=0 modes.

      dennl(0:kxmax,0:kymax) = cExB(0:kxmax,0:kymax)
      phinl(0:kxmax,0:kymax) = cPol(0:kxmax,0:kymax)
      denavgnl(0:kxmax,0:kymax) = cConvt(0:kxmax,0:kymax)
      philength(0:kxmax,0:kymax) = kvdv(0:kxmax,0:kymax)

!      phisrc(0:kxmax,0:kymax) = cphisource(0:kxmax,0:kymax)

!C   The kx >= 0, ky < 0 modes.

      dennl(0:kxmax,-kymax:-1) = cExB(0:kxmax,nycmplx-kymax:nycmplx-1)
      phinl(0:kxmax,-kymax:-1) = cPol(0:kxmax,nycmplx-kymax:nycmplx-1)
      denavgnl(0:kxmax,-kymax:-1) = cConvt(0:kxmax,nycmplx-kymax:nycmplx-1)
      philength(0:kxmax,-kymax:-1) = kvdv(0:kxmax,nycmplx-kymax:nycmplx-1)

!      phisrc(0:kxmax,-kymax:-1) = cphisource(0:kxmax,nycmplx-kymax:nycmplx-1)

!C  FFTW is not normalized
      dennl = dennl/REAL(nxreal*nyreal)
      phinl = phinl/REAL(nxreal*nyreal)
      denavgnl = denavgnl/REAL(nxreal*nyreal)
      philength = philength/REAL(nxreal*nyreal)
      
!      phisrc = phisrc/REAL(nxreal*nyreal)

! limit the feedback to low-k from mean field to the fluctuations, apply a low pass filter to the linear drive term
!      denlength(9:kxmax,-kymax:kymax) = czero  ! testing with a certain number of modes
!      denlength(1:kxmax,9:kymax) = czero
!      denlength(1:kxmax,-kymax:-9) = czero

!      philength(9:kxmax,-kymax:kymax) = czero
!      philength(1:kxmax,9:kymax) = czero
!      philength(1:kxmax,-kymax:-9) = czero

      ! take only the poloidal averaged modes for the mean field convolution
!      denavgnl(0:kxmax,1:kymax) = czero
!      denavgnl(0:kxmax,-kymax:-1) = czero

      linkdrive = philength  ! store variable for output

END SUBROUTINE convolve


SUBROUTINE makesource
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   SUBROUTINE makesource: creates an external source term for the mean field
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE fft_variables, ONLY: kxmax,kymax,nxcmplx,nycmplx,nxreal,nyreal,kxm2, plan2
      USE variables, ONLY: densource,kxnorm,rdensource,samp,sscale,czero
      IMPLICIT NONE

      INTEGER :: i,j
      COMPLEX :: source1, sink1
      REAL :: x1, x2, amp, scale
!      COMPLEX, DIMENSION(0:nxcmplx-1,0:nycmplx-1) :: kdensource

      densource = czero; source1 = czero; sink1 = czero
      x1 = 0.25d0; x2 = 0.75d0
      amp = samp; scale = sscale*REAL(kxmax)

      DO i = 0, kxmax
         CALL getgauss(source1,kxnorm*REAL(i),1.d0,x1,scale)
         CALL getgauss(sink1,kxnorm*REAL(i),-1.d0,x2,scale)
         densource(i,0) = source1 + sink1
      END DO

      densource = 0.5d0*amp*densource

!      kdensource = CMPLX(0.0,0.0)
!      kdensource(0:kxmax,0) = densource(0:kxmax,0)

!      CALL rfftwnd_f77_one_complex_to_real(plan2,kdensource,rdensource)

      ! replicate for all ky modes
!!$      DO j = -kymax, -1
!!$         densource(:,j) = densource(:,0)
!!$      END DO
!!$      DO j = 1, kymax
!!$         densource(:,j) = densource(:,0)
!!$      END DO

END SUBROUTINE makesource

SUBROUTINE getgauss(sorrow,k,mag,x0,sigma)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   SUBROUTINE get_gauss: returns value of shifted and scaled gaussian
!C                         for a given k
!C     -- parameter inputs: mag   => magnitude
!C                          x0    => shift in REAL space
!C                          sigma => scale parameter 
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE variables, ONLY: pi
      IMPLICIT NONE
      
      REAL, INTENT(in) :: k,mag,x0,sigma
      COMPLEX, INTENT(out) :: sorrow
      REAL :: temp_fact

!      sorrow = mag*SQRT(1/(2*pi*sigma)) * CMPLX(COS(2*pi*k*x0),-1.0*SIN(2*pi*k*x0)) * EXP(-0.25*k*k/sigma)
!      sorrow = mag*SQRT(pi/sigma) * CMPLX(COS(2.0d0*pi*k*x0),-1.0d0*SIN(2.0d0*pi*k*x0)) * EXP(-0.25d0*k*k/sigma)
      sorrow = mag / SQRT(pi*sigma) * CMPLX(COS(2.0d0*pi*k*x0),-1.0d0*SIN(2.0d0*pi*k*x0)) * EXP(-0.25d0*k*k/sigma)
!C                                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ phase due to shift of gaussian curve

END SUBROUTINE getgauss


SUBROUTINE make_extflow
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   SUBROUTINE make_extflow: create an external flow profile
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE fft_variables, ONLY: kxmax,kymax
      USE variables, ONLY: phi0,phi0amp,czero,pi,ic
      IMPLICIT NONE

      INTEGER :: kx,ky
      REAL :: x1,x2,coeff
      
      phi0 = czero

!  These parameters give flow in the first mode of ky
!         phi0(0,1) = phi0amp*CMPLX(1.d0,0.d0)     
      phi0(1,0) = phi0amp*CMPLX(0.d0,-1.d0)  ! negative sine function
!      phi0(1,0) = phi0amp*CMPLX(1.d0,0.d0)  ! cosine function
!      phi0(0,1) = phi0amp*CMPLX(1.d0,0.d0)  ! for circular flow


!  Step function flow
!      x1 = 2.d0*pi*0.25; x2 = 2.d0*pi*0.75
!      coeff = -2.d0*phi0amp/SQRT(pi)
!      DO kx = 1, kxmax   ! implementing a step function for velocity
!         phi0(kx,0) = coeff*CMPLX(0.d0,SIN(REAL(kx)*x1)/(REAL(kx)*x1))
!      END DO

END SUBROUTINE make_extflow


SUBROUTINE FINALIZE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine final                                            12/14/94
!C
!C   Close all files
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE visit
      USE variables
      USE fft_variables
      USE fileunits
      USE prog_flags
!      USE convolve_module
      IMPLICIT NONE
!C
!C   A little closure....
!C 
      
      DEALLOCATE(x_visit, y_visit)
      DEALLOCATE(kx_visit, ky_visit)
      
      WRITE(6,100)
100   FORMAT(4x,' das ende')

      CLOSE( u_energy ); !CLOSE( 48 ); close( 49 ); CLOSE( 95 )
      CLOSE( u_inputscopy )
      CLOSE( u_denpolavg ); CLOSE( u_phipolavg ); CLOSE( u_denavgpolavg ); CLOSE( u_Ln ); CLOSE( u_denavg1d )

      CALL rfftwnd_f77_destroy_plan(plan1); CALL rfftwnd_f77_destroy_plan(plan2)
      CALL rfftwnd_f77_destroy_plan(plan1_1d); CALL rfftwnd_f77_destroy_plan(plan2_1d)

      ! optional outputs
      IF (l_specout .EQV. .TRUE.) CLOSE( u_specout )
      IF (l_nltrans .EQV. .TRUE.) CLOSE( u_nltrans )
      IF (l_tracers .EQV. .TRUE.) CALL tracers_final

      CALL FCVFREE
 
END SUBROUTINE FINALIZE



SUBROUTINE writerestart
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C       *********************************************************
!C       *                                                       *
!C       *  Writes restart file                  (DEN 6/17/04)   *
!C       *  and writes particle restart file	(DEN 7/27/05)   *
!C       *  now accepte jet values for constant jet - SD         *
!C       *********************************************************
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE variables, ONLY: den,phi,denavg,ncl
      IMPLICIT NONE

      INTEGER :: i,j
      character(LEN=10) :: itstr

      CALL get_iterstr(itstr,ncl,LEN(itstr))
      OPEN(20,file='restart.'//TRIM(itstr)//'.bin',status='replace',form='unformatted')  ! this is a binary file 
      WRITE(20) den,phi,denavg
!      CALL write_spec(den,20)
      CLOSE(20)

!      OPEN(20,file='restart_phi.'//TRIM(itstr),status='replace') 
!      CALL write_spec(phi,20)
!      CLOSE(20)

END SUBROUTINE writerestart

SUBROUTINE writeinit
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  subroutine: writeinit
!C    Writes initspec binary files for continuation only. Replaces the old files.
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE variables, ONLY: den,phi,denavg,denekold,denenold,phiekold,phiek1old,phienold
      IMPLICIT NONE
      OPEN(unit=47,file='initspec.bin',status='replace',form='unformatted')
      WRITE(47) denekold,denenold,phiekold,phiek1old,phienold
      WRITE(47) den,phi,denavg
      CLOSE(47)

END SUBROUTINE writeinit

SUBROUTINE write_spec(specin,funit)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  SUBROUTINE write_spec: writes spectrum files, CMPLX type
!C    ** only half the spectrum is retained **
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE fft_variables
      USE variables

      IMPLICIT NONE
      INTEGER :: i, j , n
      INTEGER, INTENT(in) :: funit
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax), INTENT(in) :: specin
!      COMPLEX :: uout(kxm1,kym1)

!      uout(kxm2:kxm1,1:kym1) = u(0:kxmax,-kymax:kymax)

!!$!C  Fill in the conjugate modes...
!!$ 
!!$!C    The kx < 0, ky <= 0 modes.
!!$      uout(1:kxmax,1:kym2) = CONJG( uout(kxm1:kxm2+1:-1,kym1:kym2:-1) )
!!$
!!$!C   The kx < 0, ky > 0 modes.
!!$      uout(1:kxmax,kym2+1:kym1) = CONJG( uout(kxm1:kxm2+1:-1,kymax:1:-1) )
!!$
!!$      DO i = 1, kxm1
!!$         WRITE(file_unit,*) (uout(i,j),j=1,kym1)
!!$      END DO

!      WRITE(6,*) "File: ",file_unit," written."

      WRITE(funit) specin

END SUBROUTINE write_spec


SUBROUTINE writeit(fileunit)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                                                                      C
!C   subroutine writeit                                                 C
!C                                                                      C
!C      Writes the current values of the parameters                     C
!C      ncl, t, ek, en, gamk, gamen to the file "energy".               C
!C                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE variables
      USE prog_flags, ONLY: lcont
      IMPLICIT NONE

      INTEGER, INTENT(in) :: fileunit

      WRITE(fileunit,200) ncl,t,denek,denen,dengamek,dengamen,phiek,phiek1,phien,phigamek,phigamen,denavgek,ek
      !formatting altered to include more decimal places-SD(03/16/2010)
!200   FORMAT(i6,1x,e16.8,1x,e16.8,1x,e16.8,1x,e16.8,1x,e16.8)
200   FORMAT(I10,1x,12(1pe16.9,1x))
 
END SUBROUTINE writeit


SUBROUTINE specoutinit
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   SUBROUTINE specoutinit
!C
!C   Checks for continuation and open files.
!C              
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE variables, ONLY: nspecmodesx, nspecmodesy, specmodesx, specmodesy
      USE fileunits, ONLY: u_specout
      USE prog_flags, ONLY: lcont
      IMPLICIT NONE

      ! check for continuation
      IF (lcont .EQV. .FALSE.) THEN  ! a fresh run
         OPEN(unit=u_specout,file='specout.bin',status='replace',form='unformatted')
         WRITE(u_specout) nspecmodesx,nspecmodesy
         WRITE(u_specout) specmodesx,specmodesy
      ELSE
         OPEN(unit=u_specout,file='specout.bin',status='old', position='append',form='unformatted')
      END IF


END SUBROUTINE specoutinit



SUBROUTINE write_specout
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   SUBROUTINE write_specout
!C
!C   Outputs spectral energy of a given mode to corresponding file.
!C              
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE variables
      USE fileunits, ONLY: u_specout
      IMPLICIT NONE

      INTEGER :: i,j,n
      COMPLEX, DIMENSION(nspecmodesx*nspecmodesy) :: denspecout,phispecout,denavgspecout

      DO j = 1, nspecmodesy
         DO i = 1, nspecmodesx
            n = i + (j-1)*nspecmodesy  ! flattened index
            denspecout(n) = den(specmodesx(i),specmodesy(j))
            phispecout(n) = phi(specmodesx(i),specmodesy(j))
            denavgspecout(n) = denavg(specmodesx(i),specmodesy(j))
         END DO
      END DO

      WRITE(u_specout) ncl,t,denspecout,phispecout,denavgspecout
!      WRITE(6,*) ncl,t,denavgspecout

END SUBROUTINE write_specout



SUBROUTINE NLTRANSINIT
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   SUBROUTINE NLTRANSINIT
!C
!C   Checks for continuation and open files.
!C              
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE variables, ONLY: nspecmodesx, nspecmodesy, specmodesx, specmodesy
      USE fileunits, ONLY: u_nltrans
      USE prog_flags, ONLY: lcont
      IMPLICIT NONE

      ! check for continuation
      IF (lcont .EQV. .FALSE.) THEN  ! a fresh run
         OPEN(unit=u_nltrans,file='nltrans.bin',status='replace',form='unformatted')
         WRITE(u_nltrans) nspecmodesx,nspecmodesy
         WRITE(u_nltrans) specmodesx,specmodesy
      ELSE
         OPEN(unit=u_nltrans,file='nltrans.bin',status='old', position='append',form='unformatted')
      END IF

END SUBROUTINE NLTRANSINIT


SUBROUTINE WRITE_NLTRANS
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  SUBROUTINE NLTRANS
!C   
!C  Subroutine to calculate the nonlinear transfer terms. This is similar
!C  to the WRITE_SPECOUT routine
!C  *** this section is not complete!
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE fft_variables
      USE variables
      USE precmod
      USE fileunits, ONLY: u_nltrans
      IMPLICIT NONE

      INTEGER :: i,j,n,kxp,kyp,kydiff,kxdiff
      COMPLEX :: temp_polnltrans, temp_exbnltrans, phi_conjg, temp_den, temp_phi
      REAL :: kcrosskp, kkpsqdiff
      REAL, DIMENSION(nspecmodesx*nspecmodesy) :: polnltrans, exbnltrans
      

      ! cycle through the output modes
      DO j = 1, nspecmodesy
         DO i = 1, nspecmodesx
            n = i + (j-1)*nspecmodesy  ! flattened output index0

            ! cycle through all modes to get the sum
            temp_polnltrans = czero; temp_exbnltrans = czero
            phi_conjg = CONJG(phi(specmodesx(i),specmodesy(j)))

            ! set the difference indexing
            IF ((specmodesy(j) - kymax) .LT. -kymax) THEN
               kydiff = kymax
            ELSE
               kydiff = 0
            END IF
               
!!$OMP PARALLEL DO PRIVATE(temp_den,temp_phi) REDUCTION(+:temp_polnltrans,temp_exbnltrans)
            DO kyp = -kymax, kymax
               DO kxp = 0, kxmax
                  
!                  WRITE(*,*) "kxp = ", kxp, " kyp = ", kyp

                  ! calculate some intermediate terms
!                  kdiff = specmodesx(i) - kxp 
                  kcrosskp = REAL(kyp*specmodesx(i)-kxp*specmodesy(j))
                  kkpsqdiff = REAL( (specmodesx(i)-kxp)*(specmodesx(i)-kxp) &
                       + (specmodesy(j)-kyp)*(specmodesy(j)-kyp) - kxp*kxp - kyp*kyp)

                  ! use reality condition to recover data, this is for kx direction only
                  IF ((kxp-specmodesx(i)) .LT. 0) THEN
                     temp_den = CONJG(den(specmodesx(i)-kxp,kyp-specmodesy(j)))
                     temp_phi = CONJG(phi(specmodesx(i)-kxp,kyp-specmodesy(j)))
                  ELSE
                     temp_den = den(kxp-specmodesx(i),kyp-specmodesy(j))
                     temp_phi = phi(kxp-specmodesx(i),kyp-specmodesy(j))
                  END IF

                  temp_polnltrans = temp_polnltrans &
                       + kcrosskp*phi(kxp,kyp)*phi_conjg*temp_den
                  temp_exbnltrans = temp_exbnltrans &
                       + kcrosskp*kkpsqdiff*phi(kxp,kyp)*phi_conjg*temp_phi

               END DO
            END DO
!!$OMP END PARALLEL DO  

            polnltrans(n) = lp*dcoeff*REAL(temp_polnltrans)
            exbnltrans(n) = le*dcoeff*REAL(temp_exbnltrans)

         END DO
      END DO

      WRITE(u_nltrans) ncl,t,exbnltrans,polnltrans

END SUBROUTINE WRITE_NLTRANS


SUBROUTINE INIT_VISIT
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   SUBROUTINE INIT_VISIT: initialize outputs to Silo for ViSit
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE fft_variables, ONLY: nxreal,nyreal,kxm2,kym1,kxmax,kymax
      USE visit
      IMPLICIT NONE

      INTEGER :: j
      REAL :: dk

      ndims = 2
      nx_visit = nxreal; ny_visit = nyreal
      nkx_visit = kxm2; nky_visit = kym1
      ALLOCATE(kdims(ndims), xdims(ndims))

      ALLOCATE(x_visit(nx_visit), y_visit(ny_visit))
      ALLOCATE(kx_visit(nkx_visit), ky_visit(nky_visit))
!
!     K-arrays
!
      kdims(1) = nkx_visit; kdims(2) = nky_visit 
      dk = kxmax/REAL(nkx_visit - 1)
      DO j = 1, nkx_visit
         kx_visit(j) = REAL(j - 1)*dk
      ENDDO 
      
      dk = 2.0*kymax/REAL(nky_visit - 1)
!  swapped to coincide with data storage <-- didn't work?
      DO j = 1, nky_visit
!      DO j = nky_visit, 1, -1
         ky_visit(j) = -kymax + REAL(j - 1)*dk
      ENDDO
!
!     Spatial arrays
!
      xdims(1) = nx_visit; xdims(2) = ny_visit 
      DO j = 1, nx_visit
         x_visit(j) = REAL(j-1)/REAL(nx_visit - 1)
      ENDDO
      DO j = 1, ny_visit
         y_visit(j) = REAL(j-1)/REAL(ny_visit - 1)
      ENDDO
      
!      ivisit = 0 ! <-- this has been moved to "check_cont" subroutine
      
END SUBROUTINE INIT_VISIT



SUBROUTINE WRITE_VISIT
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   SUBROUTINE WRITE_VISIT: write outputs to Silo for ViSit
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE visit
      USE variables
      USE fft_variables
      USE prog_flags
      USE omp_variables
      USE prog_flags, ONLY:l_tracers_visit
      IMPLICIT NONE

      INTEGER :: ierr, err, optlistid, dbfile
      CHARACTER(LEN=6):: citer
      CHARACTER(LEN=10):: file2
      CHARACTER(LEN=16):: filename   ! the length of this string must coincide with the string length for the database
      INTEGER :: i, j
      REAL ak, xk, yk, ksq

      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax) :: cxcorr_visit
      COMPLEX, DIMENSION(0:nxcmplx-1,0:nycmplx-1) :: kden,kvort,kphi,kdenavg,kden2,klindrive
      REAL, DIMENSION(nxreal,nyreal) :: den_visit,vort_visit,phi_visit,denavg_visit, &
           den2_visit,lindrive_visit,tpos_visit
      REAL, DIMENSION(0:kxmax,-kymax:kymax) :: &
           denekk_visit, phiekk_visit, denavgekk_visit, phiekkphase_visit, denekkphase_visit, &
           xcorrphase_visit, linkdrive_visit
      REAL :: kx(0:nxcmplx-1),ky(0:nycmplx-1)

!C  Specify file name
      CALL get_iterstr(citer,ivisit,LEN(citer))

      filename = 'beta.'//citer//'.silo'

      ierr = DBCREATE(filename, 16, DB_CLOBBER, DB_LOCAL, "BETA_RUN", 8, DB_HDF5, dbfile)

!C  Write density and vorticity

      ierr = DBMKOPTLIST(4, optlistid)
      ierr = DBADDIOPT(optlistid, DBOPT_CYCLE, ncl)
      ierr = DBADDDOPT(optlistid, DBOPT_DTIME, t)
      ierr = DBADDCOPT(optlistid, DBOPT_XLABEL, "X", 1)
      ierr = DBADDCOPT(optlistid, DBOPT_YLABEL, "Y", 1)
       
!C    Initialize
      denekk_visit = 0.d0; phiekk_visit = 0.d0; denavgekk_visit = 0.d0
      phiekkphase_visit = 0.d0; denekkphase_visit = 0.d0
      xcorrphase_visit = 0.d0; linkdrive_visit = 0.d0
      tpos_visit = 0.d0

      kden = czero; kphi = czero; kvort = czero; kdenavg = czero; kden2 = czero
      klindrive = czero
      kx = 0.0; ky = 0.0

!C    Get real components
      kden(0:kxmax,0:kymax) = den(0:kxmax,0:kymax)
      kden(0:kxmax,nycmplx-kymax:nycmplx-1) = den(0:kxmax,-kymax:-1)
      kphi(0:kxmax,0:kymax) = phi(0:kxmax,0:kymax)
      kphi(0:kxmax,nycmplx-kymax:nycmplx-1) = phi(0:kxmax,-kymax:-1)
      kdenavg(0:kxmax,0:kymax) = denavg(0:kxmax,0:kymax)
      kdenavg(0:kxmax,nycmplx-kymax:nycmplx-1) = denavg(0:kxmax,-kymax:-1)
      klindrive(0:kxmax,0:kymax) = linkdrive(0:kxmax,0:kymax)
      klindrive(0:kxmax,nycmplx-kymax:nycmplx-1) = linkdrive(0:kxmax,-kymax:-1)

      kx(0:kxmax) = kxnorm*(/ (i,i=0,kxmax) /)
      ky(0:kymax) = kynorm*(/ (j,j=0,kymax) /)
      ky(nycmplx-kymax:nycmplx-1) = kynorm*(/ (j,j=-kymax,-1) /)

!!$OMP PARALLEL DO PRIVATE(ksq)
      DO j = 0, nycmplx-1
         DO i = 0, nxcmplx-1
            ksq = kx(i)*kx(i) + ky(j)*ky(j)
            kvort(i,j) = -1.0 * ksq * kphi(i,j)
            kden2(i,j) = -1.0 * ksq * kden(i,j)
         END DO
      END DO
!!$OMP END PARALLEL DO
       

!$      IF (l_omp .EQV. .TRUE.) THEN ! switch to OpenMP commands
!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads,plan2,kden,den_visit)
!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads,plan2,kphi,phi_visit)
!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads,plan2,kvort,vort_visit)
!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads,plan2,kden2,den2_visit)
!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads,plan2,kdenavg,denavg_visit)
!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads,plan2,klindrive,lindrive_visit)
!$      ELSE
         CALL rfftwnd_f77_one_complex_to_real(plan2,kden,den_visit)
         CALL rfftwnd_f77_one_complex_to_real(plan2,kphi,phi_visit)
         CALL rfftwnd_f77_one_complex_to_real(plan2,kvort,vort_visit)
         CALL rfftwnd_f77_one_complex_to_real(plan2,kden2,den2_visit)
         CALL rfftwnd_f77_one_complex_to_real(plan2,kdenavg,denavg_visit)
         CALL rfftwnd_f77_one_complex_to_real(plan2,klindrive,lindrive_visit)
!$      END IF

         IF (l_tracers_visit) CALL get_tpos_visit(tpos_visit,nxreal,nyreal,ncl,dt)

!C    Write to file
      ierr = DBPUTQM(dbfile, "spacemesh", 9, "x", 1, "y", 1, "z", 1, &
           x_visit, y_visit, DB_F77NULL, xdims, ndims, DB_DOUBLE, DB_COLLINEAR, optlistid, err)

      ierr = DBPUTQV1(dbfile, "density", 7, "spacemesh", 9, &
           den_visit, xdims, ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, err)

      ierr = DBPUTQV1(dbfile, "potential", 9, "spacemesh", 9, &
!                                         ^^ string length of name!!
           phi_visit, xdims, ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, err)

      ierr = DBPUTQV1(dbfile, "vorticity", 9, "spacemesh", 9, &
           vort_visit, xdims, ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, err)

      ierr = DBPUTQV1(dbfile, "dencurl", 7, "spacemesh", 9, &
           den2_visit, xdims, ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, err)

      ierr = DBPUTQV1(dbfile, "denavg", 6, "spacemesh", 9, &
           denavg_visit, xdims, ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, err)

      ierr = DBPUTQV1(dbfile, "lindrive", 8, "spacemesh", 9, &
           lindrive_visit, xdims, ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, err)

      ierr = DBPUTQV1(dbfile, "tracers", 7, "spacemesh", 9, &
           tpos_visit, xdims, ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, err)


!C  Write k-wave spectrum energy      

      ierr = DBFREEOPTLIST(optlistid)
      ierr = DBMKOPTLIST(2, optlistid)
      ierr = DBADDCOPT(optlistid, DBOPT_XLABEL, "KX", 2)
      ierr = DBADDCOPT(optlistid, DBOPT_YLABEL, "KY", 2)

      ierr = DBPUTQM(dbfile, "kmesh", 5, "kx", 2, "ky", 2, "z", &
           1, kx_visit, ky_visit, DB_F77NULL, kdims, ndims, DB_DOUBLE, DB_COLLINEAR, optlistid, err)

!!$OMP PARALLEL DO
      DO j = -kymax, kymax
         DO i = 0, kxmax
!            xk = kxnorm*i                                      !!! x-wavenumber
!            yk = kynorm*j                                      !!! y-wavenumber
!            ak = xk*xk + yk*yk                                 !!! Total k squared
            denekk_visit(i,j) = 2.d0*den(i,j)*CONJG(den(i,j))    !!! Energy in a mode
            phiekk_visit(i,j) = 2.d0*phi(i,j)*CONJG(phi(i,j))
            denavgekk_visit(i,j) = 2.d0*denavg(i,j)*CONJG(denavg(i,j))
            phiekkphase_visit(i,j) = ATAN2(REAL(phi(i,j)),AIMAG(phi(i,j))) !!! phi phase
            denekkphase_visit(i,j) = ATAN2(REAL(den(i,j)),AIMAG(den(i,j))) !!! den phase
            cxcorr_visit(i,j) = den(i,j)*CONJG(phi(i,j))       !!! phase cross-correlation between fields
            xcorrphase_visit(i,j) = ATAN2(REAL(cxcorr_visit(i,j)),AIMAG(cxcorr_visit(i,j)))
            linkdrive_visit(i,j) = 2.d0*linkdrive(i,j)*CONJG(linkdrive(i,j))  !!! linear drive terms for fluctuation fields
         END DO
      END DO
!!$OMP END PARALLEL DO
      denekk_visit(0,0) = denekk_visit(1,0)                        ! the (0,0) mode
      phiekk_visit(0,0) = phiekk_visit(1,0) 
      denavgekk_visit(0,0) = denavgekk_visit(1,0)

!      cxcorr_visit = cxcorr_visit/ABS(den*CONJG(phi))              ! normalize the cross correlation function

      ierr = DBPUTQV1(dbfile, "denkspectrum", 12, "kmesh", 5, &
           denekk_visit, kdims, ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, err)

      ierr = DBPUTQV1(dbfile, "phikspectrum", 12, "kmesh", 5, &
           phiekk_visit, kdims, ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, err)

      ierr = DBPUTQV1(dbfile, "denavgkspectrum", 15, "kmesh", 5, &
           denavgekk_visit, kdims, ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, err)

      ierr = DBPUTQV1(dbfile, "phikphase", 9, "kmesh", 5, &
           phiekkphase_visit, kdims, ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, err)

      ierr = DBPUTQV1(dbfile, "denkphase", 9, "kmesh", 5, &
           denekkphase_visit, kdims, ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, err)

!      ierr = DBPUTQV1(dbfile, "xcorr", 5, "kmesh", 5, &
!           cxcorr_visit, kdims, ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, err)

      ierr = DBPUTQV1(dbfile, "xcorrphase", 10, "kmesh", 5, &
           xcorrphase_visit, kdims, ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, err)

      ierr = DBPUTQV1(dbfile, "linkdrive", 9, "kmesh", 5, &
           linkdrive_visit, kdims, ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, err)

      ierr = DBFREEOPTLIST(optlistid)
      ierr = DBCLOSE(dbfile)

      ivisit = ivisit + 1

      WRITE(6,*) " SILO ==> "//filename," written at t = ", t

END SUBROUTINE WRITE_VISIT


SUBROUTINE get_tpos_visit(tposarr,nx,ny,ncl,dt)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  SUBROUTINE get_tpos_visit: construct an array with 1's and 0's
!c    for VISIT
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE tracers
      USE prog_flags, ONLY: nvisit
      IMPLICIT NONE
      
      INTEGER :: i,j,n
      INTEGER, INTENT(IN) :: nx, ny, ncl
      REAL, INTENT(IN) :: dt
      REAL, PARAMETER :: rfac = 1.25, rfac0 = 1.0      ! radial size of tracer point in SILO file
      REAL, DIMENSION(nx,ny), INTENT(OUT) :: tposarr
      INTEGER, DIMENSION(nmax,2,nvisit) :: tpx_tail, tpy_tail   ! temporary arrays
      REAL, DIMENSION(nmax,2,nvisit) :: px_tail_visit, py_tail_visit

      tposarr = 0.d0

      ! need to transfer coordinates to nearest "spacemesh" coordinates
      tpx_tail = INT(ANINT(px_tail*REAL(nx-1)/boxlen))   ! rounding numbers to the positions to nearest integer of mesh location
      tpy_tail = INT(ANINT(py_tail*REAL(ny-1)/boxlen))
      px_tail_visit = REAL(tpx_tail-1)/REAL(nx-1)  ! transfer back to fractions
      py_tail_visit = REAL(tpy_tail-1)/REAL(ny-1)

    
      DO i = 1, nvisit
!         DO j = 1, 2
            DO n = 1, nmax
!               tposarr(tpx_tail(n,j,i)+1,tpy_tail(n,j,i)+1) = 1.0    ! mark the locations where the tracers are
               tposarr(tpx_tail(n,1,i)+1,tpy_tail(n,1,i)+1) = 1.0 
            END DO
!         END DO
      END DO

      IF (MOD(ncl,INT(10/dt)) .EQ. 0) tposarr_old = 0 ! clear old tracers' trajectories
      tposarr = tposarr + tposarr_old
      tposarr_old = tposarr              ! save trajectories


END SUBROUTINE get_tpos_visit


SUBROUTINE write_polavg(denunit,phiunit,denavgunit,Lnunit,denavg1dunit)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  SUBROUTINE write_polavg: writes poloidal averaged density,
!C    potential, mean filed, and scale length
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE fft_variables
      USE variables, ONLY: denavg,den,phi,rLn,rdenavg1d,ncl,czero
      IMPLICIT NONE

      INTEGER :: i
      INTEGER, INTENT(in) :: denunit, phiunit, denavgunit, Lnunit, denavg1dunit
      COMPLEX, DIMENSION(0:nxcmplx-1) :: kdenavg,kden,kphi
      REAL, DIMENSION(nxreal) :: rdenavg,rden,rphi

      kdenavg = czero; kden = czero; kphi = czero

      ! extract the poloidal average, this corresponds to the ky = 0 modes
      kdenavg(0:kxmax) = denavg(0:kxmax,0)
      kden(0:kxmax) = den(0:kxmax,0)
      kphi(0:kxmax) = phi(0:kxmax,0)

      CALL rfftwnd_f77_one_complex_to_real(plan2_1d,kdenavg,rdenavg)
      CALL rfftwnd_f77_one_complex_to_real(plan2_1d,kden,rden)
      CALL rfftwnd_f77_one_complex_to_real(plan2_1d,kphi,rphi)

      WRITE(denunit,*) ncl, ( rden(i),i=1,nxreal )
      WRITE(phiunit,*) ncl, ( rphi(i),i=1,nxreal )
      WRITE(denavgunit,*) ncl, ( rdenavg(i),i=1,nxreal )
      WRITE(Lnunit,*) ncl, ( rLn(i),i=1,nxreal )
      WRITE(denavg1dunit,*) ncl, ( rdenavg1d(i),i=1,nxreal )

END SUBROUTINE write_polavg


SUBROUTINE get_iterstr(itstr,itnum,numd)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  SUBROUTINE get_iterstr: returns a string to the number of iteration
!C    with zero padding according to the number of digits, "numd"
!C    Currently, the number of places is set to 99.
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE

      INTEGER, INTENT(in) :: itnum
      INTEGER, INTENT(in) :: numd
      CHARACTER(LEN=numd), INTENT(out) :: itstr
      CHARACTER(LEN=2) :: numdstr

      IF ((numd .GE. 0) .AND. (numd .LT. 10)) THEN
         WRITE(numdstr,'(I1)') numd
      ELSE IF ((numd .GE. 10) .AND. (numd .LT. 100)) THEN
         WRITE(numdstr,'(I2)') numd
      ELSE
         WRITE(6,*) "Error in GET_ITERSTR: Number of digits is not in range."
      END IF

      WRITE(itstr,'(I0.'//TRIM(numdstr)//')') itnum

END SUBROUTINE get_iterstr


SUBROUTINE SECOND0(STIME)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  SUBROUTINE SECOND0(STIME): timing subroutine
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      IMPLICIT NONE

      REAL :: stime
      INTEGER, DIMENSION(8):: values
      CALL DATE_AND_TIME(VALUES=values)
      stime = 60.*(60.*values(5) + values(6)) + values(7)

END SUBROUTINE SECOND0



!C BEGIN TRACER SECTION ********************************************************************************************************************************************************************
!C subroutines in the MAIN level begin with the name "tracers_"

SUBROUTINE tracers_init
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C     subroutine tracers_init                RS: 06/29/05
!C
!C     initialize tracers positions at iteration NITER (modified)
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE tracers
      USE fileunits
      USE fft_variables
      USE variables, ONLY: t
      USE prog_flags, ONLY: lcont, nvisit
      IMPLICIT NONE
      
!      itrack = 0         ! counts the number entries written for tracer trajectories (so don't have to determine record length later on)
!     itrack has moved to "check_cont" subroutine
      nmax = mbox*nbox   ! number of particles
!      ixmax = 2*kxm1 - 2; jymax = 2*kym1 - 2

!      ALLOCATE( px0(nmax),py0(nmax),px(nmax),py(nmax) )
!      ALLOCATE( pvx0(nmax),pvx(nmax),pvy(nmax),pvy0(nmax) )
!      ALLOCATE( px0(mbox,nbox),py0(mbox,nbox),px(mbox,nbox),py(mbox,nbox) )  ! tracers ID are denoted by (i,j)
!      ALLOCATE( pvx0(mbox,nbox),pvx(mbox,nbox),pvy(mbox,nbox),pvy0(mbox,nbox) )
      ALLOCATE( px0(nmax,2),py0(nmax,2),px(nmax,2),py(nmax,2) )  ! tracers are paired
      ALLOCATE( px_tail(nmax,2,nvisit),py_tail(nmax,2,nvisit) )  ! paired tracers with histories
      ALLOCATE( tposarr_old(nxreal,nyreal) )                     ! marked positions for tracer output to VISIT
      ALLOCATE( pvx0(nmax,2),pvx(nmax,2),pvy(nmax,2),pvy0(nmax,2) )
      ALLOCATE( lypmax(nmax),pdistinit(nmax) ) ! tracking for maximum lypunov exponent
      ALLOCATE( deltax_old(nmax),deltay_old(nmax),timex_old(nmax),timey_old(nmax),xpos_old(nmax),ypos_old(nmax) )  ! for recording flights
      
      lypmax = 0.0; pdistinit = 0.0; lypcounter = 0
      tposarr_old = 0.d0

      IF (lcont .EQV. .FALSE.) THEN

         told = 0.0

         OPEN(unit=u_tracers,file='tposvel.bin',status='replace',form='unformatted')
         OPEN(unit=u_tflightsx,file='tflightsx.bin',status='replace',form='unformatted')
         OPEN(unit=u_tflightsy,file='tflightsy.bin',status='replace',form='unformatted')

!!$         OPEN(u_tposx,file='tposx.out',status='replace',form='formatted')
!!$         OPEN(u_tposy,file='tposy.out',status='replace',form='formatted')
!!$         OPEN(u_tvelx,file='tvelx.out',status='replace',form='formatted')
!!$         OPEN(u_tvely,file='tvely.out',status='replace',form='formatted')
!!$         OPEN(u_tflightsx,file='tflightsx.out',status='replace',form='formatted')
!!$         OPEN(u_tflightsy,file='tflightsy.out',status='replace',form='formatted')

      ELSE

         told = t

         OPEN(unit=u_tracers,file='tposvel.bin',status='old',position='append',form='unformatted')
         OPEN(unit=u_tflightsx,file='tflightsx.bin',status='old',position='append',form='unformatted')
         OPEN(unit=u_tflightsy,file='tflightsy.bin',status='old',position='append',form='unformatted')

!!$         OPEN(u_tposx,file='tposx.out',status='old',position='append',form='formatted')
!!$         OPEN(u_tposy,file='tposy.out',status='old',position='append',form='formatted')
!!$         OPEN(u_tvelx,file='tvelx.out',status='old',position='append',form='formatted')
!!$         OPEN(u_tvely,file='tvely.out',status='append',form='formatted')
!!$         OPEN(u_tflightsx,file='tflightsx.out',status='append',form='formatted')
!!$         OPEN(u_tflightsy,file='tflightsy.out',status='append',form='formatted')

      END IF

!      WRITE(u_tposx,212) nmax

      CALL pinit   ! initialize tracers' positions

END SUBROUTINE tracers_init



SUBROUTINE tracers_final
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C     subroutine tracers_final                RS: 06/29/05
!C
!C     close tracer files
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE tracers
      USE fileunits
      USE variables, ONLY: t
      IMPLICIT NONE

      INTEGER :: n

      ! write the lypunov exponents to file
      OPEN(u_tlypmax,file='tlypmax.out',status='replace',form='formatted')
      DO n = 1, nmax
         WRITE(u_tlypmax,'(I8,1x,1pe16.9,1x,1pe16.9)') n,lypmax(n)/REAL(lypcounter),pdistinit(n)
      END DO
      CLOSE(u_tlypmax)

      CALL pwritecont(51)   ! write continuation file

!!$      ! write lypunov information for continuation
!!$      OPEN(unit=78,file='tlypmax.cont',status='replace',form='formatted')
!!$      WRITE(78,1114) lypcounter
!!$      WRITE(78,1115) (pdistinit(n),n=1,nmax)
!!$      WRITE(78,1115) (lypmax(n),n=1,nmax)
!!$      CLOSE(78)
      
!!$1114  FORMAT(I15)
!!$1115  FORMAT(20(1pe16.9))

!!$      ! write final tracers' positions for continuation
!!$      OPEN(unit=71,file='tposx.cont',status='replace',form='formatted')
!!$      OPEN(unit=72,file='tposy.cont',status='replace',form='formatted')
!!$      OPEN(unit-73,file='tvelx.cont',status='replace',form='formatted')
!!$      OPEN(unit=74,file='tvely.cont',status='replace',form='formatted')
!!$      CALL pwrite(t,71,72,73,74)
!!$      CLOSE(71); CLOSE(72); CLOSE(73); CLOSE(74)

      DEALLOCATE(px0,py0,px,py,pvx,pvy,pvx0,pvy0)
      DEALLOCATE(px_tail,py_tail)
      DEALLOCATE(lypmax)
      DEALLOCATE(deltax_old,deltay_old,timex_old,timey_old)
!      CLOSE(u_tposx); CLOSE(u_tposy); CLOSE(u_tvelx); CLOSE(u_tvely)
      CLOSE(u_tracers)
      CLOSE(u_tflightsx); CLOSE(u_tflightsy)

END SUBROUTINE tracers_final


SUBROUTINE pinit 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  Particles                                                    10/16/95
!C
!C   All the particle subroutines.  Note that these 
!C   have been changed for the x-y initialization. 
!C
!C  subroutine PINIT                        
!C
!C  Initialize the particles according to NVPART
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE tracers
      USE variables, ONLY: irand, t
      USE fileunits
      USE prog_flags, ONLY: lcont
!      USE fft_variables

      IMPLICIT NONE

      LOGICAL :: lexist
      INTEGER :: i,n
      REAL, EXTERNAL :: ran2
      REAL, PARAMETER :: small = TINY(0.0d0)
      REAL, DIMENSION(nmax) :: dx0,dy0
      REAL :: tdump
      CHARACTER(LEN=128) :: tcontfile_str

      IF (lcont .EQV. .FALSE.) THEN

         SELECT CASE (nvpart)
         CASE (1)   ! Initializes tracers at random positions

            WRITE(6,*) '  Tracers equally spaced in a box. Number of tracers: ', nmax
            WRITE(u_inputscopy,*) '  Tracers equally spaced in a box. Number of tracers: ', nmax

            DO n = 1, nmax        ! loop over number of tracers
               px0(n,1) = ran2(irand+n)*boxlen
               py0(n,1) = ran2(irand+n*n)*boxlen
               px0(n,2) = px0(n,1) + smalldist*boxlen ! starting another tracer nearby, test tracers
               py0(n,2) = py0(n,1)             ! same y-position
            END DO

            px = px0; py = py0 ! assign initial positions

         CASE (2)   ! Read from input file
            
            WRITE(6,*) '  Tracers read in from "tcont.bin" file'

            ! a small error checking loop
            lexist = .FALSE.
            DO
               WRITE(6,'(/a)',advance='no') "Specify the tcont.bin file (enter 'q' to quit): "
               READ(5,'(a)') tcontfile_str
               IF (TRIM(ADJUSTL(tcontfile_str)) == "q") THEN
                  WRITE(6,'(//a/)') "User quits from specifying restart file. BETA code aborted."
                  STOP
               END IF
               INQUIRE(file=TRIM(ADJUSTL(tcontfile_str)), exist=lexist)
               IF (lexist .EQV. .TRUE.) EXIT ! break from loop 
            END DO
            
            WRITE(6,'(/a)', advance='no') "Specified tcont.bin file: "
            WRITE(6,*) TRIM(ADJUSTL(tcontfile_str))
            OPEN(unit=47, file=TRIM(ADJUSTL(tcontfile_str)), status='old', form='unformatted')
            READ(47) px0,py0
            CLOSE(47)
         
            px0(1:nmax,2) = px0(1:nmax,1) + smalldist*boxlen
            py0(1:nmax,2) = py0(1:nmax,1)
            px = px0; py = py0

         CASE DEFAULT

            STOP "WARNING: NVPART needs to defined properly!"

         END SELECT

         dx0 = px0(1:nmax,1) - px0(1:nmax,2); dy0 = py0(1:nmax,1) - py0(1:nmax,2)
         pdistinit = SQRT( dx0*dx0 + dy0*dy0 )

         ! flights initializations
         deltax_old = 0.0; deltay_old = 0.0
         timex_old = t; timey_old = t
         xpos_old = px0(1:nmax,1); ypos_old = py0(1:nmax,1)

      ELSE

         WRITE(6,*) "  Read in tracers' positions from continuation files."
         CALL preadcont(51)
         
      END IF

      pvx = 0.0; pvy = 0.0  ! no initial velocities
      px_tail = 0.0; py_tail = 0.0
      px_tail(1:nmax,1:2,1) = px(1:nmax,1:2)  ! add to time history
      py_tail(1:nmax,1:2,1) = py(1:nmax,1:2)

END SUBROUTINE pinit


SUBROUTINE padvance(phifield,dt)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  subroutine padvance
!C
!C   Advance the particles.
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE fft_variables
      USE tracers
      USE visit
      USE prog_flags, ONLY: nvisit
      IMPLICIT NONE

      INTEGER :: n, j
      REAL :: vxt, vyt
      REAL, INTENT(in) :: dt
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax), INTENT(in) :: phifield

      DO j = 1, 2         ! loop over pairs
         DO n = 1, nmax   ! loops over number of tracers
            
            vxt = 0.d0; vyt = 0.d0

            CALL pvelocity(vxt,vyt,px(n,j),py(n,j),phifield)             ! find velocity at tracer position
!            px(n,j) = px(n,j) + vxt*REAL(ixmax-1)*dt/REAL(nx_visit - 1)  ! normalize to visit's box
!            py(n,j) = py(n,j) + vyt*REAL(jymax-1)*dt/REAL(ny_visit - 1)
            px(n,j) = px(n,j) + vxt*boxlen*dt
            py(n,j) = py(n,j) + vyt*boxlen*dt
            pvx(n,j) = vxt; pvy(n,j) = vyt
         
            IF (j .EQ. 1) CALL pflights(px(n,1),py(n,1),n) ! do only for tracers, not the test tracers

         END DO
      END DO



END SUBROUTINE padvance
      
SUBROUTINE update_tracer_tail
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  SUBROUTINE update_tracer_tail
!C
!C  Updates tail of the tracers according to time step
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE tracers
      USE prog_flags, ONLY: nvisit
      IMPLICIT NONE

      INTEGER :: i,j,n
      REAL, DIMENSION(nmax,2) :: pxt, pyt      ! temporary arrays for normalized position to the box

      DO j = 1,2
         DO n = 1, nmax

            ! normalize to the box for VISIT output
            IF (px(n,j) .GT. 0.0) pxt(n,j) = MOD(px(n,j),boxlen)        
            IF (px(n,j) .LE. 0.0) pxt(n,j) = boxlen + MOD(px(n,j),boxlen)
            IF (py(n,j) .GT. 0.0) pyt(n,j) = MOD(py(n,j),boxlen)
            IF (py(n,j) .LE. 0.0) pyt(n,j) = boxlen + MOD(py(n,j),boxlen)

         END DO
      END DO

      ! update the time histories by shifting data by one index (the first index will store the most up-to-date tracers' positions)
      px_tail(1:nmax,1:2,2:nvisit) = px_tail(1:nmax,1:2,1:nvisit-1)
      py_tail(1:nmax,1:2,2:nvisit) = py_tail(1:nmax,1:2,1:nvisit-1)

      px_tail(1:nmax,1:2,1) = pxt(1:nmax,1:2)
      py_tail(1:nmax,1:2,1) = pyt(1:nmax,1:2)

END SUBROUTINE update_tracer_tail

SUBROUTINE pvelocity(vxt,vyt,xt,yt,phifield)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  subroutine pvelocity
!C
!C   Find velocity of at a particular spatial point.
!C   Accepts a position (x0,y0) and returns the velocity
!C     at that location (vx,vy) by summing over all Fourier modes.
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE fft_variables
      USE variables, ONLY: twopi,ic,kxnorm,kynorm
      USE tracers, ONLY: ixmax, jymax, boxlen
      IMPLICIT NONE

      INTEGER :: kx, ky
      REAL :: xk, yk, spx, spy
      REAL, INTENT(in) :: xt,yt
      REAL, INTENT(out) :: vxt,vyt
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax), INTENT(in) :: phifield
!      REAL, DIMENSION(0:kxmax) :: kx_vect
!      REAL, DIMENSION(-kymax:kymax) :: ky_vect

!C  Initialization
 
      vxt = 0.0; vyt = 0.0
!      kx_vect = kxnorm*(/ (kx,kx=0,kxmax) /); ky_vect = kynorm*(/ (ky,ky=-kymax,kymax) /)

!C   Normalize particle infinite domain to periodic domain of the field

!      IF(xt .GT. 0.0) spx = MOD(xt,REAL(ixmax-1))                    ! Note that SPX(Y) contains the value of PX(Y) truncated
!      IF(xt .LE. 0.0) spx = REAL(ixmax-1) + MOD(xt,REAL(ixmax-1))    ! onto the physical space.  (??)
!      IF(yt .GT. 0.0) spy = MOD(yt,REAL(jymax-1))
!      IF(yt .LE. 0.0) spy = REAL(jymax-1) + MOD(yt,REAL(jymax-1))

!  Not sure if this is implemented correctly
!      IF (xt .GT. 0.0) spx = MOD(xt,twopi/kxnorm)                    ! Note that SPX(Y) contains the value of PX(Y) truncated
!      IF (xt .LE. 0.0) spx = twopi/kxnorm + MOD(xt,twopi/kxnorm)    ! onto the physical space.  (??)
!      IF (yt .GT. 0.0) spy = MOD(yt,twopi/kynorm)
!      IF (yt .LE. 0.0) spy = twopi/kynorm + MOD(yt,twopi/kynorm)

      IF (xt .GT. 0.0) spx = MOD(xt,boxlen)                           ! modified implementation
      IF (xt .LE. 0.0) spx = boxlen + MOD(xt,boxlen)
      IF (yt .GT. 0.0) spy = MOD(yt,boxlen)
      IF (yt .LE. 0.0) spy = boxlen + MOD(yt,boxlen)
      spx = twopi*spx; spy = twopi*spy

!      spx = twopi*spx/REAL(ixmax-1)                                  ! Normalize position to a box of size 1
!      spy = twopi*spy/REAL(jymax-1)

!C  Calculate velocities by adding over the Fourier series

!!$OMP PARALLEL DO REDUCTION(+:vxt,vyt) PRIVATE(yk,xk)
      DO ky = -kymax, kymax
         yk = kynorm*REAL(ky)

         DO kx = 0, kxmax
            xk = kxnorm*REAL(kx)

            vxt = vxt - REAL(ic*yk*phifield(kx,ky)*EXP(ic*(xk*spx+yk*spy)))
            vyt = vyt + REAL(ic*xk*phifield(kx,ky)*EXP(ic*(xk*spx+yk*spy)))

!            vxt = vxt - REAL(ic*yk*phifield(kx,ky)*EXP(ic*(xk*spx+yk*spy)))
!            vyt = vyt + REAL(ic*xk*phifield(kx,ky)*EXP(ic*(xk*spx+yk*spy)))

         END DO
      END DO
!!$OMP END PARALLEL DO
!      vxt = 2.0d0*vxt/(nxreal*nyreal); vyt = 2.0d0*vyt/(nxreal*nyreal)                ! Two factor to account for symmetric part of spectrum
      vxt = 2.0d0*vxt; vyt = 2.0d0*vyt   ! omitting the division by the total number of grid points since the factor is simply a scaling factor

END SUBROUTINE pvelocity


SUBROUTINE pwrite(nentry,t,ftracers)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine PWRITE
!C
!C   Writes the particle positions at time t.
!C   Must also enter file units. This is a binary file.
!C      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE tracers
      IMPLICIT NONE

      INTEGER, INTENT(in) :: nentry
      REAL, INTENT(in) :: t
      INTEGER, INTENT(in) :: ftracers

      WRITE(ftracers) nentry,t,px,py,pvx,pvy

!!$      WRITE(fposx,1114,advance='no') t
!!$      WRITE(fposy,1114,advance='no') t
!!$      WRITE(fvelx,1114,advance='no') t
!!$      WRITE(fvely,1114,advance='no') t
!!$
!!$      DO n = 1, nmax-1
!!$         WRITE(fposx,1115,advance='no') px(n,1),px(n,2)    ! write x-position pairs
!!$         WRITE(fposy,1115,advance='no') py(n,1),py(n,2)    ! write y-position pairs
!!$         WRITE(fvelx,1115,advance='no') pvx(n,1),pvx(n,2)  ! write x-velocity pairs
!!$         WRITE(fvely,1115,advance='no') pvy(n,1),pvy(n,2)  ! write y-velocity pairs
!!$      END DO
!!$
!!$      WRITE(fposx,1115) px(nmax,1),px(nmax,2) 
!!$      WRITE(fposy,1115) py(nmax,1),py(nmax,2)
!!$      WRITE(fvelx,1115) pvx(nmax,1),pvx(nmax,2)
!!$      WRITE(fvely,1115) pvy(nmax,1),pvy(nmax,2)
!!$
!!$1114  FORMAT(1pe16.9)
!!$1115  FORMAT(1x,1pe16.9,1x,1pe16.9)

END SUBROUTINE pwrite

SUBROUTINE pwritecont(funit)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine pwritecont
!C
!C   Writes continuation parameters. This is a binary file.
!C      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE tracers
      IMPLICIT NONE

      INTEGER, INTENT(in) :: funit
      INTEGER :: n

      OPEN(unit=funit,file="tcont.bin",status="replace",form="unformatted")
      WRITE(funit) px0,py0,px,py                       ! write particle positions
      WRITE(funit) lypcounter, pdistinit, lypmax       ! write lypunov information
      WRITE(funit) xpos_old, ypos_old, deltax_old, deltay_old, timex_old, timey_old       ! write flights information
      CLOSE(funit)

END SUBROUTINE pwritecont

SUBROUTINE preadcont(funit)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine preadcont
!C
!C   Reads from continuation tracer file. This is a binary file.
!C      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE tracers
      IMPLICIT NONE

      INTEGER, INTENT(in) :: funit
      INTEGER :: n

      OPEN(unit=funit,file="tcont.bin",status="old",form="unformatted")
      READ(funit) px0,py0,px,py                       ! read particle positions
      READ(funit) lypcounter, pdistinit, lypmax       ! read lypunov information
      READ(funit) xpos_old, ypos_old, deltax_old, deltay_old, timex_old, timey_old       ! read flights information
      CLOSE(funit)


END SUBROUTINE preadcont

SUBROUTINE plypmax
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine plypmax
!C
!C   Calculates the maximum lypunov exponent
!C     and resets the tracers' position
!C      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE tracers
      USE prog_flags, ONLY: l_tracers_renew
      IMPLICIT NONE

      INTEGER :: n
      REAL, DIMENSION(nmax) :: dx,dy

      dx = px(1:nmax,1) - px(1:nmax,2)   ! difference in positions for each tracer pair
      dy = py(1:nmax,1) - py(1:nmax,2)

      lypmax = lypmax + LOG( SQRT( dx*dx + dy*dy ) / pdistinit ) ! track sum for lypunov exponent, this is an elementwise operation

      IF (l_tracers_renew) THEN
         px(1:nmax,2) = px(1:nmax,1) + smalldist                    ! reset the paired test tracer
         py(1:nmax,2) = py(1:nmax,1)
      END IF

      lypcounter = lypcounter + 1
!      WRITE(6,*) "lypmax = ", lypmax," lypcounter = ", lypcounter

END SUBROUTINE plypmax


SUBROUTINE pflights(xpos,ypos,nid)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine plypmax
!C
!C   Determines the flights and flight times.
!C   Records to file if direction changes.
!C      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE variables, ONLY: t
      USE tracers
      USE fileunits
      IMPLICIT NONE

      INTEGER, INTENT(in) :: nid      ! tracer ID
      REAL, INTENT(in) :: xpos, ypos  ! tracer's positions
      REAL :: deltax, deltay

      deltax = xpos - xpos_old(nid) ! xpos_old is a global variable specified in the TRACERS module
      deltay = ypos - ypos_old(nid)

      ! if directions change
      IF (deltax*deltax_old(nid) .LT. 0.0) THEN     ! x-component
!         WRITE(u_tflightsx,1115) nid,deltax,t-timex_old(nid)
         WRITE(u_tflightsx) nid,deltax_old(nid),t-timex_old(nid)
         xpos_old(nid) = xpos
         deltax_old(nid) = deltax
         timex_old(nid) = t
         iflightsx = iflightsx + 1

      ELSE
         deltax_old(nid) = deltax_old(nid) + deltax

      END IF
      IF (deltay*deltay_old(nid) .LT. 0.0) THEN     ! y-component
!         WRITE(u_tflightsy,1115) nid,deltay,t-timey_old(nid)
         WRITE(u_tflightsy) nid,deltay,t-timey_old(nid)
         ypos_old(nid) = ypos
         deltay_old(nid) = deltay
         timey_old(nid) = t
         iflightsy = iflightsy + 1

      ELSE
         deltay_old(nid) = deltay_old(nid) + deltay

      END IF
      
!!$1115  FORMAT(I10,2(1x,1pe16.9))

END SUBROUTINE pflights

!C END TRACER SECTION ********************************************************************************************************************************************************************


SUBROUTINE make_namelist
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine make_namelist
!C
!C   Create a namelist for postprocessing
!C      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE fft_variables, ONLY: kxmax, kymax, nxreal, nyreal
      USE variables, ONLY: t, dt, ncl
      USE tracers, ONLY: nmax, lypmaxdtstep, itrack, iflightsx, iflightsy
      USE visit, ONLY: ivisit
      IMPLICIT NONE

      OPEN(UNIT = 602, FILE = "params_out", STATUS = "replace")
      
      ! for continuing the simulation
      WRITE(602,*) "&CONTP"
      WRITE(602,*) "TIME0 = ", t, " ,"
      WRITE(602,*) "DT0 = ", dt, " ,"
      WRITE(602,*) "NCL0 = ", ncl, " ,"
      WRITE(602,*) "IVISIT0 = ", ivisit, " /"
      WRITE(602,*)

      ! for tracers
      WRITE(602,*) "&TRACERP"
      WRITE(602,*) "NMAX0 = ", nmax, " ,"
      WRITE(602,*) "ITRACK0 = ", itrack, " ,"
      WRITE(602,*) "IFLIGHTSX0 = ", iflightsx, " ,"
      WRITE(602,*) "IFLIGHTSY0 = ", iflightsy, " ,"
      WRITE(602,*) "LYPMAXDTSTEP0 = ", lypmaxdtstep, " /"
      WRITE(602,*)

      ! extract number of modes
      WRITE(602,*) "&FFTP"
      WRITE(602,*) "KXMAX = ", kxmax, " ,"
      WRITE(602,*) "KYMAX = ", kymax, " ,"
      WRITE(602,*) "NXREAL = ", nxreal, " ,"
      WRITE(602,*) "NYREAL = ", nyreal, " /"

      CLOSE(602)

END SUBROUTINE make_namelist


SUBROUTINE check_cont(lcont)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   Subroutine: check_cont (check for continuation)
!C
!C   Checks if run is a continuation or a new run by looking for file "params_out"
!C   If this is a continuation, then load some counters from previous run,
!C     and open files with "APPEND" option.
!C   If not, open files with "REPLACE".
!C   Returns logical parameter "lcont" to denote status of continuation.
!C      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE precmod
      USE variables, ONLY: t,ncl,ncl_in,denekold,denenold,phiekold,phiek1old,phienold
      USE tracers, ONLY: itrack,iflightsx,iflightsy
      USE fileunits
      USE visit, ONLY: ivisit
      IMPLICIT NONE

      LOGICAL, INTENT(out) :: lcont

      ! temporary variables for read in
      REAL :: time0,dt0
      INTEGER :: ncl0,ivisit0,nmax0,itrack0,iflightsx0,iflightsy0,lypmaxdtstep0

      ! read in NAMELIST
      NAMELIST /contp/ time0,dt0,ncl0,ivisit0
      NAMELIST /tracerp/ nmax0,itrack0,iflightsx0,iflightsy0,lypmaxdtstep0

      ! initialize energy variables
      denekold = 0.0; denenold = 0.0
      phiekold = 0.0; phiek1old = 0.0; phienold = 0.0

      INQUIRE(file='params_out',exist=lcont)
      IF (lcont .EQV. .FALSE.) THEN   ! open new files if this is a fresh run

!C  OPEN files for READ and WRITE --------------------------------------
         OPEN(u_energy,file='energy',status='replace')
         WRITE(u_energy,100)     ! writes energy file header
100      FORMAT(4x,'ncl',10x,'t',10x,'denek',10x,'denen',10x,'dengamek',10x,'dengamen',10x, &
              'phiek',10x,'phiek1',10x,'phien',10x,'phigamek',10x,'phigamen',10x,'denavgek',10x,'ek')
         OPEN(u_inputscopy,file='inputs0_copy',status='replace')

         ! mean field files
         OPEN(u_denpolavg,file='denpolavg',status='replace')
         OPEN(u_phipolavg,file='phipolavg',status='replace')
         OPEN(u_denavgpolavg,file='denavgpolavg',status='replace')
         OPEN(u_Ln,file='Lnpolavg',status='replace')
         OPEN(u_denavg1d,file='denavg1d',status='replace')

         ! initialize variables
         t = 0.0; ncl = 0; ncl_in = 0; itrack = 0; ivisit = 0
         iflightsx = 0; iflightsy = 0

      ELSE

         ! read continuation parameters
         OPEN(21,file='params_out',status='old')
         READ(21,contp)
         READ(21,tracerp)
         CLOSE(21)

         t = time0; ncl = 0; ncl_in = ncl0; itrack = itrack0; ivisit = ivisit0
         iflightsx = iflightsx0; iflightsy = iflightsy0

         WRITE(6,'(/a)') " ** This is a continuation run. **"
         WRITE(6,*) "time = ", t
         WRITE(6,*) "ncl_in = ", ncl_in, " itrack = ", itrack, " ivisit = ", ivisit

         OPEN(u_energy,file='energy',status='old',position='append')
         
         OPEN(u_denpolavg,file='denpolavg',status='old',position='append')
         OPEN(u_phipolavg,file='phipolavg',status='old',position='append')
         OPEN(u_denavgpolavg,file='denavgpolavg',status='old',position='append')
         OPEN(u_Ln,file='Lnpolavg',status='old',position='append')
         OPEN(u_denavg1d,file='denavg1d',status='old',position='append')

      END IF

!C  End OPEN files -----------------------------------------------------

!C      OPEN(9,file='run.log',status='unknown')      ! <-- log files will be "tee" from script
      OPEN(95,file='dvodpk_log',status='replace')   ! <-- to track solver's performance

END SUBROUTINE check_cont

SUBROUTINE que_hora_es(file_unit)
! ----------------------------------------------------------------------------------------
! Subroutine to display time and date
! -- file_unit = 6 for "stdout" (default screen display)
! -- file_unit = 0 for "stderr" (error display)
! ----------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER :: file_unit
      CHARACTER(LEN=10) :: timew,datew,zone_h

      CALL DATE_AND_TIME(datew,timew,zone_h)
      WRITE(UNIT=file_unit, FMT='(/" time = ",a2,":",a2,":",a2,"       date = ",a2,"/",a2,"/",a4/)') &
           timew(1:2),timew(3:4),timew(5:6),datew(7:8),datew(5:6),datew(1:4)

END SUBROUTINE que_hora_es


SUBROUTINE check_meanfield(uin,uout)
! ----------------------------------------------------------------------------------------
! before using denavg, need to check if the spatial vallues are real
! hence, inverse transform first, check values, then transform back.
! implementing a floor
! ----------------------------------------------------------------------------------------
      USE precmod
      USE fft_variables
      USE variables, ONLY: czero
      USE omp_variables
      USE prog_flags, ONLY: l_omp
      IMPLICIT NONE

      INTEGER :: i,j
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax), INTENT(IN) :: uin
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax), INTENT(OUT) :: uout
      COMPLEX, DIMENSION(0:nxcmplx-1,0:nycmplx-1) :: cuin,cuout
      REAL, DIMENSION(nxreal,nyreal) :: ru
      REAL :: rsmallval
      LOGICAL :: lval

      ! initialize
      cuin = czero; ru = 0.d0; uout = czero
      rsmallval = TINY(0.d0)
!      rsmallval = 1.0E-8
      lval = .FALSE.

      ! shuffle data into extended array
      cuin(0:kxmax,0:kymax) = uin(0:kxmax,0:kymax)
      cuin(0:kxmax,nycmplx-kymax:nycmplx-1) = uin(0:kxmax,-kymax:-1)

!$      IF (l_omp .EQV. .TRUE.) THEN
!$         CALL rfftwnd_f77_threads_one_complex_to_real(nthreads, plan2, cuin, ru)
!$      ELSE
         CALL rfftwnd_f77_one_complex_to_real(plan2, cuin, ru)
!$      END IF
      DO j = 1, nyreal
         DO i = 1, nxreal
            IF (ru(i,j) .LT. 0.d0) THEN
               ru(i,j) = rsmallval
               WRITE(*,'(/a)') "Unphysical mean field. Value changed."
               WRITE(*,*) "Value changed at: i = ",i, " j = ",j
               lval = .TRUE.
            END IF
         END DO
      END DO
      cuout = czero
!$      IF (l_omp .EQV. .TRUE.) THEN
!$         CALL rfftwnd_f77_threads_one_real_to_complex(nthreads,plan1,ru,cuout)
!$      ELSE
         CALL rfftwnd_f77_one_real_to_complex(plan1, ru, cuout)
!$      END IF

      ! unshuffle data
      uout(0:kxmax,0:kymax) = cuout(0:kxmax,0:kymax)
      uout(0:kxmax,-kymax:-1) = cuout(0:kxmax,nycmplx-kymax:nycmplx-1)

      uout = uout/REAL(nxreal*nyreal)

      ! catch errors
      IF ( (lval .EQV. .TRUE.) .AND. &
           (SUM(ABS(uout-uin)) .EQ. 0.0 ) ) THEN

         CALL FINALIZE
         STOP "Values changed, but fields have same value!"

      END IF

END SUBROUTINE check_meanfield
