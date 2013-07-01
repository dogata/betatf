!
!
! --> BETA two field MPI verson <--
!
! -- port over the structure of the one-field MPI code
! -- data storge also changed in this version to coincide with FFT requirement
! -- this version uses CVODE instead VODPK or DVODPK
! -- language is updated to F90
! -- check for SILO.INC file format also (might still be in F77 format)
! -- This is MPI version. Data structure is organized differently than serial version.
!    The extended array is used instead, which means there are zero modes in some processors.
! -- Referencing of array elements is different from the serial version, since the extended arrays are used.
!    For instance, the density arrary is defined locally as den(0:nxcmplx-1,yinit:yend).
!      nxcmplx ==> number of complex elements in kx-direction covering [0,1,...,kxmax,kxmax+1,...,nxcmplx-2,nxcmplx-1]
!      yinit   ==> locally starting index in ky-direction
!      yend    ==> locally ending index in ky-direction
! -- Node in the middle with modes that straddles positive and negative values (i.e. modes nxcmplx/2-1 and -nxcmplx/2)
!    is assigned as the node to write files. It is designated as "blacksheep".
! -- CVODE requires INTEGER*8 for arguments except for err code, which requires INTEGER*4
! -- Array size must be defined as divided through FFTW and cannot be specified externally as in the serial version.
! -- ** All processors must pass through CVODE ** So, two group communicators are not used.




MODULE fft_variables
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  MODULE fft_variables: accounts for variables required for FFT (global)
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
!           nn = 6*kxm2*kym1                               ! MPI code stores the full half

!      INTEGER, PARAMETER :: nn_serial = 6*(kym1*kxmax + kym2)
      INTEGER, PARAMETER :: nn_serial = 6*kxm2*kym1

      INTEGER, PARAMETER :: nfield = 3                   ! number of fields to solve

      INTEGER, PARAMETER :: &                            ! parameters for FFT - accounting for aliasing (require 3*(kxmax + 1) points)
!           nxreal = 128, nyreal = 128, &
           nxreal = 256, nyreal = 256, &                 ! number of points in REAL-space (number of modes: 625 to use radix 5), 256 for 41 modes
!           nxreal = 512, nyreal = 512, &                ! for 128x128 modes
!           nxreal = 729, nyreal = 729, &                ! using radix 3 at N = 729 for kxmax = kymax = 242 modes
!           nxreal=585,nyreal=585, &                     ! dealiasing requires 193*3 = 579 modes, radix 5 is used for optimized size
           nxcmplx = nxreal/2 + 1, nycmplx = nyreal      ! number of points in k-space (this is twice the one-field version)           

!C constant definitions needed by FFTW ------------------------------------------------
      INTEGER, PARAMETER :: FFTW_FORWARD = 1,FFTW_BACKWARD = -1
      INTEGER, PARAMETER :: FFTW_REAL_TO_COMPLEX = -1,FFTW_COMPLEX_TO_REAL = 1
      INTEGER, PARAMETER :: FFTW_ESTIMATE = 0,FFTW_MEASURE = 1
      INTEGER, PARAMETER :: FFTW_OUT_OF_PLACE = 0,FFTW_IN_PLACE = 8,FFTW_USE_WISDOM = 16
      INTEGER, PARAMETER :: FFTW_THREADSAFE = 128
      INTEGER, PARAMETER :: FFTW_TRANSPOSED_ORDER = 1,FFTW_NORMAL_ORDER = 0
!C ------------------------------------------------------------------------------------

      INTEGER(kind=8) :: planb, planf, planb_serial, planf_serial
      
END MODULE fft_variables


MODULE variables
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  MODULE variables: accounts for main variables in code
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE

!C
!C      Constants
!C
      REAL, PARAMETER :: pi = 4.d0*ATAN(1.d0)
      REAL, PARAMETER :: twopi = 2.d0*pi
      COMPLEX, PARAMETER :: ic = CMPLX(0.d0,1.d0), czero = CMPLX(0.d0,0.d0)
      REAL :: epsi, rhos, emass, nue, rtepsi, rhos3
 
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
!C      denekmax = maximum energy
!C      denekmin = minimum energy
!C      denenmax = maximum enstrophy
!C      denenmin = minimum enstrophy 

      REAL :: kxnorm,kynorm,v1,gama1,gama2,gama3,amp
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

      INTEGER :: ku,kl,kg1sq,kg3sq,kd1sq,kd3sq,kinit
    
!C	den 	= the array which contains the density field in k-space.
!C      phi     = array contains the potential field in k-space
!C      denavg  = density mean field (x = radial direction, y = poloidal direction), d<n>/dy = 0 only ky = 0 modes are valid
!C      phi0    = this field is for an external field
!C      denout  = array with all modes, but (0,0) mode is centered in the array, index runs from [-kxmax -> kxmax] and [-kymax, kymax]
!C      phitr   = phi-field for tracer propagation

      COMPLEX, DIMENSION(:,:), ALLOCATABLE :: den, phi, phiout, denout
      COMPLEX, DIMENSION(:,:), ALLOCATABLE :: denavg, denavgout, phi0, phitr

!C      phi0amp = amplitude of an external flow
      REAL :: phi0amp 

!C      rLn       = length scale stored to write to file
      
      REAL, DIMENSION(:), ALLOCATABLE :: rLn

!C       philinarr     = array of linear terms
!C       damparr       = array for damping terms for potential field
!C       dampdenavgarr = array for damping terms on the mean field
!C       drivearr      = array of drive terms
!C       densource     = array for source terms
!C       rdensource    = real-space array for source term

      COMPLEX, DIMENSION(:,:), ALLOCATABLE :: philinarr, damplinarr, drivelinarr, &
           dampdenavglinarr, philindrivearr
      COMPLEX, DIMENSION(:,:), ALLOCATABLE :: linkdrive        ! stored array for output of the drive for fluctuation fields
      COMPLEX, DIMENSION(:,:), ALLOCATABLE :: densource
      REAL, DIMENSION(:,:), ALLOCATABLE :: rdensource

!C	ntstep	= number of time steps
!C	irand	= number to feed random number generator

      INTEGER :: ntstep
      INTEGER(kind=8):: irand  ! kind=8 is important!!

!C    Theses terms pertain to the source and the mean field.
!C       samp    = source amplitude
!C       sscale  = source horizontal scale
      REAL :: samp, sscale

!C       tola the tolerance vector passed to VODPK matches element wise 
!C       with the Y(nn) vector

      REAL,DIMENSION(:), ALLOCATABLE :: tola, tola_serial

END MODULE variables


MODULE prog_flags
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  MODULE prog_flags: accounts for all the code flags
!C     (e.g. fft routine, save big files, etc.)
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE

!C       lp      = dummy constant (lp is a control on Pol nonlinear terms) (pol drift NL or NS)
!C       le      = dummy constant (le is a control on ExB nonlinear terms) (ExB NL)
!C       lm      = dummy constant (lm is a control on mean field convection)
      REAL :: lp,le,lm

!C	lpout	= turn on(=1) or off(=0) nonlinearity at high k's
!C	lpmid	= turn on(=1) or off(=0) nonlinearity at mid k's
      INTEGER :: lpout,lpmid

!C	init	= flag for which initial conditions to use
      INTEGER :: init

!C	nwrt	= how often to write ENERGY files
!C	nwrs	= how often to write RESTART files for spectrum and particles
!C      nvisit  = how often to write to silo files - SD
!C      nfft    = which fft subroutine to use. fftw if nfft=1, matfft otherwise (matfft has been completely disabled!)
      INTEGER :: nwrt,nwrs,nvisit,nfft

!C	To save or not to save the BIG files: stream(3), vorticity(34), 
!C       ensemble(12), energ_vs_ak(14), energyspec(11)
      LOGICAL :: l_save

      CHARACTER(LEN=5):: BETA_version = 'x.x0'

      LOGICAL:: l_tracers, l_tracers_renew                              ! L_TRACERS = .T. uses tracers, L_TRACERS_RENEW = reinit tracers every NPT iter

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

      INTEGER(kind=8), DIMENSION(21):: iout     ! <-- notice that this is KIND=8!!
      REAL, DIMENSION(6):: rout
      INTEGER :: ipar
      REAL :: rpar

END MODULE cvode_variables


MODULE visit
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  MODULE prog_flags: accounts for ViSit outputs
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE FFT_variables, ONLY: kxm1, kym1, kxm2, kxmax, kymax
      IMPLICIT NONE

      INCLUDE "silo_f90.inc"

      INTEGER :: nx_visit, ny_visit, nkx_visit, nky_visit, ivisit
      INTEGER :: ndims 
      INTEGER, DIMENSION(:), ALLOCATABLE:: xdims, kdims
      REAL, DIMENSION(:), ALLOCATABLE:: x_visit, y_visit, kx_visit, ky_visit


END MODULE visit


MODULE file_units
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  MODULE file_units: define external file units
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

END MODULE file_units



MODULE convolve_module
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  MODULE convolve_module: declares pointer values to use
!C                          in CONVOLVE subroutine
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE fft_variables
      IMPLICIT NONE

!C  To compensate for de-aliasing, dimensions are doubled for kx and ky dimensions.
!C  although "2/3" rule applies to compensate for de-aliasing, FFTW transforms efficiently for multiples of "2";
!C  hence, this is motivation for doubling the array size. Need to only transform half on modes, since the
!C  others are just conjugates. The modes must be placed with zero modes on the boundaries.
!C

!C  Declare arrays for FFT data. Notice that these are 1D arrays as required by FFTW.

      COMPLEX, DIMENSION(:), ALLOCATABLE :: &
           cvex_flat,cvey_flat, &          ! components of velocity
           cnky_flat,cnkx_flat, &          ! components of gradients of dn/dy
           cwx_flat,cwy_flat, &            ! vorticity gradient
           cnavgkx_flat, cnavgky_flat, &   ! components of gradients of mean field
           cnavgkx1d_flat, cnavg1d_flat, & ! mean field quantities for gradient scale length
           cnavg_flat                      ! local mean field quantity

! FFTW MPI transforms are in-place transforms, so extra arrays are not required
      REAL, DIMENSION(:), ALLOCATABLE :: Pol_flat, ExB_flat, denavgnl_flat, drivenl_flat

END MODULE convolve_module


MODULE mpi_variables
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  MODULE mpi_variables: declares values for MPI
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE

      LOGICAL :: l_ihavedata                                      ! specifies if processor has local data to calculate
      INTEGER, PARAMETER :: masternode = 0                        ! masternode is set to 0 (this could be changed to any node)
      INTEGER :: blacksheep                                       ! node that straddles negative and positive (comes from extended array for dealiasing)
      INTEGER :: npes,pid,mpierr                                  ! common MPI varaibles
      INTEGER :: nylocal,yinit,yend,nxlocal_after_transform, &    ! parameters extracted from FFTW's LOCAL_SIZES routine
           xinit_after_transform,local_size  
      INTEGER :: kylmin,kylmax, &                                 ! bounding modes for local data
           nn_local, &                                            ! total local modes to push through ODE solver
           kylmin_solv,kylmax_solv, &                             ! bounding modes for ODE solver
           kymodes, &                                             ! number of non-zero ky modes
           kylmin_conjg,kylmax_conjg, &                           ! bounding modes for conjugate modes
           ncmplxfft, &                                           ! total number of CMPLX elements for FFTs
           nytransfer                                             ! number of elements to transfer, should be at least the largest nylocal
      INTEGER, DIMENSION(2) :: localindex                         ! stores local values {nylocal,yinit}
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: localindex_dir      ! directory of local indices to be stored on master node and "blacksheep" node
      INTEGER, DIMENSION(:), ALLOCATABLE :: kxext, kyext, yext    ! stored indices to map between different reference schemes
      INTEGER :: mpi_comm_solver, mpi_comm_nosolver               ! MPI group communicators for two processor groups
      INTEGER, DIMENSION(:), ALLOCATABLE :: grp_solv              ! stored indices for number of processors in the solver group and their corresponding indices
      REAL, DIMENSION(:), ALLOCATABLE :: workarr                  ! Work array for FFTs

    CONTAINS

      LOGICAL FUNCTION estoy_bien(kmode,kmin,kmax)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  LOGICAL FUNCTION estoy_bien: returns .TRUE. if kmode is in the range
!C    [kmin, kmax]
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        IMPLICIT NONE

        INTEGER, INTENT(in) :: kmode, kmin, kmax

        IF ((kmode .GE. kmin) .AND. (kmode .LE. kmax)) THEN
           estoy_bien = .TRUE.
        ELSE
           estoy_bien = .FALSE.
        END IF

      END FUNCTION estoy_bien

END MODULE mpi_variables


MODULE tracers
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  Tracer parameters
!C
!C     Total number of particles = nmax
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE

      INTEGER :: nmax, nmax_local, nvpart
      INTEGER :: itrack                                                           ! ** "itrack" counts every time iteration minus the first iteration **
      INTEGER, DIMENSION(:), ALLOCATABLE :: nmax_index                            ! index for number of tracers per processor
      REAL, PARAMETER :: boxlen = 1                                               ! box length
      REAL, DIMENSION(:,:), ALLOCATABLE :: px0, py0, px, pvx0, pvy0, py, pvx, pvy ! position and velocities for pairs
      REAL, DIMENSION(:), ALLOCATABLE :: pdistinit                                ! initial tracers' distance
      REAL :: told                                                                ! stores time of previous step, used to advance tracers
      REAL :: smalldist                                                           ! this value is specified through the input file

END MODULE tracers

PROGRAM BETA
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  PROGRAM BETA: main program
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE mpi
      USE prog_flags
      USE variables, ONLY:ntstep,ncl,t,ncl_in
      USE fft_variables, ONLY: nfield
      USE mpi_variables
      IMPLICIT NONE

      INTEGER :: i,j
      REAL :: start_time, end_time
      REAL :: time_init, time_fin

!C  Banner variables
      CHARACTER(LEN=50), PARAMETER :: banner = 'THIS IS BETA . VERSION: '
      CHARACTER(LEN=30) :: computer
      INTEGER :: imon
      CHARACTER(LEN=10) :: date0, time0, zone0
      CHARACTER(LEN=3), DIMENSION(12), PARAMETER :: &
           months = (/'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'/)
      CHARACTER(LEN=40) :: dateloc


!     >>>> Start MPI <<<<
      CALL MPI_INIT(mpierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,npes,mpierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,pid,mpierr)


!C >>>>>> Disable this before submit to PACMAN! <<<<<<<<<<<<<<<
!      WRITE(*,*) "Hit enter to start BETA!"
!      READ(*,*)
!C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      CALL check_cont(lcont)       ! check for continuation file


!C  Initialize...
      CALL initial
      CALL GATHER_DATA(blacksheep,0)  ! collects data on "blacksheep" node

      IF (pid .EQ. blacksheep) THEN

         computer = ' Darwin OSX'

!C  Print banner -------------------------------------------------------
         WRITE(6,'(//a)') ' ========================================================================================='
         CALL DATE_AND_TIME(date0,time0,zone0)
         READ(date0(5:6),'(i2)') imon
         WRITE(dateloc,100) months(imon),date0(7:8),date0(1:4),time0(1:2),time0(3:4),time0(5:6)
100      FORMAT(' DATE = ',a3,' ',a2,',',a4,' ',' TIME= ',2(a2,':'),a2)
         WRITE(6,'(1x,2a,/,1x,2a)') banner, BETA_version, computer, dateloc
         WRITE(6,'(a/)') ' ========================================================================================='
         WRITE(6,'(/4x,a,4x/)') "*** MPI VERSION ***"
!C  END printing banner ------------------------------------------------

         WRITE(10,'(1x,2a,/,1x,2a)') banner, BETA_version, computer, dateloc

         CALL INIT_VISIT   ! Declare SILO files for ViSit
         CALL WROUT        ! Write output files

      END IF


!C  Iterate over time steps...
      DO ncl = ncl_in+1, ntstep+ncl_in

         CALL timead                                  ! advance in time
         
         IF ( MOD(ncl,nvisit) .EQ. 0.0 ) THEN          
            CALL GATHER_DATA(blacksheep,0)             ! collect data on "blacksheep" node
            IF (pid .EQ. blacksheep) CALL WROUT        ! write output files
         END IF

         IF (l_tracers) CALL tracers_write(t)

      END DO
     
      ncl = ncl - 1  ! counter is automatically advanced at end of DO loop

!C  Write "restart" file
      CALL GATHER_DATA(blacksheep,0)               ! collects data on "blacksheep" node
      IF (pid .EQ. blacksheep) CALL writeinit      ! writes last "restart" file

!C  Write continuation parameter files
      CALL make_namelist

!C  Terminate all processes
      CALL FINALIZE
      CALL MPI_FINALIZE(mpierr)

END PROGRAM BETA



SUBROUTINE initial
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                                                                      C
!C   subroutine initial                                                 C
!C                                                                      C
!C   Calls the initializing subroutines                                 C
!C   and writes the initail energy spectrum.                            C
!C                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE prog_flags
      USE variables
      USE fft_variables
      USE mpi_variables
      USE convolve_module
      IMPLICIT NONE

      INTEGER :: itime, ierr

      ! initialize energy variables
      denekold = 0.0; denenold = 0.0
      phiekold = 0.0; phiek1old = 0.0; phienold = 0.0

      IF (pid .EQ. masternode) THEN
         WRITE(6,*) "Precision of REAL type is: ", PRECISION(pi)  ! print out machine precision in decimal places
         WRITE(6,*) "Range of REAL type is: ", RANGE(pi)
         WRITE(6,*) "Maximum exponent of REAL type is: ", MAXEXPONENT(pi)
         WRITE(6,*) "Minimum exponent of REAL type is: ", MINEXPONENT(pi)
      END IF

      CALL SYSTEM_CLOCK(itime)
      irand = -INT(itime)
      WRITE(10,*) 'IRAND=', irand

     ! every processors need to have same input parameters
      CALL readit                 !!! Read in parameters
      rtepsi = SQRT(epsi); rhos3 = rhos*rhos*rhos  ! define constants from read-in parameters

      CALL fftw_init              !!! Sets up MPI data structure in conjunction with FFTW
      CALL ALLOC_DATA             !!! allocates arrays

      IF (l_ihavedata) THEN

         CALL SetToleranceVector                        !!! Set the tolerance vector
         IF (samp .NE. 0.d0) CALL make_source           !!! Declares source terms
         IF (phi0amp .NE. 0.d0) CALL make_extflow       !!! Construct external flow
         CALL linterm_construct                         !!! Constructs linear terms
         
      END IF

      CALL pert                                       !!! Set the inital stream function
      CALL solver_init                                !!! Start CVODE, all processors must pass through solver
      IF (l_tracers) CALL tracers_init                !!! Initialize tracer parameters, all processors share the calculation

END SUBROUTINE initial


SUBROUTINE fftw_init
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C     subroutine fftw_init
!C
!C     Construct local data structure according to separated dimensions
!C     given by FFTW. Information of local modes is gathered on master node.
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE mpi
      USE variables
      USE prog_flags, ONLY: lp,le,lm,l_tracers
      USE fft_variables
      USE mpi_variables
      IMPLICIT NONE

      INTEGER :: i,j,n,kount,kount_conjg,ierr
      INTEGER :: kyinit, kyend

      IF (pid .EQ. masternode) THEN
         WRITE(6,'(/a,I5)') "Number of processors is: ", npes
         WRITE(6,'(2(a,I4,1x))') "nxreal = ", nxreal, "nyreal = ", nyreal
         WRITE(6,'(2(a,I4,1x)/)') "nxcmplx = ", nxcmplx, "nycmplx = ", nycmplx
      END IF

      ! declare FFTW plans according to REAL-type array
      CALL RFFTW2D_F77_MPI_CREATE_PLAN(planb,MPI_COMM_WORLD,nxreal,nyreal,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE)
      CALL RFFTW2D_F77_MPI_CREATE_PLAN(planf,MPI_COMM_WORLD,nxreal,nyreal,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE)

      IF (pid .EQ. masternode) THEN
         WRITE(6,*) "planb: ", planb
         WRITE(6,*) "planf: ", planf
      END IF

      ! Divide data -- these are REAL-type size, data is divided on the extended array
      !  Important variables are "nylocal", "yinit", and "local_size"
      !  -- number of elements of FFT array are divided according to "nylocal"
      !  -- hence, total components to solve depends on this allocation
      !  -- notice that the indices are swapped contrasting to FFTW c-implementation
      CALL RFFTWND_F77_MPI_LOCAL_SIZES(planb,nylocal,yinit, nxlocal_after_transform,xinit_after_transform,local_size)

      ncmplxfft = local_size/2              ! total number of CMPLX elements for FFT

      ! allocate workarray
      ALLOCATE( workarr(local_size) )

      ! Construct a mapping from smaller array to extended array, storing the indices in "kxext" and "kyext".
      ! Mapping is from normal indexing by number of element of [0, 1, 2, ..., nycmplx-2, nycmplx-1]
      ! to ky indexing of [0, 1, 2, ..., kymax, ..., nycmplx/2-1, -nycmplx/2, -nycmplx/2+1, ..., -kymax, ..., -2, -1]
      ! All processors will have this map as reference to which portion of the extended array it stores the data.
      ! The kx-direction will remain the same (1 to 1 mapping). Only the ky-direction has a different designation.
      ALLOCATE(kxext(0:kxmax)); ALLOCATE(kyext(0:nycmplx-1))
      ALLOCATE(yext(-nycmplx/2:INT(CEILING(REAL(nycmplx)/2.0))-1))  ! upperbound specification to deal with odd numbers
      kxext = 0; kyext = 0; yext = 0
      kxext(0:kxmax) = (/ (i,i=0,kxmax) /)                          ! defines the indices for kx array

      DO j = 0, nycmplx-1
         kyext(j) = MOD(j+nycmplx/2,nycmplx) - nycmplx/2      ! apply mapping function from j -> ky
         yext(kyext(j)) = j                                   ! defines mapping from ky -> j
!         IF (pid .EQ. masternode) WRITE(*,*) kyext(j), yext(kyext(j))
!         WRITE(6,*) "kyext(",j,") = ", kyext(j)
!         WRITE(6,*) "yext(", kyext(j),") = ", j
      END DO

      yend = yinit+nylocal-1
      kylmin = kyext(yinit); kylmax = kyext(yend)             ! bounds for locally stored modes


      ! Gather directory on master node. Specifies local data on each processor. Specifies the "blacksheep" node.
      ALLOCATE( localindex_dir(2,0:npes-1) )                  ! all nodes are required to have the directory, allocated this way to ease the MPI_GATHER command
      localindex_dir = 0                                      ! data strucure of index directory: {nylocal yinit}
      localindex = (/ nylocal,yinit /)

      CALL MPI_GATHER(localindex,2,MPI_INTEGER, &             ! gathers local indices on master node first
           localindex_dir,2,MPI_INTEGER, &
           masternode,MPI_COMM_WORLD,mpierr)  

      blacksheep = -1                                         ! initial value, which should be defined

      IF (pid .EQ. masternode) THEN                        ! master node determines the "blacksheep" node
         WRITE(6,'(/a)') "FFTW MPI data division..."
         DO i = 0, npes-1
            WRITE(6,'(a,I2,1x,4(a,I4,1x),/,7x,a,I8)', advance='no') &
                 "ID: ",i,"nylocal = ",localindex_dir(1,i),"yinit = ",localindex_dir(2,i), &
                 "kyinit = ",kyext(localindex_dir(2,i)), &
                 "kyend = ",kyext(localindex_dir(1,i)+localindex_dir(2,i)-1), &
                 "local_size = ", local_size
            IF ( kyext(localindex_dir(2,i)+localindex_dir(1,i)-1) .LT. &
                 kyext(localindex_dir(2,i)) ) THEN    
               blacksheep = i                              ! the black sheep node is the one that straddles positive and negative modes
               WRITE(6,'(a)') " ** black sheep **"
            ELSE
               WRITE(6,'(a)')
            END IF

         END DO

         IF (blacksheep .EQ. -1) THEN                      ! if none of the processor straddles the middle portion
            blacksheep = npes/2                            ! assign to one of the middle nodes
            WRITE(6,'(a,I3)') "** black sheep ** is ", blacksheep
         END IF

      END IF

      CALL MPI_BCAST(blacksheep,1,MPI_INTEGER,masternode,MPI_COMM_WORLD,mpierr)          ! master node broadcasts the blacksheep node to all processors
      CALL MPI_BCAST(localindex_dir,2*npes,MPI_INTEGER,masternode,MPI_COMM_WORLD,mpierr) ! master node broadcasts index directory

      IF (pid .EQ. masternode) THEN                        ! sends index directory to black sheep
         CALL MPI_SEND(localindex_dir, 2*npes, MPI_INTEGER, blacksheep, 665, MPI_COMM_WORLD, mpierr)

      ELSE IF (pid .EQ. blacksheep) THEN                   ! receives index directory from master node
            
         CALL MPI_RECV(localindex_dir, 2*npes, MPI_INTEGER, masternode, 665, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierr)
         WRITE(6,'(/a/)') "** Black sheep received the index directory. **"

      END IF

!      CALL split_comm_group       !!! Splits MPI comminications, defines mpi_comm_solver & mpi_comm_nosolver


      ! Determine total number of modes to push through ODE solver
      ! -- Sets the maximum and minimum ky modes to solve
      ! -- Consecutively add to "nn_local" depending on the ky mode specified. Conjugate modes are not pushed through the solver.
      ! -- Instead of scanning through all ky-modes, this could also be implemented with logic statements
      !    comparing upper and lower bounds to -kymax and kymax.
      kount = 0; nn_local = 0; kylmin_solv = 0; kylmax_solv = 0
      DO j = yinit, yend                                         ! scans though number ky-elements
         IF (estoy_bien(kyext(j),-kymax,kymax)) THEN             ! if ky-element is in the range [-kymax,kymax]
            IF (kount .EQ. 0) THEN                               ! if this is first loop, then stores starting values
               kylmin_solv = kyext(j)
               kylmax_solv = kyext(j)
            ELSE                                                 ! tests if the current ky-element is the smallest or largest
               kylmin_solv = MIN(kyext(j),kylmin_solv)
               kylmax_solv = MAX(kyext(j),kylmax_solv)
            END IF
            kount = kount + 1
         END IF

         ! dealing with conjugate modes
         IF (estoy_bien(kyext(j),-kymax,0)) THEN    ! if the maximum ky mode is NOT in the conjugate modes range  (conjugate mode range: 1 <= ky <= kymax)
            nn_local = nn_local + kxm2

            IF (kyext(j) .NE. 0) THEN               ! skips the ky = 0 modes
               IF (kount_conjg .EQ. 0) THEN         ! find the required modes to be conjugated
                  kylmin_conjg = kyext(j)
                  kylmax_conjg = kyext(j)
               ELSE
                  kylmin_conjg = MIN(kyext(j),kylmin_conjg)
                  kylmax_conjg = MAX(kyext(j),kylmax_conjg)
               END IF
               kount_conjg = kount_conjg + 1
            END IF

         ELSEIF (estoy_bien(kyext(j),1,kymax)) THEN ! if the maximum ky mode is in the conjugate modes range
            nn_local = nn_local + kxmax
         END IF


      END DO

      kymodes = kount                    ! number of non-zero ky modes
      IF (kymodes .LT. 0) STOP "Something wrong with calculating 'kymodes'."


      ! Set processor flag if it contains data. Nodes without data remain idle until FFTW calls. Only "blacksheep" writes output files.
      IF (kount .EQ. 0) THEN
         WRITE(6,'(a)') ""
         WRITE(6,300) "ID: ",pid," does not have any local data."
         WRITE(6,'(a)') ""
         l_ihavedata = .FALSE.
      ELSE
         WRITE(6,'(a)') ""
         WRITE(6,300) "ID: ",pid," has local data from ", kylmin_solv," to ",kylmax_solv, &
              ", a total of ", kount," modes."
         WRITE(6,300) "ID: ",pid," has conjugate modes from ", kylmin_conjg," to ",kylmax_conjg, &
              ", a total of ", kount_conjg," modes."
         WRITE(6,300) "ID: ",pid," index starts from yinit = ",yinit, &
              " to yend = ",yend, &
              " a total of ",yend-yinit+1," elements."
         WRITE(6,'(a)') ""
         l_ihavedata = .TRUE.
      END IF
300   FORMAT(a,4(I4,a))

      nn_local = 2*nn_local        ! total number of terms including complex part per field
      nn_local = 3*nn_local        ! total number of modes for three fields

      nytransfer = INT(nycmplx/npes) + 1      ! number of y-elements to transfer, larger size to eliminate data truncation

      WRITE(6,*) "Number of elements to solve = ", nn_local, " Global length = ", nn

      ! setup FFTW on the node "blacksheep" for output
      IF (pid .EQ. blacksheep) THEN
         CALL RFFTW2D_F77_CREATE_PLAN(planb_serial,nxreal,nyreal,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE)
         CALL RFFTW2D_F77_CREATE_PLAN(planf_serial,nxreal,nyreal,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE)  ! <-- serial forward transforms are not used
      END IF

      CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)

END SUBROUTINE fftw_init


SUBROUTINE readit
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C     subroutine READIT
!C
!C     Reads in the parameter values.
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE fft_variables
      USE variables
      USE prog_flags
      USE mpi_variables
      USE file_units
      USE tracers
      IMPLICIT NONE
 
!C   Namelist input.
      NAMELIST /runp/ dt,ntstep,nwrt,nwrs,nvisit,tol,kxnorm,kynorm
      NAMELIST /physp/ v1,rhos,nue,emass,epsi,gama1,gama2,gama3,din,dout,dmid,dcoeff,dmeanf
      NAMELIST /physp2/ lp,le,lm
      NAMELIST /intilp/ amp,ku,kl,kg1sq,kg3sq,ratio
      NAMELIST /linzones/ kg1sq,kg3sq,kd1sq,kd3sq,kinit
      NAMELIST /pertin/ init,psi0,flreal,flima, samp, sscale, phi0amp
      NAMELIST/parts/ nvpart,nmax
      NAMELIST /flags/ l_tracers,drivedamp

      OPEN(4,file='inputs_mpi',status='old')

      ! all nodes read
      READ (4,NML = runp)
      READ (4,NML = physp)
      READ (4,NML = physp2)
      READ (4,NML = intilp)
      READ (4,NML = linzones)
      READ (4,NML = pertin)
      READ (4,NML = parts)
      READ (4,NML = flags)

      ! master node writes ONLY
      IF (pid .EQ. blacksheep) THEN
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

      CLOSE(4)
 
END SUBROUTINE readit


SUBROUTINE SetToleranceVector
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C     
!C     subroutine SetToleranceVector
!C  
!C	  This routine will set the values of a tolerance vector 
!C        to be passed to VODPK. The layout follow the format of the
!C        Y(nn) vector, where each entry is inversly proportional
!C        to the wavenumber squared  (i.e. ATOL(i) = TOL/(kx^2+ky^2))
!C	  where k^2 = (kx^2+ky^2).  This yields tighter absolute tolerance
!C        in higher modes.
!C
!C     routine called for nodes that have data only
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE variables, ONLY: tol,tola
      USE mpi_variables
      USE fft_variables, ONLY: kxmax,kymax
      IMPLICIT NONE

      INTEGER :: i,j
      INTEGER :: kx, ky, kount, kountm, kountd, kounts

      kount = 0
      kountm = nn_local/(2*3)    ! last element for REAL-type data for one field
      kountd = 2*kountm          ! last element for IMG-type data for one field
!      kountdm = 3*kountm         ! last element for REAL-type data for two fields
      kounts = 4*kountm          ! last element for IMG-type data for two fields

!C  Modes arrangement to coincide with serial version when running with only one CPU
      ! do for modes -kymax <= ky <= -1, kx = 0
      DO j = yinit, yend
         IF (estoy_bien(kyext(j),-kymax,-1)) THEN
            kount = kount + 1
            tola(kount) = tol/REAL(kyext(j)*kyext(j))

         END IF
      END DO

      ! do for zeroth mode
      IF (kyext(yinit) .EQ. 0) THEN
         kount = kount + 1
         tola(kount) = tol
      END IF

      ! do for -kymax <= ky <= -1, 1 <= kx <= kxmax
      DO j = yinit,yend
         IF (estoy_bien(kyext(j),-kymax,-1)) THEN
            DO i = 1, kxmax
               kount = kount + 1
               tola(kount) = tol/REAL(kyext(j)*kyext(j) + kxext(i)*kxext(i))

            END DO
         END IF
      END DO

      ! do for 0 <= ky <= kymax, 1 <= kx <= kxmax
      DO j = yinit, yend
         IF (estoy_bien(kyext(j),0,kymax)) THEN
            DO i = 1, kxmax
               kount = kount + 1
               tola(kount) = tol/REAL(kyext(j)*kyext(j) + kxext(i)*kxext(i))

            END DO
         END IF
      END DO

      tola(kountm+1:kountd) = tola(1:kountm)   ! apply to imaginary part

      IF (kount .NE. kountm) STOP "Tolerance vector elements not matched."

      tola(kountd+1:kounts) = tola(1:kountd)   ! apply the same to second field
      tola(kounts+1:nn_local) = tola(1:kountd) ! apply the same to third field

END SUBROUTINE SetToleranceVector


SUBROUTINE pert
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  subroutine pert                                              12/16/94
!C
!C  Sets the initial perturbation.
!C 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE mpi
      USE fft_variables
      USE variables
      USE prog_flags, ONLY: init
      USE mpi_variables
      IMPLICIT NONE
 
      REAL, EXTERNAL :: RAN2

      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax) :: dent,phit,denavgt
      INTEGER :: i, j
      REAL :: realdummy, xk, yk, ak, a1, a2, a3
      REAL :: zmodeavg
      CHARACTER(LEN=256) :: restartfile_str
 
!C   Initialization.
      dent = czero; phit = czero; denavgt = czero
      zmodeavg = 100.d0

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
!C   4               Read in amplitudes from file initspec.
!C
!C   5               Peaked initial profile with random phases,
!C                   plus a jet or vortex read from file initspec.
!C
!C   6               Peaked initial profile with random phases,
!C                   no energy in the k < kinit modes, 
!C                   plus a jet or vortex read from file initspec.
 
      SELECT CASE (init)
      CASE (3)
         ! For random initial phases, only the "blacksheep" node will specify the modes and write to file.
         ! Other nodes will read from this file to ensure that complex conjugates are enforced.
         IF (pid .EQ. blacksheep) THEN
            WRITE(6,'(a)') "Peaked initial profile with random phases."
            WRITE(10,'(a)') "Peaked initial profile with random phases."

            DO j = -kymax, kymax
               yk = kynorm*REAL(j)
               DO i = 0, kxmax
                  xk = kxnorm*REAL(i)
                  ak = SQRT(xk*xk + yk*yk)
                  a1 = psi0/(1.d0 + ak**ratio)
                  a2 = a1*2.d0*SIN(2.d0*pi*RAN2(irand))
                  a3 = SQRT(ABS(2.d0*a1*a1 - a2*a2))
                  dent(i,j) = CMPLX(a3,a2)
               END DO
            END DO

            dent(0,0) = czero   ! <-- no average perturbation
            phit = dent

            denavgt(0,0) = CMPLX(zmodeavg,0.d0)

            dent(0,1:kymax) = CONJG(dent(0,-1:-kymax:-1))  ! enforce conjugacy
            phit(0,1:kymax) = CONJG(phit(0,-1:-kymax:-1))
            denavgt(0,1:kymax) = CONJG(denavgt(0,-1:-kymax:-1))

            ! write initial random spectrum
            OPEN(20,file='restart_init.bin',status='replace',form='unformatted')  ! this is a binary file
            WRITE(20) dent,phit,denavgt
            CLOSE(20)            
!            WRITE(*,*) "Finished writing 'restart_init.bin' file."
!            OPEN(20,file='restart_init.out',status='replace',form='formatted')  ! this is a binary file
!            WRITE(20,*) dent,phit,denavgt
!            CLOSE(20)   
         END IF
         
!         WRITE(*,*) pid, "blacksheep =", blacksheep
!         CALL MPI_BCAST(dent,kxm1*kym2,MPI_DOUBLE_COMPLEX,blacksheep,MPI_COMM_WORLD,mpierr) ! all processors wait for data
!         CALL MPI_BCAST(phit,kxm1*kym2,MPI_DOUBLE_COMPLEX,blacksheep,MPI_COMM_WORLD,mpierr)
!         CALL MPI_BCAST(denavgt,kxm1*kym2,MPI_DOUBLE_COMPLEX,blacksheep,MPI_COMM_WORLD,mpierr)
         CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)

         IF (l_ihavedata) THEN
            OPEN(48,file='restart_init.bin',status='old',form='unformatted') ! all processors designated with data read from file
            READ(48) dent,phit,denavgt
            CLOSE(48)

            dent(0,1:kymax) = CONJG(den(0,-1:-kymax:-1))  ! enforce conjugacy
            phit(0,1:kymax) = CONJG(phi(0,-1:-kymax:-1))
            denavgt(0,1:kymax) = CONJG(denavg(0,-1:-kymax:-1))

            DO j = -kymax,kymax                                             ! sorts data
               IF (estoy_bien(j,kylmin_solv,kylmax_solv)) THEN              ! restricts to modes within modes pushed through ODE
                  den(0:kxmax,yext(j)) = dent(0:kxmax,j)
                  phi(0:kxmax,yext(j)) = phit(0:kxmax,j)
                  denavg(0:kxmax,yext(j)) = denavgt(0:kxmax,j)

               END IF
            END DO

         END IF


      CASE (4)
         IF (pid .EQ. blacksheep) THEN
            WRITE(6,*) 'Read in amplitudes from restart.bin file'
            WRITE(10,*)'Read in amplitudes from restart.bin file'
         END IF

!         restartfile_str = "restart_betatf40sat.bin" ! for testing only!!
         restartfile_str = "restart_betatf40.bin"
!         restartfile_str = "restart_betatf128_cont.bin"
!               restartfile_str = "restart_betatf10.bin"

         IF (pid .EQ. blacksheep) THEN
            WRITE(6,'(/a)', advance='no') "Specified restart.bin file: "
            WRITE(6,*) TRIM(ADJUSTL(restartfile_str))
         END IF

         IF (l_ihavedata) THEN
            OPEN(48,file=TRIM(ADJUSTL(restartfile_str)),status='old',form='unformatted') ! all processors designated with data read from file
!            READ(48) realdummy,realdummy,realdummy,realdummy,realdummy
            READ(48) dent,phit,denavgt
            CLOSE(48)

            dent(0,1:kymax) = CONJG(dent(0,-1:-kymax:-1))  ! enforce conjugacy
            phit(0,1:kymax) = CONJG(phit(0,-1:-kymax:-1))
            denavgt(0,1:kymax) = CONJG(denavgt(0,-1:-kymax:-1))

            DO j = -kymax,kymax                                             ! sorts data
               IF (estoy_bien(j,kylmin_solv,kylmax_solv)) THEN              ! restricts to modes within modes pushed through ODE
                  den(0:kxmax,yext(j)) = dent(0:kxmax,j)
                  phi(0:kxmax,yext(j)) = phit(0:kxmax,j)
                  denavg(0:kxmax,yext(j)) = denavgt(0:kxmax,j)

               END IF
            END DO

         END IF

      CASE DEFAULT
         CALL MPI_FINALIZE(mpierr)
         STOP 'WARNING: Must have INIT == 4'
      END SELECT

      IF (pid .EQ. blacksheep) WRITE(6,'(a)') "Initial perturbation specified."

END SUBROUTINE pert


SUBROUTINE solver_init
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine solver_init                                      01/13/12
!C
!C   Initializes the ODE solver
!C 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      USE mpi
      USE fft_variables
      USE cvode_variables
      use mpi_variables
      USE variables
      IMPLICIT NONE

      INTEGER :: ierr
      REAL, DIMENSION(nn_local) :: y

      CALL FNVINITP(MPI_COMM_WORLD,1,INT(nn_local,kind=8),INT(nn,kind=8),ierr)    ! parallel call
      IF (IERR .NE. 0) THEN
         WRITE(6,20) IERR
20       FORMAT(///' SUNDIALS_ERROR: FNVINITS returned IER = ', I5)
         STOP
      ENDIF

!C  Initialize CVODE solver
!C   Convert initial data to vector form
      CALL convert(den,phi,denavg,y)

!C   Allocate memory before solving
!C     CALL FCVMALLOC(T0=t, Y0=y, METH=2, ITMETH=1, IATOL=2, RTOL=tol,
!C    & ATOL=tola, IOUT=iout, ROUT=rout, IPAR=ipar, RPAR=rpar, IER=ierr)
      CALL FCVMALLOC(t, y, 2, 1, 2, tol, tola, iout, rout, ipar, rpar, ierr)
      IF (IERR .NE. 0) THEN
         WRITE(6,30) IERR
30       FORMAT(///' SUNDIALS_ERROR: FCVMALLOC returned IER = ', I5)
         STOP
      ENDIF


!   CALL FCVSPGMR(IPRETYPE, IGSTYPE, MAXL, DELT, IER)
!     IPRETYPE = preconditioner, 0 for no preconditioner
!     IGSTYPE  = Gram-Schimdt process, 1 for modified, 2 for classic
!     MAXL     = maximum Krylov subspace dimension
!     DELT     = linear convergence, 0.0 for default
      CALL FCVSPGMR(0,2,5,0.0,IERR)
!      CALL FCVSPTFQMR(0,2,5,0.0,IERR)
      IF (IERR .NE. 0) THEN
         WRITE(6,40) IERR
40       FORMAT(///' SUNDIALS_ERROR: FCVDIAG returned IER = ', I5)
         STOP
      ENDIF

!!$C   Use SPGMR treatment (might try other solving routines also)
!!$C     IPRETYPE = 2         : right-only preconditioning
!!$C     IGSTYPE  = 1         : modified Gram-Schmidt process
!!$C     MAXL     = MIN(5,nn) : Maximum Krylov subspace dimension, imported from dvodpk.f specifications
!!$C     DELT     = DEFAULT   : linear convergence tolerance factor, DEFAULT = 0.0
!!$      CALL FCVSPGMR(IPRETYPE=2, IGSTYPE=1, MAXL=5, DELT=0.0, IER=ierr)
!!$      IF (IERR .NE. 0) THEN
!!$         WRITE(6,50) IERR
!!$ 50      FORMAT(///' SUNDIALS_ERROR: FCVSPGMR returned IER = ', I5)
!!$         CALL FCVFREE
!!$         STOP
!!$      ENDIF

      CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)
      IF (pid .EQ. blacksheep) WRITE(6,'(a)') "Solver started."

END SUBROUTINE solver_init




SUBROUTINE timead
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C                                                                      C
!C  subroutine timead                                                   C
!C                                                                      C
!C      Calls the ODE solver vodpk.                                     C
!C                                                                      C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE variables
      USE cvode_variables
      USE fft_variables
      use mpi_variables
      IMPLICIT NONE

      INTEGER :: ierr,j
      REAL :: tf
      REAL, DIMENSION(nn_local) :: y
      COMPLEX, DIMENSION(0:kxmax,yinit:yend) :: dentempout,junkout,dendiff

      tf = t + dt

      CALL convert(den,phi,denavg,y)   ! flattens data into 1D array, all processors need to call this

!C   Call CVODE to solve for a timestep
!C     CALL FCVODE(TOUT=tf, T=t, Y=y, ITASK=1, IER=ierr)
      CALL FCVODE(tf, t, y, 1, ierr)

!C   Check for problems (refer to reference for more outputs)
      IF (pid .EQ. blacksheep) THEN
         WRITE(95,50) T, IOUT(3), IOUT(4), IOUT(9), ROUT(2)
50       FORMAT(/' t = ', E11.3, 3X, 'nsteps = ', I5,'  nfuncevals = ',I4,'  q = ', I2, '  h = ', E14.6)
      END IF

      IF (IERR .NE. 0) THEN
         WRITE(6,60) IERR, IOUT(15)
60       FORMAT(///' SUNDIALS_ERROR: FCVODE returned IER = ', I5, /, &
              '                 Linear returned IER = ', I5)
         CALL FCVFREE
         STOP
      ENDIF
         
!C   Convert y back to den
 
      CALL invert(y,den,phi,denavg)    ! expands data into 2D array and updates density array, all processors need to call this

END SUBROUTINE timead


SUBROUTINE convert(den_in,phi_in,denavg_in,y)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  subroutine convert                                           12/14/94
!C
!C   Puts a complex two-dimensional array, u(kx,ky), 
!C   into a one-dimensional array, y.
!C 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE fft_variables
      USE mpi_variables
      IMPLICIT NONE
 
      INTEGER :: j
      INTEGER :: kxarg, kount, kountm,  kountd, kountdm, kounts, kountsm
      COMPLEX,DIMENSION(0:kxmax,yinit:yend), INTENT(in):: den_in,phi_in,denavg_in
      REAL, DIMENSION(nn_local), INTENT(out):: y

      kount = 0
      kountm = nn_local/(2*3)               ! half REAL-type and half CMPLX-type
      kountd = 2*kountm                     ! starting index for real parts of potential terms
      kountdm = 3*kountm                    ! starting index for imaginary parts of potential terms
      kounts = 4*kountm                     ! starting index for mean field
      kountsm = 5*kountm                    ! starting index for imaginary parts of mean filed terms

      y = 0.d0

!C  Modes pushed through the solver are separated into two sections to deal with conjugate modes.
!C  The arrangement is made such that this coincides with the serial case on one CPU
      ! first portion coincides with the kx = 0 modes
      DO j = yinit, yend
         IF (estoy_bien(kyext(j),-kymax,-1)) THEN
            kount = kount + 1
            y(kount) = REAL( den_in(0,j) )
            y(kountm+kount) = AIMAG( den_in(0,j) )
            y(kountd+kount) = REAL( phi_in(0,j) )
            y(kountdm+kount) = AIMAG( phi_in(0,j) )
            y(kounts+kount) = REAL( denavg_in(0,j) )
            y(kountsm+kount) = AIMAG( denavg_in(0,j) )
         END IF
      END DO

      ! zeroth mode
      IF (kyext(yinit) .EQ. 0) THEN
         kount = kount + 1
         y(kount) = REAL( den_in(0,yinit) )
         y(kountm+kount) = AIMAG( den_in(0,yinit) )
         y(kountd+kount) = REAL( phi_in(0,yinit) )
         y(kountdm+kount) = AIMAG( phi_in(0,yinit) )
         y(kounts+kount) = REAL( denavg_in(0,yinit) )
         y(kountsm+kount) = AIMAG( denavg_in(0,yinit) )
      END IF

      ! do modes -kymax <= ky <= -1, 1 <= kx <= kxmax
      DO j = yinit, yend
         IF (estoy_bien(kyext(j),-kymax,-1)) THEN
            DO kxarg = 1, kxmax
               kount = kount + 1
               y(kount) = REAL( den_in(kxarg,j) )
               y(kountm+kount) = AIMAG( den_in(kxarg,j) )
               y(kountd+kount) = REAL( phi_in(kxarg,j) )
               y(kountdm+kount) = AIMAG( phi_in(kxarg,j) )
               y(kounts+kount) = REAL( denavg_in(kxarg,j) )
               y(kountsm+kount) = AIMAG( denavg_in(kxarg,j) )
            END DO
         END IF
      END DO

      ! do modes 0 <= ky <= kymax, 1 <= kx <= kxmax
      DO j = yinit, yend
         IF (estoy_bien(kyext(j),0,kymax)) THEN
            DO kxarg = 1, kxmax
               kount = kount + 1
               y(kount) = REAL( den_in(kxarg,j) )
               y(kountm+kount) = AIMAG( den_in(kxarg,j) )
               y(kountd+kount) = REAL( phi_in(kxarg,j) )
               y(kountdm+kount) = AIMAG( phi_in(kxarg,j) )
               y(kounts+kount) = REAL( denavg_in(kxarg,j) )
               y(kountsm+kount) = AIMAG( denavg_in(kxarg,j) )
            END DO
         END IF
      END DO

      IF (kount .NE. kountm) STOP "CONVERT routine: Number of data elements does not match."

END SUBROUTINE convert



SUBROUTINE invert(y,den_out,phi_out,denavg_out)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine invert                                          12/16/94
!C
!C   Puts a one-dimensional array into the 2-D complex array, u(kx,ky).
!C 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE fft_variables
      USE mpi_variables
      USE variables, ONLY: czero
      IMPLICIT NONE

      INTEGER :: j
      INTEGER :: kxarg, kount, kountm, kountd, kountdm, kounts, kountsm
      COMPLEX, DIMENSION(0:kxmax,yinit:yend), INTENT(out):: den_out,phi_out,denavg_out
      REAL, DIMENSION(nn_local), INTENT(in):: y

      kount = 0
      kountm = nn_local/(2*3)
      kountd = 2*kountm
      kountdm = 3*kountm
      kounts = 4*kountm
      kountsm = 5*kountm

      den_out = czero; phi_out = czero; denavg_out = czero    ! initialize to zero

!C  Replicate data arrangement such that it coincides with serial version for one cpu
      ! extract from first portion first
      DO j = yinit, yend
         IF (estoy_bien(kyext(j),-kymax,-1)) THEN
            kount = kount + 1
            den_out(0,j) = CMPLX( y(kount), y(kountm+kount) )
            phi_out(0,j) = CMPLX( y(kountd+kount), y(kountdm+kount) )
            denavg_out(0,j) = CMPLX( y(kounts+kount), y(kountsm+kount) )
         END IF
      END DO

      ! zeroth mode
      IF (kyext(yinit) .EQ. 0) THEN
         kount = kount + 1
         den_out(0,yinit) = CMPLX( y(kount), y(kountm+kount) )
         phi_out(0,yinit) = CMPLX( y(kountd+kount), y(kountdm+kount) )
         denavg_out(0,yinit) = CMPLX( y(kounts+kount), y(kountsm+kount) )
      END IF

      ! do for -kymax <= ky <= -1, 1 <= kx <= kxmax
      DO j = yinit, yend
         IF (estoy_bien(kyext(j),-kymax,-1)) THEN
            DO kxarg = 1, kxmax
               kount = kount + 1
               den_out(kxarg,j) = CMPLX( y(kount), y(kountm+kount) )
               phi_out(kxarg,j) = CMPLX( y(kountd+kount), y(kountdm+kount) )
               denavg_out(kxarg,j) = CMPLX( y(kounts+kount), y(kountsm+kount) )
            END DO
         END IF
      END DO

      ! do for 0 <= ky <= kymax, 1 <= kx <= kxmax
      DO j = yinit, yend
         IF (estoy_bien(kyext(j),0,kymax)) THEN
            DO kxarg = 1, kxmax
               kount = kount + 1
               den_out(kxarg,j) = CMPLX( y(kount), y(kountm+kount) )
               phi_out(kxarg,j) = CMPLX( y(kountd+kount), y(kountdm+kount) )
               denavg_out(kxarg,j) = CMPLX( y(kounts+kount), y(kountsm+kount) )
            END DO
         END IF
      END DO

!      WRITE(6,*) "INVERT called."
      IF (kount .NE. kountm) STOP "INVERT routine: Number of data elements does not match."

      CALL GATHER_CONJG(den_out,den_out)  ! share conjugate modes
      CALL GATHER_CONJG(phi_out,phi_out)
      CALL GATHER_CONJG(denavg_out,denavg_out)

END SUBROUTINE invert



SUBROUTINE FCVFUN(T, Y, YDOT, IPAR, RPAR, IER)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   SUBROUTINE FCVFUN: specifies right-hand side of diff.eq.
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE fft_variables
      USE mpi_variables
      IMPLICIT NONE

      REAL, INTENT(in):: t
      INTEGER(kind=8), DIMENSION(:), INTENT(in):: ipar
      REAL, DIMENSION(:), INTENT(in):: rpar
      REAL, DIMENSION(nn_local), INTENT(in):: y
      REAL, DIMENSION(nn_local), INTENT(out):: ydot
      INTEGER, INTENT(out):: IER
      COMPLEX, DIMENSION(0:kxmax,yinit:yend):: den_in,phi_in,denavg_in, &
           denkp, phikp, denavgkp

!C  Convert y to den.
 
      CALL invert(y,den_in,phi_in,denavg_in)
 
!C  Get y-prime.  (The right hand side)
 
      ! all processors need to call this due to the solver MPI call
      CALL fden(t,den_in,phi_in,denavg_in,denkp,phikp,denavgkp)
 
!C  Convert denkp back to yprime.
         
      CALL convert(denkp,phikp,denavgkp,ydot)

!C  Return success
      IER = 0


END SUBROUTINE FCVFUN


SUBROUTINE fden(tnow,den_in,phi_in,denavg_in,denkp,phikp,denavgkp)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine fden                                             12/15/94
!C
!C   This calculates the derivative and also handles the particles.
!C 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE variables
      USE prog_flags, ONLY: lp,le,lm,l_tracers
      USE fft_variables
      USE mpi_variables
      USE tracers, ONLY: told
      IMPLICIT NONE
 
      INTEGER :: i,j,n
      REAL, INTENT(IN) :: tnow
      COMPLEX, DIMENSION(0:kxmax,yinit:yend), INTENT(IN) :: den_in,phi_in,denavg_in
      COMPLEX, DIMENSION(0:kxmax,yinit:yend), INTENT(OUT) :: denkp,phikp,denavgkp
      COMPLEX, DIMENSION(0:kxmax,yinit:yend) :: dennl,phinl,denavgnl,drivenl,phi_temp
      COMPLEX, DIMENSION(0:kxmax,yinit:yend) :: junkin,junkout

      denkp = czero; phikp = czero; denavgkp = czero; junkin = czero
      dennl = czero; phinl = czero; denavgnl = czero; drivenl = czero

!C   Calculate df/dt=denkp at kx, ky where f is the stream function.

      IF (l_ihavedata) THEN 
         phi_temp = phi_in + phi0
         IF ((lp .NE. 0.0) .OR. (le .NE. 0.0) .OR. (lm .NE. 0.0)) THEN 
            CALL convolve(den_in,phi_temp,denavg_in,dennl,phinl,denavgnl,drivenl)   ! check flags
         END IF
         IF (lm .EQ. 0.0) drivenl = philindrivearr*phi_temp    ! constant gradient if mean field is not present
!         drivenl = philindrivearr*phi_in  ! <-- over ride to fixed gradient drive

         phikp = (-1.d0*(1.d0 - rtepsi*dcoeff)*v1*drivenl - damplinarr*phi_temp &
              + nue*rtepsi*(drivelinarr*den_in - phi_in) + 0.5d0*lp*dcoeff*phinl)*philinarr
         denkp = -1.d0*v1*dcoeff*drivenl + nue*(phi_in - den_in) - 0.5d0*le*dcoeff*dennl

         IF (lm .NE. 0.0) denavgkp = densource - 0.5d0*lm*denavgnl - dampdenavglinarr*denavg_in

      ELSE    ! idle processors simply wait for FFTW calls

         IF ((lp .NE. 0.0) .OR. (le .NE. 0.0) .OR. (lm .NE. 0.0)) THEN 
            CALL convolve(junkin,junkin,junkin,junkout,junkout,junkout,junkout)
         END IF

      END IF

!C   Advance tracers with current field -- apply to all processors
      IF (l_tracers) THEN
         CALL tracers_advance(tnow - told)
         told = tnow
      END IF

END SUBROUTINE fden


SUBROUTINE linterm_construct
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine linterm_construct                                22/12/11
!C
!C   Construct linear portions of the propagator
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE fft_variables
      USE variables
      use mpi_variables
      IMPLICIT NONE
 
      logical :: temp_logic
      INTEGER :: i,j,kxarg,kyarg
      REAL :: xk, yk, ak, ak1
      REAL :: kd1, kd3, kg1, kg3

      philinarr = czero; damplinarr = czero; drivelinarr = czero; dampdenavglinarr = czero
      philindrivearr = czero

      kg1 = SQRT(REAL(kg1sq)); kg3 = SQRT(REAL(kg3sq))  ! this is used for square drive and damp regions

!C   Start the main loops over the modes.

      DO j = yinit, yend
         
         IF (estoy_bien(kyext(j),-kymax,kymax)) THEN
        
            kyarg = kyext(j)

            DO kxarg = 0, kxmax
 
!C   Defining normalized components of the k-vector.
 
               xk = kxnorm*REAL(kxarg)
               yk = kynorm*REAL(kyarg)
 
!C   The magnitudes of k^2 and k^4
 
               ak = yk*yk + xk*xk
               ak1 = (1.0d0 - rtepsi + ak*rhos*rhos)

               IF ((kxarg .EQ. 0) .AND. (kyarg .EQ. 0)) THEN        ! Enforce no drive nor damp on (0,0) mode
                  damplinarr(kxarg,j) = czero; drivelinarr(kxarg,j) = czero
!                  WRITE(6,*) "ID:",pid,"(0,0) mode evaluation!"
               ELSE


!C   Form the linear damping terms
!C
!C       gama1 : damping proportional to k^0.
!C       gama2 : damping proportional to k^2.
!C       gama3 : damping proportional to k^4.
!C
                  IF ((ABS(kyarg) .LE. kg1) .AND. (ABS(kxarg) .LE. kg1)) THEN ! (a square)
                     damplinarr(kxarg,j) = gama1
                  ELSE IF ((ABS(kyarg) .GT. kg3) .OR. (ABS(kxarg) .GT. kg3)) THEN
                     damplinarr(kxarg,j) = gama3*ak*ak
                  ELSE
                     damplinarr(kxarg,j) = gama2*ak
                  ENDIF

!C   Low-k driving proportional to din
!C   High-k driving proportional to dout
!C   Mid-k driving proportional to dmid

                  ! Specify drive regions
                  IF (ak .LT. REAL(kd1sq)) THEN
                     drivelinarr(kxarg,j) = din
                  ELSE IF (ak .GE. REAL(kd3sq)) THEN
                     drivelinarr(kxarg,j) = dout
                  ELSE
                     drivelinarr(kxarg,j) = dmid
                  END IF

                  ! linear damping for the mean field
                  dampdenavglinarr(kxarg,j) = dmeanf * ak

               END IF

                  ! linear division array for potential (phi)
                  philinarr(kxarg,j) = 1.d0/ak1

                  ! linear diamagnetic drift term
                  philindrivearr(kxarg,j) = - ic * yk

            END DO

         END IF

      END DO


CONTAINS
      LOGICAL FUNCTION write_me( some_arr,file_unit )
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C     LOGICAL FUNCTION write_me :
!C        -- writes linear array
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
            USE fft_variables, ONLY: kxmax, kymax
            use mpi_variables
            IMPLICIT NONE
            INTEGER :: i,j
            INTEGER, INTENT(in) :: file_unit
            COMPLEX, DIMENSION(0:kxmax,yinit:yend), INTENT(IN) :: some_arr
            
            DO j = yinit, yend
               WRITE(file_unit,*) (REAL(some_arr(i,j)*CONJG(some_arr(i,j))), i=0,kxmax)
!               WRITE(file_unit,*) (some_arr(i,j), i=0,kxmax)
            END DO

            write_me = .TRUE.

      END FUNCTION write_me

END SUBROUTINE linterm_construct


SUBROUTINE convolve(den_in,phi_in,denavg_in,dennl,phinl,denavgnl,drivenl)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine convolve                                        12/14/94
!C
!C     This subroutine handles the nonlinear convolutions.
!C 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      USE variables
      USE prog_flags, ONLY: lp,le,lm
      USE fft_variables
      USE mpi_variables
      USE convolve_module
      IMPLICIT NONE


!C  Argument declarations.

      COMPLEX, DIMENSION(0:kxmax,yinit:yend), INTENT(in) :: den_in, phi_in, denavg_in
      COMPLEX, DIMENSION(0:kxmax,yinit:yend), INTENT(out) :: dennl, phinl, denavgnl, drivenl
      REAL, DIMENSION(local_size) :: rwy_flat  ! testing variables
      COMPLEX, DIMENSION(ncmplxfft) :: c2wy_flat

!C  Misc...
 
      INTEGER :: i, j, jstart
      REAL :: kx,ky,ksq
      
!C  Initialize arrays to zero

      dennl = czero; phinl = czero; denavgnl = czero; drivenl = czero
      cvex_flat = czero; cvey_flat = czero;
      IF (lp .NE. 0.0) THEN 
         cwx_flat = czero; cwy_flat = czero; Pol_flat = 0.0
      END IF
      IF (le .NE. 0.0) THEN 
         cnkx_flat = czero; cnky_flat = czero; ExB_flat = 0.0
      END IF
      IF (lm .NE. 0.0) THEN
         cnavgkx_flat = czero; cnavgky_flat = czero; drivenl_flat = 0.0
         cnavg_flat = czero;
         cnavgkx1d_flat = czero; cnavg1d_flat = czero
      END IF


!      WRITE(*,*) "ID: ", pid, " CONVOLVE called! at ncl = ",ncl

!C  Calculates stuff to be transformed to x-space . . .

      IF (l_ihavedata) THEN                                       ! restrict to processors with data
         DO j = yinit, yend                                       ! restrict to data range only

            IF (estoy_bien(kyext(j),-kymax,kymax)) THEN

               DO i = 0, kxmax
                  kx = kxnorm*REAL(kxext(i)); ky = kynorm*REAL(kyext(j))
                  ksq = kx*kx + ky*ky

                  cvex_flat((j-yinit)*nxcmplx+1+i) =         ic * ky * phi_in(i,j)
                  cvey_flat((j-yinit)*nxcmplx+1+i) = -1.d0 * ic * kx * phi_in(i,j)

                  IF (lp .NE. 0.0) THEN
                     cwx_flat((j-yinit)*nxcmplx+1+i) = ic * kx * ksq * phi_in(i,j)
                     cwy_flat((j-yinit)*nxcmplx+1+i) = ic * ky * ksq * phi_in(i,j)
                  END IF
                  IF (le .NE. 0.0) THEN
                     cnkx_flat((j-yinit)*nxcmplx+1+i) = -1.d0 * ic * kx * den_in(i,j)
                     cnky_flat((j-yinit)*nxcmplx+1+i) = -1.d0 * ic * ky * den_in(i,j)
                  END IF
                  IF (lm .NE. 0.0) THEN
                     cnavgkx_flat((j-yinit)*nxcmplx+1+i) = -1.d0 * ic * kx * denavg_in(i,j)
                     cnavgky_flat((j-yinit)*nxcmplx+1+i) = -1.d0 * ic * ky * denavg_in(i,j)
                     cnavg_flat((j-yinit)*nxcmplx+1+i) = denavg_in(i,j)
                     IF (kyext(j) .EQ. 0) THEN             ! output only ky = 0 modes for gradient scale length. However, this is a 2D transform to maintain the form.
                        cnavgkx1d_flat((j-yinit)*nxcmplx+1+i) = -1.d0 * ic * kx * denavg_in(i,j)
                        cnavg1d_flat((j-yinit)*nxcmplx+1+i) = denavg_in(i,j)
                     END IF
                  END IF
               

               END DO

            END IF

         END DO
      END IF

!C  Off to REAL space . . .

      CALL RFFTWND_F77_MPI(planb,1,cvex_flat,workarr,1,FFTW_NORMAL_ORDER)       ! normal order gives the correct answer
      CALL RFFTWND_F77_MPI(planb,1,cvey_flat,workarr,1,FFTW_NORMAL_ORDER)


      ! Polarization non-linearity
      IF (lp .NE. 0.0) THEN

         CALL RFFTWND_F77_MPI(planb,1,cwx_flat,workarr,1,FFTW_NORMAL_ORDER)
         CALL RFFTWND_F77_MPI(planb,1,cwy_flat,workarr,1,FFTW_NORMAL_ORDER)

         Pol_flat = TRANSFER(cvex_flat,(/ 0.d0 /)) * TRANSFER(cwx_flat,(/ 0.d0 /)) &       ! Calculate the nonlinear terms
              + TRANSFER(cvey_flat,(/ 0.d0 /)) * TRANSFER(cwy_flat,(/ 0.d0 /))                             

         CALL RFFTWND_F77_MPI(planf,1,Pol_flat,workarr,1,FFTW_NORMAL_ORDER)              ! Back to Fourier space . . .

         IF (l_ihavedata) THEN
            ! initial index strides by "nxcmplx" but end index strides by "kxm2 = kxmax + 1", only [0,kxmax] modes are retained
            DO j = yinit,yend
               jstart = (j-yinit)*2*nxcmplx
               phinl(0:kxmax,j) = TRANSFER(Pol_flat( jstart+1:jstart+2*kxm2 ),(/ (0.d0,0.d0) /))
            END DO
         END IF

         phinl = phinl/REAL(nxreal*nyreal)

      END IF

      ! ExB non-linearity
      IF (le .NE. 0.0) THEN
         CALL RFFTWND_F77_MPI(planb,1,cnkx_flat,workarr,1,FFTW_NORMAL_ORDER)
         CALL RFFTWND_F77_MPI(planb,1,cnky_flat,workarr,1,FFTW_NORMAL_ORDER)

         ExB_flat = TRANSFER(cvex_flat,(/ 0.d0 /)) * TRANSFER(cnkx_flat,(/ 0.d0 /)) &      ! Calculate the nonlinear terms
              + TRANSFER(cvey_flat,(/ 0.d0 /)) * TRANSFER(cnky_flat,(/ 0.d0 /))                               

         CALL RFFTWND_F77_MPI(planf,1,ExB_flat,workarr,1,FFTW_NORMAL_ORDER)              ! Back to Fourier space . . .

         IF (l_ihavedata) THEN
            DO j = yinit,yend
               jstart = (j-yinit)*2*nxcmplx
               dennl(0:kxmax,j) =  TRANSFER(ExB_flat( jstart+1:jstart + 2*kxm2 ),(/ (0.d0,0.d0) /))
            END DO
         END IF

         dennl = dennl/REAL(nxreal*nyreal)   ! FFTW is not normalized

      END IF

      ! ExB non-linearity on the mean field and nonlinear drive term
      ! The flag to cover the case without the mean field is already declared as a linear phi drive term.
      IF (lm .NE. 0.0) THEN
         CALL RFFTWND_F77_MPI(planb,1,cnavgkx_flat,workarr,1,FFTW_NORMAL_ORDER)
         CALL RFFTWND_F77_MPI(planb,1,cnavgky_flat,workarr,1,FFTW_NORMAL_ORDER)
         CALL RFFTWND_F77_MPI(planb,1,cnavg_flat,workarr,1,FFTW_NORMAL_ORDER)
         CALL RFFTWND_F77_MPI(planb,1,cnavgkx1d_flat,workarr,1,FFTW_NORMAL_ORDER)
         CALL RFFTWND_F77_MPI(planb,1,cnavg1d_flat,workarr,1,FFTW_NORMAL_ORDER)
         
         denavgnl_flat = TRANSFER(cvex_flat,(/ 0.d0 /)) * TRANSFER(cnavgkx_flat,(/ 0.d0 /)) &      ! Calculate the nonlinear terms
              + TRANSFER(cvey_flat,(/ 0.d0 /)) * TRANSFER(cnavgky_flat,(/ 0.d0 /))       
         drivenl_flat = -1.d0*TRANSFER(cvex_flat,(/ 0.d0 /)) * TRANSFER(cnavgkx1d_flat,(/ 0.d0 /)) &
              / TRANSFER(cnavg1d_flat,(/ 0.d0 /))
!         drivenl_flat = -1.d0*TRANSFER(cvex_flat,(/ 0.d0 /)) * TRANSFER(cnavgkx_flat,(/ 0.d0 /)) &  ! nonlinear drive with local gradients
!              / TRANSFER(cnavg_flat,(/ 0.d0 /))
         
         CALL RFFTWND_F77_MPI(planf,1,denavgnl_flat,workarr,1,FFTW_NORMAL_ORDER)              ! Back to Fourier space . . .
         CALL RFFTWND_F77_MPI(planf,1,drivenl_flat,workarr,1,FFTW_NORMAL_ORDER)

         IF (l_ihavedata) THEN
            DO j = yinit,yend
               jstart = (j-yinit)*2*nxcmplx
               denavgnl(0:kxmax,j) = TRANSFER(denavgnl_flat( jstart+1:jstart + 2*kxm2 ),(/ (0.d0,0.d0) /))
               drivenl(0:kxmax,j) = TRANSFER(drivenl_flat( jstart+1:jstart + 2*kxm2 ),(/ (0.d0,0.d0) /))
            END DO
         END IF

         denavgnl = denavgnl/REAL(nxreal*nyreal)
         drivenl = drivenl/REAL(nxreal*nyreal)

      END IF

END SUBROUTINE convolve



SUBROUTINE FINALIZE
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine FINALIZE                                            
!C
!C   Close all files
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE visit
      USE variables
      USE fft_variables
      USE mpi_variables
      USE convolve_module
      USE prog_flags, ONLY: l_tracers
      IMPLICIT NONE
!C
!C   A little closure....
!C 

      DEALLOCATE(kxext); DEALLOCATE(kyext)
      IF (ALLOCATED(cvex_flat)) DEALLOCATE(cvex_flat)
      IF (ALLOCATED(cvey_flat)) DEALLOCATE(cvey_flat)
      IF (ALLOCATED(cnky_flat)) DEALLOCATE(cnky_flat)
      IF (ALLOCATED(cnkx_flat)) DEALLOCATE(cnkx_flat)
      IF (ALLOCATED(cwx_flat)) DEALLOCATE(cwx_flat)
      IF (ALLOCATED(cwy_flat)) DEALLOCATE(cwy_flat)
      IF (ALLOCATED(Pol_flat)) DEALLOCATE(Pol_flat)
      IF (ALLOCATED(ExB_flat)) DEALLOCATE(ExB_flat)

      IF (ALLOCATED(workarr)) DEALLOCATE(workarr)

      IF (l_ihavedata) THEN
         DEALLOCATE(den); DEALLOCATE(phi); DEALLOCATE(denavg)
         DEALLOCATE(philinarr); DEALLOCATE(drivelinarr); DEALLOCATE(damplinarr)
         DEALLOCATE(dampdenavglinarr); DEALLOCATE(philindrivearr)
         DEALLOCATE(densource); DEALLOCATE(tola)
         IF (ALLOCATED(phi0)) DEALLOCATE(phi0)
!         DEALLOCATE(temparr); DEALLOCATE(recvarr)
         CALL FCVFREE
      END IF
      
      IF (pid .EQ. blacksheep) THEN 
         WRITE(6,'(///4x,a///)') '  das ende'
         CLOSE( 95 ); CLOSE( 10 )

!         DEALLOCATE(kden)
         DEALLOCATE(denout); DEALLOCATE(phiout); DEALLOCATE(denavgout)
         DEALLOCATE(x_visit, y_visit); DEALLOCATE(kx_visit, ky_visit)

         CALL RFFTWND_F77_DESTROY_PLAN(planb_serial)
         CALL RFFTWND_F77_DESTROY_PLAN(planf_serial)

      END IF

      CALL RFFTWND_F77_MPI_DESTROY_PLAN(planb)
      CALL RFFTWND_F77_MPI_DESTROY_PLAN(planf)
      WRITE(6,'(a,I4,a)') "ID: ", pid, " FFTW plans cleared."
 
      IF (l_tracers) CALL tracers_final

END SUBROUTINE FINALIZE


SUBROUTINE energy
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine energy                                           12/14/94
!C
!C   Calculates the total energy and enstrophy.
!C 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
      use fft_variables
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

!      denekmin = 0.0; denenmin = 0.0    ! assign small values for initial comparisons    
!      denekmax = 0.0; denenmax = 0.0

!C  Reduced loops -----------------------------------------------
!C  -- first part calculates ky modes for kx = 0
!C  -- second part calculates the rest of modes

      xk = 0.0
      DO j = -kymax, -1
         yk = kynorm*j
         ak = xk*xk + yk*yk

         ! calculated different from one-field model
         denektemp = denout(0,j)*CONJG(denout(0,j))*2.0d0
         denentemp = denout(0,j)*CONJG(denout(0,j))*ak*2.0d0
         phiektemp = phiout(0,j)*CONJG(phiout(0,j))*2.0d0
         phiek1temp = phiout(0,j)*CONJG(phiout(0,j))*2.d0*(1.d0-rtepsi+ak*rhos*rhos)  ! <-- when need to calculate conserved quantity
         phientemp = phiout(0,j)*CONJG(phiout(0,j))*ak*2.0d0
         denavgektemp = denavgout(0,j)*CONJG(denavgout(0,j))*2.0d0

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

            denektemp = denout(i,j)*CONJG(denout(i,j))*2.0d0
            denentemp = denout(i,j)*CONJG(denout(i,j))*ak*2.0d0
            phiektemp = phiout(i,j)*CONJG(phiout(i,j))*2.0d0
            phiek1temp = phiout(i,j)*CONJG(phiout(i,j))*2.d0*(1.d0-rtepsi+ak*rhos*rhos)  ! <-- when need to calculate conserved quantity
            phientemp = phiout(i,j)*CONJG(phiout(i,j))*ak*2.0d0
            denavgektemp = denavgout(i,j)*CONJG(denavgout(i,j))*2.0d0

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
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: fileunit

      WRITE(fileunit,200) ncl,t,denek,denen,dengamek,dengamen,phiek,phiek1,phien,phigamek,phigamen,denavgek,ek
200   FORMAT(I10,1x,12(1pe16.9,1x))
 
END SUBROUTINE writeit


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


SUBROUTINE writerestart
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C       *********************************************************
!C       *                                                       *
!C       *  Writes restart file                  (DEN 6/17/04)   *
!C       *  and writes particle restart file	(DEN 7/27/05)   *
!C       *  now accepte jet values for constant jet - SD         *
!C       *********************************************************
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE variables, ONLY: denout,phiout,denavgout,ncl
      IMPLICIT NONE

      INTEGER :: i,j
      character(LEN=10) :: itstr

      CALL get_iterstr(itstr,ncl,LEN(itstr))
      OPEN(20,file='restart.'//TRIM(itstr)//'.bin',status='replace',form='unformatted')  ! this is a binary file 
      WRITE(20) denout,phiout,denavgout
      CLOSE(20)

END SUBROUTINE writerestart



SUBROUTINE INIT_VISIT
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   SUBROUTINE INIT_VISIT: initialize outputs to Silo for ViSit
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE fft_variables
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
      END DO 
      
      dk = 2.0*kymax/REAL(nky_visit - 1)
      DO j = 1, nky_visit
         ky_visit(j) = -kymax + REAL(j - 1)*dk
      END DO
!
!     Spatial arrays
!
      xdims(1) = nx_visit; xdims(2) = ny_visit 
      DO j = 1, nx_visit
         x_visit(j) = REAL(j-1)/REAL(nx_visit - 1)
      END DO
      DO j = 1, ny_visit
         y_visit(j) = REAL(j-1)/REAL(ny_visit - 1)
      END DO
      
END SUBROUTINE INIT_VISIT



SUBROUTINE WRITE_VISIT
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   SUBROUTINE WRITE_VISIT: write outputs to Silo for ViSit
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE visit
      USE variables, ONLY: denout,phiout,denavgout,t,ncl,kxnorm,kynorm,czero
      USE fft_variables
      use mpi_variables
      IMPLICIT NONE

      INTEGER :: ierr, err, optlistid, dbfile
      CHARACTER(LEN=6) :: citer
      CHARACTER(LEN=16) :: filename
      INTEGER :: i, j
      REAL :: ak, xk, yk, ksq
     
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax) :: cxcorr_visit
      COMPLEX, DIMENSION(0:nxcmplx-1,0:nycmplx-1) :: kden,kvort,kphi,kdenavg,kden2,klindrive
      REAL, DIMENSION(nxreal,nyreal) :: den_visit,vort_visit,phi_visit,denavg_visit,den2_visit,lindrive_visit
      REAL, DIMENSION(0:kxmax,-kymax:kymax) :: &
           denekk_visit, phiekk_visit, denavgekk_visit, phiekkphase_visit, denekkphase_visit, &
           xcorrphase_visit, linkdrive_visit
      REAL :: kx(0:nxcmplx-1),ky(0:nycmplx-1)

!      WRITE(*,*) "filename = ", filename, " ivisit = ", ivisit

      CALL get_iterstr(citer,ivisit,LEN(citer))  ! convert iteration number to character string with padded zeroes

      filename = 'beta.'//citer//'.silo'

      ierr = DBCREATE(filename, 16, DB_CLOBBER, DB_LOCAL, "BETA_RUN", 8, DB_HDF5, dbfile)

!      WRITE(*,*) "filename = ", filename, " Silo error = ", ierr

      ierr = DBMKOPTLIST(4, optlistid)
      ierr = DBADDIOPT(optlistid, DBOPT_CYCLE, ncl)
      ierr = DBADDDOPT(optlistid, DBOPT_DTIME, t)
      ierr = DBADDCOPT(optlistid, DBOPT_XLABEL, "X", 1)
      ierr = DBADDCOPT(optlistid, DBOPT_YLABEL, "Y", 1)

!C  Write density and vorticity
       
!C    Initialize

      kden = czero; kphi = czero; kvort = czero; kdenavg = czero; kden2 = czero
      klindrive = czero
      kx = 0.0; ky = 0.0

!C    Get real components
      kden(0:kxmax,0:kymax) = denout(0:kxmax,0:kymax)
      kden(0:kxmax,nycmplx-kymax:nycmplx-1) = denout(0:kxmax,-kymax:-1)
      kphi(0:kxmax,0:kymax) = phiout(0:kxmax,0:kymax)
      kphi(0:kxmax,nycmplx-kymax:nycmplx-1) = phiout(0:kxmax,-kymax:-1)
      kdenavg(0:kxmax,0:kymax) = denavgout(0:kxmax,0:kymax)
      kdenavg(0:kxmax,nycmplx-kymax:nycmplx-1) = denavgout(0:kxmax,-kymax:-1)

      kx(0:kxmax) = kxnorm*(/ (i,i=0,kxmax) /)
      ky(0:kymax) = kynorm*(/ (j,j=0,kymax) /)
      ky(nycmplx-kymax:nycmplx-1) = kynorm*(/ (j,j=-kymax,-1) /)

      ! calculate different sectors
      DO j = 0, kymax     ! calculates vorticity
         DO i = 0, kxmax
            ksq = kxnorm*kxnorm*kxext(i)*kxext(i) + kynorm*kynorm*kyext(j)*kyext(j)
            kvort(i,j) = -1.0 * ksq * kphi(i,j)
            kden2(i,j) = -1.0 * ksq * kden2(i,j)
         END DO
      END DO

!      WRITE(*,*) "kvort constructed at t = ", t

      DO j = nycmplx-kymax, nycmplx-1
         DO i = 0, kxmax
            ksq = kxnorm*kxnorm*REAL( kxext(i)*kxext(i) ) &
                 + kynorm*kynorm*REAL( kyext(j)*kyext(j) )
            kvort(i,j) = -1.0 * ksq * kphi(i,j)
            kden2(i,j) = -1.0 * ksq * kden2(i,j)
         END DO
      END DO
       

      ! FFTW on one node only
      CALL RFFTWND_F77_ONE_COMPLEX_TO_REAL(planb_serial,kden,den_visit)
      CALL RFFTWND_F77_ONE_COMPLEX_TO_REAL(planb_serial,kvort,vort_visit)
      CALL RFFTWND_F77_ONE_COMPLEX_TO_REAL(planb_serial,kphi,phi_visit)
      CALL RFFTWND_F77_ONE_COMPLEX_TO_REAL(planb_serial,kden2,den2_visit)
      CALL RFFTWND_F77_ONE_COMPLEX_TO_REAL(planb_serial,kdenavg,denavg_visit)

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


!C  Write k-wave spectrum energy      

      ierr = DBFREEOPTLIST(optlistid)
      ierr = DBMKOPTLIST(2, optlistid)
      ierr = DBADDCOPT(optlistid, DBOPT_XLABEL, "KX", 2)
      ierr = DBADDCOPT(optlistid, DBOPT_YLABEL, "KY", 2)

      ierr = DBPUTQM(dbfile, "kmesh", 5, "kx", 2, "ky", 2, "z", &
           1, kx_visit, ky_visit, DB_F77NULL, kdims, ndims, DB_DOUBLE, DB_COLLINEAR, optlistid, err)

      DO j = -kymax, kymax
         DO i = 0, kxmax
            denekk_visit(i,j) = 2.d0*denout(i,j)*CONJG(denout(i,j))    !!! Energy in a mode
            phiekk_visit(i,j) = 2.d0*phiout(i,j)*CONJG(phiout(i,j))
            denavgekk_visit(i,j) = 2.d0*denavgout(i,j)*CONJG(denavgout(i,j))
            phiekkphase_visit(i,j) = ATAN2(REAL(phiout(i,j)),AIMAG(phiout(i,j))) !!! phi phase
            denekkphase_visit(i,j) = ATAN2(REAL(denout(i,j)),AIMAG(denout(i,j))) !!! den phase
            cxcorr_visit(i,j) = denout(i,j)*CONJG(phiout(i,j))       !!! phase cross-correlation between fields
            xcorrphase_visit(i,j) = ATAN2(REAL(cxcorr_visit(i,j)),AIMAG(cxcorr_visit(i,j)))
         END DO
      END DO


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

      ierr = DBPUTQV1(dbfile, "xcorrphase", 10, "kmesh", 5, &
           xcorrphase_visit, kdims, ndims, DB_F77NULL, 0, DB_DOUBLE, DB_NODECENT, DB_F77NULL, err)

      ierr = DBFREEOPTLIST(optlistid)
      ierr = DBCLOSE(dbfile)

      ivisit = ivisit + 1
       
      WRITE(6,*) " SILO ==> "//filename//" written at t = ",t," ncl = ",ncl


END SUBROUTINE WRITE_VISIT


SUBROUTINE SECOND0(STIME)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  SUBROUTINE SECOND0(STIME): timing subroutine
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE

      REAL :: stime
      INTEGER, DIMENSION(8):: values
      CALL DATE_AND_TIME(VALUES=values)
      stime = 60.*(60.*values(5) + values(6)) + values(7)

END SUBROUTINE SECOND0

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



SUBROUTINE SPLIT_COMM_GROUP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   SUBROUTINE split_comm_group: Creates two new communicator groups;
!C     one for passing through the solver routine (mpi_comm_solver),
!C     the other with no solver routine (mpi_comm_nosolver).
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE mpi
      USE fft_variables
      USE mpi_variables
      IMPLICIT NONE

      INTEGER :: i,j
      INTEGER :: yinit_temp,yend_temp
      INTEGER :: kount_nodata, kount_data, numsolvprocs, numnosolvprocs
      INTEGER :: group_world, group_solver, group_nosolver
      INTEGER, DIMENSION(:), ALLOCATABLE :: solvgrpbuf,nosolvgrpbuf

      ! array containing number of processors in the solver group and their ranks
      ! format of the array is [no_procs_for_solver rank_ids...]
      ALLOCATE( solvgrpbuf(npes+1) ); solvgrpbuf = -1
      ALLOCATE( nosolvgrpbuf(npes+1) ); nosolvgrpbuf = -1

      kount_data = 0; kount_nodata = 0

      IF (pid .EQ. masternode) THEN                                      ! rootid sorts who is part of the solver
         DO i = 0,npes-1
            yinit_temp = localindex_dir(2,i)                             ! extract bounding indices for each processor
            yend_temp = yinit_temp + localindex_dir(1,i) - 1

            IF ( estoy_bien(kyext(yinit_temp),-kymax,kymax) .OR. &       ! check if bounding modes are within range
                 estoy_bien(kyext(yend_temp),-kymax,kymax) ) THEN 
               kount_data = kount_data + 1
               solvgrpbuf(kount_data+1) = i
            ELSE                                                         ! if not solving
               kount_nodata = kount_nodata + 1
               nosolvgrpbuf(kount_nodata+1) = i
            END IF

         END DO

         solvgrpbuf(1) = kount_data; nosolvgrpbuf(1) = kount_nodata

         WRITE(6,'(/2(a,I2,1x)/)') "Number of processors solving: ",solvgrpbuf(1)," not solving: ",nosolvgrpbuf(1)

      END IF

      !  Sends data to every processors
      CALL MPI_BCAST (solvgrpbuf,npes+1,MPI_INTEGER,masternode,MPI_COMM_WORLD,mpierr)
      CALL MPI_BCAST (nosolvgrpbuf,npes+1,MPI_INTEGER,masternode,MPI_COMM_WORLD,mpierr)

      !  Extract the original group handle
      CALL MPI_COMM_GROUP(MPI_COMM_WORLD,group_world,mpierr)

      !  Creates a new group based on existing group
      numsolvprocs = solvgrpbuf(1); numnosolvprocs = nosolvgrpbuf(1)
      WRITE(*,*) "ID: ",pid,"solvgrpbuf = ", solvgrpbuf(2:1+numsolvprocs)
      WRITE(*,*) "ID: ",pid,"nosolvgrpbuf = ", nosolvgrpbuf(2:1+numnosolvprocs)
      CALL MPI_GROUP_INCL(group_world,numsolvprocs,solvgrpbuf(2:1+numsolvprocs),group_solver,mpierr)
      CALL MPI_GROUP_INCL(group_world,numnosolvprocs,nosolvgrpbuf(2:1+numnosolvprocs),group_nosolver,mpierr)

      !  Creates new communicators
      CALL MPI_COMM_CREATE(MPI_COMM_WORLD,group_solver,mpi_comm_solver,mpierr)
      CALL MPI_COMM_CREATE(MPI_COMM_WORLD,group_nosolver,mpi_comm_nosolver,mpierr)

!      WRITE(6,*) "MPI_COMM_SOLVER = ", mpi_comm_solver

      ! output the solver group table
      ! -- an index table containing [ flag if in solver group, number of processors with non-zero ky modes, their processor IDs ]
      ALLOCATE( grp_solv(1+numsolvprocs) )
!!$      IF ( (kylmax-kylmin+1) .NE. 0 ) THEN       ! check if the processor is part of the solver group
!!$         grp_solv(1) = 1
!!$      ELSE
!!$         grp_solv(1) = -1
!!$      END IF
      grp_solv = solvgrpbuf


END SUBROUTINE SPLIT_COMM_GROUP


SUBROUTINE make_source
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   SUBROUTINE make_source: returns "densource"
!C      -- currently the source is defined only for ky = 0 modes
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE fft_variables, ONLY: kxmax,kymax,nxcmplx,nycmplx,nxreal,nyreal,kxm2
      USE variables, ONLY: densource,kxnorm,samp,sscale,czero
      USE mpi_variables
      IMPLICIT NONE

      INTEGER :: i,j
      COMPLEX :: source1, sink1
      REAL :: x1, x2, amp, scale

      densource = czero; source1 = czero; sink1 = czero
      x1 = 0.25d0; x2 = 0.75d0
      amp = samp; scale = sscale*REAL(kxmax)

      IF (estoy_bien(0,kylmin_solv,kylmax_solv)) THEN           ! check if ky = 0 mode is in the solving range

         DO i = 0, kxmax
            CALL getgauss(source1,kxnorm*REAL(i),1.d0,x1,scale)
            CALL getgauss(sink1,kxnorm*REAL(i),-1.d0,x2,scale)
            densource(i,0) = source1 + sink1
         END DO

      END IF

      densource = 0.5d0*amp*densource

END SUBROUTINE make_source

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
      USE variables, ONLY: pi
      IMPLICIT NONE
      
      REAL, INTENT(in) :: k,mag,x0,sigma
      COMPLEX, INTENT(out) :: sorrow
      REAL :: temp_fact

      sorrow = mag / SQRT(pi*sigma) * CMPLX(COS(2.0d0*pi*k*x0),-1.0d0*SIN(2.0d0*pi*k*x0)) * EXP(-0.25d0*k*k/sigma)
!C                                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ phase due to shift of gaussian curve

END SUBROUTINE getgauss


SUBROUTINE make_extflow
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   SUBROUTINE make_extflow: create an external flow profile
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE mpi_variables, ONLY: kylmin_solv, kylmax_solv, estoy_bien
      USE fft_variables, ONLY: kxmax,kymax
      USE variables, ONLY: phi0,phi0amp,czero,pi,ic
      IMPLICIT NONE

      INTEGER :: kx,ky
      REAL :: x1,x2,coeff
      
      phi0 = czero

      IF (estoy_bien(0,kylmin_solv,kylmax_solv)) THEN
         phi0(1,0) = phi0amp*CMPLX(0.d0,-1.d0)  ! negative sine function
      END IF


END SUBROUTINE make_extflow


SUBROUTINE GATHER_CONJG(fieldin,fieldout)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   SUBROUTINE GATHER_CONJG: gathers conjugate modes on required
!C                      processors.
!C     This is a similar method to GATHER_FIELD except for the number of data.
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE mpi
      USE fft_variables
      USE mpi_variables
!      USE variables
      IMPLICIT NONE

      COMPLEX, DIMENSION(0:kxmax,yinit:yend), INTENT(IN) :: fieldin
      COMPLEX, DIMENSION(0:kxmax,yinit:yend), INTENT(OUT) :: fieldout
      INTEGER :: i,j,ky, kount
      COMPLEX, DIMENSION(nytransfer) :: sendbuf                 ! send same number of data
      COMPLEX, DIMENSION(npes*nytransfer) :: recvbuf            ! receive data and place according to rank
      COMPLEX, DIMENSION(0:nycmplx-1) :: kdat                   ! total number kx = 0 modes
      COMPLEX, dimension(kymax) :: conjg_modes

      sendbuf = 0.d0; recvbuf = 0.d0

      IF (l_ihavedata) THEN      ! fill send buffer with data
         sendbuf(1:nylocal) = fieldin(0,yinit:yend)
      END IF

      ! Gather data from all processors on all processors (perhaps there's a way to split to just communicate within the solver group)
      ! -- Received data is arranged according to the rank of the processor.
      CALL MPI_ALLGATHER(sendbuf,nytransfer,MPI_DOUBLE_COMPLEX, &
           recvbuf,nytransfer,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,mpierr)

      IF (l_ihavedata) THEN
         ! Truncate undesirable modes (1<= ky <= kymax), and store the desired modes (-kymax <= ky <= -1).
         ! Sorts out the data according to index directory
         kount = 0
         DO i = 0, npes-1                       ! iterate through number of processors
            DO j = 1, localindex_dir(1,i)       ! move starting position by amount of nytransfer and end at number of elements on each processor (nylocal)
               kdat(kount) = recvbuf(i*nytransfer+j)
               kount = kount + 1
            END DO
         END DO

         IF (kount .NE. nycmplx) THEN
            WRITE(*,*) "kount = ", kount, " nycmplx = ", nycmplx
            WRITE(6,'(/a,I4)') "ID: ", pid
            STOP "Forced stop in GATHER_CONJG. Data accumulation failure."
         END IF

         ! Calculate the conjugate modes.
         ! Format of "kdata" => [0,1,...,kymax,...,-kymax,...,-1]. The order is reversed to calculated conjugate modes.
         conjg_modes(1:kymax) = CONJG(kdat(nycmplx-1:nycmplx-kymax:-1))

         ! Fill in conjugate modes for processors holding the conjugate modes (1 <= ky <= kymax)
         DO ky = 1, kymax
            IF (estoy_bien(ky,kylmin_solv,kylmax_solv)) THEN
               fieldout(0,yext(ky)) = conjg_modes(ky)
            END IF
         
         END DO

      END IF
      
END SUBROUTINE GATHER_CONJG


SUBROUTINE GATHER_FIELD(fieldin,fieldout,rootid,gfflag)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   SUBROUTINE GATHER_FIELD: collect segmented data and broadcast the complete field
!C     This is a similar method to GATHER_DATA except the broadcasted result.
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE mpi
      USE fft_variables
      USE mpi_variables
      IMPLICIT NONE

      INTEGER :: i, j, kount
      INTEGER, INTENT(IN) :: rootid, gfflag
      COMPLEX, DIMENSION(0:kxmax,yinit:yend), INTENT(IN) :: fieldin
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax), INTENT(OUT) :: fieldout
      COMPLEX, DIMENSION(kxm2,npes*nytransfer) :: recvbuf   ! <-- use to receive data
      COMPLEX, DIMENSION(0:kxmax,0:nycmplx-1) :: kdat       ! <-- use to sort data
      COMPLEX, DIMENSION(kxm2,nytransfer) :: sendbuf        ! <-- use to transfer data
      COMPLEX, PARAMETER :: czero = CMPLX(0.0,0.0)

      recvbuf = czero; kdat = czero; sendbuf = czero

      IF (l_ihavedata) THEN      ! fill the send array with data
         sendbuf(1:kxm2,1:nylocal) = fieldin(0:kxmax,yinit:yend)
      END IF

      SELECT CASE (gfflag)
      CASE (0)
         ! Collects data on to node "rootid"
         CALL MPI_GATHER(sendbuf,kxm2*nytransfer,MPI_DOUBLE_COMPLEX, &
              recvbuf,kxm2*nytransfer,MPI_DOUBLE_COMPLEX,rootid,MPI_COMM_WORLD,mpierr)
         
      CASE (1)
         ! Gather data from all processors on all processors
         ! -- Received data is arranged according to the rank of the processor.
         CALL MPI_ALLGATHER(sendbuf,kxm2*nytransfer,MPI_DOUBLE_COMPLEX, &
              recvbuf,kxm2*nytransfer,MPI_DOUBLE_COMPLEX,MPI_COMM_WORLD,mpierr)

      END SELECT

      ! sorts out the data according to index directory
      kount = 0
      DO i = 0, npes-1
         DO j = 1, localindex_dir(1,i)       ! move starting position by amount of nytransfer and end at number of elements on each processor
            kount = kount + 1
            kdat(0:kxmax,kount-1) = recvbuf(1:kxm2,i*nytransfer+j)
         END DO
      END DO

      IF (kount .NE. nycmplx) STOP "Forced stop in GATHER_DATA. Data accumulation failure."

      fieldout = czero
         
      fieldout(0:kxmax,-kymax:-1) = kdat(0:kxmax,nycmplx-kymax:nycmplx-1)   ! obtain [0 -> kxmax] and [-kymax -> -1] modes
      fieldout(0:kxmax,0:kymax) = kdat(0:kxmax,0:kymax)                     ! obtain [0 -> kxmax] and [0 -> kymax] modes

      ! Fill in the conjugate modes...
      fieldout(0,1:kymax) = CONJG( fieldout(0,-1:-kymax:-1) )

END SUBROUTINE GATHER_FIELD


SUBROUTINE GATHER_DATA(rootid,gdflag)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   SUBROUTINE GATHER_DATA: collect segmented data and broadcast the complete field
!C     depending on parameter "mflag"
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE mpi_variables, ONLY: pid,blacksheep
      USE fft_variables, ONLY: kxmax,kymax
      USE variables, ONLY: den,phi,denavg,denout,phiout,denavgout,phitr
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: rootid, gdflag
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax) :: cjunk

      SELECT CASE (gdflag)
      CASE (0)
         ! Collects all fields on to the "rootid"
         IF (pid .EQ. blacksheep) THEN
            CALL GATHER_FIELD(den,denout,rootid,0)
            CALL GATHER_FIELD(phi,phiout,rootid,0)
            CALL GATHER_FIELD(denavg,denavgout,rootid,0)
         ELSE
            CALL GATHER_FIELD(den,cjunk,rootid,0)
            CALL GATHER_FIELD(phi,cjunk,rootid,0)
            CALL GATHER_FIELD(denavg,cjunk,rootid,0)
         END IF
      CASE (1)
         ! Collects and broadcast "phi" field
         CALL GATHER_FIELD(phi,phitr,-1,1)   ! "-1" for ID is to introduce an error
      END SELECT

END SUBROUTINE GATHER_DATA


SUBROUTINE WROUT
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   SUBROUTINE WROUT: calls other subroutines to write outputs
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT NONE

      CALL energy                             ! calculates energy
      CALL writeit(10)                        ! write to energy file
      CALL writerestart                       ! write "restart" file
      CALL write_visit                        ! write output to ViSit

END SUBROUTINE WROUT


SUBROUTINE ALLOC_DATA
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   SUBROUTINE ALLOC_DATA: allocates arrays
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE variables
      USE fft_variables
      USE mpi_variables
      USE prog_flags, ONLY: lp,le,lm,l_tracers
      USE convolve_module
      IMPLICIT NONE
      
      ! all processors need to have a subset of data
      ALLOCATE( den(0:kxmax,yinit:yend) ); den = czero  ! subset modes of the complete extended array
      ALLOCATE( phi(0:kxmax,yinit:yend) ); phi = czero
      ALLOCATE( denavg(0:kxmax,yinit:yend) ); denavg = czero

      ! allocate dependent data -- allocate for only processors that contain local data
      IF (l_ihavedata) THEN
         ALLOCATE( philinarr(0:kxmax,yinit:yend) ); philinarr = czero     ! linear arrays are allocated to match size of data
         ALLOCATE( philindrivearr(0:kxmax,yinit:yend) ); philindrivearr = czero
         ALLOCATE( damplinarr(0:kxmax,yinit:yend) ); damplinarr = czero
         ALLOCATE( drivelinarr(0:kxmax,yinit:yend) ); drivelinarr = czero
         ALLOCATE( dampdenavglinarr(0:kxmax,yinit:yend) ); dampdenavglinarr = czero
         ALLOCATE( densource(0:kxmax,yinit:yend) ); densource = czero
         ALLOCATE( phi0(0:kxmax,yinit:yend) ); phi0 = czero

         ALLOCATE( tola(nn_local) ); tola = 0.d0                     ! tolerance for total modes being solved
         ALLOCATE( tola_serial(nn_serial) ); tola_serial = 0.d0

      END IF

      ! allocate arrays for output
      IF (pid .EQ. blacksheep) THEN 
         ALLOCATE( denout(0:kxmax,-kymax:kymax) )  ! stores half of complete array
         ALLOCATE( phiout(0:kxmax,-kymax:kymax) )
         ALLOCATE( denavgout(0:kxmax,-kymax:kymax) )
      END IF

      ! allocate arrays for the total phi field to compute tracer propagation
      IF (l_tracers) THEN
         ALLOCATE( phitr(0:kxmax,-kymax:kymax) )
      END IF

      ! allocate arrays for FFTW -- all processors must have locally allocated space for FFT
      IF ( (lp .NE. 0.0) .OR. (le .NE. 0.0)  ) THEN
         ALLOCATE( cvex_flat(ncmplxfft) ); ALLOCATE( cvey_flat(ncmplxfft) )

         IF (lp .NE. 0.0) THEN   ! if Polarization non-linearity flag is on
            ALLOCATE( cwx_flat(ncmplxfft) ); ALLOCATE( cwy_flat(ncmplxfft) )
            ALLOCATE(Pol_flat(local_size))
            IF (pid .EQ. blacksheep) WRITE(6,'(/a/)') "Polarization non-linearity enabled."
         END IF
         IF (le .NE. 0.0) THEN   ! if ExB non-linearity flag is on
            ALLOCATE( cnky_flat(ncmplxfft) ); ALLOCATE( cnkx_flat(ncmplxfft) )
            ALLOCATE(ExB_flat(local_size))
            IF (pid .EQ. blacksheep) WRITE(6,'(/a/)') "ExB non-linearity enabled."
         END IF
         IF (lm .NE. 0.0) THEN   ! if mean field flag is on
            ALLOCATE( cnavgky_flat(ncmplxfft) ); ALLOCATE( cnavgkx_flat(ncmplxfft) )
            ALLOCATE( cnavgkx1d_flat(ncmplxfft) ); ALLOCATE( cnavg1d_flat(ncmplxfft) )
            ALLOCATE( cnavg_flat(ncmplxfft) )
            ALLOCATE(denavgnl_flat(local_size)); ALLOCATE(drivenl_flat(local_size))
            IF (pid .EQ. blacksheep) WRITE(6,'(/a/)') "ExB non-linearity on mean field enabled."
         END IF

      ELSE

         IF (pid .EQ. blacksheep) WRITE(6,'(/a/)') "Running without non-linear terms."

      END IF

END SUBROUTINE ALLOC_DATA


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
      USE mpi_variables, ONLY: pid, blacksheep
      USE variables, ONLY: t,ncl,ncl_in,denekold,denenold,phiekold,phiek1old,phienold
      USE tracers, ONLY: itrack
      USE file_units
      USE visit, ONLY: ivisit
      IMPLICIT NONE

      LOGICAL, INTENT(OUT) :: lcont

      ! temporary variables for read in
      REAL :: time0,dt0
      INTEGER :: ncl0,ivisit0,nmax0,itrack0,iflightsx0,iflightsy0,lypmaxdtstep0

      ! read in NAMELIST
      NAMELIST /contp/ time0,dt0,ncl0,ivisit0
      NAMELIST /tracerp/ nmax0,itrack0

      ! initialize energy variables
      denekold = 0.0; denenold = 0.0
      phiekold = 0.0; phiek1old = 0.0; phienold = 0.0

      INQUIRE(file='params_out',exist=lcont)
      IF (lcont .EQV. .FALSE.) THEN   ! open new files if this is a fresh run

!C  OPEN files for READ and WRITE --------------------------------------

         IF (pid .EQ. blacksheep) THEN

            OPEN(u_energy,file='energy',status='replace')
            WRITE(u_energy,100)     ! writes energy file header
100         FORMAT(4x,'ncl',10x,'t',10x,'denek',10x,'denen',10x,'dengamek',10x,'dengamen',10x, &
                 'phiek',10x,'phiek1',10x,'phien',10x,'phigamek',10x,'phigamen',10x,'denavgek',10x,'ek')
            OPEN(u_inputscopy,file='inputs0_copy',status='replace')

         END IF

         ! initialize variables
         t = 0.0; ncl = 0; ncl_in = 0; itrack = 0; ivisit = 0

      ELSE

         ! read continuation parameters
         OPEN(21,file='params_out',status='old')
         READ(21,contp)
         READ(21,tracerp)
         CLOSE(21)

         t = time0; ncl = 0; ncl_in = ncl0; itrack = itrack0; ivisit = ivisit0

         IF (pid .EQ. blacksheep) THEN

            WRITE(6,'(/a)') " ** This is a continuation run. **"
            WRITE(6,*) "time = ", t
            WRITE(6,*) "ncl_in = ", ncl_in, " itrack = ", itrack, " ivisit = ", ivisit

            OPEN(u_energy,file='energy',status='old',position='append')

         END IF

      END IF

!C  End OPEN files -----------------------------------------------------

      ! blacksheep prints solver's performance
      IF (pid .EQ. blacksheep) OPEN(95,file='dvodpk_log',status='replace')   ! <-- to track solver's performance

END SUBROUTINE check_cont

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
      USE tracers, ONLY: nmax, itrack
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
      WRITE(602,*) "ITRACK0 = ", itrack, " /"
      WRITE(602,*)

      ! extract number of modes
      WRITE(602,*) "&FFTP"
      WRITE(602,*) "KXMAX = ", kxmax, " ,"
      WRITE(602,*) "KYMAX = ", kymax, " ,"
      WRITE(602,*) "NXREAL = ", nxreal, " ,"
      WRITE(602,*) "NYREAL = ", nyreal, " /"

      CLOSE(602)

END SUBROUTINE make_namelist

!C BEGIN TRACER SECTION ********************************************************************************************************************************************************************

SUBROUTINE tracers_init
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C     subroutine tracers_init
!C
!C     initialize tracers positions
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE variables, ONLY: t
      USE mpi
      USE tracers
      USE mpi_variables
      USE prog_flags, ONLY: lcont
      USE file_units
      IMPLICIT NONE

      CHARACTER(LEN=4) :: pid_str

      !  Divide the number of tracers to processors that do not have data
      IF ( nmax .LT. npes) THEN
         nmax_local = 1                    ! force to have one particle per procesor
         IF (pid .EQ. blacksheep) THEN 
            WRITE(6,'(/a/)') "Forcing number of particles per node = 1 since NMAX < NPES."
         END IF
      ELSE
         IF ( (pid .EQ. (npes-1)) .AND. (MOD(nmax,npes) .NE. 0) ) THEN
            nmax_local = INT(nmax/npes) + MOD(nmax,npes)
         ELSE
            nmax_local = INT(nmax/npes)
         END IF
         WRITE(6,'(2(a,I4))') "ID: ",pid," number of tracers = ", nmax_local
      END IF

      ALLOCATE( px0(nmax_local,2),py0(nmax_local,2),px(nmax_local,2),py(nmax_local,2) )         ! tracers are paired
      ALLOCATE( pvx0(nmax_local,2),pvx(nmax_local,2),pvy(nmax_local,2),pvy0(nmax_local,2) )
      ALLOCATE( pdistinit(nmax_local) )

      CALL get_iterstr(pid_str,pid,LEN(pid_str))
      IF (lcont) THEN
         told = t
         OPEN(unit=u_tracers,file='tposvel'//TRIM(pid_str)//'.bin',status='old',position='append',form='unformatted')
      ELSE
         told = 0.0
         OPEN(unit=u_tracers,file='tposvel.'//TRIM(pid_str)//'bin',status='replace',form='unformatted')
      END IF

      CALL pinit

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
      USE file_units
      USE variables, ONLY: t
      IMPLICIT NONE

      CALL pwritecont(51)   ! write continuation file

      DEALLOCATE(px0,py0,px,py,pvx,pvy,pvx0,pvy0)
      CLOSE(u_tracers)

END SUBROUTINE tracers_final


SUBROUTINE pinit 
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  subroutine PINIT                        
!C
!C  Initialize the particles according to NVPART
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE tracers
      USE variables, ONLY: irand, t
      USE prog_flags, ONLY: lcont
      USE file_units

      IMPLICIT NONE

      LOGICAL :: lexist
      INTEGER :: i,n
      REAL, EXTERNAL :: ran2
      REAL, PARAMETER :: small = TINY(0.0d0)
      REAL, DIMENSION(nmax_local) :: dx0,dy0
      REAL :: tdump
      CHARACTER(LEN=128) :: tcontfile_str

      IF (lcont .EQV. .FALSE.) THEN

         SELECT CASE (nvpart)
         CASE (1)   ! Initializes tracers at random positions

            WRITE(6,*) '  Tracers equally spaced in a box. Number of tracers: ', nmax
            WRITE(u_inputscopy,*) '  Tracers equally spaced in a box. Number of tracers: ', nmax

            DO n = 1, nmax_local        ! loop over number of tracers
               px0(n,1) = ran2(irand+n)*boxlen
               py0(n,1) = ran2(irand+n*n)*boxlen
               px0(n,2) = px0(n,1) + smalldist*boxlen ! starting another tracer nearby, test tracers
               py0(n,2) = py0(n,1)             ! same y-position
            END DO

            px = px0; py = py0 ! assign initial positions
            pvx = 0.0; pvy = 0.0  ! no initial velocities

         CASE DEFAULT

            STOP "WARNING: NVPART needs to defined properly!"

         END SELECT

         dx0 = px0(1:nmax_local,1) - px0(1:nmax_local,2); dy0 = py0(1:nmax_local,1) - py0(1:nmax_local,2)
         pdistinit = SQRT( dx0*dx0 + dy0*dy0 )

      ELSE

         WRITE(6,*) "  Read in tracers' positions from continuation files."
         CALL preadcont(51)
         
      END IF

END SUBROUTINE pinit

SUBROUTINE tracers_advance(dt)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C  subroutine padvance
!C
!C   Advance the tracers.
!C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE variables, ONLY: phi
      USE fft_variables
      USE tracers
      USE prog_flags, ONLY: nvisit
      IMPLICIT NONE

      INTEGER :: n, j
      REAL :: vxt, vyt
      REAL, INTENT(in) :: dt
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax) :: phifield

!      WRITE(*,*) "Send phi field for tracers."
      CALL GATHER_FIELD(phi,phifield,-1,1)

      DO j = 1, 2         ! loop over pairs
         DO n = 1, nmax_local   ! loops over number of tracers
            
            vxt = 0.d0; vyt = 0.d0

            CALL pvelocity(vxt,vyt,px(n,j),py(n,j),phifield)             ! find velocity at tracer position
            px(n,j) = px(n,j) + vxt*boxlen*dt
            py(n,j) = py(n,j) + vyt*boxlen*dt
            pvx(n,j) = vxt; pvy(n,j) = vyt
         
         END DO
      END DO

END SUBROUTINE tracers_advance

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
      USE fft_variables
      USE variables, ONLY: twopi,ic,kxnorm,kynorm
      USE tracers, ONLY: boxlen
      IMPLICIT NONE

      INTEGER :: kx, ky
      REAL :: xk, yk, spx, spy
      REAL, INTENT(in) :: xt,yt
      REAL, INTENT(out) :: vxt,vyt
      COMPLEX, DIMENSION(0:kxmax,-kymax:kymax), INTENT(in) :: phifield

!C  Initialization
 
      vxt = 0.0; vyt = 0.0

!C   Normalize particle infinite domain to periodic domain of the field

      IF (xt .GT. 0.0) spx = MOD(xt,boxlen)                           ! modified implementation
      IF (xt .LE. 0.0) spx = boxlen + MOD(xt,boxlen)
      IF (yt .GT. 0.0) spy = MOD(yt,boxlen)
      IF (yt .LE. 0.0) spy = boxlen + MOD(yt,boxlen)
      spx = twopi*spx; spy = twopi*spy

!C  Calculate velocities by adding over the Fourier series

      DO ky = -kymax, kymax
         yk = kynorm*REAL(ky)

         DO kx = 0, kxmax
            xk = kxnorm*REAL(kx)

            vxt = vxt - REAL(ic*yk*phifield(kx,ky)*EXP(ic*(xk*spx+yk*spy)))
            vyt = vyt + REAL(ic*xk*phifield(kx,ky)*EXP(ic*(xk*spx+yk*spy)))

         END DO
      END DO

      vxt = 2.0d0*vxt; vyt = 2.0d0*vyt   ! omitting the division by the total number of grid points since the factor is simply a scaling factor

END SUBROUTINE pvelocity

SUBROUTINE tracers_write(t)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine TRACERS_WRITE
!C
!C   Writes the particle positions at time t.
!C   Must also enter file units. This is a binary file.
!C      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE tracers
      USE file_units
      IMPLICIT NONE

      REAL, INTENT(IN) :: t

      itrack = itrack + 1      ! specify as one entry
      WRITE(u_tracers) itrack,t,px,py,pvx,pvy

END SUBROUTINE tracers_write

SUBROUTINE preadcont(funit)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine preadcont
!C
!C   Reads from continuation tracer file. This is a binary file.
!C      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE tracers
      USE mpi_variables, ONLY: pid
      IMPLICIT NONE

      INTEGER, INTENT(in) :: funit
      CHARACTER(LEN=4) :: pid_str

      CALL get_iterstr(pid_str,pid,LEN(pid_str))
      OPEN(unit=funit,file="tcont"//TRIM(pid_str)//".bin",status="old",form="unformatted")
      READ(funit) px0,py0,px,py                       ! read particle positions
      CLOSE(funit)

END SUBROUTINE preadcont

SUBROUTINE pwritecont(funit)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!C
!C   subroutine pwritecont
!C
!C   Writes continuation parameters. This is a binary file.
!C      
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      USE tracers
      USE mpi_variables, ONLY: pid
      IMPLICIT NONE

      INTEGER, INTENT(in) :: funit
      INTEGER :: n
      CHARACTER(LEN=4) :: pid_str

      CALL get_iterstr(pid_str,pid,LEN(pid_str))
      OPEN(unit=funit,file="tcont"//TRIM(pid_str)//".bin",status="replace",form="unformatted")
      WRITE(funit) px0,py0,px,py                       ! write particle positions
      CLOSE(funit)

END SUBROUTINE pwritecont

!C END TRACER SECTION ********************************************************************************************************************************************************************
