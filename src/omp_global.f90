!==========================================================================================
! Modification notes
! ------------------
! Since executing the runs for the ICB DCU paper:
! 26/08/09  Changed DC_DCprox from 1.0 to 1.5 (DCs were too close together)
! 27/08/09
! try turning off use_DCflux, and set fix_vascularity, to test steadystate case
! added correction to exit_prob in steadystate case
! try turn_off_chemotaxis with TC_TO_DC = 0, TC_COGNATE_FRACTION = 0 and log_results = .false.
! (has the effect of pushing NTcells to 32500, i.e. increasing Tres by 1.3 -> 31 hrs)
! Problem: with chemotaxis I'm seeing Tres = 43 hrs when it should be 12!!!
! Fixed
!
! 01/10/09
! Error: noticed that TC_STIM_HALFLIFE was not being used!!! All previous results were in
! fact based on an infinite halflife of TCR stimulation (no decay).
! Fixed
!
! 02/10/09
! Error: noticed that TC_STIM_WEIGHT can be < 1 when use_cytokines is false (no IL-2 effect)
! This needs to be fixed: set TC_STIM_WEIGHT = 1 when .not.use_cytokines.
! Fixed
! All cases will need to be rerun!
!
! 14/10/09
! For expts to simulate Henrickson2008, make T cells live longer:
! life_median1 = life_median2 = 96
!
! 19/10/2009
! Changed DC%density to represent the number of pMHC/DC.
! The threshold level for any TCR signalling was set DC_pMHC_THRESHOLD = 30
! on the basis of Henrickson2008.
! When TCR signalling is turned on, the rate of signalling is k.A.D
! where A is the TCR avidity and D is the pMHC count.
! Normalize avidity based on M-peptide as having a median value of 1, then the rate
! constant k must be adjusted (start at 1/30).
!
! 20/10/2009
! Added density_min = DC_pMHC_THRESHOLD/2 to subroutine updater.
! The idea is to reduce stimulation rate as the threshold is approached
!==========================================================================================

module omp_global

use omp_lib
use omp_CD69
use omp_IL2
use omp_IL7
use omp_IL_dummy
use par_zig_mod
use winsock

implicit none

INTEGER,  PARAMETER  ::  DP=SELECTED_REAL_KIND( 12, 60 )
integer, parameter :: REAL_KIND = 4

! General parameters
integer, parameter :: BIG_INT = 2**30
integer, parameter :: NEUMANN_MODEL = 1
integer, parameter :: MOORE18_MODEL = 2
integer, parameter :: MOORE26_MODEL = 3
integer, parameter :: MODEL = MOORE18_MODEL
integer, parameter :: MAXRELDIR = 26
integer, parameter :: TRAFFIC_MODE_1 = 1
integer, parameter :: TRAFFIC_MODE_2 = 2
integer, parameter :: EXIT_EVERYWHERE = 1
integer, parameter :: EXIT_LOWERHALF = 2
!integer, parameter :: EXIT_CHEMOTAXIS = 3
integer, parameter :: EXIT_BLOB_PORTALS = 3
integer, parameter :: EXIT_SURFACE_PORTALS = 4
integer, parameter :: EXIT_NONE = 5
integer, parameter :: TASK_TAG = 1
integer, parameter :: CELL_TAG = 2
integer, parameter :: LOC_TAG  = 3
integer, parameter :: OCC_TAG  = 4
integer, parameter :: OCNT_TAG = 5
integer, parameter :: CCNT_TAG = 6
integer, parameter :: DC_TAG   = 7
integer, parameter :: ADD_TAG   = 8
integer, parameter :: CYT_L2R_TAG   = 10
integer, parameter :: CYT_R2L_TAG   = 11
integer, parameter :: GLOBAL1_TAG   = 12
integer, parameter :: GLOBAL2_TAG   = 13
integer, parameter :: DC1_TAG   = 14
integer, parameter :: DC2_TAG   = 15
integer, parameter :: DC3_TAG   = 16
integer, parameter :: DC4_TAG   = 17
integer, parameter :: DC5_TAG   = 18
integer, parameter :: WX_TAG   = 19
integer, parameter :: CYT_INFO_TAG   = 20
integer, parameter :: CYT_DATA_TAG   = 21
integer, parameter :: RES1_TAG   = 22
integer, parameter :: RES2_TAG   = 23
integer, parameter :: RES3_TAG   = 24
integer, parameter :: RES4_TAG   = 25
integer, parameter :: RES5_TAG   = 26
integer, parameter :: VIS1_TAG   = 27
integer, parameter :: VIS2_TAG   = 28
integer, parameter :: MOL1_TAG   = 29
integer, parameter :: MOL2_TAG   = 30

! Boundary types
integer, parameter :: CONC_BC = 1
integer, parameter :: RATE_BC = 2

! T cell region
integer, parameter :: LYMPHNODE = 1
integer, parameter :: PERIPHERY = LYMPHNODE + 1

! Site identifiers for proximity limits
integer, parameter :: BDRY_SITE = 1
integer, parameter :: EXIT_SITE = 2
integer, parameter :: HEV_SITE = 3
integer, parameter :: DC_SITE = 4

! T cell activation stage
integer, parameter :: NAIVE      = 1
integer, parameter :: TRANSIENT  = 2
integer, parameter :: CLUSTERS   = 3
integer, parameter :: SWARMS     = 4
integer, parameter :: DIVIDING   = 5
integer, parameter :: FINISHED   = 6
integer, parameter :: STAGELIMIT = 6
!integer, parameter :: ACTIVATED  = 5   ! NOT USED NOW

integer, parameter :: TCP_PORT_0 = 5000		! main communication port (logging) 
integer, parameter :: TCP_PORT_1 = 5001		! data transfer port (plotting)
integer, parameter :: TCP_PORT_2 = 5002
integer, parameter :: TCP_PORT_3 = 5003

integer, parameter :: CT_CONSTANT = 1
integer, parameter :: CT_HILL = 2
integer, parameter :: CT_HENRICKSON = 3

integer, parameter :: STAGE_BYTE = 1
integer, parameter :: GENERATION_BYTE = 2
integer, parameter :: SLOT_NUM1 = 1
integer, parameter :: SLOT_NUM2 = 2
integer, parameter :: BOTH = 3
integer, parameter :: IL2_TAG  = 1
integer, parameter :: IL4_TAG  = 2
integer, parameter :: IL7_TAG  = 3
integer, parameter :: IL9_TAG  = 4
integer, parameter :: IL15_TAG = 5
integer, parameter :: IL21_TAG = 6
integer, parameter :: MAX_CYT = 6
character*(5), parameter :: cyt_name(MAX_CYT) = (/ 'IL-2 ','IL-4 ','IL-7 ','IL-9 ','IL-15','IL-21' /)

integer, parameter :: NCTYPES = 4
integer, parameter :: NONCOG_TYPE_TAG  = 1
integer, parameter :: COG_TYPE_TAG  = 2
integer, parameter :: COG_CD4_TAG  = 2
integer, parameter :: COG_CD8_TAG  = 3
integer, parameter :: COG_TREG_TAG  = 4

! Revised system
integer, parameter :: CD4 = 1
integer, parameter :: CD8 = 2

integer, parameter :: TAGGED_CELL = 100
integer, parameter :: RES_TAGGED_CELL = 101
integer, parameter :: CHEMO_TAGGED_CELL = 102
integer, parameter :: OUTSIDE_TAG = -BIG_INT + 1

integer, parameter :: neumann(3,6) = reshape((/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /), (/3,6/))
integer, parameter :: jumpvec2D(3,8) = reshape((/ 1,0,0, 1,1,0, 0,1,0, -1,1,0, -1,0,0, -1,-1,0, 0,-1,0, 1,-1,0 /), (/3,8/))
integer, parameter :: nfcell = 10, nfout = 11, nfvec = 12, nfpath = 13, nfres = 14, nfdcbind = 15, nftraffic = 16, nfrun = 17, &
					  nftravel=18, nfcmgui=19, nfpos=20, nflog=21, nfchemo=22, nfDCinfo=23
real, parameter :: BIG_TIME = 100000
real, parameter :: BALANCER_INTERVAL = 10
integer, parameter :: SCANNER_INTERVAL = 100
real, parameter :: DELTA_T = 0.25       ! minutes
logical, parameter :: use_add_count = .true.    ! keep count of sites to add/remove, do the adjustment at regular intervals 
logical, parameter :: save_input = .true.

!logical, parameter :: IN_VITRO = .true.
integer, parameter :: NZ_IN_VITRO = 3
integer, parameter :: NDCcore_IN_VITRO = 8
!real, parameter :: VITRO_FRACTION = 0.3		! temporary measure, for testing 
real, parameter :: STACKRAD_2 = 3.5
real, parameter :: STACKRAD_3 = 2.3


! DC parameters
logical, parameter :: DC_motion = .false.
logical, parameter :: RANDOM_DCFLUX = .false.
integer, parameter :: DCDIM = 4         ! MUST be an even number
integer, parameter :: cDCDIM = 4        ! MUST be an even number
integer, parameter :: NDCcore = 7      ! In fact DC vol = 1400 = 5.6*250
real, parameter :: DC_DCprox = 1.3      ! closest placement of DCs, units DC_RADIUS (WAS 1.0 for ICB DCU paper)
real, parameter :: bdry_DCprox = 0.5	! closest placement of DC to bdry, units DC_RADIUS
real, parameter :: exit_DCprox = 4.0    ! closest placement of DC to exit, units sites
real, parameter :: exit_prox = 1.0      ! closest placement of exit to exit, units chemo_radius
! closest placement of exit to bdry?
logical, parameter :: incapable_DC_dies = .false.	! don't want to remove low-antigen DCs prematurely
logical, parameter :: reuse_DC_index = .false.	! index reuse conflicts with the GUI animation code

! Diffusion parameters
logical, parameter :: use_cytokines = .false.       ! to use IL-2
logical, parameter :: use_diffusion = .false.		! For cytokines
integer, parameter :: NDIFFSTEPS = 6    ! divisions of DELTA_T for diffusion computation

! T cell parameters
integer, parameter :: traffic_mode = TRAFFIC_MODE_2
logical, parameter :: use_blob = .true.
logical, parameter :: use_HEV_portals = .true.
logical, parameter :: random_cognate = .true.	! number of cognate seed cells is random or determined
integer, parameter :: MMAX_GEN = 20     ! max number of generations (for array dimension only)
integer, parameter :: NGEN_EXIT = 8     ! minimum non-NAIVE T cell generation permitted to exit (exit_rule = 1)
real, parameter :: CHEMO_MIN = 0.05		! minimum level of exit chemotactic influence (at r = chemo_radius)
integer :: exit_rule = 1                ! 1 = use NGEN_EXIT, 2 = use EXIT_THRESHOLD, 3 = use S1P1
logical :: USE_S1P = .true.				! this is the default
logical :: COMPUTE_OUTFLOW = .false.
real :: CTYPE_FRACTION(2)

! Chemokine/receptor parameters
integer, parameter :: MAX_CHEMO = 1
integer, parameter :: MAX_RECEPTOR = 1
integer, parameter :: CCL3 = 1
integer, parameter :: CCR1 = 1
logical, parameter :: USE_CHEMOKINE_GRADIENT = .true.	! with use_DC_chemotaxis
logical, parameter :: USE_CELL_SITES = .true.		! specify concentrations at the DC sites
logical, parameter :: USE_GENERAL_CODE = .true.	! this code covers dynamic and SS solvers, conc and secretion b.c.s
logical, parameter :: use_ODE_diffusion = .false.	! for general case, use ODE system (dynamic solver)
logical, parameter :: USE_DC_SECRETION = .false.
logical, parameter :: USE_ORIGINAL_CODE = .not.USE_GENERAL_CODE	! simple case, fixed secretion into DC neighborhood

!specify bdry concentrations, not secretion rates

! Data above this line almost never change
!==============================================================================================================

! Run parameters

! These parameters are used only if exit_region = EXIT_LOWERHALF
logical, parameter :: constant_efactor = .true.
real, parameter :: etheta = 1./25000.
real, parameter :: ZIN_FRACTION = 1.0   ! fraction of blob radius for traffic inflow in upper hemisphere
!real, parameter :: ZOUT_FRACTION = 1.0  ! not used
! This method of excess adjustment is no longer used
!real, parameter :: excess_factor0 = 11.495
!real, parameter :: excess_factor1 = 0.0158e-3
!real, parameter :: excess_factor2 = -8.750


!logical, parameter :: vary_vascularity = .true.
!logical, parameter :: use_chemotaxis = .false. ! now based on exit_region == EXIT_CHEMOTAXIS
!logical, parameter :: fix_avidity = .true.
!integer, parameter :: avidity_nlevels = 8
!!logical, parameter :: avidity_logscale = .true.
!!real, parameter :: avidity_min = -1    ! if avidity_logscale then avidity_min is log10(actual min)
!!real, parameter :: avidity_step = 0.185861   ! and log10(actual max) = log10(actual min) + (nlevels-1)*avidity_step
!logical, parameter :: avidity_logscale = .false.
!real, parameter :: avidity_min = 0.3
!real, parameter :: avidity_step = 0.2

!logical, parameter :: use_DC_chemotaxis = .true.
integer, parameter :: MAX_EX_CHEMO = 4, MAX_DC_CHEMO = 6	! max allowed influences on a cell
!real, parameter :: chemo_K_DC = 2.0
real, parameter :: BASE_DCchemo = 0.0

! Parameters and switches for calibration
logical, parameter :: calibrate_motility = .false.
logical, parameter :: motility_param_range = .false.
logical, parameter :: motility_save_paths = .false.
logical, parameter :: calibrate_diffusion = .false.
logical, parameter :: compute_travel_time = .false.
integer, parameter :: n_multiple_runs = 1

! Parameters and switches for testing
logical, parameter :: test_vascular = .false.
logical, parameter :: turn_off_chemotaxis = .false.		! to test the chemotaxis model when cells are not attracted to exits
logical, parameter :: L_selectin = .false.				! T cell inflow is suppressed - to simulate Franca's experiment

! Debugging parameters
!logical, parameter :: dbug = .false.
logical, parameter :: dbug_cog = .false.
logical, parameter :: avid_debug = .false.
logical, parameter :: debug_DCchemotaxis = .false.
integer, parameter :: CHECKING = 0
integer :: idbug = 0

real, parameter :: TAGGED_CHEMO_ACTIVITY = 1.0

! To investigate the effect of chemotaxis on residence time.
! The situation to be simulated is one in which most cells are not subject to exit chemotaxis,
! but a small fraction of tagged cells are.  The question is: what is the effect on the residence time
! of the tagged cells?
logical, parameter :: TAGGED_EXIT_CHEMOTAXIS = .false.	! ==> evaluate_residence_time = .true.
logical, parameter :: evaluate_residence_time = .false.
integer, parameter :: istep_res1 = 4*60*24*3			! 3 days (was 5000)
integer, parameter :: istep_res2 = istep_res1 + 4*60*24	! 1 day of tagging

! To investigate the effect of chemotaxis on DC contacts.
! The situation to be simulated is one in which most cells are not subject to DC chemotaxis,
! but a small fraction of tagged cells are.  The question is: what is the effect on the DC visit frequency
! of the tagged cells?
! Modify this to read the parameters from an input file.
logical :: track_DCvisits = .true.
logical :: TAGGED_DC_CHEMOTAXIS = .true.
integer :: istep_DCvisits = 24*60*4	! 24 hours delay before counting visits (when use_traffic)
real :: t_log_DCvisits = 2*24*60		! log DC visits of retained cells for 2 days
logical :: USE_DC_COGNATE = .true.	! if only a fraction (0.5) of DCs bear cognate antigen, secrete chemokine
real :: DC_COGNATE_FRACTION = 0.5
logical :: retain_tagged_cells = .true.
real :: RETAIN_TAGLIMIT_FRACTION = 0.5
logical :: DC_CHEMO_NOTRAFFIC = .true.	! cells are tagged initially, no traffic
real :: TAGGED_CHEMO_FRACTION = 0.1		! was 0.1 for exit chemotaxis
real :: DC_CHEMO_FRACTION = 0.2		! fraction of cells tagged when DC_CHEMO_NOTRAFFIC
real :: HI_CHEMO_FRACTION = 0.5		! fraction of tagged cells with full chemotactic susceptibility
real :: HI_CHEMO = 1.0				! chemotactic susceptibility of the HI_CHEMO_FRACTION cells
real :: LO_CHEMO = 0.0				! chemotactic susceptibility of the other tagged cells

! To investigate how chemotaxis influences T cell paths.
! Set up a single exit at the centre of the blob.  Tag cells starting at a specified distance from the exit,
! and log their paths until they exit.  Cells are either chemotactic or not.
logical, parameter :: TAGGED_LOG_PATHS = .false.
integer, parameter :: MAX_LOG_PATHS = 100		! number of paths to track for non-chemo (1) and chemo (2) cells
integer, parameter :: MAX_LOG_PATH_SITES = 1000		! max number of tracking start locations
!integer :: log_path_list(MAX_LOG_PATHS,2)		! lists of IDs of tracked cells
integer :: log_path_site(3,MAX_LOG_PATH_SITES)		! list of tracking start locations
integer :: n_log_path_sites						! number of tracking start locations
integer :: n_log_path(2)						! current numbers of cells being tracked
real, parameter :: R_LOG_PATH = 7				! distance of start location from exit
real, parameter :: LOG_PATH_FACTOR = 0.11		! inflow factor for a single exit (depends on blob size, 18.1 -> 0.11)
integer, parameter :: MAX_PATH_STEPS = 10000	! maximum number of steps to log on a path
type path_type
	integer :: kcell					! cell number
	integer :: id						! cell ID
	integer :: kstep					! number of last location logged (kstep = 1 is the start location)
	logical :: in						! flag for cell in the blob
	integer :: pos(3,MAX_PATH_STEPS)	! sequence of cell locations
end type
type (path_type) :: log_path(MAX_LOG_PATHS,2)

! Parameters for controlling data capture for graphical purposes etc.
logical, parameter :: save_pos_cmgui = .false.          ! To make movies
integer, parameter :: save_interval_hours = 48
real, parameter :: save_length_hours = 0.5      ! 30 minutes
logical, parameter :: generate_exnode = .false.
logical, parameter :: evaluate_stim_dist = .false.
integer, parameter :: ntaglimit_base = 200000	
logical, parameter :: save_DCbinding = .false.
logical, parameter :: log_firstDCcontact = .true.
integer, parameter :: ntres = 60    ! 15 min
logical, parameter :: log_results = .false.
logical, parameter :: log_traffic = .true.
integer, parameter :: TCR_nlevels = 10
real, parameter :: TCR_limit = 2000
integer, parameter :: MAX_AVID_LEVELS = 30

!-------------------------------------------------------------
! Cytokine section
!integer, parameter :: N_CYT = 2
!integer, parameter :: CYT_TAG(1:N_CYT) = (/IL2_TAG, IL7_TAG/)
!integer, parameter :: CYT_NP = IL2_NP + IL7_NP
! diffusion calibration
integer, parameter :: N_CYT = 1
integer, parameter :: CYT_TAG(1:N_CYT) = (/IL2_TAG/)
integer, parameter :: CYT_NP = 0    !IL2_NP
logical, parameter :: CD25_SWITCH = .true.
!-------------------------------------------------------------

! Result type
type result_type
    integer :: dN_EffCogTC(NCTYPES)
    integer :: dN_EffCogTCGen(MMAX_GEN)
    integer :: N_EffCogTC(NCTYPES)
    integer :: N_EffCogTCGen(MMAX_GEN)
    integer :: dN_Dead
    integer :: N_Dead
end type

! Type definitions
!type global_type
	integer :: NHEV
    integer :: NTcells0
    integer :: NTcells
    integer :: NTcellsPer
    integer :: NDC
    integer :: NDCalive
    integer :: NDCcapable
    integer :: Nsites
    integer :: Nexits
    integer :: lastexit
    real :: Radius
    real :: Radius0
!    real :: DCactivity
	real :: ave_residence_time
    real :: InflowTotal
    real :: OutflowTotal
    real :: VEGF
    real :: Vascularity
    real :: dVdt
    real :: c_vegf
!!    real, allocatable :: Inflow(:)
!!    real, allocatable :: Outflow(:)
!end type

type cog_type
    sequence
	real :: avidity			! level of TCR avidity with DC
	real :: stimulation		! TCR stimulation level
!    real :: entrytime       ! time that the cell entered the paracortex (by HEV or cell division)
	real :: dietime			! time that the cell dies
	real :: dividetime		! time that the cell divides
	real :: stagetime		! time that a cell can pass to next stage
	real :: stimrate        ! rate of TCR stimulation
	real :: CD69            ! level of CD69 expression
	real :: S1P1            ! level of S1P1 expression
	real :: CCR7            ! level of CCR7 expression
    real :: IL_state(CYT_NP)    ! receptor model state variable values
    real :: IL_statep(CYT_NP)   ! receptor model state variable time derivative values
    integer :: status       ! holds data in bytes: 1=stage, 2=generation
	integer :: cogID		! index in the list of cognate cells
	integer :: cnt(2)		! count of timesteps 
							! (1) sharing a DC with activated CD4/CD8
							! (2) in close proximity to an activated CD4/CD8
end type

type cell_type
    sequence
    integer :: ID
    integer :: site(3)
    integer :: step
!	real :: DCchemo			! level of susceptibility to DC chemotaxis
	real :: receptor_level(MAX_RECEPTOR)
	real :: receptor_saturation_time(MAX_RECEPTOR)
    integer :: tag
    integer(2) :: ctype
	integer(2) :: lastdir
    integer(2) :: DCbound(2)     ! DCbound(k) = bound DC, allow binding to MAX_DC_BIND DC, MAX_DC_BIND <= 2
    real :: entrytime       ! time that the cell entered the paracortex (by HEV or cell division)
	real :: unbindtime(2)
!    type(cog_type),    pointer :: cptr => NULL()    ! pointer to cognate cell data
    type(cog_type),    pointer :: cptr    ! because NULL is used by winsock (from ifwinty).  NULLIFY() instead.
    ! For DC scanning statistics
	integer :: visits(2),revisits(2),ndclist
	integer(2), allocatable :: dclist(:)	! list of DC visited by a T cell
!	integer(2) :: dclist(2)	! list of DC visited by a T cell
end type

type DC_type
    sequence
    integer :: ID               ! unique ID number
    integer :: site(3)          ! DC location
    integer :: nsites
    integer :: nbound           ! current number of bound T cells
    integer :: ncogbound        ! current number of bound cognate T cells
    integer, allocatable :: cogbound(:)	! list of current bound cognate T cells
    real :: density             ! current antigen density
    real :: dietime             ! time DC will die
    real :: stimulation         ! amount of stimulation provided by a DC
    logical :: capable          ! can DC deliver TCR stimulation?
    logical :: alive            ! is DC alive?
    logical :: cognate			! is the DC bearing antigen cognate for the T cells?
    real :: secretion			! chemokine secretion rate
    real :: conc(MAX_CHEMO)		! chemokine concentration
end type

type occupancy_type
    sequence
    integer(2) :: DC(0:DCDIM-1)     ! DC(0) = number of near DCs, DC(k) = ID of kth near DC
    integer(2) :: cDC(0:cDCDIM-1)   ! cDC(0) = number of DCs within chemo_radius, cDC(k) = ID of kth DC
    integer :: indx(2)
    integer :: exitnum          ! > 0 => index of closest exit, = 0 => no close exit, < 0 => an exit
    real :: chemo_conc				! DC chemokine concentration
    real :: chemo_grad(3)			! DC chemokine gradient vector
    integer :: DC_nbdry
    integer :: isrc					! index into source list
end type

type exit_type
    sequence
    integer :: ID               ! unique ID number
    integer :: site(3)          ! DC location
end type

type dist_type
	integer :: class
	real :: p1, p2, p3
end type

type counter_type
    logical :: logscale
    integer :: period
    integer :: nbins
    real :: binmin, binstep
    integer, allocatable :: ndist(:)
    real :: total
end type

type HEV_type
	integer :: site(3)
end type

type source_type
	integer :: bc
	integer :: site(3)
	real :: level(MAX_CHEMO)
end type

!---------------------------------------------------
! Parameters to read from cell parameters input file
!---------------------------------------------------
real :: TC_AVIDITY_MEAN = 1.0
real :: TC_AVIDITY_SHAPE = 1.2			! shape -> 1 gives normal dist with small variance
real :: TC_CD8_FRACTION					! fraction of all T cells that are CD8
real :: TC_COGNATE_FRACTION(2)			! fraction of T cells that are cognate initially (CD4 and CD8)
real :: TC_CONVERSION_TIME = 48			! time (in hours) over which T cells become cognate
real :: TC_STIM_RATE_CONSTANT = 0.83	! rate const for TCR stimulation (-> molecules/min)
real :: TC_STIM_WEIGHT = 0.1			! contribution of stimulation to act level
real :: TC_STIM_HALFLIFE = 24			! hours
integer :: TC_MAX_GEN = 20              ! maximum number of TC generations

real :: DC_ANTIGEN_MEAN = 10			! mean DC antigen density
real :: DC_ANTIGEN_SHAPE = 1.2			! DC antigen density shape param
real :: DC_LIFETIME_MEAN = 3.5			! days
real :: DC_LIFETIME_SHAPE  = 1.2		! days
real :: DC_ACTIV_TAPER = 12				! time (hours) over which DC activity decays to zero
real :: DC_BIND_DELAY = 2.5				! delay after unbinding before next binding (mins)
!real :: DC_BIND_ALFA = 0.95				! binding prob parameter
real :: DC_MULTIBIND_PROB = 0.0			! reducing factor to bind prob for each current DC binding
real :: DC_DENS_BY_STIM = 0.0002        ! rate of reduction of density by TCR stimulation
real :: DC_DENS_HALFLIFE                ! half-life of DC activity (hours)

integer :: optionA
integer :: optionB
integer :: optionC
real :: IL2_PRODUCTION_TIME			    ! duration of IL2/CD25 production (hrs)
real :: IL2_THRESHOLD			        ! TCR stimulation needed to initiate IL-2/CD25 production
real :: ACTIVATION_THRESHOLD    		! combined stimulation needed for activation
real :: FIRST_DIVISION_THRESHOLD(2)		! activation level needed for first division
real :: DIVISION_THRESHOLD(2)			! activation level needed for subsequent division
real :: EXIT_THRESHOLD(2)               ! Activation/CD69 level S below which exit is permitted
real :: STIMULATION_LIMIT				! maximum activation level
real :: CD25_DIVISION_THRESHOLD         ! CD25 store level needed for division of activated cell
real :: CD25_SURVIVAL_THRESHOLD         ! CD25 store level needed for survival of activated cell

type(dist_type) :: divide_dist1
type(dist_type) :: divide_dist2
real :: CD8_DIVTIME_FACTOR = 1.5		! CD8 divide time as multiple of CD4 divide time

real :: TC_FRACTION
real :: TC_RADIUS
integer :: DC_RADIUS
real :: BLOB_RADIUS
real :: FLUID_FRACTION

real(DP) :: GAMMA                           ! controls crowding
real(DP) :: BETA                            ! speed: 0 < beta < 1
real(DP) :: RHO                             ! persistence: 0 < rho < 1

integer :: TC_TO_DC                     ! number of T cells for every DC
real :: DC_FACTOR                       ! multiplying factor for DC number (initial and influx)
!integer :: NTC_DC                      ! number of T cells that can bind to a DC
integer :: MAX_TC_BIND                  ! number of T cells that can bind to a DC
integer :: MAX_DC_BIND                  ! number of DC that a cell can bind to simultaneously
integer :: MAX_COG_BIND                 ! number of cognate T cells that can bind to a DC simultaneously
real :: DCrate_100k                     ! DC influx rate corresponding to NTcells0 = 100k
real :: T_DC1                           ! Duration of constant DC influx (hours)
real :: T_DC2                           ! Time of cessation of DC influx (hours)
logical :: DC_INJECTION					! DCs were injected into experimental animals?
real :: T_DC_INJECTION					! Time of DC injection
!real :: DC_FRACTION						! Fraction of DCs that are bearing antigen
logical :: use_traffic = .false.
logical :: use_exit_chemotaxis
logical :: use_DC_chemotaxis
logical :: computed_outflow
logical :: use_single_DC = .false.		! To evaluate chemokine field from a single DC

real :: RESIDENCE_TIME(2)                  ! T cell residence time in hours -> inflow rate
! Vascularity parameters
real :: Inflammation_days1 = 4          ! Days of plateau level - parameters for VEGF_MODEL = 1
real :: Inflammation_days2 = 5          ! End of inflammation
real :: Inflammation_level = 1.0		! This is the level of inflammation (scaled later by NTcells0)
integer :: exit_region                   ! determines blob region for cell exits
real :: efactor                         ! If constant_efactor = true, this is the factor for the p correction
integer :: VEGF_MODEL                   ! 1 = VEGF signal from inflammation, 2 = VEGF signal from DCactivity
real :: chemo_radius					! radius of chemotactic influence (um)
real :: chemo_K_exit                    ! level of chemotactic influence towards exits
real :: chemo_K_DC                      ! level of chemotactic influence towards DCs

logical :: fix_avidity                  ! true if avidity takes discrete values, false if a distribution is used
logical :: avidity_logscale             ! true if actual avidity = 10^(avidity_min + i*avidity_step)
integer :: avidity_nlevels              ! If fix_avidity, number of discrete avidity values
real :: avidity_min                     ! minimum value
real :: avidity_step                    ! step between equi-spaced values
real :: days                            ! number of days to simulate
integer :: seed(2)                      ! seed vector for the RNGs
integer :: NT_GUI_OUT					! interval between GUI outputs (timesteps)
integer :: SPECIES						! animal species source of T cells
logical :: IN_VITRO = .false.			! select in vivo or in vitro simulation
real :: IV_WELL_DIAMETER				! diameter of in vitro well (mm)
integer :: IV_NTCELLS					! initial T cell population in vitro
real :: IV_COGNATE_FRACTION				! fraction of in vitro cells that are cognate for DC antigen
logical :: IV_SHOW_NONCOGNATE = .false.	! display non-cognate T cells
character*(128) :: fixedfile

!---------------------------------------------------
! end of parameters to read from input file
!---------------------------------------------------

!---------------------------------------------------
! More input parameters
!---------------------------------------------------
integer :: NX

! T cell parameters
logical :: TCR_splitting = .false.      ! enable sharing of integrated TCR signal between progeny cells
real :: transient_stagetime				! stagetime for TRANSIENT (Stage 1)
real :: clusters_stagetime				! stagetime for CLUSTERS (Stage 2)
real :: transient_bindtime				! bindtime for TRANSIENT (Stage 1)
real :: clusters_bindtime				! bindtime for CLUSTERS (Stage 2)
real :: swarms_bindtime					! bindtime for SWARMS (Stage 3)
real :: TC_life_median1					! median lifetime of naive T cells
real :: TC_life_median2					! median lifetime of activated T cells
real :: TC_life_shape					! shape parameter for lifetime of T cells
integer :: NTC_LN = 3.0e07				! number of T cells in a LN
integer :: NTC_BODY = 1.6e09			! number of circulating T cells in the whole body
integer :: NLN_RESPONSE					! number of LNs in the response
real :: K1_S1P1 = 0.005					! S1P1/CD69 system parameters
real :: K2_S1P1 = 0.05
real :: K1_CD69 = 0.04
real :: K2_CD69 = 0.01
! Parameters added to control HEV and DC placement (for DCvisits simulations)
real :: R_HEV_min
real :: R_HEV_max
real :: R_DC_min
real :: R_DC_max

! DC parameters
integer :: NDCsites						! Number of lattice sites occupied by the DC core (soma)
logical :: use_DCflux = .true.
real :: DC_pMHC_THRESHOLD = 10          ! DC pMHC limit for TCR stimulation capability
real :: DC_STIM_THRESHOLD = 300		    ! TCR stimulation Hill function parameter
integer :: N_STIM = 1					! TCR stimulation Hill function exponent
integer :: CONTACT_RULE = CT_HENRICKSON ! rule for determining the duration of T cell - DC contact
real :: ABIND1 = 0.4, ABIND2 = 0.8      ! binding to a DC

! Egress parameters
real :: exit_fraction = 1.0/1000.       ! number of exits as a fraction of T cell population
real :: Ksurfaceportal = 40				! calibration factor for number of surface portals
logical :: suppress_egress = .false.	! transient suppression of egress (see EGRESS_SUPPRESSION_TIME1, 2)

!---------------------------------------------------
! end of more parameters to be read from input file
!---------------------------------------------------

! Cytokine data
real, allocatable :: cyt(:,:,:,:)
real, allocatable :: cyt_constit(:), cyt_mols(:), dcyt_mols(:)
real :: cyt0(MAX_CYT) = (/3.0,1.0,1.0,1.0,1.0,1.0/)
real :: K_diff(MAX_CYT), delta_diff(MAX_CYT)     ! = D.dt/(dx*dx)
integer :: Ncytokines, cytokines(MAX_CYT), cyt_seq(MAX_CYT), NP_offset(MAX_CYT+1)
real :: cyt_init(MAX_CYT), cyt_mean(MAX_CYT)

! Geometry data
integer :: NY, NZ
integer, allocatable :: xdomain(:),xoffset(:),zdomain(:),zoffset(:)
integer :: blobrange(3,2)
integer, allocatable :: nz_sites(:), nz_totsites(:), nz_cells(:), nz_excess(:)
type(HEV_type), allocatable :: HEVlist(:)
real :: proximity_limit(4,4)
real :: R_limit(4,2)
real :: DELTA_X, PI
real :: TagRadius
real :: x0,y0,z0   ! centre in global coordinates (units = grids)
real :: Centre(3)
real :: Vc, Ve

! Motility data
integer :: nreldir, njumpdirs
integer :: jumpvec(3,27)    ! 14 is no-jump case (0,0,0)
!integer :: jumpvec2D(3,8)
integer :: reldir(6,MAXRELDIR)
integer :: DCoffset(3,max(NDCcore,NDCcore_IN_VITRO))
real(DP) :: dirprob(0:MAXRELDIR)

integer :: nreldir2D, njumpdirs2D
integer :: reldir2D(8,8)
real(DP) :: dirprob2D(0:8)
logical :: diagonal_jumps
real :: ep_factor      ! (25k) 2.4 when K1_S1P1 = 0.01, 2.8 when K1_S1P1 = 0.001  ! based on no-DC case!
                       ! (50k) 2.3 when K1_S1P1 = 0.01, chemo_K_exit = 0.3

! Chemotaxis data
integer :: chemo_N
real :: chemo_exp, DC_CHEMO_DELAY	! gets set = DC_BIND_DELAY
real :: t_taglimit

real, allocatable :: chemo_r(:,:,:)
real, allocatable :: chemo_p(:,:,:,:)

! Cell data
type(occupancy_type), allocatable :: occupancy(:,:,:)
type(cell_type), allocatable, target :: cellist(:)
type(DC_type), allocatable :: DClist(:)
!type(DC_type), allocatable :: FDClist(:)
!integer :: NFDC
type(source_type), allocatable :: sourcelist(:)
integer :: nsources
integer, allocatable :: cognate_list(:)
integer, allocatable :: DCdeadlist(:)
integer, allocatable :: DCinjected(:)
integer, allocatable :: gaplist(:)
integer, allocatable :: zrange2D(:,:,:)	! for 2D case
integer :: DCstack(-4:4,-4:4,3)			! for 2D case
integer :: lastID, MAX_COG, lastcogID, nlist, n2Dsites, ngaps, ntaglimit, ntagged=0, ntagged_left, ID_offset, ncogseed(2)
integer :: lastNTcells, k_nonrandom(2)
integer :: max_nlist, max_ngaps
integer :: nadd_sites, MAX_DC, ndeadDC, NDCtotal, ndivisions
real :: excess_factor, lastbalancetime
real :: scale_factor	! scaling from model to one (or more) whole LNs
real :: Fcognate		! fraction of T cells in circulation that are cognate
logical :: use_DC, use_cognate

! Result data
type(result_type) :: localres, totalres
integer :: nvisits(2,2), nrevisits(2,2)
integer, allocatable :: DCvisits(:,:,:,:)
integer, allocatable :: DCtotvisits(:,:,:,:)
!integer :: tot_nvisits, tot_nrevisits, tot_DCvisits(0:MAX_DC+1)
integer :: noutflow_tag, ninflow_tag
real :: restime_tot
real, allocatable :: Tres_dist(:)
type(counter_type) :: avid_count,avid_count_total
integer :: dcbind(0:50)
logical :: firstSummary

! Travel time computation
integer :: ntravel
integer :: N_TRAVEL_COG, N_TRAVEL_DC, N_TRAVEL_DIST
integer :: k_travel_cog, k_travel_dc
integer, allocatable :: travel_dc(:)
real, allocatable :: travel_cog(:), travel_dist(:,:,:)

! Time to first DC contact computation
integer :: firstDC_n(2,2)
real :: firstDC_tot(2,2)
real :: firstDC_dist(2,2,10000)

! Miscellaneous data
type(exit_type), allocatable :: exitlist(:)
!type(global_type) :: globalvar
integer :: max_exits        ! size of the exitlist(:) array
integer :: nbindmax, nbind1, nbind2   ! these parameters control the prob of a T cell
real :: min_transit_time = 60		! minimum time a T cell spends in the DCU (min)
real :: CD69_threshold				! level of CD69 below which egress can occur
real :: last_portal_update_time		! time that the number of exit portals was last updated
logical :: initialized, steadystate
integer :: navid = 0
integer ::  Nsteps, nsteps_per_min, istep
integer :: Mnodes
integer :: IDtest
integer :: total_in = 0, total_out = 0
integer :: nIL2thresh = 0           ! to store times to IL2 threshold
real :: tIL2thresh = 0
integer :: ndivided(MMAX_GEN) = 0   ! to store times between divisions
real :: tdivided(MMAX_GEN) = 0
real :: DCdecayrate                 ! base rate of depletion of DC activity (/min) 
real :: TCRdecayrate                ! rate of decay of integrated TCR stimulation (/min)
real :: max_TCR = 0
real :: avidity_level(MAX_AVID_LEVELS)  ! discrete avidity levels for use with fix_avidity
logical :: vary_vascularity = .true.     ! to allow inflammation to change vascularity (false if VEGF_MODEL = 0)

integer :: check_egress(1000)		! for checking traffic
integer :: check_inflow

! Vascularity parameters
real :: VEGF_alpha = 4.0e-7         ! rate constant for dependence on inflammation (/min) (alpha_G in hev.m) (was 5.0e-7)
real :: VEGF_beta = 5.0e-8			! rate constant for basal VEGF production (beta_G in hev.m) (was 4.0e-8)
real :: VEGF_decayrate = 0.002      ! VEGF decay rate (/min)	(was 0.002)
!real :: vasc_maxrate = 0.0006       ! max rate constant for vascularity growth (/min)
real :: vasc_maxrate = 0.001       ! max rate constant for vascularity growth (/min)  (was 0.003)
!real :: vasc_beta = 1.5				! Hill function parameter
real :: vasc_beta = 2.0				! Hill function parameter
integer :: vasc_n = 2               ! Hill function exponent
! NOTE: all the parameter changes I've been trying do not get the peak N/N0 above 4, with 
! inflam = 1, days 3,4.

real :: vasc_decayrate				! vascularity decay rate (/min) (deduced)
real :: VEGF_baserate				! base rate of production of VEGF (> 0 for VEGF-vascularity model VEGF_MODEL = 1)
real :: c_vegf_0					! steady-state VEGF concentration (VEGF_MODEL = 1)
real :: Kflow1                          ! Inflow dependence on expansion
real :: Kflow2                          ! Outflow dependence on expansion
real :: Kflow3                          ! Maximum amplification of inflow from DC activity
real :: Bflow                           ! Hill function parameter for amplification
integer :: Nflow                        ! Hill function power coefficient

real :: TC_AVIDITY_MEDIAN				! median T cell avidity
real :: DC_ANTIGEN_MEDIAN				! median DC antigen density
real :: DC_LIFETIME_MEDIAN				! days

character*(128) :: inputfile
character*(128) :: outputfile
!character*(128) :: resultfile
character*(128) :: DCinjectionfile
character*(2048) :: logmsg
TYPE(winsockport) :: awp_0, awp_1, awp_2, awp_3
logical :: use_TCP = .true.         ! turned off in para_main()
logical :: use_CPORT1 = .false.
logical :: stopped, clear_to_send, simulation_start, par_zig_init
logical :: dbug = .false.
logical :: USE_PORTAL_EGRESS			! use fixed exit portals rather than random exit points
logical :: BLOB_PORTALS					! egress to the sinus at portals throughout the blob
logical :: SURFACE_PORTALS				! egress to the sinus at portals on the blob surface
logical :: FIXED_NEXITS = .false.		! the number of exit portals is held fixed
real :: XFOLLICLE = 0.6					! normalized x boundary of follicular interface "cap"
real :: EGRESS_SUPPRESSION_TIME1 = 12	! hours
real :: EGRESS_SUPPRESSION_TIME2 = 24	! hours
real :: EGRESS_SUPPRESSION_RAMP = 6		! hours
real :: INLET_R_FRACTION = 0.7			! fraction of blob radius within which ingress occurs
logical :: RELAX_INLET_EXIT_PROXIMITY = .false.	! override INLET_R_FRACTION, allow inlet closer to exits
real :: INLET_LAYER_THICKNESS = 5		! if RELAX_INLET_EXIT_PROXIMITY, inlets are within this distance of the blob boundary
real :: INLET_EXIT_LIMIT = 5			! if RELAX_INLET_EXIT_PROXIMITY, this determines how close an inlet point can be to an exit portal.
real :: CHEMO_K_RISETIME = 120			! if RELAX_INLET_EXIT_PROXIMITY, this the the time for chemotaxis to reach full strength (mins) 

! PERIPHERY parameters
logical, parameter :: SIMULATE_PERIPHERY = .true.
integer, parameter :: PERI_GENERATION = 2
real, parameter :: PERI_PROBFACTOR = 10

!DEC$ ATTRIBUTES DLLEXPORT :: ntravel, N_TRAVEL_COG, N_TRAVEL_DC, N_TRAVEL_DIST, k_travel_cog, k_travel_dc
!DEC$ ATTRIBUTES DLLEXPORT :: travel_dc, travel_cog, travel_dist
!DEC$ ATTRIBUTES DLLEXPORT :: nsteps	!istep
contains

!---------------------------------------------------------------------
! Uniform randomly generates an integer I: n1 <= I <= n2
!---------------------------------------------------------------------
integer function random_int(n1,n2,kpar)
!use IFPORT
integer :: n1,n2,kpar
integer :: k,R

!k = irand()     ! intrinsic
if (n1 == n2) then
    random_int = n1
elseif (n1 > n2) then
    write(logmsg,*) 'ERROR: random_int: n1 > n2: ',n1,n2
    call logger(logmsg)
    stop
endif
R = par_shr3(kpar)
k = abs(R)
random_int = n1 + mod(k,(n2-n1+1))

end function

!--------------------------------------------------------------------------------
! Returns a permutation of the elements of a()
!--------------------------------------------------------------------------------
subroutine permute(a,n,kpar)
integer :: a(*),n,kpar
integer :: i,k,tmp

do i = 1,n
    k = random_int(1,n,kpar)
	tmp = a(i)
	a(i) = a(k)
	a(k) = tmp
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine waste_time(n,dummy)
integer :: k, n
real :: dummy
real(DP) :: rsum,R
integer :: kpar=0

rsum = 0
do k = 1,n
!    call random_number(R)
    R = par_uni(kpar)
    rsum = rsum + R
enddo
dummy = rsum
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real function norm(r)
real :: r(3)

norm = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real function norm2(r)
real :: r(3)

norm2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine normalize(r)
real :: r(3)

r = r/norm(r)
end subroutine

!-----------------------------------------------------------------------------------------
! Distance from the blob centre (units = grids)
!-----------------------------------------------------------------------------------------
real function cdistance(site)
integer :: site(3)
real :: r(3)

r = site - Centre
cdistance = norm(r)
end function

!-----------------------------------------------------------------------------------------
! Distance from the iexit exit (units = grids)
!-----------------------------------------------------------------------------------------
real function ExitDistance(site,iexit)
integer :: site(3),iexit
real :: r(3)

r = site - exitlist(iexit)%site
ExitDistance = norm(r)
end function

!--------------------------------------------------------------------------------
! A site is taggable if it is less than a specified distance TagRadius from
! the centre.
!--------------------------------------------------------------------------------
logical function taggable(site)
integer :: site(3)
real :: d, r(3)

r = site - Centre
d = norm(r)
if (d <= TagRadius) then
    taggable = .true.
else
    taggable = .false.
endif
end function

!-----------------------------------------------------------------------------------------
! Squeeze gaps out of cellist array, adjusting occupancy array.
!-----------------------------------------------------------------------------------------
subroutine squeezer(force)
logical :: force
integer :: last, k, site(3), indx(2), i, n, region

!write(*,*) 'squeezer'
if (ngaps == 0) return
if (.not.force .and. (ngaps < max_ngaps/2)) return
if (dbug) write(nflog,*) 'squeezer: ',ngaps,max_ngaps,nlist

n = 0
do k = 1,nlist
    if (cellist(k)%ID == 0) then    ! a gap
        n = n+1
!        write(*,*) 'gap at : ',k
    endif
enddo

last = nlist
k = 0
n = 0
do
    k = k+1
    if (cellist(k)%ID == 0) then    ! a gap
        if (k == last) exit
        do
            if (last == 0) then
                write(nflog,*) 'last = 0: k: ',k
                stop
            endif
            if (cellist(last)%ID == 0) then
                last = last-1
                n = n+1
                if (n == ngaps) exit
            else
                exit
            endif
        enddo
        if (n == ngaps) exit
        call copycell2cell(cellist(last),cellist(k),k)
!        cellist(k) = cellist(last)
		if (associated(cellist(last)%cptr)) then
			call get_region(cellist(last)%cptr,region)
		else
			region = LYMPHNODE
		endif
		if (region == LYMPHNODE) then
	        site = cellist(last)%site
	        indx = occupancy(site(1),site(2),site(3))%indx
	        do i = 1,2
	            if (indx(i) == last) indx(i) = k
	        enddo
	        occupancy(site(1),site(2),site(3))%indx = indx
	    endif
        last = last-1
        n = n+1
    endif
    if (n == ngaps) exit
enddo
nlist = nlist - ngaps
ngaps = 0
if (dbug) write(nflog,*) 'squeezed: ',n,nlist

end subroutine

!-----------------------------------------------------------------------------------------
! Copy the contents of cellist(kfrom) to the entry cellist(kto)
! Need to fix up the corresponding cognate_list() entry as well.
! Note that:
!   cell = cellist(kcell)
!   kcog = cell%cptr%cogID
!   k = cognate_list(kcog)
! =>
!   k = kcell
!-----------------------------------------------------------------------------------------
subroutine copycell2cell(cell_from,cell_to,kcell)
integer :: kcell
type(cell_type) :: cell_from, cell_to
integer :: ctype, stype, kcog

ctype = cell_from%ctype
!stype = struct_type(ctype)
if (associated(cell_from%cptr)) then
	stype = COG_TYPE_TAG
else
	stype = NONCOG_TYPE_TAG
endif
if (stype == NONCOG_TYPE_TAG .and. associated(cell_to%cptr)) then
    deallocate(cell_to%cptr)
endif
if (stype == COG_TYPE_TAG) then
    if (.not.associated(cell_to%cptr)) then
        allocate(cell_to%cptr)
    endif
    cell_to%cptr = cell_from%cptr
    kcog = cell_to%cptr%cogID
    cognate_list(kcog) = kcell
elseif (stype /= NONCOG_TYPE_TAG) then
    write(*,*) 'ERROR: copycell2cell: istep, ID, ctype, stype: ',istep,cell_from%ID,ctype,stype
    stop
endif
cell_to%ID = cell_from%ID
cell_to%site = cell_from%site
!cell_to%DCchemo = cell_from%DCchemo
cell_to%receptor_level = cell_from%receptor_level
cell_to%receptor_saturation_time = cell_from%receptor_saturation_time
cell_to%tag = cell_from%tag
cell_to%ctype = cell_from%ctype
cell_to%lastdir = cell_from%lastdir
cell_to%DCbound = cell_from%DCbound
cell_to%entrytime = cell_from%entrytime
cell_to%unbindtime = cell_from%unbindtime
cell_to%visits = cell_from%visits
cell_to%revisits = cell_from%revisits
cell_to%ndclist = cell_from%ndclist
allocate(cell_from%dclist(max(1,NDC)))
cell_to%dclist = cell_from%dclist

if (cell_from%ctype == 0) then
    write(*,*) 'ERROR: copycell2cell: ctype = 0'
    stop
endif

end subroutine

!--------------------------------------------------------------------------------
! Returns:
! 0 if no slots are occupied
! 1 if slot 1 is occupied
! 2 if slot 2 is occupied
! 3 if both slots are occupied
!--------------------------------------------------------------------------------
integer function getslots(site)
integer :: site(3)
integer :: k

getslots = 0
do k = 1,2
    if (occupancy(site(1),site(2),site(3))%indx(k) > 0) then
        getslots = getslots + k
    endif
enddo
end function


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine initial_binding
integer :: kcell, site(3), dc(0:3), k, idc
real(DP) :: R
real :: bindtime
type(cell_type) :: tcell
integer :: kpar=0
logical :: cognate

do kcell = 1,nlist
    tcell = cellist(kcell)
    cognate = (associated(tcell%cptr))
    site = tcell%site
    dc = occupancy(site(1),site(2),site(3))%DC
    if (dc(0) > 0) then
        do k = 1,dc(0)
            idc = dc(k)
            if (.not.DClist(idc)%capable) cycle
            if (cognate .and. DClist(idc)%ncogbound == MAX_COG_BIND) cycle
            if (bindDC(idc,kpar)) then
!                call random_number(R)
                R = par_uni(kpar)
                bindtime = 1 + R*5.     ! mean = 3.5 min
!                call random_number(R)
                R = par_uni(kpar)
                cellist(kcell)%DCbound(1) = idc
                cellist(kcell)%unbindtime(1) = bindtime*R
                DClist(idc)%nbound = DClist(idc)%nbound + 1
                if (cognate) then
                    DClist(idc)%ncogbound = DClist(idc)%ncogbound + 1
                    call AddCogbound(idc,kcell)
                endif
                exit
            endif
        enddo
    endif
enddo
!write(*,*) 'DC occupancy:'
!write(*,'(8i6)') DClist(1:NDC)%nbound
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_DCproximity
integer :: x,y,z,k,cnt(0:DCDIM-1)
integer :: dc(DCDIM-1)

cnt = 0
do x = 1,NX
    do y = 1,NY
        do z = 1,NZ
            k = occupancy(x,y,z)%DC(0)
            if (k >= 0) then
                cnt(k) = cnt(k) + 1
                dc = occupancy(x,y,z)%DC(1:DCDIM-1)
                if (k >= 2) then
                    call randomizeDC(dc,k)
                    occupancy(x,y,z)%DC(1:DCDIM-1) = dc
                endif
            endif
        enddo
    enddo
enddo
!write(*,*) 'DC proximity counts: ',cnt

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine randomizeDC(dc,k)
integer :: k, dc(:)
integer :: i,p(DCDIM-1),tmp(DCDIM-1)
integer :: kpar=0

do i = 1,k
    p(i) = i
enddo
call permute(p,k,kpar)
tmp(1:k) = dc(1:k)
do i = 1,k
    dc(i) = tmp(p(i))
enddo
!write(*,*) tmp(1:k),dc(1:k)
end subroutine

!-----------------------------------------------------------------------------------------
! Is site near a DC?
! The criterion for a DC site might be different from an exit site.
! prox = DC_DCprox*DC_RADIUS for DC - DC
!-----------------------------------------------------------------------------------------
logical function tooNearDC(site,prox)
integer :: site(3)
real :: prox
integer :: idc
real :: r(3), d

if (NDCalive == 0) then
    tooNearDC = .false.
    return
endif
do idc = 1,NDC
    if (.not.DClist(idc)%alive) cycle
    r = site - DClist(idc)%site
    d = norm(r)     ! units sites
    if (d < prox) then
        tooNearDC = .true.
        return
    endif
enddo
tooNearDC = .false.
end function

!-----------------------------------------------------------------------------------------
! prox is the minimum distance from the boundary (sites)
!-----------------------------------------------------------------------------------------
logical function tooNearBdry(site,prox)
integer :: site(3)
real :: prox
real :: d, dmin

!dmin = proxfac*DC_RADIUS
if (use_blob) then
    d = cdistance(site)
    if (Radius - d < prox) then
        tooNearBdry = .true.
        return
    endif
!    write(*,*) 'tooNearBdry: ',Radius,d,dmin
else
    if (site(1) < dmin .or. (NX - site(1)) < dmin) then
        tooNearBdry = .true.
        return
    endif
    if (site(2) < dmin .or. (NY - site(2)) < dmin) then
        tooNearBdry = .true.
        return
    endif
    if (site(3) < dmin .or. (NZ - site(3)) < dmin) then
        tooNearBdry = .true.
        return
    endif
endif
tooNearBdry = .false.
end function

!-----------------------------------------------------------------------------------------
! Is site near an exit?  prox is the minimum separation (sites)
!-----------------------------------------------------------------------------------------
logical function tooNearExit(site,prox)
integer :: site(3)
real :: prox
integer :: iexit
real :: r(3), d

if (Nexits == 0) then
    tooNearExit = .false.
    return
endif
if (.not.allocated(exitlist)) then
	call logger('exitlist not allocated')
	stop
endif
do iexit = 1,lastexit
    if (exitlist(iexit)%ID == 0) cycle  ! exit no longer exists
    r = site - exitlist(iexit)%site
    d = norm(r)
    if (d < prox) then
        tooNearExit = .true.
        return
    endif
enddo
tooNearExit = .false.
end function

!-----------------------------------------------------------------------------------------
! Is site near an HEV?  prox is the minimum separation (sites)
!-----------------------------------------------------------------------------------------
logical function tooNearHEV(site,prox)
integer :: site(3)
real :: prox
integer :: ihev
real :: r(3), d

if (NHEV == 0) then
    tooNearHEV = .false.
    return
endif
do ihev = 1,NHEV
    r = site - HEVlist(ihev)%site
    d = norm(r)
    if (d < prox) then
        tooNearHEV = .true.
        return
    endif
enddo
tooNearHEV = .false.
end function

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
integer function neighbourhoodCount(site)
integer :: site(3)
integer :: site2(3), k, count

count = 0
do k = 1,27
	site2 = site + jumpvec(:,k)
	if (site2(1) < 1 .or. site2(1) > NX) cycle
	if (site2(2) < 1 .or. site2(2) > NY) cycle
	if (site2(3) < 1 .or. site2(3) > NZ) cycle
	if (occupancy(site2(1),site2(2),site2(3))%indx(1) >= 0) then
		count = count + 1
	endif
enddo
neighbourhoodCount = count
end function

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
logical function bindDC(idc,kpar)
integer :: idc,kpar
integer :: n
real(DP) :: p, R

n = DClist(idc)%nbound
if (n >= nbind2) then
    bindDC = .false.
elseif (n <= nbind1) then
    bindDC = .true.
else
    p = real(nbind2 - n)/(nbind2 - nbind1)
!    call random_number(R)
    R = par_uni(kpar)
    if (R < p) then
        bindDC = .true.
    else
        bindDC = .false.
    endif
endif
end function

!-----------------------------------------------------------------------------------------
! Add a cognate T cell index to the list of bound cognate cells for DC idc
!-----------------------------------------------------------------------------------------
subroutine AddCogbound(idc,kcell)
integer :: idc, kcell
integer :: i

do i = 1,MAX_COG_BIND
	if (DClist(idc)%cogbound(i) == 0) then
		DClist(idc)%cogbound(i) = kcell
		return
	endif
enddo
write(logmsg,*) 'Error: AddCogbound: ',idc,kcell
call logger(logmsg)
stop
end subroutine

!-----------------------------------------------------------------------------------------
! Remove a cognate T cell index from the list of bound cognate cells for DC idc
!-----------------------------------------------------------------------------------------
subroutine RemoveCogbound(idc,kcell)
integer :: idc, kcell
integer :: i

do i = 1,MAX_COG_BIND
	if (DClist(idc)%cogbound(i) == kcell) then
		DClist(idc)%cogbound(i) = 0
		return
	endif
enddo
write(logmsg,*) 'Error: RemoveCogbound: ',idc,kcell
call logger(logmsg)
stop
end subroutine


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!subroutine showsite(node,site)
!integer :: node, site(3)
!
!write(*,*) 'showsite: ',site
!write(*,*) 'indx: ',occupancy(site(1),site(2),site(3))%indx
!end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
integer function cell_count()
integer :: kcell, ntot
ntot = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle             ! skip gaps in the list
    ntot = ntot + 1
enddo
cell_count = ntot
end function

!-----------------------------------------------------------------------------------------
! The type (cognate or non-cognate) of a T cell can be random or strictly determined.
!-----------------------------------------------------------------------------------------
!integer function select_cell_type(kpar)
subroutine select_cell_type(ctype,cognate,kpar)
integer :: ctype, kpar
logical :: cognate
!integer, save :: k = 0
integer :: nratio

ctype = select_CD4_CD8()
if (random_cognate) then
!    select_cell_type = random_cell_type(kpar)
    if (par_uni(kpar) < TC_COGNATE_FRACTION(ctype)) then
		cognate = .true.
	else
		cognate = .false.
	endif
else
    ! Needs to be fixed to account for two cognate populations
!	if (mod(istep,6*4*60) == 1) then	! update Fcognate every 6 hours - otherwise mod(k_nonrandom,nratio) messed up
!		Fcognate = TC_COGNATE_FRACTION - (scale_factor*Ncogseed)/NTC_BODY
!	endif
!    nratio = 1./Fcognate
    nratio = 1./TC_COGNATE_FRACTION(ctype)
    k_nonrandom(ctype) = k_nonrandom(ctype) + 1
    if (mod(k_nonrandom(ctype),nratio) == 0) then
!        select_cell_type = COG_TYPE_TAG
		cognate = .true.
    else
!        select_cell_type = NONCOG_TYPE_TAG
		cognate = .false.
    endif
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
integer function select_CD4_CD8
integer :: kpar = 0

if (par_uni(kpar) < CTYPE_FRACTION(CD8)) then
	select_CD4_CD8 = CD8
else
	select_CD4_CD8 = CD4
endif
end function

!-----------------------------------------------------------------------------------------
! The global variable values are computed from the global data:
!    NTcells
!    NDC
!    Radius
!    LNactivation
!    Inflow
!    Outflow
!    FlowFraction(:)
!-----------------------------------------------------------------------------------------
subroutine set_globalvar
real :: inflow0

if (TAGGED_LOG_PATHS) then
	InflowTotal = LOG_PATH_FACTOR*NTcells0*DELTA_T/(ave_residence_time*60)
	OutflowTotal = InflowTotal
	return
endif
if (IN_VITRO) then
	InflowTotal = 0
	OutflowTotal = 0
	return
endif
if (use_traffic) then
    inflow0 = NTcells0*DELTA_T/(ave_residence_time*60)
else
    inflow0 = 0
endif

if (.not.steadystate) then     ! surrogate for modeling an immune response
    call generate_traffic(inflow0)
else
    InflowTotal = inflow0
    OutflowTotal = inflow0
endif
if (istep == 1 .and. .not.use_TCP) then
	write(*,'(a,i8,4f8.2)') 'NTcells,Inflow,Outflow: ',NTcells0,InflowTotal,OutflowTotal
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Total T cell inflow and outflow are generated from the vascularity and baseline
! inflow, inflow0.
!-----------------------------------------------------------------------------------------
subroutine generate_traffic(inflow0)
real :: inflow0
real :: act, expansion, actfactor, tnow
real :: inflow, outflow
!real, parameter :: T1 = 8*60, T2 = 16*60, T3 = 24*60

!if (.not.use_vascularity) then
!    act = get_DCactivity()
!    DCactivity = act
!    InflowTotal = inflow0
!    OutflowTotal = inflow0
!    steadystate = .true.
!    return
!endif
if (traffic_mode == TRAFFIC_MODE_1) then    ! naive
    write(*,*) 'generate_traffic: do not use TRAFFIC_MODE_1'
    stop
    act = get_DCactivity()
    act = act*100/NTcells ! to make act into a concentration, and Bflow approx 1.0
    expansion = real(NTcells)/NTcells0 - 1
    actfactor = 1 + (Kflow3-1)*act**Nflow/(act**Nflow + Bflow**Nflow)
    write(*,'(a,4f6.2)') 'act,actfactor,expansion,radius: ',act,actfactor,expansion,Radius
    inflow = inflow0*actfactor*(1 + Kflow1*expansion)
    inflow = max(inflow,inflow0)
    outflow = inflow0*(1 + Kflow2*expansion)
    outflow = max(outflow,inflow0)
else        !traffic_mode == TRAFFIC_MODE_2 or TRAFFIC_MODE_3
	! Note: if inflammation signal = 0 the vascularity (and inflow) should be constant
    tnow = istep*DELTA_T
    inflow = inflow0*Vascularity   ! level of vascularity (1 = steady-state)
!    if (tnow <= T1) then
!        outflow = ((T1-tnow)/T1)*inflow0
!    elseif (tnow <= T2) then
!        outflow = 0
!    elseif (tnow <= T3) then
!        outflow = ((tnow-T2)/(T3-T2))*NTcells*DELTA_T/(RESIDENCE_TIME*60)
!    else
        outflow = NTcells*DELTA_T/(ave_RESIDENCE_TIME*60)
!    endif
endif

if (L_selectin) then
	if (mod(istep,1000) == 0) then
		call logger('Using L_selectin!')
	endif
    OutflowTotal = NTcells*DELTA_T/(ave_RESIDENCE_TIME*60)
    InflowTotal = 0
    return
endif

!DCactivity = act
InflowTotal = inflow
! This is a kludge to induce a return to steady-state maintenance when NTcells drops
! back to "close enough" to the steady-state value.
! I think this is no longer needed.
! Try removing it
!if (use_exit_chemotaxis .and. NTcells < 0.99*NTcells0) then
!    OutflowTotal = InflowTotal      ! => steady-state with chemotaxis
!    steadystate = .true.
!else
    OutflowTotal = outflow
!endif
!if (mod(istep,240) == 0) then
!	write(logmsg,*) 'generate_traffic: inflow: ',inflow0,Vascularity,InflowTotal
!	call logger(logmsg)
!endif
end subroutine

!-----------------------------------------------------------------------------------------
! Pre-programmed inflammation signal is flat at Inflammation_100k (for NTcells0 = 100k)
! for Inflammation_days1 then ends at Inflammation_days2
!-----------------------------------------------------------------------------------------
real function get_inflammation()
real :: tnow, plateau

tnow = istep*DELTA_T    ! mins
tnow = tnow/(24*60)     ! days
plateau = Inflammation_level*NTcells0
if (tnow < Inflammation_days1) then
    get_inflammation = plateau
elseif (tnow > Inflammation_days2) then
    get_inflammation = 0
else
    get_inflammation = plateau*(Inflammation_days2 - tnow)/(Inflammation_days2 - Inflammation_days1)
endif
end function

!-----------------------------------------------------------------------------------------
! Computes the total DC activity summed over all DC that are capable.
! For this to work there must be a model for:
! (a) influx of cognate DC during the response
! (b) the decline of DC stimulating ability (%density) either as a function of time,
!     or in proportion to TCR stimulation delivered.  This is provided by update_DCstate()
!-----------------------------------------------------------------------------------------
real function get_DCactivity()
integer :: idc, ncap
real :: a, tnow, td

if (test_vascular) then
    ! A simple variation in DC activity to test vascular development
    tnow = istep*DELTA_T    ! mins
    td = tnow/(60*24)       ! days
    if (td < 2) then
        a = 1
    elseif (td < 5) then
        a = 1 - (td-2)/3
    else
        a = 0
    endif
    get_DCactivity = a*NDC
    return
endif

ncap = 0
a = 0
do idc = 1,NDC
    if (DClist(idc)%capable) then
        ncap = ncap + 1
        a = a + DClist(idc)%density
    endif
enddo
get_DCactivity = a
end function

!--------------------------------------------------------------------------------------
! Interim
!--------------------------------------------------------------------------------------
subroutine set_stage(p,stage)
type(cog_type), pointer :: p
integer :: stage
integer :: status, oldstage, region
integer(1) :: statusbyte(4)
equivalence (status,statusbyte)

call get_stage(p,oldstage,region)
status = p%status
!statusbyte(STAGE_BYTE) = stage
statusbyte(STAGE_BYTE) = stage + (region-1)*STAGELIMIT
p%status = status
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine set_stage_region(p,stage,region)
type(cog_type), pointer :: p
integer :: stage, region
integer :: status
integer(1) :: statusbyte(4)
equivalence (status,statusbyte)

status = p%status
statusbyte(STAGE_BYTE) = stage + (region-1)*STAGELIMIT
p%status = status
end subroutine

!----------------------------------------------------------------------------------------
! Interim
!----------------------------------------------------------------------------------------
integer function old_get_stage(p)
type(cog_type), pointer :: p
integer :: status
integer(1) :: statusbyte(4)
equivalence (status,statusbyte)

status = p%status
old_get_stage = statusbyte(STAGE_BYTE)
end function

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine get_stage(p,stage,region)
type(cog_type), pointer :: p
integer :: stage, region
integer :: status
integer(1) :: statusbyte(4)
equivalence (status,statusbyte)

status = p%status
stage = statusbyte(STAGE_BYTE)
region = LYMPHNODE
if (stage > STAGELIMIT) then
	stage = stage - STAGELIMIT
	region = PERIPHERY
endif
!write(logmsg,*) 'get_stage: ',stage,region
!call logger(logmsg)
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine get_region(p,region)
type(cog_type), pointer :: p
integer :: region
integer :: stage, status
integer(1) :: statusbyte(4)
equivalence (status,statusbyte)

status = p%status
stage = statusbyte(STAGE_BYTE)
region = LYMPHNODE
if (stage > STAGELIMIT) then
	region = PERIPHERY
endif
end subroutine

!--------------------------------------------------------------------------------------
! Interim
!--------------------------------------------------------------------------------------
subroutine set_generation(p,gen)
type(cog_type), pointer :: p
integer :: gen
integer :: status
integer(1) :: statusbyte(4)
equivalence (status,statusbyte)

status = p%status
statusbyte(GENERATION_BYTE) = gen
p%status = status
end subroutine

!----------------------------------------------------------------------------------------
! Interim
!----------------------------------------------------------------------------------------
integer function get_generation(p)
type(cog_type), pointer :: p
integer :: status
integer(1) :: statusbyte(4)
equivalence (status,statusbyte)

status = p%status
get_generation = statusbyte(GENERATION_BYTE)
end function

!---------------------------------------------------------------------
!---------------------------------------------------------------------
real function get_IL2store(p)
type(cog_type), pointer :: p
integer :: IL2store_index

if (.not.use_cytokines) then
    get_IL2store = 0
    return
endif
IL2store_index = NP_offset(cyt_seq(IL2_TAG)) + IL2_Store
get_IL2store = p%IL_state(IL2store_index)
end function

!-----------------------------------------------------------------------------------------
! Display the properties of a cognate cell.  The calling program must ascertain that kcell
! is cognate.
!-----------------------------------------------------------------------------------------
subroutine show_cognate_cell(kcell)
integer :: kcell
type (cog_type), pointer :: p
!integer :: cogID
type(cell_type) :: tcell
integer :: gen, stage, region

tcell = cellist(kcell)
if (.not.associated(tcell%cptr)) then
    write(*,*) 'ERROR: show_cognate_cell: cptr not associated: ',kcell
    stop
endif
p => tcell%cptr
write(*,*) 'Cognate cell: ',p%cogID,kcell,cellist(kcell)%ID
write(*,'(a,i10,a,3i4,a,i2)') '  ID: ',tcell%ID,' site: ',tcell%site,' ctype: ',tcell%ctype
write(*,'(a,i2,a,2i6,a,2f10.2)') '  lastdir: ',tcell%lastdir,' DCbound: ',tcell%DCbound,' unbindtime: ',tcell%unbindtime
gen = get_generation(p)
!stage = get_stage(p)
call get_stage(p,stage,region)
write(*,'(a,i8,a,i2,a,i2)') '   cogID: ',p%cogID,' gen: ',gen,' stage: ', stage
write(*,'(a,4f10.2)') '   times: entry,die,div,stage: ',tcell%entrytime,p%dietime,p%dividetime,p%stagetime
write(*,'(a,3f8.2)') 'avidity, stimulation, IL2 store:: ', p%avidity,p%stimulation,get_IL2store(p)

end subroutine

!-----------------------------------------------------------------------------------------
! Called whenever balancer carries out add_sites or removeSites.
! cognate_list(k) = 0 when the kth cognate cell has gone (left or died).
! If Mnodes = 1 this is called only once, after PlaceCells.  After that the list
! is maintained directly when cognate cells arrive or leave.
!-----------------------------------------------------------------------------------------
subroutine make_cognate_list(ok)
logical :: ok
integer :: kcell, ctype, stype, cogID
type (cog_type), pointer :: p

!write(*,*) 'make_cognate_list: ', lastcogID
ok = .true.
cognate_list(1:lastcogID) = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    ctype = cellist(kcell)%ctype
!    stype = struct_type(ctype)
	if (associated(cellist(kcell)%cptr)) then
		stype = COG_TYPE_TAG
	else
		stype = NONCOG_TYPE_TAG
	endif
    if (stype == COG_TYPE_TAG) then
        p => cellist(kcell)%cptr
        cogID = p%cogID
        if (cogID == 0) then
            lastcogID = lastcogID + 1
            if (lastcogID > MAX_COG) then
                write(logmsg,'(a,i6)') 'Error: make_cognate_list: cognate_list dimension exceeded: ',MAX_COG
                call logger(logmsg)
                ok = .false.
                return
            endif
            cogID = lastcogID
            p%cogID = cogID
        endif
        cognate_list(cogID) = kcell
    endif
enddo
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine save_pos
!integer :: k, kcell, site(3)
character*(8) :: msg
integer :: error

!write(nfout,*) 'STEP: ',istep,lastcogID
!do k = 1,lastcogID
!    kcell = cognate_list(k)
!    if (kcell > 0) then
!        site = cellist(kcell)%site
!        gen = get_generation(cellist(kcell)%cptr)
!    else
!        site = 0
!        gen = 0
!    endif
!    write(nfout,'(i8,3i4)') kcell,site,gen
!enddo
if (save_pos_cmgui) then
	call save_exnode
elseif (use_TCP) then
	call save_cell_positions
	msg = 'VTK'
	clear_to_send = .false.
    call winsock_send(awp_1,msg,len_trim(msg),error)
endif
end subroutine

!--------------------------------------------------------------------------------
! We display only T cells that are still in the lymphnode
!--------------------------------------------------------------------------------
subroutine save_cell_positions
!!!use ifport
integer :: k, kcell, site(3), j, idc, dcsite(3)
!integer :: dcstate = 1
real :: dcstate
integer :: itcstate, stype, ctype, stage, region
real :: Tcell_diam = 0.9
real :: DC_diam = 1.8
!real :: spectrum_max = 10, spectrum_freefraction = 0.9
integer :: gen, bnd(2)
logical :: ex
character*(12) :: fname = 'cell_pos.dat'
character*(9) :: removefile = 'TO_REMOVE'

if (simulation_start) then
	inquire(file=fname,exist=ex)
	if (ex) then
		call unlink(fname)
	endif
	inquire(file=removefile,exist=ex)
	if (ex) then
		call unlink(removefile)
	endif
endif
simulation_start = .false.

if (.not.clear_to_send) then
	! wait until the file called removefile exists, then remove it
	inquire(file=removefile,exist=ex)
	if (.not.ex) then
!		call logger('wait')
		do
			inquire(file=removefile,exist=ex)
			if (.not.ex) then
!				call millisleep(10) ! no good at all
			else
				exit
			endif
		enddo
	endif
	call unlink(removefile)
	clear_to_send = .true.
endif


!inquire(file=fname,exist=ex)
!if (ex) then
!	ex = .false.
!!	write(*,*) 'waiting for file removal ...'
!	call logger('waiting for file removal ...')
!	do
!		inquire(file=fname,exist=ex)
!		if (ex) then
!!!!			call sleepqq(100)
!			call logger('millisleep(10)')
!			call millisleep(10)
!		else
!			exit
!		endif
!	enddo
!!	write(*,*) 'file removed'
!	call logger('file removed')
!endif

open(nfpos,file=fname,status='new')
! DC section
if (NDC > 0) then
    do k = 1,NDC
        if (DClist(k)%alive) then
            site = DClist(k)%site
            dcstate = min(1.0,DClist(k)%density/DC_ANTIGEN_MEDIAN)
!            dcstate = 1
            ! Need dcstate to convey antigen density level (normalized to 0-1)
            write(nfpos,'(a2,i4,3i4,f4.1,f5.2)') 'D ',k-1, site, DC_diam, dcstate
!            write(logmsg,'(a2,i4,3i4,f4.1,f5.2)') 'D ',k-1, site, DC_diam, dcstate
!            call logger(logmsg)
        endif
    enddo
endif

if (.not.IV_SHOW_NONCOGNATE) then
	! T cell section
	do k = 1,lastcogID
		kcell = cognate_list(k)
		if (kcell > 0) then
			call get_stage(cellist(kcell)%cptr,stage,region)
			if (region /= LYMPHNODE) cycle
			site = cellist(kcell)%site
	!        tcstate = mod(kcell,2) + 1
			gen = get_generation(cellist(kcell)%cptr)
			bnd = cellist(kcell)%DCbound
	!        if (bnd(1) == 0 .and. bnd(2) == 0) then
	!            tcbound = 0
	!        else
	!            tcbound = 1
	!        endif
!			if (get_stage(cellist(kcell)%cptr) == NAIVE) then
			if (stage == NAIVE) then
				itcstate = 0
			else
				if (bnd(1) == 0 .and. bnd(2) == 0) then
	!				tcstate = (gen-1.0)/(TC_MAX_GEN-1.0)*spectrum_max*spectrum_freefraction
	!				tcstate = (gen/TC_MAX_GEN)*spectrum_max*spectrum_freefraction
					itcstate = gen
				else
	!				tcstate = spectrum_max
					itcstate = 99
				endif
			endif
			! Need tcstate to convey non-activated status, i.e. 0 = non-activated
			write(nfpos,'(a2,i6,3i4,f4.1,i3)') 'T ',k-1, site, Tcell_diam, itcstate
		endif
	enddo
	! Bond section
	do k = 1,lastcogID
		kcell = cognate_list(k)
		if (kcell > 0) then
			call get_stage(cellist(kcell)%cptr,stage,region)
			if (region /= LYMPHNODE) cycle
			site = cellist(kcell)%site
			do j = 1,2
				idc = cellist(kcell)%DCbound(j)
				if (idc /= 0) then
					if (DClist(idc)%capable) then
						dcsite = DClist(idc)%site
	!                    dcstate = mod(idc,2) + 1
						dcstate = 1
	!                    write(nfpos,'(a,i8,3i4,4i4)') 'Node: ',nd, site, dcsite-site, dcstate
						write(nfpos,'(a2,2i5)') 'B ',k-1,idc-1
					endif
				endif
			enddo
		endif
	enddo
else
	! T cell section
	do kcell = 1,nlist
		if (cellist(kcell)%ID == 0) cycle  ! gap
		site = cellist(kcell)%site
		bnd = cellist(kcell)%DCbound
		ctype = cellist(kcell)%ctype
!		stype = struct_type(ctype)
		if (associated(cellist(kcell)%cptr)) then
			stype = COG_TYPE_TAG
		else
			stype = NONCOG_TYPE_TAG
		endif
		if (stype == NONCOG_TYPE_TAG) then
			itcstate = -1
		else
			call get_stage(cellist(kcell)%cptr,stage,region)
			if (region /= LYMPHNODE) cycle
			gen = get_generation(cellist(kcell)%cptr)
!			if (get_stage(cellist(kcell)%cptr) == NAIVE) then
			if (stage == NAIVE) then
				itcstate = 0
			else
				if (bnd(1) == 0 .and. bnd(2) == 0) then
					itcstate = gen
				else
					itcstate = 99
				endif
			endif
		endif
		! Need tcstate to convey non-activated status, i.e. 0 = non-activated
		write(nfpos,'(a2,i6,3i4,f4.1,i3)') 'T ',kcell-1, site, Tcell_diam, itcstate
	enddo
	! Bond section
	do kcell = 1,nlist
		if (cellist(kcell)%ID == 0) cycle  ! gap
		if (associated(cellist(kcell)%cptr)) then
			call get_stage(cellist(kcell)%cptr,stage,region)
			if (region /= LYMPHNODE) cycle
		endif
		site = cellist(kcell)%site
		do j = 1,2
			idc = cellist(kcell)%DCbound(j)
			if (idc /= 0) then
				if (DClist(idc)%capable) then
					dcsite = DClist(idc)%site
!                    dcstate = mod(idc,2) + 1
					dcstate = 1
!                    write(nfpos,'(a,i8,3i4,4i4)') 'Node: ',nd, site, dcsite-site, dcstate
					write(nfpos,'(a2,2i5)') 'B ',kcell-1,idc-1
				endif
			endif
		enddo
	enddo
endif
write(nfpos,'(a2,i6)') 'E ',istep
close(nfpos)

end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine save_exnode
integer :: k, kcell, site(3), nd, j, idc, dcsite(3)
integer :: dcstate = 1
real :: tcstate
real :: Tcell_diam = 0.9
real :: DC_diam = 1.8
real :: spectrum_max = 10, spectrum_freefraction = 0.6
integer :: gen, stage, region, bnd(2)
character*(64) :: fname = '\CMGUI\DCU\dcu_00000.exnode'

write(fname(16:20),'(i5.5)') istep
write(*,*) fname
open(nfcmgui,file=fname,status='replace')

nd = 0

! DC section
if (NDC > 0) then
    write(nfcmgui,*) 'Group name : DC'
    write(nfcmgui,*) '  #Fields=3'
    write(nfcmgui,*) '  1) coordinates, coordinate, rectangular cartesian, #Components=3'
    write(nfcmgui,*) '    x.  Value index= 1, #Derivatives= 0'
    write(nfcmgui,*) '    y.  Value index= 2, #Derivatives= 0'
    write(nfcmgui,*) '    z.  Value index= 3, #Derivatives= 0'
    write(nfcmgui,*) '  2) diameter, field, real, #Components=1'
    write(nfcmgui,*) '    value.'
    write(nfcmgui,*) '  3) dcstate, field, real, #Components=1'
    write(nfcmgui,*) '    value.'

    do k = 1,NDC
        if (DClist(k)%alive) then
            nd = nd+1
            site = DClist(k)%site
            dcstate = 1
            write(nfcmgui,'(a,i8,3i4,f4.1,i3)') 'Node: ',nd, site, DC_diam, dcstate
        endif
    enddo
endif

! T cell section
write(nfcmgui,*) 'Group name : Tcell'
write(nfcmgui,*) '  #Fields=3'
write(nfcmgui,*) '  1) coordinates, coordinate, rectangular cartesian, #Components=3'
write(nfcmgui,*) '    x.  Value index= 1, #Derivatives= 0'
write(nfcmgui,*) '    y.  Value index= 2, #Derivatives= 0'
write(nfcmgui,*) '    z.  Value index= 3, #Derivatives= 0'
write(nfcmgui,*) '  2) diameter, field, real, #Components=1'
write(nfcmgui,*) '    value.'
write(nfcmgui,*) '  3) tcstate, field, real, #Components=1'
write(nfcmgui,*) '    value.'
!write(nfcmgui,*) '  3) tcgen, field, integer, #Components=1'
!write(nfcmgui,*) '    value.'
!write(nfcmgui,*) '  4) tcbound, field, integer, #Components=1'
!write(nfcmgui,*) '    value.'

do k = 1,lastcogID
    kcell = cognate_list(k)
    if (kcell > 0) then
		call get_stage(cellist(kcell)%cptr,stage,region)
		if (region /= LYMPHNODE) cycle
        nd = nd+1
        site = cellist(kcell)%site
!        tcstate = mod(kcell,2) + 1
        gen = get_generation(cellist(kcell)%cptr)
        bnd = cellist(kcell)%DCbound
!        if (bnd(1) == 0 .and. bnd(2) == 0) then
!            tcbound = 0
!        else
!            tcbound = 1
!        endif
        if (bnd(1) == 0 .and. bnd(2) == 0) then
            tcstate = (gen-1.0)/(TC_MAX_GEN-1.0)*spectrum_max*spectrum_freefraction
        else
            tcstate = spectrum_max
        endif
!        write(nfcmgui,'(a,i8,3i4,f4.1,i3,i2)') 'Node: ',nd, site, Tcell_diam, gen, tcbound
        write(nfcmgui,'(a,i8,3i4,f4.1,f6.2)') 'Node: ',nd, site, Tcell_diam, tcstate
    endif
enddo

! Bond section
write(nfcmgui,*) 'Group name : Bond'
write(nfcmgui,*) '  #Fields=3'
write(nfcmgui,*) '  1) coordinates, coordinate, rectangular cartesian, #Components=3'
write(nfcmgui,*) '    x.  Value index= 1, #Derivatives= 0'
write(nfcmgui,*) '    y.  Value index= 2, #Derivatives= 0'
write(nfcmgui,*) '    z.  Value index= 3, #Derivatives= 0'
write(nfcmgui,*) '  2) vector, coordinate, rectangular cartesian, #Components=3'
write(nfcmgui,*) '    x.  Value index= 1, #Derivatives= 0'
write(nfcmgui,*) '    y.  Value index= 2, #Derivatives= 0'
write(nfcmgui,*) '    z.  Value index= 3, #Derivatives= 0'
write(nfcmgui,*) '  3) dcstate, field, real, #Components=1'
write(nfcmgui,*) '    value.'

do k = 1,lastcogID
    kcell = cognate_list(k)
    if (kcell > 0) then
        site = cellist(kcell)%site
        do j = 1,2
            idc = cellist(kcell)%DCbound(j)
            if (idc /= 0) then
                if (DClist(idc)%capable) then
                    nd = nd+1
                    dcsite = DClist(idc)%site
!                    dcstate = mod(idc,2) + 1
                    dcstate = 1
                    write(nfcmgui,'(a,i8,3i4,4i4)') 'Node: ',nd, site, dcsite-site, dcstate
                endif
            endif
        enddo
    endif
enddo

close(nfcmgui)
end subroutine

!--------------------------------------------------------------------------------
! The ability of a DC to provide TCR stimulation is conveyed by the current
! value of antigen density, %density.  The level decays over time, and possibly
! is reduced by every episode of TCR stimulation, at the rate DC_DENS_BY_STIM.
! When t > %dietime - DC_ACTIV_TAPER the density reduces each time step by a factor
! such that when %dietime is reached the product of the factors = 0.1
! tfactor**n = 0.1 where n = DC_ACTIV_TAPER/DELTA_T
! tfactor = (0.1)**(1/n)
! A DC loses its TCR stimulating capability when %density falls below the
! limiting value DC_pMHC_THRESHOLD.
! A DC is scheduled to die when t > %dietime.  At this point %capability is set
! to .false., but the DC waits until %nbound = 0 before it dies.
!--------------------------------------------------------------------------------
subroutine update_DCstate(ok)
logical :: ok
real :: tnow, tfactor, decay_factor
integer :: idc, nalive, nbound, ncbound, ncapable

if (NDCalive == 0) then
	NDCcapable = 0
	ok = .true.
	return
endif
ncapable = 0
tfactor = 0.1**(DELTA_T/(60*DC_ACTIV_TAPER))
tnow = istep*DELTA_T
decay_factor = 1 - DCdecayrate*DELTA_T
nalive = 0
do idc = 1,NDC
    if (DClist(idc)%alive) then
        if (.not.IN_VITRO .and. DC_outside(idc)) then	! A kludge to remove DCs stranded outside the blob
            DClist(idc)%dietime = min(DClist(idc)%dietime,tnow)
		endif
        nbound = DClist(idc)%nbound
        ncbound = DClist(idc)%ncogbound
        if (save_DCbinding) then
            dcbind(ncbound) = dcbind(ncbound) + 1
        endif
        if (tnow > DClist(idc)%dietime) then
            DClist(idc)%capable = .false.
            write(*,*) 'update_DCstate: not capable: ',idc,nbound
            if (nbound /= 0) then
            elseif (nbound == 0) then
				write(logmsg,'(a,i4,f8.2,i6)') 'DC dies: ',idc,DClist(idc)%dietime,NDCalive
!				write(nflog,'(a,i4,f8.2)') 'DC dies: ',idc,DClist(idc)%dietime
				call logger(logmsg)
                DClist(idc)%alive = .false.
                ndeadDC = ndeadDC + 1
                DCdeadlist(ndeadDC) = idc
                ! later we need to adjust occupancy to reflect the death of DC idc
            else
                write(logmsg,'(a,2i6)') 'Error: update_DCstatus: nbound < 0: ',idc,nbound
                call logger(logmsg)
				ok = .false.
				return
            endif
        endif
        if (DClist(idc)%capable) then
            DClist(idc)%density = DClist(idc)%density*decay_factor - DC_DENS_BY_STIM*DClist(idc)%stimulation
            DClist(idc)%stimulation = 0     ! reset the tally of stimulation delivered to zero (OMP)
            if (tnow > DClist(idc)%dietime - 60*DC_ACTIV_TAPER) then	! antigen density decays over DC_ACTIV_TAPER
                DClist(idc)%density = DClist(idc)%density*tfactor
            endif
            if (DClist(idc)%density < DC_pMHC_THRESHOLD) then
                DClist(idc)%capable = .false.
                if (incapable_DC_dies) then
                    ! A non-capable DC might as well die - but this would eliminate low-antigen DCs
                    DClist(idc)%dietime = tnow + 60     ! give it 1 hour to die
                endif
!                write(*,*) 'DC incapable: ',idc
            endif
        endif
    endif
    if (DClist(idc)%alive) nalive = nalive + 1
    if (DClist(idc)%capable) ncapable = ncapable + 1
enddo
NDCalive = nalive
NDCcapable = ncapable
!write(nflog,*) 'update_DCstate: live DCs: ',nalive
if (nalive == 0) then
    write(logmsg,*) 'No live DC'
    call logger(logmsg)
endif
!if (ncapable == 0) then
!	call logger("No capable DC")
!	return
!endif
ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
! If a DC is left outside the blob as it contracts, the simplest thing to do is to make
! it die.  Better if the DC moves.
!-----------------------------------------------------------------------------------------
logical function DC_outside(idc)
integer :: idc
real :: d

d = cdistance(DClist(idc)%site)
if (d - Radius > 1) then
	DC_outside = .true.
else
	DC_outside = .false.
endif
end function

!-----------------------------------------------------------------------------------------
! When a DC dies it is labelled as dead in the DC list (%alive = .false.) and the count
! of DCs that have died is incremented.
! Later, after gather_data(), the DCs that have died are cleared out of occupancy()%indx().
! The reassignment of occupancy()%DC values is done all at once in
! reassign_DC(), called from balancer().
!-----------------------------------------------------------------------------------------
subroutine clearDC(idc)
integer :: idc
integer :: k, site0(3), site(3)

site0 = DClist(idc)%site
occupancy(site0(1),site0(2),site0(3))%indx = 0
occupancy(site0(1),site0(2),site0(3))%DC = 0
if (NDCsites > 1) then  ! need to clear non-central sites occupied by the DC
    do k = 2,NDCsites
        site = site0 + DCoffset(:,k)
        if (occupancy(site(1),site(2),site(3))%indx(1) == -idc) then
            occupancy(site(1),site(2),site(3))%indx = 0
            occupancy(site(1),site(2),site(3))%DC = 0
        endif
    enddo
endif
DClist(idc)%ID = 0
DClist(idc)%density = 0
DClist(idc)%stimulation = 0
DClist(idc)%dietime = 0
DClist(idc)%site = 0
DClist(idc)%nbound = 0
DClist(idc)%ncogbound = 0
deallocate(DClist(idc)%cogbound)
if (IN_VITRO) then	! Need to adjust zrange2D()
	call resetDCzrange(site0(1),site0(2))
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Set the zrange2D() values for sites near a DC at (xdc,ydc).
! For now, just use DCstack().
! Note that zrange2D was initialized to (1,1) for all interior sites before DC placement.
! Since the ranges of influence of DCs may overlap, both lower and upper z limits can
! only increase those values set for previously encountered DCs.
!-----------------------------------------------------------------------------------------
subroutine setDCzrange(xdc,ydc)
integer :: xdc, ydc
integer :: x,y,z,dx,dy,zval(3)

do dx = -4,4
	do dy = -4,4
		x = xdc + dx
		y = ydc + dy
		zval = DCstack(dx,dy,:)
		if (zval(1) == 0) cycle
		do z = 1,3
			if (zval(z) > 0) then
				zrange2D(x,y,1) = max(z,zrange2D(x,y,1))
				exit
			endif
		enddo
		do z = 3,1,-1
			if (zval(z) > 0) then
				zrange2D(x,y,2) = max(z,zrange2D(x,y,2))
				exit
			endif
		enddo
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! When a DC dies the zrange2D() values in its vicinity must be adjusted.
!-----------------------------------------------------------------------------------------
subroutine resetDCzrange(xdc,ydc)
integer :: xdc, ydc
integer :: x,y,z,dx,dy,zval(3)

do dx = -4,4
	do dy = -4,4
		x = xdc + dx
		y = ydc + dy
		zval = DCstack(dx,dy,:)
		if (zval(1) == 0) cycle
		if (zrange2D(x,y,2) == 1) cycle
		do z = 3,1,-1
			if (zval(z) > 1) then
				if (zval(z) == zrange2D(x,y,2)) then
					zrange2D(x,y,:) = 1
					write(nflog,*) 'resetDCzrange: ',xdc,ydc,dx,dy,zval(z)
					exit
				endif
			endif
		enddo
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
! Determine the number of DC that enter the paracortex between times t1 and t2.
! For now just use a rate of arrivals/min., declining to zero between times
! tdc1, tdc2 (hours).
! This rate must be scaled to be proportional to the steady-state number of
! T cells, i.e. NTcells0
! A rate of 0.1 is not bad for NTcells0 = 100k
! It also must be scaled by the specified rate of occurrence of DCs, i.e. to
! take account of NT_TO_DC.  The value NT_TO_DC = 200 can be used as a baseline.
! THIS IS HARD_WIRED, NEEDS TO DEPEND ON STATE OF INFECTION IN THE TISSUE
! In the case of DCs injected into experimental animals, the schedule of DCs
! estimated to enter the paracortex in each hour is read from an input file.
! The number of DCs entering is precomputed for each hour, scaled by the initial
! T cell population.
!--------------------------------------------------------------------------------
integer function DCinflux(t1,t2,kpar)
integer :: kpar
integer :: ih1, ih2
real :: t1, t2
real :: rate, temp, dn
real :: DCrate0
real(DP) :: R
real, save :: dn_last = 0

!write(logmsg,*) 'DCinflux: dn_last ',dn_last
!call logger(logmsg)
if (DC_INJECTION) then
	ih1 = t1/60
	ih2 = t2/60
	if (ih1 < ih2) then
		DCinflux = DCinjected(ih2)
	else
		DCinflux = 0
	endif
else
	DCrate0 = DC_FACTOR*DCrate_100k*(NTcells0/1.0e5)     ! DC influx is scaled by the initial T cell population
	!DCrate0 = DCrate0*(200./TC_TO_DC)                              ! and scaled by 1/TC_TO_DC (TRY CANCELLING THIS)
	if (t1 > T_DC2*60) then
		DCinflux = 0
		return
	elseif (t1 < T_DC1*60) then
		rate = DCrate0      ! rate /min
	else
		rate = DCrate0*(T_DC2 - t1/60)/(T_DC2 - T_DC1)
	endif
	if (RANDOM_DCFLUX) then
		temp = rate*(t2-t1)
		DCinflux = temp
		dn = temp - DCinflux
		R = par_uni(kpar)
		if (R < dn) then
			DCinflux = DCinflux + 1
		endif
	else
		! Try making the DC influx deterministic
		temp = rate*(t2-t1) + dn_last
		DCinflux = temp
		dn_last = temp - DCinflux
	endif
endif
!write(logmsg,*) 'DC rate: ',rate,DCinflux,dn_last
!call logger(logmsg)
end function

!--------------------------------------------------------------------------------
! Make a random choice of an integer from 1 - N on the basis of probabilities in
! the array p(:) (assumed to be normalized).
!--------------------------------------------------------------------------------
integer function random_choice(p,N,kpar)
integer :: N,kpar
real(DP) :: p(:)
integer :: k
real(DP) :: R, psum

!call random_number(R)
R = par_uni(kpar)
psum = 0
do k = 1,N
    psum = psum + p(k)
    if (R <= psum) then
        random_choice = k
        return
    endif
enddo
write(logmsg,*) 'ERROR: random_choice: ',N,p
call logger(logmsg)
stop
end function

!-----------------------------------------------------------------------------------------
! Given c (0 < c < 0.5), the fraction of the sphere volume that is cut off, compute
! the distance of the cut from the centre, d, as x = d/R, where R = sphere radius.
!-----------------------------------------------------------------------------------------
subroutine get_slice(c,x)
real :: c,x,xp
real, parameter :: epsilon = 0.0001
integer :: k

x = 0.5
xp = x
k = 0
do
    k = k+1
    if (k > 100) then
        write(logmsg,*) 'ERROR: get_slice: failed to converge'
		call logger(logmsg)
        stop
    endif
    x = 1 - 4*c/(2-(1+x)*x**3)
    if (abs(x-xp) < epsilon) exit
    xp = x
enddo

end subroutine

!--------------------------------------------------------------------------------
! The chemotactic step probabilities for all possible sites within an exit's
! SOI are precomputed.  The array is indexed by the offset of each site from
! the attractant site, given by (x,y,z).
! CHANGED now array is indexed by the offset of the attractant site from the site.
! This makes more sense when the vector v is the chemokine concentration gradient.
! Note that the weight given to a jump direction is found from the cosine^2 of
! the angle between the jump direction and the vector from the site to the
! attractant site.
!--------------------------------------------------------------------------------
subroutine chemo_setup
integer :: x, y, z, k, r(3), s(3)
real :: r2, s2, rmod,smod,cosa
real, allocatable :: w(:)

write(logmsg,*) 'chemo_setup'
call logger(logmsg)
allocate(chemo_r(0:chemo_N,0:chemo_N,0:chemo_N))
allocate(chemo_p(-chemo_N:chemo_N,-chemo_N:chemo_N,-chemo_N:chemo_N,njumpdirs))
allocate(w(njumpdirs))

do x = -chemo_N,chemo_N
    do y = -chemo_N,chemo_N
        do z = -chemo_N,chemo_N
            r = (/x,y,z/)
            r2 = dot_product(r,r)
            rmod = sqrt(r2)
            w = 0
            do k = 1,njumpdirs
                if (k == 14) cycle
                s = jumpvec(:,k)
                s2 = dot_product(s,s)
                smod = sqrt(s2)
                cosa = dot_product(r,s)/(rmod*smod)
!                if (cosa < 0) then
                if (cosa > 0) then
                    w(k) = cosa*cosa/smod
                endif
            enddo
            w = w/sum(w)
            chemo_p(x,y,z,:) = w	! Note that these probs sum to 1
            if (x >= 0 .and. y >= 0 .and. z >= 0) then
                chemo_r(x,y,z) = sqrt(r2)
            endif
        enddo
    enddo
enddo
deallocate(w)
end subroutine

!--------------------------------------------------------------------------------
! Computes the jump probabilities (absolute directions) accounting for chemotaxis
! towards exit (or DC).  On input:
!   p(:) holds the jump probabilities not accounting for  chemotaxis
!   v(:) is the site offset relative to the exit (or DC)
! CHANGED
!   v(:) is offset the exit (or DC) offset relative to the site
! OR
!   v(:) is the concentration gradient vector
!   f is the amount of chemotactic influence
!   c is the amount of CCR7 ligand influence (Cyster).
! Note that f incorporates both the distance from the exit (or DC) and the cell's
! susceptibility to chemotaxis, which may be the S1P1 level of the T cell.
! On return p(:) holds the modified jump probabilities.
! Note: should have p remaining unchanged as chemo_K -> 0
! Note: when njumpdirs = 27, jump 14 corresponds to (0,0,0) - unused.
! Note: code modifications now base v and f on net chemotactic attraction of 
! multiple exits and DCs.
!--------------------------------------------------------------------------------
subroutine chemo_probs(p,v,f,c)
real(DP) :: p(:)
integer :: v(:)
real(DP) :: f, c
integer :: k
real(DP) :: pc(MAXRELDIR+1)

if (f == 0 .and. c == 1) then
    return
endif
if (f > 1) then
	write(logmsg,*) 'ERROR: chemo_probs: f > 0: ',f
	call logger(logmsg)
	return
endif
p = p/sum(p)
if (f > 0) then
    pc(1:njumpdirs) = chemo_p(v(1),v(2),v(3),:)
else
    pc = 0
endif
do k = 1,njumpdirs
    if (p(k) > 0) then      ! prevents jumps in disallowed directions
        p(k) = (1-f)*c*p(k) + f*pc(k)
    endif
enddo
end subroutine

!--------------------------------------------------------------------------------
! Computes the jump probabilities (absolute directions) accounting for chemotaxis
! On input:
!   p(:) holds the jump probabilities not accounting for  chemotaxis
!   v(:) was the site offset relative to the exit, used only to approximate direction
!        by allowing a discrete number of directions (chemo_N determines these, in fact
!        only a small subset of the array positions are used - those roughly falling
!        on a sphere of radius chemo_N)
!   f is the amount of chemotactic influence
! Note that f incorporates both the magnitude of chemokine gradients and the cell's
! susceptibility to chemotaxis.
! On return p(:) holds the modified jump probabilities.
! Note: when njumpdirs = 27, jump 14 corresponds to (0,0,0) - unused.
! Note: code modifications now base v and f on net chemotactic attraction of 
! multiple exits and DCs.
!--------------------------------------------------------------------------------
subroutine chemo_probs_pre(p,v,f)
real(DP) :: p(:)
integer :: v(:)
real :: f
integer :: k
real(DP) :: pc(MAXRELDIR+1)

if (f == 0) then
    return
endif
if (f > 1) then
	write(logmsg,*) 'ERROR: chemo_probs: f > 0: ',f
	call logger(logmsg)
	return
endif
p = p/sum(p)
pc(1:njumpdirs) = chemo_p(v(1),v(2),v(3),:)
do k = 1,njumpdirs
    if (p(k) > 0) then      ! prevents jumps in disallowed directions
        p(k) = (1-f)*p(k) + f*pc(k)
    endif
enddo
end subroutine

!--------------------------------------------------------------------------------
! Computes the functional dependence of chemotactic effect on distance r from
! the exit.  r is in units of site spacing.
!--------------------------------------------------------------------------------
real function chemo_g(r)
real :: r

chemo_g = min(1.0,(1.0/r)**chemo_exp)
end function

!--------------------------------------------------------------------------------
! Determines the degree to which a cell is subject to chemotaxis.
! (Determined by CD69 level, or S1P1 level.)
! If use_exit_chemotaxis is true, i.e. chemotaxis is used to control cell exit, any cell
! that gets close enough to the exit will leave the paracortex.  Is this
! acceptable?
! For a noncognate cell, should rise from 0 to 1 in an hour or so.
!--------------------------------------------------------------------------------
real function chemo_active_exit(cell)
type(cell_type), pointer :: cell
real :: tnow, t

if (turn_off_chemotaxis) then
    chemo_active_exit = 0
    return
endif

if (RELAX_INLET_EXIT_PROXIMITY) then
    tnow = istep*DELTA_T
    t = tnow - cell%entrytime
    if (t > CHEMO_K_RISETIME) then
		chemo_active_exit = 1
	else
		chemo_active_exit = t/CHEMO_K_RISETIME
	endif
	return
endif

if (TAGGED_EXIT_CHEMOTAXIS) then
    tnow = istep*DELTA_T
    t = tnow - cell%entrytime
    if (cell%tag == RES_TAGGED_CELL) then     ! testing effect of S1P1
        chemo_active_exit = (1 - exp(-K1_S1P1*t))*TAGGED_CHEMO_ACTIVITY
    else
		chemo_active_exit = 0
	endif
	return
endif

if (TAGGED_LOG_PATHS) then
	if (cell%tag == CHEMO_TAGGED_CELL) then
		chemo_active_exit = CHEMO_K_EXIT
	else
		chemo_active_exit = 0
	endif
	return
endif

if (associated(cell%cptr)) then     ! cognate cell
    chemo_active_exit = cell%cptr%S1P1
!    if (CD69 < CD69_threshold) then
!        chemo_active = 1 - CD69/CD69_threshold
!    else
!        chemo_active = 0
!    endif
else
    tnow = istep*DELTA_T
    t = tnow - cell%entrytime
!    if (cell%tag == RES_TAGGED_CELL) then     ! testing effect of S1P1
!        chemo_active = 1 - exp(-K1_S1P1*0.1*t)
!        chemo_active = 0
!    else
        chemo_active_exit = 1 - exp(-K1_S1P1*t)
!    endif
endif
end function

!--------------------------------------------------------------------------------
! Returns the level of CCR7 ligand (i.e. CCL19/21) at a distance r sites from an
! exit site, within the exit SOI.  The value must be in the range (0,1).
! The parameters CCR7_R1 and CCR7_R2 specify the range of r over which the
! ligand level ranges linearly from 0 to 1.  This is a simple way to program a
! decreased level of CCR7 near an exit, following Cyster's ideas.
! Note: currently this always returns 1
!--------------------------------------------------------------------------------
real function CCR7_ligand(r)
real :: r
real, parameter :: CCR7_R1 = 0, CCR7_R2 = 0

if (r < CCR7_R1) then
    CCR7_ligand = 0
elseif (r < CCR7_R2) then
    CCR7_ligand = (r - CCR7_R1)/(CCR7_R2 - CCR7_R1)
else
    CCR7_ligand = 1
endif
end function

!--------------------------------------------------------------------------------
! Returns cell's susceptibility to DC chemotaxis
!--------------------------------------------------------------------------------
real function chemo_active_DC(cell)
type(cell_type), pointer :: cell

if (turn_off_chemotaxis) then
    chemo_active_DC = 0
    return
endif
!if (associated(cell%cptr)) then     ! cognate cell
!    chemo_active_DC = cell%cptr%DCchemo
!else
!	chemo_active_DC = BASE_DCchemo
!endif
!chemo_active_DC = cell%DCchemo
chemo_active_DC = cell%receptor_level(CCR1)
end function

!--------------------------------------------------------------------------------
! Determines e(:), the location of the nearest exit, if there is a close one.
! Otherwise returns e(:) = 0.
! Note: currently only one exit number is stored in occupancy()
!--------------------------------------------------------------------------------
subroutine nearest_exit(site,in_exit_SOI,e)
integer :: site(3), e(3)
logical :: in_exit_SOI
integer :: iexit

iexit = occupancy(site(1),site(2),site(3))%exitnum
if (iexit == 0) then
    in_exit_SOI = .false.
    e = (/0,0,0/)
    return
endif
in_exit_SOI = .true.
e = exitlist(abs(iexit))%site
!write(*,*) 'exit site: ',e
end subroutine

!--------------------------------------------------------------------------------
! The criterion for a near exit is based on chemo_radius
!--------------------------------------------------------------------------------
subroutine near_exits(site,ne,ee)
integer :: site(3), ne, ee(3,*)
integer :: iexit

iexit = occupancy(site(1),site(2),site(3))%exitnum
if (iexit == 0) then
	ne = 0
    return
endif
ne = 1
ee(:,1) = exitlist(abs(iexit))%site
!write(*,*) 'exit site: ',e
end subroutine

!--------------------------------------------------------------------------------
! This determines which DCs (if any) the site is within chemotactic range of.
! The criterion for a near DC is based on first sites within DC_RADIUS (%DC)
! then sites not in the %DC list but within chemo_radius.
!--------------------------------------------------------------------------------
subroutine near_DCs(site,nd,near)
integer :: site(3), nd, near(*)
integer :: i
integer(2) :: DC(0:max(DCDIM,cDCDIM)-1)     ! DC(0) = number of near DCs, DC(k) = ID of kth near DC

nd = 0
DC(0:DCDIM-1) = occupancy(site(1),site(2),site(3))%DC
if (DC(0) /= 0) then
	do i = 1,DC(0)
		if (DClist(DC(i))%alive) then
			nd = nd + 1
			near(nd) = DC(i)
			if (nd == MAX_DC_CHEMO) return
		endif
	enddo
endif
DC(0:cDCDIM-1) = occupancy(site(1),site(2),site(3))%cDC
if (DC(0) /= 0) then
	do i = 1,DC(0)
		if (DClist(DC(i))%alive) then
			nd = nd + 1
			near(nd) = DC(i)
			if (nd == MAX_DC_CHEMO) return
		endif
	enddo
endif
end subroutine

!----------------------------------------------------------------------------------------
! Convert a half life in hours to a decay rate /min
!----------------------------------------------------------------------------------------
real function DecayRate(halflife)
real :: halflife

DecayRate = log(2.0)/(halflife*60)    ! rate/min
end function


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine logger(msg)
character*(*) :: msg
integer :: error
logical :: isopen
character*(1) :: LF = char(94)

error = 0
if (use_TCP) then
    if (awp_0%is_open) then
        call winsock_send(awp_0,trim(msg)//LF,len_trim(msg)+1,error)
    endif
else
	write(*,*) trim(msg)
endif
inquire(unit=nflog,OPENED=isopen)
if (isopen) then
	write(nflog,*) 'msg: ',trim(msg)
	if (error /= 0) then
	    write(nflog,'(a,i4)') 'winsock_send error: ',error
	    close(nflog)
	endif
endif
if (error /= 0) stop
end subroutine

!--------------------------------------------------------------------------------
!     NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTH'S PASCAL
!     BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.
!     SINGLE PRECISION, ALSO CHANGES THE ORDER OF THE ASSOCIATED ARRAY T.
!--------------------------------------------------------------------------------
SUBROUTINE qsort(a, n, t)
IMPLICIT NONE

INTEGER, INTENT(IN)    :: n
REAL, INTENT(INOUT)    :: a(n)
INTEGER, INTENT(INOUT) :: t(n)

!     Local Variables

INTEGER                :: i, j, k, l, r, s, stackl(15), stackr(15), ww
REAL                   :: w, x

s = 1
stackl(1) = 1
stackr(1) = n

!     KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.

10 CONTINUE
l = stackl(s)
r = stackr(s)
s = s - 1

!     KEEP SPLITTING A(L), ... , A(R) UNTIL L >= R.

20 CONTINUE
i = l
j = r
k = (l+r) / 2
x = a(k)

!     REPEAT UNTIL I > J.

DO
  DO
    IF (a(i).LT.x) THEN                ! Search from lower end
      i = i + 1
      CYCLE
    ELSE
      EXIT
    END IF
  END DO

  DO
    IF (x.LT.a(j)) THEN                ! Search from upper end
      j = j - 1
      CYCLE
    ELSE
      EXIT
    END IF
  END DO

  IF (i.LE.j) THEN                     ! Swap positions i & j
    w = a(i)
    ww = t(i)
    a(i) = a(j)
    t(i) = t(j)
    a(j) = w
    t(j) = ww
    i = i + 1
    j = j - 1
    IF (i.GT.j) EXIT
  ELSE
    EXIT
  END IF
END DO

IF (j-l.GE.r-i) THEN
  IF (l.LT.j) THEN
    s = s + 1
    stackl(s) = l
    stackr(s) = j
  END IF
  l = i
ELSE
  IF (i.LT.r) THEN
    s = s + 1
    stackl(s) = i
    stackr(s) = r
  END IF
  r = j
END IF

IF (l.LT.r) GO TO 20
IF (s.NE.0) GO TO 10

RETURN
END SUBROUTINE qsort

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_xyz(k)
integer :: k
integer :: kcell, xyzsum(3)

xyzsum = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    xyzsum = xyzsum + cellist(kcell)%site
enddo
write(nfres,'(2i6,3i12)') istep,k,xyzsum
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine checkExits(msg)
character*(*) :: msg
integer :: iexit,site(3)
logical :: ok = .true.

write(logmsg,'(a,i8,a,a)') 'checkExits: ',istep,'  ',msg
call logger(logmsg)
do iexit = 1,lastexit
    if (exitlist(iexit)%ID == 0) cycle
	site = exitlist(iexit)%site
	if (occupancy(site(1),site(2),site(3))%exitnum /= -exitlist(iexit)%ID) then
		write(logmsg,*) 'checkExits: ',iexit,site,occupancy(site(1),site(2),site(3))%exitnum
		call logger(logmsg)
		ok = .false.
	endif
enddo
if (.not.ok) stop
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine check_exit(iexit)
integer :: iexit
integer :: x,y,z,a,nc,nctot(0:10),site(3)
integer :: n = 5

site = exitlist(iexit)%site
nctot(0) = sitecells(site,0,0,0)
do a = 1,n
    nc = 0
    do z = -a,a,2*a
        do x = -a,a
            do y = -a,a
                nc = nc + sitecells(site,x,y,z)
            enddo
        enddo
    enddo
    do z = -a+1,a-1
        y = -a
        do x = -a+1,a
            nc = nc + sitecells(site,x,y,z)
        enddo
        x = a
        do y = -a+1,a
            nc = nc + sitecells(site,x,y,z)
        enddo
        y = a
        do x = -a,a-1
            nc = nc + sitecells(site,x,y,z)
        enddo
        x = -a
        do y = -a,a-1
            nc = nc + sitecells(site,x,y,z)
        enddo
    enddo
    nctot(a) = nc
enddo
write(nfout,'(10i6)') nctot(0:n),sum(nctot(0:n))
end subroutine

!--------------------------------------------------------------------------------
! Returns the number of cells on the site offset by (x,y,z) from site(:)
!--------------------------------------------------------------------------------
integer function sitecells(site,x,y,z)
integer :: site(3),x,y,z
integer :: k,indx(2)

indx = occupancy(site(1)+x,site(2)+y,site(3)+z)%indx
sitecells = 0
do k = 1,2
    if (indx(k) > 0) sitecells = sitecells + 1
enddo
end function

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine check_zdistribution
real :: tot(NZ)
integer :: kcell, site(3)

tot = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    site = cellist(kcell)%site
    tot(site(3)) = tot(site(3)) + 1
enddo
tot = tot/nlist
!write(*,'(a,i8,2f6.3)') 'z distribution: ',nlist,tot(NZ/2+10),tot(NZ/2-10)
!write(*,'(10f6.3)') tot
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine check_bindings(str)
character*(*) :: str
integer :: kcell, k, bnd(2), idc
integer :: n26, cells(1000)
integer, allocatable :: nbnd(:)

!call MPI_BARRIER ( MPI_COMM_WORLD, ierr )
write(*,*) 'check_bindings: ',str,' ',istep
allocate(nbnd(NDC))
nbnd = 0
n26 = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    bnd = cellist(kcell)%DCbound
    if (bnd(1) == 0 .and. bnd(2) /= 0) then
        write(*,*) 'Bad DCbound order: ',kcell,bnd
        stop
    endif
    do k = 1,MAX_DC_BIND
        idc = bnd(k)
        if (idc > 0) then
            nbnd(idc) = nbnd(idc) + 1
            if (idc == 26) then
                n26 = n26 + 1
                cells(n26) = kcell
            endif
        endif
    enddo
enddo
deallocate(nbnd)
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine big_check_bindings
integer :: kcell, k, bnd(2), idc
integer, allocatable :: nbnd(:)

allocate(nbnd(NDC))
nbnd = 0
do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    bnd = cellist(kcell)%DCbound
    do k = 1,MAX_DC_BIND
        idc = bnd(k)
        if (idc > 0) then
            nbnd(idc) = nbnd(idc) + 1
        endif
    enddo
enddo
do idc = 1,NDC
    if (.not.DClist(idc)%alive) cycle
    if (nbnd(idc) /= DClist(idc)%nbound) then
        write(*,*) 'big_check_bindings: inconsistent nbnd: ',idc,nbnd(idc),DClist(idc)%nbound
        stop
    endif
enddo
deallocate(nbnd)
end subroutine

!-----------------------------------------------------------------------------------------
! For a location xyz, check that the occupancy() info is consistent with the info in
! the cellist() entries corresponding to occupancy()%indx (if any).  There could be
! 0, 1, or 2 cellist() entries.
!-----------------------------------------------------------------------------------------
subroutine checkslots(msg,xyz)
character*(*) :: msg
integer :: xyz(3)
integer :: indx(2), k, cells, slots, site(3)
logical :: occ(2)

occ = .false.
cells = 0
slots = getslots(xyz)
indx = occupancy(xyz(1),xyz(2),xyz(3))%indx
do k = 1,2
    if (indx(k) > 0) then
        cells = cells + k
        occ(k) = .true.
        site = cellist(indx(k))%site
        if (xyz(1) /= site(1) .or. xyz(2) /= site(2) .or. xyz(3) /= site(3)) then
            write(*,'(a,a,8i6)') msg,' checkslots: site error: ',k,xyz,site
            stop
        endif
    elseif (indx(1) < 0) then
        write(*,*) msg,' checkslots: indx: ',xyz,indx
        stop
    endif
enddo
if (slots /= cells) then
    write(*,'(a,a,6i4,2L2)') msg,' checkslots: mismatch: ',xyz,slots,cells,occ
    stop
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Checks global cell locations for tag validity.
!-----------------------------------------------------------------------------------------
subroutine check_tagged
integer :: k,d2,n,site(3)
type(cell_type) :: cell

do k = 1,nlist
    cell = cellist(k)
    if (cell%tag == TAGGED_CELL) then
        n = n+1
        site = cell%site
        if (.not.taggable(site)) then
            write(*,*) 'Bad tagged cell: ',site,d2,Radius*Radius
        endif
    endif
enddo
write(*,*) 'did check_tagged'
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check(msg,x,y,z)
character*(*) :: msg
integer :: x,y,z,slots,k,indx(2),site(3)
logical :: occ(2)

site = (/x,y,z/)
slots = getslots(site)
occ = .false.
indx = occupancy(x,y,z)%indx
do k = 1,2
    if (cellist(indx(k))%ctype > 0) occ(k) = .true.
enddo
if (slots == 1 .and. (.not.occ(1) .or. occ(2))) then
    write(*,*) 'Bad check: ',msg,x,y,z,slots,occ
    stop
endif
if (slots == 2 .and. (.not.occ(2) .or. occ(1))) then
    write(*,*) 'Bad check: ',msg,x,y,z,slots,occ
    stop
endif
if (slots == 3 .and. .not.(occ(1) .and. occ(2))) then
    write(*,*) 'Bad check: ',msg,x,y,z,slots,occ
    stop
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine checker
integer :: ic,x,y,z,slots,slot,cells,k,indx(2),site(3)
integer, allocatable :: tot(:)
logical :: occ(2)

allocate(tot(NX))
do x = 1,NX
    tot(x) = 0
    do y = 1,NY
        do z = 1,NZ
            site = (/x,y,z/)
            slots = getslots(site)
            cells = 0
            occ = .false.
            indx = occupancy(x,y,z)%indx
            do k = 1,2
                if (indx(k) > 0) then
                    tot(x) = tot(x) + 1
                    cells = cells + k
                    occ(k) = .true.
                    site = cellist(indx(k))%site
                    if (x /= site(1) .or. y /= site(2) .or. z /= site(3)) then
                        write(*,'(a,2i2,2i7,6i4)') 'checker: site error: ',k,indx(k),cellist(indx(k))%ID,x,y,z,site
                        stop
                    endif
                endif
            enddo
            if (slots /= cells) then
                write(*,'(a,6i4,2L2)') 'checker: mismatch: ',x,y,z,slots,cells,occ
                stop
            endif
        enddo
    enddo
enddo

do ic = 1,nlist
    if (cellist(ic)%ID == 0) cycle  ! gap
    site = cellist(ic)%site
    indx = occupancy(site(1),site(2),site(3))%indx
    if (ic == indx(1)) then
        slot = 1
    elseif (ic == indx(2)) then
        slot = 2
    else
        write(*,'(a,7i6)') 'ERROR: checker: bad indx: ',ic,site,indx
        stop
    endif
enddo
deallocate(tot)
write(*,*) 'checked OK: ',' nlist: ',nlist
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine checkcellsite(kcell)
integer :: kcell
integer :: id,site(3)

site = cellist(kcell)%site
id = cellist(kcell)%ID
write(*,'(a,2i8,3i4,4i8)') 'cell: site,indx: ',kcell,id,site,occupancy(site(1),site(2),site(3))%indx
site = cellist(kcell)%site
id = cellist(kcell)%ID
write(*,'(a,2i8,3i4,4i8)') 'big_: site,indx: ',kcell,id,site,occupancy(site(1),site(2),site(3))%indx
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_big_occ
integer :: x,y,z,k,kcell,indx(2)
logical :: OK

OK = .true.
do x = 1,NX
    do y = 1,NY
        do z = 1,NZ
            indx = occupancy(x,y,z)%indx
            do k = 1,2
                kcell = indx(k)
                if (kcell < 0) then
                    if (kcell /= OUTSIDE_TAG .and. -kcell > NDC) then
                        OK = .false.
                        write(*,*) x,y,z,k,kcell
                    endif
                endif
            enddo
        enddo
    enddo
enddo
if (.not.OK) stop
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine checkcell(str)
character*(*) :: str
integer :: ictest = 12

write(*,*) 'checkcell: ',str,'  ',ictest
!write(*,*) cellist(ictest)%ID,cellist(ictest)%PTR1%entrytime
end subroutine

!-----------------------------------------------------------------------------------------
! Checks the two sites (before and after jump) to see if the one free slot is always #2
! HAPPENS ALL THE TIME
!-----------------------------------------------------------------------------------------
subroutine check_site_indx(site1,site2)
integer :: site1(3),site2(3)
integer :: indx(2)

indx = occupancy(site1(1),site1(2),site1(3))%indx
if (indx(1) == 0 .and. indx(2) /= 0) write(*,*) 'check_site_indx: ',site1,indx
indx = occupancy(site2(1),site2(2),site2(3))%indx
if (indx(1) == 0 .and. indx(2) /= 0) write(*,*) 'check_site_indx: ',site1,indx
end subroutine

end module
