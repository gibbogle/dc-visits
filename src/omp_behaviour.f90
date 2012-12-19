! T cell behaviour
module omp_behaviour

use omp_global
use omp_motility

implicit none

!INTEGER,  PARAMETER  ::  DP=SELECTED_REAL_KIND( 12, 60 )


integer, parameter :: POST_DIVISION = SWARMS    ! on division, cells enter SWARMS

integer, parameter :: NORMAL_DIST      = 1
integer, parameter :: LOGNORMAL_DIST   = 2
integer, parameter :: EXPONENTIAL_DIST = 3
integer, parameter :: CONSTANT_DIST    = 4

integer, parameter :: Nbindtime_nc = 10
!real, parameter :: CYTOKINESIS_TIME = 20.0

type(dist_type), allocatable :: stage_dist(:,:,:), life_dist(:), divide_dist(:)

! Non-cognate bindtime distribution
real :: dc_bindtime_nc(10) = (/ 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5 /)
!double apc_contact_prob[] = {6,19,18,14,10.5,7.5,5,3,1.5,1};							! average = 3.4 min
real :: dc_bindprob_nc(10) = (/ 6.0, 26.0, 22.0, 15.0, 9.0, 6.0, 3.5, 2.0, 1.0, 0.5 /)	! average = 3.0 min
real :: dc_cummul_prob_nc(10), dc_initial_cummul_prob_nc(10)

! Mean bindtime for cognate cells at the different stages (CD4 and CD8)
! (Not used for TRANSIENT and CLUSTERS when CONTACT_RULE = CT_HENRICKSON)
!real :: dc_mean_bindtime_c(6,2) = reshape((/ 3.0,60.0,180.0,60.0,11.0,11.0, 3.0,60.0,180.0,60.0,11.0,11.0 /), (/6,2/))
real :: dc_mean_bindtime_c(5,2) = reshape((/ 3.0,10.0,180.0,20.0,0.0, 3.0,10.0,180.0,20.0,0.0 /), (/5,2/))

logical, parameter :: USE_STAGETIME(FINISHED) = &
	(/ .false.,  .true.,    .true.,    .true.,    .true.,   .true. /)
!       NAIVE   TRANSIENT   CLUSTERS   SWARMS   DIVIDING  FINISHED

! Stage times for generation 1 cells (hours)
!real :: mean_stagetime_1(6,2) = reshape((/ 0.0,1.0,12.0,2.0,2.0,2.0, 0.0,1.0,12.0,2.0,2.0,2.0 /), (/6,2/))
real :: mean_stagetime_1(5,2) = reshape((/ 0.0,1.0,12.0,2.0,0.33, 0.0,1.0,12.0,2.0,0.33 /), (/5,2/))
! Stage times for generation n > 1 cells (hours)
!real :: mean_stagetime_n(6,2) = reshape((/ 0.0,1.0,3.0,2.0,2.0,2.0, 0.0,1.0,3.0,2.0,2.0,2.0 /), (/6,2/))
real :: mean_stagetime_n(5,2) = reshape((/ 0.0,1.0,3.0,2.0,0.33, 0.0,1.0,3.0,2.0,0.33 /), (/5,2/))

logical, parameter :: revised_staging = .true.

contains

!--------------------------------------------------------------------------------------
! Compute parameters for probability distributions of lifetime and time-to-divide
! In the case of NORMAL distributions, p1 = mean, p2 = std dev
! For X LOGNORMAL distributed, p1 and p2 are mean and std dev of Y, and X = exp(Y)
! Lifetime and dividetime parameters are given in hours, then converted to minutes
! when used.  Program time is in minutes.
!--------------------------------------------------------------------------------------
subroutine setup_dists
!real :: divide_median1, divide_median2
!real :: divide_shape1, divide_shape2
!real :: life_tc1 = 15, life_tc2 = 18
!real :: life_median1 = 48, life_median2 = 24
!real :: life_median1 = 196, life_median2 = 196      !<------------- Note: hard-coded values
!real :: life_shape1 = 1.5, life_shape2 = 1.4
integer :: i

!life_dist(1)%class = EXPONENTIAL_DIST
!life_dist(i)%p1 = life_tc1
!do i = 1,TC_MAX_GEN
!	life_dist(i)%class = EXPONENTIAL_DIST
!	life_dist(i)%p1 = life_tc2
!enddo

allocate(life_dist(TC_MAX_GEN))        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate(divide_dist(TC_MAX_GEN))      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

life_dist(1)%class = LOGNORMAL_DIST
life_dist(1)%p1 = log(60*TC_life_median1)
life_dist(1)%p2 = log(TC_life_shape)
do i = 2,TC_MAX_GEN
	life_dist(i)%class = LOGNORMAL_DIST
	life_dist(i)%p1 = log(60*TC_life_median2)
	life_dist(i)%p2 = log(TC_life_shape)
enddo

!life_dist(1)%class = NORMAL_DIST
!life_dist(1)%p1 = log(60*life_median1)
!life_dist(1)%p2 = log(life_shape1)
!do i = 2,TC_MAX_GEN
!	life_dist(i)%class = NORMAL_DIST
!	life_dist(i)%p1 = log(60*life_median1)
!	life_dist(i)%p2 = log(life_shape1)
!enddo


!divide_dist(1)%class = LOGNORMAL_DIST
!divide_dist(1)%class = CONSTANT_DIST
!if (divide_dist(1)%class == LOGNORMAL_DIST) then
!	divide_dist(1)%p1 = log(60*divide_median1)
!	divide_dist(1)%p2 = log(divide_shape1)
!	do i = 2,TC_MAX_GEN
!		divide_dist(i)%class = LOGNORMAL_DIST
!		divide_dist(i)%p1 = log(60*divide_median2)
!		divide_dist(i)%p2 = log(divide_shape2)
!	enddo
!elseif (divide_dist(1)%class == CONSTANT_DIST) then
!	do i = 1,TC_MAX_GEN
!		divide_dist(i)%class = CONSTANT_DIST
!		divide_dist(i)%p1 = 60*divide_median2
!		divide_dist(i)%p2 = 0
!	enddo
!endif


divide_dist(1)%class = divide_dist1%class
divide_dist(1)%p1 = divide_dist1%p1
divide_dist(1)%p2 = divide_dist1%p2
divide_dist(2)%class = divide_dist2%class
divide_dist(2)%p1 = divide_dist2%p1
divide_dist(2)%p2 = divide_dist2%p2
do i = 3,TC_MAX_GEN
	divide_dist(i)%class = divide_dist(2)%class
	divide_dist(i)%p1 = divide_dist(2)%p1
	divide_dist(i)%p2 = divide_dist(2)%p2
enddo

!do i = 1,TC_MAX_GEN
!	if (divide_dist(i)%class == NORMAL_DIST) then
!		divide_dist(i)%p1 = 60*divide_dist(i)%p1
!		divide_dist(i)%p2 = 60*divide_dist(i)%p2
!	elseif (divide_dist(i)%class == LOGNORMAL_DIST) then
!		divide_dist(i)%p1 = log(60*divide_dist(i)%p1)
!		divide_dist(i)%p2 = log(divide_dist(i)%p2)
!	elseif (divide_dist(i)%class == CONSTANT_DIST) then
!		divide_dist(i)%p1 = 60*divide_dist(i)%p1
!		divide_dist(i)%p2 = 0
!	endif
!enddo

call make_bind_probs

end subroutine

!--------------------------------------------------------------------------------------
! If a cognate cell is NAIVE it is long-lived - effectively it does not die.
! When a cognate cell has contact with a cognate DC its lifetime is immediately limited.
!--------------------------------------------------------------------------------------
real function TClifetime(ptr)
type (cog_type), pointer :: ptr
integer :: gen, stage, region
real :: p1, p2
integer :: kpar = 0

TClifetime = 0
!stage = get_stage(ptr)
call get_stage(ptr,stage,region)
if (stage == NAIVE) then
    TClifetime = BIG_TIME
    return
endif
gen = get_generation(ptr)
if (gen < 1 .or. gen > TC_MAX_GEN) then
    write(*,*) 'TClifetime: bad gen: ',gen
!    stop
endif
p1 = life_dist(gen)%p1
p2 = life_dist(gen)%p2
select case (life_dist(gen)%class)
case (NORMAL_DIST)
	TClifetime = rv_normal(p1,p2,kpar)
case (LOGNORMAL_DIST)
	TClifetime = rv_lognormal(p1,p2,kpar)
case (CONSTANT_DIST)
	TClifetime = p1
end select
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real function dividetime(gen,ictype)
integer :: gen,ictype
real :: p1, p2
integer :: kpar = 0

dividetime = 0
p1 = divide_dist(gen)%p1
p2 = divide_dist(gen)%p2
!write(*,*) 'dividetime: ',istep,gen,divide_dist(gen)%p1,divide_dist(gen)%p2
select case (divide_dist(gen)%class)
case (NORMAL_DIST)
	dividetime = rv_normal(p1,p2,kpar)
case (LOGNORMAL_DIST)
!    write(*,*) 'dividetime: ',istep,gen,p1,p2
	dividetime = rv_lognormal(p1,p2,kpar)
!	write(*,*) 'dividetime: ',dividetime
case (CONSTANT_DIST)
	dividetime = p1
end select
!if (ictype == CD8_CELL) then
!	dividetime = dividetime*CD8_DIVTIME_FACTOR
!endif
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real function DClifetime(kpar)
integer :: kpar
real :: p1,p2

p1 = log(DC_LIFETIME_MEDIAN*24*60)
p2 = log(DC_LIFETIME_SHAPE)
DClifetime = rv_lognormal(p1,p2,kpar)
write(*,*) 'DClifetime: ',DClifetime 
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real function DCdensity(kpar)
integer :: kpar
real :: p1,p2

p1 = log(DC_ANTIGEN_MEDIAN)
p2 = log(DC_ANTIGEN_SHAPE)
DCdensity = rv_lognormal(p1,p2,kpar)
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real(DP) function rv_normal(p1,p2,kpar)
integer :: kpar
real :: p1,p2
real(DP) :: R

R = par_rnor(kpar)
rv_normal = p1+R*p2
end function

!--------------------------------------------------------------------------------------
! When Y is normal N(p1,p2) then X = exp(Y) is lognormal with
!   median = m = exp(p1)
!   shape  = s = exp(p2)
! Also median = m = mean/(s^2/2)
! kpar = parallel process number
!--------------------------------------------------------------------------------------
real(DP) function rv_lognormal(p1,p2,kpar)
integer :: kpar
real :: p1,p2
real(DP) :: R,z

R = par_rnor(kpar)
z = p1 + R*p2
rv_lognormal = exp(z)
end function

!--------------------------------------------------------------------------------------
! For testing.
!--------------------------------------------------------------------------------------
real(DP) function my_rnor()
real(DP) :: sum, R
integer :: k
integer :: kpar=0

sum = 0
do k = 1,12
!    call random_number(R)
    R = par_uni(kpar)
    sum = sum + R
enddo
my_rnor = sum - 6.0
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real(DP) function rv_exponential(p1)
real :: p1
real(DP) :: r
integer :: kpar = 0

r = par_rexp(kpar)
rv_exponential = p1*r
end function

!--------------------------------------------------------------------------------------
! Cumulative probability distribution for a lognormal variate with median m, shape s
! Computes Pr(X < a) where X = exp(Y) and Y is N(p1,p2), p1 = log(m), p2 = log(s)
! Pr(X < a) = Pr(Y < log(a)) = Pr(p1 + p2*R < log(a)) = Pr(R < (log(a)-p1)/p2)
! where R is N(0,1)
!--------------------------------------------------------------------------------------
real(DP) function cum_prob_lognormal(a,p1,p2)
real :: a, p1, p2
real(DP) :: b, prob

!p1 = log(m)
!p2 = log(s)
b = (log(a) - p1)/p2
prob = 0.5 + 0.5*erf(b/sqrt(2.0))
cum_prob_lognormal = prob
end function

!--------------------------------------------------------------------------------
! The probability histogram for a (supposedly) lognormal variable X is given
! by the probability values p(i), i=1,..,n and the associated interval
! delimiting values x(i), i=1,..,n+1
! where p(i) = Pr{x(i) <= X < x(i+1)}
! The cumulative probability for a lognormal variate with median m, shape s is
! Pr(X < a) where X = exp(Y) and Y is N(p1,p2), p1 = log(m), p2 = log(s)
! Pr(X < a) = Pr(Y < log(a)) = Pr(p1 + p2*R < log(a)) = Pr(R < (log(a)-p1)/p2)
! where R is N(0,1)
!--------------------------------------------------------------------------------
subroutine lognfit(n,x,p,mrange,srange,mean,shape)
integer :: n
real :: x(*), p(*), mrange(2), srange(2)
real :: mean, shape
integer :: i
real :: alo, ahi, p1, p2, c, err

p1 = log(mrange(1))
p2 = log(srange(1))
err = 0
do i = 1,n
	alo = x(i)
	ahi = x(i+1)
	c = cum_prob_lognormal(ahi,p1,p2) - cum_prob_lognormal(alo,p1,p2)
	err = err + (c - p(i))**2
enddo
end subroutine

!--------------------------------------------------------------------------------------
! Note: _nc refers to non-cognate T cells
!--------------------------------------------------------------------------------------
subroutine make_bind_probs
integer :: i
real :: psum, ptsum
real, allocatable :: pt(:)

allocate(pt(Nbindtime_nc))
psum = 0
ptsum = 0
do i = 1,Nbindtime_nc
	psum = psum + dc_bindprob_nc(i)
	pt(i) = dc_bindprob_nc(i)*dc_bindtime_nc(i)
	ptsum = ptsum + pt(i)
enddo
do i = 1,Nbindtime_nc
	dc_bindprob_nc(i) = dc_bindprob_nc(i)/psum
	pt(i) = pt(i)/ptsum
enddo
psum = 0
ptsum = 0
do i = 1,Nbindtime_nc
	psum = psum + dc_bindprob_nc(i)
	dc_cummul_prob_nc(i) = psum
	ptsum = ptsum + pt(i)
	dc_initial_cummul_prob_nc(i) = ptsum
enddo
!write(*,*) 'dc_cummul_prob_nc:'
!write(*,'(10f5.2)') dc_cummul_prob_nc
!write(*,*) 'dc_initial_cummul_prob_nc:'
!write(*,'(10f5.2)') dc_initial_cummul_prob_nc 
!write(*,*)

deallocate(pt)

end subroutine

!----------------------------------------------------------------------------------------
! These parameters are not yet in an input file, but they will be ...
!----------------------------------------------------------------------------------------
subroutine old_unread_cell_params

!NX = 100

! T cell parameters
!use_traffic = .false.
!TCR_splitting = .false.  ! enable sharing of integrated TCR signal between progeny cells
!transient_stagetime = 1.0
!clusters_stagetime = 13.0
!transient_bindtime = 10.0
!clusters_bindtime = 180.0
!swarms_bindtime = 20.0
!TC_life_median1 = 196				! median lifetime of naive T cells
!TC_life_median2 = 196				! median lifetime of activated T cells
!TC_life_shape= 1.4					! shape parameter for lifetime of T cells

! DC parameters
!use_DCflux = .false.
!DC_pMHC_THRESHOLD = 10   ! DC density limit for TCR stimulation capability (0.05)
!DC_STIM_THRESHOLD = 300
!CONTACT_RULE = CT_HENRICKSON
!ABIND1 = 0.4    ! binding to a DC
!ABIND2 = 0.8

! Vascularity parameters
!Inflammation_days1  ! Days of plateau level - parameters for VEGF_MODEL = 1
!Inflammation_days2  ! End of inflammation
!Inflammation_level	! This is the level of inflammation

! Egress parameters
!chemo_radius = 5  ! radius of chemotactic influence (sites)
!exit_fraction = 1.0/1000.    ! number of exits as a fraction of T cell population

end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine read_fixed_params(ok)
logical :: ok
logical :: ext
real :: TC_CD8_FRACTIONzzz
logical :: use_DCvisits_params
character*(128) :: DCparamfile

inquire(file=fixedfile,exist=ext)
if (.not.ext) then
	call logger("Fixed parameter input file not found")
	ok = .false.
	return
endif
open(nfcell,file=fixedfile,status='old',err=99)
ok = .true.
! T cell parameters
!read(nfcell,*) use_traffic			!= .false.
read(nfcell,*) TCR_splitting		!= .false.  ! enable sharing of integrated TCR signal between progeny cells
read(nfcell,*) transient_stagetime	!= 1.0
read(nfcell,*) clusters_stagetime	!= 13.0
read(nfcell,*) transient_bindtime	!= 10.0
read(nfcell,*) clusters_bindtime	!= 180.0
read(nfcell,*) swarms_bindtime		!= 20.0
read(nfcell,*) TC_life_median1		!= 196				! median lifetime of naive T cells
read(nfcell,*) TC_life_median2		!= 196				! median lifetime of activated T cells
read(nfcell,*) TC_life_shape		!= 1.4				! shape parameter for lifetime of T cells
read(nfcell,*) NTC_LN				!= 3.0e07			! number of T cells in a LN
read(nfcell,*) NTC_BODY				!= 1.6e09			! number of circulating T cells in the whole body
read(nfcell,*) NLN_RESPONSE								! number of LNs in the response
!write(*,*) 'NTC_LN, NTC_BODY: ',NTC_LN, NTC_BODY

! DC parameters
read(nfcell,*) use_DCflux			!= .false.
read(nfcell,*) DC_pMHC_THRESHOLD	!= 10   ! DC density limit for TCR stimulation capability (0.05)
read(nfcell,*) DC_STIM_THRESHOLD	!= 300
read(nfcell,*) N_STIM
read(nfcell,*) CONTACT_RULE			!= CT_HENRICKSON
read(nfcell,*) ABIND1				!= 0.4    ! binding to a DC
read(nfcell,*) ABIND2				!= 0.8

! Egress parameters
read(nfcell,*) exit_fraction		!= 1.0/1000.    ! number of exits as a fraction of T cell population

!CD25/IL2
read(nfcell,*) optionA                      ! 1,2
read(nfcell,*) optionB                      ! 1,2,3
read(nfcell,*) optionC                      ! 1,2
read(nfcell,*) IL2_PRODUCTION_TIME			! duration of IL2/CD25 production
read(nfcell,*) CD25_DIVISION_THRESHOLD	    ! CD25 store level needed for division of activated cell (optionA = 1)
read(nfcell,*) CD25_SURVIVAL_THRESHOLD		! CD25 store level needed for survival of activated cell (optionC = 1)

read(nfcell,*) TC_RADIUS					! radius of T cell (um)
read(nfcell,*) TC_STIM_WEIGHT				! contribution of stimulation to act level
read(nfcell,*) TC_MAX_GEN				    ! maximum T cell generation
read(nfcell,*) TC_CD8_FRACTIONzzz			! fraction of all T cells that are CD8 (now in the GUI)
read(nfcell,*) CD8_DIVTIME_FACTOR			! CD8 divide time as multiple of CD4 divide time
read(nfcell,*) DC_MULTIBIND_PROB			! reducing factor to bind prob for each current DC binding
read(nfcell,*) MAX_DC_BIND					! number of DC that a cell can bind to simultaneously
read(nfcell,*) divide_dist1%class
read(nfcell,*) divide_dist2%class
read(nfcell,*) GAMMA						! controls crowding

read(nfcell,*) DC_ACTIV_TAPER				! time (hours) over which DC activity decays to zero
read(nfcell,*) DC_DENS_BY_STIM              ! rate of reduction of density by TCR stimulation
read(nfcell,*) DC_FACTOR                    ! multiplying factor for DC number (initial and influx)

read(nfcell,*) efactor                      ! If constant_efactor = true, this is the factor for the p correction
read(nfcell,*) VEGF_MODEL                   ! 1 = VEGF signal from inflammation, 2 = VEGF signal from DCactivity

read(nfcell,*) fix_avidity                  ! true if avidity takes discrete values, false if a distribution is used
read(nfcell,*) avidity_logscale             ! true if actual avidity = 10^(avidity_min + i*avidity_step)
read(nfcell,*) avidity_nlevels              ! If fix_avidity, number of discrete avidity values
read(nfcell,*) avidity_min                  ! minimum value
read(nfcell,*) avidity_step                 ! step between equi-spaced values

read(nfcell,*) K1_S1P1
read(nfcell,*) K2_S1P1
read(nfcell,*) K1_CD69
read(nfcell,*) K2_CD69
! Parameters added to control HEV and DC placement
read(nfcell,*) R_HEV_min
read(nfcell,*) R_HEV_max
read(nfcell,*) R_DC_min
read(nfcell,*) R_DC_max
! Parameters added for DCvisits simulations
!read(nfcell,'(L)') use_DCvisits_params
!if (use_DCvisits_params) then
!	track_DCvisits = .true.
!	read(nfcell,'(a)') DCparamfile
!	close(nfcell)
!	open(nfcell,file=DCparamfile,status='old')
!	read(nfcell,'(L)') TAGGED_DC_CHEMOTAXIS		! DC visits are logged for tagged T cells only
!	read(nfcell,*) istep_DCvisits	! 24*60*4	! delay before counting visits (when use_traffic) (steps)
!	read(nfcell,*) t_log_DCvisits	! 2*24*60	! time to log DC visits of retained cells (mins)
!	read(nfcell,'(L)') USE_DC_COGNATE 		! if only a fraction (0.5) of DCs bear cognate antigen, secrete chemokine
!	read(nfcell,*) DC_COGNATE_FRACTION		! fraction of DCs that bear cognate antigen
!	read(nfcell,*) RETAIN_TAGLIMIT_FRACTION ! 0.5 (when use_traffic)
!	read(nfcell,'(L)') DC_CHEMO_NOTRAFFIC   ! cells are tagged initially, no use_traffic
!	read(nfcell,*) DC_CHEMO_FRACTION 		! fraction of cells tagged when DC_CHEMO_NOTRAFFIC
!	read(nfcell,*) HI_CHEMO_FRACTION 		! fraction of tagged cells with full chemotactic susceptibility
!	read(nfcell,*) HI_CHEMO 				! chemotactic susceptibility of the HI_CHEMO_FRACTION cells
!	read(nfcell,*) LO_CHEMO 				! chemotactic susceptibility of the other tagged cells
!endif
MAX_DC_BIND = 1
close(nfcell)
return
99	continue
ok = .false.
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine read_cell_params(ok)
logical :: ok
real :: sigma, divide_mean1, divide_shape1, divide_mean2, divide_shape2, real_DCradius
integer :: i, invitro, shownoncog, ncpu_dummy, dcsinjected, ispecial
integer :: usetraffic, useexitchemo, useDCchemo, computedoutflow, useCCL3secretion
character(64) :: specialfile
character(4) :: logstr
logical, parameter :: use_chemo = .false.

ok = .false.
!write(*,*) 'Read cell parameter file: ',inputfile
!call logger('Reading cell parameter file')
!call logger(inputfile)
open(nfcell,file=inputfile,status='old')
read(nfcell,*) TC_AVIDITY_MEAN              ! mean of avidity distribution (only if fix_avidity = false)
read(nfcell,*) TC_AVIDITY_SHAPE			    ! shape -> 1 gives normal dist with small variance
read(nfcell,*) TC_CD8_FRACTION				! fraction of all T cells that are CD8
read(nfcell,*) TC_COGNATE_FRACTION(1)		! fraction of CD4 T cells that are cognate
read(nfcell,*) TC_COGNATE_FRACTION(2)		! fraction of CD8 T cells that are cognate
!read(nfcell,*) TC_CONVERSION_TIME			! time (in hours) over which T cells become cognate
read(nfcell,*) TC_STIM_RATE_CONSTANT		! rate const for TCR stimulation (-> molecules/min)
!read(nfcell,*) TC_STIM_WEIGHT				! contribution of stimulation to act level
read(nfcell,*) TC_STIM_HALFLIFE				! halflife of T cell stimulation (hours)
!read(nfcell,*) TC_MAX_GEN				    ! maximum T cell generation
!read(nfcell,*) divide_dist1%p1
!read(nfcell,*) divide_dist1%p2
!read(nfcell,*) divide_dist2%p1
!read(nfcell,*) divide_dist2%p2
read(nfcell,*) divide_mean1
read(nfcell,*) divide_shape1
read(nfcell,*) divide_mean2
read(nfcell,*) divide_shape2
read(nfcell,*) BETA							! speed: 0 < beta < 1		(0.65)
read(nfcell,*) RHO							! persistence: 0 < rho < 1	(0.95)

read(nfcell,*) DC_ANTIGEN_MEAN				! mean DC antigen density
read(nfcell,*) DC_ANTIGEN_SHAPE				! DC antigen density shape param
read(nfcell,*) DC_LIFETIME_MEAN				! days
read(nfcell,*) DC_LIFETIME_SHAPE 			! days
!read(nfcell,*) DC_ACTIV_TAPER				! time (hours) over which DC activity decays to zero
read(nfcell,*) DC_BIND_DELAY				! delay after unbinding before next binding (min)
!read(nfcell,*) DC_BIND_ALFA					! binding prob parameter
!read(nfcell,*) DC_MULTIBIND_PROB			! reducing factor to bind prob for each current DC binding
!read(nfcell,*) DC_DENS_BY_STIM              ! rate of reduction of density by TCR stimulation
read(nfcell,*) DC_DENS_HALFLIFE             ! base half-life of DC activity (hours)
read(nfcell,*) MAX_TC_BIND						! number of T cells that can bind to a DC
read(nfcell,*) MAX_COG_BIND					! number of cognate T cells that can bind to a DC simultaneously

!read(nfcell,*) optionA                      ! 1,2
!read(nfcell,*) optionB                      ! 1,2,3
!read(nfcell,*) optionC                      ! 1,2
!read(nfcell,*) IL2_PRODUCTION_TIME			! duration of IL2/CD25 production

read(nfcell,*) IL2_THRESHOLD			    ! stimulation needed to initiate IL-2/CD25 production
read(nfcell,*) ACTIVATION_THRESHOLD			! stimulation needed for activation
read(nfcell,*) FIRST_DIVISION_THRESHOLD(1)	! activation level needed for first division
read(nfcell,*) DIVISION_THRESHOLD(1)		! activation level needed for subsequent division
read(nfcell,*) EXIT_THRESHOLD(1)			! activation level below which exit is permitted
read(nfcell,*) STIMULATION_LIMIT			! maximum activation level
!read(nfcell,*) CD25_DIVISION_THRESHOLD	    ! CD25 store level needed for division of activated cell (optionA = 1)
!read(nfcell,*) CD25_SURVIVAL_THRESHOLD		! CD25 store level needed for survival of activated cell (optionC = 1)

!read(nfcell,*) divide_dist1%class
!read(nfcell,*) divide_dist2%class
!read(nfcell,*) CD8_DIVTIME_FACTOR			! CD8 divide time as multiple of CD4 divide time

read(nfcell,*) NX							! rule of thumb: about 4*BLOB_RADIUS
read(nfcell,*) BLOB_RADIUS					! initial T cell blob size (sites)
read(nfcell,*) TC_FRACTION					! T cell fraction of paracortex volume
!read(nfcell,*) TC_RADIUS						! radius of T cell (um)
read(nfcell,*) FLUID_FRACTION				! fraction of paracortex that is fluid
read(nfcell,*) real_DCradius						! radius of DC sphere of influence (um)

!read(nfcell,*) GAMMA						! controls crowding

read(nfcell,*) TC_TO_DC						! number of T cells for every DC
!read(nfcell,*) DC_FACTOR                    ! multiplying factor for DC number (initial and influx)
!read(nfcell,*) MAX_DC_BIND					! number of DC that a cell can bind to simultaneously
read(nfcell,*) DCrate_100k                  ! DC influx rate corresponding to NTcells0 = 100k
read(nfcell,*) T_DC1                        ! Duration of constant DC influx (hours)
read(nfcell,*) T_DC2                        ! Time of cessation of DC influx (hours)
read(nfcell,*) DCsinjected					! select DC injection into experimental animals
read(nfcell,*) DCinjectionfile				! file with schedule of DC injection
!read(nfcell,*) T_DC_INJECTION				! Time of DC injection
!read(nfcell,*) DC_FRACTION					! Fraction of DCs that are bearing antigen

read(nfcell,*) usetraffic					! use T cell trafficking
read(nfcell,*) useexitchemo                 ! use exit chemotaxis
read(nfcell,*) useDCchemo					! use DC chemotaxis
read(nfcell,*) computedoutflow				! compute outflow (with inflow)
read(nfcell,*) RESIDENCE_TIME(CD4)          ! CD4 T cell residence time in hours -> inflow rate
read(nfcell,*) RESIDENCE_TIME(CD8)          ! CD8 T cell residence time in hours -> inflow rate
! Vascularity parameters
read(nfcell,*) Inflammation_days1	        ! Days of plateau level - parameters for VEGF_MODEL = 1
read(nfcell,*) Inflammation_days2	        ! End of inflammation
read(nfcell,*) Inflammation_level	        ! This is the level of inflammation
read(nfcell,*) exit_region                  ! blob region for cell exits (1 = all z, 2 = lower half, 3 = blob portals, 4 = surface portals)
!read(nfcell,*) efactor                     ! If constant_efactor = true, this is the factor for the p correction
!read(nfcell,*) VEGF_MODEL                  ! 1 = VEGF signal from inflammation, 2 = VEGF signal from DCactivity
read(nfcell,*) chemo_radius			        ! radius of chemotactic influence (um)
read(nfcell,*) chemo_K_exit                 ! level of chemotactic influence towards exits
read(nfcell,*) chemo_K_DC                  ! level of chemotactic influence towards DCs

!!read(nfcell,*) receptor(CCR1)%strength      ! relative strength of CCL3-CCR1 chemotactic influence (chemo_K_DC)
    read(nfcell,*) receptor(CCR1)%level(1)
    read(nfcell,*) receptor(CCR1)%level(2)
    read(nfcell,*) receptor(CCR1)%level(3)
    read(nfcell,*) receptor(CCR1)%saturation_threshold	! CCL3 signal level to saturate CCR1 receptor
    read(nfcell,*) receptor(CCR1)%refractory_time	! Time for CCR1 receptor to recover sensitivity after desensitization
    read(nfcell,*) chemo(CCL3)%diff_coef
    read(nfcell,*) chemo(CCL3)%halflife
    read(nfcell,*) useCCL3secretion
    read(nfcell,*) chemo(CCL3)%bdry_rate
    read(nfcell,*) chemo(CCL3)%bdry_conc
    read(nfcell,*) chemo(CCL3)%radius		! Sites within this radius of the DC receive CCL3 concentration (+ all DC sites)

!read(nfcell,*) fix_avidity                  ! true if avidity takes discrete values, false if a distribution is used
!read(nfcell,*) avidity_logscale             ! true if actual avidity = 10^(avidity_min + i*avidity_step)
!read(nfcell,*) avidity_nlevels              ! If fix_avidity, number of discrete avidity values
!read(nfcell,*) avidity_min                  ! minimum value
!read(nfcell,*) avidity_step                 ! step between equi-spaced values

read(nfcell,*) days							! number of days to simulate
read(nfcell,*) seed(1)						! seed vector(1) for the RNGs
read(nfcell,*) seed(2)						! seed vector(2) for the RNGs
read(nfcell,*) ncpu_dummy					! just a placeholder for ncpu, not used currently
read(nfcell,*) NT_GUI_OUT					! interval between GUI outputs (timesteps)
read(nfcell,*) SPECIES						! animal species source of T cells
read(nfcell,*) invitro 						! select in vivo or in vitro simulation
read(nfcell,*) IV_WELL_DIAMETER				! diameter of in vitro well (mm)
read(nfcell,*) IV_NTCELLS					! initial T cell population in vitro
read(nfcell,*) IV_COGNATE_FRACTION			! fraction of in vitro cells that are cognate for DC antigen
read(nfcell,*) shownoncog			        ! display non-cognate T cells
read(nfcell,*) ispecial						! special case
read(nfcell,*) specialfile					! special case input data file
read(nfcell,*) fixedfile					! file with "fixed" parameter values
close(nfcell)

!call logger('Finished reading cell parameter file')

DC_RADIUS = real_DCradius                       ! convert real -> integer
DCrate_100k = DCrate_100k*Inflammation_level	! DC influx is scaled by the inflammation level
T_DC1 = T_DC1*24								! convert days -> hours
T_DC2 = T_DC2*24
CTYPE_FRACTION(CD8) = TC_CD8_FRACTION
CTYPE_FRACTION(CD4) = 1 - CTYPE_FRACTION(CD8)
ave_residence_time = CTYPE_FRACTION(CD4)*residence_time(CD4) + CTYPE_FRACTION(CD8)*residence_time(CD8)
if (useCCL3secretion) then
	chemo(CCL3)%use_secretion = .true.
else
	chemo(CCL3)%use_secretion = .false.
endif
if (shownoncog == 0) then
    IV_SHOW_NONCOGNATE = .false.
else
    IV_SHOW_NONCOGNATE = .true.
endif
if (invitro == 1) then
	IN_VITRO = .true.
else
	IN_VITRO = .false.
	IV_SHOW_NONCOGNATE = .false.
endif
if (DCsinjected == 1) then
	DC_INJECTION = .true.
else
	DC_INJECTION = .false.
endif
if (usetraffic == 1) then
	USE_TRAFFIC = .true.
else
	USE_TRAFFIC = .false.
endif
if (useexitchemo == 1) then
	use_exit_chemotaxis = .true.
else
	use_exit_chemotaxis = .false.
	chemo_K_exit = 0
endif

receptor(CCR1)%strength = chemo_K_DC
if (useDCchemo == 1 .and. receptor(CCR1)%strength > 0) then
	use_DC_chemotaxis = .true.
endif
if (ispecial /= 0) then
	call read_specialcase(ispecial,specialfile,ok)
	if (.not.ok) return
endif

if (use_DC_chemotaxis .and. receptor(CCR1)%strength > 0) then
	receptor(CCR1)%used = .true.
	chemo(CCL3)%used = .true.
	chemo_K_DC = receptor(CCR1)%strength
	DC_CHEMO_DELAY = DC_BIND_DELAY
else
	use_DC_chemotaxis = .false.
	receptor(CCR1)%used = .false.
	chemo(CCL3)%used = .false.
	chemo_K_DC = 0
endif
if (useCCL3secretion) then
	chemo(CCL3)%use_secretion = .true.
else
	chemo(CCL3)%use_secretion = .false.
endif

if (computedoutflow == 1) then
	computed_outflow = .true.
else
	computed_outflow = .false.
endif
if (track_DCvisits) then
	TC_COGNATE_FRACTION = 0
endif
BLOB_PORTALS = .false.
SURFACE_PORTALS = .false.
if (exit_region == EXIT_BLOB_PORTALS) then 
	use_portal_egress = .true.
	BLOB_PORTALS = .true.
elseif (exit_region == EXIT_SURFACE_PORTALS) then
	use_portal_egress = .true.
	SURFACE_PORTALS = .true.
	call logger('use_portal_egress')
else
	use_portal_egress = .false.
endif

if ((exit_region == EXIT_EVERYWHERE .or. exit_region == EXIT_LOWERHALF) .and.  &
	use_exit_chemotaxis) then
	write(logmsg,*) 'Error: chemotaxis requires portal egress'
	call logger(logmsg)
	ok = .false.
	return
endif
if (DC_INJECTION .and. .not.L_selectin) then
	write(logmsg,*) 'Error: DC_INJECTION requires L_selectin'
	call logger(logmsg)
	ok = .false.
	return
endif
call read_fixed_params(ok)
if (.not.ok) then
	write(logmsg,'(a,a)') 'Error reading fixed input data file: ',fixedfile
	call logger(logmsg)
	return
endif
if (DCrate_100k == 0) then
	use_DCflux = .false.
endif

mean_stagetime_1(2,:) = transient_stagetime		! Note that these times are not needed for gen>1,
mean_stagetime_1(3,:) = clusters_stagetime		! because the stages are not revisited
dc_mean_bindtime_c(2,:) = transient_bindtime	! Not used for CT_HENRICKSON
dc_mean_bindtime_c(3,:) = clusters_bindtime		! Not used for CT_HENRICKSON
dc_mean_bindtime_c(4,:) = swarms_bindtime		! Always used

if (mod(NX,2) /= 0) NX = NX+1					! ensure that NX is even
if (TC_MAX_GEN > MMAX_GEN) then
    write(logmsg,*) 'TC_MAX_GEN > MMAX_GEN: ',TC_MAX_GEN,MMAX_GEN
    call logger(logmsg)
    stop
endif
PI = 4*atan(1.0)
DELTA_X = (4*PI/(3*TC_FRACTION))**0.33333*TC_RADIUS

NDCsites = NDCcore
! Override parameter settings for IN_VITRO case
if (IN_VITRO) then
	! Note that now BLOB_RADIUS is actually petri dish radius
	DELTA_X = 1.5*TC_RADIUS		! needs to be adjusted
	BLOB_RADIUS = (1000*IV_WELL_DIAMETER/2)/DELTA_X		! mm -> sites
	NDCsites = NDCcore_IN_VITRO
	use_DCflux = .false.
	use_traffic = .false.
	exit_region = EXIT_NONE
!	TC_COGNATE_FRACTION = IV_COGNATE_FRACTION
	call logger("In vitro simulation  ")
endif

DC_RADIUS = DC_RADIUS/DELTA_X				    ! convert from um to lattice grids
chemo_radius = chemo_radius/DELTA_X				! convert from um to lattice grids
!write(*,*) 'DC_RADIUS, chemo_radius: ',DC_RADIUS,chemo_radius
if (chemo_radius < DC_RADIUS) then
	call logger('Currently it is required that chemo_radius >= DC_RADIUS')
	stop
endif
chemo_N = max(3,int(chemo_radius + 0.5))	! convert from um to lattice grids
!write(*,*) 'chemo_N: ',chemo_N
! Note that currently exit chemotaxis and DC chemotaxis are treated in the same way - same decay
chemo_exp = log(1/CHEMO_MIN)/log(chemo_radius)
FIRST_DIVISION_THRESHOLD(2) = FIRST_DIVISION_THRESHOLD(1)
DIVISION_THRESHOLD(2) = DIVISION_THRESHOLD(1)
EXIT_THRESHOLD(2) = EXIT_THRESHOLD(1)
sigma = log(TC_AVIDITY_SHAPE)
TC_AVIDITY_MEDIAN = TC_AVIDITY_MEAN/exp(sigma*sigma/2)
sigma = log(DC_ANTIGEN_SHAPE)
DC_ANTIGEN_MEDIAN = DC_ANTIGEN_MEAN/exp(sigma*sigma/2)
sigma = log(DC_LIFETIME_SHAPE)
DC_LIFETIME_MEDIAN = DC_LIFETIME_MEAN/exp(sigma*sigma/2)

sigma = log(divide_shape1)
divide_dist1%p1 = log(60*divide_mean1/exp(sigma*sigma/2))
divide_dist1%p2 = sigma
sigma = log(divide_shape2)
divide_dist2%p1 = log(60*divide_mean2/exp(sigma*sigma/2))
divide_dist2%p2 = sigma
!write(*,*) 'divide_dist2: ',divide_dist2

if (TC_TO_DC > 0 .or. DC_INJECTION) then
    use_DC = .true.
else
    use_DC = .false.
endif

if ((TC_COGNATE_FRACTION(1) == 0 .and. TC_COGNATE_FRACTION(2) == 0) .or. .not.use_DC) then
    use_cognate = .false.
else
    use_cognate = .true.
endif

!if (exit_region == EXIT_CHEMOTAXIS) then
!    use_chemotaxis = .true.
!    if (exit_rule /= 3) then
!        call logger('EXIT_CHEMOTAXIS needs exit_rule = 3')
!        stop
!    endif
!else
!    use_chemotaxis = .false.
!endif

if (chemo_K_exit == 0.0) then
    ep_factor = 2.4
elseif (chemo_K_exit <= 0.2) then
    ep_factor = 2.3
elseif (chemo_K_exit <= 0.4) then
    ep_factor = 2.2
elseif (chemo_K_exit <= 0.6) then
    ep_factor = 2.15
elseif (chemo_K_exit <= 0.8) then
    ep_factor = 2.1
else
    ep_factor = 2.1
endif

if (DC_DENS_HALFLIFE > 0) then
!    DCdecayrate = log(2.0)/(DC_DENS_HALFLIFE*60)    ! rate/min
	DCdecayrate = DecayRate(DC_DENS_HALFLIFE)
else
    DCdecayrate = 0
endif

if (TC_STIM_HALFLIFE > 0) then
!    TCRdecayrate = log(2.0)/(TC_STIM_HALFLIFE*60)    ! rate/min
    TCRdecayrate = DecayRate(TC_STIM_HALFLIFE)    ! rate/min
else
    TCRdecayrate = 0
endif
if (.not.use_cytokines) then
    TC_STIM_WEIGHT = 1
endif

if (VEGF_MODEL == 0) then
    vary_vascularity = .false.
else
    vary_vascularity = .true.
endif

if (.not.vary_vascularity .and. .not.use_cognate) then
    steadystate = .true.
else
    steadystate = .false.
endif

!if (exit_rule == 3) then
!    CD69_threshold = EXIT_THRESHOLD(1)
!endif


Ve = DELTA_X*DELTA_X*DELTA_X
Vc = FLUID_FRACTION*Ve

Nsteps = days*60*24/DELTA_T
call make_outputfilename
open(nfout,file=outputfile,status='replace')
write(nfout,*) 'TC_COGNATE_FRACTION: ',TC_COGNATE_FRACTION
write(logmsg,*) 'Opened nfout: ',outputfile
call logger(logmsg)
!if (save_input) then
!    call save_inputfile(inputfile)
!    call save_parameters
!	call save_inputfile(fixedfile)
!endif

call setup_dists
!write(*,*) 'open resultfile: ',resultfile
!open(nfres,file=resultfile,status='replace')
!write(nfres,'(i3,a)') ntres,' ntres'
!write(nfres,'(i2,a)') TC_MAX_GEN, ' TC_MAX_GEN'
!write(nfres,'(i2,f6.0,a)') TCR_nlevels, TCR_limit, ' TCR_nlevels,TCR_limit'
if (avidity_logscale) then
    logstr = ' log'
else
    logstr = ' lin'
endif
if (fix_avidity) then
!    write(nfres,'(i2,2f10.6,a,a)') avidity_nlevels,avidity_min,avidity_step,logstr,'  avidity_nlevels,avidity_min,avidity_step'
    do i = 1,avidity_nlevels
        if (avidity_logscale) then
            avidity_level(i) = 10**(avidity_min + (i-1)*avidity_step)
        else
            avidity_level(i) = avidity_min + (i-1)*avidity_step
        endif
    enddo
!    write(*,*) 'avidity_level'
!    write(*,'(10f8.4)') avidity_level
else
!    write(nfres,'(i2,2f6.2,a,a)') 0,0.,0.,logstr,'  avid_nlevels,avid_min,avid_step'
endif

if (save_DCbinding) then
    open(nfdcbind,file='DCbinding.out',status='replace')
    dcbind = 0
endif
ok = .true.
end subroutine

!----------------------------------------------------------------------------------------
! Special cases:
! 1 
! To quantify DC visits for an artificial situation of no trafficking.
! The T cell population does not change, a small fraction of T cells are tagged.  These 
! tagged cells are assigned either HI or LO chemotactic susceptibility.  There is no TCR
! stimulation, no T cell activation.  DC visits are counted.
! 2
! To quantify DC visits with trafficking.
! 3
! DC chemokine secretion is linked to contact with cognate CD4  cells.
!----------------------------------------------------------------------------------------
subroutine read_specialcase(icase,casefile,ok)
integer :: icase
character*(*) :: casefile
logical :: ok

!if (use_DCvisits_params) then
ok = .true.
write(logmsg,*) 'read_specialcase: icase, casefile: ',icase,'  ',casefile
call logger(logmsg)
if (icase == 1) then
	if (receptor(CCR1)%strength == 0) then
		call logger('Error: need to specify chemo_K_DC for special case: TAGGED_DC_CHEMOTAXIS')
		ok = .false.
		return
	endif
	track_DCvisits = .true.
	use_traffic = .false.
	open(nfcell,file=casefile,status='old')
	read(nfcell,'(L)') USE_DC_CHEMOTAXIS		! DC chemotaxis
	read(nfcell,'(L)') TAGGED_DC_CHEMOTAXIS		! DC visits are logged for tagged T cells only
	read(nfcell,*) istep_DCvisits	! 24*60*4	! delay before counting visits (when use_traffic) (steps)
	read(nfcell,*) t_log_DCvisits	! 2*24*60	! time to log DC visits of retained cells (mins)
	read(nfcell,'(L)') USE_DC_COGNATE 		! if only a fraction (0.5) of DCs bear cognate antigen, secrete chemokine
	read(nfcell,*) DC_COGNATE_FRACTION		! fraction of DCs that bear cognate antigen
	read(nfcell,*) RETAIN_TAGGED_CELLS		! applies with use_traffic
	read(nfcell,'(L)') DC_CHEMO_NOTRAFFIC   ! cells are tagged initially, no use_traffic
	read(nfcell,*) DC_CHEMO_FRACTION 		! fraction of cells tagged when DC_CHEMO_NOTRAFFIC
	read(nfcell,*) HI_CHEMO_FRACTION 		! fraction of tagged cells with full chemotactic susceptibility
	read(nfcell,*) HI_CHEMO 				! chemotactic susceptibility of the HI_CHEMO_FRACTION cells
	read(nfcell,*) LO_CHEMO 				! chemotactic susceptibility of the other tagged cells
	close(nfcell)
elseif (icase == 2) then
	track_DCvisits = .true.
	use_traffic = .true.
	open(nfcell,file=casefile,status='old')
	read(nfcell,'(L)') USE_DC_CHEMOTAXIS		! DC chemotaxis
	read(nfcell,'(L)') TAGGED_DC_CHEMOTAXIS		! DC visits are logged for tagged T cells only
	read(nfcell,*) istep_DCvisits	! 24*60*4	! delay before counting visits (when use_traffic) (steps)
	read(nfcell,*) t_log_DCvisits	! 2*24*60	! time to log DC visits of retained cells (mins)
	read(nfcell,'(L)') USE_DC_COGNATE 		! if only a fraction (0.5) of DCs bear cognate antigen, secrete chemokine
	read(nfcell,*) DC_COGNATE_FRACTION		! fraction of DCs that bear cognate antigen
	read(nfcell,*) RETAIN_TAGGED_CELLS		! applies with use_traffic
	read(nfcell,*) RETAIN_TAGLIMIT_FRACTION ! 0.5 (when use_traffic)
	read(nfcell,'(L)') DC_CHEMO_NOTRAFFIC   ! cells are tagged initially, no use_traffic
	read(nfcell,*) TAGGED_CHEMO_FRACTION 		! fraction of cells tagged when DC_CHEMO_NOTRAFFIC
	read(nfcell,*) HI_CHEMO_FRACTION 		! fraction of tagged cells with full chemotactic susceptibility
	read(nfcell,*) HI_CHEMO 				! chemotactic susceptibility of the HI_CHEMO_FRACTION cells
	read(nfcell,*) LO_CHEMO 				! chemotactic susceptibility of the other tagged cells
	close(nfcell)

endif
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine make_outputfilename

if (track_DCvisits) then
	! Filename reporting parameters determining HEV and DC regions
	!write(outputfile,'(a,a,f3.1,a,f3.1,a,f3.1,a,f3.1,a)') 'V_','RHV',R_HEV_min,'_',R_HEV_max, &
	!	'RDC',R_DC_min,'_',R_DC_max,'.out'
	! Filename reporting parameters determining chemotaxis for DC visits, USE_DC_COGNATE = false
	!write(outputfile,'(a,f3.1,a,f3.1,a,f3.1,a,f3.1,a)') 'notraffic_DCF',DC_CHEMO_FRACTION,'_HCF',HI_CHEMO_FRACTION, &
	!	'_LO',LO_CHEMO,'_HI',HI_CHEMO,'.out'
	! Filename reporting parameters determining chemotaxis for DC visits, USE_DC_COGNATE = true
	write(outputfile,'(a,f3.1,a,f3.1,a,f3.1,a,f3.1,a,f3.1,a)') 'DC_CogF',DC_COGNATE_FRACTION,'_TF',DC_CHEMO_FRACTION, &
		'_HCF',HI_CHEMO_FRACTION,'_LO',LO_CHEMO,'_HI',HI_CHEMO,'.out'
else
	outputfile = 'para.out'
endif
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine save_inputfile(cellfile)
character(LEN=*) :: cellfile
character*(128) :: line

open(nfcell,file=cellfile,status='old')
do
    read(nfcell,'(a)',end=999) line
    write(nfout,'(a)') line
enddo
999 continue
close(nfcell)
end subroutine

!----------------------------------------------------------------------------------------
! Save parameters (values hard-coded, not yet in input file) 
!----------------------------------------------------------------------------------------
subroutine save_parameters

write(nfout,*) 'PARAMETER data'
write(nfout,*) '--------------'
write(nfout,*) 'DELTA_X: ',DELTA_X
write(nfout,*) 'BALANCER_INTERVAL: ',BALANCER_INTERVAL
write(nfout,*) 'use_add_count: ',use_add_count
write(nfout,*) 'use_blob: ',use_blob
write(nfout,*) 'use_traffic: ',use_traffic
write(nfout,*) 'use_exit_chemotaxis: ',use_exit_chemotaxis
write(nfout,*) 'use_DC_chemotaxis: ',use_DC_chemotaxis
write(nfout,*) 'computed_outflow: ',computed_outflow
write(nfout,*) 'use_cytokines: ',use_cytokines
write(nfout,*) 'use_cognate: ',use_cognate
write(nfout,*) 'random_cognate: ',random_cognate
write(nfout,*) 'use_diffusion: ',use_diffusion
!write(nfout,*) 'NX: ',NX
write(nfout,*) 'NGEN_EXIT: ',NGEN_EXIT
write(nfout,*) 'constant_efactor: ',constant_efactor
write(nfout,*) 'efactor: ',efactor
write(nfout,*) 'ZIN_FRACTION: ',ZIN_FRACTION
!write(nfout,*) 'ZOUT_FRACTION: ',ZOUT_FRACTION
write(nfout,*) 'use_DC: ',use_DC
write(nfout,*) 'use_DCflux: ',use_DCflux
write(nfout,*) 'DCDIM: ',DCDIM
write(nfout,*) 'DC_DCprox: ',DC_DCprox
write(nfout,*) 'DC_pMHC_THRESHOLD: ',DC_pMHC_THRESHOLD
write(nfout,*) 'ABIND1: ',ABIND1
write(nfout,*) 'ABIND2: ',ABIND2
write(nfout,*) 'NDIFFSTEPS: ',NDIFFSTEPS
write(nfout,*) 'VEGF_MODEL: ',VEGF_MODEL
write(nfout,*) 'VEGF_alpha: ',VEGF_alpha
write(nfout,*) 'VEGF_beta: ',VEGF_beta
write(nfout,*) 'VEGF_decayrate: ',VEGF_decayrate
write(nfout,*) 'vasc_maxrate: ',vasc_maxrate
write(nfout,*) 'vasc_decayrate: ',vasc_decayrate
write(nfout,*) 'vasc_beta: ',vasc_beta
write(nfout,*) 'vasc_n: ',vasc_n
write(nfout,*) 'fix_avidity: ',fix_avidity
write(nfout,*) 'avidity_nlevels: ',avidity_nlevels
write(nfout,*) 'avidity_logscale: ',avidity_logscale
write(nfout,*) 'avidity_min: ',avidity_min
write(nfout,*) 'avidity_step: ',avidity_step

end subroutine

!-----------------------------------------------------------------------------------------
! Read the DC injection schedule file, determine the total number of DCs entering the LN
! in the experiment.  The total for the simulation run is determined by scaling the
! experiment total by the ratio of the initial T cell population numbers for the simulation
! and the experiment.  The total NDCrun must then allocated to the hourly intervals in
! the appropriate proportions.
!-----------------------------------------------------------------------------------------
subroutine setup_DCinjected
integer :: nthours, NTCexpt, NDCexpt, ntimes, i, ihour, n, NDCrun, nsum, ndeficit
real :: ratio
integer, allocatable :: DCexpt(:)

open(nfDCinfo,file=DCinjectionfile,status='old')
nthours = days*24 + 1
allocate(DCinjected(nthours))
allocate(DCexpt(nthours))
DCinjected = 0
DCexpt = 0
read(nfDCinfo,*) NTCexpt
read(nfDCinfo,*) T_DC_INJECTION
read(nfDCinfo,*) ntimes
NDCexpt = 0
do i = 1,ntimes
	read(nfDCinfo,*) ihour, n
	if (ihour > nthours) exit
	DCexpt(ihour) = n
	NDCexpt = NDCexpt + n
enddo
close(nfDCinfo)
ratio = real(NTcells0)/NTCexpt
NDCrun = ratio*NDCexpt + 0.5
write(logmsg,*) 'NTcells0, NTCexpt, ratio, NDCrun: ',NTcells0,NTCexpt,ratio,NDCrun
call logger(logmsg)
nsum = 0
do ihour = 1,nthours
	DCinjected(ihour) = ratio*DCexpt(ihour) + 0.5
	nsum = nsum + DCinjected(ihour)
enddo
ndeficit = NDCrun - nsum
do ihour = 1,nthours
	if (ndeficit == 0) exit
	if (DCinjected(ihour) > 0) then
		if (ndeficit > 0) then
			DCinjected(ihour) = DCinjected(ihour) + 1
			ndeficit = ndeficit - 1
		elseif (ndeficit < 0) then
			DCinjected(ihour) = DCinjected(ihour) - 1
			ndeficit = ndeficit + 1
		endif
	endif
enddo
deallocate(DCexpt)
do ihour = 1,nthours
	write(logmsg,*) ihour,DCinjected(ihour)
	call logger(logmsg)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! When a T cell dies it is removed from the cell list (%ID -> 0)
! and removed from occupancy()%indx()
! The count of sites to add is decremented, for later adjustment of the blob size.
!-----------------------------------------------------------------------------------------
subroutine Tcell_death(kcell)
integer :: kcell
integer :: k, idc, site(3), indx(2), ctype, stype, region
logical :: cognate

!write(logmsg,*) 'Tcell_death: ',kcell
!call logger(logmsg)
cognate = (associated(cellist(kcell)%cptr))
if (cognate) then
	call get_region(cellist(kcell)%cptr,region)
endif
! Is the cell bound?
do k = 1,MAX_DC_BIND
    idc = cellist(kcell)%DCbound(k)
    if (idc /= 0) then
        DClist(idc)%nbound = DClist(idc)%nbound - 1
        if (DClist(idc)%nbound < 0) then
            write(*,*) 'ERROR: Tcell_death: DC nbound < 0: ',kcell,k,idc,cellist(k)%DCbound
            stop
        endif
        if (cognate) then
            DClist(idc)%ncogbound = DClist(idc)%ncogbound - 1
            call RemoveCogbound(idc,kcell)
!            if (idc == idbug) then
!                write(*,*) 'Tcell_death: decremented ncogbound to: cell: ',DClist(idc)%ncogbound,kcell
!            endif
            if (DClist(idc)%ncogbound < 0) then
                write(*,*) 'ERROR: Tcell_death: DC ncogbound < 0: ',kcell,k,idc,cellist(k)%DCbound
                stop
            endif
        endif
    endif
enddo
cellist(kcell)%DCbound = 0
cellist(kcell)%ID = 0
totalres%dN_Dead = totalres%dN_Dead + 1
totalres%N_Dead = totalres%N_Dead + 1
ctype = cellist(kcell)%ctype
!stype = struct_type(ctype)
if (associated(cellist(kcell)%cptr)) then
	stype = COG_TYPE_TAG
else
	stype = NONCOG_TYPE_TAG
endif
if (stype == COG_TYPE_TAG) then
    cognate_list(cellist(kcell)%cptr%cogID) = 0
endif
ngaps = ngaps + 1
if (ngaps > max_ngaps) then
    write(*,*) 'Tcell_death: ngaps > max_ngaps'
    stop
endif
gaplist(ngaps) = kcell

if (region == PERIPHERY) then
	NTcellsPer = NTcellsPer - 1
	return
else
	NTcells = NTcells - 1
endif
site = cellist(kcell)%site
indx = occupancy(site(1),site(2),site(3))%indx
if (indx(1) == kcell) then
    occupancy(site(1),site(2),site(3))%indx(1) = indx(2)
    occupancy(site(1),site(2),site(3))%indx(2) = 0
elseif (indx(2) == kcell) then
    occupancy(site(1),site(2),site(3))%indx(2) = 0
else
    write(logmsg,*) 'ERROR: Tcell_death: cell not at site: ',kcell,site,indx
    call logger(logmsg)
    stop
endif

! Now we need to make a site unavailable
! Remove (make OUTSIDE) a boundary site near to the specified site.
! The criteria for selection of a site to remove are that it is on the blob boundary,
! preferably "sticking out", and that it is vacant.  A site may be made vacant by 
! moving a cell to a nearby available site.
! One way to handle this is to maintain a count of the number of sites to be added(removed).
! At regular intervals the counts can be aggregated and if the total is big enough (+ or -)
! they can all be processed within a gather-scatter cycle.
if (use_add_count) then
    nadd_sites = nadd_sites - 1
else
    write(*,*) 'Tcell_death: No site removal code, use add count'
    stop
endif
end subroutine

!-----------------------------------------------------------------------------------------
! A cognate cell divides.  It has already been determined that there is space for the extra cell.
! One cell retains the cogID, the other gets cogID = 0 for now, to be allocated
! sequentially later, after gather_data() (unless Mnodes = 1, in which case the
! allocation of cogID is done in create_Tcell.)
! The other cell is placed in freeslot of site2 (which may be the same site as kcell).
! If the cell is in the periphery the site, slot etc is irrelevant.
!-----------------------------------------------------------------------------------------
subroutine cell_division(kcell,site2,freeslot,ok)
integer :: kcell, site2(3), freeslot
logical :: ok
integer :: icnew, ctype, gen, region, site(3), indx(2)
integer :: iseq, tag, kfrom, kto
real :: tnow
type(cog_type), pointer :: p1, p2
!real :: IL_state(CYT_NP)
logical :: cognate
integer :: kpar = 0

ok = .true.
tnow = istep*DELTA_T
!call show_cognate_cell(kcell)
p1 => cellist(kcell)%cptr
cognate = .true.
gen = get_generation(p1)
call get_region(p1,region)
!write(*,*) 'cell_division: ',kcell,gen,istep,tnow
if (gen == TC_MAX_GEN) then
    write(logmsg,*) 'cell_division: reached maximum generation: ',kcell
	call logger(logmsg)
    return
endif
if (gen == 1) then
    write(logmsg,*) 'First division at hour: ',tnow/60
    call logger(logmsg)
else
!    write(*,*) 'Time for next division: ',kcell,gen,(tnow-cellist(kcell)%entrytime)/60
endif
ndivided(gen) = ndivided(gen) + 1
tdivided(gen) = tdivided(gen) + (tnow - p1%dividetime)

if (ngaps > 0) then
    icnew = gaplist(ngaps)
    ngaps = ngaps - 1
else
    nlist = nlist + 1
    if (nlist > max_nlist) then
		write(logmsg,*) 'Error: cell_division: cell list full: ',nlist
		call logger(logmsg)
		ok = .false.
		return
	endif
    icnew = nlist
endif
cellist(kcell)%lastdir = random_int(1,6,kpar)
gen = gen + 1
call set_generation(p1,gen)
call set_stage(p1,POST_DIVISION)
if (TCR_splitting) then
    p1%stimulation = p1%stimulation/2
endif
ctype = cellist(kcell)%ctype
tag = 0
p1%dietime = tnow + TClifetime(p1)
p1%dividetime = tnow
if (revised_staging) then
    p1%stagetime = tnow + dividetime(gen,ctype)
else
    p1%stagetime = tnow
endif
if (region == LYMPHNODE) then
	site = cellist(kcell)%site
	indx = occupancy(site(1),site(2),site(3))%indx
endif
call create_Tcell(icnew,cellist(icnew),site2,ctype,cognate,gen,tag,POST_DIVISION,region,ok)
if (.not.ok) return

p2 => cellist(icnew)%cptr
!if (p2%avidity /= p1%avidity) then
!    navid = navid - 1   ! reset counter for setting avidity in create_Tcell
!    p2%avidity = p1%avidity
!endif

p2%stimulation = p1%stimulation
p2%status = p1%status
cellist(icnew)%entrytime = tnow
if (revised_staging) then
    p2%stagetime = tnow + dividetime(gen,ctype)
else
    p2%stagetime = tnow
endif

if (use_cytokines) then
    do iseq = 1,Ncytokines
        tag = cyt_tag(iseq)
        kfrom = NP_offset(iseq)+1
        kto = NP_offset(iseq+1)
        select case(tag)
        case(IL2_TAG)
            call IL2_divide(p1%IL_state(kfrom:kto),p2%IL_state(kfrom:kto),p1%IL_statep(kfrom:kto),p2%IL_statep(kfrom:kto))
        case(IL4_TAG)
            call IL4_divide(p1%IL_state(kfrom:kto),p2%IL_state(kfrom:kto),p1%IL_statep(kfrom:kto),p2%IL_statep(kfrom:kto))
        case(IL7_TAG)
            call IL7_divide(p1%IL_state(kfrom:kto),p2%IL_state(kfrom:kto),p1%IL_statep(kfrom:kto),p2%IL_statep(kfrom:kto))
        case(IL9_TAG)
            call IL9_divide(p1%IL_state(kfrom:kto),p2%IL_state(kfrom:kto),p1%IL_statep(kfrom:kto),p2%IL_statep(kfrom:kto))
        case(IL15_TAG)
            call IL15_divide(p1%IL_state(kfrom:kto),p2%IL_state(kfrom:kto),p1%IL_statep(kfrom:kto),p2%IL_statep(kfrom:kto))
        case(IL21_TAG)
            call IL21_divide(p1%IL_state(kfrom:kto),p2%IL_state(kfrom:kto),p1%IL_statep(kfrom:kto),p2%IL_statep(kfrom:kto))
        end select
    enddo
    !p2%IL_state = p1%IL_state
endif
p2%CD69 = p1%CD69       ! for now just assume replication of the CD69 and S1P1 expression
p2%S1P1 = p1%S1P1

ndivisions = ndivisions + 1
if (region /= LYMPHNODE) then
	NTcellsPer = NTcellsPer + 1
	return
endif

occupancy(site2(1),site2(2),site2(3))%indx(freeslot) = icnew
if (use_add_count) then
    nadd_sites = nadd_sites + 1
else
    write(logmsg,*) 'cell_division: No site removal code, use add count'
	call logger(logmsg)
	ok = .false.
    return
endif
NTcells = NTcells + 1
if (IN_VITRO .and. NTcells > 0.95*n2Dsites) then
	ok = .false.
	write(logmsg,'(a,2i6)')'More than 95% of the available space is occupied by T cells: ',NTcells,n2Dsites
	call logger(logmsg)
endif
end subroutine

!--------------------------------------------------------------------------------
! Place n DCs
! Note that we may have NDCsites > 1 !!!
! The procedure for placing DCs is as follows:
! The central DC is placed on a site that satisfies these conditions:
! (1) The site is not OUTSIDE_TAG or a DC site
! (2) There are not two T cells at this site
! (3) The site is not too near another DC, i.e. it is more than DC_DCprox*DC_RADIUS
!     from any DC.
! (4) The site is not too near an exit, i.e. it is more than exit_DCprox from any exit.
! (5) The site is not too near the blob boundary, i.e. it is more than
!     bdry_DCprox*DC_RADIUS from the sphere with radius = Radius
! (6) Either the site is free, or it is occupied by a T cell that can be moved
!     to a neighbouring site.
! A site meeting these requirements is selected for the DC, then as many as
! possible of the neighbouring sites are also allocated to the NDCsites-1 other
! sites occupied by a DC, by subroutine addDCsite().  The count of DC sites
! allocated is stored in DC%nsites.
!--------------------------------------------------------------------------------
subroutine place_DCs(n,nadded)
integer :: n, nadded
integer :: i, x, y, z, site1(3), site2(3), freeslot, err
integer :: indx(2), jslot, idc, k, kcell
integer :: xmin, xmax, ymin, ymax, zmin, zmax
integer :: kpar = 0
real(DP) :: R
real :: tnow, dist, rvec(3), prox, tmins
logical :: OK
type(DC_type) :: DC
integer, parameter :: kmax = 100000

!ndbug = DClist(idbug)%ncogbound
!write(*,*) 'place_DCs (1): ',DClist(idbug)%ncogbound
tnow = istep*DELTA_T
xmin = x0 - Radius
xmax = x0 + Radius + 1
ymin = y0 - Radius
ymax = y0 + Radius + 1
zmin = z0 - Radius
zmax = z0 + Radius + 1
xmin = max(1,xmin)
xmax = min(NX,xmax)
ymin = max(1,ymin)
ymax = min(NY,ymax)
zmin = max(1,zmin)
zmax = min(NZ,zmax)
nadded = 0
do i = 1,n
	OK = .false.
	k = 0
    do
		k = k+1
		if (k > kmax) then
			write(logmsg,*) 'Error: place_DCs: unable to find space for a DC'
			call logger(logmsg)
			OK = .false.
			stop
		endif
        R = par_uni(kpar)
        x = xmin + R*(xmax-xmin)
        R = par_uni(kpar)
        y = ymin + R*(ymax-ymin)
        R = par_uni(kpar)
        z = zmin + R*(zmax-zmin)
        site1 = (/x,y,z/)   ! global location
        indx = occupancy(x,y,z)%indx
        if (indx(1) < 0) cycle                          ! OUTSIDE_TAG or DC
        if (indx(1) /= 0 .and. indx(2) /= 0) cycle      ! two T cells at this site
!        prox = proximity_limit(DC_SITE,DC_SITE)
!        if (tooNearDC(site1,prox)) cycle
!        prox = proximity_limit(DC_SITE,EXIT_SITE)
!        if (tooNearExit(site1,prox)) cycle
!!        prox = 0.5*DC_DCprox*DC_RADIUS
!!        prox = bdry_DCprox*DC_RADIUS
!        prox = proximity_limit(DC_SITE,BDRY_SITE)
!        if (.not.tooNearBdry(site1,prox)) then
		call CheckSite(DC_SITE,site1,ok)
		if (ok) then
            jslot = 0
            if (indx(1) /= 0) then
                jslot = 1
            elseif (indx(2) /= 0) then
                jslot = 2
            endif
            if (jslot == 0) then ! free site, use it
                exit
            else  ! one T cell here in jslot, it must be moved
                ! Can the T cell be bumped to a neighbour site?
                call get_free_slot(occupancy,NX,site1,site2,freeslot)
                if (freeslot == 0) cycle    ! cannot be bumped, try again
               ! Move the cell in site1/jslot to site2/freeslot
                kcell = indx(jslot)
!                write(*,*) 'place_DCs: bump cell: ',kcell,site1,site2
                occupancy(x,y,z)%indx = 0
                occupancy(site2(1),site2(2),site2(3))%indx(freeslot) = kcell
                cellist(kcell)%site = site2
                ! Now site1 is free to use
                exit
            endif
        endif
    enddo
!    if (DClist(idbug)%ncogbound /= ndbug) then
!        write(*,*) 'place_DCs (a): ndbug changed: ',i,DClist(idbug)%ncogbound
!        stop
!    endif
    if (.not.OK) exit
    idc = 0
    if (reuse_DC_index) then
	    do k = 1,NDC
		    if (.not.DClist(k)%alive) then
		        idc = k
		        exit
		    endif
		enddo
	endif
    if (idc == 0) then    ! If there isn't a free spot in DClist()
        NDC = NDC + 1
        if (NDC > max_DC) then
			write(logmsg,'(a,i6)') 'Error: place_DCs: number of DCs exceeds limit: ',max_DC
			call logger(logmsg)
			ok = .false.
			return
		endif
        idc = NDC
    endif
    DC%ID = idc
    DC%alive = .true.
    DC%capable = .true.
    DC%site = site1
    DC%nsites = 1
    DC%stimulation = 0
    DC%nbound = 0
    DC%ncogbound = 0
    allocate(DC%cogbound(MAX_COG_BIND))
    DC%density = DCdensity(kpar)
    if (DC_INJECTION) then
        ! If the DCs were injected into an experimental animal, we need to account for the antigen decay since
        ! the time of injection, which was at T_DC_INJECTION hours
        tmins = -T_DC_INJECTION*60 + istep*DELTA_T
        DC%density = DC%density*exp(-DCdecayrate*tmins)
    endif
    DC%dietime = tnow + DClifetime(kpar)
    occupancy(site1(1),site1(2),site1(3))%indx(1) = -idc
    do k = 2,NDCsites
        call addDCsite(idc,site1,k,err)
!        if (DClist(idbug)%ncogbound /= ndbug) then
!            write(*,*) 'place_DCs (b): ndbug changed: ',i,k,DClist(idbug)%ncogbound
!            stop
!        endif
        if (err == 0) then
            DC%nsites = DC%nsites + 1
        else
!            write(*,*) 'addDCsite: idc,k,err: ',idc,k,err
        endif
    enddo
    DClist(idc) = DC
    nadded = nadded + 1
!    write(*,*) 'Added DC at: ',site1,' with nsites: ',DC%nsites
    ! now the DC proximity data in occupancy()%DC must be updated
    ! this is done by reassign_DC() called from balancer()
!    if (DClist(idbug)%ncogbound /= ndbug) then
!        write(*,*) 'place_DCs (c): ndbug changed: ',i,DClist(idbug)%ncogbound
!        stop
!    endif
enddo
NDCtotal = NDCtotal + nadded

! Now check the DC locations wrt the blob perimeter
do k = 1,NDC
    if (DClist(k)%alive) then
        rvec = DClist(k)%site - (/x0,y0,z0/)
        dist = norm(rvec)
        if (dist > Radius) then
            write(logmsg,*) 'Place_DCs: warning: DC distance: ',k,dist,Radius
			call logger(logmsg)
        endif
    endif
enddo
!write(*,*) 'place_DCs (2): ',DClist(idbug)%ncogbound
end subroutine

!-----------------------------------------------------------------------------------------
! Convert one of the sites near a DC into a DC peripheral site for DC idc.
! The index k indicates which peripheral site of 2:NDCsites to convert.
! A site is a candidate for a DC site provided:
! (1) It is not outside the blob, or already a DC site
! (2) It does not hold two T cells
! (3) It is not too close the the blob boundary
! If a candidate site is free, it is used.  If it holds a single T cell, the cell
! is moved if possible to a neighbouring site.
! Note that when a T cell is bumped it retains its binding (if any) to a DC, even though
! it may have been moved to a site that is - strictly speaking - not within the SOI of
! the DC.
!-----------------------------------------------------------------------------------------
subroutine addDCsite(idc,site0,k,err)
integer :: idc, site0(3), k, err
integer :: indx(2), jslot, kcell, freeslot, site1(3), site2(3)
real :: prox
logical :: OK

site1 = site0 + DCoffset(:,k)
indx = occupancy(site1(1),site1(2),site1(3))%indx
if (indx(1) < 0) then                       ! OUTSIDE_TAG or DC
    err = 1
    return
endif
if (indx(1) /= 0 .and. indx(2) /= 0) then   ! two T cells at this site
    err = 2
    return
endif
!if (tooNearDC(site1,NDC)) then    ! too close to another DC
!    err = 3
!    return
!endif
!prox = 0.5*DC_DCprox*DC_RADIUS
prox = bdry_DCprox*DC_RADIUS
if (tooNearBdry(site1,prox)) then                ! too close to the blob boundary
    err = 4
    return
endif

OK = .false.
jslot = 0
if (indx(1) /= 0) then
    jslot = 1
elseif (indx(2) /= 0) then
    jslot = 2
endif
if (jslot == 0) then ! free site, use it
    OK = .true.
else  ! one T cell here in jslot, it must be moved
    ! Can the T cell be bumped to a neighbour site?
    call get_free_slot(occupancy,NX,site1,site2,freeslot)
    if (freeslot == 0) then    ! cannot be bumped
        err = 5
        return
    endif
   ! Move the cell in site1/jslot to site2/freeslot
    kcell = indx(jslot)
    occupancy(site1(1),site1(2),site1(3))%indx = 0
    occupancy(site2(1),site2(2),site2(3))%indx(freeslot) = kcell
    cellist(kcell)%site = site2
    ! Now site1 is free to use
    OK = .true.
endif
!if (OK) then
    err = 0
    occupancy(site1(1),site1(2),site1(3))%indx = -idc
!else
!    err = 6
!endif
end subroutine

!-----------------------------------------------------------------------------------------
! Locate a free slot in a site adjacent to site1: site2 (freeslot)
! Returns freeslot = 0 if there is no free space in an adjacent site.
! The occupancy array occ() can be either occupancy() or big_occupancy(), hence
! the need for xlim (= NXX or NX)
!-----------------------------------------------------------------------------------------
subroutine get_free_slot(occ,xlim,site1,site2,freeslot)
type(occupancy_type) :: occ(:,:,:)
integer :: xlim, site1(3), site2(3), freeslot
logical :: Lwall, Rwall
integer :: i, indx2(2)

if (site1(1) == 1) then
    Lwall = .true.
else
    Lwall = .false.
endif
if (site1(1) == xlim) then
    Rwall = .true.
else
    Rwall = .false.
endif

if (site1(1) < 1) then
    write(logmsg,*) 'get_free_slot: bad site1: ',site1
	call logger(logmsg)
    stop
endif
do i = 1,27
    if (i == 14) cycle       ! i = 14 corresponds to the no-jump case
	site2 = site1 + jumpvec(:,i)
    if (Lwall .and. site2(1) < 1) then
        cycle
    endif
    if (Rwall .and. site2(1) > xlim) then
        cycle
    endif
    if (site2(2) < 1 .or. site2(2) > NY .or. site2(3) < 1 .or. site2(3) > NZ) cycle
	indx2 = occ(site2(1),site2(2),site2(3))%indx
	if (indx2(1) >= 0) then         ! not OUTSIDE_TAG or DC
	    if (indx2(1) == 0) then     ! slot 1 is free
	        freeslot = 1
            return
	    elseif (indx2(2) == 0) then ! slot 2 is free
	        freeslot = 2
            return
        endif
    endif
enddo
freeslot = 0
end subroutine

!-----------------------------------------------------------------------------------------
! All associations of sites with DC are recomputed.
! That is for every site (x,y,z), create the list of DC that are near this site:
! occupancy(x,y,z)%DC(1:3), where occupancy(x,y,z)%DC(0) = number of nearby DC (if >= 0)
! To avoid having to explicitly select the closest DC (when there are more than DCDIM-1 near
! a site), the order of scanning the DClist is randomized.
!-----------------------------------------------------------------------------------------
subroutine reassign_DC(kpar,ok)
integer :: kpar
logical :: ok
integer, allocatable :: perm(:)
integer :: xdc, ydc, zdc, xmin, xmax, ymin, ymax, zmin, zmax
integer :: idc, kdc, k, x, y, z, y2, z2, d2, site(3), nassigned
logical :: added

!write(*,*) 'reassign_DC'

occupancy(:,:,:)%DC(0) = 0
occupancy(:,:,:)%cDC(0) = 0
allocate(perm(MAX_DC))
do k = 1, NDC
    perm(k) = k
enddo

call permute(perm,NDC,kpar)
NDCalive = 0
do k = 1,NDC
    idc = perm(k)
    if (.not.DClist(idc)%alive) cycle
!    write(*,*) 'idc: ',idc
!	write(*,*) '(4) cell 16243: ',cellist(16243)%site,occupancy(70,63,88)%indx
	if (DClist(idc)%nsites < NDCsites) then
		call assignDCsites(idc,nassigned,ok)
		if (.not.ok) return
!		write(*,*) 'assigned DC sites for: ',idc,nassigned
		DClist(idc)%nsites = nassigned
	endif
    NDCalive = NDCalive + 1
    site = DClist(idc)%site
    xdc = site(1)
    ydc = site(2)
    zdc = site(3)
    xmin = xdc - DC_RADIUS
    xmax = xdc + DC_RADIUS
    ymin = ydc - DC_RADIUS
    ymax = ydc + DC_RADIUS
    zmin = zdc - DC_RADIUS
    zmax = zdc + DC_RADIUS
    xmin = max(1,xmin)
    xmax = min(NX,xmax)
    ymin = max(1,ymin)
    ymax = min(NY,ymax)
    zmin = max(1,zmin)
    zmax = min(NZ,zmax)
    do z = zmin,zmax
        z2 = (z-zdc)*(z-zdc)
        do y = ymin,ymax
            y2 = (y-ydc)*(y-ydc)
            do x = xmin,xmax
                d2 = (x-xdc)*(x-xdc) + y2 + z2
                added = .false.
                ! The following procedure is correct only if DC_RADIUS <= chemo_radius
	            if (d2 <= DC_RADIUS*DC_RADIUS) then
	                kdc = occupancy(x,y,z)%DC(0)
	                if (kdc < 0) cycle			! this is a DC site
		            if (kdc < DCDIM-1) then		! not all possible DC assigned
						kdc = kdc + 1
						occupancy(x,y,z)%DC(kdc) = idc
						occupancy(x,y,z)%DC(0) = kdc
						added = .true.
					endif
				endif
	            if (.not.added .and. (d2 <= chemo_radius*chemo_radius)) then
	                kdc = occupancy(x,y,z)%cDC(0)
		            if (kdc == cDCDIM-1) cycle   ! all possible DC assigned
		            kdc = kdc + 1
		            occupancy(x,y,z)%cDC(kdc) = idc
	                occupancy(x,y,z)%cDC(0) = kdc
	            endif
            enddo
        enddo
    enddo
enddo
deallocate(perm)
ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
! Try to expand any DCs that are not yet occupying their full allocation of sites.
!-----------------------------------------------------------------------------------------
subroutine growDC
integer :: k, idc, err, site0(3), site(3)

do idc = 1,NDC
    if (.not.DClist(idc)%alive) cycle
    if (DClist(idc)%nsites < NDCsites) then
        site0 = DClist(idc)%site
        do k = 2,NDCsites
            site = DClist(idc)%site + DCoffset(:,k)
            if (occupancy(site(1),site(2),site(3))%indx(1) /= -idc) then
                call addDCsite(idc,site0,k,err)
                if (err == 0) then
                    DClist(idc)%nsites = DClist(idc)%nsites + 1
                endif
            endif
        enddo
!        write(*,*) 'growDC: ',idc,DClist(idc)%nsites
    endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Test moveDC() with a randomly selected DC, random step direction.
!-----------------------------------------------------------------------------------------
subroutine test_moveDC
integer :: idc, dir
integer :: kpar=0

write(*,*) 'test_moveDC'
if (NDCsites /= 7) then
    write(*,*) 'Valid only for NDCsites = 7'
    stop
endif
do
    idc = random_int(1,NDC,kpar)
    if (DClist(idc)%alive) exit
enddo
dir = random_int(1,6,kpar)
call moveDC(idc,dir)
end subroutine

!-----------------------------------------------------------------------------------------
! A DC is allowed to move only if it has grown to its full extent,
! i.e. if DClist()%nsites = NDCsites
! For now allow only moves in the 6 principal directions (Neumann).
! Currently valid only for NDCsites = 7.
! Cells in the 5 sites in the path of the DC step are moved.
! For 4 sites, the shift is by one site, for the one site in line with the DC centre the
! shift is 3 sites.
!-----------------------------------------------------------------------------------------
subroutine moveDC(idc,dir)
integer :: idc, dir
integer :: k, i, kcell, indx(2), site0(3), site1(3), site2(3), step(3)

if (DClist(idc)%nsites /= NDCsites) return
step = neumann(:,dir)
site0 = DClist(idc)%site
!write(*,*) 'moveDC: ',idc,dir,'  ',site0,'  ',step
do k = 2,NDCsites
    if (all(DCoffset(:,k) == step)) then       ! this is in the direction of step
        ! move site contents by 3 sites in opposite direction to step
        site1 = site0 - step       ! old DC site
        site2 = site0 + 2*step     ! new DC site
!        write(*,'(a,7i4)') 'step dir: ',k,site1,site2
        indx = occupancy(site2(1),site2(2),site2(3))%indx
        if (any(indx == OUTSIDE_TAG)) then
            write(*,*) 'moveDC: site2 outside!'
            stop
        endif
        occupancy(site1(1),site1(2),site1(3))%DC = occupancy(site2(1),site2(2),site2(3))%DC
        occupancy(site1(1),site1(2),site1(3))%indx = indx
        do i = 1,2
            kcell = indx(i)
            if (kcell /= 0) then
                cellist(kcell)%site = site1
            endif
        enddo
        occupancy(site2(1),site2(2),site2(3))%indx = -idc
        occupancy(site0(1),site0(2),site0(3))%indx = -idc
        site2 = site0  + step
        occupancy(site2(1),site2(2),site2(3))%indx = -idc
    elseif (all(DCoffset(:,k) == -step)) then   ! this is the direction opposite to step
        ! do nothing
    else                                        ! one of the other 4 directions
        ! move site contents by 1 site in opposite direction to step
        site1 = site0 + DCoffset(:,k)       ! old DC site, new T cell site
        site2 = site1 + step                ! new DC site, old T cell site
!        write(*,'(a,7i4)') 'other dir: ',k,site1,site2
        indx = occupancy(site2(1),site2(2),site2(3))%indx
        if (any(indx == OUTSIDE_TAG)) then
            write(*,*) 'moveDC: site2 outside!'
            stop
        endif
        occupancy(site1(1),site1(2),site1(3))%DC = occupancy(site2(1),site2(2),site2(3))%DC
        occupancy(site1(1),site1(2),site1(3))%indx = indx
        do i = 1,2
            kcell = indx(i)
            if (kcell /= 0) then
                cellist(kcell)%site = site1
            endif
        enddo
        occupancy(site2(1),site2(2),site2(3))%indx = -idc
    endif
enddo
DClist(idc)%site = site0 + step
end subroutine

!-----------------------------------------------------------------------------------------
! Create a new T cell.
! We could give a progeny cell (the result of cell division) the same ID as its parent.
!-----------------------------------------------------------------------------------------
subroutine create_Tcell(kcell,cell,site,ctype,cognate,gen,tag,stage,region,ok)
type(cell_type) :: cell
integer :: kcell, site(3), ctype, gen, tag, stage, region
logical :: cognate, ok
integer :: stype, cogID, i
real :: tnow, param1, param2, R
integer :: kpar = 0

ok = .true.
tnow = istep*DELTA_T
!stype = struct_type(ctype)
!write(*,*) 'create_Tcell: ',kcell,ctype,stype
cell%entrytime = tnow
if (.not.cognate) then
    if (dbug) then
        write(*,*) 'create_Tcell: ',kcell,site,associated(cell%cptr)
    endif
    if (associated(cell%cptr)) then
        deallocate(cell%cptr)
    endif
else
    if (.not.associated(cell%cptr)) then
        allocate(cell%cptr)
    endif
    param1 = log(TC_AVIDITY_MEDIAN)
    param2 = log(TC_AVIDITY_SHAPE)
    cell%cptr%status = 0
    call set_generation(cell%cptr,gen)
    call set_stage_region(cell%cptr,stage,region)
    if (fix_avidity) then
        i = mod(navid,avidity_nlevels)
        navid = navid + 1
        cell%cptr%avidity = avidity_level(i+1)
!        if (avidity_logscale) then
!            cell%cptr%avidity = 10**(avidity_min + i*avidity_step)
!        else
!            cell%cptr%avidity = avidity_min + i*avidity_step
!        endif
!        write(nfout,'(a,4i7,2f8.4)') 'create_Tcell: add: ',i,navid,kcell,lastcogID + 1, &
!            cell%cptr%avidity,cellist(kcell)%cptr%avidity
!        if (i == 3 .and. istep >= 4500) then
!            avid_debug = .true.
!        endif
    else
        cell%cptr%avidity = rv_lognormal(param1,param2,kpar)
    endif
    cell%cptr%stimulation = 0
!    cell%cptr%IL_state = 0
!    cell%cptr%IL_statep = 0
    if (use_cytokines) then
        call IL2_init_state(cell%cptr%IL_state,cell%cptr%IL_statep)
    endif
    ! What should the initial CD69 level be?  If 0, can a cognate cell exit immediately?
    ! We would prefer not to impose a time or generation constraint on the exit of
    ! cognate T cells, but otherwise if CD69 is initially 0, a cell will be susceptible
    ! to chemotaxis and exit until it has received enough TCR signal to drive CD69 high.
    cell%cptr%CD69 = 0
    cell%cptr%S1P1 = 0
!	cell%cptr%DCchemo = BASE_DCchemo
    cell%cptr%dietime = tnow + TClifetime(cell%cptr)
    cell%cptr%dividetime = tnow
    cell%cptr%stagetime = BIG_TIME
	cell%cptr%cnt = 0

! Maintain cognate_list at start or if we are running on a single node
! Otherwise cogID and cognate_list is maintained by make_cognate_list
!    if (istep == 0 .or. Mnodes == 1) then
        lastcogID = lastcogID + 1
        if (lastcogID > MAX_COG) then
            write(logmsg,'(a,i6)') 'Error: create_Tcell: cognate_list dimension exceeded: ',MAX_COG
            call logger(logmsg)
            ok = .false.
            return
        endif
        cogID = lastcogID
        cell%cptr%cogID = cogID
        cognate_list(cogID) = kcell
!    else
!        cell%cptr%cogID = 0
!    endif
endif
lastID = lastID + 1     ! Each node makes its own numbers, with staggered offset
cell%ID = lastID
cell%site = site
cell%ctype = ctype
cell%tag = tag
cell%step = 0
cell%DCbound = 0
cell%unbindtime = -1000
cell%lastdir = random_int(1,6,kpar)
!cell%DCchemo = BASE_DCchemo
cell%receptor_level(CCR1) = BASE_DCchemo
cell%receptor_saturation_time = 0
! For recording DC scanning statistics
if (track_DCvisits .and. tag == TAGGED_CELL) then
    ntagged = ntagged + 1
    if (.not.allocated(cell%dclist)) then
        allocate(cell%dclist(max(1,NDC)))
    endif
    cell%dclist = 0
    cell%ndclist = 0
    cell%visits = 0
    cell%revisits = 0
!    if (DC_CHEMO_NOTRAFFIC) then	! also needed for use_traffic case
		R = par_uni(kpar)
		if (R < HI_CHEMO_FRACTION) then
!			cell%DCchemo = HI_CHEMO
			cell%receptor_level(CCR1) = HI_CHEMO
		else
!			cell%DCchemo = LO_CHEMO
			cell%receptor_level(CCR1) = LO_CHEMO
!		endif
	endif
endif

end subroutine

!--------------------------------------------------------------------------------
! Add a cell (kcell) with characteristics (ctype, gen, stage) at site.
!--------------------------------------------------------------------------------
subroutine add_Tcell(site,ctype,cognate,gen,tag,stage,region,kcell,ok)
integer :: site(3), ctype, gen, tag, stage, region, kcell
logical :: cognate, ok
integer :: indx(2)

ok = .true.
if (ngaps > 0) then
    kcell = gaplist(ngaps)
    ngaps = ngaps - 1
else
    nlist = nlist + 1
    if (nlist > max_nlist) then
		write(logmsg,*) 'Error: add_Tcell: cell list full: ',nlist
		call logger(logmsg)
		ok = .false.
		return
	endif
    kcell = nlist
endif
if (dbug) then
    write(*,'(a,9i7,L2)') 'add_Tcell: ',istep,kcell,site,ctype,gen,stage,region,cognate
endif
call create_Tcell(kcell,cellist(kcell),site,ctype,cognate,gen,tag,stage,region,ok)
if (.not.ok) return

indx = occupancy(site(1),site(2),site(3))%indx
if (indx(1) == 0) then
    indx(1) = kcell
elseif (indx(2) == 0) then
    indx(2) = kcell
else
    write(*,*) 'ERROR: add_Tcell: no free slot: ',site,indx
    stop
endif
occupancy(site(1),site(2),site(3))%indx = indx

end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine add_site_local(site)
integer :: site(3)

if (use_add_count) then
    nadd_sites = nadd_sites - 1
else
    write(*,*) 'add_site_local: no code to add site, use add count'
    stop
endif
end subroutine

!--------------------------------------------------------------------------------
! Add a vacant site at a boundary to account for the T cell added at site
! In the case of an end node (me = 0 or me = Mnodes-1) the added site can be
! at a different x value, but for internal nodes the added site  must go at
! (or near) the same x value.
! For the spherical blob case, try to place the site at a point near the bdry
! on a line drawn through site from the centre.
! NOT USED
!--------------------------------------------------------------------------------
subroutine add_vacant_site(site,kpar)
integer :: site(3),kpar
integer :: k, site0(3),newsite(3)
real(DP) :: R
real :: dxyz(3)
logical :: redo

if (dbug) write(*,'(a,4i6)') 'add_vacant_site: ',site
site0 = site
dxyz = real(site0) - Centre
do k = 2,3
!    call random_number(R)
    R = par_uni(kpar)
    dxyz(k) = dxyz(k) + (R - 0.5)
enddo
call normalize(dxyz)
redo = .false.
k = 0
do
    k = k+1
    newsite = site0 + k*0.5*dxyz
    if (newsite(1) < 1 .or. newsite(1) > NX) then
        site0 = site0 + (k-1)*0.5*dxyz
        dxyz(1) = 0
        call normalize(dxyz)
        redo = .true.
        exit
    endif
    if (newsite(2) < 1 .or. newsite(2) > NY .or. newsite(3) < 1 .or. newsite(3) > NZ) then
        write(*,*) 'ERROR: add_vacant_site: reached grid limits (a): ',k,site,dxyz
        stop
    endif
    if (occupancy(newsite(1),newsite(2),newsite(3))%indx(1) == OUTSIDE_TAG) then
        exit
    endif
enddo
if (redo) then
    if (dbug) write(*,*) 'redo: ', site0
    k = 0
    do
        k = k+1
        newsite = site0 + k*0.5*dxyz
        if (newsite(2) < 1 .or. newsite(2) > NY .or. newsite(3) < 1 .or. newsite(3) > NZ) then
            write(*,*) 'ERROR: add_vacant_site: reached grid limits (b): ',k,site,dxyz
            newsite = site0 + (k-1)*0.5*dxyz
            write(*,*) newsite,occupancy(newsite(1),newsite(2),newsite(3))%indx(1)
            stop
        endif
        if (occupancy(newsite(1),newsite(2),newsite(3))%indx(1) == OUTSIDE_TAG) then
            exit
        endif
    enddo
endif
if (dbug) write(*,'(a,4i6)') 'newsite: ',newsite
occupancy(newsite(1),newsite(2),newsite(3))%indx = 0
if (dbug) write(*,'(a,7i6)') 'site, vacant site: ',site,newsite

end subroutine



!--------------------------------------------------------------------------------
! nb is the current number of T cells bound to the DC.  Probability of binding
! varies from 1 at nb <= BIND_ALFA*MAX_TC_BIND to 0 at nb >= MAX_TC_BIND
!--------------------------------------------------------------------------------
!real function bindprob(nb)
!integer(2) :: nb
!
!if (nb >= MAX_TC_BIND) then
!	bindprob = 0.0
!elseif (nb < DC_BIND_ALFA*MAX_TC_BIND) then
!	bindprob = 1.0
!else
!	bindprob = 1.0 - (nb - DC_BIND_ALFA*MAX_TC_BIND)/((1.0-DC_BIND_ALFA)*MAX_TC_BIND)
!endif
!end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
subroutine test_exp
integer :: i, N=100
real(DP) :: R, sum
integer :: kpar = 0

sum = 0
write(*,*) 'start'
do i = 1,N
	R = par_rexp(kpar)
	write(*,'(i6,f12.8)') i,R
	sum = sum + R
enddo
write(*,*) 'mean = ',sum/N
sum = 0
write(*,*) 'start'
do i = 1,N
!	call random_number(R)
	R = par_uni(kpar)
	R = -log(R)
	write(*,'(i6,f12.8)') i,R
	sum = sum + R
enddo
write(*,*) 'mean = ',sum/N
end subroutine

!-----------------------------------------------------------------------------------------
! Initial cell position data is loaded into occupancy() and cellist().
!-----------------------------------------------------------------------------------------
subroutine PlaceCells(ok)
logical :: ok
integer :: id, cogid, x, y, z, site(3), ctype
integer :: idc, kdc, k, x2, y2, z2, gen, tag, stage, region, ncogDC
integer :: xdc, ydc, zdc, xmin, xmax, ymin, ymax, zmin, zmax, nzlim, nassigned, NDCrequired
integer :: cnt(2)
real(DP) :: R
real :: d, d2, p1, p2, tnow, prox, iv_fraction=0, tmins
integer, allocatable :: permc(:)
logical :: added, done, cognate
integer :: kpar = 0

!write(logmsg,*) 'PlaceCells: Radius: ',Radius
!call logger(logmsg)
ok = .false.

cnt = 0
IDTEST = 0
tnow = 0
nlist = 0
n2Dsites = 0
do x = 1,NX
    do y = 1,NY
	    do z = 1,NZ
            occupancy(x,y,z)%DC = 0
            occupancy(x,y,z)%cDC = 0
            occupancy(x,y,z)%indx = 0
            occupancy(x,y,z)%exitnum = 0
            occupancy(x,y,z)%isrc = 0
            site = (/x,y,z/)
            d = cdistance(site)
            if (IN_VITRO) then
				if (d > Radius) then
					occupancy(x,y,z)%indx = OUTSIDE_TAG
					zrange2D(x,y,:) = 0
				elseif (z == 1) then	! count only layer 1 sites
					n2Dsites = n2Dsites + 1
					zrange2D(x,y,:) = 1
				endif
            else
				if (use_blob .and. d > Radius) then
					occupancy(x,y,z)%indx = OUTSIDE_TAG
				else
					nlist = nlist + 1
				endif
            endif
        enddo
    enddo
enddo
! At this stage n2Dsites is the number of sites within the 2D dish.

if (IN_VITRO) then
	iv_fraction = IV_NTCELLS/real(n2Dsites)
	if (iv_fraction > 0.9) then
		call logger("Error: PlaceCells: Too many cells, reduce IV_NTCELLS")
		ok = .false.
		return
	endif
	nlist = IV_NTCELLS
endif

!
! Note: Need to account for DCs with cognate and non-cognate antigen.
! In the case of DC_INJECTION, there needs to be a distinction between the DCs initially
! resident in the LN, and injected DCs.
!
NDC = 0
NDCalive = 0
if (TC_TO_DC > 0) then   ! DC placement needs to be checked to account for SOI overlap
!if (use_DC) then   ! DC placement needs to be checked to account for SOI overlap
    NDCrequired = DC_FACTOR*nlist/TC_TO_DC
    if (use_single_DC) then
		NDCrequired = 1
	endif
    NDCtotal = NDCrequired
    if (NDCrequired > MAX_DC) then
        write(logmsg,'(a,i6)') 'Error: PlaceCells: NDC exceeds MAX_DC: ',MAX_DC
        call logger(logmsg)
		ok = .false.
        return
    endif
    p1 = log(DC_ANTIGEN_MEDIAN)
    p2 = log(DC_ANTIGEN_SHAPE)
    do idc = 1,NDCrequired
        do
            R = par_uni(kpar)
            x = 1 + R*NX
            R = par_uni(kpar)
            y = 1 + R*NY
			if (IN_VITRO) then
				z = 1
			else
	            R = par_uni(kpar)
		        z = 1 + R*NZ
		    endif
		    if (use_single_DC) then
				x = NX/2
				y = NY/2
				z = NZ/2
			endif
            site = (/x,y,z/)
            if (occupancy(x,y,z)%indx(1) < 0) cycle     ! OUTSIDE_TAG or DC
			call CheckSite(DC_SITE,site,ok)
			if (ok) exit
        enddo
        DClist(idc)%ID = idc
        DClist(idc)%site = site
        DClist(idc)%nsites = 1
        NDC = NDC + 1
		NDCalive = NDC
! Changed the treatment of injected DCs.  Now assume that all DCs originally in the paracortex are non-cognate,
! while the schedule of arriving DCs with peptide is read from an input file.
!        if (DC_INJECTION) then
!            ! If the DCs were injected into an experimental animal, we need to account for the antigen decay since
!            ! the time of injection, which was at T_DC_INJECTION hours
!            tmins = -T_DC_INJECTION*60
!            DClist(idc)%density = DClist(idc)%density*exp(-DCdecayrate*tmins)
!        endif
	    DClist(idc)%dietime = tnow + DClifetime(kpar)
	    DClist(idc)%alive = .true.
	    DClist(idc)%stimulation = 0
	    DClist(idc)%nbound = 0
	    DClist(idc)%ncogbound = 0
	    allocate(DClist(idc)%cogbound(MAX_COG_BIND))
		DClist(idc)%cognate = .true.
        if (DC_INJECTION) then
		    DClist(idc)%capable = .false.
	        DClist(idc)%density = 0
        else
		    DClist(idc)%capable = .true.
	        DClist(idc)%density = DCdensity(kpar)
		endif
    enddo
    do idc = 1,NDC
		call assignDCsites(idc,nassigned,ok)
		if (.not.ok) return
        DClist(idc)%nsites = nassigned
        nlist = nlist - nassigned
		site = DClist(idc)%site
		xdc = site(1)
		ydc = site(2)
		zdc = site(3)
		xmin = xdc - DC_RADIUS
        xmax = xdc + DC_RADIUS
        ymin = ydc - DC_RADIUS
        ymax = ydc + DC_RADIUS
        zmin = zdc - DC_RADIUS
        zmax = zdc + DC_RADIUS
        xmin = max(1,xmin)
        xmax = min(NX,xmax)
        ymin = max(1,ymin)
        ymax = min(NY,ymax)
        zmin = max(1,zmin)
        zmax = min(NZ,zmax)
        do x = xmin,xmax
	        x2 = (x-xdc)*(x-xdc)
	        do y = ymin,ymax
		        y2 = (y-ydc)*(y-ydc)
		        do z = zmin,zmax
			        z2 = (z-zdc)*(z-zdc)
			        d2 = x2 + y2 + z2
			        added = .false.
					! The following procedure is correct only if DC_RADIUS <= chemo_radius
			        if (d2 <= DC_RADIUS*DC_RADIUS) then
			            kdc = occupancy(x,y,z)%DC(0)
				        if (kdc < 0) cycle   ! don't touch a DC site
			            if (kdc < DCDIM-1) then     ! can add to the list
			                kdc = kdc + 1
			                occupancy(x,y,z)%DC(0) = kdc
			                occupancy(x,y,z)%DC(kdc) = idc
			                added = .true.
				        endif
				    endif
			        if (.not.added .and. (d2 <= chemo_radius*chemo_radius)) then
			            kdc = occupancy(x,y,z)%cDC(0)
			            if (kdc == cDCDIM-1) cycle     ! can add no more to the list
		                kdc = kdc + 1
		                occupancy(x,y,z)%cDC(0) = kdc
		                occupancy(x,y,z)%cDC(kdc) = idc
			        endif
		        enddo
	        enddo
        enddo
        if (IN_VITRO) then
			call setDCzrange(xdc,ydc)
		endif
    enddo
    call check_DCproximity
endif
write(logmsg,*) 'Number of DCs: ',NDC
call logger(logmsg)
if (USE_DC_COGNATE) then
	ncogDC = DC_COGNATE_FRACTION*NDC
	do idc = ncogDC+1,NDC
		DClist(idc)%cognate = .false.
	enddo
	write(logmsg,*) 'Number of DCs bearing cognate antigen: ',ncogDC
	call logger(logmsg)
endif
if (IN_VITRO) then
	nzlim = 1
	n2Dsites = n2Dsites - NDC*NDCsites
else
	nzlim = NZ
	allocate(permc(nlist))
	do k = 1,nlist
	    permc(k) = k
	enddo
	call permute(permc,nlist,kpar)
endif
istep = 0
ntagged  = 0
ntagged_left = 0
id = 0
cogid = 0
done = .false.
do while (.not.done)
do x = 1,NX
    do y = 1,NY
	    do z = 1,nzlim
	        if (occupancy(x,y,z)%indx(1) == 0) then ! vacant site, not OUTSIDE_TAG or DC
				if (IN_VITRO) then	! populate a specified fraction of the dish area
					if (id == IV_NTCELLS) then
						done = .true.
						exit
					endif
		            R = par_uni(kpar)
		            if (R > iv_fraction) cycle
				endif
                id = id+1
                lastID = id
                site = (/x,y,z/)
                gen = 1
                stage = NAIVE
                region = LYMPHNODE
                call select_cell_type(ctype,cognate,kpar)
		        tag = 0
                cnt(ctype) = cnt(ctype) + 1
                if (evaluate_residence_time) then
	                cognate = .false.
                elseif (calibrate_motility) then
	                cognate = .false.
                    if (taggable(site)) then
						tag = TAGGED_CELL
	                    ntagged = ntagged + 1
					endif
			    elseif (track_DCvisits) then
	                cognate = .false.
					if (DC_CHEMO_NOTRAFFIC) then	! cells are tagged initially
			            R = par_uni(kpar)
			            if (R < DC_CHEMO_FRACTION) then
							tag = TAGGED_CELL
						endif
				    endif
			    else
                    if (cognate) then
                        ncogseed(ctype) = ncogseed(ctype) + 1
                    endif
                endif
                if (IN_VITRO) then
					k = id
				else
	                k = permc(id)
	            endif
                call create_Tcell(k,cellist(k),site,ctype,cognate,gen,tag,stage,region,ok)
                if (.not.ok) return
                occupancy(x,y,z)%indx(1) = k
! Redundant - in create_Tcell()
!                if (track_DCvisits .and. DC_CHEMO_NOTRAFFIC .and. tag == TAGGED_CELL) then
!					R = par_uni(kpar)
!					if (R < HI_CHEMO_FRACTION) then
!						cellist(k)%DCchemo = HI_CHEMO
!					else
!						cellist(k)%DCchemo = LO_CHEMO
!					endif
!				endif
            endif
	    enddo
	    if (done) exit
    enddo
    if (done) exit
enddo
if (.not.IN_VITRO) done = .true.
enddo
if (.not.IN_VITRO) then
	deallocate(permc)
endif
call make_cognate_list(ok)
if (.not.ok) return

!write(*,*) 'nlist,id,RESIDENCE_TIME: ',nlist,id,RESIDENCE_TIME
write(nfout,*) 'nlist,RESIDENCE_TIME: ',nlist,RESIDENCE_TIME
nlist = id	! this is already the case for 3D blob
NTcells = nlist
Nsites = NTcells + NDC*NDCsites	! not relevant for IN_VITRO
NTcells0 = NTcells
Radius0 = Radius
scale_factor = real(NTC_LN)*NLN_RESPONSE/NTcells0
!write(*,*) 'scale_factor: ',scale_factor

if (use_DC .and. .not.IN_VITRO) call initial_binding
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine PlaceHEVPortals(ok)
logical :: ok
real(DP) :: R
integer :: kpar = 0
integer :: ihev, k, x, y, z, site(3)
integer :: nhev_required = 50	! testing

allocate(HEVlist(nhev_required))
do ihev = 1,nhev_required
	k = 0
	do
		k = k+1
		if (k > 10000) then
			call logger('PlaceHEVPortals: too many iterations')
			ok = .false.
			return
		endif
        R = par_uni(kpar)
        x = 1 + R*NX
        R = par_uni(kpar)
        y = 1 + R*NY
        R = par_uni(kpar)
        z = 1 + R*NZ
        site = (/x,y,z/)
		call CheckSite(HEV_SITE,site,ok)
		if (ok) exit
	enddo
	HEVlist(ihev)%site = site
	write(*,*) 'HEV: ',ihev,site
enddo
NHEV = nhev_required
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine assignDCsites(idc,nassigned,ok)
integer :: idc, nassigned
logical :: ok
integer :: k, site(3), site1(3)

ok = .false.
site = DClist(idc)%site
occupancy(site(1),site(2),site(3))%indx = -idc     ! This site holds a DC (centre)
occupancy(site(1),site(2),site(3))%DC(0) = -idc    ! This site holds a DC
nassigned = 1
!write(*,*) 'assignDCsites: ',site
do k = 2,NDCsites
    site1 = site + DCoffset(:,k)
    if (occupancy(site1(1),site1(2),site1(3))%indx(1) == -idc) then
		nassigned = nassigned + 1
		cycle
	endif
    ! Reduce available blob sites only if the site is within the blob.
    ! A DC site should never fall outside, I think.  Check this.
    if (.not.inside_xyz(site1)) then
		write(*,*) 'Error: assignDCsites: site outside grid: '
		write(*,*) 'DC site: ',site
		write(*,*) 'DCoffset: ',k,DCoffset(:,k)
		write(*,*) 'site1: ',site1
		call logger("Error: assignDCsites: site outside grid")
		return
	endif
    if (occupancy(site1(1),site1(2),site1(3))%indx(1) == OUTSIDE_TAG) then
		if (.not.IN_VITRO .or. (IN_VITRO .and. site1(3) == 1)) then
			call logger('assignDCsites: DC hits boundary')
			return
		endif
	endif
    if (IN_VITRO) then
		if (site1(3) == 1 .and. occupancy(site1(1),site1(2),site1(3))%indx(1) /= OUTSIDE_TAG) then
			nassigned = nassigned + 1
		endif
	else
		if (occupancy(site1(1),site1(2),site1(3))%indx(1) /= 0) cycle
		if (occupancy(site1(1),site1(2),site1(3))%indx(2) /= 0) cycle
		nassigned = nassigned + 1
	endif
    occupancy(site1(1),site1(2),site1(3))%indx = -idc     ! This site holds a DC (soma)
    occupancy(site1(1),site1(2),site1(3))%DC(0) = -idc    ! This site holds a DC
enddo
ok = .true.
end subroutine

!-----------------------------------------------------------------------------
! Number of exit portals required for the current cell population, ncells. 
! For SURFACE_PORTALS, linearly proportional to surface area, i.e. to
! ncells^(2/3), factor = exit_fraction
! Modified to use an adjustment that was computed with:
!	Tres = 12
!	chemotaxis factor = 1
!	chemotaxis radius = 100 um
! exit_fraction = mN + c
! where
!	N = ncells
!	m = 1.607E-5
!	c = 0.00602
! 
! Better is a quadratic fit (from Excel, steady_chemo_0_1.xls):
! exit_fraction = a.x^2 + b.x + c
! where x = Ncells/1000
!
! The number of exit portals is based on the residence time for CD4 cells
!-----------------------------------------------------------------------------
integer function requiredExitPortals(ncells)
integer :: ncells
real :: a, b, c, x, Fe
real, parameter :: pow = 2./3.
!real, parameter :: m = 1.607E-8, c = 0.00602
! parameters for chemotaxis
!real, parameter :: a_chemo_12h = -2.228E-08
!real, parameter :: b_chemo_12h = 2.624E-05
!real, parameter :: c_chemo_12h = 5.510E-03

! parameters for no chemotaxis, INLET_R_FRACTION = 1.0
! Best fit:
! y = -3.873E-08x2 + 5.150E-05x + 1.029E-02
!real, parameter :: a_nochemo_24h = -3.873E-08
!real, parameter :: b_nochemo_24h = 5.150E-05
!real, parameter :: c_nochemo_24h = 1.029E-02

! parameters for no chemotaxis, INLET_R_FRACTION = 0.7
! Best fit:
!y = -4.058E-08x2 + 5.590E-05x + 1.006E-02
real, parameter :: a_nochemo_24h = -4.058E-08
real, parameter :: b_nochemo_24h = 5.590E-05
real, parameter :: c_nochemo_24h = 1.006E-02

! This is to adjust Ne when there is general exit chemotaxis, to generate steady-state
real, parameter :: K_Ke = 1.0		! 0.68 for Ke = 1.0, 1.0 for Ke = 0.0 (Pe = 0.02)

! TESTING
!if (use_exit_chemotaxis) then 
	a = a_nochemo_24h
	b = b_nochemo_24h
	c = c_nochemo_24h
! try making it constant
!	a = 0
!	b = 0
! Calibration...
!	c = 12.5E-03*24/residence_time	! for R=23 (50k)
!	c = 15.5E-03*24/residence_time	! for R=29 (100k)
!	c = 19.3E-03*24/residence_time	! for R=36.3 (200k)
!	c = 22.0E-03*24/residence_time	! for R=41.5 (300k)
!	c = 24.5E-03*24/residence_time	! for R=45.7 (400k)
!	c = 26.5E-03*24/residence_time	! for R=49.2 (500k)
!else
!	a = a_nochemo_12h
!	b = b_nochemo_12h
!	c = c_nochemo_12h
!endif
if (TAGGED_LOG_PATHS) then
	requiredExitPortals = 1
elseif (FIXED_NEXITS) then
	x = NTcells0/1000
	exit_fraction = (a*x**2 + b*x + c)*24/residence_time(CD4)
	requiredExitPortals = exit_fraction*NTcells0**pow + 0.5
elseif (RELAX_INLET_EXIT_PROXIMITY) then	! just for ncells = 100k
	x = ncells/1000
	Fe = 2.50E-03*x**3.595E-01		! power law fit (8 points, 51k - 1.1m cells)
	exit_fraction = Fe*24.0/residence_time(CD4)
	requiredExitPortals = K_Ke*exit_fraction*ncells**pow + 0.5 
!	write(logmsg,'(a,f7.4)') 'exit_fraction: ',exit_fraction 
!	call logger(logmsg)
else
	x = ncells/1000
!	Fe = (a*x**2 + b*x + c)
!	Fe = 3.028E-03*x**3.571E-01		! power law fit (7 points)
	Fe = 2.990E-03*x**3.595E-01		! power law fit (8 points, 51k - 1.1m cells)
	exit_fraction = Fe*24.0/residence_time(CD4)
	requiredExitPortals = K_Ke*exit_fraction*ncells**pow + 0.5 
!	write(logmsg,'(a,f7.4)') 'exit_fraction: ',exit_fraction 
!	call logger(logmsg)
endif
end function

!---------------------------------------------------------------------
! Currently, exit locations are distributed randomly through the blob. 
! Sites within the SOI of an exit are labelled with %exitnum, unless
! the site is also within the SOI of another exit, in which case the
! site is labelled with the exit index of the nearest exit, or in
! the case of a tie the choice is made randomly.
! When a site is an exit portal location (centre), %exitnum = -(the exit index)
! Note:
! When the blob grows, more exits will be added, and when the blob
! shrinks again exits must be removed.  When an exit is removed from
! the blob we set exitlist(iexit)%ID = 0, and all sites within the SOI
! must have %exitnum either set to 0 or, if within the SOI of another
! exit, set to the ID of the nearest exit.
! This is analogous to the treatment of DCs
!
! PlaceExitPortals > AddExitPortal > ChoosePortalSite > fPlaceExitPortal
!---------------------------------------------------------------------
subroutine PlaceExitPortals(ok)
integer :: Nex, iexit, site(3)
logical :: ok
logical :: testing = .false.

!if (exit_region /= EXIT_CHEMOTAXIS) then
!    write(*,*) 'Error: placeExits: not EXIT_CHEMOTAXIS'
!    stop
!endif
call logger('PlaceExitPortals')
ok = .true.
if (TAGGED_LOG_PATHS) then
	Nex = 1
    max_exits = 10*Nex
    allocate(exitlist(max_exits))       ! Set the array size to 10* the initial number of exits
    Nexits = 1
    site = Centre
	call PlaceExitPortal(1,site)
	return
endif
if (testing) then
    lastexit = 1
    Nexits = lastexit
    allocate(exitlist(lastexit))
    exitlist(1)%ID = 1
    exitlist(1)%site = (/NX/2,NY/2,NZ/2/)
    return
endif
lastexit = 0
Nexits = 0
if (use_traffic) then
	Nex = requiredExitPortals(NTcells0)
    max_exits = 10*Nex
    write(logmsg,*) 'NTcells0, Nex, max_exits: ',NTcells0,Nex,max_exits
    call logger(logmsg)
    allocate(exitlist(max_exits))       ! Set the array size to 10* the initial number of exits
else
    Nex = 0
    write(logmsg,*) 'No exits'
    call logger(logmsg)
    return
endif
do iexit = 1,Nex
	call AddExitPortal
!	write(logmsg,*) 'Placed exit portal: ',iexit,Nex
!	call logger(logmsg)
enddo
write(logmsg,*) 'Number of exit portals: ',Nex
call logger(logmsg)
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine getExitNum(iexit)
integer :: iexit
integer :: i

iexit = 0
! First look for an unused index
do i = 1,lastexit
	if (exitlist(i)%ID == 0) then
		iexit = i
		exit
	endif
enddo
if (iexit == 0) then
	lastexit = lastexit + 1
	iexit = lastexit
endif
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine AddExitPortal
integer :: site(3)
integer :: iexit

call ChoosePortalSite(site)
call getExitNum(iexit)
Nexits = Nexits + 1
if (lastexit > max_exits) then
	write(logmsg,*) 'Error: AddExitPortal: too many exits: need to increase max_exits: ',max_exits
	call logger(logmsg)
	stop
endif
call PlaceExitPortal(iexit,site)
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
integer function blobNeighbours(x,y,z)
integer :: x,y,z
integer :: nb, k, xx, yy, zz

nb = 0
blobNeighbours = 0
if (x<1 .or. x>NX) return
if (y<1 .or. y>NY) return
if (z<1 .or. z>NZ) return
if (occupancy(x,y,z)%indx(1) < 0) return
do k = 1,27
	if (k == 14) cycle
	xx = x + jumpvec(1,k)
	yy = y + jumpvec(2,k)
	zz = z + jumpvec(3,k)
	if (occupancy(xx,yy,zz)%indx(1) >= 0) nb = nb + 1
enddo
blobNeighbours = nb
end function

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine getBoundarySite(u,site,ok)
real :: u(3)
integer :: site(3)
logical :: ok
integer :: k, x, y, z, nb
real :: r0, r

ok = .false.		
r0 = radius
do k = -10,10
	r = r0 + k*0.5
	x = x0 + u(1)*r
	y = y0 + u(2)*r
	z = z0 + u(3)*r
	nb = blobNeighbours(x,y,z)
	if (nb == 0) then
		exit
	elseif (nb <= 17) then
		ok = .true.		
		exit
	endif
enddo
site = (/x,y,z/)
end subroutine

!-----------------------------------------------------------------------------
! Exit sites (portals) are either on the blob surface, or within the blob.
! If this is the initial placement (istep = 0) then only blob boundary and
! other exit portals need be considered.
! Otherwise HEV and DC locations need to be avoided as well.
!-----------------------------------------------------------------------------
subroutine ChoosePortalSite(site)
integer :: site(3)
integer :: xex, yex, zex
real(DP) :: R
real :: prox, u(3), theta, rx
integer :: kpar = 0
logical :: ok

!call logger('ChoosePortalSite')
if (SURFACE_PORTALS) then
	! randomly choose direction in 3D, then locate sites near this vector on the blob boundary
	! Need to restrict the sinus interface by removing the follicle interface.
	! Assume that the follicle interface is a cap from (normalized) x = XFOLLICLE (< 1)
	do
		R = par_uni(kpar)
		u(1) = 2*R - 1
		if (u(1) > XFOLLICLE) cycle
		rx = sqrt(1-u(1)*u(1))
		R = par_uni(kpar)
		theta = 2*PI*R
		u(2) = rx*cos(theta)
		u(3) = rx*sin(theta)
		call getBoundarySite(u,site,ok)
		if (.not.ok) cycle
!		xex = site(1)
!		yex = site(2)
!		zex = site(3)
		call CheckSite(EXIT_SITE,site,ok)
		if (.not.ok) cycle
!		if (use_DC) then
!			if (tooNearDC(site,exit_DCprox)) then		! exit_DCprox is min distance in sites
!				call logger('tooNearDC')
!				cycle
!			endif
!		endif
!		prox = exit_prox*chemo_N				! chemo_N is chemo_radius in units of sites
!		prox = 0.5*prox
!!		if (USE_PORTAL_EGRESS) then
!!			prox = 0.3*prox
!!		endif
!		if (tooNearExit(site,prox)) then	! too near another exit
!!			call logger('tooNearExit')
!			cycle
!		endif	
		exit
	enddo
else	! blob portals
	do
		R = par_uni(kpar)
		xex = 1 + R*NX
		R = par_uni(kpar)
		yex = 1 + R*NY
		R = par_uni(kpar)
		zex = 1 + R*NZ
		site = (/xex,yex,zex/)
		if (occupancy(xex,yex,zex)%indx(1) < 0) cycle     ! OUTSIDE_TAG or DC
		if (use_DC) then
	!        write(*,*) 'use_DC?'
	!        stop
			if (tooNearDC(site,exit_DCprox)) cycle    ! exit_DCprox is min distance in sites
		endif
		prox = exit_prox*chemo_N				! chemo_N is chemo_radius in units of sites
		if (tooNearExit(site,prox)) cycle       ! too near another exit
		prox = 0.5*exit_prox*chemo_N
		if (tooNearBdry(site,prox)) cycle       ! too near the blob boundary
		exit
	enddo
endif
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine PlaceExitPortal(iexit,site)
integer :: iexit, site(3)
integer :: kexit, x, y, z, x2, y2, z2, ek(3), vk(3), ns
integer :: xex, yex, zex, xmin, xmax, ymin, ymax, zmin, zmax
real(DP) :: R
real :: d2, d2k
integer :: kpar = 0

if (iexit == lastexit + 1) then
	lastexit = iexit
endif
ns = 0
exitlist(iexit)%ID = iexit      ! when an exit site is lost (because the blob retracted) set %ID = 0
exitlist(iexit)%site = site
xex = site(1)
yex = site(2)
zex = site(3)
occupancy(xex,yex,zex)%exitnum = -iexit     ! This site holds an exit 
!write(logmsg,*) 'PlaceExitPortal: site,exitnum: ',xex,yex,zex,-iexit
!call logger(logmsg)

xmin = xex - chemo_N
xmax = xex + chemo_N
ymin = yex - chemo_N
ymax = yex + chemo_N
zmin = zex - chemo_N
zmax = zex + chemo_N
xmin = max(1,xmin)
xmax = min(NX,xmax)
ymin = max(1,ymin)
ymax = min(NY,ymax)
zmin = max(1,zmin)
zmax = min(NZ,zmax)
do x = xmin,xmax
    x2 = (x-xex)*(x-xex)
    do y = ymin,ymax
        y2 = (y-yex)*(y-yex)
        do z = zmin,zmax
	        z2 = (z-zex)*(z-zex)
	        d2 = x2 + y2 + z2
	        if (d2 <= chemo_N*chemo_N) then
		        if (d2 > 0) then   ! don't touch exit site
     ! NOTE!!!  If this site is already marked as within another exit's SOI, we need to
     ! determine which exit is closer, and set %exitnum to this exit's ID
                    kexit = occupancy(x,y,z)%exitnum
                    if (kexit > 0) then
                        ek = exitlist(kexit)%site
                        vk = (/x,y,z/) - ek
                        d2k = dot_product(vk,vk)    ! distance from the other exit (kexit)
!                            write(*,*) 'other exit: ',iexit,kexit,d2,d2k
                        if (d2k < d2) then
!	                            write(*,*) 'closer exit: ',iexit,kexit,d2,d2k
                            cycle
                        endif
                        if (d2k == d2) then     ! toss a coin
                            R = par_uni(kpar)
                            if (R < 0.5) cycle
                        endif
                        occupancy(x,y,z)%exitnum = iexit  ! the site is closer to iexit than to kexit
                        ns = ns + 1
                    elseif (kexit == 0) then
                        occupancy(x,y,z)%exitnum = iexit    ! this site is closest to exit iexit
                        ns = ns + 1
                    endif
                endif
	        endif
        enddo
    enddo
enddo
!write(*,*) 'near sites: ',iexit,ns
end subroutine

!---------------------------------------------------------------------
! Remove the last nremex exits in the list.
! Perhaps it would make more sense to remove exits that are close to 
! others, or to adjust the locations of remaining exits.
! NOT USED
!---------------------------------------------------------------------
subroutine removeExits1(nremex)
integer :: nremex
integer :: Nex, k, iexit, site(3)

!write(*,*) 'removeExits: ',nremex
Nex = lastexit
k = 0
do iexit = Nex,1,-1
    if (exitlist(iexit)%ID == 0) cycle
    k = k+1
    site = exitlist(iexit)%site
!    write(logmsg,'(a,4i6)') 'removeExits: removeExitPortal: ',iexit,site
!    call logger(logmsg)
    call RemoveExitPortal(site)
!    write(logmsg,'(a,4i6)') 'did removeExitPortal'
!    call logger(logmsg)
    if (k == nremex) exit
enddo
!call checkExits("in removeExits")
end subroutine

!---------------------------------------------------------------------
! Remove exits that are close to others.
!---------------------------------------------------------------------
subroutine RemoveExitPortals(nremex)
integer :: nremex
integer :: Nex, k, iexit, kexit, site1(3), site2(3)
real :: r(3), sum, d2, cmin
real, allocatable :: closeness(:)

!write(*,*) 'RemoveExitPortals: ',nremex
Nex = lastexit
allocate(closeness(Nex))
closeness = 0
do iexit = 1,Nex
	sum = 0
    if (exitlist(iexit)%ID == 0) cycle
    site1 = exitlist(iexit)%site
    do kexit = 1,Nex
		if (exitlist(kexit)%ID == 0) cycle
		if (iexit == kexit) cycle
	    site2 = exitlist(kexit)%site
	    r = site1 - site2
	    d2 = norm2(r)
	    sum = sum + 1/d2
	enddo
	closeness(iexit) = sum
enddo

do k = 1,nremex
	cmin = 1.0e10
	do kexit = 1,Nex
		if (closeness(kexit) == 0) cycle
		if (closeness(kexit) < cmin) then
			cmin = closeness(kexit)
			iexit = kexit
		endif
	enddo
	closeness(iexit) = 0
    site1 = exitlist(iexit)%site
    call RemoveExitPortal(site1)
enddo
!call checkExits("in RemoveExitPortals")
deallocate(closeness)
end subroutine

!---------------------------------------------------------------------
! Remove the exit portal at site.
!---------------------------------------------------------------------
subroutine RemoveExitPortal(site)
integer :: site(3)
integer :: iexit,xex,yex,zex,xmin,xmax,ymin,ymax,zmin,zmax,x,y,z,ek(3),vk(3),k,kmin,i
real :: d2k, d2min

xex = site(1)
yex = site(2)
zex = site(3)
iexit = -occupancy(xex,yex,zex)%exitnum
if (iexit <= 0) then
	write(logmsg,*) 'Error: RemoveExitPortal: no exit at: ',site,iexit
	call logger(logmsg)
	do i = 1,lastexit
		site = exitlist(i)%site
		write(logmsg,*) 'exit: ',i,site,occupancy(site(1),site(2),site(3))%exitnum
		call logger(logmsg)
	enddo
	stop
endif
exitlist(iexit)%ID = 0
occupancy(xex,yex,zex)%exitnum = iexit      ! to ensure that the site is processed in the next section
											! effectively it is treated like a site within the portal's SOI
!write(*,*) 'RemoveExitPortal: site,exitnum: ',site,iexit

xmin = xex - chemo_N
xmax = xex + chemo_N
ymin = yex - chemo_N
ymax = yex + chemo_N
zmin = zex - chemo_N
zmax = zex + chemo_N
xmin = max(1,xmin)
xmax = min(NX,xmax)
ymin = max(1,ymin)
ymax = min(NY,ymax)
zmin = max(1,zmin)
zmax = min(NZ,zmax)
do x = xmin,xmax
    do y = ymin,ymax
        do z = zmin,zmax
            if (occupancy(x,y,z)%exitnum == iexit) then
                occupancy(x,y,z)%exitnum = 0
                ! at this point we need to check to see if (x,y,z) is within the SOI of another exit
                d2min = 9999
                kmin = 0
                do k = 1,lastexit
                    if (exitlist(k)%ID == 0) cycle
                    ek = exitlist(k)%site
                    vk = (/x,y,z/) - ek
                    d2k = dot_product(vk,vk)    ! distance from the other exit (kexit)
		            if (d2k <= chemo_N*chemo_N) then
		                if (d2k < d2min) then
		                    d2min = d2k
		                    kmin = k
		                endif
		            endif
		        enddo
		        if (kmin > 0) then
		            occupancy(x,y,z)%exitnum = kmin
!					write(*,*) '  RemoveExitPortal: site,exitnum: ',x,y,z,kmin
		        endif
            endif
        enddo
    enddo
enddo
if (iexit == lastexit) then
	lastexit = lastexit - 1
endif
Nexits = Nexits - 1
end subroutine

!---------------------------------------------------------------------
! Move the exit portal to a site closer to the blob centre.
!---------------------------------------------------------------------
subroutine MoveExitPortalInwards(site0)
integer :: site0(3)
integer :: site(3), site1(3), iexit0, iexit1, dx, dy, dz
real :: d0, dc, d2, d2min, r(3)
real, parameter :: dlim = 1.95

iexit0 = -occupancy(site0(1),site0(2),site0(3))%exitnum
! First, remove the exit
call RemoveExitPortal(site0)
!write(logmsg,*) 'Removed exit within moveExitPortalInwards: ',iexit0,site0
!call logger(logmsg)
!call checkexits("in moveExitPortalInwards")
! Now add the exit back at a closer site.  Choose a new site at
! least a distance of dlim closer to blob centre.  Choose the
! site closest to site0 from among the neighbours within the
! specified distance from the centre.
! NOTE: need to check that there is not a DC or another exit portal nearby.  
! If there is then a completely new portal site will need to be chosen, 
! using ChoosePortalSite()
d0 = cdistance(site0)
d2min = 1.0e10
site1 = 0
do dx = -3,3
	do dy = -3,3
		do dz = -3,3
			site = site0 + (/dx,dy,dz/)
			dc = cdistance(site)
			if (d0-dc >= dlim) then
				r = site - site0
				d2 = norm2(r)
				if (d2 < d2min) then
					d2min = d2
					site1 = site
				endif
			endif
		enddo
	enddo
enddo
if (site1(1) == 0) then
	write(logmsg,*) 'Error: moveExitPortalInwards: no site to move to: ',iexit0,site0
	call logger(logmsg)
	stop
endif
call getExitNum(iexit1)
!write(logmsg,*) 'moveExitPortalInwards: ',iexit0,' from: ',site0
!call logger(logmsg)
!write(logmsg,*) '                new #: ',iexit1,'   at: ',site1
!call logger(logmsg)
Nexits = Nexits + 1
call PlaceExitPortal(iexit1,site1)
if (exitlist(iexit1)%site(1) /= site1(1) .or. &
	exitlist(iexit1)%site(2) /= site1(2) .or. &
	exitlist(iexit1)%site(3) /= site1(3)) then
	write(logmsg,*) 'Error: moveExitPortalInwards: bad site: ',iexit1,exitlist(iexit1)%site,site1
	call logger(logmsg)
	stop
endif
!write(logmsg,*) 'now: ',iexit1,' is at: ',exitlist(iexit1)%site
!call logger(logmsg)
!write(logmsg,*) 'and exit #10 is at: ',exitlist(10)%site
!call logger(logmsg)
!write(logmsg,'(a,6i6)') 'moveExitPortalInwards: PlaceExitPortal: ',istep,iexit,site1,Nexits
!call logger(logmsg)
end subroutine

!---------------------------------------------------------------------
! The exit portal iexit needs to be moved by one site in the direction 
! closest to that given by the unit vector v(:).
!---------------------------------------------------------------------
subroutine MoveExitPortal(iexit0,v)
integer :: iexit0
real :: v(3)
integer :: iexit1, site0(3), site1(3), k, kmax
real :: proj, pmax, jump(3)

pmax = 0
do k = 1,27
	if (k == 14) cycle
	jump = jumpvec(:,k)
	proj = dot_product(jump,v)/norm(jump)
	if (proj > pmax) then
		pmax = proj
		kmax = k
	endif
enddo
site0 = exitlist(iexit0)%site
site1 = site0 + jumpvec(:,kmax)
write(logmsg,'(a,i4,a,3i4,a,3i4)') 'moveExitPortal: ',iexit0,' from: ',site0,' to: ',site1
call logger(logmsg)
call RemoveExitPortal(site0)
call getExitNum(iexit1)
Nexits = Nexits + 1
call PlaceExitPortal(iexit1,site1)
end subroutine

!---------------------------------------------------------------------
! Check the suitability of site(:) for locating sitetype:
! EXIT_SITE, HEV_SITE, or DC_SITE
!---------------------------------------------------------------------
subroutine CheckSite(sitetype,site,ok)
integer :: sitetype, site(3)
logical :: ok
real :: proxlimit, r

select case(sitetype)
case(EXIT_SITE)
	proxlimit = proximity_limit(EXIT_SITE,EXIT_SITE)
	if (tooNearExit(site,proxlimit)) then	! too near another exit
		ok = .false.
		return
	endif
	proxlimit = proximity_limit(EXIT_SITE,DC_SITE)
	if (tooNearDC(site,proxlimit)) then	! too near a DC
		ok = .false.
		return
	endif
	proxlimit = proximity_limit(EXIT_SITE,HEV_SITE)
	if (tooNearHEV(site,proxlimit)) then	! too near an HEV
		ok = .false.
		return
	endif
case(HEV_SITE)
	r = cdistance(site)
	if (r/radius < R_limit(HEV_SITE,1) .or. r/radius > R_limit(HEV_SITE,2)) then
		ok = .false.
		return
	endif
   proxlimit = proximity_limit(HEV_SITE,HEV_SITE)
    if (tooNearHEV(site,proxlimit)) then
		ok = .false.
		return
	endif
    proxlimit = proximity_limit(HEV_SITE,EXIT_SITE)
    if (tooNearExit(site,proxlimit))  then
		ok = .false.
		return
	endif
    proxlimit = proximity_limit(HEV_SITE,BDRY_SITE)
    if (tooNearBdry(site,proxlimit))  then
		ok = .false.
		return
	endif
    proxlimit = proximity_limit(HEV_SITE,DC_SITE)
	if (tooNearDC(site,proxlimit)) then	! too near an HEV
		ok = .false.
		return
	endif
case(DC_SITE)
	r = cdistance(site)
	if (r/radius < R_limit(DC_SITE,1) .or. r/radius > R_limit(DC_SITE,2)) then
		ok = .false.
		return
	endif
    proxlimit = proximity_limit(DC_SITE,DC_SITE)
    if (tooNearDC(site,proxlimit)) then
		ok = .false.
		return
	endif
    proxlimit = proximity_limit(DC_SITE,EXIT_SITE)
    if (tooNearExit(site,proxlimit))  then
		ok = .false.
		return
	endif
    proxlimit = proximity_limit(DC_SITE,BDRY_SITE)
    if (tooNearBdry(site,proxlimit))  then
		ok = .false.
		return
	endif
    proxlimit = proximity_limit(DC_SITE,HEV_SITE)
	if (tooNearHEV(site,proxlimit)) then	! too near an HEV
		ok = .false.
		return
	endif
end select
ok = .true.
end subroutine

!---------------------------------------------------------------------
! Display exit portal location info
!---------------------------------------------------------------------
subroutine displayExits
integer :: Nex, k, iexit, kexit, site1(3), site2(3), count
real :: r(3), sum, d2

Nex = lastexit
do iexit = 1,Nex
	sum = 0
    if (exitlist(iexit)%ID == 0) cycle
    site1 = exitlist(iexit)%site
    do kexit = 1,Nex
		if (exitlist(kexit)%ID == 0) cycle
		if (iexit == kexit) cycle
	    site2 = exitlist(kexit)%site
	    r = site1 - site2
	    d2 = norm2(r)
	    sum = sum + 1/d2
	enddo
	count = neighbourhoodCount(site1)
	write(logmsg,'(5i4,3f10.3)') iexit,site1,count,cdistance(site1),(PI/180)*atan2(site1(2)-Centre(2),site1(3)-Centre(3)),1000*sum
	call logger(logmsg)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! need to ensure that exits do not get too close together
!-----------------------------------------------------------------------------------------
subroutine CheckExitSpacing(msg)
character*(*) :: msg
integer :: iexit1, iexit2, site1(3), site2(3), im1=0, im2=0
real :: r(3), d, dmin, f(3), df(3)
logical :: ok = .true.
real :: dlim

dlim = exit_prox*chemo_radius
dmin = 1.0e10
do iexit1 = 1,lastexit
	if (exitlist(iexit1)%ID == 0) cycle
	site1 = exitlist(iexit1)%site
	f = 0
	do iexit2 = 1,lastexit
		if (iexit1 == iexit2) cycle
		if (exitlist(iexit2)%ID == 0) cycle
		site2 = exitlist(iexit2)%site
		r = site1 - site2
		d = norm(r)
		if (d < dmin) then
			dmin = d
			im1 = iexit1
			im2 = iexit2
		endif
		if (d < dlim) then
			write(logmsg,'(2i4,f6.2)') iexit1,iexit2,d
			call logger(logmsg)
			df = r/d		! unit vector in direction of force
			df = df/d		! scale by inverse distance
			f = f + df
		endif
	enddo
	if (f(1) /= 0 .or. f(2) /= 0 .or. f(3) /= 0) then
		f = f/norm(f)
		call MoveExitPortal(iexit1,f)
	endif
enddo
write(logmsg,*) msg,': Min exit spacing: ',im1,im2,dmin
call logger(logmsg)
end subroutine


!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine init_cytokine
integer :: x, y, z, icyt, iseq, tag, offset
integer :: np(MAX_CYT)
real :: dt, dx, mols_pM

Ncytokines = N_CYT
cyt_seq = 0
do iseq = 1,Ncytokines
    tag = cyt_tag(iseq)
    cyt_seq(tag) = iseq
enddo
do iseq = 1,Ncytokines
    tag = cyt_tag(iseq)
    select case(tag)
    case (IL2_TAG)
        np(iseq) = IL2_NP
    case (IL4_TAG)
        np(iseq) = IL4_NP
    case (IL7_TAG)
        np(iseq) = IL7_NP
    case (IL9_TAG)
        np(iseq) = IL9_NP
    case (IL15_TAG)
        np(iseq) = IL15_NP
    case (IL21_TAG)
        np(iseq) = IL21_NP
    end select
enddo
offset = 0
do iseq = 1,Ncytokines
    NP_offset(iseq) = offset
    offset = offset + np(iseq)
    NP_offset(iseq+1) = offset
enddo

allocate(cyt_constit(Ncytokines))
if (calibrate_diffusion) then
    cyt_constit = 0
else
    cyt_constit = 0
    iseq = cyt_seq(IL2_TAG)
    cyt_constit(iseq) = IL2_CONSTIT_RATE    ! units molecules/um^3/min
endif

! For testing
K_diff = 2.0E-12                ! m^2/s
dt = DELTA_T*60.0/NDIFFSTEPS    ! s
dx = DELTA_X*1.0E-6             ! m
do icyt = 1,Ncytokines
    delta_diff(icyt) = K_diff(icyt)*dt/(dx*dx)
    if (6*delta_diff(icyt) >= 1.0) then     ! was 6*
        write(*,*) 'Bad delta_diff: need to increase NDIFFSTEPS: ',delta_diff(icyt)
!        stop
    endif
    if (calibrate_diffusion) then
        cyt_init(icyt) = 1.0
    else
        cyt_init(icyt) = cyt0(icyt)
    endif
enddo

allocate(cyt(NX,NY,NZ,Ncytokines))
allocate(cyt_mols(Ncytokines))
allocate(dcyt_mols(Ncytokines))
dcyt_mols = 0
write(*,*) 'Initial cytokine conc., mols'
do icyt = 1,Ncytokines
    if (calibrate_diffusion) then
        do x = 1,NX/2
            cyt(x,:,:,icyt) = cyt_init(icyt)
        enddo
        do x = NX/2,NX
            cyt(x,:,:,icyt) = 0
        enddo
    else
        do z = 1,NZ
            do y = 1,NY
                do x = 1,NX
                    if (occupancy(x,y,z)%indx(1) >= 0) then
                        cyt(x,y,z,icyt) = cyt_init(icyt)
                    else
                        cyt(x,y,z,icyt) = 0
                    endif
                enddo
            enddo
        enddo
        ! need to make # of mols correspond to cyt_init
        mols_pM = L_um3*M_pM/(NTcells*Vc*Navo)
        cyt_mols(icyt) = cyt_init(icyt)/mols_pM   ! -> total number of molecules
        write(*,'(a,f8.3,f12.0)') cyt_name(cyt_tag(icyt)),cyt_init(icyt),cyt_mols(icyt)
    endif
enddo

end subroutine

!--------------------------------------------------------------------------------------
! Motility behaviour
! There are potentially 4 motility states:
! (1) Naive cell, very short non-cognate DC contacts 3 min (motility level 1)
! (2) Short cognate DC contacts 11 min (motility level 2)
! (3) Long cognate DC contacts > 1 hr clusters (motility level 3)
! (4) Short cognate DC contacts 18 min  swarms (motility level 4)
! It isn't clear what level of motility occurs at each stage. Possibly
! the contact duration is related to the motility, since reduced motility
! leads to a higher probability of rebinding.
! A simple way to vary motility is to keep the persistence parameter rho fixed
! and vary alpha (plongjump), the probability of 2 steps, and possibly beta,
! the probability of moving at all.  In any case the jump parameters for each case
! are stored as pjump(2) and dirprob(0:6).
! In principle: speed & Cm -> rho, beta, delta -> dirprob(), pjump().
!
! Mark Miller says that their motility measurements didn't exclude periods when
! T cells were in contact with DC, therefore we have no info about motility in
! various stages of activation.  For now it is safest to use the same motility
! parameters for all stages.
!--------------------------------------------------------------------------------------

!---------------------------------------------------------------------
! TCR stimulation rate depends on the DC's pMHC count and the T cell's
! TCR avidity for the antigen.
! New formulation of TCR stimulation and bind-time.
! The DC antigen density is a measure of the number of pMHC complexes on the cell surface
! For now, this is set initially proportional to the concentration of peptide (pM) used to "pulse" DCs
! Henrickson2008 measured 2.5x10^4 = 25000 pMHC/DC after pulsing 5T33 DCs for 3 hr with 10 uM M- or C-peptide
! The half-life of the pMHC complex was estimated at 6.01 hrs for M-peptide and 2.36 hrs for C-peptide
! Therefore the initial pMHC counts after 18 hours were 25000/8 = 3125 for M-peptide and 130 for C-peptide
! Since 100 pM of M-peptide induced only about 10% proliferation, while 200 pM induced about 80%, it seems
! that the threshold for TCR signalling is about 30 pMHC (assuming that the initial pMHC count is linear
! with antigen concentration).
! Base stimulation rate r is a Hill function of a = pMHC*avidity, max value 1
! Actual rate of change of S is this rate scaled by TC_STIM_RATE_CONSTANT
! Duration of binding is also determined from r.  Currently the specified bind time is
! treated as the maximum, and the actual value depends linearly on r up to this max.
!---------------------------------------------------------------------
real function stimulation_rate(pMHC,avidity)
real :: pMHC, avidity
real :: a
a = (pMHC - DC_pMHC_THRESHOLD)*avidity
stimulation_rate = a**N_STIM/(a**N_STIM + DC_STIM_THRESHOLD**N_STIM)
end function

!---------------------------------------------------------------------
! Updates T cell state, and DC %stimulation (for use in modifying
! %density in update_DCstate).
!---------------------------------------------------------------------
subroutine updater(ok)
logical :: ok
integer :: kcell, ctype, stype, region, iseq, tag, kfrom, kto, k, ncog, ntot
integer :: site(3), site2(3), freeslot, indx(2), status, DC(2), idc
real :: C(N_CYT), mrate(N_CYT), tnow, dstim, S, cyt_conc, mols_pM, Ctemp, dstimrate, stimrate
logical :: divide_flag, producing, first, dbg, unbound, flag, flag1
!logical, save :: first = .true.
type (cog_type), pointer :: p

!write(*,*) 'updater: ',me
ok = .true.
dbg = .false.
flag = .false.
flag1 = .false.
ntot = 0
ncog = 0
tnow = istep*DELTA_T
! Scaling factor to convert total number of molecules in the region to conc in pM
mols_pM = L_um3*M_pM/(NTcells*Vc*Navo)

do kcell = 1,nlist
    if (cellist(kcell)%ID == 0) cycle
    ntot = ntot + 1
    if (dbg) write(*,*) 'kcell: ',kcell
    ctype = cellist(kcell)%ctype
!    stype = struct_type(ctype)
	if (associated(cellist(kcell)%cptr)) then
		stype = COG_TYPE_TAG
	else
		stype = NONCOG_TYPE_TAG
	endif
    if (stype /= COG_TYPE_TAG) cycle
    ! Only cognate cells considered
    ncog = ncog + 1
    p => cellist(kcell)%cptr
    if (.not.associated(p)) then
        write(*,*) 'ERROR: updater: p not associated: ',kcell
        stop
    endif
    call get_region(p,region)
	! Cell death
    if (tnow > p%dietime) then
        call Tcell_death(kcell)
        cycle
    endif
	! TCR stimulation decay
    p%stimulation = p%stimulation*(1 - TCRdecayrate*DELTA_T)

	if (region == LYMPHNODE) then
	! TCR stimulation
		unbound = .false.
		DC = cellist(kcell)%DCbound
		stimrate = 0
		do k = 1,MAX_DC_BIND
			idc = DC(k)
			if (idc /= 0) then
				if (DClist(idc)%capable) then
	!                dstimrate = TC_STIM_RATE_CONSTANT*DClist(idc)%density*p%avidity
					dstimrate = TC_STIM_RATE_CONSTANT*stimulation_rate(DClist(idc)%density,p%avidity)
					stimrate = stimrate + dstimrate
					dstim = dstimrate*DELTA_T
					p%stimulation = p%stimulation + dstim
					p%stimulation = min(p%stimulation, STIMULATION_LIMIT)
					DClist(idc)%stimulation = DClist(idc)%stimulation + dstim
				else    ! unbind T cell from incapable DC
					idc = cellist(kcell)%DCbound(k)
					DClist(idc)%nbound = DClist(idc)%nbound - 1
					DClist(idc)%ncogbound = DClist(idc)%ncogbound - 1
                    call RemoveCogbound(idc,kcell)
					cellist(kcell)%DCbound(k) = 0
					unbound = .true.
				endif
			endif
		enddo
		p%stimrate = stimrate
		if (unbound .and. cellist(kcell)%DCbound(1) == 0 .and. cellist(kcell)%DCbound(2) /= 0) then
			cellist(kcell)%DCbound(1) = cellist(kcell)%DCbound(2)
			cellist(kcell)%DCbound(2) = 0
			cellist(kcell)%unbindtime(1) = cellist(kcell)%unbindtime(2)
			cellist(kcell)%unbindtime(2) = tnow
		endif

		site = cellist(kcell)%site
		if (use_cytokines) then
			! IL receptor stimulation
			status = p%status
			if (use_diffusion) then
				C = cyt(site(1),site(2),site(3),:)
			endif
			S = p%stimulation
			do iseq = 1,Ncytokines
				tag = cyt_tag(iseq)
				kfrom = NP_offset(iseq)+1
				kto = NP_offset(iseq+1)
				if (use_diffusion) then
					cyt_conc = C(iseq)
				else
					cyt_conc = cyt_mols(iseq)*mols_pM   ! -> conc in pM
				endif
				ctemp = cyt_conc
				select case(tag)
				case(IL2_TAG)
					producing = IL2_production_status(p,tnow)
					if (p%IL_state(kfrom) == 0) then    ! temporary measure to detect first call of IL2_update
						first = .true.
					else
						first = .false.
					endif
					call IL2_update(p%cogID,ctype,tnow,S,p%IL_state(kfrom:kto),p%IL_statep(kfrom:kto), &
						first,producing,cyt_conc,Vc,DELTA_T,mrate(iseq),flag1)
				case(IL4_TAG)
					call IL4_update(p%IL_state(kfrom:kto))
				case(IL7_TAG)
					call IL7_update(p%IL_state(kfrom:kto))
				case(IL9_TAG)
					call IL9_update(p%IL_state(kfrom:kto))
				case(IL15_TAG)
					call IL15_update(p%IL_state(kfrom:kto))
				case(IL21_TAG)
					call IL21_update(p%IL_state(kfrom:kto))
				end select

				if (use_diffusion) then
	! Concentration units
	! mrate = mass rate of flow in molecules/min
	! mrate*L_um3/Vc = molecules/L/min
	! mrate*L_um3*DELTA_T/Vc = molecules/L
	! mrate*L_um3*DELTA_T/Vc/Navo = moles/L = M
	! mrate*L_um3*DELTA_T*M_pM/Vc/Navo = pM
	! To convert total number of molecules in the region to conc in pM
	! mols * L_um3*M_pM/(NTcells*Vc*Navo)
					C(iseq) = cyt_conc
	!                C(iseq) = C(iseq) + (mrate(iseq)/Vc)*DELTA_T*M_pM*L_um3/Navo    ! Vc/L_um3 ->free vol in L
					if (C(iseq) < 0) then
						write(*,'(a,6i6,4f8.4)') 'WARNING: cyt < 0: ',kcell,p%cogID,iseq,site,Ctemp,C(iseq)
						C(iseq) = 0
					endif
				else
					! increment total number of molecules of this cytokine
					dcyt_mols(iseq) = dcyt_mols(iseq) + (cyt_conc - ctemp)/(mols_pM*NTcells)
	!                write(*,*) cyt_mols(iseq),dcyt_mols(iseq),ctemp,cyt_conc
				endif
			enddo
			if (use_diffusion) then
				cyt(site(1),site(2),site(3),:) = C
			endif
		endif

		if (use_S1P) then
			! Note that stimrate is normalized with TC_STIM_RATE_CONSTANT for consistency
			! with the calibration of S1P1/CD69 dynamics, which was carried out with
			! normalized stimulation of approx. 1  This is easier than adjusting K1_CD69
			call S1P1_update(p%CD69,p%S1P1,p%stimrate/TC_STIM_RATE_CONSTANT,DELTA_T)
	!        call S1P1_update(p%CD69,p%S1P1,p%stimrate,DELTA_T)
		endif
	endif

! Stage transition
    call updatestage(kcell, tnow, divide_flag)

! Cell division
    if (divide_flag) then
		if (region == LYMPHNODE) then
			indx = occupancy(site(1),site(2),site(3))%indx
			freeslot = 0
			if (indx(1) == 0) then
				site2 = site
				freeslot = 1
			elseif (indx(2) == 0) then
				site2 = site
				freeslot = 2
			else
				call get_free_slot(occupancy,NX,site,site2,freeslot)
			endif
		else
			site2 = 0
			freeslot = -1
!			call logger("T cell divides in the periphery")
		endif
        if (freeslot /= 0) then     ! there is a free slot at site2 (which may be = site)
            call cell_division(kcell,site2,freeslot,ok)
            if (.not.ok) return
        endif
    endif

    if (flag) then
        call show_cognate_cell(kcell)
    endif

enddo
call UpdateHelp
end subroutine

!--------------------------------------------------------------------------------
! The idea is to record the time spent by a cognate CD8 cell bound to a DC at the
! same time as a cognate activated CD4 cell.
! At this time it is not clear how or when help is effective.  For a start we can
! record total timesteps.
!--------------------------------------------------------------------------------
subroutine UpdateHelp
integer :: idc, ncb, k1, k2, n1, n2, kcell1, kcell2
integer :: ctype1, ctype2
logical :: help

do idc = 1,NDC
	ncb = DClist(idc)%ncogbound
	if (ncb == 0) cycle
	n1 = 0
	do k1 = 1,MAX_COG_BIND
		kcell1 = DClist(idc)%cogbound(k1)
		if (kcell1 == 0) cycle
		n1 = n1+1
		ctype1 = cellist(kcell1)%ctype
		help = .false.
		n2 = 0
		do k2 = 1,MAX_COG_BIND
			kcell2 = DClist(idc)%cogbound(k2)
			if (kcell2 == 0) cycle
			n2 = n2+1
			ctype2 = cellist(kcell2)%ctype
			if (ctype1 == CD4 .and. ctype2 == CD8) then
				help = .true.
				exit
			elseif (ctype1 == CD8 .and. ctype2 == CD4) then
				help = .true.
				exit
			endif
			if (n2 == ncb) exit
		enddo
		if (help) then
			cellist(kcell1)%cptr%cnt(1) = cellist(kcell1)%cptr%cnt(1) + 1
		endif
		if (n1 == ncb) exit
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
! Note: USE_STAGETIME(icstage) means use stagetime for the transition from icstage.
! p%stagetime = time that the cell is expected to make transition to the next stage.
!--------------------------------------------------------------------------------
subroutine updatestage1(kcell,tnow,divide_flag)
integer :: kcell
logical :: divide_flag
real :: tnow
integer :: stage, region, gen, ctype
real :: stagetime
type(cog_type), pointer :: p

divide_flag = .false.
p => cellist(kcell)%cptr
!stage = get_stage(p)
call get_stage(p,stage,region)
if (stage >= CLUSTERS) then
    if (.not.cansurvive(p)) then
        write(logmsg,*) 'cell IL2 store too low: ',kcell,p%cogID
        call logger(logmsg)
        p%dietime = tnow
        stage = FINISHED
        call set_stage(p,stage)
    endif
endif
if (stage == FINISHED) return
ctype = cellist(kcell)%ctype
stagetime = p%stagetime
if (tnow > stagetime) then		! time constraint to move to next stage is met
	! May be possible to make the transition from stage
	select case(stage)
	case (NAIVE)        ! possible transition from NAIVE to TRANSIENT
	    if (p%stimulation > 0) then
	        gen = get_generation(p)
            call set_stage(p,TRANSIENT)
            p%dietime = tnow + TClifetime(p)
		    if (USE_STAGETIME(TRANSIENT)) then
			    p%stagetime = tnow + get_stagetime(p,ctype)
		    else
			    p%stagetime = 0		! time is not criterion for next transition
            endif
        endif
	case (TRANSIENT)    ! possible transition from TRANSIENT to CLUSTERS
	    if (reached_IL2_threshold(p)) then
	        nIL2thresh = nIL2thresh + 1
	        tIL2thresh = tIL2thresh + (tnow - cellist(kcell)%entrytime)
!	        write(*,'(a,i6,2f8.1)') '========= Reached IL2 threshold: ',kcell,(tnow - cellist(kcell)%entrytime),p%stimulation
            call set_stage(p,CLUSTERS)
		    if (USE_STAGETIME(CLUSTERS)) then
			    p%stagetime = tnow + get_stagetime(p,ctype)
		    else
			    p%stagetime = 0		! time is not criterion for next transition
		    endif
        endif
	case (CLUSTERS)     ! possible transition from CLUSTERS to SWARMS
	    if (reached_act_threshold(p)) then
            call set_stage(p,SWARMS)
		    if (USE_STAGETIME(SWARMS)) then     ! Must use stagetime(SWARMS) to get to ACTIVATED
			    p%stagetime = tnow + get_stagetime(p,ctype)
		    else
			    p%stagetime = 0
		    endif
        endif
! Code commented because we are now using a revised staging without ACTIVATED
!	case (SWARMS)       ! possible transition from SWARMS to ACTIVATED
!        call set_stage(p,ACTIVATED)
!		if (USE_STAGETIME(ACTIVATED)) then      ! Must use stagetime(ACTIVATED) to get to ACTIVATED
!            gen = get_generation(p)
!            ctype = cellist(kcell)%ctype
!			if (gen == 1) then
!				p%stagetime = tnow
!			else
!				p%stagetime = tnow + dividetime(gen,ctype)
!			endif
!		else
!			p%stagetime = 0
!		endif
!	case (ACTIVATED)     ! possible transition from ACTIVATED to DIVIDING
!        gen = get_generation(p)
!        ctype = cellist(kcell)%ctype
!		if (gen == TC_MAX_GEN) then
!            call set_stage(p,FINISHED)
!			p%stagetime = BIG_TIME
!		elseif (candivide(p,ctype)) then
!            call set_stage(p,DIVIDING)
!			p%stagetime = tnow + CYTOKINESIS_TIME
!		endif
	case (DIVIDING)
	    divide_flag = .true.
	case (FINISHED)

    end select
endif
end subroutine

!--------------------------------------------------------------------------------
! Note: USE_STAGETIME(icstage) means use stagetime for the transition from icstage.
! p%stagetime = time that the cell is expected to make transition to the next stage.
! In this simplified version, the sequence is (revised_staging = .true.):
! NAIVE -> TRANSIENT -> CLUSTERS -> SWARMS -> DIVIDING -> SWARMS
! On cell division, the stage is set to SWARMS, and the time for the next stage
! transition (to DIVIDING, provided candivide()) is the division time.
!--------------------------------------------------------------------------------
subroutine updatestage(kcell,tnow,divide_flag)
integer :: kcell
logical :: divide_flag
real :: tnow
integer :: stage, region, gen, ctype
real :: stagetime
type(cog_type), pointer :: p

divide_flag = .false.
p => cellist(kcell)%cptr
!stage = get_stage(p)
call get_stage(p,stage,region)
if (stage >= CLUSTERS) then
    if (.not.cansurvive(p)) then
        write(logmsg,*) 'cell IL2 store too low: ',kcell,p%cogID
        call logger(logmsg)
        p%dietime = tnow
        stage = FINISHED
        call set_stage(p,stage)
    endif
endif
if (stage == FINISHED) return
ctype = cellist(kcell)%ctype
stagetime = p%stagetime
if (tnow > stagetime) then		! time constraint to move to next stage is met
	! May be possible to make the transition from stage
	select case(stage)
	case (NAIVE)        ! possible transition from NAIVE to TRANSIENT
	    if (p%stimulation > 0) then
	        gen = get_generation(p)
            call set_stage(p,TRANSIENT)
            p%dietime = tnow + TClifetime(p)
		    if (USE_STAGETIME(TRANSIENT)) then
			    p%stagetime = tnow + get_stagetime(p,ctype)
		    else
			    p%stagetime = 0		! time is not criterion for next transition
            endif
        endif
	case (TRANSIENT)    ! possible transition from TRANSIENT to CLUSTERS
	    if (reached_IL2_threshold(p)) then
	        nIL2thresh = nIL2thresh + 1
	        tIL2thresh = tIL2thresh + (tnow - cellist(kcell)%entrytime)
!	        write(*,'(a,i6,2f8.1)') '========= Reached IL2 threshold: ',kcell,(tnow - cellist(kcell)%entrytime),p%stimulation
            call set_stage(p,CLUSTERS)
		    if (USE_STAGETIME(CLUSTERS)) then
			    p%stagetime = tnow + get_stagetime(p,ctype)
		    else
			    p%stagetime = 0		! time is not criterion for next transition
		    endif
        endif
	case (CLUSTERS)     ! possible transition from CLUSTERS to SWARMS
	    if (reached_act_threshold(p)) then
            call set_stage(p,SWARMS)
            gen = get_generation(p)
            ctype = cellist(kcell)%ctype
		    if (USE_STAGETIME(SWARMS)) then
				p%stagetime = tnow + dividetime(gen,ctype)
		    else
			    p%stagetime = 0		! time is not criterion for next transition
		    endif
        endif
	case (SWARMS)       ! possible transition from SWARMS to DIVIDING
        gen = get_generation(p)
        ctype = cellist(kcell)%ctype
		if (gen == TC_MAX_GEN) then
            call set_stage(p,FINISHED)
			p%stagetime = BIG_TIME
!			write(logmsg,*) 'updatestage: division limit reached for cell: ',kcell
!			call logger(logmsg)
		elseif (candivide(p,ctype)) then
            call set_stage(p,DIVIDING)
		    if (USE_STAGETIME(DIVIDING)) then
!				p%stagetime = tnow + CYTOKINESIS_TIME
				p%stagetime = tnow + get_stagetime(p,ctype)
		    else
			    p%stagetime = 0		! time is not criterion for next transition
		    endif
		endif
	case (DIVIDING)
	    divide_flag = .true.
	case (FINISHED)

    end select
endif
end subroutine

!--------------------------------------------------------------------------------------
! The stage time is the time that a cell will spend in a stage (after naive stage).
! May possibly on cell type (CD4/CD8).  Applies only to cognate cells.
! Note that the concept of time in a stage is just a surrogate for accurate biological
! mechanisms that determine stage transitions.
! We might be able to override time in stage 4 (SWARMS) by use of activation threshold.
!--------------------------------------------------------------------------------------
real function get_stagetime(p,ctype)
type(cog_type), pointer :: p
integer :: ctype
integer :: stage, region, gen, cd4_8

!stage = get_stage(p)
call get_stage(p,stage,region)
cd4_8 = ctype - 1       ! converts ctype -> 1=CD4, 2=CD8
gen = get_generation(p)
if (gen == 1) then
	if (stage == SWARMS) then	! use 1st division time (only if revised_stagetime = .false.)
		get_stagetime = dividetime(gen,cd4_8)
	else
		get_stagetime = 60*mean_stagetime_1(stage,cd4_8)		! hr -> min
	endif
else
	get_stagetime = 60*mean_stagetime_n(stage,cd4_8)		! hr -> min
endif
end function

!--------------------------------------------------------------------------------
! Returns .true. if the cognate cell is permitted to bind to a DC.
!--------------------------------------------------------------------------------
logical function bindable(p)
type(cog_type), pointer :: p
integer :: stage, region

bindable = .true.
!if (get_stage(p) == FINISHED) then
call get_stage(p,stage,region)
if (stage == FINISHED) then
    bindable = .false.
endif
end function

!--------------------------------------------------------------------------------------
! The bind time depends on the cell stage, and possibly on cell type (CD4/CD8).
! Rules currently apply only during initial activation, Phases 1 & 2 (TRANSIENT, CLUSTERS)
! CT_CONSTANT
! Just use mean bind times for cognate cells.
! CT_HILL
! Now try making bind time depend on the TCR stimulation rate of the binding
! CT_HENRICKSON
! The bind time is drawn from a lognormal distn, fixed shape, median increases with
! the level of TCR stimulation, up to a limiting value.
! HYPOTHESIS: bind time should also increase with level of stimulation rate,
! i.e. CT_HILL should be combined with CT_HENRICKSON - but how?
!--------------------------------------------------------------------------------------
real function get_bindtime(p,cognate,ctype,pMHC,kpar)
type(cog_type), pointer :: p
logical :: cognate
integer :: ctype, kpar
real :: pMHC
integer :: stage, region, i, cd4_8
real(DP) :: R
real :: stimrate, btime, p1, median
real, parameter :: CT_shape = 2.0, CT_max_factor = 4.0
real, parameter :: p2 = log(CT_shape)

get_bindtime = 0
if (cognate) then
    cd4_8 = ctype - 1
!	stage = get_stage(p)
	call get_stage(p,stage,region)
	btime = dc_mean_bindtime_c(stage,cd4_8)
	if (stage > NAIVE .and. stage < SWARMS) then
	    select case (CONTACT_RULE)
	    case (CT_CONSTANT)
	        get_bindtime = btime
	    case (CT_HILL)
	        stimrate = stimulation_rate(pMHC,p%avidity)
	        get_bindtime = btime*(0.10 + stimrate)
	    case (CT_HENRICKSON)
	        median = CT_median(p,pMHC)
	        p1 = log(median)
	        get_bindtime = min(CT_max_factor*median,rv_lognormal(p1,p2,kpar))
	    end select
	else
	    get_bindtime = btime
    endif
else
    R = par_uni(kpar)
	do i = 1,Nbindtime_nc
		if (R < dc_cummul_prob_nc(i)) then
			get_bindtime = dc_bindtime_nc(i)
			return
		endif
	enddo
	get_bindtime = dc_bindtime_nc(Nbindtime_nc)
endif
end function

!----------------------------------------------------------------------------------------
! Estimates the median DC contact time on the basis of the level of TCR stimulation.
! The time varies linearly with stimulation, from a minimum value at zero stimulation,
! bounded by a maximum value.
! Use a ramp function for now.
! Min = 5, max = 70
! The slope is set by requiring that the median when S = IL2_THRESHOLD, i.e. the
! Phase 1/Phase 2 transition, takes a specified value T12_median
! (Currently pMHC is not used)
!----------------------------------------------------------------------------------------
real function CT_median(p,pMHC)
type(cog_type), pointer :: p
real :: pMHC
real, parameter :: min_median = 5, max_median = 70, T12_median = 25		!<------ hard-wired
real :: slope

slope = (T12_median - min_median)/IL2_THRESHOLD
CT_median = min(max_median,min_median + slope*p%stimulation)
end function

!----------------------------------------------------------------------------------------
! A cell gets permission to divide when it is in the ACTIVATED stage and the activation
! (weighted sum of TCR stimulation and CD25/IL2 signal) exceeds a threshold.
! Note that when optionA = 1, the CD25 signal alone must also exceed a threshold for division,
! and if CD25_SWITCH is true and gen = 1 failure to reach the thresholds cancels division
! permanently.
! In this revised version, the use of a weighted sum of signals for activation (with the
! parameter TC_STIM_WEIGHT) has been separated from the optionA cases.
!----------------------------------------------------------------------------------------
logical function candivide(p,ctype)
type(cog_type), pointer :: p
integer :: ctype
integer :: gen, cd4_8
real :: div_thresh, act, CD25signal
!
! NOTE: Need to better account for first division time, and need to check CD4/CD8
!
candivide = .false.
!p => cellist(kcell)%cptr
!cd4_8 = cellist(kcell)%ctype - 1
cd4_8 = ctype - 1
gen = get_generation(p)
if (gen == 1) then			! naive cell
	div_thresh = FIRST_DIVISION_THRESHOLD(cd4_8)
else											! clone, lower threshold for division
	div_thresh = DIVISION_THRESHOLD(cd4_8)
endif
act = get_activation(p)     ! weighted sum of TCR signal and CD25/IL2 signal
if (optionA == 1) then
    CD25signal = get_IL2store(p)
    if (act > div_thresh .and. CD25signal > CD25_DIVISION_THRESHOLD) then ! allow division
	    if (gen > TC_MAX_GEN) then
		    write(*,*) 'candivide: bad gen: ',gen
		    stop
	    endif
	    candivide = .true.
    elseif (gen == 1 .and. CD25_SWITCH) then
    !	Tcell(icell)%cog%IL2_pass = .false.		! not enough IL2, division is cancelled
    !	Tcell(icell)%cog%IL2status = ibclr(Tcell(icell)%cog%IL2status,IL2_FLAG3_BIT)
    !	write(*,*) 'Division cancelled for: ',kcell
        call set_stage(p,FINISHED)
	    p%stagetime = BIG_TIME
	endif
elseif (optionA == 2) then
    if (act > div_thresh) then ! allow division
	    if (gen > TC_MAX_GEN) then
		    write(*,*) 'candivide: bad gen: ',gen
		    stop
	    endif
	    candivide = .true.
	endif
endif
end function

!----------------------------------------------------------------------------------------
! For now base this decision on the cytokine production threshold.
! This must happen only once for a given cell.(?)
!----------------------------------------------------------------------------------------
logical function reached_IL2_threshold(p)
type(cog_type), pointer :: p

!if (p%stimulation > cytokine(IL2_CYT)%cyt_threshold) then
if (p%stimulation > IL2_THRESHOLD) then
	reached_IL2_threshold = .true.
else
	reached_IL2_threshold = .false.
endif
end function

!----------------------------------------------------------------------------------------
! This must happen only once for a given cell.(?)
!----------------------------------------------------------------------------------------
logical function reached_act_threshold(p)
type(cog_type), pointer :: p
real :: activation

activation = get_activation(p)
if (activation > ACTIVATION_THRESHOLD) then
	reached_act_threshold = .true.
else
	reached_act_threshold = .false.
endif
end function

!----------------------------------------------------------------------------------------
! T cell activation is determined as the weighted sum of TCR stimulation and IL2_Store,
! which is the current level of integrated CD25 signal.
!----------------------------------------------------------------------------------------
real function get_activation(p)
type(cog_type), pointer :: p

!if (optionA == 1) then
    get_activation = TC_STIM_WEIGHT*p%stimulation + (1-TC_STIM_WEIGHT)*get_IL2store(p)
!else
!    get_activation = p%stimulation
!endif
end function

!----------------------------------------------------------------------------------------
! IL2_production_status is .true. if IL-2 and CD25 are being produced, else .false.
! OptionB:
! = 1  IL2 is produced for a maximum period of IL2_PRODUCTION_TIME in gen = 1 only
! = 2  IL2 is produced in gen = 1 only
! = 3  IL2 is produced always
!----------------------------------------------------------------------------------------
logical function IL2_production_status(p,t)
type(cog_type), pointer :: p
real :: t
integer :: gen, stage, region

IL2_production_status = .false.
!stage = get_stage(p)
call get_stage(p,stage,region)
if (stage == FINISHED) then
    return
endif
gen = get_generation(p)
select case (optionB)
case (1)
    if (t < IL2_PRODUCTION_TIME*60 .and. gen == 1) then
        IL2_production_status = .true.
    else
        IL2_production_status = .false.
    endif
case (2)
    if (gen == 1) then
        IL2_production_status = .true.
    else
        IL2_production_status = .false.
    endif
case (3)
    IL2_production_status = .true.
end select
end function

!----------------------------------------------------------------------------------------
! If optionC = 1, an activated cell cannot survive if IL2 store drops below
! CD25_SURVIVAL_THRESHOLD
!----------------------------------------------------------------------------------------
logical function cansurvive(p)
type(cog_type), pointer :: p
real :: CD25signal

cansurvive = .true.
if (optionC == 1) then
    CD25signal = get_IL2Store(p)
    if (CD25signal < CD25_SURVIVAL_THRESHOLD) then
        cansurvive = .false.
    endif
endif
end function

end module
