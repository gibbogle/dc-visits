!-----------------------------------------------------------------------------------------
subroutine chemo_jumper(kcell,indx1,kslot1,go,kpar)
integer :: kpar,kcell,indx1(2),kslot1
logical :: go
type (cell_type), pointer :: cell
integer :: fullslots1,fullslots2,site1(3),site2(3),kslot2,stype
integer :: irel,dir1,lastdir1,indx2(2),k,kr,rv(3),id, ichemo, nfull, nrest, nout
integer :: savesite2a(3,MAXRELDIR+1), saveslots2a(MAXRELDIR+1)
real(DP) :: p(MAXRELDIR+1),psum, R, pR, psumm, stay_prob,  psave(MAXRELDIR+1)
real :: tnow, v(3), vsum(3), f
logical :: ischemo, cognate

tnow = istep*DELTA_T
cell => cellist(kcell)
cognate = associated(cell%cptr)
	
id = cell%id
site1 = cell%site
if (site1(1) < 1) then
    write(logmsg,*) 'chemo_jumper: bad site1: ',site1
    call logger(logmsg)
    stop
endif
fullslots1 = 0
do k = 1,2
    if (indx1(k) > 0) then
        fullslots1 = fullslots1 + k
    endif
enddo
stay_prob = dirprob(0)

ischemo = .false.
do kr = 1,MAX_RECEPTOR
    if (cell%receptor_saturation_time(kr) /= 0) then
!        if (tnow > cell%receptor_saturation_time(kr) + T_RECEPTOR_REFRACTORY) then
        if (tnow > cell%receptor_saturation_time(kr) + receptor(kr)%refractory_time) then
            cell%receptor_saturation_time(kr) = 0
        endif
    endif
    if (receptor(kr)%used .and. (cell%receptor_level(kr) > 0) .and. (cell%receptor_saturation_time(kr) == 0)) then
        ischemo = .true.
        exit
    endif
enddo

vsum = 0
if (ischemo) then
    do kr = 1,MAX_RECEPTOR
        if (receptor(kr)%used .and. (cell%receptor_saturation_time(kr) == 0)) then
			ichemo = receptor(kr)%chemokine
			f = receptor(kr)%sign*cell%receptor_level(kr)*receptor(kr)%strength
			v = chemo(ichemo)%grad(:,site1(1),site1(2),site1(3))
			if (receptor_saturation(kr,site1,f,v)) then
			    cell%receptor_saturation_time(kr) = tnow
			    f = f/2
			endif
	    	vsum = vsum + f*v
	    endif
	enddo
	! For exit chemotaxis:
	! Need to create estimate of v() that corresponds to the direction of vsum,
	! nearest discrete location on the 3D lattice (for chemo_p(x,y,z))
	! This is an approximation to increase speed - it enables a table lookup.
	! Note that we use only the direction of v (magnitude is insignificant at this stage,
	! since it is accounted for in f)
    if (norm(vsum) > 0) then
		f = min(1.0,norm(vsum))     ! Note: f is in (0-1)
		v = vsum/norm(vsum)
		rv = chemo_N*v
	else
		f = 0
		ischemo = .false.
	endif
    stay_prob = dirprob(0)
    stay_prob = (1-f)*stay_prob
else
    stay_prob = dirprob(0)
endif

if (fullslots1 /= BOTH) then
    R = par_uni(kpar)
    if (R <= stay_prob) then    ! case of no jump
	    go = .false.
        return
    endif
endif
! Now we must jump (if possible)

stype = struct_type(int(cell%ctype))     ! COG_TYPE_TAG or NONCOG_TYPE_TAG

! Compute jump probabilities in the absence of chemotaxis
site1 = cell%site
lastdir1 = cell%lastdir
p = 0
savesite2a = 0
saveslots2a = 0
nfull = 0
nrest = 0
nout = 0
do irel = 1,nreldir
	dir1 = reldir(lastdir1,irel)
	site2 = site1 + jumpvec(:,dir1)
	if (inside_xyz(site2)) then
	    indx2 = occupancy(site2(1),site2(2),site2(3))%indx
		if (indx2(1) >= 0) then     ! not OUTSIDE_TAG or DC
            fullslots2 = 0
            do k = 1,2
                if (indx2(k) > 0) then
                    fullslots2 = fullslots2 + k
                endif
            enddo
            if (fullslots2 == BOTH) then
                nfull = nfull + 1
                cycle
            elseif (fullslots2 /= 0) then
                nrest = nrest + 1
                p(dir1) = dirprob(irel)*GAMMA
            else
                nrest = nrest + 1
                p(dir1) = dirprob(irel)
            endif
            saveslots2a(dir1) = fullslots2
        else
	        nout = nout + 1
		endif
	endif
	savesite2a(:,dir1) = site2
enddo
if (sum(p) == 0) then
    go = .false.
    return
endif

if (ischemo) then
	psave = p
!	call chemo_probs(p,v,f)
	call chemo_probs_pre(p,rv,f)     ! this is the precomputed version
endif
psum = sum(p)

if (psum == 0) then
	go = .false.
	return
else
    go = .true.
endif

! Now choose a direction on the basis of these probs p()
R = par_uni(kpar)
pR = psum*R
psumm = 0
do dir1 = 1,njumpdirs
   	psumm = psumm + p(dir1)
   	if (pR <= psumm) then
   		exit
   	endif
enddo
if (dir1 > njumpdirs) then
    dir1 = 0
    do k = 1,njumpdirs
        if (p(k) > 0) then
            dir1 = k
            exit
        endif
    enddo
endif
site2 = savesite2a(:,dir1)

fullslots2 = saveslots2a(dir1)
! new code
if (diagonal_jumps) then
	dir1 = fix_lastdir(dir1,kpar)
elseif (dir1 == 0) then
	dir1 = random_int(1,6,kpar)
endif

if (fullslots2 == 0) then       ! randomly select a slot
    R = par_uni(kpar)
    if (R <= 0.5) then
        kslot2 = SLOT_NUM1
    else
        kslot2 = SLOT_NUM2
    endif
elseif (fullslots2 == SLOT_NUM1) then
    kslot2 = SLOT_NUM2
elseif (fullslots2 == SLOT_NUM2) then
    kslot2 = SLOT_NUM1
else
    write(logmsg,*) 'ERROR in jumper: jump to crowded site'
	call logger(logmsg)
    stop
endif
cell%site = site2
cell%lastdir = dir1
occupancy(site2(1),site2(2),site2(3))%indx(kslot2) = kcell
occupancy(site1(1),site1(2),site1(3))%indx(kslot1) = 0
end subroutine

!-----------------------------------------------------------------------------------------
! Determine whether a receptor at site1(:) is saturated.
! Currently the decision is made from:
! f = total strength, which is cell%receptor_level(kr)*receptor(kr)%strength
! v = chemokine gradient
!-----------------------------------------------------------------------------------------
logical function receptor_saturation(kr,site,f,v)
integer :: kr, site(3)
real :: f, v(3), a(3)
real :: magnitude
!real :: saturation_threshold = 0.8

a = f*v
magnitude = norm(a)
if (magnitude > receptor(kr)%saturation_threshold) then
    receptor_saturation = .true.
else
    receptor_saturation = .false.
endif
end function

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