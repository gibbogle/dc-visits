module omp_motility

!use omp_transfer
use omp_global
use chemokine

implicit none

contains

!--------------------------------------------------------------------------------
! The jump site is selected from all neighbouring sites on the basis of jump
! direction probs for both this cell and for the neighbours.
! If go = .false. the cell doesn't move
! (Try adding a drift probability in the case of traffic)
!--------------------------------------------------------------------------------
subroutine jumper(kcell,indx1,kslot1,go,kpar)
integer :: kpar,kcell,indx1(2),kslot1
logical :: go
type (cell_type), pointer :: cell
integer :: fullslots1,fullslots2,site1(3),site2(3),kslot2,stype
integer :: irel,dir1,lastdir1,indx2(2),k,z
integer :: savesite2(3,MAXRELDIR), saveslots2(MAXRELDIR)
integer :: savesite2a(3,MAXRELDIR+1), saveslots2a(MAXRELDIR+1)
real(DP) :: psum, p(MAXRELDIR+1), R, pR, psumm 

cell => cellist(kcell)
site1 = cell%site
if (site1(1) < 1) then
    write(*,*) 'jumper: bad site1: ',site1
    stop
endif
fullslots1 = 0
do k = 1,2
    if (indx1(k) > 0) then
        fullslots1 = fullslots1 + k
    endif
enddo
if (fullslots1 /= BOTH) then
    R = par_uni(kpar)
    if (R <= dirprob(0)) then    ! case of no jump
	    go = .false.
        return
    endif
endif
! Now we must jump (if possible)

!stype = struct_type(int(cell%ctype))     ! COG_TYPE_TAG or NONCOG_TYPE_TAG
if (associated(cell%cptr)) then
	stype = COG_TYPE_TAG
else
	stype = NONCOG_TYPE_TAG
endif

z = site1(3)
lastdir1 = cell%lastdir
psum = 0
p = 0
savesite2a = 0
saveslots2a = 0
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
                cycle
            elseif (fullslots2 /= 0) then
!                p(irel) = dirprob(irel)*GAMMA
                p(dir1) = dirprob(irel)*GAMMA
            else
!                p(irel) = dirprob(irel)
                p(dir1) = dirprob(irel)
            endif
            if (exit_region == EXIT_LOWERHALF) then
                if (evaluate_residence_time .or. stype == NONCOG_TYPE_TAG) then
                    if (site2(3) < z) then
                        if (nz_excess(z) > 0) then
!                            p(irel) = p(irel)*(1 + excess_factor*nz_excess(z)/real(nz_sites(z)))
                            p(dir1) = p(dir1)*(1 + excess_factor*nz_excess(z)/real(nz_sites(z)))
                        endif
                    endif
                endif
            endif
!            psum = psum + p(irel)
!            saveslots2(irel) = fullslots2
            saveslots2a(dir1) = fullslots2
		endif
	endif
!	savesite2(:,irel) = site2
	savesite2a(:,dir1) = site2
enddo
psum = sum(p)
if (psum == 0) then
	go = .false.
	return
else
    go = .true.
endif

! Now choose a direction on the basis of these probs p()
R = 0
R = par_uni(kpar)
pR = psum*R
psumm = 0
!do irel = 1,nreldir
do dir1 = 1,njumpdirs
   	psumm = psumm + p(dir1)
   	if (pR <= psumm) then
   		exit
   	endif
enddo
!if (irel > nreldir) irel = nreldir
!dir1 = reldir(lastdir1,irel)
if (dir1 > njumpdirs) then
    dir1 = 0
    do k = 1,njumpdirs
        if (p(k) > 0) then
            dir1 = k
            exit
        endif
    enddo
    if (dir1 == 0) then
        write(*,*) 'jumper: bad dir1: ',dir1
        write(*,*) 'R, psum, psumm: ',R,psum,psumm
        write(*,*) p
        stop
    endif
endif
!site2 = savesite2(:,irel)
site2 = savesite2a(:,dir1)
!fullslots2 = saveslots2(irel)
fullslots2 = saveslots2a(dir1)
! new code
if (diagonal_jumps) then
	dir1 = fix_lastdir(dir1,kpar)
elseif (dir1 == 0) then
	dir1 = random_int(1,6,kpar)
endif

if (fullslots2 == 0) then       ! randomly select a slot
    R = par_uni(kpar)
    if (R < 0.5) then
        kslot2 = SLOT_NUM1
    else
        kslot2 = SLOT_NUM2
    endif
elseif (fullslots2 == SLOT_NUM1) then
    kslot2 = SLOT_NUM2
elseif (fullslots2 == SLOT_NUM2) then
    kslot2 = SLOT_NUM1
else
    write(*,*) 'ERROR in jumper: jump to crowded site'
    stop
endif
cell%site = site2
cell%lastdir = dir1
occupancy(site2(1),site2(2),site2(3))%indx(kslot2) = kcell
occupancy(site1(1),site1(2),site1(3))%indx(kslot1) = 0
end subroutine

!-----------------------------------------------------------------------------------------
! The cell chemotactic susceptibility is conveyed by 'chemo', which may be the S1P1 level.
!-----------------------------------------------------------------------------------------
subroutine chemo_jumper1(kcell,indx1,kslot1,go,kpar)
integer :: kpar,kcell,indx1(2),kslot1
logical :: go
type (cell_type), pointer :: cell
integer :: fullslots1,fullslots2,site1(3),site2(3),kslot2,stype
integer :: irel,dir1,lastdir1,indx2(2),k,z,e(3),v(3)
!integer :: ne, ee(3,MAX_EX_CHEMO), nd, dd(3,MAX_DC_CHEMO)
integer :: savesite2a(3,MAXRELDIR+1), saveslots2a(MAXRELDIR+1)
real(DP) :: psum, R, pR, psumm, p(MAXRELDIR+1), stay_prob, f, c
real :: rad, chemo_exit
logical :: in_exit_SOI	!, in_DC_SOI

cell => cellist(kcell)
site1 = cell%site
if (site1(1) < 1) then
    write(*,*) 'jumper: bad site1: ',site1
    stop
endif
fullslots1 = 0
do k = 1,2
    if (indx1(k) > 0) then
        fullslots1 = fullslots1 + k
    endif
enddo
! Need first to determine if the cell is subject to chemotaxis
in_exit_SOI = .false.
if (Nexits > 0) then
    chemo_exit = chemo_active_exit(cell)    ! the degree of chemotactic activity, possibly based on S1P1 as surrogate
    if (chemo_exit > 0) then
        ! Then need to determine if the cell is within the SOI of an exit.
        call nearest_exit(site1,in_exit_SOI,e)
    endif
endif
if (in_exit_SOI) then
	v = site1 - e
	rad = chemo_r(abs(v(1)),abs(v(2)),abs(v(3)))
	if (rad == 0) then
		go = .false.
		return
	endif
	f = chemo_K_exit*chemo_exit*chemo_g(rad)	! Note: inserting chemo_K_exit here (was missing) will mean need to
												! change ep_factor in portal_traffic()
    stay_prob = dirprob(0)
    c = 1
!    c = CCR7_ligand(rad)
!    stay_prob = (1-f)*(1 - c + c*stay_prob)
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

!stype = struct_type(int(cell%ctype))     ! COG_TYPE_TAG or NONCOG_TYPE_TAG
if (associated(cell%cptr)) then
	stype = COG_TYPE_TAG
else
	stype = NONCOG_TYPE_TAG
endif

z = site1(3)
lastdir1 = cell%lastdir
p = 0
savesite2a = 0
saveslots2a = 0
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
                cycle
            elseif (fullslots2 /= 0) then
!                p(irel) = dirprob(irel)*GAMMA
                p(dir1) = dirprob(irel)*GAMMA
            else
!                p(irel) = dirprob(irel)
                p(dir1) = dirprob(irel)
            endif
!            saveslots2(irel) = fullslots2
            saveslots2a(dir1) = fullslots2
		endif
	endif
!	savesite2(:,irel) = site2
	savesite2a(:,dir1) = site2
enddo
if (sum(p) == 0) then
    go = .false.
    return
endif

if (in_exit_SOI) then
    if (v(1) == 0 .and. v(2) == 0 .and. v(3) == 0) then
        write(*,*) 'cell at exit: ',kcell
        go = .false.
        return
    endif
    call chemo_probs(p,v,f,c)
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
!do irel = 1,nreldir
do dir1 = 1,njumpdirs
   	psumm = psumm + p(dir1)
   	if (pR <= psumm) then
   		exit
   	endif
enddo
!if (irel > nreldir) irel = nreldir
!dir1 = reldir(lastdir1,irel)
if (dir1 > njumpdirs) then
    dir1 = 0
    do k = 1,njumpdirs
        if (p(k) > 0) then
            dir1 = k
            exit
        endif
    enddo
    if (dir1 == 0) then
        write(*,*) 'chemo_jumper: bad dir1: ',dir1
        write(*,*) 'R, psum, psumm: ',R,psum,psumm
        write(*,*) p
        stop
    endif
endif
!site2 = savesite2(:,irel)
site2 = savesite2a(:,dir1)
!fullslots2 = saveslots2(irel)
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
    write(*,*) 'ERROR in jumper: jump to crowded site'
    stop
endif
cell%site = site2
cell%lastdir = dir1
occupancy(site2(1),site2(2),site2(3))%indx(kslot2) = kcell
occupancy(site1(1),site1(2),site1(3))%indx(kslot1) = 0

end subroutine

!-----------------------------------------------------------------------------------------
! The cell chemotactic susceptibility is conveyed by 'chemo', which may be the S1P1 level.
! In the limit of CHEMO_K_EXIT = 0, this should be the same as jumper(), i.e. the same
! as setting USE_EXIT_CHEMOTAXIS = false.
!-----------------------------------------------------------------------------------------
subroutine chemo_jumper2(kcell,indx1,kslot1,go,kpar)
integer :: kpar,kcell,indx1(2),kslot1
logical :: go
type (cell_type), pointer :: cell
integer :: fullslots1,fullslots2,site1(3),site2(3),kslot2,stype
integer :: irel,dir1,lastdir1,indx2(2),k,z,e(3),v(3)
integer :: ne, ee(3,MAX_EX_CHEMO), nd, neardc(MAX_DC_CHEMO), kk, ned, idc
integer :: savesite2a(3,MAXRELDIR+1), saveslots2a(MAXRELDIR+1)
real(DP) :: p(MAXRELDIR+1),psum, R, pR, psumm, stay_prob, f, c
real :: tnow, rad, chemo_exit=0, chemo_DC=0
real :: ff(MAX_EX_CHEMO+MAX_DC_CHEMO), vv(3,MAX_EX_CHEMO+MAX_DC_CHEMO), vsum(3)

!if (kcell == 310 .and. istep == 4801) then
!	dbug = .true.
!else
	dbug = .false.
!endif
tnow = istep*DELTA_T
cell => cellist(kcell)
site1 = cell%site
if (site1(1) < 1) then
    write(*,*) 'chemo_jumper: bad site1: ',site1
    stop
endif
fullslots1 = 0
do k = 1,2
    if (indx1(k) > 0) then
        fullslots1 = fullslots1 + k
    endif
enddo
stay_prob = dirprob(0)
ned = 0

! Need first to determine if the cell is subject to chemotaxis
ne = 0
if (Nexits > 0 .and. chemo_K_exit > 0) then
    chemo_exit = chemo_active_exit(cell)    ! the degree of chemotactic activity, possibly based on S1P1 as surrogate
    if (chemo_exit > 0) then
        ! Then need to determine if the cell is within the SOI of any exits.
        call near_exits(site1,ne,ee) ! new
        if (dbug) then
			write(*,*) 'Near exits: ',ne
			do k = 1,ne
				write(*,*) 'ee: ',k,ee(:,k)
			enddo
		endif
    endif
endif
nd = 0
chemo_DC = 0
if (use_DC_chemotaxis .and. chemo_K_DC > 0 .and. NDC > 0) then
	if (TAGGED_DC_CHEMOTAXIS) then
		if (cell%tag == TAGGED_CELL) then
			if (tnow-cell%unbindtime(2) > DC_CHEMO_DELAY) then
				chemo_DC = chemo_active_DC(cell)    ! the degree of chemotactic activity
			endif
		endif
	elseif (tnow-cell%unbindtime(2) > DC_CHEMO_DELAY) then
		chemo_DC = chemo_active_DC(cell)    ! the degree of chemotactic activity
	endif
    if (chemo_DC > 0) then
		if (USE_CHEMOKINE_GRADIENT) then
			nd = 1
			if (ne /= 0) then
				write(logmsg,*) 'Combination of DC gradient chemotaxis and exit chemotaxis not yet implemented'
				call logger(logmsg)
				stop
			endif
		else
	        ! Need to determine if the cell is within the SOI of any DCs.
		    call near_DCs(site1,nd,neardc) ! new
		endif
    endif
endif

vsum = 0
ff = 0
if (ne > 0 .and. chemo_K_exit > 0) then
	do k = 1,ne
		e = ee(:,k)
!		v = site1 - e	! Note: vector v directed AWAY from exit 
! CHANGED
		v = e - site1	! v now directed towards attractant site - like chemokine gradient vector
		rad = chemo_r(abs(v(1)),abs(v(2)),abs(v(3)))
		if (rad == 0) then
			f = 0
		else
			f = chemo_K_exit*chemo_exit*chemo_g(rad)  
			vsum = vsum + (f/norm(real(v)))*v
		endif
	enddo
endif

if (nd > 0) then
	! augment vsum with DC chemotaxis
	if (USE_CHEMOKINE_GRADIENT) then
!		vsum = -occupancy(site1(1),site1(2),site1(3))%chemo_grad	! - sign for consistency with radius method
		vsum = chemo(1)%grad(:,site1(1),site1(2),site1(3))	
		! sign change no longer needed - chemo_p(:) changed 
		vsum = CHEMO_K_DC*vsum	
	else
		do k = 1,nd
			kk = ne + k
			idc = neardc(k)
			e = DClist(idc)%site
			v = site1 - e	! Note: vector v directed AWAY from DC
	!		vv(:,kk) = v ! new
			rad = chemo_r(abs(v(1)),abs(v(2)),abs(v(3)))
			if (rad == 0) then
				go = .false.
				return
			endif
			f = chemo_K_DC*chemo_DC*chemo_g(rad)
			vsum = vsum + (f/norm(real(v)))*v
		enddo
	endif
endif
ned = ne + nd

if (ned > 0) then
	f = min(1.0,norm(vsum))
	! Need to create estimate of v() that corresponds to the direction of vsum,
	! nearest discrete location on the 3D lattice (for chemo_p(x,y,z))
	! This is an approximation to increase speed - it enables a table lookup.
	! Note that we use only the direction of v (magnitude is insignificant at this stage,
	! since it has been accounted for in f)
	v = chemo_N*vsum/norm(vsum)	! this is a quick approximation, needs checking! OK.
    stay_prob = dirprob(0)
    c = 1
!    c = CCR7_ligand(rad)	! no unique rad now - this needs to be changed
!    stay_prob = (1-f)*(1 - c + c*stay_prob)
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

!stype = struct_type(int(cell%ctype))     ! COG_TYPE_TAG or NONCOG_TYPE_TAG
if (associated(cell%cptr)) then
	stype = COG_TYPE_TAG
else
	stype = NONCOG_TYPE_TAG
endif

! Compute jump probabilities in the absence of chemotaxis
site1 = cell%site
z = site1(3)
lastdir1 = cell%lastdir
p = 0
savesite2a = 0
saveslots2a = 0
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
                cycle
            elseif (fullslots2 /= 0) then
!                p(irel) = dirprob(irel)*GAMMA
                p(dir1) = dirprob(irel)*GAMMA
            else
!                p(irel) = dirprob(irel)
                p(dir1) = dirprob(irel)
            endif
!            saveslots2(irel) = fullslots2
            saveslots2a(dir1) = fullslots2
		endif
	endif
!	savesite2(:,irel) = site2
	savesite2a(:,dir1) = site2
enddo
if (sum(p) == 0) then
    go = .false.
    return
endif

if (nd == 1 .and. kcell == idbug) then
	psum = sum(p)
	idc = neardc(1)
	write(nfchemo,'(a,i4,a,i6,3i4,a,i6,3i4)') 'istep: ',istep,' cell: ',kcell,site1,' DC: ',idc,DClist(idc)%site
	write(nfchemo,'(5x,10f6.3)') p/psum
endif
if (ne*CHEMO_K_EXIT > 0 .or. nd*CHEMO_K_DC > 0) then
	call chemo_probs(p,v,f,c)
endif
psum = sum(p)
if (nd == 1 .and. kcell == idbug) then
	write(nfchemo,'(a,f6.3,a,3f7.3,a,3i2)') 'f: ',f,'vsum: ',vsum,'  v: ',v
	idc = neardc(1)
	write(nfchemo,'(5x,10f6.3)') p/psum
endif

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
    if (dir1 == 0) then
        write(*,*) 'chemo_jumper: bad dir1: ',dir1,kcell,istep
        write(*,*) 'R, psum, psumm: ',R,psum,psumm
        write(*,*) p
        stop
    endif
endif
!site2 = savesite2(:,irel)
site2 = savesite2a(:,dir1)
!write(*,*) 'dir1,site2: ',dir1,site2

!fullslots2 = saveslots2(irel)
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
    write(*,*) 'ERROR in jumper: jump to crowded site'
    stop
endif
cell%site = site2
cell%lastdir = dir1
occupancy(site2(1),site2(2),site2(3))%indx(kslot2) = kcell
occupancy(site1(1),site1(2),site1(3))%indx(kslot1) = 0

end subroutine

!-----------------------------------------------------------------------------------------
! This is the version from bcell-abm
!-----------------------------------------------------------------------------------------
subroutine chemo_jumper(kcell,indx1,kslot1,go,kpar)
integer :: kpar,kcell,indx1(2),kslot1
logical :: go
type (cell_type), pointer :: cell
integer :: fullslots1,fullslots2,site1(3),site2(3),kslot2
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
!			call logger('Receptor resensitized')
        endif
    endif
	f = cell%receptor_level(kr)*receptor(kr)%strength
    if (receptor(kr)%used .and. (f > 0) .and. (cell%receptor_saturation_time(kr) == 0)) then
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
!			    call logger('Receptor saturation')
			endif
	    	vsum = vsum + f*v
!	    	! TRY THIS  - not needed now, chemo_p(:) has been changed
!	    	vsum = -vsum
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

!-----------------------------------------------------------------------------------------
! In this version the parallel section has been made into a subroutine.
! Try varying the sweep order, i.e. 0,1 and 1,0
!-----------------------------------------------------------------------------------------
subroutine mover(ok)
logical :: ok
integer :: kpar=0, sweep, nsweeps, sweep1, sweep2, dsweep

if (Mnodes > 1) then
!DEC$ IF .NOT. DEFINED (_OPENMP)
!    stop
!DEC$ ENDIF
endif

if (Mnodes == 1) then
    nsweeps = 1
else
    nsweeps = 2
endif

if (mod(istep,2) == 0) then
	sweep1 = 0
	sweep2 = nsweeps-1
	dsweep = 1
else
	sweep2 = 0
	sweep1 = nsweeps-1
	dsweep = -1
endif
!do sweep = 0,nsweeps-1
do sweep = sweep1,sweep2,dsweep
	if (IN_VITRO) then
		if (Mnodes > 1) then
		!$omp parallel do
			do kpar = 0,Mnodes-1
				call par_mover2D(sweep,kpar,ok)
			enddo
		!$omp end parallel do
		else
			call par_mover2D(sweep,kpar,ok)
		endif
		if (.not.ok) return
	else
		if (Mnodes > 1) then
		!$omp parallel do
			do kpar = 0,Mnodes-1
				call par_mover(sweep,kpar)
			enddo
		!$omp end parallel do
		else
			call par_mover(sweep,kpar)
		endif
	endif
enddo
ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine par_mover(sweep,kpar)
integer :: sweep
integer :: kpar
integer :: site1(3), kcell, indx(2), slot, z, slice, site(3), stage, region
logical :: go
type(cell_type), pointer :: cell
integer :: z_lo,z_hi

slice = sweep + 2*kpar
do kcell = 1,nlist
	if (Mnodes == 1) then
	    z_lo = 1
	    z_hi = NZ
	else
	    z_lo = zoffset(slice) + 1
	    z_hi = zoffset(slice+1)
	endif
    cell => cellist(kcell)
    if (cell%ID == 0) cycle             ! skip gaps in the list
    if (associated(cell%cptr)) then
		call get_stage(cell%cptr,stage,region)
		if (region /= LYMPHNODE) cycle
	endif
    if (cell%step == istep) cycle
    if (cell%DCbound(1) /= 0 .or. cell%DCbound(2) /= 0) then
        cycle     ! skip cells bound to DC
    endif
    site1 = cell%site
    z = site1(3)
    if (zdomain(z) /= slice) cycle      ! not in the slice for this processor
    indx = occupancy(site1(1),site1(2),site1(3))%indx
    if (indx(1) < 0) then
        write(*,*) 'Error: par_mover: OUTSIDE_TAG or DC: ',kcell,site1,indx
        stop
    endif
    if (kcell == indx(1)) then
        slot = 1
    elseif (kcell == indx(2)) then
        slot = 2
    else
        write(*,'(a,6i8)') 'Error: par_mover: bad indx: ',kcell,site1,indx
        stop
    endif
    ! TESTING ---------- to track down why case use_exit_chemotaxis + chemo_K_exit = 0
    !                    is not the same as not use_exit_chemotaxis
	if (use_exit_chemotaxis .or. use_DC_chemotaxis) then
		call chemo_jumper(kcell,indx,slot,go,kpar)
	else
		call jumper(kcell,indx,slot,go,kpar)
	endif
    cell%step = istep
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! For (mainly) 2D motion, i.e. IN_VITRO simulation.
!-----------------------------------------------------------------------------------------
subroutine par_mover2D(sweep,kpar,ok)
integer :: sweep
integer :: kpar
logical :: ok
integer :: site1(3), kcell, indx(2), slot, x, slice, site(3)
logical :: go
type(cell_type), pointer :: cell
integer :: x_lo,x_hi

slice = sweep + 2*kpar
do kcell = 1,nlist
	if (Mnodes == 1) then
	    x_lo = 1
	    x_hi = NX
	else
	    x_lo = xoffset(slice) + 1
	    x_hi = xoffset(slice+1)
	endif
    cell => cellist(kcell)

    if (cell%ID == 0) cycle             ! skip gaps in the list
    if (cell%step == istep) cycle
    if (cell%DCbound(1) /= 0 .or. cell%DCbound(2) /= 0) then
        cycle     ! skip cells bound to DC
    endif
    site1 = cell%site
    x = site1(1)
    if (xdomain(x) /= slice) cycle      ! not in the slice for this processor
    indx = occupancy(site1(1),site1(2),site1(3))%indx
    if (indx(1) < 0) then
        write(logmsg,'(a,6i6)') 'Error: par_mover2D: OUTSIDE_TAG or DC: ',kcell,site1,indx
        call logger(logmsg)
        ok = .false.
        return
    endif
    if (kcell == indx(1)) then
        slot = 1
    elseif (kcell == indx(2)) then
        slot = 2
    else
        write(logmsg,'(a,6i8)') 'Error: par_mover2D: bad indx: ',kcell,site1,indx
        call logger(logmsg)
        ok = .false.
        return
    endif
    if (use_exit_chemotaxis) then
        call chemo_jumper2D(kcell,indx,slot,go,kpar)
        site = cellist(kcell)%site
    else
        call jumper2D(kcell,indx,slot,go,kpar)
    endif
    cell%step = istep
enddo
ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
! In this version the parallel section is in this subroutine, requiring the use of
! a PRIVATE clause.
! Note that if cell was not declared as a pointer, the FIRSTPRIVATE(cell) clause is
! needed to avoid an error when cell is assigned.
!-----------------------------------------------------------------------------------------
subroutine mover2
integer :: kpar
integer :: site1(3), kcell, indx(2), slot, x, xlocal
logical :: go
type(cell_type),pointer :: cell
integer :: x_lo,x_hi,sweep
integer :: i, nslice, sum, sump, nsweeps, cnt
integer, save :: xlim(0:8)
integer :: xtotal(0:8)
integer, allocatable :: xcount(:)

if (istep == 1) then    ! must be executed when the blob size changes
    allocate(xcount(NX))
    xcount = 0
    do kcell = 1,nlist
        cell => cellist(kcell)
        if (cell%ID == 0) cycle             ! skip gaps in the list
        i = cellist(kcell)%site(1)
        xcount(i) = xcount(i) + 1
    enddo
    nslice = NTcells/(2*Mnodes)
    i = 0
    sum = 0
    sump = 0
    i = 1
    do x = 1,NX
        sum = sum + xcount(x)
        if (sum > i*nslice) then
            xlim(i) = x
            xtotal(i) = sum - sump
            sump = sum
            write(*,*) 'i,sum: ',i,sum,xtotal(i)
            i = i+1
            if (i == 2*Mnodes) exit
        endif
    enddo
    xlim(2*Mnodes) = NX
    xtotal(2*Mnodes) = NTcells - sum
    write(*,*) 'i,sum: ',i,sum,xtotal(i)
    deallocate(xcount)
endif

if (Mnodes > 1) then
!DEC$ IF .NOT. DEFINED (_OPENMP)
    stop
!DEC$ ENDIF
endif

if (Mnodes == 1) then
    nsweeps = 1
else
    nsweeps = 2
endif

!write(*,*) 'start sweep loop'
do sweep = 0,nsweeps-1
!write(*,*) 'sweep: ',sweep

!$omp parallel PRIVATE(kcell,kpar,cell,x_lo,x_hi,site1,xlocal,indx,slot,go,cnt)
kpar = omp_get_thread_num()
!write(*,*) 'kpar: ',kpar
do kcell = 1,nlist
	if (Mnodes == 1) then
	    x_lo = 1
	    x_hi = NX
	else
	    x_lo = xlim(sweep+2*kpar) + 1
	    x_hi = xlim(sweep+2*kpar+1)
	endif
    cell => cellist(kcell)
    if (cell%ID == 0) cycle             ! skip gaps in the list
    if (cell%step == istep) cycle
    if (cell%DCbound(1) /= 0 .or. cell%DCbound(2) /= 0) cycle     ! skip cells bound to DC
    site1 = cell%site
    xlocal = site1(1)
    if (xlocal < x_lo .or. xlocal > x_hi) cycle      ! not in the slice for this processor
    indx = occupancy(site1(1),site1(2),site1(3))%indx
    if (indx(1) < 0) then
        write(*,*) 'stage1: OUTSIDE_TAG or DC: ',kcell,site1,indx
        stop
    endif
    if (kcell == indx(1)) then
        slot = 1
    elseif (kcell == indx(2)) then
        slot = 2
    else
        write(*,'(a,7i6)') 'ERROR: stage1: bad indx: ',kcell,site1,indx
        stop
    endif
    call jumper(kcell,indx,slot,go,kpar)
    cell%step = istep
enddo
!$omp end parallel
enddo
end subroutine

!--------------------------------------------------------------------------------
! 2D motility
! In this case only one cell can occupy a site - no slot 2 - and no passing.
! This is valid because cells are not crammed together in the petri dish.
! We now wish to relax this constraint.  A site might be crowded, as a result
! of cell division, but in jumper2D cells can move only to vacant sites.
! Only near a DC are there sites available with z > 1.
!--------------------------------------------------------------------------------
subroutine jumper2D_old(kcell,indx1,kslot1,go,kpar)
integer :: kpar,kcell,indx1(2),kslot1
logical :: go
type (cell_type), pointer :: cell
integer :: site1(3),site2(3),dcsite(3)
integer :: irel,dir1,lastdir1,indx2(2),k,z,dx,dy,irel_z,idc,nslots
integer :: savesite2(3,26)
real(DP) :: psum, p(26), R, pfac
real :: dcdist, del(3)
integer :: DC(0:DCDIM-1)     ! DC(0) = number of near DCs, DC(k) = ID of kth near DC
logical :: DCstacking, ok
real, parameter :: upfactor(2) = (/0.5, 0.3/)		! reducing prob of jumps up to z=2, z=3

cell => cellist(kcell)
site1 = cell%site
if (indx1(1) > 0 .and. indx1(2) > 0) then
	nslots = 2
else
	nslots = 1
endif
if (nslots == 1) then
	R = par_uni(kpar)
	if (R <= dirprob2D(0)) then    ! case of no jump
		go = .false.
		return
	endif
endif

! Now we must jump (if possible)
DCstacking = .false.
!if (site1(3) > 1) then		! cell is already stacked (still need to determine DC location)
!	DCstacking = .true.
!endif
DCsite = 0
DC = occupancy(site1(1),site1(2),site1(3))%DC
if (DC(0) > 0) then		! there is a DC nearby
	do idc = 1,DC(0)
		dcsite = DClist(DC(idc))%site
		del = dcsite-site1
		dcdist = norm(del)
		if (dcdist < STACKRAD_2 + 1.5) then
			DCstacking = .true.
			! Close enough to a DC for stacking to be possible.
			! This is an attribute of site1, therefore can be precomputed, and needs
			! to be recomputed only when the DC population changes (or DCs move).
			! In fact it depends only on (x,y), i.e. on site1(1) and site1(2).
			exit
		endif
	enddo
endif

lastdir1 = cell%lastdir
irel_z = nreldir2D
p = 0
psum = 0
if (site1(3) > 1) then
	if (occupancy(site1(1),site1(2),site1(3)-1)%indx(1) == 0) then
!		write(*,*) 'in space: ',kcell,site1
		irel_z = irel_z + 1		! drop cell into vacant site below with high probability
		p(irel_z) = 999
		psum = psum + p(irel_z)
		savesite2(:,irel_z) = (/ site1(1),site1(2),site1(3)-1 /)
	endif
endif
do irel = 1,nreldir2D
    p(irel) = 0
	dir1 = reldir2D(lastdir1,irel)
	site2 = site1 + jumpvec2D(:,dir1)
	if (inside_xy(site2)) then
		if (DCstacking) then	! need to look at other sites in the stack
			dx = dcsite(1) - site2(1)
			dy = dcsite(2) - site2(2)
			do k = 1,3
				z = DCstack(dx,dy,k)
				ok = .false.
				pfac = 1
				if (z > 0 .and. z == site2(3)-1) then
					if (occupancy(site2(1),site2(2),z)%indx(1) == 0) then
						ok = .true.
						pfac = upfactor(site2(3)-1)
					endif
				endif
				if (z == site2(3)+1) then
					if (occupancy(site2(1),site2(2),z)%indx(1) == 0 .and. &		! site below occupied
						occupancy(site2(1),site2(2),z-1)%indx(1) /= 0) then
						ok = .true.
					endif
				endif
				if (ok) then
					irel_z = irel_z + 1
					p(irel_z) = pfac*dirprob2D(irel)
					psum = psum + p(irel_z)
					savesite2(:,irel_z) = (/ site2(1),site2(2),z /)
				endif
			enddo
		endif
		indx2 = occupancy(site2(1),site2(2),site2(3))%indx
		if (indx2(1) == 0) then     ! not OUTSIDE_TAG or DC
			if (site2(3) > 1) then
				if (occupancy(site2(1),site2(2),site2(3)-1)%indx(1) == 0) then	! site below vacant
					cycle
				endif
			endif
			p(irel) = dirprob2D(irel)
			psum = psum + p(irel)
			savesite2(:,irel) = site2
		endif
	endif
enddo
if (psum == 0) then
	go = .false.
	cell%lastdir = random_int(1,8,kpar)
	return
else
    go = .true.
endif

! Now choose a direction on the basis of these probs p()
R = par_uni(kpar)
R = psum*R
!if (kcell == 59) then
!	write(*,*) 'R: ',R
!	write(*,'(10f6.2)') p(1:irel_z)
!endif
psum = 0
do irel = 1,irel_z
   	psum = psum + p(irel)
   	if (R <= psum) then
   		exit
   	endif
enddo
if (irel > irel_z) then
	if (irel_z == nreldir2D) then
		go = .false.
		return
	else
		irel = irel_z
	endif
endif
site2 = savesite2(:,irel)
if (irel <= nreldir2D) then
	dir1 = reldir2D(lastdir1,irel)
else	! note that the direction of a stack jump is lost here.
	dir1 = 0
endif
if (dir1 == 0) then
	dir1 = random_int(1,8,kpar)
endif

cell%site = site2
cell%lastdir = dir1
occupancy(site2(1),site2(2),site2(3))%indx(1) = kcell
if (nslots == 1) then
	if (indx1(2) > 0) then
		write(*,*) 'Error: jumper2D: cell in slot 2 only: ',site1,indx1
		stop
	endif
	occupancy(site1(1),site1(2),site1(3))%indx(1) = 0
elseif (indx1(2) == kcell) then
	occupancy(site1(1),site1(2),site1(3))%indx(2) = 0
else
	occupancy(site1(1),site1(2),site1(3))%indx(1) = indx1(2)
	occupancy(site1(1),site1(2),site1(3))%indx(2) = 0
endif
!if (site2(3) == 3) then
!	write(*,*) 'hi: ',kcell,site2,irel,occupancy(site2(1),site2(2),site2(3)-1)%indx(1)
!endif
end subroutine

!--------------------------------------------------------------------------------
! 2D motility
! In this case only one cell can occupy a site - no slot 2 - and no passing.
! This is valid because cells are not crammed together in the petri dish.
! We now wish to relax this constraint.  A site might be crowded, as a result
! of cell division, but in jumper2D cells can move only to vacant sites.
! Only near a DC are there sites available with z > 1.
!--------------------------------------------------------------------------------
subroutine jumper2D(kcell,indx1,kslot1,go,kpar)
integer :: kpar,kcell,indx1(2),kslot1
logical :: go
type (cell_type), pointer :: cell
integer :: site1(3),site2(3)
integer :: irel,dir1,lastdir1,z,irel_z,nslots
integer :: savesite2(3,26)
real(DP) :: psum, p(26), R, pfac
real, parameter :: upfactor(2) = (/0.5, 0.3/)		! reducing prob of jumps up by 1, 2 sites
real, parameter :: downfactor(2) = (/1.3, 1.5/)		! increasing prob of jumps down by 1, 2 sites

cell => cellist(kcell)
site1 = cell%site
if (indx1(1) > 0 .and. indx1(2) > 0) then
	nslots = 2
else
	nslots = 1
endif
if (nslots == 1) then
	R = par_uni(kpar)
	if (R <= dirprob2D(0)) then    ! case of no jump
		go = .false.
		return
	endif
endif

! Now we must jump (if possible)

lastdir1 = cell%lastdir
irel_z = nreldir2D
p = 0
psum = 0
if (site1(3) > 1) then
	z = site1(3) - 1
	if (occupancy(site1(1),site1(2),z)%indx(1) == 0) then
		if (zrange2D(site1(1),site1(2),1) <= z) then
!			write(*,*) 'in space: ',kcell,site1
			irel_z = irel_z + 1		! drop cell into vacant site below with high probability
			p(irel_z) = 999
			psum = psum + p(irel_z)
			savesite2(:,irel_z) = (/ site1(1),site1(2),z /)
		endif
	endif
endif
do irel = 1,nreldir2D
    p(irel) = 0
	dir1 = reldir2D(lastdir1,irel)
	site2 = site1 + jumpvec2D(:,dir1)
	if (inside_xy(site2)) then
		! Look for the lowest free site
		if (zrange2D(site2(1),site2(2),1) == 0) cycle	! outside site
		do z = zrange2D(site2(1),site2(2),1),zrange2D(site2(1),site2(2),2)
			if (occupancy(site2(1),site2(2),z)%indx(1) == 0) then
				pfac = 1.0
				if (z > site1(3)) then
					pfac = upfactor(z-site1(3))
				elseif (z < site1(3)) then
					pfac = downfactor(site1(3)-z)
				endif
				p(irel) = pfac*dirprob2D(irel)
				psum = psum + p(irel)
				savesite2(:,irel) = (/ site2(1),site2(2),z /)
				exit
			endif
		enddo
	endif
enddo
if (psum == 0) then
	go = .false.
	cell%lastdir = random_int(1,8,kpar)
	return
else
    go = .true.
endif

! Now choose a direction on the basis of these probs p()
R = par_uni(kpar)
R = psum*R
!if (kcell == 59) then
!	write(*,*) 'R: ',R
!	write(*,'(10f6.2)') p(1:irel_z)
!endif
psum = 0
do irel = 1,irel_z
   	psum = psum + p(irel)
   	if (R <= psum) then
   		exit
   	endif
enddo
if (irel > irel_z) then
	if (irel_z == nreldir2D) then
		go = .false.
		return
	else
		irel = irel_z
	endif
endif
site2 = savesite2(:,irel)
if (irel <= nreldir2D) then
	dir1 = reldir2D(lastdir1,irel)
else	! note that the direction of a stack jump is lost here.
	dir1 = 0
endif
if (dir1 == 0) then
	dir1 = random_int(1,8,kpar)
endif

cell%site = site2
cell%lastdir = dir1
occupancy(site2(1),site2(2),site2(3))%indx(1) = kcell
if (nslots == 1) then
	if (indx1(2) > 0) then
		write(*,*) 'Error: jumper2D: cell in slot 2 only: ',site1,indx1
		stop
	endif
	occupancy(site1(1),site1(2),site1(3))%indx(1) = 0
elseif (indx1(2) == kcell) then
	occupancy(site1(1),site1(2),site1(3))%indx(2) = 0
else
	occupancy(site1(1),site1(2),site1(3))%indx(1) = indx1(2)
	occupancy(site1(1),site1(2),site1(3))%indx(2) = 0
endif
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine chemo_jumper2D(kcell,indx1,kslot1,go,kpar)
integer :: kpar,kcell,indx1(2),kslot1
logical :: go
end subroutine

!--------------------------------------------------------------------------------
! Check that the DC binding info for a T cell is correct.
!--------------------------------------------------------------------------------
subroutine check_DCbound(DCbound)
integer(2) :: DCbound(2)
integer :: i, idc
logical :: err
type(DC_type) :: DC

err = .false.
do i = 1,2
    idc = DCbound(i)
    if (idc == 0) cycle
    if (NDCalive == 0) then
        write(*,*) 'check_DCbound: NDCalive = 0, T cell bound to DC: ',idc
        err = .true.
    endif
    DC = DClist(idc)
    if (.not.DC%alive) then
        write(*,*) 'check_DCbound: T cell bound to dead DC: ',idc
        err = .true.
    endif
    if (DC%nbound == 0) then
        write(*,*) 'check_DCbound: T cell bound to DC with %nbound = 0: ',idc
        err = .true.
    endif
enddo
if (err) stop
end subroutine

!--------------------------------------------------------------------------------
! Crude check for a site inside, just looking at the allowable ranges of x, y and z.
!--------------------------------------------------------------------------------
logical function inside_xyz(site)
integer :: site(3)

if (site(1) < 1 .or. site(1) > NX .or. site(2) < 1 .or. site(2) > NY .or. site(3) < 1 .or. site(3) > NZ) then
    inside_xyz = .false.
else
    inside_xyz = .true.
endif
end function

!--------------------------------------------------------------------------------
! Crude check for a site inside 2D, just looking at the allowable ranges of x, y.
!--------------------------------------------------------------------------------
logical function inside_xy(site)
integer :: site(3)

if (site(1) < 1 .or. site(1) > NX .or. site(2) < 1 .or. site(2) > NY) then
    inside_xy = .false.
else
    inside_xy = .true.
endif
end function
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine compute_Cm(tagsite,ntagged,nvar0,nvar,dt,Cm)
integer :: tagsite(3,ntagged,0:nvar)
integer :: ntagged,nvar0,nvar
real :: dt,Cm
integer :: k,j,d(3),d2,d2sum
real, allocatable :: r2mean(:)

allocate(r2mean(nvar))

do k = 1,nvar
    d2sum = 0
    do j = 1,ntagged
        d = tagsite(:,j,k) - tagsite(:,j,0)
        d2 = d(1)*d(1) + d(2)*d(2) + d(3)*d(3)
        d2sum = d2sum + d2
    enddo
    r2mean(k) = d2sum/ntagged
enddo
write(*,'(10f6.1)') r2mean
call bestfit(r2mean,nvar0,nvar,dt,Cm)
Cm = Cm*DELTA_X*DELTA_X
!write(*,*) 'Cm: ',Cm

deallocate(r2mean)

end subroutine

!--------------------------------------------------------------------------------------
! Use least squares to find the best straight line fit to r2mean() from n1 to n2.
! In fact we might need to reduce the range of points.
! We need to use a range within which the slope (dr2/dt) doesn't vary too much from
! an "average" value.  Note that the larger the persistence parameter rho the longer
! it takes for the r2 plot to become linear.
!--------------------------------------------------------------------------------------
subroutine bestfit(r2mean,n1,n2,dt,Cm)
real :: r2mean(:),dt,Cm
integer :: n1,n2
integer :: n,i
real :: sumt,sumt2,sumy,sumyt,t,y,a,b

n = n2 - n1 + 1		! number of data points

sumt = 0
sumt2 = 0
sumy = 0
sumyt = 0
do i = n1,n2
	t = i*dt
	y = r2mean(i)
	sumt = sumt + t
	sumt2 = sumt2 + t*t
	sumy = sumy + y
	sumyt = sumyt + y*t
enddo
a = (n*sumyt - sumy*sumt)/(n*sumt2 - sumt*sumt)		! slope of line
b = (sumy - a*sumt)/n								! intercept
if (IN_VITRO) then
	Cm = a/4
else
	Cm = a/6
endif
write(*,*) 'a,b: ',a,b

end subroutine

!---------------------------------------------------------------------
! Need to precompute array reldir(:,:)
! FOR NRELDIR = 6 (diagonal_jumps = .false.)
! The directions 1 - 6 are defined according to the neumann array,
! i.e. 1 is -x, 2 is +x, 3 is -y, 4 is +y, 5 is -z, 6 is +z
! The relative directions are numbered as follows:
! irel = 1 gives the same direction as lastdir
! irel = 6 gives the opposite direction to lastdir
! irel = 2,3,4,5 cover the remaining directions - order not important
! at the moment since all directions normal to lastdir are equally likely
! FOR NRELDIR = 17 (diagonal_jumps = .true.)
! The reldir(:,:) array uses the same set of 6 previous jump directions,
! restricted to the axes.  Diagonal previous jumps are accommodated by
! making a random selection of an axis direction from either 2 or 3
! possibilities (depending on whether it was a D2 or D3 jump).
! The possible jumps, now including diagonal moves, are stored in jumpvec(:)
! The jumpvec(:) entries are ordered with z varying fastest, then y, then x,
! each ranging -1, 0, +1.
! For dirprob(:) calculation:
! There are 3 groups of jumps, with sets of probs
! Group 1: Q(1) = P(1), D1 jump in the same as last direction
! Group 2: Q(2) = P(2)+P(3)+P(4)+P(5) (D2 jumps) + P(6)+P(7)+P(8)+P(9) (D3 jumps)
! Group 3: Q(3) = P(10)+P(11)+P(12)+P(13) (D1 jumps) + P(14)+P(15)+P(16)+P(17) (D2 jumps)
! i.e. Q(1) = prob of no direction change, Q(2) = prob of direction change < 90 deg,
! Q(3) = prob of direction change = 90 deg.
! Setting D1 jumps (von Neumann) to 1, the D2 jumps have length L2 = sqrt(2)
! and the D3 jumps have length L3 = sqrt(3)
! Therefore jump probs must be scaled appropriately to give symmetry within a set
! In other words:
!		P(i)/P(j) = L2/L3 for j=6,7,8,9 j=2,3,4,5
! and	P(i)/P(j) = 1/L2  for i=14,15,16,17 j=10,11,12,13
! Then choosing P(2) to represent D2 jump prob in Group 2, and P(10) to
! represent D1 jump prob in Group 3:
!		Q(2) = 4(1 + L2/L3).P(2)
!		Q(3) = 4(1 + 1/L2).P(10)
! It is always the case that Q(1) + Q(2) + Q(3) = BETA
!
! When RHO = 1, Q(1) = BETA = p1, Q(2) + Q(3) = 0
!
! When RHO = 0, P(10) = P(1), and P(2) = P(1)/L2
! therefore:	Q(3) = 4(1 + 1/L2).P(1)
! and			Q(2) = 4(1 + L2/L3).P(1)/L2
! giving:		P(1).[1 + 4(1 + L2/L3)/L2 + 4(1 + 1/L2)] = BETA
!				P(1) = BETA/[1 + 4(1 + L2/L3)/L2 + 4(1 + 1/L2)] = p0
! and:			Q(3)/Q(2) = (1 + L2)/(1 + L2/L3) = q320
!
! Now for 0 < RHO < 1
! first set Q(1) = P(1) to vary linearly with RHO from
! the value at 0 (p0) to the value at 1 (p1)
!	P(1) = p0 + RHO*(p1-p0)
! Now we know that Q(2) + Q(3) = BETA - Q(1)
! and the third parameter ALPHA is used to determine how Q(3)/Q(2)
! varies with RHO, by making it change linearly from q23 at RHO=0
! to ALPHA*q23 at RHO = 1
! Now we have all the info to solve for Q(1), Q(2), Q(3) and all P(:)
!
! Revised Moore model.
! Now allow reverse directions, but disallow D3 diagonal jumps, giving
! 18 possible directions.
! In the preferred direction the surface is a prolate spheroid, with
! a = 1, b = 1-RHO.
! In the reverse direction the surface is either a sphere (a = b) or
! an oblate spheroid (a = b^2).
!
! IN_VITRO case - 2D motion
! This case is simple, because the 8 jump directions can be ordered
! anti-clockwise (i.e. increasing theta), with the x-axis dir = 1.
! The relative directions are then found by starting at a different
! point in the sequence 1,2,..,8.
!---------------------------------------------------------------------
subroutine make_reldir
integer :: lastdir,irel,k,ix,iy,iz,i,site(3)
integer :: reldir18(6,18),reldir26(6,26)
real(DP) :: qsum,E, theta, p(8), psum, dd

if (IN_VITRO) then
	njumpdirs2D = 8
	nreldir2D = 8
	do lastdir = 1,nreldir2D
		do k = 1,njumpdirs2D
			irel = lastdir + k-1
			if (irel > 8) irel = irel - 8
			reldir2D(lastdir,k) = irel
		enddo
	enddo
	! Consider an ellipse with eccentricity = e, with major radius = a, minor radius = b,
	! and e = sqrt((a^2 - b^2)/a^2)  The distances from the negative focus at -ea to
	! the ellipse, in the jump directions, are proportional to the jump probabilities,
	! suitably weighted to account for the jump distances (1 or sqrt(2)).
	! Let p(i) = prob. of jump i
	!     d(i) = length of jump i
	!     r(i) = distance from focus to the ellipse in direction i = a(1-e^2)/(1-e.cos(theta))
	! then p(i) = k.r(i)/d(i)

	BETA = 0.95d0
	RHO = 0.95d0
	dirprob2D(0) = 1 - BETA	! this is the prob of no jump
	e = RHO
	psum = 0
	do k = 1,8
		theta = (k-1)*PI/4
		p(k) = (1 - e*e)/(1 - e*cos(theta))
		if (mod(k,2) == 0) then
			p(k) = p(k)/sqrt(2.0d0)
		endif
		psum = psum + p(k)
	enddo
	dirprob2D(1:8) = BETA*p/psum
	write(*,'(a,9f7.4)') 'dirprob2D: ',dirprob2D
	! Now set up DCstack.  This must be consistent with DCoffset()
	DCstack = 0
	! Set all DC soma sites to -1
	do k = 1,NDCsites
		site = DCoffset(:,k)
		DCstack(site(1),site(2),site(3)+1) = -1
	enddo
	do ix = -4,4
		do iy = -4,4
			dd = sqrt(real(ix)*ix + iy*iy)
			if (dd < STACKRAD_2) then
				do iz = 1,2
					if (DCstack(ix,iy,iz) == 0) then
						DCstack(ix,iy,iz) = iz
					endif
				enddo
				if (dd < STACKRAD_3) then
					if (DCstack(ix,iy,3) == 0) then
						DCstack(ix,iy,3) = 3
					endif
				endif
			endif
		enddo
	enddo
	write(*,*) 'DCstack:'
	write(nflog,*) 'DCstack:'
	do ix = -4,4
		write(*,'(i2,9(1x,3i2))') ix,((DCstack(ix,iy,iz),iz=1,3),iy=-4,4)
		write(nflog,'(i2,9(1x,3i2))') ix,((DCstack(ix,iy,iz),iz=1,3),iy=-4,4)
	enddo
	return
endif

if (MODEL == NEUMANN_MODEL) then
	diagonal_jumps = .false.
	njumpdirs = 6
	nreldir = 6
	reldir = 0
	do lastdir = 1,6
		reldir(lastdir,1) = lastdir
		if (mod(lastdir,2) == 0) then
			reldir(lastdir,6) = lastdir - 1
		else
			reldir(lastdir,6) = lastdir + 1
		endif
		irel = 2
		do k = 1,6
			if (k /= reldir(lastdir,1) .and. k /= reldir(lastdir,6)) then
				reldir(lastdir,irel) = k
				irel = irel + 1
			endif
		enddo
	enddo
    	write(*,*) 'reldir'
	    do lastdir = 1,6
		    write(*,'(6i6)') (reldir(lastdir,k),k=1,6)
	    enddo
	jumpvec(:,1:6) = neumann(:,1:6)

else
	if (MODEL == MOORE18_MODEL) then
    	nreldir = 18
	elseif (MODEL == MOORE26_MODEL) then
    	nreldir = 26
    endif
    njumpdirs = 27
	diagonal_jumps = .true.
	k = 0
	do ix = -1,1
		do iy = -1,1
			do iz = -1,1
				k = k+1
				jumpvec(:,k) = (/ix,iy,iz/)
			enddo
		enddo
	enddo

! Data for revised Moore model.  D3 jumps are excluded, but reverse directions are allowed.
!	reldir(1,:) = (/  5,  2, 4, 6, 8, 11,13,15,17, 10,12,16,18, 22,24,26,20, 23 /)	! -x
!	reldir(2,:) = (/ 23, 20,22,24,26, 11,13,15,17, 10,12,16,18,  2, 4, 6, 8,  5 /)	! +x
!	reldir(3,:) = (/ 11,  2,10,12,20,  5,13,15,23,  4, 6,22,24,  8,16,18,26, 17 /)	! -y
!	reldir(4,:) = (/ 17,  8,16,18,26,  5,13,15,23,  4, 6,22,24,  2,10,12,20, 11/)	! +y
!	reldir(5,:) = (/ 13,  4,10,16,22,  5,11,17,23,  2, 8,20,26,  6,12,18,24, 15 /)	! -z
!	reldir(6,:) = (/ 15,  6,12,18,24,  5,11,17,23,  2, 8,20,26,  4,10,16,22, 13 /)	! +z

! Data for revised Moore18 model.  D3 jumps are excluded, but reverse directions are allowed.
	reldir18(1,:) = (/  5,  2, 4, 6, 8, 11,13,15,17, 10,12,16,18, 22,24,26,20, 23 /)	! -x
	reldir18(2,:) = (/ 23, 20,22,24,26, 11,13,15,17, 10,12,16,18,  2, 4, 6, 8,  5 /)	! +x
	reldir18(3,:) = (/ 11,  2,10,12,20,  5,13,15,23,  4, 6,22,24,  8,16,18,26, 17 /)	! -y
	reldir18(4,:) = (/ 17,  8,16,18,26,  5,13,15,23,  4, 6,22,24,  2,10,12,20, 11/)	! +y
	reldir18(5,:) = (/ 13,  4,10,16,22,  5,11,17,23,  2, 8,20,26,  6,12,18,24, 15 /)	! -z
	reldir18(6,:) = (/ 15,  6,12,18,24,  5,11,17,23,  2, 8,20,26,  4,10,16,22, 13 /)	! +z

! Data for revised Moore26 model.  D3 jumps are excluded, but reverse directions are allowed.
	reldir26(1,:) = (/  5,  2, 4, 6, 8,  1, 3, 7, 9, 11,13,15,17, 10,12,16,18, 19,21,27,25, 20,22,24,26, 23 /)	! -x
	reldir26(2,:) = (/ 23, 20,22,24,26, 19,21,27,25, 11,13,15,17, 10,12,16,18,  1, 3, 7, 9,  2, 4, 6, 8,  5 /)	! +x
	reldir26(3,:) = (/ 11,  2,10,12,20,  1, 3,19,21,  5,13,15,23,  4, 6,22,24,  7, 9,25,27,  8,16,18,26, 17 /)	! -y
	reldir26(4,:) = (/ 17,  8,16,18,26,  7, 9,25,27,  5,13,15,23,  4, 6,22,24,  1, 3,19,21,  2,10,12,20, 11 /)	! +y
	reldir26(5,:) = (/ 13,  4,10,16,22,  1, 7,19,25,  5,11,17,23,  2, 8,20,26,  3, 9,21,27,  6,12,18,24, 15 /)	! -z
	reldir26(6,:) = (/ 15,  6,12,18,24,  3, 9,21,27,  5,11,17,23,  2, 8,20,26,  1, 7,19,25,  4,10,16,22, 13 /)	! +z

    if (MODEL == MOORE18_MODEL) then
        reldir(1:6,1:18) = reldir18
    elseif (MODEL == MOORE26_MODEL) then
        reldir = reldir26
    endif

endif
call compute_dirprobs
!write(logmsg,'(a,7f6.3)') dirprob(0:nreldir)
!call logger(logmsg)
qsum = 0
do i = 0,nreldir
    qsum = qsum + dirprob(i)
enddo
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine compute_dirprobs
real(DP) :: L2,L3,b,c,e,d(MAXRELDIR),R(MAXRELDIR)

if (IN_VITRO) then

	return
endif
if (MODEL == NEUMANN_MODEL) then
	diagonal_jumps = .false.
	nreldir = 6
! Old von Neumann model (no reverses)
!	dirprob(0) = 1 - beta	! this is the prob of no jump
!	dirprob(1) = (1 + 4*rho)*beta/5
!	do i = 2,5
!		dirprob(i) = (beta - dirprob(1))/4
!	enddo
!	dirprob(6) = 0			! try zero prob of a reversal of direction
! New Neumann model (with reverses)
	b = 1-RHO
	b = b*b
	e = BETA/(1 + 4*b + b**2)
	dirprob(0) = 1 - BETA	! this is the prob of no jump
	dirprob(1) = e
	dirprob(2:5) = e*b
	dirprob(6) = e*b**2
else
    ! New prolate spheroid approach
	diagonal_jumps = .true.
    if (MODEL == MOORE18_MODEL) then
	    nreldir = 18
	    L2 = sqrt(2.0d0)
	    b = 1 - RHO
	    b = b*b		! revised version
	    c = b*b		! oblate spheroid on reverse side
	    d(1) = 1
	    R(1) = 1
	    d(2:5) = L2
	    R(2:5) = L2*b/sqrt(b*b+1)
	    d(6:9) = 1
	    R(6:9) = b
	    d(10:13) = L2
	    R(10:13) = b
	    d(14:17) = L2
	    R(14:17) = L2*b*c/sqrt(b*b+c*c)
	    d(18) = 1
	    R(18) = c
	    e = BETA/(1 + 4*b/sqrt(b*b+1) + 4*b + 4*b/L2 + 4*b*c/sqrt(b*b+c*c) + c)
	    dirprob(0) = 1 - BETA
	    dirprob(1:nreldir) = e*R(1:nreldir)/d(1:nreldir)
    elseif (MODEL == MOORE26_MODEL) then
    	nreldir = 26
	    L2 = sqrt(2.0d0)
	    L3 = sqrt(3.0d0)
	    b = 1 - RHO
	    b = b*b		! revised version
	    c = b*b		! oblate spheroid on reverse side
	    d(1) = 1
	    R(1) = 1
	    d(2:5) = L2
	    R(2:5) = L2*b/sqrt(b*b+1)
	    d(6:9) = L3
	    R(6:9) = L2*b/sqrt(b*b+1)
	    d(10:13) = 1
	    R(10:13) = b
	    d(14:17) = L2
	    R(14:17) = b
	    d(18:21) = L3
	    R(18:21) = L2*b*c/sqrt(b*b+c*c)
	    d(22:25) = L2
	    R(22:25) = L2*b*c/sqrt(b*b+c*c)
	    d(26) = 1
	    R(26) = c
	    e = BETA/(1 + 4*(1+L2/L3)*b/sqrt(b*b+1) + 4*b + 4*b/L2 + 4*(1+L2/L3)*b*c/sqrt(b*b+c*c) + c)
	    dirprob(0) = 1 - BETA
	    dirprob(1:nreldir) = e*R(1:nreldir)/d(1:nreldir)
    endif
endif
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine make_probvectors
integer :: icase, i, k, dir, jump(3)
real :: vec(3), veclen, scale
real :: xang = 0.3, yang = -0.15, zang = -0.5  ! radians
real :: scale_M18 = 10  ! for M18
real :: scale_N = 6     ! for N
real :: arrow_head = 0.07

write(*,*) 'make_probvectors'
call make_reldir
if (MODEL == NEUMANN_MODEL) then
    scale = scale_N
elseif (MODEL == MOORE18_MODEL) then
    scale = scale_M18
endif
open(nfvec,file='d:\immunesystem\lattice\motilitypaper\jump_vectors\cmgui\prob.exnode',status='replace')
beta = 1.0
k = 0
do icase = 0,4
    rho = (icase-1)*0.2
    call compute_dirprobs
    write(nfvec,'(a,i1)') 'Group name : probability',icase
    write(nfvec,'(a)') ' #Fields=1'
    write(nfvec,'(a)') ' 1) vector, field, rectangular cartesian, #Components=3'
    write(nfvec,'(a)') ' x.  Value index= 1, #Derivatives= 0'
    write(nfvec,'(a)') ' y.  Value index= 2, #Derivatives= 0'
    write(nfvec,'(a)') ' z.  Value index= 3, #Derivatives= 0'
    if (icase > 0) then
        do i = 1,nreldir
            k = k+1
            dir = reldir(2,i)
            jump = jumpvec(:,dir)
            vec = jump*dirprob(i)
            call rotate(vec,xang,yang,zang)
            vec = scale*vec
            veclen = norm(vec)
            if (veclen > arrow_head) then
                vec = vec*(veclen-arrow_head)/veclen
            else
                vec = 0
            endif
            write(nfvec,'(a,i3)') 'Node: ',k
            write(nfvec,'(3f8.4)') vec
        enddo
    else
        k = k+1
        vec = (/0.1,0.0,0.0/)
        call rotate(vec,xang,yang,zang)
        vec = scale*vec
        veclen = norm(vec)
        if (veclen > arrow_head) then
            vec = vec*(veclen-arrow_head)/veclen
        else
            vec = 0
        endif
         write(nfvec,'(a,i3)') 'Node: ',k
        write(nfvec,'(3f8.4)') vec
    endif
enddo
close(nfvec)
end subroutine

!---------------------------------------------------------------------
! Rotate the vector about the X-axis by xang, then about the Y-axis by yang
!---------------------------------------------------------------------
subroutine rotate(vec,xang,yang,zang)
real :: vec(3),xang,yang,zang
real :: new(3), cosa, sina

cosa = cos(xang)
sina = sin(xang)
new(1) = vec(1)
new(2) = cosa*vec(2) - sina*vec(3)
new(3) = sina*vec(2) + cosa*vec(3)
vec = new
cosa = cos(zang)
sina = sin(zang)
new(3) = vec(3)
new(1) = cosa*vec(1) - sina*vec(2)
new(2) = sina*vec(1) + cosa*vec(2)
vec = new
cosa = cos(yang)
sina = sin(yang)
new(2) = vec(2)
new(1) = cosa*vec(1) - sina*vec(3)
new(3) = sina*vec(1) + cosa*vec(3)
vec = new

end subroutine

!---------------------------------------------------------------------
! Choose one of the 6 axis directions to allocate to the direction of
! previous step given by jump, which takes values from 0 - 27.
! This could be already on an axis, or could be a D2 or D3 diagonal,
! in which case the choice is random.
!---------------------------------------------------------------------
integer function fix_lastdir(jump,kpar)
integer :: jump, kpar
integer :: k,nax
!                            0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
integer :: naxes(0:27) = (/  0, 3, 2, 3, 2, 1, 2, 3, 2, 3, 2, 1, 2, 1, 0, 1, 2, 1, 2, 3, 2, 3, 2, 1, 2, 3, 2, 3 /)
integer :: axes(3,27) = reshape( (/ &
!   1      2      3      4      5      6      7      8      9     10     11     12     13     14
1,3,5, 1,3,0, 1,3,6, 1,5,0, 1,0,0, 1,6,0, 1,4,5, 1,4,0, 1,4,6, 3,5,0, 3,0,0, 3,6,0, 5,0,0, 0,0,0, &
!  15     16     17     18     19     20     21     22     23     24     25     26     27
6,0,0, 4,5,0, 4,0,0, 4,6,0, 2,3,5, 2,3,0, 2,3,6, 2,5,0, 2,0,0, 2,6,0, 2,4,5, 2,4,0, 2,4,6 /), (/3,27/) )

if (diagonal_jumps) then
	nax = naxes(jump)
!	if (dbug) write(nfres,*) 'fix_lastdir: jump,nax: ',nax,jump
	if (nax == 0) then
	    write(*,*) 'Should not get here: fix_lastdir: nax=0'
	    stop
		fix_lastdir = random_int(1,6,kpar)
	elseif (nax == 1) then
		fix_lastdir = axes(1,jump)
	else
		k = random_int(1,nax,kpar)
!		if (dbug) write(nfres,*) 'random_int: ',k
		fix_lastdir = axes(k,jump)
	endif
else
	write(*,*) 'fix_lastdir: Not MOORE model: ',jump
	stop
endif
end function

end module
