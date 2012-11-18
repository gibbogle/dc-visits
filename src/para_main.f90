!-----------------------------------------------------------------------------------------
! Runs omp_para by calls to libpara32.dll or libpara64.dll 
!-----------------------------------------------------------------------------------------
#define COMPILE_THIS

PROGRAM omp_para
use omp_main_mod
use omp_global
integer :: ncpu, res, summarydata(100)
character*(128) :: infile,outfile
character*(64) :: travelfile = 'travel_time_dist.out'
integer :: status, nlen, cnt, i, inbuflen, outbuflen
integer :: jstep, hour, ntot, ncog(2), inflow, exits
character*(128) :: b, c, progname

call process_command_line(ncpu,infile,outfile)

#ifdef COMPILE_THIS

outfile = 'para.out'
!resfile = 'para.res'
inbuflen = len(infile)
outbuflen = len(outfile)

call get_command (b, nlen, status)
if (status .ne. 0) then
    write (*,*) 'get_command failed with status = ', status
    stop
end if
!write (*,*) 'command line = ', b(1:len)
call get_command_argument (0, c, nlen, status)
if (status .ne. 0) then
    write (*,*) 'Getting command name failed with status = ', status
    stop
end if
!write (*,*) 'command name = ', c(1:len)
progname = c(1:nlen)
cnt = command_argument_count ()
!write (*,*) 'number of command arguments = ', cnt
if (cnt < 2) then
    write(*,*) 'Use: ',trim(progname),' num_cpu input_file'
    stop
endif

do i = 1, cnt
    call get_command_argument (i, c, nlen, status)
    if (status .ne. 0) then
        write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', i
        stop
    end if
    if (i == 1) then
!        read(c(1:len),'(i)') ncpu
        read(c(1:nlen),*) ncpu															! --> ncpu
        write(*,*) 'Requested threads: ',ncpu
    elseif (i == 2) then
        infile = c(1:nlen)																! --> infile
        write(*,*) 'Input file: ',infile
    elseif (i == 3) then
        outfile = c(1:nlen)																! --> outfile
        write(*,*) 'Output file: ',outfile
!    elseif (i == 4) then
!        resfile = c(1:len)																! --> resfile
!        write(*,*) 'Result file: ',resfile
    endif
end do

!ncpu = 1
!infile = 'expta3.inp'

#endif

!runfile = 'running.out'
if (compute_travel_time) then
	N_TRAVEL_COG = 1
	N_TRAVEL_DC = 20
	N_TRAVEL_DIST = 50
	allocate(travel_cog(N_TRAVEL_COG))
	allocate(travel_dc(N_TRAVEL_DC))
	allocate(travel_dist(N_TRAVEL_COG,N_TRAVEL_DC,N_TRAVEL_DIST))
! Test values
	travel_cog(1) = 0.0001
!	travel_cog(2) = 0.001
	travel_dc(1) = 200
	travel_dc(2) = 300
	travel_dc(3) = 400
	travel_dc(4) = 500
	travel_dc(5) = 600
	travel_dc(6) = 700
	travel_dc(7) = 800
	travel_dc(8) = 900
	travel_dc(9) = 1000
	travel_dc(10) = 1100
	travel_dc(11) = 1200
	travel_dc(12) = 1300
	travel_dc(13) = 1400
	travel_dc(14) = 1500
	travel_dc(15) = 1600
	travel_dc(16) = 1700
	travel_dc(17) = 1800
	travel_dc(18) = 1900
	travel_dc(19) = 2000
	travel_dc(20) = 2100

	do k_travel_cog = 1, N_TRAVEL_COG
		do k_travel_dc = 1,N_TRAVEL_DC
			call execute(ncpu,infile,inbuflen,outfile,outbuflen) 
			if (ntravel > 0) then
				travel_dist(k_travel_cog,k_travel_dc,:) = travel_dist(k_travel_cog,k_travel_dc,:)/ntravel
			endif
		enddo
	enddo
	open(nftravel,file=travelfile,status='replace')
	do k_travel_cog = 1, N_TRAVEL_COG
		do k_travel_dc = 1,N_TRAVEL_DC
			write(nftravel,'(2i4,50f7.4)') k_travel_cog,k_travel_dc,travel_dist(k_travel_cog,k_travel_dc,:)
		enddo
	enddo
	close(nftravel)
	stop
endif

write(*,*) 'call execute'
call execute(ncpu,infile,inbuflen,outfile,outbuflen)
!call get_dimensions(NX,NY,NZ,Nsteps)
do jstep = 1,Nsteps
	call simulate_step(res)
	if (res < 0) then
		write(*,*) 'Error exit'
		stop
	endif
	if (res > 0) then
		write(*,*) 'Successful execution'
		exit
	endif
	if (mod(jstep,240) == 0) then
		call get_summary(summarydata)
!summaryData(1:13) = (/int(tnow/60),istep,globalvar%NDCalive,nact,ntot,ncogseed,ncog,Ndead,	nbnd,int(globalvar%InflowTotal),globalvar%Nexits, teffgen/)
		hour = summaryData(1)
		ntot = summaryData(5)
		ncog(1:2) = summaryData(7:8)
		inflow = summaryData(11)
		exits = summaryData(12)
		write(*,'(5(a,i6))') 'Hour: ',hour,' ncells: ',ntot,' ncog: ',ncog(1),' inflow: ',inflow,' nexits: ', exits		
	endif
enddo
call terminate_run(0)
end

