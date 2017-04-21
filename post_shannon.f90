program shannon
!******************************************************************************80
! SE1 is Shannon entropy computed by field value as it is.
! SE2 is Shannon entropy computed by normalized field value. Every time step is
! normalized by maximum of that time step.
!******************************************************************************80

	use global_variable

	implicit none

	integer*4, parameter :: nScale=256
	integer*4 :: i1, i2, i3, i4
	integer*4 :: rowLen, fid
	real*8 :: min, max, min_old, max_old
	real*8, dimension(0:nScale) :: ck1, ck2
	integer*4 :: check_sum1, check_sum2
	integer*4, dimension(:,:), allocatable :: hk1, hk2
	real*8 :: dScale1, dScale2
	integer*4 :: turn_on1, turn_on2
	real*8, dimension(:), allocatable :: u_node, u_normal
	real*8, dimension(:), allocatable :: SE1, SE2
	character(20) :: filename, dirname
	character(30) :: filename1, filename2
	logical :: L_exist

	write(*,'(1x,a)',advance='no') '               forward=1, TR=2: '
	read(*,*) F_or_T
	write(*,'(1x,a)',advance='no') '             simulation number: '
	read(*,*) src_sim
	write(*,'(1x,a)',advance='no') 'response(1=disp 2=vel 3=accel): '
	read(*,*) response
	print *

	if( F_or_T==FORWARD ) then
		write(dirname,'(a8,i3.3,a1)') 'forward/', src_sim, '/'
		call read_src_sim_fw()
	elseif( F_or_T==TR ) then
		write(dirname,'(a3,i3.3,a1)') 'TR/', src_sim, '/'
		call read_src_sim_TR()
	endif

	call read_model()

	allocate( u_node(nNode) )
	allocate( u_normal(nNode) )
	inquire( iolength=rowLen ) u_node

	! Find out the minimum and maximum values.
	min_old = 0d0
	max_old = 0d0
	do i1 = 1, nTstep

		if( response==DISP ) then
			write(filename,'(a9,i4.4,a4)') 'disp_step', i1, '.dat'
		elseif( response==VEL ) then
			write(filename,'(a8,i4.4,a4)') 'vel_step', i1, '.dat'
		elseif( response==ACCEL ) then
			write(filename,'(a10,i4.4,a4)') 'accel_step', i1, '.dat'
		endif

		open(newunit=fid,file=trim(adjustl(dirname))//trim(adjustl(filename)),status='old',action='read',form='unformatted',access='direct',recl=rowLen)

		read(fid,rec=1) u_node

		min = minval(u_node)
		max = maxval(u_node)

		if( min_old>min ) min_old = min
		if( max_old<max ) max_old = max

		close(fid)
	enddo

	print '(1x,a,es12.4)', 'minimum = ', min_old
	print '(1x,a,es12.4)', 'maximum = ', max_old
	print *

	! scale level 'ck'
	dScale1 = (max_old-min_old)/nScale
	dScale2 = 1d0/nScale
	ck1(0) = min_old
	ck2(0) = 0d0
	do i1 = 1, nScale
		ck1(i1) = ck1(0) + i1*dScale1
		ck2(i1) = ck2(0) + i1*dScale2
	enddo

	print '(1x,a,es12.4)', '           scale from ', ck1(0)
	print '(1x,a,es12.4)', '           scale to   ', ck1(nScale)
	print '(1x,a,es12.4)', 'normalized scale from ', ck2(0)
	print '(1x,a,es12.4)', 'normalized scale to   ', ck2(nScale)
	print *

	! histogram 'hk's
	allocate( hk1(nScale,nTstep) )
	allocate( hk2(nScale,nTstep) )
	hk1 = 0
	hk2 = 0
	do i1 = 1, nTstep

		if( mod(i1,20)==1 ) call system_clock( count1, count_rate, count_max )
		
		if( response==DISP ) then
			write(filename,'(a9,i4.4,a4)') 'disp_step', i1, '.dat'
		elseif( response==VEL ) then
			write(filename,'(a8,i4.4,a4)') 'vel_step', i1, '.dat'
		elseif( response==ACCEL ) then
			write(filename,'(a10,i4.4,a4)') 'accel_step', i1, '.dat'
		endif

		open(newunit=fid,file=trim(adjustl(dirname))//trim(adjustl(filename)),status='old',action='read',form='unformatted',access='direct',recl=rowLen)

		read(fid,rec=1) u_node
		u_normal = u_node/maxval(u_node)

		check_sum1 = 0
		check_sum2 = 0

		do i2 = 1, nNode

			turn_on1 = 0
			turn_on2 = 0
			
			do i3 = 1, nScale
				if( u_node(i2)<=ck1(i3) .and. turn_on1==0 ) then
					hk1(i3,i1) = hk1(i3,i1)+1
					check_sum1 = check_sum1+1
					turn_on1 = 1
				endif
				if( u_normal(i2)<=ck2(i3) .and. turn_on2==0 ) then
					hk2(i3,i1) = hk2(i3,i1)+1
					check_sum2 = check_sum2+1
					turn_on2 = 1
				endif
				if( turn_on1==1 .and. turn_on2==1 ) exit
			enddo
		enddo

		if( check_sum1/=nNode .or. check_sum2/=nNode ) then
			print '(1x,a,i4)', 'error in counting histogram, time step: ', i1
			print '(1x,a,i4)', ''
			stop
		endif

		close(fid)

		if( mod(i1,20)==0 ) then
			call system_clock( count2 )
			write(*,'(1x,a5,i4,a1,i4,a9,f7.4,a4)') 'step ', i1, '/', nTstep, ', elapse=', dble(count2-count1)/dble(count_rate), 'sec.'
		endif

	enddo

	! Shannon entropy 'SE'
	allocate( SE1(nTstep) )
	allocate( SE2(nTstep) )
	SE1 = 0d0
	SE2 = 0d0

	do i1 = 1, nTstep
		do i2 = 1, nScale
			if( hk1(i2,i1)/=0 ) then
				SE1(i1) = SE1(i1) - dble(hk1(i2,i1))/dble(nNode)*log(dble(hk1(i2,i1))/dble(nNode))/log(2d0)
			endif
			if( hk2(i2,i1)/=0 ) then
				SE2(i1) = SE2(i1) - dble(hk2(i2,i1))/dble(nNode)*log(dble(hk2(i2,i1))/dble(nNode))/log(2d0)
			endif
		enddo
	enddo

	! Print Shannon entropy.
	inquire( iolength=rowLen ) SE1

	if( response==DISP ) then
		write(filename1,'(a26)') 'Shannon_Entropy_disp_1.dat'
		write(filename2,'(a26)') 'Shannon_Entropy_disp_2.dat'
	elseif( response==VEL ) then
		write(filename1,'(a25)') 'Shannon_Entropy_vel_1.dat'
		write(filename2,'(a25)') 'Shannon_Entropy_vel_2.dat'
	elseif( response==ACCEL ) then
		write(filename1,'(a27)') 'Shannon_Entropy_accel_1.dat'
		write(filename2,'(a27)') 'Shannon_Entropy_accel_2.dat'
	endif

	open(newunit=fid,file=trim(adjustl(dirname))//trim(adjustl(filename1)),status='replace',action='write',form='unformatted',access='direct',recl=rowLen)
	write(fid,rec=1) SE1
	close(fid)

	inquire( iolength=rowLen ) SE2
	open(newunit=fid,file=trim(adjustl(dirname))//trim(adjustl(filename2)),status='replace',action='write',form='unformatted',access='direct',recl=rowLen)
	write(fid,rec=1) SE2
	close(fid)

	deallocate( u_node, u_normal )
	deallocate( hk1, hk2 )
	deallocate( SE1, SE2 )

contains

subroutine read_src_sim_fw()

	implicit none

	integer*4 :: i1, fid
	character(30) :: filename

	write(filename,'(a8,i3.3,a19)') 'forward/', src_sim, '/simulation_log.txt'
	open(newunit=fid,file=filename,status='old',action='read')

	do i1 = 1, 3; read(fid,*); enddo
	read(fid,'(43x,i12)') model_num

	do i1 = 1, 6; read(fid,*); enddo
	read(fid,'(43x,f12.4)') dt
	read(fid,'(43x,i12)') nTstep

	close(fid)

	print '(1x,a)', '===== Simulation Info ====='
	print '(4x,a14,i12)', 'model number: ', model_num
	print '(4x,a14,f12.4)', 'dt: ', dt
	print '(4x,a14,i12)', 'nTstep: ', nTstep
	print *

end subroutine read_src_sim_fw

subroutine read_src_sim_TR()

	implicit none

	integer*4 :: i1, fid
	character(25) :: filename

	write(filename,'(a3,i3.3,a19)') 'TR/', src_sim, '/simulation_log.txt'
	open(newunit=fid,file=filename,status='old',action='read')

	do i1 = 1, 4; read(fid,*); enddo
	read(fid,'(43x,i12)') model_num

	do i1 = 1, 6; read(fid,*); enddo
	read(fid,'(43x,f12.4)') dt
	read(fid,'(43x,i12)') nTstep

	close(fid)

	print '(1x,a)', '===== Simulation Info ====='
	print '(4x,a14,i12)', 'model number: ', model_num
	print '(4x,a14,f12.4)', 'dt: ', dt
	print '(4x,a14,i12)', 'nTstep: ', nTstep
	print *

end subroutine read_src_sim_TR

end program