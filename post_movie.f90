program movie

	use global_variable

	implicit none

	integer*4 :: fid1, fid2, rowLen
	integer*4 :: i1, i2, i3
	real*8, allocatable, dimension(:) :: u, u_node
	character(20) :: filename, dirname

	! Read simulation and model info.
	write(*,'(1x,a)',advance='no') '               forward=1, TR=2: '
	read(*,*) F_or_T
	write(*,'(1x,a)',advance='no') '                directory name: '
	read(*,*) dirname
	write(*,'(1x,a)',advance='no') 'response(1=disp 2=vel 3=accel): '
	read(*,*) response

	if( F_or_T==FORWARD ) then
		dirname = 'forward/'//trim(dirname)//'/'
		call read_src_sim_fw()
	elseif( F_or_T==TR ) then
		dirname = 'TR/'//trim(dirname)//'/'
		call read_src_sim_TR()
	endif

	call read_model()

	! Open output file.
	allocate( u(nDOF_disp) )
	inquire( iolength=rowLen ) u
	if( response==DISP ) then
		open(newunit=fid1,file=trim(dirname)//'disp.dat',status='old',action='read',form='unformatted',access='direct',recl=rowLen)
	elseif( response==VEL ) then
		open(newunit=fid1,file=trim(dirname)//'vel.dat',status='old',action='read',form='unformatted',access='direct',recl=rowLen)
	elseif( response==ACCEL ) then
		open(newunit=fid1,file=trim(dirname)//'accel.dat',status='old',action='read',form='unformatted',access='direct',recl=rowLen)
	endif

	! Start making movie.
	allocate( u_node(nNode) )
	inquire( iolength=rowLen ) u_node

	! Simulation info.
	print '(1x,a)', '*** simulation info ***'
	print '(1x,a12,i7)', '   nTstep = ', nTstep
	print '(1x,a,f7.3)', '       dt = ', dt
	print '(1x,a12,i7)', '    nNode = ', nNode
	print '(1x,a12,i7)', 'nDOF_disp = ', nDOF_disp
	print *

	do i2 = 1, nTstep

		read(fid1,rec=i2) u

		if( response==DISP ) then
			write(filename,'(a9,i4.4,a4)'), 'disp_step', i2, '.dat'
		elseif( response==VEL ) then
			write(filename,'(a9,i4.4,a4)'), 'vel_step', i2, '.dat'
		elseif( response==ACCEL ) then
			write(filename,'(a10,i4.4,a4)'), 'accel_step', i2, '.dat'
		endif

		open(newunit=fid2,file=trim(dirname)//trim(adjustl(filename)),status='replace',action='write',form='unformatted',access='direct',recl=rowLen)

		u_node = 0d0

		do i1 = 1, nDOF_disp
			i3 = DOF_disp(i1)
			u_node( DOF2node(i3) ) = u_node( DOF2node(i3) ) + u(i1)**2
		enddo

		do i1 = 1, nNode
			u_node(i1) = sqrt(u_node(i1))
		enddo

		write(fid2,rec=1) u_node

		close( fid2 )

	enddo

	close( fid1 )

	deallocate( u )
	deallocate( u_node )

contains

subroutine read_src_sim_fw()

	implicit none

	integer*4 :: i1, fid
	character(40) :: filename

	filename = trim(dirname)//'simulation_log.txt'
	open(newunit=fid,file=trim(filename),status='old',action='read')

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
	character(40) :: filename

	filename = trim(dirname)//'simulation_log.txt'
	open(newunit=fid,file=trim(filename),status='old',action='read')

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

end program movie