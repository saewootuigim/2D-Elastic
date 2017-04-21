program post_fomel

	implicit none

	integer*4 :: i1, i2
	integer*4 :: fid1, fid2, row_len
	integer*4 :: F_or_T
	integer*4 :: model_num
	integer*4 :: nDOF_disp
	integer*4 :: nTstep
	real*8    :: numerator(2), denominator(2)
	real*8, dimension(:), allocatable :: u, focus1, focus2
	character(4) :: sim_num
	character(40) :: dirname, model_dir

	write(*,'(1x,a19)',advance='no') '  1=forward, 2=TR: '
	read(*,*) F_or_T
	write(*,'(1x,a19)',advance='no') 'simulation number: '
	read(*,*) sim_num

	if( F_or_T==1 ) then
		dirname = 'forward/'//trim(sim_num)//'/'
	elseif( F_or_T==2 ) then
		dirname = 'TR/'//trim(sim_num)//'/'
	endif

	print '(1x,a,a)', 'reading simulation info from ', trim(dirname)
	open(newunit=fid1,file=trim(dirname)//'constants.txt',status='old',action='read')
	if( F_or_T==1 ) then
		read(fid1,*) model_num
		read(fid1,*);read(fid1,*) nTstep
	elseif( F_or_T==2 ) then
		read(fid1,*);read(fid1,*) model_num
		read(fid1,*);read(fid1,*) nTstep
	endif
	close(fid1)

	write(model_dir,'(a6,i1,a1)') 'model/', model_num, '/'
	print '(1x,a,a)', 'reading model info from ', trim(model_dir)
	open(newunit=fid1,file=trim(model_dir)//'constants.txt',status='old',action='read')
	do i1 = 1, 7; read(fid1,*); enddo; read(fid1,*) nDOF_disp
	close(fid1)

	allocate( u(nDOF_disp) )
	inquire( iolength=row_len) u
	open(newunit=fid1,file=trim(dirname)//'disp.dat',status='old',action='read',form='unformatted',access='direct',recl=row_len)

	print '(1x,a,i8)', '               model number: ', model_num
	print '(1x,a,i8)', '                 time steps: ', nTstep
	print '(1x,a,i8)', 'number of displacement DOFs: ', nDOF_disp

    allocate( focus1(nTstep) )
	allocate( focus2(nTstep) )

	do i1 = 1, nTstep
		read(fid1,rec=i1) u
		numerator = 0d0
		denominator = 0d0
		do i2 = 1, nDOF_disp
			numerator(1)   = numerator(1)   + abs(u(i2))
			denominator(1) = denominator(1) + u(i2)**2
			numerator(2)   = numerator(2)   + u(i2)**4
			denominator(2) = denominator(2) + u(i2)**2
		enddo
		focus1(i1) = (numerator(1)/sqrt(denominator(2))-1d0)/(sqrt(dble(nDOF_disp))-1d0)
		focus2(i1) = (nDOF_disp*numerator(2)/denominator(2)**2-1d0)/(dble(nDOF_disp)-1d0)
	enddo

	close(fid1)

	inquire( iolength=row_len) focus1
	open(newunit=fid1,file=trim(dirname)//'focus_Fomel1.dat',status='replace',action='write',form='unformatted',access='direct',recl=row_len)
	open(newunit=fid2,file=trim(dirname)//'focus_Fomel2.dat',status='replace',action='write',form='unformatted',access='direct',recl=row_len)
	write(fid1,rec=1) focus1
	write(fid2,rec=1) focus2
	close(fid1)
	close(fid2)

	deallocate( u, focus1, focus2 )

end program