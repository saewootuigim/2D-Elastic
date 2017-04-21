subroutine create_output_directory()

    use global_variable
    
    implicit none

    integer*4 :: i1
    character(1) :: asdf
    character(4) :: simulation_number
    character(12) :: dirname
    logical :: L_exist

    i1 = 0

    ! Get simulation number from user.
    write(*,'(1x,a)') 'Name this simulation.'
    if( F_or_T==FORWARD ) then
        write(*,'(4x,a8)',advance='no') 'forward/'
        read(*,'(a4)') simulation_number
        dirname = 'forward/'//trim(simulation_number)//'/'
    elseif( F_or_T==TR ) then
        write(*,'(4x,a3)',advance='no') 'TR/'
        read(*,'(a4)') simulation_number
        dirname = 'TR/'//trim(simulation_number)//'/'
    endif

    ! Check if that simulation number already exists.
    inquire( directory=trim(dirname), exist=L_exist )
    if( L_exist==.true. ) then
        write(*,'(4x,a,a,a)',advance='no') 'directory ', trim(dirname), ' already exists. Do you want to remove it?(y/n)'
        read(*,*) asdf
        if( asdf=='y' .or. asdf=='Y' ) then
            call system( 'rm -r '//trim(dirname) )
        else
            print '(4x,a)', 'program stops.'
            stop
        endif
    endif
    
    ! Create directory.
    call system( 'mkdir -p '//trim(dirname) )

    dir_sim_out = trim(dirname)

end subroutine create_output_directory




subroutine write_fw_simulation_info()

    use global_variable
    
    implicit none

    integer*4 :: i1, fid
    character(8)  :: date
    character(10) :: time

    call DATE_AND_TIME( date, time )
    call system_clock( count2 )

    open(newunit=fid,file=trim(dir_sim_out)//'simulation_log.txt',status='replace',action='write')

    write(fid,101) '********************************************************************************'
    write(fid,101) '                                                                         GENERAL'
    write(fid,101) '********************************************************************************'
    write(fid,103) '             model: ', model_num
    write(fid,104) '            elapse: ', dble(count2-count1)/dble(count_rate), ' sec.'
    write(fid,105) '       finished at: ', date(1:4)//'/'//date(5:6)//'/'//date(7:8), time(1:2)//':'//time(3:4)//':'//time(5:6)
    write(fid,*)
    write(fid,101) '********************************************************************************'
    write(fid,101) '                                                                            TIME'
    write(fid,101) '********************************************************************************'
    write(fid,102) '                dt: ', dt
    write(fid,103) '            nTstep: ', nTstep
    write(fid,102) ' simulation length: ', Tsim
    write(fid,*)
    write(fid,101) '********************************************************************************'
    write(fid,101) '                                                                            LOAD'
    write(fid,101) '********************************************************************************'
    if( fw_src_sgnl==RICKER ) then
        write(fid,106) '     source signal: ', 'Ricker'
        write(fid,102) ' central frequency: ', f0
        write(fid,102) '  loading duration: ', Tload
    elseif( fw_src_sgnl==GAUSS ) then
        write(fid,106) '     source signal: ', 'Gauss'
        write(fid,102) '           average: ', 0.1
        write(fid,102) '  loading duration: ', Tload
    elseif( fw_src_sgnl==RECT ) then
        write(fid,106) '     source signal: ', 'Rectangle'
        write(fid,102) '            offset: ', f0
        write(fid,102) '  loading duration: ', Tload
    elseif( fw_src_sgnl==DGAUSS ) then
        write(fid,106) '     source signal: ', 'derivative of Gaussian'
        write(fid,102) '           average: ', 0.1
        write(fid,102) '  loading duration: ', Tload
    endif
    write(fid,103) ' number of loading nodes: ', fw_nPntLd
    do i1 = 1, fw_nPntLd
        write(fid,103) '       node: ', fw_load_node(i1)
    enddo
    do i1 = 1, fw_nPntLd
        write(fid,107) ' coordinate: ', fw_load_coord(i1,:)
    enddo
    do i1 = 1, fw_nPntLd
        write(fid,108) '        DOF: ', fw_load_DOF(i1,:)
    enddo

    close(fid)
    
    open(newunit=fid,file=trim(dir_sim_out)//'constants.txt',status='replace',action='write')
    write(fid,'(i12)') model_num
    write(fid,'(f12.4)') dt
    write(fid,'(i12)') nTstep
    write(fid,'(f12.4)') Tsim
    write(fid,'(i12)') fw_nPntLd
    do i1 = 1, fw_nPntLd
        write(fid,'(2f12.4)') fw_load_coord(i1,:)
    enddo
    close(1)

101 format(a80)
102 format(a43,f12.4)
103 format(a43,i12)
104 format(a43,f12.4,a6)
105 format(a43,a10,1x,a8)
106 format(a43,a12)
107 format(a43,2f12.4)
108 format(a43,2i12)

end subroutine write_fw_simulation_info




subroutine write_TR_simulation_info()

    use global_variable
    
    implicit none

    integer*4 :: i1, fid
    character(8)  :: date
    character(10) :: time

    call DATE_AND_TIME( date, time )
    call system_clock( count2 )

    open(newunit=fid,file=trim(dir_sim_out)//'simulation_log.txt',status='replace',action='write')

    write(fid,101) '********************************************************************************'
    write(fid,101) '                                                                         GENERAL'
    write(fid,101) '********************************************************************************'
    write(fid,103) '        source simulation: ', src_sim
    write(fid,103) '                    model: ', model_num
    write(fid,104) '                   elapse: ', dble(count2-count1)/dble(count_rate), ' sec.'
    write(fid,105) '              finished at: ', date(1:4)//'/'//date(5:6)//'/'//date(7:8), time(1:2)//':'//time(3:4)//':'//time(5:6)
    write(fid,*)
    write(fid,101) '********************************************************************************'
    write(fid,101) '                                                                            TIME'
    write(fid,101) '********************************************************************************'
    write(fid,102) '                       dt: ', dt
    write(fid,103) '                   nTstep: ', nTstep
    write(fid,102) '        simulation length: ', Tsim
    write(fid,*)
    write(fid,101) '********************************************************************************'
    write(fid,101) '                                                                            LOAD'
    write(fid,101) '********************************************************************************'
    if( TR_surf_trac_resp==DISP ) then
        write(fid,106) 'surface response used: ', 'displacement'
    elseif( TR_surf_trac_resp==VEL ) then
        write(fid,106) 'surface response used: ', '    velocity'
    elseif( TR_surf_trac_resp==ACCEL ) then
        write(fid,106) 'surface response used: ', 'acceleration'
    endif
    if( TR_apply_sink==YES ) then
        write(fid,106) 'was load applied at sink?: ', 'yes'
    elseif( TR_apply_sink==NO ) then
        write(fid,106) 'was load applied at sink?: ', 'no'
    endif
    if( fw_src_sgnl==RICKER ) then
        write(fid,106) 'What was the load?: ', 'Ricker'
    elseif( fw_src_sgnl==GAUSS ) then
        write(fid,106) 'What was the load?: ', 'Gauss'
    endif
    if( TR_apply_surf_type==1 ) then
        write(fid,106) 'was surface traction applied?: ', 'no'
    elseif( TR_apply_surf_type==2 ) then
        write(fid,106) 'was surface traction applied?: ', 'yes, normalized'
    elseif( TR_apply_surf_type==3 ) then
        write(fid,106) 'was surface traction applied?: ', 'yes, rescaled'
    endif
    if( TR_DtN==NO ) then
        write(fid,106) 'boundary condition: ', 'Dirichlet'
    elseif( TR_DtN==YES ) then
        write(fid,106) 'boundary condition: ', 'Nuemann'
    endif

    close(fid)

    open(newunit=fid,file=trim(dir_sim_out)//'constants.txt',status='replace',action='write')
    write(fid,'(i12)') src_sim
    write(fid,'(i12)') model_num
    write(fid,'(f12.4)') dt
    write(fid,'(i12)') nTstep
    write(fid,'(f12.4)') Tsim
    close(fid)

101 format(a80)
102 format(a43,f12.4)
103 format(a43,i12)
104 format(a43,f12.4,a5)
105 format(a43,a10,1x,a8)
106 format(a43,a15)

end subroutine write_TR_simulation_info




subroutine normalize_surface_output()

    use global_variable

    implicit none

    integer*4 :: i1, fid, rowLen
    real*8 :: u_max, u_max_old
    real*8, dimension(nDOF_disp_surf) :: u

    print '(1x,a)', 'normalizing surface response'
    print *

    u_max_old = 0d0

    inquire( iolength=rowLen ) u

    if( response==DISP ) then
        open(newunit=fid,file=trim(adjustl(dir_sim_out))//'disp_surf.dat', status='old',action='readwrite',form='unformatted',access='direct',recl=rowLen)
    elseif( response==VEL ) then
        open(newunit=fid,file=trim(adjustl(dir_sim_out))//'vel_surf.dat',  status='old',action='readwrite',form='unformatted',access='direct',recl=rowLen)
    elseif( response==ACCEL ) then
        open(newunit=fid,file=trim(adjustl(dir_sim_out))//'accel_surf.dat',status='old',action='readwrite',form='unformatted',access='direct',recl=rowLen)
    endif

    do i1 = 1, nTstep
        read(fid,rec=i1) u
        u_max = maxval( abs(u) )
        if( u_max>u_max_old ) u_max_old = u_max
    enddo

    print '(1x,a)', 'Surface output is normalized.'
    if( response==DISP ) then
        print '(1x,a,es12.4)', 'displacement surface response MAX ', u_max_old
    elseif( response==VEL ) then
        print '(1x,a,es12.4)', 'velocity surface response MAX ', u_max_old
    elseif( response==ACCEL ) then
        print '(1x,a,es12.4)', 'acceleration surface response MAX ', u_max_old
    endif

    do i1 = 1, nTstep
        read(fid,rec=i1) u
        u = u/u_max_old
        write(fid,rec=i1) u
    enddo

    close( fid )

    print *

end subroutine normalize_surface_output