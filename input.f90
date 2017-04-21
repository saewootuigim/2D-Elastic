subroutine read_input_file()

    use global_variable

    implicit none

    integer*4 :: i1, fid
    logical :: L_exist

    inquire( file='input.txt', exist=L_exist )
    if( L_exist==.false. ) then
        print '(a)', 'ERROR *******************************************************************'
        print '(a)', 'Input file does not exist.'
        stop
    endif

    open(newunit=fid,file='input.txt',status='old')

    read(fid,*)
    read(fid,*) NQuad
    allocate( xi( NQuad ) )
    allocate( wi( NQuad ) )

    read(fid,*)
    read(fid,*)
    read(fid,*) quad_meth

    read(fid,*)
    read(fid,*)
    read(fid,*) model_num

    read(fid,*)
    read(fid,*)
    read(fid,*) F_or_T

    read(fid,*)
    read(fid,*)
    read(fid,*) dt

    read(fid,*)
    read(fid,*)
    read(fid,*) Tsim

    nTstep = nint(Tsim/dt)

    read(fid,*)
    read(fid,*)
    read(fid,*) response

    read(fid,*)
    read(fid,*)
    read(fid,*) fw_nPntLd

    allocate( fw_load_node ( fw_nPntLd   ) )
    allocate( fw_load_DOF  ( fw_nPntLd,2 ) )
    allocate( fw_load_coord( fw_nPntLd,2 ) )

    read(fid,*)
    read(fid,*)
    read(fid,*) fw_PntLd_dir

    read(fid,*)
    read(fid,*)
    do i1 = 1, fw_nPntLd
        read(fid,*) fw_load_coord(i1,:)
    enddo

    read(fid,*)
    read(fid,*)
    read(fid,*) fw_src_sgnl

    read(fid,*)
    read(fid,*)
    read(fid,*) f0

    read(fid,*)
    read(fid,*)
    read(fid,*) Tload
    
    read(fid,*)
    read(fid,*)
    read(fid,*) TR_DtN

    if( TR_DtN==NO ) then
        write(*,'(a10)') 'DtN: NO, Dirichlet boundary condition will apply.'
    elseif( TR_DtN==YES ) then
        write(*,'(a10)') 'DtN: YES, Neumann boundary condition will apply.'
    endif

    read(fid,*)
    read(fid,*)
    read(fid,*) src_sim

    write(dir_sim_in,'(a8,i3.3,a1)') 'forward/', src_sim, '/'

    read(fid,*)
    read(fid,*)
    read(fid,*) TR_surf_trac_case

    read(fid,*)
    read(fid,*)
    read(fid,*) TR_surf_trac_dir

    read(fid,*)
    read(fid,*)
    read(fid,*) TR_surf_trac_resp

    read(fid,*)
    read(fid,*)
    read(fid,*) a0, b0

    print '(1x,a,f8.4,a,f8.4)', 'PML constants: a0=', a0, ', b0=', b0

    read(fid,*)
    read(fid,*)
    read(fid,*) beta, gamma

    print '(1x,a,f8.4,a,f8.4)', 'Newmark constants: beta=', beta, ', gamma=', gamma

    read(fid,*)
    read(fid,*)
    read(fid,*) fw_load_intensity

    read(fid,*);read(fid,*);read(fid,*)
    read(fid,*) TR_apply_sink

    read(fid,*);read(fid,*);read(fid,*)
    read(fid,*) TR_apply_surf_type
    
    close(fid)

end subroutine read_input_file




subroutine read_model()

    use global_variable

    implicit none

    integer*4 :: i1, fid
    integer*4 :: fileLen
    character(len=8) :: dirname

    print *
    print '(1x,a,i1)', 'reading model ', model_num

    write(dirname,'(a6,i1,a1)') 'model/', model_num, '/'

    ! Read model constants.
    open(newunit=fid,file=dirname//'constants.txt',status='old',action='read')

    read(fid,'(i12)') model_num
    read(fid,'(i12)') nNode
    read(fid,'(i12)') nElem
    read(fid,'(f12.4)') pml_l
    read(fid,'(f12.4)') rdwd
    read(fid,'(f12.4)') rdht
    read(fid,'(i12)') nDOF
    read(fid,'(i12)') nDOF_disp
    read(fid,'(i12)') nDOF_disp_surf
    read(fid,'(i12)') nDOF_disp_inc
    read(fid,'(i12)') nDOF_disp_RD
    read(fid,'(i12)') nDOF_disp_PML
    read(fid,'(i12)') nDOF_stress

    close( fid )
    
    ! Read material property.
    open(newunit=fid,file=dirname//'parameters_PMK.txt',status='old',action='read')

    do i1 = 1, 49; read(fid,*); enddo
    read(fid,*) nPhyProp
    
    allocate(   E( nPhyProp ) )
    allocate( rho( nPhyProp ) )
    allocate(  nu( nPhyProp ) )

    read(fid,*);read(fid,*)
    do i1 = 1, nPhyProp
        read(fid,*) E(i1)
    enddo

    read(fid,*);read(fid,*)
    do i1 = 1, nPhyProp
        read(fid,*) rho(i1)
    enddo

    read(fid,*);read(fid,*)
    do i1 = 1, nPhyProp
        read(fid,*) nu(i1)
    enddo

    close(fid)

    ! Read other info.
    allocate(   elem2node( nElem, 8 ) )
    allocate(    node2DOF( nNode, 5 ) )
    allocate(     node2xy( nNode, 2 ) )
    allocate(    DOF2node( nDOF ) )
    allocate(    DOF_surf( nDOF_disp_surf ) )
    allocate(    DOF_disp( nDOF_disp ) )
    allocate( node2region( nNode ) )
    allocate( elem2region( nElem ) )
    allocate( elem2phyprp( nElem ) )

    inquire(iolength=fileLen) elem2node
    open(newunit=fid,file=dirname//'elem2node.dat',status='old',action='read',form='unformatted',access='direct',recl=fileLen)
    read(fid,rec=1) elem2node
    close(fid)

    inquire(iolength=fileLen) node2DOF
    open(newunit=fid,file=dirname//'node2DOF.dat',status='old',action='read',form='unformatted',access='direct',recl=fileLen)
    read(fid,rec=1) node2DOF
    close(fid)

    inquire(iolength=fileLen) node2xy
    open(newunit=fid,file=dirname//'node2xy.dat',status='old',action='read',form='unformatted',access='direct',recl=fileLen)
    read(fid,rec=1) node2xy
    close(fid)

    inquire(iolength=fileLen) DOF2node
    open(newunit=fid,file=dirname//'DOF2node.dat',status='old',action='read',form='unformatted',access='direct',recl=fileLen)
    read(fid,rec=1) DOF2node
    close(fid)

    inquire(iolength=fileLen) DOF_surf
    open(newunit=fid,file=dirname//'DOF_surf.dat',status='old',action='read',form='unformatted',access='direct',recl=fileLen)
    read(fid,rec=1) DOF_surf
    close(fid)

    inquire(iolength=fileLen) DOF_disp
    open(newunit=fid,file=dirname//'DOF_disp.dat',status='old',action='read',form='unformatted',access='direct',recl=fileLen)
    read(fid,rec=1) DOF_disp
    close(fid)

    inquire(iolength=fileLen) node2region
    open(newunit=fid,file=dirname//'node2region.dat',status='old',action='read',form='unformatted',access='direct',recl=fileLen)
    read(fid,rec=1) node2region
    close(fid)

    open(newunit=fid,file=trim(dirname)//'elem2region.txt',status='old',action='read')
    read(fid,'(<nElem>i2)') elem2region
    close(fid)

    open(newunit=fid,file=trim(dirname)//'elem2phyprp.txt',status='old',action='read')
    read(fid,'(<nElem>i2)') elem2phyprp
    close(fid)

    print *

end subroutine read_model




subroutine read_src_sim()

    use global_variable

    implicit none

    integer*4 :: fid, i1

    open(newunit=fid,file=trim(adjustl(dir_sim_in))//'simulation_log.txt',status='old',action='read')

    do i1 = 1, 3; read(fid,*); enddo
    read(fid,'(43x,i12)') model_num

    do i1 = 1, 6; read(fid,*); enddo
    read(fid,'(43x,f12.4)') dt
    read(fid,'(43x,i12)') nTstep_src

    do i1 = 1, 8; read(fid,*); enddo
    read(fid,'(43x,i12)') fw_nPntLd

    if( allocated(fw_load_node)==.true. ) then
        deallocate( fw_load_node )
        allocate( fw_load_node(fw_nPntLd) )
    endif
    if( allocated(fw_load_coord)==.true. ) then
        deallocate( fw_load_coord )
        allocate( fw_load_coord(fw_nPntLd,2) )
    endif
    if( allocated(fw_load_DOF)==.true. ) then
        deallocate( fw_load_DOF )
        allocate( fw_load_DOF(fw_nPntLd,2) )
    endif

    do i1 = 1, fw_nPntLd
        read(fid,'(43x,i12)') fw_load_node(i1)
    enddo
    do i1 = 1, fw_nPntLd
        read(fid,'(43x,2f12.4)') fw_load_coord(i1,:)
    enddo
    do i1 = 1, fw_nPntLd
        read(fid,'(43x,2i12)') fw_load_DOF(i1,:)
    enddo

    print '(1x,a)', 'source simulation info'
    print '(1x,a23,i12)', 'simulation number: ', src_sim
    print '(1x,a23,i12)', 'model number: ', model_num
    print '(1x,a23,i12)', 'total time steps: ', nTstep_src
    print '(1x,a23,f12.4)', 'dt: ', dt
    print '(1x,a23,i12)', 'number of point loads: ', fw_nPntLd
    do i1 = 1, fw_nPntLd
        print '(1x,a23,i12)', 'node: ', fw_load_node(i1)
    enddo
    do i1 = 1, fw_nPntLd
        print '(1x,a23,2f12.4)', 'coordinate: ', fw_load_coord(i1,:)
    enddo
    do i1 = 1, fw_nPntLd
        print '(1x,a23,2i12)', 'DOF: ', fw_load_DOF(i1,:)
    enddo
    print *

    close(fid)

end subroutine read_src_sim




subroutine open_surface_input_file()

    use global_variable

    implicit none

    integer*4 :: rowLen
    real*8, dimension(nDOF_disp_surf) :: u

    inquire( iolength=rowLen ) u

    if( TR_surf_trac_resp==DISP ) then
        open(newunit=fid_in_surf,file=trim(adjustl(dir_sim_in))//'disp_surf.dat' ,status='old',action='read',form='unformatted',access='direct',recl=rowLen)
    elseif( TR_surf_trac_resp==VEL ) then
        open(newunit=fid_in_surf,file=trim(adjustl(dir_sim_in))//'vel_surf.dat'  ,status='old',action='read',form='unformatted',access='direct',recl=rowLen)
    elseif( TR_surf_trac_resp==ACCEL ) then
        open(newunit=fid_in_surf,file=trim(adjustl(dir_sim_in))//'accel_surf.dat',status='old',action='read',form='unformatted',access='direct',recl=rowLen)
    endif

end subroutine open_surface_input_file




subroutine close_surface_input_file()

    use global_variable
    
    implicit none

    close( fid_in_surf )

end subroutine close_surface_input_file