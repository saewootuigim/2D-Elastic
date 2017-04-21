program main

    use global_variable
    
    implicit none

    character(8)  :: date
    character(10) :: time

integer*4 :: fid, i1

    ! Start PETSc.
    call PetscInitialize( petsc_null_character, ierr )

    ! starting time
    call system_clock( count1, count_rate, count_max )

    ! starting time stamp
    call date_and_time( date, time ) 
    print *
    print '(1x,a13,a4,a1,a2,a1,a2,1x,a2,a1,a2,a1,a2)', 'Staing time: ', &
        date(1:4), '/', date(5:6), '/', date(7:8), time(1:2), ':', time(3:4), ':', time(5:6)
        
    ! Read input file.
    call read_input_file()

    ! Read source simulation info.
    if( F_or_T==TR ) call read_src_sim()
    
    ! Create output directory.
    call create_output_directory()
    
    ! Gauss quadrature nodes and weights.
    if( quad_meth==G_L ) then
        call Gauss_Legendre( NQuad, xi, wi )
    elseif( quad_meth==G_K ) then
        call Gauss_Kronrod ( NQuad, xi, wi )
    elseif( quad_meth==G_Lbt ) then
        NQuad = 5
        if(allocated(xi)) deallocate(xi)
        allocate(xi(NQuad))
        xi = [-1D0, -.6546536707079771437983D0, 0D0, .6546536707079771437983D0, 1D0]
        if(allocated(wi)) deallocate(wi)
        allocate(wi(NQuad))
        wi = [.1D0, .5444444444444444444444D0, .7111111111111111111111D0, .5444444444444444444444D0, .1D0]
    endif
    
    ! Read model.
    call read_model()

    ! For DtN TR simulation, adjust DOF.
    if( F_or_T==TR .and. TR_DtN==NO ) call adjust_DOF_for_DtN()
    
    ! Create PETSc objects.
    call MatCreateAIJ( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nDOF, nDOF, 150, PETSC_NULL_INTEGER, 150, PETSC_NULL_INTEGER, M, ierr )
    call MatCreateAIJ( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nDOF, nDOF, 150, PETSC_NULL_INTEGER, 150, PETSC_NULL_INTEGER, C, ierr )
    call MatCreateAIJ( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nDOF, nDOF, 150, PETSC_NULL_INTEGER, 150, PETSC_NULL_INTEGER, K, ierr )
    call VecCreateMPI( PETSC_COMM_WORLD, PETSC_DECIDE, nDOF, F, ierr )
    call VecSetOption( F, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE, ierr )
    call VecSet( F, 0D0, ierr )

    if( F_or_T==TR .and. TR_DtN==NO ) then
        ! call MatCreateAIJ( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nDOF, nDOF_disp_surf, nDOF_disp_surf, PETSC_NULL_INTEGER, 20, PETSC_NULL_INTEGER, Ms, ierr )
        call MatCreateAIJ( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nDOF, nDOF_disp_surf, nDOF_disp_surf, PETSC_NULL_INTEGER, 20, PETSC_NULL_INTEGER, Ks, ierr )
    endif

    ! Assemble M, C, K in regular domain.
    call assemble_regular()

    ! Assemble M, C, K in PML domain.
    call assemble_PML()

    ! Assemble F.
    if( F_or_T==FORWARD ) then
        call find_loading_node_fw()
        call assemble_load()
    elseif( F_or_T==TR ) then
        call find_loading_node_TR()
        if( TR_apply_sink==YES ) call assemble_load()
    endif
    
    ! Assemble PETSc objects.
    call MatAssemblyBegin( M, MAT_FINAL_ASSEMBLY, ierr )
    call MatAssemblyEnd  ( M, MAT_FINAL_ASSEMBLY, ierr )
    call MatAssemblyBegin( C, MAT_FINAL_ASSEMBLY, ierr )
    call MatAssemblyEnd  ( C, MAT_FINAL_ASSEMBLY, ierr )
    call MatAssemblyBegin( K, MAT_FINAL_ASSEMBLY, ierr )
    call MatAssemblyEnd  ( K, MAT_FINAL_ASSEMBLY, ierr )
    call VecAssemblyBegin( F, ierr )
    call VecAssemblyEnd  ( F, ierr )
    if( F_or_T==TR .and. TR_DtN==NO ) then
        ! call MatAssemblyBegin( Ms, MAT_FINAL_ASSEMBLY, ierr )
        ! call MatAssemblyEnd  ( Ms, MAT_FINAL_ASSEMBLY, ierr )
        call MatAssemblyBegin( Ks, MAT_FINAL_ASSEMBLY, ierr )
        call MatAssemblyEnd  ( Ks, MAT_FINAL_ASSEMBLY, ierr )
    endif

    ! Solve.
    call Newmark()

    ! Normalize surface response.
    if( F_or_T==FORWARD ) call normalize_surface_output()
    
    ! Write simulation info.
    if( F_or_T==FORWARD ) then
        call write_fw_simulation_info()
    elseif( F_or_T==TR ) then
        call write_TR_simulation_info()
    endif

    ! Destroy PETSc objects.
    call MatDestroy( M, ierr )
    call MatDestroy( C, ierr )
    call MatDestroy( K, ierr )
    call VecDestroy( F, ierr )
    if( F_or_T==TR .and. TR_DtN==NO ) then
        ! call MatDestroy( Ms, ierr )
        call MatDestroy( Ks, ierr )
    endif

    ! starting time stamp
    call date_and_time( date, time )
    call system_clock( count2 )
    if( count2<count1 ) count2 = count2+count_max
    print '(1x,a13,a4,a1,a2,a1,a2,1x,a2,a1,a2,a1,a2)', 'Ending time: ', &
        date(1:4), '/', date(5:6), '/', date(7:8), time(1:2), ':', time(3:4), ':', time(5:6)
    print '(1x,a8,f12.4,a4)', 'Elapse: ', dble(count2-count1)/dble(count_rate), 'sec.'

    ! Terminate PETSc.
    call PetscFinalize( ierr )

contains

    subroutine adjust_DOF_for_DtN()

        implicit none

        integer*4 :: i1, i2
        integer*4, allocatable, dimension(:) :: arr_adjust

        print *, 'DOF numbers adjusting...'

        ! Adjust DOF2node.
        allocate(arr_adjust(nDOF-nDOF_disp_surf))

        do i1 = nDOF_disp_surf+1, nDOF
            arr_adjust(i1-nDOF_disp_surf) = DOF2node(i1)
        enddo

        deallocate(DOF2node)
        allocate(DOF2node(nDOF-nDOF_disp_surf))
        DOF2node = arr_adjust

        deallocate(arr_adjust)

        ! Adjust DOF_disp
        allocate(arr_adjust(nDOF_disp-nDOF_disp_surf))
        
        i2 = 0
        do i1 = 1, nDOF_disp
            if( DOF_disp(i1)>nDOF_disp_surf ) then
                i2 = i2+1
                arr_adjust(i2) = DOF_disp(i1)-nDOF_disp_surf
            endif
        enddo

        deallocate(DOF_disp)
        allocate(DOF_disp(nDOF_disp-nDOF_disp_surf))
        DOF_disp = arr_adjust

        deallocate(arr_adjust)

        ! Adjust node2DOF.
        do i1 = 1, 5
            do i2 = 1, nNode
                if( node2DOF(i2,i1)==0 ) cycle
                if( node2region(i2)==0 ) then
                    node2DOF(i2,i1) = -node2DOF(i2,i1)
                else
                    node2DOF(i2,i1) =  node2DOF(i2,i1)-nDOF_disp_surf
                endif
            enddo
        enddo

        ! Adjust number of DOF.
        nDOF         = nDOF         - nDOF_disp_surf
        nDOF_disp    = nDOF_disp    - nDOF_disp_surf
        nDOF_disp_RD = nDOF_disp_RD - nDOF_disp_surf
        fw_load_DOF  = fw_load_DOF  - nDOF_disp_surf

        write(*,'(10x,2a10)') 'before', 'after'
        write(*,'(1x,a9,2i10)') '     nDOF', nDOF+nDOF_disp_surf, nDOF
        write(*,'(1x,a9,2i10)') 'nDOF_disp', nDOF_disp+nDOF_disp_surf, nDOF_disp
        write(*,*)


    end subroutine adjust_DOF_for_DtN

end program main
