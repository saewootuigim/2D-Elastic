subroutine assemble_regular()

    use global_variable
    
    implicit none

    integer*4                   :: elem, i1, i2
    integer*4, dimension(16)    :: address, address_row, address_col
    real*8,    dimension(16,16) :: m_elem, c_elem, k_elem

    call system_clock( count3 )
    print *
    print '(1x,a)', 'assembling regular domain M, C and K'

  ! Start assembly.
    do elem = 1, nElem

      ! If the element is not it the regular domain, cycle.
        if( elem2region(elem)==2 .or. elem2region(elem)==3 .or. elem2region(elem)==4 .or. &
            elem2region(elem)==5 .or. elem2region(elem)==6 ) cycle

      ! Generate element matrices for the mass and stiffness.
        call element_m( elem, m_elem, 'rd' )
        call element_k( elem, k_elem )

      ! Obtain the address which maps each element matrix entry to the global matrix.
      ! 'address' is subtracted by 1 in order to comply the 'C' indexing, which begins with 0, not 1.
        do i1 = 1, 2
            do i2 = 1, 8
                address(i2+(i1-1)*8) = node2DOF(elem2node(elem,i2),i1)-1
            end do
        end do
      
      ! Assemble
        c_elem = 0d0
        
        call MatSetValues( M, 16, address, 16, address, transpose(m_elem), ADD_VALUES, ierr )
        call MatSetValues( C, 16, address, 16, address, transpose(c_elem), ADD_VALUES, ierr )
        call MatSetValues( K, 16, address, 16, address, transpose(k_elem), ADD_VALUES, ierr )

      ! For DtN TR, form force vector.
        if( F_or_T==TR .and. TR_DtN==NO ) then

            if( maxval(node2xy(elem2node(elem,:),2))<-.002D0 ) cycle

            do i1 = 1, 2
                do i2 = 1, 8
                    address_row(i2+(i1-1)*8) =  node2DOF(elem2node(elem,i2),i1)-1
                    address_col(i2+(i1-1)*8) = -node2DOF(elem2node(elem,i2),i1)-1
                enddo
            enddo
            
            ! call MatSetValues( Ms, 16, address_row, 16, address_col, transpose(m_elem), ADD_VALUES, ierr )
            call MatSetValues( Ks, 16, address_row, 16, address_col, transpose(k_elem), ADD_VALUES, ierr )
        endif
    end do

    call system_clock( count4 )

    print '(1x,a,f8.4,a)', 'elapse: ', dble(count4-count3)/dble(count_rate), 'sec.'
    print *
    
end subroutine assemble_regular




subroutine assemble_PML()

    use global_variable
    
    implicit none
    
    integer*4                   :: i1, elem
    integer*4, dimension(16)    :: address_disp
    integer*4, dimension(24)    :: address_stress
    real*8,    dimension(16,16) :: ma_elem, mb_elem, mc_elem
    real*8,    dimension(24,24) :: na_elem, nb_elem, nc_elem
    real*8,    dimension(16,24) :: aeu_elem, apu_elem, ael_elem_trans, apl_elem_trans
    real*8,    dimension(24,16) :: ael_elem, apl_elem

    call system_clock( count3 )
    print *
    print '(1x,a)', 'assembling PML domain M, C and K'

  ! Start assembly.
    do elem = 1, nElem

      ! Cycle if it is not PML region.
        if( elem2region(elem)==1 .or. elem2region(elem)==7 ) cycle

      ! Generate element matrices for the mass and stiffness.
        call element_m( elem, ma_elem, 'a' )
        call element_m( elem, mb_elem, 'b' )
        call element_m( elem, mc_elem, 'c' )

        call element_n( elem, na_elem, 'a' )
        call element_n( elem, nb_elem, 'b' )
        call element_n( elem, nc_elem, 'c' )

        call element_a( elem, aeu_elem, 'eu' )
        call element_a( elem, apu_elem, 'pu' )
        call element_a( elem, ael_elem_trans, 'el' )
        call element_a( elem, apl_elem_trans, 'pl' )

      ! Transpose ael_elem_trans and store in ael_elem.
        ael_elem = transpose(ael_elem_trans)
        apl_elem = transpose(apl_elem_trans)

      ! Obtain the addresses which map each element matrix entry to the global matrix.
        do i1 = 1, 8
            address_disp  (i1   ) = node2DOF(elem2node(elem,i1),1)-1
            address_disp  (i1+8 ) = node2DOF(elem2node(elem,i1),2)-1
            address_stress(i1   ) = node2DOF(elem2node(elem,i1),3)-1
            address_stress(i1+8 ) = node2DOF(elem2node(elem,i1),4)-1
            address_stress(i1+16) = node2DOF(elem2node(elem,i1),5)-1
        end do

      ! Assemble M matrix or diag_M vector depend on what kind of solver is used.
        call MatSetValues( M, 16,   address_disp, 16,   address_disp, transpose(ma_elem),  ADD_VALUES, ierr )
        call MatSetValues( M, 24, address_stress, 24, address_stress, transpose(na_elem),  ADD_VALUES, ierr )

      ! Assemble C matrix.
        call MatSetValues( C, 16,   address_disp, 16,   address_disp, transpose(mb_elem),  ADD_VALUES, ierr )
        call MatSetValues( C, 24, address_stress, 24, address_stress, transpose(nb_elem),  ADD_VALUES, ierr )
        call MatSetValues( C, 16,   address_disp, 24, address_stress, transpose(aeu_elem), ADD_VALUES, ierr )
        call MatSetValues( C, 24, address_stress, 16,   address_disp, transpose(ael_elem), ADD_VALUES, ierr )
        
      ! Assemble K marix.
        call MatSetValues( K, 16,   address_disp, 16,   address_disp, transpose(mc_elem),  ADD_VALUES, ierr )
        call MatSetValues( K, 24, address_stress, 24, address_stress, transpose(nc_elem),  ADD_VALUES, ierr )
        call MatSetValues( K, 16,   address_disp, 24, address_stress, transpose(apu_elem), ADD_VALUES, ierr )
        call MatSetValues( K, 24, address_stress, 16,   address_disp, transpose(apl_elem), ADD_VALUES, ierr )

    end do

    call system_clock( count4 )

    print '(1x,a,f8.4,a)', 'elapse: ', dble(count4-count3)/dble(count_rate), 'sec.'
    print *
    
end subroutine assemble_PML



! This subroutine must be refined later.
subroutine assemble_load()

    use global_variable
    
    implicit none

	integer*4 :: i1
	real*8,    allocatable, dimension(:,:) :: y

	allocate( y( fw_nPntLd,2 ) )
	y = 0d0

	! loading direction
	if( fw_PntLd_dir==X_DIR ) then
		if( fw_load_intensity==EQ_FORCE ) then
			y(:,1) = 1d0
		elseif( fw_load_intensity==EQ_RESPONSE ) then
! model 6, loading at ( -60, -190) and ( 60,  260)
			y(1,1) = 1.046373717084540d0
			y(2,1) = 1d0
! model 6, loading at ( -40, -260) and ( 40, -260)
!           y(1,1) = 0.999762899804392d0
!           y(2,1) = 1d0
! model 6, loading at ( -20, -220) and ( 20, -240)
!         	y(1,1) = 1d0
!          	y(2,1) = 0.879201056936112d0
        endif
	elseif( fw_PntLd_dir==Y_DIR ) then
		y(:,2) = 1d0
	elseif( fw_PntLd_dir==XY_DIR ) then
		y      = 1d0
	endif

	! Assemble
	do i1 = 1, fw_nPntLd
		call VecSetValues( F, 2, fw_load_DOF(i1,:)-1, y(i1,:), INSERT_VALUES, ierr )
	enddo

	deallocate( y )

end subroutine assemble_load