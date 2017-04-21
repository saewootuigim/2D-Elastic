subroutine find_loading_node_fw()

    use global_variable
    
    implicit none

    integer*4 :: i1, i2
    real*8    :: compare_to_this
    real*8, allocatable, dimension(:) :: dist, dist_old

    print *
    print '(1x,a)', 'finding loading point'

    allocate( dist_old( fw_nPntLd ) )
    allocate(     dist( fw_nPntLd ) )

    dist_old = 2d0
    dist     = 2d0

    compare_to_this = sqrt(1.5d0**2+1.5d0**2)

    ! Find the nearest node from the 
    do i1 = 1, nNode
        do i2 = 1, fw_nPntLd

            if( abs(node2xy(i1,1)-fw_load_coord(i2,1))<compare_to_this .and. abs(node2xy(i1,2)-fw_load_coord(i2,2))<compare_to_this ) then
                dist(i2) = sqrt( (node2xy(i1,1)-fw_load_coord(i2,1))**2+(node2xy(i1,2)-fw_load_coord(i2,2))**2 )
                if( dist(i2)<dist_old(i2) ) then
                    fw_load_node(i2) = i1
                    dist_old(i2) = dist(i2)
                endif
            endif

        enddo
    enddo

    ! Update coordinate. Create fw_load_DOF.
    do i1 = 1, fw_nPntLd
        fw_load_DOF  (i1,:) = node2DOF(fw_load_node(i1),1:2)
        fw_load_coord(i1,:) = node2xy (fw_load_node(i1),1:2)
    enddo

    deallocate( dist, dist_old )

    ! Print the loading situation.
    do i1 = 1, fw_nPntLd
        write(*,101) i1, fw_load_node(i1), fw_load_coord(i1,:), fw_load_DOF(i1,:)
    enddo
    print *

101 format(1X,'loading node ',I2,', node #: ',I6,', coordinate: (',F9.4,',',F9.4,'), DOF: (', I6,',',I6,')')

end subroutine find_loading_node_fw




subroutine find_loading_node_TR()

    use global_variable
    
    implicit none

    integer*4 :: fid
    integer*4 :: i1, i2
    integer*4 :: DOF_to_turn_off
    real*8 :: gap
    real*8, allocatable, dimension(:) :: surf_trac_x_coord

    if( TR_surf_trac_case==0 ) then
        do i1 = 1, nDOF_disp_surf
            if( abs(node2xy(DOF2node(DOF_surf(i1)),1))>=150d0 ) then
                DOF_surf(i1) = 0
            elseif( abs(node2xy(DOF2node(DOF_surf(i1)),1))<150d0 ) then
                if( TR_surf_trac_dir==X_DIR ) then
                    DOF_to_turn_off = node2DOF(DOF2node(DOF_surf(i1)),2)
                    if( DOF_to_turn_off==DOF_surf(i1) ) DOF_surf(i1) = 0
                elseif( TR_surf_trac_dir==Y_DIR ) then
                    DOF_to_turn_off = node2DOF(DOF2node(DOF_surf(i1)),1)
                    if( DOF_to_turn_off==DOF_surf(i1) ) DOF_surf(i1) = 0
                endif
            endif
        enddo

    elseif( TR_surf_trac_case==NARROW ) then
        do i1 = 1, nDOF_disp_surf
            if( abs(node2xy(DOF2node(DOF_surf(i1)),1))>=40d0 ) then
                DOF_surf(i1) = 0
            elseif( abs(node2xy(DOF2node(DOF_surf(i1)),1))<40d0 ) then
                if( TR_surf_trac_dir==X_DIR ) then
                    DOF_to_turn_off = node2DOF(DOF2node(DOF_surf(i1)),2)
                    if( DOF_to_turn_off==DOF_surf(i1) ) DOF_surf(i1) = 0
                elseif( TR_surf_trac_dir==Y_DIR ) then
                    DOF_to_turn_off = node2DOF(DOF2node(DOF_surf(i1)),1)
                    if( DOF_to_turn_off==DOF_surf(i1) ) DOF_surf(i1) = 0
                endif
            endif
        enddo

    elseif( TR_surf_trac_case==ALL ) then
        if( TR_DtN==YES ) then
            do i1 = 1, nDOF_disp_surf
                if( TR_surf_trac_dir==X_DIR ) then
                    DOF_to_turn_off = node2DOF(DOF2node(DOF_surf(i1)),2)
                    if( DOF_to_turn_off==DOF_surf(i1) ) DOF_surf(i1) = 0
                elseif( TR_surf_trac_dir==Y_DIR ) then
                    DOF_to_turn_off = node2DOF(DOF2node(DOF_surf(i1)),1)
                    if( DOF_to_turn_off==DOF_surf(i1) ) DOF_surf(i1) = 0
                endif
            enddo
        elseif( TR_DtN==NO ) then
            if( TR_surf_trac_dir==X_DIR ) then
                where( DOF_surf>nDOF_disp_surf/2 ) DOF_surf = 0
            elseif( TR_surf_trac_dir==Y_DIR ) then
                where( DOF_surf<nDOF_disp_surf/2+1 ) DOF_surf = 0
            endif
        endif

    else
        ! x-coordinates of loading nodes
        gap = 2*rdwd/(TR_surf_trac_case+1d0)
        allocate( surf_trac_x_coord(TR_surf_trac_case) )
        do i1 = 1, TR_surf_trac_case
            surf_trac_x_coord(i1) = -rdwd + i1*gap
        enddo

        ! Search which surface DOF does it correspond.
        do i2 = 1, nDOF_disp_surf
            if( TR_surf_trac_dir==X_DIR ) then
                DOF_to_turn_off = node2DOF(DOF2node(DOF_surf(i2)),2)
                if( DOF_to_turn_off==DOF_surf(i2) ) then
                    DOF_surf(i2) = 0
                else
                    do i1 = 1, TR_surf_trac_case
                        if( abs(node2xy(DOF2node(DOF_surf(i2)),1)-surf_trac_x_coord(i1))<.45d0 ) then
                            exit
                        elseif( i1==TR_surf_trac_case ) then
                            DOF_surf(i2) = 0
                        endif
                    enddo
                endif
            elseif( TR_surf_trac_dir==Y_DIR ) then
                DOF_to_turn_off = node2DOF(DOF2node(DOF_surf(i2)),1)
                if( DOF_to_turn_off==DOF_surf(i2) ) then
                    DOF_surf(i2) = 0
                else
                    do i1 = 1, TR_surf_trac_case
                        if( abs(node2xy(DOF2node(DOF_surf(i2)),1)-surf_trac_x_coord(i1))<.45d0 ) then
                            exit
                        elseif( i1==TR_surf_trac_case ) then
                            DOF_surf(i2) = 0
                        endif
                    enddo
                endif
            endif
        enddo

        print '(1x,a)', 'surface traction applied on follwing nodes'
        print '(1x,4a12)', 'DOF', 'node', 'x-coord', 'y-coord'
        do i1 = 1, nDOF_disp_surf
            if( DOF_surf(i1)/=0 ) then
                print '(1x,2i12,2f12.4)', DOF_surf(i1), DOF2node(DOF_surf(i1)), node2xy(DOF2node(DOF_surf(i1)),1:2)
            endif
        enddo

    endif

!   open(newunit=fid,file='DOF_surf.txt',status='replace',action='write')
!   write(fid,'(a,i1)') 'time reverse case: ', TR_surf_trac_case
!   write(fid,'(a,i1)') 'loading direction: ', TR_surf_trac_dir
!   do i1 = 1, nDOF_disp_surf
!       if( DOF_surf(i1)/=0 ) then
!           write(fid,'(i6,2f10.4)') DOF_surf(i1), node2xy(DOF2node(DOF_surf(i1)),1:2)
!       endif
!   enddo
!   close(fid)

end subroutine find_loading_node_TR
