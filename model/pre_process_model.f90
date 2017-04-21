! This file is to read Pranav's input file.
!==============================================================================
!    + Reads FEM data file created by ANSYS (BE AWARE OF THE SPECIFIC FORMAT!)
!    + Updates associated module variables
!    + Generates degrees of freedom matrix
!    + Creates mapping between unpartitioned and partitioned unknowns
!==============================================================================

program pre_process_model

    integer*4 :: model_num
    integer*4 :: nNode, nElem
    integer*4 :: nPhyProp
    real*8    :: pml_l, rdwd, rdht

    integer*4 :: nDOF
    integer*4 :: nDOF_disp
    integer*4 :: nDOF_disp_srf
    integer*4 :: nDOF_disp_inc
    integer*4 :: nDOF_disp_RD
    integer*4 :: nDOF_disp_PML
    integer*4 :: nDOF_stress

    integer*4, allocatable, dimension(:)   :: DOF2node
    integer*4, allocatable, dimension(:,:) :: node2DOF
    integer*4, allocatable, dimension(:)   :: DOF_surf
    integer*4, allocatable, dimension(:)   :: DOF_disp

    integer*4, allocatable, dimension(:)   :: node2region
    real*8,    allocatable, dimension(:,:) :: node2xy

    integer*4, allocatable, dimension(:)   :: elem2region
    integer*4, allocatable, dimension(:)   :: elem2phyprp
    integer*4, allocatable, dimension(:,:) :: elem2node

    real*8,    allocatable, dimension(:)   :: E
    real*8,    allocatable, dimension(:)   :: rho
    real*8,    allocatable, dimension(:)   :: nu

    character(40) :: dirname

    print *
    write(*,'(1x,a14)',advance='no') 'model number: '
    read(*,*) model_num

    call read_parameter_PMK()
    call read_model( model_num )

    print *, 'Model info is stored in ''model_log.txt'''
    print *

    deallocate( DOF2node, node2DOF, DOF_surf, node2region, node2xy, elem2region, elem2phyprp, elem2node, E, rho, nu )

contains

    subroutine read_parameter_PMK()

        implicit none

        integer*4 :: i
        character(len=2) dirname
        logical :: L_exist

        write(dirname,'(i1,a1)') model_num, '/'
        inquire( file=dirname//'parameters_PMK.txt', exist=L_exist )
        if( L_exist==.false. ) then
            print '(1x,a)', 'ERROR **************************************************************************'
            print '(1x,a)', '      parameters_PMK.txt does no exist.'
            print '(1x,a)', '      This file contains the material properties.'
            stop
        endif

        open(unit=1,file=dirname//'parameters_PMK.txt',status='old',action='read')
        do i=1,49; read(1,*); enddo

        read(1,*) nPhyProp

        allocate( E  ( nPhyProp ) )
        allocate( rho( nPhyProp ) )
        allocate( nu ( nPhyProp ) )

        read(1,*)
        read(1,*)
        do i = 1, nPhyProp
            read(1,*) E(i)
        enddo

        read(1,*)
        read(1,*)
        do i = 1, nPhyProp
            read(1,*) rho(i)
        enddo

        read(1,*)
        read(1,*)
        do i = 1, nPhyProp
            read(1,*) nu(i)
        enddo

        close(1)

    end subroutine read_parameter_PMK

    subroutine read_model( model_num )

        implicit none

        integer*4, intent(in) :: model_num

        integer*4 :: i, j, k, i3, i4, fileLen, stat
        integer*4 :: iElem, iNode
        integer*4 :: locID, locPR
        integer*4 :: n_fixed_node, n_pml_node, n_reg_node, n_srf_node
        integer*4 :: elems(8)
        real*8 :: ndcoors(2)

        character(5) :: dummy

        select case( model_num )
            case(1)
                open(unit=1,file='1/mesh_50_100_4.txt',status='old',action='read',iostat=stat)
            case(2)
                open(unit=1,file='2/mesh_50_100_4_bottom_fixed.txt',status='old',action='read')
        end select
        if( stat/=0 ) then
            print *, 'file is not found.'
            stop
        endif

        ! read FEM geometry
        read(1,'(1X,D20.8)') rdwd          ! regular domain half width
        read(1,'(1X,D20.8)') rdht          ! regular domain height
        read(1,'(1X,D20.8)') pml_l

        ! read nodes and elements statistics
        read(1, '(2(1X,I12))') nNode, nElem

        !......................................
        allocate( node2xy(nNode,2) )        ! node ids and location flags (unsorted)
        allocate(node2DOF(nNode,5) )
        allocate(node2region(nNode))        ! node map (sorted)
        allocate(elem2region(nElem))        ! element map (sorted)
        allocate(elem2phyprp(nElem))
        allocate(elem2node(nElem,8))        ! element data (unsorted)
        !......................................

        ! read node coordinate
        n_fixed_node = 0
        n_pml_node   = 0
        n_reg_node   = 0
        n_srf_node   = 0
        read(1,'(1x,a5)') dummy

        do i = 1, nNode
            read(1,'(1x,i12,1x,2e20.8)') iNode, ndcoors
            node2xy(iNode,:) = ndcoors

            ! lateral fixed boundary
            if( abs(abs(ndcoors(1))-(rdwd+pml_l))<.001d0 ) then
                node2DOF(iNode,1:2) = 0                 ! Turn off displacement.
                node2DOF(iNode,3:5) = 1                 ! Turn on stress.
                node2region(iNode) = 4                  ! 4 = PML fixed node
                n_fixed_node = n_fixed_node + 1
            ! bottom boundary of PML domain
            elseif( abs(abs(ndcoors(2))-rdht)<.001d0 .and. abs(ndcoors(1))>rdwd-.001d0 ) then
                node2DOF(iNode,1:2) = 0                 ! Turn off displacement.
                node2DOF(iNode,3:5) = 1                 ! Turn on stress.
                node2region(iNode) = 4                  ! 4 = PML fixed node
                n_fixed_node = n_fixed_node + 1
            ! bottom boundary of regular domain
            elseif( abs(abs(ndcoors(2))-rdht)<.001d0 ) then
                node2DOF(iNode,1:5) = 0                 ! Turn off displacement and stress.
                node2region(iNode) = 3                  ! 3 = RD fixed node
                n_fixed_node = n_fixed_node + 1
            ! PML
            elseif( abs(ndcoors(1))>rdwd-.001d0 ) then
                node2DOF(iNode,1:5) = 1                 ! Turn on displacement and stress.
                node2region(iNode) = 2                  ! 2 = PML
                n_pml_node = n_pml_node + 1
            ! regular domain
            else
                node2DOF(iNode,1:2) = 1                 ! Turn on displacement.
                node2DOF(iNode,3:5) = 0                 ! Turn off stress.
                if( abs(ndcoors(2))>.001d0 ) then
                    node2region(iNode) = 1              ! 1 = RD
                    n_reg_node = n_reg_node + 1
                else
                    node2region(iNode) = 0              ! 0 = RD surface
                    n_srf_node = n_srf_node + 1
                endif
            endif
        enddo

        print '(1x,a17,i7)', 'fixed node=', n_fixed_node
        print '(1x,a17,i7)', 'PML node=', n_pml_node
        print '(1x,a17,i7)', 'RD node=', n_reg_node
        print '(1x,a17,i7)', 'surface node=', n_srf_node
        print *
        print '(1x,a17,i7)', 'sum of all=', n_fixed_node+n_pml_node+n_reg_node+n_srf_node
        print '(1x,a17,i7)', 'nNode=', nNode
        print *
        if( n_fixed_node+n_pml_node+n_reg_node+n_srf_node/=nNode ) then
            print *, 'ERROR: number of nodes doesn''t match.'
            stop
        endif

        ! read element block
        READ (1, '(1X,A5)') dummy
        DO i = 1, nElem
            READ (1, '(11(1X,I12))') iElem, elems, locID, locPR
            elem2node(iElem,:) = elems
            elem2phyprp(iElem) = locPR
            elem2region(iElem) = locID
        END DO

        CLOSE(1)

        ! "node2DOF"
        allocate(DOF2node(sum(node2DOF)))

        nDOF = 0
        DOF2node = 0
        i3 = 0

        do i = 0, 4
            do j = 1, 5
                do k = 1, nNode

                    ! surface nodes -> regular domain nodes -> PML domain nodes -> fixed boundary nodes
                    if( node2region(k)==i ) then
                        if( node2DOF(k,j)==1 ) then
                            nDOF = nDOF + 1
                            node2DOF(k,j) = nDOF
                            DOF2node(nDOF) = k
                            if( node2region(k)==3 .or. (node2region(k)==4 .and. (j==1 .or. j==2)) ) then
                                print *, 'DOF of fixed boundary has turn on.'
                                print *, 'node=', k, ', node region=', i, ', DOF type=', j
                                print *
                                print *, '*** legend ******************************************************************'
                                print *, '* node region * 0=surface, 1=RD, 2=PML, 3=RD fixed nodes, 4=PML fixed nodes *'
                                print *, '*    DOF type * 1=ux,      2=uy, 3=Sxx, 4=Syy,            5=Sxy             *'
                                print *, '*****************************************************************************'
                                stop
                            endif
                        endif
                    endif

                enddo

                if( j==2 ) then
                    if( i==0 ) then                 ! surface DOFs
                        nDOF_disp_srf = nDOF-i3
                    elseif( i==1 ) then             ! regular domain DOFs
                        nDOF_disp_RD  = nDOF-i3
                    elseif( i==2 ) then             ! PML domain displacement DOFs
                        nDOF_disp_PML = nDOF-i3
                    endif
                endif

            enddo

            i3 = nDOF

        enddo

        nDOF_disp = nDOF_disp_srf+nDOF_disp_RD+nDOF_disp_PML

        print '(1x,a17,i7)', 'nDOF_disp_srf=', nDOF_disp_srf
        print '(1x,a17,i7)', 'nDOF_disp_RD =', nDOF_disp_RD
        print '(1x,a17,i7)', 'nDOF_disp_PML=', nDOF_disp_PML
        print '(1x,a17,i7)', 'nDOF_disp    =', nDOF_disp
        print *

        !**********************************************************************************************************************
        ! Check if elements have correct domain numbering.
        !**********************************************************************************************************************
        if( model_num==2 ) then
            i3 = 0; i4 = 0
            do i = 1, nElem
                if( abs(sum(node2xy(elem2node(i,:),1))/8d0)<=rdwd ) then
                    i3 = i3 + 1
                    elem2region(i) = 1
                endif
                if( abs(sum(node2xy(elem2node(i,:),1))/8d0)>=rdwd ) then
                    i4 = i4 + 1
                    elem2region(i) = 2
                endif
            enddo
        endif
        !**********************************************************************************************************************

        !**********************************************************************************************************************
        ! Check if node region numbering is correct.
        !**********************************************************************************************************************
        if( model_num==2 ) then
            i3 = 0; i4 = 0
            do i = 1, nNode
                if( node2region(i)==3 .or. node2region(i)==4 ) then
                    if( abs(abs(node2xy(i,2))-rdht)>.001d0 .and. abs(abs(node2xy(i,1))-rdwd-pml_l)>.001d0 ) then
                        print '(1x,a5,i6,a,2f23.16,a1)', 'node ', i, ' is marked as fixed node.', node2xy(i,:)
                        i4 = i4 + 1
                        print *, 'accumulated = ', i4
                    endif
                endif
            enddo
        endif
        !**********************************************************************************************************************

        if( size(DOF2node)/=nDOF ) then
            print '(1x,a)', 'size of DOF2node and nDOF are different.'
            stop
        endif

        ! Collect the surface DOFs.
        allocate( DOF_surf( nDOF_disp_srf ) )
        DOF_surf = 0

        i3 = 0
        do i = 1, 2
            do j = 1, nNode
                if( node2region(j)==0 ) then
                    i3 = i3 + 1
                    DOF_surf(i3) = node2DOF(j,i)
                endif
            enddo
        enddo

        if( nDOF_disp_srf/=i3 ) then
            print '(1x,a)', 'number of surface DOFs is different from what has counted'
            stop
        endif

        ! Collect displacement DOFs.
        allocate( DOF_disp( nDOF_disp ) )
        DOF_disp = 0

        i3 = 0
        do i = 1, 2
            do j = 1, nNode
                if( node2region(j)/=3 .and. node2region(j)/=4 ) then
                    i3 = i3 + 1
                    DOF_disp(i3) = node2DOF(j,i)
                endif
            enddo
        enddo

        if( nDOF_disp==maxval(DOF_disp) ) then
            print *, '1:nDOF_disp is good to use. DOF_disp is not needed.'
            print *
        endif

        ! Print all results.
        write(dirname,'(i1,a1)') model_num, '/'

        inquire(iolength=fileLen) elem2node
        open(unit=1,file=trim(dirname)//'elem2node.dat',status='replace',action='write',form='unformatted',access='direct',recl=fileLen)
        open(unit=2,file=trim(dirname)//'elem2node.txt',status='replace',action='write')
        write(1,rec=1) elem2node
        do i = 1, nElem
            write(2,'(8i10)') elem2node(i,:)
        enddo
        close(1)
        close(2)

        inquire(iolength=fileLen) node2DOF
        open(unit=1,file=trim(dirname)//'node2DOF.dat',status='replace',action='write',form='unformatted',access='direct',recl=fileLen)
        open(unit=2,file=trim(dirname)//'node2DOF.txt',status='replace',action='write')
        write(1,rec=1) node2DOF
        do i = 1, nNode
            write(2,'(5i10)') node2DOF(i,:)
        enddo
        close(1)
        close(2)

        inquire(iolength=fileLen) node2xy
        open(unit=1,file=trim(dirname)//'node2xy.dat',status='replace',action='write',form='unformatted',access='direct',recl=fileLen)
        open(unit=2,file=trim(dirname)//'node2xy.txt',status='replace',action='write')
        write(1,rec=1) node2xy
        do i = 1, nNode
            write(2,'(2f14.8)') node2xy(i,:)
        enddo
        close(1)
        close(2)

        inquire(iolength=fileLen) DOF2node
        open(unit=1,file=trim(dirname)//'DOF2node.dat',status='replace',action='write',form='unformatted',access='direct',recl=fileLen)
        open(unit=2,file=trim(dirname)//'DOF2node.txt',status='replace',action='write')
        write(1,rec=1) DOF2node
        do i = 1, nDOF
            write(2,'(i10)') DOF2node(i)
        enddo
        close(1)
        close(2)

        inquire(iolength=fileLen) node2region
        open(unit=1,file=trim(dirname)//'node2region.dat',status='replace',action='write',form='unformatted',access='direct',recl=fileLen)
        open(unit=2,file=trim(dirname)//'node2region.txt',status='replace',action='write')
        write(1,rec=1) node2region
        do i = 1, nNode
            write(2,'(i10)') node2region(i)
        enddo
        close(1)
        close(2)

        inquire(iolength=fileLen) DOF_surf
        open(unit=1,file=trim(dirname)//'DOF_surf.dat',status='replace',action='write',form='unformatted',access='direct',recl=fileLen)
        open(unit=2,file=trim(dirname)//'DOF_surf.txt',status='replace',action='write')
        write(1,rec=1) DOF_surf
        write(2,'(i10)') DOF_surf
        close(1)
        close(2)

        inquire(iolength=fileLen) DOF_disp
        open(unit=1,file=trim(dirname)//'DOF_disp.dat',status='replace',action='write',form='unformatted',access='direct',recl=fileLen)
        open(unit=2,file=trim(dirname)//'DOF_disp.txt',status='replace',action='write')
        write(1,rec=1) DOF_disp
        write(2,'(i10)') DOF_disp
        close(1)
        close(2)

        open(unit=2,file=trim(dirname)//'elem2region.txt',status='replace',action='write')
        write(2,'(<nElem>i2)') elem2region
        close(2)

        open(unit=2,file=trim(dirname)//'elem2phyprp.txt',status='replace',action='write')
        write(2,'(<nElem>i2)') elem2phyprp
        close(2)

        open(unit=1,file=trim(dirname)//'model_log.txt',status='replace',action='write')

        write(1,101) '********************************************************************************'
        write(1,101) '                                                                       DIMENSION'
        write(1,101) '********************************************************************************'
        write(1,102) '                               PML length: ', pml_l
        write(1,102) '                regular domain half width: ', rdwd
        write(1,102) '                    regular domain height: ', rdht
        write(1,*)
        write(1,101) '********************************************************************************'
        write(1,101) '                                                                       MESH INFO'
        write(1,101) '********************************************************************************'
        write(1,103) '                          number of nodes: ', nNode
        write(1,103) '                       number of elements: ', nElem
        write(1,*)
        write(1,101) '********************************************************************************'
        write(1,101) '                                                          DEGREE OF FREEDOM INFO'
        write(1,101) '********************************************************************************'
        write(1,103) '                            number of DOF: ', nDOF
        write(1,103) '                    number of surface DOF: ', nDOF_disp_srf
        write(1,103) 'number of regular domain displacement DOF: ', nDOF_disp_RD
        write(1,103) '    number of PML domain displacement DOF: ', nDOF_disp_PML
        write(1,103) '     number of inclusion displacement DOF: ', nDOF_disp_inc
        write(1,103) '                 number of PML stress DOF: ', nDOF_stress
        write(1,*)
        write(1,101) '********************************************************************************'
        write(1,101) '                                                             MATERIAL PROPERTIES'
        write(1,101) '********************************************************************************'
        write(1,103) '            number of material properties: ', nPhyProp
        write(1,104) '                   elastic modulus(N/m^2): ', E
        write(1,104) '                     mass density(kg/m^3): ', rho
        write(1,104) '                  Poisson''s ratio(N/m^2): ', nu
        close(1)

    101 format(a80)
    102 format(a43,f11.4)
    103 format(a43,i11)
    104 format(a43,<nPhyProp>(es11.4/43x))

        open(unit=1,file=trim(dirname)//'constants.txt',status='replace',action='write')
        write(1,'(i12)') model_num
        write(1,'(i12)') nNode
        write(1,'(i12)') nElem
        write(1,'(f12.4)') pml_l
        write(1,'(f12.4)') rdwd
        write(1,'(f12.4)') rdht
        write(1,'(i12)') nDOF
        write(1,'(i12)') nDOF_disp
        write(1,'(i12)') nDOF_disp_srf
        write(1,'(i12)') nDOF_disp_inc
        write(1,'(i12)') nDOF_disp_RD
        write(1,'(i12)') nDOF_disp_PML
        write(1,'(i12)') nDOF_stress
        close(1)

    end subroutine read_model

end program pre_process_model