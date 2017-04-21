subroutine element_k( elem, k_elem )

    use global_variable
    
    implicit none

	integer*4,                 intent(in) :: elem
	real*8, dimension(16,16), intent(out) :: k_elem

	integer*4              :: i1, i2, i3, i4
	real*8                 :: det_Jacobi, lambda, mu
	real*8, dimension(8)   :: N
	real*8, dimension(2,2) :: Jacobi, inv_Jacobi
	real*8, dimension(2,8) :: coord_trans
	real*8, dimension(8,2) :: dN
	real*8, dimension(8,8) :: NxNx, NyNy, NxNy, NyNx, TopLft, TopRgt, BotLft, BotRgt

	! Initialize k_elem
	k_elem = 0d0

	! Compute elasticity constants.
	lambda = E(elem2phyprp(elem))*nu(elem2phyprp(elem))/((1d0+nu(elem2phyprp(elem)))*(1d0-2d0*nu(elem2phyprp(elem))))
	mu = E(elem2phyprp(elem))/(2d0*(1d0+nu(elem2phyprp(elem))))

	! Obtain transformed coordinate matrix.
	do i1 = 1, 8
		do i2 = 1, 2
			coord_trans(i2,i1) = node2xy(elem2node(elem,i1),i2)
		end do
	end do

	! Gauss integration
	do i1 = 1, NQuad
		do i2 = 1, NQuad

			! Obtain the shape function and the derivative of shape function's function values at sample point.
			! Coordinate system is in xi-eta.
			call shape8( xi(i2), xi(i1), N, dN )

			! Find Jacobian.
			Jacobi = matmul( coord_trans, dN )

			! The determinant of Jacobian
			det_Jacobi = Jacobi(1,1)*Jacobi(2,2) - Jacobi(1,2)*Jacobi(2,1)

			if(det_Jacobi<0) then
				print '(a,es11.4,a,i6)', 'Jacobian is negative(', det_Jacobi, ') for element ', elem
				print '(a,a)', 'in subroutine "elem_k".'
				print '(a)', 'The coordinate of this point is'
				do i3=1,2; print '(9es11.4)', (coord_trans(i3,i4),i4=1,8); enddo
				print '(a,i1,a,i1)', 'The loop indices are i1=', i1, ' i2=', i2
				stop
			end if

			! the inverse of Jacobian
			inv_Jacobi(1,1) =  Jacobi(2,2)/det_Jacobi
			inv_Jacobi(2,2) =  Jacobi(1,1)/det_Jacobi
			inv_Jacobi(1,2) = -Jacobi(1,2)/det_Jacobi
			inv_Jacobi(2,1) = -Jacobi(2,1)/det_Jacobi

			! transformed 'derivative of shape function values' in x-y coordinate.
			dN = matmul( dN, inv_Jacobi )

			! Form the each block of K matrix.
			do i3 = 1, 8
				do i4 = 1, 8
					NxNx(i3,i4) = dN(i3,1) * dN(i4,1) ! Nx*Nx'
					NyNy(i3,i4) = dN(i3,2) * dN(i4,2) ! Ny*Ny'
					NxNy(i3,i4) = dN(i3,1) * dN(i4,2) ! Nx*Ny'
					NyNx(i3,i4) = dN(i3,2) * dN(i4,1) ! Ny*Nx'
				end do
			end do

			TopLft = (lambda+2d0*mu)*NxNx +             mu *NyNy ! Left Top
			TopRgt =         lambda *NxNy +             mu *NyNx ! Right Top
			BotLft =             mu *NxNy +         lambda *NyNx ! Left Bottom
			BotRgt =             mu *NxNx + (lambda+2d0*mu)*NyNy ! Right Bottom

			! (lambda+2mu)(Nx*Nx')+mu(Ny*Ny')
			k_elem( 1: 8, 1: 8) = k_elem( 1: 8, 1: 8) + TopLft*wi(i1)*wi(i2)*det_Jacobi
			! (lambda    )(Nx*Ny')+mu(Ny*Nx')
			k_elem( 1: 8, 9:16) = k_elem( 1: 8, 9:16) + TopRgt*wi(i1)*wi(i2)*det_Jacobi
			! (lambda    )(Ny*Nx')+mu(Nx*Ny')
			k_elem( 9:16, 1: 8) = k_elem( 9:16, 1: 8) + BotLft*wi(i1)*wi(i2)*det_Jacobi
			! (lambda+2mu)(Nx*Nx')+mu(Nx*Nx')
			k_elem( 9:16, 9:16) = k_elem( 9:16, 9:16) + BotRgt*wi(i1)*wi(i2)*det_Jacobi
		end do
	end do

end subroutine element_k




subroutine element_m( elem, m_elem, flag )

    use global_variable
    
    implicit none

	integer*4,                 intent(in) :: elem
	character(*),              intent(in) :: flag
	real*8, dimension(16,16), intent(out) :: m_elem

	integer*4              :: i1, i2, i3, i4
	real*8                 :: det_Jacobi, fac1
	real*8                 :: ax, ay, bx, by, axdy, aydx, bxdy, bydx, aa, bb, cc
	real*8, dimension(2)   :: location
	real*8, dimension(8)   :: N
	real*8, dimension(2,2) :: Jacobi, inv_Jacobi
	real*8, dimension(2,8) :: coord_trans
	real*8, dimension(8,2) :: dN
	real*8, dimension(8,8) :: NN

	! Initialize m_elem
	m_elem = 0d0

	! Obtain transformed coordinate matrix.
	do i1 = 1, 8
		do i2 = 1, 2
			coord_trans(i2,i1) = node2xy(elem2node(elem,i1),i2)
		end do
	end do

	! Gauss integration
	do i1 = 1, NQuad
		do i2 = 1, NQuad

			! Obtain shape function and derivative of shape function.
			call shape8( xi(i2), xi(i1), N, dN )

			! Jacobian
			Jacobi = matmul( coord_trans, dN )

			! The determinant of Jacobian
			det_Jacobi = Jacobi(1,1)*Jacobi(2,2) - Jacobi(1,2)*Jacobi(2,1)

			if( det_Jacobi<0 ) then
				print '(1x,a)', '-- ERROR --'
				print '(4x,a,i6,a,a,a)', 'Jacobian is negative in element ', elem, ' in subroutine "elem_m". "flag" was ', flag, '.'
				print '(4x,a)', 'The coordinate of this point is'
				do i3 = 1, 2
					print '(9f8.2)', coord_trans(i3,:)
				end do
				print '(4x,a34,f10.6,a1,f10.6,a1)', 'The quadrature point is (xi,eta)=(', xi(i1), ',', xi(i2), ')'
				print '(4x,a)', 'shape functions and derivatives of shape functions are: N, N_xi, N_eta'
				do i3 = 1, 8
					print '(3f10.6)', N(i3), dN(i3,:)
				end do
				print '(4x,a)', 'Jacobi matrix of this quadrature point is'
				do i3 = 1, 2
					print '(2f8.4)', Jacobi(i3,:)
				end do
				print '(4x,a,e16.6e3,a)', 'Jacobian of this quadrature point is', det_Jacobi, '.'
				stop
			end if

			! the inverse of Jacobian
			inv_Jacobi(1,1) =  Jacobi(2,2)/det_Jacobi
			inv_Jacobi(2,2) =  Jacobi(1,1)/det_Jacobi
			inv_Jacobi(1,2) = -Jacobi(1,2)/det_Jacobi
			inv_Jacobi(2,1) = -Jacobi(2,1)/det_Jacobi

			! transformed 'derivative of shape function values' in x-y coordinate.
			dN = matmul( dN, inv_Jacobi )

			! Obtain PML constants if needed.
			if( flag/='rd' ) then
				location = matmul( coord_trans, N )
				call PML_coefficient( location, a0, b0, rdwd, rdht, pml_l, ax, ay, bx, by, axdy, aydx, bxdy, bydx, aa, bb, cc )
			end if

			! calculate fac1.
			select case( flag )
				case('rd')
					fac1 = rho( elem2phyprp( elem ) )
				case('a') ! Ma
					fac1 = rho( elem2phyprp( elem ) ) * aa
				case('b') ! Mb
					fac1 = rho( elem2phyprp( elem ) ) * bb
				case('c') ! Mc
					fac1 = rho( elem2phyprp( elem ) ) * cc
			end select

			! Form the block of matrices, the temps.
			do i3 = 1, 8
				do i4 = 1, 8
					NN(i4,i3) = N(i4)*N(i3)
				end do
			end do

			! Assign to the element M matrix.
			m_elem( 1: 8, 1: 8) = m_elem( 1: 8, 1: 8) + NN*fac1*wi(i1)*wi(i2)*det_Jacobi
			m_elem( 9:16, 9:16) = m_elem( 9:16, 9:16) + NN*fac1*wi(i1)*wi(i2)*det_Jacobi
		end do
	end do

end subroutine element_m




subroutine element_n( elem, n_elem, flag )

    use global_variable
    
    implicit none

	integer*4,                 intent(in) :: elem
	character(*),              intent(in) :: flag
	real*8, dimension(24,24), intent(out) :: n_elem

	integer*4              :: i1, i2, i3, i4
	real*8                 :: det_Jacobi, fac1
	real*8                 :: ax, ay, bx, by, axdy, aydx, bxdy, bydx, aa, bb, cc
	real*8, dimension(2)   :: location
	real*8, dimension(8)   :: N
	real*8, dimension(2,2) :: Jacobi, inv_Jacobi
	real*8, dimension(2,8) :: coord_trans
	real*8, dimension(8,2) :: dN
	real*8, dimension(8,8) :: NN
	!------------------------------------------------------------------------------------------------------------------------------

	! Initialize n_elem
	n_elem = 0d0

	! Obtain transformed coordinate matrix.
	do i1 = 1, 8
		do i2 = 1, 2
			coord_trans(i2,i1) = node2xy(elem2node(elem,i1),i2)
		end do
	end do

	! Gauss integration
	do i1 = 1, NQuad
		do i2 = 1, NQuad

			! Obtain shape function and derivative of shape function.
			call shape8( xi(i2), xi(i1), N, dN )

			! Jacobian
			Jacobi = matmul( coord_trans, dN )

			! The determinant of Jacobian
			det_Jacobi = Jacobi(1,1)*Jacobi(2,2) - Jacobi(1,2)*Jacobi(2,1)

			if( det_Jacobi<0 ) then
				print '(a,es11.4,a,i6)', 'Jacobian is negative(', det_Jacobi, ') for element ', elem
				print '(a,a)', 'in subroutine "elem_n". flag was ', flag
				print '(a)', 'The coordinate of this point is'
				do i3=1,2; print '(9es11.4)', (coord_trans(i3,i4),i4=1,8); enddo
				print '(a,i1,a,i1)', 'The loop indices are i1=', i1, ' i2=', i2
				stop
			end if

			! the inverse of Jacobian
			inv_Jacobi(1,1) =  Jacobi(2,2)/det_Jacobi
			inv_Jacobi(2,2) =  Jacobi(1,1)/det_Jacobi
			inv_Jacobi(1,2) = -Jacobi(1,2)/det_Jacobi
			inv_Jacobi(2,1) = -Jacobi(2,1)/det_Jacobi

			! transformed 'derivative of shape function values' in x-y coordinate.
			dN = matmul( dN, inv_Jacobi )

			! Obtain PML constants.
			location = matmul( coord_trans, N )
			call PML_coefficient( location, a0, b0, rdwd, rdht, pml_l, ax, ay, bx, by, axdy, aydx, bxdy, bydx, aa, bb, cc )

			! Compute facs.
			select case(flag)
				case('a') ! Na
					fac1 = aa
				case('b') ! Nb
					fac1 = bb
				case('c') ! Nc
					fac1 = cc
			end select

			! Form the block.
			do i3 = 1, 8
				do i4 = 1, 8
					NN(i3,i4) = N(i3)*N(i4)
				end do
			end do

			! normal stress
			n_elem( 1: 8, 1: 8) = n_elem( 1: 8, 1: 8) + NN*fac1*wi(i1)*wi(i2)*det_Jacobi
			n_elem( 9:16, 9:16) = n_elem( 9:16, 9:16) + NN*fac1*wi(i1)*wi(i2)*det_Jacobi
			n_elem(17:24,17:24) = n_elem(17:24,17:24) + NN*fac1*wi(i1)*wi(i2)*det_Jacobi*2d0
		end do
	end do

end subroutine element_n




subroutine element_a( elem, a_elem, flag )

    use global_variable
    
    implicit none

	integer*4,                 intent(in) :: elem
	character(*),              intent(in) :: flag
	real*8, dimension(16,24), intent(out) :: a_elem

	integer*4              :: i1, i2, i3, i4
	real*8                 :: lambda, mu
	real*8                 :: det_Jacobi
	real*8                 :: ax, ay, bx, by, axdy, aydx, bxdy, bydx, aa, bb, cc
	real*8, dimension(2)   :: location
	real*8, dimension(8)   :: N
	real*8, dimension(2,2) :: Jacobi, inv_Jacobi
	real*8, dimension(2,8) :: coord_trans
	real*8, dimension(8,2) :: dN
	real*8, dimension(8,8) :: NxPsi, NyPsi, NPsi

	! Initialize a_elem
	a_elem = 0d0

	! Compute elasticity constants.
	lambda = E(elem2phyprp(elem))*nu(elem2phyprp(elem))/((1d0+nu(elem2phyprp(elem)))*(1d0-2d0*nu(elem2phyprp(elem))))
	mu = E(elem2phyprp(elem))/(2d0*(1d0+nu(elem2phyprp(elem))))

	! Obtain transformed coordinate matrix.
	do i1 = 1, 8
		do i2 = 1, 2
			coord_trans(i2,i1) = node2xy(elem2node(elem,i1),i2)
		end do
	end do

	! Gauss integration
	do i1 = 1, NQuad
		do i2 = 1, NQuad

			! Obtain shape function and derivative of shape function.
			call shape8( xi(i2), xi(i1), N, dN )

			! Find Jacobian.
			Jacobi = matmul( coord_trans, dN )

			! The determinant of Jacobian
			det_Jacobi = Jacobi(1,1)*Jacobi(2,2) - Jacobi(1,2)*Jacobi(2,1)

			if( det_Jacobi<0 ) then
				print '(a,es11.4,a,i6)', 'Jacobian is negative(', det_Jacobi, ') for element ', elem
				print '(a,a)', 'in subroutine "elem_a". flag was ', flag
				print '(a)', 'The coordinate of this point is'
				do i3=1,2; print '(9es11.4)', (coord_trans(i3,i4),i4=1,8); enddo
				print '(a,i1,a,i1)', 'The loop indices are i1=', i1, ' i2=', i2
				stop
			end if

			! the inverse of Jacobian
			inv_Jacobi(1,1) =  Jacobi(2,2) / det_Jacobi
			inv_Jacobi(2,2) =  Jacobi(1,1) / det_Jacobi
			inv_Jacobi(1,2) = -Jacobi(1,2) / det_Jacobi
			inv_Jacobi(2,1) = -Jacobi(2,1) / det_Jacobi

			! function values of 'derivative of shape function' in x-y coordinate.
			dN = matmul( dN, inv_Jacobi )

			! Obtain PML constants.
			location = matmul( coord_trans, N )
			call PML_coefficient( location, a0, b0, rdwd, rdht, pml_l, ax, ay, bx, by, axdy, aydx, bxdy, bydx, aa, bb, cc )

			! Form matric block.
			do i3 = 1, 8
				do i4 = 1, 8
					NxPsi(i3,i4) = dN(i3,1)*N(i4)
					NyPsi(i3,i4) = dN(i3,2)*N(i4)
					NPsi(i3,i4)  = N(i3)*N(i4)
				end do
			end do

			! Assign to a_elem.
			select case(flag)
				case('eu') !    A_{e}  + A_{Me}
					a_elem( 1: 8, 1: 8) = a_elem( 1: 8, 1: 8) +     (ay*NxPsi+aydx*NPsi)*wi(i1)*wi(i2)*det_Jacobi
					a_elem( 9:16, 9:16) = a_elem( 9:16, 9:16) +     (ax*NyPsi+axdy*NPsi)*wi(i1)*wi(i2)*det_Jacobi
					a_elem( 1: 8,17:24) = a_elem( 1: 8,17:24) +     (ax*NyPsi+axdy*NPsi)*wi(i1)*wi(i2)*det_Jacobi
					a_elem( 9:16,17:24) = a_elem( 9:16,17:24) +     (ay*NxPsi+aydx*NPsi)*wi(i1)*wi(i2)*det_Jacobi
				case('el') ! - (A_{mu e} + A_{lambda e})
					a_elem( 1: 8, 1: 8) = a_elem( 1: 8, 1: 8) - (2d0*mu+lambda)*ay*NxPsi*wi(i1)*wi(i2)*det_Jacobi
					a_elem( 1: 8, 9:16) = a_elem( 1: 8, 9:16) -         lambda *ay*NxPsi*wi(i1)*wi(i2)*det_Jacobi
					a_elem( 1: 8,17:24) = a_elem( 1: 8,17:24) -  2d0*mu        *ax*NyPsi*wi(i1)*wi(i2)*det_Jacobi
					a_elem( 9:16, 1: 8) = a_elem( 9:16, 1: 8) -         lambda *ax*NyPsi*wi(i1)*wi(i2)*det_Jacobi
					a_elem( 9:16, 9:16) = a_elem( 9:16, 9:16) - (2d0*mu+lambda)*ax*NyPsi*wi(i1)*wi(i2)*det_Jacobi
					a_elem( 9:16,17:24) = a_elem( 9:16,17:24) -  2d0*mu        *ay*NxPsi*wi(i1)*wi(i2)*det_Jacobi
				case('pu') !    A_{p}  + A_{Mp}
					a_elem( 1: 8, 1: 8) = a_elem( 1: 8, 1: 8) +     (by*NxPsi+bydx*NPsi)*wi(i1)*wi(i2)*det_Jacobi
					a_elem( 9:16, 9:16) = a_elem( 9:16, 9:16) +     (bx*NyPsi+bxdy*NPsi)*wi(i1)*wi(i2)*det_Jacobi
					a_elem( 1: 8,17:24) = a_elem( 1: 8,17:24) +     (bx*NyPsi+bxdy*NPsi)*wi(i1)*wi(i2)*det_Jacobi
					a_elem( 9:16,17:24) = a_elem( 9:16,17:24) +     (by*NxPsi+bydx*NPsi)*wi(i1)*wi(i2)*det_Jacobi
				case('pl') ! - (A_{mu p} + A_{lambda p})
					a_elem( 1: 8, 1: 8) = a_elem( 1: 8, 1: 8) - (2d0*mu+lambda)*by*NxPsi*wi(i1)*wi(i2)*det_Jacobi
					a_elem( 1: 8, 9:16) = a_elem( 1: 8, 9:16) -         lambda *by*NxPsi*wi(i1)*wi(i2)*det_Jacobi
					a_elem( 1: 8,17:24) = a_elem( 1: 8,17:24) -  2d0*mu        *bx*NyPsi*wi(i1)*wi(i2)*det_Jacobi
					a_elem( 9:16, 1: 8) = a_elem( 9:16, 1: 8) -         lambda *bx*NyPsi*wi(i1)*wi(i2)*det_Jacobi
					a_elem( 9:16, 9:16) = a_elem( 9:16, 9:16) - (2d0*mu+lambda)*bx*NyPsi*wi(i1)*wi(i2)*det_Jacobi
					a_elem( 9:16,17:24) = a_elem( 9:16,17:24) -  2d0*mu        *by*NxPsi*wi(i1)*wi(i2)*det_Jacobi
			end select
		end do
	end do

end subroutine element_a




subroutine PML_coefficient( location, a0, b0, hw, h, pml_l, ax, ay, bx, by, axdy, aydx, bxdy, bydx, a, b, c )
!###############################################################################################################################
! This subroutine returns \alpha and \beta of PML region at given 'location'
! as well as the derivatives of \alpha and \beta with respect to x and y at
! given 'location'.
! Revision information
! Arash Fathi
! Wed. OCT. 12, 2011
! Sat. Jan. 28, 2012 - MPML
! Fri. Feb. 17, 2012
! Seungbum Koo
! Mon. Mar. 26, 2012 - conversion to f90
! Mon. Jun. 09, 2014
!###############################################################################################################################

	implicit none

	!------ input ------------------------------------------------------------------------------------------------------------------
	real*8,               intent(in) :: a0, b0, hw, h, pml_l
	real*8, dimension(2), intent(in) :: location

	!------ output -----------------------------------------------------------------------------------------------------------------
	real*8, intent(out) :: ax, ay, bx, by, axdy, aydx, bxdy, bydx, a, b, c
	!-------------------------------------------------------------------------------------------------------------------------------
	! ax: alpha(x) or a(x) below,      a = a(x) * a(y)
	! ay: a(y)                         b = a(x)b(y) + a(y)b(x)
	! bx: beta(x) or b(x) below,       c = b(x) * b(y)
	! by: b(y)
	! axdy: da(x)/dy
	! aydx: da(y)/dx
	! bxdy: db(x)/dy
	! bydx: db(y)/dx
	!-------------------------------------------------------------------------------------------------------------------------------

	!------ local variables --------------------------------------------------------------------------------------------------------
	integer*4, parameter :: m=2
	real*8,    parameter :: ratio=.1d0
	real*8               :: alpha01, alpha02, alpha04
	real*8               :: beta01, beta02, beta04
	real*8               :: L_PML_1, L_PML_2, L_PML_4
	real*8               :: s01, s02, s04
	!-------------------------------------------------------------------------------------------------------------------------------
	! coordinate extention terms in PML are
	! a(s) = 1+a_0*((s-s_0)*n_s/L_PML)^m
	! b(s) =   b_0*((s-s_0)*n_s/L_PML)^m
	! MPLM terms are a(s) or b(s) times 'ratio'. n_s is +1 for 1 and 3 faces and
	! -1 for 2 and 4 faces.
	!
	! alpha0x: a_0 for face x. x varies from 1 to 4 which stand for left, right
	!       top and bottom respectively.
	! beta0x: b_0 for face x
	! L_PML_x: PML length of face x.
	! s0x: starting point location of PML of face x
	! m: polynomial degree. We choose 2 here.
	! ratio: This variable is multiplied to a(s) and b(s) to make MPML terms.
	!-------------------------------------------------------------------------------------------------------------------------------

	!------ load PML parameters ----------------------------------------------------------------------------------------------------
	!       beta0
	alpha01 = a0 ! right face
	alpha02 = a0 ! left face
	alpha04 = a0 ! bottom face

	!       beta0
	beta01  = b0
	beta02  = b0
	beta04  = b0

	!       length
	L_PML_1 = pml_l
	L_PML_2 = pml_l
	L_PML_4 = pml_l

	!       PML starting point location	  
	s01     =  hw
	s02     = -hw
	s04     = -h

	!------ figure out where does current location is in ---------------------------------------------------------------------------
	!       --------------
	!       | 5 | RD | 1 |    The numbering of PML regions are as the figure.
	!       -------------- s04
	!       | 4 | 3  | 2 |
	!       --------------
	!          s02  s01
	!-------------------------------------------------------------------------------------------------------------------------------
	! region 1
	if( location(1)>=s01 .and. location(2)>s04 ) then
		ax = 1.d0 + alpha01 * ((location(1) - s01) / L_PML_1) ** m
		ay = 1.d0 + alpha01 * ((location(1) - s01) / L_PML_1) ** m * ratio
		bx = beta01 * ((location(1) - s01) / L_PML_1) ** m
		by = bx * ratio
		axdy = .0d0
		aydx = alpha01 * m * ((location(1) - s01) / L_PML_1) ** (m-1) * ratio / L_PML_1
		bxdy = .0d0
		bydx = beta01 * m * ((location(1) - s01) / L_PML_1) ** (m-1) * ratio / L_PML_1

	! region 2
	elseif( location(1)>=s01 .and. location(2)<=s04 ) then
		ax = 1.d0 + alpha01 * ((location(1) - s01) / L_PML_1) ** m + alpha04 * ((s04 - location(2)) / L_PML_4) ** m * ratio
		ay = 1.d0 + alpha04 * ((s04 - location(2)) / L_PML_4) ** m + alpha01 * ((location(1) - s01) / L_PML_1) ** m * ratio
		bx = beta01 * ((location(1) - s01) / L_PML_1) ** m + beta04 * ((s04 - location(2)) / L_PML_4) ** m * ratio
		by = beta04 * ((s04 - location(2)) / L_PML_4) ** m + beta01 * ((location(1) - s01) / L_PML_1) ** m * ratio
		axdy = -alpha04 * m * ((s04 - location(2)) / L_PML_4) ** (m-1) * ratio / L_PML_4
		aydx =  alpha01 * m * ((location(1) - s01) / L_PML_1) ** (m-1) * ratio / L_PML_1
		bxdy =  -beta04 * m * ((s04 - location(2)) / L_PML_4) ** (m-1) * ratio / L_PML_4
		bydx =   beta01 * m * ((location(1) - s01) / L_PML_1) ** (m-1) * ratio / L_PML_1

	! region 3
	elseif( location(1)<s01 .and. location(1)>s02 .and. location(2)<=s04 ) then
		ay = 1.d0 + alpha04 * ((s04 - location(2)) / L_PML_4) ** m
		ax = 1.d0 + alpha04 * ((s04 - location(2)) / L_PML_4) ** m * ratio
		by = beta04 * ((s04 - location(2)) / L_PML_4) ** m
		bx = by * ratio
		aydx = .0d0
		axdy = -alpha04 * m * ((s04 - location(2)) / L_PML_4) ** (m-1) * ratio / L_PML_4
		bydx = .0d0
		bxdy = -beta04 * m * ((s04 - location(2)) / L_PML_4) ** (m-1) * ratio / L_PML_4

	! region 4
	elseif( location(1)<=s02 .and. location(2)<=s04 ) then
		ax = 1.d0 + alpha02 * ((s02 - location(1)) / L_PML_2)**m + alpha04 * ((s04 - location(2)) / L_PML_4)**m * ratio
		ay = 1.d0 + alpha04 * ((s04 - location(2)) / L_PML_4)**m + alpha02 * ((s02 - location(1)) / L_PML_2)**m * ratio
		bx = beta02 * ((s02 - location(1)) / L_PML_2)**m + beta04 * ((s04 - location(2)) / L_PML_4)**m * ratio
		by = beta04 * ((s04 - location(2)) / L_PML_4)**m + beta02 * ((s02 - location(1)) / L_PML_2)**m * ratio
		axdy = -alpha04 * m * ((s04 - location(2)) / L_PML_4)**(m-1) * ratio / L_PML_4
		aydx = -alpha02 * m * ((s02 - location(1)) / L_PML_2)**(m-1) * ratio / L_PML_2
		bxdy = -beta04  * m * ((s04 - location(2)) / L_PML_4)**(m-1) * ratio / L_PML_4
		bydx = -beta02  * m * ((s02 - location(1)) / L_PML_2)**(m-1) * ratio / L_PML_2

	! region 5
	elseif( location(1)<=s02 .and. location(2)>s04 ) then
		ax = 1.d0 + alpha02 * ((s02 - location(1)) / L_PML_2) ** m
		ay = 1.d0 + alpha02 * ((s02 - location(1)) / L_PML_2) ** m * ratio
		bx = beta02 * ((s02 - location(1)) / L_PML_2 ) ** m
		by = beta02 * ((s02 - location(1)) / L_PML_2 ) ** m * ratio
		axdy = .0d0
		aydx = -alpha02 * m * ((s02 - location(1)) / L_PML_2) ** (m-1) * ratio / L_PML_2
		bxdy = .0d0
		bydx = -beta02 * m * ((s02 - location(1)) / L_PML_2) ** (m-1) * ratio / L_PML_2

	! wrong location
	else
		write(*,'(4x,a85)') 'Integration point is not in PML region while calling the "PML_alpha_beta" subroutine.'
		write(*,'(4x,a18,a1,f10.4,a1,f10.4,a1)') 'current location: ', '(', location(1), ',', location(2), ')'
		stop
	end if

	a = ax * ay
	b = ax * by + ay * bx
	c = bx * by

end subroutine PML_coefficient