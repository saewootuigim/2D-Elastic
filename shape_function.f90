subroutine shape4(xi, eta, N, dN)

    implicit none
    
!------ input -----------------------------------------------------------------------------------------------------------------
    real(8), intent(in) :: xi, eta
    
!------ output ----------------------------------------------------------------------------------------------------------------
    real(8), dimension(4),   intent (out) :: N
    real(8), dimension(4,2), intent (out) :: dN
    
!------ shape functions -------------------------------------------------------------------------------------------------------

    N(1) = -.25d0 * (xi+1.d0) * (eta-1.d0)
    N(2) =  .25d0 * (xi+1.d0) * (eta+1.d0)
    N(3) = -.25d0 * (xi-1.d0) * (eta+1.d0)
    N(4) =  .25d0 * (xi-1.d0) * (eta-1.d0)

!------ derivatives of shape functions with respect to local coordinate--------------------------------------------------------
! The first column contains the derivatives with respect to xi and
! the second column contains the derivatives with respect to eta.
    dN(1,1) = -.25d0 * (eta-1.d0)
    dN(2,1) =  .25d0 * (eta+1.d0) 
    dN(3,1) = -.25d0 * (eta+1.d0)
    dN(4,1) =  .25d0 * (eta-1.d0)

    dN(1,2) = -.25d0 * (xi+1.d0)
    dN(2,2) =  .25d0 * (xi+1.d0)
    dN(3,2) = -.25d0 * (xi-1.d0)
    dN(4,2) =  .25d0 * (xi-1.d0)

end subroutine shape4




subroutine shape8(xi, eta, N, dN)
    
    implicit none
    
!------ input -----------------------------------------------------------------------------------------------------------------
    real(8), intent(in) :: xi, eta
    
!------ output ----------------------------------------------------------------------------------------------------------------
    real(8), dimension(8),   intent (out) :: N
    real(8), dimension(8,2), intent (out) :: dN
    
!------ shape functions -------------------------------------------------------------------------------------------------------

    N(1) = 0.25d0 * ( 1.0d0+xi ) * (1.0d0-eta) * (-1.0d0+xi-eta)
    N(2) = 0.25d0 * ( 1.0d0+xi ) * (1.0d0+eta) * (-1.0d0+xi+eta)
    N(3) = 0.25d0 * ( 1.0d0-xi ) * (1.0d0+eta) * (-1.0d0-xi+eta)
    N(4) = 0.25d0 * ( 1.0d0-xi ) * (1.0d0-eta) * (-1.0d0-xi-eta)
    N(5) = 0.50d0 * (1.0d0-eta*eta) * (1.0d0+xi)
    N(6) = 0.50d0 * (1.0d0-xi*xi) * (1.0d0+eta)
    N(7) = 0.50d0 * (1.0d0-eta*eta) * (1.0d0-xi)
    N(8) = 0.50d0 * (1.0d0-xi*xi) * (1.0d0-eta)

!------ derivatives of shape functions with respect to local coordinate--------------------------------------------------------
! The first column contains the derivatives with respect to xi and
! the second column contains the derivatives with respect to eta.
    dN(1,1) =  0.25d0 * (1.0d0-eta) * (2.0d0*xi-eta) 
    dN(2,1) =  0.25d0 * (1.0d0+eta) * (2.0d0*xi+eta) 
    dN(3,1) =  0.25d0 * (1.0d0+eta) * (2.0d0*xi-eta) 
    dN(4,1) =  0.25d0 * (1.0d0-eta) * (2.0d0*xi+eta) 
    dN(5,1) =  0.50d0 * (1.0d0-eta*eta)
    dN(6,1) = -1.00d0 * (1.0d0+eta)*xi 
    dN(7,1) = -0.50d0 * (1.0d0-eta*eta)
    dN(8,1) = -1.00d0 * (1.0d0-eta)*xi

    dN(1,2) =  0.25d0 * (1.0d0+xi) * (-xi+2.0d0*eta) 
    dN(2,2) =  0.25d0 * (1.0d0+xi) * ( xi+2.0d0*eta) 
    dN(3,2) =  0.25d0 * (1.0d0-xi) * (-xi+2.0d0*eta) 
    dN(4,2) =  0.25d0 * (1.0d0-xi) * ( xi+2.0d0*eta) 
    dN(5,2) = -1.00d0 * (1.0d0+xi) * eta
    dN(6,2) =  0.50d0 * (1.0d0-xi*xi)
    dN(7,2) = -1.00d0 * (1.0d0-xi)*eta
    dN(8,2) = -0.50d0 * (1.0d0-xi*xi)

end subroutine shape8




subroutine shape9(xi, eta, N, dN)

    implicit none
    
!------ input -----------------------------------------------------------------------------------------------------------------
    real(8), intent(in) :: xi, eta
    
!------ output ----------------------------------------------------------------------------------------------------------------
    real(8), dimension(9),   intent (out) :: N
    real(8), dimension(9,2), intent (out) :: dN
    
!------ shape functions -------------------------------------------------------------------------------------------------------
    N(1) = 0.25d0 *  xi * -eta * (1.0d0+xi) * (1.0d0-eta)
    N(2) = 0.25d0 *  xi *  eta * (1.0d0+xi) * (1.0d0+eta)
    N(3) = 0.25d0 * -xi *  eta * (1.0d0-xi) * (1.0d0+eta)
    N(4) = 0.25d0 * -xi * -eta * (1.0d0-xi) * (1.0d0-eta)
    N(5) = 0.50d0 *  xi *        (1.0d0+xi) *               (1.0d0-eta*eta)
    N(6) = 0.50d0 *        eta *              (1.0d0+eta) * (1.0d0-xi*xi)
    N(7) = 0.50d0 * -xi *        (1.0d0-xi) *               (1.0d0-eta*eta)
    N(8) = 0.50d0 *       -eta *              (1.0d0-eta) * (1.0d0-xi*xi) 
    N(9) = (1.0d0-xi*xi) * (1.0d0-eta*eta)

!------ derivatives of shape functions with respect to local coordinate--------------------------------------------------------
! The first column contains the derivatives with respect to xi and
! the second column contains the derivatives with respect to eta.
    dN(1,1) = 0.25d0 * -eta * (1.0d0-eta) * ( 1.0d0 + 2.0d0 * xi)  
    dN(2,1) = 0.25d0 *  eta * (1.0d0+eta) * ( 1.0d0 + 2.0d0 * xi)
    dN(3,1) = 0.25d0 *  eta * (1.0d0+eta) * (-1.0d0 + 2.0d0 * xi)
    dN(4,1) = 0.25d0 * -eta * (1.0d0-eta) * (-1.0d0 + 2.0d0 * xi)
    dN(5,1) = 0.50d0 * (1.0d0-eta*eta) * ( 1.0d0 + 2.0d0 * xi)
    dN(6,1) =  eta * (1.0d0+eta) * -xi 
    dN(7,1) = 0.50d0 * (1.0d0-eta*eta) * (-1.0d0 + 2.0d0 * xi)
    dN(8,1) = -eta * (1.0d0-eta) * -xi
    dN(9,1) = -2.0d0 * xi * (1.0d0-eta*eta) 

    dN(1,2) = 0.25d0 *  xi * (1.0d0+xi) * (-1.0d0 + 2.0d0 * eta)
    dN(2,2) = 0.25d0 *  xi * (1.0d0+xi) * ( 1.0d0 + 2.0d0 * eta)
    dN(3,2) = 0.25d0 * -xi * (1.0d0-xi) * ( 1.0d0 + 2.0d0 * eta)
    dN(4,2) = 0.25d0 * -xi * (1.0d0-xi) * (-1.0d0 + 2.0d0 * eta)
    dN(5,2) =  xi * (1.0d0+xi) * -eta
    dN(6,2) = 0.50d0 * (1.0d0-xi*xi) * ( 1.0d0 + 2.0d0 * eta)
    dN(7,2) = -xi * (1.0d0-xi) * -eta
    dN(8,2) = 0.50d0 * (1.0d0-xi*xi) * (-1.0d0 + 2.0d0 * eta)
    dN(9,2) = -2.0d0 * eta * (1.0d0-xi*xi)

end subroutine shape9