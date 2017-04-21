subroutine newmark()

    use global_variable
    
    implicit none

!-- Include PETSc header files. -----------------------------------------------80
#include "petsc/finclude/petscsysdef.h"
#include "petsc/finclude/petscvecdef.h"
#include "petsc/finclude/petscmatdef.h"
#include "petsc/finclude/petscvec.h90"
#include "petsc/finclude/petscksp.h"
#include "petsc/finclude/petscpc.h"
#include "petsc/finclude/petscis.h"
!------------------------------------------------------------------------------80

!-- Local variables -----------------------------------------------------------80
    integer*4 :: i1, n
    integer*4 :: fid1, fid2, fid3, fid4, rowLen
    integer*4 :: dont_print
    real*8    :: a1, a2, a3, a4, a5, a6, t, fac0
    real*8    :: T_src_sim, T_cur_sim
    
    integer*4, dimension(fw_nPntLd*2) :: DOF_src
    real*8, dimension(nDOF_disp_surf) :: surf_source
    real*8, dimension(nDOF_disp)      :: output_disp
    real*8, dimension(fw_nPntLd*2)    :: output_src
    PetscScalar, pointer :: u1_disp(:), u2_vel(:), u3_accel(:)
PetscInt high, low
    Mat :: K_bar
    Vec :: F_eff, u1, u2, u3, u1_old, u2_old, u3_old, temp, u_srf_k
    IS  :: perm
    MatFactorInfo :: info(mat_factorinfo_size)
!------------------------------------------------------------------------------80

    ! constants needed in this subroutine
    if( F_or_T==TR ) then
        T_src_sim = nTstep_src*dt
        T_cur_sim = nTstep*dt
    endif
    dont_print = 0
    fac0 = 0d0

    print *
    print '(1x,a)', 'begin solving'

    ! Create PETSc objects.
    call VecCreateMPI( PETSC_COMM_WORLD, PETSC_DECIDE, nDOF, F_eff,  ierr )
    call VecCreateMPI( PETSC_COMM_WORLD, PETSC_DECIDE, nDOF, u1,     ierr )
    call VecCreateMPI( PETSC_COMM_WORLD, PETSC_DECIDE, nDOF, u2,     ierr )
    call VecCreateMPI( PETSC_COMM_WORLD, PETSC_DECIDE, nDOF, u3,     ierr )
    call VecCreateMPI( PETSC_COMM_WORLD, PETSC_DECIDE, nDOF, temp,   ierr )
    call VecCreateMPI( PETSC_COMM_WORLD, PETSC_DECIDE, nDOF, u1_old, ierr )
    call VecCreateMPI( PETSC_COMM_WORLD, PETSC_DECIDE, nDOF, u2_old, ierr )
    call VecCreateMPI( PETSC_COMM_WORLD, PETSC_DECIDE, nDOF, u3_old, ierr )
    call MatCreateAIJ( PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nDOF, &
                       nDOF, 100, PETSC_NULL_INTEGER, 100, PETSC_NULL_INTEGER, K_bar, ierr )

    call VecSetOption( F_eff, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE, ierr )

    ! Do pre-processing to apply Dirichlet boundary condition.
    if( F_or_T==TR .and. TR_DtN==NO ) then
        call VecCreateMPI( PETSC_COMM_WORLD, PETSC_DECIDE, nDOF_disp_surf, u_srf_k, ierr )
        call VecSetOption( u_srf_k, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE, ierr )
        ! call VecDuplicate( u_srf_k, u_srf_m, ierr )
    endif

    ! Apply the initial contidion. (u, v and a are initially all at rest.)
    call VecSet( u1, 0d0, ierr )
    call VecSet( u2, 0d0, ierr )
    call VecSet( u3, 0d0, ierr )

    ! Define the Newmark method coefficients.
    a1 = 1d0 / (beta*dt**2)
    a2 = 1d0 / (beta*dt)
    a3 = 1d0 / (2d0*beta) - 1d0
    a4 = gamma / (beta*dt)
    a5 = gamma / beta - 1d0
    a6 = dt * (gamma/(2d0*beta)-1d0)

    ! Form the effective stiffness, say, K = K + b1 * M + b4 * C.
    call MatAXPY( K, a1, M, DIFFERENT_NONZERO_PATTERN, ierr )
    call MatAXPY( K, a4, C, DIFFERENT_NONZERO_PATTERN, ierr )

    ! Perform LU factorization on K.
    print '(1x,a)', 'LU factorization'
    call system_clock( count3 )

    call MatGetFactor(K, 'superlu_dist', MAT_FACTOR_LU, K_bar, ierr)
    call MatGetSize(K, n, PETSC_NULL_INTEGER, ierr)
    call ISCreateStride(PETSC_COMM_WORLD, n, 0, 1, perm, ierr)
    call MatLUFactorSymbolic(K_bar, K, perm, perm, info, ierr)
    call MatLUFactorNumeric(K_bar, K, info, ierr)

    call system_clock( count4 )
    print '(1x,a,f8.4,a)', 'elapse: ', dble(count4-count3)/dble(count_rate), 'sec.'
    print *

    ! OUTPUT FILES
    if( response==DISP ) then
        ! if( F_or_T==TR .and. TR_DtN==NO ) then
        !     inquire( iolength=rowLen ) surf_source, output_disp
        ! else
        !     inquire( iolength=rowLen ) output_disp
        ! endif
        inquire( iolength=rowLen ) output_disp
        open(newunit=fid1,file=trim(adjustl(dir_sim_out))//'disp.dat',status='replace',action='write',form='unformatted',access='direct',recl=rowLen)
        inquire( iolength=rowLen ) output_src
        open(newunit=fid2,file=trim(adjustl(dir_sim_out))//'disp_src.dat',status='replace',action='write',form='unformatted',access='direct',recl=rowLen)
        print '(4x,3a,i8)', '       output file opened. ', trim(adjustl(dir_sim_out))//'disp.dat', ', file id=', fid1
        print '(4x,3a,i8)', 'output source file opened. ', trim(adjustl(dir_sim_out))//'disp_src.dat', ', file id=', fid2

    elseif( response==VEL ) then
        inquire( iolength=rowLen ) output_disp
        open(newunit=fid1,file=trim(adjustl(dir_sim_out))//'vel.dat',status='replace',action='write',form='unformatted',access='direct',recl=rowLen)
        inquire( iolength=rowLen ) output_src
        open(newunit=fid2,file=trim(adjustl(dir_sim_out))//'vel_src.dat',status='replace',action='write',form='unformatted',access='direct',recl=rowLen)
        print '(4x,3a,i8)', '       output file opened. ', trim(adjustl(dir_sim_out))//'vel.dat', ', file id=', fid1
        print '(4x,3a,i8)', 'output source file opened. ', trim(adjustl(dir_sim_out))//'vel_src.dat', ', file id=', fid2

    elseif( response==ACCEL ) then
        inquire( iolength=rowLen ) output_disp
        open(newunit=fid1,file=trim(adjustl(dir_sim_out))//'accel.dat',status='replace',action='write',form='unformatted',access='direct',recl=rowLen)
        inquire( iolength=rowLen ) output_src
        open(newunit=fid2,file=trim(adjustl(dir_sim_out))//'accel_src.dat',status='replace',action='write',form='unformatted',access='direct',recl=rowLen)
        print '(4x,3a,i8)', '        output file opened. ', trim(adjustl(dir_sim_out))//'accel.dat', ', file id=', fid1
        print '(4x,3a,i8)', ' output source file opened. ', trim(adjustl(dir_sim_out))//'accel_src.dat', ', file id=', fid2
    endif

    ! SURFACE FILES
    inquire( iolength=rowLen ) surf_source
    if( F_or_T==FORWARD ) then
        if( response==DISP ) then
            open(newunit=fid3,file=trim(adjustl(dir_sim_out))//'disp_surf.dat',status='replace',action='write',form='unformatted',access='direct',recl=rowLen)
            print '(3x,3a,i8)', 'output surface file opened. ', trim(adjustl(dir_sim_out))//'disp_surf.dat', ', file id=', fid3

        elseif( response==VEL ) then
            open(newunit=fid3,file=trim(adjustl(dir_sim_out))//'vel_surf.dat',status='replace',action='write',form='unformatted',access='direct',recl=rowLen)
            print '(3x,3a,i8)', 'output surface file opened. ', trim(adjustl(dir_sim_out))//'vel_surf.dat', ', file id=', fid3

        elseif( response==ACCEL ) then
            open(newunit=fid3,file=trim(adjustl(dir_sim_out))//'accel_surf.dat',status='replace',action='write',form='unformatted',access='direct',recl=rowLen)
            print '(3x,3a,i8)', 'output surface file opened. ', trim(adjustl(dir_sim_out))//'accel_surf.dat', ', file id=', fid3

        endif

    elseif( F_or_T==TR ) then
        if( TR_surf_trac_resp==DISP ) then
            open(newunit=fid3,file=trim(adjustl(dir_sim_in))//'disp_surf.dat',status='old',action='read',form='unformatted',access='direct',recl=rowLen)
            print '(4x,3a,i8)', 'input surface file opened. ', trim(adjustl(dir_sim_in))//'disp_surf.dat', ', file id=', fid3
            open(newunit=fid4,file=trim(adjustl(dir_sim_out))//'disp_surf.dat',status='replace',action='write',form='unformatted',access='direct',recl=rowLen)

        elseif( TR_surf_trac_resp==VEL ) then
            open(newunit=fid3,file=trim(adjustl(dir_sim_in))//'vel_surf.dat',status='old',action='read',form='unformatted',access='direct',recl=rowLen)
            print '(4x,3a,i8)', 'input surface file opened. ', trim(adjustl(dir_sim_in))//'vel_surf.dat', ', file id=', fid3
            open(newunit=fid4,file=trim(adjustl(dir_sim_out))//'vel_surf.dat',status='replace',action='write',form='unformatted',access='direct',recl=rowLen)

        elseif( TR_surf_trac_resp==ACCEL ) then
            open(newunit=fid3,file=trim(adjustl(dir_sim_in))//'accel_surf.dat',status='old',action='read',form='unformatted',access='direct',recl=rowLen)
            print '(4x,3a,i8)', 'input surface file opened. ', trim(adjustl(dir_sim_in))//'accel_surf.dat', ', file id=', fid3
            open(newunit=fid4,file=trim(adjustl(dir_sim_out))//'accel_surf.dat',status='replace',action='write',form='unformatted',access='direct',recl=rowLen)

        endif

    endif

    ! DOFs of source location
    do i1 = 1, fw_nPntLd
        DOF_src(2*i1-1:2*i1) = fw_load_DOF(i1,:)
    enddo
    print '(1x,a,<2*fw_nPntLd>(1x,i7))', 'DOF source: ', DOF_src

    ! Start Newmark time integration.
    print '(1x,a)', 'starting time integration'
    if( F_or_T==FORWARD ) then
        print '(4x,a,f6.3,a5)', 'simulation length: ', Tsim, ' sec.'
    elseif( F_or_T==TR ) then
        print '(4x,a,f6.3,a5)', 'simulation length of the source simulation: ', T_src_sim, ' sec.'
        print '(4x,a,f6.3,a5)', 'simulation length of current TR simulation: ', T_cur_sim, ' sec.'
    endif

    do i1 = 1, nTstep

        if( mod(i1,10)==1 ) call system_clock( count3 )

      ! u1, u2 and u3 are copied into u1_old, u2_old and u3_old.
      ! These old terms are needed when retrieving u2_new and u3_new.
        call VecCopy( u1, u1_old, ierr )
        call VecCopy( u2, u2_old, ierr )
        call VecCopy( u3, u3_old, ierr )

      ! Calculate the time of current time step.
        t = dble(i1)*dt

      ! Initialize the effective force vector.
        call VecSet( F_eff, 0d0, ierr )

      ! Compute load for the forward problem.
        if( F_or_T==FORWARD ) then

            if( fw_src_sgnl==RICKER ) then
                if( t<=Tload ) then
                    fac0 = (1d0-2d0*pi**2*(f0*t-1d0)**2)*exp(-pi**2*(f0*t-1d0)**2)
                else
                    fac0 = 0d0
                endif
            elseif( fw_src_sgnl==GAUSS ) then
                if( t<=Tload ) then
                    fac0 = (exp(-(t-.1d0)**2/.0001d0)-exp(-.01d0/.0001d0))/(1d0-exp(-.01d0/.0001d0))
                else
                    fac0 = 0d0
                endif
            elseif( fw_src_sgnl==RECT ) then
                if( t>f0 .and. t<=f0+Tload ) then
                    fac0 = 1D0
                else
                    fac0 = 0D0
                endif
            elseif( fw_src_sgnl==DGAUSS ) then
                if( t<=Tload ) then
                    fac0 = -233.188D0*exp(-1D4*(-0.1D0+t)**2)*(-0.1D0+t)
                else
                    fac0 = 0D0
                endif
            endif

            call VecAXPY( F_eff, fac0, F, ierr )

      ! Deal with load for the TR problem.
        elseif( F_or_T==TR ) then

! =============================================================================
! To count the effect of evanescent wave, apply the same load at the source.
! =============================================================================

            ! if there is sink
            if( TR_apply_sink==YES ) then
                
                if( dont_print==0 ) then
                    write(*,'(1x,a)',advance='no') 'Sink is activated.'
                    dont_print = 1
                endif

                if( t<T_src_sim-Tload ) then
                    fac0 = 0d0
                elseif( t>T_src_sim ) then
                    if( dont_print==3 ) then
                        write(*,'(f7.3)') 'sink deactivated, t=', t
                        dont_print = 4
                    endif
                    fac0 = 0d0
                else
                    if( dont_print==2 ) then
                        write(*,'(a3,f7.3)',advance='no') ' t=', t
                        dont_print = 3
                    endif
                    if( fw_src_sgnl==RICKER ) then
                        fac0 = (1d0-2d0*pi**2*(f0*(T_src_sim-t)-1d0)**2)*exp(-pi**2*(f0*(T_src_sim-t)-1d0)**2)
                    elseif( fw_src_sgnl==GAUSS ) then
                        fac0 = (exp(-((T_src_sim-t)-.1d0)**2/.0001d0)-exp(-.01d0/.0001d0))/(1d0-exp(-.01d0/.0001d0))
                    elseif( fw_src_sgnl==RECT ) then
                        fac0 = 1d0
                    endif
                endif
    
                ! F_eff = F_eff+fac0*F
                if( dont_print==1 ) then
                    if( fw_src_sgnl==RICKER ) then
                        print '(1x,a)', 'as Ricker.'
                    elseif( fw_src_sgnl==GAUSS ) then
                        print '(1x,a)', 'as Gauss.'
                    elseif( fw_src_sgnl==RECT ) then
                        print '(1x,a)', 'as rectagular.'
                    elseif( fw_src_sgnl==DGAUSS ) then
                        print '(1x,a)', 'as derivative of Gauss.'
                    endif
                    dont_print = 2
                endif
                call VecAXPY( F_eff, fac0, F, ierr )

            endif
! =============================================================================

            ! Read surface displacement in time reversed order.
            if( TR_apply_surf_type/=1 ) then

                ! Apply Dirichlet boundary condition.
                if( TR_DtN==NO ) then
                    if( dont_print/=1 ) print '(1x,a)', 'Surface Dirichlet datum is applied.'
                    if( i1<=nTstep_src ) then
                        read(fid3,rec=nTstep_src-i1+1) surf_source
                    elseif( i1>nTstep_src ) then
                        surf_source = 0d0
                    endif

                    call VecSetValues( u_srf_k, nDOF_disp_surf, DOF_surf-1, -surf_source, INSERT_VALUES, ierr )
                    call VecAssemblyBegin( u_srf_k, ierr )
                    call VecAssemblyEnd  ( u_srf_k, ierr )
                    call MatMult( Ks, u_srf_k, F_eff, ierr )

                    if( dont_print/=1 ) dont_print = 1

                ! Apply Neumann boundary condition.
                elseif( TR_DtN==YES ) then
                    if( dont_print/=1 ) print '(1x,a)', 'Surface traction is applied.'
                    if( i1<=nTstep_src ) then
                        read(fid3,rec=nTstep_src-i1+1) surf_source
                        if( TR_apply_surf_type==3 ) then ! yes and rescale
                            if( dont_print/=1 ) then
                                print '(4x,a)', 'Surface traction is rescaled to induce the same source location displacement response'
                                print '(4x,a)', 'intensity as the source simulation.'
                                dont_print = 1
                            endif
                            if( src_sim==4 ) then
                                surf_source = surf_source*0.011583337393109d0
                            elseif( src_sim==5 ) then
                                surf_source = surf_source*0.009706135131560d0
                            elseif( src_sim==6 ) then
                                surf_source = surf_source*0.01308432666d0
                            elseif( src_sim==7 ) then
                                surf_source = surf_source*0.014203630499858d0
                            else
                                print *, 'ERROR: newmark.F90, line 258'
                                stop
                            endif
! =============================================================================
! The number below is multiplied to "surf_source" in order to make the source
! location displacement response to have the same maximum amplitude to the 
! forward simulation.
! 
! for those simulations use F004 as their source simulation: 0.011583337393109d0
! for those simulations use F005 as their source simulation: 0.009706135131560d0
! =============================================================================
                        endif
                    elseif( i1>nTstep_src ) then
                        surf_source = 0d0
                    endif

                    ! Assign it into force vector F. Only DOFs in x-direction are assigned.
                    call VecSetValues( F_eff, nDOF_disp_surf, DOF_surf-1, surf_source, INSERT_VALUES, ierr )
                    call VecAssemblyBegin( F_eff, ierr )
                    call VecAssemblyEnd  ( F_eff, ierr )

                    if( dont_print/=1 ) dont_print = 1
                endif
            endif

        end if

      ! Calculate effective load vector.
      ! F_eff = F + M (a1 * u1 + a2 * u2 + a3 * u3) + C (a4 * u1 + a5 * u2 + a6 * u3)
        call VecSet( temp, 0d0, ierr )
        call VecAXPY( temp, a1, u1, ierr )
        call VecAXPY( temp, a2, u2, ierr )
        call VecAXPY( temp, a3, u3, ierr )
        call MatMultAdd( M, temp, F_eff, F_eff, ierr )

        call VecSet( temp, 0d0, ierr )
        call VecAXPY( temp, a4, u1, ierr )
        call VecAXPY( temp, a5, u2, ierr )
        call VecAXPY( temp, a6, u3, ierr )
        call MatMultAdd( C, temp, F_eff, F_eff, ierr )

      ! Solve it to obtain the displacement of current time step.
        call MatSolve( K_bar, F_eff, u1, ierr )

      ! Retrieve the velocity at the current step (i-th step).
      ! u2_new = a4 * (u1_new - u1_old) - a5 * u2_old - a6 * u3_old
        call VecSet( temp, 0d0, ierr )
        call VecAXPY( temp, a4, u1, ierr )
        call VecAXPY( temp, -a4, u1_old, ierr )
        call VecAXPY( temp, -a5, u2_old, ierr )
        call VecWAXPY( u2, -a6, u3_old, temp, ierr )

      ! Retrieve the acceleration at the current step (i-th step).
      ! u3_new = a1 * (u1_new - u1_old) - a2 * u2_old - a3 * u3_old
        call VecSet(temp, 0d0, ierr)
        call VecAXPY(temp, a1, u1, ierr)
        call VecAXPY(temp, -a1, u1_old, ierr)
        call VecAXPY(temp, -a2, u2_old, ierr)
        call VecWAXPY(u3, -a3, u3_old, temp, ierr)

      ! Print results.
        if( response==DISP ) then
            call VecGetArrayF90( u1, u1_disp, ierr )
            write(fid1,rec=i1) u1_disp(DOF_disp)
            write(fid2,rec=i1) u1_disp(DOF_src)
            if( F_or_T==FORWARD ) then
                write(fid3,rec=i1) u1_disp(DOF_surf)
            elseif( F_or_T==TR .and. TR_DtN==NO ) then
                write(fid4,rec=i1)  surf_source
            else
                write(fid4,rec=i1) u1_disp(DOF_surf)
            endif
            call VecRestoreArrayF90( u1, u1_disp, ierr )
        elseif( response==VEL ) then
            call VecGetArrayF90( u2, u2_vel, ierr )
            write(fid1,rec=i1) u2_vel(DOF_disp)
            write(fid2,rec=i1) u2_vel(DOF_src)
            if( F_or_T==FORWARD ) then
                write(fid3,rec=i1) u1_disp(DOF_surf)
            elseif( F_or_T==TR .and. TR_DtN==NO ) then
                write(fid4,rec=i1)  surf_source
            else
                write(fid4,rec=i1) u2_vel(DOF_surf)
            endif
            call VecRestoreArrayF90( u2, u2_vel, ierr )
        elseif( response==ACCEL ) then
            call VecGetArrayF90( u3, u3_accel, ierr )
            write(fid1,rec=i1) u3_accel(DOF_disp)
            write(fid2,rec=i1) u3_accel(DOF_src)
            if( F_or_T==FORWARD ) then
                write(fid3,rec=i1) u1_disp(DOF_surf)
            elseif( F_or_T==TR .and. TR_DtN==NO ) then
                write(fid4,rec=i1)  surf_source
            else
                write(fid4,rec=i1) u3_accel(DOF_surf)
            endif
            call VecRestoreArrayF90( u3, u3_accel, ierr )
        endif

      ! Printing progress.
        if( mod(i1,10)==0 ) then
            call system_clock( count4 )
            if( F_or_T==forward .or. ( F_or_T==TR .and. TR_apply_sink==YES ) ) then
                print '(a,i4,a,i4,a,f6.2,a,es12.4)', 'step ', i1, '/', nTstep, ', elapse=', dble(count4-count3)/dble(count_rate), ' sec., fac0=', fac0
            elseif( F_or_T==TR .and. TR_apply_sink==NO ) then
                print '(a,i4,a,i4,a,f6.2,a)', 'step ', i1, '/', nTstep, ', elapse=', dble(count4-count3)/dble(count_rate), ' sec.'
            endif
        endif

    end do

    ! Close output files.
    close( fid1 )
    close( fid2 )
    close( fid3 )
    close( fid4 )

  ! Destroy all the PETSc objects used in this subroutine.
    call MatDestroy(  K_bar, ierr )
    call VecDestroy(  F_eff, ierr )
    call VecDestroy(     u1, ierr )
    call VecDestroy(     u2, ierr )
    call VecDestroy(     u3, ierr )
    call VecDestroy( u1_old, ierr )
    call VecDestroy( u2_old, ierr )
    call VecDestroy( u3_old, ierr )
    call VecDestroy(   temp, ierr )
    call ISDestroy(    perm, ierr )

    if( F_or_T==TR .and. TR_DtN==NO ) then
        ! call VecDestroy( u_srf_m, ierr )
        call VecDestroy( u_srf_k, ierr )
    endif

end subroutine newmark