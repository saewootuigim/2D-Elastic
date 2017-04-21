module global_variable

!-- Include PETSc header files. -----------------------------------------------80
#include "petsc/finclude/petscsys.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"
!------------------------------------------------------------------------------80

! PETSc variables
    integer*4 :: ierr
    Mat       :: M, C, K, Ms, Ks
    Vec       :: F

! math constants
    real*8,    parameter :: pi=acos(-1d0)

! strings to control the simulation
    integer*4, parameter :: NULL=0
    integer*4, parameter :: X_DIR=1,   Y_DIR=2, XY_DIR=3
    integer*4, parameter :: RICKER=1,  GAUSS=2, RECT=3, DGAUSS=4
    integer*4, parameter :: NARROW=1,  ALL=2
    integer*4, parameter :: G_L=1,     G_K=2,     G_Lbt=3
    integer*4, parameter :: FORWARD=1, TR=2
    integer*4, parameter :: DISP=1,    VEL=2,   ACCEL=3
    integer*4, parameter :: EQ_FORCE=1, EQ_RESPONSE=2
    integer*4, parameter :: YES=1, NO=2

! for timing
    integer*4 :: count1, count2, count3, count4, count_rate, count_max

! quadrature nodes and weights
    integer*4 :: NQuad
    integer*4 :: quad_meth
    real*8, allocatable, dimension(:) :: xi, wi

! model related
    integer*4 :: model_num
    integer*4 :: nNode, nElem
    integer*4 :: nPhyProp
    real*8    :: pml_l, rdwd, rdht

    integer*4 :: nDOF
    integer*4 :: nDOF_disp
    integer*4 :: nDOF_disp_surf
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

! source simulation related
    integer*4 :: src_sim
    integer*4 :: nTstep_src

    character(12) dir_sim_in
    character(20) dir_sim_out

! simulation related
    integer*4 :: F_or_T

    integer*4                              :: fw_src_sgnl
    integer*4                              :: fw_nPntLd
    integer*4                              :: fw_PntLd_dir
    integer*4                              :: fw_load_intensity
    integer*4, allocatable, dimension(:)   :: fw_load_node
    integer*4, allocatable, dimension(:,:) :: fw_load_DOF
    real*8,    allocatable, dimension(:,:) :: fw_load_coord

    real*8    :: dt
    real*8    :: Tsim
    integer*4 :: nTstep

    integer*4 :: response
    real*8    :: f0
    real*8    :: Tload

    real*8    :: a0, b0

    integer*4 :: TR_DtN
    integer*4 :: TR_surf_trac_case
    integer*4 :: TR_surf_trac_dir
    integer*4 :: TR_surf_trac_resp
    integer*4 :: TR_apply_sink
    integer*4 :: TR_apply_surf_type

! Newmark solver related
    real*8    :: beta, gamma

! IO related
    integer*4 :: fid_out_disp,      fid_out_vel,      fid_out_accel
    integer*4 :: fid_out_disp_surf, fid_out_vel_surf, fid_out_accel_surf
    integer*4 :: fid_in_surf
    integer*4 :: fid_out_src

end module global_variable
