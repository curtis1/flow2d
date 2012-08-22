!********************************************************************
!*                                                                  *
!*                                                          MPXLIB  *
!*                                                                  *
!********************************************************************
!* Author: Nathaniel R. Morgan (2004)                               *
!*                                                                  *
!* MPXLIB is a (M)ulti(P)hase dynamics(X) (LIB)ary for modeling     *
!* 2-phase flow problems.  MPXLIB solves the full transient         *
!* Navier-Stokes equations and has the ability to capture the       *
!* interface between two fluids of different properties using the   *
!* Level Set (LS) method.                                           *    
!*                                                                  *
!*                                                                  *
!* The staggered grid is shown below:                               *
!*                                                                  *
!* (bottom left corner)                                             *
!*                                                                  *
!* v----o----v----o  y(2)                                           *
!*      |         |                                                 *
!* P    u    P    u                                                 *
!*      |         |                                                 *
!* v----o----v----o  y(1)                                           *
!*            |         |                                           *
!* P    u    P    u                                                 *
!*                                                                  *
!*     x(1)     x(2)                                                *
!*                                                                  *
!********************************************************************

! NOTES:
! Modified by Jessica Sanders Summer, 2008
! Modified by Curtis Lee, Summer 2010
PROGRAM MPXLIB

        IMPLICIT NONE
        !
        ! INTEGERS 
        ! JDS 7/17/08 Removing temperature variables
        integer :: Nx,Ny,i,j,step,r,l,k                             ! loop integers
        integer :: ios                                              ! input file
        integer :: file_count,count,dumpcount                       ! counting integers
        integer :: iter,itermax,iterprint                           ! iteration integers, not used
        integer :: MBDCND,N,M,NBDCND,IDIMF,IERROR                   ! for fish pack solver
        integer :: solvertype,precond
        integer :: IC
        integer :: xdirgo,ydirgo
        integer :: BCnorth,BCsouth,BCeast,BCwest, &
                           BCPnorth,BCPsouth,BCPeast,BCPwest, &
                           BCphisouth,BCphinorth,BCphiwest,&
                           BCphieast                                ! BC types - wall,flow,etc
	    
	    ! Coefficient arrays for GALS routines
	    integer, dimension(5,6) :: GA
	    integer :: width
	    real(kind=8), dimension(6) :: scal,herm
		
		real(kind=8) :: volume, volume0								! volume tracking for GALS
                           
        ! FOR RIGID BODY
        integer :: rb_shape                                         ! Shape of the rigid body
                                                                    ! 1 = circle
                                                                    ! 2 = ellipse
                                                                    ! 3 = square
        integer :: LS_type                                          ! Level Set Method Type
                                                                    ! 1 = Standard LS
                                                                    ! 2 = CLSVOF
                                                                    ! 3 = Gradient-Augmented L.S.
        integer :: int_type                                         ! Type of internal component(s)
                                                                    ! 0 = none
                                                                    ! 1 = torsional disc
                                                                    ! 2 = 2D trans
                                                                    ! 3 = 1D trans
        integer :: RB_lock                                          ! Restrain motion of R
                                                                    ! 0 = No restraints
                                                                    ! 1 = Vertical motion only

        ! LOGICAL AND CHARACTER
        logical :: startup,trackphase,immersedbc                    ! booleans  
        logical :: rigidbody                                        ! Include rigid body
        logical :: wedge,beach,drop,reef                            ! Wedge
        logical :: dambreak, wave, slosh                            ! Benchmarks
        logical :: inside                                           ! dummy variable
        character(40) :: filename,variables,startfile               ! strings
        !
        !
        ! REALS
        real(kind=8) :: pi                                          ! pi
        real(kind=8) :: Lx,Ly                                       ! domain size
        real(kind=8) :: umax,vmax,dtvisc,dtadv,&
                        dtg,dtsigma,dtinternal                      ! time step variables
        real(kind=8) :: deltaT,hx,hy,deltaTau               ! incremenit
        real(kind=8) :: time,time_max, &
                        time_print,time_count                       ! time variables
        real(kind=8) :: freq,wedgeT									! wedge frequency in Hz (Ex: 0.75)
        real(kind=8) :: stability,reintime
        real(kind=8) :: rho_one,rho_two,rho_three,mu_one,&        
                        mu_two,mu_three,sigma                       ! properties
        real(kind=8) :: kappa_n,kappa_s,kappa_w,kappa_e             ! interface geometry
        real(kind=8) :: gx,gy,Fsigma                                ! body forces
        real(kind=8) :: Tolerance                                   ! iteration variables
        real(kind=8) :: unorth,usouth,uwest,ueast,&
                        vnorth,vsouth,vwest,veast                                
        real(kind=8) :: Pnorth,Psouth,Pwest,Peast                   ! dirclet pressure BC's
        real(kind=8) :: contC,contA,contR,contM                     ! contact line BC's
        real(kind=8) :: u_E,u_W,u_N,u_S,&
                        u_ne,u_se,u_nw,u_sw,u_O, &
                        v_E,v_W,v_N,v_S, &
                        v_ne,v_nw,v_se,v_sw,v_O                     ! velocities
        real(kind=8) :: u_xe,u_xw,u_yn,u_ys,&
                        v_xn,v_xs,u_xO,u_yO                         ! velocity derivatives
        real(kind=8) :: u_ye,u_yw,v_yn,v_ys,&
                        v_xe,v_xw,v_xO,v_yO                         ! velocity derivatives
        real(kind=8) :: mu_n,mu_s,mu_e,mu_w, &                        
                        rho_O,rho_n,rho_s,rho_e,rho_w               ! properties
        real(kind=8) :: u_plus,v_plus                               ! upwinding terms
        real(kind=8) :: u_minus,v_minus                             ! upwinding terms
        real(kind=8) :: dudx_plus,dudx_minus                        ! upwinding terms
        real(kind=8) :: dudy_plus,dudy_minus                        ! upwinding terms
        real(kind=8) :: dvdx_plus,dvdx_minus                        ! upwinding terms
        real(kind=8) :: dvdy_plus,dvdy_minus                        ! upwinding terms
        real(kind=8) :: Visc                                        ! viscous or diffusion term
        real(kind=8) :: DiffX, DiffY                                ! derivatives in the respective directions
        real(kind=8) :: A,B,C,D,PERTRB,ELMBDA                       ! for fishpac solver
        real(kind=8) :: an,as,aw,ae,aO,omega                        ! for SOR
        real(kind=8) :: Pavg,PavgOld,PRes,maxPRes                   ! pressure variables (SOR)
        real(kind=8) :: deltau,deltav                               ! change in u and v (immersed boundary)                                                        
        real(kind=8) :: xcenter,ycenter                             ! dummy variables for ICs
        real(kind=8) :: a_phiLS,b_phiLS,c_phiLS,d_phiLS,&
                        e_phiLS,f_phiLS,v_phiLS
        real(kind=8) :: waveh,tw,sw
        real(kind=8) :: MassInt
        !
        real(kind=8) :: alpha_1,alpha_2,beta_1,beta_2
        real(kind=8) :: Coef(4)
        real(kind=8) :: normal(2)
        
        ! Benchmarking Variables
        real(kind=8), DIMENSION(2) :: testa, testb, tempa, tempb
        integer :: maxi, maxj
        
        ! For Calculating Energy/Momentum
        real(kind=8) :: total_energy                                ! Total system energy
        real(kind=8) :: lin_momx                                    ! Momentum in x-direction
        real(kind=8) :: lin_momy                                    ! Momentum in y-direction
        real(kind=8) :: net_mom                                     ! Net momentum
        real(kind=8) :: gradnorm

        ! Rigid Body Variables
        real(kind=8) :: a_rbLS, b_rbLS, c_rbLS, d_rbLS              ! Rigid body geometry
        
        real(kind=8) :: rb1_mass,buoy_J                              ! Rigid body mass
        real(kind=8) :: rb_mass                                     ! TOTAL mass (RB + Internals)
        
        real(kind=8) :: n_rb_dx,n_rb_dy                             ! Rigid body displacement        
        real(kind=8) :: n1_rb_dx,n1_rb_dy
        real(kind=8) :: rx,ry
        
        real(kind=8) :: n_w_dy                                      ! Wedge velocities
        real(kind=8) :: n1_w_dy        
        real(kind=8) :: n12_w_vy
        real(kind=8) :: n1_w_vy        
        
        real(kind=8) :: n_rb_vx,n_rb_vy                             ! Rigid body velocities
        real(kind=8) :: n12_rb_vx,n12_rb_vy
        real(kind=8) :: n1_rb_vx,n1_rb_vy
        
        real(kind=8) :: n_rb_ax,n_rb_ay                             ! Rigid body acceleration        
        real(kind=8) :: n1_rb_ax,n1_rb_ay

        real(kind=8) :: big_theta       
        real(kind=8) :: th_change,th_domain,th_old,th_counter       ! Change in RB's angle (from original value)
        real(kind=8) :: n_rb_vom,n12_rb_vom                         ! Rotational velocities
        real(kind=8) :: n1_rb_vom
        real(kind=8) :: n_rb_aom,n1_rb_aom                          ! Rotational acceleration
        real(kind=8) :: xc,yc                                       ! Center of mass (Actually, geometric center... -CL)
        real(kind=8) :: u_12,v_12                                   ! Intermediate variables for setting
                                                                    ! fluid velocity to RB velocity
        real(kind=8) :: adj_factor

        ! Internal Disc Variables               
        real(kind=8) :: small_theta_n, small_theta_n1        
        real(kind=8) :: n_rb2_vom,n12_rb2_vom                       ! Rotational velocities
        real(kind=8) :: n1_rb2_vom
        real(kind=8) :: n_rb2_aom,n1_rb2_aom                        ! Rotational acceleration
        
        real(kind=8) :: delta_theta                                 ! Angle difference between disc and RB
        real(kind=8) :: delta_vom                                   ! Angular velocity difference
        real(kind=8) :: int_torque                                  ! Internal torque 

        real(kind=8) :: disk_mass,small_J                            ! Mass/moment of disc
        real(kind=8) :: k12                                         ! Torsional spring constant
        real(kind=8) :: c12                                         ! Torsional damping constant
      
        ! Slider Variables
        real(kind=8), DIMENSION(2) :: m_d                           ! Displacement
        real(kind=8), DIMENSION(2) :: m12_d  
        real(kind=8), DIMENSION(2) :: m1_d                               
        
        real(kind=8), DIMENSION(2) :: m_v                           ! Velocities
        real(kind=8), DIMENSION(2) :: m12_v
        real(kind=8), DIMENSION(2) :: m1_v
        
        real(kind=8), DIMENSION(2) :: m_a                           ! Acceleration        
        real(kind=8), DIMENSION(2) :: m1_a
                 
        real(kind=8), DIMENSION(2) :: Fi
        real(kind=8), DIMENSION(2,2) :: Ki
        real(kind=8), DIMENSION(2,2) :: Ci
        
        real(kind=8), DIMENSION(2,2) :: MAT1                        ! Dummy Matrices
        real(kind=8), DIMENSION(2,2) :: INV                         ! Inverse Matrix
        real(kind=8), DIMENSION(2) :: Fk                            ! Force on buoy due to mass
        real(kind=8) :: int_forcex, int_forcey
        real(kind=8) :: dmt                                         ! Determinant

        real(kind=8) :: slider_mass                                 
        real(kind=8) :: scgx,scgy                                   ! Slider Center of Mass
        real(kind=8) :: kx12            
        real(kind=8) :: x_eqm,y_eqm,y_diff                          ! Equilibrium position                                    
        real(kind=8) :: ky12, kyadj, kscale
        real(kind=8) :: ybound                                      ! Maximum spring displacement
        real(kind=8) :: cy12

        real(kind=8) :: max_n,max_s                     ! Spring constant adapts to keep mass in buoy
        real(kind=8) :: kx_adj,ky_adj                               ! Adjusted spring constants

        ! Integral terms for the RHS of rigid body equations
        real(kind=8) :: n_xint_conv, n_xint_grav, n_xint_inert      ! For the first step the convective and gravity terms
                                                                    ! For the correction step, the inertial terms
        real(kind=8) :: n_yint_conv, n_yint_grav, n_yint_inert      ! For the first step the convective and gravity terms
                                                                    ! For the correction step the inertial term
        ! Integral terms for the RHS of rigid body equations
        real(kind=8) :: n1_xint_conv, n1_xint_grav                  ! For the first step the convective and gravity terms
        real(kind=8) :: n1_xint_inert                               ! For the correction step, the inertial terms
        real(kind=8) :: n1_yint_conv, n1_yint_grav                  ! For the first step the convective and gravity terms
        real(kind=8) :: n1_yint_inert                               ! For the correction step, the intertial terms
        
        ! Integral terms for calculating fluid forces on wedge
        
        real(kind=8) :: n_w_conv, n_w_grav                  
        real(kind=8) :: n_w_inert                     
        real(kind=8) :: n1_w_conv, n1_w_grav                
        real(kind=8) :: n1_w_inert
        real(kind=8) :: wedge_force
        real(kind=8) :: wedge_work        
        real(kind=8) :: int_work

        ! Integral terms for the RHS of rigid body equations
        real(kind=8) :: r_int_conv, r_int_grav                      ! For the first step the convective and gravity terms
        real(kind=8) :: r_int_inert                                 ! For the correction step, the inertial terms
        real(kind=8) :: r1_int_conv, r1_int_grav                    ! For the first step the convective and gravity terms
        real(kind=8) :: r1_int_inert                                ! For the correction step, the inertial terms
                                                                                                                
        real(kind=8) :: ustar_int, vstar_int                        ! n+1 convective dummy variables for integral claculations
        
        real(kind=8) :: aellipse, bellipse                          ! ellipse variables
        real(kind=8) :: xo, yo, th, xn, yn, tho
        real(kind=8) :: c1,c2        
        
        real(kind=8) :: x1,x2,x3,x4                                 ! square variables
        real(kind=8) :: y1,y2,y3,y4
        real(kind=8) :: nx1,nx2,nx3,nx4                             ! square variables
        real(kind=8) :: ny1,ny2,ny3,ny4
        real(kind=8) :: cgx,cgy,ncgx,ncgy                           ! center of mass
        real(kind=8), DIMENSION(100) :: rb_data                     ! RB data - used for graphics

        real(kind=8) :: cost,sint
        real(kind=8) :: Det,a_det,b_det,c_det,temp        
        real(kind=8) :: c_xo,c_yo,c_tho
        !
        ! ALLOCATABLE ARRAYS
        real(kind=8), ALLOCATABLE :: x(:),y(:)                       ! coordinates
        real(kind=8), ALLOCATABLE :: phi(:,:),phiLS(:,:),phiLSn(:,:)             ! interface variables
        real(kind=8), ALLOCATABLE :: Fsig(:,:)
        integer, allocatable :: prox(:,:)
        !!
        !! For the third phase...(solid)
        real(kind=8), ALLOCATABLE :: s_phiLS(:,:)                    ! interface variables
        real(kind=8), allocatable :: s_phiLSn(:,:)
        real(kind=8), ALLOCATABLE :: s_H(:,:),s_HOLD(:,:)            ! heavieside field
        !!
        real(kind=8), allocatable :: w_phiLS(:,:)
        real(kind=8), ALLOCATABLE :: w_H(:,:)                        ! heavieside field
        !!
        real(kind=8), ALLOCATABLE :: H(:,:),HOLD(:,:)                ! heavieside field
        real(kind=8), ALLOCATABLE :: u(:,:),v(:,:),&
                                     P(:,:),Pold(:,:),&
                                     uINT(:,:),vINT(:,:)             ! variables
        real(kind=8), ALLOCATABLE :: F(:,:),Avof(:,:),&
                                     Bvof(:,:),Cvof(:,:)             ! CLSVOF variables
        real(kind=8), ALLOCATABLE :: ustar(:,:)                      ! temporary u-comp velocity
        real(kind=8), ALLOCATABLE :: vstar(:,:)                      ! temporary v-comp velocity
        real(kind=8), ALLOCATABLE :: u_old(:,:)                      ! An array to hold a variable
        real(kind=8), ALLOCATABLE :: v_old(:,:)                      ! An array to hold a variable   
        real(kind=8), ALLOCATABLE :: phiLSstar(:,:)                  ! temporary level set field
        real(kind=8), ALLOCATABLE :: s_phiLSstar(:,:)                ! temporary level set field        
        real(kind=8), ALLOCATABLE :: PPERHS(:,:)                     ! Poission right-hand-side
        real(kind=8), ALLOCATABLE :: BDA(:),BDB(:)                   ! BC arrays for fishpac solver
        real(kind=8), ALLOCATABLE :: BDC(:),BDD(:)                   ! BC arrays for fishpac solver
        real(kind=8), dimension(10000) :: WORK                       ! work array for fishpac solution
        real(kind=8), ALLOCATABLE :: xwhole(:),ywhole(:)             ! whole node locations
        
        real(kind=8), ALLOCATABLE :: n_rbdelom(:,:)                  ! Rotation tensor
        real(kind=8), ALLOCATABLE :: n1_rbdelom(:,:)                 ! Rotation tensor
        real(kind=8), ALLOCATABLE :: exp_big_theta(:,:) 
        
        real(kind=8), ALLOCATABLE :: x_v(:),y_v(:),z_v(:)            ! square vairables
        real(kind=8), ALLOCATABLE :: nx_v(:),ny_v(:)
        real(kind=8), ALLOCATABLE :: pt(:),rot_mat(:,:),t(:,:)

        ! Arrays for Gradient-Augmented Level Set Method (Nave 2010)
        real(kind=8), ALLOCATABLE :: phi_x(:,:),phi_xn(:,:)    ! x-direction gradient of phi
        real(kind=8), ALLOCATABLE :: phi_y(:,:),phi_yn(:,:)    ! y-direction gradient of phi
        real(kind=8), ALLOCATABLE :: phi_xy(:,:)                     ! cross-derivative of phi
        real(kind=8), ALLOCATABLE :: u_half(:,:)                     ! x-velocity at half time step
        real(kind=8), ALLOCATABLE :: v_half(:,:)                     ! y-velocity at half time step
        real(kind=8), ALLOCATABLE :: sox1(:),sox2(:),sox3(:)         ! position updates for SO-RK3
        real(kind=8), ALLOCATABLE :: sox2grad(:) 
        real(kind=8), ALLOCATABLE :: MX1(:,:),MX2(:,:),MX3(:,:)      ! matrix updates for SO-RK3
        real(kind=8), ALLOCATABLE :: LP01(:),LP10(:),LP25(:)         ! dummy variables for SO-RK3
        real(kind=8), ALLOCATABLE :: LPM01(:,:),LPM10(:,:),LPM25(:,:)! 1st number -> sox#
                                                                     ! 2nd number -> time change
        real(kind=8), ALLOCATABLE :: MM1(:,:),MM2(:,:),HC(:)         ! matrix multiplication dummy variables
        real(kind=8), ALLOCATABLE :: node(:)                
        real(kind=8), ALLOCATABLE :: ID(:,:)                         ! Identity matrix            
        !
        ! EXTERNAL FUNCTION DECLARATION
        REAL(kind=8), EXTERNAL :: LSHEAVY,LSHEAVY2
        REAL(kind=8), EXTERNAL :: LSCURVATURE
        REAL(kind=8), EXTERNAL :: GALSCURV
        REAL(kind=8), EXTERNAL :: EXPINT, ERF
        REAL(kind=8), EXTERNAL :: NORM
        REAL(kind=8), EXTERNAL :: INTERP
        REAL(kind=8), EXTERNAL :: LAGRPOLY
        !REAL(kind=8), EXTERNAL :: HERMCUBE
        REAL(kind=8), EXTERNAL :: HERMBILINEAR
        REAL(kind=8), EXTERNAL :: volfraction
        !
        !
        !
        ! //////////////////////////// INPUTS \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\        
        filename = 'inputfile_new.txt'                                   ! The name of the input file
        
        pi = 3.1415926535897932384626433832795

        ! OPENING DATA FILE
        OPEN (UNIT=2,FILE=filename,STATUS="OLD",IOSTAT=ios)
        
        ! READ THE INPUT FILE
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        !
        ! ** GEOMETRY **
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=60) Nx                                      ! Nx corresponds to the number of interior x whole nodes
        READ (UNIT=2,FMT=60) Ny                                      ! Ny corresponds to the number of interior y whole nodes
        READ (UNIT=2,FMT=70) Lx                                      ! total x length 
        READ (UNIT=2,FMT=70) Ly                                      ! total y length
        !
        !
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        !
        !
        ! ** SOLVERS **
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=60) solvertype                              ! solver for pressure poisson equation
                                                                     ! (1=FFT,2=PCG,3=SOR)
                                                                     ! NOTE: FFT only works for single-phase flow and square grids
        READ (UNIT=2,FMT=60) itermax                                 ! maximum iterations for pressure convergence -> int
        READ (UNIT=2,FMT=70) tolerance                               ! tolerance to stop pressure iteration
        READ (UNIT=2,FMT=60) iterprint                               ! print pressure residual every iterprint interations -> int
        READ (UNIT=2,FMT=60) precond                                 ! preconditioner (PCG solver), 0=none, 1=diagonal -> int
        READ (UNIT=2,FMT=70) omega                                   ! SOR weight, omega=1 corresponds to GS
        !
        !
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        !
        !
        ! ** TIME **
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=70) time_print                              ! print every time_print seconds  -> real(8)
        READ (UNIT=2,FMT=70) time_max                                ! max time -> real(8)
        READ (UNIT=2,FMT=70) stability                               ! coefficient to reduce time step such that code is stable.
        !
        !
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        !
        !
        ! ** INTERFACE **
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=80) trackphase                              ! track the interface
        READ (UNIT=2,FMT=60) LS_type                                 ! LS Method Type
        READ (UNIT=2,FMT=60) width                                   ! GALS Narrow band width
        READ (UNIT=2,FMT=70) reintime                                ! For standard L.S. only
        READ (UNIT=2,FMT=80) immersedbc                              ! use immersedbc
        !
        !
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        READ (UNIT=2,FMT=80) slosh									 ! slosh benchmark
        READ (UNIT=2,FMT=80) dambreak                                ! dambreak benchmark
        READ (UNIT=2,FMT=80) wave                                    ! wave benchmark
        READ (UNIT=2,FMT=70) waveh                                   ! wave height        
        READ (UNIT=2,FMT=80) reef                                    ! reef for wave benchmark
        READ (UNIT=2,FMT=80) rigidbody                               ! add rigid body
        READ (UNIT=2,FMT=80) wedge                                   ! use wedge
        READ (UNIT=2,FMT=80) beach                                   ! use beach (with wedge)
        READ (UNIT=2,FMT=70) freq									 ! wedge frequency (Hz)
        READ (UNIT=2,FMT=70) wedgeT									 ! wedge time limit
        !
        !
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        !
        !
        !** PROPERTIES ** 
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=60) rb_shape                                ! shape of the rigid body
                                                                     !   1 = circle
                                                                     !   2 = ellipse
        READ (UNIT=2,FMT=60) int_type                                ! type of internal component(s)
                                                                     !   0 = none
                                                                     !   1 = disc
                                                                     !   2 = 2D trans
                                                                     !   3 = 1D trans
        READ (UNIT=2,FMT=60) RB_lock                                 ! restrict motion of RB
                                                                     !   0 = No restraints
                                                                     !   1 = Vertical motion only
        READ (UNIT=2,FMT=70) rb1_mass                                ! mass of the rigid body
        READ (UNIT=2,FMT=70) buoy_J                                  ! rotational moment of inertia of the rigid body
        READ (UNIT=2,FMT=70) rho_one                                 ! density of first fluid, H=1, phiLS>0        (liquid)
        READ (UNIT=2,FMT=70) rho_two                                 ! density of second fluid, H=0, phiLS<0        (vapor)
        READ (UNIT=2,FMT=70) mu_one                                  ! viscosity of first fluid
        READ (UNIT=2,FMT=70) mu_two                                  ! viscosity of second fluid
        READ (UNIT=2,FMT=70) sigma                                   ! surface tension        
        READ (UNIT=2,FMT=70) gx                                      ! gravity constant (positive is down)
        READ (UNIT=2,FMT=70) gy                                      ! gravity constant (positive is down)
        !
        !
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        !
        !
        !** BOUNDARY TYPES **        
        ! velocity BC types - (1) dirclet wall, (2) neumman wall (also symmetry),   
        !                     (3) free surface, (4) pressure, (5) flow, (6) contact line
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=60) BCnorth 
        READ (UNIT=2,FMT=60) BCsouth  
        READ (UNIT=2,FMT=60) BCeast
        READ (UNIT=2,FMT=60) BCwest
        !
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        !
        !** BOUNDARY CONDITIONS **
        ! velocity boundary conditions [applicable for BC types (1) and (4)]
        READ (UNIT=2,FMT=50) variables                               ! skip text line        
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=70) unorth 
        READ (UNIT=2,FMT=70) usouth 
        READ (UNIT=2,FMT=70) ueast
        READ (UNIT=2,FMT=70) uwest
        !
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        !        
        READ (UNIT=2,FMT=70) vnorth 
        READ (UNIT=2,FMT=70) vsouth 
        READ (UNIT=2,FMT=70) veast 
        READ (UNIT=2,FMT=70) vwest
        !
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        !
        ! pressure dirclet BC's
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=70) Pnorth 
        READ (UNIT=2,FMT=70) Psouth 
        READ (UNIT=2,FMT=70) Peast
        READ (UNIT=2,FMT=70) Pwest 
        !
        READ (UNIT=2,FMT=50) variables                               ! skip blank line
        !
        ! contact line boundary conditions 
        ! Ucl = contC*(angle - contA)**contM IF angle(phi) > contA
        ! Ucl = -contC*(contR - angle)**contM IF angle(phi) < contR
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=50) variables                               ! skip text line
        READ (UNIT=2,FMT=70) contC                                   ! constant        
        READ (UNIT=2,FMT=70) contA                                   ! static advancing contact angle
        READ (UNIT=2,FMT=70) contR                                   ! static receding contact angle
        READ (UNIT=2,FMT=70) contM                                   ! mobility exponent
        !        
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        !
        !** INITIAL CONDITIONS FOR FREE SURFACE**        
        ! phiLS = a*sqrt( (x-b)^2 +(y-c)^2 ) + d + ex + fy
        READ (UNIT=2,FMT=50) variables  
        READ (UNIT=2,FMT=50) variables 
        READ (UNIT=2,FMT=70) a_phiLS
        READ (UNIT=2,FMT=70) b_phiLS
        READ (UNIT=2,FMT=70) c_phiLS
        READ (UNIT=2,FMT=70) d_phiLS
        READ (UNIT=2,FMT=70) e_phiLS
        READ (UNIT=2,FMT=70) f_phiLS
        READ (UNIT=2,FMT=70) v_phiLS
        !
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        !
        !** INITIAL CONDITIONS FOR CIRCULAR RIGID BODY**        
        ! s_phiLS = a*sqrt( (x-b)^2 +(y-c)^2 ) + d
        READ (UNIT=2,FMT=50) variables  
        READ (UNIT=2,FMT=50) variables 
        READ (UNIT=2,FMT=70) a_rbLS
        READ (UNIT=2,FMT=70) b_rbLS
        READ (UNIT=2,FMT=70) c_rbLS
        READ (UNIT=2,FMT=70) d_rbLS
        !
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        !
        !** INITIAL CONDITIONS FOR ELLIPTICAL RIGID BODY**        
        ! s_phiLS = 
        READ (UNIT=2,FMT=50) variables  
        READ (UNIT=2,FMT=50) variables 
        READ (UNIT=2,FMT=70) aellipse
        READ (UNIT=2,FMT=70) bellipse
        READ (UNIT=2,FMT=70) xo
        READ (UNIT=2,FMT=70) yo
        READ (UNIT=2,FMT=70) tho        
        !
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        !
        !** INITIAL CONDITIONS FOR SQUARE RIGID BODY**        
        ! s_phiLS = 
        READ (UNIT=2,FMT=50) variables  
        READ (UNIT=2,FMT=50) variables 
        READ (UNIT=2,FMT=70) c_xo
        READ (UNIT=2,FMT=70) c_yo
        READ (UNIT=2,FMT=70) x1
        READ (UNIT=2,FMT=70) y1
        READ (UNIT=2,FMT=70) x2
        READ (UNIT=2,FMT=70) y2
        READ (UNIT=2,FMT=70) x3
        READ (UNIT=2,FMT=70) y3
        READ (UNIT=2,FMT=70) x4
        READ (UNIT=2,FMT=70) y4
        READ (UNIT=2,FMT=70) c_tho
        READ (UNIT=2,FMT=70) cgx
        READ (UNIT=2,FMT=70) cgy
        !        
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        !
        !** TORSIONAL DISC PARAMETERS (int_type = 1) ***        
        READ (UNIT=2,FMT=50) variables   
        READ (UNIT=2,FMT=70) disk_mass                                ! mass of the torsional disc
        READ (UNIT=2,FMT=70) small_J                                 ! rotational moment of inertia of the disc
        READ (UNIT=2,FMT=70) k12                                     ! torsional spring constant
        READ (UNIT=2,FMT=70) c12                                     ! torsional damping constant
        !        
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        !
        !** SLIDER PARAMETERS (int_type = 2 OR 3) ***
        READ (UNIT=2,FMT=50) variables   
        READ (UNIT=2,FMT=70) slider_mass                             ! mass of the slider
        READ (UNIT=2,FMT=70) kx12                                    ! slider spring constant (x-direction)
        READ (UNIT=2,FMT=70) y_diff                                  ! slider damping constant (x-direction)
        READ (UNIT=2,FMT=70) ky12                                    ! slider spring constant (y-direction)
        READ (UNIT=2,FMT=70) kscale                                  ! slider spring ramps up near boundary
        READ (UNIT=2,FMT=70) cy12                                    ! slider damping constant (y-direction)
        !
        READ (UNIT=2,FMT=50) variables  ! skip blank line
        !
        CLOSE (UNIT=2)                                               ! close the input file
        !
        ! ///////////////////// INITIALIZE VARIABLES \\\\\\\\\\\\\\\\\\\\\\
        ! allocating the size of each arrays        
        ALLOCATE (x(Nx+1))
        ALLOCATE (y(Ny+1))
        ALLOCATE (xwhole(Nx+2))
        ALLOCATE (ywhole(Ny+2))
        !        
        ALLOCATE (phi(Nx+2,Ny+2))
        ALLOCATE (phiLS(Nx+2,Ny+2))
        ALLOCATE (phiLSn(Nx+2,Ny+2))
        ALLOCATE (prox(Nx+2,Ny+2))
        AllOCATE (H(Nx+2,Ny+2))
        AllOCATE (HOLD(Nx+2,Ny+2))
        !
        ALLOCATE (phi_x(Nx+2,Ny+2))                  ! G.A. LS Method
        ALLOCATE (phi_xn(Nx+2,Ny+2))
        ALLOCATE (phi_y(Nx+2,Ny+2))
        ALLOCATE (phi_yn(Nx+2,Ny+2))
        ALLOCATE (phi_xy(Nx+2,Ny+2))
 
        ALLOCATE (Fsig(Nx+2,Ny+2))         
        AllOCATE (u_half(Nx+1,Ny+2))
        AllOCATE (v_half(Nx+2,Ny+1))        
        ALLOCATE (sox1(2))
        ALLOCATE (sox2(2))
        ALLOCATE (sox2grad(2))
        ALLOCATE (sox3(2))
        ALLOCATE (MX1(2,2))
        ALLOCATE (MX2(2,2))
        ALLOCATE (MX3(2,2))
        ALLOCATE (LP01(2))
        ALLOCATE (LP10(2))
        ALLOCATE (LP25(2))
        ALLOCATE (LPM01(2,2))
        ALLOCATE (LPM10(2,2))
        ALLOCATE (LPM25(2,2))
        ALLOCATE (MM1(2,2))
        ALLOCATE (MM2(2,2))     
        ALLOCATE (HC(2))
        ALLOCATE (node(2))
        ALLOCATE (ID(2,2))
        !
        ALLOCATE (s_phiLS(Nx+2,Ny+2))
        ALLOCATE (s_phiLSn(Nx+2,Ny+2))
        AllOCATE (s_H(Nx+2,Ny+2))
        AllOCATE (s_HOLD(Nx+2,Ny+2))
        !
        ALLOCATE (w_phiLS(Nx+2,Ny+2))
        AllOCATE (w_H(Nx+2,Ny+2))
        !
        ALLOCATE (F(Nx+2,Ny+2))
        ALLOCATE (Avof(Nx+2,Ny+2))
        ALLOCATE (Bvof(Nx+2,Ny+2))
        ALLOCATE (Cvof(Nx+2,Ny+2))
        !
        AllOCATE (u(Nx+1,Ny+2))
        AllOCATE (v(Nx+2,Ny+1))
        AllOCATE (u_old(Nx+1,Ny+2))
        AllOCATE (v_old(Nx+2,Ny+1))
        
        AllOCATE (P(Nx+2,Ny+2))
        ALLOCATE (Pold(Nx+2,Ny+2))
        AllOCATE (uINT(Nx+1,Ny+2))
        AllOCATE (vINT(Nx+2,Ny+1))
        
        !For the square
        AllOCATE (x_v(3))
        AllOCATE (y_v(3))
        AllOCATE (nx_v(3))
        AllOCATE (ny_v(3))
        AllOCATE (z_v(3))
        AllOCATE (pt(2))
        ALLOCATE (rot_mat(2,2))
        AllOCATE (t(2,3))
        !
        AllOCATE (ustar(Nx+1,Ny+2))
        AllOCATE (vstar(Nx+2,Ny+1)) 
        AllOCATE (phiLSstar(Nx+2,Ny+2))
        AllOCATE (s_phiLSstar(Nx+2,Ny+2))
        
        ALLOCATE (n_rbdelom(2,2))
        ALLOCATE (n1_rbdelom(2,2))
        ALLOCATE (exp_big_theta(2,2))
        !        
        ! PPE RHS
        AllOCATE (PPERHS(Nx,Ny))
        !
        ! BC arrays for fishpac solver
        AllOCATE (BDA(Ny))
        AllOCATE (BDB(Ny))   
        AllOCATE (BDC(Nx))
        AllOCATE (BDD(Nx))
        !
        !
        ! time and counter initialization
        time = 0.0                                                   ! time is zero at start
        file_count = 0                                               ! print file counter 
        time_count = 0.0                                             ! counter for when to print
        dumpcount = 0                                                ! dump file counter
        count = 0
        !
        ! Identity matrix
        ID(1,1) = 1.0
        ID(1,2) = 0.0
        ID(2,1) = 0.0
        ID(2,2) = 1.0
        !
        ! RK4 coefficients
        Coef(1) = 1.0/4.0
        Coef(2) = 1.0/3.0
        Coef(3) = 1.0/2.0
        Coef(4) = 1.0
        
        ! CLSVOF coefficients
        xdirgo = 1
        ydirgo = 2

        ! Initialize rigid body variables
        
        n_rb_dx = 0.0
        n_rb_dy = 0.0
        n1_rb_dx = 0.0
        n1_rb_dy = 0.0
        n_w_dy = 0.0
        n1_w_dy = 0.0
        n12_w_vy = 0.0
        n1_w_vy = 0.0
        
        n_rb_vx = 0.0
        n_rb_vy = 0.0
        n12_rb_vx = 0.0
        n12_rb_vy = 0.0
        n1_rb_vx = 0.0
        n1_rb_vy = 0.0

        n_rb_ax = 0.0
        n_rb_ay = 0.0
        n1_rb_ax = 0.0
        n1_rb_ay = 0.0
        
        ! Rotation variables
        
        c1 = 0.0
        c2 = 0.0
        
        rx = 0.0
        ry = 0.0
        
        n_rbdelom(1,1) = 1.0
        n_rbdelom(1,2) = 0.0
        n_rbdelom(2,1) = 0.0
        n_rbdelom(2,2) = 1.0
        
        n1_rbdelom(1,1) = 1.0
        n1_rbdelom(1,2) = 0.0
        n1_rbdelom(2,1) = 0.0
        n1_rbdelom(2,2) = 1.0
        
        ! The exponential map
        exp_big_theta(1,1) = 0.0
        exp_big_theta(1,2) = 0.0
        exp_big_theta(2,1) = 0.0
        exp_big_theta(2,2) = 0.0
        
        big_theta = 0.0
        
        n_rb_vom = 0.0
        n12_rb_vom = 0.0
        n1_rb_vom = 0.0
        
        n_rb_aom = 0.0
        n1_rb_aom = 0.0

        th_change = 0.0
        th_domain = 0.0
        th_old = 0.0
        th_counter = 0.0
        
        Fsig = 0.0
        
        ! //////// STORING RIGID BODY INFORMATION FOR GRAPHICS POSTPROCESSING \\\\\\\\

        rb_data = 0.0
        rb_data(1) =    rb_shape
        rb_data(2) =    int_type
        rb_data(3) =    RB_lock
        rb_data(4) =    rb1_mass
        rb_data(5) =    buoy_J
        rb_data(6) =    a_rbLS
        rb_data(7) =    b_rbLS
        rb_data(8) =    c_rbLS
        rb_data(9) =    d_rbLS
        rb_data(10) =   aellipse
        rb_data(11) =   bellipse
        rb_data(12) =   xo
        rb_data(13) =   yo
        rb_data(14) =   tho
        rb_data(15) =   c_xo
        rb_data(16) =   c_yo
        rb_data(17) =   x1
        rb_data(18) =   y1
        rb_data(19) =   x2
        rb_data(20) =   y2
        rb_data(21) =   x3
        rb_data(22) =   y3
        rb_data(23) =   x4
        rb_data(24) =   y4
        rb_data(26) =   c_tho
        rb_data(27) =   cgx
        rb_data(28) =   cgy
        rb_data(29) =   disk_mass
        rb_data(30) =   small_J
        rb_data(31) =   k12
        rb_data(32) =   c12
        rb_data(33) =   slider_mass
        rb_data(34) =   kx12
        rb_data(35) =   y_diff
        rb_data(36) =   ky12
        rb_data(37) =   cy12
        IF (wedge .EQV. .TRUE.) THEN
            rb_data(38) =   1
        ELSE
            rb_data(38) = 0
        ENDIF
        rb_data(39) =   n_w_dy

        ! Max spring displacement

        ybound = -sin(c_tho)*(cgx - (x3+x4)/2.0) +cos(c_tho)*(cgy - (y3+y4)/2.0) 

        ! RB/Internals components
           
        dtinternal = 10.0

        ncgx = cgx
        ncgy = cgy

        IF (int_type == 0) THEN
             
            rb_mass = rb1_mass

        ELSEIF (int_type == 1) THEN
             
            rb_mass = rb1_mass + disk_mass
             
                IF ((k12 /= 0.0) .and. (c12 /= 0.0)) THEN
                    dtinternal = min( sqrt(min(buoy_J,small_J)/k12), sqrt(min(buoy_J,small_J)/c12) )
                ENDIF
           
        ELSEIF (int_type == 2) THEN

            rb_mass = rb1_mass
            x_eqm = cgx - sin(c_tho)*y_diff
            y_eqm = cgy + cos(c_tho)*y_diff 
            scgx = x_eqm - 0.0*sin(c_tho)*cos(c_tho)*slider_mass*(-gy)/ky12                             
            scgy = y_eqm + 0.0*cos(c_tho)*cos(c_tho)*slider_mass*(-gy)/ky12
             
        end if
        
        ! ---- FLUID VELOCITY ----
        u = 0.0                
        v = 0.0
        u_old = 0.0
        v_old = 0.0
        u_half = 0.0
        v_half = 0.0        
                         
        total_energy = 0.0
        lin_momx = 0.0
        lin_momy = 0.0

        ! Disc
        small_theta_n = 0.0
        small_theta_n1 = 0.0
 
        n_rb2_vom = 0.0
        n12_rb2_vom = 0.0
        n1_rb2_vom = 0.0
        
        n_rb2_aom = 0.0
        n1_rb2_aom = 0.0
        delta_theta = 0.0
        delta_vom = 0.0
        int_torque = 0.0
        int_work = 0.0

        ! Slider
        m_d(1) = 0.0
        m_d(2) = 0.0*cos(c_tho)*slider_mass*(-gy)/ky12          ! Displacement    
        m12_d = 0.0
        m1_d = 0.0
        
        m_v = 0.0                            ! Velocities
        m12_v = 0.0
        m1_v = 0.0
        
        m_a = 0.0                            ! Acceleration  
        m1_a = 0.0
           
        Fi = 0.0
        Ki = 0.0
        Ci = 0.0

        MAT1 = 0.0
        INV = 0.0
        Fk = 0.0
        int_forcex = 0.0
        int_forcey = 0.0
           

        ! Initialize intergral variables two different directions
        n_xint_conv  = 0.0
        n_xint_grav  = 0.0
        n_xint_inert = 0.0
        
        n_yint_conv  = 0.0
        n_yint_grav  = 0.0
        n_yint_inert = 0.0
        
        ! Initialize intergral variables two different directions
        n1_xint_conv  = 0.0
        n1_xint_grav  = 0.0
        n1_xint_inert = 0.0
        
        n1_yint_conv  = 0.0
        n1_yint_grav  = 0.0
        n1_yint_inert = 0.0

        ! Initialize integral variables for wedge (in one direction)
        
        n_w_conv  = 0.0
        n_w_grav  = 0.0
        n_w_inert = 0.0
        
        n1_w_conv  = 0.0
        n1_w_grav  = 0.0
        n1_w_inert = 0.0

        wedge_force = 0.0
        wedge_work = 0.0

        ! Initialize rotation integral variables
        r_int_conv  = 0.0
        r_int_grav  = 0.0
        r_int_inert = 0.0
        
        r1_int_conv  = 0.0
        r1_int_grav  = 0.0
        r1_int_inert = 0.0
        
        ustar_int = 0.0
        vstar_int = 0.0
        
        xn = 0.0
        yn = 0.0
        th = 0.0
        xc = 0.0
        yc = 0.0
    
        ! 1 = neumman, 2 = dirclet, 3 = extrapolation
        BCphisouth = 1        ! 3
        BCphinorth = 1        ! 3
        BCphiwest = 1        ! 1
        BCphieast = 1        ! 1
        
        ! ////////////////////// COORDINATE GENERATION \\\\\\\\\\\\\\\\\\\\\\
        hx = Lx/(real(Nx))
        hy = Ly/(real(Ny))
        x(1) = 0.0
        y(1) = 0.0
        x(Nx+1) = Lx
        y(Ny+1) = Ly
        
        ! Coefficient arrays for GALS routines
	    GA(1,1) = 1
	    GA(1,2) = 2
	    GA(1,3) = 1
	    GA(1,4) = 3
	    GA(1,5) = 1
	    GA(1,6) = 2
	    
	    GA(2,1) = 1
	    GA(2,2) = 1
	    GA(2,3) = 2
	    GA(2,4) = 1
	    GA(2,5) = 3
	    GA(2,6) = 2
	    
	    GA(3,1) = 1
	    GA(3,2) = 1
	    GA(3,3) = -1
	    GA(3,4) = 1
	    GA(3,5) = 1
	    GA(3,6) = -1
	    
	    GA(4,1) = 1
	    GA(4,2) = -1
	    GA(4,3) = 1
	    GA(4,4) = 1
	    GA(4,5) = 1
	    GA(4,6) = -1
	    
	    GA(5,1) = 1
	    GA(5,2) = -1
	    GA(5,3) = -1
	    GA(5,4) = 1
	    GA(5,5) = 1
	    GA(5,6) = 1
	    
    	scal(1) = 1.0
    	scal(2) = 1.0/hx
    	scal(3) = 1.0/hy
    	scal(4) = 1.0/(hx**2.0)
    	scal(5) = 1.0/(hy**2.0)
    	scal(6) = 1.0/(hx*hy)
        
        rb_data(40) = hx
        rb_data(41) = hy

        DO i=2,Nx
            x(i) = hx*real(i-1)
        ENDDO

        DO j=2,Ny
            y(j) = hy*real(j-1)
        ENDDO

        ! whole node locations
        DO i=1,Nx+1
            xwhole(i) = x(i)-0.5*hx
        ENDDO

        xwhole(Nx+2)=x(Nx+1)+0.5*hx

        DO j=1,Ny+1
            ywhole(j) = y(j)-0.5*hy
        ENDDO

        ywhole(Ny+2)=y(Ny+1)+0.5*hy

        
        ! ///////////////// INITIAL CONDITIONS \\\\\\\\\\\\\\\\\\\

        ! First Time Step
        deltaT = min(1.0/( max(mu_one/rho_one,mu_two/rho_two)*&
                         (2.0/(hx**2.0) + 2.0/(hy**2.0)) ),dtinternal)

		print *, "deltaT = ",deltaT

        ! ---- IMMERSED BOUNDARY ----        
        
        s_phiLS = 0.0                                                ! Rigid Body
        phiLS = 0.0                                                  ! Free surface
        phi = 0.0                                                    ! Immersed Boundary
 		prox = 1.0
         
        IF (immersedbc .EQV. .TRUE.) THEN
            DO i=1,Nx+2
                DO j=1,Ny+2
                    ! phi>0 is solid, phi<0 is fluid
                    phi(i,j) = -sqrt( (xwhole(i)-0.5)**2.0 +&
                                      (ywhole(j)-0.5)**2.0 ) + 0.2 
                ENDDO
            ENDDO
        ENDIF
                
        ! ---- PRESSURE        ----        
        P = 0.0                                                      ! initializing array to zero

        ! ---------------------------------- !
        ! ---- LEVEL SET INITIALIZATION ---- !
        ! ---------------------------------- !
        
        IF (trackphase .EQV. .TRUE.) THEN
        
            DO i=1,Nx+2
                DO j=1,Ny+2

                    ! DAMBREAK LEVEL SET
                    ! Use domain dimensions 5 x 1.25
                    IF (dambreak .eqv. .true.) THEN

                        a_phiLS = 0.05715

                        IF (xwhole(i) >= ywhole(j)) THEN
                            phiLS(i,j) = a_phiLS - 1.0*xwhole(i)
                            phi_x(i,j) = -1.0
                            phi_y(i,j) = 0.0
                        ELSE
                            phiLS(i,j) = a_phiLS - 1.0*ywhole(j)
                            phi_x(i,j) = 0.0
                            phi_y(i,j) = -1.0
                        ENDIF

                            phi_xy(i,j) = 0.0

                    ! WAVE LEVEL SET
                    ! Use domain dimensions 20*d_phils x (>=2*d_phiLS)
                    ELSEIF (wave .eqv. .true.) THEN
                        
                        ! FIRST-ORDER WAVE PROFILE
                        !d_phiLS = Ly/2.0
                        !phiLS(i,j) = waveh/(cosh(sqrt(3.0*waveh)/2.0*xwhole(i))**2.0) + d_phiLS - ywhole(j)
						!
						!phi_x(i,j) = (phiLS(i,j)-d_phiLS+ywhole(j))*(-sqrt(3.0*waveh)/2.0)* &
						!			 tanh(sqrt(3.0*waveh)/2.0*xwhole(i))
						!phi_y(i,j) = -1.0
                        
                        ! SECOND-ORDER WAVE PROFILE
                        !phiLS(i,j) = 1.0 + waveh/(cosh(d_phiLS*(xwhole(i)-10.0))**2.0) &
                        !             - 0.75*(waveh**2.0)/(cosh(d_phiLS*(xwhole(i)-10.0))**2.0) &
                        !             * (1.0 - 1.0/(cosh(d_phiLS*(xwhole(i)-10.0))**2.0)) - ywhole(j)
                        !phi_x(i,j) = -2.0*d_phiLS*(waveh-0.75*waveh*waveh) &
                        !             *tanh(d_phiLS*(xwhole(i)-10.0))/(cosh(d_phiLS*(xwhole(i)-10.0))**2.0)
                        !phi_y(i,j) = -1.0

                        ! THIRD-ORDER WAVE PROFILE
                        d_phiLS = sqrt(3.0*waveh/4.0)*(1.0-5.0*waveh/8.0+71.0*waveh*waveh/128.0)
                        sw = 1/cosh(d_phiLS*(xwhole(i)-10.0))
                        tw = tanh(d_phiLS*(xwhole(i)-10.0))
                        phiLS(i,j) = 1.0 + waveh*sw*sw - 0.75*waveh*waveh*sw*sw*tw*tw &
                                         + waveh**3.0*(5.0/8.0*sw*sw*tw*tw - 101.0/80.0*(sw**4.0)*tw*tw) - ywhole(j)
                        phi_x(i,j) = -1/40.0*d_phiLS*waveh*tw*sw*sw*(50.0*waveh**2.0*tw*tw + 101.0*waveh**2.0*sw**4.0 &
                                     - 50.0*waveh**2.0*sw*sw - 202.0*waveh*waveh*tw*tw*sw*sw - 60.0*waveh*tw*tw &
                                     + 60.0*waveh*sw*sw + 80.0)
                        phi_y(i,j) = -1.0

                        !if (i < Nx+2 .AND. phiLS(i,j) > -0.20) then
                             
                             ! THIRD-ORDER VELOCITY (worse than 2nd-order?)
                             !u(i,j) = sqrt(gy)*(waveh*sw*sw &
                             !         - waveh**2.0*(-0.25*sw*sw + sw**4.0 + ywhole(j)**2.0*(1.5*sw*sw-2.25*(sw**4.0))) &
                             !         - waveh**3.0*(19.0/40.0*sw*sw + 0.2*sw**4.0 - 1.2*sw**6.0 &
                             !                    + ywhole(j)**2.0*(-1.5*sw*sw-15.0/4.0*sw**4.0+7.5*sw**6.0) &
                             !                    + ywhole(j)**4.0*(-3.0/8.0*sw*sw + 45.0/16.0*sw**4.0 - 45.0/16.0*sw**6.0)))
                             ! SECOND-ORDER VELOCITY
                             u(i,j) = sqrt(gy)*(waveh*(1.0 - 1.25*waveh - 1.5*waveh*(2.0*(ywhole(j)-1.0)+(ywhole(j)-1.0)**2.0)) &
                                      * 1.0/(cosh(d_phiLS*(xwhole(i)-10.0))**2.0) &
                                      + waveh*waveh*(1.25+2.25*(2.0*(ywhole(j)-1.0)+(ywhole(j)-1.0)**2.0)) &
                                      * 1.0/(cosh(d_phiLS*(xwhole(i)-10.0))**4.0)) 
                                     
                        !endif
                        !if (j < Ny+2 .AND. phiLS(i,j) > -0.20) then
                             ! THIRD-ORDER VELOCITY
                             !v(i,j) = - sqrt(3.0*gy*waveh)*ywhole(j)*tw*(-waveh*sw*sw &
                             !         + waveh**2.0*(3.0/8.0*sw*sw + 2.0*sw**4.0 + ywhole(j)**2.0*(0.5*sw*sw - 1.5*sw**4.0)) &
                             !         + waveh**3.0*(49.0/640.0*sw*sw - 17.0/20.0*sw**4.0 - 18.0/5.0*sw**6.0 &
                             !                       + ywhole(j)**2.0*(-13.0/16.0*sw*sw - 25.0/16.0*sw**4.0 + 7.5*sw**6.0) &
                             !                       + ywhole(j)**4.0*(-3.0/40.0*sw*sw + 9.0/8.0*sw**4.0 - 27.0/16.0*sw**6.0)))
                             ! SECOND-ORDER VELOCITY
                             v(i,j) = sqrt(3.0*gy)*(waveh**(3.0/2.0))*ywhole(j) &
                                    * tanh(d_phiLS*(xwhole(i)-10.0))/(cosh(d_phiLS*(xwhole(i)-10.0))**2.0) &
                                    * (1.0-7.0/8.0*waveh-0.5*waveh*(2.0*(ywhole(j)-1.0)+(ywhole(j)-1.0)**2.0) &
                                    - 0.5*waveh*(1.0-6.0*(ywhole(j)-1.0)-3.0*(ywhole(j)-1.0)**2.0) &
                                    * 1.0/(cosh(d_phiLS*(xwhole(i)-10.0))**2.0))

                        !endif

                    ! STANDARD LEVEL SET
                    ELSE

                        phiLS(i,j) = a_phiLS*sqrt( (xwhole(i)-b_phiLS)**2.0+&
                                                   (ywhole(j)-c_phiLS)**2.0 )+&
                                                    d_phiLS+&
                                                    e_phiLS*xwhole(i)+&
                                                    f_phiLS*ywhole(j)
                                
                        phi_x(i,j) = e_phiLS  ! a_phiLS*(1.0 / sqrt((ywhole(i)-b_phiLS)**2.0+&
                                              !                     (ywhole(j)-c_phiLS)**2.0 ))*&
                                              !                     (xwhole(i)-b_phiLS)+&
                                              !                      e_phiLS
                        phi_y(i,j) = f_phiLS  ! a_phiLS*(1.0 / sqrt((xwhole(i)-b_phiLS)**2.0+&
                                              !                     (ywhole(j)-c_phiLS)**2.0 ))*&
                                              !                     (ywhole(j)-c_phiLS)+&
                                              !                      f_phiLS
                        phi_xy(i,j) = 0       ! - a_phiLS/2.0 * (((xwhole(i)-b_phiLS)**2.0+&
                                              !         (ywhole(j)-c_phiLS)**2.0)**(-1.5))*&
                                              !         (xwhole(i)-b_phiLS)*&
                                              !         (ywhole(j)-c_phiLS) + &
                                              
                        !phiLS(i,j) =   exp(-(xwhole(i)-0.5)**2.0-(ywhole(j)-0.75)**2.0) - exp(-0.15**2.0)
                        !phi_x(i,j) = - 2.0*(xwhole(i) - 0.50)*exp(-(xwhole(i)-0.5)**2.0-(ywhole(j)-0.75)**2.0)
                        !phi_y(i,j) = - 2.0*(ywhole(j) - 0.75)*exp(-(xwhole(i)-0.5)**2.0-(ywhole(j)-0.75)**2.0)
                        !phi_xy(i,j) = 0

                        !drop = .TRUE.
                        drop = .FALSE.
                        if (drop) then
        
                            !phiLS(i,j) =   exp(-(xwhole(i)-0.004)**2.0-(ywhole(j)-0.006)**2.0) - exp(-0.001**2.0)
                            !phi_x(i,j) = - 2.0*(xwhole(i) - 0.004)*exp(-(xwhole(i)-0.004)**2.0-(ywhole(j)-0.006)**2.0)
                            !phi_y(i,j) = - 2.0*(ywhole(j) - 0.006)*exp(-(xwhole(i)-0.004)**2.0-(ywhole(j)-0.006)**2.0)

                            phiLS(i,j) = -sqrt((xwhole(i)-0.003125)**2.0+(ywhole(j)-0.01)**2.0) + 0.00125
                            phi_x(i,j) = (xwhole(i)-0.003125)/(phiLS(i,j)-0.00125)
                            phi_y(i,j) = (ywhole(j)-0.01)/(phiLS(i,j)-0.00125)

                            if (ywhole(j) < 0.004) then

                                phiLS(i,j) = -ywhole(j)+0.002
                                phi_x(i,j) = 0.0
                                phi_y(i,j) = -1.0
                                phi_xy(i,j) = 0.0
                            
                            endif

                        endif 
                                                                                    
                    ENDIF
                    
                ENDDO
            ENDDO
        ENDIF
        
    	! ------------------------------ !
        ! ----  SOME GALS ROUTINES  ---- !
        ! ------------------------------ !
        
        if (LS_type == 3) then

    		CALL CROSSDERV(Nx,Ny,hx,hy,phiLS,phi_x,phi_y,phi_xy,prox,width)   
			CALL REBUILD(Nx,Ny,hx,hy,xwhole,ywhole,GA,scal,phiLS,phi_x,phi_y,phi_xy,prox,width)
			
			! Volume routines
			volume = 0.0
			CALL VOLUMESOLVE(Nx,Ny,hx,hy,xwhole,ywhole,GA,scal,phiLS,phi_x,phi_y,phi_xy,prox,volume)
			! Initial volume
			volume0 = volume
			
		ENDIF
		
        ! ------------------------------ !
        ! ---- WEDGE INITIALIZATION ---- !
        ! ------------------------------ !
        
        IF (wedge .EQV. .TRUE.) THEN
            DO i=1,Nx+2
                DO j=1,Ny+2        
        
                    IF (ywhole(j) > 3.82402*xwhole(i)-0.340938) THEN
                        w_phiLS(i,j) = -(-xwhole(i) + 0.114)
                    else
                        w_phiLS(i,j) = -(-0.46345*ywhole(j) + 0.88612*xwhole(i) - 0.05699)
                    ENDIF
                        
                ENDDO
            ENDDO
        ENDIF

        ! ----------------------------------- !
        ! ---- RIGID BODY INITIALIZATION ---- !
        ! ----------------------------------- !
        
        IF (rigidbody .EQV. .TRUE.) THEN

            ! CIRCLE        
            IF (rb_shape == 1) THEN
                    
                th = 0.0
                th_old = th
                xc = b_rbLS
                yc = c_rbLS
                      
                rb_data(60) = xc            ! geometric center & angle
                rb_data(61) = yc
                rb_data(62) = tho  

                DO i=1,Nx+2
                    DO j=1,Ny+2

                        s_phiLS(i,j) = a_rbLS*sqrt( (xwhole(i)-b_rbLS)**2.0+&
                                                    (ywhole(j)-c_rbLS)**2.0 )+&
                                                     d_rbLS
                                                        
                    ENDDO
                ENDDO
                
            ! ELLIPSE
            ELSEIF (rb_shape == 2) THEN
        
                xn = xo
                yn = yo
                xc = xn
                yc = yn
                th = pi/tho
                th_old = th
                
                rb_data(60) = xc            ! geometric center & angle
                rb_data(61) = yc
                rb_data(62) = tho  

                DO i=1,Nx+2
                    DO j=1,Ny+2
                                
                        s_phiLS(i,j) = -100*&
                                        ( bellipse**2*( xwhole(i)*cos(th)+&
                                                        ywhole(j)*sin(th)-&
                                                        (xo*cos(th) + yo*sin(th)))**2 + &
                                          aellipse**2*(-xwhole(i)*sin(th) + ywhole(j)*cos(th) - &
                                                        (-xo*sin(th) + yo*cos(th)))**2 - &
                                          aellipse**2*bellipse**2)

                    ENDDO
                ENDDO
        
            ! SQUARE
            ELSEIF (rb_shape == 3) THEN
        
                xc = c_xo  ! Geometric Center
                yc = c_yo
                
                !th = pi/c_tho
                th = c_tho
                th_old = th
                cost = cos(c_tho)
                sint = sin(c_tho)
                rot_mat(1,1) = cost
                rot_mat(1,2) = -sint
                rot_mat(2,1) = sint
                rot_mat(2,2) = cost
                
                rb_data(60) = xc            ! geometric center & angle
                rb_data(61) = yc
                rb_data(62) = th  
                        
                ! x1-4,y-4 are ordered counterclockwise
                
                DO r=1,4
                
                    IF (r == 1) THEN
                        
                        ! Region 1
                        x_v(1) = xc;
                        x_v(2) = x1;
                        x_v(3) = x2;
                                
                        y_v(1) = yc;
                        y_v(2) = y1;
                        y_v(3) = y2;
                                
                    ELSEIF (r == 2) THEN
                        
                        ! Region 2
                        x_v(1) = xc;
                        x_v(2) = x2;
                        x_v(3) = x3;
                                
                        y_v(1) = yc;
                        y_v(2) = y2;
                        y_v(3) = y3;
                                
                    ELSEIF (r == 3) THEN
                        
                        ! Region 3
                        x_v(1) = xc;
                        x_v(2) = x3;
                        x_v(3) = x4;
                                
                        y_v(1) = yc;
                        y_v(2) = y3;
                        y_v(3) = y4;
                                
                    ELSEIF (r == 4) THEN
                        
                        ! Region 4
                        x_v(1) = xc;
                        x_v(2) = x4;
                        x_v(3) = x1;
                                
                        y_v(1) = yc;
                        y_v(2) = y4;
                        y_v(3) = y1;
                                
                    ENDIF

                    nx_v(1) = rot_mat(1,1)*(xc-cgx) + rot_mat(1,2)*(yc-cgy)
                    ny_v(1) = rot_mat(2,1)*(xc-cgx) + rot_mat(2,2)*(yc-cgy)                        
                    nx_v(2) = rot_mat(1,1)*(x_v(2)-cgx) + rot_mat(1,2)*(y_v(2)-cgy)
                    ny_v(2) = rot_mat(2,1)*(x_v(2)-cgx) + rot_mat(2,2)*(y_v(2)-cgy)
                    nx_v(3) = rot_mat(1,1)*(x_v(3)-cgx) + rot_mat(1,2)*(y_v(3)-cgy)
                    ny_v(3) = rot_mat(2,1)*(x_v(3)-cgx) + rot_mat(2,2)*(y_v(3)-cgy)

                    x_v(1) = nx_v(1) + cgx
                    y_v(1) = ny_v(1) + cgy                        
                    x_v(2) = nx_v(2) + cgx
                    y_v(2) = ny_v(2) + cgy
                    x_v(3) = nx_v(3) + cgx
                    y_v(3) = ny_v(3) + cgy
                                    
                    z_v(1) = 1.0;
                    z_v(2) = 0.0;
                    z_v(3) = 0.0;
                                        
                    Det = x_v(1)*(y_v(2)*z_v(3)-y_v(3)*z_v(2)) - &
                          y_v(1)*(x_v(2)*z_v(3)-x_v(3)*z_v(2)) + &
                          z_v(1)*(x_v(2)*y_v(3)-x_v(3)*y_v(2))
                                          
                    a_det = y_v(2)*z_v(3)-y_v(3)*z_v(2) - &
                            y_v(1)*(z_v(3)-z_v(2)) + &
                            z_v(1)*(y_v(3)-y_v(2))
                                          
                    b_det = x_v(1)*(z_v(3)-z_v(2)) - &
                            x_v(2)*z_v(3)-x_v(3)*z_v(2) + &
                            z_v(1)*(x_v(2)-x_v(3))
                                          
                    c_det = x_v(1)*(y_v(2)-y_v(3)) - &
                            y_v(1)*(x_v(2)-x_v(3)) + &
                            x_v(2)*y_v(3)-x_v(3)*y_v(2)
                                                
                    a_rbLS = a_det/Det
                    b_rbLS = b_det/Det
                    c_rbLS = c_det/Det
                                
                    DO i=1,Nx+2
                        DO j=1,Ny+2
                        
                            ! Check on the location of the grid point.
                            ! Is it in the correct region?
                                
                            pt(1) = xwhole(i)
                            pt(2) = ywhole(j)
                                
                            t(1,1) = x_v(1);
                            t(2,1) = y_v(1);
                            t(1,2) = x_v(2);
                            t(2,2) = y_v(2);
                            t(1,3) = x_v(3);
                            t(2,3) = y_v(3);
                                
                            DO m = 1,3,2
                                
                                k = mod(m,3) + 1
                                        
                                temp = ( pt(1) - t(1,m) ) * ( t(2,k) - t(2,m) ) - &
                                ( pt(2) - t(2,m) ) * ( t(1,k) - t(1,m) )
                                                   
                                IF (0.0 < temp) THEN
                                    inside = .FALSE.
                                    exit
                                ENDIF
                  
                                inside = .TRUE.
                            ENDDO
                                
                            IF (inside .EQV. .TRUE.) THEN
                                
                                s_phiLS(i,j) = (-a_rbLS*xwhole(i)-&
                                                b_rbLS*ywhole(j) + 1)/c_rbLS
                                                                   
                            ENDIF
                        ENDDO                
                    ENDDO
                ENDDO
            ENDIF
        ENDIF                       


        IF (trackphase .EQV. .TRUE.) THEN

            ! --------------------------------------- !
            ! ---- LEVEL SET BOUNDARY CONDITIONS ---- !
            ! --------------------------------------- !

            IF (LS_Type /= 3) THEN
            
                ! Reinitialize as a signed distance function & enforce boundary conditions
                CALL LSREIN(Nx,Ny,hx,hy,phiLS,H,&
                            BCphisouth,BCphinorth,BCphiwest,BCphieast,ywhole,d_phiLS)
        
            ENDIF

            ! ----------------------------- !
            ! ---- HEAVIESIDE FUNCTION ---- !
            ! ------------------------------!

            volume = 0.0
            DO i=1,Nx+2
                DO j=1,Ny+2
                    ! Heavieside function
                    H(i,j) = LSHEAVY(Nx,Ny,i,j,hx,hy,phiLS)
                    volume = volume + H(i,j)*hx*hy
                ENDDO
            ENDDO
            if (LS_type == 1) then
                volume0 = volume           
            endif
 
            ! --------------------------!                  
            ! ---- VOLUME FRACTION ---- !
            ! --------------------------!

            DO i=2,Nx+1
                DO j=2,Ny+1 
                    Avof(i,j) = phiLS(i,j)
                    Bvof(i,j) = (phiLS(i+1,j)-phiLS(i-1,j))/(2.0*hx)
                    Cvof(i,j) = (phiLS(i,j+1)-phiLS(i,j-1))/(2.0*hy)
                    
                    alpha_1 = xwhole(i)-hx/2.0
                    alpha_2 = xwhole(i)+hx/2.0
                    beta_1 = ywhole(j)-hy/2.0
                    beta_2 = ywhole(j)+hy/2.0
    
                    !F(i,j) = volfraction(Avof(i,j),Bvof(i,j),Cvof(i,j), &
                    !                     alpha_1,alpha_2,beta_1,beta_2,hx,hy,xwhole(i),ywhole(j))
                ENDDO
            ENDDO
            Avof(:,1)=Avof(:,2)
            Avof(:,Ny+2)=Avof(:,Ny+1)
            Avof(1,:)=Avof(2,:)
            Avof(Nx+2,:)=Avof(Nx+1,:)
            Bvof(:,1)=Bvof(:,2)
            Bvof(:,Ny+2)=Bvof(:,Ny+1)
            Bvof(1,:)=Bvof(2,:)
            Bvof(Nx+2,:)=Bvof(Nx+1,:)
            Cvof(:,1)=Cvof(:,2)
            Cvof(:,Ny+2)=Cvof(:,Ny+1)
            Cvof(1,:)=Cvof(2,:)
            Cvof(Nx+2,:)=Cvof(Nx+1,:)
            F(:,1)=F(:,2)
            F(:,Ny+2)=F(:,Ny+1)
            F(1,:)=F(2,:)
            F(Nx+2,:)=F(Nx+1,:)

            ! -----------------------------------------------!
            ! ---- INITIAL FLUID MASS (TO BE CONSERVED) ---- !
            ! -----------------------------------------------!
        
            MassInt = 0.0
            DO i=2,Nx+1
                DO j=2,Ny+1
                    
                    MassInt = (rho_one*H(i,j) + rho_two*(1.0-H(i,j)))*hx*hy + MassInt
                
                ENDDO 
            ENDDO

        else
            
            H = 1.0    
        
        ENDIF        

        ! -------------------------------------- !
        ! ---- HEAVISIDE FUNCTION OF SOLIDS ---- !
        ! -------------------------------------- !
        
        ! The rigid body
        IF (rigidbody .EQV. .TRUE.) THEN
            DO i=1,Nx+2
                DO j=1,Ny+2
                    ! Heavieside function
                    s_H(i,j) = LSHEAVY(Nx,Ny,i,j,hx,hy,s_phiLS)
                    !s_H(i,j) = 0.0
                ENDDO
            ENDDO
        else
          
            s_H = 1.0      
        
        ENDIF        

        ! The wedge
        IF (wedge .EQV. .TRUE.) THEN
            DO i=1,Nx+2
                DO j=1,Ny+2
                    ! Heavieside function
                    w_H(i,j) = LSHEAVY(Nx,Ny,i,j,hx,hy,w_phiLS)
                    !s_H(i,j) = 0.0
                ENDDO
            ENDDO
        ENDIF
        
        ! ----------------------------------- !
        ! ---- INITIAL WORK CALCULATIONS ---- !
        ! ----------------------------------- !  

        IF (wedge .EQV. .TRUE.) THEN 
          
            DO i=2,Nx+1
                DO j=2,Ny     
              
                    ! density at the half node
                    rho_O = 0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1)))+&
                            0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))

                    IF (0.5*(w_H(i,j)+w_H(i,j+1)) > 0) THEN                                      
                
                        vstar_int = 0.0
                
                        n_w_inert = n1_w_inert + &
                                    0.5*(w_H(i,j)+w_H(i,j+1))*&
                                    rho_O/DeltaT*(v(i,j) - v_old(i,j))*&
                                    hx*hy
                                                                                        
                        n_w_conv = n1_w_conv + &
                                   0.5*(w_H(i,j)+w_H(i,j+1))*&
                                   rho_O*vstar_int*hx*hy
                                
                        n_w_grav = n1_w_grav + 0.5*&
                                   (w_H(i,j)+w_H(i,j+1))*rho_O*hx*hy*gy

                    ENDIF
              
                    wedge_force = n_w_inert + n_w_conv + n_w_grav
                    wedge_work = wedge_work - (wedge_force * n1_w_vy * deltaT)

                ENDDO
            ENDDO
        ENDIF
             
        ! Fluid Energy Contribution
        IF (rigidbody .EQV. .FALSE.) THEN
            DO i=2,Nx+1
                DO j=2,Ny+1
  
                    IF (i < Nx+1) THEN
                        rho_O = 0.5*(rho_one*H(i+1,j) + rho_two*(1.0-H(i+1,j)))+&
                                0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))
                                                 
                        total_energy = total_energy + 0.5*rho_O*(u(i,j)**2)*hx*hy  ! kinetic (in x)
                        lin_momx = lin_momx + rho_O*u(i,j)*hx*hy
                    ENDIF
              
                    IF (j < Ny+1) THEN
                        rho_O = 0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1)))+&
                                0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))                          

                        total_energy = total_energy + 0.5*rho_O*(v(i,j)**2)*hx*hy ! kinetic (in y)
                        lin_momy = lin_momy + rho_O*v(i,j)*hx*hy
                    ENDIF

                    IF ((i < Nx+1) .AND. (j < Ny+1)) THEN 
                        
                        rho_O = rho_one*H(i,j) + rho_two*(1.0-H(i,j))

                        total_energy = total_energy + &
                                       ( -P(i,j) + rho_O*gy*(j-1)*hy ) * &
                                       hx*hy
              
                    ENDIF

                ENDDO
            ENDDO
        ENDIF
        
        ! ------------------------------------------ !
        ! ---- INITIAL SET OF RIGID BODY FORCES ---- !
        ! ------------------------------------------ !   

        IF (rigidbody .EQV. .TRUE.) THEN
          
            ! Collect some integrals
                
            DO i=2,Nx
                DO j=2,Ny+1
                        
                    ! density at the half node
                    rho_O = 0.5*(rho_one*H(i+1,j) + rho_two*(1.0-H(i+1,j)))+&
                            0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))
                                                 
                    total_energy = total_energy + 0.5*rho_O*(u(i,j)**2)*hx*hy  ! kinetic (in x)
                    lin_momx = lin_momx + rho_O*u(i,j)*hx*hy

                    IF (0.5*(s_H(i,j)+s_H(i+1,j)) > 0) THEN
                                                
                        ! dummy variable for the convective opperator
                        ustar_int= 0.0
                                                
                        n1_xint_inert = n1_xint_inert+&
                                        0.5*(s_H(i,j)+s_H(i+1,j))*&
                                        rho_O/DeltaT*(u(i,j) - u_old(i,j))*&
                                        hx*hy
                                                
                        n1_xint_conv = n1_xint_conv+&
                                       0.5*(s_H(i,j)+s_H(i+1,j))*&
                                       rho_O*ustar_int*hx*hy
                                
                        n1_xint_grav = n1_xint_grav+&
                                       0.5*(s_H(i,j)+s_H(i+1,j))*&
                                       rho_O*hx*hy*gx
                                                
                        ! Rotations
                        ! rs is the s-distance from the location at i,j to the centroid of the object.
                        ry = ywhole(j) - ncgy
                                                
                        r1_int_inert = r1_int_inert - &
                                       ry*0.5*(s_H(i,j)+s_H(i+1,j))*&
                                       rho_O/deltaT*(u(i,j) - u_old(i,j))*&
                                       hx*hy
                                                
                        r1_int_conv = r1_int_conv - &
                                      ry*0.5*(s_H(i,j)+s_H(i+1,j))*&
                                      rho_O*ustar_int*hx*hy
                                
                        r1_int_grav = r1_int_grav - &
                                      ry*0.5*(s_H(i,j)+s_H(i+1,j))*&
                                      rho_O*hx*hy*gx
                                
                    ENDIF

                ENDDO
            ENDDO


            DO i=2,Nx+1
                DO j=2,Ny
                        
                    ! density at the half node
                    rho_O = 0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1)))+&
                            0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))                          

                    total_energy = total_energy + 0.5*rho_O*(v(i,j)**2)*hx*hy ! kinetic (in y)
                    lin_momy = lin_momy + rho_O*v(i,j)*hx*hy

                    IF (0.5*(s_H(i,j)+s_H(i,j+1)) > 0) THEN
                                        
                        ! dummy variable for the convective operator
                        vstar_int = 0.0
                                                        
                        n_yint_inert = n1_yint_inert +&
                                       0.5*(s_H(i,j)+s_H(i,j+1))*&
                                       rho_O/DeltaT*(v(i,j) - v_old(i,j))*&
                                       hx*hy
                                                                                        
                        n_yint_conv = n1_yint_conv + &
                                      0.5*(s_H(i,j)+s_H(i,j+1))*&
                                      rho_O*vstar_int*hx*hy
                                
                        n_yint_grav = n1_yint_grav + &
                                      0.5*(s_H(i,j)+s_H(i,j+1))*&
                                      rho_O*hx*hy*gy
                                                
                        ! Rotations
                        ! rs is the s-distance from the location at i,j to the centroid of the object.
                        rx = xwhole(i) - ncgx
                                                
                        r_int_inert = r1_int_inert + &
                                      rx*0.5*(s_H(i,j)+s_H(i,j+1))*&
                                      rho_O/deltaT*(v(i,j) - v_old(i,j))*&
                                      hx*hy
                                                                                        
                        r_int_conv = r1_int_conv + &
                                     rx*0.5*(s_H(i,j)+s_H(i,j+1))*&
                                     rho_O*vstar_int*hx*hy
                                
                        r_int_grav = r1_int_grav + &
                                     rx*0.5*(s_H(i,j)+s_H(i,j+1))*&
                                     rho_O*hx*hy*gy
                    ENDIF                                  
                
                ENDDO
            ENDDO

            ! Pressure and gravitational contributions to total energy
            DO i=2,Nx
                DO j=2,Ny     
              
                    ! density at the whole node
                    rho_O = rho_one*H(i,j) + rho_two*(1.0-H(i,j))

                    total_energy = total_energy + &
                                   ( -P(i,j) + rho_O*gy*(j-1)*hy ) * &
                                   hx*hy

                ENDDO
            ENDDO
            
            ! ---------------------------------------- !
            ! ---- INITIAL SET OF INTERNAL FORCES ---- !
            ! ---------------------------------------- !   
          
            ! Initial internal forces          
            IF (int_type == 1) THEN                               ! Torque exerted between disc and RB
              
                small_theta_n = th
                delta_theta = small_theta_n - th
                delta_vom = n_rb2_vom - n_rb_vom
                int_torque = k12*delta_theta + c12*delta_vom 
                rb_data(50) = cgx                                 ! Save data for post-processing
                rb_data(51) = cgy
                rb_data(52) = small_theta_n
                
                ! Acceleration
                n_rb2_aom = -int_torque/small_J
    
                ! Energy
                int_work = int_work + c12*(delta_vom**2)*deltaT
                total_energy = total_energy + &
                               0.5*small_J*(n_rb2_vom**2) + &
                               0.5*k12*(delta_theta**2)
     
            ELSEIF (int_type == 2) THEN

                kyadj = ky12 !*exp(log(kscale)*abs(m_d(2))/ybound)

                ! Gravitational and inertial forces acting on internal mass
                Fi(1) = slider_mass*( cos(th)*(- n_rb_ax - 0.0)+sin(th)*(- n_rb_ay - 0.0))
                Fi(2) = slider_mass*(-sin(th)*(- n_rb_ax - 0.0)+cos(th)*(- n_rb_ay - 0.0))
              
                Ki(1,1) = kx12
                Ki(1,2) = 0.0
                Ki(2,1) = 0.0
                Ki(2,2) = kyadj + slider_mass*(buoy_J*n1_rb_vom)**2
              
                Ci(1,1) = 0.0
                Ci(1,2) = 2.0*slider_mass*buoy_J*n1_rb_vom
                Ci(2,1) = 0.0
                Ci(2,2) = cy12
                
                ! Initial Acceleration, m
                m_a(1) = (Fi(1) - Ci(1,1)*m_v(1) - Ci(1,2)*m_v(2) - Ki(1,1)*m_d(1) - Ki(1,2)*m_d(2))/slider_mass
                m_a(2) = (Fi(2) - Ci(2,1)*m_v(1) - Ci(2,2)*m_v(2) - Ki(2,1)*m_d(1) - Ki(2,2)*m_d(2))/slider_mass
                
                ! Intermediate Displacement, m+1/2
                m12_d(1) = m_d(1) + deltaT*m_v(1) + deltaT**2.0/4.0*m_a(1)
                m12_d(2) = m_d(2) + deltaT*m_v(2) + deltaT**2.0/4.0*m_a(2)
                
                ! Intermediate Velocity, m+1/2
                m12_v(1) = m_v(1) + deltaT/2.0*m_a(1)
                m12_v(2) = m_v(2) + deltaT/2.0*m_a(2)
                
                ! Inverse Matrix
                MAT1(1,1) = slider_mass*ID(1,1) + deltaT/2.0*Ci(1,1) + deltaT**2.0/4.0*Ki(1,1)
                MAT1(1,2) = slider_mass*ID(1,2) + deltaT/2.0*Ci(1,2) + deltaT**2.0/4.0*Ki(1,2)
                MAT1(2,1) = slider_mass*ID(2,1) + deltaT/2.0*Ci(2,1) + deltaT**2.0/4.0*Ki(2,1)
                MAT1(2,2) = slider_mass*ID(2,2) + deltaT/2.0*Ci(2,2) + deltaT**2.0/4.0*Ki(2,2)
              
                dmt = MAT1(1,1)*MAT1(2,2)
              
                INV(1,1) = MAT1(2,2)/dmt
                INV(1,2) = -MAT1(1,2)/dmt
                INV(2,1) = -MAT1(2,1)/dmt
                INV(2,2) = MAT1(1,1)/dmt 
                
                ! New Acceleration, m+1
                m1_a(1) = INV(1,1)*(Fi(1) - Ci(1,1)*m12_v(1) - Ci(1,2)*m12_v(2) - Ki(1,1)*m12_d(1) - Ki(1,2)*m12_d(2)) + &
                          INV(1,2)*(Fi(2) - Ci(2,1)*m12_v(1) - Ci(2,2)*m12_v(2) - Ki(2,1)*m12_d(1) - Ki(2,2)*m12_d(2))
                m1_a(2) = INV(2,1)*(Fi(1) - Ci(1,1)*m12_v(1) - Ci(1,2)*m12_v(2) - Ki(1,1)*m12_d(1) - Ki(1,2)*m12_d(2)) + &
                          INV(2,2)*(Fi(2) - Ci(2,1)*m12_v(1) - Ci(2,2)*m12_v(2) - Ki(2,1)*m12_d(1) - Ki(2,2)*m12_d(2))
                        
                ! New Displacement, m+1
                m1_d(1) = m12_d(1) + deltaT**2.0/4.0*m1_a(1)
                m1_d(2) = m12_d(2) + deltaT**2.0/4.0*m1_a(2)
              
                ! New Velocity, m+1
                m1_v(1) = m12_v(1) + deltaT/2.0*m1_a(1)
                m1_v(2) = m12_v(2) + deltaT/2.0*m1_a(2)              
                
                m_d(1) = m1_d(1)
                m_d(2) = m1_d(2)
                m_v(1) = m1_v(1)
                m_v(2) = m1_v(2)
                m_a(1) = m1_a(1)
                m_a(2) = m1_a(2)
                
                Fk(1) = kx12*m_d(1)
                Fk(2) = kyadj*m_d(2) + cy12*m_v(1)
                
                int_forcex = cos(th)*Fk(1) - sin(th)*Fk(2)
                int_forcey = sin(th)*Fk(1) + cos(th)*Fk(2)
                int_torque = m_d(2)*Fk(1)
                
                scgx = x_eqm + cos(th)*m_d(1) - sin(th)*m_d(2)
                scgy = y_eqm + sin(th)*m_d(1) + cos(th)*m_d(2)
                
                rb_data(50) = ncgx
                rb_data(51) = ncgy
                rb_data(52) = scgx
                rb_data(53) = scgy 
                
                print *, ''
                print *, th
                print *, ncgx
                print *, ncgy
                print *, m_d(1)
                print *, m_d(2)
                
                  
            ENDIF
      
            ! --------------------------------- !
            ! ---- RIGID BODY ACCELERATION ---- !
            ! --------------------------------- !   
      
            ! First set of rigid body accelerations
            n_rb_ax = (n_xint_inert + n_xint_conv + n_xint_grav + int_forcex)/rb_mass - gx        
            n_rb_ay = (n_yint_inert + n_yint_conv + n_yint_grav + int_forcey)/rb_mass - gy        
            n_rb_aom = (r_int_inert + r_int_conv + r_int_grav + int_torque)/buoy_J

            ! Restrict motion
            IF (RB_lock == 1) THEN
                n_rb_ax = 0.0
                n_rb_aom = 0.0
                n1_rb_ax = 0.0
                n1_rb_aom = 0.0
            ENDIF

            ! Energy Update
            total_energy = total_energy + &
                           rb_mass*gy*cgy + &
                           0.5*rb_mass*((n_rb_vx**2)+(n_rb_vy**2)) + &
                           0.5*buoy_J*(n_rb_vom**2)
            lin_momx = lin_momx + rb_mass*n_rb_vx
            lin_momy = lin_momy + rb_mass*n_rb_vy
          

        ENDIF

        ! ---- SAVE INITIAL DATA ----
        !call SAVEDATA(Nx,Ny,Lx,Ly,time,x,y,phiLS,s_phiLS,w_phiLS,H,s_H,w_H,P,u,v,file_count,rb_data)
        call SAVEDATA2(Nx,Ny,Lx,Ly,time,xwhole,ywhole,phiLS,s_phiLS,w_phiLS,H,s_H,w_H,P,u,v, &
                       file_count,rb_data,phi_x,phi_y,phi_xy,Fsig)
        
    !**************************************************************************************************!
    !****************************************** SOLVE EQUATIONS ***************************************!
    !**************************************************************************************************!
    !**************************************************************************************************!
    !****************************************                     *************************************!
    !****************************************   START TIME LOOP   *************************************!
    !****************************************                     *************************************!
    !**************************************************************************************************!
    !**************************************************************************************************!
    !**************************************************************************************************!    
                               

    DO  ! time loop
        
        !---------------------------!
        !---- BENCHMARKING DATA ----!
        !---------------------------!
        
        ! Locate wave peak
        IF (wave .eqv. .true.) THEN
        
            maxi = 1
            maxj = 1
            DO i = 2,Nx+1
                DO j = 2,Ny+1
                    
                    IF ((phiLS(i,j) >= 0.0) .AND. (j > maxj)) THEN
                        maxi = i
                        maxj = j
                    ENDIF
                                            
                ENDDO
            ENDDO
            
    	    ! Approximate interface location
            tempa(1) = -phiLS(maxi,maxj)*hy/(phiLS(maxi,maxj+1) - phiLS(maxi,maxj))
            tempa(2) = ywhole(maxj)+tempa(1)

		! Locate free surface position peak
        ELSEIF (wedge .eqv. .true.) THEN
        
            maxi = Nx/2
            maxj = 1

        	DO j = 2,Ny+1
                    
                    IF ((phiLS(i,j) >= 0.0) .AND. (j > maxj)) THEN
                        maxj = j
                    ENDIF
                                            
            ENDDO
            
    	    ! Approximate interface location
            tempa(1) = -phiLS(maxi,maxj)*hy/(phiLS(maxi,maxj+1) - phiLS(maxi,maxj))
            tempa(2) = ywhole(maxj)+tempa(1)
                        
        ! Slosh benchmark data
        ELSEIF (slosh .eqv. .true.) THEN
        
            maxi = 2
            maxj = 1
            DO j = 2,Ny+1
                IF ((phiLS(i,j) >= 0.0) .AND. (j > maxj)) THEN
                    maxj = j
                ENDIF
            ENDDO

            testa(1) = ywhole(maxj-1)
            testa(2) = phiLS(maxi,maxj-1)
            testb(1) = ywhole(maxj+2)
            testb(2) = phiLS(maxi,maxj+2)
 
            DO j=1,25
        
                tempa(1) = xwhole(maxi)
                tempa(2) = (testa(1)+testb(1))/2.0

                CALL HERMCUBE(Nx,Ny,hx,hy,xwhole,ywhole,phiLS,phi_x,phi_y,phi_xy,herm, &
                              GA,scal,tempa,maxi,maxj,1)

                if (herm(1)*testa(2) < 0.0) then
                    testb(1) = tempa(2)
                    testb(2) = herm(1)
                elseif (herm(1)*testb(2) < 0.0) then
                    testa(1) = tempa(2)
                    testa(2) = herm(1)
                elseif (abs(herm(1)) < 1e-16) then
                    continue
                endif

            ENDDO

 
    	    ! Approximate interface location
            !tempa(1) = -phiLS(maxi,maxj)*hy/(phiLS(maxi,maxj+1) - phiLS(maxi,maxj))
            !tempa(2) = ywhole(maxj)+tempa(1)
            tempa(1) = tempa(2)
            tempa(2) = herm(1)
            
        ! Locate column front position
        ELSEIF (dambreak .eqv. .true.) THEN
        
            maxi = 1
            maxj = 1
            DO j = 2,Ny+1
            	DO i = 2,Nx+1
            	    IF ((phiLS(i,j) >= 0.0) .AND. (i > maxi)) THEN
            	        maxi = i
            	        maxj = j
            	    ENDIF
            	ENDDO
            ENDDO
            
            if (maxi > Nx) then
            	tempa(1) = 0
            	tempa(2) = 0
            else    
    	        ! Approximate interface location
            	tempa(1) = -phiLS(maxi,maxj)*hx/(phiLS(maxi+1,maxj) - phiLS(maxi,maxj))
            	tempa(2) = xwhole(maxi)+tempa(1)
            endif
            
        ENDIF
        
        !Calculate net momentum
        net_mom = sqrt((lin_momx**2) + (lin_momy**2))
        
        !-------------------------------!
        !---- WRITE RB DATA TO FILE ----!
        !-------------------------------!

        ! CREATING DATA FILE
        OPEN (UNIT=8,FILE='centers.dat', POSITION ='append', IOSTAT=ios)
        WRITE (UNIT=8,FMT=90) ncgx, ncgy, n_rb_vx, n_rb_vy, &
                              n_yint_grav, n_yint_conv, n_yint_inert
        CLOSE (UNIT=8)
                
        OPEN (UNIT=7,FILE='bench.dat', POSITION ='append', IOSTAT=ios)
        WRITE (UNIT=7,FMT=90) time, real(maxi), tempa(1), tempa(2), volume 
        CLOSE (UNIT=7)
          
        OPEN (UNIT=7,FILE='time_step.dat', POSITION ='append', IOSTAT=ios)
        WRITE (UNIT=7,FMT=90) deltaT, stability 
        CLOSE (UNIT=7)

        IF ((rigidbody .EQV. .TRUE.) .and. (int_type == 1)) THEN
            OPEN (UNIT=7,FILE='disk.dat', POSITION ='append', IOSTAT=ios)
            WRITE (UNIT=7,FMT=90) deltaT, th, small_theta_n, n_rb_aom, int_torque
            CLOSE (UNIT=7)

        ELSEIF ((rigidbody .EQV. .TRUE.) .and. (int_type == 2)) THEN
            OPEN (UNIT=7,FILE='trans.dat', POSITION ='append', IOSTAT=ios)
            WRITE (UNIT=7,FMT=90) deltaT, m_d(1), m_d(2), m_v(1), m_v(2), m_a(1), m_a(2)
            CLOSE (UNIT=7)
          
        ENDIF

        !---------------------------!
        !---- FIND MAX VELOCITY ----!
        !---------------------------!

        umax = max(unorth,usouth)
        vmax = max(vwest,veast)
        DO i=1,Nx+1
            DO j=1,Ny+2
                umax = max(abs(u(i,j)),umax) ! largest u velocity
            ENDDO
        ENDDO
        DO i=1,Nx+2 
            DO j=1,Ny+1
                vmax = max(abs(v(i,j)),vmax) ! largest velocity
            ENDDO
        ENDDO                


        !-----------------------------!
        !---- CALCULATE TIME STEP ----!  
        !-----------------------------!

        ! IF solution blows up or wrinkles form, THEN reduce stability parameter
        dtvisc=1.0/( max(mu_one/rho_one,mu_two/rho_two)*&
                     (2.0/(hx**2.0) + 2.0/(hy**2.0)) )
        dtadv = 1.0/(umax/hx + vmax/hy)
        dtg=sqrt( min(hx,hy)/max(gx,gy,1.0e-16) )
        dtsigma=sqrt( (rho_one+rho_two)*min(hx,hy)**3.0/&
                      (4.0*pi*sigma+1.0e-16) )
        deltaT=min(dtvisc,dtadv,dtg,dtsigma,dtinternal)
        
        IF (LS_type == 3) THEN
            deltaT = min(deltaT,hx/umax,hy/vmax)
        ENDIF
        
        deltaT=stability*deltaT  ! stability must be less than 0.9
        deltaT = min(deltaT,time_print)
        
        print *, "deltaT = ",deltaT
        
        !deltaT = hx

        !-------------------------------!
        !---- RIGID BODY KINEMATICS ----!
        !-------------------------------!
        !---- Velocity: n --> n + 1/2   !  
        !---- Disp:     n --> n + 1     !  
        !-------------------------------!
          
        ! Calculate Rigid Body motion
        IF (rigidbody .EQV. .TRUE.) THEN
            
            ! Check for motion restrictions      
            IF (RB_lock /= 1) THEN

                ! Compute Relative rotation vector 
                big_theta = deltaT*n_rb_vom + 0.5*deltaT**2*n_rb_aom
                
                 c1 = (sin(abs(big_theta)))/abs(big_theta)
                 c2 = 2*(sin(0.5*abs(big_theta))**2)/(big_theta**2)        
                
                 ! The exponential map
                 exp_big_theta(1,1) = 1.0 - c2*big_theta**2
                 exp_big_theta(1,2) = -c1*big_theta
                 exp_big_theta(2,1) =  c1*big_theta
                 exp_big_theta(2,2) = 1.0 - c2*big_theta**2
                
                 ! Rotational velocity at (n+1/2)
                 n12_rb_vom = n_rb_vom + 0.5*deltaT/buoy_J*&
                              (r_int_inert + r_int_conv + r_int_grav)
                
            ENDIF

            ! Calculate the rigid body velocity at v(n+1/2)
            n12_rb_vx = n_rb_vx + 0.5*deltaT*n_rb_ax
            n12_rb_vy = n_rb_vy + 0.5*deltaT*n_rb_ay
            
            ! Calculate d(n+1)                
            n1_rb_dx = n_rb_dx + deltaT*n12_rb_vx 
            n1_rb_dy = n_rb_dy + deltaT*n12_rb_vy
            
            ! Check for motion restrictions   
            IF (RB_lock /= 1) THEN 
                n1_rbdelom(1,1) = n_rbdelom(1,1)*exp_big_theta(1,1) +&
                                  n_rbdelom(1,2)*exp_big_theta(2,1)
                n1_rbdelom(1,2) = n_rbdelom(1,1)*exp_big_theta(1,2) +&
                                  n_rbdelom(1,2)*exp_big_theta(2,2)
                n1_rbdelom(2,1) = n_rbdelom(2,1)*exp_big_theta(1,1) +&
                                  n_rbdelom(2,2)*exp_big_theta(2,1)
                n1_rbdelom(2,2) = n_rbdelom(2,1)*exp_big_theta(1,2) +&
                                  n_rbdelom(2,2)*exp_big_theta(2,2)
            ENDIF
            
            ! Recover change in rigid body's angle
            IF (n1_rbdelom(1,1) < 0) THEN                 ! Expand domain of asin to +/-pi
                th_domain = sign(pi,n1_rbdelom(2,1)) - asin(n1_rbdelom(2,1))
            else    
                th_domain = asin(n1_rbdelom(2,1))            
            ENDIF
            
            IF ((abs(th_old) > pi/2.0) .and. &            ! Track number of revolutions 
               (sign(dble(1),th_old) /= sign(dble(1),th_domain))) THEN
                 
                th_counter = th_counter + sign(dble(1),th_old)            
            ENDIF

            th_change = th_domain + 2.0*pi*th_counter
            th_old = th_domain

            !-------------------------------!
            !----- INTERNAL KINEMATICS -----!
            !-------------------------------!
            !---- Disc: Same as RB          !  
            !---- Slider: Implicit, n->n+1  !  
            !-------------------------------!

            IF (int_type == 1) THEN

                ! Internal rotational velocity at (n+1/2)
                n12_rb2_vom = n_rb2_vom + 0.5*deltaT*n_rb2_aom
                
                small_theta_n1 = small_theta_n + deltaT*n12_rb2_vom 
                th = c_tho + th_change
                delta_theta = small_theta_n1 - th
                delta_vom = n12_rb2_vom - n12_rb_vom  
                int_torque = k12*delta_theta + c12*delta_vom 
                rb_data(50) = ncgx                                 ! Save data for post-processing
                rb_data(51) = ncgy
                rb_data(52) = small_theta_n1

                ! Acceleration
                n1_rb2_aom = -int_torque/small_J

                ! Update velocity
                n1_rb2_vom = n12_rb2_vom + 0.5*deltaT*n1_rb2_aom
           
                ! Energy update
                int_work = int_work + c12*(delta_vom**2)*deltaT
                total_energy = total_energy + &
                               0.5*small_J*(n1_rb2_vom**2) + &
                               0.5*k12*(delta_theta**2)

            ELSEIF (int_type == 2) THEN
              
                kyadj = ky12 !*exp(log(kscale)*abs(m_d(2))/ybound)
                
                ! Gravitational and inertial forces acting on internal mass
                Fi(1) = slider_mass*( cos(th)*(- n_rb_ax - 0.0)+sin(th)*(- n_rb_ay - 0.0))
                Fi(2) = slider_mass*(-sin(th)*(- n_rb_ax - 0.0)+cos(th)*(- n_rb_ay - 0.0))
              
                Ki(1,1) = kx12
                Ki(1,2) = 0.0
                Ki(2,1) = 0.0
                Ki(2,2) = kyadj + slider_mass*(buoy_J*n1_rb_vom)**2.0
              
                Ci(1,1) = 0.0
                Ci(1,2) = 2.0*slider_mass*buoy_J*n1_rb_vom
                Ci(2,1) = 0.0
                Ci(2,2) = cy12
              
                ! Initial Acceleration, m
                m_a(1) = (Fi(1) - Ci(1,1)*m_v(1) - Ci(1,2)*m_v(2) - Ki(1,1)*m_d(1) - Ki(1,2)*m_d(2))/slider_mass
                m_a(2) = (Fi(2) - Ci(2,1)*m_v(1) - Ci(2,2)*m_v(2) - Ki(2,1)*m_d(1) - Ki(2,2)*m_d(2))/slider_mass
                
                ! Intermediate Displacement, m+1/2
                m12_d(1) = m_d(1) + deltaT*m_v(1) + deltaT**2.0/4.0*m_a(1)
                m12_d(2) = m_d(2) + deltaT*m_v(2) + deltaT**2.0/4.0*m_a(2)
              
                ! Intermediate Velocity, m+1/2
                m12_v(1) = m_v(1) + deltaT/2.0*m_a(1)
                m12_v(2) = m_v(2) + deltaT/2.0*m_a(2)
              
                ! Inverse Matrix
                MAT1(1,1) = slider_mass*ID(1,1) + deltaT/2.0*Ci(1,1) + deltaT**2.0/4.0*Ki(1,1)
                MAT1(1,2) = slider_mass*ID(1,2) + deltaT/2.0*Ci(1,2) + deltaT**2.0/4.0*Ki(1,2)
                MAT1(2,1) = slider_mass*ID(2,1) + deltaT/2.0*Ci(2,1) + deltaT**2.0/4.0*Ki(2,1)
                MAT1(2,2) = slider_mass*ID(2,2) + deltaT/2.0*Ci(2,2) + deltaT**2.0/4.0*Ki(2,2)
                
                dmt = MAT1(1,1)*MAT1(2,2)
                
                INV(1,1) = MAT1(2,2)/dmt
                INV(1,2) = -MAT1(1,2)/dmt
                INV(2,1) = -MAT1(2,1)/dmt
                INV(2,2) = MAT1(1,1)/dmt 
                
                ! New Acceleration, m+1
                m1_a(1) = INV(1,1)*(Fi(1) - Ci(1,1)*m12_v(1) - Ci(1,2)*m12_v(2) - Ki(1,1)*m12_d(1) - Ki(1,2)*m12_d(2)) + &
                          INV(1,2)*(Fi(2) - Ci(2,1)*m12_v(1) - Ci(2,2)*m12_v(2) - Ki(2,1)*m12_d(1) - Ki(2,2)*m12_d(2))
                m1_a(2) = INV(2,1)*(Fi(1) - Ci(1,1)*m12_v(1) - Ci(1,2)*m12_v(2) - Ki(1,1)*m12_d(1) - Ki(1,2)*m12_d(2)) + &
                          INV(2,2)*(Fi(2) - Ci(2,1)*m12_v(1) - Ci(2,2)*m12_v(2) - Ki(2,1)*m12_d(1) - Ki(2,2)*m12_d(2))
                
                ! New Displacement, m+1
                m1_d(1) = m12_d(1) + deltaT**2.0/4.0*m1_a(1)
                m1_d(2) = m12_d(2) + deltaT**2.0/4.0*m1_a(2)
              
                ! New Velocity, m+1
                m1_v(1) = m12_v(1) + deltaT/2.0*m1_a(1)
                m1_v(2) = m12_v(2) + deltaT/2.0*m1_a(2)              
              
                m_d(1) = m1_d(1)
                m_d(2) = m1_d(2)
                m_v(1) = m1_v(1)
                m_v(2) = m1_v(2)
                m_a(1) = m1_a(1)
                m_a(2) = m1_a(2)
                
                Fk(1) = kx12*m_d(1)
                Fk(2) = kyadj*m_d(2) + cy12*m_v(1)
              
                int_forcex = cos(th)*Fk(1) - sin(th)*Fk(2)
                int_forcey = sin(th)*Fk(1) + cos(th)*Fk(2)
                int_torque = m_d(2)*Fk(1)
              
                x_eqm = ncgx - sin(th)*y_diff 
                y_eqm = ncgy + cos(th)*y_diff
 
                scgx = x_eqm + cos(th)*m_d(1) - sin(th)*m_d(2)
                scgy = y_eqm + sin(th)*m_d(1) + cos(th)*m_d(2)
              
                rb_data(50) = ncgx
                rb_data(51) = ncgy
                rb_data(52) = scgx
                rb_data(53) = scgy 
              
                !print *, ''
                !print *, th
                !print *, ncgx
                !print *, ncgy
                !print *, m_d(1)
                !print *, m_d(2)
            
            ENDIF     

        ENDIF
                
        !----------------------------------------------------------------------------!
        !------------------------------ VELOCITY UPDATE -----------------------------!
        !----------------------------------------------------------------------------!          

        ! apply advection BC's
        call BCVELOCITY(Nx,Ny,u,v,BCnorth,BCsouth,BCeast,BCwest,&
                        unorth,usouth,ueast,uwest,&
                        vnorth,vsouth,veast,vwest,&
                        hx,hy,phiLS,contC,contA,contR,contM,1)
                                                
        !-------------------------------------------!          
        ! ---- FORCE THE FLUID WITH THE SOLIDS ---- !
        !-------------------------------------------!          

        IF (rigidbody .EQV. .TRUE.) THEN          
            
            ! Enforce the rigid body boundary condition on fluid velocity                 
            DO i=1,Nx+1
                DO j=1,Ny+2
                                                
                    IF ( 0.5*(s_phiLS(i+1,j) + s_phiLS(i,j)) >= 0) THEN

                        ry = ywhole(j) - ncgy

                        ! new velocity, n+1/2 time level
                        u(i,j) = n12_rb_vx - n12_rb_vom*ry

                    ENDIF
                ENDDO
            ENDDO

            DO i=1,Nx+2
                DO j=1,Ny+1
                                
                    IF ( 0.5*(s_phiLS(i,j) + s_phiLS(i,j+1)) >= 0) THEN

                        rx = xwhole(i) - ncgx

                        ! new velocity, n+1/2 time level
                        v(i,j) = n12_rb_vy + n12_rb_vom*rx

                    ENDIF
                ENDDO
            ENDDO
            
        ENDIF

        IF (wedge .EQV. .TRUE.) THEN
        
            DO i=1,Nx+1
                DO j=1,Ny+2
                                                
                    IF ( 0.5*(w_phiLS(i+1,j) + w_phiLS(i,j)) >= 0) THEN

                    ! new velocity, n+1/2 time level
                    u(i,j) = 0.0

                    ENDIF
                ENDDO
            ENDDO

            DO i=1,Nx+2
                DO j=1,Ny+1
                                
                    IF ( 0.5*(w_phiLS(i,j) + w_phiLS(i,j+1)) >= 0) THEN

                    ! new velocity, n+1/2 time level
                    v(i,j) = n1_w_vy

                    ENDIF
                ENDDO
            ENDDO
            
        ENDIF

        ! ---- BEACH ----
        ! Requires dimensions: 9.114 x ~0.7595 (Actual Tank Dimensions)
        ! Remains stationary so no need for its own Heaviside function
        if (beach) then
            DO i=1,Nx+1
                DO j=1,Ny+2

                    if (xwhole(i) .GE. 8.942 .AND. ywhole(j) .LE. 0.520) then
                        u(i,j) = 0.0
                    elseif (xwhole(i) < 8.942 .AND. ywhole(j) .LE. 0.5019*xwhole(i) - 3.9683) then
                        u(i,j) = 0.0
                    endif

                ENDDO
            ENDDO
            DO i=1,Nx+2
                DO j=1,Ny+1

                    if (xwhole(i) .GE. 8.942 .AND. ywhole(j) .LE. 0.520) THEN
                        v(i,j) = 0.0
                    elseif (xwhole(i) < 8.942 .AND. ywhole(j) .LE. 0.5019*xwhole(i) - 3.9683) then
                        v(i,j) = 0.0
                    endif

                ENDDO
            ENDDO
        endif

        ! ---- REEF ----
        ! Requires dimensions: 40.0 x 2.0
        ! Remains stationary so no need for its own Heaviside function
        ! Intended for use with wave option
        if (reef) then
            DO i=1,Nx+1
                DO j=1,Ny+2

                    if (ywhole(j) .LE. 0.85 .AND. xwhole(i) .GE. 20.0) then
                        u(i,j) = 0.0
                    endif

                ENDDO
            ENDDO
            DO i=1,Nx+2
                DO j=1,Ny+1

                    if (ywhole(j) .LE. 0.85 .AND. xwhole(i) .GE. 20.0) then
                        v(i,j) = 0.0
                    endif

                ENDDO
            ENDDO
        endif

        !----- ENSURE THAT INTERMEDIATE VELOCITY IS DEFINED ON BOUNDARIES ----!
        ustar = u
        vstar = v 
                        
        ! ---- CONVECTIVE TERMS WITH NO UPWINDING ---- !
        ! ---- Unstable for very convective problems ---- !
        ! calculating the convection term for the x-momentum equation
        !DO j=2,Ny+1                ! j=Ny+2 and j=1 are boundary nodes for u
        !        DO i=2,Nx        ! i=1 and i=Nx+1 are boundary nodes for u
        !
        !                ! dummy variable for the convective opperator
        !                ustar(i,j)=((u(i+1,j)+u(i,j))**2.0-(u(i,j)+u(i-1,j))**2.0)/hx*0.25 + &
        !                                   ((u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j)) - &
        !                                   (u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1)))/hy*0.25
        !                ! i corresponds to half nodes
        !                ! j corresponds to whole nodes        
        !        ENDDO
        !ENDDO
              
        ! calculating the convection term for the y-momentum equation
        !DO j=2,Ny                ! j=1 and j=Ny+1 are boundary nodes for v
        !        DO i=2,Nx+1        ! i=1 and i=Nx+2 are boundary nodes for v
        !
        !                ! dummy variable for the convective opperator
        !                vstar(i,j)=((u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j)) - &
        !                                   (u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j)))/hx*0.25 + &
        !                                   ((v(i,j+1)+v(i,j))**2.0-(v(i,j)+v(i,j-1))**2.0)/hy*0.25
        !                ! i corresponds to whole nodes
        !                ! j corresponds to half nodes
        !
        !        ENDDO
        !ENDDO

        ! ---- CONVECTIVE TERMS WITH FIRST ORDER UPWINDING ---- !
        ! ---- Diffusive for high velocity gradients ---- !
        ! First order upwinded convective terms
        ! calculating the convection term for the x-momentum equation
        !DO j=2,Ny+1                ! j=Ny+2 and j=1 are boundary nodes for u
        !        DO i=2,Nx        ! i=1 and i=Nx+1 are boundary nodes for u
        !
        !                u_plus = 0.5*(u(i,j) + abs(u(i,j)))
        !                u_minus = 0.5*(u(i,j) - abs(u(i,j)))
        !
        !                v_plus = 0.5*(v(i,j) + abs(v(i,j)))
        !                v_minus = 0.5*(v(i,j) - abs(v(i,j)))
        !                
        !                dudx_minus = (u(i,j) - u(i-1,j))/hx
        !                dudx_plus  = (u(i+1,j) - u(i,j))/hx
        !
        !                dudy_minus = (u(i,j) - u(i,j-1))/hy
        !                dudy_plus  = (u(i,j+1) - u(i,j))/hy
        !
        !                ! dummy variable for the convective opperator
        !                ustar(i,j) = u_plus*dudx_minus + u_minus*dudx_plus + &
        !                             v_plus*dudy_minus + v_minus*dudy_plus
        !                
        !                ! i corresponds to half nodes
        !                ! j corresponds to whole nodes        
        !        ENDDO
        !ENDDO

        ! calculating the convection term for the y-momentum equation
        !DO j=2,Ny                ! j=1 and j=Ny+1 are boundary nodes for v
        !        DO i=2,Nx+1        ! i=1 and i=Nx+2 are boundary nodes for v
        !
        !
        !                u_plus = 0.5*(u(i,j) + abs(u(i,j)))
        !                u_minus = 0.5*(u(i,j) - abs(u(i,j)))
        !
        !                v_plus = 0.5*(v(i,j) + abs(v(i,j)))
        !                v_minus = 0.5*(v(i,j) - abs(v(i,j)))
        !                
        !                dvdx_minus = (v(i,j) - v(i-1,j))/hx
        !                dvdx_plus  = (v(i+1,j) - v(i,j))/hx
        !
        !                dvdy_minus = (v(i,j) - v(i,j-1))/hy
        !                dvdy_plus  = (v(i,j+1) - v(i,j))/hy
        !
        !                ! dummy variable for the convective opperator
        !                vstar(i,j) = u_plus*dvdx_minus + u_minus*dvdx_plus + &
        !                             v_plus*dvdy_minus + v_minus*dvdy_plus
        !                
        !                ! i corresponds to half nodes
        !                ! j corresponds to whole nodes        
        !        ENDDO
        !ENDDO

        ! ---- CONVECTIVE TERMS WITH SECOND ORDER UPWINDING ---- !
        ! First order near the boundaries
        ! Diffusive for high velocity gradients
        ! calculating the convection term for the x-momentum equation
        DO j=2,Ny+1                ! j=Ny+2 and j=1 are boundary nodes for u
            DO i=2,Nx        ! i=1 and i=Nx+1 are boundary nodes for u
                
                u_plus = 0.5*(u(i,j) + abs(u(i,j)))
                u_minus = 0.5*(u(i,j) - abs(u(i,j)))

                v_plus = 0.5*(v(i,j) + abs(v(i,j)))
                v_minus = 0.5*(v(i,j) - abs(v(i,j)))
                                                        
                IF (i == 2) THEN
                    dudx_minus = (u(i,j) - u(i-1,j))/hx
                    dudx_plus  = (u(i+1,j) - u(i,j))/hx
                ELSEIF (i == Nx) THEN
                    dudx_minus = (u(i,j) - u(i-1,j))/hx
                    dudx_plus  = (u(i+1,j) - u(i,j))/hx
                else
                    dudx_minus = (3*u(i,j) - 4*u(i-1,j) + u(i-2,j))/(2*hx)
                    dudx_plus  = (-u(i+2,j) + 4*u(i+1,j) - 3*u(i,j))/(2*hx)
                ENDIF                                   

                IF (j == 2) THEN
                    dudy_minus = (u(i,j) - u(i,j-1))/hy
                    dudy_plus  = (u(i,j+1) - u(i,j))/hy
                ELSEIF (j == Ny+1) THEN
                    dudy_minus = (u(i,j) - u(i,j-1))/hy
                    dudy_plus  = (u(i,j+1) - u(i,j))/hy
                else
                    dudy_minus = (3*u(i,j) - 4*u(i,j-1) + u(i,j-2))/(2*hy)
                    dudy_plus  = (-u(i,j+2) + 4*u(i,j+1) - 3*u(i,j))/(2*hy)
                ENDIF
 
                ! dummy variable for the convective opperator
                ustar(i,j) = u_plus*dudx_minus + u_minus*dudx_plus + &
                             v_plus*dudy_minus + v_minus*dudy_plus
                                
                ! i corresponds to half nodes
                ! j corresponds to whole nodes        
            ENDDO
        ENDDO

        ! calculating the convection term for the y-momentum equation
        DO j=2,Ny                ! j=1 and j=Ny+1 are boundary nodes for v
            DO i=2,Nx+1        ! i=1 and i=Nx+2 are boundary nodes for v
                
                
                u_plus = 0.5*(u(i,j) + abs(u(i,j)))
                u_minus = 0.5*(u(i,j) - abs(u(i,j)))

                v_plus = 0.5*(v(i,j) + abs(v(i,j)))
                v_minus = 0.5*(v(i,j) - abs(v(i,j)))
                               
                IF (i == 2) THEN
                    dvdx_minus = (v(i,j) - v(i-1,j))/hx
                    dvdx_plus  = (v(i+1,j) - v(i,j))/hx
                ELSEIF (i == Nx+1) THEN
                    dvdx_minus = (v(i,j) - v(i-1,j))/hx
                    dvdx_plus  = (v(i+1,j) - v(i,j))/hx
                else
                    dvdx_minus = (3*v(i,j) - 4*v(i-1,j) + v(i-2,j))/(2*hx)
                    dvdx_plus  = (-v(i+2,j) + 4*v(i+1,j) - 3*v(i,j))/(2*hx)
                ENDIF                                   

                IF (j == 2) THEN
                    dvdy_minus = (v(i,j) - v(i,j-1))/hy
                    dvdy_plus  = (v(i,j+1) - v(i,j))/hy
                ELSEIF (j == Ny) THEN
                    dvdy_minus = (v(i,j) - v(i,j-1))/hy
                    dvdy_plus  = (v(i,j+1) - v(i,j))/hy
                else
                    dvdy_minus = (3*v(i,j) - 4*v(i,j-1) + v(i,j-2))/(2*hy)
                    dvdy_plus  = (-v(i,j+2) + 4*v(i,j+1) - 3*v(i,j))/(2*hy)
                ENDIF

                ! dummy variable for the convective opperator
                vstar(i,j) = u_plus*dvdx_minus + u_minus*dvdx_plus + &
                             v_plus*dvdy_minus + v_minus*dvdy_plus
                                
                ! i corresponds to half nodes
                ! j corresponds to whole nodes        
            ENDDO
        ENDDO
                
        ! apply viscous BC's
        call BCVELOCITY(Nx,Ny,u,v,BCnorth,BCsouth,BCeast,BCwest,&
                        unorth,usouth,ueast,uwest,&
                        vnorth,vsouth,veast,vwest,&
                        hx,hy,phiLS,contC,contA,contR,contM,2)

        ! finding the intermediate velocity field for the x-momentum equation
        !DO i = 2,Nx
        !        DO j = 2,Ny+1
        DO j=2,Ny+1                ! j=Ny+2 and j=1 are boundary nodes for u
            DO i=2,Nx        ! i=1 and i=Nx+1 are boundary nodes for u
                        
                ! finding the viscous term
                ! dummy variables
                u_E = u(i+1,j)                ! note, i corresponds to half nodes
                u_O = u(i,j)                ! "
                u_W = u(i-1,j)                ! "
                u_N = u(i,j+1)                ! "
                u_S = u(i,j-1)                ! "
    
                v_ne = v(i+1,j)                ! note, j corresponds to half nodes
                v_nw = v(i,j)                ! "
                v_se = v(i+1,j-1)        ! "
                v_sw = v(i,j-1)                ! "
    
                ! u derivatives at the respective locations, n,s,e,w
                u_xe = (u_E - u_O)/hx
                u_xw = (u_O - u_W)/hx
                u_yn = (u_N - u_O)/hy
                u_ys = (u_O - u_S)/hy
    
                !    v derivatives derivatives at the respective locations, n,s
                v_xn = (v_ne - v_nw)/hx
                v_xs = (v_se - v_sw)/hx
    
                ! viscosity at the respective locations, n,s,e,w
                mu_n = 0.25*(mu_one*H(i,j) + mu_two*(1.0-H(i,j)) + &
                             mu_one*H(i+1,j) + mu_two*(1.0-H(i+1,j)) + &
                             mu_one*H(i,j+1) + mu_two*(1.0-H(i,j+1)) + &
                             mu_one*H(i+1,j+1) + mu_two*(1.0-H(i+1,j+1)))
                mu_s = 0.25*(mu_one*H(i,j) + mu_two*(1.0-H(i,j)) + & 
                             mu_one*H(i+1,j) + mu_two*(1.0-H(i+1,j)) + &
                             mu_one*H(i,j-1) + mu_two*(1.0-H(i,j-1)) + &
                           mu_one*H(i+1,j-1) + mu_two*(1.0-H(i+1,j-1)))
                mu_e = mu_one*H(i+1,j) + mu_two*(1.0-H(i+1,j))
                mu_w = mu_one*H(i,j) + mu_two*(1.0-H(i,j))                                                         

                                
                ! the viscous term
                Visc = (2.0*mu_e*u_xe - 2.0*mu_w*u_xw)/hx + &
                       (mu_n*(u_yn + v_xn) - mu_s*(u_ys + v_xs))/hy                        
                       ! note, i corresponds to ihalf
                                
                ! density at the half node
                rho_O = 0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j))) + &
                        0.5*(rho_one*H(i+1,j) + rho_two*(1.0-H(i+1,j)))
                                
                ! surface tension - only for fluid/vapor interface
                IF (trackphase .EQV. .TRUE.) THEN 
                    ! get the curvature 
                    !if (LS_type == 1) then
                        kappa_e = LSCURVATURE(Nx,Ny,i+1,j,hx,hy,phiLS)
                        kappa_w = LSCURVATURE(Nx,Ny,i,j,hx,hy,phiLS)
                        ! body force
                        Fsigma = sigma*0.5*(kappa_e+kappa_w)*(H(i+1,j)-H(i,j))/hx
                    !elseif (LS_type == 3) then
                    !    kappa_e = GALSCURV(Nx,Ny,hx,hy,xwhole,ywhole,phiLS,phi_x,phi_y,phi_xy,prox,i,j,&
                    !                       GA,scal,1)
                    !    !Fsigma = kappa_e                   
                    !    Fsigma = sigma*kappa_e*(H(i+1,j)-H(i,j))/hx
                    !endif
                else
                    Fsigma = 0.0
                ENDIF
                                                              
                ! collect integral terms, IF the node is underneath
                ! the rigid body
                                
                ! ustar on left is intermediate velocity, ustar on right is the conv term
                ustar(i,j) = u(i,j) + deltaT*(-ustar(i,j) + Visc/rho_O - gx - Fsigma/rho_O)        
                                        
            ENDDO
        ENDDO ! end of calculating x-momentum's viscous term

		Fsig = 0.0

        ! finding the intermediate velocity field for the y-momentum equation
        DO j=2,Ny                        ! j=1 and j=Ny+1 are boundary nodes for v
            DO i=2,Nx+1                ! i=1 and i=Nx+2 are boundary nodes for v
                
                
                ! dummy variables
                u_ne = u(i,j+1)                ! note, i corresponds to ihalf
                u_se = u(i,j)                ! "
                u_nw = u(i-1,j+1)        ! "
                u_sw = u(i-1,j)                ! "
    
                v_E = v(i+1,j)                ! note, j corresponds to jhalf
                v_O = v(i,j)                ! "
                v_W = v(i-1,j)                ! "
                v_N = v(i,j+1)                ! "
                v_S = v(i,j-1)                ! "

                ! u derivatives at the respective locations, e,w
                u_ye = (u_ne - u_se)/hy
                u_yw = (u_nw - u_sw)/hy
    
                ! v derivatives at the respective locations, n,s,e,w
                v_yn = (v_N - v_O)/hy
                v_ys = (v_O - v_S)/hy
                v_xe = (v_E - v_O)/hx
                v_xw = (v_O - v_W)/hx
                          
                ! viscosity at the respective locations, n,s,e,w
                mu_n = mu_one*H(i,j+1) + mu_two*(1.0-H(i,j+1))
                mu_s = mu_one*H(i,j) + mu_two*(1.0-H(i,j))
                mu_e = 0.25*(mu_one*H(i,j) + mu_two*(1.0-H(i,j)) + & 
                             mu_one*H(i+1,j) + mu_two*(1.0-H(i+1,j)) + &
                             mu_one*H(i+1,j+1) + mu_two*(1.0-H(i+1,j+1)) + &
                             mu_one*H(i,j+1) + mu_two*(1.0-H(i,j+1)))
                mu_w = 0.25*(mu_one*H(i,j) + mu_two*(1.0-H(i,j)) + &
                             mu_one*H(i-1,j) + mu_two*(1.0-H(i-1,j)) + &
                             mu_one*H(i-1,j+1) + mu_two*(1.0-H(i-1,j+1)) + &
                             mu_one*H(i,j+1) + mu_two*(1.0-H(i,j+1)) )                                                         
    
                Visc = (mu_e*(u_ye + v_xe) - mu_w*(u_yw + v_xw))/hx + &
                       (2.0*mu_n*v_yn - 2.0*mu_s*v_ys)/hy                                        
                       ! note, j corresponds to jhalf
    
                ! density at the half node
                rho_O = 0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j))) + &
                        0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1)))
                        
                ! surface tension
                IF (trackphase .EQV. .TRUE.) THEN
                    ! get the curvature
                    !if (LS_type == 1) then
                        kappa_n = LSCURVATURE(Nx,Ny,i,j+1,hx,hy,phiLS)
                        kappa_s = LSCURVATURE(Nx,Ny,i,j,hx,hy,phiLS)
                        ! body force
                        Fsigma = sigma*0.5*(kappa_n+kappa_s)*(H(i,j+1)-H(i,j))/hy
                    !elseif (LS_type == 3) then
                    !    kappa_n = GALSCURV(Nx,Ny,hx,hy,xwhole,ywhole,phiLS,phi_x,phi_y,phi_xy,prox,i,j,&
                    !                       GA,scal,2)
                    !    !Fsigma = kappa_n
                    !    Fsigma = sigma*kappa_n*(H(i,j+1)-H(i,j))/hy
                    !endif
                else
                    Fsigma = 0.0
                ENDIF
                
                !Fsig(i,j) = -deltaT*Fsigma/rho_O
                if (prox(i,j) > 0 .and. prox(i,j) < width-1) then
                	Fsig(i,j) = kappa_n
            	else
            		Fsig(i,j) = 0.0
            	endif
                !print *, "kappa_n"
                !print *, kappa_n
                !print *, "kappa_s"
                !print *, kappa_s
    
                !rho_O = rho_one
                !beta_O = beta_one
                                
                !print *, "convective vstarij"
                !print *, vstar(i,j)
                                
                ! vstar on left is intermediate velocity, vstar on right is the conv term
                vstar(i,j) = v(i,j) + deltaT*(-vstar(i,j) + Visc/rho_O - gy - Fsigma/rho_O)        
                                
                !print *, "vij"
                !print *, v(i,j)
                                
                !print *, "deltaT"
                !print *, deltaT
                                
                !print *, "Visc"
                !print *, Visc
                                
                !print *, "Fsigma"
                !print *, Fsigma
                                
                !print *, "intermediate vstarij"
                !print *, vstar(i,j)
                                        
            ENDDO
        ENDDO
    
        ! -------------------------------------------------------------------------------------- !
        ! ----------------------------------- PRESSURE SOLVE ----------------------------------- !
        ! -------------------------------------------------------------------------------------- !          
          
        IF (solvertype==1) THEN ! solver 1 only works for single phase
                        
            ! solving for pressure at n+1 implicitly so all the variables on the RHS of
            ! the PPE must be at n+1 or in the case of velocity they must be at star

            ! this is the RHS of the pressure poission equation (PPE)
            PPERHS(1:Nx,1:Ny) = rho_one* &
                                ((ustar(2:Nx+1,2:Ny+1) - ustar(1:Nx,2:Ny+1))/hx + &
                                 (vstar(2:Nx+1,2:Ny+1) - vstar(2:Nx+1,1:Ny))/hy)/deltaT
                                                                
                                                                
            A=0.0                ! X(0)
            B=Lx                ! X(XL)
            M=Nx
            MBDCND=3        ! 3=DERIVATIVE SPECIFIED AT BOTH ENDS
            BDA=0.0
            BDB=0.0
            C=0.0                ! Y(0)
            D=Ly                ! Y(YL)
            N=Ny
            NBDCND=3        ! 3=DERIVATIVE SPECIFIED AT BOTH ENDS, 2=reference p=0 at y=0
            BDC=0.0
            BDD=0.0
            ELMBDA=0.0        ! LAMBDA COEFFICIENT IN HELMHOLTZ EQ. (ALWAYS ZERO FOR THIS PPE)
            IDIMF=Nx

            ! FISHPAK ROUTINES FOR PPE SOLUTION
            !CALL HSTCRT (A,B,M,MBDCND,BDA,BDB,C,D,N,NBDCND,BDC,BDD, &
            !             ELMBDA,PPERHS,IDIMF,PERTRB,IERROR,WORK)

            IF (IERROR.NE.0) THEN
                WRITE(*,*) 'ERROR IN FISHPACK IN HSTCRT, IERROR= ',IERROR
                !PAUSE
            ENDIF

            IF (WORK(1).GE.100000.0) THEN
                WRITE(*,*) 'WORK ARRAY IN FISHPACK UNDERDIMENSIONED'
                !PAUSE
            ENDIF

            P(2:nx+1,2:ny+1) = PPERHS(1:nx,1:ny)
                        
            ! Neumann BC's for pressure
            P(1:nx+2,1) = P(1:nx+2,2)
            P(1:nx+2,ny+2) = P(1:nx+2,ny+1)
            P(1,1:ny+2) = P(2,1:ny+2)
            P(nx+2,1:ny+2) = P(nx+1,1:ny+2)
            ! ---- END OF PPE SOLUTION ----
                
        ENDIF
                
        IF (solvertype==2) THEN
                        
            !--- Specifing the pressure BC type ---
            ! north
            IF (BCnorth==3 .or. BCnorth==4) THEN
                ! Dirclet Pressure BC
                BCPnorth = 2
            else
                ! Neumann Pressure BC
                BCPnorth = 1
            ENDIF
            ! south
            IF (BCsouth==3 .or. BCsouth==4) THEN
                ! Dirclet Pressure BC
                BCPsouth = 2
            else
                ! Neumann Pressure BC
                BCPsouth = 1
            ENDIF 
            ! east
            IF (BCeast==3 .or. BCeast==4) THEN
                ! Dirclet Pressure BC
                BCPeast = 2
            else
                ! Neumann Pressure BC
                BCPeast = 1
            ENDIF
            ! west
            IF (BCwest==3 .or. BCwest==4) THEN
                ! Dirclet Pressure BC
                BCPwest = 2
            else
                ! Neumann Pressure BC
                BCPwest = 1
            ENDIF
                        
                        
            ! this is the RHS of the pressure poission equation (PPE) 
            DO i=1,Nx  ! looping over all interior nodes
                DO j=1,Ny
                                        
                    ! ---- PHASE TRANSITION ----
                    ! get the normal vector
                    ! call LSNORMAL(Nx,Ny,i+1,j+1,hx,hy,phiLS,normal)
        
                    !print *, i
                    !print *, j

                    !---- PPE RHS ----
                    deltau = ustar(i+1,j+1)-ustar(i,j+1)
                    deltav = vstar(i+1,j+1)-vstar(i+1,j)
                    !print *, "delta_u"
                    !print *, deltau
                    !print *, "delta_v"
                    !print *, deltav
                    ! JDS 6/17/2008 Removed masstrans term from PPERHS
                    PPERHS(i,j) = ( deltau/hx + deltav/hy )/deltaT
                    ! note, i in ustar corresponds to ihalf
                    ! note, j in vstar corresponds to jhalf
                    ! note, delstar is at whole nodes
                    !print *, "hx"
                    !print *, hx
                    !print *, "hy"
                    !print *, hy
                    !print *, "deltaT"
                    !print *, deltaT
                    !print *, "PPERHSij"
                    !print *, PPERHS(i,j)
                    !call ArrtoVec(Nx,Ny,PPERHS,junk)
                    !junk2 = norm(Nx,Ny,junk)
                    !print *, "norm of PPERHS"
                    !print *, junk2
                ENDDO
            ENDDO ! end of the PPERHS                        
                        
            Pold = P
            call PCG(Nx,Ny,hx,hy,P(2:Nx+1,2:Ny+1),PPERHS,H,s_H,&
                     rho_one,rho_two,rho_three,&
                     tolerance,itermax,iterprint,precond,&
                     BCPnorth,BCPsouth,BCPeast,BCPwest,Pnorth,Psouth,Peast,Pwest)
                        

            ! Neumann BC's for pressure
            IF (BCPsouth==1) THEN
                P(1:nx+2,1) = P(1:nx+2,2)
            ENDIF
            IF (BCPnorth==1) THEN
                P(1:nx+2,ny+2) = P(1:nx+2,ny+1)
            ENDIF
            IF (BCPwest==1) THEN
                P(1,1:ny+2) = P(2,1:ny+2)
            ENDIF
            IF (BCPeast==1) THEN
                P(nx+2,1:ny+2) = P(nx+1,1:ny+2)
            ENDIF
                        
            ! Dirclet BC's for pressure
            IF (BCPsouth==2) THEN
                P(1:nx+2,1) = Psouth
            ENDIF
            IF (BCPnorth==2) THEN
                P(1:nx+2,ny+2) = Pnorth
            ENDIF
            IF (BCPwest==2) THEN
                P(1,1:ny+2) = Pwest
            ENDIF
            IF (BCPeast==2) THEN
                P(nx+2,1:ny+2) = Peast
            ENDIF
            ! ---- END OF PPE SOLUTION ----
        ENDIF
  
        ! ---- SOR -----
        IF (solvertype==3) THEN 
                        
            ! this is the RHS of the pressure poission equation (PPE) rho(2:Nx+1,2:Ny+1)
            PPERHS(1:Nx,1:Ny) = ((ustar(2:Nx+1,2:Ny+1) - ustar(1:Nx,2:Ny+1))/hx + &
                                 (vstar(2:Nx+1,2:Ny+1) - vstar(2:Nx+1,1:Ny))/hy)/deltaT
                                 ! note, i in ustar corresponds to ihalf
                                 ! note, j in vstar corresponds to jhalf
                                 ! note, delstar is at whole nodes
                        
            DO iter = 1,itermax

                ! ---- APPLY BC's ----
                ! Neumann BC's for pressure
                P(1:nx+2,1) = P(1:nx+2,2)
                P(1:nx+2,ny+2) = P(1:nx+2,ny+1)
                P(1,1:ny+2) = P(2,1:ny+2)
                P(nx+2,1:ny+2) = P(nx+1,1:ny+2)
    
                ! corners (not used in derivative stencil)
                P(1,1) = 0.5*(P(2,1)+P(1,2))
                P(1,Ny+2)=0.5*(P(2,Ny+2)+P(1,Ny-1))
                P(Nx+2,1)=0.5*(P(Nx+1,1)+P(Nx+2,2))
                P(Nx+2,Ny+2)=0.5*(P(Nx+1,Ny+2)+P(Nx+2,Ny+1))
                                                   
                P(floor(real(Nx)/2.0)+1,1) = 0.0

                Pavg = 0.0        ! initializing back to zero
                maxPRes = 0.0 ! initializing back to zero
                DO i=1,Nx
                    DO j=1,Ny
                        ! density 
                            rho_e = 0.5*(rho_one*H(i+2,j+1) + rho_two*(1.0-H(i+2,j+1))) + &
                                    0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1)))
                            rho_w = 0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1))) + &
                                    0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1)))
                            rho_n = 0.5*(rho_one*H(i+1,j+2) + rho_two*(1.0-H(i+1,j+2))) + &
                                    0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1)))
                            rho_s = 0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1))) + &
                                    0.5*(rho_one*H(i+1,j) + rho_two*(1.0-H(i+1,j)))
                                                
                            aw = 1.0/(rho_w*hx**2.0) !1.0/hx**2.0
                            ae = 1.0/(rho_e*hx**2.0) !1.0/hx**2.0
                            an = 1.0/(rho_n*hy**2.0) !1.0/hx**2.0
                            as = 1.0/(rho_s*hy**2.0) !1.0/hx**2.0
                            aO = -aw-ae-an-as                         !-4.0/hx**2.0

                            PRes = PPERHS(i,j)-ae*P(i+2,j+1)-aw*P(i,j+1) - &
                                   an*P(i+1,j+2)-as*P(i+1,j)-ao*P(i+1,j+1)
                            P(i+1,j+1) = P(i+1,j+1) + omega/aO*PRes
                                                
                            maxPRes = max(abs(PRes),maxPRes)  ! largest error
                            Pavg = P(i+1,j+1) + Pavg        ! summing up the pressure
                    ENDDO
                ENDDO
                Pavg = Pavg/(Nx*Ny)

                ! for Neumann BC's
                !DO i=2,Nx+2
                !        DO j=2,Ny+2
                !                P(i,j) = P(i,j)-Pavg  ! keeping the solution bounded
                !        ENDDO
                !ENDDO
    
                ! ---- check for convergence ----
                IF (iter>1) THEN
                    ! convergence criteria
                    IF (maxPRes<=tolerance) THEN
                        exit
                    ENDIF
                ENDIF

                IF (mod(iter,iterprint)==0) THEN
                    print *, 'Pressure Residual =', maxPRes
                ENDIF

                IF (iter+1==itermax) THEN
                    print *, 'ERROR-> Pressure Failed to Converge!'
                ENDIF
                                        
                PavgOld = Pavg

            ENDDO  ! end of iteration
                
            ! ---- END OF PPE SOLUTION ----
        ENDIF

        !----------------------------------------------------------------------------!
        !-------------------------------- NEW VELOCITY ------------------------------!
        !----------------------------------------------------------------------------!   
        
        DO i=2,Nx
            DO j=2,Ny+1
                ! density at the half node
                rho_O = 0.5*(rho_one*H(i+1,j) + rho_two*(1.0-H(i+1,j))) + &
                        0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))

                ! Hold old variables
                u_old(i,j) = u(i,j)
                !u(i,j) = - 2.0*cos(pi*time/8.0)*(sin(pi*x(i))**2.0)*sin(pi*ywhole(j))*cos(pi*ywhole(j))

                ! new velocity, n+1 time level
                u(i,j) = ustar(i,j) - deltaT/rho_O*(P(i+1,j)-P(i,j))/hx  ! note, i corresponds to ihalf
           	
            ENDDO
        ENDDO        
        
        DO i=2,Nx+1
            DO j=2,Ny
                ! density at the half node
                rho_O = 0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1))) + &
                        0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))
                                
                ! Hold old variable
                v_old(i,j) = v(i,j)
                !v(i,j) = 2.0*cos(pi*time/8.0)*(sin(pi*y(j))**2.0)*sin(pi*xwhole(i))*cos(pi*xwhole(i))
                
                ! new velocity, n+1 time level
                v(i,j) = vstar(i,j) - deltaT/rho_O*(P(i,j+1)-P(i,j))/hy  ! note, j corresponds to jhalf
            	
            ENDDO
        ENDDO
                
        !----------------------------------------------------------------------------!
        !----------------------------- EVOLVE LEVEL SET -----------------------------!
        !----------------------------------------------------------------------------!   
        
        IF ((trackphase .EQV. .TRUE.)) THEN
            ! All variables,(phiLS, u, v, P, etc) are at time level n, 
            ! but after the new level-set field is calculated the properties 
            ! will be at the n+1 time level.
            
            IF (LS_type == 1) THEN
                

                !---- STANDARD LEVEL SET METHOD ----!
                !---- USING RK4 ----
                ! PDE => d(phiLS)/d(t) = -V*grad(phiLS)
                DO i=2,Nx+1
                    DO j=2,Ny+1
						phiLSn(i,j) = phiLS(i,j)                                                       
                    ENDDO
                ENDDO                
            	DO step=1,4                
                    ! getting the RHS                
                    DO i=2,Nx+1
                        DO j=2,Ny+1
                            ! calculate spatial derivatives
                            IF (i<4 .or. j<4) THEN
                                CALL ENODIFF(Nx,Ny,i,j,hx,hy,u,v,phiLS,DiffX,DiffY)
                            ELSEIF (i>Nx-2 .or. j>Ny-2) THEN
                                CALL ENODIFF(Nx,Ny,i,j,hx,hy,u,v,phiLS,DiffX,DiffY)
                            else
                                !CALL ENODIFF(Nx,Ny,i,j,hx,hy,u,v,phiLS,DiffX,DiffY)
                                                          
                                ! use WENO on interior (up to 5th order accurate)
                                CALL WENODIFF(Nx,Ny,i,j,hx,hy,u,v,phiLS,DiffX,DiffY)
                            ENDIF
                                                
                                                
                            ! ---- PHASE TRANSITION ----
                            ! get the normal vector
                            !call LSNORMAL(Nx,Ny,i,j,hx,hy,phiLS,normal)
                                                  
                            ! density at the O node
                            !rho_O = rho_one*H(i,j) + rho_two*(1.0-H(i,j))
  
                            !  JDS 7/16/2008 Took umflux and vmflux out of equation                                                
                            ! phiLSstar is a dummy variable for advection
                            phiLSstar(i,j) = -( 0.5*(u(i,j)+u(i-1,j)))*DiffX - &
                                              ( 0.5*(v(i,j)+v(i,j-1)))*DiffY                

                        ENDDO
                    ENDDO
                                
                    ! get next value
                    DO i=2,Nx+1
                        DO j=2,Ny+1
                            ! note, Coef is positive

                            phiLS(i,j) = phiLSn(i,j) + Coef(step)*deltaT*phiLSstar(i,j)      
                                                        
                        ENDDO
                    ENDDO

                    ! apply boundary conditions (extrapolation)
                    CALL BCLEVELSET(Nx,Ny,hx,hy,ywhole,phiLS,&
                              BCphisouth,BCphinorth,BCphiwest,BCphieast,d_phiLS)
                              
                ENDDO
                ! phiLS is now at n+1
                                
                !---- REINITIALIZATION ----
                ! reinitialize the new LS field to a signed distance function
                !CALL LSREIN(Nx,Ny,hx,hy,phiLS,H,&
                !            BCphisouth,BCphinorth,BCphiwest,BCphieast,ywhole,d_phiLS)
            


            ELSEIF (LS_type == 2) THEN      !*********** CLSVOF phase change doesn't work yet ************
                
                uINT = u
                vINT = v
                
                !---------------------START OF CLSVOF---------------------
                !call CLSVOF(Nx,Ny,u,v,uINT,vINT,phiLS,F,Avof,Bvof,Cvof,hx,hy,deltaT,&
                !            xwhole,ywhole,xdirgo,ydirgo,tolerance,rho_one,rho_two)  ! rho_one, rho_two  advect density version
                !call CLSVOF(Nx,Ny,uINT,vINT,phiLS,F,Avof,Bvof,Cvof,hx,hy,deltaT,&
                !           xwhole,ywhole,xdirgo,ydirgo,tolerance,rho_one,rho_two)

                ! get new order
                xdirgo=xdirgo+1
                ydirgo=ydirgo+1
                IF (xdirgo>2) THEN
                    xdirgo=1
                ELSEIF (ydirgo>2) THEN
                    ydirgo=1
                ENDIF
                ! boundary conditions
                Avof(:,1)=Avof(:,2)
                Avof(:,Ny+2)=Avof(:,Ny+1)
                Avof(1,:)=Avof(2,:)
                Avof(Nx+2,:)=Avof(Nx+1,:)
                Bvof(:,1)=Bvof(:,2)
                Bvof(:,Ny+2)=Bvof(:,Ny+1)
                Bvof(1,:)=Bvof(2,:)
                Bvof(Nx+2,:)=Bvof(Nx+1,:)
                Cvof(:,1)=Cvof(:,2)
                Cvof(:,Ny+2)=Cvof(:,Ny+1)
                Cvof(1,:)=Cvof(2,:)
                Cvof(Nx+2,:)=Cvof(Nx+1,:)

                !---------------------END OF CLSVOF---------------------
        

            ELSEIF (LS_type == 3) THEN
              
            !---- GRADIENT-AUGMENTED LEVEL SET METHOD ----!         
    
				CALL VOLUMESOLVE(Nx,Ny,hx,hy,xwhole,ywhole,GA,&
								 scal,phiLS,phi_x,phi_y,phi_xy,prox,volume)
								 
				if (1.0-abs(volume/volume0) > 0.001) then
					CALL VOLUMEFIX(Nx,Ny,hx,hy,phiLS,prox,volume,volume0)
				endif
								 
								 
				print *, "VOLUME FRACTION = ",volume/volume0
        
                call BCVELOCITY(Nx,Ny,u,v,BCnorth,BCsouth,BCeast,BCwest,&
                                unorth,usouth,ueast,uwest,&
                                vnorth,vsouth,veast,vwest,&
                                hx,hy,phiLS,contC,contA,contR,contM,1)
           		
           		CALL GALS(Nx,Ny,hx,hy,deltaT,xwhole,ywhole,GA,scal,phiLS,phi_x, &
           		          phi_y,phi_xy,prox,x,y,u,v,u_old,v_old,volume,volume0,time,width)
            
            ENDIF
                                    
            !---------------------END OF LSM---------------------


            ! ------------------------------------------------ !
            ! ----------- NEW HEAVISIDE FUNCTION ------------- !
            ! ------------------------------------------------ !              
                                       

             

            IF (LS_type == 3) THEN
                DO i=1,Nx+2
                    DO j=1,Ny+2
                         ! Heavieside function
                         H(i,j) = LSHEAVY(Nx,Ny,i,j,hx,hy,phiLS)
                         !H(i,j) = LSHEAVY2(Nx,Ny,i,j,hx,hy,phiLS,phi_x,phi_y)
                    ENDDO
                ENDDO
            ELSE
                volume = 0.0
                DO i=1,Nx+2
                    DO j=1,Ny+2
                        ! Heavieside function
                        H(i,j) = LSHEAVY(Nx,Ny,i,j,hx,hy,phiLS)
                        volume = volume + H(i,j)*hx*hy
                    ENDDO
                ENDDO
                print *, "VOLUME FRACTION = ",volume/volume0
            ENDIF
        
        ENDIF  ! end of trackphase


        ! ------------------------------------------------ !
        ! -------------- ADVANCE THE WEDGE --------------- !
        ! ------------------------------------------------ !                         
          
        IF (wedge .EQV. .TRUE.) THEN      
            
            IF (time < wedgeT) THEN      
                n12_w_vy = (2*pi)*(0.081)*(freq)*sin(2*pi*(freq)*(time + 0.5*deltaT))
                
            else
                n12_w_vy = 0
            ENDIF

            ! Calculate d(n+1)                
            !n1_rb_dx = n_rb_dx + deltaT*n12_rb_vx 
            n1_w_dy = n_w_dy + deltaT*n12_w_vy
                
                
            DO i=1,Nx+2
                DO j=1,Ny+2
                                
                    IF (ywhole(j) > 3.82402*xwhole(i) - 0.340938 + n1_w_dy) THEN
                        w_phiLS(i,j) = -(-xwhole(i) + 0.114)
                    else
                        w_phiLS(i,j) = -(-0.46345*(ywhole(j)-n1_w_dy) + 0.88612*xwhole(i) - 0.05699)
                    ENDIF
                
                    ! ---- NEW HEAVIESIDE FCN ----
                    ! finding the new heavieside function at each node
                    w_H(i,j) = LSHEAVY(Nx,Ny,i,j,hx,hy,w_phiLS)                                
                                                
                ENDDO
            ENDDO        

            DO i=2,Nx+1
                DO j=2,Ny
                        
                    ! density at the half node
                    rho_O = 0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1))) + &
                            0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))


                    IF (0.5*(w_H(i,j)+w_H(i,j+1)) > 0) THEN
                                        
                        ! density at the whole node
                        !rho_O = rho_one*H(i,j) + rho_two*(1.0-H(i,j))                        
                        !IF (s_phiLS(i,j) >= 0) THEN
                                              
                        u_plus = 0.5*(u(i,j) + abs(u(i,j)))
                        u_minus = 0.5*(u(i,j) - abs(u(i,j)))

                        v_plus = 0.5*(v(i,j) + abs(v(i,j)))
                        v_minus = 0.5*(v(i,j) - abs(v(i,j)))
                               
                        IF (i == 2) THEN
                            dvdx_minus = (v(i,j) - v(i-1,j))/hx
                            dvdx_plus  = (v(i+1,j) - v(i,j))/hx
                        ELSEIF (i == Nx+1) THEN
                            dvdx_minus = (v(i,j) - v(i-1,j))/hx
                            dvdx_plus  = (v(i+1,j) - v(i,j))/hx
                        else
                            dvdx_minus = (3*v(i,j) - 4*v(i-1,j) + v(i-2,j))/(2*hx)
                            dvdx_plus  = (-v(i+2,j) + 4*v(i+1,j) - 3*v(i,j))/(2*hx)
                        ENDIF                                   

                        IF (j == 2) THEN
                            dvdy_minus = (v(i,j) - v(i,j-1))/hy
                            dvdy_plus  = (v(i,j+1) - v(i,j))/hy
                        ELSEIF (j == Ny) THEN
                            dvdy_minus = (v(i,j) - v(i,j-1))/hy
                            dvdy_plus  = (v(i,j+1) - v(i,j))/hy
                        else
                            dvdy_minus = (3*v(i,j) - 4*v(i,j-1) + v(i,j-2))/(2*hy)
                            dvdy_plus  = (-v(i,j+2) + 4*v(i,j+1) - 3*v(i,j))/(2*hy)
                        ENDIF

                        ! dummy variable for the convective operator
                        vstar_int = u_plus*dvdx_minus + u_minus*dvdx_plus + &
                                    v_plus*dvdy_minus + v_minus*dvdy_plus

                        IF ((wedge .EQV. .TRUE.) .and. (0.5*(w_H(i,j)+w_H(i,j+1)) > 0)) THEN                                      
                   
                            n1_w_inert = n1_w_inert + &
                                         0.5*(w_H(i,j)+w_H(i,j+1))*&
                                         rho_O/DeltaT*(v(i,j) - v_old(i,j))*&
                                         hx*hy
                                                                                        
                            n1_w_conv = n1_w_conv + &
                                        0.5*(w_H(i,j)+w_H(i,j+1))*&
                                        rho_O*vstar_int*hx*hy
                                
                            n1_w_grav = n1_w_grav + 0.5*&
                                        (w_H(i,j)+w_H(i,j+1))*rho_O*hx*hy*gy

                        ENDIF                                
                    ENDIF
                ENDDO
            ENDDO

            wedge_force = n_w_inert + n_w_conv + n_w_grav
            wedge_work = wedge_work - (wedge_force * n1_w_vy * deltaT)
            
            IF (time < wedgeT) THEN
                n1_w_vy = (2*pi)*(0.081)*(freq)*sin(2*pi*(freq)*time)            
            else
                n1_w_vy = 0      
            ENDIF

        ENDIF
          
        ! Fluid Energy Contribution
        IF (rigidbody .EQV. .FALSE.) THEN
            
            total_energy = 0.0
            lin_momx = 0.0
            lin_momy = 0.0
    
            DO i=2,Nx+1
                DO j=2,Ny+1
    
                    IF (i < Nx+1) THEN
                        rho_O = 0.5*(rho_one*H(i+1,j) + rho_two*(1.0-H(i+1,j)))+&
                                0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))
                                                 
                        total_energy = total_energy + 0.5*rho_O*(u(i,j)**2)*hx*hy  ! kinetic (in x)
                        lin_momx = lin_momx + rho_O*u(i,j)*hx*hy
                    ENDIF
              
                    IF (j < Ny+1) THEN
                        rho_O = 0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1)))+&
                                0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))                          

                        total_energy = total_energy + 0.5*rho_O*(v(i,j)**2)*hx*hy  ! kinetic (in y)
                        lin_momy = lin_momy + rho_O*v(i,j)*hx*hy
                    ENDIF
 
                    IF ((i < Nx+1) .AND. (j < Ny+1)) THEN 
                        rho_O = rho_one*H(i,j) + rho_two*(1.0-H(i,j))

                        total_energy = total_energy + &                        ! Pressure and gravitational potential
                                       ( -P(i,j) + rho_O*gy*(j-1)*hy ) * &
                                       hx*hy
                    ENDIF

                ENDDO
            ENDDO
        ENDIF

        !-------------------------------------------!
        !---- RECONSTRUCT RIGID BODY LEVEL SETS ----!
        !-------------------------------------------!  
        
        ! Reconstruct the level set
        IF (rigidbody .EQV. .TRUE.) THEN
            
            ! Level sets

            IF (rb_shape == 1) THEN
            
                th = 0.0 + th_change 
  
                xc = b_rbLS + n1_rb_dx
                yc = c_rbLS + n1_rb_dy                        
                            
                DO i=1,Nx+2
                    DO j=1,Ny+2
                
                        ! Circle
                        s_phiLS(i,j) = a_rbLS*sqrt( (xwhole(i)-(b_rbLS + n1_rb_dx))**2.0 + &
                                       (ywhole(j)-(c_rbLS + n1_rb_dy))**2.0 ) + d_rbLS             
                
                    ENDDO
                ENDDO

            ELSEIF (rb_shape == 2) THEN  
            
                xn = xo + n1_rb_dx
                yn = yo + n1_rb_dy
                xc = xn
                yc = yn

                th = pi/tho + th_change   
                   
                DO i=1,Nx+2
                    DO j=1,Ny+2                      
                                                        
                        s_phiLS(i,j) = -100*&
                                       (bellipse**2*(xwhole(i)*cos(th) + ywhole(j)*sin(th) - &
                                                    (xn*cos(th) + yn*sin(th)))**2 + &
                                        aellipse**2*(-xwhole(i)*sin(th) + ywhole(j)*cos(th) - &
                                                    (-xn*sin(th) + yn*cos(th)))**2 - &
                                        aellipse**2*bellipse**2)        
                    ENDDO
                ENDDO

            ELSEIF (rb_shape == 3) THEN
                        
                xc = c_xo + n1_rb_dx
                yc = c_yo + n1_rb_dy
                nx1 = x1 + n1_rb_dx
                ny1 = y1 + n1_rb_dy
                nx2 = x2 + n1_rb_dx
                ny2 = y2 + n1_rb_dy
                nx3 = x3 + n1_rb_dx
                ny3 = y3 + n1_rb_dy
                nx4 = x4 + n1_rb_dx
                ny4 = y4 + n1_rb_dy
                ncgx = cgx + n1_rb_dx
                ncgy = cgy + n1_rb_dy
                                  
                th = c_tho + th_change
    
                cost = cos(th)
                sint = sin(th)
                rot_mat(1,1) = cost
                rot_mat(1,2) = -sint
                rot_mat(2,1) = sint
                rot_mat(2,2) = cost
                                          
                ! x1-4,y-4 are ordered counterclockwise
                                
                DO r=1,4
                                
                    IF (r == 1) THEN
                                        
                        ! Region 1
                        x_v(1) = xc;
                        x_v(2) = nx1;
                        x_v(3) = nx2;
                                                      
                        y_v(1) = yc;
                        y_v(2) = ny1;
                        y_v(3) = ny2;
                                                
                    ELSEIF (r == 2) THEN
                                        
                        ! Region 2
                        x_v(1) = xc;
                        x_v(2) = nx2;
                        x_v(3) = nx3;
                        
                        y_v(1) = yc;
                        y_v(2) = ny2;
                        y_v(3) = ny3;
                                                
                    ELSEIF (r == 3) THEN
                                        
                        ! Region 3
                        x_v(1) = xc;
                        x_v(2) = nx3;
                        x_v(3) = nx4;
                                                
                        y_v(1) = yc;
                        y_v(2) = ny3;
                        y_v(3) = ny4;
                                                
                    ELSEIF (r == 4) THEN
                                        
                        ! Region 4
                        x_v(1) = xc;
                        x_v(2) = nx4;
                        x_v(3) = nx1;
                                                
                        y_v(1) = yc;
                        y_v(2) = ny4;
                        y_v(3) = ny1;
                                                
                    ENDIF
                                        
                    nx_v(1) = rot_mat(1,1)*(xc-ncgx) + rot_mat(1,2)*(yc-ncgy)
                    ny_v(1) = rot_mat(2,1)*(xc-ncgx) + rot_mat(2,2)*(yc-ncgy)                        
                    nx_v(2) = rot_mat(1,1)*(x_v(2)-ncgx) + rot_mat(1,2)*(y_v(2)-ncgy)
                    ny_v(2) = rot_mat(2,1)*(x_v(2)-ncgx) + rot_mat(2,2)*(y_v(2)-ncgy)
                    nx_v(3) = rot_mat(1,1)*(x_v(3)-ncgx) + rot_mat(1,2)*(y_v(3)-ncgy)
                    ny_v(3) = rot_mat(2,1)*(x_v(3)-ncgx) + rot_mat(2,2)*(y_v(3)-ncgy)
    
                    x_v(1) = nx_v(1) + ncgx
                    y_v(1) = ny_v(1) + ncgy                        
                    x_v(2) = nx_v(2) + ncgx
                    y_v(2) = ny_v(2) + ncgy
                    x_v(3) = nx_v(3) + ncgx
                    y_v(3) = ny_v(3) + ncgy
                                                    
                    z_v(1) = 1.0;
                    z_v(2) = 0.0;
                    z_v(3) = 0.0;
                            
                    Det = x_v(1)*(y_v(2)*z_v(3)-y_v(3)*z_v(2)) - &
                          y_v(1)*(x_v(2)*z_v(3)-x_v(3)*z_v(2)) + &
                          z_v(1)*(x_v(2)*y_v(3)-x_v(3)*y_v(2))
                                                          
                    a_det = y_v(2)*z_v(3)-y_v(3)*z_v(2) - &
                            y_v(1)*(z_v(3)-z_v(2)) + &
                            z_v(1)*(y_v(3)-y_v(2))
                                                          
                    b_det = x_v(1)*(z_v(3)-z_v(2)) - &
                            x_v(2)*z_v(3)-x_v(3)*z_v(2) + &
                            z_v(1)*(x_v(2)-x_v(3))
                                                          
                    c_det = x_v(1)*(y_v(2)-y_v(3)) - &
                            y_v(1)*(x_v(2)-x_v(3)) + &
                            x_v(2)*y_v(3)-x_v(3)*y_v(2)
                                                                
                    a_rbLS = a_det/Det
                    b_rbLS = b_det/Det
                    c_rbLS = c_det/Det
                                                
                    DO i=1,Nx+2
                        DO j=1,Ny+2
                                        
                            ! Check on the location of the grid point.
                            ! Is it in the correct region?
                                                    
                            pt(1) = xwhole(i)
                            pt(2) = ywhole(j)
            
                            t(1,1) = x_v(1);
                            t(2,1) = y_v(1);
                            t(1,2) = x_v(2);
                            t(2,2) = y_v(2);
                            t(1,3) = x_v(3);
                            t(2,3) = y_v(3);
                
                            DO m = 1,3,2
                                                
                                k = mod(m,3) + 1
                                                        
                                temp = ( pt(1) - t(1,m) ) * ( t(2,k) - t(2,m) ) - &
                                       ( pt(2) - t(2,m) ) * ( t(1,k) - t(1,m) )
                                                                   
                                IF (0.0 < temp) THEN
                                    inside = .FALSE.
                                    exit
                                ENDIF
                                inside = .TRUE.
                            ENDDO
                                                
                            IF (inside .EQV. .TRUE.) THEN
                                                
                                s_phiLS(i,j) = (-a_rbLS*xwhole(i) - b_rbLS*ywhole(j) + 1)/c_rbLS
                                                                                   
                            ENDIF
                        ENDDO                
                    ENDDO
                ENDDO
            ENDIF
   
            ! Save RB location & orientation
            rb_data(60) = xc
            rb_data(61) = yc
            rb_data(62) = th     
                                        
            ! ---- NEW HEAVIESIDE FCN ----
            ! finding the new heavieside function at each node
            DO i=1,Nx+2
              DO j=1,Ny+2
                s_H(i,j) = LSHEAVY(Nx,Ny,i,j,hx,hy,s_phiLS)
              ENDDO
            ENDDO
                
                
            ! ------------------------------------- !
            ! ---- CALCULATE RIGID BODY FORCES ---- !
            ! ------------------------------------- !   
                      
            ! Collect some integrals
            
            total_energy = 0.0
            lin_momx = 0.0
            lin_momy = 0.0

            DO i=2,Nx
                DO j=2,Ny+1
                        
                    ! density at the half node
                    rho_O = 0.5*(rho_one*H(i+1,j) + rho_two*(1.0-H(i+1,j))) + &
                            0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))

                    total_energy = total_energy + 0.5*rho_O*(u(i,j)**2)*hx*hy  ! kinetic (in x)                                
                    lin_momx = lin_momx + rho_O*u(i,j)*hx*hy               
                     
                    IF (0.5*(s_H(i,j)+s_H(i+1,j)) > 0)  THEN
                                              
                        u_plus = 0.5*(u(i,j) + abs(u(i,j)))
                        u_minus = 0.5*(u(i,j) - abs(u(i,j)))

                        v_plus = 0.5*(v(i,j) + abs(v(i,j)))
                        v_minus = 0.5*(v(i,j) - abs(v(i,j)))
                                                        
                        IF (i == 2) THEN
                            dudx_minus = (u(i,j) - u(i-1,j))/hx
                            dudx_plus  = (u(i+1,j) - u(i,j))/hx
                        ELSEIF (i == Nx) THEN
                            dudx_minus = (u(i,j) - u(i-1,j))/hx
                            dudx_plus  = (u(i+1,j) - u(i,j))/hx
                        else
                            dudx_minus = (3*u(i,j) - 4*u(i-1,j) + u(i-2,j))/(2*hx)
                            dudx_plus  = (-u(i+2,j) + 4*u(i+1,j) - 3*u(i,j))/(2*hx)
                        ENDIF                                   

                        IF (j == 2) THEN
                            dudy_minus = (u(i,j) - u(i,j-1))/hy
                            dudy_plus  = (u(i,j+1) - u(i,j))/hy
                        ELSEIF (j == Ny+1) THEN
                            dudy_minus = (u(i,j) - u(i,j-1))/hy
                            dudy_plus  = (u(i,j+1) - u(i,j))/hy
                        else
                            dudy_minus = (3*u(i,j) - 4*u(i,j-1) + u(i,j-2))/(2*hy)
                            dudy_plus  = (-u(i,j+2) + 4*u(i,j+1) - 3*u(i,j))/(2*hy)
                        ENDIF
 
                        ! dummy variable for the convective operator
                        ustar_int = u_plus*dudx_minus + u_minus*dudx_plus + &
                                    v_plus*dudy_minus + v_minus*dudy_plus
                
                        n1_xint_inert = n1_xint_inert +&
                                        0.5*(s_H(i,j)+s_H(i+1,j))*&
                                        rho_O/DeltaT*(u(i,j) - u_old(i,j))*&
                                        hx*hy
                                                
                        n1_xint_conv = n1_xint_conv +&
                                        0.5*(s_H(i,j)+s_H(i+1,j))*&
                                        rho_O*ustar_int*hx*hy
                                
                        n1_xint_grav = n1_xint_grav + &
                                        0.5*(s_H(i,j)+s_H(i+1,j))*&
                                        rho_O*hx*hy*gx
                                                
                        ! Rotations
                        ! rs is the s-distance from the location at i,j to the centroid of the object.
                        ry = ywhole(j) - ncgy
                                                
                        r1_int_inert = r1_int_inert - &
                                       ry*0.5*(s_H(i,j)+s_H(i+1,j))*&
                                       rho_O/deltaT*(u(i,j) - u_old(i,j))*&
                                       hx*hy
                                                
                        r1_int_conv = r1_int_conv - &
                                      ry*0.5*(s_H(i,j)+s_H(i+1,j))*&
                                      rho_O*ustar_int*hx*hy
                                
                        r1_int_grav = r1_int_grav - &
                                      ry*0.5*(s_H(i,j)+s_H(i+1,j))*&
                                      rho_O*hx*hy*gx
                                
                    ENDIF
                ENDDO
            ENDDO

            DO i=2,Nx+1
                DO j=2,Ny
                        
                    ! density at the half node
                    rho_O = 0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1))) + &
                            0.5*(rho_one*H(i,j) + rho_two*(1.0-H(i,j)))
                                                
                    total_energy = total_energy + 0.5*rho_O*(v(i,j)**2)*hx*hy ! kinetic (in y)
                    lin_momy = lin_momy + rho_O*v(i,j)*hx*hy

                    IF (0.5*(s_H(i,j)+s_H(i,j+1)) > 0)then
                                        
                        ! density at the whole node
                        !rho_O = rho_one*H(i,j) + rho_two*(1.0-H(i,j))                        
                        !IF (s_phiLS(i,j) >= 0) THEN
                                        
                        u_plus = 0.5*(u(i,j) + abs(u(i,j)))
                        u_minus = 0.5*(u(i,j) - abs(u(i,j)))

                        v_plus = 0.5*(v(i,j) + abs(v(i,j)))
                        v_minus = 0.5*(v(i,j) - abs(v(i,j)))
                               
                        IF (i == 2) THEN
                            dvdx_minus = (v(i,j) - v(i-1,j))/hx
                            dvdx_plus  = (v(i+1,j) - v(i,j))/hx
                        ELSEIF (i == Nx+1) THEN
                            dvdx_minus = (v(i,j) - v(i-1,j))/hx
                            dvdx_plus  = (v(i+1,j) - v(i,j))/hx
                        else
                            dvdx_minus = (3*v(i,j) - 4*v(i-1,j) + v(i-2,j))/(2*hx)
                            dvdx_plus  = (-v(i+2,j) + 4*v(i+1,j) - 3*v(i,j))/(2*hx)
                        ENDIF                                   

                        IF (j == 2) THEN
                            dvdy_minus = (v(i,j) - v(i,j-1))/hy
                            dvdy_plus  = (v(i,j+1) - v(i,j))/hy
                        ELSEIF (j == Ny) THEN
                            dvdy_minus = (v(i,j) - v(i,j-1))/hy
                            dvdy_plus  = (v(i,j+1) - v(i,j))/hy
                        else
                            dvdy_minus = (3*v(i,j) - 4*v(i,j-1) + v(i,j-2))/(2*hy)
                            dvdy_plus  = (-v(i,j+2) + 4*v(i,j+1) - 3*v(i,j))/(2*hy)
                        ENDIF

                        ! dummy variable for the convective opperator
                        vstar_int = u_plus*dvdx_minus + u_minus*dvdx_plus + &
                                    v_plus*dvdy_minus + v_minus*dvdy_plus           

                        n1_yint_inert = n1_yint_inert + &
                                        0.5*(s_H(i,j)+s_H(i,j+1))*&
                                        rho_O/DeltaT*(v(i,j) - v_old(i,j))*&
                                        hx*hy
                                                                                        
                        n1_yint_conv = n1_yint_conv + &
                                       0.5*(s_H(i,j)+s_H(i,j+1))*&
                                       rho_O*vstar_int*hx*hy
                                
                        n1_yint_grav = n1_yint_grav + 0.5*&
                                       (s_H(i,j)+s_H(i,j+1))*rho_O*hx*hy*gy
                                                
                        ! Rotations
                        ! rs is the s-distance from the location at i,j to the centroid of the object.
                        rx = xwhole(i) - ncgx
                                                
                        r1_int_inert = r1_int_inert + &
                                       rx*0.5*(s_H(i,j)+s_H(i,j+1))*&
                                       rho_O/deltaT*(v(i,j) - v_old(i,j))*&
                                       hx*hy
                                                                                           
                        r1_int_conv = r1_int_conv + &
                                      rx*0.5*(s_H(i,j)+s_H(i,j+1))*&
                                      rho_O*vstar_int*hx*hy
                                
                        r1_int_grav = r1_int_grav + &
                                      rx*0.5*(s_H(i,j)+s_H(i,j+1))*&
                                      rho_O*hx*hy*gy
                                            
                    ENDIF
                ENDDO
            ENDDO
                
            ! Pressure and gravitational contributions to total energy
            DO i=2,Nx
                DO j=2,Ny     
              
                    ! density at the whole node
                    rho_O = rho_one*H(i,j) + rho_two*(1.0-H(i,j))

                    total_energy = total_energy + &
                                   ( -P(i,j) + rho_O*gy*(j-1)*hy ) * &
                                   hx*hy
                ENDDO
            ENDDO
            
            ! Upddate the convective angular velocity
            n1_rb_vom = n_rb_vom + 0.5*deltaT/buoy_J*&
                                   (r_int_inert + r_int_conv + r_int_grav) &
                                 + 0.5*deltaT/buoy_J*&
                                   (r1_int_inert + r1_int_conv + r1_int_grav)
                
                
            ! --------------------------------- !
            ! ---- RIGID BODY ACCELERATION ---- !
            ! --------------------------------- !   
                
            ! Compute updated acceleration ::                
            n1_rb_ax = (n1_xint_inert + n1_xint_conv + n1_xint_grav + int_forcex)/rb_mass - gx
            n1_rb_ay = (n1_yint_inert + n1_yint_conv + n1_yint_grav + int_forcey)/rb_mass - gy
            n1_rb_aom = 2/deltaT*(n1_rb_vom - n_rb_vom) - n_rb_aom + (int_torque/buoy_J)

            ! Restrict motion
            IF (RB_lock == 1) THEN
                n1_rb_ax = 0.0
                n1_rb_aom = 0.0      
            ENDIF 

            ! Calculate final velocity ::         
            n1_rb_vx = n12_rb_vx + 0.5*deltaT*n1_rb_ax                
            n1_rb_vy = n12_rb_vy + 0.5*deltaT*n1_rb_ay

            ! Calculate Energy ::
            total_energy = total_energy + &
                           rb_mass*gy*ncgy + &
                           0.5*rb_mass*((n1_rb_vx**2)+(n1_rb_vy**2)) + &
                           0.5*buoy_J*(n1_rb_vom**2)
            lin_momx = lin_momx + rb_mass*n1_rb_vx
            lin_momy = lin_momy + rb_mass*n1_rb_vy

        ENDIF   ! End rigid body    
                
        ! PRINT diveregent errors
        !DO i=2,Nx
        !        DO j=2,Ny
        !                IF ( abs((u(i,j)-u(i-1,j))/hx+(v(i,j)-v(i,j-1))/hy) >= 1.0e-10 ) THEN
        !                        print *, (u(i,j)-u(i-1,j))/hx+(v(i,j)-v(i,j-1))/hy 
        !                ENDIF
        !        ENDDO
        !ENDDO

        ! -------------------------------- !
        ! ---- ADVANCE VALUES IN TIME ---- !
        ! -------------------------------- !             
        ! n+1 -> n
                
        n_xint_conv  = n1_xint_conv
        n_xint_grav  = n1_xint_grav
        n_xint_inert = n1_xint_inert
        
        n_yint_conv  = n1_yint_conv
        n_yint_grav  = n1_yint_grav
        n_yint_inert = n1_yint_inert
                
        r_int_conv  = r1_int_conv
        r_int_grav  = r1_int_grav
        r_int_inert = r1_int_inert
          
        n_w_conv  = n1_w_conv
        n_w_grav  = n1_w_grav
        n_w_inert = n1_w_inert
                      
        ! Reinitialize intergral variables two different directions
                
        n1_xint_conv  = 0.0
        n1_xint_grav  = 0.0
        n1_xint_inert = 0.0
        
        n1_yint_conv  = 0.0
        n1_yint_grav  = 0.0
        n1_yint_inert = 0.0
                
        r1_int_conv  = 0.0
        r1_int_grav  = 0.0
        r1_int_inert = 0.0

        n1_w_conv  = 0.0
        n1_w_grav  = 0.0
        n1_w_inert = 0.0

        n_w_dy = n1_w_dy
        rb_data(39) = n_w_dy
              
        n_rb_dx = n1_rb_dx
        n_rb_dy = n1_rb_dy
        
        n_rb_vx = n1_rb_vx
        n_rb_vy = n1_rb_vy

        n_rb_ax = n1_rb_ax
        n_rb_ay = n1_rb_ay
                
        n_rbdelom(1,1) = n1_rbdelom(1,1)
        n_rbdelom(1,2) = n1_rbdelom(1,2)
        n_rbdelom(2,1) = n1_rbdelom(2,1)
        n_rbdelom(2,2) = n1_rbdelom(2,2)

        n_rb_vom = n1_rb_vom
        n_rb_aom = n1_rb_aom

        n_rb2_vom = n1_rb2_vom
        n_rb2_aom = n1_rb2_aom

        small_theta_n = small_theta_n1
             
        time = time + deltaT  ! increment the time to n+1
        time_count = time_count + deltaT  ! increment the print clock
        count = count + 1    ! Total iterations through the DO loop

        ! ------------------------------------- !
        ! ---- REINITIALIZE GALS LEVEL SET ---- !
        ! ------------------------------------- !   

        IF ((LS_type == 3) .AND. (time  > 0.0) .AND. (mod(time,reintime) < deltaT)) THEN
			CALL REBUILD(Nx,Ny,hx,hy,xwhole,ywhole,GA,scal,phiLS,phi_x,phi_y,phi_xy,prox,width)
        ENDIF

        IF ((LS_type /= 3) .AND. (time  > 0.0) .AND. (mod(time,reintime) < deltaT)) THEN
            CALL LSREIN(Nx,Ny,hx,hy,phiLS,H,&
                        BCphisouth,BCphinorth,BCphiwest,BCphieast,ywhole,d_phiLS)
        ENDIF

        !---------------------!
        !---- SAVING DATA ----!
        !---------------------!
        
        IF (time>=time_max .or. mod(time,time_print)<deltaT) THEN
            time_count = 0.0
            file_count = file_count + 1
            
            ! ----------------------------------------- !
            ! ---- REINITIALIZE STANDARD LEVEL SET ---- !
            ! ----------------------------------------- ! 
            
            !IF (LS_type == 1) THEN                                    
            !    CALL LSREIN(Nx,Ny,hx,hy,phiLS,H,&
            !                BCphisouth,BCphinorth,BCphiwest,BCphieast,ywhole,d_phiLS)
                !CALL LSVOLREIN(Nx,Ny,hx,hy,phiLS,H,&
                !               BCphisouth,BCphinorth,BCphiwest,BCphieast,&
                !              ywhole,d_phiLS,rho_one,rho_two,MassInt)
            !ENDIF
            ! apply level set BC's prior to saving
            !call BCLEVELSET(Nx,Ny,phiLS,2,2,1,1)
            !CALL BCLEVELSET(Nx,Ny,hx,hy,ywhole,phiLS,&
            !                BCphisouth,BCphinorth,BCphiwest,BCphieast,d_phiLS)
            !CALL BCLEVELSET(Nx,Ny,hx,hy,ywhole,s_phiLS,&
            !                BCphisouth,BCphinorth,BCphiwest,BCphieast,d_phiLS)
                                                           
            !call SAVEDATA(Nx,Ny,Lx,Ly,time,x,y,phiLS,s_phiLS,w_phiLS,H,s_H,w_H,P,u,v,file_count,rb_data)
            call SAVEDATA2(Nx,Ny,Lx,Ly,time,xwhole,ywhole,phiLS,s_phiLS,w_phiLS,H,s_H,w_H,P,u,v, &
                           file_count,rb_data,phi_x,phi_y,phi_xy,Fsig)
                                                      
        ENDIF        ! end of data saving
                
        !---- EXIT ----
        IF (time>=time_max) exit  ! exit IF t is greater than tmax

    ENDDO  ! end of time loop


    ! FORMATS for input file
    50  FORMAT (A)
    60  FORMAT (16X, I10)
    70  FORMAT (16X, E10.4)
    80  FORMAT (16X, L)
        
    !90  FORMAT (10F13.6)  !9F10.4
    90  FORMAT (10F13.9)  !9F10.4
    95  FORMAT (2I5)

END PROGRAM MPXLIB
!********************************************************************
!*                                                                                                                                        *
!*                                                                ENODIFF                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is for advecting scalar quantities T, PhiLS, etc. Note,                *
!* this program is only used near the boundaries when advecting the        *
!* the level set field because the WENODIFF is used in the interior.*
!*                                                                                                                                        *
!********************************************************************
SUBROUTINE ENODIFF(Nx,Ny,i,j,hx,hy,u,v,VAR,DiffX,DiffY)
        IMPLICIT NONE
        
        INTEGER, INTENT(IN):: Nx,Ny,i,j
        
        REAL(kind=8), INTENT(IN) :: hx,hy
        REAL(kind=8), INTENT(INOUT) :: DiffX,DiffY
        REAL(kind=8), DIMENSION(Nx+1,Ny+2), INTENT(IN) :: u
        REAL(kind=8), DIMENSION(Nx+2,Ny+1), INTENT(IN) :: v
        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: VAR

        REAL(kind=8) :: Dnegx,Dnegy,Dnegz,Dposx,Dposy,Dposz
        REAL(kind=8) :: DiffRight,DiffLeft,DiffTop,DiffBottom
        REAL(kind=8) :: Dnegnegx,Dposnegx,Dposposx,Dnegnegy,Dposnegy,Dposposy


        ! EXTERNAL FUNCTION DECLARATION
        REAL(kind=8), EXTERNAL :: MFCN

        DiffX = 0.0
        DiffY = 0.0

        ! difference operators
        Dnegx = (VAR(i,j) - VAR(i-1,j))/hx ! Back diff at i
        Dposx = (VAR(i+1,j) - VAR(i,j))/hx ! Forward diff at i
        Dnegy = (VAR(i,j) - VAR(i,j-1))/hy ! Back diff at j
        Dposy = (VAR(i,j+1) - VAR(i,j))/hy ! Forward diff at j

        IF ((i==2 .OR. i==Nx+1) .OR. (j==2 .OR. j==Ny+1)) THEN
                ! near edges I use 1st order differences
                DiffLeft = Dnegx
                DiffRight = Dposx
                DiffBottom = Dnegy
                DiffTop = Dposy
        ELSE
                Dnegnegx = (VAR(i-2,j) - 2.0*VAR(i-1,j) + VAR(i,j))/hx**2.0  ! Central 2diff at i-1
                Dposnegx = (VAR(i-1,j) - 2.0*VAR(i,j) + VAR(i+1,j))/hx**2.0  ! Central 2diff at i
                Dposposx = (VAR(i,j) - 2.0*VAR(i+1,j) + VAR(i+2,j))/hx**2.0  ! Central 2diff at i+1
                Dnegnegy = (VAR(i,j-2) - 2.0*VAR(i,j-1) + VAR(i,j))/hy**2.0  ! Central 2diff at j-1
                Dposnegy = (VAR(i,j-1) - 2.0*VAR(i,j) + VAR(i,j+1))/hy**2.0  ! Central 2diff at j
                Dposposy = (VAR(i,j) - 2.0*VAR(i,j+1) + VAR(i,j+2))/hy**2.0  ! Central 2diff at j+1        

                DiffLeft = Dnegx + hx/2.0*MFCN(Dposnegx, Dnegnegx)    ! Diff on left face
                DiffRight = Dposx - hx/2.0*MFCN(Dposnegx, Dposposx)   ! Diff on right face
                DiffBottom = Dnegy + hy/2.0*MFCN(Dposnegy, Dnegnegy)  ! Diff on bottom face
                DiffTop = Dposy - hy/2.0*MFCN(Dposnegy, Dposposy)     ! Diff on top face
        ENDIF

        IF (0.5*(u(i,j)+u(i-1,j)) > 0.0) THEN
                DiffX = DiffLeft
        ELSE
                DiffX = DiffRight
        END IF

        IF (0.5*(v(i,j)+v(i,j-1)) > 0.0) THEN
                DiffY = DiffBottom
        ELSE
                DiffY = DiffTop
        END IF


END SUBROUTINE ENODIFF
!********************************************************************
!*                                                                                                                                        *
!*                                                                MFNC                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is the switch function of the 2nd order ENO schemes                        *
!*                                                                                                                                        *
!********************************************************************
REAL(kind=8) FUNCTION MFCN(a,b)
        IMPLICIT NONE
        REAL(kind=8), INTENT(IN) :: a,b

        IF (a*b > 0.0) THEN
                IF (ABS(a) <= ABS(b)) THEN
                        MFCN = a;
                ELSEIF (ABS(a) > ABS(b)) THEN
                        MFCN = b;
                END IF
        ELSE
                MFCN = 0.0;
        END IF

END FUNCTION MFCN
!********************************************************************
!*                                                                                                                                        *
!*                                                          WENODIFF                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is for advecting the smooth scalar fields away from the                *
!* boundaries.  The scheme is up to 5th        order accurate and will                * 
!* always be more the 3rd order        accurate; the accuracy depends on        * 
!* the smoothness of the field.        For details see, Osher S., Fedkiw,  *
!* R. "Level Set Methods and Dynamic Implicit Surfaces", pp.33-37,        *
!* 2003.                                                                                                                        *
!*                                                                                                                                        *
!********************************************************************
SUBROUTINE WENODIFF(Nx,Ny,i,j,hx,hy,u,v,VAR,DiffX,DiffY)

        IMPLICIT NONE

        INTEGER, INTENT(IN):: Nx,Ny,i,j
        
        REAL(kind=8), INTENT(IN) :: hx,hy
        REAL(kind=8), INTENT(OUT) :: DiffX,DiffY
        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: VAR
        REAL(kind=8), DIMENSION(Nx+1,Ny+2), INTENT(IN) :: u
        REAL(kind=8), DIMENSION(Nx+2,Ny+1), INTENT(IN) :: v

        REAL(kind=8) :: vxone,vxtwo,vxthree,vxfour,vxfive
        REAL(kind=8) :: vyone,vytwo,vythree,vyfour,vyfive
        REAL(kind=8) :: varxone,varxtwo,varxthree,varyone,varytwo,varythree 
        REAL(kind=8) :: sxone,sxtwo,sxthree,syone,sytwo,sythree
        REAL(kind=8) :: axone,axtwo,axthree,ayone,aytwo,aythree,e
        REAL(kind=8) :: wone,wtwo,wthree


        ! X-direction
        
                IF (0.5*(u(i,j)+u(i-1,j))<0.0) THEN
                        ! forward differences, u<0
                        vxone = (VAR(i+3,j) - VAR(i+2,j))/hx
                        vxtwo = (VAR(i+2,j) - VAR(i+1,j))/hx
                        vxthree = (VAR(i+1,j) - VAR(i,j))/hx
                        vxfour = (VAR(i,j) - VAR(i-1,j))/hx
                        vxfive = (VAR(i-1,j) - VAR(i-2,j))/hx
                else
                        ! backward differences, u>0
                        vxone = (VAR(i-2,j) - VAR(i-3,j))/hx
                        vxtwo = (VAR(i-1,j) - VAR(i-2,j))/hx
                        vxthree = (VAR(i,j) - VAR(i-1,j))/hx
                        vxfour = (VAR(i+1,j) - VAR(i,j))/hx
                        vxfive = (VAR(i+2,j) - VAR(i+1,j))/hx
                ENDIF

                varxone = vxone/3.0 - 7.0/6.0*vxtwo + 11.0/6.0*vxthree
                varxtwo = -vxtwo/6.0 + 5.0/6.0*vxthree + vxfour/3.0
                varxthree = vxthree/3.0 + 5.0/6.0*vxfour - vxfive/6.0

                ! smoothness parameters
                sxone = 13.0/12.0*(vxone - 2.0*vxtwo + vxthree)**2.0 + &
                                         0.25*(vxone - 4.0*vxtwo + 3.0*vxthree)**2.0
                sxtwo = 13.0/12.0*(vxtwo - 2.0*vxthree + vxfour)**2.0 + 0.25*(vxtwo - vxfour)**2.0
                sxthree = 13.0/12.0*(vxthree - 2.0*vxfour + vxfive)**2.0 + &
                                           0.25*(3.0*vxthree - 4.0*vxfour + 3.0*vxfive)**2.0

                ! find smoothness ratios
                e = 0.000001*MAX(vxone**2.0,vxtwo**2.0,vxthree**2.0,vxfour**2.0,vxfive**2.0) + 1.0e-15
                axone = 0.1/(sxone + e)**2.0
                axtwo = 0.6/(sxtwo + e)**2.0
                axthree = 0.3/(sxthree + e)**2.0

                ! find weights, THEN find derivative
                wone = axone/(axone + axtwo + axthree)
                wtwo = axtwo/(axone + axtwo + axthree)
                wthree = axthree/(axone + axtwo + axthree)
                DiffX = wone*varxone + wtwo*varxtwo + wthree*varxthree

        
        ! Y-direction
        
                IF (0.5*(v(i,j)+v(i,j-1))<0.0) THEN
                        ! forward differences, v<0
                        vyone = (VAR(i,j+3) - VAR(i,j+2))/hy
                        vytwo = (VAR(i,j+2) - VAR(i,j+1))/hy
                        vythree = (VAR(i,j+1) - VAR(i,j))/hy
                        vyfour = (VAR(i,j) - VAR(i,j-1))/hy
                        vyfive = (VAR(i,j-1) - VAR(i,j-2))/hy
                else
                        ! backward differences, v>0
                        vyone = (VAR(i,j-2) - VAR(i,j-3))/hy
                        vytwo = (VAR(i,j-1) - VAR(i,j-2))/hy
                        vythree = (VAR(i,j) - VAR(i,j-1))/hy
                        vyfour = (VAR(i,j+1) - VAR(i,j))/hy
                        vyfive = (VAR(i,j+2) - VAR(i,j+1))/hy
                ENDIF

                varyone = vyone/3.0 - 7.0/6.0*vytwo + 11.0/6.0*vythree
                varytwo = -vytwo/6.0 + 5.0/6.0*vythree + vyfour/3.0
                varythree = vythree/3.0 + 5.0/6.0*vyfour - vyfive/6.0

                ! smoothness parameters
                syone = 13.0/12.0*(vyone - 2.0*vytwo + vythree)**2.0 + &
                                         0.25*(vyone - 4.0*vytwo + 3.0*vythree)**2.0
                sytwo = 13.0/12.0*(vytwo - 2.0*vythree + vyfour)**2.0 + 0.25*(vytwo - vyfour)**2.0
                sythree = 13.0/12.0*(vythree - 2.0*vyfour + vyfive)**2.0 + &
                                           0.25*(3.0*vythree - 4.0*vyfour + 3.0*vyfive)**2.0

                ! find smoothness ratios
                e = 0.000001*MAX(vyone**2.0,vytwo**2.0,vythree**2.0,vyfour**2.0,vyfive**2.0) + 1.0e-15
                ayone = 0.1/(syone + e)**2.0
                aytwo = 0.6/(sytwo + e)**2.0
                aythree = 0.3/(sythree + e)**2.0

                ! find weights, THEN find derivative
                wone = ayone/(ayone + aytwo + aythree)
                wtwo = aytwo/(ayone + aytwo + aythree)
                wthree = aythree/(ayone + aytwo + aythree)
                DiffY = wone*varyone + wtwo*varytwo + wthree*varythree

END SUBROUTINE WENODIFF
!********************************************************************
!*                                                                                                                                        *
!*                                                            LSREIN                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                          *
!*                                                                                                                                        *
!* This is a program for reinitializing to a sign distance function        *
!* using the technique suggested by, Sussman, M., Smereka, P. and        *
!* Osher, S.J., "A Level Set Method for Computing Solutions to                *
!* Solutions to Incompressible Two-Phase Flow, J. Computational                *
!* Physics, Vol. 114, pp.146-159, 1994.                                                                * 
!*                                                                                                                                        *
!********************************************************************
SUBROUTINE LSREIN(Nx,Ny,hx,hy,phiLS,H,&
                                  BCphisouth,BCphinorth,BCphiwest,BCphieast,ywhole,num)
        IMPLICIT NONE

        INTEGER, INTENT(IN):: Nx,Ny,BCphisouth,BCphinorth,BCphiwest,BCphieast
        
        INTEGER :: i,j,k,step,iter

        REAL(kind=8) :: deltaTau,u

        REAL(kind=8), INTENT(IN) :: hx,hy,num
        REAL(kind=8), DIMENSION(Ny+2), INTENT(IN) :: ywhole
        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(INOUT) :: H
        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(INOUT) :: phiLS
        REAL(kind=8), DIMENSION(Nx+2,Ny+2) :: phiLSn,RHS
        REAL(kind=8), DIMENSION(4) :: Coef
        
        ! EXTERANL FUNCTION DECLARATION
        REAL(kind=8), EXTERNAL :: DISTRHS
        REAL(kind=8), EXTERNAL :: LSHEAVY


        deltaTau = hx/2.0
        
        Coef(1) = 1.0/4.0
        Coef(2) = 1.0/3.0
        Coef(3) = 1.0/2.0
        Coef(4) = 1.0
        
        ! I need the heavieside fcn for reinitialization
        ! finding the new heavieside function at each node
        DO i=1,Nx+2
                DO j=1,Ny+2
                        ! heavieside function
                        H(i,j) = LSHEAVY(Nx,Ny,i,j,hx,hy,phiLS)
                ENDDO
        ENDDO

        ! iter is the number of artificial time steps        
        DO iter=1,4

                !////////////////// RK4 //////////////////////////
                ! PDE => d(phi)/d(Tau) = u*(1 - ||grad(phi)||)
                ! note: u corresponds to the characteristic velocity
                phiLSn = phiLS
                DO step=1,4
                                                
                        DO i=2,Nx+1
                                DO j=2,Ny+1
                                        RHS(i,j) = DISTRHS(Nx,Ny,i,j,hx,hy,phiLS,H)                
                                ENDDO
                        ENDDO
                        
                        DO i=2,Nx+1
                                DO j=2,Ny+1
                                        ! note: Coef is positive
                                        phiLS(i,j) = phiLSn(i,j) + Coef(step)*deltaTau*RHS(i,j)                                
                                ENDDO
                        ENDDO
                        
                        ! apply boundary conditions
                        !call BCLEVELSET(Nx,Ny,phiLS,2,2,1,1)
                        call BCLEVELSET(Nx,Ny,hx,hy,ywhole,phiLS,&
                                                    BCphisouth,BCphinorth,BCphiwest,BCphieast,num)
                ENDDO

        ENDDO ! end artificial time loop


END SUBROUTINE LSREIN
!********************************************************************
!*                                                                  *
!*                            LSVOLREIN                             *
!*                                                                  *
!********************************************************************
!* Author: Nathaniel R. Morgan                                      *
!*                                                                  *
!* This is a program for conserving the volume and reinitializing   *
!* to a sign distance function using the technique suggested by,    *
!* Sussman, M., Smereka, P. and Osher, S.J., "A Level Set Method    *
!* for Computing Solutions to Solutions to Incompressible Two-Phase *
!* Flow, J. Computational Physics, Vol. 114, pp.146-159, 1994       * 
!*                                                                  *
!* Edited by C. Lee to remove distance reinitialization (This is    *
!* accomplished after LSVOLREIN is finished).                       *
!*                                                                  *
!********************************************************************
SUBROUTINE LSVOLREIN(Nx,Ny,hx,hy,phiLS,H,&
                     BCphisouth,BCphinorth,BCphiwest,BCphieast,&
                     ywhole,num,rho_one,rho_two,MassInt)
    IMPLICIT NONE

    INTEGER, INTENT(IN):: Nx,Ny,BCphisouth,BCphinorth,BCphiwest,BCphieast
    
    INTEGER :: i,j,k,step,iter

    REAL(kind=8) :: deltaTau,u,Mass

    REAL(kind=8), INTENT(IN) :: hx,hy,num,rho_one,rho_two,MassInt
    REAL(kind=8), DIMENSION(Ny+2), INTENT(IN) :: ywhole
    REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(INOUT) :: H
    REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(INOUT) :: phiLS
    REAL(kind=8), DIMENSION(Nx+2,Ny+2) :: phiLSn,RHS
    
    ! EXTERANL FUNCTION DECLARATION
    REAL(kind=8), EXTERNAL :: LSHEAVY
    REAL(kind=8), EXTERNAL :: DISTRHS
    REAL(kind=8), EXTERNAL :: VOLRHS


    deltaTau = hx/2.0/2.0
    
    
    ! iter is the number of artificial time steps   
    DO iter=1,140

        ! get new heavieside function and Mass
        Mass = 0.0
        DO i=2,Nx+1
            DO j=2,Ny+1
                H(i,j) = LSHEAVY(Nx,Ny,i,j,hx,hy,phiLS)
                Mass = (rho_one*H(i,j) + rho_two*(1.0-H(i,j)))*hx*hy + Mass
            ENDDO
        ENDDO

        print *, abs(MassInt-Mass), MassInt, Mass
        IF (abs(MassInt-Mass)/MassInt<=1.0e-4) EXIT
        

        ! FIRST STEP -> reinitializing to a distance function
        ! PDE => d(phi)/d(Tau) = u*(1 - ||grad(phi)||)
        ! note: u corresponds to the characteristic velocity
        !phiLSn = phiLS                      
        !DO i=2,Nx+1
        !    DO j=2,Ny+1
        !        RHS(i,j) = DISTRHS(Nx,Ny,i,j,hx,hy,phiLS,H)     
        !    ENDDO
        !ENDDO   
        !DO i=2,Nx+1
        !    DO j=2,Ny+1
        !        phiLS(i,j) = phiLSn(i,j) + deltaTau*RHS(i,j)                
        !    ENDDO
        !ENDDO

        
        ! SECOND STEP -> reinitializing to satifies volume conservation
        ! PDE => d(phi)/d(Tau) = -(Massint-Mass)(kappa-1)*||graph(phi)||
        phiLSn = phiLS
        DO i=2,Nx+1
            DO j=2,Ny+1
                RHS(i,j) = VOLRHS(Nx,Ny,i,j,hx,hy,phiLS,MassInt,Mass)
            ENDDO
        ENDDO
        DO i=2,Nx+1
            DO j=2,Ny+1
                phiLS(i,j) = phiLSn(i,j) + deltaTau*RHS(i,j)                
            ENDDO
        ENDDO
        
        ! APPLY BOUNDARY CONDITIONS
        call BCLEVELSET(Nx,Ny,hx,hy,ywhole,phiLS,&
                        BCphisouth,BCphinorth,BCphiwest,BCphieast,num)

    ENDDO ! end artificial time loop


END SUBROUTINE LSVOLREIN
!********************************************************************
!*                                                                  *
!*                            VOLRHS                                *
!*                                                                  *
!********************************************************************
!* Author: Nathaniel R. Morgan                                      *
!*                                                                  *
!* This is a function used in reinitializing the level set field to *
!* a signed distance function.                                      *
!*                                                                  *
!* Edited by C. Lee to remove influence of curvature (which is      *
!* unstable - it shouldn't be though. To be revisited...)           *
!*                                                                  *
!********************************************************************
REAL(kind=8) FUNCTION VOLRHS(Nx,Ny,i,j,hx,hy,phiLS,MassInt,Mass)
    
    IMPLICIT NONE

    REAL(kind=8) :: u, NormGrad,kappa

    INTEGER, INTENT(IN):: Nx,Ny,i,j
    REAL(kind=8), INTENT(IN) :: hx,hy,MassInt,Mass
    REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: phiLS
    
    ! EXTERANL FUNCTION DECLARATION
    REAL(kind=8), EXTERNAL :: LSHEAVY, NORMDIFFPHI, LSCURVATURE
    
    
    kappa = LSCURVATURE(Nx,Ny,i,j,hx,hy,phiLS)

    print *, 'k ',kappa

    u = - (MassInt - Mass)*(-1.0+kappa)/(Nx*Ny)  ! mass error = characteristic velocity
    !u = (MassInt - Mass)/(Nx*Ny)
    
    NormGrad = NORMDIFFPHI(Nx,Ny,i,j,hx,hy,u,phiLS)
    
    VOLRHS = u*NormGrad  ! u*||grad(phi)||
    !VOLRHS = 0.0
END FUNCTION VOLRHS
!********************************************************************
!*                                                                                                                                        *
!*                                                          DISTRHS                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a function used in reinitializing the level set field to        *
!* a signed distance function.                                                                                *
!*                                                                                                                                        *
!********************************************************************
REAL(kind=8) FUNCTION DISTRHS(Nx,Ny,i,j,hx,hy,phiLS,H)
        
        IMPLICIT NONE

        INTEGER, INTENT(IN):: Nx,Ny,i,j
        REAL(kind=8), INTENT(IN) :: hx,hy
        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: phiLS,H
        REAL(kind=8) :: u, NormGrad

        ! EXTERANL FUNCTION DECLARATION
        REAL(kind=8), EXTERNAL :: LSHEAVY, NORMDIFFPHI

        
        u = 2.0*(H(i,j) - 1.0/2.0)  ! smoothed signum function
        NormGrad = NORMDIFFPHI(Nx,Ny,i,j,hx,hy,u,phiLS)
        
        DISTRHS = u*(1.0 - NormGrad)  ! u*(1 - ||grad(phi)||)

END FUNCTION DISTRHS
!********************************************************************
!*                                                                                                                                        *
!*                                                        NORMDIFFPHI                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a program for calculating the norm of the gradient of        *
!* phi using the method outlined in, Sethian, J.A., "Level Set                *
!* Methods and Fast Marching Methods", p 66, 2002                                        *
!*                                                                                                                                        *        
!********************************************************************
REAL(kind=8) FUNCTION NORMDIFFPHI(Nx,Ny,i,j,hx,hy,u,phi)
        
        IMPLICIT NONE

        INTEGER, INTENT(IN):: Nx,Ny,i,j

        REAL(kind=8), INTENT(IN) :: hx,hy,u
        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: phi
        REAL(kind=8) :: Dnegx,Dnegy,Dposx,Dposy, &
                                        DiffRight,DiffLeft,DiffTop,DiffBottom, &
                                        Dnegnegx,Dposnegx,Dposposx, &
                                        Dnegnegy,Dposnegy,Dposposy

        ! EXTERNAL FUNCTION DECLARATION
        REAL(kind=8), EXTERNAL :: MFCN

        ! difference operators
        Dnegx = (phi(i,j) - phi(i-1,j))/hx ! Back diff at i
        Dposx = (phi(i+1,j) - phi(i,j))/hx ! Forward diff at i
        Dnegy = (phi(i,j) - phi(i,j-1))/hy ! Back diff at j
        Dposy = (phi(i,j+1) - phi(i,j))/hy ! Forward diff at j

        IF ((i==2 .OR. i==Nx+1) .OR. (j==2 .OR. j==Ny+1)) THEN
                ! near edges I use 1st order differences
                DiffLeft = Dnegx
                DiffRight = Dposx
                DiffBottom = Dnegy
                DiffTop = Dposy

        ELSE
                Dnegnegx = (phi(i-2,j) - 2.0*phi(i-1,j) + phi(i,j))/hx**2.0  ! Central 2diff at i-1
                Dposnegx = (phi(i-1,j) - 2.0*phi(i,j) + phi(i+1,j))/hx**2.0  ! Central 2diff at i
                Dposposx = (phi(i,j) - 2.0*phi(i+1,j) + phi(i+2,j))/hx**2.0  ! Central 2diff at i+1
                Dnegnegy = (phi(i,j-2) - 2.0*phi(i,j-1) + phi(i,j))/hy**2.0  ! Central 2diff at j-1
                Dposnegy = (phi(i,j-1) - 2.0*phi(i,j) + phi(i,j+1))/hy**2.0  ! Central 2diff at j
                Dposposy = (phi(i,j) - 2.0*phi(i,j+1) + phi(i,j+2))/hy**2.0  ! Central 2diff at j+1        

                ! function M is defined in Advection
                DiffLeft = Dnegx + hx/2.0*MFCN(Dposnegx, Dnegnegx)    ! Diff on left face
                DiffRight = Dposx - hx/2.0*MFCN(Dposnegx, Dposposx)   ! Diff on right face
                
                DiffBottom = Dnegy + hy/2.0*MFCN(Dposnegy, Dnegnegy)  ! Diff on bottom face
                DiffTop = Dposy - hy/2.0*MFCN(Dposnegy, Dposposy)     ! Diff on top face
                
        ENDIF

        IF (u >= 0.0) THEN
                NORMDIFFPHI = SQRT( (MAX(DiffLeft,0.0))**2.0 + &
                                                        (MIN(DiffRight,0.0))**2.0 + &
                                                        (MAX(DiffBottom,0.0))**2.0 + &
                                                        (MIN(DiffTop,0.0))**2.0 )
        ELSE
                NORMDIFFPHI = SQRT( (MIN(DiffLeft,0.0))**2.0 + &
                                                        (MAX(DiffRight,0.0))**2.0 + &
                                                        (MIN(DiffBottom,0.0))**2.0 + &
                                                        (MAX(DiffTop,0.0))**2.0 )
        ENDIF
        
END FUNCTION NORMDIFFPHI
!********************************************************************
!*                                                                                                                                        *
!*                                                          LSHEAVY                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a program for calculating the LS heavieside function                *
!*                                                                                                                                        *        
!********************************************************************
REAL(kind=8) FUNCTION LSHEAVY(Nx,Ny,i,j,hx,hy,phiLS)
        
        IMPLICIT NONE

        INTEGER, INTENT(IN):: Nx,Ny,i,j

        REAL(kind=8), INTENT(IN) :: hx,hy

        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: phiLS
        
        REAL(kind=8) :: pi,epsilon

        PARAMETER(pi=3.141592653589)

        ! smearing the heavieside function over 1.5 cells in each direction 
        epsilon = 1.5*max(hx,hy)

        IF (phiLS(i,j) >= epsilon) THEN
                LSHEAVY = 1.0
        ELSEIF (phiLS(i,j) <= -epsilon) THEN
                LSHEAVY = 0.0
        ELSE
                LSHEAVY = 0.5 + phiLS(i,j)/(2.0*epsilon) + &
                                  sin(pi*phiLS(i,j)/epsilon)/(2.0*pi)
        ENDIF

END FUNCTION LSHEAVY
!********************************************************************
!*                                                                                                                                        *
!*                                                          LSHEAVY2                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a program for calculating the LS heavieside function                *
!*                                                                                                                                        *        
!********************************************************************
REAL(kind=8) FUNCTION LSHEAVY2(Nx,Ny,i,j,hx,hy,phiLS,phi_x,phi_y)
        
        IMPLICIT NONE

        INTEGER, INTENT(IN):: Nx,Ny,i,j

        REAL(kind=8), INTENT(IN) :: hx,hy

        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: phiLS,phi_x,phi_y
        
        REAL(kind=8) :: pi,epsilon,gradmag,phiadj

        PARAMETER(pi=3.141592653589)

        ! smearing the heavieside function over 1.5 cells in each direction 
        epsilon = 1.5*max(hx,hy)

        gradmag = sqrt((phi_x(i,j)**2 + phi_y(i,j)**2))

        IF ((abs(gradmag) > 1e-5) .AND. (i>1) .AND. (j>1) .AND. (i<Nx+2) .AND. (j<Ny+2)) THEN
            phiadj = phiLS(i,j)/gradmag
        ELSE
            phiadj = phiLS(i,j)
        ENDIF

        IF (phiadj >= epsilon) THEN
                LSHEAVY2 = 1.0
        ELSEIF (phiadj <= -epsilon) THEN
                LSHEAVY2 = 0.0
        ELSE
                LSHEAVY2 = 0.5 + phiadj/(2.0*epsilon) + &
                                  sin(pi*phiadj/epsilon)/(2.0*pi)
        ENDIF

END FUNCTION LSHEAVY2
!********************************************************************
!*                                                                                                                                        *
!*                                                          LSDELTA                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a program for calculating the LS delta function                        *
!*                                                                                                                                        *        
!********************************************************************
REAL(kind=8) FUNCTION LSDELTA(Nx,Ny,i,j,hx,hy,phiLS)
        
        IMPLICIT NONE

        INTEGER, INTENT(IN):: Nx,Ny,i,j

        REAL(kind=8), INTENT(IN) :: hx,hy
        real(kind=8) :: epsilon

        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: phiLS
        
        REAL(kind=8) :: pi

        PARAMETER(pi=3.141592653589)

        ! smearing the delta function over 1.5 cells in each direction 
        epsilon = 1.5*max(hx,hy)

        IF (abs(phiLS(i,j)) > epsilon) THEN
                LSDELTA = 0.0
        ELSE
                LSDELTA = 1.0/(2.0*epsilon)*( 1.0 + cos(pi*phiLS(i,j)/epsilon) )
        ENDIF

END FUNCTION LSDELTA
!********************************************************************
!*                                                                                                                                        *
!*                                                        LSNORMAL                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a program for calculating the normal of the level set        *
!* field.                                                                                                                        *
!*                                                                                                                                        *
!********************************************************************
SUBROUTINE LSNORMAL(Nx,Ny,i,j,hx,hy,phiLS,normal)

        IMPLICIT NONE
        
        INTEGER, INTENT(IN):: Nx,Ny,i,j 

        REAL(kind=8) :: phi_x,phi_y,norm
        
        REAL(kind=8), INTENT(IN) :: hx, hy
        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: phiLS
        REAL(kind=8), DIMENSION(2), INTENT(OUT) :: normal

        ! the i,j values correspond to the whole nodes

        ! first derivatives
        phi_x = ( phiLS(i+1,j) - phiLS(i-1,j) )/(2.0*hx)
        phi_y = ( phiLS(i,j+1) - phiLS(i,j-1) )/(2.0*hy)

        ! magnitude of the the gradient of phi
        norm = (phi_x**2.0 + phi_y**2.0)**0.5

        ! normal vector
        normal(1) = phi_x/norm
        normal(2) = phi_y/norm

        !normal(1) = 0.0                        !??????
        !normal(2) = 1.0                        !??????

ENDSUBROUTINE LSNORMAL
!********************************************************************
!*                                                                                                                                        *
!*                                                        LSCURVATURE                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a program for calculating the curvature of the level set        *
!* field.                                                                                                                        *
!*                                                                                                                                        *
!********************************************************************
REAL(kind=8) FUNCTION LSCURVATURE(Nx,Ny,i,j,hx,hy,phiLS)

        IMPLICIT NONE
        
        INTEGER, INTENT(IN):: Nx,Ny,i,j 

        REAL(kind=8) :: phi_x,phi_y,phi_xy,phi_xx,phi_yy,norm
        !REAL(kind=8), DIMENSION(2) :: normal
        
        REAL(kind=8), INTENT(IN) :: hx, hy
        REAL(kind=8), DIMENSION(Nx+2,Ny+2), INTENT(IN) :: phiLS
        

        ! the i,j values correspond to the whole nodes

        ! first derivatives
        phi_x = ( phiLS(i+1,j) - phiLS(i-1,j) )/(2.0*hx)
        phi_y = ( phiLS(i,j+1) - phiLS(i,j-1) )/(2.0*hy)
        
        ! second derivatives
        phi_xx = ( phiLS(i+1,j) - 2.0*phiLS(i,j) + phiLS(i-1,j) )/(hx**2.0)
        phi_yy = ( phiLS(i,j+1) - 2.0*phiLS(i,j) + phiLS(i,j-1) )/(hy**2.0)

        ! mixed derivative
        phi_xy = (phiLS(i+1,j+1) - phiLS(i+1,j-1))/(4.0*hx*hy) - &
                         (phiLS(i-1,j+1) - phiLS(i-1,j-1))/(4.0*hx*hy)

        ! magnitude of the the gradient of phi
        norm = (phi_x**2.0 + phi_y**2.0)**0.5

        ! normal vector
        !normal(1) = phi_x/norm
        !normal(2) = phi_y/norm

        ! curvature
		IF (norm > 1e-16) then
        	LSCURVATURE = (phi_xx*phi_y**2.0 - 2.0*phi_y*phi_x*phi_xy + phi_yy*phi_x**2.0)/norm**3.0
        ELSE
        	LSCURVATURE = 0.0
        ENDIF
            
        ! bounding the curvature to the grid scale, (-1/hx<=kappa<=1/hx)
        IF ( abs(LSCURVATURE) >= 1.0/min(hx,hy) ) THEN
                IF (LSCURVATURE < 0.0) THEN
                        LSCURVATURE = -1.0/min(hx,hy)
                else
                        LSCURVATURE = 1.0/min(hx,hy)
                ENDIF
        ENDIF
     
END FUNCTION LSCURVATURE
!********************************************************************
!*                                                                                                                                        *
!*                                                        BCVELOCITY                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a program for applying boundary conditions to velocity        *                                                                                                                        
!*                                                                                                                                        *
!********************************************************************
SUBROUTINE BCVELOCITY(Nx,Ny,u,v,BCnorth,BCsouth,BCeast,BCwest,&
                                          unorth,usouth,ueast,uwest,&
                                          vnorth,vsouth,veast,vwest,&
                                          hx,hy,phiLS,contC,contA,contR,contM,par)
        implicit none

        integer :: i,j
        integer, intent(in) :: Nx,Ny,par
        integer, intent(in) :: BCnorth,BCsouth,BCeast,BCwest

        real(kind=8), intent(in) :: unorth,usouth,ueast,uwest,&
                                                                vnorth,vsouth,veast,vwest
        real(kind=8), intent(in) :: hx,hy,contC,contA,contR,contM

        real(kind=8), dimension(Nx+1,Ny+2), intent(inout) :: u
        real(kind=8), dimension(Nx+2,Ny+1), intent(inout) :: v
        real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phiLS

        real(kind=8), external :: CONTACTSPEED

        ! par is a parameter to distinguish between viscous and advection BC's
        ! par = 1 corresponds to advection
        ! par = 2 corresponds to viscous

        IF (par==1) THEN
                ! apply advection BC's

                ! ---- NORTH BOUNDARY ----
                IF (BCnorth == 1) THEN                ! dirclet wall
                        v(2:Nx+1,Ny+1) = 0.0                                                        ! no flow through wall
                        u(1:Nx+1,Ny+2) = 2.0*unorth - u(1:Nx+1,Ny+1)        ! ghost node horizontal velocity
                ELSEIF (BCnorth == 2) THEN        ! neumman wall
                        v(2:Nx+1,Ny+1) = 0.0                                                        ! no flow through wall
                        u(1:Nx+1,Ny+2) = u(1:Nx+1,Ny+1)                                        ! horizontal slip
                ELSEIF (BCnorth == 3) THEN        ! free surface
                        v(2:Nx+1,Ny+1) = v(2:Nx+1,Ny)                                        ! constant velocity
                        u(1:Nx+1,Ny+2) = u(1:Nx+1,Ny+1)                                        ! constant velocity
                ELSEIF (BCnorth == 4) THEN        ! pressure
                        v(2:Nx+1,Ny+1) = v(2:Nx+1,Ny) !2.0*v(2:Nx+1,Ny)-v(2:Nx+1,Ny-1) ! extrapolation
                        u(1:Nx+1,Ny+2) = u(1:Nx+1,Ny+1)        !2.0*u(1:Nx+1,Ny+1)-u(1:Nx+1,Ny) ! extrapolation
                ELSEIF (BCnorth == 5) THEN !flow
                        v(2:Nx+1,Ny+1) = vnorth                                                        ! vertical velocity
                        u(1:Nx+1,Ny+2) = 2.0*unorth - u(1:Nx+1,Ny+1)        ! ghost node horizonal velocity
                ELSEIF (BCnorth == 6) THEN        ! contact line
                        v(2:Nx+1,Ny+1) = 0.0                                                        ! no flow through wall
                        j=Ny+2
                        DO i=1,Nx+1
                                u(i,j)=CONTACTSPEED(Nx,Ny,i,j,hx,hy,phiLS,contC,contA,contR,contM)
                        ENDDO
                        
                ENDIF


                ! ---- SOUTH BOUNDARY ----
                IF (BCsouth == 1) THEN                ! dirclet wall
                        v(2:Nx+1,1) = 0.0                                                                ! no flow through wall
                        u(1:Nx+1,1) = 2.0*usouth - u(1:Nx+1,2)                        ! ghost node horizontal velocity
                ELSEIF (BCsouth == 2) THEN        ! neumman wall
                        v(2:Nx+1,1) = 0.0                                                                ! no flow through wall
                        u(1:Nx+1,1) = u(1:Nx+1,2)                                                ! horizontal slip
                ELSEIF (BCsouth == 3) THEN        ! free surface
                        v(2:Nx+1,1) = v(2:Nx+1,2)                                                 ! constant velocity                                                                
                        u(1:Nx+1,1) = u(1:Nx+1,2)                                                ! constant velocity
                ELSEIF (BCsouth == 4) THEN        ! pressure
                        ! 2.0*v(2:Nx+1,2)-v(2:Nx+1,3)
                ELSEIF (BCsouth == 5) THEN        ! flow boundary
                        v(2:Nx+1,1) = vsouth                                                        ! vertical velocity
                        u(1:Nx+1,1) = 2.0*usouth - u(1:Nx+1,2)                        ! ghost node horizontal velocity
                ELSEIF (BCsouth == 6) THEN        ! contact line
                        v(2:Nx+1,1) = 0.0                                                                ! no flow through wall
                        j=1
                        DO i=1,Nx+1
                                u(i,j)=CONTACTSPEED(Nx,Ny,i,j,hx,hy,phiLS,contC,contA,contR,contM)
                        ENDDO
                        
                ENDIF
                

                ! ---- EAST BOUNDARY ----
                IF (BCeast == 1) THEN                ! dirclet wall        
                        u(Nx+1,2:Ny+1) = 0.0                                                        ! no flow through wall
                        v(Nx+2,1:Ny+1) = 2.0*veast - v(Nx+1,1:Ny+1)                ! ghost node vertical velocity
                ELSEIF (BCeast == 2) THEN        ! neumman wall
                        u(Nx+1,2:Ny+1) = 0.0                                                        ! no flow through wall
                        v(Nx+2,1:Ny+1) = v(nx+1,1:ny+1)                                        ! vertical slip
                ELSEIF (BCeast == 3) THEN        ! free surface
                        u(Nx+1,2:Ny+1) = u(Nx,2:Ny+1)                                        ! constant velocity                                                        
                        v(Nx+2,1:Ny+1) = v(Nx+1,1:ny+1)                                        ! constant velocity
                ELSEIF (BCeast == 4) THEN        ! pressure
                        ! 2.0*u(Nx,2:Ny+1)-u(Nx-1,2:Ny+1)
                ELSEIF (BCeast == 5) THEN        ! flow boundary
                        u(Nx+1,2:Ny+1) = ueast                                                        ! horizontal velocity
                        v(Nx+2,1:Ny+1) = 2.0*veast - v(nx+1,1:ny+1)                ! ghost node vertical velocity
                ELSEIF (BCeast == 6) THEN        ! contact line
                        u(Nx+1,2:Ny+1) = 0.0                                                        ! no flow through wall
                        i=Nx+2
                        DO j=1,Ny+1
                                v(i,j)=CONTACTSPEED(Nx,Ny,i,j,hx,hy,phiLS,contC,contA,contR,contM)
                        ENDDO
                                
                ENDIF
                

                ! ---- WEST BOUNDARY ----        
                IF (BCwest == 1) THEN                ! dirclet wall
                        u(1,2:Ny+1) = 0.0                                                                ! no flow through wall
                        v(1,1:Ny+1) = 2.0*vwest - v(2,1:Ny+1)                        ! ghost node vertical velocity
                ELSEIF (BCwest == 2) THEN        ! neumman wall
                        u(1,2:Ny+1) = 0.0                                                                ! no flow through wall
                        v(1,1:Ny+1) = v(2,1:Ny+1)                                                ! vertical slip
                ELSEIF (BCwest == 3) THEN        ! free surface
                        u(1,2:Ny+1) = u(2,2:Ny+1)                                                 ! constant velocity                                                                
                        v(1,1:Ny+1) = v(2,1:Ny+1)                                                ! constant velocity
                ELSEIF (BCeast == 4) THEN        ! pressure
                        ! 2.0*u(2,2:Ny+1)-u(3,2:Ny+1)
                ELSEIF (BCwest == 5) THEN        !flow boundary
                        u(1,2:Ny+1) = uwest                                                                ! horizontal velocity
                        v(1,1:Ny+1) = 2.0*vwest - v(2,1:Ny+1)                        ! ghost node vertical velocity
                ELSEIF (BCwest == 6) THEN        ! contact line
                        u(1,2:Ny+1) = 0.0                                                                ! no flow through wall
                        i=1
                        DO j=2,Ny+1
                                v(i,j)=CONTACTSPEED(Nx,Ny,i,j,hx,hy,phiLS,contC,contA,contR,contM)
                        ENDDO
                        
                ENDIF

        ELSEIF (par==2) THEN
                ! apply viscous BCs

                ! ---- NORTH BOUNDARY ----
                IF (BCnorth == 1 .or. BCnorth == 5) THEN        ! dirclet wall or flow boundary
                        ! higher-order interpolation for ghost node
                        u(1:nx+1,ny+2)=8.0*unorth/3.0 - 2.0*u(1:nx+1,ny+1) + u(1:nx+1,ny)/3.0
                ENDIF

                ! ---- SOUTH BOUNDARY ----
                IF (BCsouth == 1 .or. BCsouth == 5) THEN        ! dirclet wall or flow boundary
                        ! higher-order interpolation for ghost node
                        u(1:nx+1,1)=8.0*usouth/3.0 - 2.0*u(1:nx+1,2) + u(1:nx+1,3)/3.0
                ENDIF
                
                ! ---- EAST BOUNDARY ----
                IF (BCeast == 1 .or. BCeast == 5) THEN                ! dirclet wall or flow boundary
                        ! higher-order interpolation for ghost node
                        v(nx+2,1:ny+1)=8.0*veast/3.0 - 2.0*v(nx+1,1:ny+1) + v(nx,1:ny+1)/3.0
                ENDIF
                
                ! ---- WEST BOUNDARY ----
                IF (BCwest == 1 .or. BCwest == 5) THEN                ! dirclet wall or flow boundary
                        ! higher-order interpolation for ghost node
                        v(1,1:ny+1)=8.0*vwest/3.0 - 2.0*v(2,1:ny+1) + v(3,1:ny+1)/3.0
                ENDIF
        
        ENDIF

ENDSUBROUTINE BCVELOCITY
!********************************************************************
!*                                                                                                                                        *
!*                                                        BCLEVELSET                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a program for apply BC's to the level-set field                        *                                                                                                                        
!*                                                                                                                                        *
!********************************************************************
SUBROUTINE BCLEVELSET(Nx,Ny,hx,hy,ywhole,phiLS,&
                                          BCphisouth,BCphinorth,BCphiwest,BCphieast,num)
        implicit none
        
        ! local integers
        integer :: i,j

        ! integers
        integer, intent(in) :: Nx,Ny,BCphisouth,BCphinorth,BCphiwest,BCphieast
        
        ! local reals
        real(kind=8) :: deltafunc

        ! reals
        real(kind=8), intent(in) :: hx,hy,num
        
        ! real arrays
        real(kind=8), dimension(Ny+2), intent(in) :: ywhole
        real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phiLS
        
        ! external functions
        real(kind=8), external :: LSDELTA

        

        ! south wall
        IF (BCphisouth == 1) THEN
                phiLS(2:Nx+1,1) = phiLS(2:Nx+1,2)                        ! neumman 
        else
                ! extrapolation
                phiLS(2:Nx+1,1) = 2.0*phiLS(2:Nx+1,2) - phiLS(2:Nx+1,3)
        ENDIF
        
        ! north wall
        IF (BCphinorth == 1) THEN
                phiLS(2:Nx+1,Ny+2) = phiLS(2:Nx+1,Ny+1)                ! neumman 
        else
                ! extrapolation
                phiLS(2:Nx+1,Ny+2) = 2.0*phiLS(2:Nx+1,Ny+1) - phiLS(2:Nx+1,Ny)
        ENDIF

        ! west wall
        IF (BCphiwest == 1) THEN
                phiLS(1,2:Ny+1) = phiLS(2,2:Ny+1)                        ! neumman
        ELSEIF (BCphiwest == 2) THEN
                ! dirclet -> making phiLS = 0 at a point and neumman everywhere else
                DO j=2,Ny+1
                        deltafunc = LSDELTA(Nx,Ny,1,j,hx,hy,phiLS)
                        IF (deltafunc > 0.0) THEN
                                ! I am fixing the interface at ywhole=num
                                phiLS(1,j) = ywhole(j) + num        
                                phiLS(2,j) = ywhole(j) + num
                        else
                                phiLS(1,j) = phiLS(2,j)
                        ENDIF        
                ENDDO                                                                        
        else
                ! extrapolation
                phiLS(1,2:Ny+1) = 2.0*phiLS(2,2:Ny+1) - phiLS(3,2:Ny+1)
        ENDIF

        ! east wall
        IF (BCphieast == 1) THEN
                phiLS(Nx+2,2:Ny+1) = phiLS(Nx+1,2:Ny+1)                ! neumman
        ELSEIF (BCphieast == 2) THEN
                ! dirclet
                DO j=2,Ny+1
                        deltafunc = LSDELTA(Nx,Ny,Nx+2,j,hx,hy,phiLS)
                        IF (deltafunc > 0.0) THEN
                                ! I am fixing the interface at ywhole=num
                                phiLS(Nx+2,j) = ywhole(j) + num
                                phiLS(Nx+1,j) = ywhole(j) + num
                        else
                                phiLS(Nx+2,j) = phiLS(Nx+1,j)
                        ENDIF
                ENDDO
        else
                ! extrapolation
                phiLS(Nx+2,2:Ny+1) = 2.0*phiLS(Nx+1,2:Ny+1) - phiLS(Nx,2:Ny+1)
        ENDIF

        ! corner BC's
        ! phiLS(1,1) = 0.5*(phiLS(2,1) + phiLS(1,2))
        ! phiLS(Nx+2,1) = 0.5*(phiLS(Nx+1,1) + phiLS(Nx+2,2))
        ! phiLS(1,Ny+2) = 0.5*(phiLS(1,Ny+1) + phiLS(2,Ny+2))
        ! phiLS(Nx+2,Ny+2) = 0.5*(phiLS(Nx+1,Ny+2) + phiLS(Nx+2,Ny+1))

        ! MODIFIED BY C. LEE (6/2011)
        phiLS(1,1) = phiLS(2,1) + phiLS(1,2) - phiLS(2,2)
        phiLS(Nx+2,1) = phiLS(Nx+1,1) + phiLS(Nx+2,2) - phiLS(Nx+1,2)
        phiLS(1,Ny+2) = phiLS(1,Ny+1) + phiLS(2,Ny+2) - phiLS(2,Ny+1)
        phiLS(Nx+2,Ny+2) = phiLS(Nx+1,Ny+2) + phiLS(Nx+2,Ny+1) - phiLS(Nx+1,Ny+1)

ENDSUBROUTINE BCLEVELSET
!********************************************************************
!*                                                                                                                                        *
!*                                                        LSCONTACT                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a program for moving the contact line along a boundary        *                                                                                                        
!*                                                                                                                                        *
!********************************************************************
REAL(kind=8) FUNCTION CONTACTSPEED(Nx,Ny,i,j,hx,hy,phiLS,&
                                                                   contC,contA,contR,contM)

        !possible error- check the sign of the angle.

        implicit none
        
        ! integers 
        integer, intent(in):: Nx,Ny,i,j

        real(kind=8) :: phi_x,phi_y,phi_xy,phi_xx,phi_yy,norm
        real(kind=8) :: pi,angle,Ucl,deltafcn,epsilon

        real(kind=8), intent(in) :: hx, hy
        real(kind=8), intent(in) :: contC,contA,contR,contM
        real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phiLS
        
        ! local parameters
        parameter (pi = 3.141592653589)
        
        ! external functions
        real(kind=8), external :: LSDELTA
        !
        !
        !
        epsilon = 1.5*max(hx,hy)  ! smearing the deltafcn 

        IF (j == Ny+2) THEN                        ! north wall
                phi_x = ( phiLS(i+1,j) - phiLS(i-1,j) )/(2.0*hx) ! central difference
                phi_y = ( phiLS(i,j) - phiLS(i,j-1) )/hy                 ! backward difference
        ELSEIF (j == 1) THEN                ! south wall
                phi_x = ( phiLS(i+1,j) - phiLS(i-1,j) )/(2.0*hx) ! central difference
                phi_y = ( phiLS(i,j+1) - phiLS(i,j) )/hy                 ! forward difference
        ELSEIF (i == Nx+2) THEN                ! east wall
                phi_x = ( phiLS(i,j) - phiLS(i-1,j) )/hx                  ! backward difference
                phi_y = ( phiLS(i,j+1) - phiLS(i,j-1) )/(2.0*hy)  ! central difference
        ELSEIF (i == 1) THEN                ! west wall
                phi_x = ( phiLS(i+1,j) - phiLS(i,j) )/hx                  ! forward difference
                phi_y = ( phiLS(i,j+1) - phiLS(i,j-1) )/(2.0*hy)  ! central difference
        ENDIF

        ! normal vector: nx = phi_x/norm, ny = phi_y/norm
        angle = pi + ATAN2(phi_x,phi_y)  ! angle of contact line
                        
        ! delta function
        deltafcn = LSDELTA(Nx,Ny,i,j,hx,hy,phiLS)
                
        if (abs(phi_x) > 2.0 .OR. abs(phi_y) > 2.0) then
            deltafcn = 0.0
        endif
                
        ! contact line speed
        IF (angle > contA) THEN
                CONTACTSPEED = -deltafcn*contC*(angle-contA)**contM  
        ELSEIF (angle < contR) THEN
                CONTACTSPEED = deltafcn*contC*(contR-angle)**contM  
        else
                CONTACTSPEED = 0.0
        ENDIF
                        
ENDFUNCTION CONTACTSPEED
!********************************************************************
!*                                                                                                                                        *
!*                                                          EXPINT                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a function used in find the exponential integral function*
!* EXPINT(x) = integral from x to infinity of exp(-t)/t for x>0                *
!*                                                                                                                                        *
!********************************************************************
REAL(kind=8) FUNCTION EXPINT(x)
        implicit none
        
        integer :: i,Num

        real(kind=8), intent(in) :: x
        real(kind=8) :: Fone, Ftwo,sum,t,delta

        ! using the trapezoidal rule to integrate the function
        Num = 2000
        delta = (10.0-x)/real(Num)
        IF (delta<0.0) THEN
                EXPINT = 0.0
        else
                sum = 0.0
                t=x
                DO i=1,Num
                        Fone = exp(-t)/t
                        Ftwo = exp(-(t+delta))/(t+delta)
                        sum = sum + 0.5*delta*(Fone + Ftwo)
                        t=t+delta
                ENDDO
                EXPINT = sum
        ENDIF
        IF (x==0.0) THEN
                EXPINT = 1.0e16
        ENDIF
        IF (x>10.0) THEN
                EXPINT = 0.0
        ENDIF
ENDFUNCTION EXPINT
!********************************************************************
!*                                                                                                                                        *
!*                                                                  ERF                                                                *
!*                                                                                                                                        *
!********************************************************************
!* Author: Nathaniel R. Morgan                                                                                *
!*                                                                                                                                        *
!* This is a function used in find the error function                                *
!* ERF(x) = 2/sqrt(pi) * integral from 0 to x of exp(-t^2)                         *
!*                                                                                                                                        *
!********************************************************************
REAL(kind=8) FUNCTION ERF(x)
        implicit none
        
        integer :: i,Num

        real(kind=8), intent(in) :: x
        real(kind=8) :: Fone, Ftwo,sum,t,delta,pi
        
        pi = 3.1415926535897932384626433832795

        ! using the trapezoidal rule to integrate the function
        Num = 2000
        delta = x/real(Num)
        IF (delta<=0.0) THEN
                ERF = 0.0
        else
                sum = 0.0
                t=0.0
                DO i=1,Num
                        Fone = 2.0/sqrt(pi)*exp(-t**2.0)
                        Ftwo = 2.0/sqrt(pi)*exp(-(t+delta)**2.0)
                        sum = sum + 0.5*delta*(Fone + Ftwo)
                        t=t+delta
                ENDDO
                ERF = sum
        ENDIF
        
ENDFUNCTION ERF
