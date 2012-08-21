!GALS.F90
! Implements the narrow-band gradient-augmented level set method 
! for advection of the fluid free surface.
! Author: Curtis Lee

SUBROUTINE GALS(Nx,Ny,hx,hy,deltaT,xwhole,ywhole,GA,scal, &
				phi,phi_x,phi_y,phi_xy,prox,x,y,u,v,u_old,v_old,volume,volume0,time,width)

    implicit none
    integer, intent(in) :: Nx, Ny, width
    real(kind=8), intent(in) :: hx, hy, deltaT, time
    real(kind=8), intent(in) :: volume,volume0
    integer, dimension(5,6), intent(in) :: GA
    real(kind=8), dimension(6), intent(in) :: scal
    real(kind=8), dimension(Nx+1), intent(in) :: x
    real(kind=8), dimension(Ny+1), intent(in) :: y
    real(kind=8), dimension(Nx+2), intent(in) :: xwhole
    real(kind=8), dimension(Ny+2), intent(in) :: ywhole
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi       ! level set value
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi_x     ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi_y     ! gradient y-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi_xy    ! cross-derivative of phi
    integer, dimension(Nx+2,Ny+2), intent(inout) :: prox
    real(kind=8), dimension(Nx+1,Ny+2), intent(in) :: u            ! x-vel, ns
    real(kind=8), dimension(Nx+2,Ny+1), intent(in) :: v            ! y-vel, n
    real(kind=8), dimension(Nx+1,Ny+2), intent(in) :: u_old        ! x-vel, n-1
    real(kind=8), dimension(Nx+2,Ny+1), intent(in) :: v_old        ! y-vel, n-1
    
    ! LOCAL VARIABLES
    integer :: i,j
    real(kind=8), dimension(6) :: herm
    real(kind=8), dimension(6) :: vel
    real(kind=8), dimension(Nx+1,Ny+2) :: u_half                   ! x-vel, n-1/2
    real(kind=8), dimension(Nx+2,Ny+1) :: v_half                   ! y-vel, n-1/2
    real(kind=8), dimension(Nx+2,Ny+2) :: phin                     ! temporary values
    real(kind=8), dimension(Nx+2,Ny+2) :: phi_xn    
    real(kind=8), dimension(Nx+2,Ny+2) :: phi_yn     
    real(kind=8), dimension(6) :: vel1,vel2,vel3      ! interpolated velocities
	real(kind=8), dimension(2) :: node,X1,X2,X3
	real(kind=8), dimension(2,2) :: DX1,DX2,DX3
    
    ! Store Old Level Set Values
    DO i=1,Nx+2
        DO j=1,Ny+2
            phin(i,j) = phi(i,j)
            phi_xn(i,j) = phi_x(i,j)
            phi_yn(i,j) = phi_y(i,j)
        ENDDO
    ENDDO
    
	! Calculate Intermediate Velocities
    DO i=1,Nx+1
        DO j=1,Ny+2
            u_half(i,j) = (u(i,j) + u_old(i,j))/2.0
        ENDDO
    ENDDO
    DO i=1,Nx+2
        DO j=1,Ny+1
            v_half(i,j) = (v(i,j) + v_old(i,j))/2.0
        ENDDO
    ENDDO

    ! Superconsistent Shu-Osher RK3 Scheme
    DO i=2,Nx+1
        DO j=2,Ny+1
        
        	if (prox(i,j) > 0 .AND. prox(i,j) <= width-1) then
        
                node(1) = xwhole(i)
                node(2) = ywhole(j)
                  
                ! STEP ONE
                CALL VELO(Nx,Ny,hx,hy,x,y,xwhole,ywhole,u,v,vel1,node,time)
                  
                X1(1) = node(1) - deltaT*vel1(1)
                X1(2) = node(2) - deltaT*vel1(4)
                  
                DX1(1,1) = 1.0 - deltaT*vel1(2)
                DX1(1,2) = 0.0 - deltaT*vel1(5)
                DX1(2,1) = 0.0 - deltaT*vel1(3)
                DX1(2,2) = 1.0 - deltaT*vel1(6)
                                                                                           
                ! STEP TWO
                CALL VELO(Nx,Ny,hx,hy,x,y,xwhole,ywhole,u_old,v_old,vel2,X1,time-deltaT)
                
                X2(1) = node(1) - (deltaT/4.0)*(vel1(1) + vel2(1))
                X2(2) = node(2) - (deltaT/4.0)*(vel1(4) + vel2(4))
                      
                DX2(1,1) = 1.0 - (deltaT/4.0)*(vel1(2) + DX1(1,1)*vel2(2) + DX1(1,2)*vel2(3) )
                DX2(1,2) = 0.0 - (deltaT/4.0)*(vel1(5) + DX1(1,1)*vel2(5) + DX1(1,2)*vel2(6) )
                DX2(2,1) = 0.0 - (deltaT/4.0)*(vel1(3) + DX1(2,1)*vel2(2) + DX1(2,2)*vel2(3) )
                DX2(2,2) = 1.0 - (deltaT/4.0)*(vel1(6) + DX1(2,1)*vel2(5) + DX1(2,2)*vel2(6) )
            
                ! STEP THREE
                CALL VELO(Nx,Ny,hx,hy,x,y,xwhole,ywhole,u_half,v_half,vel3,X2,time-deltaT/2.0)
                      
                X3(1) = node(1) - (deltaT/6.0)*(vel1(1) + vel2(1) + 4.0*vel3(1))
                X3(2) = node(2) - (deltaT/6.0)*(vel1(4) + vel2(4) + 4.0*vel3(4))
                
                DX3(1,1) = 1.0 - (deltaT/6.0)*(vel1(2) + (DX1(1,1)*vel2(2) + DX1(1,2)*vel2(3)) &
                                                 + 4.0 * (DX2(1,1)*vel3(2) + DX2(1,2)*vel3(3)) )
                DX3(1,2) = 0.0 - (deltaT/6.0)*(vel1(5) + (DX1(1,1)*vel2(5) + DX1(1,2)*vel2(6)) &
                                                 + 4.0 * (DX2(1,1)*vel3(5) + DX2(1,2)*vel3(6)) )
                DX3(2,1) = 0.0 - (deltaT/6.0)*(vel1(3) + (DX1(2,1)*vel2(2) + DX1(2,2)*vel2(3)) &
                                                 + 4.0 * (DX2(2,1)*vel3(2) + DX2(2,2)*vel3(3)) )
                DX3(2,2) = 1.0 - (deltaT/6.0)*(vel1(6) + (DX1(2,1)*vel2(5) + DX1(2,2)*vel2(6)) &
                                                 + 4.0 * (DX2(2,1)*vel3(5) + DX2(2,2)*vel3(6)) )

                CALL HERMCUBE(Nx,Ny,hx,hy,xwhole,ywhole,phin,phi_xn,phi_yn,phi_xy,herm, &
                              GA,scal,X3,i,j,3)

                phi(i,j) =   herm(1)
                phi_x(i,j) = DX3(1,1)*herm(2) + DX3(1,2)*herm(3)
                phi_y(i,j) = DX3(2,1)*herm(2) + DX3(2,2)*herm(3)
                                         
            endif
                              
        ENDDO
    ENDDO

	CALL EXTRAP(Nx,Ny,hx,hy,phi,phi_x,phi_y,phi_xy,prox,width)     
    CALL BCGRAD(Nx,Ny,hx,hy,phi,phi_x,phi_y,prox)
    CALL CROSSDERV(Nx,Ny,hx,hy,phi,phi_x,phi_y,phi_xy,prox,width)     

    CALL BCCHECK(Nx,Ny,hx,hy,xwhole,ywhole,GA,scal,phi,phi_x,phi_y,phi_xy,prox,width)    
        
END SUBROUTINE GALS

SUBROUTINE HERMCUBE(Nx,Ny,hx,hy,xwhole,ywhole,phi,phi_x,phi_y,phi_xy, &
                    herm,GA,scal,point,i,j,mode)

    implicit none
    integer, intent(in) :: Nx,Ny,i,j,mode
    integer, dimension(5,6), intent(in) :: GA
    real(kind=8), intent(in) :: hx, hy
    real(kind=8), dimension(2), intent(in) :: point
    real(kind=8), dimension(6), intent(in) :: scal
    real(kind=8), dimension(Nx+2), intent(in) :: xwhole
    real(kind=8), dimension(Ny+2), intent(in) :: ywhole
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi            ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi_x          ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi_y          ! gradient y-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi_xy         ! cross-derivative of phi
    real(kind=8), dimension(6),intent(out) :: herm
    
    ! Local Variables
    integer :: a,b,k
    real(kind=8) :: chi,eta
    real(kind=8), dimension(3,2,2) :: f,g
    
    a = i
    b = j
    IF (point(1) < xwhole(i)) THEN
        a = i - 1
    ENDIF
    IF (point(2) < ywhole(j)) THEN
        b = j - 1    
    ENDIF
    
    a = max(a,1)   
    b = max(b,1)
 
    chi = (point(1) - xwhole(a))/hx
    eta = (point(2) - ywhole(b))/hy
    
    ! BASIS FUNCTION KEY:
    !      {(1,2,3 ),(1,2),(1  ,2  )}
    ! (f,g){(-,p,pp),(0,1),(chi,eta)}
    
    f(1,1,1) = 1.0 - 3.0*     chi **2.0 + 2.0*     chi **3.0
    f(1,2,1) = 1.0 - 3.0*(1.0-chi)**2.0 + 2.0*(1.0-chi)**3.0
    f(1,1,2) = 1.0 - 3.0*     eta **2.0 + 2.0*     eta **3.0
    f(1,2,2) = 1.0 - 3.0*(1.0-eta)**2.0 + 2.0*(1.0-eta)**3.0

    g(1,1,1) =      chi *(1.0 -      chi )**2.0 
    g(1,2,1) = (1.0-chi)*(1.0 - (1.0-chi))**2.0
    g(1,1,2) =      eta *(1.0 -      eta )**2.0
    g(1,2,2) = (1.0-eta)*(1.0 - (1.0-eta))**2.0 

    f(2,1,1) = 6.0*(     chi **2.0 -      chi )
    f(2,2,1) = 6.0*((1.0-chi)**2.0 - (1.0-chi))
    f(2,1,2) = 6.0*(     eta **2.0 -      eta )
    f(2,2,2) = 6.0*((1.0-eta)**2.0 - (1.0-eta))

    g(2,1,1) = 1.0 - 4.0*     chi  + 3.0*     chi **2.0
    g(2,2,1) = 1.0 - 4.0*(1.0-chi) + 3.0*(1.0-chi)**2.0
    g(2,1,2) = 1.0 - 4.0*     eta  + 3.0*     eta **2.0
    g(2,2,2) = 1.0 - 4.0*(1.0-eta) + 3.0*(1.0-eta)**2.0
    
    f(3,1,1) = -6.0 + 12.0*     chi
    f(3,2,1) = -6.0 + 12.0*(1.0-chi)
    f(3,1,2) = -6.0 + 12.0*     eta
    f(3,2,2) = -6.0 + 12.0*(1.0-eta)

    g(3,1,1) = -4.0 + 6.0*     chi
    g(3,2,1) = -4.0 + 6.0*(1.0-chi)
    g(3,1,2) = -4.0 + 6.0*     eta
    g(3,2,2) = -4.0 + 6.0*(1.0-eta)
    
    DO k = 1,mode
    
        herm(k) =      scal(k)*((          phi(a,b)       *f(GA(1,k),1,1)*f(GA(2,k),1,2)  &
      	              	         + GA(3,k)*phi(a,b+1)     *f(GA(1,k),1,1)*f(GA(2,k),2,2)  &
        	                     + GA(4,k)*phi(a+1,b)     *f(GA(1,k),2,1)*f(GA(2,k),1,2)  &
            	                 + GA(5,k)*phi(a+1,b+1)   *f(GA(1,k),2,1)*f(GA(2,k),2,2)) &
                	       + hy*(          phi_y(a,b)     *f(GA(1,k),1,1)*g(GA(2,k),1,2)  &
                    	         - GA(3,k)*phi_y(a,b+1)   *f(GA(1,k),1,1)*g(GA(2,k),2,2)  &
                                 + GA(4,k)*phi_y(a+1,b)   *f(GA(1,k),2,1)*g(GA(2,k),1,2)  &
    	                         - GA(5,k)*phi_y(a+1,b+1) *f(GA(1,k),2,1)*g(GA(2,k),2,2)) &
        	               + hx*(          phi_x(a,b)     *g(GA(1,k),1,1)*f(GA(2,k),1,2)  & 
            	                 + GA(3,k)*phi_x(a,b+1)   *g(GA(1,k),1,1)*f(GA(2,k),2,2)  &
                	             - GA(4,k)*phi_x(a+1,b)   *g(GA(1,k),2,1)*f(GA(2,k),1,2)  &
                    	         - GA(5,k)*phi_x(a+1,b+1) *g(GA(1,k),2,1)*f(GA(2,k),2,2)) &
    	                + hx*hy*(          phi_xy(a,b)    *g(GA(1,k),1,1)*g(GA(2,k),1,2)  &
        	                     - GA(3,k)*phi_xy(a,b+1)  *g(GA(1,k),1,1)*g(GA(2,k),2,2)  &
            	                 - GA(4,k)*phi_xy(a+1,b)  *g(GA(1,k),2,1)*g(GA(2,k),1,2)  &
                	             + GA(5,k)*phi_xy(a+1,b+1)*g(GA(1,k),2,1)*g(GA(2,k),2,2)))
                
    ENDDO       

    
END SUBROUTINE HERMCUBE

SUBROUTINE VELO(Nx,Ny,hx,hy,x,y,xwhole,ywhole,u,v,vel,point,time)

    implicit none
    integer, intent(in) :: Nx, Ny
    real(kind=8), intent(in) :: hx, hy, time
    real(kind=8), dimension(6), intent(out) :: vel
    real(kind=8), dimension(2), intent(in) :: point
    real(kind=8), dimension(Nx+1), intent(in) :: x
    real(kind=8), dimension(Ny+1), intent(in) :: y
    real(kind=8), dimension(Nx+2), intent(in) :: xwhole
    real(kind=8), dimension(Ny+2), intent(in) :: ywhole
    real(kind=8), dimension(Nx+1,Ny+2), intent(in) :: u
    real(kind=8), dimension(Nx+2,Ny+1), intent(in) :: v
	
    ! Local Variables
    real(kind=8) :: pi                                         
	     
    ! Local Variables
    integer :: a, b
    real(kind=8), dimension(2) :: L1, L2
    
    pi = 3.1415926535897932384626433832795
    
    ! Locate bottom-left node for u
    a = NINT((point(1)-MOD(point(1),hx))/hx) + 1
    b = NINT((point(2)+(hy/2.0)-MOD(point(2)+(hy/2.0),hy))/hy) + 1
    IF (a < 1) THEN
        a = 1
    ELSEIF (a > Nx+1) THEN
        a = Nx+1
    ENDIF
    IF (b < 2) THEN
        b = 2
    ELSEIF (b > Ny) THEN
        b = Ny
    ENDIF
                     
    L1(1) =  (point(1)-x(a+1))/hx
    L1(2) = -(point(1)-x(a))/hx
    L2(1) =  (point(2)-ywhole(b+1))/hy
    L2(2) = -(point(2)-ywhole(b))/hy        
    vel(1)=   u(a,b)     * L1(1) * L2(1) &
            + u(a+1,b)   * L1(2) * L2(1) &
            + u(a,b+1)   * L1(1) * L2(2) &
            + u(a+1,b+1) * L1(2) * L2(2)
            
    !print *, 'uL1: ', -L1(2)
    !print *, 'uL2: ', -L2(2)
                     
    L1(1) =  1.0/hx
    L1(2) = -1.0/hx
    L2(1) =  (point(2)-ywhole(b+1))/hy
    L2(2) = -(point(2)-ywhole(b))/hy        
    vel(2)=   u(a,b)     * L1(1) * L2(1) &
            + u(a+1,b)   * L1(2) * L2(1) &
            + u(a,b+1)   * L1(1) * L2(2) &
            + u(a+1,b+1) * L1(2) * L2(2)
                     
    L1(1) =  (point(1)-x(a+1))/hx
    L1(2) = -(point(1)-x(a))/hx
    L2(1) =  1.0/hy
    L2(2) = -1.0/hy              
	vel(3)=   u(a,b)     * L1(1) * L2(1) &
            + u(a+1,b)   * L1(2) * L2(1) &
            + u(a,b+1)   * L1(1) * L2(2) &
            + u(a+1,b+1) * L1(2) * L2(2)
            
                     
    ! Locate bottom-left node for v         
    a = NINT((point(1)+(hx/2.0)-MOD(point(1)+(hx/2.0),hx))/hx) + 1
    b = NINT((point(2)-MOD(point(2),hy))/hy) + 1
    IF (a < 2) THEN
        a = 2
    ELSEIF (a > Nx) THEN
        a = Nx
    ENDIF       
    IF (b < 1) THEN
        b = 1
    ELSEIF (b > Ny+1) THEN
        b = Ny+1
    ENDIF
    
    L1(1) =  (point(1)-xwhole(a+1))/hx
    L1(2) = -(point(1)-xwhole(a))/hx
    L2(1) =  (point(2)-y(b+1))/hy
    L2(2) = -(point(2)-y(b))/hy        
    vel(4)=   v(a,b)     * L1(1) * L2(1) &
            + v(a+1,b)   * L1(2) * L2(1) &
            + v(a,b+1)   * L1(1) * L2(2) &
            + v(a+1,b+1) * L1(2) * L2(2) 
            
            
    !print *, 'vL1: ', -L1(2)
    !print *, 'vL2: ', -L2(2)
                     
    L1(1) =  1.0/hx
    L1(2) = -1.0/hx
    L2(1) =  (point(2)-y(b+1))/hy
    L2(2) = -(point(2)-y(b))/hy
    vel(5) =  v(a,b)     * L1(1) * L2(1) &
            + v(a+1,b)   * L1(2) * L2(1) &
            + v(a,b+1)   * L1(1) * L2(2) &
            + v(a+1,b+1) * L1(2) * L2(2) 
                     
    L1(1) =  (point(1)-xwhole(a+1))/hx
    L1(2) = -(point(1)-xwhole(a))/hx
    L2(1) =  1.0/hy
    L2(2) = -1.0/hy
    vel(6) =  v(a,b)     * L1(1) * L2(1) &
            + v(a+1,b)   * L1(2) * L2(1) &
            + v(a,b+1)   * L1(1) * L2(2) &
            + v(a+1,b+1) * L1(2) * L2(2) 

     
    !vel(1) = -        cos(pi*time/2.0)*(sin(pi*point(1))**2.0)*sin(2.0*pi*point(2))
    !vel(2) = -     pi*cos(pi*time/2.0)*sin(2.0*pi*point(1))*sin(2.0*pi*point(2))
    !vel(3) = - 2.0*pi*cos(pi*time/2.0)*(sin(pi*point(1))**2.0)*cos(2.0*pi*point(2))       
    !vel(4) =        cos(pi*time/2.0)*(sin(pi*point(2))**2.0)*sin(2.0*pi*point(1))
    !vel(5) = 2.0*pi*cos(pi*time/2.0)*(sin(pi*point(2))**2.0)*cos(2.0*pi*point(1))
    !vel(6) =     pi*cos(pi*time/2.0)*sin(2.0*pi*point(1))*sin(2.0*pi*point(2))
    
END SUBROUTINE VELO

SUBROUTINE VELO2(Nx,Ny,hx,hy,x,y,xwhole,ywhole,u,v,vel,point,time)

    implicit none
    integer, intent(in) :: Nx, Ny
    real(kind=8), intent(in) :: hx, hy, time
    real(kind=8), dimension(6), intent(out) :: vel
    real(kind=8), dimension(2), intent(in) :: point
    real(kind=8), dimension(Nx+1), intent(in) :: x
    real(kind=8), dimension(Ny+1), intent(in) :: y
    real(kind=8), dimension(Nx+2), intent(in) :: xwhole
    real(kind=8), dimension(Ny+2), intent(in) :: ywhole
    real(kind=8), dimension(Nx+1,Ny+2), intent(in) :: u
    real(kind=8), dimension(Nx+2,Ny+1), intent(in) :: v
    
    ! Local Variables
    integer :: a, b
    real(kind=8), dimension(3) :: L1, L2, L3
    
    ! Locate bottom-left node for u
    a = NINT((point(1)-MOD(point(1),hx))/hx) + 1
    b = NINT((point(2)+(hy/2.0)-MOD(point(2)+(hy/2.0),hy))/hy) + 1
    IF (a < 1) THEN
        a = 1
    ELSEIF (a > Nx-1) THEN
        a = Nx-1
    ENDIF
    IF (b < 2) THEN
        b = 2
    ELSEIF (b > Ny-2) THEN
        b = Ny-2
    ENDIF
                     
    L1(1) =   (1.0/(2.0*hx**2))*((point(1)-x(a+1))*(point(1)-x(a+2)))
    L1(2) =  -(1.0/(    hx**2))*((point(1)-x(a))*(point(1)-x(a+2)))
    L1(3) =   (1.0/(2.0*hx**2))*((point(1)-x(a))*(point(1)-x(a+1))) 
    L2(1) =   (1.0/(2.0*hy**2))*((point(2)-ywhole(b+1))*(point(2)-ywhole(b+2)))
    L2(2) =  -(1.0/(    hy**2))*((point(2)-ywhole(b))*(point(2)-ywhole(b+2)))
    L2(3) =   (1.0/(2.0*hy**2))*((point(2)-ywhole(b))*(point(2)-ywhole(b+1)))      
    vel(1)=   u(a,  b  ) * L1(1) * L2(1) &
            + u(a+1,b  ) * L1(2) * L2(1) &
            + u(a+2,b  ) * L1(3) * L2(1) &
			+ u(a,  b+1) * L1(1) * L2(2) &
            + u(a+1,b+1) * L1(2) * L2(2) &
			+ u(a+2,b+1) * L1(3) * L2(2) &
            + u(a,  b+2) * L1(1) * L2(3) &
			+ u(a+1,b+2) * L1(2) * L2(3) &
            + u(a+2,b+2) * L1(3) * L2(3) 
                     
    L1(1) =   (1.0/(2.0*hx**2))*((point(1)-x(a+1))+(point(1)-x(a+2)))
    L1(2) =  -(1.0/(    hx**2))*((point(1)-x(  a))+(point(1)-x(a+2)))
    L1(3) =   (1.0/(2.0*hx**2))*((point(1)-x(  a))+(point(1)-x(a+1)))
    L2(1) =   (1.0/(2.0*hy**2))*((point(2)-ywhole(b+1))*(point(2)-ywhole(b+2)))
    L2(2) =  -(1.0/(    hy**2))*((point(2)-ywhole(b))*(point(2)-ywhole(b+2)))
    L2(3) =   (1.0/(2.0*hy**2))*((point(2)-ywhole(b))*(point(2)-ywhole(b+1)))     
    vel(2)=   u(a,  b  ) * L1(1) * L2(1) &
            + u(a+1,b  ) * L1(2) * L2(1) &
            + u(a+2,b  ) * L1(3) * L2(1) &
            + u(a,  b+1) * L1(1) * L2(2) &
			+ u(a+1,b+1) * L1(2) * L2(2) &
            + u(a+2,b+1) * L1(3) * L2(2) &
			+ u(a,  b+2) * L1(1) * L2(3) &
            + u(a+1,b+2) * L1(2) * L2(3) &
			+ u(a+2,b+2) * L1(3) * L2(3) 
                     
    L1(1) =   (1.0/(2.0*hx**2))*((point(1)-x(a+1))*(point(1)-x(a+2)))
    L1(2) =  -(1.0/(    hx**2))*((point(1)-x(a))*(point(1)-x(a+2)))
    L1(3) =   (1.0/(2.0*hx**2))*((point(1)-x(a))*(point(1)-x(a+1)))		  
    L2(1) =   (1.0/(2.0*hy**2))*((point(2)-ywhole(b+1))+(point(2)-ywhole(b+2)))
    L2(2) =  -(1.0/(    hy**2))*((point(2)-ywhole(  b))+(point(2)-ywhole(b+2)))
    L2(3) =   (1.0/(2.0*hy**2))*((point(2)-ywhole(  b))+(point(2)-ywhole(b+1)))       
	vel(3)=   u(a,  b  ) * L1(1) * L2(1) &
            + u(a+1,b  ) * L1(2) * L2(1) &
            + u(a+2,b  ) * L1(3) * L2(1) &
            + u(a,  b+1) * L1(1) * L2(2) &
			+ u(a+1,b+1) * L1(2) * L2(2) &
            + u(a+2,b+1) * L1(3) * L2(2) &
			+ u(a,  b+2) * L1(1) * L2(3) &
            + u(a+1,b+2) * L1(2) * L2(3) &
			+ u(a+2,b+2) * L1(3) * L2(3) 
                     
    ! Locate bottom-left node for v         
    a = NINT((point(1)+(hx/2.0)-MOD(point(1)+(hx/2.0),hx))/hx) + 1
    b = NINT((point(2)-MOD(point(2),hy))/hy) + 1
    IF (a < 2) THEN
        a = 2
    ELSEIF (a > Nx-2) THEN
        a = Nx-2
    ENDIF       
    IF (b < 1) THEN
        b = 1
    ELSEIF (b > Ny-1) THEN
        b = Ny-1
    ENDIF
    
    L1(1) =   (1.0/(2.0*hx**2))*((point(1)-xwhole(a+1))*(point(1)-xwhole(a+2)))
    L1(2) =  -(1.0/(    hx**2))*((point(1)-xwhole(a))*(point(1)-xwhole(a+2)))
    L1(3) =   (1.0/(2.0*hx**2))*((point(1)-xwhole(a))*(point(1)-xwhole(a+1)))			  
    L2(1) =   (1.0/(2.0*hy**2))*((point(2)-y(b+1))*(point(2)-y(b+2)))
    L2(2) =  -(1.0/(    hy**2))*((point(2)-y(b))*(point(2)-y(b+2)))
    L2(3) =   (1.0/(2.0*hy**2))*((point(2)-y(b))*(point(2)-y(b+1)))      
    vel(4)=   v(a,  b  ) * L1(1) * L2(1) &
                        + v(a+1,b  ) * L1(2) * L2(1) &
                        + v(a+2,b  ) * L1(3) * L2(1) &
			            + v(a,  b+1) * L1(1) * L2(2) &
			            + v(a+1,b+1) * L1(2) * L2(2) &
			            + v(a+2,b+1) * L1(3) * L2(2) &
			            + v(a,  b+2) * L1(1) * L2(3) &
			            + v(a+1,b+2) * L1(2) * L2(3) &
			            + v(a+2,b+2) * L1(3) * L2(3) 
                     
    L1(1) =   (1.0/(2.0*hx**2))*((point(1)-xwhole(a+1))+(point(1)-xwhole(a+2)))
    L1(2) =  -(1.0/(    hx**2))*((point(1)-xwhole(  a))+(point(1)-xwhole(a+2)))
    L1(3) =   (1.0/(2.0*hx**2))*((point(1)-xwhole(  a))+(point(1)-xwhole(a+1)))
    L2(1) =   (1.0/(2.0*hy**2))*((point(2)-y(b+1))*(point(2)-y(b+2)))
    L2(2) =  -(1.0/(    hy**2))*((point(2)-y(b))*(point(2)-y(b+2)))
    L2(3) =   (1.0/(2.0*hy**2))*((point(2)-y(b))*(point(2)-y(b+1)))
    vel(5) =  v(a,  b  ) * L1(1) * L2(1) &
                        + v(a+1,b  ) * L1(2) * L2(1) &
                        + v(a+2,b  ) * L1(3) * L2(1) &
			            + v(a,  b+1) * L1(1) * L2(2) &
			            + v(a+1,b+1) * L1(2) * L2(2) &
			            + v(a+2,b+1) * L1(3) * L2(2) &
			            + v(a,  b+2) * L1(1) * L2(3) &
			            + v(a+1,b+2) * L1(2) * L2(3) &
			            + v(a+2,b+2) * L1(3) * L2(3) 
                     
    L1(1) =   (1.0/(2.0*hx**2))*((point(1)-xwhole(a+1))*(point(1)-xwhole(a+2)))
    L1(2) =  -(1.0/(    hx**2))*((point(1)-xwhole(a))*(point(1)-xwhole(a+2)))
    L1(3) =   (1.0/(2.0*hx**2))*((point(1)-xwhole(a))*(point(1)-xwhole(a+1)))  
    L2(1) =   (1.0/(2.0*hy**2))*((point(2)-y(b+1))+(point(2)-y(b+2)))
    L2(2) =  -(1.0/(    hy**2))*((point(2)-y(  b))+(point(2)-y(b+2)))
    L2(3) =   (1.0/(2.0*hy**2))*((point(2)-y(  b))+(point(2)-y(b+1)))
    vel(6) =  v(a,  b  ) * L1(1) * L2(1) &
                        + v(a+1,b  ) * L1(2) * L2(1) &
                        + v(a+2,b  ) * L1(3) * L2(1) &
			            + v(a,  b+1) * L1(1) * L2(2) &
			            + v(a+1,b+1) * L1(2) * L2(2) &
			            + v(a+2,b+1) * L1(3) * L2(2) &
			            + v(a,  b+2) * L1(1) * L2(3) &
			            + v(a+1,b+2) * L1(2) * L2(3) &
			            + v(a+2,b+2) * L1(3) * L2(3) 
    
END SUBROUTINE VELO2

SUBROUTINE VELOEXACT(Nx,Ny,hx,hy,x,y,xwhole,ywhole,u,v,vel,point,time)

    implicit none
    integer, intent(in) :: Nx, Ny
    real(kind=8), intent(in) :: hx, hy, time
    real(kind=8), dimension(6), intent(out) :: vel
    real(kind=8), dimension(2), intent(in) :: point
    real(kind=8), dimension(Nx+1), intent(in) :: x
    real(kind=8), dimension(Ny+1), intent(in) :: y
    real(kind=8), dimension(Nx+2), intent(in) :: xwhole
    real(kind=8), dimension(Ny+2), intent(in) :: ywhole
    real(kind=8), dimension(Nx+1,Ny+2), intent(in) :: u
    real(kind=8), dimension(Nx+2,Ny+1), intent(in) :: v
	
    ! Local Variables
    real(kind=8) :: pi                                         
	
    pi = 3.1415926535897932384626433832795
                     
    vel(1) = -        cos(pi*time/2.0)*(sin(pi*point(1))**2.0)*sin(2.0*pi*point(2))
    vel(2) = -     pi*cos(pi*time/2.0)*sin(2.0*pi*point(1))*sin(2.0*pi*point(2))
    vel(3) = - 2.0*pi*cos(pi*time/2.0)*(sin(pi*point(1))**2.0)*cos(2.0*pi*point(2))
            
    vel(4) =        cos(pi*time/2.0)*(sin(pi*point(2))**2.0)*sin(2.0*pi*point(1))
    vel(5) = 2.0*pi*cos(pi*time/2.0)*(sin(pi*point(2))**2.0)*cos(2.0*pi*point(1))
    vel(6) =     pi*cos(pi*time/2.0)*sin(2.0*pi*point(1))*sin(2.0*pi*point(2))
    
END SUBROUTINE VELOEXACT

SUBROUTINE CROSSDERV(Nx,Ny,hx,hy,phi,phi_x,phi_y,phi_xy,prox,width)

	IMPLICIT NONE
	
    integer, intent(in) :: Nx,Ny,width
    real(kind=8), intent(in) :: hx,hy
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi        ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi_x      ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi_y      ! gradient y-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi_xy  ! cross-derivative of phi
    integer, dimension(Nx+2,Ny+2), intent(in) :: prox
    
    integer :: i,j
    
    ! Solve for Cross-Derivatives
    DO i = 2,Nx+1
        DO j = 2,Ny+1
        	if (prox(i,j) > 0 .AND. prox(i,j) < width) then
                phi_xy(i,j) =   (phi(i+1,j+1) - phi(i+1,j-1) - &
                				 phi(i-1,j+1) + phi(i-1,j-1))/(4.0*hx*hy)
                !phi_xy(i,j) =   ((phi_x(i,j+1) - phi_x(i,j-1))/(2.0*hy) + &
                !				(phi_y(i+1,j) - phi_y(i-1,j))/(2.0*hx))/2.0
            else
            	phi_xy(i,j) = 0.0
            endif
        ENDDO
    ENDDO    
        
    phi_xy(1,:) = phi_xy(2,:)     
	phi_xy(Nx+2,:) = phi_xy(Nx+1,:) 
    phi_xy(:,1) = phi_xy(:,2)        
    phi_xy(:,Ny+2) = phi_xy(:,Ny+1)       
        
END SUBROUTINE CROSSDERV

SUBROUTINE EXTRAP(Nx,Ny,hx,hy,phi,phi_x,phi_y,phi_xy,prox,width)

	IMPLICIT NONE
	
    integer, intent(in) :: Nx,Ny,width
    real(kind=8), intent(in) :: hx,hy
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi        ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi_x      ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi_y      ! gradient y-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi_xy        ! cross-derivative of phi
    integer, dimension(Nx+2,Ny+2), intent(in) :: prox
    
    integer :: i,j,k,l,count
    real(kind=8) :: phisum,phixsum,phiysum,phixysum
     
    DO i = 2,Nx+1
        DO j = 2,Ny+1
			if (prox(i,j) == width) then
							
				count = 0
				phisum = 0.0
				phixsum = 0.0
				phiysum = 0.0
				phixysum = 0.0
							
				DO k = -1,1
					DO l = -1,1
								
						if (prox(i+k,j+l) > 0 .AND. prox(i+k,j+l) < width .AND. &
						    (i+k > 1) .AND. (i+k < Nx+2) .AND. &
						    (j+l > 1) .AND. (j+l < Ny+2)) then
						    
							count = count + 1
							phisum = phisum + 2.0*phi(i+k,j+l) - phi(i+2*k,j+2*l)  
							phixsum = phixsum + 2.0*phi_x(i+k,j+l) - phi_x(i+2*k,j+2*l) 
							phiysum = phiysum + 2.0*phi_y(i+k,j+l) - phi_y(i+2*k,j+2*l)
							phixysum = phixysum + 2.0*phi_xy(i+k,j+l) - phi_xy(i+2*k,j+2*l) 
							
						endif
						
					ENDDO
				ENDDO
				
				phi(i,j) = phisum/count
				phi_x(i,j) = phixsum/count
				phi_y(i,j) = phiysum/count
				phi_xy(i,j) = phixysum/count
										
			endif
        ENDDO
    ENDDO         
     
END SUBROUTINE EXTRAP

SUBROUTINE BCCHECK(Nx,Ny,hx,hy,xwhole,ywhole,GA,scal,phi,phi_x,phi_y,phi_xy,prox,width)
    implicit none
    integer, intent(in) :: Nx,Ny,width
    integer, dimension(5,6), intent(in) :: GA
    real(kind=8), dimension(6), intent(in) :: scal
    real(kind=8), intent(in) :: hx,hy
    real(kind=8), dimension(Nx+2), intent(in) :: xwhole
    real(kind=8), dimension(Ny+2), intent(in) :: ywhole
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi            ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi_x          ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi_y          ! gradient y-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi_xy         ! cross-derivative of phi
    integer, dimension(Nx+2,Ny+2), intent(inout) :: prox

	! Local
	integer :: i,j,k
	real(kind=8), dimension(2) :: trial
    integer, dimension(Nx+2,Ny+2) :: prox_old
    real(kind=8), dimension(6) :: herm1,herm2,herm3,herm4,herm5,herm6,herm7,herm8
    logical :: trigger

    trigger = .FALSE.

    DO i=2,Nx+1
        DO j=2,Ny+1

            if (prox(i,j) == width-1 .AND. phi(i,j) < max(hx,hy)) then 

                trigger = .TRUE.
                continue

            endif
            
        ENDDO
    ENDDO

    if (trigger) then
        CALL REBUILD(Nx,Ny,hx,hy,xwhole,ywhole,GA,scal,phi,phi_x,phi_y,phi_xy,prox,width)
    endif

END SUBROUTINE BCCHECK

SUBROUTINE REBUILD(Nx,Ny,hx,hy,xwhole,ywhole,GA,scal,phi,phi_x,phi_y,phi_xy,prox,width)
    implicit none
    integer, intent(in) :: Nx,Ny,width
    integer, dimension(5,6), intent(in) :: GA
    real(kind=8), dimension(6), intent(in) :: scal
    real(kind=8), intent(in) :: hx,hy
    real(kind=8), dimension(Nx+2), intent(in) :: xwhole
    real(kind=8), dimension(Ny+2), intent(in) :: ywhole
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi            ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi_x          ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi_y          ! gradient y-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi_xy         ! cross-derivative of phi
    integer, dimension(Nx+2,Ny+2), intent(inout) :: prox

	! Local
    real(kind=8), dimension(Nx+2,Ny+2) :: phin            ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2) :: phi_xn          ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2) :: phi_yn          ! gradient y-comp
	integer :: i,j
	
	CALL PROXY(Nx,Ny,hx,hy,xwhole,ywhole,GA,scal,phi,phi_x,phi_y,phi_xy,prox,width)	

    CALL REINIT(Nx,Ny,hx,hy,xwhole,ywhole,GA,scal,phi,phi_x,phi_y,phi_xy,prox,width)
    CALL BCGRAD(Nx,Ny,hx,hy,phi,phi_x,phi_y,prox)
    CALL CROSSDERV(Nx,Ny,hx,hy,phi,phi_x,phi_y,phi_xy,prox,width)     
				
	DO i=1,Nx+2
		DO j=1,Ny+2
        	if(prox(i,j) == -1) then
        		phi(i,j) = 1.0
        		phi_x(i,j) = 0.0
        		phi_y(i,j) = 0.0
        	elseif(prox(i,j) == 0) then
        		phi(i,j) = -1.0
        		phi_x(i,j) = 0.0
        		phi_y(i,j) = 0.0
        	endif
        ENDDO
    ENDDO

	print *, '************************************************'
	print *, '************ REBUILDING NARROW BAND ************'
	print *, '************************************************'

END SUBROUTINE REBUILD

SUBROUTINE PROXY(Nx,Ny,hx,hy,xwhole,ywhole,GA,scal,phi,phi_x,phi_y,phi_xy,prox,width)
    implicit none
    integer, intent(in) :: Nx,Ny,width
    integer, dimension(5,6), intent(in) :: GA
    real(kind=8), dimension(6), intent(in) :: scal
    real(kind=8), intent(in) :: hx,hy
    real(kind=8), dimension(Nx+2), intent(in) :: xwhole
    real(kind=8), dimension(Ny+2), intent(in) :: ywhole
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi            ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi_x          ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi_y          ! gradient y-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi_xy         ! cross-derivative of phi
    integer, dimension(Nx+2,Ny+2), intent(inout) :: prox

	! Local
	integer :: i,j,k
	real(kind=8), dimension(2) :: trial
    integer, dimension(Nx+2,Ny+2) :: prox_old
    real(kind=8), dimension(6) :: herm1,herm2,herm3,herm4,herm5,herm6,herm7,herm8

	prox_old = prox
    prox = 0
    
    ! Nodes nearest to interface (1)
    DO i = 2,Nx+1
        DO j = 2,Ny+1
                            
			if (abs(phi(i,j)) < max(hx,hy)) then
                prox(i,j) = 1
            endif
                        
    	ENDDO
    ENDDO

    ! Nodes in the narrow band (1)
	DO k=2,width
    	DO i = 2,Nx+1
    	    DO j = 2,Ny+1
            
    	        if (prox(i,j) == 0 .AND. &
    	           (prox(i,  j+1) == k-1 .OR. &
    	            prox(i,  j-1) == k-1 .OR. &
            	    prox(i+1,j+1) == k-1 .OR. &
            	    prox(i+1,j)   == k-1 .OR. &
            	    prox(i+1,j-1) == k-1 .OR. &
            	    prox(i-1,j+1) == k-1 .OR. &
            	    prox(i-1,j)   == k-1 .OR. &
                	prox(i-1,j-1) == k-1)) then
           
                	prox(i,j) = k
           
            	endif
            
        	ENDDO
    	ENDDO
    ENDDO

    ! Far-away nodes having positive Level Set Value (-1)
    DO i = 1,Nx+2
        DO j = 1,Ny+2
            
            if (prox(i,j) == 0 .AND. (phi(i,j) > 0.0)) then
                prox(i,j) = -1
            endif
            
        ENDDO
    ENDDO
    
    ! Copy to boundaries appropriately
    DO i = 2,Nx+1
    	if (prox(i-1,2) > 0 .AND. prox(i-1,2) < width .OR. &
    	    prox(i  ,2) > 0 .AND. prox(i  ,2) < width .OR. &
    	    prox(i+1,2) > 0 .AND. prox(i+1,2) < width) then
        		prox(i,1) = prox(i,2)
        endif
    	if (prox(i-1,Ny+1) > 0 .AND. prox(i-1,Ny+1) < width .OR. &
    	    prox(i  ,Ny+1) > 0 .AND. prox(i  ,Ny+1) < width .OR. &
    	    prox(i+1,Ny+1) > 0 .AND. prox(i+1,Ny+1) < width) then
        		prox(i,Ny+2) = prox(i,Ny+1)
        endif      
    ENDDO
    
    DO j = 2,Ny+1
    	if (prox(2,j-1) > 0 .AND. prox(2,j-1) < width .OR. &
    	    prox(2,j  ) > 0 .AND. prox(2,j  ) < width .OR. &
    	    prox(2,j+1) > 0 .AND. prox(2,j+1) < width) then
        		prox(1,j) = prox(2,j)
        endif
    	if (prox(Nx+1,j-1) > 0 .AND. prox(Nx+1,j-1) < width .OR. &
    	    prox(Nx+1,j  ) > 0 .AND. prox(Nx+1,j  ) < width .OR. &
    	    prox(Nx+1,j+1) > 0 .AND. prox(Nx+1,j+1) < width) then
        		prox(Nx+2,j) = prox(Nx+1,j) 
        endif
    ENDDO
        
END SUBROUTINE PROXY

SUBROUTINE VOLUMESOLVE(Nx,Ny,hx,hy,xwhole,ywhole,GA,scal,phi,phi_x,phi_y,phi_xy,prox,volume)
    implicit none

    integer, intent(in) :: Nx,Ny
    integer, dimension(5,6), intent(in) :: GA
    real(kind=8), dimension(6), intent(in) :: scal
    real(kind=8), intent(in) :: hx,hy
    real(kind=8), dimension(Nx+2), intent(in) :: xwhole
    real(kind=8), dimension(Ny+2), intent(in) :: ywhole
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi        ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi_x      ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi_y      ! gradient y-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi_xy        ! cross-derivative of phi
    integer, dimension(Nx+2,Ny+2), intent(in) :: prox
    real(kind=8), intent(out) :: volume

    integer :: i,j,k,l
    integer :: k0,k1,l0,l1
    real(kind=8), dimension(2) :: X
    real(kind=8), dimension(6) :: herm
    
    volume = 0.0
    
    ! Whole cells
    DO i = 1,Nx+1
        DO j = 1,Ny+1           
            if (phi(i,j)>=0.0 .AND. phi(i+1,j)>=0.0 .AND. phi(i,j+1)>=0.0 .AND. phi(i+1,j+1)>=0.0) then
                
                ! Whole cell
                if ((i > 1) .AND. (j > 1) .AND. (i < Nx+1) .AND. (j < Ny+1)) then
                    volume = volume + hx*hy
                ! Edge cell
                elseif (((i > 1) .AND. (i < Nx+1)) .OR. ((j > 1) .AND. (j < Ny+1))) then
                    volume = volume + (hx*hy)/2.0
                ! Corner cell
                else
                    volume = volume + (hx*hy)/4.0
                endif
                
            elseif (phi(i,j)>=0.0 .OR. phi(i+1,j)>=0.0 .OR. phi(i,j+1)>=0.0 .OR. phi(i+1,j+1)>=0.0) then
                
                if (i > 1) then
                    k0 = 1
                else
                    k0 = 6
                endif
                if (j > 1) then
                    l0 = 1
                else
                    l0 = 6
                endif
                
                if (i < Nx+1) then
                    k1 = 10
                else
                    k1 = 5
                endif
                if (j < Ny+1) then
                    l1 = 10
                else
                    l1 = 5
                endif

                DO k = k0,k1
                    DO l = l0,l1
                        X(1) = xwhole(i)+(2.0*k-1.0)*0.05*hx
                        X(2) = ywhole(j)+(2.0*l-1.0)*0.05*hy
                        CALL HERMCUBE(Nx,Ny,hx,hy,xwhole,ywhole,phi,phi_x,phi_y,phi_xy,herm, &
                                      GA,scal,X,i,j,1)
                        if (herm(1) >= 0.0) then
                            volume = volume + (hx*hy)/100.0
                        endif
                    ENDDO
                ENDDO
            endif
        ENDDO
    ENDDO
    
END SUBROUTINE VOLUMESOLVE

SUBROUTINE VOLUMEFIX(Nx,Ny,hx,hy,phi,prox,volume,volume0)
    implicit none

    integer, intent(in) :: Nx,Ny
    real(kind=8), intent(in) :: hx,hy
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi   
    integer, dimension(Nx+2,Ny+2), intent(in) :: prox
    real(kind=8), intent(in) :: volume,volume0

    integer :: i,j,k,l
    real(kind=8) :: kappa

    kappa = 0.05*sqrt(hx*hx+hy*hy)

    DO i=1,Nx+2
        DO j=1,Ny+2

            phi(i,j) = phi(i,j) + kappa*(volume0 - volume)
            
        ENDDO
    ENDDO

END SUBROUTINE VOLUMEFIX

!********************************************************************
!*                                                                                                                                        *
!*                                                             BCGRAD                                                               *
!*                                                                                                                                        *
!********************************************************************
!* Author: Curtis Lee                                                                                *
!*                                                                                                                                        *
!* This is a program for applying BC's to the level set and its 
!* gradients (G.A.L.S.)                                                                                                 
!*                                                                                                                                        *
!********************************************************************
SUBROUTINE BCGRAD(Nx,Ny,hx,hy,phiLS,phi_x,phi_y,prox)

      implicit none
        
      integer :: i,j
      integer, intent(in) :: Nx,Ny
      real(kind=8), intent(in) :: hx,hy
      real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phiLS,phi_x,phi_y
      integer, dimension(Nx+2,Ny+2), intent(in) :: prox
            
      ! Extrapolate level set values for dambreak/contact-angle problems!
            
        DO i = 2,Nx+1
        	if (prox(i,1) > 0) then
          		!phiLS(i,1) = 3*phiLS(i,2) - 3*phiLS(i,3) + phiLS(i,4)
          		phiLS(i,1) = phiLS(i,2)
          		phi_x(i,1) = phi_x(i,2)
          		phi_y(i,1) = 0.0
          		!phiLS(i,1) = -4*phiLS(i,2) + 5*phiLS(i,3) - hy*(4*phi_y(i,2) + 2*phi_y(i,3)) 
          	endif
        	if (prox(i,Ny+2) > 0) then
          		!phiLS(i,Ny+2) = 3*phiLS(i,Ny+1) - 3*phiLS(i,Ny) + phiLS(i,Ny-1)
          		phiLS(i,Ny+2) = phiLS(i,Ny+1) 
          		phi_x(i,Ny+2) = phi_x(i,Ny+1)
          		phi_y(i,Ny+2) = 0.0
          		!phiLS(i,Ny+2) = 5*phiLS(i,Ny) - 4*phiLS(i,Ny+1) + hy*(2*phi_y(i,Ny) + 4*phi_y(i,Ny+1))
          	endif
        ENDDO
        DO j = 2,Ny+1
        	if (prox(1,j) > 0) then
          		!phiLS(1,j) = 3*phiLS(2,j) - 3*phiLS(3,j) + phiLS(4,j)
          		phiLS(1,j) = phiLS(2,j) 
          		phi_x(1,j) = 0.0
          		phi_y(1,j) = phi_y(2,j)
          		!phiLS(1,j) = phiLS(2,j) - hx*phi_x(2,j)
          		!phiLS(1,j) = -4*phiLS(2,j) + 5*phiLS(3,j) - hy*(4*phi_y(2,j) + 2*phi_y(3,j)) 
          	endif
        	if (prox(Nx+2,j) > 0) then
          		!phiLS(Nx+2,j) = 3*phiLS(Nx+1,j) - 3*phiLS(Nx,j) + phiLS(Nx-1,j)
          		phiLS(Nx+2,j) = phiLS(Nx+1,j)
          		phi_x(Nx+2,j) = 0.0
          		phi_y(Nx+2,j) = phi_y(Nx+1,j)
          		!phiLS(Nx+2,j) = phiLS(Nx+1,j) + hx*phi_x(Nx+1,j)
          		!phiLS(Nx+2,j) = 5*phiLS(Nx,j) - 4*phiLS(Nx+1,j) + hy*(2*phi_y(Nx,j) + 4*phi_y(Nx+1,j))
          	endif
        ENDDO
                    
        
        phiLS(1,1) = phiLS(2,2)
        phiLS(Nx+2,1) = phiLS(Nx+1,2)
        phiLS(1,Ny+2) = phiLS(2,Ny+1)
        phiLS(Nx+2,Ny+2) = phiLS(Nx+1,Ny+1)
        phi_x(1,1) = (phi_x(2,1) + phi_x(1,2))/2.0
        phi_x(Nx+2,1) = (phi_x(Nx+1,1) + phi_x(Nx+2,2))/2.0
        phi_x(1,Ny+2) = (phi_x(1,Ny+1) + phi_x(2,Ny+2))/2.0
        phi_x(Nx+2,Ny+2) = (phi_x(Nx+1,Ny+2) + phi_x(Nx+2,Ny+1))/2.0
        phi_y(1,1) = (phi_y(2,1) + phi_y(1,2))/2.0
        phi_y(Nx+2,1) = (phi_y(Nx+1,1) + phi_y(Nx+2,2))/2.0
        phi_y(1,Ny+2) = (phi_y(1,Ny+1) + phi_y(2,Ny+2))/2.0
        phi_y(Nx+2,Ny+2) = (phi_y(Nx+1,Ny+2) + phi_y(Nx+2,Ny+1))/2.0
        

ENDSUBROUTINE BCGRAD


SUBROUTINE REINIT(Nx,Ny,hx,hy,xwhole,ywhole,GA,scal,phi,phi_x,phi_y,phi_xy,prox,width)
    implicit none

    integer, intent(in) :: Nx,Ny,width
    integer, dimension(5,6), intent(in) :: GA
    real(kind=8), dimension(6), intent(in) :: scal
    real(kind=8), intent(in) :: hx,hy
    real(kind=8), dimension(Nx+2), intent(in) :: xwhole
    real(kind=8), dimension(Ny+2), intent(in) :: ywhole
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi        ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi_x      ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(inout) :: phi_y      ! gradient y-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi_xy        ! cross-derivative of phi
    integer, dimension(Nx+2,Ny+2), intent(in) :: prox

	! Local Variables
    integer :: max_iter,i,j,k,l,m,n,count,badNear
    real(kind=8) :: tol,alpha,gradmag2,step_x,step_y,step_dist
    real(kind=8) :: AA,BB,CC,DD,F1,F2
    real(kind=8), dimension(2) :: guess,prev,mid
    real(kind=8), dimension(6) :: herm
    real(kind=8), dimension(Nx+2,Ny+2) :: phin          ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2) :: phi_xn        ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2) :: phi_yn        ! gradient y-comp

    phin = phi
    phi_xn = phi_x
    phi_yn = phi_y

	prev(1) = hx
	prev(2) = hx

	! New level set 
    DO i = 2,Nx+1
        DO j = 2,Ny+1
        	if (prox(i,j) > 0) then
        		phi(i,j) = 1000.0
        	endif
        ENDDO
    ENDDO
        	
    max_iter = 100
    tol = 0.100E-15    
    alpha = 0.9
    badNear = 0
    
    DO i = 2,Nx+1
        DO j = 2,Ny+1
            
            ! For each (1) node, find the nearest interface point
            if (prox(i,j) == 1) then
                    
                herm(1) = phin(i,j)
                herm(2) = phi_xn(i,j)
                herm(3) = phi_yn(i,j)

                gradmag2 = phi_xn(i,j)**2 + phi_yn(i,j)**2
                step_x = -herm(1)*phi_xn(i,j)/gradmag2
                step_y = -herm(1)*phi_yn(i,j)/gradmag2
                
                guess(1) = xwhole(i)
                guess(2) = ywhole(j)
                
                !F1 = herm(1)
                !F2 = herm(2)*(ywhole(j)-guess(2)) - herm(3)*(xwhole(i)-guess(1))
                
                !CALL LINESEARCH(Nx,Ny,hx,hy,xwhole,ywhole,GA,scal,phin,phi_xn,phi_yn,phi_xy,i,j,guess,F1,F2,step_x,step_y,alpha)

                guess(1) = guess(1) + alpha*step_x
                guess(2) = guess(2) + alpha*step_y
                    
                ! Ensure that projected point is in the domain
                guess(1) = MIN(guess(1),xwhole(Nx+2))
                guess(2) = MIN(guess(2),ywhole(Ny+2))
                guess(1) = MAX(guess(1),xwhole(1))
                guess(2) = MAX(guess(2),ywhole(1))

                ! Find the index of the nearest node (for hermite interpolation).
                m = NINT((guess(1)+hx/2 - mod(guess(1)+hx/2.0,hx))/(hx)) + 1.0
                n = NINT((guess(2)+hy/2 - mod(guess(2)+hy/2.0,hy))/(hy)) + 1.0
                m = MAX(1,m);
                n = MAX(1,n);
                m = MIN(m,Nx+1);
                n = MIN(n,Ny+1);            
                
                ! Interpolate
            	CALL HERMCUBE(Nx,Ny,hx,hy,xwhole,ywhole,phin,phi_xn,phi_yn,phi_xy,herm,GA,scal,guess,m,n,6)
                
                F1 = herm(1)
                F2 = herm(2)*(ywhole(j)-guess(2)) - herm(3)*(xwhole(i)-guess(1))
                
                DO count = 1,max_iter

                    ! Calculate step size
                    AA = herm(2)
                    BB = herm(3)
                    CC = herm(4)*(ywhole(j)-guess(2)) + herm(3) - herm(6)*(xwhole(i)-guess(1))
                    DD = herm(6)*(ywhole(j)-guess(2)) - herm(2) - herm(5)*(xwhole(i)-guess(1))
                   
                    step_x = -1/(AA*DD-BB*CC)*( DD*F1-BB*F2)
                    step_y = -1/(AA*DD-BB*CC)*(-CC*F1+AA*F2)
                    
                    !CALL LINESEARCH(Nx,Ny,hx,hy,xwhole,ywhole,GA,scal,phin,phi_xn,phi_yn,phi_xy,i,j,guess,F1,F2,step_x,step_y,alpha)
                    
                    ! Take a step
                    guess(1) = guess(1) + alpha*step_x
                    guess(2) = guess(2) + alpha*step_y           
                    
                	! Ensure that projected point is in the domain
                    guess(1) = min(guess(1),xwhole(Nx+2))
                    guess(2) = min(guess(2),ywhole(Ny+2))
                    guess(1) = max(guess(1),xwhole(1))
                    guess(2) = max(guess(2),ywhole(1))

                    ! Find the index of the nearest node (for hermite interpolation).
                    m = NINT((guess(1)+hx/2 - mod(guess(1)+hx/2,hx))/(hx)) + 1
                    n = NINT((guess(2)+hy/2 - mod(guess(2)+hy/2,hy))/(hy)) + 1
                    m = MAX(1,m)
                    n = MAX(1,n)
                    m = MIN(m,Nx+1)
                    n = MIN(n,Ny+1)

                    ! Interpolate 
            		CALL HERMCUBE(Nx,Ny,hx,hy,xwhole,ywhole,phin,phi_xn,phi_yn, &
            					  phi_xy,herm,GA,scal,guess,m,n,6)
                    
                    F1 = herm(1)
                    F2 = herm(2)*(ywhole(j)-guess(2)) - herm(3)*(xwhole(i)-guess(1))
                    
                    if (ABS(F1) <= tol .AND. ABS(F2) <= tol) exit
                
                enddo

				if (abs(F1) <= tol) then
								
					! Compare CONVERGED interface point to all nearby points
					DO k = max(i-width,2),min(i+width,Nx+1)
						DO l = max(j-width,2),min(j+width,Ny+1)
					
							gradmag2 = SQRT((guess(1)-xwhole(k))**2 + (guess(2)-ywhole(l))**2)
							
							if (gradmag2 < abs(phi(k,l)) .AND. prox(k,l) > 0) then
								phi(k,l) = SIGN(1.0,real(phin(k,l)))*gradmag2
                				phi_x(k,l) = (xwhole(k)-guess(1))/phi(k,l)
                				phi_y(k,l) = (ywhole(l)-guess(2))/phi(k,l)
                			endif
                								
						ENDDO
					ENDDO
                    
                endif
            
            endif        
        enddo
    enddo 
    
    DO i = 2,Nx+1
        DO j = 2,Ny+1
            
            if (prox(i,j) > 1 .OR. prox(i,j) .LE. width) then
                    
                herm(1) = phi(i,j)
                herm(2) = phi_x(i,j)
                herm(3) = phi_y(i,j)

                gradmag2 = phi_x(i,j)**2 + phi_y(i,j)**2
                step_x = -herm(1)*phi_x(i,j)/gradmag2
                step_y = -herm(1)*phi_y(i,j)/gradmag2
                
                guess(1) = xwhole(i)
                guess(2) = ywhole(j)
                
                !F1 = herm(1)
                !F2 = herm(2)*(ywhole(j)-guess(2)) - herm(3)*(xwhole(i)-guess(1))
                
                !CALL LINESEARCH(Nx,Ny,hx,hy,xwhole,ywhole,GA,scal,phin,phi_xn,phi_yn,phi_xy,i,j,guess,F1,F2,step_x,step_y,alpha)

                guess(1) = guess(1) + alpha*step_x
                guess(2) = guess(2) + alpha*step_y
                    
                ! Ensure that projected point is in the domain
                guess(1) = MIN(guess(1),xwhole(Nx+2))
                guess(2) = MIN(guess(2),ywhole(Ny+2))
                guess(1) = MAX(guess(1),xwhole(1))
                guess(2) = MAX(guess(2),ywhole(1))

                ! Find the index of the nearest node (for hermite interpolation).
                m = NINT((guess(1)+hx/2 - mod(guess(1)+hx/2.0,hx))/(hx)) + 1.0
                n = NINT((guess(2)+hy/2 - mod(guess(2)+hy/2.0,hy))/(hy)) + 1.0
                m = MAX(1,m);
                n = MAX(1,n);
                m = MIN(m,Nx+1);
                n = MIN(n,Ny+1);            
                
                ! Interpolate
            	CALL HERMCUBE(Nx,Ny,hx,hy,xwhole,ywhole,phin,phi_xn,phi_yn,phi_xy,herm,GA,scal,guess,m,n,6)
                
                F1 = herm(1)
                F2 = herm(2)*(ywhole(j)-guess(2)) - herm(3)*(xwhole(i)-guess(1))
                
                DO count = 1,max_iter

                    ! Calculate step size
                    AA = herm(2)
                    BB = herm(3)
                    CC = herm(4)*(ywhole(j)-guess(2)) + herm(3) - herm(6)*(xwhole(i)-guess(1))
                    DD = herm(6)*(ywhole(j)-guess(2)) - herm(2) - herm(5)*(xwhole(i)-guess(1))
                   
                    step_x = -1/(AA*DD-BB*CC)*( DD*F1-BB*F2)
                    step_y = -1/(AA*DD-BB*CC)*(-CC*F1+AA*F2)
                    
                    !CALL LINESEARCH(Nx,Ny,hx,hy,xwhole,ywhole,GA,scal,phin,phi_xn,phi_yn,phi_xy,i,j,guess,F1,F2,step_x,step_y,alpha)
                    
                    ! Take a step
                    guess(1) = guess(1) + alpha*step_x
                    guess(2) = guess(2) + alpha*step_y           
                    
                	! Ensure that projected point is in the domain
                    guess(1) = min(guess(1),xwhole(Nx+2))
                    guess(2) = min(guess(2),ywhole(Ny+2))
                    guess(1) = max(guess(1),xwhole(1))
                    guess(2) = max(guess(2),ywhole(1))

                    ! Find the index of the nearest node (for hermite interpolation).
                    m = NINT((guess(1)+hx/2 - mod(guess(1)+hx/2,hx))/(hx)) + 1
                    n = NINT((guess(2)+hy/2 - mod(guess(2)+hy/2,hy))/(hy)) + 1
                    m = MAX(1,m)
                    n = MAX(1,n)
                    m = MIN(m,Nx+1)
                    n = MIN(n,Ny+1)

                    ! Interpolate 
            		CALL HERMCUBE(Nx,Ny,hx,hy,xwhole,ywhole,phin,phi_xn,phi_yn, &
            					  phi_xy,herm,GA,scal,guess,m,n,6)
                    
                    F1 = herm(1)
                    F2 = herm(2)*(ywhole(j)-guess(2)) - herm(3)*(xwhole(i)-guess(1))
                    
                    if (ABS(F1) <= tol .AND. ABS(F2) <= tol) exit
                
                enddo

				if (abs(F1) <= tol) then
								
					! Compare CONVERGED interface point to all nearby points
					DO k = max(i-width,2),min(i+width,Nx+1)
						DO l = max(j-width,2),min(j+width,Ny+1)
					
							gradmag2 = SQRT((guess(1)-xwhole(k))**2 + (guess(2)-ywhole(l))**2)
							
							if (gradmag2 < abs(phi(k,l)) .AND. prox(k,l) > 0) then
								phi(k,l) = SIGN(1.0,real(phin(k,l)))*gradmag2
                				phi_x(k,l) = (xwhole(k)-guess(1))/phi(k,l)
                				phi_y(k,l) = (ywhole(l)-guess(2))/phi(k,l)
                			endif
                								
						ENDDO
					ENDDO
                    
                endif    
           endif        
        enddo
    enddo    
    
    ! Adjust level set values near small features
    DO i = 2,Nx+1
        DO j = 2,Ny+1
    		if (phi(i,j)*phi(i,j+1) < 0.0 .AND. phi(i,j)*phi(i,j-1) < 0.0) then
    			phi_y(i,j) = (phi(i,j+1) - phi(i,j-1))/(2.0*hy)
    		endif
    		if (phi(i,j)*phi(i+1,j) < 0.0 .AND. phi(i,j)*phi(i-1,j) < 0.0) then
    			phi_x(i,j) = (phi(i+1,j) - phi(i-1,j))/(2.0*hx)
    		endif
    		if (phi(i,j)*phi(i+1,j) < 0.0 .AND. phi(i,j)*phi(i-1,j) < 0.0 .OR. &
    		    phi(i,j)*phi(i,j+1) < 0.0 .AND. phi(i,j)*phi(i,j-1) < 0.0) then
    			phi(i,j) = 0.5*phi(i,j)
    		endif
    	ENDDO
    ENDDO
    
END SUBROUTINE REINIT

SUBROUTINE LINESEARCH(Nx,Ny,hx,hy,xwhole,ywhole,GA,scal,phi,phi_x,phi_y,phi_xy,i,j,guess,F1,F2,step_x,step_y,alpha)
    implicit none

    integer, intent(in) :: Nx,Ny,i,j
    integer, dimension(5,6), intent(in) :: GA
    real(kind=8), dimension(6), intent(in) :: scal
    real(kind=8), intent(in) :: hx,hy,F1,F2
    real(kind=8), dimension(Nx+2), intent(in) :: xwhole
    real(kind=8), dimension(Ny+2), intent(in) :: ywhole
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi        ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi_x      ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi_y      ! gradient y-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi_xy        ! cross-derivative of phi
    real(kind=8), dimension(2),intent(in) :: guess
    real(kind=8), intent(in) :: step_x,step_y
    real(kind=8), intent(out) :: alpha

	integer :: k
    real(kind=8), dimension(6) :: herm2
    real(kind=8), dimension(2) :: X
	real(kind=8), dimension(5) :: trial
	real(kind=8) :: F1b,F2b
	 
	trial(1) = 1.0
	trial(2) = 0.5
	trial(3) = 0.25
	trial(4) = 0.1
	trial(5) = 0.01

	DO k=1,5
	
		alpha = trial(k)
	
		X(1) = guess(1)+alpha*step_x
		X(2) = guess(2)+alpha*step_y
	
    	CALL HERMCUBE(Nx,Ny,hx,hy,xwhole,ywhole,phi,phi_x,phi_y,phi_xy,herm2, &
                      GA,scal,X,i,j,3)
                      
        F1b = herm2(1)
        F2b = herm2(2)*(ywhole(j)-X(2)) - herm2(3)*(xwhole(i)-X(1))     
                      
        if ( ABS(F1b)+ABS(F2b) < ABS(F1)+ABS(F2) ) exit
	
	ENDDO
	
	!print *, 'alpha:',  alpha

END SUBROUTINE LINESEARCH

REAL(kind=8) FUNCTION GALSCURV(Nx,Ny,hx,hy,xwhole,ywhole,phi,phi_x,phi_y,phi_xy,prox,i,j,&
                               GA,scal,mode)

    IMPLICIT NONE
    
    integer, intent(in) :: Nx,Ny,i,j,mode
    real(kind=8), intent(in) :: hx,hy
    integer, dimension(5,6), intent(in) :: GA
    real(kind=8), dimension(6), intent(in) :: scal
    real(kind=8), dimension(Nx+2), intent(in) :: xwhole
    real(kind=8), dimension(Ny+2), intent(in) :: ywhole
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi        ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi_x      ! gradient x-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi_y      ! gradient y-comp
    real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: phi_xy        ! cross-derivative of phi
    integer, dimension(Nx+2,Ny+2), intent(in) :: prox

	REAL(kind=8), dimension(2) :: XX
    REAL(kind=8), dimension(6) :: herm
    REAL(kind=8) :: norm
    REAL(kind=8) :: smear,eps

    !if (prox(i,j) > 0 .AND. prox(i,j) < 3) then

    if (mode == 1) then 
        XX(1) = xwhole(i) + hx/2.0
        XX(2) = ywhole(j)
    elseif (mode == 2) then
        XX(1) = xwhole(i)
        XX(2) = ywhole(j) + hy/2.0
    endif            

    CALL HERMCUBE(Nx,Ny,hx,hy,xwhole,ywhole,phi,phi_x,phi_y,phi_xy,herm, &
                  GA,scal,XX,i,j,6)
        
    ! magnitude of the the gradient of phi
    norm = (herm(2)**2.0 + herm(3)**2.0)**0.5

    ! curvature
    if (norm > 1e-16) then
        GALSCURV = (herm(4)*herm(3)**2.0 - 2.0*herm(3)*herm(2)*herm(6) + herm(5)*herm(2)**2.0)/norm**3.0
    else
        GALSCURV = 0.0
    endif
            
    ! bounding the curvature to the grid scale, (-1/hx<=kappa<=1/hx)
    !IF ( abs(GALSCURV) >= 1.0/min(hx,hy) ) THEN
    !        IF (GALSCURV < 0.0) THEN
    !                GALSCURV = -1.0/min(hx,hy)
    !        else
    !                GALSCURV = 1.0/min(hx,hy)
    !        ENDIF
    !ENDIF

    ! Curvature times the normal vector
    !IF (mode == 1) then
    !    GALSCURV = GALSCURV*herm(2)/norm
    !ELSEIF (mode == 2) then
    !    GALSCURV = GALSCURV*herm(3)/norm
    !ENDIF

    ! Smeared delta function (in the manner of Engquist 2005)
    !if (norm > 1e-16) then    
    !    eps = 1.5*max(hx,hy)*(abs(herm(2))+abs(herm(3)))/norm
    !else
    !    eps = 1.5*max(hx,hy)
    !endif
    !if (abs(herm(1)/eps) > 1.0) then
    !    smear = 0.0
    !else
    !    smear = (1.0/eps)*(1.0-abs(herm(1)/eps))
    !endif
    !GALSCURV = GALSCURV*smear

    !else
    !    GALSCURV = 0.0
    !endif

END FUNCTION GALSCURV
