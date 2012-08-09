!********************************************************************
!*																	*
!*							    PCG									*
!*																	*
!********************************************************************
!* Author: Nathaniel R. Morgan	(June, 2004)						*
!*																	*
!* This is a program for solving the pressure poisson equation		*
!* (PPE) with variable density on a staggered grid using the		*
!* preconditioned conjugate gradient method (PCG).					*
!*																	*
!* INPUT VARIABLES:													*
!* Nx			The number of interior nodes in the x-direction		*
!* Ny			The number of interior nodes in the y-direction		*
!* deltaX		The step size in the x-direction					*
!* deltaY		The step size in the y-direction					*
!* Xarray		The initial guess of the pressure field				*
!* Barray		The RHS of the PPE, (note: AX=B)					*
!* H			The heavieside function								*
!* rho_one		The density of the first field	(H=1)				*
!* rho_two		The density of the second field (H=0)				*
!* tolerance	The tolerance										*
!* kmax			The max number of iterations						*
!* iterprint	The number of iterations between printing residual	*
!* precond		The type of preconditioning, 0=none & 1=diagonal	*
!* BCPnorth		The north BC,1=Neumman, 2=Dirclet					*
!* BCPsouth		The south BC,	"" ""								*
!* BCPeast		The east BC,	"" ""								*
!* BCPwest		The west BC,	"" ""								*
!* Pnorth		The north pressure, corresponds to the Dirclet BC,2	*
!* Psouth		The south pressure,		"" ""						*
!* Peast		The east pressure,		"" ""						*
!* Pwest		The west pressure,		"" ""						*
!*																	*
!* OUTPUT VARIABLES:												*
!* Xarray		The new pressure field								*
!*																	*
!* Warning: No guarantees are implied concerning any algorythm		*
!********************************************************************
SUBROUTINE PCG(Nx,Ny,deltaX,deltaY,Xarray,Barray,H,s_H,&
			   rho_one,rho_two,rho_three,&
			   tolerance,kmax,iterprint,precond,&
			   BCPnorth,BCPsouth,BCPeast,BCPwest,Pnorth,Psouth,Peast,Pwest)
	
	implicit none
	
	! local integers
	integer :: k,iterprint

	! integers
	integer, intent(in) :: Nx,Ny,precond,kmax
	integer, intent(in) :: BCPnorth,BCPsouth,BCPeast,BCPwest
	
	! local reals
	real(kind=8) :: dummy,err,terminate,alpha,beta,gamma,gammaold

	! reals
	real(kind=8), intent(in) :: deltaX,deltaY
	real(kind=8), intent(in) :: tolerance
	real(kind=8), intent(in) :: Pnorth,Psouth,Peast,Pwest
	real(kind=8), intent(in) :: rho_one,rho_two,rho_three

	! local real arrays
	real(kind=8), dimension(Nx*Ny) :: X,B,AV,p,r,z
	
	! real arrays
	real(kind=8), dimension(Nx,Ny), intent(in) :: Barray
	real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: H
	real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: s_H
	real(kind=8), dimension(Nx,Ny), intent(inout) :: Xarray

	! external functions
	REAL(kind=8), EXTERNAL :: NORM,VECDOTVEC

	
	! convert arrays to vectors
	call ArrtoVec(Nx,Ny,Barray,B)
	call ArrtoVec(Nx,Ny,Xarray,X)

	! I need to switch the sign on B because I multiply the A matrix by -1
	! to make it positive definite
	B=-B

	! apply BC's to B vector
	call BCB(Nx,Ny,deltaX,deltaY,BCPnorth,BCPsouth,BCPeast,BCPwest,&
			 Pnorth,Psouth,Peast,Pwest,B,H,s_H,rho_one,rho_two,rho_three)

	! get first residual, r=B-AX
	call AtimesV(Nx,Ny,deltaX,deltaY,BCPnorth,BCPsouth,BCPeast,BCPwest,&
				 X,AV,H,s_H,rho_one,rho_two,rho_three)
	r = B-AV
	

	! start the minimization process
	k=0
	dummy = norm(Nx,Ny,r)
	!print *, "norm of r"
	!print *, dummy
	err = dummy**0.5
	dummy = norm(Nx,Ny,B)
	!print *, "norm of RHS"
	!print *, dummy

	terminate = tolerance*sqrt(dummy)
	
	do
		if (err<=terminate .or. k>kmax) exit
		k=k+1
    
		! preconditioning, solving Mz=r
		if (precond == 0) then
			z=r										! no-preconditioning
		else
			! solving Pz=r
			call PCv2(Nx,Ny,deltaX,deltaY,r,z,H,s_H,&
					 rho_one,rho_two,rho_three)		! diagonal preconditioning
		endif
		
		gamma = VECDOTVEC(Nx,Ny,r,z)				! gamma is a substitutional variable
		! get new p vector
		if (k==1) then
			p=z  ! note: MP=r
		else
			! beta=r'*z/(rold'*zold)
			beta=gamma/gammaold
			p=z+beta*p
		endif
		    	
		call AtimesV(Nx,Ny,deltaX,deltaY,BCPnorth,BCPsouth,BCPeast,BCPwest,&
					 p,AV,H,s_H,rho_one,rho_two,rho_three)
		! alpha = (r'*z)/(p'*Av)
		alpha = VECDOTVEC(Nx,Ny,p,AV)     ! alpha is used here as a substitutional variable
		alpha = gamma/alpha
		
		X=X+alpha*p         ! get new solution
		r=r-alpha*AV        ! get new residual
		
		dummy = norm(Nx,Ny,r)
		err=dummy**0.5	! get new error
	
		! store old value
		gammaold=gamma

		if (mod(k,iterprint)==0) then
			print *, 'Pressure Residual =', err
		endif
	enddo

	!print *, 'Pressure Residual =', err
	if (k>kmax) then
		print *, 'ERROR: Pressure failed to converge'
	endif
	
	print *, "CG Iterations = ",k

	! convert back to array
	call VectoArr(Nx,Ny,Xarray,X)

ENDSUBROUTINE PCG
!********************************************************************
!*																	*
!*							VECDOTVEC								*
!*																	*
!********************************************************************
!* Author: Nathaniel R. Morgan										*
!*																	*
!* This function takes the 2-norm of a vector						*
!*																	*
!********************************************************************
REAL(kind=8) FUNCTION NORM(Nx,Ny,vec)
	implicit none
	
	integer :: m
	integer, intent(in) :: Nx,Ny
	real(kind=8), dimension(Nx*Ny), intent(in) :: vec
	
	NORM = 0.0
	do m=1,Nx*Ny
		NORM = NORM + vec(m)**2.0
	enddo

	NORM = sqrt(NORM) 

ENDFUNCTION
!********************************************************************
!*																	*
!*							VECDOTVEC								*
!*																	*
!********************************************************************
!* Author: Nathaniel R. Morgan										*
!*																	*
!* This function takes the dot product of two vectors				*
!*																	*
!********************************************************************
REAL(kind=8) FUNCTION VECDOTVEC(Nx,Ny,vecone,vectwo)
	implicit none
	
	integer :: m
	integer, intent(in) :: Nx,Ny
	real(kind=8), dimension(Nx*Ny), intent(in) :: vecone,vectwo
	
	VECDOTVEC = 0.0
	do m=1,Nx*Ny
		VECDOTVEC = VECDOTVEC + vecone(m)*vectwo(m)
	enddo
ENDFUNCTION VECDOTVEC
!********************************************************************
!*																	*
!*							ArrtoVec								*
!*																	*
!********************************************************************
!* Author: Nathaniel R. Morgan										*
!*																	*
!* This is a program for converting an array into a vector			*													
!*																	*
!********************************************************************
subroutine ArrtoVec(Nx,Ny,Array,vec)

	implicit none
	
	integer :: mn,ms,me,mw,mo,i,j

	integer, intent(in) :: Nx,Ny
	real(kind=8), dimension(Nx,Ny), intent(in) :: Array
	real(kind=8), dimension(Nx*Ny), intent(out) :: vec

	do j=1,Ny
		do i=1,Nx
			! m is the matrix index, in 3D m=(j-1)Nx+(k-1)*Nx*Ny+i
			mn = j*Nx + i          ! i,j+1
			ms = (j-2)*Nx + i      ! i,j-1
			me = (j-1)*Nx + i+1    ! i+1,j
			mw = (j-1)*Nx + i-1    ! i-1,j
			mo = (j-1)*Nx + i      ! i,j
        
			if (j/=Ny) then
				vec(mn) = Array(i,j+1)
			endif
			if (j/=1) then
				vec(ms) = Array(i,j-1)
			endif
			if (i/=Nx) then
				vec(me) = Array(i+1,j)
			endif
			if (i/=1) then
				vec(mw) = Array(i-1,j)
			endif
			vec(mo) = Array(i,j)
		enddo
	enddo

endsubroutine ArrtoVec
!********************************************************************
!*																	*
!*							VectoArray								*
!*																	*
!********************************************************************
!* Author: Nathaniel R. Morgan										*
!*																	*
!* This is a program for converting a vector into an array			*													
!*																	*
!********************************************************************
subroutine VectoArr(Nx,Ny,Array,Vec)
	implicit none

	integer :: i,j,m
	integer, intent(in) :: Nx,Ny
	real(kind=8), dimension(Nx,Ny), intent(out) :: Array
	real(kind=8), dimension(Nx*Ny), intent(in) :: vec

	do m=1,Nx*Ny
		! note: in 3D k=floor((m-1)/(Nx*Ny))+1
		j = floor( real(m-1)/real(Nx) )+1  
		i = m-Nx*(j-1)
    
		! storing the answer in RHS array
		Array(i,j) = Vec(m)
	enddo

ENDSUBROUTINE VectoArr
!********************************************************************
!*																	*
!*							AtimesV									*
!*																	*
!********************************************************************
!* Author: Nathaniel R. Morgan										*
!*																	*
!* This is a program for multipling the A matrix by the Krylov		*
!* vector, V.  The A matrix is NEVER stored because I only need the	*													
!* product of A times V, AV.										*
!*																	*
!********************************************************************
SUBROUTINE AtimesV(Nx,Ny,deltaX,deltaY,BCPnorth,BCPsouth,BCPeast,BCPwest, &
				   vec,AV,H,s_H,rho_one,rho_two,rho_three)
	
	implicit none

	integer :: mn,ms,me,mw,mo,i,j

	integer, intent(in) :: Nx,Ny,BCPnorth,BCPsouth,BCPeast,BCPwest

	real(kind=8) :: an,as,ae,aw,ao,Axn,Axs,Axe,Axw,Axo,rho_e,rho_w,rho_n,rho_s

	real(kind=8), intent(in) :: deltaX,deltaY,rho_one,rho_two,rho_three
	
	real(kind=8), dimension(Nx*Ny), intent(in) :: vec
	real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: H
	real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: s_H
	real(kind=8), dimension(Nx*Ny), intent(out) :: AV


	! remember X has the m index
	do j=1,Ny
		do i=1,Nx
			! m is the matrix index, in 3D m=(j-1)Nx+(k-1)*Nx*Ny+i
			mn = j*Nx + i          ! i,j+1
			ms = (j-2)*Nx + i      ! i,j-1
			me = (j-1)*Nx + i+1    ! i+1,j
			mw = (j-1)*Nx + i-1    ! i-1,j
			mo = (j-1)*Nx + i      ! i,j

			! density
			rho_e = 0.5*(rho_one*H(i+2,j+1) + rho_two*(1.0-H(i+2,j+1))) + &
					0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1)))
			rho_w = 0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1))) + &
					0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1)))
			rho_n = 0.5*(rho_one*H(i+1,j+2) + rho_two*(1.0-H(i+1,j+2))) + &
					0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1)))
			rho_s = 0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1))) + &
					0.5*(rho_one*H(i+1,j) + rho_two*(1.0-H(i+1,j)))
        
			! matrix coefficients (note: multiplied by -1)
			an = -1.0/(rho_n*deltaY**2.0)
			as = -1.0/(rho_s*deltaY**2.0)
			ae = -1.0/(rho_e*deltaX**2.0)
			aw = -1.0/(rho_w*deltaX**2.0)
			ao = -an-as-ae-aw
        
        
			! north   
			if (j/=Ny) then
				!A(mo,mn) = an 
				Axn = an*vec(mn)
			else
				Axn = 0.0
			endif
        
			! south
			if (j/=1) then
				!A(mo,ms) = as
				Axs = as*vec(ms)
			else
				Axs = 0.0
			endif
         
			! east
			if (i/=Nx) then
				!A(mo,me) = ae
				Axe = ae*vec(me)
			else
				Axe = 0.0
			endif
        
			! west
			if (i/=1) then
				!A(mo,mw) = aw
				Axw = aw*vec(mw)
			else
				Axw = 0.0
			endif
			

			!--- NEUMMAN BC ---	
			if (j==Ny) then
				! neumman BC
				if (BCPnorth==1) then
					ao = ao+an 
				endif
			endif
			if (j==1) then
				! neumman BC
				if (BCPsouth==1) then
					if (BCPnorth==1 .and. BCPwest==1 .and. BCPeast==1) then
						! if neumman on all boundaries then
						! I need a dirclet BC.  I am making
						! pressure zero at i=Nx/2,j=1
						if (i/=floor(real(Nx/2))) then
							ao = ao+as 
						endif
					else
						ao = ao+as
					endif
				endif
			endif
			if (i==Nx) then
				! neumman BC
				if (BCPeast==1) then
					ao = ao+ae 
				endif
			endif
			if (i==1) then
				! neumman BC
				if (BCPwest==1) then
					ao = ao+aw 
	            endif
			endif


			! origin
			!A(mo,mo) = ao
			Axo = ao*vec(mo)
        
			Av(mo) = Axo + Axn + Axs + Axe + Axw
		enddo
	enddo

ENDSUBROUTINE AtimesV
!********************************************************************
!*																	*
!*							 PCv									*
!*																	*
!********************************************************************
!* Author: Nathaniel R. Morgan										*
!*																	*
!* This is a program for applying a diagonal preconditioner	to the	*
!* A matrix. The diagonal preconditioner works by dividing the A	*
!* matrix by the diagonal elements (similar to the Jacobi method).	*
!*																	*
!********************************************************************
SUBROUTINE PCv(Nx,Ny,deltaX,deltaY,r,vec,H,s_H,rho_one,rho_two,rho_three)

	implicit none

	integer :: mo,i,j

	integer, intent(in) :: Nx,Ny

	real(kind=8) :: an,as,ae,aw,ao,rho_e,rho_w,rho_n,rho_s

	real(kind=8), intent(in) :: deltaX,deltaY,rho_one,rho_two,rho_three
	real(kind=8), dimension(Nx*Ny), intent(out) :: vec
	real(kind=8), dimension(Nx*Ny), intent(in) :: r
	real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: H
	real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: s_H

	! remember r has the m index
	do j=1,Ny
		do i=1,Nx
			! m is the matrix index, in 3D m=(j-1)Nx+(k-1)*Nx*Ny+i
			mo = (j-1)*Nx + i		! i,j
        
			! density
			rho_e = 0.5*(rho_one*H(i+2,j+1) + rho_two*(1.0-H(i+2,j+1))) + &
					0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1)))
			rho_w = 0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1))) + &
					0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1)))
			rho_n = 0.5*(rho_one*H(i+1,j+2) + rho_two*(1.0-H(i+1,j+2))) + &
					0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1)))
			rho_s = 0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1))) + &
					0.5*(rho_one*H(i+1,j) + rho_two*(1.0-H(i+1,j)))


			! matrix coefficients
			an = -1.0/(rho_n*deltaY**2.0)
			as = -1.0/(rho_s*deltaY**2.0)
			ae = -1.0/(rho_e*deltaX**2.0)
			aw = -1.0/(rho_w*deltaX**2.0)
			ao = -an-as-ae-aw		! note: diagonal coefficient
        
			vec(mo) = 1.0/ao*r(mo)
		enddo
	enddo

ENDSUBROUTINE PCv
!********************************************************************
!*																	*
!*							 PCv2									*
!*																	*
!********************************************************************
!* Author: Curtis Lee												*
!*																	*
!* SSOR Preconditioner for much faster convergence than PCv			*
!********************************************************************
SUBROUTINE PCv2(Nx,Ny,deltaX,deltaY,r,vec,H,s_H,rho_one,rho_two,rho_three)

	implicit none

	integer :: mn,ms,me,mw,mo,i,j,a,b

	integer, intent(in) :: Nx,Ny

	real(kind=8) :: an,as,ae,aw,ao,rho_e,rho_w,rho_n,rho_s

	real(kind=8) :: omega

	real(kind=8), intent(in) :: deltaX,deltaY,rho_one,rho_two,rho_three
	real(kind=8), dimension(Nx*Ny), intent(out) :: vec
	real(kind=8), dimension(Nx*Ny), intent(in) :: r
	real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: H
	real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: s_H
	
	! Convergence speed is HIGHLY sensitive to omega value
	omega = 1.75

	! FORWARD SWEEP (i -->, j ^^^)
	do j=1,Ny
		do i=1,Nx
			! m is the matrix index, in 3D m=(j-1)Nx+(k-1)*Nx*Ny+i

			ms = (j-2)*Nx + i     ! i,j-1
			mw = (j-1)*Nx + i-1    ! i-1,j
			mo = (j-1)*Nx + i      ! i,j        
			
			! density
			rho_e = 0.5*(rho_one*H(i+2,j+1) + rho_two*(1.0-H(i+2,j+1))) + &
					0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1)))
			rho_w = 0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1))) + &
					0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1)))
			rho_n = 0.5*(rho_one*H(i+1,j+2) + rho_two*(1.0-H(i+1,j+2))) + &
					0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1)))
			rho_s = 0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1))) + &
					0.5*(rho_one*H(i+1,j) + rho_two*(1.0-H(i+1,j)))

			! matrix coefficients
			an = 1.0/(rho_n*deltaY**2.0)
			as = 1.0/(rho_s*deltaY**2.0)
			ae = 1.0/(rho_e*deltaX**2.0)
			aw = 1.0/(rho_w*deltaX**2.0)
			ao = an+as+ae+aw		! note: diagonal coefficient
			
        	if (i > 1 .AND. j > 1) then
				vec(mo) = omega/ao*(r(mo)+aw*vec(mw)+as*vec(ms))
        	elseif (j > 1) then
				vec(mo) = omega/ao*(r(mo)+as*vec(ms))
        	elseif (i > 1) then
				vec(mo) = omega/ao*(r(mo)+aw*vec(mw))
        	else
				vec(mo) = omega/ao*(r(mo))
			endif
			
		enddo
	enddo
	
	! STATIC SWEEP
	do j=1,Ny
		do i=1,Nx
			
			mo = (j-1)*Nx + i      ! i,j    
			
			! density
			rho_e = 0.5*(rho_one*H(i+2,j+1) + rho_two*(1.0-H(i+2,j+1))) + &
					0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1)))
			rho_w = 0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1))) + &
					0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1)))
			rho_n = 0.5*(rho_one*H(i+1,j+2) + rho_two*(1.0-H(i+1,j+2))) + &
					0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1)))
			rho_s = 0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1))) + &
					0.5*(rho_one*H(i+1,j) + rho_two*(1.0-H(i+1,j)))

			! matrix coefficients
			an = 1.0/(rho_n*deltaY**2.0)
			as = 1.0/(rho_s*deltaY**2.0)
			ae = 1.0/(rho_e*deltaX**2.0)
			aw = 1.0/(rho_w*deltaX**2.0)
			ao = an+as+ae+aw		! note: diagonal coefficient
			
			vec(mo) = ((2.0-omega)/omega)*ao*vec(mo)
			
		enddo
	enddo
	
	! BACKWARDS SWEEP (i <--, j vvv)
	do j=Ny,1,-1
		do i=Nx,1,-1
		
			! m is the matrix index, in 3D m=(j-1)Nx+(k-1)*Nx*Ny+i
			mn = j*Nx + i          ! i,j+1
			me = (j-1)*Nx + i+1    ! i+1,j
			mo = (j-1)*Nx + i      ! i,j        
			
			! density
			rho_e = 0.5*(rho_one*H(i+2,j+1) + rho_two*(1.0-H(i+2,j+1))) + &
					0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1)))
			rho_w = 0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1))) + &
					0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1)))
			rho_n = 0.5*(rho_one*H(i+1,j+2) + rho_two*(1.0-H(i+1,j+2))) + &
					0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1)))
			rho_s = 0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1))) + &
					0.5*(rho_one*H(i+1,j) + rho_two*(1.0-H(i+1,j)))

			! matrix coefficients
			an = 1.0/(rho_n*deltaY**2.0)
			as = 1.0/(rho_s*deltaY**2.0)
			ae = 1.0/(rho_e*deltaX**2.0)
			aw = 1.0/(rho_w*deltaX**2.0)
			ao = an+as+ae+aw		! note: diagonal coefficient
        
        	if (i < Nx .AND. j < Ny) then
				vec(mo) = omega/ao*(vec(mo)+ae*vec(me)+an*vec(mn))
        	elseif (j < Ny) then
				vec(mo) = omega/ao*(vec(mo)+an*vec(mn))
        	elseif (i < Nx) then
				vec(mo) = omega/ao*(vec(mo)+ae*vec(me))
        	else
				vec(mo) = omega/ao*(vec(mo))
			endif
			
		enddo
	enddo
	

ENDSUBROUTINE PCv2
!********************************************************************
!*																	*
!*							 BCB									*
!*																	*
!********************************************************************
!* Author: Nathaniel R. Morgan										*
!*																	*
!* This is a program for applying the boundary conditions to the	*
!* RHS vector, B.  This program is only used for dirclet boundary	*
!* conditions.  Note, when all the boundary conditions are neumman,	*
!* I make one point on the bottom boundary a dirclet condition, ie. *
!* I set the value to zero.											*
!*																	*
!********************************************************************
SUBROUTINE BCB(Nx,Ny,deltaX,deltaY,BCPnorth,BCPsouth,BCPeast,BCPwest,&
			   Pnorth,Psouth,Peast,Pwest,B,H,s_H,rho_one,rho_two,rho_three)

	implicit none

	integer :: mn,ms,me,mw,mo,i,j

	integer, intent(in) :: Nx,Ny,BCPnorth,BCPsouth,BCPeast,BCPwest

	real(kind=8) :: an,as,ae,aw,ao,rho_e,rho_w,rho_n,rho_s

	real(kind=8), intent(in) :: deltaX,deltaY,Pnorth,Psouth,Peast,Pwest
	real(kind=8), intent(in) :: rho_one,rho_two,rho_three
	
	real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: H
	real(kind=8), dimension(Nx+2,Ny+2), intent(in) :: s_H
	real(kind=8), dimension(Nx*Ny), intent(inout) :: B



	! remember b has the m index
	do j=1,Ny
		do i=1,Nx
			! m is the matrix index, in 3D m=(j-1)Nx+(k-1)*Nx*Ny+i
			mn = j*Nx + i          ! i,j+1
			ms = (j-2)*Nx + i      ! i,j-1
			me = (j-1)*Nx + i+1    ! i+1,j
			mw = (j-1)*Nx + i-1    ! i-1,j
			mo = (j-1)*Nx + i      ! i,j

			! density
			rho_e = 0.5*(rho_one*H(i+2,j+1) + rho_two*(1.0-H(i+2,j+1))) + &
					0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1)))
			rho_w = 0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1))) + &
					0.5*(rho_one*H(i,j+1) + rho_two*(1.0-H(i,j+1)))
			rho_n = 0.5*(rho_one*H(i+1,j+2) + rho_two*(1.0-H(i+1,j+2))) + &
					0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1)))
			rho_s = 0.5*(rho_one*H(i+1,j+1) + rho_two*(1.0-H(i+1,j+1))) + &
					0.5*(rho_one*H(i+1,j) + rho_two*(1.0-H(i+1,j)))

        
			! matrix coefficients
			an = -1.0/(rho_n*deltaY**2.0)
			as = -1.0/(rho_s*deltaY**2.0)
			ae = -1.0/(rho_e*deltaX**2.0)
			aw = -1.0/(rho_w*deltaX**2.0)    
		
			
			! --- DIRCLET BC'S ---
			! north   
			if (j==Ny) then
				! dirclet BC
				if (BCPnorth==2) then
					B(mo) = B(mo)-an*Pnorth
				endif
			endif
			! south    
			if (j==1) then
				! dirclet BC
				if (BCPsouth==2) then
					B(mo) = B(mo)-as*Psouth
				endif
			endif
			! east
			if (i==Nx) then
				! dirclet BC
				if (BCPeast==2) then
					B(mo) = B(mo)-ae*Peast
				endif
			endif
			! west
			if (i==1) then
				! dirclet BC
				if (BCPwest==2) then
					B(mo) = B(mo)-aw*Pwest
				endif
			endif


			! --- NEUMMAN BC'S ---
			! note, I must have 1 dirclet condition at a point if all boundaries are neumman
			if (j==1) then
				! making the center value on the bottom edge zero 
				! if all the boundaries are neumman
				if (BCPsouth==1 .and. BCPnorth==1 .and. BCPeast==1 .and. BCPwest==1) then
					if (i==floor(real(Nx/2))+1) then
						B(mo) = B(mo)-as*0.0 
					endif
				endif
			endif

		enddo
	enddo

ENDSUBROUTINE BCB