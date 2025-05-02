program spectral
! usage: 
!        ./spectral  <rad> <omega> <Ndpt>  <Nmom> <idir>
!
!
! Find the polarization tensor alpha[1,1] or alpha[2,2] as a function of the variable s=(eps1-eps2)/eps2
!
!   alpha[k,k] (s) = s*[nu[0]-nu[1]*s+ nu[2]*s^2 - .....+ nu[n]*(-s)^n + ...]
!
! Method:
!
! set rho[0] = x1 if IDIR=1, rho[0] = x2 if IDIR=2
!
! iterate: rho[k] = int_S K(x,y) rho[k-1](y)  ds(y),     k.GE.1
!
! with K(x,y)= n(y) . grad_y G(x-y) ,   G(x-y)= -1/(2*Pi)* log|x-y|
!      n(y)= (n1(y), n2(y)) unit outer normal
!
! The inclusion has permittivity eps1 in a matrix of permittivity eps2
! of shape S given in polar coordinates:
!  
!                           r= 1.d0 + rad*dcos(omega*theta),  0 < theta < 2*pi
!  
! Inputs on command line:
!            rad: radial geometric parameter
!            omega: number of waves of star
!            Ndpt: Number of boundary points /unit length (integer)
!            Nmom: Number of moments (integer)
!            IDIR: 1 or 2
!
! Output: coefficients of series expansion
!
!          nu[0], ..., nu[Nmom],    nu[k] = int_S n[IDIR]* rho[k]*ds
!  
!
! Note #1: the series must be summed by Pade approximants to accelerate the series.
! use: maple < spect.mw  (the input data file is "numoments.dat")
! Note #2: Shape of S is found is shape.dat
!
! Author: R.V. Roy, 07/26/2023
!
!....................................................................
!  
implicit double precision (a-h,o-z)
double precision, parameter::pi=3.141592653589793238d0
double precision, allocatable :: theta(:),x(:),y(:),w(:),anormx(:),anormy(:),rho0(:),rho1(:),anu(:),amat(:,:),ak(:)
CHARACTER(100) :: arg,num1char,num2char,num3char,num4char,num5char


! output file for the moments

open(unit=7,file='numoments.dat',status='unknown')
open(unit=8,file='shape.dat',status='unknown')

!
! Step 1: read parameters/data from terminal
!

narg = COMMAND_ARGUMENT_COUNT()

  if (narg < 5) then
     write(*,*) ""
     write(*,*) "Goal: find polarization tensor alpha of 2d shape"
     write(*,*) ""
     write(*,*) "Usage:"
     write(*,*) "./spectral  <rad> <omega> <Ndpt>  <Nmom> <idir>"
     write(*,*) " Choose omega integer and 0<rad<1  "
     write(*,*) "     rad/omega: geometric parameters of star shape"
     write(*,*) "     Ndpt : Number of boundary points per unit length"
     write(*,*) "     Nmom: Number of moments (even integer)"
     write(*,*) "     idir: set to i=idir for alpha[i,i]"
     write(*,*) " "
     stop
  end if
  
CALL GET_COMMAND_ARGUMENT(1,num1char)   !first, read in the 3 input values
CALL GET_COMMAND_ARGUMENT(2,num2char)
CALL GET_COMMAND_ARGUMENT(3,num3char)
CALL GET_COMMAND_ARGUMENT(4,num4char)
CALL GET_COMMAND_ARGUMENT(5,num5char)

READ(num1char,*)rad     ! geometry of array
READ(num2char,*)omega     ! geometry of array
READ(num3char,*)Ndpt    ! number of boundary points /unit length
READ(num4char,*)Nmom    ! number of "moments"
READ(num5char,*)idir    ! direction of external field



wtime = omp_get_wtime ( )  ! wall time

twopi=2.d0*pi

! need to find area enclosed by the curve and its perimeter

area=0.d0
per = 0.d0
Np=1000
dt= twopi/dfloat(Np) 
do i=1,Np
   t= dt*dfloat(i)
   call dradial(rad,omega,t,r,dr)
   per = per +dsqrt(dr**2 +r**2)
   area=area + r**2
end do
per=per*dt ! perimeter
area =0.5d0*area*dt

Npt= Ndpt*per ! number of boundary points on l-th inclusion

print *,'Npt = ',Npt,' boundary points'
print *,'area = ',area
print *,'perimeter = ',per

! allocate arrays

allocate (theta(Npt),x(Npt),y(Npt),w(Npt),anormx(Npt),anormy(Npt),rho0(Npt),rho1(Npt),amat(Npt,Npt),ak(Npt))
allocate (anu(0:Nmom))

!
!     generate boundary points (x(k),y(k)) k=1..Npt
!
dte= twopi/dfloat(Npt)

do k=1,Npt ! loop over # of boundary points
   theta(k) = dte*dfloat(k-1) ! angles from 0 to 2*pi
   call d2radial(rad,omega,theta(k),r,dr,ddr)
   csp= dcos(theta(k))
   snp= dsin(theta(k))
   denom=dsqrt(dr**2 + r**2)
   anormx(k)=(dr*snp+r*csp)/denom ! x-normal at (x(k),y(k))
   anormy(k)= (-dr*csp+r*snp)/denom  ! y-normal at (x(k),y(k))
   !
   x(k)=r*csp
   y(k)=r*snp
   w(k) =   denom*dte
   ak(k)=  (2.d0*dr**2-r*ddr+r**2)/denom**3 ! curvature      
   if(idir.eq.1)then
      rho0(k) =  x(k) ! set rho0= x (for M[1,1] )
   end if
   if(idir.eq.2)then
      rho0(k) =  y(k) ! set rho0= y (for M[2,2])
   end if
   write(8,*)x(k),y(k)
end do


! 
!

print *,'series coefficients for M[IDIR,IDIR]'
print *,'idir:=',idir,';'

anu(0)=0.d0
do i=1,Npt
      anu(0)=anu(0)+ x(i)*anormx(i)*w(i)
end do

anu(0)= anu(0)/area
print *,'nu[0]:=',anu(0),';'


!
! compute interaction matrix from gradient of Green's function
!

do i=1,Npt   
   do j=1,Npt
      if(j.eq.i)amat(i,i)= -ak(i)/(4.d0*pi)
      if(j.ne.i)then
         xij=x(i)-x(j)
         yij=y(i)-y(j)
         dGdx= xij/(xij**2+yij**2)
         dGdy= yij/(xij**2+yij**2)
         amat(i,j)= 0.5d0*(dGdx*anormx(j)+dGdy*anormy(j))/pi ! G(x,y)= (1/2*pi)*n_y .(x-y)/|x-y|^2
      end if
   end do
end do

! main loop

do m=1,Nmom 

   do i=1,Npt
      rho1(i)=0.0
      do j=1,Npt
         rho1(i)=rho1(i)+amat(i,j)*w(j)*rho0(j)
      end do
      rho1(i)= rho1(i)+0.5d0*rho0(i)
   end do
!!!!!!!!!!!!!!!!!!!
! find mth moment
!!!!!!!!!!!!!!!!!!!
   anu(m)=0.0 
   do i=1,Npt
      if(idir.eq.1)then
         anu(m)=anu(m)+ anormx(i)*rho1(i)*w(i)
      end if
      if(idir.eq.2)then
         anu(m)=anu(m)+ anormy(i)*rho1(i)*w(i)
      end if
   end do
   anu(m)=anu(m)/area
   write(*,*)'nu[',m,']:=',anu(m),';'
!
! set rho0 <- rho1
!     
   do i=1,Npt
      rho0(i)=rho1(i)
   end do
!         
end do                    ! end of loop
!
!     print moments
!
write(7,*)'M := ',Nmom,';'
write(7,*)'Npt := ',Npt,';'

do k=0,Nmom
   write(7,*)'nu[',k,'] := ',anu(k),';'
end do

wtime = omp_get_wtime ( ) - wtime
write(*,*)' '
write(*,*)'Exec time  = ',wtime,'seconds'   
end program spectral

subroutine GradGreen(dx,dy,dGdx,dGdy)
  implicit double precision (a-h,o-z)
  dGdx= dx/(dx**2+dy**2)
  dGdy= dy/(dx**2+dy**2)
end subroutine GradGreen

subroutine dradial(a,omega,t,r,dr)
  implicit double precision (a-h,o-z)
  r= 1.d0 + a*dcos(omega*t)
  dr= -a*omega*dsin(omega*t)
end subroutine dradial


subroutine d2radial(a,omega,t,r,dr,ddr)
  implicit double precision (a-h,o-z)
  r= 1.d0 + a*dcos(omega*t)
  dr= -a*omega*dsin(omega*t)
  ddr=-a*omega**2*dcos(omega*t)
end subroutine d2radial


