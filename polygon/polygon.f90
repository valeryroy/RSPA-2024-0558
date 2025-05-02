program polygon
! usage: 
!        ./polygon  <Ngon> <Npan> <Nsub> <Nleg>  <Nmom>
!  
! Ngon: nbr of sides of polygon
! Npan: number of panels
! Nsub: number of subpanels
! Nleg: order of Gauss-Legendre quadrature on each subpanel (8, or 16)
! Nmom: Number of moments (even integer)
! 
!
! Goal: find the polarization coefficients for the polygon geometry in the form
!
! alpha[1,1] = lambda*(1- mu1[1]*lambda + mu1[2]*lambda^2+ ...)
! alpha[2,2] = lambda*(1- mu2[1]*lambda + mu2[2]*lambda^2+ ....)
!  
! where
! 
!       lambda= (epsilon-1)/(epsilon+1) (epsilon is the inclusion's permittivity normalized to that of the matrix)
!
! This code finds the coefficients mu1[1],..., mu1[Mmom]
!                                  mu2[1],..., mu2[Mmom]
!
! Output:
!         1. the coefficients mu1[k], mu2[k] in moments.dat in maple format
!         2. the boundary of the inclusion in shape.dat
!
! Author: R.V. Roy, 11/24/2023
!
!....................................................................
!  
implicit double precision (a-h,o-z)
double precision, parameter::pi=3.141592653589793238d0
double precision, allocatable :: x(:),y(:),w(:),anormx(:),anormy(:),rho0(:),rho1(:),amu(:),amat(:,:)
double precision, allocatable :: xgon(:),ygon(:)
double precision, allocatable :: Xpt(:),WW(:),Xpan(:)
CHARACTER(100) :: arg,num1char,num2char,num3char,num4char,num5char,num6char
! output file for the moments

open(unit=7,file='moments.dat',status='unknown')
open(unit=8,file='shape.dat',status='unknown')

!
! Step 1: read parameters/data from terminal
!

narg = COMMAND_ARGUMENT_COUNT()

  if (narg < 5) then
     write(*,*) ""
     write(*,*) "Polarization of N-gon"
     write(*,*) ""
     write(*,*) "Usage:"
     write(*,*) "./polygon <Ngon> <Npan> <Nsub> <Nleg>  <Nmom> "
     write(*,*) " Ngon: number of polygon sides  "
     write(*,*) " Npan: number of panels"
     write(*,*) " Nsub: number of subpanels"
     write(*,*) " Nleg: order of Gauss-Legendre quadrature on each subpanel (8, or 16)"
     write(*,*) " Nmom: Number of moments (even integer)"

     stop
  end if
  
CALL GET_COMMAND_ARGUMENT(1,num1char)   !first, read in the 6 input values
CALL GET_COMMAND_ARGUMENT(2,num2char)
CALL GET_COMMAND_ARGUMENT(3,num3char)
CALL GET_COMMAND_ARGUMENT(4,num4char)
CALL GET_COMMAND_ARGUMENT(5,num5char)

READ(num1char,*)Ngon
READ(num2char,*)Npan
READ(num3char,*)Nsub
READ(num4char,*)Nleg
READ(num5char,*)Nmom

! call makepanels to create quadrature nodes/weights
! according to paneling scheme
! Xpt: array of size Npt containing the coordinates of the quadrature nodes in [-1:1]
! W:   array of size Npt for the corresponding weights
!
Npanfin = Npan + 2*Nsub
Npt=Npanfin*Nleg

allocate(Xpt(Npt),WW(Npt),Xpan(Npanfin+1))

call makepanels(Npan,Nsub,Nleg,Xpt,Xpan,WW,Npt,Npanfin)


wtime = omp_get_wtime ( )  ! wall time

Npts=Ngon*Npt

write(7,*)'Npts := ',Npts,'; # total boundary points'

print *,'Npts = ',Npts,' total nbr boundary points'
allocate (x(Npts),y(Npts),w(Npts),anormx(Npts),anormy(Npts),rho0(Npts),rho1(Npts),amat(Npts,Npts))
allocate (xgon(0:Ngon),ygon(0:Ngon))

! allocate arrays

allocate (amu(0:Nmom))

! generate corners of N-gon
do n=0,Ngon
   phi= 2*pi*dfloat(n)/dfloat(Ngon)
   xgon(n)=cos(phi)
   ygon(n)=sin(phi)
end do

area= 0.5*Ngon*sin(2*pi/dfloat(Ngon))
print *,'area= ',area
!
!     generate boundary points (x(k),y(k)) k=1..Npt for each side n=0,Ngon-1 of !     the polygon. Repeat for idir=1,2
!
do idir=1,2
!
!
   kk=0

   do n=0,Ngon-1
         phi0=2*pi*dfloat(n)/dfloat(Ngon)
         phi1=2*pi*dfloat(n+1)/dfloat(Ngon)
         x0=cos(phi0)
         x1=cos(phi1)
         y0=sin(phi0)
         y1=sin(phi1)
      do k=1,NPT ! loop over side #1 of boundary points
         kk=kk+1
         anormx(kk)= cos(0.5*(phi0+phi1))! x-normal at (x(kk),y(kk))
         anormy(kk)= sin(0.5*(phi0+phi1))! y-normal at (x(kk),y(kk))
         x(kk) =  (x0+x1)/2.d0+ (x1-x0)/2.d0*XPT(k) ! x-boundary point
         y(kk) =  (y0+y1)/2.d0+ (y1-y0)/2.d0*XPT(k) ! y-boundary point
         w(kk) =   sin(pi/dfloat(Ngon))*WW(k) !  weights
         if(idir.eq.1)then
            rho0(kk) =  anormx(kk) ! set rho0= nx
            write(8,*)x(kk),y(kk)
         else
            rho0(kk) =  anormy(kk) ! set rho0= ny
         end if
      end do
   end do

! 
! set amu[0]= 1
!
amu(0)=0.d0 

if(idir.eq.1)then
   do i=1,NPTS
      amu(0)=amu(0)+x(i)*rho0(i)*w(i)
   end do
else
   do i=1,NPTS
      amu(0)=amu(0)+y(i)*rho0(i)*w(i)
   end do
endif

amu(0)=amu(0)/area
m=0
print *," "

if(idir.eq.1)then
   write(*,*)'mu1[',m,']:=',amu(m),';'
else
   write(*,*)'mu2[',m,']:=',amu(m),';'
end if


do i=1,Npts   
   do j=1,Npts
      if(j.eq.i)then
         amat(i,i)= 0.d0
      else
         xij=(x(j)-x(i))
         yij=(y(j)-y(i))

         call GradGreen(xij,yij,dGdx,dGdy)
         amat(i,j)=-(dGdx*anormx(i)+dGdy*anormy(i))/pi 
      end if
   end do
end do

! main loop

do m=1,Nmom 

   do i=1,Npts
      rho1(i)=0.0
      do j=1,Npts
         rho1(i)=rho1(i)+amat(i,j)*w(j)*rho0(j)
      end do
   end do
   
!!!!!!!!!!!!!!!!!!!
! find mth moment
!!!!!!!!!!!!!!!!!!!
   amu(m)=0.0
   if(idir.eq.1)then 
      do i=1,Npts
         amu(m)=amu(m)+ x(i)*rho1(i)*w(i)
      end do
   else
      do i=1,Npts
         amu(m)=amu(m)+ y(i)*rho1(i)*w(i)
      end do
   endif
   amu(m)= amu(m)/area
   if(idir.eq.1)then
      write(*,*)'mu1[',m,']:=',amu(m),';'
   else
      write(*,*)'mu2[',m,']:=',amu(m),';'
   end if
!
! set rho0 <- rho1
!     
   do i=1,Npts
      rho0(i)=rho1(i)
   end do
!         
end do                    ! end of loop
!
!     print moments
!
write(7,*)'M := ',Nmom,';'
write(7,*)'Npts := ',Npts,';'
write(7,*)'Ngon := ',Ngon,';'

if(idir.eq.1)then
   do k=0,Nmom
      write(7,*)'mu1[',k,'] := ',amu(k),';'
   end do
else
   do k=0,Nmom
      write(7,*)'mu2[',k,'] := ',amu(k),';'
   end do   
endif
!
!
end do
!
wtime = omp_get_wtime ( ) - wtime
write(*,*)' '
write(*,*)'Exec time  = ',wtime,'seconds'   
end program polygon

subroutine GradGreen(dx,dy,dGdx,dGdy)
  implicit double precision (a-h,o-z)
  dx2=dx*dx
  dy2=dy*dy
  dGdx= dx/(dx2+dy2)
  dGdy= dy/(dx2+dy2)
end subroutine GradGreen

subroutine makepanels(Npan,Nsub,N,Xpt,Xpan,W,Npts,Npanfin)
! Npan: number of panels
! Nsub: number of subpanels
! N: order of Gauss-Legendre quadrature on each subpanel (8, or 16)
! Xpt: the nodes in the interval [x0,x1]
  ! We: the weights of the quadrature
  !
  ! Note #1: Npanfin=Npan +2*Nsub
  !          Npts=Npanfin*N
  ! Note #2: Here x0=-1 and x1=1
  !
  ! Output:
  !        Xpt: array of size Npts the nodes coordinates in interval [x0,x1]
  !        Xpan: array of size Npanfin+1 containing the coordinates of the
  !              panels intervals [Xpan(i), Xpan(i+1)] i=1, Npanfin
  !        W: array of size Npts for the weights of the quadrature.
  !
implicit real*8(a-h,o-z)
double precision Xpt(Npts),W(Npts)
double precision Xpan(Npanfin+1),sinter(Npanfin+1),sinterdiff(Npanfin)
double precision WW(N),T(N)

do k=1,Npanfin+1
	sinter(k)= 0.d0
end do

do k=1,Npan+1
	sinter(k)= dfloat(k-1)/dfloat(Npan)
end do
do k=1,Npanfin
	sinterdiff(k)= 1.d0/dfloat(Npan)
end do
do i=1,Nsub
	do k=Npanfin+1,3,-1
		sinter(k)= sinter(k-1)
	end do
	sinter(2)=(sinter(1)+sinter(2))/2.d0
	do k=Npanfin,3,-1
		sinterdiff(k)= sinterdiff(k-1)
	end do
	sinterdiff(2)=sinterdiff(1)/2.d0
	sinterdiff(1)=sinterdiff(1)/2.d0
end do

do k=1,Nsub+1
	sinter(Npanfin-Nsub+k) = 1.d0- sinter(Nsub+2-k)
end do
do k=1,Nsub+2
	sinterdiff(Npanfin-Nsub-2+k) = sinterdiff(Nsub+1-k)
end do

call legendre(N,T,WW)

x0=-1.d0
x1=1.d0
delta=x1-x0
ii=0

do i=1,Npanfin+1
   xpan(i) = x0+ delta*sinter(i)
end do   

do i=0,Npanfin-1
   xA = x0+ delta*sinter(i+1)
   xB = x0+ delta*sinter(i+2)
   do j=1,N
      ii=ii+1
        Xpt(ii)= 0.5d0*(xA+xB) + 0.5d0*T(j)*(xB-xA)
        W(ii)=  0.5d0*WW(j)*(xB-xA)
	end do
end do
end subroutine makepanels
!
subroutine legendre(N,T,W)
!
! N: order of Gauss-Legendre quadrature on each subpanel (8, or 16)
! T: the nodes in the interval [-1,1]
! W: the weights of the quadrature
!
implicit real*8(a-h,o-z)
double precision T(N),W(N)
if(N.eq.8)then
   T(1)=-0.96028985649753623168
   T(2)=-0.79666647741362673959
   T(3)=-0.52553240991632898582
   T(4)=-0.18343464249564980494
   T(5)=0.18343464249564980494
   T(6)=0.52553240991632898582
   T(7)=0.79666647741362673959
   T(8)=0.96028985649753623168
   W(1)=0.10122853629037625915
   W(2)=0.22238103445337447028
   W(3)=0.31370664587788728726
   W(4)=0.36268378337836198294
   W(5)=0.36268378337836198294
   W(6)=0.31370664587788728726
   W(7)=0.22238103445337447028
   W(8)=0.10122853629037625915
else
   if(N.eq.16)then
      T(1)=-0.989400934991649932596154173450332627
      T(2)=-0.944575023073232576077988415534608345
      T(3)=-0.865631202387831743880467897712393132
      T(4)=-0.755404408355003033895101194847442268
      T(5)=-0.617876244402643748446671764048791019
      T(6)=-0.458016777657227386342419442983577574
      T(7)=-0.281603550779258913230460501460496106
      T(8)=-0.095012509837637440185319335424958063
      T(9)= 0.095012509837637440185319335424958063
      T(10)= 0.281603550779258913230460501460496106
      T(11)= 0.458016777657227386342419442983577574
      T(12)= 0.617876244402643748446671764048791019
      T(13)= 0.755404408355003033895101194847442268
      T(14)= 0.865631202387831743880467897712393132
      T(15)= 0.944575023073232576077988415534608345
      T(16)= 0.989400934991649932596154173450332627

      W(1)= 0.027152459411754094851780572456018104
      W(2)= 0.062253523938647892862843836994377694
      W(3)= 0.095158511682492784809925107602246226
      W(4)= 0.124628971255533872052476282192016420
      W(5)= 0.149595988816576732081501730547478549
      W(6)= 0.169156519395002538189312079030359962
      W(7)= 0.182603415044923588866763667969219939
      W(8)= 0.189450610455068496285396723208283105
      W(9)= 0.189450610455068496285396723208283105
      W(10)= 0.182603415044923588866763667969219939
      W(11)= 0.169156519395002538189312079030359962
      W(12)= 0.149595988816576732081501730547478549
      W(13)= 0.124628971255533872052476282192016420
      W(14)= 0.095158511682492784809925107602246226
      W(15)= 0.062253523938647892862843836994377694
      W(16)= 0.027152459411754094851780572456018104
   endif
endif
end subroutine legendre
!
!
subroutine cardinal(x,Card,N) 
  ! Compute the N  Legendre cardinal function at a point
  ! x of interval [-1:1]: the jth cardinal function satisfies
  !              C_j (x_i) = delta_ij
  ! where x_j are the poles of the Legendre polynomial of degree N
  ! The value of f(x) is given by  sum_{i=1..N} f[i] C[i](x)
  ! where f[i] = f(x_i)
  !
  implicit real*8(a-h,o-z)
  double precision Card(N)
  if(N.eq.16)then
     Card(1) = -0.1551805027118433D3*(0.1D1*x+0.9445750230732326D0)*(x+0.8656312023878317D0)*&
          (x+0.755404408355003D0)*(x+0.6178762444026437D0)*(x+0.4580167776572274D0)*&
          (x+0.2816035507792589D0)*(x+0.9501250983763744D-1)*&
          (x-0.9501250983763744D-1)*(x-0.2816035507792589D0)*(x-0.4580167776572274D0)*&
          (x-0.6178762444026437D0)*(x-0.755404408355003D0)*(x-0.8656312023878317D0)*&
          (x-0.9445750230732326D0)*(x-0.9894009349916499D0)
     Card(2) = 0.5312325937338991D3*(0.1D1*x+0.9894009349916499D0)*(x+0.8656312023878317D0)*&
          (x+0.755404408355003D0)*(x+0.6178762444026437D0)*(x+0.4580167776572274D0)*(x+0.2816035507792589D0)&
          *(x+0.9501250983763744D-1)*(x-0.9501250983763744D-1)*(x-0.2816035507792589D0)*&
          (x-0.4580167776572274D0)*(x-0.6178762444026437D0)*(x-0.755404408355003D0)*&
          (x-0.8656312023878317D0)*(x-0.9445750230732326D0)*(x-0.9894009349916499D0)
     Card(3) = -0.1001667459328063D4*(0.1D1*x+0.9894009349916499D0)*(x+0.9445750230732326D0)*&
          (x+0.755404408355003D0)*(x+0.6178762444026437D0)*(x+0.4580167776572274D0)*(x+0.2816035507792589D0)&
          *(x+0.9501250983763744D-1)*(x-0.9501250983763744D-1)*(x-0.2816035507792589D0)*(x-0.4580167776572274D0)&
          *(x-0.6178762444026437D0)*(x-0.755404408355003D0)*(x-0.8656312023878317D0)*(x-0.9445750230732326D0)&
          *(x-0.9894009349916499D0)
     Card(4) = 0.1500237583985985D4*(0.1D1*x+0.9894009349916499D0)*(x+0.9445750230732326D0)*&
          (x+0.8656312023878317D0)*(x+0.6178762444026437D0)*(x+0.4580167776572274D0)*(x+0.2816035507792589D0)&
          *(x+0.9501250983763744D-1)*(x-0.9501250983763744D-1)*(x-0.2816035507792589D0)*(x-0.4580167776572274D0)&
          *(x-0.6178762444026437D0)*(x-0.755404408355003D0)*(x-0.8656312023878317D0)*(x-0.9445750230732326D0)*&
          (x-0.9894009349916499D0)
     Card(5) = -0.1972296814457157D4*(0.1D1*x+0.9894009349916499D0)*(x+0.9445750230732326D0)*&
          (x+0.8656312023878317D0)*(x+0.755404408355003D0)*(x+0.4580167776572274D0)*(x+0.2816035507792589D0)*&
          (x+0.9501250983763744D-1)*(x-0.9501250983763744D-1)*(x-0.2816035507792589D0)*(x-0.4580167776572274D0)*&
          (x-0.6178762444026437D0)*(x-0.755404408355003D0)*(x-0.8656312023878317D0)*(x-0.9445750230732326D0)*&
          (x-0.9894009349916499D0)
     Card(6) = 0.2371134466321811D4*(0.1D1*x+0.9894009349916499D0)*(x+0.9445750230732326D0)*(x+0.8656312023878317D0)&
          *(x+0.755404408355003D0)*(x+0.6178762444026437D0)*(x+0.2816035507792589D0)*(x+0.9501250983763744D-1)*&
          (x-0.9501250983763744D-1)*(x-0.2816035507792589D0)*(x-0.4580167776572274D0)*(x-0.6178762444026437D0)*&
          (x-0.755404408355003D0)*(x-0.8656312023878317D0)*(x-0.9445750230732326D0)*(x-0.9894009349916499D0)
     Card(7) = -0.2659200099349322D4*(0.1D1*x+0.9894009349916499D0)*(x+0.9445750230732326D0)*(x+0.8656312023878317D0)&
          *(x+0.755404408355003D0)*(x+0.6178762444026437D0)*(x+0.4580167776572274D0)*(x+0.9501250983763744D-1)*&
          (x-0.9501250983763744D-1)*(x-0.2816035507792589D0)*(x-0.4580167776572274D0)*(x-0.6178762444026437D0)*&
          (x-0.755404408355003D0)*(x-0.8656312023878317D0)*(x-0.9445750230732326D0)*(x-0.9894009349916499D0)
     Card(8) = 0.2810065635073964D4*(0.1D1*x+0.9894009349916499D0)*(x+0.9445750230732326D0)*(x+0.8656312023878317D0)&
          *(x+0.755404408355003D0)*(x+0.6178762444026437D0)*(x+0.4580167776572274D0)*(x+0.2816035507792589D0)*&
          (x-0.9501250983763744D-1)*(x-0.2816035507792589D0)*(x-0.4580167776572274D0)*(x-0.6178762444026437D0)*&
          (x-0.755404408355003D0)*(x-0.8656312023878317D0)*(x-0.9445750230732326D0)*(x-0.9894009349916499D0)
     Card(9) = -0.2810065635073964D4*(0.1D1*x+0.9894009349916499D0)*(x+0.9445750230732326D0)*(x+0.8656312023878317D0)&
          *(x+0.755404408355003D0)*(x+0.6178762444026437D0)*(x+0.4580167776572274D0)*(x+0.2816035507792589D0)*&
          (x+0.9501250983763744D-1)*(x-0.2816035507792589D0)*(x-0.4580167776572274D0)*(x-0.6178762444026437D0)*&
          (x-0.755404408355003D0)*(x-0.8656312023878317D0)*(x-0.9445750230732326D0)*(x-0.9894009349916499D0)
     Card(10) = 0.2659200099349322D4*(0.1D1*x+0.9894009349916499D0)*(x+0.9445750230732326D0)*(x+0.8656312023878317D0)&
          *(x+0.755404408355003D0)*(x+0.6178762444026437D0)*(x+0.4580167776572274D0)*(x+0.2816035507792589D0)*&
          (x+0.9501250983763744D-1)*(x-0.9501250983763744D-1)*(x-0.4580167776572274D0)*(x-0.6178762444026437D0)*&
          (x-0.755404408355003D0)*(x-0.8656312023878317D0)*(x-0.9445750230732326D0)*(x-0.9894009349916499D0)
     Card(11) = -0.2371134466321811D4*(0.1D1*x+0.9894009349916499D0)*(x+0.9445750230732326D0)*&
          (x+0.8656312023878317D0)*(x+0.755404408355003D0)*(x+0.6178762444026437D0)*(x+0.4580167776572274D0)*&
          (x+0.2816035507792589D0)*(x+0.9501250983763744D-1)*(x-0.9501250983763744D-1)*(x-0.2816035507792589D0)*&
          (x-0.6178762444026437D0)*(x-0.755404408355003D0)*(x-0.8656312023878317D0)*(x-0.9445750230732326D0)*&
          (x-0.9894009349916499D0)
     Card(12) = 0.1972296814457157D4*(0.1D1*x+0.9894009349916499D0)*(x+0.9445750230732326D0)*(x+0.8656312023878317D0)&
          *(x+0.755404408355003D0)*(x+0.6178762444026437D0)*(x+0.4580167776572274D0)*(x+0.2816035507792589D0)*&
          (x+0.9501250983763744D-1)*(x-0.9501250983763744D-1)*(x-0.2816035507792589D0)*(x-0.4580167776572274D0)*&
          (x-0.755404408355003D0)*(x-0.8656312023878317D0)*(x-0.9445750230732326D0)*(x-0.9894009349916499D0)
     Card(13) = -0.1500237583985985D4*(0.1D1*x+0.9894009349916499D0)*(x+0.9445750230732326D0)*&
          (x+0.8656312023878317D0)*(x+0.755404408355003D0)*(x+0.6178762444026437D0)*(x+0.4580167776572274D0)*&
          (x+0.2816035507792589D0)*(x+0.9501250983763744D-1)*(x-0.9501250983763744D-1)*(x-0.2816035507792589D0)*&
          (x-0.4580167776572274D0)*(x-0.6178762444026437D0)*(x-0.8656312023878317D0)*(x-0.9445750230732326D0)*&
          (x-0.9894009349916499D0)
     Card(14) = 0.1001667459328063D4*(0.1D1*x+0.9894009349916499D0)*(x+0.9445750230732326D0)*(x+0.8656312023878317D0)&
          *(x+0.755404408355003D0)*(x+0.6178762444026437D0)*(x+0.4580167776572274D0)*(x+0.2816035507792589D0)*&
          (x+0.9501250983763744D-1)*(x-0.9501250983763744D-1)*(x-0.2816035507792589D0)*(x-0.4580167776572274D0)*&
          (x-0.6178762444026437D0)*(x-0.755404408355003D0)*(x-0.9445750230732326D0)*(x-0.9894009349916499D0)
     Card(15) = -0.5312325937338991D3*(0.1D1*x+0.9894009349916499D0)*(x+0.9445750230732326D0)*&
          (x+0.8656312023878317D0)*(x+0.755404408355003D0)*(x+0.6178762444026437D0)*(x+0.4580167776572274D0)*&
          (x+0.2816035507792589D0)*(x+0.9501250983763744D-1)*(x-0.9501250983763744D-1)*(x-0.2816035507792589D0)*&
          (x-0.4580167776572274D0)*(x-0.6178762444026437D0)*(x-0.755404408355003D0)*(x-0.8656312023878317D0)*&
          (x-0.9894009349916499D0)
     Card(16) = 0.1551805027118433D3*(0.1D1*x+0.9894009349916499D0)*(x+0.9445750230732326D0)*(x+0.8656312023878317D0)&
          *(x+0.755404408355003D0)*(x+0.6178762444026437D0)*(x+0.4580167776572274D0)*(x+0.2816035507792589D0)*&
          (x+0.9501250983763744D-1)*(x-0.9501250983763744D-1)*(x-0.2816035507792589D0)*(x-0.4580167776572274D0)*&
          (x-0.6178762444026437D0)*(x-0.755404408355003D0)*(x-0.8656312023878317D0)*(x-0.9445750230732326D0)
end if
!
!
if(N.eq.8)then
   Card(1) = -0.315562898850432D1*(0.1D1*x+0.7966664774136267D0)*(x+0.525532409916329D0)*(x+0.1834346424956498D0)*&
   (x-0.1834346424956498D0)*(x-0.525532409916329D0)*(x-0.7966664774136267D0)*(x-0.9602898564975362D0)
   Card(2) = 0.1013236162558063D2*(0.1D1*x+0.9602898564975362D0)*(x+0.525532409916329D0)*(x+0.1834346424956498D0)*&
        (x-0.1834346424956498D0)*(x-0.525532409916329D0)*(x-0.7966664774136267D0)*(x-0.9602898564975362D0)
   Card(3) = -0.1693945520550332D2*(0.1D1*x+0.9602898564975362D0)*(x+0.7966664774136267D0)*(x+0.1834346424956498D0)&
        *(x-0.1834346424956498D0)*(x-0.525532409916329D0)*(x-0.7966664774136267D0)*(x-0.9602898564975362D0)
   Card(4) = 0.2104530708427862D2*(0.1D1*x+0.9602898564975362D0)*(x+0.7966664774136267D0)*(x+0.525532409916329D0)&
        *(x-0.1834346424956498D0)*(x-0.525532409916329D0)*(x-0.7966664774136267D0)*(x-0.9602898564975362D0)
   Card(5) = -0.2104530708427862D2*(0.1D1*x+0.9602898564975362D0)*(x+0.7966664774136267D0)*(x+0.525532409916329D0)*&
        (x+0.1834346424956498D0)*(x-0.525532409916329D0)*(x-0.7966664774136267D0)*(x-0.9602898564975362D0)
   Card(6) = 0.1693945520550332D2*(0.1D1*x+0.9602898564975362D0)*(x+0.7966664774136267D0)*(x+0.525532409916329D0)*&
        (x+0.1834346424956498D0)*(x-0.1834346424956498D0)*(x-0.7966664774136267D0)*(x-0.9602898564975362D0)   
   Card(7) = -0.1013236162558063D2*(0.1D1*x+0.9602898564975362D0)*(x+0.7966664774136267D0)*(x+0.525532409916329D0)*&
        (x+0.1834346424956498D0)*(x-0.1834346424956498D0)*(x-0.525532409916329D0)*(x-0.9602898564975362D0)
   Card(8) = 0.315562898850432D1*(0.1D1*x+0.9602898564975362D0)*(x+0.7966664774136267D0)*(x+0.525532409916329D0)*&
        (x+0.1834346424956498D0)*(x-0.1834346424956498D0)*(x-0.525532409916329D0)*(x-0.7966664774136267D0)
end if
end subroutine cardinal
