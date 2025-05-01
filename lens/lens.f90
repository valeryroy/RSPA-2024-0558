program lens
  ! usage: 
  !        ./lens  <rad> <psi> <Npan> <Nsub> <Nleg>  <Nmom> "
  ! parameters:
  !    rad: radius of lens
  !    ipsi: half-angle of lens in radians 
  !    Npan: number of panels"
  !    Nsub: number of subpanels"
  !    Nleg: order of Gauss-Legendre quadrature on each subpanel (8, or 16)"
  !    Nmom: Number of moments (even integer)"
  !    Idir: set to 1 for sigma*11, 2 for sigma*22
  ! 
  !
  ! Goal: find the polarizations alpha[1,1] and alpha[2,2] of lens-shaped inclusion in the form
  !
  !            2*lambda*(1+ mu[1]*lambda + mu[2]*lambda[2]+ ....)
  ! where
  !       lambda= (eps-1)/(eps+1) (eps is the inclusion's permmitivity normalized to that of the matrix)
  !
  ! This code finds the coefficients mu[1],..., mu[Mmom]
  !
  ! The shape is described as y= \pm l1 ( -cos(psi) + sqrt(1 -(x*sin(psi))**2) ),  l1= l/sin(psi) )
  !
  ! or: x1= R.sin(theta)
  !
  !     x2= d- R.cos(theta),    -theta0 < theta <theta0    (theta0= pi/2 -psi)
  !
  ! Author: R.V. Roy, 06/08/2023
  !
  ! use maple pade.mw to sum the 2 series: 2*lambda*(1+ mu[1]*lambda + mu[2]*lambda[2]+ ....) 
  !....................................................................
!  
implicit double precision (a-h,o-z)
double precision, parameter::pi=3.141592653589793238d0
double precision, allocatable :: theta(:),x(:),y(:),w(:),anormx(:),anormy(:),rho0(:),rho1(:),amu(:),amat(:,:),ak(:)
double precision, allocatable :: Xpt(:),WW(:),Xpan(:)
CHARACTER(100) :: arg,num1char,num2char,num3char,num4char,num5char,num6char
! output file for the moments

open(unit=7,file='moments.dat',status='unknown')
open(unit=8,file='shape.dat',status='unknown')
open(unit=9,file='modes.dat',status='unknown')

!
! Step 1: read parameters/data from terminal
!

narg = COMMAND_ARGUMENT_COUNT()

  if (narg < 6) then
     write(*,*) "lambda formulation: find the coefficients of expansion"
     write(*,*) " "
     write(*,*) "            2*phi*lambda*(1+ mu[1]*lambda + mu[2]*lambda[2]+ ....)"
     write(*,*) ""
     write(*,*) "Usage:"
     write(*,*) "./lens  <rad> <psi> <Npan> <Nsub> <Nleg> <Nmom>"
     write(*,*) " "
     write(*,*) " rad: radius lens"
     write(*,*) " psi: half-angle of lens (in degrees)"
     write(*,*) " Npan: number of panels"
     write(*,*) " Nsub: number of subpanels"
     write(*,*) " Nleg: order of Gauss-Legendre quadrature on each subpanel (8, or 16)"
     write(*,*) " Nmom: Number of moments (even integer)"
     write(*,*) " "
     write(*,*) " call: maple  polarization.mw"
     stop
  end if
  
  CALL GET_COMMAND_ARGUMENT(1,num1char)   !first, read in the 6 input values
  CALL GET_COMMAND_ARGUMENT(2,num2char)
  CALL GET_COMMAND_ARGUMENT(3,num3char)
  CALL GET_COMMAND_ARGUMENT(4,num4char)
  CALL GET_COMMAND_ARGUMENT(5,num5char)
  CALL GET_COMMAND_ARGUMENT(6,num6char)

  READ(num1char,*)rad
  READ(num2char,*)psi
  READ(num3char,*)Npan
  READ(num4char,*)Nsub
  READ(num5char,*)Nleg
  READ(num6char,*)Nmom ! number of "moments"

! call makepanels to create quadrature nodes/weights
! according to paneling scheme
! Xpt: array of size Npt containing the coordinates of the quadrature nodes in [-1:1]
! W:   array of size Npt for the corresponding weights
!
  Npanfin = Npan + 2*Nsub
  Npt=Npanfin*Nleg

  allocate(Xpt(Npt),WW(Npt),Xpan(Npanfin+1))

  call makepanels(Npan,Nsub,Nleg,Xpt,Xpan,WW,Npt,Npanfin)

  !print *,Xpt
  !print *,WW
  kout= 100
  
  psi=(psi)*pi/180.d0
  dd= rad*cos(psi)
  thetamax=  psi
  area= rad**2*(2.d0*psi-sin(2.d0*psi))

  wtime = omp_get_wtime ( )  ! wall time

  Npts=2*Npt

  print *,'Npts = ',Npts,' total nbr boundary points'
  print *,'area= ',area
  
! allocate arrays

  allocate (amu(0:Nmom))
!
!     generate boundary points (x(k),y(k)) k=1..NPTS
!

  allocate (theta(Npts),x(Npts),y(Npts),w(Npts),anormx(Npts),anormy(Npts),rho0(Npts),rho1(Npts),ak(Npts),amat(Npts,Npts))

write(7,*)'RAD := ',rad,';'
write(7,*)'PHI := ',ipsi,';'
write(7,*)'NPAN := ',npan,';'
write(7,*)'NSUB := ',nsub,';'
write(7,*)'NLEG := ',nleg,';'
write(7,*)'M := ',Nmom,';'
write(7,*)'Npts := ',Npts,'; # total boundary points'

do Idir=1,2 ! loop over direction
   kk=0
  do k=1,NPT ! loop over side #1 of boundary points
     kk=kk+1
     theta(kk)=thetamax*XPT(k)
     anormx(kk)= sin(theta(kk)) ! x-normal at (x(kk),y(kk))
     anormy(kk)= cos(theta(kk))  ! y-normal at (x(kk),y(kk))
     x(kk) =  rad*sin(theta(kk))  ! x-boundary point
     y(kk) =  rad*cos(theta(kk))  ! y-boundary point
     w(kk) =   rad*thetamax*WW(k) !  weight
     if(idir.eq.1)then
        rho0(kk) =  anormx(kk) ! set rho0= n_x
     else
        rho0(kk) =  anormy(kk) ! set rho0= n_y
     endif
     ak(kk)=-1.d0/rad
     write(8,*)x(kk),y(kk)
  end do

  do k=1,NPT ! loop over side #2 of boundary points
     kk=kk+1
     theta(kk)=thetamax*XPT(NPT-k+1)
     anormx(kk)=sin(theta(kk)) ! x-normal at (x(kk),y(kk))
     anormy(kk)= -cos(theta(kk)) ! y-normal at (x(kk),y(kk))
     x(kk) =  rad*sin(theta(kk))    ! x-boundary point
     y(kk) = 2.d0*dd-rad*cos(theta(kk))    ! y-boundary point
     w(kk) =   rad*thetamax*WW(NPT-k+1) !  weight
     if(idir.eq.1)then
        rho0(kk) =  anormx(kk) ! set rho0= n_x
     else
        rho0(kk) =  anormy(kk) ! set rho0= n_y
     endif
     ak(kk)=-1.d0/rad
     write(8,*)x(kk),y(kk)
  end do
  
write(9,*)rho0(kout)


  NPTS=kk
  print *,'Npts = ',Npts,' total boundary points'
! 
! set amu[0]= 1
!
  amu(0)=0.d0
  m=0
  
  if(idir.eq.1)then
     do i=1,Npts
        amu(0)=amu(0)+ x(i)*rho0(i)*w(i)
     end do
  else
     do i=1,Npts
        amu(0)=amu(0)+ y(i)*rho0(i)*w(i)
     end do
   endif
   amu(0)= amu(0)/area
   if(idir.eq.1)then     
      write(*,*)'mu1[',m,']:=',amu(m),';'
   else
      write(*,*)'mu2[',m,']:=',amu(m),';'
   end if

  do i=1,Npts   
     do j=1,Npts
        if(j.eq.i)then
           amat(i,i)= -ak(i)/(2.0*pi)  ! factor of two
        else
           xij=(x(j)-x(i))
           yij=(y(j)-y(i))
!
           call GradGreen(xij,yij,dGdx,dGdy)
           amat(i,j)=-(dGdx*anormx(i)+dGdy*anormy(i))/pi ! 2*G(x,y)= (1/pi)*n_x .(y-x)/|x-y|^2
        end if
     end do
  end do

! main loop
  id=1

  do m=1,Nmom 

     do i=1,Npts
        rho1(i)=0.d0
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
   amu(m)=dfloat(id)*amu(m)/area
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
   write(9,*)rho0(kout)

end do                    ! end of loop
!
!     print moments
!
write(7,*)'IDIR := ',idir,';'

if(idir.eq.1)then
   do k=0,Nmom
      write(7,*)'mu1[',k,'] := ',amu(k),';'
   end do
else
   do k=0,Nmom
      write(7,*)'mu2[',k,'] := ',amu(k),';'
   end do
end if
!
   write(9,*)
end do ! idir loop

wtime = omp_get_wtime ( ) - wtime
write(*,*)' '
write(*,*)'Exec time  = ',wtime,'seconds'   
end program lens

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
