program test

  use sphharm
  use lebedev
  
  implicit none

  integer                       :: m,l,mp,lp
  integer                       :: npnts,n,i
  double precision              :: theta,phi,phibar
  double precision              :: sph_real  
  double precision, allocatable :: x(:),y(:),z(:),w(:)
  double precision              :: overlap
  double precision              :: pi

  pi=2.0d0*acos(0.0d0)
  
!----------------------------------------------------------------------
! Test: Check the orthonormality of the computed real spherical
!       harmonics by numerically calculating their overlaps using
!       Lebedev quadrature.
!----------------------------------------------------------------------
  !
  ! Set the number of quadrature points.
  ! This can only take one of the following valeues:
  ! 6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,590
  !
  npnts=194

  !
  ! Allocate the quadrature points and weights arrays
  !
  allocate(x(npnts),y(npnts),z(npnts),w(npnts))

  !
  ! Get the quadrature points and weights
  !
  call ld0194(x,y,z,w,n)

  !
  ! Compute the overlap integrals
  !
  ! Initialisation
  overlap=0.0d0

  ! m and l values
  l=7
  m=5
  lp=7
  mp=5
  
  ! Loop over the Lebedev quadrature points
  do i=1,n

     ! Theta
     theta=acos(z(i))

     ! Phi
     ! Note that there is some messing about here because
     ! the atan2 function returns values in (-pi,pi], whereas
     ! we need phi to be in the interval [0,pi)
     phibar=atan2(y(i),x(i))
     if (y(i)>=0.0d0) then
        phi=phibar+pi/2.0d0
     else if (y(i)<0.0d0.and.x(i)>=0.0d0) then
        phi=phibar+pi/2.0d0
     else
        phi=phibar+5.0d0*pi/2.0d0
     endif
     
     ! Conributaion to the overlap
     overlap=overlap+&
          w(i)*Ylm_real(l,m,theta,phi)*Ylm_real(lp,mp,theta,phi)
     
  enddo

  ! 4pi prefactor
  overlap=overlap*4.0d0*pi
  
  print*,''
  print*,overlap
  
end program test
