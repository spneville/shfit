subroutine fit_sphharm(lmax,npnts,distribution,coeff,rmsd)

  implicit none

  ! Input arguments
  integer, intent(in)  :: lmax
  integer, intent(in)  :: npnts
  real(8), intent(in)  :: distribution(npnts)
  real(8), intent(out) :: coeff(0:lmax,-lmax:lmax)
  real(8), intent(out) :: rmsd
  
  ! Everything else
  integer              :: n,i,l,m
  real(8), allocatable :: x(:),y(:),z(:),w(:)
  real(8), allocatable :: Ylm_vals(:,:,:)
  real(8)              :: theta,phi
  real(8)              :: modval
  
  ! External functions
  real(8) :: Ylm_real

  ! pi
  real(8) :: pi
  pi=2.0d0*acos(0.0d0)
  
!----------------------------------------------------------------------
! Allocate arrays
!----------------------------------------------------------------------
  ! Quadrature points and weights arrays
  allocate(x(npnts),y(npnts),z(npnts),w(npnts))
  x=0.0d0
  y=0.0d0
  z=0.0d0
  w=0.0d0

  ! Values of the real spherical harmonic basis functions at the
  ! quadrature points
  allocate(Ylm_vals(0:lmax,-lmax:lmax,npnts))
  Ylm_vals=0.0d0
  
!----------------------------------------------------------------------
! Get the Lebedev quadrature points and weights
! Note that currently the number of points can only be up to 590
!----------------------------------------------------------------------
  ! Get the quadrature points and weights
  select case(npnts)
  case(6)
     call ld0006(x,y,z,w,n)
  case(14)
     call ld0014(x,y,z,w,n)
  case(26)
     call ld0026(x,y,z,w,n)
  case(38)
     call ld0038(x,y,z,w,n)
  case(50)
     call ld0050(x,y,z,w,n)
  case(74)
     call ld0074(x,y,z,w,n)
  case(86)
     call ld0086(x,y,z,w,n)
  case(110)
     call ld0110(x,y,z,w,n)
  case(146)
     call ld0146(x,y,z,w,n)
  case(170)
     call ld0170(x,y,z,w,n)
  case(194)
     call ld0194(x,y,z,w,n)
  case(230)
     call ld0230(x,y,z,w,n)
  case(266)
     call ld0266(x,y,z,w,n)
  case(302)
     call ld0302(x,y,z,w,n)
  case(350)
     call ld0350(x,y,z,w,n)
  case(434)
     call ld0434(x,y,z,w,n)
  case(590)
     call ld0590(x,y,z,w,n)
  case default
     write(6,'(/,a,1x,i0,/)') 'Error in fit_sphharm.'&
          //' Invalid no. Lebedev quadrature points:',npnts
     stop
  end select

!----------------------------------------------------------------------
! Pre-calculate the values of the real spherical harmonics at the
! Lebedev quadrature points
!----------------------------------------------------------------------
  do i=1,npnts
     call xyz2ang(x(i),y(i),z(i),theta,phi)
     do l=0,lmax
        do m=-l,l
           Ylm_vals(l,m,i)=Ylm_real(l,m,theta,phi)
        enddo
     enddo
  enddo
     
!----------------------------------------------------------------------
! Calculate the coefficients for the expansion of the distribution
! function in terms of the real spherical harmonics
!----------------------------------------------------------------------
  ! Initialisation
  coeff=0.0d0
  
  ! Loop over real spherical harmonics Y_lm
  do l=0,lmax
     do m=-l,l

        ! Projection of the distribution function onto Y_lm
        do i=1,npnts
           coeff(l,m)=coeff(l,m)+&
                w(i)*distribution(i)*Ylm_vals(l,m,i)
        enddo
        coeff(l,m)=coeff(l,m)*4.0d0*pi
        
     enddo
  enddo

!----------------------------------------------------------------------
! Calculate the RMSD
!----------------------------------------------------------------------
  rmsd=0.0d0
  do i=1,npnts
     modval=0.0d0
     do l=0,lmax
        do m=-l,l
           modval=modval+coeff(l,m)*Ylm_vals(l,m,i)
        enddo
     enddo
     rmsd=rmsd+(distribution(i)-modval)**2
  enddo
  rmsd=rmsd/npnts
  rmsd=sqrt(rmsd)
  
!----------------------------------------------------------------------
! Deallocate arrays
!----------------------------------------------------------------------
  deallocate(x,y,z,w,Ylm_vals)
  
  return
  
end subroutine fit_sphharm
