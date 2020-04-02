
!######################################################################
! Ylm_real: Computes the normalised real spherical harmonics.
!      Here, as usual, theta in [0,pi] and phi in [0,2pi] are the
!      polar and azimuthal angles, respectively.
!######################################################################
  function Ylm_real(l,m,theta,phi)

    implicit none

    integer, intent(in)          :: m,l
    double precision, intent(in) :: theta,phi
    double precision             :: Ylm_real
    integer                      :: mabs
    double precision             :: Plm,prefac
    double precision             :: pi

    integer          :: factorial
    double precision :: plgndr
    
    pi=2.0d0*acos(0.0d0)
    
    !
    ! Check on the m and l values
    !
    if (m>abs(l)) then
       write(6,'(/,a,1x,i0,1x,i0,/)') 'Error in Yml.'&
            //' Illegal m,l values:',m,l
       stop
    endif

    !
    ! Check on the angle values
    !
    if (theta<0.0d0.or.theta>pi.or.phi<0.0d0.or.phi>2.0d0*pi) then
       write(6,'(/,a,1x,F5.2,1x,F5.2,/)') 'Error in Yml.' &
            //' Illegal theta,phi values:',m,l
       stop
    endif

    !
    ! Compute the associated Legendre polynomial P_l|m|
    ! Note that |m| is used here
    !
    mabs=abs(m)
    Plm=plgndr(l,mabs,cos(theta))

    !
    ! Prefactor
    !
    prefac=dble(2*l+1)*dble(factorial(l-mabs))/dble(factorial(l+mabs))
    if (mabs>0) prefac=prefac*2.0d0
    prefac=prefac/(4.0d0*pi)
    prefac=sqrt(prefac)
        
    !
    ! Real speherical harmonic value
    !
    if (m>=0) then
       Ylm_real=prefac*Plm*cos(m*phi)
    else
       Ylm_real=prefac*Plm*sin(mabs*phi)
    endif

    return

  end function Ylm_real
    
!######################################################################
! Ylm: Computes the normalised spherical harmonics.
!      Here, as usual, theta in [0,pi] and phi in [0,2pi] are the
!      polar and azimuthal angles, respectively.
!######################################################################
  function Ylm(l,m,theta,phi)

    implicit none

    integer, parameter :: dp = selected_real_kind(15,50)
    
    integer, intent(in)          :: m,l
    double precision, intent(in) :: theta,phi
    complex(dp)                  :: Ylm
    double precision             :: Plm,prefac
    double precision             :: pi

    integer          :: factorial
    double precision :: plgndr
    
    pi=2.0d0*acos(0.0d0)
    
    !
    ! Check on the m and l values
    !
    if (m>abs(l)) then
       write(6,'(/,a,1x,i0,1x,i0,/)') 'Error in Yml.'&
            //' Illegal m,l values:',m,l
       stop
    endif

    !
    ! Check on the angle values
    !
    if (theta<0.0d0.or.theta>pi.or.phi<0.0d0.or.phi>2.0d0*pi) then
       write(6,'(/,a,1x,F5.2,1x,F5.2,/)') 'Error in Yml.' &
            //' Illegal theta,phi values:',m,l
       stop
    endif
    
    !
    ! Compute the associated Legendre polynomial P_lm
    !
    Plm=plgndr(l,m,cos(theta))
    
    !
    ! Prefactor
    !
    prefac=(2*l+1)*factorial(l-m)/4.0d0*pi/factorial(l+m)
    prefac=sqrt(prefac)
    
    !
    ! Spherical harmonic value
    !
    Ylm=prefac*Plm*exp(-cmplx(0.0d0,1.0d0)*m*phi)
    
    return
    
  end function Ylm
  
!######################################################################
! plgndr: Numerical recipes routine for the calculation of the
!         Legendre polynomials
!----------------------------------------------------------------------
! Computes the associated Legendre polynomial Plm(x). Here m and l
! are integers satisfying 0 ≤ m ≤ l, while x lies in the range
! −1 ≤ x ≤ 1.
!######################################################################
  function plgndr(l,m,x)

    implicit none

    integer, intent(in)          :: l,m
    double precision, intent(in) :: x
    double precision             :: plgndr
    integer                      :: i,ll
    double precision             :: fact,pll,pmm,pmmp1,somx2

    if (m.lt.0.0d0.or.m.gt.l.or.abs(x).gt.1.0d0) then
       write(6,'(a)') 'Bad arguments in plgndr' 
       stop
    endif

    pmm=1.0d0
    
    if (m.gt.0) then
       somx2=sqrt((1.0d0-x)*(1.0d0+x))
       fact=1.0d0
       do i=1,m
          pmm=-pmm*fact*somx2
          fact=fact+2.0d0
       enddo
    endif

    if (l.eq.m) then
       plgndr=pmm
    else
       pmmp1=x*(2*m+1)*pmm
       if(l.eq.m+1) then
          plgndr=pmmp1
       else
          do ll=m+2,l
             pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m)
             pmm=pmmp1
             pmmp1=pll
          enddo
          plgndr=pll
       endif
    endif
    
    return
    
  end function plgndr
    
!######################################################################

  function factorial(n)

    implicit none

    integer, intent(in) :: n
    integer             :: factorial
    integer             :: i

    factorial=1

    if (n.eq.0) return

    do i=1,n
       factorial=factorial*i
    enddo
       
    return
    
  end function factorial
  
!######################################################################
  
