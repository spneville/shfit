module lebedev
  
contains
  
!######################################################################
  
  subroutine gen_oh(code, num, x, y, z, w, a, b, v)
    implicit logical(a-z)
    double precision x(*),y(*),z(*),w(*)
    double precision a,b,v
    integer code
    integer num
    double precision c
!hvd
!hvd   This subroutine is part of a set of subroutines that generate
!hvd   Lebedev grids [1-6] for integration on a sphere. The original 
!hvd   C-code [1] was kindly provided by Dr. Dmitri N. Laikov and 
!hvd   translated into fortran by Dr. Christoph van Wuellen.
!hvd   This subroutine was translated from C to fortran77 by hand.
!hvd
!hvd   Users of this code are asked to include reference [1] in their
!hvd   publications, and in the user- and programmers-manuals 
!hvd   describing their codes.
!hvd
!hvd   This code was distributed through CCL (http://www.ccl.net/).
!hvd
!hvd   [1] V.I. Lebedev, and D.N. Laikov
!hvd       "A quadrature formula for the sphere of the 131st
!hvd        algebraic order of accuracy"
!hvd       Doklady Mathematics, Vol. 59, No. 3, 1999, pp. 477-481.
!hvd
!hvd   [2] V.I. Lebedev
!hvd       "A quadrature formula for the sphere of 59th algebraic
!hvd        order of accuracy"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 50, 1995, pp. 283-286. 
!hvd
!hvd   [3] V.I. Lebedev, and A.L. Skorokhodov
!hvd       "Quadrature formulas of orders 41, 47, and 53 for the sphere"
!hvd       Russian Acad. Sci. Dokl. Math., Vol. 45, 1992, pp. 587-592. 
!hvd
!hvd   [4] V.I. Lebedev
!hvd       "Spherical quadrature formulas exact to orders 25-29"
!hvd       Siberian Mathematical Journal, Vol. 18, 1977, pp. 99-107. 
!hvd
!hvd   [5] V.I. Lebedev
!hvd       "Quadratures on a sphere"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 16,
!hvd       1976, pp. 10-24. 
!hvd
!hvd   [6] V.I. Lebedev
!hvd       "Values of the nodes and weights of ninth to seventeenth 
!hvd        order Gauss-Markov quadrature formulae invariant under the
!hvd        octahedron group with inversion"
!hvd       Computational Mathematics and Mathematical Physics, Vol. 15,
!hvd       1975, pp. 44-51.
!hvd
!vw
!vw    Given a point on a sphere (specified by a and b), generate all
!vw    the equivalent points under Oh symmetry, making grid points with
!vw    weight v.
!vw    The variable num is increased by the number of different points
!vw    generated.
!vw
!vw    Depending on code, there are 6...48 different but equivalent
!vw    points.
!vw
!vw    code=1:   (0,0,1) etc                                (  6 points)
!vw    code=2:   (0,a,a) etc, a=1/sqrt(2)                   ( 12 points)
!vw    code=3:   (a,a,a) etc, a=1/sqrt(3)                   (  8 points)
!vw    code=4:   (a,a,b) etc, b=sqrt(1-2 a^2)               ( 24 points)
!vw    code=5:   (a,b,0) etc, b=sqrt(1-a^2), a input        ( 24 points)
!vw    code=6:   (a,b,c) etc, c=sqrt(1-a^2-b^2), a/b input  ( 48 points)
!vw
       goto (1,2,3,4,5,6) code
       write (6,*) 'Gen_Oh: Invalid Code'
       stop 
    1  continue
       a=1.0d0
       x(1) =  a
       y(1) =  0.0d0
       z(1) =  0.0d0
       w(1) =  v
       x(2) = -a
       y(2) =  0.0d0
       z(2) =  0.0d0
       w(2) =  v
       x(3) =  0.0d0
       y(3) =  a
       z(3) =  0.0d0
       w(3) =  v
       x(4) =  0.0d0
       y(4) = -a
       z(4) =  0.0d0
       w(4) =  v
       x(5) =  0.0d0
       y(5) =  0.0d0
       z(5) =  a
       w(5) =  v
       x(6) =  0.0d0
       y(6) =  0.0d0
       z(6) = -a
       w(6) =  v
       num=num+6
       return
!vw
    2  continue
       a=sqrt(0.5d0)
       x( 1) =  0d0
       y( 1) =  a
       z( 1) =  a
       w( 1) =  v
       x( 2) =  0d0
       y( 2) = -a
       z( 2) =  a
       w( 2) =  v
       x( 3) =  0d0
       y( 3) =  a
       z( 3) = -a
       w( 3) =  v
       x( 4) =  0d0
       y( 4) = -a
       z( 4) = -a
       w( 4) =  v
       x( 5) =  a
       y( 5) =  0d0
       z( 5) =  a
       w( 5) =  v
       x( 6) = -a
       y( 6) =  0d0
       z( 6) =  a
       w( 6) =  v
       x( 7) =  a
       y( 7) =  0d0
       z( 7) = -a
       w( 7) =  v
       x( 8) = -a
       y( 8) =  0d0
       z( 8) = -a
       w( 8) =  v
       x( 9) =  a
       y( 9) =  a
       z( 9) =  0d0
       w( 9) =  v
       x(10) = -a
       y(10) =  a
       z(10) =  0d0
       w(10) =  v
       x(11) =  a
       y(11) = -a
       z(11) =  0d0
       w(11) =  v
       x(12) = -a
       y(12) = -a
       z(12) =  0d0
       w(12) =  v
       num=num+12
       return
!vw
    3  continue
       a = sqrt(1d0/3d0)
       x(1) =  a
       y(1) =  a
       z(1) =  a
       w(1) =  v
       x(2) = -a
       y(2) =  a
       z(2) =  a
       w(2) =  v
       x(3) =  a
       y(3) = -a
       z(3) =  a
       w(3) =  v
       x(4) = -a
       y(4) = -a
       z(4) =  a
       w(4) =  v
       x(5) =  a
       y(5) =  a
       z(5) = -a
       w(5) =  v
       x(6) = -a
       y(6) =  a
       z(6) = -a
       w(6) =  v
       x(7) =  a
       y(7) = -a
       z(7) = -a
       w(7) =  v
       x(8) = -a
       y(8) = -a
       z(8) = -a
       w(8) =  v
       num=num+8
       return
!vw
    4  continue
       b = sqrt(1d0 - 2d0*a*a)
       x( 1) =  a
       y( 1) =  a
       z( 1) =  b
       w( 1) =  v
       x( 2) = -a
       y( 2) =  a
       z( 2) =  b
       w( 2) =  v
       x( 3) =  a
       y( 3) = -a
       z( 3) =  b
       w( 3) =  v
       x( 4) = -a
       y( 4) = -a
       z( 4) =  b
       w( 4) =  v
       x( 5) =  a
       y( 5) =  a
       z( 5) = -b
       w( 5) =  v
       x( 6) = -a
       y( 6) =  a
       z( 6) = -b
       w( 6) =  v
       x( 7) =  a
       y( 7) = -a
       z( 7) = -b
       w( 7) =  v
       x( 8) = -a
       y( 8) = -a
       z( 8) = -b
       w( 8) =  v
       x( 9) =  a
       y( 9) =  b
       z( 9) =  a
       w( 9) =  v
       x(10) = -a
       y(10) =  b
       z(10) =  a
       w(10) =  v
       x(11) =  a
       y(11) = -b
       z(11) =  a
       w(11) =  v
       x(12) = -a
       y(12) = -b
       z(12) =  a
       w(12) =  v
       x(13) =  a
       y(13) =  b
       z(13) = -a
       w(13) =  v
       x(14) = -a
       y(14) =  b
       z(14) = -a
       w(14) =  v
       x(15) =  a
       y(15) = -b
       z(15) = -a
       w(15) =  v
       x(16) = -a
       y(16) = -b
       z(16) = -a
       w(16) =  v
       x(17) =  b
       y(17) =  a
       z(17) =  a
       w(17) =  v
       x(18) = -b
       y(18) =  a
       z(18) =  a
       w(18) =  v
       x(19) =  b
       y(19) = -a
       z(19) =  a
       w(19) =  v
       x(20) = -b
       y(20) = -a
       z(20) =  a
       w(20) =  v
       x(21) =  b
       y(21) =  a
       z(21) = -a
       w(21) =  v
       x(22) = -b
       y(22) =  a
       z(22) = -a
       w(22) =  v
       x(23) =  b
       y(23) = -a
       z(23) = -a
       w(23) =  v
       x(24) = -b
       y(24) = -a
       z(24) = -a
       w(24) =  v
       num=num+24
       return
!vw
    5  continue
       b=sqrt(1d0-a*a)
       x( 1) =  a
       y( 1) =  b
       z( 1) =  0d0
       w( 1) =  v
       x( 2) = -a
       y( 2) =  b
       z( 2) =  0d0
       w( 2) =  v
       x( 3) =  a
       y( 3) = -b
       z( 3) =  0d0
       w( 3) =  v
       x( 4) = -a
       y( 4) = -b
       z( 4) =  0d0
       w( 4) =  v
       x( 5) =  b
       y( 5) =  a
       z( 5) =  0d0
       w( 5) =  v
       x( 6) = -b
       y( 6) =  a
       z( 6) =  0d0
       w( 6) =  v
       x( 7) =  b
       y( 7) = -a
       z( 7) =  0d0
       w( 7) =  v
       x( 8) = -b
       y( 8) = -a
       z( 8) =  0d0
       w( 8) =  v
       x( 9) =  a
       y( 9) =  0d0
       z( 9) =  b
       w( 9) =  v
       x(10) = -a
       y(10) =  0d0
       z(10) =  b
       w(10) =  v
       x(11) =  a
       y(11) =  0d0
       z(11) = -b
       w(11) =  v
       x(12) = -a
       y(12) =  0d0
       z(12) = -b
       w(12) =  v
       x(13) =  b
       y(13) =  0d0
       z(13) =  a
       w(13) =  v
       x(14) = -b
       y(14) =  0d0
       z(14) =  a
       w(14) =  v
       x(15) =  b
       y(15) =  0d0
       z(15) = -a
       w(15) =  v
       x(16) = -b
       y(16) =  0d0
       z(16) = -a
       w(16) =  v
       x(17) =  0d0
       y(17) =  a
       z(17) =  b
       w(17) =  v
       x(18) =  0d0
       y(18) = -a
       z(18) =  b
       w(18) =  v
       x(19) =  0d0
       y(19) =  a
       z(19) = -b
       w(19) =  v
       x(20) =  0d0
       y(20) = -a
       z(20) = -b
       w(20) =  v
       x(21) =  0d0
       y(21) =  b
       z(21) =  a
       w(21) =  v
       x(22) =  0d0
       y(22) = -b
       z(22) =  a
       w(22) =  v
       x(23) =  0d0
       y(23) =  b
       z(23) = -a
       w(23) =  v
       x(24) =  0d0
       y(24) = -b
       z(24) = -a
       w(24) =  v
       num=num+24
       return
!vw
    6  continue
       c=sqrt(1d0 - a*a - b*b)
       x( 1) =  a
       y( 1) =  b
       z( 1) =  c
       w( 1) =  v
       x( 2) = -a
       y( 2) =  b
       z( 2) =  c
       w( 2) =  v
       x( 3) =  a
       y( 3) = -b
       z( 3) =  c
       w( 3) =  v
       x( 4) = -a
       y( 4) = -b
       z( 4) =  c
       w( 4) =  v
       x( 5) =  a
       y( 5) =  b
       z( 5) = -c
       w( 5) =  v
       x( 6) = -a
       y( 6) =  b
       z( 6) = -c
       w( 6) =  v
       x( 7) =  a
       y( 7) = -b
       z( 7) = -c
       w( 7) =  v
       x( 8) = -a
       y( 8) = -b
       z( 8) = -c
       w( 8) =  v
       x( 9) =  a
       y( 9) =  c
       z( 9) =  b
       w( 9) =  v
       x(10) = -a
       y(10) =  c
       z(10) =  b
       w(10) =  v
       x(11) =  a
       y(11) = -c
       z(11) =  b
       w(11) =  v
       x(12) = -a
       y(12) = -c
       z(12) =  b
       w(12) =  v
       x(13) =  a
       y(13) =  c
       z(13) = -b
       w(13) =  v
       x(14) = -a
       y(14) =  c
       z(14) = -b
       w(14) =  v
       x(15) =  a
       y(15) = -c
       z(15) = -b
       w(15) =  v
       x(16) = -a
       y(16) = -c
       z(16) = -b
       w(16) =  v
       x(17) =  b
       y(17) =  a
       z(17) =  c
       w(17) =  v
       x(18) = -b
       y(18) =  a
       z(18) =  c
       w(18) =  v
       x(19) =  b
       y(19) = -a
       z(19) =  c
       w(19) =  v
       x(20) = -b
       y(20) = -a
       z(20) =  c
       w(20) =  v
       x(21) =  b
       y(21) =  a
       z(21) = -c
       w(21) =  v
       x(22) = -b
       y(22) =  a
       z(22) = -c
       w(22) =  v
       x(23) =  b
       y(23) = -a
       z(23) = -c
       w(23) =  v
       x(24) = -b
       y(24) = -a
       z(24) = -c
       w(24) =  v
       x(25) =  b
       y(25) =  c
       z(25) =  a
       w(25) =  v
       x(26) = -b
       y(26) =  c
       z(26) =  a
       w(26) =  v
       x(27) =  b
       y(27) = -c
       z(27) =  a
       w(27) =  v
       x(28) = -b
       y(28) = -c
       z(28) =  a
       w(28) =  v
       x(29) =  b
       y(29) =  c
       z(29) = -a
       w(29) =  v
       x(30) = -b
       y(30) =  c
       z(30) = -a
       w(30) =  v
       x(31) =  b
       y(31) = -c
       z(31) = -a
       w(31) =  v
       x(32) = -b
       y(32) = -c
       z(32) = -a
       w(32) =  v
       x(33) =  c
       y(33) =  a
       z(33) =  b
       w(33) =  v
       x(34) = -c
       y(34) =  a
       z(34) =  b
       w(34) =  v
       x(35) =  c
       y(35) = -a
       z(35) =  b
       w(35) =  v
       x(36) = -c
       y(36) = -a
       z(36) =  b
       w(36) =  v
       x(37) =  c
       y(37) =  a
       z(37) = -b
       w(37) =  v
       x(38) = -c
       y(38) =  a
       z(38) = -b
       w(38) =  v
       x(39) =  c
       y(39) = -a
       z(39) = -b
       w(39) =  v
       x(40) = -c
       y(40) = -a
       z(40) = -b
       w(40) =  v
       x(41) =  c
       y(41) =  b
       z(41) =  a
       w(41) =  v
       x(42) = -c
       y(42) =  b
       z(42) =  a
       w(42) =  v
       x(43) =  c
       y(43) = -b
       z(43) =  a
       w(43) =  v
       x(44) = -c
       y(44) = -b
       z(44) =  a
       w(44) =  v
       x(45) =  c
       y(45) =  b
       z(45) = -a
       w(45) =  v
       x(46) = -c
       y(46) =  b
       z(46) = -a
       w(46) =  v
       x(47) =  c
       y(47) = -b
       z(47) = -a
       w(47) =  v
       x(48) = -c
       y(48) = -b
       z(48) = -a
       w(48) =  v
       num=num+48
       return
     end subroutine gen_oh

!######################################################################
     
       SUBROUTINE LD0006(X,Y,Z,W,N)
       DOUBLE PRECISION X(   6)
       DOUBLE PRECISION Y(   6)
       DOUBLE PRECISION Z(   6)
       DOUBLE PRECISION W(   6)
       INTEGER N
       DOUBLE PRECISION A,B,V
!VW
!VW    LEBEDEV    6-POINT ANGULAR GRID
!VW
       N=1
       V=0.1666666666666667D+0
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
     END SUBROUTINE LD0006

!######################################################################

     SUBROUTINE LD0014(X,Y,Z,W,N)
       DOUBLE PRECISION X(  14)
       DOUBLE PRECISION Y(  14)
       DOUBLE PRECISION Z(  14)
       DOUBLE PRECISION W(  14)
       INTEGER N
       DOUBLE PRECISION A,B,V
!VW
!VW    LEBEDEV   14-POINT ANGULAR GRID
!VW
       N=1
       V=0.6666666666666667D-1
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.7500000000000000D-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
     END SUBROUTINE LD0014

!######################################################################

     SUBROUTINE LD0026(X,Y,Z,W,N)
       DOUBLE PRECISION X(  26)
       DOUBLE PRECISION Y(  26)
       DOUBLE PRECISION Z(  26)
       DOUBLE PRECISION W(  26)
       INTEGER N
       DOUBLE PRECISION A,B,V
!VW
!VW    LEBEDEV   26-POINT ANGULAR GRID
!VW
       N=1
       V=0.4761904761904762D-1
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3809523809523810D-1
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3214285714285714D-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
     END SUBROUTINE LD0026

!######################################################################

     SUBROUTINE LD0038(X,Y,Z,W,N)
       DOUBLE PRECISION X(  38)
       DOUBLE PRECISION Y(  38)
       DOUBLE PRECISION Z(  38)
       DOUBLE PRECISION W(  38)
       INTEGER N
       DOUBLE PRECISION A,B,V
!VW
!VW    LEBEDEV   38-POINT ANGULAR GRID
!VW
       N=1
       V=0.9523809523809524D-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3214285714285714D-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4597008433809831D+0
       V=0.2857142857142857D-1
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
     END SUBROUTINE LD0038

!######################################################################
     
     SUBROUTINE LD0050(X,Y,Z,W,N)
       DOUBLE PRECISION X(  50)
       DOUBLE PRECISION Y(  50)
       DOUBLE PRECISION Z(  50)
       DOUBLE PRECISION W(  50)
       INTEGER N
       DOUBLE PRECISION A,B,V
!VW
!VW    LEBEDEV   50-POINT ANGULAR GRID
!VW
       N=1
       V=0.1269841269841270D-1
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2257495590828924D-1
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2109375000000000D-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3015113445777636D+0
       V=0.2017333553791887D-1
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
     END SUBROUTINE LD0050

!######################################################################

     SUBROUTINE LD0074(X,Y,Z,W,N)
       DOUBLE PRECISION X(  74)
       DOUBLE PRECISION Y(  74)
       DOUBLE PRECISION Z(  74)
       DOUBLE PRECISION W(  74)
       INTEGER N
       DOUBLE PRECISION A,B,V
!VW
!VW    LEBEDEV   74-POINT ANGULAR GRID
!VW
       N=1
       V=0.5130671797338464D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1660406956574204D-1
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=-0.2958603896103896D-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4803844614152614D+0
       V=0.2657620708215946D-1
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3207726489807764D+0
       V=0.1652217099371571D-1
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
     END SUBROUTINE LD0074

!######################################################################

     SUBROUTINE LD0086(X,Y,Z,W,N)
       DOUBLE PRECISION X(  86)
       DOUBLE PRECISION Y(  86)
       DOUBLE PRECISION Z(  86)
       DOUBLE PRECISION W(  86)
       INTEGER N
       DOUBLE PRECISION A,B,V
!VW
!VW    LEBEDEV   86-POINT ANGULAR GRID
!VW
       N=1
       V=0.1154401154401154D-1
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1194390908585628D-1
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3696028464541502D+0
       V=0.1111055571060340D-1
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6943540066026664D+0
       V=0.1187650129453714D-1
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3742430390903412D+0
       V=0.1181230374690448D-1
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
     END SUBROUTINE LD0086

!######################################################################

     SUBROUTINE LD0110(X,Y,Z,W,N)
       DOUBLE PRECISION X( 110)
       DOUBLE PRECISION Y( 110)
       DOUBLE PRECISION Z( 110)
       DOUBLE PRECISION W( 110)
       INTEGER N
       DOUBLE PRECISION A,B,V
!VW
!VW    LEBEDEV  110-POINT ANGULAR GRID
!VW
       N=1
       V=0.3828270494937162D-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.9793737512487512D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1851156353447362D+0
       V=0.8211737283191111D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6904210483822922D+0
       V=0.9942814891178103D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3956894730559419D+0
       V=0.9595471336070963D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4783690288121502D+0
       V=0.9694996361663028D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
     END SUBROUTINE LD0110

!######################################################################

     SUBROUTINE LD0146(X,Y,Z,W,N)
       DOUBLE PRECISION X( 146)
       DOUBLE PRECISION Y( 146)
       DOUBLE PRECISION Z( 146)
       DOUBLE PRECISION W( 146)
       INTEGER N
       DOUBLE PRECISION A,B,V
!VW
!VW    LEBEDEV  146-POINT ANGULAR GRID
!VW
       N=1
       V=0.5996313688621381D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.7372999718620756D-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.7210515360144488D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6764410400114264D+0
       V=0.7116355493117555D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4174961227965453D+0
       V=0.6753829486314477D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1574676672039082D+0
       V=0.7574394159054034D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1403553811713183D+0
       B=0.4493328323269557D+0
       V=0.6991087353303262D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
     END SUBROUTINE LD0146

!######################################################################

     SUBROUTINE LD0170(X,Y,Z,W,N)
       DOUBLE PRECISION X( 170)
       DOUBLE PRECISION Y( 170)
       DOUBLE PRECISION Z( 170)
       DOUBLE PRECISION W( 170)
       INTEGER N
       DOUBLE PRECISION A,B,V
!VW
!VW    LEBEDEV  170-POINT ANGULAR GRID
!VW
       N=1
       V=0.5544842902037365D-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.6071332770670752D-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.6383674773515093D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2551252621114134D+0
       V=0.5183387587747790D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6743601460362766D+0
       V=0.6317929009813725D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4318910696719410D+0
       V=0.6201670006589077D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2613931360335988D+0
       V=0.5477143385137348D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4990453161796037D+0
       B=0.1446630744325115D+0
       V=0.5968383987681156D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
     END SUBROUTINE LD0170

!######################################################################

     SUBROUTINE LD0194(X,Y,Z,W,N)
       DOUBLE PRECISION X( 194)
       DOUBLE PRECISION Y( 194)
       DOUBLE PRECISION Z( 194)
       DOUBLE PRECISION W( 194)
       INTEGER N
       DOUBLE PRECISION A,B,V
!VW
!VW    LEBEDEV  194-POINT ANGULAR GRID
!VW
       N=1
       V=0.1782340447244611D-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.5716905949977102D-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.5573383178848738D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6712973442695226D+0
       V=0.5608704082587997D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2892465627575439D+0
       V=0.5158237711805383D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4446933178717437D+0
       V=0.5518771467273614D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1299335447650067D+0
       V=0.4106777028169394D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3457702197611283D+0
       V=0.5051846064614808D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1590417105383530D+0
       B=0.8360360154824589D+0
       V=0.5530248916233094D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
     END SUBROUTINE LD0194

!######################################################################

     SUBROUTINE LD0230(X,Y,Z,W,N)
       DOUBLE PRECISION X( 230)
       DOUBLE PRECISION Y( 230)
       DOUBLE PRECISION Z( 230)
       DOUBLE PRECISION W( 230)
       INTEGER N
       DOUBLE PRECISION A,B,V
!VW
!VW    LEBEDEV  230-POINT ANGULAR GRID
!VW
       N=1
       V=-0.5522639919727325D-1
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4450274607445226D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4492044687397611D+0
       V=0.4496841067921404D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2520419490210201D+0
       V=0.5049153450478750D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6981906658447242D+0
       V=0.3976408018051883D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6587405243460960D+0
       V=0.4401400650381014D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4038544050097660D-1
       V=0.1724544350544401D-1
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5823842309715585D+0
       V=0.4231083095357343D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3545877390518688D+0
       V=0.5198069864064399D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2272181808998187D+0
       B=0.4864661535886647D+0
       V=0.4695720972568883D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
     END SUBROUTINE LD0230

!######################################################################

     SUBROUTINE LD0266(X,Y,Z,W,N)
       DOUBLE PRECISION X( 266)
       DOUBLE PRECISION Y( 266)
       DOUBLE PRECISION Z( 266)
       DOUBLE PRECISION W( 266)
       INTEGER N
       DOUBLE PRECISION A,B,V
!VW
!VW    LEBEDEV  266-POINT ANGULAR GRID
!VW
       N=1
       V=-0.1313769127326952D-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=-0.2522728704859336D-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.4186853881700583D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7039373391585475D+0
       V=0.5315167977810885D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1012526248572414D+0
       V=0.4047142377086219D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4647448726420539D+0
       V=0.4112482394406990D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3277420654971629D+0
       V=0.3595584899758782D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6620338663699974D+0
       V=0.4256131351428158D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8506508083520399D+0
       V=0.4229582700647240D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3233484542692899D+0
       B=0.1153112011009701D+0
       V=0.4080914225780505D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2314790158712601D+0
       B=0.5244939240922365D+0
       V=0.4071467593830964D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
     END SUBROUTINE LD0266

!######################################################################

     SUBROUTINE LD0302(X,Y,Z,W,N)
       DOUBLE PRECISION X( 302)
       DOUBLE PRECISION Y( 302)
       DOUBLE PRECISION Z( 302)
       DOUBLE PRECISION W( 302)
       INTEGER N
       DOUBLE PRECISION A,B,V
!VW
!VW    LEBEDEV  302-POINT ANGULAR GRID
!VW
       N=1
       V=0.8545911725128148D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3599119285025571D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3515640345570105D+0
       V=0.3449788424305883D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6566329410219612D+0
       V=0.3604822601419882D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4729054132581005D+0
       V=0.3576729661743367D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9618308522614784D-1
       V=0.2352101413689164D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2219645236294178D+0
       V=0.3108953122413675D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7011766416089545D+0
       V=0.3650045807677255D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2644152887060663D+0
       V=0.2982344963171804D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5718955891878961D+0
       V=0.3600820932216460D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2510034751770465D+0
       B=0.8000727494073952D+0
       V=0.3571540554273387D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1233548532583327D+0
       B=0.4127724083168531D+0
       V=0.3392312205006170D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
     END SUBROUTINE LD0302

!######################################################################

     SUBROUTINE LD0350(X,Y,Z,W,N)
       DOUBLE PRECISION X( 350)
       DOUBLE PRECISION Y( 350)
       DOUBLE PRECISION Z( 350)
       DOUBLE PRECISION W( 350)
       INTEGER N
       DOUBLE PRECISION A,B,V
!VW
!VW    LEBEDEV  350-POINT ANGULAR GRID
!VW
       N=1
       V=0.3006796749453936D-2
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.3050627745650771D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7068965463912316D+0
       V=0.1621104600288991D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4794682625712025D+0
       V=0.3005701484901752D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1927533154878019D+0
       V=0.2990992529653774D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6930357961327123D+0
       V=0.2982170644107595D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3608302115520091D+0
       V=0.2721564237310992D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6498486161496169D+0
       V=0.3033513795811141D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1932945013230339D+0
       V=0.3007949555218533D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3800494919899303D+0
       V=0.2881964603055307D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2899558825499574D+0
       B=0.7934537856582316D+0
       V=0.2958357626535696D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9684121455103957D-1
       B=0.8280801506686862D+0
       V=0.3036020026407088D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1833434647041659D+0
       B=0.9074658265305127D+0
       V=0.2832187403926303D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
     END SUBROUTINE LD0350

!######################################################################

     SUBROUTINE LD0434(X,Y,Z,W,N)
       DOUBLE PRECISION X( 434)
       DOUBLE PRECISION Y( 434)
       DOUBLE PRECISION Z( 434)
       DOUBLE PRECISION W( 434)
       INTEGER N
       DOUBLE PRECISION A,B,V
!VW
!VW    LEBEDEV  434-POINT ANGULAR GRID
!VW
       N=1
       V=0.5265897968224436D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2548219972002607D-2
       Call GEN_OH( 2, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.2512317418927307D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6909346307509111D+0
       V=0.2530403801186355D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1774836054609158D+0
       V=0.2014279020918528D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4914342637784746D+0
       V=0.2501725168402936D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6456664707424256D+0
       V=0.2513267174597564D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2861289010307638D+0
       V=0.2302694782227416D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7568084367178018D-1
       V=0.1462495621594614D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3927259763368002D+0
       V=0.2445373437312980D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.8818132877794288D+0
       V=0.2417442375638981D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9776428111182649D+0
       V=0.1910951282179532D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2054823696403044D+0
       B=0.8689460322872412D+0
       V=0.2416930044324775D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5905157048925271D+0
       B=0.7999278543857286D+0
       V=0.2512236854563495D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5550152361076807D+0
       B=0.7717462626915901D+0
       V=0.2496644054553086D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.9371809858553722D+0
       B=0.3344363145343455D+0
       V=0.2236607760437849D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
     END SUBROUTINE LD0434

!######################################################################

     SUBROUTINE LD0590(X,Y,Z,W,N)
       DOUBLE PRECISION X( 590)
       DOUBLE PRECISION Y( 590)
       DOUBLE PRECISION Z( 590)
       DOUBLE PRECISION W( 590)
       INTEGER N
       DOUBLE PRECISION A,B,V
!VW
!VW    LEBEDEV  590-POINT ANGULAR GRID
!VW
       N=1
       V=0.3095121295306187D-3
       Call GEN_OH( 1, N, X(N), Y(N), Z(N), W(N), A, B, V)
       V=0.1852379698597489D-2
       Call GEN_OH( 3, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.7040954938227469D+0
       V=0.1871790639277744D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6807744066455243D+0
       V=0.1858812585438317D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6372546939258752D+0
       V=0.1852028828296213D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5044419707800358D+0
       V=0.1846715956151242D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4215761784010967D+0
       V=0.1818471778162769D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3317920736472123D+0
       V=0.1749564657281154D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2384736701421887D+0
       V=0.1617210647254411D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1459036449157763D+0
       V=0.1384737234851692D-2
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6095034115507196D-1
       V=0.9764331165051050D-3
       Call GEN_OH( 4, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.6116843442009876D+0
       V=0.1857161196774078D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3964755348199858D+0
       V=0.1705153996395864D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.1724782009907724D+0
       V=0.1300321685886048D-2
       Call GEN_OH( 5, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5610263808622060D+0
       B=0.3518280927733519D+0
       V=0.1842866472905286D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.4742392842551980D+0
       B=0.2634716655937950D+0
       V=0.1802658934377451D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5984126497885380D+0
       B=0.1816640840360209D+0
       V=0.1849830560443660D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.3791035407695563D+0
       B=0.1720795225656878D+0
       V=0.1713904507106709D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.2778673190586244D+0
       B=0.8213021581932511D-1
       V=0.1555213603396808D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       A=0.5033564271075117D+0
       B=0.8999205842074875D-1
       V=0.1802239128008525D-2
       Call GEN_OH( 6, N, X(N), Y(N), Z(N), W(N), A, B, V)
       N=N-1
       RETURN
     END SUBROUTINE LD0590

!######################################################################

   end module lebedev
