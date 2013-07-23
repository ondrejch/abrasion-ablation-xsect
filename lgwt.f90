! program prglwt
! use types_module
! implicit none
! integer,parameter :: n = 10
! real (dp) :: a = 1,b = 2
! real (dp), dimension(n) :: arrx, arrw
! integer :: i
! call lgwt(arrx,arrw,n,a,b)
! do i = 1,n
!   print *,arrx(i), arrw(i)
! enddo
! end program prglwt

subroutine lgwt(x,w,N,a,b)
! This subroutine computes definite integrals using Legendre-Gauss Quadrature. 
! It computes the Legendre-Gauss nodes and weights  on an interval
! [a,b] with truncation order N 
! Adapted from Greg von Winckel - 02/25/2004
use types_module
implicit none
integer, intent (in) :: N
real (dp),    intent (in) :: a, b
real (dp), dimension (N), intent (out) :: x, w
real (dp), parameter  :: eps = 1e-15
real (dp), parameter  :: pi = atan(1.)*4.
real (dp), dimension (N) :: kx, ky, ky0
! Legendre-Gauss Vandermonde Matrix
real (dp), dimension (N,N+1) :: K
! Derivative of LGVM
real (dp), dimension (N) :: Kp
integer   :: i, j, N1, M1, M2
real (dp) :: tmpx, tmpy
N1=N-1; M1=N; M2=N+1;
K = 0.; Kp = 0.

! matlab: kx=linspace(-1,1,M1)';
tmpx=-1.; tmpy=1.
kx = [(tmpx+ dble(i)*(tmpy-tmpx)/dble(N-1),i=0,N-1)]

! Initial guess
ky = [(cos((2.*dble(i)+1.)*pi/(2.*dble(N1)+2.))+(0.27/dble(M1))*sin(pi*kx(i+1)*dble(N1)/dble(M2)),i=0,n-1)]

! Compute the zeros of the N+1 Legendre Polynomial
! using the recursion relation and the Newton-Raphson method
ky0=2.

! Iterate until new points are uniformly within epsilon of old points
j = 0
do while (maxval(abs(ky-ky0))>eps)
  print *,'Iterations ', j, maxval(abs(ky-ky0))
  j=j+1
  K(:,1)  = 1.
!  Kp(:,1) = 0.
  K(:,2)  = ky
!  Kp(:,2) = 1.

  do i = 2, M1
    K(:,i+1) = ((2.*dble(i)-1.)*ky*K(:,i)-(i-1)*K(:,i-1) )/dble(i)
  enddo

  Kp  = M2*( K(:,M1)-ky*K(:,M2) ) / (1.-ky**2)
  ky0 = ky
  ky  = ky0 - K(:,M2) / Kp
enddo

! Linear map from[-1,1] to [a,b]
x = (a*(1.-ky)+b*(1.+ky))/2.
! Compute the weights
w = (b-a)/((1.-ky**2)*Kp**2)*(dble(M2)/dble(M1))**2

end subroutine lgwt

