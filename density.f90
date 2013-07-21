! program prgdensity
!use types_module
! implicit none
! integer :: A ! atomic number	
! real (dp)    :: C1, D1
! !print *, "Input A "; read *, A 
! do A=1,50
! call density(A, C1, D1)
! print *,'For A = ',A,' constants C1, D1: ', C1, D1
! enddo
! end program prgdensity


subroutine density(A, C1, D1)
! This subroutine calculates the nuclear single particle densities
use types_module
implicit none
integer, intent (in)  :: A
real (dp), intent (out) :: C1, D1
real (dp), parameter  :: pi = atan(1.)*4.
real (dp)::  s, D, A1, B, rho0
real (dp) :: ap, gp, beta, ta, rho02
real (dp) :: rp = 0.87
real (dp), dimension (20) :: atab, gtab 
!integer :: i
!
if (A<20)  then
! Gauss approximation to Harmonic Well Distriutions
! Values up to A=16 from Dejager and Devries A from 16-20 expratoplated
  atab(1:20) =[1.2,1.77,1.51,1.33,1.7,1.75,1.77,1.78,1.79,1.8,1.69,1.649,1.7,1.729,1.8,1.833,1.835,1.86,1.88,1.89]
  gtab(1:20) =[0.0,0.0,0.0,0.0,0.1,0.2,0.327,0.45,0.611,0.81,1.0,1.247,1.28,1.291,1.4,1.544,1.61,1.63,1.64,1.65]
  !do i=1,20; print *,i,atab(i),gtab(i); enddo
  s = ((((atab(A))**2)/4.)-((rp**2)/6.))**0.5
  D = 1./(4.*(s**2));
  A1= ((atab(A)**3)/(8.*s**3))*(1.+(3.*gtab(A)/2)-((3.*gtab(A)*(atab(A)**2))/(8.*s**2)))
  B = ((atab(A)**3)/(8.*s**3))*((gtab(A)*(atab(A)**2))/(16.*(s**4)));
  rho0 = 1./(((((pi/D)**(3./2.)))+((3.*B/2.)*(((pi**3)/(D**5))**0.5))));
  C1 = 0.66*(rho0*(atab(A))**3/(8.*s**3))* &
(1.5+(3.*(gtab(A))/2.)-(3.*gtab(A)*(atab(A))**2/(8.*s**2))+gtab(A)*(atab(A))**2/(16.*s**4)) 
  D1 = 1./(4.*(1. +0.52*gtab(A))*s**2)
else 
! Gauss approximation to Wood Saxon Distribution 
  ap = (-6e-5*(A)**2) +(0.0342*A) +2.0976;
  gp = 0.000005*(A**2)-(0.0017*A) +2.5955; 
  beta = exp(4.4*rp/(gp*sqrt(3.)));
  ta = (8.8*rp/sqrt(3.))*((log((3.*beta-1.)/(3.-beta)))**(-1));   
  rho02 = 1./(((4.*pi*(((ap**3)/3.)+ ((pi**2)*(ta**2)*ap/3.)))));
  C1 = (ap/gp)**0.3 *rho02;
  D1 = -0.0109*ap + 0.1234
endif
!
end subroutine density
!
