subroutine abrasion(SigAbr, SigAbr2, SigAbr4, Tlab, Aproj, Atarg, Zproj, Ztarg, thetarad)
! This function calculates the Abrasion cross section for the projectile
use binomial_module
use types_module
implicit none
real (dp),  intent (in) :: Tlab, thetarad
integer, intent (in)    :: Aproj, Atarg, Zproj, Ztarg
real (dp),  intent (out):: SigAbr, SigAbr2, SigAbr4
real (dp)               :: SigAbr3, SigAbr5
real (dp),    parameter :: pi = atan(1.)*4.
real (dp),    parameter :: r0 = 1.26         ! Constant for radius calculation (fm)
real (dp),    parameter :: neutronmass = 939.565378   ! Mass of Neutron MeV
real (dp),    parameter :: atomicmass  = 931.494061   ! Atomic mass unit MeV
real (dp),    parameter :: hbarc = 197.327   ! MeV fm ! http://physics.nist.gov/cgi-bin/cuu/Value?hbcmevf
integer, parameter      :: Ninterp = 1000    ! number of Gauss interpolatoion points
real (dp), dimension(Ninterp) :: x1, w1, triarg, xb, wb, argb
real (dp), dimension(Ninterp) :: Xx1, Xx2, Xx3, Xx4, Xx5
real (dp), dimension(Ninterp) :: Pb1, Pb2, Pb3, Pb4, Pb5
real (dp), dimension(Ninterp) :: Zsum1, ZCsum2, ZCsum3, ZCsum4, ZCsum5
real (dp), dimension(Ninterp) :: ZCsumTwo, ZCsumThree, ZCsumFour, ZCsumFive
real (dp), dimension(Aproj) :: Z1, Z2C, Z3C, Z4C, Z5C
integer :: i, j, n, Nproj, Ntarg
real (dp) :: Rp, Rt, Mp, Mt, Cproj, Dproj, Ctarg, Dtarg
real (dp) :: tmpa, tmpb
real (dp) :: Elab, gamma, beta, Plab, PfL, PfT, PK, K, S, B
real (dp) :: Sigpp, Sigpn, Sig, W, V, M1, N1
real (dp) :: Ksum1, Ksum2, Ksum3, Ksum4, Ksum5, binomcoeff

! Calculate other constnats
Nproj = Aproj-Zproj                  ! Number of Neutrons in the projectile            
Ntarg = Atarg-Ztarg                  ! Number of Neutrons in the Target
Rp    = r0*((Aproj)**(1./3.))        ! Projectile radius (fm)
Rt    = r0*((Atarg)**(1./3.))        ! Target Radius (fm)
Mp    = Aproj*atomicmass             ! Projectile total mass
Mt    = Atarg*atomicmass             ! Target total mass

! The values for density parameters taken from Dejager and Devries..
! Call Density function %
! Calcualte the Gauss density constants
call density(Aproj, Cproj, Dproj)
call density(Atarg, Ctarg, Dtarg)

! Call lgwt to do the Gauss Quadrature, Calculate the Gauss points for N=1000
tmpa= -1; tmpb = 1
call lgwt(x1,w1,Ninterp,tmpa,tmpb)
Mp = Aproj*atomicmass                        ! Projectile total mass
Mt = Atarg*atomicmass                        ! Target total mass
!print *, x1(300), w1(300)

! The values for density parameters taken from Dejager and Dev

! Calcualtion of energy and momentum in lab and projectile frame
Elab    = Tlab + neutronmass
gamma   = 1. + Tlab/neutronmass
beta    = sqrt(1.-(1./gamma)**2)
Plab    = sqrt(Elab**2-neutronmass**2)
PfL     = gamma*(Plab*cos(thetarad)-beta*Elab)
PfT     = Plab*sin(thetarad)
PK      = sqrt(PfL**2+PfT**2)
K       = (PK)/(2.*hbarc*sin(thetarad/2.))
!print *,Elab, gamma, beta,Plab,PfL,PfT,PK,K

! the invariant in lab
S        = ((Mp + Mt)**2)+ 2.*Mt*Tlab
! Calculate the Slope Parameter 
B        = 0.0389*(10. + 0.5*log(S/1000000.))
!print *,"S,B ", S,B

! Calculate the nucleon-nucleon parameters
if (Tlab <= 25.) then 
  Sigpp = exp(6.51*(exp(-Tlab/135.)**0.7))/10.
else ! Tlab>25
  Sigpp = ((1.+(5./Tlab))*(40.+(109.*(cos(0.199*((Tlab)**0.5)))*exp(-0.451*((Tlab-25.)**(0.258))))))/10.
endif
Sigpn   = (38.+12500.*exp(-1.187*(((Tlab-0.1)**0.35))))/10.
Sig     = ((dble(Nproj+Ntarg)/dble(Aproj+Atarg))*Sigpn) + Sigpp*(dble(Zproj*Ztarg+Nproj*Ntarg))/dble(Aproj*Atarg)
!print *,"Sigpp, Sigpn,Sig ",Sigpp, Sigpn,Sig

! Calculate the constant values for the potential
W       = Dproj + (1./(2.*B))
V       = Dproj - ((Dproj**2)/W)
M1      = dble(Aproj*Atarg)*Cproj*Ctarg*Sig*((2.*pi*B)**(-1.5))*((pi/W)**1.5)*(((pi/(Dtarg+V))**1.5))
N1      = V - ((V**2)/(Dtarg+V))
!print *,W,V,M1,N1

! Gaussian weights !
!$omp parallel do 
do i = 1, Ninterp
  triarg(i) = (pi/4.)*(x1(i)+1.)
  xb(i)     = tan(triarg(i))
  wb(i)     = (pi/4.)*w1(i)/((cos(triarg(i)))*(cos(triarg(i))))
enddo
!$omp end parallel do

!$omp parallel do 
do n = 1, Ninterp
  argb(n) = xb(n)
  Xx1(n)  = (M1/2.)*((pi/N1)**0.5)*exp(-N1*(argb(n))**2)
  Xx2(n)  = (1./(4.*K))*((M1**2))*((pi/((2.*N1)))**0.5)*((4.*N1*(argb(n))**2)+1.)*exp(-2*N1*(argb(n)**2))
  Xx3(n)  = -1.*(((M1)**3)/(12.*K**2))*((pi/(3.*N1))**0.5)*((argb(n))**2)*(-1.*(36.*(N1)**2*((argb(n))**2)))*exp(-3.*N1*(argb(n))**2)
  Xx4(n)  = (1.*(M1**4)/(48.*K**3))*((pi/(4.*N1))**0.5)*((-24.*N1*((argb(n))**2))-(192.*(N1**2)*(argb(n)**4))+(512.*(N1**3)*(argb(n)**6))-(3.))*exp((-4.*N1*(argb(n))**2))
  Xx5(n)  = -(1.*(M1**5)/(240.*K**4))*((pi/(5.*N1))**0.5)*((8000.*(N1**3)*(argb(n)**6))-(10000.*(N1**4)*(argb(n)**8)))*exp((-5.*N1*(argb(n))**2))

  Pb1(n)  = exp(-2.*Xx1(n)/Aproj)
  Pb2(n)  = exp(-2.*Xx2(n)/Aproj)!(Xx2(n)+Xx1(n))/Aproj)
  Pb3(n)  = exp(-2.*Xx3(n)/Aproj)!(Xx3(n)+Xx2(n)+Xx1(n))/Aproj)
  Pb4(n)  = exp(-2.*Xx4(n)/Aproj)!(Xx4(n)+Xx3(n)+Xx2(n)+Xx1(n))/Aproj)
  Pb5(n)  = exp(-2.*Xx5(n)/Aproj)!(Xx5(n)+Xx4(n)+Xx3(n)+Xx2(n)+Xx1(n))/Aproj)  
!  print *,n, Xx1(n),Xx2(n),Xx3(n),Xx4(n), Pb1(n),Pb2(n),Pb3(n),Pb4(n)
enddo 
!$omp end parallel do

do j = 1, Aproj     
  binomcoeff = dble(binom(Aproj,j))
  !$omp parallel do 
  do n = 1, Ninterp
    Zsum1(n)  = 2.*pi*((argb(n)*wb(n)*((1.-Pb1(n))**j)*((Pb1(n))**(Aproj-j))))
    ZCsum2(n) = 2.*pi*((argb(n)*wb(n)*((1.-Pb2(n))**j)*((Pb2(n))**(Aproj-j))))
    ZCsum3(n) = 2.*pi*((argb(n)*wb(n)*((1.-Pb3(n))**j)*((Pb3(n))**(Aproj-j))))
    ZCsum4(n) = 2.*pi*((argb(n)*wb(n)*((1.-Pb4(n))**j)*((Pb4(n))**(Aproj-j))))
    ZCsum5(n) = 2.*pi*((argb(n)*wb(n)*((1.-Pb5(n))**j)*((Pb5(n))**(Aproj-j))))    
    
    ZCsumTwo(n)   = Zsum1(n)+ZCsum2(n)+ZCsum3(n)
    ZCsumThree(n) = ZCsum3(n)!;%+ZCsum2(n)+ZCsum3(n)
    ZCsumFour(n)  = Zsum1(n)+ZCsum2(n)+ZCsum3(n)+ZCsum4(n)+ZCsum5(n)
    ZCsumFive(n)  = ZCsum5(n)!;%+ZCsum2(n)+ZCsum3(n)+ZCsum4(n)+ZCsum5(n)
  enddo
  !$omp end parallel do
  Z1(j)  = binomcoeff*sum(Zsum1)
  Z2C(j) = binomcoeff*sum(ZCsumTwo)
  Z3C(j) = binomcoeff*sum(ZCsumThree)
  Z4C(j) = binomcoeff*sum(ZCsumFour)
  Z5C(j) = binomcoeff*sum(ZCsumFive)    
enddo

Ksum1 = sum(Z1)
Ksum2 = sum(Z2C)
Ksum3 = sum(Z3C)
Ksum4 = sum(Z4C)
Ksum5 = sum(Z5C)
     
SigAbr  = Ksum1  ! Abrasion with no correction terms
SigAbr2 = Ksum2  ! Abrasion with two correction terms
SigAbr3 = Ksum3
SigAbr4 = Ksum4  ! Abrasion with four correction terms
SigAbr5 = Ksum5

end subroutine abrasion
