!program prgabrasion1
!use types_module
!implicit none
!real (dp), parameter :: pi = atan(1.)*4.
! integer :: Ap1 = 12                ! projectile A
! integer :: At1 = 12                ! target A
! integer :: Zp1 =  6                ! projectile Z
! integer :: Zt1 =  6                ! target Z
! real (dp)        :: Tlab1 = 100.
! real (dp)        :: AA = 10.*pi/180.
! 
! 
! end program prgabrasion1

subroutine abrasioncrs1(SigAbr, SigAbr2, SigAbr4, Tlab1, Ap1, At1, Zp1, Zt1, AA)
! This function calculates the Abrasion cross section for the projectile
use binomial_module
use types_module
implicit none
real (dp),    intent (in) :: Tlab1, AA
integer, intent (in) :: Ap1, At1, Zp1, Zt1
real (dp),    intent (out):: SigAbr, SigAbr2, SigAbr4
real (dp)                 :: SigAbr3, SigAbr5
real (dp),    parameter :: pi = atan(1.)*4.
real (dp),    parameter :: r0 = 1.26         ! Constant for radius calculation (fm)
real (dp),    parameter :: am = 937.57       ! Mass of Neutron MeV
real (dp),    parameter :: amt= (am-7.)      ! Atomic mass MeV
real (dp),    parameter :: hbarc = 197.326   ! MeV fm
integer, parameter :: Ninterp=1000      ! number of Gauss interpolatoion points
real (dp), dimension(Ninterp) :: x1, w1, triarg, xb, wb, argb
real (dp), dimension(Ninterp) :: Xx1, Xx2, Xx3, Xx4, Xx5
real (dp), dimension(Ninterp) :: Pb1, Pb2, Pb3, Pb4, Pb5
real (dp), dimension(Ninterp) :: Zsum1, ZCsum2, ZCsum3, ZCsum4, ZCsum5
real (dp), dimension(Ninterp) :: ZCsumTwo, ZCsumThree, ZCsumFour, ZCsumFive
real (dp), dimension(Ap1) :: Z1, Z2C, Z3C, Z4C, Z5C
integer :: i, j, n, Np, Nt
real (dp) :: Rp, Rt, Mp, Mt, cCp, cDp, cCt, cDt
real (dp) :: a1, b1
real (dp) :: Elab, Gamma1, Beta1, Plab1, Pfl1, PfT1, PK, K, S, B
real (dp) :: Sigpp, Sigpn, Sig, W, V, M1, N1
real (dp) :: Ksum1, Ksum2, Ksum3, Ksum4, Ksum5, binomcoeff

! Calculate other constnats
Np = Ap1-Zp1                    ! Number of Neutrons in the projectile            
Nt = At1-Zt1                    ! Number of Neutrons in the Target
Rp = r0*((Ap1)**(1./3.))        ! Projectile radius (fm)
Rt = r0*((At1)**(1./3.))        ! Target Radius (fm)
Mp = Ap1*amt                    ! Projectile total mass
Mt = At1*amt                    ! Target total mass

! The values for density parameters taken from Dejager and Devries..
! Call Density function %
! Calcualte the Gauss density constants
call density(Ap1, cCp, cDp)
call density(At1, cCt, cDt)

! Call lgwt to do the Gauss Quadrature, Calculate the Gauss points for N=1000
a1= -1; b1 = 1
call lgwt(x1,w1,Ninterp,a1,b1)
Mp = Ap1*amt                        ! Projectile total mass
Mt = At1*amt                        ! Target total mass
!print *, x1(300), w1(300)

! The values for density parameters taken from Dejager and Dev

! Calcualtion of energy and momentum in lab and projectile frame
Elab     = Tlab1+am
Gamma1   = 1.+Tlab1/am
Beta1    = sqrt(1.-(1./Gamma1)**2)
Plab1    = sqrt(Elab**2-am**2)
Pfl1     = Gamma1*(Plab1*cos(AA)-Beta1*Elab)
PfT1     = Plab1*sin(AA)
PK       = sqrt(Pfl1**2+PfT1**2)
K        = (PK)/(2.*hbarc*sin(AA/2.))
!print *,Elab, Gamma1, Beta1,Plab1,Pfl1,PfT1,PK,K

! the invariant in lab
S        = ((Mp + Mt)**2)+ 2.*Mt*Tlab1
! Calculate the Slope Parameter 
B        = 0.0389*(10. + 0.5*log(S/1000000.))
!print *,"S,B ", S,B

! Calculate the nucleon-nucleon parameters
if (Tlab1<=25) then 
  Sigpp = exp(6.51*(exp(-Tlab1/135.)**0.7))/10.
else ! Tlab1>25
  Sigpp = ((1.+(5./Tlab1))*(40.+(109.*(cos(0.199*((Tlab1)**0.5)))*exp(-0.451*((Tlab1-25.)**(0.258))))))/10.
endif
Sigpn   = (38.+12500.*exp(-1.187*(((Tlab1-0.1)**0.35))))/10.
Sig     = ((dble(Np+Nt)/dble(Ap1+At1))*Sigpn) + Sigpp*(dble(Zp1*Zt1+Np*Nt))/dble(Ap1*At1)
!print *,"Sigpp, Sigpn,Sig ",Sigpp, Sigpn,Sig

! Calculate the constant values for the potential
W       = cDp + (1./(2.*B))
V       = cDp - ((cDp**2)/W)
M1      = dble(Ap1*At1)*cCp*cCt*Sig*((2.*pi*B)**(-1.5))*((pi/W)**1.5)*(((pi/(cDt+V))**1.5))
N1      = V - ((V**2)/(cDt+V))
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

  Pb1(n)  = exp(-2.*Xx1(n)/Ap1)
  Pb2(n)  = exp(-2.*Xx2(n)/Ap1)!(Xx2(n)+Xx1(n))/Ap1)
  Pb3(n)  = exp(-2.*Xx3(n)/Ap1)!(Xx3(n)+Xx2(n)+Xx1(n))/Ap1)
  Pb4(n)  = exp(-2.*Xx4(n)/Ap1)!(Xx4(n)+Xx3(n)+Xx2(n)+Xx1(n))/Ap1)
  Pb5(n)  = exp(-2.*Xx5(n)/Ap1)!(Xx5(n)+Xx4(n)+Xx3(n)+Xx2(n)+Xx1(n))/Ap1)  
!  print *,n, Xx1(n),Xx2(n),Xx3(n),Xx4(n), Pb1(n),Pb2(n),Pb3(n),Pb4(n)
enddo 
!$omp end parallel do

do j = 1, Ap1     
  binomcoeff = dble(binom(Ap1,j))
  !$omp parallel do 
  do n = 1, Ninterp
    Zsum1(n)  = 2.*pi*((argb(n)*wb(n)*((1.-Pb1(n))**j)*((Pb1(n))**(Ap1-j))))
    ZCsum2(n) = 2.*pi*((argb(n)*wb(n)*((1.-Pb2(n))**j)*((Pb2(n))**(Ap1-j))))
    ZCsum3(n) = 2.*pi*((argb(n)*wb(n)*((1.-Pb3(n))**j)*((Pb3(n))**(Ap1-j))))
    ZCsum4(n) = 2.*pi*((argb(n)*wb(n)*((1.-Pb4(n))**j)*((Pb4(n))**(Ap1-j))))
    ZCsum5(n) = 2.*pi*((argb(n)*wb(n)*((1.-Pb5(n))**j)*((Pb5(n))**(Ap1-j))))    
    
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

end subroutine abrasioncrs1
