program xsect
! Calculation of Lorentz Invariant and Double differential 
! cross sections from abrasion-ablation model
use types_module
implicit none
real (dp), parameter :: pi = atan(1.)*4.
real (dp), parameter :: neutronmass =  939.565378  ! Mass of Neutron  MeV
real (dp), parameter :: atomicmass  =  931.494061  ! Atomic mass unit MeV
integer, parameter :: N = 50                 ! steps in energy
integer, parameter :: maxBinomA = 1750       ! largest A for which we can calculate binomial coefficient
real (dp), dimension(N) :: Ek, W, TL, EF, Plab, Pf, PFL, PFT, PFF, EK1, EZK
real (dp), dimension(N) :: NuMom, NuMom1, NuMom2, NuMom3
real (dp), dimension(N) :: NuMomc2, NuMom1c2, NuMom2c2, NuMom3c2
real (dp), dimension(N) :: NuMomc4, NuMom1c4, NuMom2c4, NuMom3c4
real (dp), dimension(N) :: NuMomt, NuMom1t, NuMom2t, NuMom3t
real (dp), dimension(N) :: NuMomtc2, NuMom1tc2, NuMom2tc2, NuMom3tc2
real (dp), dimension(N) :: NuMomtc4, NuMom1tc4, NuMom2tc4, NuMom3tc4
real (dp), dimension(N) :: LnInv, LnInv1, LnInv2
real (dp), dimension(N) :: LnInvt, LnInv1t, LnInv2t
real (dp), dimension(N) :: DbDiff, DbDiff1, DbDiff2
real (dp), dimension(N) :: DbDifft, DbDiff1t, DbDiff2t
real (dp) :: H, Emin, Emax, PNN, T0lab, E0lab, P0lab, Gamma, Beta
real (dp) :: SigA, SigA1, SigA2, SigAt, SigA1t, SigA2t
real (dp) :: Ebeam, Pbeam, Pn, Tn, Gamman, Betan
real (dp) :: KF, P1, P2, P3, C1, C2, C3, N0
integer   :: i

real (dp) :: Tlab       ! projectile lab energy [MeV]
real (dp) :: amp        ! total projectile mass MeV 
integer   :: Ap         ! projectile A
integer   :: At         ! target A
integer   :: Zp         ! projectile Z
integer   :: Zt         ! target Z
real (dp) :: thetadeg   ! scattering thetarad [deg]
real (dp) :: thetarad   ! scattering thetarad [radians]

print *,'Enter projectile energy in lab [MeV] '; read *,Tlab
print *,'Enter projectile mass number ';         read *,Ap
print *,'Enter projectile proton number ';       read *,Zp
print *,'Enter target mass number ';             read *,At
print *,'Enter target proton number ';           read *,Zt
print *,'Enter scattering angle [deg] ';         read *,thetadeg

! sanity checks
if (Ap<Zp) then
  print *,"Ap cannot be smaller than Zp" ; stop
endif
if (At<Zt) then
  print *,"At cannot be smaller than Zt" ; stop
endif
if ((Ap<1) .or. (At<1) .or. (Zp<1) .or. (Zt<1)) then
  print *,"Atomic and proton numbers have to be greater than 1"; stop
endif
if ((Ap>maxBinomA) .or. (At>maxBinomA)) then
  print *,"Maximum atomic number we support is ", maxBinomA; stop
endif
if ((thetadeg<0.) .or. (thetadeg>180.)) then
  print *,"Scattering angle needs to be between 0 and 180"; stop
endif

amp = dble(Ap)*atomicmass         ! total projectile mass
thetarad  = (thetadeg*pi/180.)    ! degrees -> radians
!print *,"thetarad ", thetarad

! Calcualte the energy grid for outgoing Nucleon energy
Emin = 5.00             ! MeV
Emax = 3.*Tlab          ! MeV

Ek(1)  = Emin
Ek(N)  = Emax
H      = (Emax-Emin)/dble(N-1)
W(1)   = H/2.

do i = 2, N-1
  Ek(i) = Ek(i-1) + H
  W(i)  = H
enddo

! Incoming Beam Calculations
PNN   = sqrt(((Tlab+neutronmass)**2) - neutronmass**2)
T0lab = Tlab*dble(Ap)
E0lab = T0lab + amp
P0lab = sqrt((E0lab**2) - amp**2)
Gamma = 1. + Tlab/neutronmass
Beta  = sqrt(1. - (1./Gamma)**2)
!print *,"PNN, T0lab, E0lab, P0lab, Gamma, Beta ", PNN, T0lab, E0lab, P0lab, Gamma, Beta

! incoming nucleon
Ebeam = E0lab/dble(Ap)
if(neutronmass>Ebeam) then
  print *, ' *** Negative beam momentum!';  stop
endif
Pbeam = sqrt(Ebeam**2 - neutronmass**2)
Pn    = Pbeam
Tn    = sqrt(neutronmass**2 + Pn**2) - neutronmass
Gamman= 1. + Tn/neutronmass
Betan = sqrt(1.-(1./Gamman)**2)
!print *,"Ap, Ebeam, Pbeam, Gamman, Betan ",Ap, Ebeam, Pbeam, Gamman, Betan

! Calculation for nucleon momentum in projectile rest frame
!$omp parallel do 
do i = 1,N
  TL(i)  = Ek(i)
  EF(i)  = TL(i) + neutronmass                      ! Lab Frame Nucleon energy
  Plab(i)= sqrt((EF(i))**2-neutronmass**2)          ! Lab Frame Momentum
  Pf(i)  = Plab(i)
  PFL(i) = Gamman*(Pf(i)-(Betan*EF(i)))
  PFT(i) = Pf(i)*sin(thetarad)
  PFF(i) = sqrt((PFL(i))**2 + (PFT(i))**2) ! Projectile rest frame 
  EK1(i) = sqrt((PFF(i))**2 + neutronmass**2)
  EZK(i) = sqrt((Plab(i))**2 + neutronmass**2)  
!print *,i,EF(i),EK1(i),EZK(i)
enddo
!$omp end parallel do  

! Calculte the Fermi Momentum and Momentum Distributions 
 KF = 26.*log(dble(Ap)) + 129.
 P1 = KF*sqrt(2./5.)
 C1 = 1.
 P2 = KF*sqrt(6./5.)
 C2 = 0.03
 P3 = 500.
 C3 = 0.0008
 
! Normalization Constant
N0 = 1./((C1*(2.*pi*P1**2)**1.5)+(C2*(2.*pi*P2**2)**1.5)+((C3*(2.*pi*P3**2)**1.5)))

! Call the function get the Abrasion cross section
call abrasion(SigA, SigA1, SigA2, Tlab, Ap, At, Zp, Zt, thetarad)
print *," * abr1 [SigA, SigA1, SigA2] = ", SigA, SigA1, SigA2
call abrasion(SigAt,SigA1t,SigA2t,Tlab, At, Ap, Zt, Zp, thetarad)
print *," * abr2 [SigAt, SigA1t, SigA2t] = ", SigAt, SigA1t, SigA2t

! Calulation of Lorentz Invariant and Double Differential %%
!$omp parallel do 
do i = 1,N   
  ! Projectile,  no correction terms
  NuMom1(i) = SigA*N0*C1*exp(-((PFF(i))**2)/(2.*P1**2))
  NuMom2(i) = SigA*N0*C2*exp(-((PFF(i))**2)/(2.*P2**2))
  NuMom3(i) = SigA*N0*C3*exp(-((PFF(i))**2)/(2.*P3**2))
  NuMom(i)  = NuMom1(i) + NuMom2(i) + NuMom3(i)

  ! Projectile, two correction terms 
  NuMom1c2(i) = SigA1*N0*C1*exp(-((PFF(i))**2)/(2.*P1**2))
  NuMom2c2(i) = SigA1*N0*C2*exp(-((PFF(i))**2)/(2.*P2**2))
  NuMom3c2(i) = SigA1*N0*C3*exp(-((PFF(i))**2)/(2.*P3**2))
  NuMomc2(i)  = NuMom1c2(i) + NuMom2c2(i) + NuMom3c2(i)

  ! Projectile, four corrrection terms
  NuMom1c4(i) = SigA2*N0*C1*exp(-((PFF(i))**2)/(2.*P1**2))
  NuMom2c4(i) = SigA2*N0*C2*exp(-((PFF(i))**2)/(2.*P2**2))
  NuMom3c4(i) = SigA2*N0*C3*exp(-((PFF(i))**2)/(2.*P3**2))
  NuMomc4(i)  = NuMom1c4(i) + NuMom2c4(i) + NuMom3c4(i)

  ! Target Contribution, no correction terms
  NuMom1t(i) = SigAt*N0*C1*exp(-((Plab(i))**2)/(2.*P1**2))
  NuMom2t(i) = SigAt*N0*C2*exp(-((Plab(i))**2)/(2.*P2**2))
  NuMom3t(i) = SigAt*N0*C3*exp(-((Plab(i))**2)/(2.*P3**2))
  NuMomt(i)  = NuMom1t(i) + NuMom2t(i) + NuMom3t(i)

  ! Target, two correction terms
  NuMom1tc2(i) = SigA1t*N0*C1*exp(-((Plab(i))**2)/(2.*P1**2))
  NuMom2tc2(i) = SigA1t*N0*C2*exp(-((Plab(i))**2)/(2.*P2**2))
  NuMom3tc2(i) = SigA1t*N0*C3*exp(-((Plab(i))**2)/(2.*P3**2))
  NuMomtc2(i)  = NuMom1tc2(i) + NuMom2tc2(i) + NuMom3tc2(i)

  ! Target, four correction terms
  NuMom1tc4(i) = SigA2t*N0*C1*exp(-((Plab(i))**2)/(2.*P1**2))
  NuMom2tc4(i) = SigA2t*N0*C2*exp(-((Plab(i))**2)/(2.*P2**2))
  NuMom3tc4(i) = SigA2t*N0*C3*exp(-((Plab(i))**2)/(2.*P3**2))
  NuMomtc4(i)  = NuMom1tc4(i) + NuMom2tc4(i) + NuMom3tc4(i)

  ! Lorentz Invariant for 0,2 & 4 correction terms 
  LnInv(i)   = 10.*(EK1(i)*NuMom(i)  +EZK(i)*NuMomt(i))
  LnInv1(i)  = 10.*(EK1(i)*NuMomc2(i)+EZK(i)*NuMomtc2(i))
  LnInv2(i)  = 10.*(EK1(i)*NuMomc4(i)+EZK(i)*NuMomtc4(i))

  LnInvt(i)  = 10.*(EZK(i)*NuMomt(i))
  LnInv1t(i) = 10.*(EZK(i)*NuMomtc2(i))
  LnInv2t(i) = 10.*(EZK(i)*NuMomtc4(i))

  ! Double Differential for 0,2 and 4 correction terms.
  DbDiff(i)  = Plab(i)*LnInv(i)
  DbDiff1(i) = Plab(i)*LnInv1(i)
  DbDiff2(i) = Plab(i)*LnInv2(i)

  DbDifft(i) = Plab(i)*LnInvt(i)
  DbDiff1t(i)= Plab(i)*LnInv1t(i)
  DbDiff2t(i)= Plab(i)*LnInv2t(i)
enddo
!$omp end parallel do 

! print output
do i = 1,N   
  print *, i, DbDiff(i), DbDiff1(i), DbDiff2(i), DbDifft(i), DbDiff1t(i), DbDiff2t(i)
enddo

end program xsect
