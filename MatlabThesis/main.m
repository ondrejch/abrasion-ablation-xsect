% The program calculates the Lorentz Invariant and the Double
%%%%% differential cross sections from the abrasion process.
%%%%% The cross sections from ablation are calculated using the abrasion
%%%%% cross section generated by this program and the probability functions
%%%%% from the UBERNSPEC code.
clc;
clear;
%%% Input properties%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
Tlab= input('\n Please Enter Projectile Lab. Energy : ');
Ap= input('\n Please Enter Projectile Mass Number : ');
Zp= input('\n Please Enter Projectile Charge Number : ');
At= input('\n Please Enter Target Mass Number : ');
Zt= input('\n Please Enter the Target Charge Number : ');
theta= input('\n Please Enter the Scattering Angle (Degrees) : ');
hbar=197.326;
%% Change angle to Radians%%%%
AA=(theta*pi/180);
am=937.57; %%%mass of nucleon
amp= Ap*(am-7); %% total projectile mass
%%% Calculate the energy grid for outgoing Nucleon energy%%%
Emin=5.00; %% MeV
Emax=3*Tlab; %% MeV
Ek(1)=Emin;
Ek(50)=Emax;
N=50;
NM1=N-1;
H=(Emax-Emin)/NM1;
W(1)=H/2;
for i=2:49
  Ek(i)=Ek(i-1)+H;
  W(i)=H;
end

%%%%%%Incoming Beam Calculations%%%%%%%%%%%%
PNN = sqrt(((Tlab+am)^2)-am^2);
T0lab=Tlab*Ap;
E0lab=(T0lab+amp);
P0lab=sqrt((E0lab^2)-amp^2);
Gamma = 1+ Tlab/am;
Beta= sqrt(1-(1/Gamma)^2);
%% Incoming nucleon %%%
Ebeam=E0lab/Ap;
Pbeam=sqrt((Ebeam^2)-am^2);
Pn=Pbeam;
Tn=sqrt(am^2+Pn^2)-am;
Gamman= 1+Tn/am;
Betan= sqrt(1-(1/Gamman)^2);
%%%% Calculation for nucleon momentum in projectile rest frame%%%%
for n=1:50
  TL(n)=Ek(n);
  EF(n)=TL(n)+am; %%%% Lab Frame Nucleon energy
  Plab(n)=sqrt((EF(n))^2-am^2); %%%% Lab Frame Momentum
  Pf(n)=Plab(n);
  PFL(n)=Gamman*(Pf(n)-(Betan*EF(n)));
  PFT(n)=Pf(n)*sin(AA);
  PFF(n)=sqrt(((PFL(n))^2)+((PFT(n))^2)); %%% Projectile rest frame
  EK1(n)=sqrt(((PFF(n))^2)+am^2);
  EZK(n)=sqrt((Plab(n))^2+am^2);
end
%% Calculate the Fermi Momentum and Momentum Distributions
KF=26*log(Ap)+129;
P1=KF*sqrt(2/5);
C1=1;
P2=KF*sqrt(6/5);
C2=0.03;
P3=500;
C3=0.0008;
% Normalization Constant
N0=1/((C1*(2*pi*P1^2)^1.5)+(C2*(2*pi*P2^2)^1.5)+((C3*(2*pi*P3^2)^1.5)));

%%% Call the functions get the Abrasion cross section for the projectile
%%%%% and the target
[SigA,SigA1,SigA2]   =AbrasionCrs1(Tlab,Ap,At,Zp,Zt,AA) %projectile
[SigAt,SigA1t,SigA2t]=AbrasionCrs2(Tlab,At,Ap,Zt,Zp,AA) %Target
%%%%% Calculation of Lorentz Invariant and Double Differential %%
for k=1:50
%Projectile, no correction terms
NuMom1(k)=(SigA)*N0*C1*exp(-((PFF(k))^2)/(2*P1^2));
NuMom2(k)=(SigA)*N0*C2*exp(-((PFF(k))^2)/(2*P2^2));
NuMom3(k)=(SigA)*N0*C3*exp(-((PFF(k))^2)/(2*P3^2));
NuMom(k)=NuMom1(k)+NuMom2(k)+NuMom3(k);
%Projectile, two correction terms
NuMom1c2(k)=(SigA1)*N0*C1*exp(-((PFF(k))^2)/(2*P1^2));
NuMom2c2(k)=(SigA1)*N0*C2*exp(-((PFF(k))^2)/(2*P2^2));
NuMom3c2(k)=(SigA1)*N0*C3*exp(-((PFF(k))^2)/(2*P3^2));
NuMomc2(k)=NuMom1c2(k)+NuMom2c2(k)+NuMom3c2(k);
%Projectile, four correction terms
NuMom1c4(k)=(SigA2)*N0*C1*exp(-((PFF(k))^2)/(2*P1^2));
NuMom2c4(k)=(SigA2)*N0*C2*exp(-((PFF(k))^2)/(2*P2^2));
NuMom3c4(k)=(SigA2)*N0*C3*exp(-((PFF(k))^2)/(2*P3^2));
NuMomc4(k)=NuMom1c4(k)+NuMom2c4(k)+NuMom3c4(k);
% Target Contribution, no correction terms
NuMom1t(k)=(SigAt)*N0*C1*exp(-((Plab(k))^2)/(2*P1^2));
NuMom2t(k)=(SigAt)*N0*C2*exp(-((Plab(k))^2)/(2*P2^2));
NuMom3t(k)=(SigAt)*N0*C3*exp(-((Plab(k))^2)/(2*P3^2));
NuMomt(k)=NuMom1t(k)+NuMom2t(k)+NuMom3t(k);
% Target, two correction terms
NuMom1tc2(k)=(SigA1t)*N0*C1*exp(-((Plab(k))^2)/(2*P1^2));
NuMom2tc2(k)=(SigA1t)*N0*C2*exp(-((Plab(k))^2)/(2*P2^2));
NuMom3tc2(k)=(SigA1t)*N0*C3*exp(-((Plab(k))^2)/(2*P3^2));
NuMomtc2(k)=NuMom1tc2(k)+NuMom2tc2(k)+NuMom3tc2(k);
% Target, four correction terms
NuMom1tc4(k)=(SigA2t)*N0*C1*exp(-((Plab(k))^2)/(2*P1^2));
NuMom2tc4(k)=(SigA2t)*N0*C2*exp(-((Plab(k))^2)/(2*P2^2));
NuMom3tc4(k)=(SigA2t)*N0*C3*exp(-((Plab(k))^2)/(2*P3^2));
NuMomtc4(k)=NuMom1tc4(k)+NuMom2tc4(k)+NuMom3tc4(k);

%% Lorentz Invariant for 0,2 & 4 correction terms
LnInv(k)= 10*(EK1(k)*NuMom(k)+EZK(k)*NuMomt(k));
LnInv1(k)= 10*(EK1(k)*NuMomc2(k)+EZK(k)*NuMomtc2(k));
%LnInv2(k)= 10*(EK1(k)*NuMomc4(k)+EZK(k)*NuMomtc4(k));
LnInvt(k)=10*(EZK(k)*NuMomt(k));
LnInv1t(k)=10*(EZK(k)*NuMomtc2(k));
%LnInv2t(k)=10*(EZK(k)*NuMomtc4(k));
%Double Differential for 0,2 and 4 correction terms.
DbDiff(k)=Plab(k)*LnInv(k);
DbDiff1(k)=Plab(k)*LnInv1(k);
%DbDiff2(k)=Plab(k)*LnInv2(k);
DbDifft(k)=Plab(k)*LnInvt(k);
DbDiff1t(k)=Plab(k)*LnInv1t(k);
%DbDiff2t(k)=Plab(k)*LnInv2t(k);
end
%%Display Results and plot
DbDiff'
DbDiff1'
%plot(Ek, DbDiff, Ek, DbDiff1)
