%%%% This function calculates the Abrasion cross section for the projectile%%%%
%%%% Called by MainRun.m
function[SigAbr, SigAbr2, SigAbr4]=AbrasionCrs1(Tlab1,Ap1,At1,Zp1,Zt1,AA)
%%%%% Calculate other constants%%%%%%
Np= Ap1-Zp1; %% Number of Neutrons in the projectile
Nt= At1-Zt1; %% Number of Neutrons in the Target
r0= 1.26; %% Constant for radius calculation (fm)
Rp=r0*((Ap1)^(1/3)); %% Projectile radius (fm)
Rt= r0*((At1)^(1/3)); %% Target Radius (fm)
am=937.57; %% Mass of Neutron MeV
amt= (am-7); %% Atomic mass MeV
Mp=Ap1*amt; %% Projectile total mass
Mt=At1*amt; %% Target total mass
hbarc= 197.326; %% MeV fm

%%% The values for density parameters taken from DeJager and DeVries.
%%% Call Density function %%%%%
%%% Calculate the Gauss density constants
Ap1
[Cp,Dp]=density(Ap1);
At1
[Ct,Dt]=density(At1);

%%% Call lgwt to do the Gauss Quadrature%%%
% Calculate the Gauss points for N=1000
Nint=1000;
a1= -1;
b1 = 1;
[x1,w1]=lgwt(Nint,a1,b1);

% Calculation of the energy and momentum in lab and projectile frame
Elab=Tlab1+am;
Gamma1= 1+Tlab1/am;
Beta1=sqrt(1-(1/Gamma1)^2);
Plab1=sqrt(Elab^2-am^2);
Pfl1=Gamma1*(Plab1*cos(AA)-Beta1*Elab);
PfT1=Plab1*sin(AA);
PK=sqrt(Pfl1^2+PfT1^2);
K=(PK)/(2*hbarc*sin(AA/2));
%%% The invariant in lab
S= ((Mp + Mt)^2)+ 2*Mt*Tlab1;
%% Calculate the Slope Parameter %%%
B= 0.0389*(10+ 0.5*log((((S)/1000000))));

%%% Calculate the nucleon-nucleon parameters
if Tlab1<=25
  Sigpp=exp(6.51*(exp(-Tlab1/135)^0.7))/10;
elseif Tlab1>25
  Sigpp=((1+(5/Tlab1))*(40+(109*(cos(0.199*((Tlab1)^0.5)))*exp(-0.451*((Tlab1-25)^(0.258))))))/10;
end
Sigpn=(38+12500*exp(-1.187*(((Tlab1-0.1)^0.35))))/10;
Sig=(((Np+Nt)/(Ap1+At1))*(Sigpn))+Sigpp*((Zp1*Zt1+Np*Nt))/(Ap1*At1);
Sigpp
Sigpn
Sig

%%% Calculate the constant values for the potential
W= Dp+(1/(2*B));
V= Dp-((Dp^2)/(W));
M1= Ap1*At1*Cp*Ct*Sig*((2*pi*B)^(-1.5))*((pi/W)^1.5)*(((pi/(Dt+V))^1.5));
N1= V-(((V)^2)/(Dt+V));
W
V
M1
N1

%%%% Gaussian Weights %%%
for i=1:Nint
  triarg(i)=(pi/4)*(x1(i)+1);
  xb(i)=tan(triarg(i));
  wb(i)=((pi/4)*w1(i))/((cos(triarg(i)))*(cos(triarg(i))));
end
for k=1:Ap1 ;
  for n=1:Nint
    argb(n)= xb(n);
    X1(n)= ((M1)/2)*((pi/N1)^0.5)*exp(-N1*(argb(n))^2);
    X2(n)=(1/(4*K))*((M1^2))*((pi/((2*N1)))^0.5)*((4*N1*(argb(n))^2)+1)*exp(-2*N1*(argb(n)^2));
    X3(n)= -1*(((M1)^3)/(12*K^2))*((pi/(3*N1))^0.5)*((argb(n))^2)*(-1*(36*(N1)^2*((argb(n))^2)))*exp((-3*N1*(argb(n))^2));
    X4(n)=(1*(M1^4)/(48*K^3))*((pi/(4*N1))^0.5)*((-24*N1*((argb(n))^2))-(192*(N1^2)*(argb(n)^4))+(512*(N1^3)*(argb(n)^6))-(3))*exp((-4*N1*(argb(n))^2));
    X5(n)=-(1*(M1^5)/(240*K^4))*((pi/(5*N1))^0.5)*((8000*(N1^3)*(argb(n)^6))-(10000*(N1^4)*(argb(n)^8)))*exp((-5*N1*(argb(n))^2));

    Pb1(n)=exp(-2*X1(n)/Ap1);
    Pb2(n)=exp(-2*X2(n)/Ap1);
    Pb3(n)=exp(-2*X3(n)/Ap1);
    Pb4(n)=exp(-2*X4(n)/Ap1);
    Pb5(n)=exp(-2*X5(n)/Ap1);
    Zsum1(n)= 2*pi*((argb(n)*wb(n)*((1-Pb1(n))^k)*((Pb1(n))^(Ap1-k))));
    ZCsum2(n)=2*pi*((argb(n)*wb(n)*((1-Pb2(n))^(k))*((Pb2(n))^((Ap1-k)))));
    ZCsum3(n)=2*pi*((argb(n)*wb(n)*((1-Pb3(n))^(k))*((Pb3(n))^((Ap1-k)))));
    ZCsum4(n)=2*pi*((argb(n)*wb(n)*((1-Pb4(n))^(k))*((Pb4(n))^((Ap1-k)))));
    ZCsum5(n)=2*pi*((argb(n)*wb(n)*((1-Pb5(n))^(k))*((Pb5(n))^((Ap1-k)))));
    ZCsumTwo(n)=Zsum1(n)+ZCsum2(n)+ZCsum3(n);
    ZCsumFour(n)=Zsum1(n)+ZCsum2(n)+ZCsum3(n)+ZCsum4(n)+ZCsum5(n);
  end
  Z1(k)=nchoosek(Ap1,k)*sum(Zsum1);
  Z2C(k)=nchoosek(Ap1,k)*sum(ZCsumTwo);
  Z4C(k)=nchoosek(Ap1,k)*sum(ZCsumFour);
end

Ksum1=sum(Z1);
Ksum2=sum(Z2C);
Ksum4=sum(Z4C);
SigAbr=Ksum1; % Abrasion with no correction terms
SigAbr2=Ksum2; % Abrasion with two correction terms
SigAbr4=Ksum4; % Abrasion with four correction terms
