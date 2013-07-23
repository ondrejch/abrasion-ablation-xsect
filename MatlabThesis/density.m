%%% This subroutine calculates the nuclear single particle densities
%%% Called bt AbrasionCrs1.m and AbrasionCrs2.m

function[C1,D1]=density(A)
rp= 0.87;
pi=3.1416;

%%% Gauss approximation to Harmonic Well Distriutions
if (A <20)  
%%% Values up to A=16 from Dejager and Devries A from 16-20 expratoplated
    
a = [1.2 1.77 1.51 1.33 1.7 1.75 1.77 1.78 1.79 1.8 1.69 1.649 1.7 1.729 1.8 1.833,1.835, 1.86 1.88 1.89];
g = [0.0 0.0 0.0 0.0 0.1 0.2 0.327 0.45 0.611 0.81 1.0 1.247 1.28 1.291 1.4 1.544 1.61 1.63 1.64 1.65];

s= ((((a(A))^2)/4)-((rp^2)/6))^0.5
D=1/(4*(s^2));
A1=((a(A)^3)/(8*s^3))*(1+(3*g(A)/2)-((3*g(A)*(a(A)^2))/(8*s^2)));
B= ((a(A)^3)/(8*s^3))*((g(A)*(a(A)^2))/(16*(s^4)));
rho0= 1/(((((pi/D)^(3/2)))+((3*B/2)*(((pi^3)/(D^5))^0.5))));

C1= 0.66*(rho0*(a(A))^3/(8*s^3))*(1.5+(3*(g(A))/2)-(3*g(A)*(a(A))^2/(8*s^2))+(g(A)*(a(A))^2*(1)/(16*s^4)));
D1 = 1/(4*(1 +0.52*g(A))*s^2);

%% Gauss approximation to Wood Saxon Distribution 
elseif (A>=20) 
a= (-6e-5*(A)^2)+(0.0342*A)+2.0976;
g=0.000005*(A^2)-0.0017*A+2.5955; 
beta=exp(4.4*rp/(g*sqrt(3)));
ta=(8.8*rp/sqrt(3))*((log((3*beta-1)/(3-beta)))^-1);   
rho02=(1/(((4*pi*(((a^3)/3)+ ((pi^2)*(ta^2)*a/3))))));

    C1=(a/g)^0.3*rho02;
    D1= -0.0109*a+0.1234;
end

C1
D1
