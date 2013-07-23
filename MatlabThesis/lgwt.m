
% This script is computes definite integrals using Legendre-Gauss 
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N
% Adapted from Greg von Winckel - 02/25/2004
function [x,w]=lgwt(N,a,b)
N=N-1;
M1=N+1; M2=N+2;
kx=linspace(-1,1,M1)';
% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/M1)*sin(pi*kx*N/M2);
% Legendre-Gauss Vandermonde Matrix
K=zeros(M1,M2);
% Derivative of LGVM
Kp=zeros(M1,M2);
% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method
y0=2;
% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps
    
    
    K(:,1)=1;
    Kp(:,1)=0;
    
    K(:,2)=y;
    Kp(:,2)=1; 
    
    for k=2:M1
        K(:,k+1)=( (2*k-1)*y.*K(:,k)-(k-1)*K(:,k-1) )/k;
    end
 
    Kp=(M2)*( K(:,M1)-y.*K(:,M2) )./(1-y.^2);    % ondrejch-- here Kp turns into 1D vector. Eh?
    
    y0=y;
    y=y0-K(:,M2)./Kp;
    
end

% Linear map from[-1,1] to [a,b]
x=(a*(1-y)+b*(1+y))/2;      

% Compute the weights
w=(b-a)./((1-y.^2).*Kp.^2)*(M2/M1)^2;

% x

% w
