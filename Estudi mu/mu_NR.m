% Code of the Problem 4 - Aerospace vehicles %

close all;
clear all;
clc;

%% Data

k1=1; %Nm/m
k2=k1; %Nm/m
rho=1.225; %kg/m3
b=1; %m
for q=1:100

position(q)=q;
    
m=q*pi*rho*b^2; %kg/m

TOL=1e-6;
tol=[1e-6;1e-6];
alpha=1;

%% Constants

omegah=(k1+k2)/m; %omega_h ^2
omegaalfa=(k2-k1)/m; %omega_alfa ^2

mu=q;

Ck=@(k) besselh(1,2,k)/(besselh(1,2,k)+besselh(0,2,k)*i);

%% Development

% Coefficients

% Lift
Clh=@(k) 2*pi*(2*k*Ck(k)*i-k^2);
Clalfa=@(k) 2*pi*(2*Ck(k)+k*(1+Ck(k))*i);

% Momentum
Cmah=@(k) 2*pi*k*Ck(k)*i;
Cmaalfa=@(k) pi*(0.25*k^2+2*Ck(k)+k*(Ck(k)-1)*i);

% Matrices
Q=@(k) [-Clh(k) -Clalfa(k);Cmah(k) Cmaalfa(k)];
M=[1 0;0 1/3];
K=[1 (omegaalfa/omegah);(omegaalfa/omegah) 1]; 

A=@(k) K\(M+(1/(2*pi*mu*k^2))*Q(k));
B=@(k,tau) A(k)-tau*eye(2);
d=@(k,tau) det(B(k,tau));

G=@(k,tau) -real(d(k,tau));
H=@(k,tau) imag(d(k,tau));

fun=@(x) [G(x(1),x(2)),H(x(1),x(2))];

%Initial values

tau0=0.1;
k0=0.1;
x0=[k0;tau0];

%% Newton-Raphson
options=optimset('Display','iter','TolFun',TOL);
tic;
NRsol=fsolve(fun,x0,options);
toc;
t_NR=toc;

mu_study_NR(q)=NRsol(1);

wNR=omegah/sqrt(NRsol(2));
UinfNR=(b*wNR)/NRsol(1); % m/s

end

[xx,ind] = sort(position);
yy = smooth(position,mu_study_NR,0.4,'loess');

figure(1);
plot(xx,yy(ind),'b-','linewidth',2)
grid on;
title('Convergència de la k vs mu per Newton-Raphson');
xlabel('mu');
ylabel('k');