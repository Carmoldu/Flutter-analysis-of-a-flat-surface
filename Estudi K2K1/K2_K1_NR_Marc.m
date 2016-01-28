% Code of the Problem 4 - Aerospace vehicles %

close all;
clear all;
clc;

%% Data
for q=1:100

position0(q)=q/100;    
    
k1=1; %Nm/m
k2=k1/(101-q); %Nm/m
rho=1.225; %kg/m3
b=1; %m
 mu=22;   
m=mu*pi*rho*b^2; %kg/m

TOL=1e-6;
tol=[1e-6;1e-6];
alpha=1;

%% Constants

omegah=(k1+k2)/m; %omega_h ^2
omegaalfa=(k2-k1)/m; %omega_alfa ^2



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

K2K1_study_NR0(q)=NRsol(1);

wNR=omegah/sqrt(NRsol(2));
UinfNR=(b*wNR)/NRsol(1); % m/s

end

for q=1:100

position1(q)=q;    
    
k1=1; %Nm/m
k2=q*k1; %Nm/m
rho=1.225; %kg/m3
b=1; %m
 mu=22;   
m=mu*pi*rho*b^2; %kg/m

TOL=1e-6;
tol=[1e-6;1e-6];
alpha=1;

%% Constants

omegah=(k1+k2)/m; %omega_h ^2
omegaalfa=(k2-k1)/m; %omega_alfa ^2



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

K2K1_study_NR1(q)=NRsol(1);

wNR=omegah/sqrt(NRsol(2));
UinfNR=(b*wNR)/NRsol(1); % m/s

end
position0;
position1;
K2K1_study_NR0;
K2K1_study_NR1;
p=fliplr(position0);
A=fliplr(K2K1_study_NR0);

position=[position0 position1];
K2K1_study_NR=[K2K1_study_NR0 K2K1_study_NR1];

[xx,ind] = sort(position);
yy = smooth(position,K2K1_study_NR,0.3,'loess');

figure(1);
semilogx(xx,yy(ind),'b-','linewidth',2)
grid on;
title('Convergència de la k vs K2/K1 per Newton-Raphson');
xlabel('K2/K1');
ylabel('k');