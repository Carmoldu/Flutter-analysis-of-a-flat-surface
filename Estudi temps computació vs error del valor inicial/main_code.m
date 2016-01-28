% Code of the Problem 4 - Aerospace vehicles %

close all;
clear all;
clc;

%% Data

k1=1; %Nm/m
k2=k1; %Nm/m
rho=1.225; %kg/m3
b=1; %m
m=22*pi*rho*b^2; %kg/m

TOL=1e-6;
tol=[1e-6;1e-6];
alpha=1;

%% Constants

omegah=(k1+k2)/m; %omega_h ^2
omegaalfa=(k2-k1)/m; %omega_alfa ^2

mu=m/(pi*rho*b^2);

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

%Convergència del mètode: Newton-Raphson

E_NR=1e3; % Per entrar al while
iterationsNR=1;

k_NR=NRsol(1);

while TOL<E_NR
    
    options_2=struct('MaxIter',iterationsNR);
    NRsol=fsolve(fun,x0,options_2);
    E_NR=abs(k_NR-NRsol(1));
    vec_k(iterationsNR)=NRsol(1);
    iterationsNR=iterationsNR+1;
    
end

n_NR=iterationsNR-1;
pos=1;

for i=1:(length(vec_k)-2)
    
    kNR_x(pos)=(abs(vec_k(i)-k_NR));
    kNR_y(pos)=(abs(vec_k(i+1)-k_NR));
    pos=1+pos;
    
end

%r: regression values for each of the N matrix rows
%m: slope of regression fit for each of the N matrix rows
%offset: offset of regression fit for each of the N matrix rows

[r,m,offset]=regression(log(kNR_x),log(kNR_y));
OC_NR=m;
ratio_NR=10^(offset);

%Gràfic Newton-Raphson

figure(1);
loglog(kNR_x,kNR_y,'linewidth',1.3);
grid on;
title('Convergència Mètode Newton-Raphson');
xlabel('k_n-L');
ylabel('k_{n+1}-L');

%% Valors Propis 
tic;
k_vp=k0;
n=0;
iterationsVP=10000; %numero d'iteracions
E=1;
while n<iterationsVP
    
    n=n+1;
    vp=eig(A(k_vp));
    tau_vp=vp(2);
    E=imag(tau_vp);
    
    if abs(E)<TOL
        break;
    end
    
    k_vp=k_vp+E;
    v_tau(n)=[real(tau_vp)];
    v_k(n)=[k_vp];

end
toc;
t_vp=toc;

% Convergència del mètode: Valors Propis

for i=1:(n-2)
    kvp_x(i)=v_k(i)-k_vp;
    kvp_y(i)=v_k(i+1)-k_vp;
end

[r,m,offset]=regression(log(abs(kvp_x(100:end-2))),log(abs(kvp_y(100:end-2))));
OC_vp=m;
ratio_vp=10^(offset);

% Gràfic Valors Propis
figure(2);
loglog(kvp_x,kvp_y,'linewidth',1.3);
grid on;
title('Convergència Mètode dels Valors Propis');
xlabel('k_n-L');
ylabel('k_{n+1}-L');

%% Mètode del Residu o Gradient

K1=x0;
Kn=[1e3;1e3]; % Per entrar al while
iteracionsResidu=1;
tic;
while abs(K1(1)-Kn(1))>tol(1) || abs(K1(2)-Kn(2))>tol(2)
    Kn=K1;
    MAT=[G(Kn(1),Kn(2)) ; H(Kn(1),Kn(2))];
    K1=Kn-alpha*eye(2)*MAT;
    vec_R(iteracionsResidu)=K1(1);
    iteracionsResidu=iteracionsResidu+1;
end

toc;
%k=bw/Uinf
%tau=(omegaalpha/w)^2
wres=omegah/sqrt(K1(2));
Uinfres=b*wres/K1(1);
t_Res=toc;
k_R=K1(1);

% Convergència del mètode: Residu o Gradient
for i=1:(iteracionsResidu-2)
    kR_x(i)=abs(vec_R(i)-k_R); %k estat 
    kR_y(i)=abs(vec_R(i+1)-k_R); %k estat i+1
end

kR_x=sort(kR_x);
kR_y=sort(kR_y);

[r,m,offset]=regression(log(kR_x(10:end-1)),log(kR_y(10:end-1)));
OC_R=m;
ratio_R=10^(offset);

% Gràfic Residu o Gradient
figure(3);
loglog(kR_x,kR_y,'linewidth',1.3);
grid on;
title('Convergència Mètode del Residu o Gradient');
xlabel('k_n-L');
ylabel('k_{n+1}-L');


%% Mètode Secant
%[n+1]=[n]-alpha[derivada punt inicial]*[G;H]
Ks1=x0;
Ks=[0;0];
iteracionsSec=1;

tic;
%Matriu derivades
DER=[(G(k0+0.0001,tau0)-G(k0,tau0))/0.0001 , (G(k0,tau0+0.0001)-G(k0,tau0))/0.0001  ; (H(k0+0.0001,tau0)-H(k0,tau0))/0.0001 , (H(k0,tau0+0.0001)-H(k0,tau0))/0.0001];



while abs(Ks1(1)-Ks(1))>tol(1) || abs(Ks1(2)-Ks(2))>tol(2)
    Ks=Ks1;
    MAT1=[G(Ks(1),Ks(2)) ; H(Ks(1),Ks(2))];
  Ks1=Ks-alpha*inv(DER)*MAT1;
  vec_S(iteracionsSec)=Ks1(1);
  iteracionsSec=iteracionsSec+1;
 
end

wSec=sqrt(omegah/Ks1(2));
UinfSec=b*wSec/Ks1(1);
toc;
tempsSec=toc;

k_S=Ks1(1);
% Convergència del mètode: Secant
for i=1:(iteracionsSec-2)
    kS_x(i)=abs(vec_S(i)-k_S); %k estat i
    kS_y(i)=abs(vec_S(i+1)-k_S); %k estat i+1
end



[r,m,b]=regression(log(kS_x(100:end-1)),log(kS_y(100:end-1)));

OC_S=m;
ratio_S=10^(b);

% Gràfic Secant
figure(4);
loglog(kS_x,kS_y,'linewidth',1.3);
grid on;
title('Convergència Mètode de la Secant');
xlabel('k_n-L');
ylabel('k_{n+1}-L');

%% Resultats

fprintf(' \n');
fprintf('NEWTON RAPHSON: \n');
fprintf('k =     %.4f  \n',NRsol(1));
fprintf('tau =   %.4f  \n',NRsol(2));
fprintf('Temps = %.4f s\n',t_NR);
fprintf('Ordre de convergència =  %f\n',OC_NR);
fprintf('Ratio de convergència =  %f\n',ratio_NR);

fprintf(' \n');
fprintf('VALORS PROPIS \n');
fprintf('k =     %.4f\n',k_vp);
fprintf('tau =   %.4f\n',tau_vp);
fprintf('Temps = %.4f s\n',t_vp);
fprintf('Ordre de convergència =  %.4f\n',OC_vp);
fprintf('Ratio de convergència =  %f\n',ratio_vp);

fprintf(' \n');
fprintf('RESIDU \n');
fprintf('k =    %.4f\n',K1(1));
fprintf('tau =  %.4f\n',K1(2));
fprintf('Temps =%.4f s\n',t_Res);
fprintf('Ordre de convergència =  %.4f\n',OC_R);
fprintf('Ratio de convergència =  %f\n',ratio_R);
% fprintf('VELOCITAT = %.4f\n',Uinfres);
% fprintf('FREQ =  %.4f\n',wres);

fprintf(' \n');
fprintf('SECANT \n');
fprintf('k =  %.4f\n',Ks1(1));
fprintf('tau =  %.4f\n',Ks1(2));
fprintf('VELOCITAT = %.4f\n',UinfSec);
fprintf('FREQ =  %.4f\n',wSec);
fprintf('temps =  %.4f\n',tempsSec);
fprintf('Ordre de convergència =  %.4f\n',OC_S);