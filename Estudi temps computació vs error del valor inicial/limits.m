% Limits de converg�ncia %

close all;
clear all;
clc;

%% Data
tausolucio=0.5571;
ksolucio=0.3555;

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



%% Newton-Raphson
convergencia=true;
maxIter=1000;
MAXiter=10000;
MAXtempsiter=5;
i=1;

tempsNR=zeros(MAXiter);
iteracionsNR=zeros(MAXiter);
errorInicialNR=zeros(MAXiter);

tau0=tausolucio-0.1;
k0=ksolucio-0.1;
x0=[k0;tau0];

while convergencia && i<MAXiter
    errorInicialNR(i)=abs(k0-ksolucio)*100/ksolucio;
    errorInicialNR(i)
    
    tic;
    options=optimset('Display','iter','TolFun',TOL,'MaxIter',maxIter);
    [NRsol,exitflag]=fsolve(fun,x0,options);
    toc;
    
    iteracionsNR(i)=exitflag(1);
    tempsNR(i)=toc;
    
    if iteracionsNR(i)>=maxIter
        convergencia=false;
        fprintf('Superat l�mit de iteracions \n')
    end
    if tempsNR(i)>=MAXtempsiter
        convergencia=false;
        fprintf('Superat temps de iteracio \n')
    end
    
    i=i+1;
    k0=k0-0.1;
    tau0=tau0-0.1;
    x0=[k0;tau0];
    
   
end
tempsNR=nonzeros(tempsNR);
errorInicialNR=nonzeros(errorInicialNR);
iteracionsNR=nonzeros(iteracionsNR);

figure(1)
subplot(2,1,1)
title('Temps de computaci� i n�mero de iteracions vs error inicial pel m�tode de Newton Raphson')
plot(errorInicialNR,tempsNR)
ylabel('Temps de computaci� (s)');
subplot(2,1,2)
plot(errorInicialNR,iteracionsNR,'color','r')
xlabel('Error inicial (%)');
ylabel('Iteracions');


%% Valors Propis 
convergencia=true;
i=1;
maxIter=100000;
MAXiter=10000;
MAXtempsiter=5;

tempsVp=zeros(MAXiter);
iteracionsVp=zeros(MAXiter);
errorInicialVp=zeros(MAXiter);

tau0=tausolucio;
k0=ksolucio-0.1;
x0=[k0;tau0];
while convergencia && i<MAXiter
    errorInicialVp(i)=abs(k0-ksolucio)*100/ksolucio;
    errorInicialVp(i)
    tic;
    k_vp=k0;
    n=0;
    E=1;
    while n<=maxIter && abs(E)>TOL

        n=n+1;
        vp=eig(A(k_vp));
        tau_vp=vp(2);
        E=imag(tau_vp);

        k_vp=k_vp+E;
        v_tau(n)=[real(tau_vp)];
        v_k(n)=[k_vp];

    end  
    toc;
    tempsVp(i)=toc;
    iteracionsVp(i)=n;
    if iteracionsVp(i)>=maxIter
        convergencia=false;
        fprintf('Superat l�mit de iteracions \n')
    end
    if tempsVp(i)>=MAXtempsiter
        convergencia=false;
        fprintf('Superat temps de iteracio \n')
    end
    
    i=i+1;
    k0=k0-0.1;
end

tempsVp=nonzeros(tempsVp);
errorInicialVp=nonzeros(errorInicialVp);
iteracionsVp=nonzeros(iteracionsVp);

figure(2)
subplot(2,1,1)
title('Temps de computaci� i n�mero de iteracions vs error inicial pel m�tode de valors pr�pis')
plot(errorInicialVp,tempsVp)
ylabel('Temps de computaci� (s)');
subplot(2,1,2)
plot(errorInicialVp,iteracionsVp,'color','r')
xlabel('Error inicial(%)');
ylabel('Iteracions');

%% M�tode del Residu o Gradient
convergencia=true;
i=1;
maxIter=100000;
MAXiter=10000;
MAXtempsiter=5;

tempsRes=zeros(MAXiter);
iteracionsRes=zeros(MAXiter);
errorInicialRes=zeros(MAXiter);

tau0=tausolucio-0.1;
k0=ksolucio-0.1;
x0=[k0;tau0];

while convergencia && i<MAXiter
    errorInicialRes(i)=abs(k0-ksolucio)*100/ksolucio;
    errorInicialRes(i)
    K1=x0;
    Kn=[1e3;1e3]; % Per entrar al while
    n=0;
    tic;
    while (abs(K1(1)-Kn(1))>tol(1) || abs(K1(2)-Kn(2))>tol(2)) && n<=maxIter
        n=n+1;
        Kn=K1;
        MAT=[G(Kn(1),Kn(2)) ; H(Kn(1),Kn(2))];
        K1=Kn-alpha*eye(2)*MAT;
        vec_R(n)=K1(1);

    end

    toc;
    tempsRes(i)=toc;
    iteracionsRes(i,1)=n;
    if iteracionsRes(i)>=maxIter
        convergencia=false;
        fprintf('Superat l�mit de iteracions \n')
    end
    if tempsRes(i)>=MAXtempsiter
        convergencia=false;
        fprintf('Superat temps de iteracio \n')
    end
    
    i=i+1;
    k0=k0-0.1;
    tau0=tau0-0.1;
    x0=[k0;tau0];
end
tempsRes=nonzeros(tempsRes);
errorInicialRes=nonzeros(errorInicialRes);
iteracionsRes=nonzeros(iteracionsRes);

figure(3)
subplot(2,1,1)
title('Temps de computaci� i n�mero de iteracions vs error inicial pel m�tode del gradient')
plot(errorInicialRes,tempsRes)
ylabel('Temps de computaci� (s)');
subplot(2,1,2)
plot(errorInicialRes,iteracionsRes,'color','r')
xlabel('Error inicial (%)');
ylabel('Iteracions');


%% M�tode Secant
% 
% convergencia=true;
% i=1;
% maxIter=100000;
% MAXiter=10000;
% MAXtempsiter=5;
% 
% tempsSec=zeros(MAXiter);
% iteracionsSec=zeros(MAXiter);
% errorInicialSec=zeros(MAXiter);
% 
% tau0=tausolucio-0.1;
% k0=ksolucio-0.1;
% x0=[k0;tau0];
% 
% while convergencia && i<MAXiter
%     
%     errorInicialSec(i)=abs(k0-ksolucio)*100/ksolucio;
%     errorInicialSec(i)
%     
%     Ks1=x0;
%     Ks=[0;0];
%     n=0;
% 
%     tic;
%     %Matriu derivades
%     DER=[(G(k0+0.0001,tau0)-G(k0,tau0))/0.0001 , (G(k0,tau0+0.0001)-G(k0,tau0))/0.0001  ; (H(k0+0.0001,tau0)-H(k0,tau0))/0.0001 , (H(k0,tau0+0.0001)-H(k0,tau0))/0.0001];
% 
% 
% 
%     while (abs(Ks1(1)-Ks(1))>tol(1) || abs(Ks1(2)-Ks(2))>tol(2)) && n<=maxIter
%         n=n+1;
%         Ks=Ks1;
%         MAT1=[G(Ks(1),Ks(2)) ; H(Ks(1),Ks(2))];
%       Ks1=Ks-alpha*inv(DER)*MAT1;
%       vec_S(n)=Ks1(1);
%       n=n+1;
% 
%     end
%     toc;
%     tempsSec(i)=toc;
%     
%     iteracionsSec(i,1)=n;
%     if iteracionsSec(i)>=maxIter
%         convergencia=false;
%         fprintf('Superat l�mit de iteracions \n')
%     end
%     if tempsSec(i)>=MAXtempsiter
%         convergencia=false;
%         fprintf('Superat temps de iteracio \n')
%     end
%     
%     i=i+1;
%     k0=k0-0.1;
%     tau0=tau0-0.1;
%     x0=[k0;tau0];
% 
% end
% 
% figure(4)
% subplot(2,1,1)
% title('Temps de computaci� i n�mero de iteracions vs error inicial pel m�tode de la secant')
% plot(errorInicialSec,tempsSec)
% ylabel('Temps de computaci� (s)');
% subplot(2,1,2)
% plot(errorInicialSec,iteracionsSec,'color','r')
% xlabel('Error inicial (%)');
% ylabel('Iteracions');
