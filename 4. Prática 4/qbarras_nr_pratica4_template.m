clear all, close all, clc,
%% Simulacao de mecanismo de quatro barras
outvid=0;
% Lengths of links
a=??; b=??; c=??; d=??; % Grashof Crank-Rocker
% End-effector P
AP=??; tAP=??;
% Number of steps
N=??;
% Vector of time instants
t=linspace(0,1,N)';

%% Setup of input
% Initial t2
t2i=??;
% Vector of values for t2 (for imposed t2 in constant steps)
% t2v=t2i+linspace(0,4*pi,N)';
  % Vector of values for t2 (for imposed w2 constant)
   w2=??; t2v=t2i+w2*t; w2v=w2*ones(size(t)); a2v=zeros(size(t));
    % Vector of values for t2 (for imposed a2 constant)
    % a2=8*pi; t2v=t2i+a2/2*t.^2; a2v=a2*ones(size(t)); w2v=a2*t;

%% Newton-Raphson algorithm for the evaluation of t3 and t4 given t2
tol=??;
% Initial guesses for t3 and t4
  t3=??; t4=??; % Circuito aberto
  %t3=??; t4=??; % Circuito cruzado
disp('    it2   Iterations')
for it2=1:length(t2v)
   t2=t2v(it2); B=tol+1; iconv=0;
   while norm(B)>tol
       iconv=iconv+1;
       A=[-b*sin(t3) c*sin(t4);b*cos(t3) -c*cos(t4)];
       B=[a*cos(t2)+b*cos(t3)-c*cos(t4)-d; a*sin(t2)+b*sin(t3)-c*sin(t4)];
       Dt=-A\B;
       t3=t3+Dt(1); t4=t4+Dt(2);
   end
   if iconv>2 disp([it2 iconv]), end % Show number of iterations required to converge
   t3v(it2,1)=t3; t4v(it2,1)=t4;
end

%% Show simulation
%qbarras_nr_simulation

%% Evaluate positions of points O2,A,B,O4,P
rO2=zeros(length(t2v),2);
rA=a*[cos(t2v) sin(t2v)];
rB=rA+b*[cos(t3v) sin(t3v)];
rO4=[rO2(:,1)+d rO2(:,2)];
rP=rA+AP*[cos(t3v+tAP) sin(t3v+tAP)];

%% Post-processing for velocities and accelerations
qbarras_nr_velacel

%% Post-processing for forces and torques
% Position vectors (assuming bar-like uniform links)
rCG2=(a/2)*[cos(t2v) sin(t2v)];
R12=rO2-rCG2; R12x=R12(:,1); R12y=R12(:,2);
R32=rA-rCG2; R32x=R32(:,1); R32y=R32(:,2);

rCG3=rA+(b/2)*[cos(t3v) sin(t3v)];
R23=rA-rCG3; R23x=R23(:,1); R23y=R23(:,2);
R43=rB-rCG3; R43x=R43(:,1); R43y=R43(:,2);

rCG4=rO4+(c/2)*[cos(t4v) sin(t4v)];
R34=rB-rCG4; R34x=R34(:,1); R34y=R34(:,2);
R14=rO4-rCG4; R14x=R14(:,1); R14y=R14(:,2);

RP=rP-rCG3; RPx=RP(:,1); RPy=RP(:,2);

% Resisting forces and torques
FPx=??; FPy=??; T4=??; ipeso=??;

% Center of mass accelerations (assuming bar-like uniform links)
aG2=(a/2)*(j*a2v-w2v.^2).*exp(j*t2v); aG2x=real(aG2); aG2y=imag(aG2);
aG3=aB+(b/2)*(j*a3v-w3v.^2).*exp(j*t3v); aG3x=real(aG3); aG3y=imag(aG3);
aG4=(c/2)*(j*a4v-w4v.^2).*exp(j*t4v); aG4x=real(aG4); aG4y=imag(aG4);

% Inertia properties of links (assuming bar-like uniform links)
rhoAst=??; g=9.81;
m2=rhoAst*a; IG2=m2*a^2/12;
m3=rhoAst*b; IG3=m3*b^2/12;
m4=rhoAst*c; IG4=m4*c^2/12;

% Solution of dynamic problem (ADAPTAR ? ESCRITA DAS EQUA??ES E INC?GNITAS ESCOLHIDAS)
for id=1:length(t2v)
    Af=[??];
    Bf=[??];
    xf=Af\Bf;
    F12(id,[1 2])=xf([1 2]); %%% ADAPTAR (ANIMA??O USA SA?DA F12)
    F32(id,[1 2])=xf([3 4]); %%% ADAPTAR
    F43(id,[1 2])=xf([5 6]); %%% ADAPTAR
    F14(id,[1 2])=xf([7 8]); %%% ADAPTAR (ANIMA??O USA SA?DA F14)
    T12(id,:)=xf(9); %%% ADAPTAR
end

% Plot evolution of forces and torques
figure(31), set(31,'position',[692 1 560 840])
subplot(411)
plot(t,T12,'k-'),
set(gca,'xlim',[min(t) max(t)])
xlabel('Tempo (s)'), ylabel('$T_{12}$ (Nm)'),
subplot(412)
plot(t,-F12(:,1),'r-',t,-F12(:,2),'m--'),
set(gca,'xlim',[min(t) max(t)])
xlabel('Tempo (s)'), ylabel('$F_{21}=-F_{12}$ (N)'), legend('$F_{21x}$','$F_{21y}$','Location','NorthWest')
subplot(413)
plot(t,-F14(:,1),'b-',t,-F14(:,2),'c--'),
set(gca,'xlim',[min(t) max(t)])
xlabel('Tempo (s)'), ylabel('$F_{41}=-F_{14}$ (N)'), legend('$F_{41x}$','$F_{41y}$','Location','NorthWest')
subplot(414)
plot(t,-F12(:,1)-F14(:,1),'g-',t,-F12(:,2)-F14(:,2),'b--'),
set(gca,'xlim',[min(t) max(t)])
xlabel('Tempo (s)'), ylabel('$F_s$ (N)'), legend('$F_{sx}$','$F_{sy}$','Location','NorthWest')
% COMO INCLUIR PLOT DO TORQUE NO SUPORTE???

% Print extreme values of forces and torques
disp('Valores maximos e minimos')
disp(sprintf('T12: Max=%.4f, Min=%.4f',max(T12),min(T12)))
disp(sprintf('F21: Max=%.4f, Min=%.4f',max(sqrt(F12(:,1).^2+F12(:,2).^2)),min(sqrt(F12(:,1).^2+F12(:,2).^2))))
disp(sprintf('F41: Max=%.4f, Min=%.4f',max(sqrt(F14(:,1).^2+F14(:,2).^2)),min(sqrt(F14(:,1).^2+F14(:,2).^2))))
% COMO CALCULAR M?XIMO TORQUE NO SUPORTE???

%% Show simulation with forces
qbarras_nr_simulation_forces
