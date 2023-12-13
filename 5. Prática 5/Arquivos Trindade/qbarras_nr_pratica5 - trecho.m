clear all, close all, clc,
%% Four bar linkage synthesis
% Position 2 (w/u e z/s)
%   eq1:w*cos(t)*(cos(b2)-1)-w*sin(t)*sin(b2)+z*cos(p)*(cos(a2)-1)-z*sin(p)*sin(a2)=p21*cos(d2);
%   eq2:w*sin(t)*(cos(b2)-1)+w*cos(t)*sin(b2)+z*sin(p)*(cos(a2)-1)+z*cos(p)*sin(a2)=p21*sin(d2);
% Position 3 (w/u e z/s)
%   eq3:w*cos(t)*(cos(b3)-1)-w*sin(t)*sin(b3)+z*cos(p)*(cos(a3)-1)-z*sin(p)*sin(a3)=p31*cos(d3);
%   eq4:w*sin(t)*(cos(b3)-1)+w*cos(t)*sin(b3)+z*sin(p)*(cos(a3)-1)+z*cos(p)*sin(a3)=p31*sin(d3);

% Dados (Posicoes P1,P2,P3 e Orientacoes o1,o2,o3)
P1=[ 0.000,0.000]; o1=210.0/180*pi; 
P2=[-1.236,2.138]; o2=147.5/180*pi; 
P3=[-2.500,2.931]; o3=110.2/180*pi;
p21=norm(P2-P1);
d2=atan2((P2(2)-P1(2)),(P2(1)-P1(1)));
a2=o2-o1;
p31=norm(P3-P1);
d3=atan2((P3(2)-P1(2)),(P3(1)-P1(1)));
a3=o3-o1;
b2=30/180*pi; g2=-10/180*pi; b3=60/180*pi; g3=25/180*pi; % Caso 1

% Matrix equations
% Se dados p21,d2,a2,p31,d3,a3, e b2,b3 : calcular w1x,w1y,z1x,z1y
