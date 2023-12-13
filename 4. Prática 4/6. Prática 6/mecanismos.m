% Curva Cicloidal
% Dados
clear;
beta = pi/3;
h = 1;
t = 0:0.001:6.001;
disp(length(t))
theta = pi/3 * t;

% Equacionamento
C_0 = 2*pi*(h/beta^2);
s = []; v = []; a = []; j = [];
for i=1:length(theta)/6
    s_theta = C_0 * beta / (2 * pi) * theta(i) - (C_0 * (beta^2 / (4 * pi^2))) * sin(2 * pi * theta(i) / beta);
    v_theta = -C_0*beta/(2*pi) * cos(2*pi*theta(i)/beta) + C_0*beta/(2*pi);
    a_theta = C_0*sin(2*pi*theta(i)/beta);
    j_theta = 4*pi^2*h/beta^3 * cos(2*pi*theta(i)/beta);
s = [s, s_theta];
v = [v, v_theta];
a = [a, a_theta];
j = [j, j_theta];
end

N= 27;
delta_t1 = 3-N/60;
delta_t2 = 1+N/60;

% Vetores
%t1 = 0:(2.55/1000):delta_t1; %% Repouso incial
%t2 = delta_t1:(delta_t2/1000):4; %% Subida
t1 = 0:0.001:1;
t2 = 1:0.001:4;
t3 = 4:0.001:5; %% Repouso intermediario
t4 = 5:0.001:6; %% Descida

disp(length(t1))
disp(length(t2))
disp(length(t3))
disp(length(t4))

% Repouso / subida / repouso / descida (inverso da subida − flip)
s_vetor = [zeros(size(t1)), s, ones(size(t3)), flip(s)];
v_vetor = [zeros(size(t1)), v, zeros(size(t3)), flip(v)];
a_vetor = [zeros(size(t1)), a, zeros(size(t3)), flip(a)];
j_vetor = [zeros(size(t1)), j, zeros(size(t3)), flip(j)];
disp(length(s_vetor));
disp(length(v_vetor));
disp(length(a_vetor));
disp(length(j_vetor));

% Plot SVAJ
figure(1), set(1,'position',[0 0 644 420]),
subplot(411), plot(t,s_vetor), ylabel('S (cm)'), axis tight,
subplot(412), plot(t,v_vetor), ylabel('V (cm/s)'), axis tight,
subplot(413), plot(t,a_vetor), ylabel('A (cm/s^2)'), axis tight,
subplot(414), plot(t,j_vetor), ylabel('J (cm/s^3)'), axis tight,
xlabel('Tempo (s)')

% Draw cam
figure(2), set(2,'position',[0 0 560 420])
subplot(221), Rb=h/2; polarplot(theta,Rb+s_vetor), title('R_b=h/2'),
subplot(222), Rb=h/1; polarplot(theta,Rb+s_vetor), title('R_b=h'),
subplot(223), Rb=2*h; polarplot(theta,Rb+s_vetor), title('R_b=2h'),
subplot(224), Rb=3*h; polarplot(theta,Rb+s_vetor), title('R_b=3h'),

% Polinomial
% Dados
clear;
beta = pi/3;
h = 1;
t = 0:0.001:6.001;
theta = pi/3 * t;

% Calculo das constantes
A = [
[1 1 1 1]
[4 5 6 7]
[12 20 30 42]
[24 60 120 210]
];
v = [h 0 0 0]';
x = A\v;
C_4 = x(1); C_5 = x(2); C_6 = x(3); C_7 = x(4);

% Equacionamento
s = []; v = []; a = []; j = [];
for i=1:length(theta)/6
    s_theta = C_4*(theta(i)/beta)^4 + C_5*(theta(i)/beta)^5 + C_6*(theta(i)/beta)^6 + C_7*(theta(i)/beta)^7;
    v_theta = 1/beta*(4*C_4*(theta(i)/beta)^3 + 5*C_5*(theta(i)/beta)^4 + 6*C_6*(theta(i)/beta)^5 + 7*C_7*(theta(i)/beta)^6);
    a_theta = 1/beta^2*(12*C_4*(theta(i)/beta)^2 + 20*C_5*(theta(i)/beta)^3 + 30*C_6*(theta(i)/beta)^4 + 42*C_7*(theta(i)/beta)^5);
    j_theta = 1/beta^3*(24*C_4*(theta(i)/beta)^1 + 60*C_5*(theta(i)/beta)^2 + 120*C_6*(theta(i)/beta)^3 + 210*C_7*(theta(i)/beta)^4);
s = [s, s_theta];
v = [v, v_theta];
a = [a, a_theta];
j = [j, j_theta];
end

% Vetores
t1 = 0:0.001:3; %% Repouso incial
t2 = 3:0.001:4; %% Subida
t3 = 4:0.001:5; %% Repouso intermediario
t4 = 5:0.001:6; %% Descida

% Repouso / subida / repouso / descida (inverso da subida − flip)
s_vetor = [zeros(size(t1)), s, ones(size(t3)), flip(s)];
v_vetor = [zeros(size(t1)), v, zeros(size(t3)), flip(v)];
a_vetor = [zeros(size(t1)), a, zeros(size(t3)), flip(a)];
j_vetor = [zeros(size(t1)), j, zeros(size(t3)), flip(j)];

% Plot SVAJ
figure(3), set(1,'position',[0 0 644 420]),
subplot(411), plot(t,s_vetor), ylabel('S (cm)'), axis tight,
subplot(412), plot(t,v_vetor), ylabel('V (cm/s)'), axis tight,
subplot(413), plot(t,a_vetor), ylabel('A (cm/s^2)'), axis tight,
subplot(414), plot(t,j_vetor), ylabel('J (cm/s^3)'), axis tight,
xlabel('Tempo (s)')

% Draw cam
figure(4), set(2,'position',[0 0 560 420])
subplot(221), Rb=h/2; polarplot(theta,Rb+s_vetor), title('R_b=h/2'),
subplot(222), Rb=h/1; polarplot(theta,Rb+s_vetor), title('R_b=h'),
subplot(223), Rb=2*h; polarplot(theta,Rb+s_vetor), title('R_b=2h'),
subplot(224), Rb=3*h; polarplot(theta,Rb+s_vetor), title('R_b=3h'),