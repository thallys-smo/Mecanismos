% Número USP:  11819827
% Final do número USP: 27

N = 27;
zeta = N/4;

printf(zeta);

% Comprimento do Elos (mm)

a = 40 - zeta; %% Elo L2
b = 120 + zeta; %% Elo L3
c = 80 - zeta;  %% Elo L4
d = 100 + zeta; %% Elo L1

% Verificando a Condição de Grashof

S;L;P;Q=a;b;c;d

if (S+L)>(P+Q)
    printf("Não Grashof");
else
    printf("Grashof");
endif

% Condições Iniciais

theta_2_ini = 40.0;               % Posição inicial de θ2 em graus
theta_2_ini = theta_2_ini*pi/180; % Posição inicial de θ2 em radianos
w_2_ini = 4*pi;                   % Velocidade Angular inicial de θ2 
alpha_2_ini = 0.0;                % Aceleração Angular inicial de θ2  


%% Constantes Importantes

pL = 3;    %% Densidade Linear kg/m
g = 9.81;   %% Gravidade m/s^2

%% Massas; Pesos e Momentos de Inércia de cada Elo

m2 = a*pL; m3 = b*pL; m4 =c*pL;                       % Massa de cada Elo
p2 = m2*g; p3 = m3*g; p4 = m4*g;                      % Pesos de cada Elo
I2 = m2*(a**2)/3; I3 = m3*(b**2)/3; I4 = m4*(c**2)/3; % Momentos de Inércia de cada Elo

% Análise da Posição 

t = [0:0.01:1];

theta_3 = t ;theta_4 = t
theta_2