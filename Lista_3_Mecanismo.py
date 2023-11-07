import numpy as np
import matplotlib.pyplot as plt
import cmath

# Número USP:  11819827
# Final do número usp

N = 27
zeta = N/4

print(zeta)

# Comprimento do Elos (mm)

a = 40 - zeta  # Elo L2
b = 120 + zeta # Elo L3
c = 80 - zeta  # Elo L4
d = 100 + zeta # Elo L1

# Verificando a Condição de Grashof

S,L,P,Q=a,b,c,d

if (S+L)>(P+Q):
    print("Não Grashof, pois")
else:
    print("Grashof")

# Condições Iniciais

theta_2_ini = 40.0                    # Posição inicial de θ2 em graus
theta_2_ini = np.radians(theta_2_ini) # Posição inicial de θ2 em radianos
w_2_ini = 4*np.pi                     # Velocidade Angular inicial de θ2 
alpha_2_ini = 0.0                     # Aceleração Angular inicial de θ2  

print("As condições iniciais do mecanismos são {:.4f},{:.4f},{:.4f}".format(theta_2_ini,w_2_ini,alpha_2_ini))

# Análise da Posição

t = np.linspace(0,1,1000) #Tempo em segundos

# Variáveis Circuito Aberto

theta_2, theta_3_aberto,theta_4_aberto =[],[],[] # Posição
w_3_aberto, w_4_aberto = [],[]                   # Velocidade Angular
a_3_aberto, a_4_aberto = [],[]                   # Aceleração Angular

# Variáveis Circuito Cruzado

theta_3_cruza,theta_4_fecha = [],[]   # Posição
w_3_fecha, w_4_fecha = [],[]          # Velocidade Angular
a_3_cruza, a_4_cruza = [],[]          # Aceleração Angular

''' Analise da Posição Angular, Velocodade Angular, Aceleração Angular, Velociade Linear e Aceleração Linear'''

# Variáveis Auxiliares

K1, K2, K3, K4, K5 = d/a, d/c, (a*a-b*b+c*c+d*d)/(2*a*c), d/b, (c**2-d**2-a**2-b**2)/(2*a*c)
A,B,C,D,E,F = [],[],[],[],[],[]
KA_aberto, KB_aberto, KC_aberto, KD_aberto, KE_aberto, KF_aberto = [],[],[],[],[],[]
KA_cruza , KB_cruza , KC_cruza , KD_cruza , KE_cruza , KF_cruza  = [],[],[],[],[],[]
v_A, a_A, v_B_aberto,a_B_aberto, v_B_cruza, a_B_cruza = [],[],[],[],[],[]


for i in range(0,len(t),1):
    theta_2.append(theta_2_ini + w_2_ini*t[i])

    A.append(np.cos(theta_2[i])-K1-K2*np.cos(theta_2[i]) + K3)
    B.append(-2*np.sin(theta_2[i]))
    C.append(K1 - (K2+1)*np.cos(theta_2[i]) + K3)
    D.append(np.cos(theta_2[i])-K1+K4*np.cos(theta_2[i]) + K5)
    E.append(-2*np.sin(theta_2[i]))
    F.append(K1 + (K4 -1)*np.cos(theta_2[i]) + K5)

for i in range(0,len(t),1):

    '''Circuito Aberto'''

    # Posição

    theta_3_aberto.append(2*np.arctan((-E[i]-cmath.sqrt(E[i]*E[i]-4*D[i]*F[i]))/(2*D[i])))
    theta_4_aberto.append(2*np.arctan((-B[i]-cmath.sqrt(B[i]*B[i]-4*A[i]*C[i]))/(2*A[i])))

    # Velocidade Angular

    w_3_aberto.append(a*w_2_ini*np.sin(theta_4_aberto[i]-theta_2[i])/(b*np.sin(theta_3_aberto[i]-theta_4_aberto[i])))
    w_4_aberto.append(a*w_2_ini*np.sin(theta_2[i]-theta_3_aberto[i])/(b*np.sin(theta_4_aberto[i]-theta_3_aberto[i])))

    # Aceleração Angular

    KA_aberto.append(c*np.sin(theta_4_aberto[i]))
    KB_aberto.append(b*np.sin(theta_3_aberto[i]))
    KD_aberto.append(c*np.sin(theta_4_aberto[i]))
    KE_aberto.append(b*np.sin(theta_4_aberto[i]))
    KC_aberto.append(a*alpha_2_ini*np.sin(theta_2[i])+a*w_2_ini**2*np.cos(theta_2[i])+b*w_3_aberto[i]**2*np.cos(theta_3_aberto[i])-c*w_4_aberto[i]**2*np.cos(theta_4_aberto[i]))
    KF_aberto.append(a*alpha_2_ini*np.cos(theta_2[i])-a*w_2_ini**2*np.sin(theta_2[i])-b*w_3_aberto[i]**2*np.sin(theta_3_aberto[i])+c*w_4_aberto[i]**2*np.sin(theta_4_aberto[i]))
    a_3_aberto.append((KC_aberto[i]*KD_aberto[i]-KA_aberto[i]*KE_aberto[i])/(KA_aberto[i]*KE_aberto[i]-KB_aberto[i]*KD_aberto[i]))
    a_4_aberto.append((KC_aberto[i]*KE_aberto[i]-KB_aberto[i]*KF_aberto[i])/(KA_aberto[i]*KE_aberto[i]-KB_aberto[i]*KD_aberto[i]))


    '''Circuito Cruzado'''

    # Posição

    theta_3_cruza.append(2*np.arctan((-E[i]+cmath.sqrt(E[i]*E[i]-4*D[i]*F[i]))/(2*D[i])))
    theta_4_fecha.append(2*np.arctan((-B[i]+cmath.sqrt(B[i]*B[i]-4*A[i]*C[i]))/(2*A[i])))

    # Velocidade Angular

    w_3_fecha.append(a*w_2_ini*np.sin(theta_4_fecha[i]-theta_2[i])/(b*np.sin(theta_3_cruza[i]-theta_4_fecha[i])))
    w_4_fecha.append(a*w_2_ini*np.sin(theta_2[i]-theta_3_cruza[i])/(b*np.sin(theta_4_fecha[i]-theta_3_cruza[i])))

    # Aceleração Angular

    KA_cruza.append(c*np.sin(theta_4_fecha[i]))
    KB_cruza.append(b*np.sin(theta_3_cruza[i]))
    KD_cruza.append(c*np.sin(theta_4_fecha[i]))
    KE_cruza.append(b*np.sin(theta_4_fecha[i]))
    KC_cruza.append(a*alpha_2_ini*np.sin(theta_2[i])+a*w_2_ini**2*np.cos(theta_2[i])+b*w_3_fecha[i]**2*np.cos(theta_3_cruza[i])-c*w_4_fecha[i]**2*np.cos(theta_4_fecha[i]))
    KF_cruza.append(a*alpha_2_ini*np.cos(theta_2[i])-a*w_2_ini**2*np.sin(theta_2[i])-b*w_3_fecha[i]**2*np.sin(theta_3_cruza[i])+c*w_4_fecha[i]**2*np.sin(theta_4_fecha[i]))
    a_3_cruza.append((KC_cruza[i]*KD_cruza[i]-KA_cruza[i]*KE_cruza[i])/(KA_cruza[i]*KE_cruza[i]-KB_cruza[i]*KD_cruza[i]))
    a_4_cruza.append((KC_cruza[i]*KE_cruza[i]-KB_cruza[i]*KF_cruza[i])/(KA_cruza[i]*KE_cruza[i]-KB_cruza[i]*KD_cruza[i]))


    '''Ponto A e B'''

    # Ponto A

    v_A.append(a*w_2_ini*0.001)
    a_A.append((v_A[i]**2/a)*0.001)

    #Ponto B - Circuito Aberto

    v_B_aberto.append(w_4_aberto[i]*c*0.001)
    a_B_aberto.append(cmath.sqrt((a_4_aberto[i]*c)**2 + (w_4_aberto[i]**2*c)**2)*0.001)

    #Ponto B - Circuito Cruzado

    v_B_cruza.append(w_4_fecha[i]*c*0.001)
    a_B_cruza.append(cmath.sqrt((a_4_cruza[i]*c)**2 + (w_4_fecha[i]**2*c)**2)*0.001)

plt.figure(figsize=(15,15))

plt.subplot(2,2,1)
plt.plot(t, theta_3_aberto)
plt.xlabel('t')
plt.ylabel('θ3 - Circuito Aberto')
plt.subplot(2,2,2)
plt.plot(t, theta_4_aberto)
plt.xlabel('t')
plt.ylabel('θ4 - Circuito Aberto')
plt.subplot(2,2,3)
plt.plot(t, theta_3_cruza)
plt.xlabel('t')
plt.ylabel('θ3 - Circuito Fechado')
plt.subplot(2,2,4)
plt.plot(t, theta_4_fecha)
plt.xlabel('t')
plt.ylabel('θ4 - Circuito Fechado')
plt.suptitle('Posição Angular vs Tempo')
plt.show()

plt.figure(figsize=(15,15))
plt.subplot(2,2,1)
plt.plot(t, w_3_aberto)
plt.xlabel('t')
plt.ylabel('W3 - Circuito Aberto')
plt.subplot(2,2,2)
plt.plot(t, w_4_aberto)
plt.xlabel('t')
plt.ylabel('W4 - Circuito Aberto')
plt.subplot(2,2,3)
plt.plot(t, w_3_fecha)
plt.xlabel('t')
plt.ylabel('W3 - Circuito Fechado')
plt.subplot(2,2,4)
plt.plot(t, w_4_fecha)
plt.xlabel('t')
plt.ylabel('W4 - Circuito Fechado')
plt.suptitle('Velocidade Angular vs Tempo')
plt.show()

plt.figure(figsize=(15,15))
plt.subplot(2,2,1)
plt.plot(t, a_3_aberto)
plt.xlabel('t')
plt.ylabel('α3 - Circuito Aberto')
plt.subplot(2,2,2)
plt.plot(t, a_4_aberto)
plt.xlabel('t')
plt.ylabel('α4 - Circuito Aberto')
plt.subplot(2,2,3)
plt.plot(t, a_3_cruza)
plt.xlabel('t')
plt.ylabel('α3 - Circuito Fechado')
plt.subplot(2,2,4)
plt.plot(t, a_4_cruza)
plt.xlabel('t')
plt.ylabel('α4 - Circuito Fechado')
plt.suptitle('Aceleração Angular vs Tempo')
plt.show()

plt.figure(figsize=(30,30))
plt.subplot(3,2,1)
plt.plot(t, v_A)
plt.xlabel('t')
plt.ylabel('Velocidade Linear do Ponto A')
plt.subplot(3,2,2)
plt.plot(t, a_A)
plt.xlabel('t')
plt.ylabel('Aceleração Linear do Ponto A')
plt.subplot(3,2,3)
plt.plot(t, v_B_aberto)
plt.xlabel('t')
plt.ylabel('Velocidade Linear do Ponto B - Circuito Aberto')
plt.subplot(3,2,4)
plt.plot(t, v_B_cruza)
plt.xlabel('t')
plt.ylabel('Velocidade Linear do Ponto B - Circuito Cruzado')
plt.subplot(3,2,5)
plt.plot(t, a_B_aberto)
plt.xlabel('t')
plt.ylabel('Aceleração Linear do Ponto B - Circuito Aberto')
plt.subplot(3,2,6)
plt.plot(t, a_B_cruza)
plt.xlabel('t')
plt.ylabel('Aceleração Linear do Ponto B - Circuito Cruzado')
plt.suptitle('Velocidades e Acelerações Lineares dos pontos A e B')
plt.show()



