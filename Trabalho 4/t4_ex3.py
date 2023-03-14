import numpy as np
import matplotlib.pyplot as plt
from numba import njit

@njit
def f(t,x):
    res = 0.0
    for n in range(1,1000+1):
        res+=((2.0/(n*np.pi))*(np.sin(n*np.pi*x)*np.exp(-(n**2)*(np.pi**2)*t)))
    return 1.0 - x - res

@njit
def T_i(T,x,dx):
    T[0] = 1.0
    T[N] = 0.0
    for i in range(1,N):
        x+=dx
        T[i] = 0.0
    return T

@njit
def solexata(N,T,tempo_final):
    x_aux = 0.0
    for i in range(1,N):
        x_aux+=dx
        T[i] = f(tempo_final,x_aux)
        
    return T

@njit
def btcs(R,T,T_new,lamda,dt,t_final,t,w):
    iter = 0
    while t < t_final:
        R[1] = 1.0
        while np.max(np.abs(R)) > 1e-6:
            for i in range(1, N):
                R[i] = (1.0/(1.0+2.0*lamda))*(T[i] + lamda*(T_new[i+1]+T_new[i-1])
                - (1.0+2.0*lamda)*T_new[i])
                T_new[i] = T_new[i] + w*R[i]
        iter+=1
        t+=dt
        T = np.copy(T_new)
    #print(iter)
    return T        
@njit    
def cranknicolson(R,T,T_new,lamda,dt,t_final,t,w):
    iter = 0
    while t < t_final:
        R[1] = 1.0
        while np.max(np.abs(R)) > 1e-6:
            for i in range(1, N):
                R[i] = ((1.0-lamda)/(1.0+lamda))*T[i] + ((lamda)/(2.0+2.0*lamda))*(T_new[i+1]+T_new[i-1]) + ((lamda)/(2.0+2.0*lamda))*(T[i+1]+T[i-1]) - T_new[i]
                #R[i] = (T[i] + (lamda/2.0)*(T_new[i+1] +T_new[i-1]) + (lamda/2.0)*(T[i+1]-2.0*T[i]+T[i-1]) - T_new[i])/(1+lamda)
                #R[i] = (T[i] + (lamda/2.0)*(T[i+1]-2.0*T[i]+T[i-1]) + (lamda/2.0)*(T[i+1]+T[i-1]))/(1+lamda) - T_new[i]
                #R[i] = (1.0/(2.0+2.0*lamda))*((2.0-2.0*lamda)*T[i]+ lamda*(T[i+1]+T[i-1]+T_new[i+1]+T_new[i-1]) - (2.0+2.0*lamda)*T_new[i])
                T_new[i] = T_new[i] + w*R[i]
            iter+=1
        t+=dt
        T = np.copy(T_new)
    #print(iter)
    return T        

N = 200 #numero de passos, 10,100,200
N_curvas = 3 #numero de curvas
t = 0.0
t_final = 0.2 #0.2
dx = 1.0/N
dt = 0.001*t_final
x = 0.0
lamda = dt/(dx*dx)
w = 1.0
  
T = np.zeros(N+1)
T_exata = np.zeros(N+1)
R = np.zeros(N+1)

T = T_i(T,x,dx)
T_new = np.copy(T)
       
T_exata[0] = 1.0
T_exata[N] = 0.0
    
#Inicializando o Grafico

graf = plt.figure()
ax = graf.add_subplot()

#Parâmetros do Grafico

graf.suptitle('Comparação da Resolução da Equação do calor analítica e numérica', fontsize = 10, fontweight = 'bold', usetex = False)
ax.set_ylabel(r'erro', fontsize = 12, usetex = False)
ax.set_xlabel(r't', fontsize = 12, usetex = False)
plt.xscale('log')
plt.yscale('log')

    
#Calculando a Equacao da Conducao do Calor
    
for i in range(1,N_curvas+1):
    t = 0.0
    t_new = (t_final/N_curvas)*i
    dt_new = 0.001*t_new
    lamda_new = dt_new/(dx*dx)
    T_aux1 = np.copy(T)
    T_aux2 = np.copy(T)
    T_new2 = np.copy(T_new)
    T1 = btcs(R,T_aux1,T_new,lamda_new,dt_new,t_new,t,w)
    T3 = cranknicolson(R,T_aux2,T_new2,lamda_new,dt_new,t_new,t,w)
    
    #plt.plot(T3,'or')
    #plt.plot(T1,'pb')


T2 = solexata(N,T_exata,t_final)
erro_btcs = np.abs(T1-T2)
erro_cn = np.abs(T3-T2)
#plt.plot(T2,'k', label = 'Analítica')

plt.plot(erro_btcs,'b')
plt.plot(erro_cn,'r')

plt.show()