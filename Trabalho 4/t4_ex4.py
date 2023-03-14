import numpy as np
import matplotlib.pyplot as plt
from numba import njit

@njit
def f(t,x):
    res = 0.0
    for n in range(100):
        res+= (((-1)**n)/((2*n+1)**2))*np.sin((2*n+1)*np.pi*x)*np.exp(-((2*n+1)**2)*(np.pi**2)*0.01*t)
    return (800/(np.pi**2))*res

@njit
def T_i(T,x,dx):
    T[0] = 0.0
    for i in range(1,N+1):
        x+=dx
        T[i] = 200*x
    return T

@njit
def solexata(N,T,tempo_final):
    x_aux = 0.0
    for i in range(1,N+1):
        x_aux+=dx
        T[i] = f(tempo_final,x_aux)
        
    return T

@njit
def btcs(R,T,T_new,lamda,dt,t_final,t,w):
    iter = 0
    while t < t_final:
        R[1] = 1.0
        while np.max(np.abs(R)) > 1e-8:
            for i in range(1, N+1):
                if i != N:
                    R[i] = (1.0/(1.0+0.02*lamda))*(T[i] + 0.01*lamda*(T_new[i+1]+T_new[i-1])
                    - (1.0+0.02*lamda)*T_new[i])
                    T_new[i] = T_new[i] + w*R[i]
                else:
                    T_new[i] = (4/3)*(T_new[i-1]) - (1/3)*T_new[i-2]
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
        while np.max(np.abs(R)) > 1e-8:
            for i in range(1, N+1):
                if i != N:
                    R[i] = ((1.0-0.01*lamda)/(1.0+0.01*lamda))*T[i] + ((0.005*lamda)/(1.0+0.01*lamda))*(T_new[i+1]+T_new[i-1]) + ((0.005*lamda)/(1.0+0.01*lamda))*(T[i+1]+T[i-1]) - T_new[i]
                    #R[i] = (T[i] + (lamda/2.0)*(T_new[i+1] +T_new[i-1]) + (lamda/2.0)*(T[i+1]-2.0*T[i]+T[i-1]) - T_new[i])/(1+lamda)
                    #R[i] = (T[i] + (lamda/2.0)*(T[i+1]-2.0*T[i]+T[i-1]) + (lamda/2.0)*(T[i+1]+T[i-1]))/(1+lamda) - T_new[i]
                    #R[i] = (1.0/(2.0+2.0*lamda))*((2.0-2.0*lamda)*T[i]+ lamda*(T[i+1]+T[i-1]+T_new[i+1]+T_new[i-1]) - (2.0+2.0*lamda)*T_new[i])
                    T_new[i] = T_new[i] + w*R[i]
                else:
                    T_new[i] = (4/3)*(T_new[i-1]) - (1/3)*T_new[i-2]
            iter+=1
        t+=dt
        T = np.copy(T_new)
    #print(iter)
    return T        

N = 10 #numero de passos 10,100,200
N_curvas = 3 #numero de curvas
t = 0.0
t_final = 5.0
dx = 0.5/N
dt = 0.0001*t_final
x = 0.0
lamda = dt/(dx*dx)
w = 1.0
  
T = np.zeros(N+1)
T_exata = np.zeros(N+1)
R = np.zeros(N+1)

T = T_i(T,x,dx)
T_new = np.copy(T)
           
#Inicializando o Grafico

graf = plt.figure()
ax = graf.add_subplot()

#Parâmetros do Grafico

graf.suptitle('Comparação da Resolução da Equação do calor analítica e numérica', fontsize = 10, fontweight = 'bold', usetex = False)
ax.set_ylabel(r'erro', fontsize = 12, usetex = False)
ax.set_xlabel(r't', fontsize = 12, usetex = False)
    
#Calculando a Equacao da Conducao do Calor
    
for i in range(1,N_curvas+1):
    t = 0.0
    t_new = (t_final/N_curvas)*i
    dt_new = ((dx**2)/4)
    lamda_new = dt_new/(dx*dx)
    T_aux1 = np.copy(T)
    T_aux2 = np.copy(T)
    T_new2 = np.copy(T_new)
    T1 = btcs(R,T_aux1,T_new,lamda_new,dt_new,t_new,t,w)
    T3 = cranknicolson(R,T_aux2,T_new2,lamda_new,dt_new,t_new,t,w)
    
    #plt.plot(T3,'or')
    #plt.plot(T1,'pb')
    
T2 = solexata(N,T_exata,t_final)
#plt.plot(T2,'k')

erro_btcs = np.abs(T1 - T2)
erro_cn = np.abs(T3 - T2)

plt.plot(erro_btcs,'b')
plt.plot(erro_cn,'r')
plt.show()