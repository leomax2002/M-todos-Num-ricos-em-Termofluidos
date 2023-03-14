import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def f(t):
    return np.exp(-0.06*np.pi*t)*np.sin(2.0*t-np.pi)

def df(t):
    return 0.06*np.pi*np.exp(-0.06*np.pi*t)*np.sin(2*t)-2*np.exp(-0.06*np.pi*t)*np.cos(2*t)

def df2(t):
    return 3.96447*np.exp(-0.06*np.pi*t)*np.sin(2*t) + 0.753982*np.exp(-0.06*np.pi*t)*np.cos(2*t)
    

def eqsegundaordem(N,t,dt,R,L,C):
    
    i = np.zeros(N+1)
    E = np.zeros(N+1)
    dE = np.zeros(N+1)
    dE2 = np.zeros(N+1)

    i[0] = 0.0
    
    for k in range(N):
        
        i[k+1] = i[k] + dt*((C*df2(t[k])) + (1/R)*(df(t[k]))+(1/L)*f(t[k]))
        
    return i

def rksegundaordem(N,t,dt,R,L,C):
    i = np.zeros(N+1)

    i[0] = 0.0
    for k in range(N):
          
        a = dt*((C*df2(t[k])) + (1/R)*(df(t[k]))+(1/L)*f(t[k]))
        
        b = dt*((C*df2(t[k]+dt/2)) + (1/R)*(df(t[k]+dt/2))+(1/L)*f(t[k]+dt/2))
        
        c = dt*((C*df2(t[k]+dt/2)) + (1/R)*(df(t[k]))+(1/L)*f(t[k]+dt/2))
        
        d = dt*((C*df2(t[k])) + (1/R)*(df(t[k]))+(1/L)*f(t[k]))
        
        i[k+1] = i[k] + (1/6)*(a+2*b+2*c+d)

    return i

t_inicial = 0.0
t_final = 10.0
dt = 0.001
C = 0.3
R = 1.4
L = 1.7
io = 0.0
N = int((t_final-t_inicial)/dt)
t = np.arange(t_inicial,t_final+dt,dt)

i = eqsegundaordem(N,t,dt,R,L,C)
i_rk4 = rksegundaordem(N,t,dt,R,L,C)

#Inicializando o Gráfico

graf = plt.figure()
ax = graf.add_subplot()
 
#Parâmetros do Gráfico

graf.suptitle('Corrente no Cicuito RLC a partir dos métodos de Euler Explícito e Range-Kutta de Quarta Ordem', fontsize = 10, fontweight = 'bold', usetex = False)
ax.set_ylabel(r'i', fontsize = 12, usetex = False)
ax.set_xlabel(r't', fontsize = 12, usetex = False)
blp = mpatches.Patch(color = 'b',label = 'Corrente Euler Explícito')

blr = mpatches.Patch(color = 'r',label = 'Corrente Range-Kutta de 4º Ordem')

plt.legend(handles = [blp,blr], fontsize = 8)

    
plt.plot(t,i,'b')
plt.plot(t,i_rk4,'r')
plt.show()