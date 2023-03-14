import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#funcao desejada
def f(t,y,a,b):
    return a*y-b*y**2

#funcao para o calculo de y pelo metodo de euler explicito
def euler_explicito(N,t,dt,a,b):
    y = np.zeros(N+1)
    y[0] = 100000.0
    for k in range(N):
        y[k+1] = y[k] + dt*f(t[k], y[k],a,b)
    return y

#funcao para o calculo de y pelo metodo de range-kutta de quarta ordem
def rk4(N,t,dt,a,b):
    y = np.zeros(N+1)
    y[0] = 100000.0
    for k in range(N):
          
        k1 = dt*f(t[k],y[k],a,b)
        
        k2 = dt*f(t[k]+dt/2,y[k]+k1/2,a,b)
        
        k3 = dt*f(t[k]+dt/2,y[k]+k2/2,a,b)
        
        k4 = dt*f(t[k]+dt,y[k]+k3,a,b)
        
        y[k+1] = y[k] + (1/6)*(k1+2*k2+2*k3+k4)

    return y

dt = 1.0

t_inicial = 0.0
t_final = 20.0
N = int((t_final-t_inicial)/dt)
t = np.linspace(t_inicial,t_final,N+1)
a = 0.1
b = 10.0**(-7)


y_euler_explicito = euler_explicito(N,t,dt,a,b)
y_rk4 = rk4(N,t,dt,a,b)

#Inicializando o Gráfico

graf = plt.figure()
ax = graf.add_subplot()
 
#Parâmetros do Gráfico

graf.suptitle('Evolução da população ao longo do tempo', fontsize = 10, fontweight = 'bold', usetex = False)
ax.set_ylabel(r'N(t)', fontsize = 12, usetex = False)
ax.set_xlabel(r't', fontsize = 12, usetex = False)
blp = mpatches.Patch(color = 'b',label = 'População Euler Explícito')

blr = mpatches.Patch(color = 'r',label = 'População Range-Kutta de 4º Ordem')

plt.legend(handles = [blp,blr], fontsize = 8)
 
plt.plot(t,y_euler_explicito,'b')
plt.plot(t,y_rk4,'r')
plt.show()