import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


def eqsegundaordem(N,t,dt,c1):

    w = np.zeros(N+1)
    theta = np.zeros(N+1)

    theta[0] = 0.5
    w[0] = 0.0
    theta_aprox = np.copy(theta)
    w_aprox = np.copy(w)
    
    for k in range(N):
        
        theta[k+1] = theta[k] + dt*w[k]
        w[k+1] = w[k] + dt*(c1*np.sin(theta[k]))
        theta_aprox[k+1] = theta_aprox[k] + dt*w_aprox[k]
        w_aprox[k+1] = w_aprox[k] + dt*(c1*(theta_aprox[k]))
        
    return theta,theta_aprox

t_inicial = 0.0
t_final = 10.0
dt = 0.0001
L = 0.1
g = 9.81
c1 = -g/L
theta_inicial = 0.5
theta_linha_inicial = 0.0
N = int((t_final-t_inicial)/dt)
t = np.arange(t_inicial,t_final+dt,dt)

theta,theta_aprox = eqsegundaordem(N,t,dt,c1)

#Inicializando o Gráfico

graf = plt.figure()
ax = graf.add_subplot()
 
#Parâmetros do Gráfico

graf.suptitle('Ângulo do Pêndulo em relação à sua posição de equilíbrio ao longo do tempo', fontsize = 10, fontweight = 'bold', usetex = False)
ax.set_ylabel(r'theta', fontsize = 12, usetex = False)
ax.set_xlabel(r't', fontsize = 12, usetex = False)

blp = mpatches.Patch(color = 'b',label = 'Ângulo Aproximado')

blr = mpatches.Patch(color = 'r',label = 'Ângulo Sem aproximação')


plt.legend(handles = [blp,blr], fontsize = 8)
    
plt.plot(t,theta,'r')
plt.plot(t,theta_aprox,'b')
plt.show()