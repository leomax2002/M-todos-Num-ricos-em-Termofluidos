import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#funcao desejada
def f1(t,y,a1,b1,c1,a2,b2,c2,y2):
    return y*(a1-b1*y-c1*y2)

def f2(t,y,a1,b1,c1,a2,b2,c2,y2):
    return y*(a2-b2*y-c2*y2)


#funcao para o calculo de y pelo metodo de euler explicito
def euler_explicito(N,t,dt,a1,b1,c1,a2,b2,c2):
    y1 = np.zeros(N+1)
    y2 = np.zeros(N+1)
    y1[0] = 100000.0
    y2[0] = 100000.0
    for k in range(N):
        y1[k+1] = y1[k] + dt*f1(t[k],y1[k],a1,b1,c1,a2,b2,c2,y2[k])
        y2[k+1] = y2[k] + dt*f2(t[k],y2[k],a1,b1,c1,a2,b2,c2,y1[k])
        
    return y1,y2

#funcao para o calculo de y pelo metodo de range-kutta de quarta ordem
def rk4(N,t,dt,a1,b1,c1,a2,b2,c2):
    y1 = np.zeros(N+1)
    y2 = np.zeros(N+1)
    y1[0] = 100000.0
    y2[0] = 100000.0
    
    for k in range(N):
          
        k11 = dt*f1(t[k],y1[k],a1,b1,c1,a2,b2,c2,y2[k])
        
        k12 = dt*f1(t[k]+dt/2,y1[k]+k11/2,a1,b1,c1,a2,b2,c2,y2[k])
        
        k13 = dt*f1(t[k]+dt/2,y1[k]+k12/2,a1,b1,c1,a2,b2,c2,y2[k])
        
        k14 = dt*f1(t[k]+dt,y1[k]+k13,a1,b1,c1,a2,b2,c2,y2[k])
        
        
        k21 = dt*f2(t[k],y2[k],a1,b1,c1,a2,b2,c2,y1[k])
        
        k22 = dt*f2(t[k]+dt/2,y2[k]+k21/2,a1,b1,c1,a2,b2,c2,y1[k])
        
        k23 = dt*f2(t[k]+dt/2,y2[k]+k22/2,a1,b1,c1,a2,b2,c2,y1[k])
        
        k24 = dt*f2(t[k]+dt,y2[k]+k23,a1,b1,c1,a2,b2,c2,y1[k])
        
        y1[k+1] = y1[k] + (1/6)*(k11+2*k12+2*k13+k14)
        y2[k+1] = y2[k] + (1/6)*(k21+2*k22+2*k23+k24)

    return y1,y2

dt = 1.0

t_inicial = 0.0
t_final = 10.0
N = int((t_final-t_inicial)/dt)
t = np.linspace(t_inicial,t_final,N+1)

a1 = 0.1
b1 = 8*10.0**(-7)
c1 = 10.0**(-6)

a2 = 0.1
b2 = 8*10.0**(-7)
c2 = 10.0**(-7)



y_euler_explicito1,y_euler_explicito2 = euler_explicito(N,t,dt,a1,b1,c1,a2,b2,c2)
y_rk41,y_rk42 = rk4(N,t,dt,a1,b1,c1,a2,b2,c2)
 
#Inicializando o Gráfico

graf = plt.figure()
ax = graf.add_subplot()
 
#Parâmetros do Gráfico

graf.suptitle('Evolução das duas populações em competição ao longo do tempo', fontsize = 10, fontweight = 'bold', usetex = False)
ax.set_ylabel(r'N(t)', fontsize = 12, usetex = False)
ax.set_xlabel(r't', fontsize = 12, usetex = False) 

blp = mpatches.Patch(color = 'b',label = 'Euler Explícito Primeira População')

bly = mpatches.Patch(color = 'y',label = 'Euler Explícito Segunda População')

blr = mpatches.Patch(color = 'r',label = 'Range-Kutta de 4º Ordem Primeira População')

blg = mpatches.Patch(color = 'g',label = 'Range-Kutta de 4º Ordem Segunda População')
plt.legend(handles = [blp,blg,blr,bly], fontsize = 8)
plt.plot(t,y_euler_explicito1,'b')
plt.plot(t,y_euler_explicito2,'y')
plt.plot(t,y_rk41,'r')
plt.plot(t,y_rk42,'g')
plt.show()