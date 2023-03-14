import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#funcao desejada
def f(t,y):
    return 1.0 + y/t

#funcao obtida a partir de manipulacoes algebricas para euler implicito
def f_implicito(t,y,dt):
    return (y + dt)/(1-dt/t)

#funcao exata desejada
def f_exata(t):
    return t*np.log(t) + 2*t

#funcao para o calculo de y pelo metodo de euler explicito
def euler_explicito(N,t,dt):
    y = np.zeros(N+1)
    y[0] = 2.0
    for k in range(N):
        y[k+1] = y[k] + dt*f(t[k], y[k])
    return y

#funcao para o calculo de y pelo metodo de euler implicito
def euler_implicito(N,t,dt):
    y = np.zeros(N+1)
    y[0] = 2.0
    for k in range(N):
        y[k+1] = f_implicito(t[k+1],y[k],dt)
    return y

#funcao para o calculo de y pelo metodo de range-kutta de segunda ordem    
def rk2(N,t,dt):
    y = np.zeros(N+1)
    y[0] = 2.0
    for k in range(N):
          
        a = dt*f(t[k],y[k])
        
        b = dt*f(t[k]+dt,y[k]+a)
        
        y[k+1] = y[k] + 0.5*(a+b)
    return y

#funcao para o calculo de y pelo metodo de range-kutta de quarta ordem
def rk4(N,t,dt):
    y = np.zeros(N+1)
    y[0] = 2.0
    for k in range(N):
          
        a = dt*f(t[k],y[k])
        
        b = dt*f(t[k]+dt/2,y[k]+a/2)
        
        c = dt*f(t[k]+dt/2,y[k]+b/2)
        
        d = dt*f(t[k]+dt,y[k]+c)
        
        y[k+1] = y[k] + (1/6)*(a+2*b+2*c+d)

    return y

def exata(N,t):
    y_exata = np.zeros(N+1)   
    for k in range(N+1):
        y_exata[k] = f_exata(t[k])
    return y_exata

dt = 1.0

t_inicial = 1.0
t_final = 2.0
N = int((t_final-t_inicial)/dt)
t = np.linspace(t_inicial,t_final,N+1)

y_euler_explicito = euler_explicito(N,t,dt)
y_euler_implicito = euler_implicito(N,t,dt)
y_rk2 = rk2(N,t,dt)
y_rk4 = rk4(N,t,dt)
y_exata = exata(N,t)
 
#Inicializando o Gráfico

graf = plt.figure()
ax = graf.add_subplot()
 
#Parâmetros do Gráfico

graf.suptitle('Comparação dos resultados obtidos para os 4 métodos numéricos utilizados', fontsize = 10, fontweight = 'bold', usetex = False)
ax.set_ylabel(r'y', fontsize = 12, usetex = False)
ax.set_xlabel(r't', fontsize = 12, usetex = False)
blk = mpatches.Patch(color = 'k',label = 'Solução Analítica')

blp = mpatches.Patch(color = 'b',label = 'Euler Explícito')

blg = mpatches.Patch(color = 'g',label = 'Euler Implícito')

blr = mpatches.Patch(color = 'r',label = 'Range-Kutta de 2º Ordem')

bly = mpatches.Patch(color = 'y',label = 'Range-Kutta de 4º Ordem')
plt.legend(handles = [blk,blp,blg,blr,bly], fontsize = 8)

#Traçando os Gráficos
    
plt.plot(t,y_euler_explicito,'b')
plt.plot(t,y_euler_implicito,'g')
plt.plot(t,y_rk2,'r')
plt.plot(t,y_rk4,'y')
plt.plot(t,y_exata,'k')
plt.show()