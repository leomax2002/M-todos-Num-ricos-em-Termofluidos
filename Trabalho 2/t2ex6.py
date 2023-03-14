import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

#funcao desejada
def f(t,y):
    return 2.0-2.0*t+4*t**2-4*t**3-4*t**4

#funcao obtida a partir de manipulacoes algebricas para euler implicito
def f_implicito(t,y,dt):
    return y + dt*(2.0-2.0*t+4.0*t**2-4.0*t**3-4*t**4)

#funcao exata desejada
def f_exata(t):
    return 1.0+2.0*t-t**2.0+(4.0/3.0)*t**3-t**4-(4/5)*t**5

#funcao para o calculo de y pelo metodo de euler explicito
def euler_explicito(N,t,dt):
    y = np.zeros(N+1)
    y[0] = 1.0
    for k in range(N):
        y[k+1] = y[k] + dt*f(t[k], y[k])
    return y

#funcao para o calculo de y pelo metodo de euler implicito
def euler_implicito(N,t,dt):
    y = np.zeros(N+1)
    y[0] = 1.0
    for k in range(N):
        y[k+1] = f_implicito(t[k+1],y[k],dt)
    return y

#funcao para o calculo de y pelo metodo de range-kutta de segunda ordem    
def rk2(N,t,dt):
    y = np.zeros(N+1)
    y[0] = 1.0
    for k in range(N):
          
        a = dt*f(t[k],y[k])
        
        b = dt*f(t[k]+dt,y[k]+a)
        
        y[k+1] = y[k] + 0.5*(a+b)
    return y

#funcao para o calculo de y pelo metodo de range-kutta de quarta ordem
def rk4(N,t,dt):
    y = np.zeros(N+1)
    y[0] = 1.0
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

dt = np.array([0.001,0.002,0.0025,0.005,0.01,0.02,0.025,0.05,0.1,0.2,0.25,0.5])
y_erro_euler_explicito = np.zeros(len(dt))
y_erro_euler_implicito = np.zeros(len(dt))
y_erro_rk2 = np.zeros(len(dt))
y_erro_rk4 = np.zeros(len(dt))

t_inicial = 0.0
t_final = 1.0 

for i in range(len(dt)):
    N = int((t_final-t_inicial)/dt[i])
    t = np.linspace(t_inicial,t_final,N+1)
    y_exata = exata(N,t)
    y_euler_explicito = euler_explicito(N,t,dt[i])
    y_euler_implicito = euler_implicito(N,t,dt[i])
    y_rk2 = rk2(N,t,dt[i])
    y_rk4 = rk4(N,t,dt[i])

    y_erro_euler_explicito[i] = np.abs(y_euler_explicito[len(y_euler_explicito)-1] -y_exata[len(t)-1])
    y_erro_euler_implicito[i] = np.abs(y_euler_implicito[len(y_euler_implicito)-1] -y_exata[len(t)-1])
    y_erro_rk2[i] = np.abs(y_rk2[len(y_rk2)-1] -y_exata[len(t)-1])
    y_erro_rk4[i] = np.abs(y_rk4[len(y_rk4)-1] -y_exata[len(t)-1])

#Inicializando o Gráfico

graf = plt.figure()
ax = graf.add_subplot()
 
#Parâmetros do Gráfico

graf.suptitle('Comparação do erro entre os 4 métodos numéricos utilizados', fontsize = 10, fontweight = 'bold', usetex = False)
ax.set_ylabel(r'y', fontsize = 12, usetex = False)
ax.set_xlabel(r'dt', fontsize = 12, usetex = False)

blp = mpatches.Patch(color = 'b',label = 'Erro Euler Explícito')

blg = mpatches.Patch(color = 'g',label = 'Erro Euler Implícito')

blr = mpatches.Patch(color = 'r',label = 'Erro Range-Kutta de 2º Ordem')

bly = mpatches.Patch(color = 'y',label = 'Erro Range-Kutta de 4º Ordem')

plt.legend(handles = [blp,blg,blr,bly], fontsize = 8)
    
plt.plot(dt,y_erro_euler_explicito,'b')
plt.plot(dt,y_erro_euler_implicito,'g')
plt.plot(dt,y_erro_rk2,'r')
plt.plot(dt,y_erro_rk4,'y')
plt.yscale('log')
plt.xscale('log')
plt.show()