import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def cond_iniciais(Nx,Ny,dx,dy,T):
    x_aux = 0.0
    y_aux = 0.0
    for i in range(Nx+1):
        y_aux = 0.0
        for j in range(Ny+1):
            if (i == 2 or i == 3) and (j == 2 or j == 3):
                T[i,j] = 1.0               
            else:
                T[i,j] = 0.0
            y_aux+=dy
        x_aux+=dx                         
    return T
    
def eqcondcal(Nx,Ny,dx,dt,tempo_final,T):
    tempo = 0.0
    x_aux = 0.0
    y_aux = 0.0
    intervalos = np.linspace(0.0,1.0,Nx+1)
    T_x = np.linspace(0.0,1.0,Nx+1)
    T_y = np.linspace(0.0,1.0,Ny+1)
    X,Y = np.meshgrid(T_x,T_y)
    T_new = np.copy(T)
    while tempo <= tempo_final:
        for i in range(1,Nx):
            for j in range(1,Ny):
                if (i == 2 or i == 3) and (j == 2 or j == 3):
                    T_new[i,j] = 1.0               
                else:
                    T_new[i,j] = T[i,j] + (dt/(dx*dx))*(T[i+1,j]-2.0*T[i,j]+T[i-1,j]) + (dt/(dy*dy))*(T[i,j+1]-2.0*T[i,j]+T[i,j-1])
                    
        T = np.copy(T_new)
        tempo+=dt

    ax.dist = 11
    ax.plot_surface(X,Y,np.transpose(T), cmap=cm.coolwarm,linewidth=0, antialiased=False)

Nx = 5
Ny = 5
dx = 1.0/Nx
dy = 1.0/Ny
dt = (dx*dx)/4.0
tempo_final = 0.1

T = np.zeros((Nx+1,Nx+1),float)
T = cond_iniciais(Nx,Ny,dx,dy,T)
tempo = 0.0
tempo_final = 10.0

#Inicializando o Grafico

graf = plt.figure()
ax = plt.axes(projection = '3d')

#Parâmetros do Grafico

graf.suptitle('Resolução da Equação do calor para uma Temperatura intermediária constante', fontsize = 10, fontweight = 'bold', usetex = False)
ax.set_ylabel(r'y', fontsize = 12, usetex = False)
ax.set_xlabel(r'x', fontsize = 12, usetex = False)
    
#Calculando a Equacao da Conducao do Calor
eqcondcal(Nx,Ny,dx,dt,tempo_final,T)
plt.show() #Inicializando o Grafico