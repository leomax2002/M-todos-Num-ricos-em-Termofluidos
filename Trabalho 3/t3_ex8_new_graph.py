import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def f(t,x,y):
    res = 0.0
    for n in range(1,11):
        res+=((np.sin(2*n*np.pi/3))/(n*n*np.sinh(n*np.pi)))*np.sin(n*np.pi*x)*np.sinh(n*np.pi*y)
    return (450/(np.pi*np.pi))*res

def cond_iniciais(Nx,Ny,dx,dy,T):
    #Esquerda
    y_aux = 0.0
    for j in range(Ny+1):
        T[0,j] = 0.0
        y_aux+=dy
    #Direita
    y_aux = 0.0
    for j in range(Ny+1):
        T[Nx,j] = 0.0
        y_aux+=dy
    #Baixo
    x_aux = 0.0
    for i in range(Nx+1):
        T[i,0] = 0.0
        x_aux+=dx
    #Cima
    x_aux = 0.0
    for i in range(Nx+1):
        if x_aux <= (2/3):
            T[i,Ny] = 75*x_aux
        else:
            T[i,Ny] = 150*(1-x_aux)
        x_aux+=dx
    #quinas nao sao utilizadas para os calculos
    return T

def solexata(Nx,Ny,T,tempo_final):
    x_aux = 0.0
    y_aux = 0.0
    intervalos = np.linspace(0.0,1.0,Nx+1)
    T_x = np.linspace(0.0,1.0,Nx+1)
    T_y = np.linspace(0.0,1.0,Ny+1)
    X,Y = np.meshgrid(T_x,T_y)
    for i in range(Nx+1):
        y_aux = 0.0
        for j in range(Ny+1):
            T[i,j] = f(tempo_final,x_aux,y_aux)
            y_aux+=dy
        x_aux+=dx
    ax.dist = 11
    ax.plot_surface(X,Y,np.transpose(T), cmap=cm.coolwarm,linewidth=0, antialiased=False)
    
def eqcondcal(Nx,Ny,dx,dt,tempo_final,T):
    tempo = 0.0
    intervalos = np.linspace(0.0,1.0,Nx+1)
    T_x = np.linspace(0.0,1.0,Nx+1)
    T_y = np.linspace(0.0,1.0,Ny+1)
    X,Y = np.meshgrid(T_x,T_y)
    T_new = np.copy(T)
    while tempo <= tempo_final:

        for i in range(1,Nx):
            for j in range(1,Ny):
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
tempo_final = 10.0

T = np.zeros((Nx+1,Nx+1),float)
T_exata = np.zeros((Nx+1,Nx+1),float)
T = cond_iniciais(Nx,Ny,dx,dy,T)

#Inicializando o Grafico

graf = plt.figure()
ax = plt.axes(projection = '3d')

#Parâmetros do Grafico

graf.suptitle('Comparação da Resolução da Equação do calor analítica e numérica', fontsize = 10, fontweight = 'bold', usetex = False)
ax.set_ylabel(r'y', fontsize = 12, usetex = False)
ax.set_xlabel(r'x', fontsize = 12, usetex = False)
    
#Calculando a Equacao da Conducao do Calor
    
#solexata(Nx,Ny,T_exata,tempo_final)
eqcondcal(Nx,Ny,dx,dt,tempo_final,T)
plt.show() #Inicializando o Grafico