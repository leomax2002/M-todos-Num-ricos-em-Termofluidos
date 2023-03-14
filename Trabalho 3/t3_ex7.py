import numpy as np
import matplotlib.pyplot as plt

def f(t,x,y):
    return (np.sinh(np.pi*y)*np.sin(np.pi*x))/(np.sinh(np.pi))

def cond_iniciais(Nx,Ny,dx,dy,T):
    #Esquerda
    y_aux = 0.0
    for j in range(Ny+1):
        T[0,j] = 0.0
        y_aux+=dy
    #Direita por aproximacao de segunda ordem
    y_aux = 0.0
    for j in range(Ny+1):
        T[Nx,j] = (4/3)*T[Nx-1,j] -(1/3)*T[Nx-2,j]
        y_aux+=dy
    #Baixo
    x_aux = 0.0
    for i in range(Nx+1):
        T[i,0] = 0.0
        x_aux+=dx
    #Cima
    x_aux = 0.0
    for i in range(Nx+1):
        T[i,Ny] = np.sin(np.pi*x_aux)
        x_aux+=dx
    #quinas nao sao utilizadas para os calculos
    return T

def solexata(Nx,Ny,T,tempo_final):
    x_aux = 0.0
    y_aux = 0.0
    intervalos = np.linspace(0.0,1.0,Nx+1)
    T_x = np.linspace(0.0,1.0,Nx+1)
    T_y = np.linspace(0.0,1.0,Ny+1)
    for i in range(Nx+1):
        y_aux = 0.0
        for j in range(Ny+1):
            T[i,j] = f(tempo_final,x_aux,y_aux)
            y_aux+=dy
        x_aux+=dx
    c = plt.contour(T_x,T_y,np.transpose(T),intervalos, colors = 'k')
    plt.clabel(c, inline = True, fontsize = 15, colors = 'k')
    
def eqcondcal(Nx,Ny,dx,dt,tempo_final,T):
    tempo = 0.0
    intervalos = np.linspace(0.0,1.0,Nx+1)
    T_x = np.linspace(0.0,1.0,Nx+1)
    T_y = np.linspace(0.0,1.0,Ny+1)
    T_new = np.copy(T)
    while tempo <= tempo_final:

        for i in range(1,Nx+1):
            for j in range(1,Ny):
                if i == Nx:
                    T_new[i,j] = (4/3)*T[Nx-1,j] -(1/3)*T[Nx-2,j]
                else:
                    T_new[i,j] = T[i,j] + (dt/(dx*dx))*(T[i+1,j]-2.0*T[i,j]+T[i-1,j]) + (dt/(dy*dy))*(T[i,j+1]-2.0*T[i,j]+T[i,j-1])
                    
        T = np.copy(T_new)
        tempo+=dt

    c = plt.contour(T_x,T_y,np.transpose(T),intervalos)
    plt.clabel(c, inline = True, fontsize = 15)

Nx = 5
Ny = 5
dx = 0.5/Nx
dy = 1.0/Ny
dt = (dx*dx)/4.0

T = np.zeros((Nx+1,Nx+1),float)
T_exata = np.zeros((Nx+1,Nx+1),float)
T = cond_iniciais(Nx,Ny,dx,dy,T)

tempo = 0.0
tempo_final = 10.0

#Inicializando o Grafico

graf = plt.figure()
ax = graf.add_subplot()

#Parâmetros do Grafico

graf.suptitle('Comparação da Resolução da Equação do calor analítica e numérica', fontsize = 10, fontweight = 'bold', usetex = False)
ax.set_ylabel(r'y', fontsize = 12, usetex = False)
ax.set_xlabel(r'x', fontsize = 12, usetex = False)
    
#Calculando a Equacao da Conducao do Calor
    
solexata(Nx,Ny,T_exata,tempo_final)
eqcondcal(Nx,Ny,dx,dt,tempo_final,T)
plt.show() #Inicializando o Grafico