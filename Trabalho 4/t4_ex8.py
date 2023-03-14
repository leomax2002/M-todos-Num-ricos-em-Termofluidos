import numpy as np
import matplotlib.pyplot as plt

def eqcondcal2(Nx,Ny,dx,dt,tempo_final,T):
    tempo = 0.0
    intervalos = np.linspace(0.0,1.0,Nx+1)
    T_x = np.linspace(0.0,1.0,Nx+1)
    T_y = np.linspace(0.0,1.0,Ny+1)
    T_new = np.copy(T)
    while tempo <= tempo_final:

        for i in range(1,Nx):
            for j in range(1,Ny):
                T_new[i,j] = T[i,j] + (dt/(dx*dx))*(T[i+1,j]-2.0*T[i,j]+T[i-1,j]) + (dt/(dy*dy))*(T[i,j+1]-2.0*T[i,j]+T[i,j-1])
        T = np.copy(T_new)
        tempo+=dt

    c = plt.contour(T_x,T_y,np.transpose(T),intervalos, colors = 'k')
    plt.clabel(c, inline = True, fontsize = 15)


def cond_iniciais(Nx,Ny,dx,dy,T):
    
    for j in range(0,Ny):
        T[-1,j] = -T[0,j]
        T[Nx,j] = -T[Nx-1,j]
        #T[Nx,j] = 0.0
        
    
    for i in range(0,Nx+1):
        T[i,-1] = -T[i,0]
        #T[i,-1] = 0.0
        T[i,Ny] = 2-T[i,Ny-1]
    return T
   
def eqcondcal(Nx,Ny,dx,dt,tempo_final,T):
    tempo = 0.0
    T_new = np.copy(T)
    while tempo <= tempo_final:

        for i in range(1,Nx):
            for j in range(1,Ny):
                T_new[i,j] = T[i,j] + (dt/(dx*dx))*(T[i+1,j]-2.0*T[i,j]+T[i-1,j]) + (dt/(dy*dy))*(T[i,j+1]-2.0*T[i,j]+T[i,j-1])
        
        for j in range(0,Ny):
            T_new[-1,j] = -T_new[0,j]
            T_new[Nx,j] = -T_new[Nx-1,j]
            
        for i in range(0,Nx):
            T_new[i,-1] = -T_new[i,0]
            T_new[i,Ny] = 2-T_new[i,Ny-1]
        
        T = np.copy(T_new)
        tempo+=dt
        
    return T

Nx = 5
Ny = 5
dx = 1.0/Nx
dy = 1.0/Ny
dt = (dx*dx)/4.0
tempo_final = 0.1

T_approx = np.zeros((Nx+2,Ny+2),float)
T_approx = cond_iniciais(Nx,Ny,dx,dy,T_approx)

T = np.zeros((Nx+1,Ny+1),float)
T_exata = np.copy(T)

for i in range(0,Nx+1):
    T_exata[i,Ny] = 1.0

#Inicializando o Grafico

graf = plt.figure()
ax = graf.add_subplot()

#Parâmetros do Grafico

graf.suptitle('Comparação da Resolução da Equação do calor numérica com e em malha defasada', fontsize = 8, fontweight = 'bold', usetex = False)
ax.set_ylabel(r'y', fontsize = 12, usetex = False)
ax.set_xlabel(r'x', fontsize = 12, usetex = False)
    
#Calculando a Equacao da Conducao do Calor
    
T_approx = eqcondcal(Nx,Ny,dx,dt,tempo_final,T_approx)

for i in range(1,Nx+1):
    for j in range(1,Ny+1):
        T[i,j] = 0.25*(T_approx[i-1,j-1] + T_approx[i,j-1] + T_approx[i-1,j] + T_approx[i,j])
intervalos = np.linspace(0.0,1.0,Nx+1)
T_x = np.linspace(0.0,1.0,Nx+1)
T_y = np.linspace(0.0,1.0,Ny+1)
c = plt.contour(T_x,T_y,np.transpose(T),intervalos)
plt.clabel(c, inline = True, fontsize = 15)

#eqcondcal2(Nx,Ny,dx,dt,tempo_final,T_exata)
plt.show() #I.nicializando o Grafico