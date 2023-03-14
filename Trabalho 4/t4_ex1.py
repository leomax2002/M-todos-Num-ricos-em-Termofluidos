import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from numba import njit

@njit
def eqcondcal(Nx,Ny,dx,dt,T,T_new):
    iter = 0
    while True:
        for i in range(1,Nx):
            for j in range(1,Ny):
                T_new[i,j] = T[i,j] + (dt/(dx*dx))*(T[i+1,j]-2.0*T[i,j]+T[i-1,j]) + (dt/(dy*dy))*(T[i,j+1]-2.0*T[i,j]+T[i,j-1])
        iter+=1
        if np.abs(np.max(T_new[:len(T_new)-1:,:len(T_new)-1:])-np.abs(np.max(T[:len(T)-1:,:len(T_new)-1:]))) < 1e-6:
            break
        T = np.copy(T_new)
    print(iter)
    return T
@njit
def sor(Nx,Ny,dx,dy,R,T,w):
    R[2,2] = 1.0
    iter = 0
    while np.max(np.abs(R)) > 1e-6:
        for i in range(1,Nx):
            for j in range(1,Ny):
                R[i,j] = 0.25*(T[i+1,j] + T[i-1,j]+T[i,j+1]+T[i,j-1]-4*T[i,j])
                T[i,j] = T[i,j] + w*R[i,j] 
        iter+=1
    print(iter)
    return T
Nx = 100
Ny = 100
dx = 1.0/Nx
dy = 1.0/Ny
w = 1.95
dt = (dx*dx)/4.0
T = np.zeros((Nx+1,Ny+1),float)
T[:,Ny] = 1.0
T_new = np.copy(T)
T_aux = np.copy(T)
R = np.zeros((Nx+1,Ny+1),float)

t1 = time.time()
#T1 = sor(Nx,Ny,dx,dy,R,T,w)

T2 = eqcondcal(Nx,Ny,dx,dt,T_aux,T_new)
t2 = time.time()

print(t2-t1)


#Inicializando o Grafico
graf = plt.figure()
#ax = graf.add_subplot()
ax = plt.axes(projection = '3d')

#Parâmetros do Grafico
intervalos = np.linspace(0.0,1.0,Nx+1)
T_x = np.linspace(0.0,1.0,Nx+1)
T_y = np.linspace(0.0,1.0,Ny+1)

intervalos2 = np.linspace(0.0,1.0,Nx+1)
T_x2 = np.linspace(0.0,1.0,Nx+1)
T_y2 = np.linspace(0.0,1.0,Ny+1)
X,Y = np.meshgrid(T_x,T_y)

graf.suptitle('Comparação da Resolução da Equação do calor para diferenças finitas e SOR', fontsize = 10, fontweight = 'bold', usetex = False)
ax.set_ylabel(r'y', fontsize = 12, usetex = False)
ax.set_xlabel(r'x', fontsize = 12, usetex = False)

#c2 = plt.contour(T_x2,T_y2,np.transpose(T1),intervalos, colors = 'k')
#plt.clabel(c2, inline = True, fontsize = 0)

#c = plt.contour(T_x,T_y,np.transpose(T1),intervalos)
#plt.clabel(c, inline = False, fontsize = 0)

ax.dist = 11
ax.plot_surface(X,Y,np.transpose(T2), cmap=cm.coolwarm,linewidth=0, antialiased=False)
#ax.plot_surface(X,Y,np.transpose(T2), cmap=cm.coolwarm,linewidth=0, antialiased=False)

plt.show()