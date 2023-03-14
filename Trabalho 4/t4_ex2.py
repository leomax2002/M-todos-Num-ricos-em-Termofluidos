import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import time
from numba import njit

@njit
def solexata(Nx,Ny,T):
    for i in range(Nx+1):
        for j in range(Ny+1):
            T[i,j] = np.sin(np.pi*i*dx)*np.sin(np.pi*j*dy)
    return T

@njit
def difcent(Nx,Ny,dx,dy,T,T_new):
    iter = 0
    while True:
        #print(np.abs(R))
        for i in range(1,Nx):
            for j in range(1,Ny):
                T_new[i,j] = 0.25*(T_new[i+1,j] + T_new[i-1,j] + T_new[i,j+1] + T_new[i,j-1]
                +(2*(np.pi*np.pi)*np.sin(np.pi*i*dx)*np.sin(np.pi*j*dy)*(dx*dx)))
        iter+=1
        if np.max(np.abs(T_new)-np.abs(T)) < 1e-8:
            break
        T = np.copy(T_new)
    print('Difcent', iter)
    return T
@njit
def sor(Nx,Ny,dx,dy,R,T):
    iter = 0
    while True:
        #print(np.abs(R))
        for i in range(1,Nx):
            for j in range(1,Ny):
                #R[i,j] = 0.25*(T[i+1,j] + T[i-1,j] -4*T[i,j] + T[i,j+1] + T[i,j-1]
                #+(2*(np.pi*np.pi)*np.sin(np.pi*i*dx)*np.sin(np.pi*j*dy)*(dx*dx)))
                #R[i,j] = ((dx*dx*dy*dy)/(2*(dx*dx+dy*dy)))*(((T[i+1,j] + T[i-1,j])/(dx*dx)) 
                #-((2*(dy*dy+dx*dx))/(dx*dx*dy*dy))*T[i,j] + ((T[i,j+1] + T[i,j-1])/(dy*dy))
                #+(2*(np.pi*np.pi)*np.sin(np.pi*i*dx)*np.sin(np.pi*j*dy)))
                R[i,j] = ((dx*dx*dy*dy)/(2*(dx*dx+dy*dy)))*(((T[i+1,j] + T[i-1,j])/(dx*dx)) 
                + ((T[i,j+1] + T[i,j-1])/(dy*dy))
                +(2*(np.pi*np.pi)*np.sin(np.pi*i*dx)*np.sin(np.pi*j*dy))) - T[i,j]
                T[i,j] = T[i,j] + w*R[i,j]
        iter+=1
        if np.max(np.abs(R)) < 1e-8:
            break
    print('SOR:',iter)
    return T

Nx = 200 #10,100,200 primeiro difcent depois sor
Ny = 200 #10,100,200

dx = 1.0/Nx
dy = 1.0/Ny

w = 1.4

T = np.zeros((Nx+1,Ny+1),float)
T_aux = np.copy(T)
T_aux2 = np.copy(T)
T_new = np.copy(T)

R = np.zeros((Nx+1,Ny+1),float)
R[2,2] = 1.0

T_ex = solexata(Nx,Ny,T_aux)
t1 = time.time()
#T1 = sor(Nx,Ny,dx,dy,R,T)
T2 = difcent(Nx,Ny,dx,dy,T_aux2,T_new)
t2 = time.time()
print(t2-t1)


#Inicializando o Grafico
graf = plt.figure()
#ax = graf.add_subplot()
ax = plt.axes(projection = '3d')

intervalos = np.linspace(0.0,1.0,Nx+1)
T_x = np.linspace(0.0,1.0,Nx+1)
T_y = np.linspace(0.0,1.0,Ny+1)
X,Y = np.meshgrid(T_x,T_y)



#Parâmetros do Grafico
#intervalos = np.linspace(0.0,1.0,Nx+1)
#T_x1 = np.linspace(0.0,1.0,Nx+1)
#T_y1 = np.linspace(0.0,1.0,Ny+1)

#T_x2 = np.linspace(0.0,1.0,Nx+1)
#T_y2 = np.linspace(0.0,1.0,Ny+1)

#T_x3 = np.linspace(0.0,1.0,Nx+1)
#T_y3 = np.linspace(0.0,1.0,Ny+1)


graf.suptitle('Comparação da Resolução da Equação do calor bidimensional numérica e analítica', fontsize = 10, fontweight = 'bold', usetex = False)
ax.set_ylabel(r'y', fontsize = 12, usetex = False)
ax.set_xlabel(r'x', fontsize = 12, usetex = False)

#c = plt.contour(T_x1,T_y1,np.transpose(T_ex),intervalos, colors = 'k')
#plt.clabel(c, inline = False, fontsize = 0)

#c = plt.contour(T_x2,T_y2,np.transpose(T1),intervalos)
#plt.clabel(c, inline = True, fontsize = 15)

#c = plt.contour(T_x3,T_y3,np.transpose(T1),intervalos)
#plt.clabel(c, inline = False, fontsize = 0)

ax.dist = 11
ax.plot_surface(X,Y,np.transpose(T2), cmap=cm.coolwarm,linewidth=0, antialiased=False)

plt.show()