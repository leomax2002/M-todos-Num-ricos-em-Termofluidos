import numpy as np
import matplotlib.pyplot as plt
from numba import njit

@njit
def f(t,x,y):
    return (np.sinh(np.pi*y)*np.sin(np.pi*x))/(np.sinh(np.pi))

@njit
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

@njit
def solexata(Nx,Ny,T,tempo_final):
    x_aux = 0.0
    y_aux = 0.0
    for i in range(Nx+1):
        y_aux = 0.0
        for j in range(Ny+1):
            T[i,j] = f(tempo_final,x_aux,y_aux)
            y_aux+=dy
        x_aux+=dx
    return T
 
@njit   
def btcs(Nx,Ny,ld1,ld2,R,dt,tempo_final,T,T_new,w):
    tempo = 0.0
    while tempo < tempo_final:
        R[2,2] = 1.0
        while np.max(np.abs(R)) > 1e-6:
            for i in range(1,Nx+1):
                for j in range(1,Ny):
                    d = 1.0+2.0*ld1+2.0*ld2
                    if i == Nx:
                        T_new[i,j] = (4/3)*(T[i-1,j]) - (1/3)*(T[i-2,j])
                    else:
                        R[i,j] = (1.0/d)*(T[i,j] + ld1*(T_new[i+1,j]+T_new[i-1,j]) + ld2*(T_new[i,j+1]+T_new[i,j-1])) - T_new[i,j]
                        T_new[i,j] = T_new[i,j] + w*R[i,j]
        T = np.copy(T_new)
        tempo+=dt
    return T

@njit
def cranknicolson(Nx,Ny,ld1,ld2,R,dt,tempo_final,T,T_new,w):
    tempo = 0.0
    while tempo < tempo_final:
        R[2,2] = 1.0
        while np.max(np.abs(R)) > 1e-6:
            for i in range(1,Nx+1):
                for j in range(1,Ny):
                    d = 1.0 + ld1 + ld2
                    if i == Nx:
                        T_new[i,j] = (4/3)*(T[i-1,j]) - (1/3)*(T[i-2,j])
                    else:
                        R[i,j] = ((1-ld1-ld2)/d)*T[i,j] + (0.5/(d))*(ld1*(T[i+1,j]+T[i-1,j]+T_new[i+1,j]+T_new[i-1,j]) + ld2*(T[i,j+1]+T[i,j-1]+T_new[i,j+1]+T_new[i,j-1])) - T_new[i,j]
                        T_new[i,j] = T_new[i,j] + w*R[i,j]
        T = np.copy(T_new)
        tempo+=dt
    return T


Nx = 200
Ny = 200
dx = 0.5/Nx
dy = 1.0/Ny
tempo_final = 5.0
dt = 0.0001*tempo_final
ld1 = dt/(dx*dx)
ld2 = dt/(dy*dy)
w = 1.0

T = np.zeros((Nx+1,Ny+1),float)
T_exata = np.zeros((Nx+1,Ny+1),float)
T_aux = np.copy(T)
R = np.zeros((Nx+1,Ny+1),float)
T = cond_iniciais(Nx,Ny,dx,dy,T)
T_new = np.copy(T)
T_new_aux = np.copy(T_aux)

#Inicializando o Grafico

graf = plt.figure()
ax = graf.add_subplot()

#Parâmetros do Grafico

graf.suptitle('Comparação da Resolução da Equação do calor analítica e numérica', fontsize = 10, fontweight = 'bold', usetex = False)
ax.set_ylabel(r'erro', fontsize = 12, usetex = False)
ax.set_xlabel(r'T', fontsize = 12, usetex = False)
plt.xscale('log')
plt.yscale('log')  
    
#Calculando a Equacao da Conducao do Calor
T_aux = np.copy(T)
T_n = np.copy(T_new)
T_ex = solexata(Nx,Ny,T_exata,tempo_final)
T1 = btcs(Nx,Ny,ld1,ld2,R,dt,tempo_final,T,T_new,w)
T2 = cranknicolson(Nx,Ny,ld1,ld2,R,dt,tempo_final,T_aux,T_n,w)

erro_btcs = np.zeros(Nx+1,float)
erro_cn = np.zeros(Nx+1,float)
for i in range(Nx):
    erro_btcs[i] = np.abs(T1[i,i]-T_ex[i,i])
    erro_cn[i] = np.abs(T2[i,i]-T_ex[i,i])

print(erro_btcs)
print(erro_cn)

#intervalos = np.linspace(0.0,1.0,Nx+1)
#T_x = np.linspace(0.0,1.0,Nx+1)
#T_y = np.linspace(0.0,1.0,Ny+1)

#c = plt.contour(T_x,T_y,np.transpose(T_ex),intervalos, colors = 'k')
#plt.clabel(c, inline = False, fontsize = 0, colors = 'k')
#c1 = plt.contour(T_x,T_y,np.transpose(T1),intervalos, colors = 'b')
#plt.clabel(c1, inline = False, fontsize = 0, colors = 'k')
#c2 = plt.contour(T_x,T_y,np.transpose(T2),intervalos, colors = 'r')
#plt.clabel(c2, inline = False, fontsize = 0, colors = 'k')
    
plt.plot(erro_btcs, 'b')
plt.plot(erro_cn, 'r')
plt.show() #Inicializando o Grafico