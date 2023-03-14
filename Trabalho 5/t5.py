import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from numba import njit

###########Calcular u_star##############
@njit
def u_star(u,v,Nx,Ny,dx,dy,dt,u_str,U,C1):
    
    for i in range(1,Nx):
        for j in range(0,Ny):
            C1 = 0.25*(v[i,j+1]+v[i-1,j+1]+v[i,j]+v[i-1,j])
            
            R = -dt*(u[i,j]*((u[i+1,j]-u[i-1,j])/(2*dx)) + C1*((u[i,j+1]-u[i,j-1])/(2*dy))) + (dt/Re)*((u[i+1,j]-2*u[i,j]+u[i-1,j])/(dx*dx) + (u[i,j+1]-2*u[i,j]+u[i,j-1])/(dy*dy))
            
            u_str[i,j] = u[i,j] + R
        
    for j in range(0,Ny):
        u_str[0,j] = 0.0
        u_str[Nx,j] = 0.0
    
    for i in range(0,Nx+1):
        u_str[i,-1] = -u_str[i,0]
        u_str[i,Ny] = 2.0*U[i] - u_str[i,Ny-1]
    return u_str

###########Calcular v_star##############
@njit
def v_star(u,v,Nx,Ny,dx,dy,dt,v_str,U,C2,R):
    for i in range(Nx):
        for j in range(1,Ny):
            C2 = 0.25*(u[i+1,j]+u[i,j]+u[i+1,j-1]+u[i,j-1])
            R = -dt*(C2*((v[i+1,j]-v[i-1,j])/(2*dx)) + v[i,j]*((v[i,j+1]-v[i,j-1])/(2*dy))) + (dt/Re)*(((v[i+1,j]-2*v[i,j]+v[i-1,j])/(dx*dx)) + ((v[i,j+1]-2*v[i,j]+v[i,j-1])/(dy*dy)))
            
            v_str[i,j] = v[i,j] + R
            
    for j in range(0,Ny+1):
        v_str[-1,j] = -v_str[0,j]
        v_str[Nx,j] = -v_str[Nx-1,j]
    
    for i in range(0,Nx):
        v_str[i,0] = 0.0
        v_str[i,Ny] = 0.0
    return v_str

###########Calcular a pressão##############
@njit
def calc_p(u_str,v_str,p,Nx,Ny,dx,dy,dt,R,tol):
    e = 100.0
    while e > tol:
        R_max = 0.0
        for i in range(Nx):
            for j in range(Ny):
                if i == 0 and j == 0:
                    lamda = -(((1)/(dx*dx)) + ((1)/(dy*dy)))
                    R = ((u_str[i+1,j]-u_str[i,j])/(dt*dx)) + ((v_str[i,j+1] - v_str[i,j])/(dt*dy)) -(((p[i+1,j]-p[i,j])/(dx*dx)) + ((p[i,j+1]-p[i,j])/(dy*dy)))
                    
                elif i == 0 and j == Ny-1:
                    lamda = -(((1)/(dx*dx)) + ((1)/(dy*dy)))
                    R = ((u_str[i+1,j]-u_str[i,j])/(dt*dx)) + ((v_str[i,j+1] - v_str[i,j])/(dt*dy)) -(((p[i+1,j]-p[i,j])/(dx*dx)) + ((-1.0*p[i,j]+p[i,j-1])/(dy*dy)))
                    
                elif i == Nx-1 and j == 0:
                    lamda = -(((1)/(dx*dx)) + ((1)/(dy*dy)))
                    R = ((u_str[i+1,j]-u_str[i,j])/(dt*dx)) + ((v_str[i,j+1]-v_str[i,j])/(dt*dy)) - (((-1.0*p[i,j] + p[i-1,j])/(dx*dx)) + ((p[i,j+1]-p[i,j])/(dy*dy)))
                
                elif i == Nx-1 and j == Ny-1:
                    lamda = -(((1)/(dx*dx)) + ((1)/(dy*dy)))
                    R = ((u_str[i+1,j]-u_str[i,j])/(dt*dx)) + ((v_str[i,j+1]-v_str[i,j])/(dt*dy)) - (((-1.0*p[i,j]+p[i-1,j])/(dx*dx)) + ((-1.0*p[i,j] + p[i,j-1])/(dy*dy)))
                    
                elif i == 0 and j != 0 and j != Ny-1:
                    lamda = -(((1)/(dx*dx)) + ((2)/(dy*dy)))
                    R = ((u_str[i+1,j]-u_str[i,j])/(dt*dx)) + ((v_str[i,j+1]-v_str[i,j])/(dt*dy)) - (((p[i+1,j] - p[i,j])/(dx*dx)) + ((p[i,j+1]-2*p[i,j]+p[i,j-1])/(dy*dy)))
                    
                elif i == Nx-1 and j != 0 and j != Ny-1:
                    lamda = -(((1)/(dx*dx)) + ((2)/(dy*dy)))
                    R = ((u_str[i+1,j]-u_str[i,j])/(dt*dx)) + ((v_str[i,j+1]-v_str[i,j])/(dt*dy)) - (((-1.0*p[i,j]+p[i-1,j])/(dx*dx)) + ((p[i,j+1]-2*p[i,j]+p[i,j-1])/(dy*dy)))
                    
                elif j == 0 and i != 0 and i != Nx-1:
                    lamda = -(((2)/(dx*dx)) + ((1)/(dy*dy)))
                    R = ((u_str[i+1,j]-u_str[i,j])/(dt*dx)) + ((v_str[i,j+1]-v_str[i,j])/(dt*dy)) - (((p[i+1,j]-2*p[i,j]+p[i-1,j])/(dx*dx)) + ((p[i,j+1]-p[i,j])/(dy*dy)))
                    
                elif j == Ny-1 and i != 0 and i != Nx-1:
                    lamda = -(((2)/(dx*dx)) + ((1)/(dy*dy)))
                    R = ((u_str[i+1,j] - u_str[i,j])/(dt*dx)) + ((v_str[i,j+1]-v_str[i,j])/(dt*dy)) - (((p[i+1,j] - 2*p[i,j] + p[i-1,j])/(dx*dx)) + ((-1.0*p[i,j] + p[i,j-1])/(dy*dy)))
                    
                else:
                    lamda = -(((2)/(dx*dx)) + ((2)/(dy*dy)))
                    R = ((u_str[i+1,j] - u_str[i,j])/(dt*dx)) + ((v_str[i,j+1]-v_str[i,j])/(dt*dy)) - (((p[i+1,j]-2*p[i,j] + p[i-1,j])/(dx*dx)) + ((p[i,j+1]-2*p[i,j]+p[i,j-1])/(dy*dy)))
                    
                R = R/lamda
                p[i,j] = p[i,j] + R
                
                if np.abs(R) > R_max:
                    R_max = np.abs(R)
                
                e = R_max
    
    for i in range(0,Nx):
        p[i,-1] = p[i,0]
        p[i,Ny] = p[i,Ny-1]
    
    for j in range(0,Ny):
        p[-1,j] = p[0,j]
        p[Nx,j] = p[Nx-1,j]
    
    p[-1,-1] = p[0,0]
    p[-1,Ny] = p[0,Ny-1]
    p[Nx,-1] = p[Nx-1,0]
    p[Nx,Ny] = p[Nx-1,Ny-1]
    
    return p
                    
###########Calcular novo u############## 
@njit
def u_new(u,u_str,p,Nx,Ny,dx,dt):
    for i in range(1,Nx):
        for j in range(-1,Ny+1):
            
            u[i,j] = u_str[i,j] - dt*((p[i,j] - p[i-1,j])/(dx))
            
    return u
                               
###########Calcular novo v############## 
@njit
def v_new(v,v_str,p,Nx,Ny,dy,dt):
    for i in range(-1,Nx+1):
        for j in range(1,Ny):
            
            v[i,j] = v_str[i,j] - dt*((p[i,j] - p[i,j-1])/(dy))
    
    return v

###########Calcular a funcao de fluxo##############
@njit
def calc_phi(phi,u,v,Nx,Ny,dx,dy,dt,tol):
    for i in range(0,Nx+1):
        for j in range(0,Ny+1):
            phi[i,j] = 0.0
    
    lamda = -((2)/(dx*dx) + (2)/(dy*dy))
    e = 100.0
    
    while e > tol:
        R_max = 0
        for i in range(1,Nx):
            for j in range(1,Ny):
                
                R = -((v[i,j]-v[i-1,j])/(dx)) + ((u[i,j]-u[i,j-1])/(dy)) - ((phi[i+1,j]-2*phi[i,j]+phi[i-1,j])/(dx*dx) + (phi[i,j+1]-2*phi[i,j]+phi[i,j-1])/(dy*dy))
                
                R = R/lamda
                phi[i,j] = phi[i,j] + R
                
                if np.abs(R) > R_max:
                    R_max = np.abs(R)
        e = R_max
    
    return phi
###########Definição de Variáveis##############

Nx = 10 #N = 10, Re = 1.0
Ny = 10

Re = 1.0
#dx = 0.5*(1/(np.sqrt(Re)))
#dy = 0.5*(1/(np.sqrt(Re)))
#dt = (0.125)*Re*(dx*dx)

dx = 1.0/Nx
dy = 2.0/Ny
#dt = (0.125)*Re*(dx*dx)
dt = 0.01

C1 = 0.0
C2 = 0.0
R = 0.0
tol = 1.e-7
U = np.ones(Nx+1, float)
u = np.zeros((Nx+1,Ny+2),float)
v = np.zeros((Nx+2,Ny+1),float)
p = np.zeros((Nx+2,Ny+2),float)

for i in range(0,Nx+1):
    u[i,Ny] = 2*U[i]


phi = np.zeros((Nx+1,Ny+1),float)


u_str = np.copy(u)
v_str = np.copy(v)
p_str = np.copy(p)

t = 0.0
t_final = 10.0 #t = 5.0 Re = 1.0; t = 30.0, outros valores de Re, tol = 1e-5

########### Algoritmo ##############

while t < t_final:
    t+=dt
    u_str =  u_star(u,v,Nx,Ny,dx,dy,dt,u_str,U,C1)
    v_str = v_star(u,v,Nx,Ny,dx,dy,dt,v_str,U,C2,R)
    p = calc_p(u_str,v_str,p,Nx,Ny,dx,dy,dt,R,tol)
    u = u_new(u,u_str,p,Nx,Ny,dx,dt)
    v = v_new(v,v_str,p,Nx,Ny,dy,dt)
    phi = calc_phi(phi,u,v,Nx,Ny,dx,dy,dt,tol)
    
phi = calc_phi(phi,u,v,Nx,Ny,dx,dy,dt,tol)

########### Parâmetros do Gráfico ##############
uplot = np.zeros((Nx+1,Ny+1),float)
vplot = np.zeros((Nx+1,Ny+1),float)
pplot = np.zeros((Nx+1,Ny+1),float)

uplot_unit = np.copy(uplot)
vplot_unit = np.copy(vplot)
for i in range(0,Nx+1):
    for j in range(0,Ny+1):
        uplot[i,j] = 0.5*(u[i,j]+u[i,j-1])
        vplot[i,j] = 0.5*(v[i,j] + v[i-1,j])
        pplot[i,j] = 0.25*(p[i,j] + p[i-1,j] + p[i,j-1] + p[i-1,j-1])

for i in range(0,Nx+1):
    for j in range(0,Ny+1):
        if i != 0 and j != 0 and i != Nx:
            modulo = np.sqrt(uplot[i,j]*uplot[i,j] + vplot[i,j]*vplot[i,j])
        
            if modulo != 0.0:
                uplot_unit[i,j] = uplot[i,j]/modulo
                vplot_unit[i,j] = vplot[i,j]/modulo

graf = plt.figure()
ax = graf.add_subplot()

graf.suptitle('Gráfico de Orientação do Fluxo da Cavidade em função de u e v', fontsize = 10, fontweight = 'bold', usetex = False)
ax.set_ylabel(r'y', fontsize = 12, usetex = False)
ax.set_xlabel(r'x', fontsize = 12, usetex = False)

x = np.linspace(0,1.0,Nx+1)
y= np.linspace(0,1.0,Ny+1)
plt.quiver(x,y,np.transpose(uplot_unit),np.transpose(vplot_unit))
plt.show()

graf = plt.figure()
ax = graf.add_subplot()

graf.suptitle('Gráfico de Variação da Pressão para a Cavidade', fontsize = 10, fontweight = 'bold', usetex = False)
ax.set_ylabel(r'y', fontsize = 12, usetex = False)
ax.set_xlabel(r'x', fontsize = 12, usetex = False)



plt.contourf(x,y,np.transpose(pplot))
plt.show()

graf = plt.figure()
ax = graf.add_subplot()

graf.suptitle('Gráfico de Fluxo de Corrente para a Cavidade', fontsize = 10, fontweight = 'bold', usetex = False)
ax.set_ylabel(r'y', fontsize = 12, usetex = False)
ax.set_xlabel(r'x', fontsize = 12, usetex = False)

plt.contour(x,y,np.transpose(phi), colors = 'k')
plt.show()
#print(vplot)
#print(pplot)
#print(pplot)

#for j in range(1,Ny):
#    print(f'{uplot[8,j]:.8e}')
