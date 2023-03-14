import numpy as np
import matplotlib.pyplot as plt
from numba import njit

@njit
def calor(T,T_new,N,dx,dt,t):
    while t < 0.2:
        for i in range(N):
            for j in range(N):
                T_new[i,j] = T[i,j] + (dt/(dx**2.0))*(T[i+1,j]+T[i-1,j]+T[i,j+1]+T[i,j-1]-4.0*T[i,j])

        for i in range(N):
            T_new[i,-1] = -T_new[i,0]
            T_new[i,N] = 2.0 -T_new[i,N-1]
            
        for j in range(N):
            T_new[-1,j] = T_new[0,j]
            T_new[N,j] = T_new[N-1,j]
            
        T = np.copy(T_new)     
        t+=dt
    return T
        
    

N = 10
dx = 1.0/N
dy = 1.0/N
dt = 0.1*(dx**2.0)
T = np.zeros((N+2,N+2),float)
T_new = np.copy(T)



for i in range(N):
    T[i,N] = 2.0
   
t_final = 0.1
t = 0.0 
#while t < 5.0:
 #   for i in range(N):
  #      for j in range(N):
   #         T_new[i,j] = T[i,j] + (dt/(dx**2.0))*(T[i+1,j]+T[i-1,j]+T[i,j+1]+T[i,j-1]-4.0*T[i,j])

    #for i in range(N):
     #   T_new[i,-1] = -T_new[i,0]
     #   T_new[i,N] = 2.0 -T_new[i,N-1]
        
    #for j in range(N):
     #   T_new[-1,j] = -T_new[0,j]
     #   T_new[N,j] = -T_new[N-1,j]
        
    #T = np.copy(T_new)     
    #t+=dt   
    
T_plot = np.zeros((N+1,N+1),float)
T = calor(T,T_new,N,dx,dt,t)
for i in range(0,N+1):
    for j in range(1,N):
        T_plot[i,j] = 0.25*(T[i-1,j-1] + T[i,j-1] + T[i-1,j] + T[i,j])

for i in range(1,N):
    T_plot[i,N] = 1.0    
print(T_plot)


intervalos = np.linspace(0.0,1.0,6)
T_x = np.linspace(0.0,1.0,N+1)
T_y = np.linspace(0.0,1.0,N+1)
c = plt.contour(T_x,T_y,np.transpose(T_plot),intervalos)
plt.clabel(c, inline = True, fontsize = 15)

plt.show()