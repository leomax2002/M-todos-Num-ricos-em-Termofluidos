from this import d
import numpy as np

N = 5
dx = 1.0/N
dt = (dx**2)/4.0

T = np.zeros(N+2)

T[5] = 2.0
T_new = np.copy(T)
for n in range(4):
    for i in range(0,N):
        T_new[i] = T[i] + (dt/(dx*dx))*(T[i+1]-2.0*T[i]+T[i-1])
        
    
    T_new[-1] = -T_new[0]
    T_new[N] = 2.0 - T_new[N-1]
    T = np.copy(T_new)
T_plot = np.zeros(N+1)
x = np.linspace(0.0,1.0,N+1)
print(T)