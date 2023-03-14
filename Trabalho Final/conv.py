import numpy as np
import matplotlib.pyplot as plt

N = 20
L = 5.0
dx = L/N
dt = 0.001
u = 1.0
T = np.zeros(N+1)
x = np.linspace(0.0,1.0,N+1)

for i in range(int(N/2)-10,int(N/2)+10):
    T[i] = 1.0

T_new = np.copy(T)
for n in range(10):
    for i in range(1,N+1):
        T_new[i] = T[i] - (u*dt/dx)*(T[i]-T[i-1])
    
    T = np.copy(T_new)    
plt.plot(x,T)
plt.show()