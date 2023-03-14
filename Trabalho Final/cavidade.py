import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

N_x = 16
N_y = 16
L_x = 1.0
L_y = 1.0
dx = L_x/N_x
dy = L_y/N_y

Re = 10.0

t_final = 0.001

dt = 0.001

tol = 1.e-5

u = np.zeros((N_x+1, N_y+2), float)
v = np.zeros((N_x+2, N_y+1), float)
p = np.zeros((N_x+2, N_y+2), float)

U = np.ones((N_x+1))

for i in range(0, N_x+1):
  u[i, N_y] = 2.0*U[i]

u_star = np.copy(u)
v_star = np.copy(v)

t = 0.0

while t < t_final: 

  t = t + dt

  for i in range(1, N_x):
    for j in range(0, N_y):
      C_1 = 0.25*(v[i, j+1] + v[i-1, j+1] + v[i, j] + v[i-1, j])

      R = -dt*( u[i, j]*(u[i+1, j] - u[i-1, j])/(2.0*dx) 
          + C_1*(u[i, j+1] - u[i, j-1])/(2.0*dy) )

      R += (dt/Re)*( (u[i+1, j] - 2.0*u[i, j] + u[i-1, j])/(dx**2.0) 
          + (u[i, j+1] - 2.0*u[i, j] + u[i, j-1])/(dy**2.0) )

      u_star[i, j] = u[i, j] + R

  for j in range(0, N_y):
    u_star[0, j] = 0.0
    u_star[N_x, j] = 0.0

  for i in range(0, N_x+1):
    u_star[i, -1] = -u_star[i, 0]
    u_star[i, N_y] = 2.0*U[i] - u_star[i, N_y-1]


  for i in range(0, N_x):
    for j in range(1, N_y):
      C_2 = 0.25*(u[i+1, j] + u[i, j] + u[i+1, j-1] + u[i, j-1])

      R = - dt*(C_2*(v[i+1, j] - v[i-1, j])/(2.0*dx) 
          + v[i, j]*(v[i, j+1] - v[i, j-1])/(2.0*dy))
      
      R = R + (dt/Re)*( (v[i+1, j] - 2.0*v[i, j] + v[i-1, j])/(dx**2.0)
          + (v[i, j+1] - 2.0*v[i, j] + v[i, j-1])/(dy**2.0) )

      v_star[i, j] = v[i, j] + R

  for j in range(0, N_y+1):
    v_star[-1, j] = - v_star[0, j]
    v_star[N_x, j] = - v_star[N_x - 1, j]

  for i in range(0, N_x):
    v_star[i, 0] = 0.0
    v_star[i, N_y] = 0.0

  error = 100.0

  iter = 0

  while error > tol:

    R_max = 0.0

    for i in range(0, N_x):
      for j in range(0, N_y):

        if i == 0 and j == 0:

          lamda = - (1.0/(dx**2.0) + 1.0/(dy**2.0))

          R = (u_star[i+1, j] - u_star[i, j])/(dt*dx) + (v_star[i, j+1] - v_star[i, j])/(dt*dy)

          R = R - ( (p[i+1, j] - 1.0*p[i, j] )/(dx**2.0) 
              + (p[i, j+1] - 1.0*p[i, j] )/(dy**2.0) )

        elif i == 0 and j == N_y-1:
          
          lamda = - (1.0/(dx**2.0) + 1.0/(dy**2.0))

          R = (u_star[i+1, j] - u_star[i, j])/(dt*dx) + (v_star[i, j+1] - v_star[i, j])/(dt*dy)

          R = R - ( (p[i+1, j] - 1.0*p[i, j])/(dx**2.0) 
              + (- 1.0*p[i, j] + p[i, j-1])/(dy**2.0) )        

        elif i == N_x-1 and j == 0:

          lamda = - (1.0/(dx**2.0) + 1.0/(dy**2.0))

          R = (u_star[i+1, j] - u_star[i, j])/(dt*dx) + (v_star[i, j+1] - v_star[i, j])/(dt*dy)

          R = R - ( ( -1.0*p[i, j] + p[i-1, j])/(dx**2.0) 
              + (p[i, j+1] - 1.0*p[i, j] )/(dy**2.0) )        


        elif i == N_x-1 and j == N_y-1:

          lamda = - (1.0/(dx**2.0) + 1.0/(dy**2.0))

          R = (u_star[i+1, j] - u_star[i, j])/(dt*dx) + (v_star[i, j+1] - v_star[i, j])/(dt*dy)

          R = R - ( ( - 1.0*p[i, j] + p[i-1, j])/(dx**2.0) 
              + (- 1.0*p[i, j] + p[i, j-1])/(dy**2.0) )


        elif i == 0:

          lamda = - (1.0/(dx**2.0) + 2.0/(dy**2.0))

          R = (u_star[i+1, j] - u_star[i, j])/(dt*dx) + (v_star[i, j+1] - v_star[i, j])/(dt*dy)

          R = R - ( (p[i+1, j] - 1.0*p[i, j])/(dx**2.0) 
              + (p[i, j+1] - 2.0*p[i, j] + p[i, j-1])/(dy**2.0) )     

        elif i == N_x - 1:

          lamda = - (1.0/(dx**2.0) + 2.0/(dy**2.0))

          R = (u_star[i+1, j] - u_star[i, j])/(dt*dx) + (v_star[i, j+1] - v_star[i, j])/(dt*dy)

          R = R - ( (- 1.0*p[i, j] + p[i-1, j])/(dx**2.0) 
              + (p[i, j+1] - 2.0*p[i, j] + p[i, j-1])/(dy**2.0) )        

        elif j == 0:

          lamda = - (2.0/(dx**2.0) + 1.0/(dy**2.0))

          R = (u_star[i+1, j] - u_star[i, j])/(dt*dx) + (v_star[i, j+1] - v_star[i, j])/(dt*dy)

          R = R - ( (p[i+1, j] - 2.0*p[i, j] + p[i-1, j])/(dx**2.0) 
              + (p[i, j+1] - 1.0*p[i, j])/(dy**2.0) )        

        elif j == N_y-1:

          lamda = - (2.0/(dx**2.0) + 1.0/(dy**2.0))

          R = (u_star[i+1, j] - u_star[i, j])/(dt*dx) + (v_star[i, j+1] - v_star[i, j])/(dt*dy)

          R = R - ( (p[i+1, j] - 2.0*p[i, j] + p[i-1, j])/(dx**2.0) 
              + (- 1.0*p[i, j] + p[i, j-1])/(dy**2.0) )        


        else:

          lamda = - (2.0/(dx**2.0) + 2.0/(dy**2.0))

          R = (u_star[i+1, j] - u_star[i, j])/(dt*dx) + (v_star[i, j+1] - v_star[i, j])/(dt*dy)

          R = R - ( (p[i+1, j] - 2.0*p[i, j] + p[i-1, j])/(dx**2.0) 
              + (p[i, j+1] - 2.0*p[i, j] + p[i, j-1])/(dy**2.0) )

        R = R/lamda
        p[i, j] = p[i, j] + 1.8*R

        if np.abs(R) > R_max:
          R_max = np.abs(R)

    error = R_max

    iter = iter + 1

  for i in range(0, N_x):
    p[i, -1] = p[i, 0]
    p[i, N_y] = p[i, N_y-1]

  for j in range(0, N_y):
    p[-1, j] = p[0, j]
    p[N_x, j] = p[N_x-1, j]

  p[-1, -1] = p[0, 0]
  p[-1, N_y] = p[0, N_y-1]
  p[N_x, -1] = p[N_x-1, 0]
  p[N_x, N_y] = p[N_x-1, N_y-1]

  for i in range(1, N_x):
    for j in range(-1, N_y+1):
      u[i, j] = ( u_star[i, j] 
                  - dt*(p[i, j] - p[i-1,j])/dx  )

  for i in range(-1, N_x+1):
    for j in range(1, N_y):
      v[i, j] = ( v_star[i, j]
                - dt*(p[i, j] - p[i, j-1])/dy  )


  print('Tempo:', t, 'Iter:', iter)

print(v)
u_plot = np.zeros((N_x+1, N_y+1), float)
v_plot = np.zeros((N_x+1, N_y+1), float)
p_plot = np.zeros((N_x+1, N_y+1), float)

for i in range(0, N_x+1):
  for j in range(0, N_y+1):
    u_plot[i, j] = 0.5*(u[i, j] + u[i, j-1])
    v_plot[i, j] = 0.5*(v[i, j] + v[i-1, j])
    p_plot[i, j] = 0.25*(p[i, j] + p[i-1, j] + p[i, j-1] + p[i-1, j-1])

u_plot_unit = np.zeros((N_x+1, N_y+1), float)
v_plot_unit = np.zeros((N_x+1, N_y+1), float)

for i in range(0, N_x+1):
  for j in range(0, N_y+1):
    if i != 0 and j != 0 and i != N_x:
      modulo = np.sqrt(u_plot[i, j]**2.0 + v_plot[i, j]**2.0)
      if modulo != 0.0:
        u_plot_unit[i, j] = u_plot[i, j]/modulo
        v_plot_unit[i, j] = v_plot[i, j]/modulo


x = np.linspace(0, L_x, N_x+1)
y = np.linspace(0, L_y, N_y+1)

#plt.quiver(x, y, np.transpose(u_plot_unit), np.transpose(v_plot_unit))

#plt.contourf(x, y, np.transpose(p_plot))

#plt.show()

#print(p_plot)

for j in range(0,N_y+1):
  print(u_plot[8,j])