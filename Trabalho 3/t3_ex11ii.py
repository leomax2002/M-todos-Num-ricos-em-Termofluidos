import numpy as np
import matplotlib.pyplot as plt
import time

def mat_n(n):
    m = np.repeat([np.zeros(n,int)],n, axis = 0)
    for i in range(n):
        for j in range(n):
            if j == i-1 or j == i+1:
                m[i][j] = 1
            elif j == i:
                m[i][j] = -2
            else:
                m[i][j] = 0
    return m

#Jacobi
def jacobi(A,b,atol,n):
    x = np.zeros(n,float)
    R = np.zeros(n,float)
    x_novo = np.zeros(n)
    t1 = time.time()
    while True:
        for i in range(n):
            s = 0.0
            for j in range(n):
                s = s + A[i,j]*x[j]
            R[i] = (1.0/A[i,i])*(b[i] - s)
            x_novo[i] = x[i] + R[i]
        x = np.copy(x_novo)
        if np.sqrt(np.dot(R,R)) < atol:
            break
    t2 = time.time()
    t = t2 - t1
    return x,t
#Gauss Seidel
def gauss(A,b,atol,n):
    x = np.zeros(n,float)
    R = np.zeros(n,float)  
    t1 = time.time()  
    while True:
        for i in range(n):
            s = 0.0
            for j in range(n):
                s = s + A[i,j]*x[j]
            R[i] = (1.0/A[i,i])*(b[i]-s)
            x[i] = x[i] + R[i]
           
        if np.sqrt(np.dot(R,R)) < atol:
            break
    t2 = time.time()
    t = t2-t1
    
    return x,t

#SOR
def sor(A,b,atol,n):
    x = np.zeros(n,float)
    R = np.zeros(n,float)
    w = 1.5 #alterar w para maxima convergencia
    t1 = time.time()
    while True:
        for i in range(n):
            s = 0.0
            for j in range(n):
                s = s + A[i,j]*x[j]
            R[i] = (1.0/A[i,i])*(b[i]-s)
            x[i] = x[i] + w*R[i]
            
        if np.sqrt(np.dot(R,R)) < atol:
            break
    t2 = time.time()
    t = t2-t1
    
    return x,t

#metodo gradiente conjugado
def gradconj(A,b,atol,n):
    k = 0
    x = np.zeros(n,float)
    R = np.zeros(n,float)
    aux = np.dot(A,x)
    for i in range(n):
        R[i] = b[i] - aux[i]
    p = np.copy(R)
    t1 = time.time()
    while True:
        s = np.matmul(A,p)
        a = (np.dot(R,R))/(np.dot(p,s))
        R_aux = np.copy(R)
        for i in range(n):
            x[i] = x[i] + a*p[i]
            R[i] = R[i] - a*(s[i])
            
        if np.sqrt(np.dot(R,R)) < atol:
            break
        b = (np.dot(R,R))/(np.dot(R_aux,R_aux))
        for i in range(n):
            p[i] = R[i] + b*p[i]
    t2 = time.time()
    t = t2-t1
    
    return x,t

atol = 1e-7
n = 20

t_jacobi = np.zeros(n-1,float)
t_gauss = np.zeros(n-1,float)
t_sor = np.zeros(n-1,float)
t_gradconj = np.zeros(n-1,float)
    
for i in range(2,n+1):
    A = mat_n(i)
    b = np.zeros(i)
    b[i-1] = -1
    x_aux, t_jacobi[i-2] = jacobi(A,b,atol,i)
    x_aux, t_gauss[i-2] = gauss(A,b,atol,i)
    x_aux, t_sor[i-2] = sor(A,b,atol,i)
    x_aux, t_gradconj[i-2] = gradconj(A,b,atol,i)
    
n = np.linspace(2,n,n-1)

#Inicializando o Grafico

graf = plt.figure()
ax = graf.add_subplot()

#Parâmetros do Grafico

graf.suptitle('Tamanho da matriz por tempo de execução das iterações', fontsize = 10, fontweight = 'bold', usetex = False)
ax.set_ylabel(r't', fontsize = 12, usetex = False)
ax.set_xlabel(r'n', fontsize = 12, usetex = False)

#plt.plot(n,t_jacobi)
#plt.plot(n,t_gauss)
#plt.plot(n,t_sor)
plt.plot(n,t_gradconj, 'r')
plt.show()