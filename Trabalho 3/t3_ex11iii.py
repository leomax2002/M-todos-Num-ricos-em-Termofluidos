import numpy as np
import matplotlib.pyplot as plt

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
    r = []
    k = 0
    while True:
        for i in range(n):
            s = 0.0
            for j in range(n):
                s = s + A[i,j]*x[j]
            R[i] = (1.0/A[i,i])*(b[i] - s)
            x_novo[i] = x[i] + R[i]
        r.append(max(R))
        k+=1
        x = np.copy(x_novo)
        if np.sqrt(np.dot(R,R)) < atol:
            break
    return k, np.array(r)
#Gauss Seidel
def gauss(A,b,atol,n):
    x = np.zeros(n,float)
    R = np.zeros(n,float)  
    r = []
    k = 0 
    while True:
        for i in range(n):
            s = 0.0
            for j in range(n):
                s = s + A[i,j]*x[j]
            R[i] = (1.0/A[i,i])*(b[i]-s)
            x[i] = x[i] + R[i]
        r.append(max(R))
        k+=1
           
        if np.sqrt(np.dot(R,R)) < atol:
            break
    
    return k, np.array(r)

#SOR
def sor(A,b,atol,n):
    x = np.zeros(n,float)
    R = np.zeros(n,float)
    w = 1.5 #alterar w para maxima convergencia
    r = []
    k = 0
    while True:
        for i in range(n):
            s = 0.0
            for j in range(n):
                s = s + A[i,j]*x[j]
            R[i] = (1.0/A[i,i])*(b[i]-s)
            x[i] = x[i] + w*R[i]
        r.append(max(R))
        k+=1
        if np.sqrt(np.dot(R,R)) < atol:
            break
    
    return k, np.array(r)

#metodo gradiente conjugado
def gradconj(A,b,atol,n):
    k = 0
    r = []
    x = np.zeros(n,float)
    R = np.zeros(n,float)
    aux = np.dot(A,x)
    for i in range(n):
        R[i] = b[i] - aux[i]
    p = np.copy(R)
    while True:
        s = np.matmul(A,p)
        a = (np.dot(R,R))/(np.dot(p,s))
        R_aux = np.copy(R)
        for i in range(n):
            x[i] = x[i] + a*p[i]
            R[i] = R[i] - a*(s[i])
        r.append(max(R))
        k+=1
            
        if np.sqrt(np.dot(R,R)) < atol:
            break
        b = (np.dot(R,R))/(np.dot(R_aux,R_aux))
        for i in range(n):
            p[i] = R[i] + b*p[i]
    
    return k, np.array(r)

atol = 1e-7
n = 20

r_jacobi = np.zeros(n-1,int)
r_gauss = np.zeros(n-1,int)
r_sor = np.zeros(n-1,int)
r_gradconj = np.zeros(n-1,int)
i = 20    
A = mat_n(i)
b = np.zeros(i)
b[i-1] = -1
#k_jacobi, r_jacobi = jacobi(A,b,atol,i)
#k_gauss, r_gauss = gauss(A,b,atol,i)
#k_sor, r_sor = sor(A,b,atol,i)
k_gradconj, r_gradconj = gradconj(A,b,atol,i)

#k_jacobi = np.linspace(0,k_jacobi-1,k_jacobi)
#k_gauss = np.linspace(0,k_gauss-1,k_gauss)
#k_sor = np.linspace(0,k_sor-1,k_sor)
k_gradconj = np.linspace(0,k_gradconj-1,k_gradconj)

#Inicializando o Grafico

graf = plt.figure()
ax = graf.add_subplot()

#Parâmetros do Grafico

graf.suptitle('Iterações por R máximo', fontsize = 10, fontweight = 'bold', usetex = False)
ax.set_ylabel(r'y', fontsize = 12, usetex = False)
ax.set_xlabel(r'x', fontsize = 12, usetex = False)

#plt.plot(k_jacobi,r_jacobi)
#plt.plot(k_gauss,r_gauss, 'y')
#plt.plot(k_sor,r_sor, 'g')
plt.plot(k_gradconj,r_gradconj, 'r')
plt.show()