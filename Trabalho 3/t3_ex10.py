import numpy as np
import matplotlib.pyplot as plt

#Jacobi
def jacobi(A,x,R,b,atol):
    x_novo = np.zeros(4)
    while True:
        for i in range(4):
            s = 0.0
            for j in range(4):
                s = s + A[i,j]*x[j]
            R[i] = (1.0/A[i,i])*(b[i] - s)
            x_novo[i] = x[i] + R[i]
            
        k+=1
        x = np.copy(x_novo)
        if np.sqrt(np.dot(R,R)) < atol:
            break

    return x
#Gauss Seidel
def gauss(A,x,R,b,atol):
    while True:
        for i in range(4):
            s = 0.0
            for j in range(4):
                s = s + A[i,j]*x[j]
            R[i] = (1.0/A[i,i])*(b[i]-s)
            x[i] = x[i] + R[i]
            
        if np.sqrt(np.dot(R,R)) < atol:
            break
    
    return x

#SOR
def sor(A,x,R,b,atol):
    w = 1 #alterar w para maxima convergencia
    while True:
        for i in range(4):
            s = 0.0
            for j in range(4):
                s = s + A[i,j]*x[j]
            R[i] = (1.0/A[i,i])*(b[i]-s)
            x[i] = x[i] + w*R[i]
            
        if np.sqrt(np.dot(R,R)) < atol:
            break
    
    return x
    
A = np.array([[12.0,-2.0,3.0,1.0],[-2.0,15.0,6.0,-3.0],[1.0,6.0,20.0,-4.0],[0.0,-3.0,2.0,9.0]])
b = np.array([0.0,0.0,20.0,0.0])
x = np.zeros(4,float)
R = np.zeros(4,float)
atol = 1e-8
print(jacobi(A,x,R,b,atol))
print(gauss(A,x,R,b,atol))
print(sor(A,x,R,b,atol))