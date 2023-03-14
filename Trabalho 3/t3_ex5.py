import numpy as np
import matplotlib.pyplot as plt

def f(t,x):
    res = 0.0
    for n in range(100):
        res+= (((-1)**n)/((2*n+1)**2))*np.sin((2*n+1)*np.pi*x)*np.exp(-((2*n+1)**2)*(np.pi**2)*0.01*t)
    return (800/(np.pi**2))*res

def T_i(T,x,dx):
    T[0] = 0.0
    for i in range(1,N+1):
        x+=dx
        T[i] = 200*x
    return T

def solexata(N,T,tempo_final):
    x_aux = 0.0
    for i in range(1,N+1):
        x_aux+=dx
        T[i] = f(tempo_final,x_aux)
        
    plt.plot(T,'k')
    
def eqcondcal(N,dx,dt,tempo_final,T):
    tempo = 0.0

    c = dt/(dx*dx)
    k = 0

    while tempo <= tempo_final:
        T_new = np.copy(T)
        for i in range(1,N+1):
            #ghost-point
            if i == N:
                T_new[i] = T[i] +  0.01*c*(2.0*T[i-1]-2.0*T[i])
            
            else:
                T_new[i] = T[i] + 0.01*c*(T[i+1]-2.0*T[i]+T[i-1])
            
            
        tempo+=dt
        k+=1 
        T = np.copy(T_new)   
        if k%10 == 0:
            plt.plot(T,'o')
        
    

N = 50 #numero de passos
dx = 0.5/N
dt = (dx*dx)/2
tempo_final = 50.0
k = 0

if dt > (dx*dx)/2: #condicao de estabilidade de Von Neuman
    print('Erro')
 
else:  
    T = np.zeros(N+1)
    T_exata = np.zeros(N+1)
    x = 0.0
    T = T_i(T,x,dx)
        
    T_exata[0] = 0.0
    
    #Inicializando o Grafico

    graf = plt.figure()
    ax = graf.add_subplot()

    #Parâmetros do Grafico

    graf.suptitle('Comparação da Resolução da Equação do calor analítica e numérica', fontsize = 10, fontweight = 'bold', usetex = False)
    ax.set_ylabel(r'T', fontsize = 12, usetex = False)
    ax.set_xlabel(r't', fontsize = 12, usetex = False)
    
    #Calculando a Equacao da Conducao do Calor
    
    eqcondcal(N,dx,dt,tempo_final,T)
    solexata(N,T_exata,tempo_final)
    
    plt.show()