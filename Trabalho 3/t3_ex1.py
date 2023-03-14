import numpy as np
import matplotlib.pyplot as plt

def f(t,x):
    return np.exp(((-np.pi**2)*t)/4)*np.sin((np.pi/2)*x)

def T_i(T,x,dx):
    T[0] = 0.0
    T[N] = 0.0
    for i in range(1,N):
        x+=dx
        T[i] =  np.sin((np.pi/2.0)*x)
    return T

def solexata(N,T,tempo_final):
    x_aux = 0.0
    for i in range(1,N):
        x_aux+=dx
        T[i] = f(tempo_final,x_aux)
        
    plt.plot(T,'k')

def eqcondcal(N,dx,dt,tempo_final,T):
    tempo = 0.0

    c = dt/(dx*dx)
    k = 0

    while tempo <= tempo_final:
        T_new = np.copy(T)
        for i in range(1,N):

            T_new[i] = T[i] + c*(T[i+1]-2.0*T[i]+T[i-1])
            
        tempo+=dt
        k+=1    
        T = np.copy(T_new)
        if k%1000 == 0:
            plt.plot(T,'o')
        
    

N = 100 #numero de passos
dx = 2.0/N
dt = (dx*dx)/4
tempo_final = 10.0

k = 0

if dt > (dx*dx)/2: #condicao de estabilidade de Von Neuman
    print('Erro')
 
else:  
    T = np.zeros(N+1)
    T_exata = np.zeros(N+1)
    x = 0.0
    T = T_i(T,x,dx)
        
    T_exata[0] = 0.0
    T_exata[N] = 0.0
    
    #Inicializando o Grafico

    graf = plt.figure()
    ax = graf.add_subplot()
    
    #Parâmetros do Gráfico

    graf.suptitle('Comparação da Resolução da Equação do calor analítica e numérica', fontsize = 10, fontweight = 'bold', usetex = False)
    ax.set_ylabel(r'T', fontsize = 12, usetex = False)
    ax.set_xlabel(r't', fontsize = 12, usetex = False)
    
    #Calculando a Equacao da Conducao do Calor
    
    eqcondcal(N,dx,dt,tempo_final,T)
    solexata(N,T_exata,tempo_final)
 
    plt.show()