import numpy as np

def mediaedesvio(N):
    n = len(N)
    if n > 0:
        n2 = len(N[0])
        termos = n*n2
        res = 0
        for i in range(n):
            for j in range(n2):
                res+=N[i][j]
                
        media = res/termos
        
        res2 = 0
    
        for k in range(n):
            for l in range(n2):
                res2 += (media - N[i][j])**2
        dp = np.sqrt(res2/termos)
        
        return media, dp
    
    
    
A = np.array([[1,2,3],[4,5,6],[7,8,9]])
print(mediaedesvio(A))
                
            
                