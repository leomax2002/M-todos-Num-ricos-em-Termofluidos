import numpy as np #importação de biblioteca

def mult_matrix(A,B):
    
    nl_a = len(A)
    if nl_a > 0:
        nc_a = len(A[0])
    nl_b = len(B)
    if nl_b > 0:
        nc_b = len(B[0])
    C = np.array(nl_a*[nc_b*[0]])
    print('C inicial:')
    print(C)
    aux_c = 0
    aux_l = 0
    if nc_a == nl_b:
        while aux_l < nl_a:
            aux = 0
            for j in range(nc_a):
                aux +=  A[aux_l][j]*B[j][aux_c]
                    
            C[aux_l][aux_c] = aux
            
            aux_c+=1
            if aux_c == nc_b:
                aux_l+=1
                i = 0
                aux_c = 0
    else:
        print('Não é possível realizar a multiplicação')    
    return C


A = np.array([[1,2,3],[4,5,6]])

B = np.array([[1,2,3],[4,5,6],[7,8,9]])

print(mult_matrix(A,B))