import sys
import numpy as np

n = len( sys.argv)

if n != 3:
    print( 'Args: archivo1 archivo2' )
    sys.exit(1)
    
#Leer los datos de los dos frentes de Pareto

A = np.genfromtxt( sys.argv[1] )
B = np.genfromtxt( sys.argv[2] )

#concateno los dos frentes
A = np.concatenate (A,B)

n, m = A.shape

#Son soluciones que tenemos que comparar

def dominancia( s1, s2) :
    k = len(s1)
    prueba = np.zeros(k)
    
    i = 0
    
    while i<k :
        if s1[i] < s2[i]:
            prueba[i] = 1
        elif s1[i] > s2[i]:
            prueba[i] = -1
        i += 1
        
    # Recorremos el resultado de las pruebas
    i=0

    while i<k :
        if prueba[i] == 1:
            return 1
        
        elif prueba[i] == -1:
            return -1
        
        i += 1
        
    return 0

r = dominancia( [1,1], [9,9])
print(r)  

