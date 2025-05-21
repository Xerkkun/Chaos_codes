import numpy as np
import sys
from scipy.stats import wilcoxon,levene

n = len( sys.argv )
if n != 5 :
    print( 'Args: media1 sigma1 media2 sigma2' )
    sys.exit(1)
    
m1 = float(sys.argv[1])
sigma1 = float(sys.argv[2])
m2 = float(sys.argv[3])
sigma2 = float(sys.argv[4])

v1= np.random.normal( loc=m1, scale=sigma1, size=30 )
v2 = np.random.normal( loc=m2, scale=sigma2, size=30 )

print( v1 )
print( v2 )

#Vamos a suponer que es un problema de minimización entonces el mejor es la menor media

s, p = wilcoxon( v1, v2, zero_method='wilcox', correction=True, alternative='two-sided')

if p >= 0.05 :
    print( 'Los resultados son iguales' )
else: 
    s, p = wilcoxon( v1, v2, zero_method='wilcox', correction=True,alternative='less' )
    if p >= 0.05 :
        print( 'El segundo algoritmo es mejor' )
    else:
        s, p = wilcoxon (v1, v2, zero_method='wilcox', correction=True, alternative='greater')
        if p >= 0.05 :
            print( 'El primer algoritmo es mejor')
        
        #hay que hacer las tres pruebas, si una no es mejor no quiere decir que la otra si lo sea
        #comparar también la desviación estandar
        
        
        #Programa para calcular el hipervolumen
        #Punto de nadir reportar para que sean reproducibles y el cálculo del hipervolumen depende de ese punto
        #Punto ideal es diferente al punto de nadir
        