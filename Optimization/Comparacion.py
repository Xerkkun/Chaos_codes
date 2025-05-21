import sys
import numpy.random as random
impot scipy.stats as stats

n = len( sys.argv )
if n != 5 :
    print( "Args = media1 var1 media2 var2" )
    sys.exit(1)

m1 = floar ( sys.argv[1] )
v1 = floar ( sys.argv[2] )
m2 = floar ( sys.argv[3] )
v2 = floar ( sys.argv[4] )

# Hacemos dos distribuciones gaussianas:
vA = v1*random.randn( 30 ) + m1
vB = v2*random.randn( 30 ) + m2

# Las comparamos con la prueba de Wilcoxon
v, pvalue = stats.wilcoxon( vA, vB, alternative='two-sided' )
print( v, pvalue )

if pvalue > 0.05 :
    print( "Las distribuciones son iguales" )
else:
    print( "Las distribuciones NO son iguales" )
    v, pvalue = stats.wilcoxon( vA,vB, alternative='less' )
    if pvalue < 0.05 :
        print ("La distribución A es mayor que B")
    else:
        v, pvalue = stats.wilcoxon( vA, vB, alternative='greater')
        
