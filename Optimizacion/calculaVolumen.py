from pymoo.indicators.hv import HV 
import numpy as np
import sys

n = len(sys.argv)

if n != 2:
    print( "Args: archivo_frente_Pareto" )
    sys.exit(1)
    
nombre = sys.argv[1]

A = np.loadtxt( nombre )
print( A.shape )

ref_point = np.array( [20.0, 80.0] )

# ind = HV(ref_point = ref_ point)
# print("HV", ind(A))