import numpy as np
from pymoo.visualization.pcp import PCP

X = np.loadtxt( "vars.txt" )
F = np.loadtxt( "objs.txt" )

print( X.shape )
print( F.shape )

PCP().add(F).show()