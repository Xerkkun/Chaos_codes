from pymoo.algorithms.moo.moead import MOEAD
from pymoo.optimize import minimize
from pymoo.util.ref_dirs import get_reference_directions 
import numpy as np
import NuevoProblema

problem = NuevoProblema.NuevoProblema()
#El número de vectores, "particiones" en pymoo, es igual
#al tamaño de soluciones

vectores = get_reference_directions("uniform", 2, n_partitions=30)

algorithm = MOEAD(ref_dirs=vectores, n_neighbors=3, prob_neighbor_mating=0.7)

res = minimize(problem,
                algorithm,
                ('n_gen', 20),
                seed=1,
                verbose=True)

np.savetxt( "vars.txt", res.X )
np.savetxt( "objs.txt", res.F )
