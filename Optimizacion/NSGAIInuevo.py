from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter
import numpy as np
import NuevoProblema

problem = NuevoProblema.NuevoProblema()

algorithm = NSGA2(pop_size=50)

res = minimize(problem,
                algorithm,
                ('n_gen', 20),
                seed=1,
                verbose=True)

np.savetxt( "vars.txt", res.X )
np.savetxt( "objs.txt", res.F )
np.savetxt( "feas.txt", res.CV )
np.savetxt( "const.txt", res.G )
