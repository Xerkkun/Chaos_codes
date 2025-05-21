from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.problems import get_problem
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter
from pymoo.util.plotting import plot
import numpy as np

problem = get_problem("osy")

algorithm = NSGA2(pop_size=200)

res = minimize(problem,
                algorithm,
                ('n_gen', 250),
                seed=2,
                verbose=False)

print( res.X.shape )
print( res.F.shape )
print( res.CV.shape )
print( res.G.shape )
print( res.H.shape )
print( res.pop.shape )
#print( res.pf.shape )

print("Soluciones factibles = ", res.feas )
#"osy.py" 45L, 1284C

np.savetxt( "vars.txt", res.X )
np.savetxt( "objs.txt", res.F )
np.savetxt( "feas.txt", res.CV )
np.savetxt( "const.txt", res.G )

plot = Scatter()
plot.add(problem.pareto_front(), plot_type="line",color="black", alpha=0.7)
plot.add(res.F, facecolor="none",edgecolor="red")
plot.show()

#problem = get_problem("osy")
#plot(problem.pareto_front(), no_fill=True)
