
from pymoo.algorithms.soo.nonconvex.de import DE
from pymoo.problems import get_problem
from pymoo.operators.sampling.lhs import LHS
from pymoo.optimize import minimize
import tarea1


problem = tarea1.Tarea1()


algorithm = DE(
    pop_size=16,
    sampling=LHS(),
    variant="DE/rand/1/bin",
    CR=0.6,
    dither="vector",
    jitter=False
)

res = minimize(problem,
               algorithm,
               seed=1,
               verbose=False)

#print("Best solution found: \nX = %s\nF = %s" % (res.X, res.F))
print(res.X, res.F)