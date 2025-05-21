from pymoo.algorithms.soo.nonconvex.ga import GA
import tarea1
from pymoo.optimize import minimize

problem = tarea1.Tarea1()

algorithm = GA(
    pop_size=30,
    eliminate_duplicates=True)

res = minimize(problem,
               algorithm,
               seed=1,
               verbose=False)

#print("Best solution found: \nX = %s\nF = %s" % (res.X, res.F))
print(res.X, res.F)