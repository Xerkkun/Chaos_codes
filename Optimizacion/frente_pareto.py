import matplotlib.pyplot as plt
from pymoo.factory import get_problem
from pymoo.optimize import minimize
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.termination.default import DefaultMultiObjectiveTermination
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
import numpy as np

# Definir el problema de optimización multiobjetivo
problem = get_problem("zdt1")

# Configurar el algoritmo NSGA-II
algorithm = NSGA2(pop_size=100)

# Definir el criterio de terminación
termination = DefaultMultiObjectiveTermination()

# Ejecutar la optimización
res = minimize(problem,
               algorithm,
               termination,
               seed=1,
               save_history=True,
               verbose=True)

# Extraer las soluciones obtenidas
solutions = res.F

# Calcular el frente de Pareto
pareto_front = NonDominatedSorting().do(
    solutions, only_non_dominated_front=True)

# Filtrar las soluciones no dominadas
pareto_solutions = solutions[pareto_front]

print("Frente de Pareto:")
print(pareto_solutions)


# Visualizar el frente de Pareto
plt.scatter(solutions[:, 0], solutions[:, 1], label="Soluciones")
plt.scatter(pareto_solutions[:, 0], pareto_solutions[:,
            1], color='r', label="Frente de Pareto")
plt.xlabel("Objetivo 1")
plt.ylabel("Objetivo 2")
plt.legend()
plt.title("Frente de Pareto")
plt.show()
