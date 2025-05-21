from pymoo.factory import get_problem
from pymoo.indicators.hv import HV
from pymoo.indicators.hv import Hypervolume
from pymoo.visualization.scatter import Scatter
import numpy as np

# Crear un problema de prueba multiobjetivo
problem = get_problem("zdt1")

# Supongamos que estas son las soluciones obtenidas (frente de Pareto)
# Estas soluciones deben estar en el espacio objetivo
solutions = np.array([
    [0.1, 0.9],
    [0.3, 0.6],
    [0.4, 0.5],
    [0.7, 0.3],
    [0.8, 0.2]
])

# Definir un punto de referencia adecuado para el cálculo del hipervolumen
reference_point = np.array([1.0, 1.0])


# Crear el indicador de hipervolumen
hv = Hypervolume(ref_point=reference_point)

# Calcular el hipervolumen
hypervolume_value = hv.do(solutions)

print(f"El hipervolumen calculado es: {hypervolume_value}")


# The pareto front of a scaled zdt1 problem
pf = get_problem("zdt1").pareto_front()

# The result found by an algorithm
A = pf[::10] * 1.1

# plot the result
Scatter(legend=True).add(
    pf, label="Pareto-front").add(A, label="Result").show()

ref_point = np.array([1.2, 1.2])

ind = HV(ref_point=ref_point)
print("HV", ind(A))

#En problemas de minimización, es mejor el algoritmo que tiene mayor hypervolumen