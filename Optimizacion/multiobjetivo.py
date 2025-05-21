import pymoo.gradient.toolbox as anp
import numpy as np
import math
from pymoo.optimize import minimize

from pymoo.core.problem import ElementwiseProblem


class NuevoProblema(ElementwiseProblem):
    def __init__(self):
        super().__init__(n_var=2, n_obj=2, xl=[0.0], xu=[10.0] )

    def _evaluate(self, x, out, *args, **kwargs):
        out["F"] = np.column_stack( [ 0.5*(x**2+y**2), (x-4)*(x-4) + 4 ] )
