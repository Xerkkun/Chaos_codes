import pymoo.gradient.toolbox as anp
import numpy as np
import math
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter
from pymoo.util.plotting import plot

from pymoo.core.problem import ElementwiseProblem


class NuevoProblema(ElementwiseProblem):
    def __init__(self):
        super().__init__(n_var=1, n_obj=2, xl=[0.0], xu=[10.0] )

    def _evaluate(self, x, out, *args, **kwargs):
        out["F"] = np.column_stack( [ (x-2)*(x-2) + 4, (x-4)*(x-4) + 4 ] )
