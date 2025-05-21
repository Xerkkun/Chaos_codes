import pymoo.gradient.toolbox as anp
import numpy as np
import math
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.visualization.scatter import Scatter
from pymoo.util.plotting import plot

from pymoo.core.problem import ElementwiseProblem


class Tarea1(ElementwiseProblem):
    def __init__(self, n_var=1):
        super().__init__(n_var=n_var, n_obj=1, n_ieq_constr=0, xl=[0.0], xu=[7.0], vtype=float)

    def _evaluate(self, x, out, *args, **kwargs):
        #l = []
        # for i in range(x.shape[1] - 1):
            #val = (x[:, i + 1] - x[:, i] ** 2) ** 2 + (1 - x[:, i]) ** 2
        #print(x.shape)
        val = (x-2.0)*(x-5.0) + math.sin(1.5*math.pi*x)
        #l.append(val)
        out["F"] = val
