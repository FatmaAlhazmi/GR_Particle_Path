import numpy as np
from scipy.integrate import solve_ivp
import sympy as sp

def christoffel_symbols(metric, coords):
    n = len(coords)
    christoffel = sp.MutableDenseNDimArray.zeros(n, n, n)
    inv_metric = metric.inv()
    
    for k in range(n):
        for i in range(n):
            for j in range(n):
                christoffel[i, j, k] = 0.5 * sum([inv_metric[k, l] * (sp.diff(metric[l, j], coords[i]) +
                    sp.diff(metric[l, i], coords[j]) - sp.diff(metric[i, j], coords[l])) for l in range(n)])
    return christoffel

def geodesic_equations(christoffel, coords):
    n = len(coords)
    
    def geodesic(t, Y):
        dydt = np.zeros(2 * n)
        dydt[:n] = Y[n:]
        
        for i in range(n):
            dydt[n + i] = -sum([
                float(christoffel[j, k, i].subs({coords[l]: Y[l] for l in range(n)}).evalf())
                * Y[n + j] * Y[n + k]
                for j in range(n) for k in range(n)
            ])
        return dydt

    return geodesic

def solve_geodesic(geodesic_func, initial_conditions, t_span, t_eval):
    sol = solve_ivp(geodesic_func, t_span, initial_conditions, t_eval=t_eval, method='RK45')
    return sol
