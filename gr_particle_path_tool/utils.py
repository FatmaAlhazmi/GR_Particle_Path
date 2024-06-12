# utils.py
import sympy as sp
from scipy.integrate import solve_ivp

def christoffel_symbols(metric, coords, params):
    n = len(coords)
    christoffel = sp.MutableDenseNDimArray.zeros(n, n, n)
    inv_metric = metric.inv()
    for k in range(n):
        for i in range(n):
            for j in range(n):
                christoffel[i, j, k] = 0.5 * sum([inv_metric[k, l] * (metric[l, i].diff(coords[j]) + metric[l, j].diff(coords[i]) - metric[i, j].diff(coords[l]))
                                                  for l in range(n)])
    return christoffel

def geodesic_equations(christoffel, coords, params):
    n = len(coords)
    
    def geodesic(t, Y):
        dydt = [0] * (2 * n)
        dydt[:n] = Y[n:]
        
        for i in range(n):
            sum_term = 0
            for j in range(n):
                for k in range(n):
                    sum_term += float(christoffel[j, k, i].subs({coords[l]: Y[l] for l in range(n)}).evalf()) * Y[n + j] * Y[n + k]
            dydt[n + i] = -sum_term
        
        return dydt
    
    return geodesic

def solve_geodesic(geodesic_func, initial_conditions, t_span, t_eval):
    sol = solve_ivp(geodesic_func, t_span, initial_conditions, t_eval=t_eval, method='RK45')
    return sol
