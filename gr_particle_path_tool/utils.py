import sympy as sp

def christoffel_symbols(metric, coords):
    n = len(coords)
    christoffel = sp.MutableDenseNDimArray.zeros(n, n, n)  # Use MutableDenseNDimArray for 3D array
    inv_metric = metric.inv()

    for i in range(n):
        for j in range(n):
            for k in range(n):
                christoffel[i, j, k] = 0.5 * sum(inv_metric[i, m] * (sp.diff(metric[m, j], coords[k]) +
                                                                     sp.diff(metric[m, k], coords[j]) -
                                                                     sp.diff(metric[j, k], coords[m]))
                                                 for m in range(n))
    return christoffel

def geodesic_equations(christoffel, coords):
    n = len(coords)
    geodesic_eqs = []

    for mu in range(n):
        eq = sp.diff(sp.Symbol('x%d' % mu), 'tau', 2)
        for alpha in range(n):
            for beta in range(n):
                eq += -christoffel[mu, alpha, beta] * sp.diff(sp.Symbol('x%d' % alpha), 'tau') * sp.diff(sp.Symbol('x%d' % beta), 'tau')
        geodesic_eqs.append(eq)

    return geodesic_eqs

def solve_geodesic(equations, initial_conditions, t_span, t_eval):
    from scipy.integrate import solve_ivp
    import numpy as np

    def odes(t, y):
        dydt = []
        for eq in equations:
            subs = {sp.Symbol('tau'): t}
            subs.update({sp.Symbol(f'x{i}'): y[i] for i in range(len(y)//2)})
            subs.update({sp.Symbol(f'x{i}'): y[i+len(y)//2] for i in range(len(y)//2)})
            dydt.append(y[len(y)//2:][equations.index(eq)])
            dydt.append(eq.subs(subs))
        return np.array(dydt, dtype=float)

    sol = solve_ivp(odes, t_span, initial_conditions, t_eval=t_eval)
    return sol
