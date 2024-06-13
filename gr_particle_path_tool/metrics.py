import sympy as sp

def schwarzschild_metric(M):
    t, r, theta, phi = sp.symbols('t r theta phi')
    g = sp.zeros(4)
    g[0, 0] = -(1 - 2 * M / r)
    g[1, 1] = 1 / (1 - 2 * M / r)
    g[2, 2] = r**2
    g[3, 3] = r**2 * sp.sin(theta)**2
    return g

def christoffel_symbols(g, coords):
    n = len(coords)
    Gamma = sp.MutableDenseNDimArray(sp.zeros((n, n, n)))
    g_inv = g.inv()
    for k in range(n):
        for i in range(n):
            for j in range(n):
                Gamma[k, i, j] = 0.5 * sum([g_inv[k, l] * (sp.diff(g[l, j], coords[i]) +
                                                            sp.diff(g[l, i], coords[j]) -
                                                            sp.diff(g[i, j], coords[l]))
                                            for l in range(n)])
    return Gamma

def geodesic_equations(Gamma, coords):
    n = len(coords)
    def geodesic(t, Y):
        dYdt = sp.MutableDenseNDimArray(sp.zeros(2 * n))
        for mu in range(n):
            dYdt[mu] = Y[n + mu]
            dYdt[n + mu] = -sum([Gamma[mu, i, j].subs([(coords[k], Y[k]) for k in range(n)]) * Y[n + i] * Y[n + j]
                                 for i in range(n) for j in range(n)])
        return dYdt
    return geodesic

def solve_geodesic(geodesic_func, initial_conditions, t_span, t_eval):
    n = len(initial_conditions) // 2
    def rhs(t, Y):
        return [float(expr) for expr in geodesic_func(t, Y)]
    sol = solve_ivp(rhs, t_span, initial_conditions, t_eval=t_eval)
    return sol


def kerr_metric(M, a):
    t, r, theta, phi = sp.symbols('t r theta phi')
    g = sp.zeros(4, 4)
    rho2 = r**2 + (a * sp.cos(theta))**2
    delta = r**2 - 2 * M * r + a**2
    g[0, 0] = -(1 - 2 * M * r / rho2)
    g[0, 3] = g[3, 0] = -2 * M * a * r * sp.sin(theta)**2 / rho2
    g[1, 1] = rho2 / delta
    g[2, 2] = rho2
    g[3, 3] = (r**2 + a**2 + 2 * M * a**2 * r * sp.sin(theta)**2 / rho2) * sp.sin(theta)**2
    return g

def reissner_nordstrom_metric(M, Q):
    t, r, theta, phi = sp.symbols('t r theta phi')
    g = sp.zeros(4, 4)
    g[0, 0] = -(1 - 2 * M / r + Q**2 / r**2)
    g[1, 1] = 1 / (1 - 2 * M / r + Q**2 / r**2)
    g[2, 2] = r**2
    g[3, 3] = r**2 * sp.sin(theta)**2
    return g

def kerr_newman_metric(M, a, Q):
    t, r, theta, phi = sp.symbols('t r theta phi')
    g = sp.zeros(4, 4)
    rho2 = r**2 + (a * sp.cos(theta))**2
    delta = r**2 - 2 * M * r + a**2 + Q**2
    g[0, 0] = -(1 - 2 * M * r / rho2 + Q**2 / rho2)
    g[0, 3] = g[3, 0] = -2 * M * a * r * sp.sin(theta)**2 / rho2
    g[1, 1] = rho2 / delta
    g[2, 2] = rho2
    g[3, 3] = (r**2 + a**2 + 2 * M * a**2 * r * sp.sin(theta)**2 / rho2) * sp.sin(theta)**2
    return g

def flrw_metric(a, k):
    t, r, theta, phi = sp.symbols('t r theta phi')
    g = sp.zeros(4, 4)
    g[0, 0] = -1
    g[1, 1] = a**2 / (1 - k * r**2)
    g[2, 2] = a**2 * r**2
    g[3, 3] = a**2 * r**2 * sp.sin(theta)**2
    return g

def de_sitter_metric(Lambda):
    t, r, theta, phi = sp.symbols('t r theta phi')
    g = sp.zeros(4, 4)
    g[0, 0] = -(1 - Lambda * r**2 / 3)
    g[1, 1] = 1 / (1 - Lambda * r**2 / 3)
    g[2, 2] = r**2
    g[3, 3] = r**2 * sp.sin(theta)**2
    return g

def anti_de_sitter_metric(Lambda):
    t, r, theta, phi = sp.symbols('t r theta phi')
    g = sp.zeros(4, 4)
    g[0, 0] = -(1 + Lambda * r**2 / 3)
    g[1, 1] = 1 / (1 + Lambda * r**2 / 3)
    g[2, 2] = r**2
    g[3, 3] = r**2 * sp.sin(theta)**2
    return g

def minkowski_metric():
    t, r, theta, phi = sp.symbols('t r theta phi')
    g = sp.eye(4)
    g[0, 0] = -1
    return g

def vaidya_metric(M):
    t, r, theta, phi = sp.symbols('t r theta phi')
    g = sp.zeros(4, 4)
    g[0, 0] = -(1 - 2 * M / r)
    g[0, 1] = g[1, 0] = 2 * M / r
    g[1, 1] = 1 / (1 - 2 * M / r)
    g[2, 2] = r**2
    g[3, 3] = r**2 * sp.sin(theta)**2
    return g

def gottingen_metric():
    t, r, theta, phi = sp.symbols('t r theta phi')
    g = sp.eye(4)
    return g

def bertotti_robinson_metric(r0):
    t, r, theta, phi = sp.symbols('t r theta phi')
    g = sp.zeros(4, 4)
    g[0, 0] = -(r**2 / r0**2)
    g[1, 1] = r**2 / r0**2
    g[2, 2] = r0**2
    g[3, 3] = r0**2 * sp.sin(theta)**2
    return g
