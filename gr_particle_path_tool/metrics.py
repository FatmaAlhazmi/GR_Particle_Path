import sympy as sp

def schwarzschild_metric():
    t, r, theta, phi = sp.symbols('t r theta phi')
    M = sp.symbols('M')
    g_tt = -(1 - 2 * M / r)
    g_rr = 1 / (1 - 2 * M / r)
    g_thetatheta = r**2
    g_phiphi = r**2 * sp.sin(theta)**2
    return sp.Matrix([
        [g_tt, 0, 0, 0],
        [0, g_rr, 0, 0],
        [0, 0, g_thetatheta, 0],
        [0, 0, 0, g_phiphi]
    ])

def kerr_metric():
    t, r, theta, phi = sp.symbols('t r theta phi')
    M, a = sp.symbols('M a')
    rho2 = r**2 + (a * sp.cos(theta))**2
    delta = r**2 - 2 * M * r + a**2
    sigma2 = (r**2 + a**2)**2 - a**2 * delta * sp.sin(theta)**2
    g_tt = -(1 - 2 * M * r / rho2)
    g_rr = rho2 / delta
    g_thetatheta = rho2
    g_phiphi = sigma2 / rho2 * sp.sin(theta)**2
    g_tr = 0
    g_ttheta = 0
    g_tphi = -2 * M * a * r * sp.sin(theta)**2 / rho2
    g_rtheta = 0
    g_rphi = 0
    g_thetaphi = -2 * M * a * r * sp.sin(theta)**2 / rho2

    return sp.Matrix([
        [g_tt, g_tr, g_ttheta, g_tphi],
        [g_tr, g_rr, g_rtheta, g_rphi],
        [g_ttheta, g_rtheta, g_thetatheta, g_thetaphi],
        [g_tphi, g_rphi, g_thetaphi, g_phiphi]
    ])

def reissner_nordstrom_metric():
    t, r, theta, phi = sp.symbols('t r theta phi')
    M, Q = sp.symbols('M Q')
    g = sp.diag(-(1 - 2*M/r + Q**2/r**2), 1/(1 - 2*M/r + Q**2/r**2), r**2, r**2 * sp.sin(theta)**2)
    return g

def kerr_newman_metric():
    t, r, theta, phi = sp.symbols('t r theta phi')
    M, a, Q = sp.symbols('M a Q')
    rho2 = r**2 + a**2 * sp.cos(theta)**2
    delta = r**2 - 2*M*r + a**2 + Q**2
    g = sp.zeros(4)
    g[0, 0] = -(1 - 2*M*r / rho2 + Q**2 / rho2)
    g[0, 3] = g[3, 0] = -2*M*r*a*sp.sin(theta)**2 / rho2
    g[1, 1] = rho2 / delta
    g[2, 2] = rho2
    g[3, 3] = (r**2 + a**2 + 2*M*r*a**2*sp.sin(theta)**2 / rho2 + Q**2 / rho2) * sp.sin(theta)**2
    return g

def flrw_metric():
    t, r, theta, phi = sp.symbols('t r theta phi')
    a = sp.Function('a')(t)
    g = sp.diag(-1, a**2 / (1 - k*r**2), a**2 * r**2, a**2 * r**2 * sp.sin(theta)**2)
    return g

def de_sitter_metric():
    t, r, theta, phi = sp.symbols('t r theta phi')
    Lambda = sp.Symbol('Lambda')
    g = sp.diag(-1 + Lambda*r**2/3, 1/(1 - Lambda*r**2/3), r**2, r**2 * sp.sin(theta)**2)
    return g

def anti_de_sitter_metric():
    t, r, theta, phi = sp.symbols('t r theta phi')
    L = sp.Symbol('L')
    g = sp.diag(-1 + r**2/L**2, 1/(1 - r**2/L**2), r**2, r**2 * sp.sin(theta)**2)
    return g

def minkowski_metric():
    t, x, y, z = sp.symbols('t x y z')
    g = sp.diag(-1, 1, 1, 1)
    return g

def vaidya_metric():
    v, r, theta, phi = sp.symbols('v r theta phi')
    M = sp.Function('M')(v)
    g = sp.Matrix([
        [-1 + 2*M/r, 2, 0, 0],
        [2, 0, 0, 0],
        [0, 0, r**2, 0],
        [0, 0, 0, r**2 * sp.sin(theta)**2]
    ])
    return g

def gottingen_metric():
    # Placeholder: replace with actual Göttingen metric if available
    t, r, theta, phi = sp.symbols('t r theta phi')
    g = sp.diag(-1, 1, r**2, r**2 * sp.sin(theta)**2)
    return g

def bertotti_robinson_metric():
    # Placeholder: replace with actual Bertotti-Robinson metric if available
    t, r, theta, phi = sp.symbols('t r theta phi')
    g = sp.diag(-1 + r**2, 1/(1 - r**2), r**2, r**2 * sp.sin(theta)**2)
    return g
