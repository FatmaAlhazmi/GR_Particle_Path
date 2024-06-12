import unittest
import sympy as sp
import numpy as np
from gr_particle_path_tool.metrics import schwarzschild_metric, kerr_metric
from gr_particle_path_tool.utils import christoffel_symbols, geodesic_equations, solve_geodesic

class TestMetrics(unittest.TestCase):

    def test_schwarzschild_metric(self):
        M = sp.symbols('M')
        r, theta, phi = sp.symbols('r theta phi')
        expected_metric = sp.Matrix([
            [1 - 2 * M / r, 0, 0, 0],
            [0, 1 / (1 - 2 * M / r), 0, 0],
            [0, 0, r**2, 0],
            [0, 0, 0, r**2 * sp.sin(theta)**2]
        ]).subs(M, 1)
        metric = schwarzschild_metric().subs(M, 1)
        if not metric.equals(expected_metric):
            print("Schwarzschild Metric Test Failed")
            print("Expected Metric:")
            sp.pprint(expected_metric)
            print("Calculated Metric:")
            sp.pprint(metric)
        self.assertTrue(metric.equals(expected_metric))

    def test_kerr_metric(self):
        M, a = sp.symbols('M a')
        r, theta, phi = sp.symbols('r theta phi')
        sin = sp.sin  # Ensure sin function is correctly defined
        cos = sp.cos  # Ensure cos function is correctly defined
        # Define the expected Kerr metric here
        rho2 = r**2 + (a * cos(theta))**2
        delta = r**2 - 2 * M * r + a**2
        expected_metric = sp.Matrix([
            [-1 + 2*M*r/rho2, 0, 0, -2*M*r*a*sin(theta)**2/rho2],
            [0, rho2/delta, 0, 0],
            [0, 0, rho2, 0],
            [-2*M*r*a*sin(theta)**2/rho2, 0, 0, (r**2 + a**2 + 2*M*r*a**2*sin(theta)**2/rho2) * sin(theta)**2]
        ]).subs({M: 1, a: 0.5})
        metric = kerr_metric().subs({M: 1, a: 0.5})
        if not metric.equals(expected_metric):
            print("Kerr Metric Test Failed")
            print("Expected Metric:")
            sp.pprint(expected_metric)
            print("Calculated Metric:")
            sp.pprint(metric)
        self.assertTrue(metric.equals(expected_metric))

    def test_geodesic_solution(self):
        M = sp.symbols('M')
        metric = schwarzschild_metric().subs(M, 1)
        coords = sp.symbols('t r theta phi')
        parameters = {M: 1}
        christoffel = christoffel_symbols(metric, coords, parameters)
        initial_conditions = [0, 10, sp.pi/2, 0, 0, -0.1, 0, 0.05]
        t_span = (0, 100)
        t_eval = np.linspace(0, 100, 1000)
        geodesic_func = geodesic_equations(christoffel, coords, parameters)
        solution = solve_geodesic(geodesic_func, initial_conditions, t_span, t_eval)
        self.assertIsNotNone(solution)

if __name__ == '__main__':
    unittest.main()
