# tests/test_metrics.py

import unittest
import sympy as sp
from gr_particle_path_tool.metrics import (
    schwarzschild_metric, kerr_metric, reissner_nordstrom_metric, kerr_newman_metric,
    flrw_metric, de_sitter_metric, anti_de_sitter_metric, minkowski_metric,
    vaidya_metric, gottingen_metric, bertotti_robinson_metric
)

class TestMetrics(unittest.TestCase):

    def test_schwarzschild_metric(self):
        M = 1
        metric = schwarzschild_metric(M)
        r, theta = sp.symbols('r theta')
        expected_metric = sp.Matrix([
            [-(1 - 2 * M / r), 0, 0, 0],
            [0, 1 / (1 - 2 * M / r), 0, 0],
            [0, 0, r**2, 0],
            [0, 0, 0, r**2 * sp.sin(theta)**2]
        ])
        self.assertTrue(metric.equals(expected_metric))

    def test_kerr_metric(self):
        M, a = 1, 0.5
        metric = kerr_metric(M, a)
        r, theta = sp.symbols('r theta')
        Sigma = r**2 + a**2 * sp.cos(theta)**2
        Delta = r**2 - 2 * M * r + a**2
        expected_metric = sp.Matrix([
            [-(1 - 2 * M * r / Sigma), 0, 0, -2 * M * r * a * sp.sin(theta)**2 / Sigma],
            [0, Sigma / Delta, 0, 0],
            [0, 0, Sigma, 0],
            [-2 * M * r * a * sp.sin(theta)**2 / Sigma, 0, 0, (r**2 + a**2 + 2 * M * r * a**2 * sp.sin(theta)**2 / Sigma) * sp.sin(theta)**2]
        ])
        self.assertTrue(metric.equals(expected_metric))

    def test_reissner_nordstrom_metric(self):
        M, Q = 1, 0.5
        metric = reissner_nordstrom_metric(M, Q)
        r, theta = sp.symbols('r theta')
        expected_metric = sp.Matrix([
            [-(1 - 2 * M / r + Q**2 / r**2), 0, 0, 0],
            [0, 1 / (1 - 2 * M / r + Q**2 / r**2), 0, 0],
            [0, 0, r**2, 0],
            [0, 0, 0, r**2 * sp.sin(theta)**2]
        ])
        self.assertTrue(metric.equals(expected_metric))

    def test_kerr_newman_metric(self):
        M, a, Q = 1, 0.5, 0.3
        metric = kerr_newman_metric(M, a, Q)
        r, theta = sp.symbols('r theta')
        Sigma = r**2 + a**2 * sp.cos(theta)**2
        Delta = r**2 - 2 * M * r + a**2 + Q**2
        expected_metric = sp.Matrix([
            [-(1 - 2 * M * r / Sigma + Q**2 / Sigma), 0, 0, -2 * M * r * a * sp.sin(theta)**2 / Sigma],
            [0, Sigma / Delta, 0, 0],
            [0, 0, Sigma, 0],
            [-2 * M * r * a * sp.sin(theta)**2 / Sigma, 0, 0, (r**2 + a**2 + 2 * M * r * a**2 * sp.sin(theta)**2 / Sigma) * sp.sin(theta)**2]
        ])
        self.assertTrue(metric.equals(expected_metric))

    def test_flrw_metric(self):
        a, k = 1, 0
        metric = flrw_metric(a, k)
        r, theta = sp.symbols('r theta')
        expected_metric = sp.Matrix([
            [-1, 0, 0, 0],
            [0, a**2 / (1 - k * r**2), 0, 0],
            [0, 0, a**2 * r**2, 0],
            [0, 0, 0, a**2 * r**2 * sp.sin(theta)**2]
        ])
        self.assertTrue(metric.equals(expected_metric))

    def test_de_sitter_metric(self):
        Lambda = 1
        metric = de_sitter_metric(Lambda)
        r, theta = sp.symbols('r theta')
        expected_metric = sp.Matrix([
            [-(1 - Lambda * r**2 / 3), 0, 0, 0],
            [0, 1 / (1 - Lambda * r**2 / 3), 0, 0],
            [0, 0, r**2, 0],
            [0, 0, 0, r**2 * sp.sin(theta)**2]
        ])
        self.assertTrue(metric.equals(expected_metric))

    def test_anti_de_sitter_metric(self):
        Lambda = 1
        metric = anti_de_sitter_metric(Lambda)
        r, theta = sp.symbols('r theta')
        expected_metric = sp.Matrix([
            [-(1 + Lambda * r**2 / 3), 0, 0, 0],
            [0, 1 / (1 + Lambda * r**2 / 3), 0, 0],
            [0, 0, r**2, 0],
            [0, 0, 0, r**2 * sp.sin(theta)**2]
        ])
        self.assertTrue(metric.equals(expected_metric))

    def test_minkowski_metric(self):
        metric = minkowski_metric()
        expected_metric = sp.Matrix([
            [-1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ])
        self.assertTrue(metric.equals(expected_metric))

    def test_vaidya_metric(self):
        M = 1
        metric = vaidya_metric(M)
        r, theta = sp.symbols('r theta')
        expected_metric = sp.Matrix([
            [-(1 - 2 * M / r), 2 * M / r, 0, 0],
            [2 * M / r, 1 / (1 - 2 * M / r), 0, 0],
            [0, 0, r**2, 0],
            [0, 0, 0, r**2 * sp.sin(theta)**2]
        ])
        self.assertTrue(metric.equals(expected_metric))

    def test_gottingen_metric(self):
        metric = gottingen_metric()
        r, theta = sp.symbols('r theta')
        expected_metric = sp.Matrix([
            [-1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, r**2, 0],
            [0, 0, 0, r**2 * sp.sin(theta)**2]
        ])
        self.assertTrue(metric.equals(expected_metric))

    def test_bertotti_robinson_metric(self):
        r0 = 1
        metric = bertotti_robinson_metric(r0)
        r, theta = sp.symbols('r theta')
        expected_metric = sp.Matrix([
            [-(r**2 / r0**2), 0, 0, 0],
            [0, r**2 / r0**2, 0, 0],
            [0, 0, r0**2, 0],
            [0, 0, 0, r0**2 * sp.sin(theta)**2]
        ])
        self.assertTrue(metric.equals(expected_metric))

if __name__ == '__main__':
    unittest.main()
