import unittest
import sympy as sp
from gr_particle_path_tool.utils import christoffel_symbols

class TestUtils(unittest.TestCase):

    def test_christoffel_symbols(self):
        M = sp.symbols('M')
        r, theta, phi = sp.symbols('r theta phi')
        metric = sp.Matrix([
            [1 - 2 * M / r, 0, 0, 0],
            [0, 1 / (1 - 2 * M / r), 0, 0],
            [0, 0, r**2, 0],
            [0, 0, 0, r**2 * sp.sin(theta)**2]
        ])
        coords = sp.symbols('t r theta phi')
        params = {M: 1}
        christoffel = christoffel_symbols(metric, coords, params)
        self.assertIsNotNone(christoffel)
        self.assertTrue((christoffel[1, 1, 1] - (-M / (r**2 * (1 - 2 * M / r)))).subs(params).simplify() == 0)

if __name__ == '__main__':
    unittest.main()
