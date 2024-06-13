import unittest
from gr_particle_path_tool.event_horizon import detect_event_horizon

class TestEventHorizon(unittest.TestCase):

    def test_schwarzschild_event_horizon(self):
        result = detect_event_horizon('schwarzschild', 1.0)
        self.assertEqual(result, [2.0])  # Event horizon at r = 2M

    def test_kerr_event_horizon(self):
        result = detect_event_horizon('kerr', 1.0, 0.5)
        self.assertTrue(len(result) > 0)  # Should return event horizon roots

    def test_reissner_nordstrom_event_horizon(self):
        result = detect_event_horizon('reissner_nordstrom', 1.0, 0.5)
        self.assertTrue(len(result) > 0)  # Should return event horizon roots

if __name__ == '__main__':
    unittest.main()
