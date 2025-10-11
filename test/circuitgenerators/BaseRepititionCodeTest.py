import unittest

import sys
import pathlib
sys.path.append("..")
import ddq_circuit_generators as scg
import numpy as np

class TestBaseCircuitGeneratorHelpers(unittest.TestCase):
    def setUp(self):
        self.d = 3
        _, _, _, self.c2i = scg.prepare_coords(3)
        return super().setUp()
    
    def test_data_coords(self):
        coords = scg.data_coords(self.d)
        self.assertListEqual(coords, [(1, 1), (2, 1), (3, 1), (1, 2), (2, 2), (3, 2), (1, 3), (2, 3), (3, 3)])

    def test_z_measure_coords(self):
        d=3
        coords = scg.z_measure_coords(d)
        self.assertListEqual(coords, [(0.5, 1.5), (2.5, 1.5), (1.5, 2.5), (3.5, 2.5)])

    def test_x_measure_coords(self):
        d=3
        coords = scg.x_measure_coords(d)
        self.assertListEqual(coords, [(2.5, 0.5), (1.5, 1.5), (2.5, 2.5), (1.5, 3.5)])

class TestConstantNoiseRepCodeGenerator(unittest.TestCase):
    def setUp(self):
        np.random.seed(seed=233423)
        self.uniform_string = open('test/circuitgenerators/artifacts/uniform_noise_d_3_string.txt')
        self.uniform_defect_string = open('test/circuitgenerators/artifacts/uniform_noise_defect_center_d_3_string.txt')
        self.heterogeneous_string = open('test/circuitgenerators/artifacts/heterogeneous_noise_d_3_string.txt')
        self.heterogeneous_defect_string = open('test/circuitgenerators/artifacts/heterogeneous_noise_defect_center_d_3_string.txt')
        return super().setUp()

    def test_uniform_output_string(self):
        circ_str = scg.surface_code_circuit_string(3, 2, .005)
        self.assertEqual(circ_str,self.uniform_string.read())

    def test_uniform_defect_output_string(self):
        circ_str = scg.surface_code_circuit_string(3, 2, .005, p_def=[0.1], def_coord=[(2,2)])
        self.assertEqual(circ_str,self.uniform_defect_string.read())

    def test_heterogeneous_output_string(self):
        circ_str = scg.surface_code_circuit_string_with_sigma_noise(3, 2, .005, 0.001)
        self.assertEqual(circ_str,self.heterogeneous_string.read())

    def test_heterogeneous_defect_output_string(self):
        circ_str = scg.surface_code_circuit_string_with_sigma_noise(3, 2, .005, 0.001, p_def=[0.5], def_coord=[(2,2)])
        self.assertEqual(circ_str,self.heterogeneous_defect_string.read())

    def tearDown(self):
        self.uniform_string.close()
        self.uniform_defect_string.close()
        return super().tearDown()

if __name__ == "__main__":
    unittest.main()