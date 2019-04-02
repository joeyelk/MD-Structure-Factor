#! /usr/bin/env python

import unittest
import os
import sys

path_to_tests=str(os.path.dirname(os.path.realpath(__file__)))
path_to_tests+="/../"
sys.path.insert(0,path_to_tests)

# Import local files
import load_traj

class TestLoadingTrajectoryFiles(unittest.TestCase):

    def test_load_gro(self):
        gro_file_path=str(os.path.dirname(os.path.realpath(__file__)))
        gro_file_path+="/water.gro"
        atom_names = load_traj.load_gro(gro_file_path)
        self.assertEqual(atom_names[0],"OW1")
        self.assertEqual(atom_names[1],"HW2")
        self.assertEqual(atom_names[2],"HW3")
        self.assertEqual(atom_names[3],"OW1")
        self.assertEqual(atom_names[4],"HW2")
        self.assertEqual(atom_names[5],"HW3")

if __name__ == '__main__':
    unittest.main()
