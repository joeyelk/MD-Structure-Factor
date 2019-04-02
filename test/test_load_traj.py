#! /usr/bin/env python

import unittest
import sys
sys.path.insert(0,'../')

" Import local files "
import load_traj

class TestLoadingTrajectoryFiles(unittest.TestCase):

    def test_load_gro(self):
        atom_names = load_traj.load_gro("water.gro")
        self.assertEqual(atom_names[0],"OW1")
        self.assertEqual(atom_names[1],"HW2")
        self.assertEqual(atom_names[2],"HW3")
        self.assertEqual(atom_names[3],"OW1")
        self.assertEqual(atom_names[4],"HW2")
        self.assertEqual(atom_names[5],"HW3")

if __name__ == '__main__':
    unittest.main()
