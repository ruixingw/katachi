#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import rxcclib.geometry.molecules as rxmol
import rxcclib.file.Gaussian as rxgau
import rxcclib.file.GauAmberCOM as gAmber
import unittest
from io import StringIO


class TestConnect(unittest.TestCase):
    def test_xyzfile(self):
        benzene = rxmol.Molecule("benzene")
        benfile = rxgau.GauFile("samples/bencom")
        benfile.fchk.read()
        benzene.addatomsFromLists(benfile.fchk.atomnos,benfile.fchk.atomcoords[-1])
        self.assertEqual(benzene[1].name, 'C1')
        self.assertEqual(benzene[12].name, 'H12')

    def testmmfile(self):
        mmfile = rxgau.GauFile('samples/mmfile')
        gAmber.GauAmberCOM(mmfile)
        mmfile.com.read()
        benmol = rxmol.Molecule('benzene')
        benmol.readfromxyz(mmfile.com.xyz)
        #benmol.readtypefromlist(mmfile.com.atomtypelist)
        #benmol.readchargefromlist(mmfile.com.atomchargelist)
        for atom in benmol:
            print(atom, atom.atomtype, atom.atomcharge)
        benmol.readconnectivity(mmfile.com.connectivity)
        print("Bond: ", len(benmol.bondlist.keys()))
        print("Angle: ", len(benmol.anglelist.keys()))
        print("Dihd: ", len(benmol.dihedrallist.keys()))


if __name__ == '__main__':
    unittest.main()
