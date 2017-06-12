#!/usr/bin/env python3
import rxcclib.io.Gaussian as rxgau
import rxcclib.io.mol2 as rxmol2file
import rxcclib.geometry.molecules as rxmol
from rxcclib.ff.mol2 import Mol2
import subprocess
import unittest, os, logging
import numpy as np
from io import StringIO


rxgau.GauCOM.g09rt = 'g09'
os.system(
    'rm A* q* Q* p* esout *Gaussian* samples/bencom.fchk samples/bencom.chk samples/bencom.log'
)


class TestFile(unittest.TestCase):
    def test_comfchk(self):
        file = rxgau.GauFile('samples/bencom')
        Mol2(file)
        self.assertIsInstance(file, rxgau.GauFile)
        self.assertIsInstance(file.com, rxgau.GauCOM)
        self.assertIsInstance(file.log, rxgau.GauLOG)
        self.assertIsInstance(file.fchk, rxgau.GauFCHK)
        self.assertIsInstance(file.mol2, Mol2)
        file.com.Popen()
        file.com.wait()
        file.chk.formchk()

        self.assertEqual(file.fchk.read(), True)
        self.assertEqual(file.fchk.natom, 12)
        self.assertEqual(file.fchk.mult, 1)
        self.assertEqual(file.fchk.charge, 0)
        self.assertIsInstance(file.fchk.xyz, str)

        hess = file.fchk.find33Hessian(3, 5)
        self.assertAlmostEqual(hess[0][0], -2.62909045e-2)
        self.assertAlmostEqual(hess[1][1], 3.38743754e-2)
        self.assertAlmostEqual(hess[2][2], 7.19580040e-3)

    def test_logmol2(self):
        file = rxgau.GauFile('samples/benresp')
        Mol2(file)
        args = 'antechamber -i {} -fi gout -o {} -fo mol2 -c resp'.format(
            file.log.abspath, file.mol2.abspath)
        args = args.split()
        run = subprocess.run(args)
        file.mol2.read()
        self.assertEqual(file.mol2.atomtypelist[0], None)
        self.assertEqual(file.mol2.atomtypelist[1], 'ca')
        self.assertEqual(file.mol2.atomtypelist[12], 'ha')
        self.assertEqual(file.mol2.atomchargelist[0], None)
        self.assertEqual(file.mol2.atomchargelist[1], -0.117738)
        self.assertEqual(file.mol2.atomchargelist[12], 0.117738)

    # def test_MMcom(self):
    #     mmfile=rxgau.File('samples/mmfile')
    #     mmfile.com.read()
    #     xyz=StringIO(mmfile.com.xyz)

    #        self.assertEqual()


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    unittest.main()
