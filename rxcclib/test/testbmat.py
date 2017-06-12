#!/usr/bin/env python
import rxcclib.io.Gaussian as rxgau
import rxcclib.io.GauAmberCOM as gAmber
import rxcclib.geometry.molecules as rxmol
import rxcclib.utils as utils
from rxcclib.geometry.Bmat import Bmatrix
import numpy as np
import unittest
from io import StringIO


def findmatch(itnlcordL, inputfile):
    def readitnl(fileobj):
        intcords = []
        with open(fileobj.log.abspath) as f:
            tmp = f.read()
        tmp = tmp.split('Initial command')[-1]

        with StringIO(tmp) as f:
            for line in f:
                if line.find('Initial Parameters') < 0:
                    continue
                break
            [next(f) for x in range(0, 4)]
            for line in f:
                if line.find('---') >= 0:
                    break
                line = line.split()[2][2:-1]
                line = [int(x) for x in line.split(',')]
                intcords.append(" ".join([str(x) for x in line]))
        return intcords

    gauseq = readitnl(inputfile)

    for item in itnlcordL:
        atomset = []
        if type(item) == rxmol.Improper:
            a = item[1].atomnum
            b = item[2].atomnum
            c = item[3].atomnum
            d = item[4].atomnum
            if b > c:
                atomset = [d, c, b, a]
            else:
                atomset = [a, b, c, d]
        else:
            for atom in item:
                atomset.append(atom.atomnum)
        atomset = " ".join([str(x) for x in atomset])
        item.gauseq = gauseq.index(atomset)
    return


class TestBmat(unittest.TestCase):
    def setUp(self):
        self.file = rxgau.GauFile("samples/testBmat")
        gAmber.GauAmberCOM(self.file)
        self.file.com.read()
        self.file.fchk.read()
        self.mole = rxmol.Molecule('12CH')
        self.mole.readfromxyz(self.file.com.xyz)
        self.mole.readconnectivity(self.file.com.connectivity)

    def test_xyzfile(self):
        itnlcordL = []
        itnlcordL.extend(self.mole.dihedrallist.values())
        itnlcordL.extend(self.mole.anglelist.values())
        itnlcordL.extend(self.mole.bondlist.values())


        findmatch(itnlcordL,self.file)
        self.Bmat = Bmatrix(itnlcordL)
        self.Bmat.firstdiv()


        Btinv = np.linalg.pinv(self.Bmat.Bmat.transpose())
        gx = np.array(utils.toList(self.file.fchk.grads)) / 0.5291772086
        gq = np.dot(Btinv,gx)

        gqx=[]
        for item in itnlcordL:
            gqx.append([repr(item),self.file.fchk.intforces[item.gauseq]])
        for i,item in enumerate(gq):
            gqx[i].append(item)
        gqx=np.array(gqx)

        print(gq,len(gq))
        print(gqx,len(gqx))


if __name__ == '__main__':

    unittest.main()
