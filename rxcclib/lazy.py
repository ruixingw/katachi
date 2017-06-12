import rxcclib.geometry.molecules as rxmol
import rxcclib.file.Gaussian as rxgau
import rxcclib.utils as utils
import logging
from io import StringIO


def fchkToGeom(filename):
    filename = filename.split('.fchk')[0]
    filename = filename.split('.fch')[0]
    fileobj = rxgau.GauFile(filename)
    fileobj.fchk.read()
    geom = rxmol.Molecule('mole')
    geom.addatomsFromLists(fileobj.fchk.atomnos, fileobj.fchk.atomcoords[-1])


def ConnectivityToConnectionMatrix(constring):
    lines = constring.split('\n')[:-1]
    dimension = len(lines)
    L=[]
    for i,line in enumerate(lines):
        tmp = line.split()
        if i+1 !=int(tmp[0]):
            raise
        thislinenum = dimension-i
        tmpL = [0]*thislinenum
        if len(tmp) == 1:
            L.extend(tmpL)
        else:
            for no,order in zip(tmp[1::2],tmp[2::2]):
                no=int(no)
                order = float(order)
                tmpL[no-i-1] = order
            L.extend(tmpL)
    return L


# Read connectivity to complete neighbor info
def readconnectivity(self, conntystring):
    f = StringIO(conntystring)
    for line in f:
        tmp = line.split()
        if not tmp:
            logging.debug('End of connectivity, return.')
            return
        ite = iter(tmp)
        item0 = next(ite)
        if not item0.isdigit():
            logging.debug('End of connectivity, return.')
            return
        a = int(item0)
        try:
            b = int(next(ite))
        except StopIteration:
            continue
        self.addbond(a, b)
        while True:
            try:
                b = next(ite)
                b = next(ite)
                b = int(b)
                self.addbond(a, b)
            except StopIteration:
                break
    angles = []
    dihedrals = []
    for atom1 in self:  # atom1: atom obj
        for atom2 in atom1.neighbor:  # atom2: obj
            if atom2 is atom1:
                continue
            for atom3 in atom2.neighbor:  # atom3: obj
                if (atom3 is atom2 or atom3 is atom1):
                    continue
                a = atom1.atomnum
                b = atom2.atomnum
                c = atom3.atomnum
                if a > c:
                    a, c = c, a
                angles.append(str(a) + '-' + str(b) + '-' + str(c))
                for atom4 in atom3.neighbor:
                    if (atom4 is atom3 or atom4 is atom2 or atom4 is atom1):
                        continue
                    a = atom1.atomnum
                    b = atom2.atomnum
                    c = atom3.atomnum
                    d = atom4.atomnum
                    if b > c:
                        b, c = c, b
                        a, d = d, a
                    dihedrals.append(
                        str(a) + '-' + str(b) + '-' + str(c) + '-' + str(d))
    angles = list(set(angles))
    dihedrals = list(set(dihedrals))
    for item in angles:
        tmp = [int(x) for x in item.split('-')]
        self.addangle(*tmp)
    for item in dihedrals:
        tmp = [int(x) for x in item.split('-')]
        self.adddihedral(*tmp)


def readtypefromlist(self, L):
    if len(L) != len(self._atomlist):
        print(L, self._atomlist)
        logging.error("Error when reading atomtype from list:"
                      " length is not consistent with natom")
        raise GeometryError

    ite = iter(L)
    for atom in self:
        atom.atomtype = next(ite)


def readchargefromlist(self, L):
    if len(L) != len(self._atomlist):
        logging.error("Error when reading atomcharge from list:"
                      " length is not consistent with natom")
        raise GeometryError

    ite = iter(L)
    for atom in self:
        atom.atomcharge = next(ite)
