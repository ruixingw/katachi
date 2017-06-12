# Molecules, atoms definition
import logging
import copy
from io import StringIO
import numpy as np
import rxcclib.utils as utils

periodictable = utils.PeriodicTable()


class GeometryError(Exception):
    def __init__(self):
        pass


class Molecule(object):
    def __init__(self, name):
        if not isinstance(name, str):
            logging.error('Error: Molecule name must be str')
            raise GeometryError

        self.name = name
        self._atomlist = []
        self._bondlist = {}
        self._anglelist = {}
        self._dihedrallist = {}
        self._improperlist = {}

    @property
    def atomcoords(self):
        L = []
        for atom in self.atomlist:
            L.extend(atom.coords)
        return L

    @property
    def natom(self):
        return len(self._atomlist)

    @property
    def atomlist(self):
        return self._atomlist

    @property
    def bondlist(self):
        return self._bondlist

    @property
    def anglelist(self):
        return self._anglelist

    @property
    def dihedrallist(self):
        return self._dihedrallist

    @property
    def improperlist(self):
        return self._improperlist

    # Add atom and internal coordinates

    def addatom(self, atomnoORelement, coords, unit='angstrom'):
        tmpatom = Atom(self, atomnoORelement, coords, unit)
        self._atomlist.append(tmpatom)

        return tmpatom

    def addbond(self, atomnum1, atomnum2, order=1):
        if atomnum1 > atomnum2:
            atomnum1, atomnum2 = atomnum2, atomnum1

        tmpbond = Bond(self, atomnum1, atomnum2, order)
        key = "{}-{}".format(atomnum1, atomnum2)
        self._bondlist[key] = tmpbond
        return tmpbond

    def addangle(self, atomnum1, atomnum2, atomnum3):
        if atomnum1 > atomnum3:
            atomnum1, atomnum3 = atomnum3, atomnum1

        tmpangle = Angle(self, atomnum1, atomnum2, atomnum3)
        key = "{}-{}-{}".format(atomnum1, atomnum2, atomnum3)
        self._anglelist[key] = tmpangle


        return tmpangle

    def adddihedral(self, atomnum1, atomnum2, atomnum3, atomnum4):
        if atomnum2 > atomnum3:
            atomnum2, atomnum3 = atomnum3, atomnum2
            atomnum1, atomnum4 = atomnum4, atomnum1
        elif atomnum2 == atomnum3:
            if atomnum1 > atomnum4:
                atomnum1, atomnum4 = atomnum4, atomnum1

        tmpdihedral = Dihedral(self, atomnum1, atomnum2, atomnum3, atomnum4)
        key = "{}-{}-{}-{}".format(atomnum1, atomnum2, atomnum3, atomnum4)
        self._dihedrallist[key] = tmpdihedral

        return tmpdihedral

    def addimproper(self, atomnum1, atomnum2, atomnum3, atomnum4):
        tmpimproper = Improper(self, atomnum1, atomnum2, atomnum3, atomnum4)
        key = "{}-{}-{}-{}".format(atomnum1, atomnum2, atomnum3, atomnum4)
        self.improperlist[key] = tmpimproper

        return tmpimproper

    # Get atom and internal coordiantes
    def atom(self, atomnum):
        return self._atomlist[atomnum - 1]

    def bond(self, atomnum1, atomnum2):
        if type(atomnum1) is Atom and type(atomnum2) is Atom:
            return self.bond(atomnum1.atomnum, atomnum2.atomnum)

        if atomnum1 > atomnum2:
            atomnum1, atomnum2 = atomnum2, atomnum1
        key = "{}-{}".format(atomnum1, atomnum2)
        return self._bondlist[key]

    def angle(self, atomnum1, atomnum2, atomnum3):
        if type(atomnum1) is Atom and type(atomnum2) is Atom and type(
                atomnum3) is Atom:
            return self.angle(atomnum1.atomnum, atomnum2.atomnum,
                              atomnum3.atomnum)

        if atomnum1 > atomnum3:
            atomnum1, atomnum3 = atomnum3, atomnum1
        key = "{}-{}-{}".format(atomnum1, atomnum2, atomnum3)
        return self._anglelist[key]

    def dihedral(self, atomnum1, atomnum2, atomnum3, atomnum4):
        if type(atomnum1) is Atom and type(atomnum2) is Atom and type(
                atomnum3) is Atom and type(atomnum4) is Atom:
            return self.dihedral(atomnum1.atomnum, atomnum2.atomnum,
                                 atomnum3.atomnum, atomnum4.atomnum)

        if atomnum2 > atomnum3:
            atomnum2, atomnum3 = atomnum3, atomnum2
            atomnum1, atomnum4 = atomnum4, atomnum1
        key = "{}-{}-{}-{}".format(atomnum1, atomnum2, atomnum3, atomnum4)
        return self._dihedrallist[key]

    def improper(self, atomnum1, atomnum2, atomnum3, atomnum4):
        if type(atomnum1) is Atom and type(atomnum2) is Atom and type(
                atomnum3) is Atom and type(atomnum4) is Atom:
            return self.improper(atomnum1.atomnum, atomnum2.atomnum,
                                 atomnum3.atomnum, atomnum4.atomnum)

        key = "{}-{}-{}-{}".format(atomnum1, atomnum2, atomnum3, atomnum4)
        return self._improperlist[key]

    # slice
    def __getitem__(self, key):
        tmplist = copy.copy(self._atomlist)
        tmplist.insert(0, None)
        return tmplist.__getitem__(key)

    # iterate atoms
    def __iter__(self):
        return iter(self._atomlist)

    def __len__(self):
        return self.natom

    # Read structure from xyz
    def addatomsFromXYZ(self, string):
        f = StringIO(string)
        for line in f:
            tmp = line.split()
            atom = tmp[0]
            if atom.isdigit():
                atom = int(atom)
            coords = [tmp[1], tmp[2], tmp[3]]
            self.addatom(atom, coords, unit='angstrom')
        return

    def addatomsFromLists(self, atomlist, coordslist):
        coords = utils.atomwiseList.toGroupOfThree(coordslist)
        zipped = zip(atomlist, coordslist)
        for x, y in zipped:
            self.addatom(x, y, unit='angstrom')
        return

    def readConnectionMatrix(self, conmat):
        for i in range(self.natom):
            # row
            for j in range(i+1):
                #column
                if conmat[i][j] != 0:
                    self.addbond(self[i + 1].atomnum, self[j + 1].atomnum, order=conmat[i][j])
        return


class Atom(object):
    def __init__(self, mole, atomnoORelement, coords,
                 unit='angstrom'):  # molecule object,int,[float,float,float]
        """
        An atom object should not be directly created, should use addatom method of a Molecule instance.
        """

        self.mymolecule = mole
        self.neighbors = []
        if isinstance(atomnoORelement, int):
            self._atomno = atomnoORelement
            try:
                self._element = periodictable.element[self._atomno]
            except KeyError:
                logging.critical("Error when adding atom:"
                                 "No element for atomno:" + str(self._atomno))
                raise GeometryError

        elif isinstance(atomnoORelement, str):
            self._element = atomnoORelement
            try:
                self._atomno = periodictable.number[self._element]
            except KeyError:
                logging.critical("Error when adding atom: "
                                 " No atomic symbol:" + self._element)
                raise GeometryError
        else:
            logging.critical("Error when adding atom: Expected atomic"
                             " NO(int) or element(str) for input, received a" +
                             str(type(atomnoORelement)))
            raise GeometryError

        if unit == 'bohr':
            coords = [utils.convertor(x, "bohr", "Angstrom") for x in coords]
        tmp = type(coords)
        if tmp is list or tmp is tuple:
            self.coords = np.array(coords, dtype=np.float64)
        elif tmp is np.ndarray:
            self.coords = coords

    @property
    def atomno(self):
        return self._atomno

    @property
    def element(self):
        return self._element

    @property
    def atomnum(self):
        return self.mymolecule.atomlist.index(self) + 1

    @property
    def name(self):
        return "{}{}".format(self.element, self.atomnum)

    def __repr__(self):
        return "Atom " + self.name

    __str__ = __repr__


class Bond(object):
    def __init__(self, mole, a, b, order=1):  # self, atomnum a, atomnum b
        if a > b:
            a, b = b, a
        a = mole[a]
        b = mole[b]
        self._atomlist = [a, b]
        self.mymolecule = mole
        self.order = order
        a.neighbors.append(b)
        a.neighbors = list(set(a.neighbors))
        b.neighbors.append(a)
        b.neighbors = list(set(b.neighbors))
        self.vec = self[1].coords - self[2].coords

    def __getitem__(self, key):
        if key != 1 and key != 2:
            logging.critical('key for Bond[key] must be 1 or 2.')
            raise GeometryError
        return self._atomlist[key - 1]

    # iterate atoms
    def __iter__(self):
        return iter(self._atomlist)

    @property
    def length(self):
        return np.linalg.norm(self.vec)

    def __repr__(self):
        return "Bond " + self[1].name + '-' + self[2].name

    __str__ = __repr__


class Angle(object):
    def __init__(self, mole, a, b, c):
        if a > c:
            a, c = c, a
        a = mole[a]
        b = mole[b]
        c = mole[c]
        self._atomlist = [a, b, c]
        self.mymolecule = mole

    def __getitem__(self, key):
        if isinstance(key, int):
            if key != 1 and key != 2 and key != 3:
                logging.critical('key for Angle[key] must be 1, 2 or 3.')
                raise GeometryError
        return self._atomlist[key - 1]

    # iterate atoms
    def __iter__(self):
        return iter(self._atomlist)

    @property
    def anglevalue(self):
        v1 = self[1].coords - self[2].coords
        v2 = self[2].coords - self[3].coords
        v1u = v1 / np.linalg.norm(v1)
        v2u = v2 / np.linalg.norm(v2)
        angle = 180.0 - np.arccos(np.dot(v1u, v2u)) * 180.0 / np.pi
        if np.isnan(angle):
            if (v1u == v2u).all():
                return 0.0
            else:
                return 180.0
        return angle

    @property
    def anglesin(self):
        return np.sin(self.anglevalue * np.pi / 180)

    @property
    def anglecos(self):
        return np.cos(self.anglevalue * np.pi / 180)

    def __repr__(self):
        return (
            "Angle " + self[1].name + '-' + self[2].name + '-' + self[3].name)

    __str__ = __repr__


class Dihedral(object):
    def __init__(self, mole, a, b, c, d):
        if b > c:
            b, c = c, b
            a, d = d, a
        elif b == c and a > d:
            a, d = d, a
        a = mole[a]
        b = mole[b]
        c = mole[c]
        d = mole[d]
        self._atomlist = [a, b, c, d]
        self.mymolecule = mole

    def __getitem__(self, key):
        if isinstance(key, int):
            if key != 1 and key != 2 and key != 3 and key != 4:
                logging.critical('key for Dihedral[key] must be 1, 2, 3 or 4.')
                raise GeometryError
        return self._atomlist[key - 1]

    # iterate atoms
    def __iter__(self):
        return iter(self._atomlist)

    @property
    def anglesin(self):
        return np.sin(self.anglevalue * np.pi / 180)

    @property
    def anglecos(self):
        return np.cos(self.anglevalue * np.pi / 180)

    @property
    def anglevalue(self):
        # arctan2 is better than arccos
        # see http://math.stackexchange.com/a/47084
        # and http://stackoverflow.com/a/12011762
        # if someday this code face performance problem
        # try "new_dihedral" in this post:
        # http://stackoverflow.com/a/34245697/5677043
        v1 = self[1].coords - self[2].coords
        v2 = self[2].coords - self[3].coords
        v3 = self[3].coords - self[4].coords
        v1u = v1 / np.linalg.norm(v1)
        v2u = v2 / np.linalg.norm(v2)
        v3u = v3 / np.linalg.norm(v3)
        n1 = np.cross(v1u, v2u)
        n2 = np.cross(v2u, v3u)
        n1 = n1 / np.linalg.norm(n1)
        n2 = n2 / np.linalg.norm(n2)
        m1 = np.cross(n1, v2u)
        x = np.dot(n1, n2)
        y = np.dot(m1, n2)
        dihedral = np.arctan2(y, x) * 180.0 / np.pi
        return dihedral

    def __repr__(self):
        return ("Dihedral " + self[1].name + '-' + self[2].name + '-' +
                self[3].name + '-' + self[4].name)

    __str__ = __repr__


class Improper(Dihedral):
    def __init__(self, mole, a, b, c, d):
        a = mole[a]
        b = mole[b]
        c = mole[c]
        d = mole[d]
        self._atomlist = [a, b, c, d]
        self.mymolecule = mole

    def __repr__(self):
        return ("Improper " + self[1].name + '-' + self[2].name + '-' +
                self[3].name + '-' + self[4].name)

    __str__ = __repr__
