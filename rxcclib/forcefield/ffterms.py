import numpy as np

class HrmStr1(object):
    """
    Harmonic stretch I (Amber 1): ForceC*(distance-eqvalue)^2
    """

    def __init__(self, bondobj, ForceC, eqvalue):
        """
        Float ForceC: Force constant in Kcal mol-1 Angstrom-2
        Float distance: Bond distance in Angstrom
        Float eqvalue: Equilibrium value in Angstrom
        """
        self.bond = bondobj
        self.ForceC = ForceC
        self.eqvalue = eqvalue

    @property
    def distance(self):
        return self.bond.length

    def getEnergy(self):
        energy = self.ForceC * (self.distance - self.eqvalue)**2
        return energy


class HrmBnd1(object):
    """
    Harmonic bend (Amber 1): ForceC*(anglevalue-eqvalue)^2
    """

    def __init__(self, angleobj, ForceC, eqvalue):
        """
        Float ForceC: Force constant in Kcal mol-1 rad-2
        Float anglevalue: Angle in degree
        Float eqvalue: Equilibrium value in degree
        """
        self.angle = angleobj
        self.ForceC = ForceC
        self.eqvalue = eqvalue

    @property
    def anglevalue(self):
        return self.angle.anglevalue

    def getEnergy(self):
        diff = np.radians(self.anglevalue - self.eqvalue)
        energy = self.ForceC * diff**2
        return energy


class AmbTrs(object):

    """
    Amber torsion (Amber 1)
    Master function
    Sum_{i}{Term_i}/npaths
    """

    class Term(object):
        """
        Amber torsion (Amber 1)
        Individual terms
        ForceC*[1+cos(i*anglevalue-phase_i)]
        """

        def __init__(self, masterobj, ForceC, periodicity, phase):
            """
            Float ForceC: V/2 magnitude
            Int i: Periodicity
            Float anglevalue: Dihedral angle value
            Float phase_i: Phase offset
            """
            self.master = masterobj
            self.ForceC = ForceC
            self.periodicity = periodicity
            self.phase = phase

        def getEnergy(self):
            i = self.periodicity
            phi = self.master.anglevalue
            po = self.phase
            k = self.ForceC
            energy = k * [1 + np.cos(np.radians(i * phi - po))]
            return energy

    def __init__(self, dihedralobj, npaths):
        """
        Int npaths: Divider
        """
        self.dihedral = dihedralobj
        self.terms = {}
        self.npaths = npaths

    @property
    def anglevalue(self):
        return self.dihedral.anglevalue

    def addTerm(self, ForceC, periodicity, phase):
        term = Term(ForceC, periodicity, phase)
        self.terms[periodicity] = term
        return term

    def getEnergy(self):
        tmp = [x.getEnergy() for x in self.terms.values()]
        return sum(tmp)


class ImpTrs(AmbTrs.Term):
    def __init__(self, improperobj, ForceC, periodicity, phase):
        super().__init__(improperobj, ForceC, periodicity, phase)
        self.npaths = 2
        self.improper = self.master

    @property
    def anglevalue(self):
        return self.improper.anglevalue

