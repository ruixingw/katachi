from ..geometry import molecules as rxmol
from . import ffterms


class AmberFF(object):
    def __init__(self):
        self.bondlist = []
        self.anglelist = []
        self.dihedrallist = []
        self.improperlist = []

    def addTerms(self, itnlcordL):
        for item in itnlcordL:
            thetype = type(item)
            if thetype is rxmol.Dihedral:
                npaths = (len(item[1].neighbors) - 1) * (
                    len(item[4].neighbors) - 1)
                tmp = ffterms.AmbTrs(item, npaths)
                self.dihedrallist.append(tmp)
            elif thetype is rxmol.Angle:
                tmp = ffterms.HrmBnd1(item, 0.0, item.anglevalue)
                self.anglelist.append(tmp)
            elif thetype is rxmol.Bond:
                tmp = ffterms.HrmStr1(item, 0.0, item.length)
                self.bondlist.append(tmp)
            elif thetype is rxmol.Improper:
                tmp = ffterms.ImpTrs(item, 0.0, 2, 180.0)
                self.improperlist.append(tmp)
            else:
                raise
        return

    @property
    def terms(self):
        terms = []
        terms.extend(self.bondlist)
        terms.extend(self.anglelist)
        terms.extend(self.dihedrallist)
        terms.extend(self.improperlist)
        return terms

    def getEnergy(self):
        energy = 0
        for item in terms:
            energy += item.getEnergy()
        return energy

    def getGradient(self):
        grads = 0
        for item in terms:
            grads += item.getGradient()
        return grads

    def getHessian(self):
        hessian = 0
        for item in terms:
            hessian += item.getHessian()
        return hessian
