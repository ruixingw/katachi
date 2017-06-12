#!/usr/bin/env python3
# Construct B-matrix converting Cartesian Coordinates to Internal Coordinates
import numpy as np
import rxcclib.geometry.molecules as rxmol
import rxcclib.utils as utils


class BmatrixError(Exception):
    pass


class Bmatrix(object):
    def __init__(self, itnlcordL):
        self.mymolecule = itnlcordL[0][1].mymolecule
        self.itnlcordL = itnlcordL

    def firstdiv(self):
        Bmat = []
        for itnl in self.itnlcordL:
            Bt = []
            typeof = type(itnl)
            if typeof is rxmol.Bond:
                myatoms = [itnl[1], itnl[2]]
                v21 = itnl[1].coords - itnl[2].coords
                e21 = v21 / np.linalg.norm(v21)

                for atom in self.mymolecule:
                    if atom not in myatoms:
                        Bt.extend([0.0, 0.0, 0.0])
                    else:
                        if atom is itnl[1]:
                            Bt.extend(e21)
                        elif atom is itnl[2]:
                            Bt.extend(-e21)

                Bmat.append(Bt)
            elif typeof is rxmol.Angle:
                myatoms = [*itnl]
                v31 = itnl[1].coords - itnl[2].coords
                v32 = itnl[3].coords - itnl[2].coords
                e31 = v31 / np.linalg.norm(v31)
                e32 = v32 / np.linalg.norm(v32)

                s1 = (itnl.anglecos * e31 - e32) / (np.linalg.norm(v31) *
                                                    itnl.anglesin)
                s2 = (itnl.anglecos * e32 - e31) / (np.linalg.norm(v32) *
                                                    itnl.anglesin)
                for atom in self.mymolecule:
                    if atom not in myatoms:
                        Bt.extend([0.0, 0.0, 0.0])
                    else:
                        if atom is itnl[1]:
                            s = s1
                        elif atom is itnl[3]:
                            s = s2
                        elif atom is itnl[2]:
                            s = -s1 - s2
                        Bt.extend(s)

                Bmat.append(Bt)
            elif typeof is rxmol.Dihedral or typeof is rxmol.Improper:

                myatoms = [*itnl]
                phi1 = self.mymolecule.angle(
                    *[x.atomnum for x in myatoms[:-1]]).anglevalue
                phi2 = self.mymolecule.angle(
                    *[x.atomnum for x in myatoms[1:]]).anglevalue

                e12 = itnl[2].coords - itnl[1].coords
                e23 = itnl[3].coords - itnl[2].coords
                e32 = -e23
                e34 = itnl[4].coords - itnl[3].coords
                e43 = -e34

                n1 = np.cross(e12, e23)
                n2 = np.cross(e23, e34)
                n1 = n1/np.linalg.norm(n1)
                n2 = n2/np.linalg.norm(n2)

                l12 = np.linalg.norm(e12)
                l23 = np.linalg.norm(e23)
                l34 = np.linalg.norm(e43)

                e12u = e12 / np.linalg.norm(e12)
                e23u = e23 / np.linalg.norm(e23)
                e32u = -e23u
                e34u = e34 / np.linalg.norm(e34)
                e43u = -e34u

                n1223 = np.cross(e12u, e23u)
                n1223 = n1223 / np.linalg.norm(n1223)
                n4332 = np.cross(e43u, e32u)
                n4332 = n4332 / np.linalg.norm(n4332)

                sin123 = np.sin(phi1 * np.pi / 180)
                cos123 = np.cos(phi1 * np.pi / 180)
                sin234 = np.sin(phi2 * np.pi / 180)
                cos234 = np.cos(phi2 * np.pi / 180)

                for atom in self.mymolecule:
                    if atom not in myatoms:
                        Bt.extend([0.0, 0.0, 0.0])
                    else:
                        if atom is itnl[1]:
                            s1 = -n1223 / (l12 * sin123)
                            #print(s1)
                            one = np.identity(3, dtype=int)
                            minusone = -one

                            n1pp = np.cross(minusone, e23)

                            cos = itnl.anglecos*(
                                np.dot(
                                    np.cross(n1pp, e23u), n2)
                            )
                            sin = itnl.anglesin*(
                                np.dot(n1pp, n2)
                            )
                            s1 = cos - sin
                            #print(s1)
                            s = s1
                        elif atom is itnl[2]:
                            s2 = (n1223 / (l12 * sin123) - cos123 * n1223 /
                                  (l23 * sin123) + cos234 * n4332 /
                                  (l23 * sin234))
                            s = s2
                        elif atom is itnl[3]:
                            s3 = (n4332 / (l34 * sin234) - cos234 * n4332 /
                                  (l23 * sin234) + cos123 * n1223 /
                                  (l23 * sin123))
                            s = s3
                        elif atom is itnl[4]:
                            s4 = -n4332 / (l34 * sin234)
                            s = s4
                        else:
                            # Should never happen
                            raise BmatrixError(
                                'Unexpected atom in dihd of Bmatrix.construct')
                        Bt.extend(s)
                Bmat.append(Bt)
            else:
                # Should never happen
                raise BmatrixError("Unexpected itnl type: " + str(type(itnl)))

        self.Bmat = np.array(Bmat)
        return self.Bmat

    def renew(self):
        del self.Bmat
        self.firstdiv()
        return


if __name__ == '__main__':
    initset()
    finalfuncL, itnlcordL, unkitnlL = readgeom()

    a = Bmatrix(itnlcordL)
    a.firstdiv()
    np.set_printoptions(edgeitems=3, linewidth=45)
    #a.itnlsympy()
    res = 0
    #for i in range(0,12):
    #    res += np.dot(a.Bmat[0][i], a._deltax[i])
