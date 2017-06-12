# Molecules, atoms definition
from __future__ import print_function
import numpy as np
import rxcclib.utils as cclibutils
import inspect,logging
from io import StringIO

class MoleDefError(Exception):
    def __init__(self,value):
        self.value=value
    def __str__(self):
        return repr(self.value)

class Molecule(object):
    def __init__(self,moleculename):
        if not isinstance(moleculename,str):
            logging.error('Error: Molecule name must be str')
            raise MoleDefError('Error: Molecule name must be str')

        self.name=moleculename
        self.__atomlist=[None]
        self.__bondlist={}
        self.__anglelist={}
        self.__dihdlist={}
        self.__bondfunclist={}
        self.__anglefunclist={}
        self.__dihdfunclist={}
        self.__pointer=0
    @property
    def bondfunclist(self):
        return self.__bondfunclist
    @property
    def anglefunclist(self):
        return self.__anglefunclist
    @property
    def dihdfunclist(self):
        return self.__dihdfunclist
    @property
    def natoms(self):
        return len(self.__atomlist)-1
    @property
    def bondlist(self):
        return self.__bondlist
    @property
    def anglelist(self):
        return self.__anglelist
    @property
    def dihdlist(self):
        return self.__dihdlist
    # Add atom and internal coordinates
    def addatom(self,idorsym,coords,unit='angstrom'):
        self.__atomlist.append(Atom(self,idorsym,coords,unit))

    def addbond(self,atomnum1,atomnum2):
        if atomnum1>atomnum2:
            atomnum1,atomnum2=atomnum2,atomnum1
        self.__bondlist.update({str(atomnum1)+','+str(atomnum2):Bond(self,atomnum1,atomnum2)})


    def addangle(self,atomnum1,atomnum2,atomnum3):
        if atomnum1>atomnum3:
            atomnum1,atomnum3=atomnum3,atomnum1
        self.__anglelist.update({str(atomnum1)+','+str(atomnum2)+','+str(atomnum3):Angle(self,atomnum1,atomnum2,atomnum3)})

    def adddihd(self,atomnum1,atomnum2,atomnum3,atomnum4):
        if atomnum2>atomnum3:
            atomnum2,atomnum3=atomnum3,atomnum2
            atomnum1,atomnum4=atomnum4,atomnum1
        elif atomnum2==atomnum3:
            if atomnum1>atomnum4:
                atomnum1,atomnum4=atomnum4,atomnum1
        self.__dihdlist.update({str(atomnum1)+','+str(atomnum2)+','+str(atomnum3)+','+str(atomnum4):Dihd(self,atomnum1,atomnum2,atomnum3,atomnum4)})

    # Get atom and internal coordiantes
    def atom(self,atomnum):
        return self.__atomlist[atomnum]

    def bond(self,atomnum1,atomnum2):
        if atomnum1>atomnum2:
            atomnum1,atomnum2=atomnum2,atomnum1
        return self.__bondlist[str(atomnum1)+','+str(atomnum2)]

    def angle(self,atomnum1,atomnum2,atomnum3):
        if atomnum1>atomnum3:
            atomnum1,atomnum3=atomnum3,atomnum1
        return self.__anglelist[str(atomnum1)+','+str(atomnum2)+','+str(atomnum3)]

    def dihd(self,atomnum1,atomnum2,atomnum3,atomnum4):
        if atomnum2>atomnum3:
            atomnum2,atomnum3=atomnum3,atomnum2
            atomnum1,atomnum4=atomnum4,atomnum1
        elif atomnum2==atomnum3:
            if atomnum1>atomnum4:
                atomnum1,atomnum4=atomnum4,atomnum1
        return self.__dihdlist[str(atomnum1)+','+str(atomnum2)+','+str(atomnum3)+','+str(atomnum4)]

    # Add MM functions
    def addbondfunc(self,bondobj):
        atype=bondobj[1].atomtype
        btype=bondobj[2].atomtype
        if atype>btype:
            atype,btype=btype,atype
        link=atype+' '+btype
        iffind=False
        for key,value in self.__bondfunclist.items():
            if key.link==link:
                value.append(bondobj)
                iffind=True
        if iffind==False:
            self.__bondfunclist.update({Bondfunc(self,bondobj):[bondobj]})
        bondobj.bondtype=link

    def addanglefunc(self,angleobj):
        atype=angleobj[1].atomtype
        btype=angleobj[2].atomtype
        ctype=angleobj[3].atomtype
        if atype>ctype:
            atype,ctype=ctype,atype
        link=atype+' '+btype+' '+ctype
        iffind=False
        for key,value in self.__anglefunclist.items():
            if key.link==link:
                value.append(angleobj)
                iffind=True
        if iffind==False:
            self.__anglefunclist.update({Anglefunc(self,angleobj):[angleobj]})
        angleobj.angletype=link

    def adddihdfunc(self,dihdobj):
        atype=dihdobj[1].atomtype
        btype=dihdobj[2].atomtype
        ctype=dihdobj[3].atomtype
        dtype=dihdobj[4].atomtype
        if btype>ctype:
            atype,dtype=dtype,atype
            btype,ctype=ctype,btype
        elif btype==ctype:
            if atype>dtype:
                atype,dtype=dtype,atype
        link=atype+' '+btype+' '+ctype+' '+dtype
        iffind=False
        for key,value in self.__dihdfunclist.items():
            if key.link==link:
                value.append(dihdobj)
                iffind=True
        if iffind==False:
            self.__dihdfunclist.update({Dihdfunc(self,dihdobj):[dihdobj]})
        dihdobj.dihdtype=link


    # getitem and iteration
    def __getitem__(self,key):
        if isinstance(key,int):
            return self.__atomlist[key]
        if isinstance(key,slice):
            start=key.start
            stop=key.stop
            if start is None:
                start=1
            L=[]
            for x in range(start,stop):
                L.append(self.__atomlist[x])
            return L
    def __iter__(self):
        for i in range(self.natoms):
            yield self[i+1]



    # Read structure from coords
    def readfromxyz(self,filelikeobj):
        f=filelikeobj
        for line in f:
            tmp=line.split()
            atom=tmp[0]
            if atom.isdigit():
                atom=int(atom)
            coords=np.array([tmp[1],tmp[2],tmp[3]],dtype=float)
            self.addatom(atom,coords,unit='angstrom')
    # Read connectivity
    def readconnectivity(self,conntystring):
        f=StringIO(conntystring)
        for line in f:
            try: tmp=line.split()
            except:
                logging.debug('End of connectivity, return.')
                return
            ite=iter(tmp)
            item0=next(ite)
            if not item0.isdigit():
                logging.debug('End of connectivity, return.')
                return
            a=int(item0)
            try:
                b=int(next(ite))
            except StopIteration:
                continue
            self.addbond(a,b)
            while True:
                try:
                    b=next(ite)
                    b=next(ite)
                    b=int(b)
                    self.addbond(a,b)
                except StopIteration:
                    break
        self.geninlcoords()
    def geninlcoords(self):
        angles=[]
        dihds=[]
        for atom1 in self:  # atom1: atom obj
            for atom2 in atom1.neighbor:  # atom2: atomnum
                if atom2==atom1.atomnum:
                    continue
                for atom3 in self[atom2].neighbor: #atom3: atomnum
                    if atom3==atom2 or atom3==atom1.atomnum:
                        continue
                    a=atom1.atomnum
                    b=atom2
                    c=atom3
                    if a>c:
                        a,c=c,a
                    angles.append(str(a)+'-'+str(b)+'-'+str(c))
                    for atom4 in self[atom3].neighbor:
                        if atom4==atom3 or atom4==atom2 or atom4==atom1.atomnum:
                            continue
                        a=atom1.atomnum
                        b=atom2
                        c=atom3
                        d=atom4
                        if b>c:
                            b,c=c,b
                            a,d=d,a
                        elif b==c:
                            if a>d:
                                a,d=d,a
                        dihds.append(str(a)+'-'+str(b)+'-'+str(c)+'-'+str(d))
        angles=list(set(angles))
        dihds=list(set(dihds))
        for item in angles:
            tmp=[int(x) for x in item.split('-')]
            self.addangle(*tmp)
        for item in dihds:
            tmp=[int(x) for x in item.split('-')]
            self.adddihd(*tmp)
    def readtypefromlist(self,L):
        if len(L)!=len(self.__atomlist):
            logging.error("Error when reading atomtype from list: length is not consistent with natoms")
            raise MoleDefError("Error when reading atomtype from list: length is not consistent with natoms")
        ite=iter(L)
        next(ite)
        for atom in self:
            atom.atomtype=next(ite)
    def readchargefromlist(self,L):
        if len(L)!=len(self.__atomlist):
            logging.error("Error when reading atomcharge from list: length is not consistent with natoms")
            raise MoleDefError("Error when reading atomcharge from list: length is not consistent with natoms")
        ite=iter(L)
        next(ite)
        for atom in self:
            atom.atomcharge=next(ite)



class Atom(object):
    '''
    Class "Atom" for atom object.
    Never add an atom directly. Always use molecule.addatom method instead, otherwise the molecule.__atomlist is not correct.

    >>> mole=molecule("CO")  # Define a molecule
    >>> mole.addatom(6,np.array([0.0,0.0,0.0]))
    >>> mole.addatom('O',np.array([0.0,0.0,2.0]))
    >>> mole.name
    'CO'
    >>> mole.atom(1).name
    'C1'
    >>> mole.atom(2).name
    'O2'
    >>> mole.atom(1).mymolecule.name
    'CO'
    >>> mole.atom(1).coords
    array([ 0.,  0.,  0.])
    >>> setattr(mole.atom(1),'atomtype','c2')
    >>> mole.atom(1).atomtype
    'c2'
    '''

    periotable=cclibutils.PeriodicTable()

    def __init__(self,mole,idorsym,coords,unit='angstrom'):  # molecule object,int,[float,float,float]
        # Assertion
        callername=inspect.stack()[1][3]
        assert callername=='addatom',"Atom must be added via Molecule.addatom method"
        assert isinstance(mole,Molecule),"First argument must be a molecule object!. Use Molecule.addatom method to avoid this problem."
        assert unit!='bohr' or unit!='angstrom', "Coordinate unit must be bohr or angstrom"


        self.__mymolecule=mole
        if isinstance(idorsym,int):
            self.elementid=idorsym
            try:
                self.atomsym=Atom.periotable.element[self.elementid] #str
            except KeyError:
                logging.critical("Error when adding atom: Idtosym not defined for atomic no:"+str(self.elementid))
                raise MoleDefError("Error when adding atom: Idtosym not defined for atomic no:"+str(self.elementid))
        elif isinstance(idorsym,str):
            self.atomsym=idorsym
            try:
                self.elementid=Atom.periotable.number[self.atomsym]
            except KeyError:
                logging.critical("Error when adding atom: Idtosym not defined for atomic symbol:"+self.atomsym)
                raise MoleDefError("Error when adding atom: Idtosym not defined for atomic symbol:"+self.atomsym)
        else:
            logging.critical("Error when adding atom: Expected atomic NO(int) or symbol(str) for input, received a"+str(type(idorsym)))
            raise MoleDefError("Error when adding atom: Expected atomic NO(int) or symbol(str) for input, received a"+str(type(idorsym)))



        if unit=='bohr':
            self.coords=cclibutils.convertor(coords,"bohr","Angstrom")
        elif unit=='angstrom':
            self.coords=coords

        self.atomnum=mole.natoms+1
        self.atomtype=self.name
        self.__neighbor=[]
        self.atomcharge=None

    @property
    def neighbor(self):
        return self.__neighbor
    def addneighbor(self,atomnum):
        if isinstance(atomnum,int):
            self.__neighbor.append(atomnum)
            self.__neighbor=list(set(self.__neighbor))
        else:
            logging.error("Error when adding neighbor: atomnum must be an integer.")
            raise MoleDefError("Error when adding neighbor: atomnum must be an integer.")
    def delneighbor(self,atomnum):
        if isinstance(atomnum,int):
            self.__neighbor.remove(atomnum)
        else:
            logging.error("Error when deleting neighbor: atomnum must be an integer.")
            raise MoleDefError("Error when deleting neighbor: atomnum must be an integer.")

    @property
    def name(self):
        return str(self.atomsym)+str(self.atomnum)
    @property
    def mymolecule(self):
        return self.__mymolecule
    def __str__(self):
        return "Atom object for atom "+self.name



class Bond(object):
    def __init__(self,mole,a,b): # self, atomid a, atomid b
        if a>b:
            a,b=b,a
        self.__a=mole[a]
        self.__b=mole[b]
        self.vec=self.__a.coords-self.__b.coords
        self.mytype=self.__a.name+' '+self.__b.name
        self.__a.addneighbor(b)
        self.parm=0.0
        self.__b.addneighbor(a)
    def __getitem__(self,value):
        if value==1:
            return self.__a
        elif value==2:
            return self.__b
        else:
            raise MoleDefError("Index for bond object must be 1 or 2.")
    @property
    def length(self):
        return np.linalg.norm(self.vec)
    def __str__(self):
        return "Bond object of bond "+self.__a.name+'-'+self.__b.name

class Angle(object):
    def __init__(self,mole,a,b,c):
        if a>c:
            a,c=c,a
        self.__a=mole[a]
        self.__b=mole[b]
        self.__c=mole[c]
        self.__ab=mole[a].coords-mole[b].coords
        self.__bc=mole[b].coords-mole[c].coords
        self.mytype=self.__a.name+' '+self.__b.name+' '+self.__c.name
        self.parm=0.0
    def __getitem__(self,value):
        if value==1:
            return self.__a
        elif value==2:
            return self.__b
        elif value==3:
            return self.__c
        else:
            raise MoleDefError("Index for angle object must be 1, 2 or 3.")

    @property
    def anglevalue(self):
        v1u=self.__ab/np.linalg.norm(self.__ab)
        v2u=self.__bc/np.linalg.norm(self.__bc)
        angle=180.0-np.arccos(np.dot(v1u,v2u))*180.0/np.pi
        if np.isnan(angle):
            if(v1u==v2u).all():
                return 0.0
            else:
                return 180.0
        return angle
    def __str__(self):
        return "Angle object of angle "+self.__a.name+'-'+self.__b.name+'-'+self.__c.name
class Dihd(object):
    def __init__(self,mole,a,b,c,d):
        if b>c:
            b,c=c,b
            a,d=d,a
        elif b==c:
            if a>d:
                a,d=d,a
        self.__a=mole[a]
        self.__b=mole[b]
        self.__c=mole[c]
        self.__d=mole[d]
        self.parm=0.0
        self.mytype=self.__a.name+' '+self.__b.name+' '+self.__c.name+' '+self.__d.name
    def __getitem__(self,value):
        if value==1:
            return self.__a
        elif value==2:
            return self.__b
        elif value==3:
            return self.__c
        elif value==4:
            return self.__d
        else:
            raise MoleDefError("Index for dihedral object must be 1, 2 ,3 or 4.")

    @property
    def dihdvalue(self):
        v1=self.__a.coords-self.__b.coords
        v2=self.__b.coords-self.__c.coords
        v3=self.__c.coords-self.__d.coords
        v1u=v1/np.linalg.norm(v1)
        v2u=v1/np.linalg.norm(v2)
        v3u=v1/np.linalg.norm(v3)
        n1=np.cross(v1u,v2u)
        n2=np.cross(v2u,v3u)
        dihd=np.arccos(np.dot(n1,n2))*180.0/np.pi
        if np.isnan(dihd):
            if (n1==n2).all():
                return 0.0
            else:
                return 180.0
        return dihd
    def __str__(self):
        return "Dihedral object of dihedral "+self.__a.name+'-'+self.__b.name+'-'+self.__c.name+'-'+self.__d.name

if __name__=='__main__':
    import doctest
    doctest.testmod()
