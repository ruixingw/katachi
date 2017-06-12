# From fchk read Charge, Multiplicity, Coordinates
from __future__ import print_function
import os,time,logging
import numpy as np
from io import StringIO
import molecules as rxmol
import rxcclib.utils as cclibutils

class rxccError(Exception):
    def __init(self,value):
        self.value=value
    def __repr__(self):
        return repr(self.value)
def tryfunc(func):
    def wrapper(*args,**kw):
        try:
            func(*args,**kw)
        except Exception as ex:
            raise ex
            return False
        return True
    return wrapper
class File(object):
    def __init__(self,name):
        pwd=os.path.abspath('.')
        self.name=os.path.join(pwd,name)  # Name with absolute path

        self.__com=gauCOM(self)
        self.__log=gauLOG(self)
        self.__fchk=gauFCHK(self)
        self.__ac=amberAC(self)
        self.default='fchk'
        self.__natoms=None
        self.__mlpty=None
        self.__totalcharge=None
    # subfile objects
    @property
    def comname(self):
        return self.name+'.com'
    @property
    def logname(self):
        return self.name+'.log'
    @property
    def chkname(self):
        return self.name+'.chk'
    @property
    def fchkname(self):
        return self.name+'.fchk'
    @property
    def acname(self):
        return self.name+'.ac'

    @property
    def com(self):
        return self.__com
    @property
    def log(self):
        return self.__log
    @property
    def fchk(self):
        return self.__fchk
    @property
    def ac(self):
        return self.__ac

    # molecule properties; setter should only be called by subfile methods
    @property
    def natoms(self):
        return self.__natoms
    @natoms.setter
    def natoms(self,value):
        if self.__natoms==None:
            self.__natoms=value
        else:
            if self.__natoms!=value:
                raise rxccError("Error: natoms is already read, and not consistent with new value: Current value is "+str(self.__natoms)+", New value is "+str(value))

    @property
    def multiplicity(self):
        return self.__mlpty
    @multiplicity.setter
    def multiplicity(self,value):
        if self.__mlpty==None:
            self.__mlpty=value
        else:
            if self.__mlpty!=value:
                raise rxccError("Error: multiplicity is already read, and not consistent with new value: Current value is "+str(self.__mlpty)+", New value is "+str(value))
    @property
    def totalcharge(self):
        return self.__totalcharge
    @totalcharge.setter
    def totalcharge(self,value):
        if self.__totalcharge==None:
            self.__totalcharge=value
        else:
            if self.__totalcharge!=value:
                raise rxccError("Error: totalcharge is already read, and not consistent with new value: Current value is "+str(self.__totalcharge)+", New value is "+str(value))
    @property
    def xyzfile(self):
        souc='fchk'
        if self.default!='fchk':
            souc=self.default
        if souc=='fchk':
            return self.__fchk.xyz
        elif souc=='com':
            return self.__com.xyz
    @property
    def atomtypelist(self):
        return self.__ac.atomtypelist
    @property
    def atomchargelist(self):
        return self.__ac.atomchargelist

    def runformchk(self):
        string='formchk '+self.chkname+' '+self.fchkname
        logging.info('  '+string)
        iferror=os.popen(string)
        if iferror.read().find('Error')>=0:
            raise rxccError('   Error in formatting'+self.chkname)
            iferror.close()
            return False
        iferror.close()
        return True



class gauFCHK(object):
    def __init__(self,father):
        self.__father=father
        self.filename=self.__father.fchkname
        self.coordslist=[]
        self.atomlist=[None]
        self.readstate=None
        self.totalcharge=None
        self.multiplicity=None
        self.natoms=None
        self.hessian=[]
        self.xyz=''
    @tryfunc
    def read(self):
        if self.readstate==True:
            logging.warning("fchk.read(): fchk is already read")
            return True
        logging.info('Read fchk:'+self.__father.fchkname)
        # FCHK parser
        with open(self.__father.fchkname,'r') as f:
            string=next(f)
            for string in f:
                if string.find('Charge')==0:
                    self.totalcharge=int(string.split('I')[1])
                    self.__father.totalcharge=self.totalcharge
                if string.find('Multiplicity')==0:
                    self.multiplicity=int(string.split('I')[1])
                    self.__father.multiplicity=self.multiplicity
                if string.find('Atomic numbers')==0:
                    self.natoms=int(string.split('=')[1])
                    self.__father.natoms=self.natoms
                    string=next(f)
                    while string.find('Nuclear charges')<0:
                        self.atomlist.extend([int(x) for x in string.split()])
                        string=next(f)

                if string.find('Current cartesian coordinates')==0:
                    string=next(f)
                    while string.find('Force Field')<0:
                        self.coordslist.extend([float(x) for x in string.split()])
                        string=next(f)
                # Read Hessian
                if string.find('Cartesian Force Constants')==0:
                    string=next(f)
                    while string.find('Dipole')<0:
                        self.hessian.extend([float(x) for x in string.split()])
                        string=next(f)
                #Stop
        self.coordslist=[cclibutils.convertor(x,"bohr","Angstrom") for x in self.coordslist]
        self.coordslist=np.array(self.coordslist)
        self.atomlist=np.array(self.atomlist)
        self.hessian=np.array(self.hessian)
        for i in range(0,len(self.atomlist)-1):
            tmp=str(self.atomlist[i+1])+'   '+str(self.coordslist[3*i])+'   '+str(self.coordslist[3*i+1])+'   '+str(self.coordslist[3*i+2])+'\n'
            self.xyz+=tmp
        return True

    def findHessianElement(self,i,j): #i, j: coordinate number
        if i<j:
            i,j=j,i
        num=i*(i-1)/2+j
        num=int(num)
        return self.hessian[num-1]
    def find33Hessian(self,i,j):    #i, j: atom number
        if i<j:
            i,j=j,i
        tthess=[]
        i1=3*(i-1)+1
        i2=3*(i-1)+2
        i3=3*i
        j1=3*(j-1)+1
        j2=3*(j-1)+2
        j3=3*j
        tthess.append([self.findHessianElement(i1,j1),self.findHessianElement(i1,j2),self.findHessianElement(i1,j3)])
        tthess.append([self.findHessianElement(i2,j1),self.findHessianElement(i2,j2),self.findHessianElement(i2,j3)])
        tthess.append([self.findHessianElement(i3,j1),self.findHessianElement(i3,j2),self.findHessianElement(i3,j3)])
        tthess=np.array(tthess)
        return tthess

class amberAC(object):

    def __init__(self,father):
        self.__father=father
        self.atomtypelist=[None]
        self.atomchargelist=[None]
    @tryfunc
    def read(self):
        with open(self.__father.acname,'r') as f:
            string=f.readline()
            string=f.readline()
            for string in f.readlines():
                if string.find('BOND')>=0:
                    break
                ac=string.split()
                self.atomtypelist.append(ac[len(ac)-1])
                self.atomchargelist.append(float(ac[len(ac)-2]))
        return True


class gauCOM(object):
    g09rt='g09'
    g09a2rt='g09'
    def __init__(self,father):
        self.__father=father
        self.xyzfile=''
        self.atomlist=[None]
        self.atomtypelist=[None]
        self.atomchargelist=[None]
        self.coordslist=[]
        self.connectivity=''
        self.nozomudihdfunc=[]
        self.nozomuanglefunc=[]
        self.nozomubondfunc=[]
        self.additionfunc=[]
        self.nozomuvdw=[]
        self.xyz=''
        self.commandline=''


    # Parse
    @tryfunc
    def read(self):
        self.xyzfile=''
        self.atomlist=[None]
        self.atomtypelist=[None]
        self.atomchargelist=[None]
        self.coordslist=[]
        self.connectivity=''
        self.nozomudihdfunc=[]
        self.nozomuanglefunc=[]
        self.nozomubondfunc=[]
        self.additionfunc=[]
        self.nozomuvdw=[]
        self.xyz=''
        self.commandline=''
        with open(self.__father.comname,'r') as f:
            counter=0
            ifconnect=False
            ifamber=False
            line=''
            for line in f:
                if line=='\n':
                    counter+=1
                    continue
                if counter==0:
                    while True:
                        self.commandline+=line
                        line=next(f)
                        if line=='\n':
                            counter+=1
                            break
                    if self.commandline.find('connectivity')>=0:
                        ifconnect=True
                    if self.commandline.find('amber')>=0:
                        ifamber=True
                def molespecs(line):
                    self.xyzfile+=line
                    tmp=line.split()[0]
                    if tmp.find('-')>=0:
                        self.atomlist.append(tmp.split('-')[0])
                        if tmp.count('-')==2:
                            tmp=tmp.split('-')
                            self.atomtypelist.append(tmp[1])
                            self.atomchargelist.append(tmp[2])
                        elif tmp.count('-')==3:
                            tmp=tmp.split('-')
                            self.atomtypelist.append(tmp[1])
                            self.atomchargelist.append(-float(tmp[3]))
                    else:
                        self.atomlist.append(tmp)
                    self.coordslist.extend(line.split()[1:4])

                if counter==2:
                    self.__father.multiplicity=int(line.split()[1])
                    self.__father.totalcharge=int(line.split()[0])
                    while counter==2:
                        line=next(f)
                        if line=='\n':
                            counter+=1
                            break
                        molespecs(line)
                def connect(line):
                    self.connectivity+=line
                if counter==3:
                    if ifconnect:
                        line=next(f)
                        while counter==3:
                            if line=='\n':
                                counter+=1
                                break
                            connect(line)
                            line=next(f)

                if counter==4:
                    if ifamber:
                        line=next(f)
                        while counter==4:
                            if line=='\n':
                                counter+=1
                                break
                            thisline=mmfunction(line)
                            if thisline.type=='dihd':
                                self.nozomudihdfunc.append(thisline)
                            elif thisline.type=='angle':
                                self.nozomuanglefunc.append(thisline)
                            elif thisline.type=='bond':
                                self.nozomubondfunc.append(thisline)
                            elif thisline.type=='else':
                                self.additionfunc.append(thisline)
                            elif thisline.type=='vdw':
                                self.nozomuvdw.append(thisline)
                            line=next(f)




        self.coordslist=np.array(self.coordslist)
        for i in range(0,len(self.atomlist)-1):
            tmp=str(self.atomlist[i+1])+'   '+str(self.coordslist[3*i])+'   '+str(self.coordslist[3*i+1])+'   '+str(self.coordslist[3*i+2])+'\n'
            self.xyz+=tmp

    # File operation
    def rung09(self):
        ifchk=1 # if no chk, add.
        with open(self.__father.comname,'r') as f:
            for line in f.readlines():
                if line.find('%chk')>=0:  # !!!!!! leave later
                    line='%chk='+self.__father.chkname+'\n'
                    ifchk=0
            if ifchk==1:
                f.seek(0)
                content=f.read()
        if ifchk==1:
            with open(self.__father.comname,'w') as f:
                f.write('%chk='+self.__father.chkname+'\n')
                f.write(content)
        logging.info('Run g09 : '+gauCOM.g09rt+' '+self.__father.comname)
        os.system(gauCOM.g09rt+' '+self.__father.comname)
    def rung09a2(self):
        ifchk=1 # if no chk, add.
        with open(self.__father.comname,'r') as f:
            for line in f.readlines():
                if line.find('%chk')>=0:
                    line='%chk='+self.__father.chkname+'\n'
                    ifchk=0
            if ifchk==1:
                f.seek(0)
                content=f.read()
        if ifchk==1:
            with open(self.__father.comname,'w') as f:
                f.write('%chk='+self.__father.chkname+'\n')
                f.write(content)
        logging.info('Run g09a2 : '+gauCOM.g09a2rt+' '+self.__father.comname)
        os.system(gauCOM.g09a2rt+' '+self.__father.comname)


    def isover(self):
        logging.info('Checking g09 termination for'+self.__father.comname+'...')
        while True:
            output=''
            if not os.path.isfile(self.__father.logname):
                continue
            with open(self.__father.logname,'r') as f:
                for x in f.readlines()[:-6:-1]:
                    output+=x
            if output.find('Normal termination')>=0:
                logging.info('    ..normal termination')
                return True
            if output.find('Error termination')>=0:
                logging.critical('Error termination in '+self.__father.comname)
                raise rxccError('Error termination')
                return False
            time.sleep(2)
class gauLOG(object):
    antecommand='antechamber -c resp'
    def __init__(self,father):
        self.__father=father
        self.freq=0
    def getnatoms(self):
        with open(self.__father.logname,'r') as f:
            for x in f.readlines():
                if x.find('NAtoms')>=0:
                    self.natoms=int(x.split()[1])
                    return self.natoms
                    break

    def coordslast(self): #!!!!!! leave later
        with open(self.__father.logname,'r') as f:
            for x in list(reversed(f.readlines())):
                if x.find('orientation')>=0:
                    self.orn=x
                    break
        return self.orn
    def getfreq(self):
        freq=[]
        with open(self.__father.logname,'r') as f:
            for line in f.readlines():
                if line.find('Frequencies -- ')>=0:
                    freq.extend(line.split()[2:])
            if freq!=[]:
                freq=list(map(float,freq))
                # if (len(freq)>3*self.natoms-5):
                #     freq=freq[:len(freq)/2]
            else:
                logging.warning("No Freq found")
                return ['No Freq found']
        self.freq=freq
        return self.freq
    def runantecham(self):
        logging.info('Runing antechamber: \n')
        command=gauLOG.antecommand+' -i '+self.__father.logname+' -fi gout -o '+self.__father.acname+' -fo ac'
        logging.info(command)
        os.system(command)
class mmfunction(object):
    magicnum='XXXXXX'
    def __init__(self,line):
        fun=line.split()
        self.type=None
        self.repr=None

        def newfloat(value):
            if value=='XXXXXX':
                return mmfunction.magicnum
            else:
                return float(value)
        if fun[0]=='AmbTrs':
            self.type='dihd'
            self.a=fun[1]
            self.b=fun[2]
            self.c=fun[3]
            self.d=fun[4]
            self.npaths=float(fun[13])
            for i,paras in enumerate(fun[9:13]):
                if newfloat(paras)!=0.000:
                    i+=1
                    self.value=newfloat(paras)
                    break
            self.periodicity=i
            self.phase=int(fun[4+self.periodicity])
            self.repr=self.a+' '+self.b+' '+self.c+' '+self.d
        elif fun[0]=='HrmBnd1':
            self.type='angle'
            self.a=fun[1]
            self.b=fun[2]
            self.c=fun[3]
            self.value=newfloat(fun[4])
            self.eqvalue=newfloat(fun[5])
            self.repr=self.a+' '+self.b+' '+self.c
        elif fun[0]=='HrmStr1':
            self.type='bond'
            self.a=fun[1]
            self.b=fun[2]
            self.value=newfloat(fun[3])
            self.eqvalue=newfloat(fun[4])
            self.repr=self.a+' '+self.b
        # elif fun[0]=='ImpTrs':
        #     self.type='improper'
        #     self.a=fun[1]
        #     self.b=fun[2]
        #     self.c=fun[3]
        #     self.d=fun[4]
        #     self.value=newfloat(fun[5])
        #     self.eqvalue=newfloat(fun[6])
        #     self.repr=self.type+self.a+self.b+self.c+self.d
        elif fun[0]=='VDW':
            self.type='vdw'
            self.content=line
        else:
            self.type='else'
            self.content=line
