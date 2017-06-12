# From fchk read Charge, Multiplicity, Coordinates
import os
import subprocess
import time
import logging
import shutil
from io import StringIO
import numpy as np
import rxcclib.utils as utils

table = utils.PeriodicTable()


class rxFileError(Exception):
    def __init(self, value):
        self.value = value

    def __repr__(self):
        return repr(self.value)

    __str__ = __repr__


class GauFile(object):
    """
    General control of Gaussian Files.
    """

    def __init__(self, filename):
        pwd = os.path.split(os.path.abspath(filename))
        self.pwd = pwd[0]
        self.basename = pwd[1]
        self._com = GauCOM(self)
        self._log = GauLOG(self)
        self._chk = GauCHK(self)
        self._fchk = GauFCHK(self)
        self.comname = self.basename + '.com'
        self.chkname = self.basename + '.chk'
        self.fchkname = self.basename + '.fchk'
        self.logname = self.basename + '.log'

    @property
    def com(self):
        return self._com

    @property
    def log(self):
        return self._log

    @property
    def fchk(self):
        return self._fchk

    @property
    def chk(self):
        return self._chk




class GauFCHK(object):
    def __init__(self, parent):
        self._parent = parent

    @property
    def filename(self):
        return self._parent.fchkname

    @property
    def abspath(self):
        return os.path.join(self._parent.pwd, self.filename)

    def read(self):
        logging.info('Read fchk:' + self.abspath)
        # FCHK parser
        with open(self.abspath, 'r') as f:
            for line in f:

                if line.find('Charge') == 0:
                    self.charge = int(line.split('I')[1])
                if line.find('Mult') == 0:
                    self.mult = int(line.split('I')[1])
                if line.find('Atomic numbers') == 0:
                    self.natom = int(line.split('=')[1])
                    self.atomnos = []
                    for line in f:
                        try:
                            self.atomnos.extend([int(x) for x in line.split()])
                        except ValueError:
                            assert len(self.atomnos) == self.natom, (
                                "Error: len(atomnos)" + str(len(atomnos)) +
                                " != natom ! " + str(self.natom))
                            break

                if line.find('Current cartesian coordinates') == 0:
                    coordslist = []
                    conv = lambda x: utils.convertor(x, "bohr", "Angstrom")
                    for line in f:
                        try:
                            coordslist.extend(
                                [conv(float(x)) for x in line.split()])
                        except ValueError:
                            assert (len(coordslist) == self.natom * 3)
                            break
                    self.atomcoords = []
                    self.atomcoords.append(utils.atomwiseList.toGroupOfThree(coordslist))

                # Read Gradient(Forces)
                if line.find('Cartesian Gradient') == 0:
                    gradslist = []
                    for line in f:
                        try:
                            tmp = line.split()
                            for i, item in enumerate(tmp):
                                # 1.1234567E-123 will be 1.1234567-123, so fix it
                                if item[-4] == '-':
                                    tmp[i] = item[:-4] + 'E' + item[-4:]
                            gradslist.extend([float(x) for x in tmp])
                        except (ValueError, IndexError):
                            assert (len(gradslist) == 3 * self.natom)
                            break
                    self.grads = utils.atomwiseList.toGroupOfThree(gradslist)

                # Read Hessian
                if line.find('Cartesian Force Constants') == 0:
                    self.hessian = []
                    for line in f:
                        try:
                            tmp = line.split()
                            tmpline = []
                            for item in tmp:
                                if item[-4] == '-':
                                    item = item[:-4] + 'E' + item[-4:]
                                tmpline.append(float(item))
                            self.hessian.extend([float(x) for x in tmpline])
                        except ValueError:
                            assert (len(self.hessian) == 4.5 * self.natom**2 +
                                    1.5 * self.natom)
                            break

                # Read Internal Forces
                if line.find('Internal Forces') == 0:
                    self.intforces = []
                    for line in f:
                        try:
                            self.intforces.extend(
                                [float(x) for x in line.split()])
                        except ValueError:
                            break

                # Read Internal Hessian
                if line.find('Internal Force Constants') == 0:

                    self.inthessian = []
                    for line in f:
                        try:
                            self.inthessian.extend(
                                [float(x) for x in line.split()])
                        except ValueError:
                            break

        return True

    @property
    def xyz(self):
        tmp = list(
            zip(
                map(lambda x: table.element[x], self.atomnos), self.atomcoords[
                    -1]))
        string = ''
        for item in tmp:
            string += '{}\t{}\n'.format(
                item[0], '  '.join(map('{: 9.8f}'.format, item[1])))
        return string

    def findHessianElement(self, i, j):  # i, j: coordinate number
        return self.hessian[utils.getNum(i, j)]

    def find33Hessian(self, i, j):  # i, j: atom number
        if i < j:
            i, j = j, i
        tthess = []
        i1 = 3 * i - 3
        i2 = 3 * i - 2
        i3 = 3 * i - 1
        j1 = 3 * j - 3
        j2 = 3 * j - 2

        j3 = 3 * j - 1
        tthess.append([
            self.findHessianElement(i1, j1), self.findHessianElement(i1, j2),
            self.findHessianElement(i1, j3)
        ])
        tthess.append([
            self.findHessianElement(i2, j1), self.findHessianElement(i2, j2),
            self.findHessianElement(i2, j3)
        ])
        tthess.append([
            self.findHessianElement(i3, j1), self.findHessianElement(i3, j2),
            self.findHessianElement(i3, j3)
        ])
        tthess = np.array(tthess)
        return tthess

    def findintHessianElement(self, i, j):  # i, j: coordinate number
        return self.inthessian[utils.getNum(i, j)]


class GauCOM(object):
    gausub = 'g09'

    def __init__(self, parent):
        self._parent = parent

    @property
    def filename(self):
        return self._parent.comname

    @property
    def abspath(self):
        return os.path.join(self._parent.pwd, self.filename)

    # normal xyz COM file:
    def blocksplit(self):
        with open(self.abspath, 'r') as f:
            content = f.read()
            tmp = content.split('\n')
            block = ''
            for item in tmp:
                if item.isspace():
                    item = ''
                block += item + '\n'
            block = block.split('\n\n')
            block = [x + '\n' for x in block]

        return block

    def Popen(self, addchk=True):
        if addchk is True:
            blocks = self.blocksplit()
            if blocks[0].find('%chk') < 0:
                with open(self.abspath, 'w') as f:
                    f.write('%chk={}\n'.format(self._parent.chk.abspath))
                    for item in blocks:
                        f.write(item + '\n')

        logging.info('Run g09 : ' + GauCOM.gausub + ' ' + self.abspath + ' ' + self._parent.log.abspath)
        self.run = subprocess.Popen(
            [GauCOM.gausub, self.abspath, self._parent.log.abspath])

        return True

    def poll(self):
        """ Return False if not ended, True if normal term, raise rxFileError if error.
        """

        if not os.path.isfile(self._parent.log.abspath):
            return False

        f = utils.FileReadBackwards(self._parent.log.abspath)
        fp = iter(f)
        lastlines = ''.join([next(fp) for i in range(5)])


        if lastlines.find('Normal termination') >= 0:
            logging.info('    ..normal termination')
            return True
        elif lastlines.find('Error termination') >= 0:
            logging.error('Error termination in ' + self._parent.comname)
            raise rxFileError('G09 Error termination')
        else:
            return False

    def wait(self):
        logging.info('Wait Gaussian termination for ' + self._parent.comname +
                     '...')

        self.run.wait()
        while True:
            if self.poll() is True:
                return True


class GauLOG(object):
    def __init__(self, parent):
        self._parent = parent

    @property
    def filename(self):
        return self._parent.logname

    @property
    def abspath(self):
        return os.path.join(self._parent.pwd, self.filename)


class GauCHK(object):
    def __init__(self, parent):
        self._parent = parent

    @property
    def filename(self):
        return self._parent.chkname

    @property
    def abspath(self):
        return os.path.join(self._parent.pwd, self.filename)

    def formchk(self):
        args = 'formchk -3 {} {}'.format(self.abspath, self._parent.fchk.abspath)
        args = args.split()
        logging.info('  ' + ' '.join(args))
        iferror = subprocess.run(
            args,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT)
        if iferror.stdout.decode().find('Error') >= 0:
            raise rxFileError('   Error in formatting' + self.filename)
        return iferror

