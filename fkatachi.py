#!/usr/bin/env python3
import argparse
import copy
import logging
import shutil
import os
import time
import itertools
from io import StringIO
import numpy as np
import rxcclib.utils as utils
import rxcclib.file.Gaussian as gau
import rxcclib.geometry.molecules as rxmol
from rxcclib.file.GauAmberCOM import MMFunction
from rxcclib.file.GauAmberCOM import GauAmberCOM


def readgeom(inputfile):
    # Read info from com

    mmcom = inputfile.com
    mole = rxmol.Molecule('thisgeometry')

    mmcom = GauAmberCOM(inputfile)

    mmcom.read()
    mole.readfromxyz(mmcom.xyz)
    mole.readchargefromlist(mmcom.atomchargelist)
    mole.readtypefromlist(mmcom.atomtypelist)
    mole.readconnectivity(mmcom.connectivity)
    for atom in mole:
        atom.vdwradius = float(mmcom.vdwdict[atom.atomtype][0])
        atom.vdwwelldepth = float(mmcom.vdwdict[atom.atomtype][1])

    # Store and count finalfunc
    finalfuncL = []
    finalfuncL.extend(sorted(mmcom.bondfunc, key=lambda x: repr(x)))
    finalfuncL.extend(sorted(mmcom.anglefunc, key=lambda x: repr(x)))
    finalfuncL.extend(sorted(mmcom.dihdfunc, key=lambda x: repr(x)))
    finalfuncL.extend(sorted(mmcom.improperfunc, key=lambda x: repr(x)))

    for item in mmcom.improperfunc:
        for atom3 in mole:
            if atom3.atomtype == item.c:
                permu = itertools.permutations(atom3.neighbor, 3)
                res = []
                for tu in permu:
                    a = tu[0].atomtype == item.a or item.a == '*'
                    b = tu[1].atomtype == item.b or item.b == '*'
                    c = tu[2].atomtype == item.d or item.d == '*'
                    if a and b and c:
                        res.append([
                            tu[0].atomnum, tu[1].atomnum, atom3.atomnum,
                            tu[2].atomnum
                        ])
                res = sorted(res, key=lambda x: (str(x[1]) + str(x[3])))
                res = res[0]
                mole.addimproper(*res)

    # Match itnl and finalfunc
    itnlcordL = []
    itnlcordL.extend(sorted(mole.dihedrallist.values(), key=lambda x: repr(x)))
    itnlcordL.extend(sorted(mole.anglelist.values(), key=lambda x: repr(x)))
    itnlcordL.extend(sorted(mole.bondlist.values(), key=lambda x: repr(x)))
    itnlcordL.extend(sorted(mole.improperlist.values(), key=lambda x: repr(x)))

    unkitnlL = []
    knownitnlL = []
    for item in itnlcordL:
        for func in finalfuncL:
            if matchitnlwithfinalfunc(item, func):
                item.func = func
                if type(item) is rxmol.Dihedral:
                    item.dihdfunctions = copy.deepcopy(func.dihdfunctions)
                    item.npaths = func.npaths
                elif type(item) is rxmol.Improper:
                    item.forceconst = func.forceconst
                    item.phase = func.phase
                    item.periodicity = func.periodicity
                else:
                    item.forceconst = func.forceconst
                    item.eqvalue = func.eqvalue
                break

    return finalfuncL, itnlcordL, mole


def matchitnlwithfinalfunc(item, finalfunc):
    if type(item) is rxmol.Dihedral and finalfunc.type == 'dihd':
        a = (item[1].atomtype == finalfunc.a or finalfunc.a == '*')
        b = (item[2].atomtype == finalfunc.b or finalfunc.b == '*')
        c = (item[3].atomtype == finalfunc.c or finalfunc.c == '*')
        d = (item[4].atomtype == finalfunc.d or finalfunc.d == '*')
        forward = a and b and c and d
        a = (item[1].atomtype == finalfunc.d or finalfunc.d == '*')
        b = (item[2].atomtype == finalfunc.c or finalfunc.c == '*')
        c = (item[3].atomtype == finalfunc.b or finalfunc.b == '*')
        d = (item[4].atomtype == finalfunc.a or finalfunc.a == '*')
        backward = a and b and c and d
        if forward or backward:
            return True
        else:
            return False
    elif type(item) is rxmol.Angle and finalfunc.type == 'angle':
        a = (item[1].atomtype == finalfunc.a or finalfunc.a == '*')
        b = (item[2].atomtype == finalfunc.b or finalfunc.b == '*')
        c = (item[3].atomtype == finalfunc.c or finalfunc.c == '*')
        forward = a and b and c
        a = (item[1].atomtype == finalfunc.c or finalfunc.c == '*')
        b = (item[2].atomtype == finalfunc.b or finalfunc.b == '*')
        c = (item[3].atomtype == finalfunc.a or finalfunc.a == '*')
        backward = a and b and c
        if forward or backward:
            return True
        else:
            return False
    elif type(item) is rxmol.Bond and finalfunc.type == 'bond':
        a = (item[1].atomtype == finalfunc.a or finalfunc.a == '*')
        b = (item[2].atomtype == finalfunc.b or finalfunc.b == '*')
        forward = a and b
        a = (item[1].atomtype == finalfunc.b or finalfunc.b == '*')
        b = (item[2].atomtype == finalfunc.a or finalfunc.a == '*')
        backward = a and b
        if forward or backward:
            return True
        else:
            return False
    elif type(item) is rxmol.Improper and finalfunc.type == 'improper':
        a = (item[1].atomtype == finalfunc.a or finalfunc.a == '*')
        b = (item[2].atomtype == finalfunc.b or finalfunc.b == '*')
        c = (item[3].atomtype == finalfunc.c or finalfunc.c == '*')
        d = (item[4].atomtype == finalfunc.d or finalfunc.d == '*')
        forward = a and b and c and d
        if forward:
            return True
        else:
            return False
    else:
        return False


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


def addlink1(mmfile, itnlcordL):
    #    return ''
    content = '\n--link1--\n'
    content += ('%chk=' + mmfile.chkname + '\n')
    content += ('#p geom=allcheck ')
    content += ('freq=(readfc,modredundant,intmodes) '
                'iop(4/33=3,7/33=1,99/5=5)\n\n')
    content += ('* * K\n* * * K\n* * * * K\n')
    for item in itnlcordL:
        if type(item) == rxmol.Bond:
            content += (
                str(item[1].atomnum) + ' ' + str(item[2].atomnum) + ' A\n')
        if type(item) == rxmol.Angle:
            content += (str(item[1].atomnum) + ' ' + str(item[2].atomnum) + ' '
                        + str(item[3].atomnum) + ' A\n')
        if type(item) == rxmol.Dihedral:
            content += (
                str(item[1].atomnum) + ' ' + str(item[2].atomnum) + ' ' +
                str(item[3].atomnum) + ' ' + str(item[4].atomnum) + ' A\n')
        if type(item) == rxmol.Improper:
            content += (
                str(item[1].atomnum) + ' ' + str(item[2].atomnum) + ' ' +
                str(item[3].atomnum) + ' ' + str(item[4].atomnum) + ' A\n')
    content += '\n'
    return content


def summarize(finalfuncL, itnlcordL, originalname, finalhead, method, mmcom):
    # Summarize
    logging.info('Start Summarizing')
    for func in finalfuncL:
        if func.type == 'dihd':
            pass
        elif func.type == 'improper':
                res = 0
                i = 0
                for item in itnlcordL:
                    if item.func == func:
                        res += item.forceconst
                        i += 1
                func.forceconst = res / i
        else:
            res = 0
            eqres = 0
            i = 0
            for item in itnlcordL:
                if item.func == func:
                    res += item.forceconst
                    eqres += item.eqvalue
                    i += 1
            func.forceconst = res / i
            func.eqvalue = eqres / i


    # Build tailstring
    tailstring = ''
    for dihd in mmcom.dihdfunc:
        parm = []
        phase = []
        for item in dihd.dihdfunctions:
            parm.append(item.forceconst)
            phase.append(item.phase)

        tailstring += 'AmbTrs  ' + ' '.join(
            [x.center(3, ' ') for x in dihd.repr.split()]) + '  ' + ' '.join(
                [str(x).center(3, ' ') for x in phase]) + '  ' + ' '.join(
                    ['{:>.10f}'.format(x)
                     for x in parm]) + '   ' + str(dihd.npaths) + '\n'
    for angle in mmcom.anglefunc:
        if angle.forceconst == MMFunction.unknownsign:
            parm = '0.000'
            logging.critical('Force constant is not determined for angle ' +
                             angle.repr)
            raise
        else:
            parm = angle.forceconst
        tailstring += 'HrmBnd1  ' + ' '.join(
            [x.center(3, ' ')
             for x in angle.repr.split()]) + '  ' + '{:>.10f}'.format(
                 parm) + '  ' + '{:>9.5f}'.format(angle.eqvalue) + '\n'
    for bond in mmcom.bondfunc:
        if bond.forceconst == MMFunction.unknownsign:
            parm = '0.000'
            logging.critical('Force constant is not determined for bond ' +
                             bond.repr)
            raise
        else:
            parm = bond.forceconst
        tailstring += 'HrmStr1  ' + ' '.join(
            [x.center(3, ' ')
             for x in bond.repr.split()]) + '  ' + '{:>.10f}'.format(
                 parm) + '  ' + '{:>7.5f}'.format(bond.eqvalue) + '\n'
    for improper in mmcom.improperfunc:
        if improper.forceconst == MMFunction.unknownsign:
            logging.critical('Force constant is not determined for improper ' +
                             improper.repr)
            raise
        else:
            parm = improper.forceconst
            if parm < 0:
                improper.phase = 180
                parm = -parm
        tailstring += 'ImpTrs  ' + ' '.join([
            x.center(3, ' ') for x in improper.repr.split()
        ]) + '  ' + '{:>.10f}'.format(parm) + '  ' + '{:6.2f}'.format(
            improper.phase) + '  2.0\n'
    for x in mmcom.additionfunc:
        tailstring += x.content
    for vdw in mmcom.vdw:
        tailstring += 'VDW  ' + '  ' + vdw.atomtype + \
            '  ' + vdw.radius + '  ' + vdw.welldepth + '\n'
    tailstring += '\n\n'

    finalname = method + originalname

    with open(finalname, 'w') as f:
        f.write(finalhead + tailstring)

    return finalname



def fkatachi(inputfile, loopid, convthreshold):
    global finalfuncL
    global mole
    global itnlcordL
    finalfuncL, itnlcordL, mole = readgeom(inputfile)
    shutil.copy(inputfile.comname, 'MM0.com')

    conv = convthreshold
    prebondmax = 100
    preanglemax = 100
    while True:
        logging.info('START LOOP' + str(loopid))
        currentfile = gau.GauFile('MM'+str(loopid))
        GauAmberCOM(currentfile)
        currentfile.com.read()

        with open(currentfile.comname,'r') as f:
            content = f.read()
        with open(currentfile.comname, 'w') as f:
            f.write("%chk="+ currentfile.chkname + '\n' +
                content + addlink1(currentfile, itnlcordL))

        head = content.split('AmbTrs')[0]
        if head == content:
            head = content.split('HrmBnd1')[0]
            if head == content:
                head = content.split('HrmStr1')[0]

        currentfile.com.Popen()
        currentfile.com.wait()
        currentfile.chk.formchk()
        currentfile.fchk.read()
        forces = currentfile.fchk.intforces
        findmatch(itnlcordL, currentfile)

        bondmax = 0
        anglemax = 0
        for item in itnlcordL:
            if type(item) is rxmol.Bond:
                # convert forces to Kcal mol-1 Ang-1
                if abs(forces[item.gauseq]) > bondmax:
                    bondmax = abs(forces[item.gauseq])
                delta = 0.5 * forces[item.gauseq] * 1185.8211 / item.forceconst
            elif type(item) is rxmol.Angle:
                # convert forces to Kcal mol-1 rad-1
                if abs(forces[item.gauseq]) > anglemax:
                    anglemax = abs(forces[item.gauseq])
                delta = 180* 0.5 * forces[item.gauseq] * 627.5095 / item.forceconst / np.pi
            else:
                continue
            item.eqvalue -= delta

        logging.info('Conv Threshold: '+str(conv))
        if anglemax == 0:
            anglemax = 100
        logging.info('Bondmax:'+ str(bondmax))
        logging.info('Anglemax:'+ str(anglemax))
        if bondmax < prebondmax:
            conv = convthreshold
            prebondmax = bondmax
        if anglemax < preanglemax:
            conv = convthreshold
            preanglemax = anglemax
        loopid+=1
        conv-=1
        if abs(bondmax) < 5.0e-5 and abs(anglemax) < 5.0e-5:
            logging.info('Exit for small max forces:'+str(bondmax)+','+str(anglemax))
            with open(os.path.join('..', 'fkatachi_'+inputfile.comname), 'w') as f:
                f.write(content)
            return loopid

        if conv <= 0:
            logging.info('Exit for threshold')
            with open(os.path.join('..', 'fkatachi_'+inputfile.comname), 'w') as f:
                f.write(content)
            return loopid

        summarize(finalfuncL, itnlcordL, 'MM'+str(loopid)+'.com', head, '', inputfile.com)






if __name__ == '__main__':
    start = time.perf_counter()
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfile', help='name of com file')
    parser.add_argument('loopid', default=0)
    parser.add_argument('convthreshold', default=10)
    args = parser.parse_args()

    loopid = int(args.loopid)
    convthreshold = int(args.convthreshold)
    inputfile = args.inputfile
    inputfile = inputfile.split('.com')[0]
    print(inputfile)
    try:
        shutil.rmtree('fkatachi_'+inputfile.split('_result')[0])
    except:
        pass

    os.mkdir('fkatachi_'+inputfile.split('_result')[0])
    os.chdir('fkatachi_'+inputfile.split('_result')[0])
    shutil.copy(os.path.join('..', args.inputfile), '.')

    inputfile = gau.GauFile(inputfile)
    inputcom = GauAmberCOM(inputfile)


    logging.basicConfig(
        filename=inputfile.basename + '.katachi',
        level=logging.DEBUG,
        filemode='w')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(levelname)-8s %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    id = fkatachi(inputfile, loopid, convthreshold) -1
    id = str(id)
    stop = time.perf_counter()
    t = stop - start
    os.system('echo "'+str(t)+' '+id+ '" >> ~/tfkatachi')
