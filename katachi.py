#!/usr/bin/env python3
from __future__ import print_function
import os,argparse,logging,pdb,io,copy
import numpy as np
import time
import shutil
import molecules as rxmol
import chemfiles as rxccfile


def matchdihd(dihd,func):
    a=(dihd[1].atomtype==func.a or func.a=='*')
    b=(dihd[2].atomtype==func.b or func.b=='*')
    c=(dihd[3].atomtype==func.c or func.c=='*')
    d=(dihd[4].atomtype==func.d or func.d=='*')
    forward=a and b and c and d
    a=(dihd[1].atomtype==func.d or func.d=='*')
    b=(dihd[2].atomtype==func.c or func.c=='*')
    c=(dihd[3].atomtype==func.b or func.b=='*')
    d=(dihd[4].atomtype==func.a or func.a=='*')
    backward=a and b and c and d
    if forward or backward:
        return True
    else:
        return False
def matchbond(bond,func):
    a=(bond[1].atomtype==func.a or func.a=='*')
    b=(bond[2].atomtype==func.b or func.b=='*')
    forward=a and b
    a=(bond[1].atomtype==func.b or func.b=='*')
    b=(bond[2].atomtype==func.a or func.a=='*')
    backward=a and b
    if forward or backward:
        return True
    else:
        return False
def matchangle(angle,func):
    a=(angle[1].atomtype==func.a or func.a=='*')
    b=(angle[2].atomtype==func.b or func.b=='*')
    c=(angle[3].atomtype==func.c or func.c=='*')
    forward=a and b and c
    a=(angle[1].atomtype==func.c or func.c=='*')
    b=(angle[2].atomtype==func.b or func.b=='*')
    c=(angle[3].atomtype==func.a or func.a=='*')
    backward=a and b and c
    if forward or backward:
        return True
    else:
        return False


def katachi(mmresult,loopid,opt,convthreshold): # mmresult: phf_result_mmxxx   without extension
    # Def read initial:
    foldername = 'katachi_'+mmresult.split('_result')[0]
    try:
        shutil.rmtree(foldername)
    except:
        pass
    os.mkdir(foldername)
    os.chdir(foldername)
    shutil.copy(os.path.join('..',mmresult+'.com'), '.')

    init=rxccfile.File(mmresult)
    stdgeom=rxmol.Molecule('stdgeom')
    init.com.read()
    xyz=io.StringIO(init.com.xyz)
    stdgeom.readfromxyz(xyz)
    stdgeom.readchargefromlist(init.com.atomchargelist)
    stdgeom.readtypefromlist(init.com.atomtypelist)
    stdgeom.readconnectivity(init.com.connectivity)
    L=[]
    L.extend(stdgeom.dihdlist.values())
    L.extend(stdgeom.anglelist.values())
    L.extend(stdgeom.bondlist.values())
    nozomuL=[]
    nozomuL.extend(init.com.nozomudihdfunc)
    nozomuL.extend(init.com.nozomuanglefunc)
    nozomuL.extend(init.com.nozomubondfunc)

    stdL=[]
    stdL.extend(init.com.nozomuanglefunc)
    stdL.extend(init.com.nozomubondfunc)


    # Def iteration (stdgeom,currentfile):
    def iteration(stdgeom,currentfile,opt,loopid,convthreshold):
        os.system('cp '+currentfile.comname+' MM0.com')
        with open('MM0.com','r') as f:
            dihds=''
            for line in f:
                if line.find('AmbTrs')>=0:
                    dihds+=line
        while True:
            currentfile=rxccfile.File('MM'+str(loopid))
            if opt=='opt':
                os.system('sed -i "s/#p opt=(nomicro,cartesian) /#p /g" '+currentfile.comname)
                os.system('sed -i "s/#p/#p opt=(nomicro,cartesian)/g" '+currentfile.comname)
                os.system('sed -i "/freq/d" '+currentfile.comname)
                os.system('sed -i "/chk/d" '+currentfile.comname)
            elif opt=='calcall':
                os.system('sed -i "s/#p opt=(nomicro,cartesian,tight,calcall) /#p /g" '+currentfile.comname)
                os.system('sed -i "s/#p opt=(nomicro,cartesian,tight,calcall) /#p /g" '+currentfile.comname)
                os.system('sed -i "s/#p/#p opt=(nomicro,cartesian,tight,calcall)/g" '+currentfile.comname)
                os.system('sed -i "/freq/d" '+currentfile.comname)
                os.system('sed -i "/chk/d" '+currentfile.comname)

#            if loopid>1000:
 #               raise StopIteration
            currentfile.com.read()
            try:
                currentfile.com.rung09()
                currentfile.com.isover()
                currentfile.runformchk()
                ifstop=True
                os.system('rm '+currentfile.chkname+' '+currentfile.logname)
            except:
                logging.error("Calculation failed, try again.")
                try:
                    currentfile.com.rung09()
                    currentfile.com.isover()
                    currentfile.runformchk()
                    ifstop=True
                    os.system('rm '+currentfile.chkname+' '+currentfile.logname)
                except:
                    logging.critical('Calculation still failed, continue...')
                    currentfile.runformchk()
                    os.system('rm '+currentfile.chkname+' '+currentfile.logname)
                    ifstop=False
                    # logging.info('minimum max2 is chosen from loop'+str(minmax2loop))
                    # os.system('cp MM'+str(minmax2loop)+'.com ../'+mmresult[0:3]+'amd'+mmresult[3:]+'.com')
                    # os.system('sed -i "s/#p opt=(verytight,z-matrix,calcall)/#p /g" ../*.com')
                    # os.system('sed -i "/chk/d" ../*.com')
                    # return

            currentfile.fchk.read()


            currentgeom=rxmol.Molecule('currentgeom')
            currentgeom.readfromxyz(io.StringIO(currentfile.fchk.xyz))
            currentgeom.readchargefromlist(currentfile.com.atomchargelist)
            currentgeom.readtypefromlist(currentfile.com.atomtypelist)
            currentgeom.readconnectivity(currentfile.com.connectivity)

            for angle in currentgeom.anglelist.values():
                    for anglefunc in currentfile.com.nozomuanglefunc:
                        if matchangle(angle,anglefunc):
                            angle.nozomufunc=anglefunc
            for bond in currentgeom.bondlist.values():
                    for bondfunc in currentfile.com.nozomubondfunc:
                        if matchbond(bond,bondfunc):
                            bond.nozomufunc=bondfunc

            # reassign current eqvalue
            for nozomufunc in currentfile.com.nozomuanglefunc:
                eq=0
                i=0
                for angle in currentgeom.anglelist.values():
                    if angle.nozomufunc==nozomufunc:
                        eq+=angle.anglevalue
                        i+=1
                nozomufunc.eqvalue=eq/i
            for nozomufunc in currentfile.com.nozomubondfunc:
                eq=0
                i=0
                for bond in currentgeom.bondlist.values():
                    if bond.nozomufunc==nozomufunc:
                        eq+=bond.length
                        i+=1
                nozomufunc.eqvalue=eq/i

            currentL=[]
            currentL.extend(currentfile.com.nozomuanglefunc)
            currentL.extend(currentfile.com.nozomubondfunc)
            for item1 in currentL:
                for item2 in currentL:
                    if item1.value==item2.value and item1.repr!=item2.repr:
                        item1.eqvalue=(item1.eqvalue+item2.eqvalue)/2
                        item2.eqvalue=item1.eqvalue
                        logging.debug('Averaged old eqvalue '+item1.repr+' and '+item2.repr+' '+str(item2.eqvalue))

            delta1=[x.eqvalue-y.eqvalue for x,y in zip(stdL,currentL) if x.type=='bond']
            delta2=[x.eqvalue-y.eqvalue for x,y in zip(stdL,currentL) if x.type=='angle']
            max1=sorted(delta1,key=abs,reverse=True)[0]
            max2=sorted(delta2,key=abs,reverse=True)[0]
            try:
                if abs(max2)<abs(minmax2):
                    minmax2=abs(max2)
                    minmax2loop=loopid
            except:
                minmax2=abs(max2)
                minmax2loop=loopid

            if loopid-minmax2loop>convthreshold:
                logging.info('Stopped for convergence: max Delta2 do not decrease in '+str(convthreshold)+' cycles')
                if opt=='calcall':
                    logging.info('minimum max2 is chosen from loop'+str(minmax2loop))
                    os.system('cp MM'+str(minmax2loop)+'.com ../katachi_'+mmresult+'.com')
                    os.system('sed -i "s/#p opt=(nomicro,cartesian,tight,calcall)/#p freq/g" ../*.com')
                    os.system('sed -i "/chk/d" ../*.com')

                    return loopid
                else:
                    opt='calcall'
                    minmax2=100
                    minmax2loop=loopid

            logging.info('------------------------------')
            logging.info('Loop '+str(loopid)+' keyword '+opt+': max bond delta: '+str(max1)+'  max angle delta: '+str(max2))
            logging.info('MinMax2 is '+str(minmax2)+' at loop '+str(minmax2loop))
            logging.info('------------------------------')
            if abs(max1)<0.0001 and abs(max2)<0.01:

                if opt=='opt':
                    opt='calcall'
                    logging.info('-------------------------------------------')
                    logging.info('opt converged at '+str(max1)+' '+str(max2)+' '+str(loopid))
                    logging.info('-------------------------------------------')
                    minmax2=100
                    minmax2loop=loopid+1
                elif opt=='calcall':
                    logging.info('-------------------------------------------')
                    logging.info('calcall converged at '+str(max1)+' '+str(max2)+' '+str(loopid))
                    logging.info('-------------------------------------------')
                    os.system('cp '+currentfile.comname+' ../katachi_'+mmresult+'.com')
                    os.system('sed -i "s/#p opt=(nomicro,cartesian,tight,calcall)/#p freq/g" ../*.com')
                    os.system('sed -i "/chk/d" ../*.com')
                    if ifstop:
                        return loopid

            delta=[x.eqvalue-y.eqvalue for x,y in zip(stdL,currentL)]

            try:
                type(last)
            except UnboundLocalError:
                if loopid!=0:
                    lastfile=rxccfile.File('MM'+str(loopid-1))
                    lastfile.com.read()
                    lastL=[]
                    lastL.extend(lastfile.com.nozomuanglefunc)
                    lastL.extend(lastfile.com.nozomubondfunc)
                    last=lastL
                else:
                    last=copy.deepcopy(stdL)
            for now,std,former,delt in zip(currentL,stdL,last,delta):
                if former.eqvalue+delt>0 and former.eqvalue+delt<180:
                    now.eqvalue=former.eqvalue+delt
                if now.type=='bond' and abs(delt)>0.0001:
                    logging.warning(now.repr+' '+str(delt))
                if now.type=='angle' and abs(delt)>0.01:
                    logging.warning(now.repr+' '+str(delt))
            for item1 in currentL:
                for item2 in currentL:
                    if item1.value==item2.value and item1.repr!=item2.repr:
                        item1.eqvalue=(item1.eqvalue+item2.eqvalue)/2
                        item2.eqvalue=item1.eqvalue
                        logging.debug('Averaged new eqvalue '+item1.repr+' and '+item2.repr+' '+str(item2.eqvalue))

            last=copy.deepcopy(currentL)


            finalxyz=''
            for atom in stdgeom:
                finalxyz+=atom.atomsym+'-'+atom.atomtype+'-'+'{:<9.6f}'.format(float(atom.atomcharge))+'   '+'    '.join(["{: .12f}".format(x) for x in atom.coords])+'\n'
            finalhead=currentfile.com.commandline+'\nfinal\n\n'+str(currentfile.fchk.totalcharge)+' '+str(currentfile.fchk.multiplicity)+'\n'+finalxyz+'\n'+currentfile.com.connectivity+'\n'

            finaltail=''
            finaltail+=dihds


            for item in currentfile.com.nozomuanglefunc:
                finaltail+='HrmBnd1 '+item.repr+' '
                parm="{: .3f}".format(item.value)
                finaltail+=' '+parm+' {: .4f}'.format(item.eqvalue)+'\n'
            for item in currentfile.com.nozomubondfunc:
                finaltail+='HrmStr1 '+item.repr+' '
                parm="{: .3f}".format(item.value)
                finaltail+=' '+parm+' {: .5f}'.format(item.eqvalue)+'\n'

            for addfunc in currentfile.com.additionfunc:
                finaltail+=addfunc.content
            for nozovdw in currentfile.com.nozomuvdw:
                finaltail+=nozovdw.content

            finaltail+='\n\n'

            loopid+=1
            with open('MM'+str(loopid)+'.com','w') as f:
                f.write(finalhead+finaltail)


            del currentgeom
            del currentfile


    id = iteration(stdgeom,init,opt,loopid,convthreshold)
    return id




if __name__=='__main__':
    # Parse Input
    start = time.perf_counter()
    parser=argparse.ArgumentParser()
    parser.add_argument('mmresult',help="parameterized result MM file.")
    parser.add_argument('loopid',help="loopid",default=0)
    parser.add_argument('opt',help='opt or calcall',default='opt')
    parser.add_argument('convthreshold',help='convergence threshold',default=10)
    args=parser.parse_args()
    mmresult=args.mmresult
    mmresult=mmresult[:mmresult.find('.')]
    loopid=int(args.loopid)
    convthreshold=int(args.convthreshold)
    opt=args.opt
    if loopid==0:
        logging.basicConfig(filename=args.mmresult+'.katachiout',level=logging.DEBUG,filemode='w')
    else:
        logging.basicConfig(filename=args.mmresult+'.katachiout',level=logging.DEBUG)
    console=logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter=logging.Formatter('%(levelname)-8s %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)


    id = katachi(mmresult,loopid,opt,convthreshold) 
    id = str(id)
    stop = time.perf_counter()
    t = stop - start

    os.system('echo "'+str(t)+' '+id+ '" >> ~/tkatachi')
