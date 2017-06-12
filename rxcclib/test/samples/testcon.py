import rxcclib.utils as utils
import rxcclib.lazy as lazy
import rxcclib.file.Gaussian as rxfile
import rxcclib.geometry.molecules as rxmol

with open('cnnty.com','r') as f:
    con=f.read()


con=lazy.ConnectivityToConnectionMatrix(con)

con=utils.LTMatrix.newFromUpperMat(con)

benmol=rxmol.Molecule('ben')
benfile=rxfile.GauFile('bencom')
benfile.fchk.read()

benmol.addatomsFromLists(benfile.fchk.atomnos,benfile.fchk.atomcoords[-1])


benmol.readConnectionMatrix(con.fullmat)

