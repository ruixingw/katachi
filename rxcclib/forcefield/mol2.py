import os
import rxcclib.io.mol2 as rxmol2
class Mol2(rxmol2.Mol2):
    def __init__(self, parent):
        self._parent = parent
        parent.mol2 = self
        parent.mol2name = parent.basename + '.mol2'

    @property
    def filename(self):
        return self._parent.mol2name

    @property
    def abspath(self):
        return os.path.join(self._parent.pwd, self.filename)
