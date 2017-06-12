import os

class Mol2(object):
    def __init__(self, filename):
        self._filename = filename

    @property
    def filename(self):
        return self._filename

    @property
    def abspath(self):
        return os.path.abspath(self.filename)

    def read(self):
        self.atomtypelist = [None]
        self.atomchargelist = [None]
        with open(self.abspath, 'r') as f:
            for line in f:
                if line.find('ATOM') >= 0 and line.find('@') >= 0:
                    break
            for line in f:
                if line.find('@') >= 0:
                    break
                mol2 = line.split()
                self.atomtypelist.append(mol2[-4])
                self.atomchargelist.append(float(mol2[-1]))
        assert len(self.atomtypelist) == len(self.atomchargelist)
        return True
