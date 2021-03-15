import pytest, os
from .. import PDBManipulator

cur_dir = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))

class TestPDBManipulator:
    def setup(self):
        self.manipulator = PDBManipulator(os.path.join(cur_dir, 'data/testPDBManipulator.pdb'), 'END')

    def teardown(self):
        self.manipulator = None

    def test_attributes(self):
        assert self.manipulator.patches['NTER']['ParentAtom'] == 'N'
        assert self.manipulator.atom_id[0] == 136

    def test_exceptions(self):
        with pytest.raises(KeyError):
            self.manipulator._getMass('Z')

        with pytest.raises(KeyError):
            self.manipulator.addPatch(0, 'NT')

    def test_sortAtomId(self):
        self.manipulator.sortAtomId()
        for i, atom_id in enumerate(self.manipulator.atom_id):
            assert atom_id == i

    def test_sortResId(self):
        cur_res = self.manipulator.res_id[0][0]
        cur_order = 0
        sorted_id = []
        for res_id in self.manipulator.res_id:
            if res_id[0] == cur_res:
                sorted_id.append(cur_order)
            else:
                cur_res = res_id[0]
                cur_order += 1
                sorted_id.append(cur_order)
        self.manipulator.sortResId()
        for i, j in zip(sorted_id, self.manipulator.res_id):
            assert i == j[0]

    def test_addPatch(self):
        self.manipulator.sortResId()
        self.manipulator.sortAtomId()
        self.manipulator.addPatch(0, 'NH3')
        self.manipulator.addPatch(1, 'COOH')
        self.manipulator.writeNewFile(os.path.join(cur_dir, 'output/patched.pdb'))

    def test_writeNewFile(self):
        self.manipulator.sortResId()
        self.manipulator.sortAtomId()
        self.manipulator.writeNewFile(os.path.join(cur_dir, 'output/sorted.pdb'))