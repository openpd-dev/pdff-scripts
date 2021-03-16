import pytest, os
from .. import ScriptSamplingNonBonded

cur_dir = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))
forcefield_dir = os.path.join(cur_dir, '../forcefield/amberff14sb')

class TestScriptSamplingTorsion:
    def setup(self):
        self.script = ScriptSamplingNonBonded(
            save_dir=os.path.join(cur_dir, 'output'), 
            forcefield_file=os.path.join(forcefield_dir, 'protein.ff14SB_torsion.xml'),
            pdb_file=os.path.join(cur_dir, 'data/testPDBManipulator.pdb'),
            model_name='alpha helix validation', cuda_id=0, file_name='05_sampling_nonbonded.py'
        )

    def teardown(self):
        self.script = None

    def test_attributes(self):
        pass

    def test_exceptions(self):
        pass

    def test_writeFile(self):
        self.script.writeFile()