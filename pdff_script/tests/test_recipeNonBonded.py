import pytest, os
from .. import RecipeNonBonded

cur_dir = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))
forcefield_dir = os.path.join(cur_dir, '../forcefield/amberff14sb')

class TestRecipeTorsion:
    def setup(self):
        self.recipe = RecipeNonBonded(
            save_dir=os.path.join(cur_dir, 'output/testRecipeNonBonded'),
            forcefield_file=os.path.join(forcefield_dir, 'protein.ff14SB_nonbonded.xml'),
            model_name='PDFF alpha helix validation', cuda_id=0,
            peptide1='ASN', peptide2='TYR'
        )

    def teardown(self):
        self.recipe = None

    def test_attributes(self):
        pass

    def test_exceptions(self):
        pass

    def test_createDirs(self):
        self.recipe.createDirs()

    def test_createStrFiles(self):
        self.recipe.createDirs()
        self.recipe.createStrFiles()

    def test_createScriptFiles(self):
        self.recipe.createDirs()
        self.recipe.createStrFiles()
        self.recipe.createScriptFiles()

    def test_runRecipe(self):
        self.recipe.runRecipe()