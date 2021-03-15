import pytest, os
from .. import RecipeTorsion

cur_dir = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))
forcefield_dir = os.path.join(cur_dir, '../forcefield/amberff14sb')

class TestScriptEqNVT:
    def setup(self):
        self.recipe = RecipeTorsion(
            save_dir=os.path.join(cur_dir, 'output/testRecipe')
        )

    def teardown(self):
        self.recipe = None

    def test_attributes(self):
        pass

    def test_exceptions(self):
        pass

    def test_createDirs(self):
        self.recipe.createDirs()