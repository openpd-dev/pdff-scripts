import pytest, os
from .. import ScriptShell

cur_dir = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))

class TestScriptMinmize:
    def setup(self):
        self.script = ScriptShell(
            save_dir=os.path.join(cur_dir, 'output'), 
            model_name='alpha helix validation'
        )

    def teardown(self):
        self.script = None

    def test_attributes(self):
        pass

    def test_exceptions(self):
        pass

    def test_writeFile(self):
        self.script.writeFile()