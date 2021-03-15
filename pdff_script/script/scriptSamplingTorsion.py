import datetime
import simtk.openmm.app as app
from . import Script

class ScriptSamplingTorsion(Script):
    def __init__(
        self, save_dir: str, model_name: str, forcefield_file: str, cuda_id: int,
        file_name='05_sampling.py'
    ) -> None:
        super().__init__(save_dir, model_name, forcefield_file, cuda_id=cuda_id)
        self.file_name = file_name