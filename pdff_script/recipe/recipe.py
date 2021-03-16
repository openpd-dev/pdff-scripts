import os, shutil

class Recipe:
    def __init__(
        self, save_dir: str, model_name: str, forcefield_file: str, cuda_id: int,
    ) -> None:
        self.save_dir = save_dir
        self.model_name = model_name
        self.forcefield_file = forcefield_file
        self.cuda_id = cuda_id

    def _createDirsFromFileTree(self, parent_dir, file_tree):
        for key, value in list(file_tree.items()):
            os.mkdir(os.path.join(parent_dir, key))
            if isinstance(value, dict) and value != {}:
                self._createDirsFromFileTree(os.path.join(parent_dir, key), value)
            elif isinstance(value, set):
                for i in list(value):
                    os.mkdir(os.path.join(parent_dir, key, i))
                
    def createDirs(self):
        if os.path.exists(self.save_dir):
            shutil.rmtree(self.save_dir)
        os.mkdir(self.save_dir)
        self._createDirsFromFileTree(self.save_dir, self.file_tree)

    def createSimulationFiles(self):
        for script in self.simulation_recipe:
            script.writeFile()

    def createStrFiles(self):
        raise NotImplementedError('createStrFile method has not been overloaded')

    def runRecipe(self):
        self.createDirs()
        self.createStrFiles()
        self.createSimulationFiles()