import os, pytest
from pdff_script.recipe import recipeNonBonded
from .. import Device, Job
from .. import RecipeNonBonded

cur_dir = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))

class TestJob:
    def setup(self):
        recipe = RecipeNonBonded(
            target_dir=os.path.join(cur_dir, 'output'),
            model_name='alpha helix validation', cuda_id=0,
            peptide1='ASN', peptide2='TRP'
        )
        self.device = Device(
            hostname='10.201.129.123', username='zhenyuwei', password='Wzy99714!',
            cuda_id=0, python_path='/home/zhenyuwei/.conda/envs/openmm/bin/python',
            target_dir='/home/zhenyuwei/simulation_data/openmm/pdff',
            back_dir='/home/zhenyuwei/simulation_data/openmm/pdff_back'
        )    
        self.job = Job(recipe)

    def teardown(self):
        self.device = None
        self.job = None

    def test_attributes(self):
        assert self.job.device == None
        assert self.job.is_bound == False

    def test_exceptions(self):
        with pytest.raises(ValueError):
            self.job.testBound()

    def test_bindDevice(self):
        self.job.bindDevice(self.device)
        assert self.job.is_bound == True
        assert self.job.device != None

    def test_startClients(self):
        self.job.bindDevice(self.device)
        self.job._startClients()
        _, stdout, _ = self.job.ssh.exec_command('ls')
        res = stdout.read().decode('utf-8').split('\n')
        assert 'Arduino' in res
        self.job._closeClients()

    def test_submitFiles(self):
        self.job.bindDevice(self.device)
        self.job._startClients()
        self.job._submitFiles()
        _, stdout, _ = self.job.ssh.exec_command(
            'cd ' + self.job.remote_save_dir + '/..' + 
            ' && ls'
        )
        res = stdout.read().decode('utf-8').split('\n')
        self.job._closeClients()

    def test_runSimulation(self):
        self.job.bindDevice(self.device)
        self.job._startClients()
        self.job._submitFiles()
        stdin, stdout, stderr = self.job.ssh.exec_command(
            'export python=' + self.device.python_path +
            '&& cd ' + self.job.remote_save_dir +
            '&& cd simulation' +
            '&& $python 01_min.py'
        )
        assert stderr.read().decode('utf-8') == ''
        self.job._closeClients()

    def test_collectResults(self):
        self.job.bindDevice(self.device)
        self.job._startClients()
        self.job._submitFiles()
        stdin, stdout, stderr = self.job.ssh.exec_command(
            'cd ' + os.path.join(self.job.remote_save_dir, self.job.recipe.collect_dir) + 
            '&& mkdir test' + '&& touch cv_1.txt' + '&& mkdir test/test_sub'
        )
        stdout.read()
        self.job._collectResults()
        listdir = os.listdir(os.path.join(self.job.recipe.save_dir, self.job.recipe.collect_dir))
        assert 'test' in listdir
        assert 'cv_1.txt' in listdir
        listdir = os.listdir(os.path.join(self.job.recipe.save_dir, self.job.recipe.collect_dir, 'test'))
        assert 'test_sub' in listdir
        self.job._closeClients()

    def test_backUpFiles(self):
        self.job.bindDevice(self.device)
        self.job._startClients()
        self.job._submitFiles()
        self.job._backUpFiles()
        self.job._closeClients()