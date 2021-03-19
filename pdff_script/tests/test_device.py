import pytest
from .. import Device

class TestDevice:
    def setup(self):
        self.device = Device(
            hostname='10.201.129.123', username='zhenyuwei', password='Wzy99714!',
            cuda_id=0, python_path='/home/zhenyuwei/.conda/envs/openmm/bin/python',
            target_dir='/home/zhenyuwei/simulation_data/openmm/pdff',
            back_dir='/home/zhenyuwei/simulation_data/openmm/pdff_back'
        )

    def teardown(self):
        pass

    def test_attributes(self):
        assert self.device.hostname == '10.201.129.123'
        assert self.device.username == 'zhenyuwei'

    def test_exceptions(self):
        pass

    def test_getSSHInfo(self):
        info = self.device.getSSHInfo()
        assert info['hostname'] == '10.201.129.123'
        assert info['username'] == 'zhenyuwei'
        assert info['password'] == 'Wzy99714!'

    def test_getSFTPInfo(self):
        info = self.device.getSFTPInfo()
        assert info['username'] == 'zhenyuwei'
        assert info['password'] == 'Wzy99714!'
        with pytest.raises(KeyError):
            info['hostname']