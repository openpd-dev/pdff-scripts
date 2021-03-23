#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
file: job.py
created time : 2021/03/17
last edit time : 2021/03/18
author : Zhenyu Wei 
version : 1.0
contact : zhenyuwei99@gmail.com
copyright : (C)Copyright 2021-2021, Zhenyu Wei and Southeast University
'''

import paramiko, os
from ast import literal_eval
from pdff_distribute.recipe import Recipe
from pdff_distribute.device import Device

class Job:
    def __init__(self, recipe: Recipe) -> None:
        self.recipe = recipe
        self.device = None
        self.is_bound = False

    def bindDevice(self, device: Device):
        self.device = device
        self.is_bound = True
        self.remote_save_dir = os.path.join(self.device.target_dir, self.recipe.save_dir.split('/')[-1])

    def testBound(self):
        if self.is_bound == False:
            raise ValueError('Job has not been bound')

    def run(self):
        self.testBound()
        self._startClients()
        self._submitFiles()
        self._runSimulation()
        self._collectResults()
        self._backUpFiles()
        self._closeClients()

    def _startClients(self):
        self.ssh = paramiko.SSHClient()
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        self.ssh.connect(**self.device.getSSHInfo())
        transport = paramiko.Transport((self.device.hostname, 22))
        transport.connect(**self.device.getSFTPInfo())
        self.sftp = paramiko.SFTPClient.from_transport(transport)

    def _submitFiles(self):
        self.recipe.runRecipe()
        _, stdout, _ = self.ssh.exec_command(
            'rm -rf ' + self.remote_save_dir
        )
        stdout.read()
        for dirpath, dirnames, filenames in os.walk(self.recipe.save_dir):
            remote_dir_path = os.path.join(
                self.device.target_dir, dirpath[len(self.recipe.target_dir)+1:]
            ) # Get rid of the origin absolute dir path
            # make remote directory
            self.sftp.mkdir(remote_dir_path)
            for filename in filenames:
                local_path = os.path.join(dirpath, filename)
                remote_path = os.path.join(remote_dir_path, filename)
                # put file
                self.sftp.put(local_path, remote_path)

    def _runSimulation(self):
        stdin, stdout, stderr = self.ssh.exec_command(
            'export python=' + self.device.python_path +
            '&& cd ' + self.remote_save_dir +
            '&& chmod +x run.sh' +
            '&& ./run.sh'
        )
        out = stdout.read().decode('utf-8')
        if out != '':
            print(out)
        err = stderr.read().decode('utf-8')
        if err != '':
            print(err)

    def _collectResults(self):
        stdin, stdout, stderr = self.ssh.exec_command(
            'cd ' + os.path.join(self.remote_save_dir, self.recipe.collect_dir) +
            '''&& touch get.py && echo "import os \\nprint(list(os.walk('.')))" >> get.py''' +
            '&& ' + self.device.python_path + ' get.py && rm -rf get.py'
        )
        dir_data = literal_eval(stdout.read().decode('utf-8'))
        for dirpath, dirnames, filenames in dir_data:
            remote_dir_path = os.path.join(
                self.remote_save_dir, self.recipe.collect_dir, dirpath
            )
            local_dir_path = os.path.join(
                self.recipe.save_dir, self.recipe.collect_dir, dirpath
            ) # Get rid of the origin absolute dir path
            # make local directory
            if dirpath != '.':
                os.system('rm -rf ' + local_dir_path + '&& mkdir ' + local_dir_path)
            for filename in filenames:
                if filename == 'get.py':
                    continue
                local_path = os.path.join(local_dir_path, filename)
                remote_path = os.path.join(remote_dir_path, filename)
                # get file
                self.sftp.get(remote_path, local_path)
    
    def _backUpFiles(self):
        stdin, stdout, stderr = self.ssh.exec_command(
            'rm -rf ' + os.path.join(self.device.back_dir, self.recipe.save_dir.split('/')[-1]) + 
            '&& mv ' + self.remote_save_dir + ' ' + self.device.back_dir
        )
        stdout.read()

    def _closeClients(self):
        self.ssh.close()
        self.sftp.close()