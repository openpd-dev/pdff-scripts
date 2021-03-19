#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
file: device.py
created time : 2021/03/17
last edit time : 2021/03/18
author : Zhenyu Wei 
version : 1.0
contact : zhenyuwei99@gmail.com
copyright : (C)Copyright 2021-2021, Zhenyu Wei and Southeast University
'''

class Device:
    def __init__(
        self, hostname, username, password,
        cuda_id, python_path,
        target_dir, back_dir,
        is_local=False, label=''
    ) -> None:
        self.hostname = hostname
        self.username = username
        self.password = password
        self.cuda_id = cuda_id
        self.target_dir = target_dir
        self.back_dir = back_dir
        self.python_path = python_path
        if label == '':
            self.label = 'GPU %d of ' %self.cuda_id + self.username + '@' + self.hostname
        else:
            self.label = label 
        self.is_local = is_local

    def getSSHInfo(self):
        return {
            'hostname': self.hostname,
            'username': self.username,
            'password': self.password
        }
    
    def getSFTPInfo(self):
        return {
            'username': self.username,
            'password': self.password
        }

    def submitJob(self, job):
        job.bindDevice(self)
        job.run()