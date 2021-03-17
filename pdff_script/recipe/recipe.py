#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
file: recipe.py
created time : 2021/03/16
last edit time : 2021/03/17
author : Zhenyu Wei 
version : 1.0
contact : zhenyuwei99@gmail.com
copyright : (C)Copyright 2021-2021, Zhenyu Wei and Southeast University
'''

import os, shutil

cur_dir = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))
forcefield_dir = os.path.join(cur_dir, '../forcefield')
class Recipe:
    def __init__(
        self, save_dir: str, model_name: str, forcefield_file: str, cuda_id: int,
    ) -> None:
        self.save_dir = save_dir
        self.model_name = model_name
        self.forcefield_file = forcefield_file if 'forcefield/' in forcefield_file else os.path.join(forcefield_dir, forcefield_file)
        self.cuda_id = cuda_id

    def _createDirsFromFileTree(self, parent_dir, dir_tree):
        for key, value in list(dir_tree.items()):
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
        self._createDirsFromFileTree(self.save_dir, self.dir_tree)

    def createScriptFiles(self):
        for script in self.script_recipe:
            script.writeFile()

    def createStrFiles(self):
        raise NotImplementedError('createStrFile method has not been overloaded')

    def runRecipe(self):
        self.createDirs()
        self.createStrFiles()
        self.createScriptFiles()