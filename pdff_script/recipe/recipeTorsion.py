#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
file: recipeTorsion.py
created time : 2021/03/16
last edit time : 2021/03/18
author : Zhenyu Wei 
version : 1.0
contact : zhenyuwei99@gmail.com
copyright : (C)Copyright 2021-2021, Zhenyu Wei and Southeast University
'''

import os
import numpy as np
from . import Recipe
from ..script import *
from .. import PDBManipulator, isStandardPeptide

cur_dir = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))
template_dir = os.path.join(cur_dir, '../template')

class RecipeTorsion(Recipe):
    def __init__(
        self, save_dir: str, model_name: str, cuda_id: int, peptide1: str, peptide2: str, 
        solution: str='nacl_0.15mol', parent_id=0, forcefield_file: str='amber14/torsion.xml'
    ) -> None:
        super().__init__(save_dir, model_name, forcefield_file, cuda_id)
        self.parent_id = parent_id
        
        # Dir structure
        self.dir_tree = {
            'simulation': {},
            'str': {},
            'output': {
                'log_files': {},
                'pdb_files': {},
                'meta_files': {'back'}
            }
        }

        # Structure
        _ = isStandardPeptide(peptide1, peptide2)
        self.peptide1 = PDBManipulator(os.path.join(template_dir, 'peptide', peptide1 + '.pdb'), end_label='ENDMDL')
        self.peptide2 = PDBManipulator(os.path.join(template_dir, 'peptide', peptide2 + '.pdb'), end_label='ENDMDL')
        self.solution = PDBManipulator(os.path.join(template_dir, 'solution', solution + '.pdb'), end_label='ENDMDL')

        # Simulation recipe
        self.script_recipe = [
            ScriptMinimize(
                save_dir=os.path.join(self.save_dir, 'simulation'), 
                forcefield_file=self.forcefield_file_name,
                model_name=self.model_name, cuda_id=self.cuda_id
            ),
            ScriptHeatingNVT(
                save_dir=os.path.join(self.save_dir, 'simulation'), 
                forcefield_file=self.forcefield_file_name,
                model_name=self.model_name, cuda_id=self.cuda_id
            ),
            ScriptEqNPT(
                save_dir=os.path.join(self.save_dir, 'simulation'), 
                forcefield_file=self.forcefield_file_name,
                model_name=self.model_name, cuda_id=self.cuda_id
            ),
            ScriptEqNVT(
                save_dir=os.path.join(self.save_dir, 'simulation'), 
                forcefield_file=self.forcefield_file_name,
                model_name=self.model_name, cuda_id=self.cuda_id
            ),
            ScriptSamplingTorsion(
                save_dir=os.path.join(self.save_dir, 'simulation'), 
                pdb_file=os.path.join(self.save_dir, 'str/str.pdb'),
                forcefield_file=self.forcefield_file_name,
                model_name=self.model_name, cuda_id=self.cuda_id
            ),
            ScriptShell(
                save_dir=os.path.join(self.save_dir), 
                model_name=self.model_name, parent_id=self.parent_id
            )
        ]
        
    def createStrFiles(self):
        self.peptide1.addPatch(0, 'NH3')
        self.peptide2.addPatch(0, 'COOH')
        self.peptide2.setResIdByResId(0, 1)
        coord_peptide1_n = self.peptide1.coord[[i for i, j in enumerate(self.peptide1.atom_name) if j == 'C'][0]]
        coord_peptide2_c = self.peptide1.coord[[i for i, j in enumerate(self.peptide2.atom_name) if j == 'N'][0]]
        move_vec = coord_peptide1_n - coord_peptide2_c + np.array([5, 5, 5])
        self.peptide2.moveBy(move_vec)
        self.structure = self.peptide1
        self.structure.catManipulators(self.peptide2, self.solution)
        self.structure.writeNewFile(os.path.join(self.save_dir, 'str/str.pdb'))

    