#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
file: recipeNonBonded.py
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
from .. import BACK_BONE_ATOMS, PDBManipulator, isStandardPeptide

cur_dir = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))
template_dir = os.path.join(cur_dir, '../template')

class RecipeNonBonded(Recipe):
    def __init__(
        self, target_dir: str, model_name: str, cuda_id: int, peptide1: str, peptide2: str, 
        solution: str='nacl_0.15mol', forcefield_file: str='amber14/nonbonded.xml'
    ) -> None:
        super().__init__(target_dir, model_name, forcefield_file, cuda_id)
        _ = isStandardPeptide(peptide1, peptide2)
        self.save_dir = os.path.join(self.target_dir, peptide1+'-'+peptide2)

        # Dir structure
        self.dir_tree = {
            'simulation': {},
            'str': {},
            'output': {
                'log_files',
                'pdb_files',
                'cv_files'
            }
        }

        # Structure
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
            ScriptSamplingNonBonded(
                save_dir=os.path.join(self.save_dir, 'simulation'), 
                pdb_file=os.path.join(self.save_dir, 'str/str.pdb'),
                forcefield_file=self.forcefield_file_name,
                model_name=self.model_name, cuda_id=self.cuda_id
            ),
            ScriptShell(
                save_dir=os.path.join(self.save_dir), 
                model_name=self.model_name
            )
        ]
        
    def createStrFiles(self):
        index_peptide1_sc = [i for i in range(self.peptide1.num_atoms) if not self.peptide1.atom_name[i] in BACK_BONE_ATOMS]
        coord_peptide1_com = np.zeros(3)
        mass_sc = 0
        for index in index_peptide1_sc:
            mass_sc += self.peptide1.mass[index]
            coord_peptide1_com += self.peptide1.mass[index] * self.peptide1.coord[index]
        coord_peptide1_com /= mass_sc

        index_peptide2_sc = [i for i in range(self.peptide2.num_atoms) if not self.peptide2.atom_name[i] in BACK_BONE_ATOMS]
        coord_peptide2_com = np.zeros(3)
        mass_sc = 0
        for index in index_peptide2_sc:
            mass_sc += self.peptide2.mass[index]
            coord_peptide2_com += self.peptide2.mass[index] * self.peptide2.coord[index]
        coord_peptide2_com /= mass_sc
        self.peptide2.setResIdByResId(0, 1)
        self.peptide2.setChainNameByResId(1, 'B')

        move_vec = coord_peptide1_com - coord_peptide2_com + np.array([5, 5, 5])
        self.peptide2.moveBy(move_vec)
        self.structure = self.peptide1
        self.structure.catManipulators(self.peptide2, self.solution)
        self.structure.writeNewFile(os.path.join(self.save_dir, 'str/str.pdb'))

    