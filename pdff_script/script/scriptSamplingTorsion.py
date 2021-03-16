import datetime, os
import simtk.openmm.app as app
from . import Script
from .. import BACK_BONE_ATOMS

context = """
#########################################################
## Model: {:<45}##
## Goal: Sampling system in NVT ensmbles               ##
## Author: Zhenyu Wei                                  ##
## Date: {:<46}##
######################################################### 

## Presetting
import simtk.openmm.app as app
import simtk.openmm as openmm
import simtk.unit as unit
import os, sys, datetime, shutil
import numpy as np
start_time = datetime.datetime.now()

## Setting Variables
# Collective Variables
# The collective variables is obtained through a smart way under the project construction:
# We define a custom force = theta instead of k*theta like usual, in which way the variable can be obtained easily.

# group1 side chain center of the first peptide
sample_group1 = {:<}
# group2 ca of the first peptide
sample_group2 = {:<}
# group3 ca of the second peptide
sample_group3 = {:<}
# group4 side chain center of the second peptide
sample_group4 = {:<}

compound_force = openmm.CustomCentroidBondForce(4, 'dihedral(g1, g2, g3, g4)')
compound_force.addGroup(sample_group1, [])
compound_force.addGroup(sample_group2, [])
compound_force.addGroup(sample_group3, [])
compound_force.addGroup(sample_group4, [])
compound_force.addBond([0, 1, 2, 3])
#compound_force.setUsesPeriodicBoundaryConditions(True)
compound_force.setForceGroup(31)

bias_variable = app.BiasVariable(compound_force, -np.pi, np.pi, 0.1, True)

# MetaDynamics Varibile
meta_bias_factor = 2
meta_enegy_height = 0.01 * unit.kilojoules_per_mole
meta_add_frequency = 100
meta_save_frequency = 1000
meta_back_frequency = 5000 # Backup while $meta_back_frequency files has been saved
meta_save_id = 0

# Path Variables
path_current = sys.path[0]
path_out_pdb = path_current + '/../output/pdb_files'
path_out_log = path_current + '/../output/log_files'
path_out_metaD = path_current + '/../output/meta_files'
if os.path.exists(path_out_metaD):
    shutil.rmtree(path_out_metaD)
os.mkdir(path_out_metaD)
os.mkdir(os.path.join(path_out_metaD, 'back'))

# Simulation Variables
sim_temp = 300 * unit.kelvin
sim_interval = 1 * unit.femtoseconds
sim_steps = 200000000
sim_cutoff = 12 * unit.angstroms
sim_platform = 'CUDA'

# Output Variables
out_log_interval = 10000
out_pdb_interval = 10000 
out_prefix = 'sampling'

## Loading data
pdb = app.PDBFile(path_out_pdb + '/eq_nvt_restart.pdb')

## Setting simulation
# Forcefiled
forcefield = app.ForceField('{:<}', 'amber14/tip3pfb.xml')

# Platform
platform = openmm.Platform_getPlatformByName(sim_platform)
platform.setPropertyDefaultValue('DeviceIndex', '{:<}')

# Reporter
# Log reporter
file_log = open(path_out_log + '/' + out_prefix + '.log', 'w')
sys.stdout = file_log
log_reporter = app.StateDataReporter(sys.stdout, out_log_interval, step=True,
        potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
        temperature=True, volume=True, speed=True, density=True,
        totalSteps=sim_steps, remainingTime=True,separator='\\t')
# PDB reporter
pdb_reporter = app.PDBReporter(path_out_pdb + '/' + out_prefix + '.pdb', out_pdb_interval, enforcePeriodicBox=True)

# System
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME,nonbondedCutoff=sim_cutoff, 
                    constraints=app.HAngles, ewaldErrorTolerance=0.0005)

# Integrator
integrator = openmm.LangevinIntegrator(sim_temp, 0.001/unit.femtosecond, sim_interval) 

# MetaDynamics
# This part should always beyond the Simulation part.

# metaD = app.Metadynamics(system, [bias_variable], sim_temp, meta_bias_factor, meta_enegy_height, meta_add_frequency, meta_save_frequency, path_out_metaD)
metaD = app.Metadynamics(system, [bias_variable], sim_temp, meta_bias_factor,
    meta_enegy_height, meta_add_frequency, meta_save_frequency, path_out_metaD, 
    meta_save_id, meta_back_frequency)
    
# Simulation
simulation = app.Simulation(pdb.topology, system, integrator, platform)

simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(sim_temp)

simulation.reporters.append(log_reporter)
simulation.reporters.append(pdb_reporter)

## Run
metaD.step(simulation, sim_steps)

## Restart Files
final_state = simulation.context.getState(getPositions=True, getParameters=True, enforcePeriodicBox=True)

file_restart = open(path_out_pdb + '/' + out_prefix + '_restart.pdb', 'w')
simulation.topology.setPeriodicBoxVectors(final_state.getPeriodicBoxVectors())
app.PDBFile.writeFile(simulation.topology, final_state.getPositions(), file_restart)


end_time = datetime.datetime.now()
print('Total running time:', end='\\t')
print(end_time - start_time)

file_restart.close()
file_log.close()
"""


class ScriptSamplingTorsion(Script):
    def __init__(
        self, save_dir: str, model_name: str, forcefield_file: str, cuda_id: int,
        pdb_file: str, file_name='05_sampling.py',
    ) -> None:
        super().__init__(save_dir, model_name, file_name)
        self.forcefield_file = forcefield_file
        self.cuda_id = cuda_id
        self.file_name = file_name
        self.pdb_file = pdb_file
        
    def format_context(self):
        self.pdb_file = app.PDBFile(self.pdb_file)
        atoms = list(self.pdb_file.topology.atoms())
        res0_atoms = [atom for atom in atoms if atom.residue.index==0]
        res0_ca = [atom.index for atom in res0_atoms if atom.name=='CA']
        res0_sc = [atom.index for atom in res0_atoms if not atom.name in BACK_BONE_ATOMS]

        res1_atoms = [atom for atom in atoms if atom.residue.index==1]
        res1_ca = [atom.index for atom in res1_atoms if atom.name=='CA']
        res1_sc = [atom.index for atom in res1_atoms if not atom.name in BACK_BONE_ATOMS]
        
        self.context = context.format(
            self.model_name, str(datetime.datetime.now().replace(microsecond=0)), 
            str(res0_sc), str(res0_ca), str(res1_ca), str(res1_sc),
            self.forcefield_file, self.cuda_id
        )