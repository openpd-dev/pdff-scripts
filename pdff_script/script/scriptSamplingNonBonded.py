import datetime
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
# Sampling Variables
sample_frequency = 50                       # CV output frequency
sample_k = 50*unit.kilojoule_per_mole/unit.angstrom/unit.angstrom
sample_range = np.linspace(1, 15, 200)
sample_group1 = {:<}
sample_group2 = {:<}

sample_num_runs = len(sample_range)         # Number of biased simulation

centroid_distance = openmm.CustomCentroidBondForce(2, 'distance(g1,g2)')
centroid_distance.addGroup(sample_group1, [])
centroid_distance.addGroup(sample_group2, [])
centroid_distance.addBond([0, 1], [])
centroid_distance.setUsesPeriodicBoundaryConditions(True)
centroid_distance.setForceGroup(30)

bias_potential = openmm.CustomCentroidBondForce(2, '0.5*k*(distance(g1,g2)-r0)^2')
bias_potential.addGlobalParameter('k', sample_k)
bias_potential.addPerBondParameter('r0')
bias_potential.addGroup(sample_group1, [])
bias_potential.addGroup(sample_group2, [])
bias_potential.addBond([0, 1], [0])
bias_potential.setUsesPeriodicBoundaryConditions(True)
bias_potential.setForceGroup(31)

# Path Variables
path_current = sys.path[0]
path_out_pdb = path_current + '/../output/pdb_files'
path_out_log = path_current + '/../output/log_files'
path_out_cv = path_current + '/../output/cv_files'
if os.path.exists(path_out_cv):
    shutil.rmtree(path_out_cv)
os.mkdir(path_out_cv)

# Simulation Variables
sim_temp = 300 * unit.kelvin
sim_interval = 1 * unit.femtoseconds
sim_steps = 500000000
sim_cutoff = 12 * unit.angstroms
sim_platform = 'CUDA'

# Output Variables
out_log_interval = 100000
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
thermostat = openmm.AndersenThermostat(sim_temp, 0.001/unit.femtosecond)
system.addForce(thermostat)
system.addForce(bias_potential)
system.addForce(centroid_distance)

# Integrator
integrator = openmm.LangevinIntegrator(sim_temp, 0.001/unit.femtosecond, sim_interval) 
integrator.setIntegrationForceGroups({{0, 31}})

# Simulation
simulation = app.Simulation(pdb.topology, system, integrator, platform)

simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(sim_temp)

simulation.reporters.append(log_reporter)
simulation.reporters.append(pdb_reporter)

## Umbrella sampling
# Initial run for rearrangement of cv
state = simulation.context.getState(getPositions=True)
position = state.getPositions()
pbc = state.getPeriodicBoxVectors()
r_state = simulation.context.getState(getEnergy=True, groups={{30}})
print('Initial radius: %f an' %(r_state.getPotentialEnergy()/unit.kilojoule_per_mole*10))

bias_potential.setBondParameters(0, [0, 1], [sample_range[0]*unit.angstrom])
bias_potential.updateParametersInContext(simulation.context)
simulation.step(500000) # 0.5ns

for (i, bias) in enumerate(sample_range):
    sim_single_step = int(sim_steps/sample_num_runs)
    path_cv_file = os.path.join(path_out_cv, 'cv_%d.txt' %(i))
    cv_file = open(path_cv_file, 'a')
    cv_file.write('%f\\n' %(bias))
    #cv_force.setGlobalParameterDefaultValue(1, 20*unit.angstrom)
    bias_potential.setBondParameters(0, [0, 1], [bias*unit.angstrom])
    bias_potential.updateParametersInContext(simulation.context)
    num_output = int(sim_single_step/sample_frequency)
    print('%d steps in one run' %(num_output*sample_frequency))
    for i in range(num_output):
        simulation.step(sample_frequency)
        state = simulation.context.getState(getPositions=True)
        position = state.getPositions()
        pbc = state.getPeriodicBoxVectors()
        r_state = simulation.context.getState(getEnergy=True, groups={{30}})
        cv_file.write('%f\\n' %(r_state.getPotentialEnergy()/unit.kilojoule_per_mole*10))
        if i % 50 == 0:
            cv_file.close()
            cv_file = open(path_cv_file, 'a')

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


class ScriptSamplingNonBonded(Script):
    def __init__(
        self, save_dir: str, model_name: str, forcefield_file: str, cuda_id: int,
        pdb_file: str, file_name='05_sampling.py',
    ) -> None:
        super().__init__(save_dir, model_name, forcefield_file, cuda_id=cuda_id)
        self.file_name = file_name
        self.pdb_file = pdb_file
        
    def format_context(self):
        self.pdb_file = app.PDBFile(self.pdb_file)
        atoms = list(self.pdb_file.topology.atoms())
        res0_atoms = [atom for atom in atoms if atom.residue.index==0]
        res0_sc = [atom.index for atom in res0_atoms if not atom.name in BACK_BONE_ATOMS]

        res1_atoms = [atom for atom in atoms if atom.residue.index==1]
        res1_sc = [atom.index for atom in res1_atoms if not atom.name in BACK_BONE_ATOMS]
        self.context = context.format(
            self.model_name, str(datetime.datetime.now().replace(microsecond=0)), 
            str(res0_sc), str(res1_sc), self.forcefield_file, self.cuda_id
        )