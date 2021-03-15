import datetime
from . import Script

context = """
#########################################################
## Model: {:<45}##
## Goal: Minimize system energy                        ##
## Author: Zhenyu Wei                                  ##
## Date: {:<46}##
######################################################### 

## Presetting
import simtk.openmm.app as app
import simtk.openmm as openmm
import simtk.unit as unit
import sys, datetime
start_time = datetime.datetime.now()

## Setting Variables
# Path Variables
path_current = sys.path[0]
path_input = path_current + '/../str'
path_out_pdb = path_current + '/../output/pdb_files'
path_out_log = path_current + '/../output/log_files'

# Simulation Variables
sim_temp = 300 * unit.kelvin
sim_cutoff = 10 * unit.angstroms
sim_interval = 2 * unit.femtoseconds

## Loading data
pdb = app.PDBFile(path_input + '/str.pdb')

## Seeting Simulation
# Forcefield
forcefield = app.ForceField('{:<}', 'amber14/tip3pfb.xml')

# Integrator
integrator = openmm.LangevinIntegrator(sim_temp, 0.001/unit.femtosecond, sim_interval)

# Platform
platform = openmm.Platform_getPlatformByName('CUDA')
platform.setPropertyDefaultValue('DeviceIndex', '{:<}')

# System
system = forcefield.createSystem(pdb.topology, 
            nonbondedMethod=app.PME, nonbondedCutoff=sim_cutoff, 
            constraints=app.HBonds)
            
# Simulation
simulation = app.Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)

## Output
# Only log info in min.py
file_log = open(path_out_log + '/min.log', 'w')
sys.stdout = file_log
init_state = simulation.context.getState(getEnergy=True)
print('start to minimize system energy')
print('Initial potential energy:\\t%.2f\\t kj/mol' %(init_state.getPotentialEnergy()/unit.kilojoule_per_mole))

## Run
simulation.minimizeEnergy(tolerance=1e-4 * unit.kilocalorie_per_mole, maxIterations=1000)
#simulation.minimizeEnergy()

## Restart Files
final_state = simulation.context.getState(getPositions=True, getParameters=True, enforcePeriodicBox=True, getEnergy=True)

print('Final potential energy:\\t\\t%.2f\\t kj/mol' %(final_state.getPotentialEnergy()/unit.kilojoule_per_mole))

file_restart = open(path_out_pdb + '/min_restart.pdb', 'w')
app.PDBFile.writeFile(simulation.topology, final_state.getPositions(), file_restart)

end_time = datetime.datetime.now()
print('Total running time:', end='\\t')
print(end_time - start_time)

file_log.close()
file_restart.close()
"""

class ScriptMinimize(Script):
    def __init__(
        self, save_dir: str, model_name: str, forcefield_file: str, cuda_id: int,
        file_name='01_min.py'
    ) -> None:
        super().__init__(save_dir, model_name, forcefield_file, cuda_id=cuda_id)
        self.file_name = file_name
        self.context = context.format(
            self.model_name, str(datetime.datetime.now()), 
            self.forcefield_file, self.cuda_id
        )


