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
path_out_pdb = path_current + '/../output/pdb_files'
path_out_log = path_current + '/../output/log_files'

# Simulation Variables
sim_temp = 300 * unit.kelvin
sim_pres = 1 * unit.bar
sim_interval = 1 * unit.femtoseconds
sim_steps = 4000000
sim_cutoff = 12 * unit.angstroms
sim_platform = 'CUDA'

# Output Variables
out_log_interval = 5000
out_pdb_interval = 5000 
out_prefix = 'eq_npt'

## Loading data
pdb = app.PDBFile(path_out_pdb + '/heating_nvt_restart.pdb')

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
barostat = openmm.MonteCarloBarostat(sim_pres, sim_temp, 1000)
system.addForce(barostat)

# Integrator
integrator = openmm.LangevinIntegrator(sim_temp, 0.001/unit.femtosecond, sim_interval)  

# Simulation
simulation = app.Simulation(pdb.topology, system, integrator, platform)

simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(sim_temp)

simulation.reporters.append(log_reporter)
simulation.reporters.append(pdb_reporter)

## Run
simulation.step(sim_steps)

## Restart Files
final_state = simulation.context.getState(getPositions=True, getParameters=True, enforcePeriodicBox=True)

file_restart = open(path_out_pdb + '/' + out_prefix + '_restart.pdb', 'w')
simulation.topology.setPeriodicBoxVectors(final_state.getPeriodicBoxVectors())
app.PDBFile.writeFile(simulation.topology, final_state.getPositions(), file_restart)

end_time = datetime.datetime.now()
print('Total running time:', end='\\t')
print(end_time - start_time)

file_log.close()
file_restart.close()
"""

class ScriptEqNPT(Script):
    def __init__(
        self, save_dir: str, model_name: str, forcefield_file: str, cuda_id: int,
        file_name='03_eq_npt.py'
    ) -> None:
        super().__init__(save_dir, model_name, forcefield_file, cuda_id=cuda_id)
        self.file_name = file_name
        self.context = context.format(
            self.model_name, str(datetime.datetime.now().replace(microsecond=0)), 
            self.forcefield_file, self.cuda_id
        )
