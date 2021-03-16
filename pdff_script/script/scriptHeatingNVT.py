import datetime
from . import Script

context = """
#########################################################
## Model: {:<45}##
## Goal: Heating system to ideal temperature           ##
## Author: Zhenyu Wei                                  ##
## Date: {:<46}##
######################################################### 

## Presetting
import simtk.openmm.app as app
import simtk.openmm as openmm
import simtk.unit as unit
import sys, datetime
import numpy as np
start_time = datetime.datetime.now()

## Setting Variables
# Path Variables
path_current = sys.path[0]
path_out_pdb = path_current + '/../output/pdb_files'
path_out_log = path_current + '/../output/log_files'

# Simulation Variables
# temp, without unit as they will changed as a float in the for loop
sim_temp_start = 1
sim_temp_end = 300
# time
sim_time_single = 10 * unit.picoseconds # Single equilibrium trajectory in heating paradigm
sim_interval = 1 * unit.femtoseconds
# step
sim_steps = 5000000
# other
sim_cutoff = 12 * unit.angstroms
sim_platform = 'CUDA'

# Output Variables
out_log_interval = 5000
out_pdb_interval = 5000
out_prefix = 'heating_nvt'

## Loading data
pdb = app.PDBFile(path_out_pdb + '/min_restart.pdb')

## Setting Simulation
# In the heating paradigm, system, integrator will be changed periodically, 
# which is not including in `Setting Simulation` part, rather in `Run` part
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

## Run
sim_steps_single = sim_time_single / sim_interval
sim_temp_diff = sim_temp_end + 1 - sim_temp_start
num_iters = int(np.floor(sim_steps / sim_steps_single))

print('Start heating system from %.2f K to %.2f K\\n' %(sim_temp_start, sim_temp_end))
print('%d simulations of %d steps (%.2f ps) will be performed' %(num_iters, sim_steps_single, sim_time_single/unit.picosecond))

for temp in np.arange(sim_temp_start, sim_temp_end+1, sim_temp_diff/num_iters): 
    # Commented codes are only used while heating system in NPT ensemble
    print('\\nHeating system at %.2f K. \\n' %temp)  
    # System
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME, nonbondedCutoff=sim_cutoff, 
                constraints=app.HAngles, ewaldErrorTolerance=0.0005)
    #barostat = openmm.MonteCarloBarostat(sim_pres, sim_temp_start*unit.kelvin, 1000)
    #system.addForce(barostat)
    #thermostat = openmm.AndersenThermostat(temp*unit.kelvin, 0.001/unit.femtosecond)
    #system.addForce\\(thermostat)
    
    # Integrator      
    integrator = openmm.LangevinIntegrator(temp*unit.kelvin, 0.001/unit.femtosecond, sim_interval)  
    simulation = app.Simulation(pdb.topology, system, integrator, platform)
    #init_pbc = simulation.topology.getPeriodicBoxVectors()
    #simulation.context.setPeriodicBoxVectors(init_pbc[0], init_pbc[1], init_pbc[2])
    simulation.context.setPositions(pdb.positions)
    simulation.context.setVelocitiesToTemperature(temp*unit.kelvin)
    
    simulation.reporters.append(log_reporter)
    simulation.reporters.append(pdb_reporter)

    simulation.step(sim_steps_single)

    ## Restart Files
    final_state = simulation.context.getState(getPositions=True, getParameters=True, enforcePeriodicBox=True)
    
    file_restart = open(path_out_pdb + '/' + out_prefix + '_restart.pdb', 'w')
    #simulation.topology.setPeriodicBoxVectors(final_state.getPeriodicBoxVectors())
    
    app.PDBFile.writeFile(simulation.topology, final_state.getPositions(), file_restart)
    file_restart.close()
    
    pdb = app.PDBFile(path_out_pdb + '/' + out_prefix + '_restart.pdb')

end_time = datetime.datetime.now()
print('Total running time:', end='\\t')
print(end_time - start_time)

file_log.close()
"""

class ScriptHeatingNVT(Script):
    def __init__(
        self, save_dir: str, model_name: str, forcefield_file: str, cuda_id: int,
        file_name='02_heating_nvt.py'
    ) -> None:
        super().__init__(save_dir, model_name, forcefield_file, cuda_id=cuda_id)
        self.file_name = file_name
        self.context = context.format(
            self.model_name, str(datetime.datetime.now().replace(microsecond=0)), 
            self.forcefield_file, self.cuda_id
        )