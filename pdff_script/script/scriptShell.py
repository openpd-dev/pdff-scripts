import datetime
from . import Script

context = """
#########################################################
## Model: {:<45}##
## Goal: Script for series of simulation               ##
## Author: Zhenyu Wei                                  ##
## Date: {:<46}##
#########################################################

export PID={:<}
while [ -e /proc/$PID ]
do 
    time=$(date "+%Y-%m-%d %H:%M:%S")
    echo "$time: Process: $PID is still running"
    sleep 300
done
time=$(date "+%Y-%m-%d %H:%M:%S")
echo "$time: Start Simulation"

cd simulation 
python 01_min.py
python 02_heating_nvt.py
python 03_eq_npt.py
python 04_eq_nvt.py
python 05_sampling.py

time=$(date "+%Y-%m-%d %H:%M:%S")
echo "$time: End Simulation"
"""

class ScriptShell(Script):
    def __init__(
        self, save_dir: str, model_name: str, file_name='run.sh', parent_id=0
    ) -> None:
        super().__init__(save_dir, model_name, file_name)
        self.parent_id = parent_id
        
    def format_context(self):
        self.context = context.format(
            self.model_name, str(datetime.datetime.now().replace(microsecond=0)), 
            self.parent_id
        ) 