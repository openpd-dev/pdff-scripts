__author__ = "Zhenyu Wei"
__maintainer__ = "Zhenyu Wei" 
__email__ = "zhenyuwei99@gmail.com"
__copyright__ = "Copyright 2021-2021, Southeast University and Zhenyu Wei"
__license__ = "GPLv3"

MASS_DICT = {
    "C": 12,
    "H": 1,
    "O": 16,
    "N": 14,
    "S": 32
}

BACK_BONE_ATOMS = [
    "C", "O", "OXT",
    "CA", "HA", "HA2",
    "N", "H1", "H2", "H3", "H"
]

from .script import Script
from .scriptMinimize import ScriptMinimize
from .scriptHeatingNVT import ScriptHeatingNVT
from .scriptEqNVT import ScriptEqNVT
from .scriptEqNPT import ScriptEqNPT

__all__ = [
    'ScriptMinimize',
    'ScriptHeatingNVT',
    'ScriptEqNVT',
    'ScriptEqNPT'
]