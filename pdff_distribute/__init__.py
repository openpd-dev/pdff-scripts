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

TRIPLE_LETTER_ABBREVIATION = [
    'ALA', 'ARG', 'ASN', 'ASP',
    'CYS', 'GLN', 'GLU', 'GLY',
    'HIS', 'ILE', 'LEU', 'LYS',
    'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL'
]

SINGLE_LETTER_ABBREVIATION = [
    'A', 'R', 'N', 'D',
    'C', 'Q', 'E', 'G',
    'H', 'I', 'L', 'K',
    'M', 'F', 'P', 'S',
    'T', 'W', 'Y', 'V'
]

def isStandardPeptide(*peptides):
    for peptide in peptides:
        peptide = peptide.upper()
        if len(peptide) == 3:
            if not peptide in TRIPLE_LETTER_ABBREVIATION:
                raise ValueError(
                    'Peptide type %s is not in the standard peptide list:\n %s' 
                    %(peptide, TRIPLE_LETTER_ABBREVIATION)
                )
        elif len(peptide) == 1:
            if not peptide in SINGLE_LETTER_ABBREVIATION:
                raise ValueError(
                    'Peptide type %s is not in the standard peptide list:\n %s' 
                    %(peptide, SINGLE_LETTER_ABBREVIATION)
                )
        else:
            raise ValueError(
                'Peptide type should be 1 or 3 letters, instead of %d letters'
                %(len(peptide))
            )
    return True
    
from pdff_distribute.pdbManipulator import PDBManipulator
from pdff_distribute.script import *
from pdff_distribute.recipe import *
from pdff_distribute.device import Device
from pdff_distribute.job import Job
from pdff_distribute.manager import *