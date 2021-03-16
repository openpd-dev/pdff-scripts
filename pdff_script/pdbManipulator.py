import os, codecs, json
import numpy as np
from numpy.testing._private.utils import raises 
from . import MASS_DICT

cur_dir = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))
template_dir = os.path.join(cur_dir, 'template/patch')


class PDBManipulator(object):
    def __init__(self, file_name, end_label='TER'):
        template_file = os.path.join(template_dir, 'patches.json')
        with codecs.open(template_file, 'r', 'utf-8') as f:
            template_text = f.read()
        self.patches = json.loads(template_text)
        self.file_name = file_name
        self.file = open(self.file_name, 'r')
        line = self.file.readline()
        line = self._readNoAtom(line)
        line = self._readATOM(line, end_label)
        self._setPDBInfo()
        self.sortAtomId()
        self.sortResId()
        
    def _skipBlankLine(self, line):
        while line.startswith('\n') or line.startswith('#'):
            line = self.file.readline()
        return line

    def _skipNoAtomLine(self, line):
        while not line.startswith('ATOM'):
            line = self.file.readline()
        return line
    
    def _parseLine(self, line):
        data = line.split()
        for (i, value) in enumerate(data):
            try:
                value = float(value)
            except:
                data[i] = value
            else:
                data[i] = value
        return data

    def _parseAtomLine(self, line):
        data = []
        data.append(int(line[6:11]))
        data.append(line[12:16].strip())
        data.append(line[17:20].strip())
        data.append(line[21])
        data.append(int(line[22:26]))
        data.append(float(line[30:38]))
        data.append(float(line[38:46]))
        data.append(float(line[46:54]))
        return data

    def _readNoAtom(self, line):
        line = self._skipBlankLine(line)
        self.line_NoATOM = []
        while not line.startswith('ATOM'):
            self.line_NoATOM.append(line)
            line = self.file.readline()
            #self.pbc_inv = np.linalg.inv(self.pbc_diag)
        return line

    def _readATOM(self, line, end_label):
        line = self._skipBlankLine(line)
        line = self._skipNoAtomLine(line)
        self.raw_data = []
        self.line_ATOM = []
        while not line.startswith(end_label):
            self.raw_data.append(self._parseAtomLine(line))
            self.line_ATOM.append(line)
            line = self.file.readline()
        self.line_ATOMEND = line
        return line

    def _getMass(self, atom_name):
        keys = MASS_DICT.keys()
        atom_name = atom_name.upper()
        if len(atom_name) >= 2:
            if atom_name[0:2] in keys:
                return MASS_DICT[atom_name[0:2]]
            elif atom_name[0] in keys:
                return MASS_DICT[atom_name[0]]
            else:
                raise KeyError('Atom %s or Atom %s is not in element_mass %s: ' 
                        %(atom_name[0:2], atom_name[0], keys))
        else:
            if atom_name[0] in keys:
                return MASS_DICT[atom_name[0]]
            else:
                raise KeyError('Atom %s is not in element_mass %s: ' 
                        %(atom_name[0], keys))
                    
    def _setPDBInfo(self):
        self.num_atoms = len(self.raw_data)
        self.atom_id = np.ones(self.num_atoms, int)
        self.res_id = np.ones(self.num_atoms, int)
        self.mass = np.ones(self.num_atoms)
        self.coord = np.ones([self.num_atoms, 3])
        self.atom_name = []
        self.res_name = []
        self.chain_name = []
        for (i, data) in enumerate(self.raw_data):
            self.atom_id[i] = int(data[0])
            self.atom_name.append(data[1])
            self.res_name.append(data[2])
            self.chain_name.append(data[3])
            self.res_id[i] = int(data[4])
            self.coord[i, :] = data[5:8]
            self.mass[i] = self._getMass(data[1])
            self.raw_data[i].append(self.mass[i])
        self.num_res = np.unique(self.res_id).shape[0]
        self.num_atoms = len(self.atom_name)

    def _replaceATOMLineChainName(self, line, chain_name: str):
        if len(chain_name) != 1:
            raise ValueError('The length of chain name should be 1')
        return (line[:21] + '%1s' %(chain_name) + line[22:])

    def _replaceATOMLineResName(self, line, res_name):
        return (line[:17] + '%3s' %(res_name) + line[20:])

    def _replaceATOMLineResId(self, line, res_id):
        return (line[:22] + '%4d' %(res_id) + line[26:])

    def _replaceATOMLineAtomName(self, line, atom_name):
        return (line[:12] + '%4s' %(atom_name) + line[16:])

    def _replaceATOMLineAtomId(self, line, atom_id):
        return (line[:6] + '%5d' %(atom_id) + line[11:])

    def _replaceATOMLineCoord(self, line, coord):
        return (line[:30] + '%8.2f' %coord[0] + 
            '%8.2f' %coord[1] + 
            '%8.2f' %coord[2] +  line[54:])

    def _insertATOMLine(self, index, line):
        if index >= len(self.line_ATOM):
            self.line_ATOM.append(line)
        else:
            line_new = self.line_ATOM[:index]
            line_new.append(line)
            line_new.extend(self.line_ATOM[index:])
            self.line_ATOM = line_new
    
    def _insertAtomId(self, index, atom_id):
        if index >= self.atom_id.shape[0]:
            self.atom_id = np.append(self.atom_id, atom_id)
        else:
            self.atom_id = np.insert(self.atom_id, atom_id, self.atom_id[index])
            self.atom_id[index+1:] += 1
    
    def _insertAtomName(self, index, atom_name):
        if index >= len(self.atom_name):
            self.atom_name.append(atom_name)
        else:
            atom_name_list = self.atom_name[:index]
            atom_name_list.append(atom_name)
            atom_name_list.extend(self.atom_name[index:])
            self.atom_name = atom_name_list

    def _insertResId(self, index, res_id):
        if index >= self.res_id.shape[0]:
            self.res_id = np.append(self.res_id, res_id)
        else:
            self.res_id = np.insert(self.res_id, index, res_id)

    def _insertResName(self, index, res_name):
        if index >= len(self.res_name):
            self.res_name.append(res_name)
        else:
            res_name_list = self.res_name[:index]
            res_name_list.append(res_name)
            res_name_list.extend(self.res_name[index:])
            self.res_name = res_name_list
        
    def _insertMass(self, index, mass):
        if index >= self.mass.shape[0]:
            self.mass = np.append(self.mass, mass)
        else:
            self.mass = np.insert(self.mass, index, mass)

    def _insertCoord(self, index, coord):
        if index >= self.coord.shape[0]:
            self.coord = np.append(self.coord, coord)
        else:
            self.coord = np.insert(self.coord, index, coord)
        self.coord = self.coord.reshape(-1, 3)
    
    def _insertChainName(self, index, chain_name):
        if index >= len(self.chain_name):
            self.chain_name.append(chain_name)
        else:
            chain_name_list = self.chain_name[:index]
            chain_name_list.append(chain_name)
            chain_name_list.extend(self.chain_name[index:])
            self.chain_name = chain_name_list
                
    def updateATOMLine(self):
        for (i, line) in enumerate(self.line_ATOM[:-1]):
            line = self._replaceATOMLineResName(line, self.res_name[i])
            line = self._replaceATOMLineAtomName(line, self.atom_name[i])
            line = self._replaceATOMLineResId(line, self.res_id[i])
            self.line_ATOM[i] = line

    def sortAtomId(self):
        for atom in range(self.num_atoms):
            self.atom_id[atom] = atom
            self.line_ATOM[atom] = self._replaceATOMLineAtomId(self.line_ATOM[atom], atom)

    def sortResId(self):
        cur_res = self.res_id[0]
        cur_order = 0
        for i, res_id in enumerate(self.res_id):
            if res_id == cur_res:
                self.res_id[i] = cur_order
                self.line_ATOM[i] = self._replaceATOMLineResId(self.line_ATOM[i], cur_order)
            else:
                cur_res = res_id
                cur_order += 1
                self.res_id[i] = cur_order
                self.line_ATOM[i] = self._replaceATOMLineResId(self.line_ATOM[i], cur_order)
        self.num_res = self.res_id[-1]
    
    def addPatch(self, res_id, patch_name, is_verbose=False):
        if not patch_name in self.patches:
            raise KeyError('Patch %s are not supported')

        if is_verbose:
            print('Info) Start to add patch %s' %(patch_name))
        atom_indexes = [i[0] for i in np.argwhere(self.res_id == res_id)]
        atom_name = [self.atom_name[i] for i in atom_indexes]
        res_name = self.res_name[atom_indexes[0]]
        chain_name = self.chain_name[atom_indexes[0]]
        patch = self.patches[patch_name]
        num_matched_atoms = 0
        for patch_atom in patch['PatchAtoms']:
            if patch_atom in atom_name:
                num_matched_atoms += 1
                if is_verbose:
                    print('Info) Atom %s of Patch %s has already existed in Residue %s' %(patch_atom, patch_name, res_name))
        if num_matched_atoms != 0:
            if num_matched_atoms == len(patch['PatchAtoms']):
                if is_verbose:
                    print('Info) Patch %s has already existed, all atoms are matched\n' %(patch_name))
                return
            else:
                raise KeyError('Error) Parts of atoms of patch %s are in Residue %s ID=%d, please check patches or residue' %(patch_name, res_name, res_id))
        
        parent_index = atom_indexes[atom_name.index(patch['ParentAtom'])]
        label_index = atom_indexes[atom_name.index(patch['LabelAtom'])]
        patch_res_id = res_id
        for (i, atom) in enumerate(patch['PatchAtoms']):
            if is_verbose:
                print('Info) Add Atom %s' %(atom))
            index = atom_indexes[-1] + i + 1
            coord = self.coord[parent_index] + np.random.rand(3)*2 - 1
            self._insertATOMLine(atom_indexes[-1]+i+1, patch['PatchLines'][i] 
                    %(self.atom_id[index-1]+1+i, res_name, chain_name, res_id, coord[0], coord[1], coord[2])) 
            self._insertAtomId(index, self.atom_id[index-1]+1+i)
            self._insertAtomName(index, atom)
            self._insertResId(index, res_id)
            self._insertResName(index, res_name)
            self._insertChainName(index, chain_name)
            self._insertCoord(index, coord)
            self._insertMass(index, self._getMass(atom))
            self.num_atoms += 1

        atom_indexes = [i[0] for i in np.argwhere(self.res_id == res_id)]
        atom_name = [self.atom_name[i] for i in atom_indexes]
        for atom in patch['DeleteAtoms']:
            if is_verbose:
                print('Info) Delete Atom %s' %(atom))
            index = atom_indexes[atom_name.index(atom)]
            self.atom_id = np.delete(self.atom_id, index)
            del self.atom_name[index]
            self.res_id = np.delete(self.res_id, index)
            del self.res_name[index]
            del self.chain_name[index]
            self.coord = np.delete(self.coord, [(index-1)*3, (index-1)*3+1, (index-1)*3+2]).reshape(-1, 3)
            self.mass = np.delete(self.mass, index)
            del self.raw_data[index]
            del self.line_ATOM[index]
            self.num_atoms -= 1
        if is_verbose:
            print('Info) Patch %s has been added successfully\n' %(patch_name))

    def moveBy(self, move_vec):
        move_vec = np.array(move_vec)
        for atom in range(self.num_atoms):
            self.coord[atom, :] += move_vec
            self.line_ATOM[atom] = self._replaceATOMLineCoord(self.line_ATOM[atom], self.coord[atom, :])

    def getCenterOfMass(self):
        tot_mass = 0
        center_of_mass = np.zeros(3)
        for atom in range(self.num_atoms):
            tot_mass += self.mass[atom]
            center_of_mass += self.mass[atom] * self.coord[atom, :]
        self.center_of_mass = center_of_mass/tot_mass
        return self.center_of_mass

    def setChainNameByResId(self, res_id, target_chain_name):
        indexes = [index for index, value in enumerate(self.res_id) if value == res_id]
        for index in indexes:
            self.chain_name[index] = target_chain_name
            self.line_ATOM[index] = self._replaceATOMLineChainName(self.line_ATOM[index], target_chain_name)

    def setResIdByResId(self, res_id, target_res_id):
        if target_res_id - res_id == 1 and target_res_id <= self.num_res:
            raise Warning('The distance between res_id and target_res_id is 1, which blur the distinction between neighbor residues')
        elif res_id - target_res_id == 1:
            raise Warning('The distance between res_id and target_res_id is 1, which blur the distinction between neighbor residues')
        indexes = [index for index, value in enumerate(self.res_id) if value == res_id]
        for index in indexes:
            self.res_id[index] = target_res_id
            self.line_ATOM[index] = self._replaceATOMLineResId(self.line_ATOM[index], target_res_id)

    def setAtomIdByAtomId(self, atom_id, target_atom_id):
        if target_atom_id - atom_id == 1 and target_atom_id <= self.num_atoms:
            raise Warning('The distance between atom_id and target_atom_id is 1, which blur the distinction between neighbor atoms')
        elif atom_id - target_atom_id == 1:
            raise Warning('The distance between atom_id and target_atom_id is 1, which blur the distinction between neighbor atoms')
        judge = lambda x: x==atom_id
        index = [i for i in self.atom_id if judge(i)][0]
        self.atom_id[index] = target_atom_id
        self.line_ATOM[index] = self._replaceATOMLineAtomId(self.line_ATOM[index], target_atom_id)

    def _catManipulator(self, manipulator):
        # Combine data
        self.atom_id = np.hstack((self.atom_id, manipulator.atom_id))
        self.atom_name.extend(manipulator.atom_name)
        self.res_name.extend(manipulator.res_name)
        self.chain_name.extend(manipulator.chain_name)
        self.res_id = np.hstack((self.res_id, manipulator.res_id))
        self.coord = np.vstack((self.coord, manipulator.coord))
        self.mass = np.hstack((self.mass, manipulator.mass))
        self.raw_data.extend(manipulator.raw_data)
        self.line_ATOM.extend(manipulator.line_ATOM)

        # Update Info
        self.num_atoms += manipulator.num_atoms
        self.num_res += manipulator.num_res
        self.sortAtomId()
        self.sortResId()

        if self.line_NoATOM == []:
            self.line_NoATOM.extend(manipulator.line_NoATOM)

    def catManipulators(self, *manipulators):
        for manipulator in manipulators:
            self._catManipulator(manipulator)

    def writeNewFile(self, file_name):
        with open(file_name, 'w') as io:
            for line in self.line_NoATOM:
                io.write(line)
            for line in self.line_ATOM:
                io.write(line)
            io.write(self.line_ATOMEND)
            io.close()