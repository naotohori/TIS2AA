#!/usr/bin/env python
'''
@author: Naoto Hori
'''

from coord import Coord

class Atom :
    def __init__(self) :
        self.serial = ''
        self.name = ''
        self.alt_loc = ''
        self.res_name = ''
        self.chain_id = ''
        self.res_seq = 0
        self.ins_code = ''
        self.xyz = Coord()
        self.occupancy = 0.0
        self.temp_factor = 0.0
        self.element = ''
        self.charge = 0
        
    def show(self):
        print 'serial:', self.serial
        print 'name:', self.name
        print 'alt_loc:', self.alt_loc
        print 'res_name:', self.res_name
        print 'chain_id:', self.chain_id
        print 'res_seq:', self.res_seq
        print 'ins_code', self.ins_code
        print 'xyz', self.xyz.x, self.xyz.y, self.xyz.z
        print 'occupancy', self.occupancy
        print 'temp_factor', self.temp_factor
        print 'element', self.element
        print 'charge', self.charge
        
class Residue :
    def __init__(self):
        self.atoms = []
    
    def push_atom(self, a):
        self.atoms.append(a)

    def find_Calpha_atom(self):
        for a in self.atoms:
            if a.name == ' CA ':
                return a
        return False

    def find_atom_by_name(self, name):
        for a in self.atoms:
            if a.name.strip() == name:
                return a
        return False
        
class Chain :
    def __init__(self):
        self.residues = []
        self._where_is_atomX = []
        self.num_atom = lambda : len(self._where_is_atomX)
        self.num_res = lambda : len(self.residues)
        
    def push_residue(self, r):
        residue_id = len(self.residues)
        num_atom = len(r.atoms)
        for i in xrange(num_atom):
            self._where_is_atomX.append((residue_id, i))
        self.residues.append(r)
        
    def print_all(self):
        for residue in self.residues :
            for atom in residue.atoms :
                atom.show()
                
    def get_atom(self, idx):
        (residue_id, atom_local_id) = self._where_is_atomX[idx]
        #print (residue_id, atom_local_id)
        return (self.residues[residue_id].atoms[atom_local_id])

    def reindex(self):
        self._where_is_atomX = []
        for ir, r in enumerate(self.residues):
            for ia, a in enumerate(r.atoms):
                self._where_is_atomX.append((ir,ia))

#class Model :
#    def __init__(self):
#        self.chains = []
#        
#    def push_chain(self, c):
#        self.chains.append(c)
#
