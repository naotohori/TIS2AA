#!/usr/bin/env python


import os
import glob

from coord import Coord
from pdb_elements import Chain, Residue, Atom
from pdbfile import PdbFile

''' Selection of a database '''
SUGAR_MARK = "'"
BASE_DIR = os.path.dirname(os.path.realpath(__file__)) + '/../'

mass = {'P':30.973761, 'O':15.9994, 'C':12.0107, 'N':14.0067}

cgmodel = "TISRNA"
ATOMS_P = ('P', 'OP1', 'OP2')
ATOMS_S = ('')

#cgmodel = "TISDNA"
#ATOMS_P = ('P', 'OP1', 'OP2', "O5'")

cg_residues = []

def aa2tis(aa_residues):


    for r in residues:



for path in sorted(glob.glob('%s/*_?.pdb' % LIBPDBAA)):

    pdb_in = PdbFile(path)
    pdb_in.open_to_read()
    c = pdb_in.read_all()[0]
    pdb_in.close()

    c_out = Chain()

    res_id = 0
    atom_id = 0

    for ir, r in enumerate(c.residues):
        xyz_P = Coord()
        nP = 0
        if cgmodel == 'TISDNA' and ir > 0:
            xyz_P += xyz_nextP
            nP += 1
        xyz_nextP = Coord()
        nnextP = 0
        xyz_S = Coord()
        nS = 0
        xyz_B = Coord()
        nB = 0

        for a in r.atoms:
            name = a.name.strip()

            if name[0] == 'H':
                continue

            elif name in ATOMS_P:
                xyz_P += a.xyz
                nP += 1

            elif cgmodel == 'TISDNA' and name == 'O3'+SUGAR_MARK:
                xyz_nextP += a.xyz
                nnextP += 1

            elif name.find(SUGAR_MARK) != -1:
                xyz_S += a.xyz
                nS += 1

            else:
                xyz_B += a.xyz
                nB += 1

            nt = a.res_name.strip()

        res_id += 1
        r_cg = Residue()

        # if nP != 0:
        if ir > 0:
            atom_id += 1
            a = Atom()
            a.serial = atom_id
            a.name = ' P  '
            if cgmodel == 'TISRNA':
                a.res_name = 'R%s ' % nt
            else:
                a.res_name = '%s ' % nt
            a.chain_id = 'A'
            a.res_seq = res_id
            a.xyz = xyz_P / float(nP)
            r_cg.push_atom(a)

        if nS != 0:
            atom_id += 1
            a = Atom()
            a.serial = atom_id
            a.name = ' S  '
            if cgmodel == 'TISRNA':
                a.res_name = 'R%s ' % nt
            else:
                a.res_name = '%s ' % nt
            a.chain_id = 'A'
            a.res_seq = res_id
            a.xyz = xyz_S / float(nS)
            r_cg.push_atom(a)

        if nB != 0:
            atom_id += 1
            a = Atom()
            a.serial = atom_id
            a.name = ' B  '
            if cgmodel == 'TISRNA':
                a.res_name = 'R%s ' % nt
            else:
                a.res_name = '%s ' % nt
            a.chain_id = 'A'
            a.res_seq = res_id
            a.xyz = xyz_B / float(nB)
            r_cg.push_atom(a)

        c_out.push_residue(r_cg)

    pdb_out = PdbFile(LIBPDBCG + os.path.basename(path))
    pdb_out.open_to_write()
    pdb_out.write_all([c_out, ])
    pdb_out.close()
