#!/usr/bin/env python

''' Selection of a database '''
TARGET = 'RNA09'
SUGAR_MARK = "'"
LIBPDBAA = 'RNA09_FRAG_AA/'
LISTFILE = 'RNA09.nts'

from pdbfile import PdbFile
from pdb import Chain, Residue
from mtx_coord_transform import mtx_crd_transform
from CalcROT import calcrotation
import copy
import numpy as np


""" Prepare standard base configurations from nab/[augc]_base.pdb """
p = PdbFile('nab/a_base.pdb')
p.open_to_read()
chains = p.read_all()
p.close()
res_base_A = copy.deepcopy( chains[0].residues[0] )

p = PdbFile('nab/u_base.pdb')
p.open_to_read()
chains = p.read_all()
p.close()
res_base_U = copy.deepcopy( chains[0].residues[0] )

p = PdbFile('nab/g_base.pdb')
p.open_to_read()
chains = p.read_all()
p.close()
res_base_G = copy.deepcopy( chains[0].residues[0] )

p = PdbFile('nab/c_base.pdb')
p.open_to_read()
chains = p.read_all()
p.close()
res_base_C = copy.deepcopy( chains[0].residues[0] )


def xyz_list(res, names):
    xyzs = []
    for name in names:
        xyzs.append( res.find_atom_by_name(name).xyz.get_as_list() )
    return xyzs


""" Coordinates for superpositioning """
xyzs_A = xyz_list( res_base_A, ("N9", "C8", "N7", "C6", "N6", "C5", "C4", "N3", "C2", "N1"))
xyzs_U = xyz_list( res_base_U, ("C6", "C5", "C4", "O4", "N3", "C2", "O2", "N1"))
xyzs_G = xyz_list( res_base_G, ("N9", "C8", "N7", "C6", "O6", "C5", "C4", "N3", "C2", "N2", "N1"))
xyzs_C = xyz_list( res_base_C, ("C6", "C5", "C4", "N4", "N3", "C2", "O2", "N1"))
""" Add the position of C1' """
xyzs_A.append( [7.514,4.080,7.390] )
xyzs_U.append( [7.514,4.080,7.390] )
xyzs_G.append( [7.514,4.080,7.390] )
xyzs_C.append( [7.514,4.080,7.390] )



def generate_atom(res, names, newname, chain_id, res_seq, element):
    a = False
    i = 0
    while not a:
        a = res.find_atom_by_name(names[i])
        i += 1

    if not a:
        return False
    else:
        a.name = newname
        a.chain_id = chain_id
        a.res_seq = res_seq
        a.element = element
        return a


for l in open(LISTFILE):
    lsp = l.split()
    ifrag = int(lsp[0])
    augc = lsp[7]

    pdb = PdbFile('%s/%06i.pdb' % (LIBPDBAA, ifrag))
    pdb.open_to_read()
    chains = pdb.read_all()
    pdb.close()

    r1 = chains[0].residues[0]
    r2 = chains[0].residues[1]
    r3 = chains[0].residues[2]

    c = Chain()

    ################# Backbone of Residue 1 ####################
    r = Residue()
    r.push_atom( generate_atom(r1, ("P",), " P  ", 1, 1, ' P') )
    r.push_atom( generate_atom(r1, ("O1P", "OP1" ), " OP1", 1, 1, ' O') )
    r.push_atom( generate_atom(r1, ("O2P", "OP2" ), " OP2", 1, 1, ' O') )
    r.push_atom( generate_atom(r1, ("O5'", "O5*" ), " O5'", 1, 1, ' O') )
    r.push_atom( generate_atom(r1, ("C5'", "C5*" ), " C5'", 1, 1, ' C') )
    r.push_atom( generate_atom(r1, ("C4'", "C4*" ), " C4'", 1, 1, ' C') )
    r.push_atom( generate_atom(r1, ("O4'", "O4*" ), " O4'", 1, 1, ' O') )
    r.push_atom( generate_atom(r1, ("C3'", "C3*" ), " C3'", 1, 1, ' C') )
    r.push_atom( generate_atom(r1, ("O3'", "O3*" ), " O3'", 1, 1, ' O') )
    r.push_atom( generate_atom(r1, ("C2'", "C2*" ), " C2'", 1, 1, ' C') )
    r.push_atom( generate_atom(r1, ("O2'", "O2*" ), " O2'", 1, 1, ' O') )
    r.push_atom( generate_atom(r1, ("C1'", "C1*" ), " C1'", 1, 1, ' C') )
    r.push_atom( generate_atom(r1, ("H5'", "H5*" ), "H5' ", 1, 1, ' H') )
    r.push_atom( generate_atom(r1, ("H5''","H5**"), "H5''", 1, 1, ' H') )
    r.push_atom( generate_atom(r1, ("H4'", "H4*" ), "H4' ", 1, 1, ' H') )
    r.push_atom( generate_atom(r1, ("H3'", "H3*" ), "H3' ", 1, 1, ' H') )
    r.push_atom( generate_atom(r1, ("H2'", "H2*" ), "H2' ", 1, 1, ' H') )
    r.push_atom( generate_atom(r1, ("HO2'","HO2*"), "HO2'", 1, 1, ' H') )
    r.push_atom( generate_atom(r1, ("H1'", "H1*" ), "H1' ", 1, 1, ' H') )
    c.push_residue(r)

    ################# Backbone of Residue 2 ####################
    r = Residue()
    r.push_atom( generate_atom(r2, ("P",), " P  ", 1, 2, ' P') )
    r.push_atom( generate_atom(r2, ("O1P", "OP1" ), " OP1", 1, 2, ' O') )
    r.push_atom( generate_atom(r2, ("O2P", "OP2" ), " OP2", 1, 2, ' O') )
    r.push_atom( generate_atom(r2, ("O5'", "O5*" ), " O5'", 1, 2, ' O') )
    r.push_atom( generate_atom(r2, ("C5'", "C5*" ), " C5'", 1, 2, ' C') )
    r.push_atom( generate_atom(r2, ("C4'", "C4*" ), " C4'", 1, 2, ' C') )
    r.push_atom( generate_atom(r2, ("O4'", "O4*" ), " O4'", 1, 2, ' O') )
    r.push_atom( generate_atom(r2, ("C3'", "C3*" ), " C3'", 1, 2, ' C') )
    r.push_atom( generate_atom(r2, ("O3'", "O3*" ), " O3'", 1, 2, ' O') )
    r.push_atom( generate_atom(r2, ("C2'", "C2*" ), " C2'", 1, 2, ' C') )
    r.push_atom( generate_atom(r2, ("O2'", "O2*" ), " O2'", 1, 2, ' O') )
    r.push_atom( generate_atom(r2, ("C1'", "C1*" ), " C1'", 1, 2, ' C') )
    r.push_atom( generate_atom(r2, ("H5'", "H5*" ), "H5' ", 1, 2, ' H') )
    r.push_atom( generate_atom(r2, ("H5''","H5**"), "H5''", 1, 2, ' H') )
    r.push_atom( generate_atom(r2, ("H4'", "H4*" ), "H4' ", 1, 2, ' H') )
    r.push_atom( generate_atom(r2, ("H3'", "H3*" ), "H3' ", 1, 2, ' H') )
    r.push_atom( generate_atom(r2, ("H2'", "H2*" ), "H2' ", 1, 2, ' H') )
    r.push_atom( generate_atom(r2, ("HO2'","HO2*"), "HO2'", 1, 2, ' H') )
    r.push_atom( generate_atom(r2, ("H1'", "H1*" ), "H1' ", 1, 2, ' H') )
    c.push_residue(r)

    ################# Backbone of Residue 3 ####################
    r = Residue()
    r.push_atom( generate_atom(r3, ("P",), " P  ", 1, 3, ' P') )
    r.push_atom( generate_atom(r3, ("O1P", "OP1" ), " OP1", 1, 3, ' O') )
    r.push_atom( generate_atom(r3, ("O2P", "OP2" ), " OP2", 1, 3, ' O') )
    r.push_atom( generate_atom(r3, ("O5'", "O5*" ), " O5'", 1, 3, ' O') )
    r.push_atom( generate_atom(r3, ("C5'", "C5*" ), " C5'", 1, 3, ' C') )
    r.push_atom( generate_atom(r3, ("C4'", "C4*" ), " C4'", 1, 3, ' C') )
    r.push_atom( generate_atom(r3, ("O4'", "O4*" ), " O4'", 1, 3, ' O') )
    r.push_atom( generate_atom(r3, ("C3'", "C3*" ), " C3'", 1, 3, ' C') )
    r.push_atom( generate_atom(r3, ("O3'", "O3*" ), " O3'", 1, 3, ' O') )
    r.push_atom( generate_atom(r3, ("C2'", "C2*" ), " C2'", 1, 3, ' C') )
    r.push_atom( generate_atom(r3, ("O2'", "O2*" ), " O2'", 1, 3, ' O') )
    r.push_atom( generate_atom(r3, ("C1'", "C1*" ), " C1'", 1, 3, ' C') )
    r.push_atom( generate_atom(r3, ("H5'", "H5*" ), "H5' ", 1, 3, ' H') )
    r.push_atom( generate_atom(r3, ("H5''","H5**"), "H5''", 1, 3, ' H') )
    r.push_atom( generate_atom(r3, ("H4'", "H4*" ), "H4' ", 1, 3, ' H') )
    r.push_atom( generate_atom(r3, ("H3'", "H3*" ), "H3' ", 1, 3, ' H') )
    r.push_atom( generate_atom(r3, ("H2'", "H2*" ), "H2' ", 1, 3, ' H') )
    r.push_atom( generate_atom(r3, ("HO2'","HO2*"), "HO2'", 1, 3, ' H') )
    r.push_atom( generate_atom(r3, ("H1'", "H1*" ), "H1' ", 1, 3, ' H') )
    c.push_residue(r)


    ################# Base of Residue 2 ####################

    """##################"""
    """# Super position #"""
    """##################"""

    """ Prepare list of coordinates of all haevy atoms in base and C1' as xyzs """
    """ and get the translation-rotation matrix by comparing xyzs and xyzs_[AUGC] """
    if augc == "A":
        xyzs = xyz_list( r2, ("N9", "C8", "N7", "C6", "N6", "C5", "C4", "N3", "C2", "N1", "C1%s" % SUGAR_MARK))
        rmsd, rotmat = calcrotation(np.transpose(xyzs), np.transpose(xyzs_A))
    elif augc == "U":
        xyzs = xyz_list( r2, ("C6", "C5", "C4", "O4", "N3", "C2", "O2", "N1", "C1%s" % SUGAR_MARK))
        rmsd, rotmat = calcrotation(np.transpose(xyzs), np.transpose(xyzs_U))
    elif augc == "G":
        xyzs = xyz_list( r2, ("N9", "C8", "N7", "C6", "O6", "C5", "C4", "N3", "C2", "N2", "N1", "C1%s" % SUGAR_MARK))
        rmsd, rotmat = calcrotation(np.transpose(xyzs), np.transpose(xyzs_G))
    elif augc == "C":
        xyzs = xyz_list( r2, ("C6", "C5", "C4", "N4", "N3", "C2", "O2", "N1", "C1%s" % SUGAR_MARK))
        rmsd, rotmat = calcrotation(np.transpose(xyzs), np.transpose(xyzs_C))

    mtx = mtx_crd_transform()
    mtx.mtx[:,:] = rotmat


    """ Generate four types (AUGC) of fragment """
    for char, res_target in (('A', res_base_A), 
                             ('U', res_base_U),
                             ('G', res_base_G),
                             ('C', res_base_C),):

        c_out = copy.deepcopy(c)

        """ Apply the translation-rotation matrix """
        for a in res_target.atoms:
            acopy = copy.deepcopy(a)
            c_out.residues[1].push_atom( acopy )
            c_out.residues[1].atoms[-1].xyz.put_as_list( mtx.do_to_array( acopy.xyz.get_as_list() ) )

        """ Re-index the serial ID and chain ID etc. """
        iserial = 0
        chain_id = c_out.residues[0].atoms[0].chain_id
        for ir, r in enumerate(c_out.residues):
            res_name = r.atoms[0].res_name
            for a in r.atoms:
                iserial += 1
                a.chain_id = chain_id
                a.serial = iserial
                a.res_seq = ir + 1
                if ir == 1:
                    a.res_name = '  %s' % char
                else:
                    a.res_name = res_name
                a.ins_code = ''
                a.alt_loc = ''
    
        c_out.reindex()


        """########################"""
        """# Output the fragments #"""
        """########################"""

        pdb = PdbFile('%s/%06i_%s.pdb' % (LIBPDBAA, ifrag, char))
        pdb.open_to_write()
        pdb.write_all([c_out,])
        pdb.close()
