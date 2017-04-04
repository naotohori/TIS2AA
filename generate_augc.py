#!/usr/bin/env python

TARGET = 'RNA09'
SUGAR_MARK = "'"
LIBPDBAA = 'RNA09_FRAG_AA/'
LISTFILE = 'RNA09.nts'

from cafysis.file_io.pdb import PdbFile
from cafysis.elements.pdb import Chain, Residue
from cafysis.mtx_coord_transform import mtx_crd_transform
from CalcROT import calcrotation
import copy
import numpy as np

## Coordinates of C1', N9/N1, C4/C2 of nab/[augc].pdb
base3_A = [[7.514,4.080,7.390],[6.275,3.365,7.010],[5.089,3.947,6.630]]
base3_G = [[7.514,4.080,7.390],[6.275,3.365,7.010],[5.082,3.956,6.630]]
base3_U = [[7.514,4.080,7.390],[6.275,3.365,7.010],[5.189,4.127,6.660]]
base3_C = [[7.514,4.080,7.390],[6.275,3.365,7.010],[5.166,4.124,6.660]]


## Prepare standard base configurations from nab/[augc]_base.pdb
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


def generate_atom(res, names, newname, chain_id, res_seq):
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
    r.push_atom( generate_atom(r1, ("P",), " P  ", 1, 1) )
    r.push_atom( generate_atom(r1, ("O1P", "OP1" ), " OP1", 1, 1) )
    r.push_atom( generate_atom(r1, ("O2P", "OP2" ), " OP2", 1, 1) )
    r.push_atom( generate_atom(r1, ("O5'", "O5*" ), " O5'", 1, 1) )
    r.push_atom( generate_atom(r1, ("C5'", "C5*" ), " C5'", 1, 1) )
    r.push_atom( generate_atom(r1, ("C4'", "C4*" ), " C4'", 1, 1) )
    r.push_atom( generate_atom(r1, ("O4'", "O4*" ), " O4'", 1, 1) )
    r.push_atom( generate_atom(r1, ("C3'", "C3*" ), " C3'", 1, 1) )
    r.push_atom( generate_atom(r1, ("O3'", "O3*" ), " O3'", 1, 1) )
    r.push_atom( generate_atom(r1, ("C2'", "C2*" ), " C2'", 1, 1) )
    r.push_atom( generate_atom(r1, ("O2'", "O2*" ), " O2'", 1, 1) )
    r.push_atom( generate_atom(r1, ("C1'", "C1*" ), " C1'", 1, 1) )
    r.push_atom( generate_atom(r1, ("H5'", "H5*" ), "H5' ", 1, 1) )
    r.push_atom( generate_atom(r1, ("H5''","H5**"), "H5''", 1, 1) )
    r.push_atom( generate_atom(r1, ("H4'", "H4*" ), "H4' ", 1, 1) )
    r.push_atom( generate_atom(r1, ("H3'", "H3*" ), "H3' ", 1, 1) )
    r.push_atom( generate_atom(r1, ("H2'", "H2*" ), "H2' ", 1, 1) )
    r.push_atom( generate_atom(r1, ("HO2'","HO2*"), "HO2'", 1, 1) )
    r.push_atom( generate_atom(r1, ("H1'", "H1*" ), "H1' ", 1, 1) )
    c.push_residue(r)

    ################# Backbone of Residue 2 ####################
    r = Residue()
    r.push_atom( generate_atom(r2, ("P",), " P  ", 1, 2) )
    r.push_atom( generate_atom(r2, ("O1P", "OP1" ), " OP1", 1, 2) )
    r.push_atom( generate_atom(r2, ("O2P", "OP2" ), " OP2", 1, 2) )
    r.push_atom( generate_atom(r2, ("O5'", "O5*" ), " O5'", 1, 2) )
    r.push_atom( generate_atom(r2, ("C5'", "C5*" ), " C5'", 1, 2) )
    r.push_atom( generate_atom(r2, ("C4'", "C4*" ), " C4'", 1, 2) )
    r.push_atom( generate_atom(r2, ("O4'", "O4*" ), " O4'", 1, 2) )
    r.push_atom( generate_atom(r2, ("C3'", "C3*" ), " C3'", 1, 2) )
    r.push_atom( generate_atom(r2, ("O3'", "O3*" ), " O3'", 1, 2) )
    r.push_atom( generate_atom(r2, ("C2'", "C2*" ), " C2'", 1, 2) )
    r.push_atom( generate_atom(r2, ("O2'", "O2*" ), " O2'", 1, 2) )
    r.push_atom( generate_atom(r2, ("C1'", "C1*" ), " C1'", 1, 2) )
    r.push_atom( generate_atom(r2, ("H5'", "H5*" ), "H5' ", 1, 2) )
    r.push_atom( generate_atom(r2, ("H5''","H5**"), "H5''", 1, 2) )
    r.push_atom( generate_atom(r2, ("H4'", "H4*" ), "H4' ", 1, 2) )
    r.push_atom( generate_atom(r2, ("H3'", "H3*" ), "H3' ", 1, 2) )
    r.push_atom( generate_atom(r2, ("H2'", "H2*" ), "H2' ", 1, 2) )
    r.push_atom( generate_atom(r2, ("HO2'","HO2*"), "HO2'", 1, 2) )
    r.push_atom( generate_atom(r2, ("H1'", "H1*" ), "H1' ", 1, 2) )
    c.push_residue(r)

    ################# Backbone of Residue 3 ####################
    r = Residue()
    r.push_atom( generate_atom(r3, ("P",), " P  ", 1, 3) )
    r.push_atom( generate_atom(r3, ("O1P", "OP1" ), " OP1", 1, 3) )
    r.push_atom( generate_atom(r3, ("O2P", "OP2" ), " OP2", 1, 3) )
    r.push_atom( generate_atom(r3, ("O5'", "O5*" ), " O5'", 1, 3) )
    r.push_atom( generate_atom(r3, ("C5'", "C5*" ), " C5'", 1, 3) )
    r.push_atom( generate_atom(r3, ("C4'", "C4*" ), " C4'", 1, 3) )
    r.push_atom( generate_atom(r3, ("O4'", "O4*" ), " O4'", 1, 3) )
    r.push_atom( generate_atom(r3, ("C3'", "C3*" ), " C3'", 1, 3) )
    r.push_atom( generate_atom(r3, ("O3'", "O3*" ), " O3'", 1, 3) )
    r.push_atom( generate_atom(r3, ("C2'", "C2*" ), " C2'", 1, 3) )
    r.push_atom( generate_atom(r3, ("O2'", "O2*" ), " O2'", 1, 3) )
    r.push_atom( generate_atom(r3, ("C1'", "C1*" ), " C1'", 1, 3) )
    r.push_atom( generate_atom(r3, ("H5'", "H5*" ), "H5' ", 1, 3) )
    r.push_atom( generate_atom(r3, ("H5''","H5**"), "H5''", 1, 3) )
    r.push_atom( generate_atom(r3, ("H4'", "H4*" ), "H4' ", 1, 3) )
    r.push_atom( generate_atom(r3, ("H3'", "H3*" ), "H3' ", 1, 3) )
    r.push_atom( generate_atom(r3, ("H2'", "H2*" ), "H2' ", 1, 3) )
    r.push_atom( generate_atom(r3, ("HO2'","HO2*"), "HO2'", 1, 3) )
    r.push_atom( generate_atom(r3, ("H1'", "H1*" ), "H1' ", 1, 3) )
    c.push_residue(r)


    ################# Base of Residue 2 ####################

    c_A = copy.deepcopy(c)
    c_U = copy.deepcopy(c)
    c_G = copy.deepcopy(c)
    c_C = copy.deepcopy(c)

    base3 = []
    # C1' of the second residue
    base3.append( c.residues[1].atoms[-1].xyz.get_as_list() )
    
    if augc in ("A", "G"):
        base3.append( r2.find_atom_by_name("N9").xyz.get_as_list() )
        base3.append( r2.find_atom_by_name("C4").xyz.get_as_list() )

    elif augc in ("C", "U"):
        base3.append( r2.find_atom_by_name("N1").xyz.get_as_list() )
        base3.append( r2.find_atom_by_name("C2").xyz.get_as_list() )

    if augc == "A":
        rmsd, rotmat = calcrotation(np.transpose(base3), np.transpose(base3_A))
    elif augc == "U":
        rmsd, rotmat = calcrotation(np.transpose(base3), np.transpose(base3_U))
    elif augc == "G":
        rmsd, rotmat = calcrotation(np.transpose(base3), np.transpose(base3_G))
    elif augc == "C":
        rmsd, rotmat = calcrotation(np.transpose(base3), np.transpose(base3_C))

    mtx = mtx_crd_transform()
    mtx.mtx[:,:] = rotmat


    for a in res_base_A.atoms:
        acopy = copy.deepcopy(a)
        c_A.residues[1].push_atom( acopy )
        c_A.residues[1].atoms[-1].xyz.put_as_list( mtx.do_to_array( acopy.xyz.get_as_list() ) )

    iserial = 0
    chain_id = c_A.residues[0].atoms[0].chain_id
    for r in c_A.residues:
        for a in r.atoms:
            iserial += 1
            a.chain_id = chain_id
            a.serial = iserial

    c_A.reindex()

    pdb = PdbFile('%s/%06i_A.pdb' % (LIBPDBAA, ifrag))
    pdb.open_to_write()
    pdb.write_all([c_A,])
    pdb.close()

    for a in res_base_U.atoms:
        acopy = copy.deepcopy(a)
        c_U.residues[1].push_atom( acopy )
        c_U.residues[1].atoms[-1].xyz.put_as_list( mtx.do_to_array( acopy.xyz.get_as_list() ) )

    iserial = 0
    chain_id = c_U.residues[0].atoms[0].chain_id
    for r in c_U.residues:
        for a in r.atoms:
            iserial += 1
            a.chain_id = chain_id
            a.serial = iserial

    c_U.reindex()

    pdb = PdbFile('%s/%06i_U.pdb' % (LIBPDBAA, ifrag))
    pdb.open_to_write()
    pdb.write_all([c_U,])
    pdb.close()

    for a in res_base_G.atoms:
        acopy = copy.deepcopy(a)
        c_G.residues[1].push_atom( acopy )
        c_G.residues[1].atoms[-1].xyz.put_as_list( mtx.do_to_array( acopy.xyz.get_as_list() ) )

    iserial = 0
    chain_id = c_G.residues[0].atoms[0].chain_id
    for r in c_G.residues:
        for a in r.atoms:
            iserial += 1
            a.chain_id = chain_id
            a.serial = iserial

    c_G.reindex()

    pdb = PdbFile('%s/%06i_G.pdb' % (LIBPDBAA, ifrag))
    pdb.open_to_write()
    pdb.write_all([c_G,])
    pdb.close()

    for a in res_base_C.atoms:
        acopy = copy.deepcopy(a)
        c_C.residues[1].push_atom( acopy )
        c_C.residues[1].atoms[-1].xyz.put_as_list( mtx.do_to_array( acopy.xyz.get_as_list() ) )

    iserial = 0
    chain_id = c_C.residues[0].atoms[0].chain_id
    for r in c_C.residues:
        for a in r.atoms:
            iserial += 1
            a.chain_id = chain_id
            a.serial = iserial

    c_C.reindex()

    pdb = PdbFile('%s/%06i_C.pdb' % (LIBPDBAA, ifrag))
    pdb.open_to_write()
    pdb.write_all([c_C,])
    pdb.close()

