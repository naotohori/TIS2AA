#!/usr/bin/env python

SUGAR_MARK = "'"
from pdbfile import PdbFile
from pdb import Chain, Residue
from mtx_coord_transform import mtx_crd_transform
from CalcROT import calcrotation
import copy
import numpy as np
import sys


''' INPUT '''
if len(sys.argv) != 5:
    print 'Usage: SCRIPT [PDB file] [5 or 3] [nucleotide type (a,u,g,c)] [output PDB]'
    sys.exit(2)

filename_pdb = sys.argv[1]
end = int(sys.argv[2])
nuc = sys.argv[3].lower()
filename_out = sys.argv[4]

''' Check arguments '''
if end not in (5,3):
    print 'Error: the second argument has to be either 5 or 3'
    sys.exit(2)
if nuc not in ('a','u','g','c'):
    print 'Error: the last argument has to be either one of augc'
    sys.exit(2)


''' Read the PDB '''
p = PdbFile(filename_pdb)
p.open_to_read()
chains = p.read_all()
p.close()

if len(chains) > 1:
    print 'Error: currently available only for single chain PDB'
    sys.exit(2)
    # To do: add multiple chain case. Maybe specify chain ID as an argument


if end == 5:
    nuc_ref = chains[0].residues[0].atoms[0].res_name.lower().strip()
    if nuc_ref[-1] == '5':
        nuc_ref = nuc_ref[0]
        for a in chains[0].residues[0].atoms:
            a.res_name = '   %s' % nuc_ref.upper()
    dimer_seq = '%s%s' % (nuc, nuc_ref)
else:
    nuc_ref = chains[0].residues[-1].atoms[0].res_name.lower().strip()
    if nuc_ref[-1] == '3':
        nuc_ref = nuc_ref[0]
        for a in chains[0].residues[-1].atoms:
            a.res_name = '   %s' % nuc_ref.upper()
    dimer_seq = '%s%s' % (nuc_ref, nuc)

print nuc_ref, dimer_seq


""" Prepare standard dimerconfigurations from nab """
p = PdbFile('nab/%s.pdb' % dimer_seq)
p.open_to_read()
dimer = p.read_all()[0]
p.close()

def xyz_list(res, names):
    xyzs = []
    for name in names:
        xyzs.append( res.find_atom_by_name(name).xyz.get_as_list() )
    return xyzs


if end == 5:
    dimer_r = dimer.residues[1]
else:
    dimer_r = dimer.residues[0]

if nuc_ref == 'a':
    xyz_dimer = xyz_list( dimer_r, ("O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
                                    "N9", "C8", "N7", "C6", "N6", "C5", "C4", "N3", "C2", "N1") )
elif nuc_ref == 'u':
    xyz_dimer = xyz_list( dimer_r, ("O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
                                    "C6", "C5", "C4", "O4", "N3", "C2", "O2", "N1") )
elif nuc_ref == 'g':
    xyz_dimer = xyz_list( dimer_r, ("O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
                                    "N9", "C8", "N7", "C6", "O6", "C5", "C4", "N3", "C2", "N2", "N1") )
elif nuc_ref == 'c':
    xyz_dimer = xyz_list( dimer_r, ("O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
                                    "C6", "C5", "C4", "N4", "N3", "C2", "O2", "N1"))


if end == 5:
    target_r = chains[0].residues[0]
else:
    target_r = chains[0].residues[-1]

if nuc_ref == 'a':
    xyz_target = xyz_list( target_r, ("O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
                                      "N9", "C8", "N7", "C6", "N6", "C5", "C4", "N3", "C2", "N1") )
elif nuc_ref == 'u':
    xyz_target = xyz_list( target_r, ("O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
                                      "C6", "C5", "C4", "O4", "N3", "C2", "O2", "N1") )
elif nuc_ref == 'g':
    xyz_target = xyz_list( target_r, ("O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
                                      "N9", "C8", "N7", "C6", "O6", "C5", "C4", "N3", "C2", "N2", "N1") )
elif nuc_ref == 'c':
    xyz_target = xyz_list( target_r, ("O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
                                      "C6", "C5", "C4", "N4", "N3", "C2", "O2", "N1"))


rmsd, rotmat = calcrotation(np.transpose(xyz_target), np.transpose(xyz_dimer))
print 'rmsd=',rmsd
mtx = mtx_crd_transform()
mtx.mtx[:,:] = rotmat


c_out = copy.deepcopy(chains[0])

if end == 5:
    chain_id = c_out.residues[0].atoms[0].chain_id
    res_seq = c_out.residues[0].atoms[0].res_seq

    r = Residue()
    for a in dimer.residues[0].atoms:
        acopy = copy.deepcopy(a)
        acopy.chain_id = chain_id
        acopy.res_seq = res_seq - 1
        if acopy.res_name[-1] == '5':
            acopy.res_name = '  %s' % acopy.res_name.strip()[0]
        r.push_atom( acopy )
        r.atoms[-1].xyz.put_as_list( mtx.do_to_array( acopy.xyz.get_as_list() ) )

    c_out.residues.insert(0,r)

    ''' Delete HO5' '''
    for a in c_out.residues[1].atoms:
        if a.name.strip() == 'HO5%s' % SUGAR_MARK:
            c_out.residues[1].atoms.remove(a)
            break

    ''' Push P, OP1, OP2 from the second nucoletide of the dimer '''
    for a in dimer.residues[1].atoms:
        acopy = copy.deepcopy(a)
        acopy.chain_id = chain_id
        acopy.res_seq = res_seq
        if acopy.res_name[-1] == '3':
            acopy.res_name = '  %s' % acopy.res_name.strip()[0]
        if a.name.strip() == 'OP2':
            c_out.residues[1].atoms.insert(0,acopy)
            c_out.residues[1].atoms[0].xyz.put_as_list( mtx.do_to_array( acopy.xyz.get_as_list() ) )
    for a in dimer.residues[1].atoms:
        acopy = copy.deepcopy(a)
        acopy.chain_id = chain_id
        acopy.res_seq = res_seq
        if acopy.res_name[-1] == '3':
            acopy.res_name = '  %s' % acopy.res_name.strip()[0]
        if a.name.strip() == 'OP1':
            c_out.residues[1].atoms.insert(0,acopy)
            c_out.residues[1].atoms[0].xyz.put_as_list( mtx.do_to_array( acopy.xyz.get_as_list() ) )
    for a in dimer.residues[1].atoms:
        acopy = copy.deepcopy(a)
        acopy.chain_id = chain_id
        acopy.res_seq = res_seq
        if acopy.res_name[-1] == '3':
            acopy.res_name = '  %s' % acopy.res_name.strip()[0]
        if a.name.strip() == 'P':
            c_out.residues[1].atoms.insert(0,acopy)
            c_out.residues[1].atoms[0].xyz.put_as_list( mtx.do_to_array( acopy.xyz.get_as_list() ) )


if end == 3:
    ''' Delete HO3' '''
    for a in c_out.residues[-1].atoms:
        if a.name.strip() == 'HO3%s' % SUGAR_MARK:
            c_out.residues[-1].atoms.remove(a)
            break

    chain_id = c_out.residues[-1].atoms[-1].chain_id
    res_seq = c_out.residues[-1].atoms[-1].res_seq

    r = Residue()
    for a in dimer.residues[1].atoms:
        acopy = copy.deepcopy(a)
        acopy.chain_id = chain_id
        acopy.res_seq = res_seq + 1
        if acopy.res_name[-1] == '3':
            acopy.res_name = '  %s' % acopy.res_name.strip()[0]
        r.push_atom( acopy )
        r.atoms[-1].xyz.put_as_list( mtx.do_to_array( acopy.xyz.get_as_list() ) )

    c_out.residues.append(r)

''' Asign serial numbers '''
serial = 0
for r in c_out.residues:
    for a in r.atoms:
        serial += 1
        a.serial = serial

p = PdbFile(filename_out)
p.open_to_write()
p.write_all([c_out,])
p.close()
