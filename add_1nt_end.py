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
    nuc_end = chains[0].residues[0].atoms[0].res_name.lower().strip()
    if nuc_end[-1] == '5':
        nuc_end = nuc_end[0]
        for a in chains[0].residues[0].atoms:
            a.res_name = '   %s' % nuc_end.upper()
    dimer_seq = '%s%s' % (nuc, nuc_end)

    print 'The first nucleotide of the PDB is ', nuc_end.upper(), '.'
    print 'Now ', nuc.upper(),' is being added by using a dimer', dimer_seq.upper(),'.'

else: # end == 3
    nuc_end = chains[0].residues[-1].atoms[0].res_name.lower().strip()
    if nuc_end[-1] == '3':
        nuc_end = nuc_end[0]
        for a in chains[0].residues[-1].atoms:
            a.res_name = '   %s' % nuc_end.upper()
    dimer_seq = '%s%s' % (nuc_end, nuc)

    print 'The last nucleotide of the PDB is ', nuc_end.upper(), '.'
    print 'Now ', nuc.upper(),' is being added by using a dimer', dimer_seq.upper(),'.'


""" Prepare standard dimer configurations from nab """
p = PdbFile('nab/%s.pdb' % dimer_seq)
p.open_to_read()
dimer = p.read_all()[0]
p.close()


""" a function to extract coordinates by atom names """
def xyz_list(res, names):
    xyzs = []
    for name in names:
        xyzs.append( res.find_atom_by_name(name).xyz.get_as_list() )
    return xyzs


""" Residues for superposition """
if end == 5:
    dimer_res = dimer.residues[1]  # Superposition is done by the 2nd nucleotide of the dimer
    target_res = chains[0].residues[0]
else: # end == 3
    dimer_res = dimer.residues[0]  # Superposition is done by the 1st nucleotide of the dimer '''
    target_res = chains[0].residues[-1]


""" Preparae coordinates from the dimer for superposigion """
if nuc_end == 'a':
    xyz_dimer = xyz_list( dimer_res, ("O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
                                    "N9", "C8", "N7", "C6", "N6", "C5", "C4", "N3", "C2", "N1") )
elif nuc_end == 'u':
    xyz_dimer = xyz_list( dimer_res, ("O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
                                    "C6", "C5", "C4", "O4", "N3", "C2", "O2", "N1") )
elif nuc_end == 'g':
    xyz_dimer = xyz_list( dimer_res, ("O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
                                    "N9", "C8", "N7", "C6", "O6", "C5", "C4", "N3", "C2", "N2", "N1") )
elif nuc_end == 'c':
    xyz_dimer = xyz_list( dimer_res, ("O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
                                    "C6", "C5", "C4", "N4", "N3", "C2", "O2", "N1"))


""" Preparae coordinates from the target PDB for superposigion """
if nuc_end == 'a':
    xyz_target = xyz_list( target_res, ("O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
                                      "N9", "C8", "N7", "C6", "N6", "C5", "C4", "N3", "C2", "N1") )
elif nuc_end == 'u':
    xyz_target = xyz_list( target_res, ("O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
                                      "C6", "C5", "C4", "O4", "N3", "C2", "O2", "N1") )
elif nuc_end == 'g':
    xyz_target = xyz_list( target_res, ("O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
                                      "N9", "C8", "N7", "C6", "O6", "C5", "C4", "N3", "C2", "N2", "N1") )
elif nuc_end == 'c':
    xyz_target = xyz_list( target_res, ("O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'",
                                      "C6", "C5", "C4", "N4", "N3", "C2", "O2", "N1"))


""" Superposition """
rmsd, rotmat = calcrotation(np.transpose(xyz_target), np.transpose(xyz_dimer))
print 'rmsd=',rmsd
mtx = mtx_crd_transform()
mtx.mtx[:,:] = rotmat


""" Generate new PDB """
c_out = copy.deepcopy(chains[0])

if end == 5:
    """ New 1st nucleotide """
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

    ''' Delete HO5' from the new 2nd nucleotide '''
    for a in c_out.residues[1].atoms:
        if a.name.strip() == 'HO5%s' % SUGAR_MARK:
            c_out.residues[1].atoms.remove(a)
            break

    ''' Push P, OP1, OP2 from the 2nd nucoletide of the dimer into the new 2nd nucleotide'''
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
    ''' Delete HO3' form the last nucleotide (it will be the 2nd last nucleotide) '''
    for a in c_out.residues[-1].atoms:
        if a.name.strip() == 'HO3%s' % SUGAR_MARK:
            c_out.residues[-1].atoms.remove(a)
            break

    ''' Push new last nucleotide '''
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


''' Re-asign serial numbers '''
serial = 0
for r in c_out.residues:
    for a in r.atoms:
        serial += 1
        a.serial = serial


''' Write new PDB file '''
p = PdbFile(filename_out)
p.open_to_write()
p.write_all([c_out,])
p.close()
