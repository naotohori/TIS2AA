#!/usr/bin/env python

DIR_NAB = "/home/hori/TIS2AA/nab"
SUGAR_MARK = "'"
from pdbfile import PdbFile
from pdb_elements import Chain, Residue
from mtx_coord_transform import mtx_crd_transform
from CalcROT import calcrotation
import copy
import numpy as np
import sys
import argparse

''' INPUT '''
parser = argparse.ArgumentParser(
    description='Extend terminal(s) by one nucleotide',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--5end', dest='nuc5', default=None,
                    action='store',
                    help='Nucleotide (a,u,g,c) to be added to the 5-end')

parser.add_argument('--3end', dest='nuc3', default=None,
                    action='store',
                    help='Nucleotide (a,u,g,c) to be added to the 3-end')

parser.add_argument('pdb_in',  help='Target PDB file')
parser.add_argument('pdb_out', help='Output PDB file')

args = parser.parse_args()

''' Check arguments '''
flg_5 = False
flg_3 = False

if args.nuc5 is not None:
    flg_5 = True
    if args.nuc5 not in ('a', 'u', 'g', 'c', 'A', 'U', 'G', 'C'):
        print ('Error: --end5 has to be either one of augc')
        sys.exit(2)

if args.nuc3 is not None:
    flg_3 = True
    if args.nuc3 not in ('a', 'u', 'g', 'c', 'A', 'U', 'G', 'C'):
        print ('Error: --end3 has to be either one of augc')
        sys.exit(2)

if not flg_5 and not flg_3:
    print ('Error: either or both of --3end or --5end option is needed')
    sys.exit(2)

''' Read the PDB '''
p = PdbFile(args.pdb_in)
p.open_to_read()
chains = p.read_all()
p.close()

if len(chains) > 1:
    print ('Error: currently available only for single chain PDB')
    sys.exit(2)
    # To do: add multiple chain case. Maybe specify chain ID as an argument


""" a function to extract coordinates by atom names """
def xyz_list(res, names):
    xyzs = []
    for name in names:
        xyzs.append( res.find_atom_by_name(name).xyz.get_as_list() )
    return xyzs


""" function to derive superposition matrix"""
def dimer_superposition(dimer_res, target_res, nuc_end):

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
    print ('rmsd=', rmsd)
    mtx = mtx_crd_transform()
    mtx.mtx[:,:] = rotmat

    return mtx


""" Generate new PDB """
c_out = copy.deepcopy(chains[0])

if flg_5:

    nuc_end5 = chains[0].residues[0].atoms[0].res_name.lower().strip()
    if nuc_end5[0] in ('r', ):
        nuc_end5 = nuc_end5[1:]
    if nuc_end5[-1] == '5':
        nuc_end5 = nuc_end5[0]
        for a in c_out.residues[0].atoms:
            a.res_name = '  %s' % nuc_end5.upper()
    dimer_seq5 = '%s%s' % (args.nuc5.lower(), nuc_end5)

    print ('The first nucleotide of the PDB is ', nuc_end5.upper(), '.')
    print ('Now ', args.nuc5.upper(),' is being added by using a dimer', dimer_seq5.upper(),'.')

    """ Prepare standard dimer configurations from nab """
    p = PdbFile('%s/%s.pdb' % (DIR_NAB, dimer_seq5))
    p.open_to_read()
    dimer5 = p.read_all()[0]
    p.close()

    """ Residues for superposition """
    dimer_res = dimer5.residues[1]  # Superposition is done by the 2nd nucleotide of the dimer
    target_res = chains[0].residues[0]

    mtx5 = dimer_superposition(dimer_res, target_res, nuc_end5)

    """ New 1st nucleotide """
    chain_id = c_out.residues[0].atoms[0].chain_id
    res_seq = c_out.residues[0].atoms[0].res_seq

    r = Residue()
    for a in dimer5.residues[0].atoms:
        acopy = copy.deepcopy(a)
        acopy.chain_id = chain_id
        acopy.res_seq = res_seq - 1
        if acopy.res_name[-1] == '5':
            acopy.res_name = '  %s' % acopy.res_name.strip()[0]
        r.push_atom( acopy )
        r.atoms[-1].xyz.put_as_list( mtx5.do_to_array( acopy.xyz.get_as_list() ) )

    c_out.residues.insert(0,r)

    ''' Delete HO5' from the new 2nd nucleotide '''
    for a in c_out.residues[1].atoms:
        if a.name.strip() == 'HO5%s' % SUGAR_MARK:
            c_out.residues[1].atoms.remove(a)
            break

    ''' Push P, OP1, OP2 from the 2nd nucoletide of the dimer into the new 2nd nucleotide'''
    for a in dimer5.residues[1].atoms:
        acopy = copy.deepcopy(a)
        acopy.chain_id = chain_id
        acopy.res_seq = res_seq
        if acopy.res_name[-1] == '3':
            acopy.res_name = '  %s' % acopy.res_name.strip()[0]
        if a.name.strip() == 'OP2':
            c_out.residues[1].atoms.insert(0,acopy)
            c_out.residues[1].atoms[0].xyz.put_as_list( mtx5.do_to_array( acopy.xyz.get_as_list() ) )
    for a in dimer5.residues[1].atoms:
        acopy = copy.deepcopy(a)
        acopy.chain_id = chain_id
        acopy.res_seq = res_seq
        if acopy.res_name[-1] == '3':
            acopy.res_name = '  %s' % acopy.res_name.strip()[0]
        if a.name.strip() == 'OP1':
            c_out.residues[1].atoms.insert(0,acopy)
            c_out.residues[1].atoms[0].xyz.put_as_list( mtx5.do_to_array( acopy.xyz.get_as_list() ) )
    for a in dimer5.residues[1].atoms:
        acopy = copy.deepcopy(a)
        acopy.chain_id = chain_id
        acopy.res_seq = res_seq
        if acopy.res_name[-1] == '3':
            acopy.res_name = '  %s' % acopy.res_name.strip()[0]
        if a.name.strip() == 'P':
            c_out.residues[1].atoms.insert(0,acopy)
            c_out.residues[1].atoms[0].xyz.put_as_list( mtx5.do_to_array( acopy.xyz.get_as_list() ) )

if flg_3:
    nuc_end3 = chains[0].residues[-1].atoms[0].res_name.lower().strip()
    if nuc_end3[0] in ('r', 'R'):
        nuc_end3 = nuc_end3[1:]
    if nuc_end3[-1] == '3':
        nuc_end3 = nuc_end3[0]
        for a in c_out.residues[-1].atoms:
            a.res_name = '  %s' % nuc_end3.upper()
    dimer_seq3 = '%s%s' % (nuc_end3, args.nuc3.lower())

    print ('The last nucleotide of the PDB is ', nuc_end3.upper(), '.')
    print ('Now ', args.nuc3.upper(),' is being added by using a dimer', dimer_seq3.upper(),'.')

    """ Prepare standard dimer configurations from nab """
    p = PdbFile('%s/%s.pdb' % (DIR_NAB, dimer_seq3))
    p.open_to_read()
    dimer3 = p.read_all()[0]
    p.close()

    """ Residues for superposition """
    dimer_res = dimer3.residues[0]  # Superposition is done by the 1st nucleotide of the dimer '''
    target_res = chains[0].residues[-1]

    mtx3 = dimer_superposition(dimer_res, target_res, nuc_end3)

    ''' Delete HO3' form the last nucleotide (it will be the 2nd last nucleotide) '''
    for a in c_out.residues[-1].atoms:
        if a.name.strip() == 'HO3%s' % SUGAR_MARK:
            c_out.residues[-1].atoms.remove(a)
            break

    ''' Push new last nucleotide '''
    chain_id = c_out.residues[-1].atoms[-1].chain_id
    res_seq = c_out.residues[-1].atoms[-1].res_seq

    r = Residue()
    for a in dimer3.residues[1].atoms:
        acopy = copy.deepcopy(a)
        acopy.chain_id = chain_id
        acopy.res_seq = res_seq + 1
        if acopy.res_name[-1] == '3':
            acopy.res_name = '  %s' % acopy.res_name.strip()[0]
        r.push_atom( acopy )
        r.atoms[-1].xyz.put_as_list( mtx3.do_to_array( acopy.xyz.get_as_list() ) )

    c_out.residues.append(r)


''' Re-asign serial numbers '''
serial = 0
for r in c_out.residues:
    for a in r.atoms:
        serial += 1
        a.serial = serial


''' Write new PDB file '''
p = PdbFile(args.pdb_out)
p.open_to_write()
p.write_all([c_out,])
p.close()
