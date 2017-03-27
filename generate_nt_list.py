#!/usr/bin/env python

MIN_NUM_ATOM_PER_NT = 22

TARGET = 'RNA05'
PDBGLOB = 'RNA05/pdb/*.pdb'
SUGAR_SIGN = "*"
#TARGET = 'RNAlibraries'
#PDBGLOB = 'RNAlibraries/ROT-5/*.pdb'
#SUGAR_SIGN = "'"

import math
import numpy as np
import glob
import os
from cafysis.file_io.pdb import PdbFile
from cafysis.torsion import torsion


def check_necessary_atoms_exist(res, res_name):
    '''
    A nucleotide should have all necessary atom types.
    Each atom type should appear only once.
    '''
    atomnames = [atom.name.strip() for atom in res.atoms]

    # Phosphate 
    if ( atomnames.count("P"  ) != 1 or 
         atomnames.count("O1P") != 1 or atomnames.count("O2P") != 1 ):
        return False

    # Sugar
    if ( atomnames.count("O5%s" % SUGAR_SIGN) != 1 or atomnames.count("C5%s" % SUGAR_SIGN) != 1 or
         atomnames.count("C4%s" % SUGAR_SIGN) != 1 or atomnames.count("O4%s" % SUGAR_SIGN) != 1 or
         atomnames.count("C3%s" % SUGAR_SIGN) != 1 or atomnames.count("O3%s" % SUGAR_SIGN) != 1 or
         atomnames.count("C2%s" % SUGAR_SIGN) != 1 or atomnames.count("O2%s" % SUGAR_SIGN) != 1 or
         atomnames.count("C1%s" % SUGAR_SIGN) != 1 ):
        return False

    # Base 
    if ( atomnames.count("N1" ) != 1 or atomnames.count("C2" ) != 1 or
         atomnames.count("N3" ) != 1 or atomnames.count("C4" ) != 1 or
         atomnames.count("C5" ) != 1 or atomnames.count("C6" ) != 1 ):
        return False

    if res_name in ("C","U"):
        if atomnames.count("O2") != 1:
            return False
        if res_name == "C":
            if atomnames.count("N4") != 1:
                return False
        if res_name == "U":
            if atomnames.count("O4") != 1:
                return False

    if res_name in ("A","G"):
        if atomnames.count("N7") != 1:
            return False
        if atomnames.count("C8") != 1:
            return False
        if atomnames.count("N9") != 1:
            return False
        if res_name == "G":
            if atomnames.count("N2") != 1:
                return False
            if atomnames.count("O6") != 1:
                return False
        if res_name == "A":
            if atomnames.count("N6") != 1:
                return False

    return True

'''
def check_clash(chains)

    for c in chains:
        natom = c.hum_atom()

        for i1 in range(natom):
            a1 = c.get_atom(i1)
            radi_a1 = get_atom_radius(a1.name)

            for i2 in range(i1+1,natom):
                a2 = c.get_atom(i2)
                radi_a2 = get_atom_radius(a2.name)

                d = a2.xyz.distance( a1.xyz )
                if d < (radi_a1 + radi_a2 - 0.4):
                    return False

    return True
'''

def perpendicular_distance(a, b, c):

    ba = b - a
    ca = c - a
    H = a + np.dot(ca,ba) * ba / np.dot(ba,ba)
    d = math.sqrt( np.dot(c-H,c-H) )
    return d


def calc_pucker(res, res_next):

    for a in res.atoms:
        if a.name.strip() == "C5%s" % SUGAR_SIGN:
            xyzC5 = a.xyz.get_as_ndarray()
        elif a.name.strip() == "C4%s" % SUGAR_SIGN:
            xyzC4 = a.xyz.get_as_ndarray()
        elif a.name.strip() == "C3%s" % SUGAR_SIGN:
            xyzC3 = a.xyz.get_as_ndarray()
        elif a.name.strip() == "O3%s" % SUGAR_SIGN:
            xyzO3 = a.xyz.get_as_ndarray()
        elif a.name.strip() == "C1%s" % SUGAR_SIGN:
            xyzC1 = a.xyz.get_as_ndarray()
        elif ((a.name.strip() == "N1" and a.res_name.strip() in ("U","C")) or
              (a.name.strip() == "N9" and a.res_name.strip() in ("A","G")) ):
            xyzN1N9 = a.xyz.get_as_ndarray()

    for a in res_next.atoms:
        if a.name.strip() == "P":
            xyzP_next = a.xyz.get_as_ndarray()

    delta = torsion(xyzC5, xyzC4, xyzC3, xyzO3, flg_degree=True, flg_360=True)

    d = perpendicular_distance(xyzN1N9, xyzC1, xyzP_next)

    if (54.0 <= delta) and (delta <= 114.0) and (2.9 < d):
        return 3
    elif (117.0 <= delta) and (delta <=177.0) and (d <=2.9):
        return 2
    else:
        return 1

seq = []
                  # False if ...
flg_restype = []  #   not either of AUGC
flg_natom   = []  #   number of atoms is less than 22 (to remove terminal phosphate etc.)
flg_Bfact   = []  #   any atom has B-factor more than 60
flg_noterm  = []  #   terminal nucleotide of the chain
flg_atoms   = []  #   (after above tests) any of necessary atom types does not exist
flg_clash   = []  #   (after above tests) any atom clash with vdW radii overlap > 0.4 angstrom

pdbfiles = glob.glob(PDBGLOB)

for pdbfile in pdbfiles:

    pdb = PdbFile(pdbfile)
    pdb.flg_HETATM = True   # Accept HETATM !!
    pdb.open_to_read()
    chains = pdb.read_all()

    for ic, c in enumerate(chains):

        for ir, r in enumerate(c.residues):

            flg_so_far = True
            resname = r.atoms[0].res_name.strip()

            # Number of atoms in the nt
            if len(r.atoms) < MIN_NUM_ATOM_PER_NT:  # To remove terminal phosphate and others
                flg_natom.append(False)
                flg_so_far = False
            else:
                flg_natom.append(True)
                #print pdbfile, resname

            # Sequence
            if resname in ("A", "C", "U", "G"):
                flg_restype.append(True)
                seq.append(resname)
            elif resname in ("HOH","CA","MN","BRO","ROB"):
                flg_restype.append(False)
                flg_so_far = False
                seq.append("o") # others
            else:
                flg_restype.append(False)
                flg_so_far = False
                seq.append("m") # modified

            # B-factor
            flg = True
            for a in r.atoms:
                if a.temp_factor > 60.0:
                    flg = False
                    flg_so_far = False
            flg_Bfact.append(flg)

            # Terminal or not
            if ir == 0 or ir == (c.num_res() - 1):
                flg_noterm.append( False )
                flg_so_far = False
            else:
                flg_noterm.append( True )

            # All necessary atoms
            if flg_so_far:
                flg = check_necessary_atoms_exist(r,resname)
                flg_atoms.append( flg )
                flg_so_far = flg
            else:
                flg_atoms.append( False )

            '''
            # Clash
            if flg_so_far:
                flg = check_clash()
                flg_clash.append( flg )
                flg_so_far = flg
            else:
                flg_clash.append( False )
            '''


nres = len(seq)

flg_final = []  
for i in range(nres):

    # Terminal nucleotides never pass.
    if not flg_noterm[i]:
        flg_final.append( False )
        continue

    # Itself passes all the quality test, AND both neighboring nucleotides are also qualified.
    if (flg_restype[i-1] and flg_natom[i-1] and flg_Bfact[i-1] and flg_atoms[i-1] and 
        flg_restype[i  ] and flg_natom[i  ] and flg_Bfact[i-1] and flg_atoms[i  ] and 
        flg_restype[i+1] and flg_natom[i+1] and flg_Bfact[i+1] and flg_atoms[i+1] ):
        flg_final.append( True )
    else:
        flg_final.append( False )


##################################
## Output list
##################################
f_out = open('%s.nts' % TARGET,'w')

puckers = []

iseq = -1
for pdbfile in pdbfiles:

    pdb = PdbFile(pdbfile)
    pdb.flg_HETATM = True   # Accept HETATM !!
    pdb.open_to_read()
    chains = pdb.read_all()

    for ic, c in enumerate(chains):

        for ir, r in enumerate(c.residues):

            iseq += 1

            resname = r.atoms[0].res_name.strip()
            f_out.write('%13s %2i %2s %4i %4i %4s %2i' 
                     % (os.path.basename(pdbfile)[:-4], ic+1, r.atoms[0].chain_id, 
                        ir+1, r.atoms[0].res_seq, resname, len(r.atoms)))

            # Number of atoms in the nt
            if flg_natom[iseq]:
                f_out.write(' o' )
            else:
                f_out.write(' x' )

            # Sequence
            f_out.write(' %s' % seq[iseq])

            # B-factor
            if flg_Bfact[iseq]:
                f_out.write(' o' )
            else:
                f_out.write(' x' )

            # Terminal or not
            if flg_noterm[iseq]:
                f_out.write(' o' )
            else:
                f_out.write(' x' )

            if flg_atoms[iseq]:
                f_out.write(' o' )
            else:
                f_out.write(' x' )

            # Pucker
            if flg_final[iseq]:
                puc = calc_pucker(r, c.residues[ir+1])
                if puc == 1:
                    flg_final[iseq] = False
            else:
                puc = 0
            f_out.write(' %i' % puc)

            if flg_final[iseq]:
                f_out.write(' o' )
            else:
                f_out.write(' x' )


            f_out.write('\n')


print len(pdbfiles)
print len(seq)
print sum(flg_restype)
print sum(flg_natom)
print sum(flg_Bfact)
print sum(flg_final)

