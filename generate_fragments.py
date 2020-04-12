#!/usr/bin/env python

import math
import numpy as np
import glob
import os
from pdbfile import PdbFile
from pdb import Chain
from torsion import torsion

MIN_NUM_ATOM_PER_NT = 22
flg_PDB = True  # True: output PDB for each fragment
# False: only make a list

TARGET = 'RNA09'
PDBGLOB = 'RNA09/pdb/*.pdb'
SUGAR_MARK = "'"
LIBPDBAA = 'RNA09_FRAG_AA/'


# TARGET = 'RNA05'
# PDBGLOB = 'RNA05/pdb/*.pdb'
# SUGAR_MARK = "*"
# LIBPDBAA = 'RNA05_FRAG_AA/'

# TARGET = 'RNAlibraries'
# PDBGLOB = 'RNAlibraries/ROT-5/*.pdb'
# SUGAR_MARK = "'"


def check_necessary_atoms_exist(res, res_name):
    """
    A nucleotide should have all necessary atom types.
    Each atom type should appear only once.
    """
    atomnames = [atom.name.strip() for atom in res.atoms]

    # Phosphate 
    if ((atomnames.count("P") != 1 or
         atomnames.count("O1P") != 1 or atomnames.count("O2P") != 1)
            and
            (atomnames.count("P") != 1 or
             atomnames.count("OP1") != 1 or atomnames.count("OP2") != 1)):
        return False

    # Sugar
    if (atomnames.count("O5%s" % SUGAR_MARK) != 1 or atomnames.count("C5%s" % SUGAR_MARK) != 1 or
            atomnames.count("C4%s" % SUGAR_MARK) != 1 or atomnames.count("O4%s" % SUGAR_MARK) != 1 or
            atomnames.count("C3%s" % SUGAR_MARK) != 1 or atomnames.count("O3%s" % SUGAR_MARK) != 1 or
            atomnames.count("C2%s" % SUGAR_MARK) != 1 or atomnames.count("O2%s" % SUGAR_MARK) != 1 or
            atomnames.count("C1%s" % SUGAR_MARK) != 1 or
            atomnames.count("H5%s" % SUGAR_MARK) != 1 or atomnames.count("H5%s%s" % (SUGAR_MARK, SUGAR_MARK)) != 1 or
            atomnames.count("H4%s" % SUGAR_MARK) != 1 or atomnames.count("H3%s" % SUGAR_MARK) != 1 or
            atomnames.count("H2%s" % SUGAR_MARK) != 1 or atomnames.count("HO2%s" % SUGAR_MARK) != 1 or
            atomnames.count("H1%s" % SUGAR_MARK) != 1):
        return False

    # Base 
    if (atomnames.count("N1") != 1 or atomnames.count("C2") != 1 or
            atomnames.count("N3") != 1 or atomnames.count("C4") != 1 or
            atomnames.count("C5") != 1 or atomnames.count("C6") != 1):
        return False

    if res_name in ("C", "U"):
        if atomnames.count("O2") != 1:
            return False
        if res_name == "C":
            if atomnames.count("N4") != 1:
                return False
        if res_name == "U":
            if atomnames.count("O4") != 1:
                return False

    if res_name in ("A", "G"):
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


def perpendicular_distance(a, b, c):
    """
    return: perpendicular distance from the point C 
            to the line connecting A and B
    """

    ba = b - a
    ca = c - a
    H = a + np.dot(ca, ba) * ba / np.dot(ba, ba)
    d = math.sqrt(np.dot(c - H, c - H))
    return d


def calc_pucker(res, res_next):
    """
    args: two consequtive residues
    return: 1 = neither 2'-end nor 3'-end
            2 = 2'-end
            3 = 3'-end
    """
    xyzC5 = res.find_atom_by_name("C5%s" % SUGAR_MARK).xyz.get_as_ndarray()
    xyzC4 = res.find_atom_by_name("C4%s" % SUGAR_MARK).xyz.get_as_ndarray()
    xyzC3 = res.find_atom_by_name("C3%s" % SUGAR_MARK).xyz.get_as_ndarray()
    xyzO3 = res.find_atom_by_name("O3%s" % SUGAR_MARK).xyz.get_as_ndarray()
    xyzC1 = res.find_atom_by_name("C1%s" % SUGAR_MARK).xyz.get_as_ndarray()
    if res.atoms[0].res_name.strip() in ("U", "C"):
        xyzN1N9 = res.find_atom_by_name("N1").xyz.get_as_ndarray()
    elif res.atoms[0].res_name.strip() in ("A", "G"):
        xyzN1N9 = res.find_atom_by_name("N9").xyz.get_as_ndarray()

    #    for a in res.atoms:
    #        if a.name.strip() == "C5%s" % SUGAR_MARK:
    #            xyzC5 = a.xyz.get_as_ndarray()
    #        elif a.name.strip() == "C4%s" % SUGAR_MARK:
    #            xyzC4 = a.xyz.get_as_ndarray()
    #        elif a.name.strip() == "C3%s" % SUGAR_MARK:
    #            xyzC3 = a.xyz.get_as_ndarray()
    #        elif a.name.strip() == "O3%s" % SUGAR_MARK:
    #            xyzO3 = a.xyz.get_as_ndarray()
    #        elif a.name.strip() == "C1%s" % SUGAR_MARK:
    #            xyzC1 = a.xyz.get_as_ndarray()
    #        elif ((a.name.strip() == "N1" and a.res_name.strip() in ("U","C")) or
    #              (a.name.strip() == "N9" and a.res_name.strip() in ("A","G")) ):
    #            xyzN1N9 = a.xyz.get_as_ndarray()

    xyzP_next = res_next.find_atom_by_name("P").xyz.get_as_ndarray()
    #    for a in res_next.atoms:
    #        if a.name.strip() == "P":
    #            xyzP_next = a.xyz.get_as_ndarray()

    delta = torsion(xyzC5, xyzC4, xyzC3, xyzO3, flg_degree=True, flg_360=True)

    d = perpendicular_distance(xyzN1N9, xyzC1, xyzP_next)

    if (54.0 <= delta) and (delta <= 114.0) and (2.9 < d):
        return 3
    elif (117.0 <= delta) and (delta <= 177.0) and (d <= 2.9):
        return 2
    else:
        return 1


def calc_pseudo(r1, r2, r3, phos="P", sugar="C4", flg_degree=True, flg_360=True):
    if sugar == 'C4':
        for a in r1.atoms:
            if a.name.strip() == "C4%s" % SUGAR_MARK:
                xyzS1 = a.xyz.get_as_ndarray()
        for a in r2.atoms:
            if a.name.strip() == "C4%s" % SUGAR_MARK:
                xyzS2 = a.xyz.get_as_ndarray()
        for a in r3.atoms:
            if a.name.strip() == "C4%s" % SUGAR_MARK:
                xyzS3 = a.xyz.get_as_ndarray()

    elif sugar == 'geo':
        xyzS1 = np.zeros((3,))
        n = 0
        for a in r1.atoms:
            if (a.name.strip() == "C1%s" % SUGAR_MARK or
                    a.name.strip() == "C2%s" % SUGAR_MARK or
                    a.name.strip() == "O2%s" % SUGAR_MARK or
                    a.name.strip() == "C3%s" % SUGAR_MARK or
                    a.name.strip() == "O3%s" % SUGAR_MARK or
                    a.name.strip() == "C4%s" % SUGAR_MARK or
                    a.name.strip() == "O4%s" % SUGAR_MARK or
                    a.name.strip() == "C5%s" % SUGAR_MARK or
                    a.name.strip() == "O5%s" % SUGAR_MARK):
                xyzS1 += a.xyz.get_as_ndarray()
                n += 1
        if n == 9:
            xyzS1 /= 9.0
        else:
            print('Error: n != 9 for xyzS1 in calc_pseudo')
            sys.exit(2)

        xyzS2 = np.zeros((3,))
        n = 0
        for a in r2.atoms:
            if (a.name.strip() == "C1%s" % SUGAR_MARK or
                    a.name.strip() == "C2%s" % SUGAR_MARK or
                    a.name.strip() == "O2%s" % SUGAR_MARK or
                    a.name.strip() == "C3%s" % SUGAR_MARK or
                    a.name.strip() == "O3%s" % SUGAR_MARK or
                    a.name.strip() == "C4%s" % SUGAR_MARK or
                    a.name.strip() == "O4%s" % SUGAR_MARK or
                    a.name.strip() == "C5%s" % SUGAR_MARK or
                    a.name.strip() == "O5%s" % SUGAR_MARK):
                xyzS2 += a.xyz.get_as_ndarray()
                n += 1
        if n == 9:
            xyzS2 /= 9.0
        else:
            print('Error: n != 9 for xyzS2 in calc_pseudo')
            sys.exit(2)

        xyzS3 = np.zeros((3,))
        n = 0
        for a in r3.atoms:
            if (a.name.strip() == "C1%s" % SUGAR_MARK or
                    a.name.strip() == "C2%s" % SUGAR_MARK or
                    a.name.strip() == "O2%s" % SUGAR_MARK or
                    a.name.strip() == "C3%s" % SUGAR_MARK or
                    a.name.strip() == "O3%s" % SUGAR_MARK or
                    a.name.strip() == "C4%s" % SUGAR_MARK or
                    a.name.strip() == "O4%s" % SUGAR_MARK or
                    a.name.strip() == "C5%s" % SUGAR_MARK or
                    a.name.strip() == "O5%s" % SUGAR_MARK):
                xyzS3 += a.xyz.get_as_ndarray()
                n += 1
        if n == 9:
            xyzS3 /= 9.0
        else:
            print('Error: n != 9 for xyzS2 in calc_pseudo')

    if phos == 'P':
        for a in r2.atoms:
            if a.name.strip() == "P":
                xyzP2 = a.xyz.get_as_ndarray()
        for a in r3.atoms:
            if a.name.strip() == "P":
                xyzP3 = a.xyz.get_as_ndarray()
    elif phos == 'PO2':
        xyzP2 = np.zeros((3,))
        n = 0
        for a in r2.atoms:
            if (a.name.strip() == "P" or
                    a.name.strip() == "OP1" or a.name.strip() == "OP2" or
                    a.name.strip() == "O1P" or a.name.strip() == "O2P"):
                xyzP2 += a.xyz.get_as_ndarray()
                n += 1
        if n == 3:
            xyzP2 /= 3.0
        else:
            print('Error: n != 3 for xyzP2 in calc_pseudo')

        xyzP3 = np.zeros((3,))
        n = 0
        for a in r3.atoms:
            if (a.name.strip() == "P" or
                    a.name.strip() == "OP1" or a.name.strip() == "OP2" or
                    a.name.strip() == "O1P" or a.name.strip() == "O2P"):
                xyzP3 += a.xyz.get_as_ndarray()
                n += 1
        if n == 3:
            xyzP3 /= 3.0
        else:
            print('Error: n != 3 for xyzP3 in calc_pseudo')

    eta = torsion(xyzS1, xyzP2, xyzS2, xyzP3, flg_degree=flg_degree, flg_360=flg_360)
    theta = torsion(xyzP2, xyzS2, xyzP3, xyzS3, flg_degree=flg_degree, flg_360=flg_360)

    return eta, theta


if __name__ == "__main__":

    seq = []
    # False if ...
    flg_natom = []  # Number of atoms is less than 22 (to remove terminal phosphate etc.)
    #    (this is for checking consistency to Humphris-Narayanan & Pyle)
    flg_restype = []  # Neither of AUGC
    flg_Bfact = []  # Any atom has B-factor more than 60
    flg_noterm = []  # Terminal nucleotide of the chain
    #   (after above tests)
    flg_atoms = []  # Any of necessary atom types does not exist
    # flg_clash   = []  #   Any atom clash with vdW radii overlap > 0.4 angstrom

    pdbfiles = glob.glob(PDBGLOB)

    for pdbfile in pdbfiles:

        pdb = PdbFile(pdbfile)
        pdb.flg_HETATM = True  # Accept HETATM !!
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
                    # print pdbfile, resname

                # Sequence
                if resname in ("A", "C", "U", "G"):
                    flg_restype.append(True)
                    seq.append(resname)
                elif resname in ("HOH", "CA", "MN", "BRO", "ROB"):
                    flg_restype.append(False)
                    flg_so_far = False
                    seq.append("o")  # others
                else:
                    flg_restype.append(False)
                    flg_so_far = False
                    seq.append("m")  # modified

                # B-factor
                flg = True
                for a in r.atoms:
                    if a.temp_factor > 60.0:
                        flg = False
                        flg_so_far = False
                flg_Bfact.append(flg)

                # Terminal or not
                if ir == 0 or ir == (c.num_res() - 1):
                    flg_noterm.append(False)
                    flg_so_far = False
                else:
                    flg_noterm.append(True)

                # All necessary atoms
                if flg_so_far:
                    flg = check_necessary_atoms_exist(r, resname)
                    flg_atoms.append(flg)
                    flg_so_far = flg
                else:
                    flg_atoms.append(False)

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
            flg_final.append(False)
            continue

        # Itself passes all the quality test, AND both neighboring nucleotides are also qualified.
        if (flg_restype[i - 1] and flg_natom[i - 1] and flg_Bfact[i - 1] and flg_atoms[i - 1] and
                flg_restype[i] and flg_natom[i] and flg_Bfact[i - 1] and flg_atoms[i] and
                flg_restype[i + 1] and flg_natom[i + 1] and flg_Bfact[i + 1] and flg_atoms[i + 1]):
            flg_final.append(True)
        else:
            flg_final.append(False)

    # #################################
    # # Output list
    # #################################
    f_out = open('%s.nts_all' % TARGET, 'w')
    f_final = open('%s.nts' % TARGET, 'w')

    puckers = []

    iseq = -1
    ifinal = 0
    for pdbfile in pdbfiles:

        pdb = PdbFile(pdbfile)
        pdb.flg_HETATM = True  # Accept HETATM !!
        pdb.open_to_read()
        chains = pdb.read_all()

        for ic, c in enumerate(chains):

            for ir, r in enumerate(c.residues):

                iseq += 1

                resname = r.atoms[0].res_name.strip()

                if resname == "HOH":  # Do not output water
                    continue

                f_out.write('%13s %2i %2s %4i %4i %4s %2i'
                            % (os.path.basename(pdbfile)[:-4], ic + 1, r.atoms[0].chain_id,
                               ir + 1, r.atoms[0].res_seq, resname, len(r.atoms)))

                # Number of atoms in the nt
                if flg_natom[iseq]:
                    f_out.write(' o')
                else:
                    f_out.write(' N')

                # Sequence
                f_out.write(' %s' % seq[iseq])

                # B-factor
                if flg_Bfact[iseq]:
                    f_out.write(' o')
                else:
                    f_out.write(' B')

                # Terminal or not
                if flg_noterm[iseq]:
                    f_out.write(' o')
                else:
                    f_out.write(' T')

                if flg_atoms[iseq]:
                    f_out.write(' o')
                else:
                    f_out.write(' A')

                # Pucker
                if flg_final[iseq]:
                    puc = calc_pucker(r, c.residues[ir + 1])
                    if puc == 1:
                        flg_final[iseq] = False
                else:
                    puc = 0
                puckers.append(puc)
                f_out.write(' %i' % puc)

                if flg_final[iseq]:
                    f_out.write(' o')
                else:
                    f_out.write(' x')

                if puc in (2, 3):
                    ifinal += 1

                    f_final.write('%6i %6i %13s %2i %2s %4i %4i %s'
                                  % (ifinal, iseq, os.path.basename(pdbfile)[:-4], ic + 1, r.atoms[0].chain_id,
                                     ir + 1, r.atoms[0].res_seq, seq[iseq]))
                    f_final.write('  %i' % puc)

                    eta1, theta1 = calc_pseudo(c.residues[ir - 1], r, c.residues[ir + 1], phos='P', sugar='C4')
                    eta2, theta2 = calc_pseudo(c.residues[ir - 1], r, c.residues[ir + 1], phos='PO2', sugar='geo')

                    f_final.write('  %5.1f %5.1f  %5.1f %5.1f' % (eta1, theta1, eta2, theta2))
                    f_final.write('\n')

                    if flg_PDB:
                        cout = Chain()
                        cout.push_residue(c.residues[ir - 1])
                        cout.push_residue(c.residues[ir])
                        cout.push_residue(c.residues[ir + 1])

                        pdb_out = PdbFile('%s/%06i.pdb' % (LIBPDBAA, ifinal))
                        pdb_out.open_to_write()
                        pdb_out.write_all([cout, ])
                        pdb_out.close()

                f_out.write('\n')

    f_out.write('# PDB files: %i\n' % len(pdbfiles))
    f_out.write('# residues : %i\n' % len(seq))
    f_out.write('# restype filter: %i\n' % sum(flg_restype))
    f_out.write('# number-of-atoms filter: %i\n' % sum(flg_natom))
    f_out.write('# B-factor filter: %i\n' % sum(flg_Bfact))
    f_out.write('# Terminal filter: %i\n' % sum(flg_noterm))
    f_out.write('# atomtype filter: %i\n' % sum(flg_atoms))
    f_out.write('# Final: %i\n' % sum(flg_final))
    f_out.write('#     2: %i\n' % puckers.count(2))
    f_out.write('#     3: %i\n' % puckers.count(3))
