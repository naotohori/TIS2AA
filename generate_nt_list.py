#!/usr/bin/env python

MIN_NUM_ATOM_PER_NT = 22

import glob
import os
from cafysis.file_io.pdb import PdbFile

TARGET = 'RNA05'
pdbfiles = glob.glob('RNA05/pdb/*.pdb')

f_out = open('%s.nts' % TARGET,'w')


seq = []

flg_restype = []
flg_natom   = []
flg_Bfact   = []  # True: OK,  False: more than 60
flg_noterm  = []


for pdbfile in pdbfiles:

    pdb = PdbFile(pdbfile)
    pdb.flg_HETATM = True   # Accept HETATM
    pdb.open_to_read()
    chains = pdb.read_all()

    for ic, c in enumerate(chains):

        for ir, r in enumerate(c.residues):

            resname = r.atoms[0].res_name.strip()
            f_out.write('%10s %2i %2s %4i %4i %4s %2i' 
                     % (os.path.basename(pdbfile)[:-4], ic+1, r.atoms[0].chain_id, 
                        ir+1, r.atoms[0].res_seq, resname, len(r.atoms)))

            # Number of atoms in the nt
            if len(r.atoms) < MIN_NUM_ATOM_PER_NT:  # To remove terminal phosphate and others
                flg_natom.append(False)
                f_out.write(' x' )
            else:
                flg_natom.append(True)
                f_out.write(' o' )
                #print pdbfile, resname

            # Sequence
            if resname in ('A', 'C', 'U', 'G'):
                flg_restype.append(True)
                f_out.write(' o')
                seq.append(resname)
            elif resname in ('HOH','CA','MN','BRO','ROB'):
                flg_restype.append(False)
                f_out.write(' x' )
                seq.append('o') # others
            else:
                flg_restype.append(False)
                f_out.write(' x' )
                seq.append('m') # modified

            # B-factor
            flg = True
            for a in r.atoms:
                if a.temp_factor > 60.0:
                    flg = False
            if flg:
                f_out.write(' o' )
            else:
                f_out.write(' x' )
            flg_Bfact.append(flg)

            # Terminal or not
            if ir == 0 or ir == c.num_res()-1:
                flg_noterm.append( False )
                f_out.write(' x' )
            else:
                flg_noterm.append( True )
                f_out.write(' o' )

            f_out.write('\n')



nres = len(seq)
flg_final = [] # Both ends are continus to qualified nucleotides.
for i in range(nres):
    if not flg_noterm[i]:
        flg_final.append( False )
    else:
        if (flg_restype[i-1] and flg_natom[i-1] and flg_Bfact[i-1] and
            flg_restype[i  ] and flg_natom[i  ] and flg_Bfact[i-1] and
            flg_restype[i+1] and flg_natom[i+1] and flg_Bfact[i+1] ):
            flg_final.append( True )
        else:
            flg_final.append( False )


print len(pdbfiles)
print len(seq)
print sum(flg_restype)
print sum(flg_natom)
print sum(flg_Bfact)
print sum(flg_final)

