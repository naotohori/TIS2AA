#!/usr/bin/env python

''' Selection of a database '''
SUGAR_MARK = "'"
LIBPDBAA = 'RNA09_FRAG_AA/'  # input
LIBPDBCG = 'RNA09_FRAG_CG/'  # output
LISTFILE = 'RNA09.nts'

PSEUDO_LIMIT = 2.5   ## +/- 2.5  ===>  bins of 5 degree
RMSD_ACCEPT = 2.0

from cafysis.elements.pdb import Chain, Residue, Atom
from cafysis.file_io.pdb import PdbFile
from cafysis.torsion import torsion
from CalcROT import calcrotation
import sys
import numpy as np

#p = PdbFile(sys.argv[1])
p = PdbFile('f024_run019_2000.pdb')
p.open_to_read()
chains_cg = p.read_all()
p.close()

f_log = open('tis2aa.log','w')

#lib_augc = []
#lib_pucker = []
lib_pseudo = []
for l in open(LISTFILE):
    lsp = l.split()
    #lib_augc.append( lsp[7] )
    #lib_pucker.append( lsp[8] )
    lib_pseudo.append( (float(lsp[9]), float(lsp[10])) )


def angle_diff(ang1, ang2):
    d = abs(ang1 - ang2)
    if d > 180.0:
        d = d - 180.0
    return d

for c_cg in chains_cg:

    seq = []
    xyzs_P = []
    xyzs_S = []
    xyzs_B = []

    for ir, r in enumerate(c_cg.residues):

        flg_P = False; flg_S = False; flg_B = False;
        for a in r.atoms:
            if a.name.strip() == 'P':
                xyzs_P.append( a.xyz.get_as_list() )
                flg_P = True
            elif a.name.strip() == 'S':
                xyzs_S.append( a.xyz.get_as_list() )
                flg_S = True
            elif a.name.strip()[1] == 'b':
                xyzs_B.append( a.xyz.get_as_list() )
                flg_B = True
        if not flg_P:
            xyzs_P.append( False )
        if not flg_S:
            print 'no S in cg model, ir=', ir
        if not flg_B:
            print 'no B in cg model, ir=', ir

        res_name = r.atoms[0].res_name.strip()[-1]
        if res_name in ('A','U','G','C'):
            seq.append( res_name )
        else:
            print 'res_name is not valid:', r.atoms[0].res_name

    nres = len(c_cg.residues)
    c_aa = Chain()
    
    # for the first nucleotide
    c_aa.push_residue( Residue() )

    for ir in range(1, nres-1):

        xyzS1 = np.array( xyzs_S[ir-1] )
        xyzP2 = np.array( xyzs_P[ir] )
        xyzS2 = np.array( xyzs_S[ir] )
        xyzB2 = np.array( xyzs_B[ir] )
        xyzP3 = np.array( xyzs_P[ir+1] )
        xyzS3 = np.array( xyzs_S[ir+1] )

        eta    = torsion(xyzS1, xyzP2, xyzS2, xyzP3, flg_degree=True, flg_360=True)
        theta  = torsion(xyzP2, xyzS2, xyzP3, xyzS3, flg_degree=True, flg_360=True)
        

        limit = PSEUDO_LIMIT
        
        flg_found = False
        while not flg_found:

            lib_cand = []
            angle_cand = []
            for ilib, pseudo in enumerate(lib_pseudo):
                d_eta   = angle_diff(eta,   pseudo[0])
                d_theta = angle_diff(theta, pseudo[1])

                if d_eta <= limit and d_theta <= limit:
                    lib_cand.append(ilib+1)
                    angle_cand.append( (d_eta, d_theta) )
    
            xyzPSBP_ref = [ xyzP2, xyzS2, xyzB2, xyzP3 ]
    
            best_lib = 0
            best_rmsd = 999999.9
            best_mtx = None
            best_angle = None
            
            for lib, angle in zip(lib_cand, angle_cand):
    
                p = PdbFile('%s/%06i_%s.pdb' % (LIBPDBCG, lib, seq[ir]))
                p.open_to_read()
                c = p.read_all()[0]
                xyzPSBP = [ c.residues[1].atoms[0].xyz.get_as_list(), # P2
                            c.residues[1].atoms[1].xyz.get_as_list(), # S2
                            c.residues[1].atoms[2].xyz.get_as_list(), # B2
                            c.residues[2].atoms[0].xyz.get_as_list()] # P3
    
                rmsd, mtx = calcrotation(np.transpose(xyzPSBP_ref), np.transpose(xyzPSBP))
    
                if rmsd < best_rmsd:
                    best_rmsd = rmsd
                    best_lib = lib
                    best_mtx = mtx
                    best_angle = angle
    
                f_log.write("#%3i %06i %f  %f %f\n" % (ir, lib, rmsd, angle[0], angle[1]))

            if best_rmsd <= RMSD_ACCEPT:
                flg_found = True
            else:
                limit = limit + PSEUDO_LIMIT

        f_log.write('%3i %6i %f  %f %f\n' % (ir, len(lib_cand), best_rmsd, best_angle[0], best_angle[1]))

f_log.close()
