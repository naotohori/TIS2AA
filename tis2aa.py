#!/usr/bin/env python

""" Selection of a database """
SUGAR_MARK = "'"
LIBPDBAA = 'RNA09_FRAG_AA/'  # input
LIBPDBCG = 'RNA09_FRAG_CG/'  # output
LISTFILE = 'RNA09.nts'

""" Parameters to search library """
PSEUDO_LIMIT = 2.0   ## +/- 2.5  ===>  bins of 5 degree
PSEUDO_MAX = 10.0
RMSD_ACCEPT = 1.0
RMSD_MAX = 6.0
flg_adjust_P_position = False

from cafysis.elements.pdb import Chain, Residue, Atom
from cafysis.elements.coord import Coord
from cafysis.file_io.pdb import PdbFile
from cafysis.torsion import torsion
from cafysis.mtx_coord_transform import mtx_crd_transform
from CalcROT import calcrotation
import sys
import numpy as np
import copy


def angle_diff(ang1, ang2):
    """ Calculate deviation of angles defined between 0 and 360 degree """
    d = abs(ang1 - ang2)
    if d > 180.0:
        d = d - 180.0
    return d


if __name__ == "__main__":

    if not len(sys.argv) in (3,4):
        print 'Usage: SCRIPT [input CG PDB] [output log file] [output AA PDB]'
        print '  or : SCRIPT [input CG PDB] [output AA PDB]  (log file will be named tis2aa.log)'
        sys.exit(2)
    
    if len(sys.argv) == 4:
        f_log = open(sys.argv[2],'w')
    else:
        f_log = open('tis2aa.log','w')


    """ Read CG PDB """
    p = PdbFile(sys.argv[1])
    p.open_to_read()
    chains_cg = p.read_all()
    p.close()
    
    """ Open AA PDB """
    p = PdbFile(sys.argv[-1])
    p.open_to_write()

    
    """ Read pseudo angles from library list (.nts) file """
    #lib_pucker = []
    lib_pseudo = []
    for l in open(LISTFILE):
        lsp = l.split()
        #lib_pucker.append( lsp[8] )
        lib_pseudo.append( (float(lsp[9]), float(lsp[10])) )
    
    
    chains_aa = []
    for c_cg in chains_cg:
            
        seq = []
        xyzs_P = []
        xyzs_S = []
        xyzs_B = []
    
        """ Store sequence and coordinates of CG model """
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

        
        """ For the first residue """
        #c_aa.push_residue( Residue() )
    
        """ From the second to the second last residues """
        for ir in range(1, nres-1):
    
            """ Calcualte pseudo angles """
            xyzS1 = np.array( xyzs_S[ir-1] )
            xyzP2 = np.array( xyzs_P[ir] )
            xyzS2 = np.array( xyzs_S[ir] )
            xyzB2 = np.array( xyzs_B[ir] )
            xyzP3 = np.array( xyzs_P[ir+1] )
            xyzS3 = np.array( xyzs_S[ir+1] )
    
            eta    = torsion(xyzS1, xyzP2, xyzS2, xyzP3, flg_degree=True, flg_360=True)
            theta  = torsion(xyzP2, xyzS2, xyzP3, xyzS3, flg_degree=True, flg_360=True)


            """ For superposition """
            xyzSPS_ref = [ xyzS1, xyzP2, xyzS2 ]
            xyzSPSB_ref = [ xyzS1, xyzP2, xyzS2, xyzB2 ]
            
    
            """ Search best fragment """
            ilimit = 0
            flg_found = False
            while not flg_found:
    
                """ Psuedo-angle range """
                ilimit += 1
                limit = PSEUDO_LIMIT * ilimit
    
                if limit > PSEUDO_MAX:
                    if best_rmsd <= RMSD_MAX:
                        break
                    else:
                        print 'limit > PSEUDO_MAX; could not find library for ir=',ir
                        sys.exit(2)
    

                """ List up candiate fragments by pseudo angles """
                cand_lib = []
                cand_angle = []
                for ilib, pseudo in enumerate(lib_pseudo):
                    d_eta   = angle_diff(eta,   pseudo[0])
                    d_theta = angle_diff(theta, pseudo[1])
    
                    if d_eta <= limit and d_theta <= limit:
                        cand_lib.append(ilib+1)
                        cand_angle.append( (d_eta, d_theta) )
        
        
                """ Evaluate all candidate fragments by RMSD """
                best_lib = 0
                best_rmsd = 9999.9
                best_mtx = None
                best_angle = None
                
                for lib, angle in zip(cand_lib, cand_angle):
        
                    """ Read fragment coordinates """
                    p = PdbFile('%s/%06i_%s.pdb' % (LIBPDBCG, lib, seq[ir]))
                    p.open_to_read()
                    c = p.read_all()[0]
                    xyzSPSB = [ c.residues[0].atoms[1].xyz.get_as_list(), # S1
                                c.residues[1].atoms[0].xyz.get_as_list(), # P2
                                c.residues[1].atoms[1].xyz.get_as_list(), # S2
                                c.residues[1].atoms[2].xyz.get_as_list()] # B2
        
                    rmsd, mat = calcrotation(np.transpose(xyzSPSB_ref), np.transpose(xyzSPSB))
                    p.close()
        
                    if rmsd < best_rmsd:
                        best_rmsd = rmsd
                        best_lib = lib
                        best_angle = angle
                        #xyzSPS = [ c.residues[0].atoms[1].xyz.get_as_list(), # S1
                        #           c.residues[1].atoms[0].xyz.get_as_list(), # P2
                        #           c.residues[1].atoms[1].xyz.get_as_list()] # S2
                        #rmsd, mat = calcrotation(np.transpose(xyzSPS_ref), np.transpose(xyzSPS))
                        best_mtx = mat
        
                    f_log.write("#%3i %06i %f  %f %f\n" % (ir, lib, rmsd, angle[0], angle[1]))
    
                if best_rmsd <= RMSD_ACCEPT:
                    flg_found = True
    
            f_log.write('%3i %6i %f  %f %f\n' % (ir, len(cand_lib), best_rmsd, best_angle[0], best_angle[1]))


            p = PdbFile('%s/%06i_%s.pdb' % (LIBPDBAA, best_lib, seq[ir]))
            p.open_to_read()
            c = p.read_all()[0]
            frag_res = c.residues[1]
            p.close()
    
            mtx = mtx_crd_transform()
            mtx.mtx[:,:] = best_mtx

            r_aa = Residue()
            for a in frag_res.atoms:
                a_new = copy.deepcopy(a)
                a_new.xyz.put_as_list( mtx.do_to_array( a.xyz.get_as_list() ) )
                r_aa.push_atom( a_new )
            
            """ Adjust the position of P to be exactly identical """
            if flg_adjust_P_position:
                aaP = Coord()
                nP = 0
                for a in r_aa.atoms:
                    if a.name.strip() in ("P","OP1", "OP2", "O1P", "O2P"):
                        aaP = aaP + a.xyz
                        nP += 1
                if nP != 3:
                    print 'Error: nP != 3, ir=,',ir
                    sys.exit(2)
                aaP = aaP / float(3)
                trans = xyzP2 - aaP.get_as_ndarray()
                mtx.reset()
                mtx.translation( trans[0], trans[1], trans[2] )
                for a in r_aa.atoms:
                    a.xyz.put_as_list( mtx.do_to_array( a.xyz.get_as_list() ) )

            c_aa.push_residue( r_aa )
    
        chains_aa.append( c_aa )
    
    
    """ Re-numbering serial and residue IDs """
    serial = 0
    res_seq = 0
    for c in chains_aa:
        for r in c.residues:
            res_seq += 1
            for a in r.atoms:
                serial += 1
                a.serial = serial
                a.res_seq = res_seq
                res_name = a.res_name.strip()[-1]
                a.res_name = "  %s" % res_name
        c.reindex()
    

    """ Write AA PDB """
    p = PdbFile(sys.argv[-1])
    p.open_to_write()
    p.write_all( chains_aa )
    p.close()
    
    f_log.close()
    
