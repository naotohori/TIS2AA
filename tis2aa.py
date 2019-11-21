#!/usr/bin/env python
''' Conversion from TIS model to All-atomistic model'''
__author__ = "Naoto Hori"

import os 
import argparse
from pdb import Chain, Residue, Atom
from coord import Coord
from pdbfile import PdbFile
from torsion import torsion
from mtx_coord_transform import mtx_crd_transform
from CalcROT import calcrotation
import sys
import numpy as np
import copy

################################################################################
""" Selection of a database """
SUGAR_MARK = "'"
#BASEDIR = os.path.expanduser("~/TIS2AA/")
BASEDIR = '%s/TIS2AA/' % os.environ['HOME']

""" Parameters to search library """
PSEUDO_BIN = 2.0   ## +/- 2.5  ===>  bins of 5 degree
PSEUDO_MAX = 20.0
RMSD_ACCEPT = 1.0
RMSD_MAX = 6.0


################################################################################
def angle_diff(ang1, ang2):
    """ Calculate deviation of angles defined between 0 and 360 degree """
    d = abs(ang1 - ang2)
    if d > 180.0:
        d = 360.0 - d
    return d


################################################################################
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
             description='Structural conversion from TIS model to All-atomistic model',
             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--binpseudo', dest='pseudo_bin', default=PSEUDO_BIN,
                        action='store', type=float, 
                        help='Half of bin size for pseudo angle') 

    parser.add_argument('--maxpseudo', dest='pseudo_max', default=PSEUDO_MAX,
                        action='store', type=float, 
                        help='Maximum pseudo angle allowed to be deviate') 

    parser.add_argument('--rmsd', dest='rmsd_accept', default=RMSD_ACCEPT,
                        action='store', type=float, 
                        help='RMSD to accept a fragment') 

    parser.add_argument('--maxrmsd', dest='rmsd_max', default=RMSD_MAX,
                        action='store', type=float, 
                        help='Maximum rmsd for searching a fragment')

    parser.add_argument('--Pexact', dest='flg_Pexact', default=False,
                        action='store_true', 
                        help='Align the position of P exactly')

    parser.add_argument('--basedir', dest='basedir', default=BASEDIR,
                        action='store', 
                        help='Path for the directory where database exist')

    parser.add_argument('--log', dest='logfilename', default='tis2aa.log',
                        action='store', 
                        help='Log filename')
    
    parser.add_argument('--circ', dest='circ', default=False,
                        action='store_true', 
                        help='circular RNA/DNA')
    
    parser.add_argument('pdb_in',  help='Target CG PDB file')
    parser.add_argument('pdb_out', help='Output AA PDB file')

    args = parser.parse_args()

    LIBPDBAA = args.basedir + 'RNA09_FRAG_AA/'  # input
    LIBPDBCG = args.basedir + 'RNA09_FRAG_CG/'  # output
    LISTFILE = args.basedir + 'RNA09.nts'
    
    ################################################################################
    """ Open log file """
    f_log = open(args.logfilename,'w')

    """ Read CG PDB """
    p = PdbFile(args.pdb_in)
    p.open_to_read()
    chains_cg = p.read_all()
    p.close()
    
    
    """ Read pseudo angles from library list (.nts) file """
    #lib_pucker = []
    lib_pseudo = []
    for l in open(LISTFILE):
        lsp = l.split()
        #lib_pucker.append( lsp[8] )
        lib_pseudo.append( (float(lsp[9]), float(lsp[10])) )
    
    nlib = len(lib_pseudo)
    
    chains_aa = []
    for c_cg in chains_cg:
            
        seq = []
        xyzs_P = []
        xyzs_S = []
        xyzs_B = []
    
        """ Store sequence and coordinates of CG model """
        for ir, r_cg in enumerate(c_cg.residues):
    
            flg_P = False; flg_S = False; flg_B = False;
            for a in r_cg.atoms:
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
                f_log.write('no S in cg model, ir=%i' % ir)
            if not flg_B:
                print 'no B in cg model, ir=', ir
                f_log.write('no B in cg model, ir=%i' % ir)
    
            res_name = r_cg.atoms[0].res_name.strip()[-1]
            if res_name in ('A','U','G','C'):
                seq.append( res_name )
            else:
                print 'res_name is not valid:', r_cg.atoms[0].res_name
                f_log.write('Warning: res_name is not valid: %s' % r_cg.atoms[0].res_name)
    
        nres = len(c_cg.residues)
        c_aa = Chain()

        
        """ The first residue """
        if not args.circ:
            xyzSBPS_ref = np.array( [xyzs_S[0], xyzs_B[0], xyzs_P[1], xyzs_S[1]] )
            best_rmsd = 9999.9
            best_lib = 0
            best_mtx = None
            for ilib in range(nlib):
                p = PdbFile('%s/%06i_%s.pdb' % (LIBPDBCG, ilib+1, seq[0]))
                p.open_to_read()
                c = p.read_all()[0]
                xyzSBPS = [ c.residues[1].atoms[1].xyz.get_as_list(), # S2
                            c.residues[1].atoms[2].xyz.get_as_list(), # B2
                            c.residues[2].atoms[0].xyz.get_as_list(), # P3
                            c.residues[2].atoms[1].xyz.get_as_list()] # S3
                rmsd, mat = calcrotation(np.transpose(xyzSBPS_ref), np.transpose(xyzSBPS))
    
                if rmsd < best_rmsd:
                    best_rmsd = rmsd
                    best_lib = ilib + 1
                    best_mtx = mat
                    if rmsd <= args.rmsd_accept:
                        break
    
                p.close()
                f_log.write("#%3i %06i %f\n" % (1, ilib+1, rmsd))
    
            p = PdbFile('%s/%06i_%s.pdb' % (LIBPDBAA, best_lib, seq[0]))
            p.open_to_read()
            c = p.read_all()[0]
            frag_res = c.residues[1]
            p.close()
         
            mtx = mtx_crd_transform()
            mtx.mtx[:,:] = best_mtx
        
            r_aa = Residue()
            for a in frag_res.atoms:
                if a.name.strip() in ("P","OP1","OP2","O1P", "O2P"):  # No P particle at the first nucleotide
                   continue
                a_new = copy.deepcopy(a)
                a_new.xyz.put_as_list( mtx.do_to_array( a.xyz.get_as_list() ) )
                r_aa.push_atom( a_new )
    
            c_aa.push_residue( r_aa )
            f_log.write('%3i %6i  %6i %f\n' % (1, nlib, best_lib, best_rmsd))

    
        """ From the second to the second last residues """
        if args.circ:
            ir_start = 0
            ir_end = nres-1
        else:
            ir_start = 1
            ir_end = nres-2

        for ir in range(ir_start, ir_end+1):
    
            ir0 = ir - 1
            if ir0 < 0:
                ir0 += nres
            ir2 = ir + 1
            if ir2 > nres-1:
                ir2 -= nres

            """ Calcualte pseudo angles """
            xyzS1 = np.array( xyzs_S[ir0] )
            xyzP2 = np.array( xyzs_P[ir] )
            xyzS2 = np.array( xyzs_S[ir] )
            xyzB2 = np.array( xyzs_B[ir] )
            xyzP3 = np.array( xyzs_P[ir2] )
            xyzS3 = np.array( xyzs_S[ir2] )
    
            eta    = torsion(xyzS1, xyzP2, xyzS2, xyzP3, flg_degree=True, flg_360=True)
            theta  = torsion(xyzP2, xyzS2, xyzP3, xyzS3, flg_degree=True, flg_360=True)


            """ For superposition """
            #xyzSPS_ref  = [ xyzS1, xyzP2, xyzS2 ]
            xyzSPSB_ref = [ xyzS1, xyzP2, xyzS2, xyzB2 ]
            
    
            """ Search best fragment """
            ilimit = 0
            flg_found = False
            while not flg_found:
    
                """ Psuedo-angle range """
                ilimit += 1
                limit = args.pseudo_bin * ilimit
    
                if limit > args.pseudo_max:
                    if best_rmsd <= args.rmsd_max:
                        break
                    else:
                        print 'limit > PSEUDO_MAX; could not find library for ir=',ir
                        f_log.write('Error: limit > PSEUDO_MAX; could not find library for ir=%i' % ir)
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
        
                    f_log.write("#%3i %06i %f  %f %f\n" % (ir+1, lib, rmsd, angle[0], angle[1]))
    
                if best_rmsd <= args.rmsd_accept:
                    flg_found = True
    
            f_log.write('%3i %6i  %6i %f  %f %f\n' % (ir+1, len(cand_lib), best_lib, best_rmsd, best_angle[0], best_angle[1]))


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
            if args.flg_Pexact:
                aaP = Coord()
                nP = 0
                for a in r_aa.atoms:
                    if a.name.strip() in ("P","OP1", "OP2", "O1P", "O2P"):
                        aaP = aaP + a.xyz
                        nP += 1
                if nP != 3:
                    print 'Error: nP != 3, ir=,',ir
                    f_log.write('Error: nP != 3, ir=%i' % ir)
                    sys.exit(2)
                aaP = aaP / float(3)
                trans = xyzP2 - aaP.get_as_ndarray()
                mtx.reset()
                mtx.translation( trans[0], trans[1], trans[2] )
                for a in r_aa.atoms:
                    a.xyz.put_as_list( mtx.do_to_array( a.xyz.get_as_list() ) )

            c_aa.push_residue( r_aa )
    

        """ The last residue """
        if not args.circ:
            xyzSPSB_ref = np.array( [xyzs_S[-2], xyzs_P[-1], xyzs_S[-1], xyzs_B[-1]] )
            best_rmsd = 9999.9
            best_lib = 0
            best_mtx = None
            for ilib in range(nlib):
                p = PdbFile('%s/%06i_%s.pdb' % (LIBPDBCG, ilib+1, seq[-1]))
                p.open_to_read()
                c = p.read_all()[0]
                xyzSPSB = [ c.residues[0].atoms[1].xyz.get_as_list(), # S-2
                            c.residues[1].atoms[0].xyz.get_as_list(), # P-1
                            c.residues[1].atoms[1].xyz.get_as_list(), # S-1
                            c.residues[1].atoms[2].xyz.get_as_list()] # B-1
                rmsd, mat = calcrotation(np.transpose(xyzSPSB_ref), np.transpose(xyzSPSB))
    
                if rmsd < best_rmsd:
                    best_rmsd = rmsd
                    best_lib = ilib + 1
                    best_mtx = mat
                    if rmsd <= args.rmsd_accept:
                        break
    
                p.close()
                f_log.write("#%3i %06i %f\n" % (nres, ilib+1, rmsd))
    
            p = PdbFile('%s/%06i_%s.pdb' % (LIBPDBAA, best_lib, seq[-1]))
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
    
            c_aa.push_residue( r_aa )
            f_log.write('%3i %6i  %6i %f\n' % (nres, nlib, best_lib, best_rmsd))
     
        """ Add the chain """
        chains_aa.append( c_aa )
    
    
    """ Re-numbering serial and residue IDs """
    serial = 0
    res_seq = 0
    for c_aa in chains_aa:
        for r_aa in c_aa.residues:
            res_seq += 1
            for a in r_aa.atoms:
                serial += 1
                a.serial = serial
                a.res_seq = res_seq
                res_name = a.res_name.strip()[-1]
                a.res_name = "  %s" % res_name
        c.reindex()
    

    """ Write AA PDB """
    p = PdbFile(args.pdb_out)
    p.open_to_write()
    p.write_all( chains_aa )
    p.close()
    
    f_log.close()
    
