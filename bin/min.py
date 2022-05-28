#!/usr/bin/env python
''' Wrapper of tleap & sander '''
__author__ = "Naoto Hori"

import os
import subprocess
import argparse
import shutil

parser = argparse.ArgumentParser(description='Minimize given structure using sander (jsut a wrapper for sander!)',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--noH', dest='flg_noH', default=False,
                    action='store_true', help='Remove hydrogen atoms from output')
parser.add_argument('--rest', dest='restraint', default=50.0,
                    action='store', type=float, 
                    help='[sander parameter] restraint_wt ; Put a negative value if no restraint desired.') 
parser.add_argument('--maxcyc', dest='maxcyc', default=500,
                    action='store', type=int, help='[sander parameter] maxcyc') 
parser.add_argument('--ncyc', dest='ncyc', default=200,
                    action='store', type=int, help='[sander parameter] ncyc') 
parser.add_argument('--prefix', dest='prefix', default='min',
                    action='store', help='prefix for sander output files')
parser.add_argument('--circ', dest='circ', default=False,
                    action='store_true', help='circular RNA/DNA')

parser.add_argument('pdb_in',  help='target PDB file')
parser.add_argument('pdb_out', help='output minimized PDB file')

args = parser.parse_args()

################################################################################

file_pdb = args.pdb_in
file_minpdb = args.pdb_out

if args.circ:
    file_leap_src = 'leaprc.circRNA.OL3'
else:
    file_leap_src = 'leaprc.RNA.OL3'

file_min = args.prefix + '.in'
file_out = args.prefix + '.out'
file_mdinfo = args.prefix + '.mdinfo'
file_leap = args.prefix + '.leap.in'
file_prmtop = args.prefix + '.prmtop'
file_inpcrd = args.prefix + '.inpcrd'
file_rst = args.prefix + '.rst'

################################################################################
'''
Leap
'''
""" Get the first and last residue ID if circ """
if args.circ:
    ir_ini = 1
    ll = ''
    for l in open(file_pdb,'r'):
        if l[0:6] == 'ATOM  ':
            ll = l
    ir_end = int(ll[22:26])

f_out = open(file_leap, 'w')
f_out.write('source %s\n' % file_leap_src)
f_out.write('x = loadpdb %s\n' % file_pdb)
if args.circ:
    f_out.write("bond x.%i.O3' x.%i.P\n" % (ir_end, ir_ini))
f_out.write('saveamberparm x %s %s\n' % (file_prmtop, file_inpcrd))
f_out.write('quit\n')
f_out.close()

#tleap -f leap.in
FNULL = open(os.devnull, 'w')
subprocess.call( ['tleap','-f',file_leap], stdout=FNULL )


'''
sander
'''
f_out = open(file_min,'w')
f_out.write("initial minimization whole system\n")
f_out.write("&cntrl\n")
f_out.write("  imin   = 1\n")
f_out.write("  ncyc = %i, maxcyc = %i\n" % (args.ncyc, args.maxcyc))
f_out.write("  ntb    = 0\n")
f_out.write("  cut    = 999.0\n")
if not args.restraint < 0.0:
    f_out.write("  ntr = 1\n")
    f_out.write("  restraint_wt = %f\n" % args.restraint)
    f_out.write('  restraintmask="@P"\n')
f_out.write("/\n")
f_out.close()
#sander -O -i min.in -o mini.out -p prmtop -c inpcrd -ref inpcrd -r mini.rst
subprocess.call( ['sander','-O', '-i',file_min, '-o', file_out, '-inf', file_mdinfo, '-p', file_prmtop,
                  '-c', file_inpcrd, '-ref', file_inpcrd, '-r', file_rst] )


'''
ambpdb
'''
#ambpdb -p prmtop -c mini.rst > VPK_0_34.aa.min.pdb
f_minpdb = open(file_minpdb, 'w')
subprocess.call( ['ambpdb', '-p', file_prmtop, '-c', file_rst], stdout=f_minpdb )
f_minpdb.close()


'''
Remove hydrogens (optional)
'''
if args.flg_noH:

    pdb = open(file_minpdb+'.noH','w') 

    serial = 0
    for l in open(file_minpdb):
        if l[0:4] == 'ATOM':
            if l[12:16].strip()[0] == 'H':
                continue
            else:
                serial += 1
                l = l[0:6] + '%5i'%serial + l[11:]
                pdb.write(l)

    pdb.close()
    shutil.move( file_minpdb+'.noH', file_minpdb)
