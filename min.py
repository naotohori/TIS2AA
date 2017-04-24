#!/usr/bin/env python

import subprocess
import argparse
import shutil

parser = argparse.ArgumentParser(description='Minimize given structure using sander (jsut a wrapper for sander!)')

parser.add_argument('--noH', dest='flg_noH', default=False,
                    action='store_true', help='Hydrogen atoms are removed from output')
parser.add_argument('--rest', dest='restraint', action='store',
                             default=50.0, help='restraint_wt in sander') 
parser.add_argument('--ncyc', dest='ncyc', action='store',
                             default=200, help='ncyc in sander') 
parser.add_argument('--maxcyc', dest='maxcyc', action='store',
                             default=500, help='maxcyc in sander') 
parser.add_argument('--prefix', dest='prefix', action='store',
                            default='min', help='prefix for output files')

parser.add_argument('pdb_in', help='target PDB file')
parser.add_argument('pdb_out', help='output minimized PDB file')

args = parser.parse_args()

flg_noH = args.flg_noH
restraint = args.restraint
ncyc = args.ncyc
maxcyc = args.maxcyc
prefix = args.prefix

file_pdb = args.pdb_in
file_minpdb = args.pdb_out

file_min = prefix + '.in'
file_out = prefix + '.out'
file_leap = prefix + '.leap.in'
file_prmtop = prefix + '.prmtop'
file_inpcrd = prefix + '.inpcrd'
file_rst = prefix + '.rst'

################################################################################

'''
Leap
'''
f_out = open(file_leap, 'w')
f_out.write('source leaprc.RNA.OL3\n')
f_out.write('x = loadpdb %s\n' % file_pdb)
f_out.write('saveamberparm x %s %s\n' % (file_prmtop, file_inpcrd))
f_out.write('quit\n')
f_out.close()

#tleap -f leap.in
subprocess.call( ['tleap','-f',file_leap] )


'''
sander
'''
f_out = open(file_min,'w')
f_out.write("initial minimization whole system\n")
f_out.write("&cntrl\n")
f_out.write("  imin   = 1\n")
f_out.write("  ncyc = %i, maxcyc = %i\n" % (ncyc, maxcyc))
f_out.write("  ntb    = 0\n")
f_out.write("  cut    = 999.0\n")
f_out.write("  ntr = 1\n")
f_out.write("  restraint_wt = %f\n" % restraint)
f_out.write('  restraintmask="@P"\n')
f_out.write("/\n")
f_out.close()
#sander -O -i min.in -o mini.out -p prmtop -c inpcrd -ref inpcrd -r mini.rst
subprocess.call( ['sander','-O', '-i',file_min, '-o', file_out, '-p', file_prmtop,
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
if flg_noH:

    pdb = open(file_minpdb+'.noH','w') 

    serial = 0
    for l in open(file_minpdb):
        if l[0:4] == 'ATOM':
            if l[12:16].strip()[0] == 'H':
                continue
            else:
                serial += 1
                l = l[0:6] + '%5i'%serial + l[12:]
                pdb.write(l)

    pdb.close()
    shutil.move( file_minpdb+'.noH', file_minpdb)
