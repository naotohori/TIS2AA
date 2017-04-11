#!/usr/bin/env python

import sys
import subprocess

if len(sys.argv) != 2:
    print 'Usage: min.py [name]'
    sys.exit(2)

name = sys.argv[1]

f_out = open('%s.leap.in' % name, 'w')
f_out.write('source leaprc.RNA.OL3\n')
f_out.write('x = loadpdb %s.aa.pdb\n' % name)
f_out.write('saveamberparm x %s.prmtop %s.inpcrd\n' % (name,name))
f_out.write('quit\n')
#for l in open('leap.in.temp'):
#    f_out.write( l.replace('##NAME##', name) )
f_out.close()

#tleap -f leap.in
subprocess.call( ['tleap','-f','%s.leap.in' % name] )

#sander -O -i min.in -o mini.out -p prmtop -c inpcrd -ref inpcrd -r mini.rst
f_out = open('min.in','w')
f_out.write("initial minimization whole system\n")
f_out.write("&cntrl\n")
f_out.write("  imin   = 1\n")
f_out.write("  ncyc = 500, maxcyc = 2000\n")
f_out.write("  ntb    = 0\n")
f_out.write("  cut    = 999.0\n")
f_out.write("  ntr = 1\n")
f_out.write("  restraint_wt = 5.0\n")
f_out.write('  restraintmask="@P"\n')
f_out.write("/\n")
f_out.close()
subprocess.call( ['sander','-O', '-i','min.in', '-o', '%s.min.out' % name, '-p', '%s.prmtop' % name,
                  '-c', '%s.inpcrd' % name, '-ref', '%s.inpcrd' % name, '-r', '%s.rst' % name] )

#ambpdb -p prmtop -c mini.rst > VPK_0_34.aa.min.pdb
f_minpdb = open('%s.aa.min.pdb' % name, 'w')
subprocess.call( ['ambpdb', '-p', '%s.prmtop' % name, '-c', '%s.rst' % name], stdout=f_minpdb )
f_minpdb.close()
