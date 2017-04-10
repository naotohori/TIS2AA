#!/usr/bin/env python

import sys
import subprocess

if len(sys.argv) != 2:
    print 'Usage: min.py [name]'
    sys.exit(2)

name = sys.argv[1]

f_out = open('%s.leap.in' % name, 'w')
for l in open('leap.in.temp'):
    l.replace('##NAME##', name)
    f_out.write(l)
f_out.close()

#tleap -f leap.in
subprocess.call( ['tleap','-f','%s.leap.in' % name] )

#sander -O -i min.in -o mini.out -p prmtop -c inpcrd -ref inpcrd -r mini.rst
subprocess.call( ['sander','-O', '-i','min.in', '-o', '%s.min.out' % name, '-p', '%s.prmtop' % name,
                  '-c', '%s.inpcrd' % name, '-ref', '%s.inpcrd' % name, '-r', '%s.rst' % name] )

#ambpdb -p prmtop -c mini.rst > VPK_0_34.aa.min.pdb
f_minpdb = open('%s.aa.min.pdb' % name, 'w')
subprocess.call( ['ambpdb', '-p', '%s.prmtop' % name, '-c', '%s.rst' % name], stdout=f_minpdb )
f_minpdb.close()
