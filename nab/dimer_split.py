#!/usr/bin/env python

import re

files = [ 'aa_uu.pdb', 'cg_cg.pdb', 'ac_gu.pdb', 'ga_uc.pdb', 'gc_gc.pdb',
          'ua_ua.pdb', 'ug_ca.pdb', 'gg_cc.pdb', 'ag_cu.pdb', 'au_au.pdb',]

filename_re = re.compile('(\S\S)_(\S\S).pdb')

for f in files:
    r = filename_re.match(f)
    dimer1 = r.group(1)
    dimer2 = r.group(2)

    f1 = open('%s.pdb' % dimer1, 'w')
    if dimer1 == dimer2:
        f2 = None
    else:
        f2 = open('%s.pdb' % dimer2, 'w')

    flg_1 = True
    for l in open(f):
        if flg_1:
            f1.write(l)
            if l[0:3] == 'TER':
                flg_1 = False
        elif f2 is not None:
            f2.write(l)

    f1.close()
    if f2 is not None:
        f2.close()

