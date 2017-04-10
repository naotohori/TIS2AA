#!/usr/bin/env python

import math
import numpy as np

def torsion(p1, p2, p3, p4, flg_degree=False, flg_360=False):
    '''
    p1, p2, p3, p4 are numpy arrays for four coordinates by which the torsion angle phi is defined.
    '''
    
    v12 = p1 - p2
    v32 = p3 - p2
    v34 = p3 - p4

    m = np.cross(v12, v32)
    n = np.cross(v32, v34)

    if np.dot(v12,n) >= 0.0: # 0 < phi < 180
        phi =   math.acos( np.dot(m,n) / math.sqrt(np.dot(m,m)*np.dot(n,n)) )

    else: # -180 < phi < 0
        phi = - math.acos( np.dot(m,n) / math.sqrt(np.dot(m,m)*np.dot(n,n)) )

        if flg_360: # convert to (180 < phi < 360)
            phi = 2 * math.pi + phi

    if flg_degree:
        return phi * 180.0 / math.pi
    else:
        return phi

def torsion_xy(p1, p2, p3, p4):
    '''
    p1, p2, p3, p4 are numpy arrays for four coordinates by which the torsion angle phi is defined.
    '''
    
    v12 = p1 - p2
    v32 = p3 - p2
    v34 = p3 - p4

    m = np.cross(v12, v32)
    n = np.cross(v32, v34)

    if np.dot(v12,n) >= 0.0: # 0 < phi < 180
        phi =   math.acos( np.dot(m,n) / math.sqrt(np.dot(m,m)*np.dot(n,n)) )

    else: # -180 < phi < 0
        phi = - math.acos( np.dot(m,n) / math.sqrt(np.dot(m,m)*np.dot(n,n)) )

    return math.cos(phi), math.sin(phi)


if __name__ == "__main__":
    print 'Test of torsion.py'
    print ('''
    ATOM  22440  C4'   G A   0      24.717  92.775  96.702  1.00 31.19           C
    ATOM  22458  P     A A   1      22.332  94.740  95.607  1.00 34.55           P
    ATOM  22463  C4'   A A   1      20.808  97.714  93.719  1.00 34.31           C
    ATOM  22480  P     G A   2      18.954  96.788  90.471  1.00 34.67           P
    ATOM  22485  C4'   G A   2      19.175  99.605  87.735  1.00 33.40           C
    ''')

    a = np.array([ 24.717 , 92.775 , 96.702])
    b = np.array([ 22.332 , 94.740 , 95.607])
    c = np.array([ 20.808 , 97.714 , 93.719])
    d = np.array([ 18.954 , 96.788 , 90.471])
    e = np.array([ 19.175 , 99.605 , 87.735])

    print "Torsion of C4'-P-C4'-P"
    print '  flg_360=False(default): %f rad, %f deg' % (torsion(a,b,c,d), torsion(a,b,c,d,True))
    print '  flg_360=True          : %f rad, %f deg' % (torsion(a,b,c,d,flg_360=True), torsion(a,b,c,d,True,flg_360=True))
    print '  cos(phi), sin(phi)    : x = %f, y = %f' % torsion_xy(a,b,c,d)

    print "Torsion of P-C4'-P-C4'"
    print '  flg_360=False(default): %f rad, %f deg' % (torsion(b,c,d,e), torsion(b,c,d,e,True))
    print '  flg_360=True          : %f rad, %f deg' % (torsion(b,c,d,e,flg_360=True), torsion(b,c,d,e,True,flg_360=True))
    print '  cos(phi), sin(phi)    : x = %f, y = %f' % torsion_xy(b,c,d,e)

