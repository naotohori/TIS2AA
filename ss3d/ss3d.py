#!/usr/bin/env python

import sys
import random
import math

from index import *
from rnaAform import *
from rnaDT15 import *

random.seed(10)

seq = 'UGGGGGUUUUCCCCCU'
ss  = '.(((((....))))).'
#seq = 'UUG'
#ss  = '...'

dt = 0.001
maxt = 1000
nnt = len(ss)

f_xyz = open('ss3d.xyz','w')
f_ene = open('ss3d.ene','w')

#print('SS: %s' % (ss,))
#print('length: %i' % (nnt,))

#RNAP = 0
#RNAS = 1
#RNAB = 2

#BASEA = 0
#BASEU = 1
#BASEG = 2
#BASEC = 3

def random_coord():
    return [random.random(),random.random(),random.random()]

nnt = len(ss)
xyz = []
mp2type = []
for i in range(nnt):
    if i > 0:
        xyz.append(random_coord())
        mp2type.append(RNAP)
    xyz.append(random_coord())
    mp2type.append(RNAS)
    xyz.append(random_coord())
    mp2type.append(RNAB)

nmp = len(xyz) 

#################### Bond ####################
bd_k = []
bd_mp = []
bd_r0 = []

def add_bond(i,j,k,r0):
    bd_mp.append((i,j))
    bd_k.append(k)
    bd_r0.append(r0)

for i in range(nnt):
    if i > 0:
        # SP
        add_bond(3*i-3,3*i-1,DT15.BL_SP,ARNA.BL_SP)
        # PS
        add_bond(3*i-1,3*i,DT15.BL_PS,ARNA.BL_PS)
    # SB
    if seq[i] == 'A':
        add_bond(3*i,3*i+1,DT15.BL_SB,ARNA.BL_SA)
    elif seq[i] == 'U':
        add_bond(3*i,3*i+1,DT15.BL_SB,ARNA.BL_SU)
    elif seq[i] == 'G':
        add_bond(3*i,3*i+1,DT15.BL_SB,ARNA.BL_SG)
    elif seq[i] == 'C':
        add_bond(3*i,3*i+1,DT15.BL_SB,ARNA.BL_SC)

def bond(frc):

    ene = 0.0
    for (i,j), k, r0 in zip(bd_mp, bd_k, bd_r0):
        
        v = [x2 - x1 for x1,x2 in zip(xyz[i],xyz[j])]
        dist = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
        ddist = dist - r0

        # calc force
        f = -2.0 * k * ddist / dist
        frc[i] = [x - f*y for x,y in zip(frc[i],v)] 
        frc[j] = [x + f*y for x,y in zip(frc[j],v)] 

        # calc energy
        ene += k * ddist**2

    return ene

#################### Angle ####################
angl_k = []
angl_mp = []
angl_r0 = []

def add_angl(i,j,k,c,r0):
    angl_mp.append((i,j,k))
    angl_k.append(c)
    angl_r0.append(r0)

for i in range(nnt):
    if i > 0:
        # BSP
        if seq[i-1] == 'A':
            add_angl(3*i-2,3*i-3,3*i-1,DT15.BA_BSP,ARNA.BA_ASP)
        elif seq[i-1] == 'U':
            add_angl(3*i-2,3*i-3,3*i-1,DT15.BA_BSP,ARNA.BA_USP)
        elif seq[i-1] == 'G':
            add_angl(3*i-2,3*i-3,3*i-1,DT15.BA_BSP,ARNA.BA_GSP)
        elif seq[i-1] == 'C':
            add_angl(3*i-2,3*i-3,3*i-1,DT15.BA_BSP,ARNA.BA_CSP)
        # SPS
        add_angl(3*i-3,3*i-1,3*i,DT15.BA_SPS,ARNA.BA_SPS)
        # PSB
        if seq[i] == 'A':
            add_angl(3*i-1,3*i,3*i+1,DT15.BA_PSB,ARNA.BA_PSA)
        elif seq[i] == 'U':
            add_angl(3*i-1,3*i,3*i+1,DT15.BA_PSB,ARNA.BA_PSU)
        elif seq[i] == 'G':
            add_angl(3*i-1,3*i,3*i+1,DT15.BA_PSB,ARNA.BA_PSG)
        elif seq[i] == 'C':
            add_angl(3*i-1,3*i,3*i+1,DT15.BA_PSB,ARNA.BA_PSC)
    if i < nnt-1:
        # PSP
        add_angl(3*i-1,3*i,3*i+2,DT15.BA_PSP,ARNA.BA_PSP)

def angl(frc):

    ene = 0.0
    for (i,j,k), c, r0 in zip(angl_mp, angl_k, angl_r0):
        
        v21 = [x2 - x1 for x1,x2 in zip(xyz[i],xyz[j])]
        v32 = [x3 - x2 for x2,x3 in zip(xyz[j],xyz[k])]

        c11 = v21[0] * v21[0] + v21[1] * v21[1] + v21[2] * v21[2]
        c22 = v32[0] * v32[0] + v32[1] * v32[1] + v32[2] * v32[2]
        c21 = v32[0] * v21[0] + v32[1] * v21[1] + v32[2] * v21[2]

        co_theta = - c21 / math.sqrt(c11 * c22)
        if co_theta > 1.0:
            co_theta = 1.0
        elif co_theta < -1.0:
            co_theta = -1.0

        dba = math.acos(co_theta) - r0

        t3 = c11 * c22 - c21**2
        if t3 <= 1.0:
            t3 = 1.0

        # calc force
        f = 2.0 * c * dba / math.sqrt(t3) 

        f21 = [f * (v1 * c21 / c11 - v2) for v1,v2 in zip(v21,v32)]
        f32 = [f * (v1 * c21 / c22 - v2) for v1,v2 in zip(v32,v21)]

        frc[i] = [x - f1      for x,f1    in zip(frc[i],f21)]
        frc[j] = [x + f1 - f2 for x,f1,f2 in zip(frc[j],f21,f32)]
        frc[k] = [x + f2      for x,f2    in zip(frc[k],f32)]

        # calc energy
        ene += c * dba**2

    return ene

#################### Base pair ####################
bp_mp = []
bp_k = []
bp_r0 = []
def add_bp(i,j,k,r0):
    bp_mp.append((i,j))
    bp_k.append(k)
    bp_r0.append(r0)

#bp_pairs = []
q = []
for i,s in enumerate(ss):
    if s == '(':
        q.append(i)
    elif s == ')':
        j = q.pop()
        #bp_pairs.append((j,i))
        add_bp(3*i+1,3*j+1,10.0,5.0)

if len(q) != 0:
    print ('Error: len(q) != 0')
    sys.exit(2)

def bp_spring(frc):

    ene = 0.0
    for (i,j), k, r0 in zip(bp_mp, bp_k, bp_r0):
        
        v = [x2 - x1 for x1,x2 in zip(xyz[i],xyz[j])]
        dist = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
        ddist = dist - r0

        # calc force
        f = -2.0 * k * ddist / dist
        frc[i] = [x - f*y for x,y in zip(frc[i],v)] 
        frc[j] = [x + f*y for x,y in zip(frc[j],v)] 

        # calc energy
        ene += k * ddist**2

    return ene

##########################################################

def write_xyz(t,xyz):
    f_xyz.write ('%i\n' % (nmp,))
    f_xyz.write ('# %i\n' % (t,))
    for i,x in enumerate(xyz):
        if (mp2type[i] == RNAP):
            f_xyz.write ('P')
        elif (mp2type[i] == RNAS):
            f_xyz.write ('S')
        elif (mp2type[i] == RNAB):
            f_xyz.write ('B')
        f_xyz.write (' %f %f %f\n' % (x[0],x[1],x[2]))

def write_ene(t,Ebond,Ebp,Eangl):
    ene = Ebond + Ebp + Eangl
    f_ene.write('%10i %f %f %f %f\n' % (t,ene, Ebp, Ebond, Eangl))

def move_to_origin(xyz):
    orig = [0.0]*3
    for x in xyz:
        orig[0] += x[0]
        orig[1] += x[1]
        orig[2] += x[2]
    fl = float(len(xyz))
    orig[0] /= fl
    orig[1] /= fl
    orig[2] /= fl
    for x in xyz:
        x[0] -= orig[0]
        x[1] -= orig[1]
        x[2] -= orig[2]

ene_pre = 0.0
for t in range(maxt):

    frc = [[0.0,0.0,0.0] for i in range(nmp)]
    Ebond = bond(frc)
    Ebp   = bp_spring(frc)
    Eangl = angl(frc) 
    for x,f in zip(xyz, frc):
        x[0] = x[0] + dt * f[0]
        x[1] = x[1] + dt * f[1]
        x[2] = x[2] + dt * f[2]

    move_to_origin(xyz)
    write_xyz(t,xyz)
    write_ene(t,Ebond,Ebp,Eangl)

    ene = Ebond + Eangl + Ebp
    if (abs(ene - ene_pre) < 0.000001):
        break
    ene_pre = ene

f_xyz.close()
f_ene.close()
