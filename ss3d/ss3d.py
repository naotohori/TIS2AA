#!/usr/bin/env python

import sys
import random
import math

from index import *
from rnaAform import *
from rnaDT15 import *

random.seed(11)

seq = 'UGGGGGUUUUCCCCCU'
ss  = '.(((((....))))).'
#seq = 'UUG'
#ss  = '...'

dt = 0.005
maxt = 10000
nnt = len(ss)
tempk = 300.0
kB = 0.00198720359

f_xyz = open('ss3d.xyz','w')
f_ene = open('ss3d.ene','w')
f_ene.write('# t  Etotal  Ebond  Eangl  Ebp  Estack  Eexv\n')

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
    #return [random.random(),random.random(),random.random()]
    return [random.uniform(-100.0,100.0),random.uniform(-100.0,100.0),random.uniform(-100.0,100.0)]

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
    angl_r0.append(r0 / 180.0 * math.pi)

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
    
#################### Base stack ####################
stack_mp = []
stack_k = []
stack_r0 = []
def add_stack(i1,i2,i3,i4,i5,i6,i7,k0,k1,k2,k3,r01,r02,r03):
    stack_mp.append((i1,i2,i3,i4,i5,i6,i7))
    stack_k.append((k0,k1,k2,k3))
    stack_r0.append((r01,r02/180.0*math.pi,r03/180.0*math.pi))
    print (i1,i2,k0,r01)

for jnt in range(1,nnt-2):
    imp1 = jnt*3 + 1
    imp2 = jnt*3 + 4
    imp3 = jnt*3 - 1
    imp4 = jnt*3
    imp5 = jnt*3 + 2
    imp6 = jnt*3 + 3
    imp7 = jnt*3 + 5

    if seq[jnt] == 'A':
        if seq[jnt+1] == 'A':
            native = ARNA.ST_AA ; type_str = 'A-A'
        elif seq[jnt+1] == 'U':
            native = ARNA.ST_AU ; type_str = 'A-U'
        elif seq[jnt+1] == 'G':
            native = ARNA.ST_AG ; type_str = 'A-G'
        elif seq[jnt+1] == 'C':
            native = ARNA.ST_AC ; type_str = 'A-C'
    elif seq[jnt] == 'U':
        if seq[jnt+1] == 'A':
            native = ARNA.ST_UA ; type_str = 'U-A'
        elif seq[jnt+1] == 'U':
            native = ARNA.ST_UU ; type_str = 'U-U'
        elif seq[jnt+1] == 'G':
            native = ARNA.ST_UG ; type_str = 'U-G'
        elif seq[jnt+1] == 'C':
            native = ARNA.ST_UC ; type_str = 'U-C'
    elif seq[jnt] == 'G':
        if seq[jnt+1] == 'A':
            native = ARNA.ST_GA ; type_str = 'G-A'
        elif seq[jnt+1] == 'U':
            native = ARNA.ST_GU ; type_str = 'G-U'
        elif seq[jnt+1] == 'G':
            native = ARNA.ST_GG ; type_str = 'G-G'
        elif seq[jnt+1] == 'C':
            native = ARNA.ST_GC ; type_str = 'G-C'
    elif seq[jnt] == 'C':
        if seq[jnt+1] == 'A':
            native = ARNA.ST_CA ; type_str = 'C-A'
        elif seq[jnt+1] == 'U':
            native = ARNA.ST_CU ; type_str = 'C-U'
        elif seq[jnt+1] == 'G':
            native = ARNA.ST_CG ; type_str = 'C-G'
        elif seq[jnt+1] == 'C':
            native = ARNA.ST_CC ; type_str = 'C-C'
    
    (h,s,Tm) = DT15.STACK_PARAM[type_str]
    u0 = -h + kB * (tempk - Tm) * s

    add_stack(imp1,imp2,imp3,imp4,imp5,imp6,imp7,
              u0, DT15.ST_DIST, DT15.ST_DIH, DT15.ST_DIH,
              native, ARNA.DIH_PSPS, ARNA.DIH_SPSP)

def stack(frc):

    ene = 0.0

    for (i1,i2,i3,i4,i5,i6,i7), (k0,k1,k2,k3), (r01,r02,r03) in zip(stack_mp, stack_k, stack_r0):

        ediv = 1.0
        fv = [[0.0,0.0,0.0] for i in range(8)]  # fv[0] is dummy
     
        #===== calc vectors =====
        v21 = [x2 - x1 for x1,x2 in zip(xyz[i1],xyz[i2])]
        v34 = [x3 - x4 for x3,x4 in zip(xyz[i3],xyz[i4])]
        v54 = [x5 - x4 for x4,x5 in zip(xyz[i4],xyz[i5])]
        v56 = [x5 - x6 for x5,x6 in zip(xyz[i5],xyz[i6])]
        v76 = [x7 - x6 for x6,x7 in zip(xyz[i6],xyz[i7])]

        #===== 1. Distance between 1 and 2 =====
        dist = math.sqrt(v21[0]**2 + v21[1]**2 + v21[2]**2)
        ddist = dist - r01
        ediv += k1 * ddist**2

        f = 2.0 * k1 * ddist / dist
        fv[1] = [- f * x  for x in v21]
        fv[2] = [  f * x  for x in v21]

        #===== 2. Dihedral angle among 3-4-5-6 =====
        m = [0.0,0.0,0.0]
        n = [0.0,0.0,0.0]
        m[0] = v34[1]*v54[2] - v34[2]*v54[1]
        m[1] = v34[2]*v54[0] - v34[0]*v54[2]
        m[2] = v34[0]*v54[1] - v34[1]*v54[0]
        n[0] = v54[1]*v56[2] - v54[2]*v56[1]
        n[1] = v54[2]*v56[0] - v54[0]*v56[2]
        n[2] = v54[0]*v56[1] - v54[1]*v56[0]

        dmm = m[0]**2 + m[1]**2 + m[2]**2
        dnn = n[0]**2 + n[1]**2 + n[2]**2
        d5454 = v54[0]**2 + v54[1]**2 + v54[2]**2
        abs54 = math.sqrt(d5454)
        d5654 = sum([x*y for x,y in zip(v56,v54)])
        d3454over5454 = sum([x*y for x,y in zip(v34,v54)]) / d5454
        d5654over5454 = d5654 / d5454

        dih = math.atan2(sum([x*y for x,y in zip(v34,n)])*math.sqrt(d5454), sum([x*y for x,y in zip(m,n)]))
        d = dih - r02
        if d > math.pi:
            d = d - 2*math.pi
        elif d < -math.pi:
            d = d + 2*math.pi

        ediv += k2 * d**2

        f_i = [+ 2.0 * k2 * d * abs54 / dmm * x for x in m]
        f_l = [- 2.0 * k2 * d * abs54 / dnn * x for x in n]

        fv[3] = f_i
        fv[4] = [(-1.0 + d3454over5454) * x - d5654over5454 * y for x,y in zip(f_i,f_l)]
        fv[5] = [(-1.0 + d5654over5454) * y - d3454over5454 * x for x,y in zip(f_i,f_l)]
        fv[6] = f_l

        #===== 3. Dihedral angle among 7-6-5-4 =====
        m[0] = v76[1]*v56[2] - v76[2]*v56[1]
        m[1] = v76[2]*v56[0] - v76[0]*v56[2]
        m[2] = v76[0]*v56[1] - v76[1]*v56[0]
        #n[0] = v56[1]*v54[2] - v56[2]*v54[1]
        #n[1] = v56[2]*v54[0] - v56[0]*v54[2]
        #n[2] = v56[0]*v54[1] - v56[1]*v54[0]
        n = [-x for x in n]

        dmm = m[0]**2 + m[1]**2 + m[2]**2
        # dnn does not change.
        d5656 = v56[0]**2 + v56[1]**2 + v56[2]**2
        abs56 = math.sqrt(d5656)
        d7656over5656 = sum([x*y for x,y in zip(v76,v56)]) / d5656
        d5456over5656 = d5654 / d5656

        dih = math.atan2(sum([x*y for x,y in zip(v76,n)])*math.sqrt(d5656), sum([x*y for x,y in zip(m,n)]))
        d = dih - r03
        if d > math.pi:
            d = d - 2*math.pi
        elif d < -math.pi:
            d = d + 2*math.pi

        ediv += k3 * d**2

        f_i = [+ 2.0 * k3 * d * abs56 / dmm * x for x in m]
        f_l = [- 2.0 * k3 * d * abs56 / dnn * x for x in n]
        fv[7] = [f0 + x for f0,x in zip(fv[7],f_i)]
        fv[6] = [f0 + (-1.0 + d7656over5656) * x - d5456over5656 * y for f0,x,y in zip(fv[6],f_i,f_l)]
        fv[5] = [f0 + (-1.0 + d5456over5656) * y - d7656over5656 * x for f0,x,y in zip(fv[5],f_i,f_l)]
        fv[4] = [f0 + x for f0,x in zip(fv[4],f_l)]

        #===== Total =====
        fv = [[k0 / ediv**2 * x for x in f0] for f0 in fv]

        frc[i1] = [x + y for x,y in zip(frc[i1],fv[1])]
        frc[i2] = [x + y for x,y in zip(frc[i2],fv[2])]
        frc[i3] = [x + y for x,y in zip(frc[i3],fv[3])]
        frc[i4] = [x + y for x,y in zip(frc[i4],fv[4])]
        frc[i5] = [x + y for x,y in zip(frc[i5],fv[5])]
        frc[i6] = [x + y for x,y in zip(frc[i6],fv[6])]
        frc[i7] = [x + y for x,y in zip(frc[i7],fv[7])]

        ene += k0 / ediv

    return ene

#################### Excluded Volume ####################
        
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

def exv(frc, factor=1.0):

    a = 1.5852
    a2 = a*a
    dij = 3.2
    coef = 0.2 * 0.2
    DE_MAX = 1000.0

    ene = 0.0

    for i in range(nmp-1):
        for j in range(i+1,nmp):

            v = [x2 - x1 for x1,x2 in zip(xyz[i],xyz[j])]
            dist = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)

            if dist > dij:
                continue

            dr = dist + a - dij
            dr2 = dr * dr

            roverdist2 = a2 / dr2
            roverdist4 = roverdist2 * roverdist2
            roverdist6 = roverdist4 * roverdist2
            roverdist8 = roverdist4 * roverdist4
            roverdist12 = roverdist8 * roverdist4
            roverdist14 = roverdist2 * roverdist4 * roverdist8

            # calc force
            f = factor * abs(12.0 * coef * (roverdist14 - roverdist8) * dr / a2 / dist)
            if f > DE_MAX:
                f = DE_MAX

            frc[i] = [x - f*y for x,y in zip(frc[i],v)] 
            frc[j] = [x + f*y for x,y in zip(frc[j],v)] 

            # calc energy
            ene += factor * coef * (roverdist12 - 2*roverdist6 + 1.0)

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

def write_ene(t,Ebond,Eangl,Ebp,Estack,Eexv):
    ene = Ebond + Eangl + Ebp + Estack + Eexv
    f_ene.write('%10i %f %f %f %f %f %f\n' % (t,ene, Ebp, Ebond, Eangl, Estack, Eexv))

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


#flg_first = True
#if flg_first:
#    frc = [[0.0,0.0,0.0] for i in range(nmp)]
#    Ebond = bond(frc)
#    Ebp   = bp_spring(frc)
#    Eangl = angl(frc) 
#
#    norm_max = max(max(frc),abs(min(frc)))
#    print (norm_max)

ene_pre = 0.0
write_xyz(0,xyz)
for t in range(1,maxt+1):

    frc = [[0.0,0.0,0.0] for i in range(nmp)]
    Ebond = bond(frc)
    Eangl = angl(frc) 
    Ebp   = bp_spring(frc)
    Estack= stack(frc)
    #Eexv  = exv(frc,1.0)
    Eexv  = 0.0

    for x,f in zip(xyz, frc):
        x[0] = x[0] + dt * f[0]
        x[1] = x[1] + dt * f[1]
        x[2] = x[2] + dt * f[2]

    move_to_origin(xyz)
    write_xyz(t,xyz)
    write_ene(t-1,Ebond,Eangl,Ebp,Estack,Eexv)

    ene = Ebond + Eangl + Ebp
    if (abs(ene - ene_pre) < 0.000001):
        break
    ene_pre = ene

frc = [[0.0,0.0,0.0] for i in range(nmp)]
Ebond = bond(frc)
Ebp   = bp_spring(frc)
Eangl = angl(frc)
Eexv  = exv(frc,1.0)
write_ene(t,Ebond,Eangl,Ebp,Estack,Eexv)

f_xyz.close()
f_ene.close()
