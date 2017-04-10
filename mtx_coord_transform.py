#!/usr/bin/env python
# coding: UTF-8
'''
Created on 2013/06/28 (based on matrix_transform_make.py)
@author: Naoto Hori
'''

from numpy import identity, dot
from math import cos, sin

class mtx_crd_transform():
    def __init__(self):
        self.mtx = identity(4)
        
    def reset(self):
        self.mtx = identity(4)
        
    def show(self):
        for i in range(4):
            print tuple(self.mtx[i])
            
    def do_to_array(self, d):
        '''配列d[x,y,z]を受け取って、mtxを施して返す。
          dの値は変更されない。
          使い方：
            d = mtx_crd_transform.do_to_array(d)'''
        return dot(self.mtx, d+[1.0,])[0:3]
    
    def do_to_data(self,d):
        for i,v in enumerate(d):
            d[i][0:3] = dot(self.mtx, v+[1.0,])[0:3]
        
    def translation(self,x,y,z):
        '''並進移動'''
        self.mtx[0,3] += x
        self.mtx[1,3] += y
        self.mtx[2,3] += z

    def rotate_by_mtx(self, mtx_rot):
        self.mtx = dot(mtx_rot, self.mtx)
    
    def rotate(self, nx, ny, nz, t):
        '''任意の単位ベクトル(nx,ny,nz)を軸として、角度tだけ回転'''
        ope = identity(4) # operation matrix
        
        ope[0,0] = nx*nx*(1.0-cos(t)) +    cos(t)
        ope[1,0] = nx*ny*(1.0-cos(t)) + nz*sin(t)
        ope[2,0] = nz*nx*(1.0-cos(t)) - ny*sin(t)
    
        ope[0,1] = nx*ny*(1.0-cos(t)) - nz*sin(t)
        ope[1,1] = ny*ny*(1.0-cos(t)) +    cos(t)
        ope[2,1] = ny*nz*(1.0-cos(t)) + nx*sin(t)
    
        ope[0,2] = nz*nx*(1.0-cos(t)) + ny*sin(t)
        ope[1,2] = ny*nz*(1.0-cos(t)) - nx*sin(t)
        ope[2,2] = nz*nz*(1.0-cos(t)) +    cos(t)
        
        self.mtx = dot(ope, self.mtx)
     
    def rotate_x(self, t):
        '''X軸まわりに、角度tだけ回転'''
        ope = identity(4) # operation matrix
        
        ope[1,1] = cos(t)
        ope[2,1] = sin(t)
    
        ope[1,2] = -sin(t)
        ope[2,2] = cos(t)
        
        self.mtx = dot(ope, self.mtx)
        
    def rotate_y(self, t):
        '''Y軸まわりに、角度tだけ回転'''
        ope = identity(4) # operation matrix
        
        ope[0,0] = cos(t)
        ope[2,0] = -sin(t)
    
        ope[0,2] = sin(t)
        ope[2,2] = cos(t)
        
        self.mtx = dot(ope, self.mtx)
        
    def rotate_z(self, t):
        '''Z軸まわりに、角度tだけ回転'''
        ope = identity(4) # operation matrix
        
        ope[0,0] = cos(t)
        ope[1,0] = sin(t)
        
        ope[0,1] = -sin(t)
        ope[1,1] = cos(t)
        
        self.mtx = dot(ope, self.mtx)
         
    def euler_zxz(self,a,b,c):
        '''Z-X-Z系のオイラー角で回転'''
        '''これは、rotate_z, rotate_x, rotate_zを連続で呼び出すのと同じ'''
        ope = identity(4) # operation matrix
        
        ope[0,0] = cos(a)*cos(c) - sin(a)*cos(b)*sin(c)
        ope[1,0] = cos(a)*sin(c) + sin(a)*cos(b)*cos(c)
        ope[2,0] = sin(a)*sin(b)
    
        ope[0,1] = - sin(a)*cos(c) - cos(a)*cos(b)*sin(c)
        ope[1,1] = - sin(a)*sin(c) + cos(a)*cos(b)*cos(c)
        ope[2,1] = cos(a)*sin(b)
    
        ope[0,2] = sin(b)*sin(c)
        ope[1,2] = - sin(b)*cos(c)
        ope[2,2] = cos(b)
        
        self.mtx = dot(ope, self.mtx)
         
if __name__ == "__main__" :
    import sys
    
    if not len(sys.argv) in (7,8):
        print ''
        print 'This script makes a homogeneous transformation matrix,'
        print 'angles of which is defined by Z-X-Z Euler angles.'
        print ''
        print 'Usage: % SCRIPT [alpha] [beta] [gamma] [x] [y] [z] [[output]]'
        print ''
        print 'When "output" is specified, the matrix will be written in the file.'
        print 'Otherwise STDOUT is used to display.'
        sys.exit(2)

    a = float(sys.argv[1])
    b = float(sys.argv[2])
    c = float(sys.argv[3])
    x = float(sys.argv[4])
    y = float(sys.argv[5])
    z = float(sys.argv[6])
    
    mtx = mtx_crd_transform()
    mtx.translation(x, y, z)
    mtx.euler_zxz(a, b, c)
    
    if len(sys.argv) == 8: # Output to a file
        file_mat = file(sys.argv[-1],'w')
        file_mat.write('#matrix\n')
        for i in range(4):
            file_mat.write('%15.10f %15.10f %15.10f %15.10f\n' % tuple(mtx.mtx[i]))
        file_mat.write('#a: %f\n#b: %f\n#c: %f\n' % (a,b,c)) 
        file_mat.write('#x: %f\n#y: %f\n#z: %f\n' % (x,y,z)) 
        file_mat.close()
    else:  # Display on STDOUT
        sys.stdout.write('#matrix\n')
        for i in range(4):
            sys.stdout.write('%15.10f %15.10f %15.10f %15.10f\n' % tuple(mtx.mtx[i]))
        sys.stdout.write('\n')
        sys.stdout.write('#a: %f\n#b: %f\n#c: %f\n' % (a,b,c)) 
        sys.stdout.write('#x: %f\n#y: %f\n#z: %f\n' % (x,y,z)) 
