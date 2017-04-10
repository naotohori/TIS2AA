#!/usr/bin/env python
#vim:fileencoding=UTF-8
'''
@author: Naoto Hori
'''

import math
from numpy import dot

class Coord(object) :
    def __init__(self, x=0.0, y=0.0, z=0.0) :
        self.x = x
        self.y = y
        self.z = z

    def get_as_tuple(self):
        return (self.x, self.y, self.z)
    
    def get_as_list(self):
        return [self.x, self.y, self.z]
    
    def get_as_ndarray(self):
        import numpy
        return numpy.array([self.x, self.y, self.z])
    
    def put_as_list(self, xyz):
        self.x = xyz[0]
        self.y = xyz[1]
        self.z = xyz[2]
    
    def distance(self, another_coord):
        return math.sqrt((another_coord.x - self.x) ** 2 
                       + (another_coord.y - self.y) ** 2
                       + (another_coord.z - self.z) ** 2)

    def dot(self, another_coord):
        return self.x * another_coord.x + self.y * another_coord.y + self.z * another_coord.z

    def cross(self, another_coord):
        ret = Coord()
        ret.x = self.y * another_coord.z - self.z * another_coord.y
        ret.y = self.z * another_coord.x - self.x * another_coord.z
        ret.z = self.x * another_coord.y - self.y * another_coord.x
        return ret

    def move(self, delta_coord):
        self.x += delta_coord.x
        self.y += delta_coord.y
        self.z += delta_coord.z

    def transform(self, mtx): 
        c = [self.x, self.y, self.z, 1.0] 
        c = dot(mtx,c)
        self.x = c[0]
        self.y = c[1]
        self.z = c[2]
    
    def __add__(self, other):
        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self
    
    def __sub__(self, other):
        rt = Coord()
        rt.x = self.x - other.x
        rt.y = self.y - other.y
        rt.z = self.z - other.z
        return rt
        
    def __div__(self, n):
        '''重心を求める際などに使う'''
        self.x /= n
        self.y /= n
        self.z /= n
        return self

    def __mul__(self, n):
        '''重心を求める際などに使う'''
        self.x *= n
        self.y *= n
        self.z *= n
        return self

    def norm(self,):
        import math
        return math.sqrt(self.dot(self))

def Angle(c1,c2,c3):
    v13 = c1 - c3
    v12 = c1 - c2
    cos_theta = v13.dot(v12) / sqrt(v13.dot(v13) * v12.dot(v12))
    return cos_theta
