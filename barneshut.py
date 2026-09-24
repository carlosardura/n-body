import numpy as np

# octant in 3D space
class Octree:
    def __init__(self, rx: float, ry: float, rz: float, size: float):
        self.xmin = rx
        self.ymin = ry
        self.zmin = rz
        self.size = size

        self.xmax = rx + size
        self.ymax = ry + size
        self.zmax = rz + size

    def NEZ(self): # north, east, zenith
        return Octree(self.xmin+self.size/2.0, self.ymin+self.size/2.0, self.zmin+self.size/2.0, self.size/2.0)

    def NWZ(self): # north, west, zenith
        return Octree(self.xmin, self.ymin+self.size/2.0, self.zmin+self.size/2.0, self.size/2.0)

    def SWZ(self): # south, west, zenith
        return Octree(self.xmin, self.ymin, self.zmin+self.size/2.0, self.size/2.0)

    def SEZ(self): # south, east, zenith
        return Octree(self.xmin+self.size/2.0, self.ymin, self.zmin+self.size/2.0, self.size/2.0)

    def NEN(self): # north, east, nadir
        return Octree(self.xmin+self.size/2.0, self.ymin+self.size/2.0, self.zmin, self.size/2.0)

    def NWN(self): # north, west, nadir
        return Octree(self.xmin, self.ymin+self.size/2.0, self.zmin, self.size/2.0)

    def SWN(self): # south, west, nadir
        return Octree(self.xmin, self.ymin, self.zmin, self.size/2.0)

    def SEN(self): # south, east, nadir
        return Octree(self.xmin+self.size/2.0, self.ymin, self.zmin, self.size/2.0)


class Node:
    def __init__(self, octree):
        self.octree = octree
        self.body = None
        self.R_cm = None
        self.M_cm = 0.0
        self.external = True

        self.NWZ = None
        self.NEZ = None
        self.SWZ = None
        self.SEZ = None
        self.NWN = None
        self.NEN = None
        self.SWN = None
        self.SEN = None

    def insertBody(self, body):
        if self.body is not None:   # non-empty nodes
            if self.external:
                self.external = False   # more than one particle                
                self._new_octant(self.body)
            self._new_octant(body)

            # updates the center of mass
            self.R_cm = (self.M_cm * self.R_cm + body.m * body.r) / (self.M_cm + body.m)
            self.M_cm += body.m
            
        else:   # empty nodes
            self.body = body
            self.R_cm = body.r
            self.M_cm = body.m
            self.external = True


    def _new_octant(self, particle):
        px, py, pz = particle.r
        mid_x = self.octree.xmin + self.octree.size / 2.0
        mid_y = self.octree.ymin + self.octree.size / 2.0
        mid_z = self.octree.zmin + self.octree.size / 2.0

        is_east = px >= mid_x
        is_north = py >= mid_y
        is_zenith = pz >= mid_z

        if is_zenith:
            if is_north:
                if is_east: # NEZ
                    if self.NEZ is None: self.NEZ = Node(self.octree.NEZ())
                    self.NEZ.insertBody(particle)
                else:       # NWZ
                    if self.NWZ is None: self.NWZ = Node(self.octree.NWZ())
                    self.NWZ.insertBody(particle)
            else:
                if is_east: # SEZ
                    if self.SEZ is None: self.SEZ = Node(self.octree.SEZ())
                    self.SEZ.insertBody(particle)
                else:       # SWZ
                    if self.SWZ is None: self.SWZ = Node(self.octree.SWZ())
                    self.SWZ.insertBody(particle)
        else:
            if is_north:
                if is_east: # NEN
                    if self.NEN is None: self.NEN = Node(self.octree.NEN())
                    self.NEN.insertBody(particle)
                else:       # NWN
                    if self.NWN is None: self.NWN = Node(self.octree.NWN())
                    self.NWN.insertBody(particle)
            else:
                if is_east: # SEN
                    if self.SEN is None: self.SEN = Node(self.octree.SEN())
                    self.SEN.insertBody(particle)
                else:       # SWN
                    if self.SWN is None: self.SWN = Node(self.octree.SWN())
                    self.SWN.insertBody(particle)