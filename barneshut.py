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


    def inOctree(self, body):
            x, y, z = body.r
            return (self.xmin <= x < self.xmax) and \
                (self.ymin <= y < self.ymax) and \
                (self.zmin <= z < self.zmax)


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

                self.NWZ = Node(self.octree.NWZ())
                self.NEZ = Node(self.octree.NEZ())
                self.SWZ = Node(self.octree.SWZ())
                self.SEZ = Node(self.octree.SEZ())
                self.NWN = Node(self.octree.NWN())
                self.NEN = Node(self.octree.NEN())
                self.SWN = Node(self.octree.SWN())
                self.SEN = Node(self.octree.SEN())
                
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
        if self.NWZ.octree.inOctree(particle):
            self.NWZ.insertBody(particle)
        elif self.NEZ.octree.inOctree(particle):
            self.NEZ.insertBody(particle)
        elif self.SWZ.octree.inOctree(particle):
            self.SWZ.insertBody(particle)
        elif self.SEZ.octree.inOctree(particle):
            self.SEZ.insertBody(particle)
        elif self.NWN.octree.inOctree(particle):
            self.NWN.insertBody(particle)
        elif self.NEN.octree.inOctree(particle):
            self.NEN.insertBody(particle)
        elif self.SWN.octree.inOctree(particle):
            self.SWN.insertBody(particle)
        elif self.SEN.octree.inOctree(particle):
            self.SEN.insertBody(particle)