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