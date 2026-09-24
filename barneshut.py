import numpy as np

# octant in 3D space
class Octree:
    def __init__(self, rx: float, ry: float, rz: float, size: float, depth: int = 0):
        self.xmin = rx
        self.ymin = ry
        self.zmin = rz
        self.size = size
        self.depth = depth

        self.mid_x = rx + size / 2.0
        self.mid_y = ry + size / 2.0
        self.mid_z = rz + size / 2.0

    def NEZ(self): # north, east, zenith
        return Octree(self.mid_x, self.mid_y, self.mid_z, self.size/2.0, self.depth + 1)

    def NWZ(self): # north, west, zenith
        return Octree(self.xmin, self.mid_y, self.mid_z, self.size/2.0, self.depth + 1)

    def SWZ(self): # south, west, zenith
        return Octree(self.xmin, self.ymin, self.mid_z, self.size/2.0, self.depth + 1)

    def SEZ(self): # south, east, zenith
        return Octree(self.mid_x, self.ymin, self.mid_z, self.size/2.0, self.depth + 1)

    def NEN(self): # north, east, nadir
        return Octree(self.mid_x, self.mid_y, self.zmin, self.size/2.0, self.depth + 1)

    def NWN(self): # north, west, nadir
        return Octree(self.xmin, self.mid_y, self.zmin, self.size/2.0, self.depth + 1)

    def SWN(self): # south, west, nadir
        return Octree(self.xmin, self.ymin, self.zmin, self.size/2.0, self.depth + 1)

    def SEN(self): # south, east, nadir
        return Octree(self.mid_x, self.ymin, self.zmin, self.size/2.0, self.depth + 1)


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
        if self.M_cm > 0:   # non-empty nodes
            if self.octree.depth < 20:   # depth limit 
                if self.external:
                    self.external = False   # more than one particle                
                    self._new_octant(self.body)
                    self.body = None

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

        is_east = px >= self.octree.mid_x
        is_north = py >= self.octree.mid_y
        is_zenith = pz >= self.octree.mid_z

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


    def multipole_acceptance_criterion(self, body, theta: float, epsilon: float, G: float):
        if self.M_cm == 0 or self.body is body:   # ignore empty nodes and self-interactions
            return 0.0, 0.0, 0.0

        dx = self.R_cm[0] - body.r[0]
        dy = self.R_cm[1] - body.r[1]
        dz = self.R_cm[2] - body.r[2]
        d2 = dx**2 + dy**2 + dz**2 + epsilon**2

        if (self.octree.size**2 / d2 < theta**2) or self.external:   # valid approximation
            a = G * self.M_cm / (d2 * np.sqrt(d2))
            return a*dx, a*dy,a*dz
        
        else:
            ax, ay, az = 0.0, 0.0, 0.0
            children = [self.NWZ, self.NEZ, self.SWZ, self.SEZ, self.NWN, self.NEN, self.SWN, self.SEN]
            
            for child in children:
                if child is not None:
                    # recursion
                    dax, day, daz = child.multipole_acceptance_criterion(body, theta, epsilon, G)
                    ax += dax
                    ay += day
                    az += daz
            return ax, ay, az