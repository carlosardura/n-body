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


class Node:
    def __init__(self, octree):
        self.octree = octree
        self.body = None
        self.R_cm = None
        self.M_cm = 0.0
        self.external = True
        self.children = [None] * 8

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
        index = (px >= self.octree.mid_x) | ((py >= self.octree.mid_y) << 1) | ((pz >= self.octree.mid_z) << 2)
        
        if self.children[index] is None:
        # initializes the octants when the first particle is added
            nx = self.octree.mid_x if (index & 1) else self.octree.xmin
            ny = self.octree.mid_y if (index & 2) else self.octree.ymin
            nz = self.octree.mid_z if (index & 4) else self.octree.zmin
            
            new_octree = Octree(nx, ny, nz, self.octree.size / 2.0, self.octree.depth + 1)
            self.children[index] = Node(new_octree)

        self.children[index].insertBody(particle)


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
            for child in self.children:
                if child is not None:
                    # recursion
                    dax, day, daz = child.multipole_acceptance_criterion(body, theta, epsilon, G)
                    ax += dax
                    ay += day
                    az += daz
            return ax, ay, az