import numpy as np

class Cell:
    __slots__ = ['xmin', 'ymin', 'zmin', 'mid_x', 'mid_y', 'mid_z', 'size', 'depth', 'body', 'R_cm', 'M_cm', 'external', 'children']

    def __init__(self, rx: float, ry: float, rz: float, size: float, depth: int = 0):
        self.xmin = rx
        self.ymin = ry
        self.zmin = rz
        
        self.mid_x = rx + size / 2.0
        self.mid_y = ry + size / 2.0
        self.mid_z = rz + size / 2.0

        self.size = size
        self.depth = depth

        self.body = None
        self.R_cm = None
        self.M_cm = 0.0
        self.external = True
        self.children = [None] * 8

    def insertBody(self, body):
        if self.M_cm > 0:   # non-empty nodes
            if self.depth < 20:   # depth limit 
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
        index = (px >= self.mid_x) | ((py >= self.mid_y) << 1) | ((pz >= self.mid_z) << 2)
        
        if self.children[index] is None:
        # initializes the octants when the first particle is added
            nx = self.mid_x if (index & 1) else self.xmin
            ny = self.mid_y if (index & 2) else self.ymin
            nz = self.mid_z if (index & 4) else self.zmin
            
            self.children[index] = Cell(nx, ny, nz, self.size / 2.0, self.depth + 1)
        self.children[index].insertBody(particle)


    def multipole_acceptance_criterion(self, body, theta: float, epsilon: float, G: float):
        if self.M_cm == 0 or self.body is body:   # ignore empty nodes and self-interactions
            return 0.0, 0.0, 0.0

        ax, ay, az = 0.0, 0.0, 0.0
        stack = [self]
        while stack:
            current = stack.pop()

            if current.M_cm == 0 or current.body is body:
                continue

            dx = current.R_cm[0] - body.r[0]
            dy = current.R_cm[1] - body.r[1]
            dz = current.R_cm[2] - body.r[2]
            d2 = dx**2 + dy**2 + dz**2 + epsilon**2

            if (current.size**2 / d2 < theta**2) or current.external:   # valid approximation
                a = G * current.M_cm / (d2 * np.sqrt(d2))
                ax += a * dx
                ay += a * dy
                az += a * dz
        
            else:
                for child in current.children:
                    if child is not None:
                        stack.append(child)

        return ax, ay, az