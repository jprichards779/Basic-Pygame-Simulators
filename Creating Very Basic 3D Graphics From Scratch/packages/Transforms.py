from packages.rotation_functions import*
from packages.array_functions import *

"""
This class is responsible for transforming the coordinates of our shapes.
I accepts inputs corresponding to coordinates, positioning, colour of each fragment/element,
and parameters for surface curvature such as parabolic.

        object_pos = [[0,0,0],     [0,0,0],        [0,0,0],[0,0,0],        [0,0,0]]
        ----------------------------------------------------------------------------------
                    World         Shape         Shape Orientation       Shape Orientation    
                    Positioning   location      a_xy, a_xz, a_yz        a_xy, a_xz, a_yz
                    wrt Camera    in World      about World's local     about Shape's own
                    x,y,z         x,y,z         origin                  origin

                    
        camera_pos = [[0,0,0],      [0,0,0]]    
        -----------------------------------------------------------------           
                    Camera        Camera orientation  
                    Positioning   a_xy, a_xz, a_yz  
                    x,y,z         Retotes the world about camera 
                                  Giving apperance of camera turning

"""

class Transform:
    def __init__(self, coordinates = [], 
                 object_pos = [[0,0,0],[0,0,0], [0,0,0], [0,0,0], [0,0,0]], 
                 camera_pos = [[0,0,0], [0,0,0]], colour=None, parabolic=False):
        R = List_Reduction()
        self.parabolic = parabolic
        self.Trans_dict={} 
        self.coordinates, self.colour = coordinates, colour
        self.object_pos, self.camera_pos = object_pos, camera_pos
        self.object_pos[2] = norm_angle_list(self.object_pos[2])
        self.object_pos[3] = norm_angle_list(self.object_pos[3])
        self.object_pos[4] = norm_angle_list(self.object_pos[4])
        coor = tuple(R.flatten_ND_list(self.coordinates))
        ob   = tuple(R.flatten_ND_list(self.object_pos))
        cam  = tuple(R.flatten_ND_list(self.camera_pos))
        self.data = (coor,ob,cam)
        self.local = []
        self.out = self.output() 
        # ^ This is then stored in the Shape class, ready to be entered into 
        #   the Screen class for rendering
    def output(self):
        if self.data in self.Trans_dict: 
            vals = self.Trans_dict.get(self.data)
            return list(vals) 
        elif self.data not in self.Trans_dict: 
            self.shape_motion()
            self.Trans_dict[self.data] = tuple(self.group_motion(self.local[0], self.local[1], self.local[2]))
            return self.group_motion(self.local[0], self.local[1], self.local[2])
    def transpose(self):
        M = [[n[0] for n in self.coordinates], 
             [n[1] for n in self.coordinates],
             [n[2] for n in self.coordinates]]
        return M
    def shape_motion(self):
        Arr = self.transpose()
        x, y, z  = Arr[0],Arr[1],Arr[2]
        x,y,z = rotate(x,y,z, angles=self.object_pos[4])
        x = [n + self.object_pos[1][0] for n in x]
        y = [n + self.object_pos[1][1] for n in y]
        z = [n + self.object_pos[1][2] for n in z]
        self.local = [x, y, z]
        if self.parabolic:
            self.local[2] = self.equation_z(self.local[0], self.local[1], self.local[2])
    def equation_z(self,x,y,z):
        z = [z[n] - 0.1*(x[n]**2+y[n]**2) for n in range(len(z))]
        return z
    def group_motion(self, x, y, z):
        x,y,z = rotate(x,y,z, angles=self.object_pos[2])
        x,y,z = rotate(x,y,z, angles=self.object_pos[3])
        X = [n + self.object_pos[0][0]-self.camera_pos[0][0] for n in x]
        Y = [n + self.object_pos[0][1]-self.camera_pos[0][1] for n in y]
        Z = [n - self.object_pos[0][2]+self.camera_pos[0][2] for n in z]
        return [X,Y,Z]
    


