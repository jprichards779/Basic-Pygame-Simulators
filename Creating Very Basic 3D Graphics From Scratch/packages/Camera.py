from packages.rotation_functions import*
from packages.array_functions import *

"""
Camera class takes Transform(Shape()) as argument.
Takes care of both the projection of rays onto a flat screen and further transformation 
due to the camera's position and orientation. 

Camera_dict is used to store key value data momentarily to optimise performance.
Storing all data simply grinds the program to a hault, once the time it takes 
to search the dictionary becomes longer than the time to recalculate the values instead.
"""

class Camera:
    front_view=False
    plan_view=False
    def __init__(self, obj=None):
        self.size, self.threshold, self.zoom = 8, 0.1, 0.2
        assert self.zoom > 0
        self.obj = obj
        self.Camera_dict = {}
        self.camera_pos, self.object_pos = self.obj.camera_pos, self.obj.object_pos
        self.camera_pos[1] = norm_angle_list(self.camera_pos[1])
        self.X, self.Y, self.Z  = self.obj.out[0], self.obj.out[1], self.obj.out[2]
        self.colour = self.obj.colour
        # dictionary data 
        R = List_Reduction()
        ob = tuple(R.flatten_ND_list(self.object_pos))
        cam = tuple(R.flatten_ND_list(self.camera_pos))
        self.data = (ob,cam, self.colour, tuple(self.X), tuple(self.Y), tuple(self.Z))
        # calculated data 
        self.object_points = []
        self.orient_camera_XY_YZ()
        self.radial_dist = self.avg_radius() # Used to orer shape data by radial distance
        self.points = []
        self.output() 

    def output(self):
        x1,y1,z1 = self.object_points
        self.points = self.projection(x1,y1,z1, front_view = Camera.front_view, plan = Camera.plan_view)
    def avg_radius(self):
        x, y, z = self.object_points
        r = ((sum(x)/len(x))**2 + (sum(y)/len(y))**2 + (sum(z)/len(z))**2)**0.5
        return float(r)
    def orient_camera_XY_YZ(self):
        x,y,z = rotate(self.X, self.Y,self.Z, 
                        [-self.camera_pos[1][0],-self.camera_pos[1][1],self.camera_pos[1][2]])
        self.object_points = x,y,z
    def projection(self, x,y,z, front_view = False, plan = False):
        # if self.data in self.Camera_dict:
            # return self.Camera_dict.get(self.data)
        if plan and front_view:   
            front_view = False 
        if self.data not in self.Camera_dict:
            X_array, Y_array, Z_array = x,y,z
            obj_distance = min(Y_array[n] for n in range(len(Y_array)))
            if plan and not front_view:         
                self.size = 1
                VectorArray = [1 for i in X_array]
                X_axis = [round(self.size*X_array[n]*self.zoom*0.1/VectorArray[n], 
                        10) for n in range(len(X_array))]
                Y_axis = [round(-self.size*Y_array[n]*self.zoom*0.1/VectorArray[n], 
                        10) for n in range(len(X_array))]
                points = pygame_array(X_axis, Y_axis, WIDTH, HEIGHT)
                self.Camera_dict[self.data]=tuple(points)
            elif front_view and not plan:
                self.size = 5
                VectorArray = [1 for i in X_array]
                X_axis = [round(0.5*self.size*X_array[n]*self.zoom*0.1/VectorArray[n], 
                        10) for n in range(len(X_array))]
                Y_axis = [round(0.5*self.size*Z_array[n]*self.zoom*0.1/VectorArray[n], 
                        10) for n in range(len(X_array))]
                points = pygame_array(X_axis, Y_axis, WIDTH, HEIGHT)
                self.Camera_dict[self.data]=tuple(points)
            elif not plan and not front_view:
                VectorArray = Y_array
                self.size = self.size
                if obj_distance >= self.threshold > 0: 
                    X_axis = [round(self.size*X_array[n]*self.zoom/(VectorArray[n]*1), 
                            10) for n in range(len(X_array)) if VectorArray[n] !=0]
                    Y_axis = [round(self.size*Z_array[n]*self.zoom/(VectorArray[n]*1), 
                            10) for n in range(len(Z_array)) if VectorArray[n] !=0]
                    points = pygame_array(X_axis, Y_axis, WIDTH, HEIGHT)
                else: points = [(0,0), (0,0), (0,0)]   
                self.Camera_dict[self.data]=tuple(points)
        return self.Camera_dict[self.data]
