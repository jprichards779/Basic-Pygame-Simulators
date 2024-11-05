"""
The following World class employs the scaling of inputs from the event loop
to manipulate the motion of different shapes
"""
from packages.Shapes import *

class World4:
    ob,cam =[],[]
    colour = (0,255,0)
    def __init__(self, app, colour=(0,100,200)):
        self.app=app
        self.front_view = app.front_view
        self.plan_view = app.plan_view
        self.background_colour = app.background_colour
        self.colour=World4.colour
        self.line_colour = app.line_colour
        self.ob, self.cam = self.app.object_pos, self.app.camera_pos
        self.objects, self.v, self.points = [],[],[]
        self.collect_elements() 
    def collect_elements(self):
        ob, cam = self.ob, self.cam  
        x0,y0,z0 = ob[0][0], ob[0][1], ob[0][2]
        x1,y1,z1 = ob[1][0], ob[1][1], ob[1][2]
        x4,y4,z4 = ob[4][0], ob[4][1], ob[4][2]
        # 1. Double angle transform to produce a 
        #    vertical sheet from the horizontal one
        # 2. Easier to create a sheet composed of Tile objects, 
        #    hard coded hard from a specially made vertical Tile 
        #    for future convenience.

        # ** 2. 
        # Vertical_Sheet([[x0,y0,z0], [0*x1, y1-5, z1-6], ob[2], ob[3], [0,0,0]], cam, colour = self.app.sheet_colour, parabolic=False),
        # Vertical_Sheet([[x0,y0,z0], [0*x1, y1+5, z1-6], ob[2], ob[3], [0,0,0]],  cam, colour = self.app.sheet_colour, parabolic=False),


        M = [Horizontal_Sheet([[x0,y0,z0], [0*x1, y1-5, z1], ob[2], ob[3], [0,0,0]], cam, colour = self.app.sheet_colour, parabolic=False),
             Horizontal_Sheet([[x0,y0,z0], [0*x1, y1+5, z1], ob[2], ob[3], [0,0,0]],  cam, colour = self.app.sheet_colour, parabolic=False),
             # * 1.
             Horizontal_Sheet([[x0,y0,z0], [0*x1-5, y1-6, z1+6], [ob[2][0], ob[2][1]+90, ob[2][2]+90], [ob[3][0], ob[3][1], ob[3][2]], [0,0,0]], cam, colour = self.app.sheet_colour, parabolic=False),
             Horizontal_Sheet([[x0,y0,z0], [0*x1+5, y1-6, z1+6], [ob[2][0], ob[2][1]+90, ob[2][2]+90], [ob[3][0], ob[3][1], ob[3][2]], [0,0,0]], cam, colour = self.app.sheet_colour, parabolic=False),
             # * 2. 
             # Commenting out the two previous elements and uncomment the following
             # before adding ** 2 will create the same result.
             Cylinder([[x0,y0,z0+1.5], [x1-50, y1-1.5, z1], [ob[2][0]+90, ob[2][1], ob[2][2]], ob[3], ob[4]],
                       cam, self.colour, front=1,back=1.8,l=0.5, hole_size=0.5, sides = "top,front"),
             Cylinder([[x0,y0,z0+1.5], [x1-50, y1+1.5, z1], [ob[2][0]+90, ob[2][1], ob[2][2]], ob[3], ob[4]],
                       cam, self.colour, front=1.8,back=1,l=0.5, hole_size=0.5, sides = "top,back"),
             Cylinder([[x0,y0,z0+1.5], [x1-50, y1, z1], [ob[2][0]+90, ob[2][1], ob[2][2]], ob[3], ob[4]],
                       cam, self.colour, front=1,back=1,l=2.5, hole_size=0, sides = "top"),
             
             Cylinder([[x0,y0,z0+1.5], [x1-53, y1-1.5, z1], [ob[2][0]+90, ob[2][1], ob[2][2]], ob[3], ob[4]],
                       cam, self.app.cylinder_colour, front=1,back=1.8,l=0.5, hole_size=0.5, sides = "top,front"),
             Cylinder([[x0,y0,z0+1.5], [x1-53, y1+1.5, z1], [ob[2][0]+90, ob[2][1], ob[2][2]], ob[3], ob[4]],
                       cam, self.app.cylinder_colour, front=1.8,back=1,l=0.5, hole_size=0.5, sides = "top,back"),
             Cylinder([[x0,y0,z0+1.5], [x1-53, y1, z1], [ob[2][0]+90, ob[2][1], ob[2][2]], ob[3], ob[4]],
                       cam, self.app.cylinder_colour, front=1,back=1,l=2.5, hole_size=0.5, sides = "top"),
                       

             Cylinder([[x0,y0,z0+1.5], [x1-70, y1-1.5, z1], [ob[2][0]+90, ob[2][1], ob[2][2]], ob[3], ob[4]],
                       cam, self.app.cylinder_colour, front=1,back=1.8,l=0.5, hole_size=0.5, sides = "top,front"),
             Cylinder([[x0,y0,z0+1.5], [x1-70, y1+1.5, z1], [ob[2][0]+90, ob[2][1], ob[2][2]], ob[3], ob[4]],
                       cam, self.app.cylinder_colour, front=1.8,back=1,l=0.5, hole_size=0.5, sides = "top,back"),
             Cylinder([[x0,y0,z0+1.5], [x1-70, y1, z1], [ob[2][0]+90, ob[2][1], ob[2][2]], ob[3], ob[4]],
                       cam, self.app.cylinder_colour, front=1,back=1,l=2.5, hole_size=0.5, sides = "top"),
             
             Cylinder([[x0,y0,z0+1.5], [x1-73, y1-1.5, z1], [ob[2][0]+90, ob[2][1], ob[2][2]], ob[3], ob[4]],
                       cam, self.colour, front=1,back=1.8,l=0.5, hole_size=0.5, sides = "top,front"),
             Cylinder([[x0,y0,z0+1.5], [x1-73, y1+1.5, z1], [ob[2][0]+90, ob[2][1], ob[2][2]], ob[3], ob[4]],
                       cam, self.colour, front=1.8,back=1,l=0.5, hole_size=0.5, sides = "top,back"),
             Cylinder([[x0,y0,z0+1.5], [x1-73, y1, z1], [ob[2][0]+90, ob[2][1], ob[2][2]], ob[3], ob[4]],
                       cam, self.colour, front=1,back=1,l=2.5, hole_size=0.5, sides = "top")]
        for i in M:
            self.objects.append(i.objects) # Each element delivered by Transfrom class
            self.v.append(i.v) # <--- vertices for each element /slice
            for n in i.points:
                self.points.append(n)
    def display(self, wireframe=False, colours=False, flatten=False):
        Screen(self).camera_feed(wireframe, colours, flatten)