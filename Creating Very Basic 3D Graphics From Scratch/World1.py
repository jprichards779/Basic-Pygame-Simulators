
"""
The following World class employs the scaling of inputs from the event loop
to manipulate the motion of different shapes
"""
from packages.Shapes import *

class World:
    ob,cam =[],[]
    def __init__(self, app, colour=Colours.DARK_BLUE):
        self.app=app
        self.front_view = app.front_view
        self.plan_view = app.plan_view
        self.background_colour = app.background_colour
        self.colour=colour
        self.line_colour=app.line_colour
        self.ob, self.cam = self.app.object_pos, self.app.camera_pos
        self.input= [] 
        self.objects, self.v, self.points = [],[],[]
        self.collect_elements() 
    def collect_elements(self):
        ob, cam = self.ob, self.cam 
        x,y,z = ob[1][0], ob[1][1], ob[1][2]
        M = [Cylinder([ob[0], [x, y-1.5, z], ob[2], ob[3], ob[4]],
                       cam,colour=(10,60,150),front=1.1,back=1.8,l=0.5, hole_size=0, sides = "top,front,back"),
             Cylinder([ob[0], [x, y+1.5, z], ob[2], ob[3], ob[4]],
                       cam,colour=(10,60,150),front=1.8,back=1.1,l=0.5, hole_size=0, sides = "top,back,front"),
             Cylinder([ob[0], [x, y, z], ob[2], ob[3], ob[4]],
                       cam,colour=(10,60,150),front=1,back=1,l=2.5, hole_size=0, sides = "top"),
             Block_Stacker([[ob[0][0]-4, ob[0][1]+4, ob[0][2]], [x, y, z], ob[2], ob[3], [ob[4][0], ob[4][1], ob[4][2]]],cam,)]
        for i in M:
            self.input.append(i)
            self.objects.append(i.objects) # Each element delivered by Transfrom class
            self.v.append(i.v)
            for n in i.points:
                self.points.append(n)
    def display(self, wireframe=False, colours=False, flatten=False):
        Screen(self).camera_feed(wireframe, colours, flatten)
        