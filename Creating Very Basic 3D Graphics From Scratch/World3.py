"""
The following World class employs the scaling of inputs from the event loop
to manipulate the motion of different shapes
"""
from packages.Shapes import *

class World3:
    ob,cam =[],[]
    def __init__(self, app, colour=(0,100,200)):
        self.app=app
        self.front_view = app.front_view
        self.plan_view = app.plan_view
        self.background_colour = app.background_colour
        self.colour=app.sheet_colour
        self.line_colour = app.line_colour
        self.ob, self.cam = self.app.object_pos, self.app.camera_pos
        self.input= [] 
        self.objects, self.v, self.points = [],[],[]
        self.collect_elements() 
    def collect_elements(self):
        ob, cam = self.ob, self.cam 
        x,y,z = ob[0][0], ob[0][1], ob[0][2]
        M = [Horizontal_Sheet([[x,y,z], ob[1], ob[2], ob[3], ob[4]], cam, parabolic=True)]
        for i in M:
            self.input.append(i)
            self.objects.append(i.objects) # Each element delivered by Transfrom class
            self.v.append(i.v) # <--- vertices for each element /slice
            for n in i.points:
                self.points.append(n)
    def display(self, wireframe=False, colours=False, flatten=False):
        Screen(self).camera_feed(wireframe, colours, flatten)