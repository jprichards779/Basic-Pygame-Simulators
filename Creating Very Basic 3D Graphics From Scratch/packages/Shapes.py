from packages.Transforms import*
from packages.Camera import *
from packages.Screen import *

"""
Primary Shape class hosts input vertices list (containing a 2D array) 
while taking care of all transformation classes needed to display a shape.

Data is fed into all instances of this class, which instantiates the Transform 
class and stores data to be fed into the screen class via World class. 

Attempted to reformat using inheritance here but it seemed to cause some ugly 
side effects, so just stayed with this format for now.

The cubic() function form shape_functions module is designed to be malleable 
so various shapes can be built out of it. 
"""


class Shape:
    def __init__(self, vertices=[], object_pos = [[0,0,0] for n in range(5)], 
                 camera_pos = [[0,0,0], [0,0,0]], colour=Colours.DARK_BLUE):
        self.colour = colour 
        self.object_pos,self.camera_pos = object_pos, camera_pos
        self.vertices = vertices
        self.obj_points,self.objects, self.points,self.outline = [],[],[],[]
    def generate_points(self, colours=[], parabolic = False):
        ob, cam = self.object_pos, self.camera_pos
        if len(colours) < len(self.vertices):
            dif = len(self.vertices)-len(colours)
            for n in range(dif):
                colours.append(self.colour)  
        for n in range(len(self.vertices)):
            T = Transform(self.vertices[n], ob, cam, colours[n], parabolic)
            [x,y,z] = T.out
            self.objects.append(T) 


"""
Secondary Shape classes, make use of the framework provided by Shape() class 
to calculate data for different shapes"""
class Cube:
    def __init__(self, object_pos = [[0,0,0], [0,0,0], [0,0,0], [0,0,0], [0,0,0]], 
                 camera_pos = [[0,0,0], [0,0,0]], colour=Colours.DARK_BLUE, scale=1):
        assert scale > 0
        vertices = [cubic(vt=[[scale for n in range(3)] for n in range(4)], 
                                 vb=[[scale for n in range(3)] for n in range(4)])]
        self.v = [i[n] for i in vertices for n in range(len(i))]
        S = Shape(self.v,object_pos,camera_pos, colour)
        S.generate_points()
        self.obj_points, self.objects, self.points = S.obj_points, S.objects, S.points


class Six_Colour_Cube:
    def __init__(self, object_pos = [[0,0,0], [0,0,0], [0,0,0], [0,0,0], [0,0,0]], 
                 camera_pos = [[0,0,0], [0,0,0]], colour=Colours.DARK_BLUE, scale=1):
        self.colour=colour
        vertices = cubic(vt=[[scale for n in range(3)] for n in range(4)], 
                                vb=[[scale for n in range(3)] for n in range(4)])
        self.v = vertices
        S = Shape(self.v,object_pos,camera_pos, colour)
        S.generate_points(colours=[(100,0,0), (0,100,0), (0,0,200), 
                             (255,0,0), (190,190,200), (0,0,50)])
        self.obj_points, self.objects, self.points = S.obj_points, S.objects, S.points


"""Joining multiple cubes and removing sides that occupy the same points in space"""
class Block_Stacker:
    def __init__(self, object_pos = [[0,0,0], [0,0,0], [0,0,0], [0,0,0], [0,0,0]], 
                 camera_pos = [[0,0,0], [0,0,0]], colour=Colours.DARK_BLUE):
        x,y,z = object_pos[1][0], object_pos[1][1], object_pos[1][2]
        self.v, custom_colours = cube_stacker(locations=[[x-1,y,z],[x,y,z],[x+1,y,z]], 
                                                     colours =[(100,0,0), (0,255,0),(0,0,100)])
        S = Shape(self.v,object_pos,camera_pos, colour)
        S.generate_points(custom_colours)
        self.obj_points, self.objects, self.points = S.obj_points, S.objects, S.points


class Cylinder_Elem:
    angle_rad = math.pi/9
    angle_deg = angle_rad*(180/math.pi)
    def __init__(self, object_pos = [[0,0,0], [0,0,0], [0,0,0], [0,0,0], [0,0,0]], 
                 camera_pos = [[0,0,0], [0,0,0]], colour=Colours.DARK_BLUE, 
                 back=2,front=1, l=1, hole_size=0, sides = "top, front"):
        """If hole_size = front = back then thickness of ring = 0"""
        self.colour=colour
        hyp = 1
        X, Y, Z = hyp*math.sin(self.angle_rad/2), l, hyp*math.cos(self.angle_rad/2) # half angle each way = whole angle
        hole = hole_size
        self.vertices = [cubic2(vt=[[hole*X, Y, -hole*Z], [hole*X, Y, -hole*Z], 
                                           [hole*X, Y, -hole*Z], [hole*X, Y, -hole*Z]], 
                                  vb=[[back*X, Y, back*Z], [front*X, Y, front*Z], 
                                      [front*X, Y, front*Z], [back*X, Y, back*Z]],
                                  sides=sides)] 
        self.v = [i[n] for i in self.vertices for n in range(len(i))]
        S = Shape(self.v,object_pos,camera_pos, colour)
        S.generate_points(colours=[])
        self.obj_points, self.objects, self.points = S.obj_points, S.objects, S.points



class Horizontal_Tile:
    def __init__(self, 
                 object_pos = [[0,0,0], [0,0,0], [0,0,0], [0,0,0], [0,0,0]], 
                 camera_pos = [[0,0,0], [0,0,0]], colour=(255,255,255), 
                 location=[0,0,0], parabolic=False):
        self.parabolic =parabolic
        vt=[[1,1,1], [1,1,1], [1,1,1], [1,1,1]]
        vertices = [cubic(L=location, vt=vt, 
                          vb=[[1,1,1], [1,1,1], [1,1,1], [1,1,1]])]
        self.vertices = [i[1] for i in vertices]
        S = Shape(self.vertices,object_pos,camera_pos, colour)
        S.generate_points(parabolic=self.parabolic)
        self.obj_points, self.objects, self.points = S.obj_points, S.objects, S.points

class Vertical_Tile:
    def __init__(self, 
                 object_pos = [[0,0,0], [0,0,0], [0,0,0], [0,0,0], [0,0,0]], 
                 camera_pos = [[0,0,0], [0,0,0]], colour=(255,255,255), 
                 location=[0,0,0], parabolic=False):
        self.parabolic =parabolic
        vt=[[1,1,1], [1,1,1], [1,1,1], [1,1,1]]
        vertices = [cubic(L=location, vt=vt, 
                          vb=[[1,1,1], [1,1,1], [1,1,1], [1,1,1]])]
        self.vertices = [i[2] for i in vertices]
        S = Shape(self.vertices,object_pos,camera_pos, colour)
        S.generate_points(parabolic=self.parabolic)
        self.obj_points, self.objects, self.points = S.obj_points, S.objects, S.points


"""
Tertiary shape classes can be called in the World clases, or on their own in the main loop for inspection.
These build more complex shapes out ouf the secondary shape classes.

These could be simplified into a general case which accepts an input 
but this project is still being worked on. 
"""

class Horizontal_Sheet:
    def __init__(self, object_pos = [[0,0,0],[0,0,0], [0,0,0], [0,0,0], [0,0,0]], 
                 camera_pos = [[0,0,0], [0,0,0]], colour=Colours.DARK_BLUE, parabolic=False):
        self.colour=colour
        self.parabolic=parabolic
        self.ob, self.cam = object_pos, camera_pos
        self.input= [] 
        self.v, self.objects, self.points = [],[], []
        self.collect_elements()    
    def collect_elements(self):
        ob, cam = self.ob, self.cam
        x,y,z = ob[1][0], ob[1][1], ob[1][2]
        M = [Horizontal_Tile([ob[0], [x+i,y+j,z], ob[2], ob[3], ob[4]], 
                  cam, colour=self.colour, parabolic=self.parabolic) for i in range(-5,6) 
                  for j in range(-5,6)]
        for i in M:
            self.v.append(i.vertices)
            self.input.append(i)
            self.objects.append(i.objects)
    """Can get rid of all display functions not invoked in main loop"""
    def display(self, wireframe=False, colours=False, flatten=False):
        Screen(self.objects).camera_feed(wireframe, colours, flatten)


class Vertical_Sheet:
    def __init__(self, object_pos = [[0,0,0],[0,0,0], [0,0,0], [0,0,0], [0,0,0]], 
                 camera_pos = [[0,0,0], [0,0,0]], colour=Colours.DARK_BLUE, parabolic=False):
        self.colour=colour
        self.parabolic=parabolic
        self.ob, self.cam = object_pos, camera_pos
        self.input= [] 
        self.v, self.objects, self.points = [],[], []
        self.collect_elements()    
    def collect_elements(self):
        ob, cam = self.ob, self.cam
        x,y,z = ob[1][0]-5, ob[1][1], ob[1][2]
        M = [Vertical_Tile([ob[0], [x,y+j,z+i], ob[2], ob[3], ob[4]], 
                  cam, colour=self.colour, parabolic=self.parabolic) for i in range(-5,6) 
                  for j in range(-5,6)]
        for i in M:
            self.v.append(i.vertices)
            self.input.append(i)
            self.objects.append(i.objects)
    """Can get rid of all display functions not invoked in main loop"""
    def display(self, wireframe=False, colours=False, flatten=False):
        Screen(self.objects).camera_feed(wireframe, colours, flatten)



class Cylinder:
    angle = Cylinder_Elem.angle_deg
    def __init__(self, object_pos = [[0,0,0],[0,0,0], [0,0,0], [0,0,0], [0,0,0]], 
                 camera_pos = [[0,0,0], [0,0,0]], colour=Colours.DARK_BLUE,
                 back=1,front=1,l=1, hole_size=0, sides="top, base, left, right, back, front"):
        self.dimensions = [back,front,l, hole_size, sides]
        self.colour=colour
        self.ob, self.cam = object_pos, camera_pos
        self.input= [] 
        self.v, self.objects, self.points= [],[], []
        self.collect_elements()    
    def collect_elements(self):
        [b,f,l,h, s]=self.dimensions
        ob, cam = self.ob, self.cam
        angle=self.angle
        times = 360//angle
        a1,a2,a3 = ob[4][0], ob[4][1]*times, ob[4][2]
        a2 += self.angle/2
        M = [Cylinder_Elem([ob[0], ob[1], ob[2], ob[3], [a1,a2+angle*i,a3]], cam, self.colour, back=b, front=f, l=l, hole_size=h, sides=s) 
                           for i in range(int(times))]
                        # This for loop means custom colours apply to all listed elements
        for i in M:
            self.v.append(i.vertices)
            self.input.append(i)
            self.objects.append(i.objects)
    """Can get rid of all display functions not invoked in main loop"""
    def display(self, wireframe=False, colours=False, flatten=False):
        Screen(self.objects).camera_feed(wireframe, colours, flatten)


class Cylinder_Elem2:
    def __init__(self, object_pos = [[0,0,0], [0,0,0], [0,0,0], [0,0,0], [0,0,0]], 
                 camera_pos = [[0,0,0], [0,0,0]], colour=Colours.DARK_BLUE, 
                 back=2,front=1, l=1, hole_size=0, sides="top, front,back,base,left,right"):
        """If hole_size = front = back then thickness of ring = 0"""
        hyp = 1
        Y = 1*l
        angle = math.pi/18 # translates to 20 degree segments
        self.colour=colour
        self.angle=angle
        hole = hole_size
        X, Z = hyp*math.sin(angle/2), hyp*math.cos(angle/2) # half angle each way = whole angle
        self.vertices = [cubic2(vt=[[hole*X, Y, -hole*Z], [hole*X, Y, -hole*Z], 
                                           [hole*X, Y, -hole*Z], [hole*X, Y, -hole*Z]], 
                                  vb=[[back*X, Y, back*Z], [front*X, Y, front*Z], 
                                      [front*X, Y, front*Z], [back*X, Y, back*Z]],
                                  sides=sides)]
        # print(self.vertices)
        self.v = [i[n] for i in self.vertices for n in range(len(i))]
        S = Shape(self.v,object_pos,camera_pos, colour)
        S.generate_points()
        self.obj_points, self.objects, self.points = S.obj_points, S.objects, S.points


class Cylinder2:
    def __init__(self, object_pos = [[0,0,0],[0,0,0], [0,0,0], [0,0,0], [0,0,0]], 
                 camera_pos = [[0,0,0], [0,0,0]], colour=Colours.DARK_BLUE, 
                 back=1,front=1,l=1, hole_size=0, sides="top, base, left, right, back, front"):
        self.dimensions = [back,front,l, hole_size, sides]
        self.colour=colour
        self.ob, self.cam = object_pos, camera_pos
        self.input= [] 
        self.v, self.objects, self.points = [],[],[]
        self.collect_elements()    
    def collect_elements(self):
        [b,f,l,h, s]=self.dimensions
        ob, cam = self.ob, self.cam
        angle = 10 # ! MUST be equal to frustum angle
        times = 360//angle
        a1,a2,a3 = ob[4][0], ob[4][1]*times, ob[4][2]
        M = [Cylinder_Elem2([ob[0], ob[1], ob[2], ob[3], [a1,a2+angle*i,a3]], cam, self.colour, 
                           back=b, front=f, l=l, hole_size=h, sides=s) 
                           for i in range(int(times))]
        for i in M:
            self.v.append(i.vertices)
            self.input.append(i)
            self.objects.append(i.objects)
            for n in i.points:
                self.points.append(n)
    """Can get rid of all display functions not invoked in main loop"""
    def display(self, wireframe=False, colours=False, flatten=False):
        Screen(self.objects).camera_feed(wireframe, colours, flatten)



