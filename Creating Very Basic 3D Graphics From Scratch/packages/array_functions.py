import math  
import pygame
from pygame.locals import* 
from packages.rotation_functions import *
from packages.shape_functions import*
pygame.init()
WIDTH = 1000
HEIGHT = 600
screen = pygame.display.set_mode((WIDTH, HEIGHT))
PI = math.pi


"""
The functions here contain matrices which dictate the traversal of 
vertices such that they appear as desired on a graph plotting interface.
Crucial rendering functions, invoked in the Screen class, can be found here.

Also contains algorithms for reducing arrays structures down 
to 1D lists and removing multiple vertices occupying the 
same coordinates. 
Some may be redundant. 
"""


def non_dupes_only(i=[], n=0):
    if n not in i:
        i.append(n)

"""Helper functions"""
class List_Reduction:
    def __init__(self):
        self.output_list = []
    def non_dupes_only(self, i=[], n=0):
        if n not in i:
            i.append(n)
    def flatten_ND_list(self, array):
        for item in array:
            if type(item) != list or item == " ":
                self.output_list.append(item)
            else:
                self.flatten_ND_list(item)
        return self.output_list
    def flatten_2D(self, L):
        m = []
        for n in L:
            if type(n)!=list:
                self.non_dupes_only(m, n)
            elif type(n) == list:
                for i in n:
                    self.non_dupes_only(m, i)
        return m
# j = [[[1,[[]],3],[4,5]], [[6,7,8],[9]]]
# print(List_Reduction().flatten_ND_list(j))

class Weight_Methods:
    def __init__(self):
        self.sides = {}
    def shed_weight2(self, M, P):
        if type(M) == tuple:
            M = list(M) 
        out = M
        for n in P:
            if n in out:
                out.remove(n)
            else: out.append(n)
        return out
    def uniqueSides_only2(self,cubes=[]):
        if len(cubes) > 1:
            prev = cubes.pop(0)
            curr = cubes[0]
            r = self.shed_weight2(prev,curr)
            cubes[0] = r
            self.uniqueSides_only2(cubes)
        return cubes[0]

def pygame_array(list1, list2, WIDTH,HEIGHT):
    scaled_coordinates = [] 
    for n in range(len(list1)):
        x = 0.5*WIDTH*list1[n] + 0.5*WIDTH
        y = 0.5*WIDTH*list2[n] + 0.5*HEIGHT
        scaled_coordinates.append((x, y))
    return scaled_coordinates


def pg_array_to_points(arr, scale=1):
    X, Y = [], []
    for n in arr:
        x = -1 + 2*n[0]/WIDTH
        y = (2*n[1]-HEIGHT)/WIDTH
        x*=scale
        y*=scale
        X.append(x)
        Y.append(y)
    return X, Y

def draw_axes(lines, colour=(100,0,50)):
    global_y = [(WIDTH/2, HEIGHT), (WIDTH/2, -HEIGHT/2)]
    global_x = [(-WIDTH, HEIGHT/2), (WIDTH, HEIGHT/2)]
    pygame.draw.lines(lines, colour, False, global_x)
    pygame.draw.lines(lines, colour, False, global_y)

def draw_object(lines, line_colour, shape_colour, points, transparent = False):
    if not transparent:
        # pygame.draw.polygon(lines, colour, points, width=10)
        pygame.draw.polygon(lines, shape_colour, points)
        pygame.draw.lines(lines, line_colour, False, points, width=1)
    else: pygame.draw.lines(lines, line_colour, False, points, width=1)


def highlight_vertices(lines, line_colour, points, radius=1):
    for n in points:
        pygame.draw.circle(lines,color=line_colour,center=n,radius=radius)

def render_object2(lines, colour, points, transparent = False):
    pygame.draw.polygon(lines, colour, points)


def shed_weight(M, P):
    if type(M) == tuple:
        M = list(M) # *
    out = M
    for n in P:
        if n in out:
                out.remove(n)
        else: out.append(n)
    return out

def uniqueSides_only(cubes=[]):
    if len(cubes) > 1:
        prev = cubes.pop(0)
        curr = cubes[0]
        r = shed_weight(prev,curr)
        cubes[0] = r
        uniqueSides_only(cubes)
    return cubes[0]


"""
WARNING: 
These stacker functions, only work if our quadrahedrons are cubes. I didn't predict this. 
"""
def cube_stacker(locations=[[0,0,0]], colours =[(0,0,0)]):
    v=[]
    for i in locations:
        if i not in v:
            v.append(i)
    locations = v
    Methods = Weight_Methods()
    assert len(colours)>=len(locations)
    dict1 = {}
    colours, colours_2 = tuple(colours), []
    v = [cubic2(L=i) for i in locations]
    v2 = Methods.uniqueSides_only2(v[:]) # <-- Must be copied so length not affected
    for i in range(len(v)):
        for j in v[i]:
                dict1[j] = colours[i]
    for i in range(len(v2)):
        if v2[i] in dict1:
                colours_2.append(dict1[v2[i]])
    return v2, colours_2
# v,c = cube_stacker(locations=[[0,0,0], [1,0,0]], colours=[(100,0,0), (0,100,0)])
# print(len(v))



