from packages.rotation_functions import*
from packages.array_functions import *
from packages.Camera import*

""" 
To colour our shapes in Pygame, we must create an instance of Pygame's surface class. 
This is reinstanciated in our World class for every itteration in the event loop. 
The effect is to create a flip chart. The Sceen class below is structured to deliver 
all the shape elements onto the same surface while ensuring they are drawn in the 
correct order - furthest objects are not drawn over nearest.
 
The camera_feed() function ensures objects nearest objects are drawn last.
    wireframe=False means shapes are opaque, 
    colours=False means all shapes are coloured black,
    flatten=False means all radial distances are calculated for each fragment.   
"""
class Screen:
    def __init__(self, world, no_lines=True): 
        self.lines=pygame.Surface((WIDTH, HEIGHT))  
        self.lines.fill(world.background_colour)
        self.line_colour = world.line_colour
        self.world=world
        self.no_lines = no_lines
        self.objects = world.objects
        self.avg_distances = []
        self.data_load = 0
    def camera_feed(self, wireframe = False, colours=True, flatten=False):
        # Dictionary to store optimal amount of calculated values 
        # minimising beed to recalculate previously visitted coordinates
        Screen_dict = {} 
        objects = List_Reduction().flatten_ND_list(self.objects)
        transparent = wireframe  
        if flatten: 
            # flatten means objects are not reordered by average radial distance of each
            # element.
            colours = False
            for n in objects: 
                CAM = Camera(n)
                Camera.front_view = self.world.front_view
                Camera.plan_view = self.world.plan_view
                points = CAM.points
                draw_object(self.lines, n.colour, points, transparent)
        elif not flatten: 
            for n in objects:
                CAM = Camera(n)
                Camera.front_view = self.world.front_view
                Camera.plan_view = self.world.plan_view
                colour = CAM.colour
                average = CAM.radial_dist
                if average in self.avg_distances: 
                    # In cases where shapes all have equal radial distance from lens
                    average += objects.index(n)/(10**10)
                    # ...this is simplistic mathematical way to ensure no cases are missed 
                self.avg_distances.append(average)
                self.avg_distances.sort(reverse=True) # reordering
                if average not in Screen_dict:
                    Screen_dict[average] = (CAM.points,CAM.colour)
            for i in self.avg_distances:
                (points,colour) = Screen_dict[i]
                if colours:
                    draw_object(self.lines, self.line_colour, colour, points, transparent)
                elif not colours: 
                    draw_object(self.lines, self.line_colour, Colours.BLACK, points, transparent) 
                    if self.no_lines:
                        size = 10/i**0.5
                        highlight_vertices(self.lines,(255,255,255), points,radius=size)
                draw_axes(self.lines)
        screen.blit(self.lines, (0, 0)) 

        