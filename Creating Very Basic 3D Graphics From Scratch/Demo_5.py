from World4 import *

"""
NESTED ROTATIONS
        >>>>>   PRESS KEY(S)...
        >>>>>   to switch between 
        >>>>>   wireframe and 
        >>>>>   opaque/coloured modes 
In this demo, more can be adjusted from the app class and makes it easier to play around 
with. Visit World4 module to see how the numerical inputs are scaled to coordinate the motion             
"""

class App: 
    # Initialise our world

    M = [[
          [0,20,-2], [0,0,0],  # World distance    # Object instance in world frame
          [0,0,15], [0,0,0],  # Constant world tilt is performed by keeping all other rotations in 3rd list
          [0,0,0]             # Final rotation transform spins allowed objects about their own axes in space
         ],           
         [
          [0,-5,0], [0,0,0]  # Camera distance  # Single orientation transform 
         ]] 
                
    object_pos = M[0]
    camera_pos = M[1]

    # Object colours
    World4.colour = (90,100,90)
    line_colour = (0,0,255)
    cylinder_colour = (10,50,220)
    sheet_colour = (0,50,110)
    background_colour = (0,0,10)

    # Display
    def __init__(self):
        self.front_view = False
        self.plan_view = False
        self.frames = 0
        self.run = True

    def wheel_move(self, ob):
        ob[1][0] += 1
    def spin_world(self, ob):
        # Horizontal rotation below, following vertical rotation 
        # operation ob[2][2] tilts the world on its axis.
        ob[3][0] += 1.5 
    def spin_wheel(self,ob):
        ob[4][1] += 2
    def camera_move(self, cam):
        cam[0][2] += 0.1
    def camera_orient(self, cam):
        cam[1][2] -= 0.1

    def main(self):
        ob, cam = self.object_pos, self.camera_pos
        wireframe, colours = False, False
        while self.run and abs(ob[3][0]) <220:
            self.frames += 1
            for event in pygame.event.get():
                if event.type == pygame.QUIT: 
                    self.run=False

                """ >>>>  PRESS KEY(S)  <<<<< """
                if event.type == KEYDOWN:
                    App.line_colour = (0,0,40)
                    colours=True
                elif event.type == KEYUP:
                    App.line_colour = (0,0,255)
                    colours=False

            World4(self).display(wireframe,colours,flatten=False)
            self.spin_world(ob)
            self.wheel_move(ob)
            self.spin_wheel(ob)

            # self.camera_move(cam)
            # self.camera_orient(cam)
            pygame.display.update() 
        pygame.quit()  
App().main()