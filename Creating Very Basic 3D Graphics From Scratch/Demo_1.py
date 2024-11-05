from World1 import *

"""
This program shows different shapes which can be built in the Shapes module. 
We can see here how they can be rotated in various ways.

We have a blue, cylindrical pillar and a cuboid comprised of smaller cubes.
Notice that the pillar is rotating about it's own axis. 
It's axis is then being rotated.

____Visit World1 module to see how the inputs are scaled to coordinate the motion___

Notice how the program speeds up after a few seconds. This is because the program stores 
visited data in dictionaries so that calculations aren't unnecessarily repeated. 

The initial slowing arises from the need to check if a data point exists or not, to avoid 
overwriting. This is might be overkill following recent modifications I made. 
The idea was to stop the globally scoped dictionaries overwriting and/or getting too large 
because this visibly slowed the frame rate, defeating their whole purpose. 

"""


class App: 
    # Initialise our world
    M = [[
          [1,3,0], [0,0,0],  # World distance    # Object instance in world frame
          [0,0,0], [0,0,0],  # Double rotation transform turns allowed objects about world origin
          [0,0,0]            # Final rotation transform spins allowed objects about their own axes in space
         ],           
         [
          [0,-2,0], [0,0,0]  # Camera distance  # Single orientation transform 
         ]] 
                
    object_pos, camera_pos = M[0], M[1]
    # Display
    front_view, plan_view = False, False
    # Object colours       
    line_colour = (200,200,200) # adjusting line colour 
    background_colour = (0,0,10)
    def __init__(self):
        self.run = True

    def main(self):
        ob, cam = self.object_pos, self.camera_pos
        # Adjust the following boolean variables as desired
        wireframe, colours = False, True
        while self.run:
            for event in pygame.event.get():
                if event.type == pygame.QUIT: 
                    self.run=False
            # >>>> PRESS A KEY 
            # >>>> on your machine:
                if event.type == KEYDOWN:
                    wireframe = True
                    colours = False
                elif event.type == KEYUP:
                    wireframe = False
                    colours = True


            # World's display function calls Screen class to render shapes
            World(self).display(wireframe,colours,flatten=False)
            # World rotation
            ob[3][0] += 1 
            # Set object spinning 
            ob[4][1] += 2 
            pygame.display.update() 
        pygame.quit()  
App().main()

