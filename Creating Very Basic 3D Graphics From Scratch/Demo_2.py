from World2 import *

"""
This program shows different shapes which can be built in the Shapes module. 
We can see here how they can be rotated in various ways.

We have a blue, cylindrical pillar and a cuboid comprised of smaller cubes.
Notice that the pillar is rotating about it's own axis. 
It's axis is then being rotated.

____Visit World2 module to see how the inputs are scaled to coordinate the motion___

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
          [0,20,0.8], [-90,0,0],  # World distance    # Object instance in world frame
          [80,0,0],   [0,0,0],    # Double rotation transform turns allowed objects about world origin
          [0,0,0]                 # Final rotation transform spins allowed objects about their own axes in space
         ],           
         [
          [0,0,0], [10,0,0]      # Camera distance  # Single orientation transform 
         ]] 
                
    object_pos, camera_pos = M[0], M[1]
    # Display
    front_view, plan_view = False, False
    # Object colours
    line_colour=(0,100,255)
    background_colour = (0,0,0)
    
    def __init__(self):
        self.run = True
    def main(self):
        ob, cam = self.object_pos, self.camera_pos
        wireframe, colours = False, True
        while self.run and ob[1][0] < 100:
            for event in pygame.event.get():
                if event.type == pygame.QUIT: 
                    self.run=False
            World2(self).display(wireframe,colours,flatten=False)
            ob[1][0] += 1
            # ob[2][0] += 1 # rotates 
            ob[4][1] += 1
            pygame.display.update() 
        pygame.quit()  
App().main()











