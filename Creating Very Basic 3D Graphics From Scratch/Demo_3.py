from World2 import *

"""
Here we can see that we can dictate whether the world is coloured in 
or not.

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
          [1,15,-1], [0,0,0],   # World distance    # Object instance in world frame
          [0,0,0], [90,0,0],    # Double rotation transform turns allowed objects about world origin
          [0,0,0]               # Final rotation transform spins allowed objects about their own axes in space
         ],           
         [
          [0,0,0], [10,0,-5]    # Camera distance  # Single orientation transform 
         ]] 
                

    object_pos, camera_pos = M[0], M[1]
    # Display
    front_view, plan_view = False, False
    # Object colours
    line_colour=(0,100,255)
    background_colour = (0,0,10)
    
    def __init__(self):
        self.run = True

    def main(self):
        ob, cam = self.object_pos, self.camera_pos
        # The following line dictates whether we want transparency
        wireframe, colours = False, False
        # ....or shapes to be coloured acoording to what we've specified. 
        # Play around witht he parameters to see what happens. 
        while self.run and ob[1][0] < 60:
            for event in pygame.event.get():
                if event.type == pygame.QUIT: 
                    self.run=False
            World2(self).display(wireframe,colours,flatten=False)
            ob[1][0] += 0.2
            ob[4][1] -= 2
            pygame.display.update() 
        pygame.quit()  
App().main()

