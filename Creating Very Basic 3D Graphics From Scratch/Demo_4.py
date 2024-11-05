from World3 import * 

"""
This program shows the sheet class which has had it's vertices 
scaled to form a 3D parabola. It's possible to fluctuate the constant used to scale 
this. Similarly, wave like behaviour can be achieved with some careful thought. 

____Visit World3 module to see how the inputs are scaled to coordinate the motion___

Notice how the program speeds up after a few seconds.
This is because previous coordinates are being revisited, as we would expect with 
periodic motion and rotations. The program stores visited data in dictionaries
so that calculations aren't unnecessarily repeated. 

The initial slowing arises from the need to check if a data point exists or not to avoid 
overwriting. This is might be overkill following recent modifications I made. 
The idea was to stop the globally scoped dictionaries overwriting and/or getting too large 
because this visibly slowed the frame rate, defeating their whole purpose. 

"""

class App: 
    # Initialise our world 
    M = [[
          [0,15,-5], [0,0,0],  # World distance    # Object instance in world frame
          [0,0,0], [0,0,-20],  # Double rotation transform turns allowed objects about world origin
          [0,0,0]              # Final rotation transform spins allowed objects about their own axes in space
         ],           
         [
          [0,0,0], [0,0,-5]    # Camera distance  # Single orientation transform 
         ]] 
                  
    object_pos = M[0]
    camera_pos = M[1]
    # Display
    front_view = False
    plan_view = False
    # Object colours
    line_colour = (0,255,0)
    sheet_colour = (0,0,150)
    background_colour = (0,0,10)

    def __init__(self):
        self.run = True
    def main(self):
        ob, cam = self.object_pos, self.camera_pos
        wireframe, colours = False, True
        # wireframe, colours = False, False
        while self.run:
            for event in pygame.event.get():
                if event.type == pygame.QUIT: 
                    self.run=False
            World3(self).display(wireframe,colours,flatten=False)
            ob[2][0] += 1
            pygame.display.update() 
        pygame.quit()  
App().main()




