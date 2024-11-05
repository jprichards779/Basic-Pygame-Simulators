import math


"""
Instances of the the following Mass class hold all the data 
concerning the their whereabouts"""

class Mass:
    id=0                    
    distance_unit = 1       
    scale = 1000           
    def __init__(self,m=0,s=[0,0],v=[0,0], colour=(255,255,255), avg_density=1000):
        self.ID = Mass.id   
        Mass.id += 1
        self.assertions(m,v)
        self.colour = colour
        self.m, self.avg_density, self.s, self.v = m,avg_density,s,v
        self.initialise_data_structures()
    
    def initialise_data_structures(self):
        self.others, self.v_mag, self.r = [],[],[]
        self.r_mag, self.g, self.gR = [],[],[]
        self.screen_points = []
        self.dot_diameter, self.real_diameter = self.calc_sphere_diam()
        
    def assertions(self,m,v):
        if m > 10**32: m=10**32
        assert (v[0] and v[1]) < 2*10**8     
        assert Mass.distance_unit > 0 

    def calc_sphere_diam(self):
        enlarge = Mass.scale
        # Diameter of the object D
        D = 2*(3*self.m/(4*math.pi*self.avg_density))**(1/3)
        # Dot diameter d
        d = enlarge*2*(3*self.m/(4*math.pi*self.avg_density))**(1/3)/Mass.distance_unit
        if d < 1: d = 1
        return d, D 