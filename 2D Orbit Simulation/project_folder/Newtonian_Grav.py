from mass import *
import helper_functions

"""
Gravitation class is responsible for calculating initial data provided by 
each Mass instance created in the Main class, then updating each mass's data 
structures for the next calculation. """

class Gravitation:
    time_step = 5000                 # Default 5000 seconds per frame 
    def __init__(self,main):
        assert Gravitation.time_step > 0 and Gravitation.time_step <= 10**4
        self.main = main
        assert self.main.TIME_LAPSE >=0 and self.main.TIME_LAPSE <= 1
        self.current_system = main.input
        self.initialise_data_structures()

    def initialise_data_structures(self):
        self.map, self.p_total = {}, []
        self.dT = Gravitation.time_step*self.main.TIME_LAPSE
        self.removed, self.rem_ids, self.new_ids, self.new_system = [],[],[],[]

    # 1. Creating a method which ientidies all mass instances surrounding the current 
    # mass. A dictionary / self.map contains {mass : surrounding masses} elements.
    def mass_network(self):
        dict = {}
        for ind in range(len(self.current_system)):
            if type(self.current_system[ind]) == Mass:
                vals = []
                for n in self.current_system:
                    if self.current_system[ind] != n:
                        vals.append(n)
                dict[self.current_system[ind]] = tuple(vals)
        self.map = dict
        
    # 2. The following method will itterate through this dictionary and update the 
    # Mass.others data structure, creating a 'gravitational network' of mass instances
    def get_neighbours(self):
        for n in self.map:
            others = []
            for i in self.map[n]:
                others.append(i)
            n.others=others

    # 3. 
    def r_vectors(self):
        for n in self.map:
            out = []
            if len(n.others)>0:
                for i in self.map[n]:
                    r = [-(n.s[j]-i.s[j]) for j in range(2)]
                    out.append(r)
            n.r = out

    # 4. 
    def R_mag(self):
        for n in self.current_system:
            result = []
            for i in n.r:
                dist = ((i[0])**2 + (i[1])**2)**0.5
                result.append(dist)
            n.r_mag = result

    # 5. 
    def g_vectors(self):
        for n in self.current_system:
            result = []
            if len(n.others) > 0:
                for k in range(len(n.r_mag)):
                    gx = round(self.main.G*n.others[k].m*n.r[k][0]/(n.r_mag[k]**3),10)
                    gy = round(self.main.G*n.others[k].m*n.r[k][1]/(n.r_mag[k]**3),10)
                    result.append([gx,gy])
            n.g = result

    # 6. 
    def resultant_g(self):
        for n in self.current_system:
            gR_x, gR_y = 0,0
            if len(n.others)>0:
                for i in n.g:
                    gR_x += i[0]
                    gR_y += i[1]
            n.gR = [gR_x, gR_y]

    # 7.
    def calc_velocity(self):
        for n in self.current_system:
            vx,vy = 0,0
            if len(n.gR) >0:
                if abs(vx) < 2*10**8:vx += n.v[0] + n.gR[0]*self.dT 
                if abs(vy) < 2*10**8: vy += n.v[1] + n.gR[1]*self.dT
                n.v = [vx,vy]
  
    # 8.
    def reposition(self):
        for n in self.current_system:
            if len(n.gR) >0: 
                sx,sy = 0,0
                sx += n.s[0] + n.v[0]*self.dT 
                sy += n.s[1] + n.v[1]*self.dT 
                n.s = [sx,sy]
                n.screen_points = helper_functions.pygame_array([n.s[0]], [n.s[1]], 
                                                                self.main.screen_width, 
                                                                self.main.screen_height)


    """
    The above methods summaries the essential vector equations underpinning orbital motion in classical 
    physics and are all invoked in method [5.] of Main. 
    But in reality we need to account for collision events as we don't just want stable orbits. 
    (It turns out pygame has methods designed to deal with collisions so this extra work was the cost of my ignorance)
    Without an effort to account for a masses internal gravity, the path equations 
    above will simply allow masses to get too close and experience assymptotic gravity, causing 
    objects to accelrate through one another / false gravitational slingshots. The following resolve this.
    """

    # 9. The following method will allow the the machine to differentiate between 
    #    collided masses and other masses. 
    def remove_collided(self):
        removed,new = [],[]
        for n in self.current_system:
            for distance in n.r_mag:
                ind = n.r_mag.index(distance)
                size_other = n.others[ind].real_diameter
                total_dist = n.real_diameter + size_other
                """
                Modification I)
                Originally calculated closing velocity magnitude without resulant gravity magnitudes.
                This modification doesn't make a great deal of apparent difference but just adds to the 
                function's mass removal capability"""
                gRx, gRy= 0,0
                if (len(n.gR) and len(n.others[ind].gR)) == 2:
                    gRx = abs(n.gR[0]) + abs(n.others[ind].gR[0])
                    gRy = abs(n.gR[1]) + abs(n.others[ind].gR[1])
                vf_x = abs(n.v[0]) + abs(n.others[ind].v[0]) + gRx*self.dT   # v = v0 + g*t
                vf_y = abs(n.v[1]) + abs(n.others[ind].v[1]) + gRy*self.dT
                vf_mag = (vf_x**2 +vf_y**2)**0.5                             # Closing velocity
                del_threshold = (0.5*n.real_diameter + 0.5*size_other)/total_dist
                LIMIT = del_threshold + vf_mag*self.dT                       # s = s0 + v*t
                if distance <= LIMIT: 
                    removed.append(n) 
            if n in removed: pass 
            else: new.append(n)   
        return removed, new
 
    # 10.
    # The physical attrubutes of the removed masses to calculate the those of a new mass
    # which will replace the removed ones ready for the next itteration. 
    # The new mass will appear on screen as the resulting outcome of the collided masses
    def assymilate(self):
        removed, new = self.remove_collided()
        rem_pos_x, rem_pos_y = [n.s[0] for n in removed], [n.s[1] for n in removed]
        if len(removed) > 0:
            masses_collected = [n.m for n in removed]
            ind = masses_collected.index(max(masses_collected))
            self.substitute_colour = removed[ind].colour
            m_final = 0
            density = 0
            px,py = 0,0
            sx,sy = 0,0
            for n in removed:
                m_final += n.m
                density += n.avg_density*n.m 
                px += n.m*n.v[0]
                py += n.m*n.v[1]
                for i in n.others:
                    if i in removed: 
                        mass_ratio = i.m/(n.m+i.m) 
                        sx += n.s[0]*(1-mass_ratio)   
                        sy += n.s[1]*(1-mass_ratio)
                p_final= [px,py]
                v_final = [n/m_final for n in p_final]  
                avg_density = density/m_final
                if sx >=min(rem_pos_x) and sx<max(rem_pos_x):
                    if sy >=min(rem_pos_y) and sy<max(rem_pos_y):
                        s_final = [sx,sy]
                        new.append(Mass(m=m_final, s=s_final, v=v_final, colour = self.substitute_colour, 
                                        avg_density = avg_density))
                        for n in removed: 
                            if n.ID == self.main.center_object_ID:
                                new[-1].ID = self.main.center_object_ID
                        self.rem_ids = [n.ID for n in removed]
                        self.new_ids = [n.ID for n in new]
                        self.current_system = new
