"""Shape Colours"""
class Colours:
    BLACK = (0,0,0)
    WHITE = (255, 255, 255)
    OFF_WHITE = (200,200,200)
    RED = (100,0,30)
    MAGENTA = (255,0,255)
    GREEN = (0,150,0)
    BLUE = (50,0,255)
    DARK_BLUE = (0,0,100)
    DARK_SEA = (0,10,60)
    GRAY = (50,50,60) 
    def __init__(self,colour):
        self.colour = colour


"""
* tuples may be preferable for optimal 
performance using dictionaries. 
"""
def cubic(L=[0,0,0], 
                 vt=[[1, 1, 1] for n in range(4)], 
                 vb=[[1, 1, 1] for n in range(4)]):
    """Vt and vb contain x,y,z scale factors"""
    u, v, w = -0.5,0.5,0.5
    V1,V2,V3,V4,V5,V6,V7,V8 = vb[0],vb[1],vb[2],vb[3],vt[0],vt[1],vt[2],vt[3]     
    p1 = [u*V1[0]+L[0], u*V1[1]+L[1], w*V1[2]+L[2]] 
    p2 = [u*V2[0]+L[0], v*V2[1]+L[1], w*V2[2]+L[2]]
    p3 = [v*V3[0]+L[0], v*V3[1]+L[1], w*V3[2]+L[2]]
    p4 = [v*V4[0]+L[0], u*V4[1]+L[1], w*V4[2]+L[2]]
    p5 = [u*V5[0]+L[0], u*V5[1]+L[1], -w*V5[2]+L[2]] 
    p6 = [u*V6[0]+L[0], v*V6[1]+L[1], -w*V6[2]+L[2]]
    p7 = [v*V7[0]+L[0], v*V7[1]+L[1], -w*V7[2]+L[2]]
    p8 = [v*V8[0]+L[0], u*V8[1]+L[1], -w*V8[2]+L[2]]
    vertices = [p1,p2,p3,p4,p5,p6,p7,p8]
    upper,lower = [p1,p2,p3,p4], [p5,p6,p7,p8]
    M = [[p1,p2,p3,p4,p1], # top
         [p5,p6,p7,p8,p5], # base
         [p1,p2,p6,p5,p1], # side - left
         [p4,p3,p7,p8,p4], # side - right *** Elements swapped to register duplicates
         [p2,p3,p7,p6,p2], # side - back
         [p1,p4,p8,p5,p1]] # side - front *** Elements swapped to register duplicates
    return M
# print(quadrahedron())
# print(len(uniqueSides_only([quadrahedron(L=[0,0,0]), quadrahedron(L=[1,0,0])])))


def square(L=[0,0,0], vt=[[1, 1, 1] for n in range(4)], 
                            vb=[[1, 1, 1] for n in range(4)]):
    return cubic(L,vt,vb)[5:6]


def normalised_angle(angle):
    if angle==0:
        return 0.0
    elif angle!=0:
        abs_val = abs(angle)
        sign = angle/abs_val
        M = abs_val//360
        if abs_val >= 360:
                abs_val-=M*360
        if abs_val != 0:
                return sign*abs_val
        else: return float(abs_val)

def norm_angle_list(l=[0,0,0]):
    return [normalised_angle(l[i]) for i in range(len(l))]


def x_y_z(M):
    x = [n[i][0] for n in M for i in range(len(n))]
    y = [n[i][1] for n in M for i in range(len(n))]
    z = [n[i][2] for n in M for i in range(len(n))]
    return x, y, z


    
def select_sides(M=cubic(), input= "top, base, left, right, back, front"):
    out = []
    l = ""
    for n in input:
        if not n.isalpha():
            l+= " "
        else: l+=n.lower()
    m = l.split()
    dict = {"top": M[0],
            "base": M[1],
            "left": M[2],
            "right": M[3],
            "back": M[4],
            "front": M[5]}
    for i in m:
        if i not in dict:
            ind = i.index(m)
            print(f"Input error @ {ind}")
        else: out.append(dict.get(i))
    return out


def cubic2(L=[0,0,0],  
                  vt=[[1, 1, 1] for n in range(4)], 
                  vb=[[1, 1, 1] for n in range(4)],
                  sides="top, base, left, right, back, front"):
    u, v, w = -0.5,0.5,0.5
    V1,V2,V3,V4,V5,V6,V7,V8 = vb[0],vb[1],vb[2],vb[3],vt[0],vt[1],vt[2],vt[3] 
    p1 = (u*V1[0]+L[0], u*V1[1]+L[1], w*V1[2]-L[2])
    p2 = (u*V2[0]+L[0], v*V2[1]+L[1], w*V2[2]-L[2])
    p3 = (v*V3[0]+L[0], v*V3[1]+L[1], w*V3[2]-L[2])
    p4 = (v*V4[0]+L[0], u*V4[1]+L[1], w*V4[2]-L[2])
    p5 = (u*V5[0]+L[0], u*V5[1]+L[1], -w*V5[2]-L[2])
    p6 = (u*V6[0]+L[0], v*V6[1]+L[1], -w*V6[2]-L[2])
    p7 = (v*V7[0]+L[0], v*V7[1]+L[1], -w*V7[2]-L[2])
    p8 = (v*V8[0]+L[0], u*V8[1]+L[1], -w*V8[2]-L[2])
    M = ((p1,p2,p3,p4,p1), # top
         (p5,p6,p7,p8,p5), # base
         (p1,p2,p6,p5,p1), # side - left
         (p4,p3,p7,p8,p4), # side - right *** Elements swapped to register duplicates
         (p2,p3,p7,p6,p2), # side - back
         (p1,p4,p8,p5,p1)) # side - front *** Elements swapped to register duplicates
    
    return select_sides(M=M, input=sides)
# print(len(uniqueSides_only([quadrahedron2(L=[0,0,0]), quadrahedron2(L=[1,0,0])])))



