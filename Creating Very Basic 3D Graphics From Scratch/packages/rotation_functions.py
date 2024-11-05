
import math
PI = math.pi
tilt_dict={}

"""
The following tilt function is designed for the reordering of rotation
opertations. 
These operations do not commute, which means that if you re-order the each 
operation you will get a different result. 

As an example, scroll down and  uncomment the print statements beneath the 
basis_rotation() function and notice whether our basis vector ends up pointing somewhere else
or not. 
This means we could design a rotation operation that allows users to rotate the shapes about 
any axis in theory. 
"""
def tilt(x=[], y=[], z=[], a=0,b=1,c=2, angle=0):
     indexing = [a,b,c]
     M = [x,y,z]
     assert len(x) == len(y) == len(z)
     data = (tuple(x), tuple(y), tuple(z), tuple(indexing), angle)
     if data not in tilt_dict:
        e1, e2, e3 = M[indexing[0]], M[indexing[1]], M[indexing[2]] 
        i = tuple(e1)
        j = tuple([round(e2[n]*math.cos(PI*angle/180) - e3[n]*math.sin(PI*angle/180),10)
                for n in range(len(e2))])
        k = tuple([round(e3[n]*math.cos(PI*angle/180) + e2[n]*math.sin(PI*angle/180), 10) 
                for n in range(len(e3))])
        M_Out = [i, j, k]
        result = M_Out[indexing[0]],M_Out[indexing[1]],M_Out[indexing[2]]
        tilt_dict[data] = tuple(result)
        return tuple(result)
     else: 
         return tilt_dict[data]
     
"""
The rotate() function below is used in the transform and camera class to carry out all transformations
Again, dictionaries are used to store data to assist with the performance of the program
abd reduce jitter. 
"""
ro_1_dict={}
def rotate(X=[],Y=[],Z=[], angles=[0,0,0]):
    angles = [round(n,2) for n in angles]
    data = (tuple(X), tuple(Y), tuple(Z), tuple(angles))
    if data not in ro_1_dict: # 
        x,y,z = tilt(X,Y,Z, 2,1,0, angles[0]) # XY
        x,y,z  = tilt(x,y,z,  1,0,2, angles[1]) # XZ  - half and half
        x,y,z  = tilt(x,y,z, 0,1,2, angles[2]) # YZ
        ro_1_dict[data] = x,y,z
    return ro_1_dict[data] 
# print(rotate1([1],[1],[1]))
# print(rotate1([1],[1],[1]))
# print(rotate1([1],[1],[1]))



def basis_rotation(element=[], angles=[[0,20,0], [0,20,0]]):
    B = []
    basis = [1,1,1]
    for n in angles:
        r = rotate([basis[0]],[basis[1]],[basis[2]],n)
        basis = [r[0][0],r[1][0],r[2][0]]
        B.append(basis)
    return B
# print(basis_rotation())
# print(basis_rotation(angles=[[30,20,0], [0,0,0]]))
# print(basis_rotation(angles=[[30,0,0], [0,20,0]])) # rotate vector sideway, rotate vetor up
# print(basis_rotation(angles=[[0,20,0], [30,0,0]]))


"""
_________________Personal Note_________________

Instead of using tranform matrices in industry standard, 4x4 format, I thought I'd leave them 
as I originally designed them. 
I have deliberately avoided help from the internet for the project besides python syntax,
for the benefits of my learning. 

I derived all the equations below myself as an exercise during lockdown and have left them in their original format to 
hopefully lend credibility to this in a world of online resources. 
"""

