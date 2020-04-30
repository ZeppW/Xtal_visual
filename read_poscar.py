#! /usr/bin/python
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np


eps = 1e-4

def norm(vec):
    norm=np.sqrt(vec[0]**2+vec[1]**2+vec[2]**2)
    return norm

def dot_product(vec1, vec2):
    dp = vec1[0]*vec2[0] + vec1[1]*vec2[1] + vec1[2]*vec2[2]
    return dp

def num_list(l):
    for i in range(len(l)):
        l[i]=float(l[i])
    return l

def GetPos():
    Input=open(pwd+"POSCAR","r")
    cont=Input.readlines()
    Input.close()
    i=0
#   Get lattice constant
    while len(cont[i].split()) != 1 or len(cont[i+1].split()) != 3:
        i+=1
    Lattice=cont[i+1:i+4]
    for a in range(len(Lattice)):
        Lattice[a]=num_list(Lattice[a].split())
#   Get atomes
    Name=cont[i+4].split()
    Num=list(map(int,cont[i+5].split()))
    i=i+7
    j=0
    tag=0
    POS=[]
    Atom=[]
#   Get positions
    for num in Num:
        POS.append([])
        for k in range(num):
            ar=num_list(cont[i+k].split())
            POS[j].append(ar)
            Atom.append(Name[tag])
        j+=1
        i+=Num[tag]
        tag+=1
    return Lattice,Atom,POS

def LatticePara(La):
    para=[]
    for vec in La:
        Norm = norm(vec)
        para.append(Norm)
    print("a = %f, b = %f, c=%f\n" %(para[0],para[1],para[2]))
    return para

def NewPos(POS,La):
    '''change position from fractal to length'''
    for i in range(len(POS)):
        for j in range(len(POS[i])):
            POS[i][j]=list(POS[i][j][0]*np.array(La[0])+POS[i][j][1]*np.array(La[1])+POS[i][j][2]*np.array(La[2]))
    return POS

def DrawStr(POS):
    f=plt.figure(1)
    color=['r','g','b','c','m','y','k']
    tag=0
    ax=f.add_subplot(111, projection='3d')
    for mat in POS:
        mat=np.array(mat)
        mat=np.transpose(mat)
        x=mat[0]
        y=mat[1]
        z=mat[2]
        ax.scatter(x,y,z,c=color[tag%len(color)])
        tag+=1
    plt.show()

def NNDist(AtoNum,POS,Atom,NN):
    '''return the name and position for NN atomes
    NN = 1: Nearest neighbor
    NN = 2: Next nearest neighbor
    '''
#   Calculate distance
    Dist = []
    POS = [item for sub in POS for item in sub]
    AtoName = Atom[AtoNum]
    AtoPos = POS[AtoNum]
    for i in range(len(POS)):
        if i != AtoNum:
            dist = np.sqrt((POS[i][0]-POS[AtoNum][0])**2+(POS[i][1]-POS[AtoNum][1])**2+(POS[i][2]-POS[AtoNum][2])**2)
            Dist.append([dist,i+1,Atom[i]])
#   Sort and select distance
    Dist = sorted(Dist, key = lambda x:x[0])
#   Remove same atomes & get min values
    min_dist = Dist[0][0]
    sec_min_dist = Dist[0][0]
    NNN_tag = 0;
    for i in range(len(Dist)):
        Name = Dist[i][2]
        if abs(sec_min_dist - min_dist) < eps:
            NNN_tag += 1;
            sec_min_dist = Dist[i][0]

#   for i in range(len(Dist)):
#       if Name == AtoName:
#           Dist.remove(Dist[i])
    sec_min_dist = Dist[NNN_tag][0]
    if NN == 1:
        NN_POS = [POS[Dist[index][1]-1] for index in range(len(Dist)) if (Dist[index][0] - min_dist) < eps]
    elif NN == 2:
        NN_POS = [POS[Dist[index][1]-1] for index in range(len(Dist)) if (Dist[index][0] - sec_min_dist) < eps]
    else:
        print("ERROR!")
    return AtoPos, NN_POS

def DipField(DipPos, DipDir, VecR):
    '''
    Calculate magnetic field created by magnetic moments
    at certain points
    '''
    g = 1
#   CGS units
    muB = 9.274e-21
    speedc = 1
    LenUnit = 1e-8
#   Hund's rule
    J = 1
    m = g*J*muB
    DipPos = np.array(DipPos)*LenUnit
    DipDir = np.array(DipDir)
    VecR = np.array(VecR)*LenUnit
    Mag = m*DipDir/norm(DipDir)
    VecL = VecR - DipPos
    NormL = norm(VecL)
#   B = 1/c*(2*r(m.r)/r^4 - m/r^2)
    VecB = 3*VecL*dot_product(Mag, VecL)/(NormL**5) - Mag/(NormL**3)
    return VecB

def field_calc(VPos, TmPOS):
    B = 0
    DipDir = [1, 0, 0]
    i = 0
    for pos in TmPOS:
        i += 1
        temp = DipField(pos, DipDir, VPos)
        B += temp
        print("{}-th atom, position: {}, magnetic field: {}\n".format(i, pos, temp))
    return B


def field_plot():
    DipPos = [0, 0, 0]
    DipDir = [1, 0, 0]
    len_unit = 1
    xx = np.linspace(-20,20,10)*len_unit
    yy = np.linspace(-20,20,10)*len_unit
    zz = np.linspace(-20,20,10)*len_unit
    xx = xx.tolist()
    yy = yy.tolist()
    zz = zz.tolist()
    VecRList = [[i, j, k] for i in xx for j in yy for k in zz]
    VecBList = []
    for v in VecRList:
        VecBList.append(DipField(DipPos, DipDir, v))
    xx, yy, zz = list(map(list, zip(*VecRList)))    
    Bx, By, Bz = list(map(list, zip(*VecBList)))
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.quiver(xx, yy, zz, Bx, By, Bz, length = 0.2)
    ax.set_xlabel('x')

    plt.show()

def main():
    Lattice,Atom,POS = GetPos()
    POS = NewPos(POS,Lattice)
#   para = LatticePara(Lattice)
    AtoNum = int(input("The index of atom you want to observe\n"))
    AtoNum = AtoNum - 1
    DipPos, NN_POS = NNDist(AtoNum,POS,Atom,1)
    DipPos, NNN_POS = NNDist(AtoNum,POS,Atom,2)
    print(NN_POS, NNN_POS)
    magB = field_calc(DipPos, POS[1])
#   magB = field_calc(DipPos, NN_POS) + field_calc(DipPos, NNN_POS)
    print("total magnetic field: {}\n".format(magB))
#   POS = NewPos(POS,Lattice)
#   DrawStr(POS)

if __name__ == "__main__":
   import sys
   import math
   import matplotlib
   from mpl_toolkits.mplot3d import axes3d
   import matplotlib.pyplot as plt
   import numpy as np
   pwd=input("Please input path of your file\n")
   main()
