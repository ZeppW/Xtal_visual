import heapq as pq
import numpy as np

class crsytal_str:
    global pwd
    pwd = input("Please input path of your file\n")
    def __init__(self):
        tag = input("POSCAR input 1;CONTCAR input 0\n")
        if tag:
            Input=open(pwd+"POSCAR","r")
        else:
            Input=open(pwd+"CONTCAR","r")

        self.cont=Input.readlines()
        Input.close()

        
        self.lattice = None
        self.laconst = []
        self.pos=[]
        self.name=[]

    #   Get lattice vectors
        rows=0
        while len(self.cont[rows].split()) != 1 or len(self.cont[rows+1].split()) != 3:
            rows+=1
        self.lattice = self.cont[rows+1:rows+4]

     #   Get lattice constant after lattice vectors are obtained       
        for i in range(3):
            self.lattice[i]=np.array(self.lattice[i].split()).astype(float)
            Norm = np.linalg.norm(self.lattice[i])
            self.laconst.append(Norm)

    #   Get atomes
        AtomName=self.cont[rows+4].split()
        AtomNum=list(map(int,self.cont[rows+5].split()))
        rows=rows+7
        level=0
    #   Get positions
        for num in AtomNum:
            for k in range(num):
                ar=list(map(float,self.cont[rows+k].split()))
                pq.heappush(self.pos, ar)
                self.name.append(AtomName[level])
            rows+=AtomNum[level]
            level+=1

#    def 
    
#    def cartpos:

#    def fractional:

#    def distance:

#    def neighbor:

#    def dipfield:

#    def drawstr:

#    def drawfield:
