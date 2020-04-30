from mpl_toolkits.mplot3d import axes3d
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import matplotlib.pyplot as plt
import numpy as np
import copy
import itertools

''' for drawing '''
class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


''' all positions are fractional '''

class atom:
    def __init__(self, pos):
        ## no one for atom, only zero
        self.pos = np.array(pos)
        self.moment = 0
        self.mdir = np.array([0,0,0])
    def set_moment(self, mvector):
        self.moment = np.linalg.norm(mvector)
        ## direction is normalized
        if self.moment == 0:
            self.mdir = np.array([0,0,0])
        else:
            self.mdir = np.array(mvector)/self.moment
        

class cell:
    def __init__(self, lattice_const):
        self.lattice_const = np.array(lattice_const)
        self.atoms = []
        self.current_num = 0

    def copy(self):
       new_cell = cell(self.lattice_const)
       new_cell.atoms = self.atoms
       return new_cell
    
    def insert_atom(self, new_atom):
        self.atoms.append(new_atom)
        self.current_num += 1

    def display_atom(self):
        for i in range(len(self.atoms)):
            print(self.atoms[i].pos)
    def plot(self, color_num = 0):
        fig = plt.figure()
        color = ['r','g','b','c','m','y','k']
        ax = fig.add_subplot(111, projection='3d')
        if self.atoms[0] is None:
            raise ValueError("NO atom inserted!")
            exit()
        else:
            print("total {} atoms".format(len(self.atoms)))
        ## generating a position matrix for vector plot
        vector_m = []
        for i in range(self.current_num):
            col = np.array([self.atoms[i].pos, self.atoms[i].mdir])
            vector_m.append(col.flatten())

        ##plot a cube frame
        r = [0, 1]
        for s, e in itertools.combinations(np.array(list(itertools.product(r, r, r))), 2):
            if np.sum(np.abs(s-e)) == r[1]-r[0]:
                ax.plot3D(*zip(s, e), color="b")
        
        ## plot moment and atom
        X, Y, Z, U, V, W = zip(*vector_m)
        ax.scatter(X, Y, Z, c=color[color_num])
        ax.quiver(X, Y, Z, U, V, W, length=.05)
        ax.set_xlim([-0.1, 1.1])
        ax.set_ylim([-0.1, 1.1])
        ax.set_zlim([-0.1, 1.1])
        plt.show()


def expand(old_cell, supercell_dim):
    supercell_dim = np.array(supercell_dim)
    supercell_size = np.prod(supercell_dim)
    new_lattice_const = old_cell.lattice_const* supercell_dim
    new_cell = cell(new_lattice_const)

    ## set up atoms in new supercell
    for old_atom in old_cell.atoms:
        old_pos = old_atom.pos/supercell_dim
        for i in range(supercell_dim[0]):
            for j in range(supercell_dim[1]):
                for k in range(supercell_dim[2]):
                    new_pos = old_pos + np.array([i,j,k])/supercell_dim
                    new_atom = atom(new_pos)
                    new_atom.set_moment(old_atom.mdir)
                    new_cell.insert_atom(new_atom)

    return new_cell

def centrosym(cell):
    ## notice the copy and deepcopy: if you change copy, you will also change reference; deepcopy does not
    centrol_cell = copy.deepcopy(cell)
    i = 0
    for i in range(len(cell.atoms)):
        aatom = cell.atoms[i]
        product = np.product(aatom.pos)
        tmp = np.zeros(3)
        if product == 0:
            zero_index = np.where(aatom.pos == 0)[0]
            N0 = len(zero_index)
            ## remove the first term with they are all zeros
            all_pos = [list(i) for i in itertools.product([0, 1], repeat=N0)][1:]
            for pos in all_pos:
                tmp = aatom.pos
                np.put(tmp, zero_index, pos)
                new_atom = atom(tmp)
                new_atom.set_moment(aatom.mdir)
                centrol_cell.insert_atom(new_atom)
    print("total {} atoms after centrosymmetry".format(len(centrol_cell.atoms)))
    return centrol_cell

def dipole_field(cell, endpoint):
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

    VecB = np.zeros(3)

    for i in range(len(cell.atoms)):
        aatom = cell.atoms[i]
        startpoint = aatom.pos
        VecR = (startpoint - endpoint)* cell.lattice_const* LenUnit ##in cm
        NormR = np.linalg.norm(VecR) 
        Moment = aatom.mdir* m ##in G

        #   B = 1/c*(3*r(m.r)/r^5 - m/r^3)
        dB = 3*VecR*np.dot(Moment, VecR)/(NormR**5) - Moment/(NormR**3)
        VecB += dB
        #print(dB)

    return VecB



tmvo4 = cell([7.07120, 7.07120, 6.26060])
tm1 = atom([0, 0, 0.5])
tm2 = atom([0.5, 0.5, 0])
tm3 = atom([0, 0.5, 0.75])
tm4 = atom([0.5, 0, 0.25])

mu = np.array([1, 0, 0])
## AFM between nearest neighbor
tm1.set_moment(mu)
tm2.set_moment(mu)
tm3.set_moment(-mu)
tm4.set_moment(-mu)

tmvo4.insert_atom(tm1)
tmvo4.insert_atom(tm2)
tmvo4.insert_atom(tm3)
tmvo4.insert_atom(tm4)
tmvo4.plot()

## V sites display
Vsite = cell([7.07120, 7.07120, 6.26060])
V1 = atom([0.5, 0.5, 0.5])
V2 = atom([0.5, 0, 0.75])
V3 = atom([0, 0.5, 0.25])
V4 = atom([0,0,0])

Vsite.insert_atom(V1)
Vsite.insert_atom(V2)
Vsite.insert_atom(V3)
Vsite.insert_atom(V4)

Vsite_sym = centrosym(Vsite)
#Vsite_sym.display_atom()
Vsite_sym.plot()

n_size = 3
n_array = np.array([1, 1, 1])* n_size
expand_cell = expand(tmvo4, n_array)
## translate sites to supercell center
target_site = [(x.pos + int((n_size - 1)/2))/n_size for x in Vsite_sym.atoms]
supercell = centrosym(expand_cell)
supercell.plot()

## result vector plot
vector_field = []
for vsite in target_site:
    vfield = dipole_field(supercell, vsite)
    print(vfield)
    col = np.array([vsite, vfield])
    vector_field.append(col.flatten())

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
cell_size = 1/(2*n_size)
##plot a cube frame
r = [0.5 - cell_size, 0.5 + cell_size]
for s, e in itertools.combinations(np.array(list(itertools.product(r, r, r))), 2):
    if np.sum(np.abs(s-e)) == r[1]-r[0]:
        ax.plot3D(*zip(s, e), color="b")

##plot axis arrows
arrow_prop_dict = dict(mutation_scale=20, arrowstyle='->',shrinkA=0, shrinkB=0)
a = Arrow3D([0.5 - cell_size, 0.5 - cell_size + 1/n_size], [0.5 - cell_size, 0.5 - cell_size], \
            [0.5 - cell_size, 0.5 - cell_size], **arrow_prop_dict, color='g')
ax.add_artist(a)
a = Arrow3D([0.5 - cell_size, 0.5 - cell_size], [0.5 - cell_size, 0.5 - cell_size + 1/n_size], \
            [0.5 - cell_size, 0.5 - cell_size], **arrow_prop_dict, color='g')
ax.add_artist(a)
a = Arrow3D([0.5 - cell_size, 0.5 - cell_size],[0.5 - cell_size, 0.5 - cell_size], \
            [0.5 - cell_size, 0.5 - cell_size + 1/n_size], **arrow_prop_dict, color='r')
ax.add_artist(a)

ax.text(0.5 - cell_size + 1.1/n_size, 0.5 - cell_size, 0.5 - cell_size, r'$a$')
ax.text(0.5 - cell_size, 0.5 - cell_size + 1.1/n_size, 0.5 - cell_size, r'$b$')
ax.text(0.5 - cell_size, 0.5 - cell_size, 0.5 - cell_size + 1.1/n_size, r'$c$')

## vector field
X, Y, Z, U, V, W = zip(*vector_field)
ax.quiver(X, Y, Z, U, V, W, length=1e-4)
ax.set_xlim([0.5 - cell_size, 0.5 + cell_size])
ax.set_ylim([0.5 - cell_size, 0.5 + cell_size])
ax.set_zlim([0.5 - cell_size, 0.5 + cell_size])
ax.axis("off")
plt.show()
    

                    
                
    
        
        
