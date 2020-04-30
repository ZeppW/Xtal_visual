from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np


''' all positions are fractional '''
class atom:
    def __init__(self, pos):
        self.pos = np.array(pos)
        self.moment = 0
        self.mdir = np.array([0,0,0])
    def set_moment(self, mvector):
        self.moment = np.linalg.norm(mvector)
        self.mdir = np.array(mvector)/self.moment
        

class cell:
    def __init__(self, lattice_const, atom_num):
        self.lattice_const = np.array(lattice_const)
        self.atom_num = atom_num
        self.atoms = np.empty(self.atom_num, dtype = object)
        self.current_num = 0
    def insert_atom(self, new_atom):
        self.atoms[self.current_num] = new_atom
        self.current_num = self.current_num + 1
    def plot(self, color_num = 0):
        fig = plt.figure()
        color = ['r','g','b','c','m','y','k']
        ax = fig.add_subplot(111, projection='3d')
        if self.atoms[0] is None:
            raise ValueError("NO atom inserted!")
            exit()
        ## generating a position matrix for vector plot
        vector_m = []
        for i in range(self.atom_num):
            col = np.array([self.atoms[i].pos, self.atoms[i].mdir])
            vector_m.append(col.flatten())

        X, Y, Z, U, V, W = zip(*vector_m)
        ax.scatter(X, Y, Z, c=color[color_num])
        ax.quiver(X, Y, Z, U, V, W, length=.05)
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])
        ax.set_zlim([0, 1])
        plt.show()


def expand(old_cell, supercell_dim):
    supercell_dim = np.array(supercell_dim)
    supercell_size = np.prod(supercell_dim)
    new_lattice_const = old_cell.lattice_const* supercell_dim
    new_atom_num = old_cell.atom_num* supercell_size
    new_cell = cell(new_lattice_const, new_atom_num)

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
    centrol_cell = cell
    for aatom in cell.atoms:
        if aatom.pos[0] == 0 or aatom.pos[0] == 1:
            new_pos = [1 - aatom.pos[0], aatom.pos[1], aatom.pos[2]]
            new_atom = atom(new_pos)
            new_atom.set_moment(aatom.mdir)
            centrol_cell.insert_atom(new_atom)


        if aatom.pos[1] == 0 or aatom.pos[1] == 1:
            new_pos = [aatom.pos[0], 1 - aatom.pos[1], aatom.pos[2]]
            new_atom = atom(new_pos)
            new_atom.set_moment(aatom.mdir)
            centrol_cell.insert_atom(new_atom)

        if aatom.pos[2] == 0 or aatom.pos[2] == 1:
            new_pos = [aatom.pos[0], aatom.pos[1], 1 - aatom.pos[2]]
            new_atom = atom(new_pos)
            new_atom.set_moment(aatom.mdir)
            centrol_cell.insert_atom(new_atom)

    return centrol_cell


test_cell = cell([10,10,10], 3)
atom1 = atom([0,0,0.5])
atom1.set_moment([0,0,1])
atom2 = atom([0,0.5,0])
atom2.set_moment([0,0,-1])
atom3 = atom([0.5,0,0])
atom3.set_moment([0,0,1])
test_cell.insert_atom(atom1)
test_cell.insert_atom(atom2)
test_cell.insert_atom(atom3)
test_cell.plot()
expand_cell = expand(test_cell, [3,3,3])
expand_cell.plot()
centrol_cell = centrosym(expand_cell)
centrol_cell.plot()
                    
                
    
        
        
