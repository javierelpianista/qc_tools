import numpy as np

class Molecule:
    def __init__(self):
        # Atomic number of each atom in the molecule
        self.Z = []
        # Atomic symbols (or custom names) of each atom in the molecule
        self.atoms = []
        # x, y and z coordinates of each atom in the molecule
        self.coords = []
        # 

    @classmethod
    def fromatoms(self, atoms, coords):
        self.atoms = np.array(atoms)
        self.coords = np.array(coords)

    @classmethod
    def from_xyz_file(self, filename, unit = 'angs'):
        read_file = open(filename, 'r')
        n_set = False
        n = 0

        if unit == 'angs':
            scalef = angs_to_bohr
        elif unit == 'bohr':
            scalef = 1
        else:
            scalef = None

        for line in read_file.readlines():
            if n_set:
                if line.strip():
                    n+=1
                    data = line.strip().split()
                    self.atoms.append(data[0])
                    coords = np.array([float(data[i])*scalef for i in range(1,4)])
                    self.coords.append(coords)
            else:
                n_read = int(line.strip())
                n_set = True

        if n != n_read:
            print("Numbers don't match")

