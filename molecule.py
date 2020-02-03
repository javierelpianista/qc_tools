import constants as co
import numpy as np

class Molecule:
    def __init__(self):
        # Atomic number of each atom in the molecule
        self.Z = []
        # Atomic symbols (or custom names) of each atom in the molecule
        self.atoms = []
        # x, y and z coordinates of each atom in the molecule
        self.coords = []
        # Charge of the molecule
        self.charge = 0
        # Multiplicity
        self.multiplicity = 0

    @classmethod
    def from_atoms(cls, atoms, coords):
        molecule = Molecule()
        molecule.atoms = np.array(atoms)
        molecule.coords = np.array(coords)

        return(molecule)

    @classmethod
    def from_xyz_file(cls, filename, unit = 'angs', extended = True):
        molecule = Molecule()

        read_file = open(filename, 'r')
        n_set = False
        n = 0
    
        if unit == 'angs':
            scalef = co.angs_to_bohr
        elif unit == 'bohr':
            scalef = 1
        else:
            scalef = None

        for i, line in enumerate(read_file.readlines()):
            if i == 0:
                n_read = int(line.strip())
            elif i == 1:
                if extended:
                    try:
                        charge, multiplicity = [int(k) for k in line.strip().split()]
                    except(ValueError):
                        charge = 0
                        multiplicity = 0
                    molecule.charge = charge
                    molecule.multiplicity = multiplicity

            else:
                n+=1
                data = line.strip().split()
                molecule.atoms.append(data[0])
                coords = np.array([float(data[i])*scalef for i in range(1,4)])
                molecule.coords.append(coords)

        if n != n_read:
            print("Numbers don't match")
            raise

        return molecule

    def natoms(self):
        return len(self.atoms)

    def xyz(self, unit = 'angs', add_header = True):
        if add_header:
            string = '{:d}\n{:d}{:d}\n'.format(self.natoms(), self.charge, self.multiplicity)
        else:
            string = ''

        if unit == 'bohr':
            scale = 1
        elif unit == 'angs':
            scale = co.bohr_to_angs
        else:
            raise('Unit {:s} not known. Available options are bohr and angs.'.format(unit))

        for n in range(self.natoms()):
            string = string + '{:8s}{:15.7f}{:15.7f}{:15.7f}\n'. \
                format(self.atoms[n], *[scale*x for x in self.coords[n]])

        return string

if __name__ == "__main__":
    ammonia = Molecule.from_xyz_file('nh3.xyz')

