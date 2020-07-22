import qc_tools.constants as co
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
        self.multiplicity = 1

    @classmethod
    def from_atoms(cls, atoms, coords):
        molecule = Molecule()
        molecule.atoms = np.array(atoms)
        molecule.coords = np.array(coords)

        return(molecule)

    @classmethod
    def from_molecules(cls, mol1, mol2):
        import warnings

        molecule = Molecule()
        molecule.atoms  = np.array(mol1.atoms + mol2.atoms) 
        molecule.coords = np.array(mol1.coords + mol2.coords) 
        molecule.Z      = np.array(mol1.Z + mol2.Z)
        molecule.charge = mol1.charge + mol2.charge

        if mol1.multiplicity == 1 and mol2.multiplicity == 1:
            molecule.multiplicity = 1
        else:
            warnings.warn('WARNING!!! I don''t know what to do to add the multiplicity of two molecules.' +
                    'I set them to 1.')

            molecule.multiplicity = 1

        return(molecule)

    @classmethod
    def from_xyz_file(cls, filename, unit = 'angs', extended = True, charge = 0, multiplicity = 1):
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

        molecule.charge       = charge
        molecule.multiplicity = multiplicity

        for i, line in enumerate(read_file.readlines()):
            if i == 0:
                n_read = int(line.strip())
            elif i == 1:
                if extended:
                    try:
                        charge, multiplicity = [int(k) for k in line.strip().split()]
                        molecule.charge = charge
                        molecule.multiplicity = multiplicity
                    except(ValueError):
                        pass
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
    
    @classmethod
    def from_cml_file(cls, filename, unit = 'angs', extended = True, charge = 0, multiplicity = 1):
        elements = []
        coords   = []

        read = False
        for line in open(filename, 'r').readlines():
            if '</atomArray>' in line:
                read = False

            if read:
                data = line.split()
                
                # Read element type
                info = data[2].split('"')
                elements.append(info[1])

                # Read coordinates
                coord = []
                for n in range(3,6):
                    info = data[n].split('"')
                    coord.append(float(info[1]))
                coords.append(coord)

            if '<atomArray>' in line: 
                read = True

        molecule = Molecule.from_atoms(elements, coords)
        return(molecule)

    def natoms(self):
        return len(self.atoms)

    def xyz(self, unit = 'angs', add_header = True):
        if add_header:
            string = '{:d}\n{:d} {:d}\n'.format(self.natoms(), self.charge, self.multiplicity)
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

    def __add__(self, other):
        if not isinstance(other, Molecule):
            raise Exception('Expected other Molecule in __add__ operator')

        result = Molecule()

        for atom, coords in zip(self.atoms + other.atoms, self.coords + other.coords):
            result.atoms.append(atom)
            result.coords.append(coords)

        result.charge = self.charge + other.charge

        if self.multiplicity != 1 or other.multiplicity != 1:
            print('One of the multiplicities is different from 1. Please set it manually.')

        result.multiplicity = 1

        return(result)

# Translate a molecule by a set amount. vec should be a list containing the x, y and z values.
def translate(mol, vec, unit = 'angs'):
    from copy import deepcopy

    if unit == 'angs':
        scalef = co.angs_to_bohr
    elif unit == 'bohr':
        scalef = 1
    else:
        scalef = None

    mol2 = deepcopy(mol)

    for n in range(mol2.natoms()):
        for m in range(3):
            mol2.coords[n][m] += vec[m]*scalef

    return(mol2)

# Return the center of mass of a molecule
def COM(mol, unit = 'angs'):
    M = 0.
    pos = np.zeros(3)

    if unit == 'angs':
        scalef = co.bohr_to_angs
    elif unit == 'bohr':
        scalef = 1
    else:
        scalef = None

    for atom, coords in zip(mol.atoms, mol.coords):
        mass = co.atomic_masses[co.atomic_symbols.index(atom[0].upper() + atom[1].lower() if len(atom) == 2 else atom[0].upper())]

        M += mass
        pos += mass*np.array(coords)

    return(scalef*pos/M)

# Give the coordinates of a midbond center located in the 1/r^6-weighted average of the atomic positions of each monomer (this should be moved to a different file later on)
def mb_r6(monoA, monoB):
    weight = np.zeros((len(monoA.atoms), len(monoB.atoms)))
    coords = np.zeros(3)

    dp_sum = float(0.0)
    for i in range(len(monoA.atoms)):
        for j in range(len(monoB.atoms)):
            rab = 0.0
            for k in range(3):
                rab = rab + (monoA.coords[i][k] - monoB.coords[j][k])**2

            weight[i][j] = 1/rab**3
            dp_sum = dp_sum + weight[i][j]

    weight = weight / dp_sum

    for k in range(3):
        for i in range(len(monoA.atoms)):
            for j in range(len(monoB.atoms)):
                coords[k] = coords[k] + weight[i][j]*0.5*(monoA.coords[i][k] + monoB.coords[j][k])
        
    return coords

