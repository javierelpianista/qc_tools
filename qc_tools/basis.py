import qc_tools.constants as co
import numpy as np

# A simple shell, containing one or more contracted primitives of the same angular momentum
class Shell:
    def __init__(self):
        self.exponents    = []
        self.coefficients = []
        self.l = None

    # Check if it's a valid shell
    def check(self):
        if len(self.exponents) != len(self.coefficients) or len(self.exponents) == 0:
            raise Exception('In Shell.check; the number of exponents should be equal to the number of coefficients')
        if not isinstance(self.l, int):
            raise Exception('In Shell.check; the angular momentum l should be an integer')
        else:
            if not self.l >= 0:
                raise Exception('In Shell.check; the angular momentum l should be an integer equal or larger than 0.')

    @classmethod
    def from_vectors(cls, expos, coefs, l):
        shell = Shell()
        shell.exponents = expos
        shell.coefficients = coefs
        shell.l = l

        return(shell)

    def gbs_string(self):
        string = ''
        string += '{:2}{:4d}{:8.3f}\n'.format(co.angular_momenta[self.l], len(self.exponents), 1.00)

        for expo, coef in zip(self.exponents, self.coefficients):
            string += '{:18.8e}{:18.8e}\n'.format(expo, coef)

        return(string)
    
class BasisSet:
    def __init__(self): 
        self.atoms  = []
        self.shells = {}

    @classmethod
    def from_gbs_file(cls, filename):
        bset = BasisSet()

        # First, preprocess the basis file: remove comments and empty lines and separate text fragments belonging to each atom
        file_fragments = []
        fragment = []
        for line in open(filename, 'r'):
            liner = line.strip().split('!')[0]
            if liner == '****':
                file_fragments.append(fragment)
                fragment = []
            else: 
                if len(liner) > 0:
                    fragment.append(liner)

        # Now process each fragment separately
        for fragment in file_fragments:
            if len(fragment) == 0:
                continue

            for n, line in enumerate(fragment):
                if n == 0:
                    atom, charge = line.split()
                    if atom in bset.atoms:
                        raise Exception('ERROR! In BasisSet.from_gbs_file. Repeated atom in file {}'.format(filename))

                    bset.atoms.append(atom)
                    #bset.shells[atom] = []
                    shell_list = []
                    read_shell = True
                    nshells = 0

                else:
                    if read_shell:
                        nshells += 1
                        data = line.split()
                        shell_type = data[0]
                        nprim = int(data[1])
                        scf = float(data[2].replace('D','E'))
                        if scf != 1.00:
                            raise Warning('WARNING!!! In BasisSet.from_gbs_file: the scaling factor for shell {} of atom {} is not equal to 1. Not sure what to do here.'.format(nshells, atom))
                        read_shell = False
                        iprim = 0
                        shell = Shell()

                    else:
                        iprim += 1
                        data = line.split()
                        expo = float(data[0].replace('D','E'))
                        coef = float(data[1].replace('D','E'))
                        shell.exponents.append(expo)
                        shell.coefficients.append(coef)
                        shell.l = co.angular_momenta.index(shell_type)

                        if iprim == nprim:
                            shell.check()
                            shell_list.append(shell)
                            read_shell = True

            bset.shells[atom] = sorted(shell_list, key = lambda shell: shell.l)

        return(bset)

    def atom_shell_gbs(self, atoms):
        if not isinstance(atoms, list):
            atoms = [atoms]

        string = ''
        for atom in atoms:
            if atom not in self.atoms:
                raise Exception('ERROR! In BasisSet.atom_shell_gbs: atom {} not included in the basis set.'.format(atom))

            string += '{:4}{}\n'.format(atom, 0)
            for shell in self.shells[atom]:
                string += shell.gbs_string()

        return(string)

    def atom_shell_dalton(self, atom):
        if atom not in self.atoms:
            raise Exception('ERROR! In BasisSet.atom_shell_daltom: atom {} not included in the basis set.'.format(atom))

        max_l = self.atom_max_l(atom)

        # Initialize one set of exponents for each l
        all_expos = []
        for l in range(max_l+1):
            all_expos.append(set())

        # Add each unique exponent to each shell
        for shell in self.shells[atom]:
            for expo in shell.exponents:
                all_expos[shell.l].add(expo)

        # Now produce the corresponding string
        string = ''
        for l in range(max_l+1):
            curr_shells = list(filter(lambda x: x.l == l, self.shells[atom]))
            n_curr_shells = len(curr_shells)

            string += '{:1}{:4}{:5}\n'.format('H', len(all_expos[l]), n_curr_shells)
            for expo in sorted(all_expos[l], reverse=True):

                coefs = np.zeros(n_curr_shells)

                for k, shell in enumerate(curr_shells):
                    for n, expo_sh in enumerate(shell.exponents):
                        if expo == expo_sh:
                            coefs[k] = shell.coefficients[n]

                if n_curr_shells > 3:
                    string += '{:20.10f}'.format(expo)
                    for n in range(3):
                        string += '{:20.10f}'.format(coefs[n])
                    string += '\n' + 20*' '
                    for n in range(3, n_curr_shells):
                        string += '{:20.10f}'.format(coefs[n])

                    string += '\n'

                else:
                    string += '{:20.10f}'.format(expo)
                    for n in range(n_curr_shells):
                        string += '{:20.10f}'.format(coefs[n])
                    string += '\n'

        return(string)

    def atom_max_l(self, atom):
        max_l = 0
        for shell in self.shells[atom]:
            if shell.l > max_l: max_l = shell.l

        return(max_l)

    def to_gbs(self, filename=None):
        string = ''
        for atom in self.atoms:
            string += self.atom_shell_gbs(atom)
            string += '****\n'

        if filename:
            open(filename, 'w').write(string)
        return(string)

    def reduced(self, level):
        other = BasisSet()

        for atom in self.atoms:
            other.atoms.append(atom)
            other.shells[atom] = []
            max_l = self.atom_max_l(atom)

            for shell in self.shells[atom]:
                if shell.l <= max(0, max_l - level):
                    other.shells[atom].append(shell)

        return(other)

