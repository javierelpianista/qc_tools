import os
import numpy as np
from qc_tools.molecule import Molecule
import qc_tools.constants as co

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

# Create a sapt-fastdf input file
def fastdf_run(opts, fastdf_options=None, monoA=None, monoB=None, mb=None):
    # Check for mandatory options

    if 'basis_set' not in opts:
        raise Exception('Error! basis_set should be specified in the options')

    if 'reduction_level' not in opts:
        reduction_level = 0
    else:
        reduction_level = opts['reduction_level']

    # The default method is saptdft
    if 'method' in opts:
        method = opts['method']
    else:
        method = 'saptdft'

    # The default functional is PBE
    if 'functional' in opts:
        functional = opts['functional']
    else:
        functional = 'PBE'

    # Check if ionization potentials are defined; if so, ask for grac correction.
    # If not, set it to false.
    if 'ipA' in opts and 'ipB' in opts:
        ipa = opts['ipA'] 
        ipb = opts['ipB']
        grac = True
    elif 'ips' in opts:
        ipa = opts['ips']
        ipb = opts['ips']
        grac = True
    else:
        grac = False

    # Handle the way monomers are taken into account
    if monoA == None:
        if 'monoA_file' in opts:
            monoA = Molecule.from_xyz_file(opts['monoA_file'])
        else:
            if os.path.isfile('A.xyz'):
                monoA = Molecule.from_xyz_file('A.xyz')
            else:
                raise Exception('Monomer A hasn\'t been defined, the monoA_file option either, and there\'s no file named A.xyz')

    if monoB == None:
        if 'monoB_file' in opts:
            monoB = Molecule.from_xyz_file(opts['monoB_file'])
        else:
            if os.path.isfile('B.xyz'):
                monoB = Molecule.from_xyz_file('B.xyz')
            else:
                raise Exception('Monomer B hasn''t been defined, the monoB_file option either, and there''s no file named B.xyz')

    # Handle the midbond functions
    midbond = mb
    if 'mb_r6' in opts:
        if mb == None:
            coords = mb_r6(monoA, monoB)
            midbond = Molecule.from_atoms([opts['mb_r6']], [coords])

    # ------------------------------------------------------------------------------
    # Write the string that will be the .sapt file
    out_str = 'method {\n'

    out_str += 4*' ' + '{}\n'.format(opts['method'].lower())

    if 'memory' in opts:
        out_str += 4*' ' + 'memory ' + opts['memory'] + '\n'

    out_str += 4*' ' + 'functional {}\n'.format(functional)

    if grac:
        out_str += 4*' ' + 'grac {} {}\n'.format(ipa, ipb)
    else:
        out_str += 4*' ' + 'grac false\n'

    if 'corrections' in opts and opts['corrections'] != None:
        out_str += 4*' ' + 'corrections'
        for corr in opts['corrections']:
            out_str += ' ' + corr

        out_str += '\n'

    out_str += '}\n\n'

    # Basis section
    out_str += 'basis {\n'
    out_str += 4*' ' + 'main_basis_set {}\n'.format(opts['basis_set'])
    if reduction_level != 0:
        out_str += 4*' ' + 'off_basis_set reduce {:d}\n'.format(reduction_level)
    out_str += '}\n\n'

    # Monomer A section
    out_str += 'monomer {\n' + 4*' '
    out_str += monoA.xyz(add_header = False).strip().replace('\n', '\n    ') + '\n'
    if monoA.charge != 0:
        out_str += 4*' ' + 'charge {}\n'.format(monoA.charge)
    out_str += 4*' ' + 'unit angs\n'
    out_str += '}\n\n'

    # Monomer B section
    out_str += 'monomer {\n' + 4*' '
    out_str += monoB.xyz(add_header = False).strip().replace('\n', '\n    ') + '\n'
    if monoB.charge != 0:
        out_str += 4*' ' + 'charge {}\n'.format(monoB.charge)
    out_str += 4*' ' + 'unit angs\n'
    out_str += '}\n\n'

    # Midbond section
    if midbond != None:
        out_str += 'midbond {\n' + 4*' '
        out_str += midbond.xyz(add_header = False).strip().replace('\n', '\n    ') + '\n'
        out_str += 4*' ' + 'unit angs\n'
        out_str += '}\n\n'

    if fastdf_options != None:
        out_str += 'options {\n'
        for key in fastdf_options:
            out_str += 4*' ' + '{} {}\n'.format(key, fastdf_options[key])
        
        out_str += '}\n'

    return(out_str)

# ------------------------------------------------------------------------------
# Here we define functions and classes for accessing output energies

# This dictionary relates short and long names for each sapt-fastdf term
terms_dict = {
        'E^{(1)}_{elst}'         :   'e1elst',
        'E^{(1)}_{exch}'         :   'e1exch',
        'E^{(2)}_{ind,r}'        :   'e2indr',
        'E^{(2)}_{ex-ind,r}'     :   'e2exindr',
        'E^{(2)}_{disp,r}'       :   'e2dispr',
        'E^{(2)}_{exch-disp,r}'  :   'e2exdispr',
        'E^{(10)}_{elst}'        :   'e10elst',
        'E^{(10)}_{exch}'        :   'e10exch',
        'E^{(20)}_{ind,r}'       :   'e20indr',
        'E^{(20)}_{ex-ind,r}'    :   'e20exindr',
        'E^{(20)}_{disp,r}'      :   'e20dispr',
        'E^{(20)}_{exch-disp,r}' :   'e20exdispr',
        '\delta(HF,r)'           :   'deltahfr',
        '\delta(HF)'             :   'deltahf',
        'E^{HF}_int'             :   'hfint'
        }

terms_dict_rev = dict((reversed(item) for item in terms_dict.items()))

class Fastdf_results:
    def __init__(self):
        self.basis_set = None
        self.functional = None
        self.ips = [None, None]
        self.dipole = [None, None]
        self.values = {}

    def read_data(self, output_file):
        read_values = False
        read_input = False

        for line in open(output_file).readlines():

            if line.strip() == 'RESULTS':
                read_values = True

            if 'sapt-fastdf input file' in line:
                read_input = True
                input_count = 0

            if read_input:
                linel = line.strip().lower()
                if 'grac' in linel:
                    if not 'false' in linel:
                        data = linel.split()
                        self.ips[0] = float(data[3])
                        self.ips[1] = float(data[4])

                if '='*80 in linel:
                    input_count += 1
                    if input_count == 2:
                        read_input = False

                if 'functional' in linel:
                    data = linel.split()
                    self.functional = data[3]

                if 'main_basis_set' in linel:
                    data = linel.split()
                    self.basis_set = data[3]

                if 'off_basis_set' in linel:
                    data = linel.split()
                    self.basis_set = data[3]

            if read_values:
                for key in terms_dict:
                    if key in line.strip():
                        data = line.strip().split()
                        key_short = terms_dict[key]
                        self.values[key_short] = float(data[1])

    def get_value(self, key, unit='mhartree'):
        if unit == 'kcalmol':
            scalef = co.hartree_to_kcalmol/1000
        elif unit == 'mhartree':
            scalef = 1
        else:
            raise Exception('unit {} not recognized'.format(unit))
        
        found = False
        for my_key in terms_dict_rev:
            if key == my_key:
                found = True

        if not found:
            raise Exception('Error! key {} not available'.format(key))

        return self.values[key]*scalef

    def get_dipole_moment(self, filename, unit = 'debye'):
        if unit == 'debye':
            string = 'Magnitude (Debye)'
        elif unit == 'bohr':
            string = 'Magnitude (a.u.)'
        else:
            raise Exception('Unit {} not recognized'.format(unit))

        entries = []

        for line in open(filename, 'r'):
            if string in line.strip():
                entries.append(line.strip())

        for i, entry in enumerate(entries[-2:]):
            data = entry.split()
            self.dipole[i] = float(data[-1])

    def print_values(self):
        for i, j in self.values.items():
            print('{:30}{:15.8f}'.format(terms_dict_rev[i], j))

    def E_saptdft(self):
        req_terms = ['e1elst', 'e1exch', 'e2indr', 'e2exindr', 'e2dispr', 'e2exdispr']

        if all(k in self.values for k in req_terms):
            req_list = [self.values[x] for x in req_terms]
            return sum(req_list)
        else:
            return None


