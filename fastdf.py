import os
import numpy as np
from molecule import Molecule


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

if __name__ == '__main__':
    options = {
            'method' : 'HFDc1', 
            'functional' : 'PBE',
            'memory' : '16gb', 

            'ips' : 12.7871,

            'basis_set' : 'aug-cc-pvtz',

            'mb_r6' : 'M1'
    }

    fastdf_options = {
            'orca_grid' : 5,
            'orca_finalgrid' : 6,
            'nprocs' : 8,
            'orca_ex_ri_method' : 'RIJK',
            'add_deltahf' : 'true'
        }

    print(fastdf_run(options, fastdf_options))
