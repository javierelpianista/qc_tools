import os
import copy
import numpy as np
from qc_tools.molecule import Molecule
from qc_tools.fastdf import mb_r6
import qc_tools.memory as mem
import qc_tools.constants as co

default_jobname = 'in'

def get_header(**kwargs):
    header = '! {} {}'.format(kwargs['method'], kwargs['basis_set'])

    if 'aux_basis_set' in kwargs:
        header += ' ' + kwargs['aux_basis_set']

    if 'conv_speed' in kwargs:
        header += ' ' + kwargs['conv_speed']

    if kwargs['method'].lower() == 'dlpno-ccsd(t)':
        header += ' ' + kwargs['aux_basis_set_c']

        if 'pno' in kwargs:
            header += ' ' + kwargs['pno']
        else:
            header += ' tightpno'

    return(header)

def get_molecule_string(molecule, ghost = None):
    molecule_string = '%coords\n    Ctyp xyz\n    Charge {}\n    Mult {}\n    coords\n'.format(
            molecule.charge, molecule.multiplicity)

    for atom, coords in zip(molecule.atoms, molecule.coords):
        molecule_string += '{:8d}{:15.8f}{:15.8f}{:15.8f}\n'.format(co.atomic_symbols.index(atom)+1,
            *[x*co.bohr_to_angs for x in coords])

    if isinstance(ghost, Molecule):
        for atom, coords in zip(ghost.atoms, ghost.coords):
            molecule_string += '{:8d}:{:15.8f}{:15.8f}{:15.8f}\n'.format(co.atomic_symbols.index(atom)+1,
                *[x*co.bohr_to_angs for x in coords])

    molecule_string += '    end\n'
    molecule_string += 'end'

    return(molecule_string)

def run(**kwargs):
    if not 'jobname' in kwargs:
        kwargs['jobname'] = default_jobname

    header = get_header(**kwargs)

    if not 'ghost' in kwargs:
        kwargs['ghost'] = None

    molecule_string = get_molecule_string(kwargs['molecule'], ghost = kwargs['ghost'])

    if not 'nprocs' in kwargs: kwargs['nprocs'] = 1
    pal_string = '%pal\n    nprocs {}\nend'.format(kwargs['nprocs'])

    if 'memory' in kwargs:
        memory = mem.str_to_int(kwargs['memory'])
        maxcore_string = '%maxcore {}\n'.format(memory/kwargs['nprocs']/1024/1024)

    if any(i in kwargs for i in ['guess', 'maxiter', 'conv_tol']):
        scf_string = '%scf\n'

        if 'guess' in kwargs:
            scf_string += '    guess {}\n'.format(kwargs['guess'])
            if kwargs['guess'] == 'moread':
                scf_string += '    moinp "{}"\n'.format(kwargs['moinp'])

        if 'maxiter' in kwargs:
            scf_string += '    maxiter {}\n'.format(kwargs['maxiter'])

        if 'conv_tol' in kwargs:
            scf_string += '    convergence {}\n'.format(kwargs['conv_tol'])

        scf_string += 'end\n\n'
    else:
        scf_string = ''

    base_string = '%base "{}"'.format(kwargs['jobname'])
    ans = header + '\n\n' + base_string + '\n\n' + pal_string + '\n\n' + maxcore_string + '\n' + scf_string + molecule_string + '\n'

    return(ans)

def ip_run(**kwargs):
    if not 'jobname' in kwargs:
        kwargs['jobname'] = default_jobname

    if not 'ion_maxiter' in kwargs:
        kwargs['ion_maxiter'] = 1000

    if not 'conv_tol' in kwargs:
        kwargs['conv_tol'] = 'tight'

    run1 = run(**kwargs)

    molecule2 = kwargs['molecule']
    molecule2.charge += 1
    molecule2.multiplicity += 1

    opts = copy.copy(kwargs)
    opts['molecule']       = molecule2
    opts['jobname']        = kwargs['jobname'] + '_ion'
    opts['guess']          = 'moread'
    opts['moinp']          = kwargs['jobname'] + '.gbw'
    opts['maxiter']        = kwargs['ion_maxiter']
    opts['conv_speed']     = 'SlowConv'

    run2 = run(**opts)

    ans = run1 + '\n$new_job\n\n' + run2
    return(ans)

def get_ip_from_output(filename):
    read_file = open(filename, 'r')

    results = []

    for line in read_file:
        if 'FINAL SINGLE POINT ENERGY' in line:
            results.append(float(line.strip().split()[4]))

    if len(results) != 2:
        raise Exception('The number of times the line FINAL SINGLE POINT ENERGY appeared in {} is {}. Expected 2.'.format(filename, len(results)))

    return(co.mhartree_to_ev*(results[1]-results[0]))

