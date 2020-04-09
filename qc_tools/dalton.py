import os
import qc_tools.molecule as mol
import qc_tools.basis as basis
import qc_tools.constants as co
import qc_tools.memory as mem

basis_set_dir = '/data/jgarcia/SAPT/SAPT2020.1-devel/basis_sets'

def mol_file(**kwargs):
    string = 'INTGRL\n\n\n'
    nat = 0

    if 'monoA' in kwargs:
        nat += len(kwargs['monoA'].atoms)

    if 'monoB' in kwargs:
        nat += len(kwargs['monoB'].atoms)

    if 'midbond' in kwargs:
        nat += len(kwargs['midbond'].atoms)

    string += '{:5}{:3}{:2}\n'.format(nat, 0, 0)

    count = {}

    for atom, coords in zip(kwargs['monoA'].atoms, kwargs['monoA'].coords):
        if kwargs['kind'] in ('A', 'MA', 'B'):
            basis_set = 'basis_set'
        elif kwargs['kind'] == 'MB':
            basis_set = 'off_basis_set'

        max_l  = kwargs[basis_set].atom_max_l(atom)

        if 'A' in kwargs['kind']:
            charge = co.atomic_symbols.index(atom)+1
        else:
            charge = 0

        string += 5*' ' + (3*'{:5d}').format(charge,1,max_l+1)
        string += ((max_l+1)*'{:5d}').format(*([1]*(max_l+1)))
        string += '\n'

        if atom in count:
            count[atom] += 1
        else:
            count[atom] = 1

        string += ('{:4}' + 3*'{:20.13f}').format(atom + str(count[atom]), *coords)
        string += '\n'
        string += kwargs[basis_set].atom_shell_dalton(atom)

    # Now midbonds (if any)
    for atom, coords in zip(kwargs['midbond'].atoms, kwargs['midbond'].coords):
        basis_set = 'mb_basis_set'

        max_l  = kwargs[basis_set].atom_max_l(atom)

        charge = 0

        string += 5*' ' + (3*'{:5d}').format(charge,1,max_l+1)
        string += ( (max_l+1)*'{:5d}').format(*([1]*(max_l+1)))
        string += '\n'

        if atom in count:
            count[atom] += 1
        else:
            count[atom] = 1

        string += ('{:4}' + 3*'{:20.13f}').format(atom + str(count[atom]), *coords)
        string += '\n'
        string += kwargs[basis_set].atom_shell_dalton(atom)

    # Now atoms from monomer B
    for atom, coords in zip(kwargs['monoB'].atoms, kwargs['monoB'].coords):
        if kwargs['kind'] in ('B', 'MB', 'A'):
            basis_set = 'basis_set'
        elif kwargs['kind'] == 'MA':
            basis_set = 'off_basis_set'

        max_l  = kwargs[basis_set].atom_max_l(atom)

        if 'B' in kwargs['kind']:
            charge = co.atomic_symbols.index(atom)+1
        else:
            charge = 0

        string += 5*' ' + (3*'{:5d}').format(charge,1,max_l+1)
        string += ( (max_l+1)*'{:5d}').format(*([1]*(max_l+1)))
        string += '\n'

        if atom in count:
            count[atom] += 1
        else:
            count[atom] = 1

        string += ('{:4}' + 3*'{:20.13f}').format(atom + str(count[atom]), *coords)
        string += '\n'
        string += kwargs[basis_set].atom_shell_dalton(atom)

    return(string)

def aux_file(**kwargs):
    string = ''
    for key in ['monoA', 'midbond', 'monoB']:
        if key in ['monoA', 'monoB']:
            basis_set = 'aux_basis_set'
        elif key == 'midbond':
            basis_set = 'aux_mb_basis_set'

        for atom, coords in zip(kwargs[key].atoms, kwargs[key].coords):
            try:
                num = co.atomic_symbols.index(atom)+1
            except:
                num = 0
            string += '0 z\n'
            string += ('{:2}' + 3*'{:20.12f}').format(num, *coords) + '\n'

            for shell in kwargs[basis_set].shells[atom]:
                string += '{:6d} {}\n'.format(len(shell.exponents), co.angular_momenta[shell.l])
                for expo, coef in zip(shell.exponents, shell.coefficients):
                    string += (2*'{:20.12f}').format(expo,coef) + '\n'

    return(string)

def dal_file(**kwargs):
    string = '**DALTON INPUT\n'
    string += '.RUN WAVE FUNCTION\n' 
    string += '**INTEGRALS\n'
    string += '.NOSUP\n.PRINT\n    1\n'

    if kwargs['kind'] in 'AB':
        string += '.NOTWO\n'

    string += '**WAVE FUNCTIONS\n.DFT\n{}\n'.format(kwargs['functional'])
    string += '.INTERFACE\n'

    if kwargs['kind'] in 'AB':
        string += '.STOP\nAFTER MO-ORTHNORMALIZATION\n'

    string += '*AUXILLIARY INPUT\n'
    string += '.NOSUPMAT\n'
    string += '*ORBITALS\n'
    string += '.NOSUPSYM\n'

    if kwargs['kind'] in ['MA', 'MB']:
        string += '*DFT INPUT\n'
        string += '.CKSAUX\n'
        string += '.DFTELS\n0.01\n'
        string += '.DFTAC\n'
        if kwargs['kind'] == 'MA':
            ip = kwargs['ipA']/co.mhartree_to_ev
        elif kwargs['kind'] == 'MB':
            ip = kwargs['ipB']/co.mhartree_to_ev
        string += '{:.4f} {} {} {}\n'.format(ip, 3, 4, 0)
        string += '.RADINT\n'
        string += '1.0E-13\n'
        string += '.ANGINT\n{}\n'.format(kwargs['cks_nang'])
        string += '*HF INPUT\n.THRESH\n5.0D-8\n'

    string += '*ORBITAL INPUT\n'

    if kwargs['kind'] in 'AB':
        string += '.MOSTART\n'
        string += 'H1DIAG\n'
    elif kwargs['kind'] in ['MA', 'MB']:
        string += '.AO DELETE\n'
        string += '1.0D-7\n'
        string += '.CMOMAX\n'
        string += '5.0D+3\n'

    string += '*END OF INPUT'

    return(string)

def basis_string(**kwargs):
    key_list = ['monoA', 'monoB']

    if 'midbond' in kwargs and kwargs['midbond']:
        key_list.insert(1, 'midbond')

    bstring = ''
    tstring = ''

    for key in key_list:
        if key == 'midbond':
            which = 'M'
            bset     = kwargs['mb_basis_set']
            off_bset = kwargs['mb_basis_set']
        if key == 'monoA':
            which = 'A'
            bset     = kwargs['basis_set']
            off_bset = kwargs['off_basis_set']
        elif key == 'monoB':
            which = 'B'
            bset     = kwargs['basis_set']
            off_bset = kwargs['off_basis_set']

        for atom, coords in zip(kwargs[key].atoms, kwargs[key].coords):
            for shell in bset.shells[atom]:
                bstring += co.angular_momenta[shell.l]
                if shell in off_bset.shells[atom]:
                    tstring += 'M'
                else:
                    tstring += which

    bstring += '+'
    tstring += '+'

    return bstring, tstring

def P_file(**kwargs):
    bstring, tstring = basis_string(**kwargs)

    memory = int(mem.str_to_int(kwargs['memory'])/8)

    string = 'File generated by qc_tools\n\n'

    string += '&TRN\n'
    string += 'ISITALCH=F, ISITG88=F, ISITG90=F, ISITHNDO=F, ISITMICR=F,\n'
    string += 'ISITATM=F, ISITCADP=F, ISITGAMS=F, ISITDALT=T, T2EL=F,\n'
    string += 'OUT=F, TOLER=15, DIMER=F, BLKMB=F, SPHG=T, MEMTRAN=\"{}\",\n'.format(memory)

    string += 'BASIS=\'{}\'\n'.format(bstring)
    string += ' TAGS=\'{}\'\n'.format(tstring)
    string += '&END\n\n'

    string += '&CCINP\n'
    string += ' CCPRINT=F, VCRIT=1.0D-10, TOLITER=1.0D-5, MEMCC="40000000"\n'
    string += '&END\n\n'

    string += '&INPUTCOR\n'
    string += ' SAPTKS=T, PRINT=F, MEMSAPT="{}"\n'.format(memory)
    string += '&END\n\n'

    string += '&DF\n'
    string += '  MEMTRAN="{}"\n'.format(memory)
    string += '&END\n'

    return(string)

def dalton_run(**kwargs):
    # Set jobname
    if 'jobname' in kwargs:
        jobname = kwargs['jobname']
    else:
        jobname = 'in'

    # Define basis sets and read them from files
    bset = basis.BasisSet.from_gbs_file(os.path.join(basis_set_dir, kwargs['basis_set'] + '.gbs'))

    if 'off_basis_set' in kwargs:
        off_bset = basis.BasisSet.from_gbs_file(os.path.join(basis_set_dir, kwargs['off_basis_set'] + '.gbs'))
    elif 'reduction' in kwargs:
        off_bset = bset.reduced(kwargs['reduction'])
    else:
        off_bset = bset

    if 'aux_basis_set' in kwargs:
        aux_bset = basis.BasisSet.from_gbs_file(os.path.join(basis_set_dir, kwargs['aux_basis_set'] + '.gbs'))
    else:
        aux_bset_name = kwargs['basis_set'] + '-ri'
        aux_bset_fname = os.path.join(basis_set_dir, aux_bset_name + '.gbs')

        if os.path.exists(os.path.join(basis_set_dir, aux_bset_fname)):
            aux_bset = basis.BasisSet.from_gbs_file(aux_bset_fname)
        else:
            raise Exception('Auxiliary basis set not specified, and the default one {} wasn\'t found in the repository.'.format(aux_bset_name))

    if 'mb_basis_set' in kwargs:
        mb_bset = basis.BasisSet.from_gbs_file(os.path.join(basis_set_dir, kwargs['mb_basis_set'] + '.gbs'))
    else:
        mb_bset = basis.BasisSet.from_gbs_file(os.path.join(basis_set_dir, 'midbond.gbs'))

    if 'aux_mb_basis_set' in kwargs:
        aux_mb_bset = basis.BasisSet.from_gbs_file(os.path.join(basis_set_dir, kwargs['aux_mb_basis_set'] + '.gbs'))
    else:
        aux_mb_bset = basis.BasisSet.from_gbs_file(os.path.join(basis_set_dir, 'aux_midbond.gbs'))

    # Now write the corresponding input files
    new_opts = kwargs

    new_opts['basis_set']        = bset
    new_opts['aux_basis_set']    = aux_bset
    new_opts['off_basis_set']    = off_bset
    new_opts['mb_basis_set']     = mb_bset
    new_opts['aux_mb_basis_set'] = aux_mb_bset

    for kind in ['MA', 'MB', 'A', 'B']:
        new_opts['kind'] = kind

        string = mol_file(**new_opts)
        open(jobname + kind + '.mol', 'w').write(string)
    
        string = dal_file(**new_opts)
        open(jobname + kind + '.dal', 'w').write(string)

    open(jobname + '.aux', 'w').write(aux_file(**new_opts))
    open(jobname + 'P.data', 'w').write(P_file(**new_opts))

