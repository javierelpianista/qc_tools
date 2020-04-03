import os

fastdf_executable = '/data/sw/SAPT2020.1/bin/sapt-fastdf'
orca_executable   = '/data/jgarcia/orca/orca_4_2_1_linux_x86-64_openmpi314/orca'

available_kinds = ['fastdf', 'orca']

default_scratch_dir = '/scratch/local/jgarcia'

def submit_script(opts):
    if 'kind' in opts:
        kind = opts['kind']
    else:
        error_str = 'ERROR!!! please state which kind of run you are preparing with the kind option. Available kinds are:' 
        for kind in available_kinds: 
            error_str += ' ' + kind
        raise Exception(error_str)
    if 'jobname' in opts:
        jobname = opts['jobname'].strip()
    else:
        raise Exception('ERROR!!! jobname is a mandatory option')

    if 'input_args' in opts:
        input_args = opts['input_args']
        if not isinstance(input_args, list):
            input_args = [input_args]
    else:
        input_args = [jobname]

    if 'localdir' in opts:
        localdir = opts['localdir']
    else:
        localdir = os.getcwd()

    if 'scratchdir' in opts:
        if opts['scratchdir'].lower() == '_localdir':
            scratchdir = localdir
        else:
            scratchdir = opts['scratchdir']
    else:
        scratchdir = default_scratch_dir + '/$SLURM_JOB_ID'

    if 'executable' in opts:
        executable = opts['executable']
    else:
        if kind == 'fastdf':
            executable = fastdf_executable
        elif kind == 'orca':
            executable = orca_executable

    if 'remove_scratch' in opts:
        remove_scratch = opts['remove_scratch']
    else:
        remove_scratch = True

# ------------------------------------------------------------------------------
    out_str = '#!/bin/bash\n\n'

    out_str += '#SBATCH --job-name={}\n'.format(opts['jobname'])

    if 'memory' in opts:
        out_str += '#SBATCH --mem={}\n'.format(opts['memory'])

    if 'nprocs' in opts:
        out_str += '#SBATCH --tasks-per-node={}\n'.format(opts['nprocs'])

    if 'time' in opts:
        out_str += '#SBATCH --time={}\n'.format(opts['time'])

    if 'output_filename' in opts:
        out_str += '#SBATCH --output={}\n'.format(opts['output_filename'])
    else:
        out_str += '#SBATCH --output=%x-%j.out\n'

    if 'exclude' in opts:
        out_str += '#SBATCH --exclude={}\n'.format(opts['exclude'])

    if 'exclusive' in opts and opts['exclusive']:
        out_str += '#SBATCH --exclusive\n'

    out_str += '\n'
    out_str += 'echo "Compute node: $HOSTNAME"\n'
    out_str += 'echo "Job ID: $SLURM_JOBID"\n'
    out_str += 'echo "Started at:"\n'
    out_str += 'date\n'
    out_str += '\n'

    if 'modules' in opts:
        out_str +='module purge\n'

        modules = opts['modules']
        out_str += 'module load'
        if isinstance(modules, list):
            for i in modules:
                out_str += ' ' + i
        else:
            out_str += ' ' + modules
        out_str += '\n'
    out_str += '\n'

    out_str += 'JOBNAME={}\n'.format(opts['jobname'])
    out_str += '\n'

    out_str += 'LOCALDIR={}\n'.format(localdir)
    out_str += 'SCRATCHDIR={}\n'.format(scratchdir)
    out_str += 'mkdir -p $SCRATCHDIR\n'
    out_str += 'cp -r $LOCALDIR/* $SCRATCHDIR/\n'
    out_str += 'cd $SCRATCHDIR\n'
    out_str += '\n'

    args_str = ''
    for i, arg in enumerate(input_args):
        if i != 0:
            args_str += ' '
        args_str += arg
    out_str += '{} {}\n'.format(executable, args_str)

    if kind == 'fastdf':
        out_str += 'cat *.timer\n'

    if remove_scratch:
        out_str += 'rm -r $SCRATCHDIR\n'

    out_str += '\n'

    out_str += 'echo "Done at:"\n'
    out_str += 'hostname\n'
    out_str += 'date\n'

    return out_str

if __name__ == '__main__':
    opts = {
            'kind'      : 'orca',
            'jobname'   : 'test1.inp',
            'filename'  : 'sub_test.sh',
            'nprocs'    : 8,
            'time'      : '1:00:00',
            'memory'    : '8gb',
            'exclusive' : True,
            'modules'   : 'gcc/8.2.0', 
            }

    print(submit_script(opts))
