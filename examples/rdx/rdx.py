import sys
sys.path.append('/home/jgarcia/qc_tools')

import qc_tools.molecule as mol
import qc_tools.sbatch as sb
import qc_tools.fastdf as fdf

fastdf_options = {
        'method' : 'SAPTDFT',
        'functional' : 'PBE', 
        'memory'     : '16gb', 

        'ips'        : 10.61244,
        'basis_set'  : 'aug-cc-pvtz',

        'mb_r6'      : 'M1'
    }

fastdf_additional_options = {
        'orca_grid' : 5,
        'orca_finalgrid' : 6,
        'nprocs' : 8,
        'orca_ex_ri_method' : 'RIJK', 
        'freeze_core' : 'false'
        }

sbatch_options = {
        'kind' : 'fastdf',
        'jobname' : 'rdx',
        'nprocs'  : 8,
        'memory'  : '18gb', 
        'time'    : '2:00:00', 
        'exclusive' : True,
        'modules'   : ['gcc/8.2.0', 'openmpi/3.1.3-intel']
        }

with open('rdx.sapt', 'w') as write_file:
    write_file.write(fdf.fastdf_run(fastdf_options, fastdf_additional_options))

with open('sub_rdx.sh', 'w') as write_file:
    write_file.write(sb.submit_script(sbatch_options))
