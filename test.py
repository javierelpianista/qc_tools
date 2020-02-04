import qc_tools.molecule as mol
import qc_tools.fastdf as fdf
import qc_tools.sbatch as sb

ammonia = mol.Molecule.from_xyz_file('A.xyz')
print(ammonia.xyz())

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

print(fdf.fastdf_run(options, fastdf_options))

opts = {
        'kind'      : 'fastdf',
        'jobname'   : 'test1',
        'filename'  : 'sub_test.sh',
        'nprocs'    : 8,
        'time'      : '1:00:00',
        'memory'    : '8gb',
        'exclusive' : True,
        'modules'   : 'gcc/8.2.0', 
        }
print(sb.submit_script(opts))
