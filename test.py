import qc_tools.molecule as mol
import qc_tools.orca as orca

opts = {
        'method'        : 'PBE0',
        'basis_set'     : 'aug-cc-pVTZ',
        'aux_basis_set' : 'aug-cc-pVTZ/JK',
        'jobname'       : '124',
        'memory'        : '100gb', 
        }

water = mol.Molecule.from_xyz_file('A.xyz', charge = 0)

print(orca.ip_run(molecule = water, **opts))
