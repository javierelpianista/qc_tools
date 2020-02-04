import sys
sys.path.append('/home/jgarcia/qc_tools')

import qc_tools.molecule as mol
import qc_tools.sbatch as sb
import qc_tools.fastdf as fdf

values = fdf.Fastdf_results()
values.read_energies('output')
values.print_values()
print(values.E_saptdft())
