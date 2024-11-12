import subprocess
import sys

import pandas as pd 

sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

# snakemake variables
func_effects_data = snakemake.input.functional_data
output_pdb_file_name = snakemake.output.output_file

data = pd.read_csv(func_effects_data)
data['reference_site'] = data['site']
data.to_csv(output_pdb_file_name)