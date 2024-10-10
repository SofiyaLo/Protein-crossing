import pandas as pd
import os
import pr_func

genome_folder = '/workspaces/Protein-crossing/genome_folder'
table_sp = '/workspaces/Protein-crossing/table_with_species.csv'
print(pr_func.gene_comb(genome_folder, table_sp))