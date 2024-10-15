import pr_func

genome_folder = r'C:\Users\Professional\PycharmProjects\Protein-crossing\genome_folder'
table_sp = r'C:\Users\Professional\PycharmProjects\Protein-crossing\table_with_species.csv'
table_gene_comb = pr_func.gene_comb(genome_folder, table_sp)
missing_list = pr_func.list_of_missing_genes(table_gene_comb)
print(missing_list)