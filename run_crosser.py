import pr_func


genome_folder = r'C:\Users\Professional\PycharmProjects\Protein-crossing\test_data_for_sonya_short_summaries'
table_sp = r'C:\Users\Professional\PycharmProjects\Protein-crossing\table_new_species.csv'
table_gene_comb = pr_func.gene_comb(genome_folder, table_sp)
missing_list = pr_func.list_of_missing_genes(table_gene_comb)


pr_func.plot_clustered_heatmap(table_gene_comb)
print(pr_func.PCA_graph(pr_func.gene_comb(genome_folder, table_sp)))