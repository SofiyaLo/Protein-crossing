import pr_func
import PoFF_analysis

genome_folder = r'C:\Users\Professional\PycharmProjects\Protein-crossing\test_data_for_sonya_short_summaries'
table_sp = r'C:\Users\Professional\PycharmProjects\Protein-crossing\table_new_species.csv'
#table_gene_comb = pr_func.gene_comb(genome_folder, table_sp)
#missing_list = pr_func.list_of_missing_genes(table_gene_comb)
pr_ortho_tsv = r'C:\Users\Professional\PycharmProjects\Protein-crossing\Prot_ortho_results\po_no_sensitive.proteinortho.tsv'
pr_ortho_tsv1 = r'C:\Users\Professional\PycharmProjects\Protein-crossing\Prot_ortho_results\very_sensitive.proteinortho.tsv'

w_ortho_df = PoFF_analysis.ortho_weight_reader(pr_ortho_tsv, 4)
unw_ortho_df = PoFF_analysis.ortho_unweight_reader(pr_ortho_tsv, 4)

#print(PoFF_analysis.PCA_graph_spic(w_ortho_df))
#print(PoFF_analysis.PCA_graph_genes(w_ortho_df))
#print(PoFF_analysis.PCA_graph_spic(unw_ortho_df))
print(PoFF_analysis.PCA_graph_genes(unw_ortho_df))