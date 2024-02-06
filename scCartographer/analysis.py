#Get possible enrichr queries
names = gp.get_library_name('Human') # default: Human "e.g GO_Biological_Process_2018"
#Define the databases to query
database='GO_Biological_Process_2018'

enr = pre_ranked_enr(Example_DE,["gene","weighted_impact"],database = database)

print("Normalised enrichment score ranges, for continous phenotype tests, a positive value indicates correlation with phenotype profile or top of the ranked list, negative values indicate inverse corelation with profile or corelation with the bottom of the list")
# print(max(enr.res2d.nes))
# print(min(enr.res2d.nes))
print(enr.res2d.shape)
enr.res2d.head(3)
#For ranked sets of individual comparisons, it may be wise to determine an ES cut-off

#enr.res2d.sort_values(by=['nes'], ascending=False)

enr = permutation_ranked_enr(input_DE = DE, cluster_1 = "BM",cluster_2 = "YS",input_DE_clust_col = "cluster",input_ranking_col = "logfc",input_gene_col = "gene",database = "GO_Biological_Process_2018")
print("Normalised enrichment score ranges, for descrete phenotype tests, positive values indicate correlation with the first phenotype and a negative value indicates correlation with the second phenotype.")
# print(max(enr.res2d.nes))
# print(min(enr.res2d.nes))
enr.res2d.head(3)