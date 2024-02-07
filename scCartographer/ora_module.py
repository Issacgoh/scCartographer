
#Functional modules
import pandas as pd
import numpy as np
import mygene
import gseapy as gp
import matplotlib.pyplot as plt
import mygene
import scipy.sparse as sparse
from sklearn.metrics.pairwise import cosine_similarity

def importlibs():
    import pandas as pd
    import numpy as np
    import mygene
    import gseapy as gp
    import matplotlib.pyplot as plt
    import mygene
    import scipy.sparse as sparse
    from sklearn.metrics.pairwise import cosine_similarity

# ENSDB-HGNC Option 1 
#from gseapy.parser import Biomart
#bm = Biomart()
## view validated marts#
#marts = bm.get_marts()
## view validated dataset
#datasets = bm.get_datasets(mart='ENSEMBL_MART_ENSEMBL')
## view validated attributes
#attrs = bm.get_attributes(dataset='hsapiens_gene_ensembl')
## view validated filters
#filters = bm.get_filters(dataset='hsapiens_gene_ensembl')
## query results
#queries = ['ENSG00000125285','ENSG00000182968'] # need to be a python list
#results = bm.query(dataset='hsapiens_gene_ensembl',
#                       attributes=['ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'go_id'],
#                       filters={'ensemble_gene_id': queries}
                      
# ENSDB-HGNC Option 2
def convert_hgnc(input_gene_list):
    import mygene
    mg = mygene.MyGeneInfo()
    geneList = input_gene_list 
    geneSyms = mg.querymany(geneList , scopes='ensembl.gene', fields='symbol', species='human')
    return(pd.DataFrame(geneSyms))
# Example use: convert_hgnc(['ENSG00000148795', 'ENSG00000165359', 'ENSG00000150676'])

# Scanpy_degs_to_long_format
def convert_scanpy_degs(input_dataframe):
    if 'concat' in locals() or 'concat' in globals():
        del(concat)
    degs = input_dataframe
    n = degs.loc[:, degs.columns.str.endswith("_n")]
    n = pd.melt(n)
    p = degs.loc[:, degs.columns.str.endswith("_p")]
    p = pd.melt(p)
    l = degs.loc[:, degs.columns.str.endswith("_l")]
    l = pd.melt(l)
    n = n.replace(regex=r'_n', value='')
    n = n.rename(columns={"variable": "cluster", "value": "gene"})
    p = (p.drop(["variable"],axis = 1)).rename(columns={ "value": "p_val"})
    l = (l.drop(["variable"],axis = 1)).rename(columns={ "value": "logfc"})
    return(pd.concat([n,p,l],axis=1))
#Usage: convert_scanpy_degs(scanpy_degs_file)

# Clean convert gene list to list
def as_gene_list(input_df,gene_col):
    gene_list = input_df[gene_col]
    glist = gene_list.squeeze().str.strip().tolist()
    return(glist)

# No ranking enrichr
def enrich_no_rank(input_gene_list,database,species="Human",description="enr_no_rank",outdir = "./enr",cutoff=0.5):
    # list, dataframe, series inputs are supported
    enr = gp.enrichr(gene_list=input_gene_list,
                     gene_sets=database,
                     organism=species, # don't forget to set organism to the one you desired! e.g. Yeast
                     #description=description,
                     outdir=outdir,
                     # no_plot=True,
                     cutoff=cutoff # test dataset, use lower value from range(0,1)
                    )
    return(enr)
    #Usge: enrich_no_rank(gene_as_list)
    
# Custom genelist test #input long format degs or dictionary of DEGS
def custom_local_GO_enrichment(input_gene_list,input_gmt,input_gmt_key_col,input_gmt_values,description="local_go",Background='hsapiens_gene_ensembl',Cutoff=0.5):
    
    #Check if GMT input is a dictionary or long-format input
    if isinstance(input_gmt, dict):
        print("input gmt is a dictionary, proceeding")
        dic = input_gmt
    else:
        print("input gmt is not a dictionary, if is pandas df,please ensure it is long-format proceeding to convert to dictionary")
        dic =  input_gmt.groupby([input_gmt_key_col])[input_gmt_values].agg(lambda grp: list(grp)).to_dict()
        
    enr_local = gp.enrichr(gene_list=input_gene_list,
                 description=description,
                 gene_sets=dic,
                 background=Background, # or the number of genes, e.g 20000
                 cutoff=Cutoff, # only used for testing.
                 verbose=True)
    return(enr_local)
    #Example_usage: custom_local_GO_enrichment(input_gene_list,input_gmt,input_gmt_key_col,input_gmt_values) #input gmt can be long-format genes and ontology name or can be dictionary of the same   

    
# Pre-ranked GS enrichment
def pre_ranked_enr(input_gene_list,gene_and_ranking_columns,database ='GO_Biological_Process_2018',permutation_num = 1000, outdir = "./enr_ranked",cutoff=0.25,min_s=5,max_s=1000):
    glist = input_gene_list[gene_and_ranking_columns]
    pre_res = gp.prerank(glist, gene_sets=database,
                     threads=4,
                     permutation_num=permutation_num, # reduce number to speed up testing
                     outdir=outdir,
                     seed=6,
                     min_size=min_s,
                     max_size=max_s)
    return(pre_res)
    #Example usage: pre_ranked_hyper_geom(DE, gene_and_ranking_columns = ["gene","logfc"],database=['KEGG_2016','GO_Biological_Process_2018'])

    
# GSEA module for permutation test of differentially regulated genes
# gene set enrichment analysis (GSEA), which scores ranked genes list (usually based on fold changes) and computes permutation test to check if a particular gene set is more present in the Up-regulated genes, 
# among the DOWN_regulated genes or not differentially regulated.
#NES = normalised enrichment scores accounting for geneset size
def permutation_ranked_enr(input_DE,cluster_1,cluster_2,input_DE_clust_col,input_ranking_col ,input_gene_col ,database = "GO_Biological_Process_2018"):
    input_DE = input_DE[input_DE[input_DE_clust_col].isin([cluster_1,cluster_2])]
    #Make set2 negative values to represent opposing condition
    input_DE[input_ranking_col].loc[input_DE[input_DE_clust_col].isin([cluster_2])] = -input_DE[input_ranking_col].loc[input_DE[input_DE_clust_col].isin([cluster_2])]
    enr_perm = pre_ranked_enr(input_DE,[input_gene_col,input_ranking_col],database,permutation_num = 100, outdir = "./enr_ranked_perm",cutoff=0.5)
    return(enr_perm)
    #Example usage:permutation_ranked_enr(input_DE = DE, cluster_1 = "BM",cluster_2 = "YS",input_DE_clust_col = "cluster",input_ranking_col = "logfc",input_gene_col = "gene",database = "GO_Biological_Process_2018")
    #input long-format list of genes and with a class for permutaion, the logfc ranking should have been derived at the same time


#Creating similarity matrix from nested gene lists
def create_sim_matrix_from_enr(input_df,nested_gene_column="Genes",seperator=";",term_col="Term"):
#    input_df[gene_col] = input_df[gene_col].astype(str)
#    input_df[gene_col] = input_df[gene_col].str.split(";")
#    uni_val = list(input_df.index.unique())
#    sim_mat = pd.DataFrame(index=uni_val, columns=uni_val)
#    exploded_df = input_df.explode(gene_col)
#    # Ugly loop for cosine gs similarity matrix (0-1)
#    for i in (uni_val):
#        row = exploded_df[exploded_df.index.isin([i])]
#        for z in (uni_val):
#            col = exploded_df[exploded_df.index.isin([z])]
#            col_ls = col[gene_col]
#            row_ls = row[gene_col]
#            sim = len(set(col_ls) & set(row_ls)) / float(len(set(col_ls) | set(row_ls)))
#            sim_mat.loc[i , z] = sim

#    Check term col in columns else, check index as it\s sometimes heree
    if not term_col in list(input_df.columns):
        input_df[term_col] = input_df.index

#    Create a similarity matrix by cosine similarity
    input_df = input_df.copy()
    gene_col = nested_gene_column #"ledge_genes"
    input_df[gene_col] = input_df[gene_col].astype(str)
    input_df[gene_col] = input_df[gene_col].str.split(seperator)
    term_vals = list(input_df[term_col].unique())
    uni_val = list(input_df[term_col].unique())
    sim_mat = pd.DataFrame(index=uni_val, columns=uni_val)
    exploded_df = input_df.explode(gene_col)
    arr = np.array(input_df[gene_col])
    vals = list(exploded_df[nested_gene_column])
    import scipy.sparse as sparse
    def binarise(sets, full_set):
        """Return sparse binary matrix of given sets."""
        return sparse.csr_matrix([[x in s for x in full_set] for s in sets])
    # Turn the matrix into a sparse boleen matrix of binarised values
    sparse_matrix = binarise(arr, vals)
    from sklearn.metrics.pairwise import cosine_similarity
    similarities = cosine_similarity(sparse_matrix)
    sim_mat = pd.DataFrame(similarities)
    sim_mat.index = uni_val
    sim_mat.columns = uni_val
    return(sim_mat)
#Example usage : sim_mat = create_sim_matrix_from_enr(enr.res2d)


#Creating similarity matrix from GO terms
def create_sim_matrix_from_term(input_df,nested_gene_column="Term",seperator=" ",term_col="Term"):

#    Check term col in columns else, check index as it\s sometimes heree
    if not term_col in list(input_df.columns):
        input_df[term_col] = input_df.index

#    Create a similarity matrix by cosine similarity
    input_df = input_df.copy()
    gene_col = nested_gene_column #"ledge_genes"
    #input_df[gene_col] = input_df[gene_col].astype(str)
    input_df[gene_col] = input_df[gene_col].str.split(seperator)
    term_vals = list(input_df[term_col].unique())
    uni_val = list(input_df.index.unique())
    sim_mat = pd.DataFrame(index=uni_val, columns=uni_val)
    exploded_df = input_df.explode(gene_col)
    arr = np.array(input_df[gene_col])
    vals = list(exploded_df[nested_gene_column])
    import scipy.sparse as sparse
    def binarise(sets, full_set):
        """Return sparse binary matrix of given sets."""
        return sparse.csr_matrix([[x in s for x in full_set] for s in sets])
    sparse_matrix = binarise(arr, vals)
    from sklearn.metrics.pairwise import cosine_similarity
    similarities = cosine_similarity(sparse_matrix)
    sim_mat = pd.DataFrame(similarities)
    sim_mat.index = uni_val
    sim_mat.columns = uni_val
    return(sim_mat)

#Creating similarity matrix from GO terms
def create_sim_matrix_from_term2(input_df,nested_gene_column="Term",seperator=" ",term_col="Term"):
#    Horrifically bad cosine similairty estimate for word frequency
#    Check term col in columns else, check index as it\s sometimes heree
    if not term_col in list(input_df.columns):
        input_df[term_col] = input_df.index
    input_df = input_df.copy()
    gene_col = nested_gene_column #"ledge_genes"
    #input_df[gene_col] = input_df[gene_col].astype(str)
    term_vals = list(input_df[term_col].unique())
    input_df[gene_col] = input_df[gene_col].str.split(seperator)
    uni_val = list(input_df.index.unique())
    sim_mat = pd.DataFrame(index=uni_val, columns=uni_val)
    exploded_df = input_df.explode(gene_col)

    nan_value = float("NaN")
    exploded_df.replace("", nan_value, inplace=True)
    exploded_df.dropna(subset = [gene_col], inplace=True)
    arr = np.array(input_df[gene_col])

    vals = list(exploded_df[nested_gene_column])

    import scipy.sparse as sparse
    def binarise(sets, full_set):
        """Return sparse binary matrix of given sets."""
        return sparse.csr_matrix([[x in s for x in full_set] for s in sets])
    sparse_matrix = binarise(arr, vals)
    from sklearn.metrics.pairwise import cosine_similarity
    similarities = cosine_similarity(sparse_matrix)
    sim_mat = pd.DataFrame(similarities)
    sim_mat.index = uni_val
    sim_mat.columns = uni_val
    return(sim_mat)
    #Example usage : sim_mat = create_sim_matrix_from_enr(enr.res2d)
    
    
def compute_weighted_impact(varm_file, top_loadings, threshold=0.05):
    """
    Computes the weighted impact of the features of a low-dimensional model.

    Parameters:
    varm_file (str): The path to the file containing the variable loadings of the model.
    top_loadings (pd.DataFrame): A dataframe containing the top loadings for each class.
    threshold (float): The p-value threshold for significance.

    Returns:
    top_loadings_lowdim (pd.DataFrame): A dataframe containing the top weighted impacts for each class.
    """
    # Load the variable loadings from the file
    model_varm = pd.read_csv(varm_file, index_col=0)

    # Map the feature names to the column names of the variable loadings
    feature_set = dict(zip(sorted(top_loadings['feature'].unique()), model_varm.columns))

    # Melt the variable loadings dataframe and add a column for p-values
    varm_melt = pd.melt(model_varm.reset_index(), id_vars='index')
    varm_melt['pvals'] = np.nan

    # Compute the p-values for each variable
    for variable in varm_melt['variable'].unique():
        varm_loadings = varm_melt[varm_melt['variable'] == variable]
        med = np.median(varm_loadings['value'])
        mad = np.median(np.abs(varm_loadings['value'] - med))
        pvals = scipy.stats.norm.sf(varm_loadings['value'], loc=med, scale=1.4826*mad)
        varm_melt.loc[varm_melt['variable'] == variable, 'pvals'] = pvals

    # Filter the variables based on the p-value threshold
    varm_sig = varm_melt[varm_melt['pvals'] < threshold]

    # Compute the weighted impact for each feature of each class
    top_loadings_lowdim = pd.DataFrame(columns=['class', 'feature', 'weighted_impact', 'e^coef_pval', 'e^coef', 'is_significant_sf'])
    top_loadings_lw = top_loadings.groupby('class').head(10)
    top_loadings_lw['feature'] = top_loadings_lw['feature'].map(feature_set)

    for classes in top_loadings_lw['class'].unique():
        for feature in top_loadings_lw.loc[top_loadings_lw['class'] == classes, ['feature', 'e^coef']].values:
            temp_varm_sig = varm_sig[varm_sig['variable'] == feature[0]]
            temp_varm_sig['weighted_impact'] = temp_varm_sig['value'] * feature[1]
            temp_varm_sig = temp_varm_sig[['index', 'weighted_impact']]
            temp_varm_sig.columns = ['feature', 'weighted_impact']
            temp_varm_sig['class'] = classes
            temp_varm_sig['e^coef_pval'] = top_loadings_lw.loc[(top_loadings_lw['class'] == classes) & (top_loadings_lw['feature'] == feature[0]), 'e^coef_pval'].values[0]
            temp_varm_sig['e^coef'] = top_loadings_lw.loc[(top_loadings_lw['class'] == classes) & (top_loadings_lw['feature'] == feature[0]), 'e^coef'].values[0]
            temp_varm_sig['is_significant_sf'] = top_loadings_lw.loc[(top_loadings_lw['class'] == classes) & (top_loadings_lw['feature'] == feature[0]), 'is_significant_sf'].values[0]
            top_loadings_lowdim = pd.concat([top_loadings_lowdim, temp_varm_sig], ignore_index=True)
    # Return the top 100 features with the highest weighted impact for each class
    top_loadings_lowdim = top_loadings_lowdim.sort_values('weighted_impact', ascending=False).groupby('class').head(100)
    return top_loadings_lowdim
