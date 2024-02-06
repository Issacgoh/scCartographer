import gseapy
names = gseapy.get_library_name()
print(names)

enr.res2d



temp_varm_sig

enr.res2d

enr.res2d

#Create seperate columns for term and GO term (Not needed for direct hypergeom results)
enr_terms = enr.res2d.copy()
# enr_terms["Term"] = enr_terms.index
enr_terms["GO"] = enr_terms["Term"].copy()
#Get substring between "(GO)" for GO terms
go_terms = enr_terms["Term"].apply(lambda st: st[st.find("(GO")+1:st.find(")")])
enr_terms["go_terms"] = go_terms
enr_terms["Term"] = enr_terms["Term"].str.replace('\(GO.*$', '')

sim_mat_gene = create_sim_matrix_from_enr(enr_terms,nested_gene_column="Lead_genes",seperator=";",term_col="Term")
#sim_mat.to_csv("./enr_sim_mat.csv")
sim_mat_gene

sim_mat_term = create_sim_matrix_from_term2(enr_terms,term_col = "Term",nested_gene_column ="Term")
print(sim_mat_term.shape)
sim_mat_term.columns = sim_mat_gene.columns
sim_mat_term.index = sim_mat_gene.index
sim_mat_term

sim_mat = (np.array(sim_mat_term)**2)* (np.array(sim_mat_gene)**2)
sim_mat = pd.DataFrame(sim_mat)
sim_mat.index = sim_mat_term.index
sim_mat.columns = sim_mat_term.columns
sim_mat

import leidenalg as la
import igraph as ig
import pandas as pd
import networkx as nx

sim_mat = pd.DataFrame(sim_mat.to_numpy(), index=sim_mat.index, columns=sim_mat.columns)
g = nx.from_pandas_adjacency(sim_mat)
#g.edges(data=True)

nx.draw(g,with_labels = True,node_color='b',node_size=500);

g

sim_mat = sim_mat_gene[:]
#Contruct graph
sim_mat_val = sim_mat.values
g = ig.Graph.Adjacency((sim_mat_val > 0).tolist())

# Add edge weights and node labels.
#g.es['weight'] = sim_mat_val[sim_mat_val.nonzero()]
#g.vs['label'] = list(sim_mat.index)  # or a.index/a.columns
#g.es['width'] = sim_mat_val[sim_mat_val.nonzero()]