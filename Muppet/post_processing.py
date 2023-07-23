from collections import Counter
from collections import defaultdict
import scanpy as sc
import pandas as pd
import pickle as pkl
import numpy as np
import scipy
import matplotlib.pyplot as plt
import re
import glob
import os
import sys
#from geosketch import gs
from numpy import cov
import scipy.cluster.hierarchy as spc
import seaborn as sns; sns.set(color_codes=True)
from sklearn.linear_model import LogisticRegression
import sklearn
from pathlib import Path
import requests
import psutil
import random
import threading
import tracemalloc
import itertools
import math
import warnings
import sklearn.metrics as metrics

def load_models(model_dict,model_run):
    if (Path(model_dict[model_run])).is_file():
        # Load data (deserialize)
        model = pkl.load(open(model_dict[model_run], "rb"))
        return model
    elif 'http' in model_dict[model_run]:
        print('Loading model from web source')
        r_get = requests.get(model_dict[model_run])
        fpath = './model_temp.sav'
        open(fpath , 'wb').write(r_get.content)
        model = pkl.load(open(fpath, "rb"))
        return model

def load_adatas(adatas_dict,data_merge, data_key_use,QC_normalise):
    if data_merge == True:
        # Read
        gene_intersect = {} # unused here
        adatas = {}
        for dataset in adatas_dict.keys():
            if 'https' in adatas_dict[dataset]:
                print('Loading anndata from web source')
                adatas[dataset] = sc.read('./temp_adata.h5ad',backup_url=adatas_dict[dataset])
            adatas[dataset] = sc.read(data[dataset])
            adatas[dataset].var_names_make_unique()
            adatas[dataset].obs['dataset_merge'] = dataset
            adatas[dataset].obs['dataset_merge'] = dataset
            gene_intersect[dataset] = list(adatas[dataset].var.index)
        adata = list(adatas.values())[0].concatenate(list(adatas.values())[1:],join='inner')
        return adatas, adata
    elif data_merge == False:
        if 'https' in adatas_dict[data_key_use]:
            print('Loading anndata from web source')
            adata = sc.read('./temp_adata.h5ad',backup_url=adatas_dict[data_key_use])
        else: 
            adata = sc.read(adatas_dict[data_key_use])
    if QC_normalise == True:
        print('option to apply standardisation to data detected, performing basic QC filtering')
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
        sc.pp.log1p(adata)
        
    return adata

# resource usage logger
class DisplayCPU(threading.Thread):
    def run(self):
        tracemalloc.start()
        starting, starting_peak = tracemalloc.get_traced_memory()
        self.running = True
        self.starting = starting
        currentProcess = psutil.Process()
        cpu_pct = []
        peak_cpu = 0
        while self.running:
            peak_cpu = 0
#           time.sleep(3)
#             print('CPU % usage = '+''+ str(currentProcess.cpu_percent(interval=1)))
#             cpu_pct.append(str(currentProcess.cpu_percent(interval=1)))
            cpu = currentProcess.cpu_percent()
        # track the peak utilization of the process
            if cpu > peak_cpu:
                peak_cpu = cpu
                peak_cpu_per_core = peak_cpu/psutil.cpu_count()
        self.peak_cpu = peak_cpu
        self.peak_cpu_per_core = peak_cpu_per_core
        
    def stop(self):
#        cpu_pct = DisplayCPU.run(self)
        self.running = False
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        return current, peak
    
# projection module
def reference_projection(adata, model, dyn_std,partial_scale):
    
    class adata_temp:
        pass
    from sklearn.preprocessing import StandardScaler
    print('Determining model flavour')
    try:
        model_lr =  model['Model']
        print('Consuming celltypist model')
    except:# hasattr(model, 'coef_'):
        print('Consuming non-celltypist model')
        model_lr =  model
    print(model_lr)
    
#     model_lr =  model['Model']

    if train_x_partition == 'X':
        print('Matching reference genes in the model')
        k_x = np.isin(list(adata.var.index), list(model_lr.features))
        if k_x.sum() == 0:
            raise ValueError(f"🛑 No features overlap with the model. Please provide gene symbols")
        print(f"🧬 {k_x.sum()} features used for prediction")
        #slicing adata
        k_x_idx = np.where(k_x)[0]
        # adata_temp = adata[:,k_x_idx]
        adata_temp.var = adata[:,k_x_idx].var
        adata_temp.X = adata[:,k_x_idx].X
        adata_temp.obs = adata[:,k_x_idx].obs
        lr_idx = pd.DataFrame(model_lr.features, columns=['features']).reset_index().set_index('features').loc[list(adata_temp.var.index)].values
        # adata_arr = adata_temp.X[:,list(lr_idexes['index'])]

        # slice and reorder model
        ni, fs, cf = model_lr.n_features_in_, model_lr.features, model_lr.coef_
        model_lr.n_features_in_ = lr_idx.size
        model_lr.features = np.array(model_lr.features)[lr_idx]
        model_lr.coef_ = np.squeeze(model_lr.coef_[:,lr_idx]) #model_lr.coef_[:, lr_idx]
        
        if partial_scale == True:
            print('scaling input data, default option is to use incremental learning and fit in mini bulks!')
            # Partial scaling alg
            scaler = StandardScaler(with_mean=False)
            n = adata_temp.X.shape[0]  # number of rows
            # set dyn scale packet size
            x_len = len(adata_temp.var)
            y_len = len(adata.obs)
            if y_len < 100000:
                dyn_pack = int(x_len/10)
                pack_size = dyn_pack
            else:
                # 10 pack for every 100,000
                dyn_pack = int((y_len/100000)*10)
                pack_size = int(x_len/dyn_pack)

            batch_size =  1000#pack_size#500  # number of rows in each call to partial_fit
            index = 0  # helper-var
            while index < n:
                partial_size = min(batch_size, n - index)  # needed because last loop is possibly incomplete
                partial_x = adata_temp.X[index:index+partial_size]
                scaler.partial_fit(partial_x)
                index += partial_size
            adata_temp.X = scaler.transform(adata_temp.X)
    
    # model projections
    print('Starting reference projection!')
    if train_x_partition == 'X':
        train_x = adata_temp.X
        pred_out = pd.DataFrame(model_lr.predict(train_x),columns = ['predicted'],index = list(adata.obs.index))
        proba =  pd.DataFrame(model_lr.predict_proba(train_x),columns = model_lr.classes_,index = list(adata.obs.index))
        pred_out = pred_out.join(proba)
        
    elif train_x_partition in list(adata.obsm.keys()): 
        print('{low_dim: this partition modality is still under development!}')
        train_x = adata.obsm[train_x_partition]
        pred_out = pd.DataFrame(model_lr.predict(train_x),columns = ['predicted'],index = list(adata.obs.index))
        proba =  pd.DataFrame(model_lr.predict_proba(train_x),columns = model_lr.classes_,index = list(adata.obs.index))
        pred_out = pred_out.join(proba)
    
    else:
        print('{this partition modality is still under development!}')
    ## insert modules for low dim below

    # Simple dynamic confidence calling
    pred_out['confident_calls'] = pred_out['predicted']
    pred_out.loc[pred_out.max(axis=1)<(pred_out.mean(axis=1) + (1*pred_out.std(axis=1))),'confident_calls'] = pred_out.loc[pred_out.max(axis=1)<(pred_out.mean(axis=1) + (1*pred_out.std(axis=1))),'confident_calls'].astype(str) + '_uncertain'
    # means_ = self.model.scaler.mean_[lr_idx] if self.model.scaler.with_mean else 0
    return(pred_out,train_x,model_lr,adata_temp)

def freq_redist_68CI(adata,clusters_reassign):
    if freq_redist != False:
        print('Frequency redistribution commencing')
        cluster_prediction = "consensus_clus_prediction"
        lr_predicted_col = 'predicted'
        pred_out[clusters_reassign] = adata.obs[clusters_reassign].astype(str)
        reassign_classes = list(pred_out[clusters_reassign].unique())
        lm = 1 # lambda value
        pred_out[cluster_prediction] = pred_out[clusters_reassign]
        for z in pred_out[clusters_reassign][pred_out[clusters_reassign].isin(reassign_classes)].unique():
            df = pred_out
            df = df[(df[clusters_reassign].isin([z]))]
            df_count = pd.DataFrame(df[lr_predicted_col].value_counts())
            # Look for classificationds > 68CI
            if len(df_count) > 1:
                df_count_temp = df_count[df_count[lr_predicted_col]>int(int(df_count.mean()) + (df_count.std()*lm))]
                if len(df_count_temp >= 1):
                    df_count = df_count_temp
            #print(df_count)     
            freq_arranged = df_count.index
            cat = freq_arranged[0]
        #Make the cluster assignment first
            pred_out[cluster_prediction] = pred_out[cluster_prediction].astype(str)
            pred_out.loc[pred_out[clusters_reassign] == z, [cluster_prediction]] = cat
        # Create assignments for any classification >68CI
            for cats in freq_arranged:
                #print(cats)
                cats_assignment = cats#.replace(data1,'') + '_clus_prediction'
                pred_out.loc[(pred_out[clusters_reassign] == z) & (pred_out[lr_predicted_col] == cats),[cluster_prediction]] = cats_assignment
        min_counts = pd.DataFrame((pred_out[cluster_prediction].value_counts()))
        reassign = list(min_counts.index[min_counts[cluster_prediction]<=2])
        pred_out[cluster_prediction] = pred_out[cluster_prediction].str.replace(str(''.join(reassign)),str(''.join(pred_out.loc[pred_out[clusters_reassign].isin(list(pred_out.loc[(pred_out[cluster_prediction].isin(reassign)),clusters_reassign])),lr_predicted_col].value_counts().head(1).index.values)))
        return pred_out

### Feature importance notes
#- If we increase the x feature one unit, then the prediction will change e to the power of its weight. We can apply this rule to the all weights to find the feature importance.
#- We will calculate the Euler number to the power of its coefficient to find the importance.
#- To sum up an increase of x feature by one unit increases the odds of being versicolor class by a factor of x[importance] when all other features remain the same.

#- For low-dim, we look at the distribution of e^coef per class, we extract the 


# class coef_extract:
#     def __init__(self, model,features, pos):
# #         self.w = list(itertools.chain(*(model.coef_[pos]).tolist())) #model.coef_[pos]
#         self.w = model.coef_[class_pred_pos]
#         self.features = features 

def long_format_features(top_loadings):
    p = top_loadings.loc[:, top_loadings.columns.str.endswith("_e^coef")]
    p = pd.melt(p)
    n = top_loadings.loc[:, top_loadings.columns.str.endswith("_feature")]
    n = pd.melt(n)
    l = top_loadings.loc[:, top_loadings.columns.str.endswith("_coef")]
    l = pd.melt(l)
    n = n.replace(regex=r'_feature', value='')
    n = n.rename(columns={"variable": "class", "value": "feature"})
    p = (p.drop(["variable"],axis = 1)).rename(columns={ "value": "e^coef"})
    l = (l.drop(["variable"],axis = 1)).rename(columns={ "value": "coef"})
    concat = pd.concat([n,p,l],axis=1)
    return concat

def model_feature_sf(long_format_feature_importance, coef_use):
        long_format_feature_importance[str(coef_use) + '_pval'] = 'NaN'
        for class_lw in long_format_feature_importance['class'].unique():
            df_loadings = long_format_feature_importance[long_format_feature_importance['class'].isin([class_lw])]
            comps = coef_use #'e^coef'
            U = np.mean(df_loadings[comps])
            std = np.std(df_loadings[comps])
            med =  np.median(df_loadings[comps])
            mad = np.median(np.absolute(df_loadings[comps] - np.median(df_loadings[comps])))
            # Survival function scaled by 1.4826 of MAD (approx norm)
            pvals = scipy.stats.norm.sf(df_loadings[comps], loc=med, scale=1.4826*mad) # 95% CI of MAD <10,000 samples
            #pvals = scipy.stats.norm.sf(df_loadings[comps], loc=U, scale=1*std)
            df_loadings[str(comps) +'_pval'] = pvals
            long_format_feature_importance.loc[long_format_feature_importance.index.isin(df_loadings.index)] = df_loadings
        long_format_feature_importance['is_significant_sf'] = False
        long_format_feature_importance.loc[long_format_feature_importance[coef_use+ '_pval']<0.05,'is_significant_sf'] = True
        return long_format_feature_importance
# Apply SF to e^coeff mat data
#         pval_mat = pd.DataFrame(columns = mat.columns)
#         for class_lw in mat.index:
#             df_loadings = mat.loc[class_lw]
#             U = np.mean(df_loadings)
#             std = np.std(df_loadings)
#             med =  np.median(df_loadings)
#             mad = np.median(np.absolute(df_loadings - np.median(df_loadings)))
#             pvals = scipy.stats.norm.sf(df_loadings, loc=med, scale=1.96*U)

class estimate_important_features: # This calculates feature effect sizes of the model
    def __init__(self, model, top_n):
        print('Estimating feature importance')
        classes =  list(model.classes_)
         # get feature names
        try:
            model_features = list(itertools.chain(*list(model.features)))
        except:
            warnings.warn('no features recorded in data, naming features by position')
            print('if low-dim lr was submitted, run linear decoding function to obtain true feature set')
            model_features = list(range(0,model.coef_.shape[1]))
            model.features = model_features
        print('Calculating the Euler number to the power of coefficients')
        impt_ = pow(math.e,model.coef_)
        try:
            self.euler_pow_mat = pd.DataFrame(impt_,columns = list(itertools.chain(*list(model.features))),index = list(model.classes_))
        except:
            self.euler_pow_mat = pd.DataFrame(impt_,columns = list(model.features),index = list(model.classes_))
        self.top_n_features = pd.DataFrame(index = list(range(0,top_n)))
        # estimate per class feature importance
        
        print('Estimating feature importance for each class')
        mat = self.euler_pow_mat
        for class_pred_pos in list(range(0,len(mat.T.columns))):
            class_pred = list(mat.T.columns)[class_pred_pos]
            #     print(class_pred)
            temp_mat =  pd.DataFrame(mat.T[class_pred])
            temp_mat['coef'] = model.coef_[class_pred_pos]
            temp_mat = temp_mat.sort_values(by = [class_pred], ascending=False)
            temp_mat = temp_mat.reset_index()
            temp_mat.columns = ['feature','e^coef','coef']
            temp_mat = temp_mat[['feature','e^coef','coef']]
            temp_mat.columns =str(class_pred)+ "_" + temp_mat.columns
            self.top_n_features = pd.concat([self.top_n_features,temp_mat.head(top_n)], join="inner",ignore_index = False, axis=1)
            self.to_n_features_long = model_feature_sf(long_format_features(self.top_n_features),'e^coef')
            
    
    # plot class-wise features
def model_class_feature_plots(top_loadings, classes, comps):
    import matplotlib.pyplot as plt
    for class_temp in classes:
        class_lw = class_temp
        long_format = top_loadings
        df_loadings = long_format[long_format['class'].isin([class_lw])]
        plt.hist(df_loadings[comps])
        for i in ((df_loadings[comps][df_loadings[str(comps) +'_pval']<0.05]).unique()):
            plt.axvline(x=i,color='red')
        med = np.median(df_loadings[comps])
        plt.axvline(x=med,color='blue')
        plt.xlabel('feature_importance', fontsize=12)
        plt.title(class_lw)
        #plt.axvline(x=med,color='pink')
        df_loadings[comps][df_loadings[str(comps) +'_pval']<0.05]
        print(len(df_loadings[comps][df_loadings[str(comps) +'_pval']<0.05]))
        #Plot feature ranking
        plot_loading = pd.DataFrame(pd.DataFrame(df_loadings[comps][df_loadings[str(comps) +'_pval']<0.05]).iloc[:,0].sort_values(ascending=False))
        table = plt.table(cellText=plot_loading.values,colWidths = [1]*len(plot_loading.columns),
        rowLabels= list(df_loadings['feature'][df_loadings.index.isin(plot_loading.index)].reindex(plot_loading.index)), #plot_loading.index,
        colLabels=plot_loading.columns,
        cellLoc = 'center', rowLoc = 'center',
        loc='right', bbox=[1.4, -0.05, 0.5,1])
        table.scale(1, 2)
        table.set_fontsize(10)
        
def report_f1(model,train_x, train_label):
    ## Report accuracy score
    from sklearn.model_selection import cross_val_score
    from sklearn.model_selection import RepeatedStratifiedKFold
    from sklearn import metrics
    import seaborn as sn
    import pandas as pd
    from sklearn.metrics import confusion_matrix
    import matplotlib.pyplot as plt
    
    # cv = RepeatedStratifiedKFold(n_splits=2, n_repeats=2, random_state=1)
    # # evaluate the model and collect the scores
    # n_scores = cross_val_score(lr, train_x, train_label, scoring='accuracy', cv=cv, n_jobs=-1)
    # # report the model performance
    # print('Mean Accuracy: %.3f (%.3f)' % (np.mean(n_scores), np.std(n_scores)))

    # Report Precision score
    metric = pd.DataFrame((metrics.classification_report(train_label, model.predict(train_x), digits=2,output_dict=True))).T
    cm = confusion_matrix(train_label, model.predict(train_x))
    #cm = confusion_matrix(train_label, model.predict_proba(train_x))
    df_cm = pd.DataFrame(cm, index = model.classes_,columns = model.classes_)
    df_cm = (df_cm / df_cm.sum(axis=0))*100
    plt.figure(figsize = (20,15))
    sn.set(font_scale=1) # for label size
    pal = sns.diverging_palette(240, 10, n=10)
    #plt.suptitle(('Mean Accuracy 5 fold: %.3f std: %.3f' % (np.mean(n_scores),  np.std(n_scores))), y=1.05, fontsize=18)
    #Plot precision recall and recall
    table = plt.table(cellText=metric.values,colWidths = [1]*len(metric.columns),
    rowLabels=metric.index,
    colLabels=metric.columns,
    cellLoc = 'center', rowLoc = 'center',
    loc='bottom', bbox=[0.25, -0.6, 0.5, 0.3])
    table.scale(1, 2)
    table.set_fontsize(10)

    sn.heatmap(df_cm, annot=True, annot_kws={"size": 16},cmap=pal) # font size
    print(metrics.classification_report(train_label, model.predict(train_x), digits=2))

def subset_top_hvgs(adata_lognorm, n_top_genes):
    dispersion_norm = adata_lognorm.var['dispersions_norm'].values.astype('float32')

    dispersion_norm = dispersion_norm[~np.isnan(dispersion_norm)]
    dispersion_norm[
                ::-1
            ].sort()  # interestingly, np.argpartition is slightly slower

    disp_cut_off = dispersion_norm[n_top_genes - 1]
    gene_subset = adata_lognorm.var['dispersions_norm'].values >= disp_cut_off
    return(adata_lognorm[:,gene_subset])

def prep_scVI(adata, 
              n_hvgs = 5000,
              remove_cc_genes = True,
              remove_tcr_bcr_genes = False
             ):
    ## Remove cell cycle genes
    if remove_cc_genes:
        adata = panfetal_utils.remove_geneset(adata,genes.cc_genes)

    ## Remove TCR/BCR genes
    if remove_tcr_bcr_genes:
        adata = panfetal_utils.remove_geneset(adata, genes.IG_genes)
        adata = panfetal_utils.remove_geneset(adata, genes.TCR_genes)
        
    ## HVG selection
    adata = subset_top_hvgs(adata, n_top_genes=n_hvgs)
    return(adata)

# Modified LR train module, does not work with low-dim by default anymore, please use low-dim adapter
def LR_train(adata, train_x, train_label, penalty='elasticnet', sparcity=0.2,max_iter=200,l1_ratio =0.2,tune_hyper_params =False,n_splits=5, n_repeats=3,l1_grid = [0.2,0.5,0.8], c_grid = [0.2,0.4,0.6]):
    if tune_hyper_params == True:
        train_labels=train_label
        results = tune_lr_model(adata, train_x_partition = train_x, random_state = 42,  train_labels = train_labels, n_splits=n_splits, n_repeats=n_repeats,l1_grid = l1_grid, c_grid = c_grid)
        print('hyper_params tuned')
        sparcity = results.best_params_['C']
        l1_ratio = results.best_params_['l1_ratio']
    
    lr = LogisticRegression(penalty = penalty, C = sparcity, max_iter =  max_iter, n_jobs=thread_num)
    if (penalty == "l1"):
        lr = LogisticRegression(penalty = penalty, C = sparcity, max_iter =  max_iter, dual = True, solver = 'liblinear',multi_class = 'ovr', n_jobs=thread_num ) # one-vs-rest
    if (penalty == "elasticnet"):
        lr = LogisticRegression(penalty = penalty, C = sparcity, max_iter =  max_iter, dual=False,solver = 'saga',l1_ratio=l1_ratio,multi_class = 'ovr', n_jobs=thread_num)
    if train_x == 'X':
        subset_train = adata.obs.index
        # Define training parameters
        train_label = adata.obs[train_label].values
#        predict_label = train_label[subset_predict]
#        train_label = train_label[subset_train]
        train_x = adata.X#[adata.obs.index.isin(list(adata.obs[subset_train].index))]
#        predict_x = adata.X[adata.obs.index.isin(list(adata.obs[subset_predict].index))]
    elif train_x in adata.obsm.keys():
        # Define training parameters
        train_label = adata.obs[train_label].values
#        predict_label = train_label[subset_predict]
#         train_label = train_label[subset_train]
        train_x = adata.obsm[train_x]
#        predict_x = train_x
#        train_x = train_x[subset_train, :]
        # Define prediction parameters
#        predict_x = predict_x[subset_predict]
#        predict_x = pd.DataFrame(predict_x)
#        predict_x.index = adata.obs[subset_predict].index
    # Train predictive model using user defined partition labels (train_x ,train_label, predict_x)
    model = lr.fit(train_x, train_label)
    model.features = np.array(adata.var.index)
    return model

def tune_lr_model(adata, train_x_partition = 'X', random_state = 42,  train_labels = None, n_splits=5, n_repeats=3,l1_grid = [0.1,0.2,0.5,0.8], c_grid = [0.1,0.2,0.4,0.6]):
    import bless as bless
    from sklearn.gaussian_process.kernels import RBF
    from numpy import arange
    from sklearn.model_selection import RepeatedKFold
    from sklearn.datasets import make_classification
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import f1_score
    from sklearn.model_selection import GridSearchCV

    # If latent rep is provided, randomly sample data in spatially aware manner for initialisation
    r = np.random.RandomState(random_state)
    if train_x_partition in adata.obsm.keys():
        lvg = bless.bless(tune_train_x, RBF(length_scale=20), lam_final = 2, qbar = 2, random_state = r, H = 10, force_cpu=True)
    #     try:
    #         import cupy
    #         lvg_2 = bless(adata.obsm[train_x_partition], RBF(length_scale=10), 10, 10, r, 10, force_cpu=False)
    #     except ImportError:
    #         print("cupy not found, defaulting to numpy")
        adata_tuning = adata[lvg.idx]
        tune_train_x = adata_tuning.obsm[train_x_partition][:]
    else:
        print('no latent representation provided, random sampling instead')
        prop = 0.1
        random_vertices = []
        n_ixs = int(len(adata.obs) * prop)
        random_vertices = random.sample(list(range(len(adata.obs))), k=n_ixs)
        adata_tuning = adata[random_vertices]
        tune_train_x = adata_tuning.X
        
    if not train_labels == None:
        tune_train_label = adata_tuning.obs[train_labels]
    elif train_labels == None:
        try:
            print('no training labels provided, defaulting to unsuperived leiden clustering, updates will change this to voronoi greedy sampling')
            sc.tl.leiden(adata_tuning)
        except:
            print('no training labels provided, no neighbors, defaulting to unsuperived leiden clustering, updates will change this to voronoi greedy sampling')
            sc.pp.neighbors(adata_hm, n_neighbors=15, n_pcs=50)
            sc.tl.leiden(adata_tuning)
        tune_train_label = adata_tuning.obs['leiden']
    ## tune regularization for multinomial logistic regression
    print('starting tuning loops')
    X = tune_train_x
    y = tune_train_label
    grid = dict()
    # define model
    cv = RepeatedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=random_state)
    #model = LogisticRegression(penalty = penalty, max_iter =  200, dual=False,solver = 'saga', multi_class = 'multinomial',)
    model = LogisticRegression(penalty = penalty, C = sparcity, max_iter =  100, n_jobs=4)
    if (penalty == "l1"):
        model = LogisticRegression(penalty = penalty, C = sparcity, max_iter =  100, dual = True, solver = 'liblinear',multi_class = 'ovr', n_jobs=4 ) # one-vs-rest
    if (penalty == "elasticnet"):
        model = LogisticRegression(penalty = penalty, C = sparcity, max_iter =  100, dual=False,solver = 'saga',l1_ratio=l1_ratio,multi_class = 'ovr', n_jobs=4) # use multinomial class if probabilities are descrete
        grid['l1_ratio'] = l1_grid
    grid['C'] = c_grid
    # define search
    search = GridSearchCV(model, grid, scoring='neg_mean_absolute_error', cv=cv, n_jobs=-1)
    # perform the search
    results = search.fit(X, y)
    # summarize
    print('MAE: %.3f' % results.best_score_)
    print('Config: %s' % results.best_params_)
    return results

def prep_training_data(adata_temp,feat_use,batch_key, model_key, batch_correction=False, var_length = 7500,penalty='elasticnet',sparcity=0.2,max_iter = 200,l1_ratio = 0.1,partial_scale=True,train_x_partition ='X',theta = 3,tune_hyper_params=False ):
    model_name = model_key + '_lr_model'
    print('performing highly variable gene selection')
    sc.pp.highly_variable_genes(adata_temp, batch_key = batch_key, subset=False)
    adata_temp = subset_top_hvgs(adata_temp,var_length)
    #scale the input data
    if partial_scale == True:
        print('scaling input data, default option is to use incremental learning and fit in mini bulks!')
        # Partial scaling alg
        #adata_temp.X = (adata_temp.X)
        scaler = StandardScaler(with_mean=False)
        n = adata_temp.X.shape[0]  # number of rows
        # set dyn scale packet size
        x_len = len(adata_temp.var)
        y_len = len(adata_temp.obs)
        if y_len < 100000:
            dyn_pack = int(x_len/10)
            pack_size = dyn_pack
        else:
            # 10 pack for every 100,000
            dyn_pack = int((y_len/100000)*10)
            pack_size = int(x_len/dyn_pack)
        batch_size =  1000#pack_size#500  # number of rows in each call to partial_fit
        index = 0  # helper-var
        while index < n:
            partial_size = min(batch_size, n - index)  # needed because last loop is possibly incomplete
            partial_x = adata_temp.X[index:index+partial_size]
            scaler.partial_fit(partial_x)
            index += partial_size
        adata_temp.X = scaler.transform(adata_temp.X)
#     else:
#         sc.pp.scale(adata_temp, zero_center=True, max_value=None, copy=False, layer=None, obsm=None)
    if (train_x_partition != 'X') & (train_x_partition in adata_temp.obsm.keys()):
        print('train partition is not in OBSM, defaulting to PCA')
        # Now compute PCA
        sc.pp.pca(adata_temp, n_comps=100, use_highly_variable=True, svd_solver='arpack')
        sc.pl.pca_variance_ratio(adata_temp, log=True,n_pcs=100)
        
        # Batch correction options
        # The script will test later which Harmony values we should use 
        if(batch_correction == "Harmony"):
            print("Commencing harmony")
            adata_temp.obs['lr_batch'] = adata_temp.obs[batch_key]
            batch_var = "lr_batch"
            # Create hm subset
            adata_hm = adata_temp[:]
            # Set harmony variables
            data_mat = np.array(adata_hm.obsm["X_pca"])
            meta_data = adata_hm.obs
            vars_use = [batch_var]
            # Run Harmony
            ho = hm.run_harmony(data_mat, meta_data, vars_use,theta=theta)
            res = (pd.DataFrame(ho.Z_corr)).T
            res.columns = ['X{}'.format(i + 1) for i in range(res.shape[1])]
            # Insert coordinates back into object
            adata_hm.obsm["X_pca_back"]= adata_hm.obsm["X_pca"][:]
            adata_hm.obsm["X_pca"] = np.array(res)
            # Run neighbours
            #sc.pp.neighbors(adata_hm, n_neighbors=15, n_pcs=50)
            adata_temp = adata_hm[:]
            del adata_hm
        elif(batch_correction == "BBKNN"):
            print("Commencing BBKNN")
            sc.external.pp.bbknn(adata_temp, batch_key=batch_var, approx=True, metric='angular', copy=False, n_pcs=50, trim=None, n_trees=10, use_faiss=True, set_op_mix_ratio=1.0, local_connectivity=15) 
        print("adata1 and adata2 are now combined and preprocessed in 'adata' obj - success!")


    # train model
#    train_x = adata_temp.X
    #train_label = adata_temp.obs[feat_use]
    print('proceeding to train model')
    model = LR_train(adata_temp, train_x = train_x_partition, train_label=feat_use, penalty=penalty, sparcity=sparcity,max_iter=max_iter,l1_ratio = l1_ratio,tune_hyper_params = tune_hyper_params)
    model.features = list(adata_temp.var.index)
    return model

def regression_results(y_true, y_pred):
    # Regression metrics
    explained_variance=metrics.explained_variance_score(y_true, y_pred)
    mean_absolute_error=metrics.mean_absolute_error(y_true, y_pred) 
    mse=metrics.mean_squared_error(y_true, y_pred) 
    mean_squared_log_error=metrics.mean_squared_log_error(y_true, y_pred)
    median_absolute_error=metrics.median_absolute_error(y_true, y_pred)
    r2=metrics.r2_score(y_true, y_pred)
    print('explained_variance: ', round(explained_variance,4))    
    print('mean_squared_log_error: ', round(mean_squared_log_error,4))
    print('r2: ', round(r2,4))
    print('MAE: ', round(mean_absolute_error,4))
    print('MSE: ', round(mse,4))
    print('RMSE: ', round(np.sqrt(mse),4))

from collections import Counter
from collections import defaultdict
import scanpy as sc
import pandas as pd
import pickle as pkl
import numpy as np
import scipy
import matplotlib.pyplot as plt
import re
import glob
import os
import sys
#from geosketch import gs
from numpy import cov
import scipy.cluster.hierarchy as spc
import seaborn as sns; sns.set(color_codes=True)
from sklearn.linear_model import LogisticRegression
import sklearn
from pathlib import Path
import requests
import psutil
import random
import threading
import tracemalloc
import itertools
import math
import warnings
import sklearn.metrics as metrics

def load_models(model_dict,model_run):
    if (Path(model_dict[model_run])).is_file():
        # Load data (deserialize)
        model = pkl.load(open(model_dict[model_run], "rb"))
        return model
    elif 'http' in model_dict[model_run]:
        print('Loading model from web source')
        r_get = requests.get(model_dict[model_run])
        fpath = './model_temp.sav'
        open(fpath , 'wb').write(r_get.content)
        model = pkl.load(open(fpath, "rb"))
        return model

def load_adatas(adatas_dict,data_merge, data_key_use,QC_normalise):
    if data_merge == True:
        # Read
        gene_intersect = {} # unused here
        adatas = {}
        for dataset in adatas_dict.keys():
            if 'https' in adatas_dict[dataset]:
                print('Loading anndata from web source')
                adatas[dataset] = sc.read('./temp_adata.h5ad',backup_url=adatas_dict[dataset])
            adatas[dataset] = sc.read(data[dataset])
            adatas[dataset].var_names_make_unique()
            adatas[dataset].obs['dataset_merge'] = dataset
            adatas[dataset].obs['dataset_merge'] = dataset
            gene_intersect[dataset] = list(adatas[dataset].var.index)
        adata = list(adatas.values())[0].concatenate(list(adatas.values())[1:],join='inner')
        return adatas, adata
    elif data_merge == False:
        if 'https' in adatas_dict[data_key_use]:
            print('Loading anndata from web source')
            adata = sc.read('./temp_adata.h5ad',backup_url=adatas_dict[data_key_use])
        else: 
            adata = sc.read(adatas_dict[data_key_use])
    if QC_normalise == True:
        print('option to apply standardisation to data detected, performing basic QC filtering')
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
        sc.pp.log1p(adata)
        
    return adata

# resource usage logger
class DisplayCPU(threading.Thread):
    def run(self):
        tracemalloc.start()
        starting, starting_peak = tracemalloc.get_traced_memory()
        self.running = True
        self.starting = starting
        currentProcess = psutil.Process()
        cpu_pct = []
        peak_cpu = 0
        while self.running:
            peak_cpu = 0
#           time.sleep(3)
#             print('CPU % usage = '+''+ str(currentProcess.cpu_percent(interval=1)))
#             cpu_pct.append(str(currentProcess.cpu_percent(interval=1)))
            cpu = currentProcess.cpu_percent()
        # track the peak utilization of the process
            if cpu > peak_cpu:
                peak_cpu = cpu
                peak_cpu_per_core = peak_cpu/psutil.cpu_count()
        self.peak_cpu = peak_cpu
        self.peak_cpu_per_core = peak_cpu_per_core
        
    def stop(self):
#        cpu_pct = DisplayCPU.run(self)
        self.running = False
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        return current, peak
    
# projection module
def reference_projection(adata, model, dyn_std,partial_scale):
    
    class adata_temp:
        pass
    from sklearn.preprocessing import StandardScaler
    print('Determining model flavour')
    try:
        model_lr =  model['Model']
        print('Consuming celltypist model')
    except:# hasattr(model, 'coef_'):
        print('Consuming non-celltypist model')
        model_lr =  model
    print(model_lr)
    
#     model_lr =  model['Model']

    if train_x_partition == 'X':
        print('Matching reference genes in the model')
        k_x = np.isin(list(adata.var.index), list(model_lr.features))
        if k_x.sum() == 0:
            raise ValueError(f"🛑 No features overlap with the model. Please provide gene symbols")
        print(f"🧬 {k_x.sum()} features used for prediction")
        #slicing adata
        k_x_idx = np.where(k_x)[0]
        # adata_temp = adata[:,k_x_idx]
        adata_temp.var = adata[:,k_x_idx].var
        adata_temp.X = adata[:,k_x_idx].X
        adata_temp.obs = adata[:,k_x_idx].obs
        lr_idx = pd.DataFrame(model_lr.features, columns=['features']).reset_index().set_index('features').loc[list(adata_temp.var.index)].values
        # adata_arr = adata_temp.X[:,list(lr_idexes['index'])]

        # slice and reorder model
        ni, fs, cf = model_lr.n_features_in_, model_lr.features, model_lr.coef_
        model_lr.n_features_in_ = lr_idx.size
        model_lr.features = np.array(model_lr.features)[lr_idx]
        model_lr.coef_ = np.squeeze(model_lr.coef_[:,lr_idx]) #model_lr.coef_[:, lr_idx]
        
        if partial_scale == True:
            print('scaling input data, default option is to use incremental learning and fit in mini bulks!')
            # Partial scaling alg
            scaler = StandardScaler(with_mean=False)
            n = adata_temp.X.shape[0]  # number of rows
            # set dyn scale packet size
            x_len = len(adata_temp.var)
            y_len = len(adata.obs)
            if y_len < 100000:
                dyn_pack = int(x_len/10)
                pack_size = dyn_pack
            else:
                # 10 pack for every 100,000
                dyn_pack = int((y_len/100000)*10)
                pack_size = int(x_len/dyn_pack)

            batch_size =  1000#pack_size#500  # number of rows in each call to partial_fit
            index = 0  # helper-var
            while index < n:
                partial_size = min(batch_size, n - index)  # needed because last loop is possibly incomplete
                partial_x = adata_temp.X[index:index+partial_size]
                scaler.partial_fit(partial_x)
                index += partial_size
            adata_temp.X = scaler.transform(adata_temp.X)
    
    # model projections
    print('Starting reference projection!')
    if train_x_partition == 'X':
        train_x = adata_temp.X
        pred_out = pd.DataFrame(model_lr.predict(train_x),columns = ['predicted'],index = list(adata.obs.index))
        proba =  pd.DataFrame(model_lr.predict_proba(train_x),columns = model_lr.classes_,index = list(adata.obs.index))
        pred_out = pred_out.join(proba)
        
    elif train_x_partition in list(adata.obsm.keys()): 
        print('{low_dim: this partition modality is still under development!}')
        train_x = adata.obsm[train_x_partition]
        pred_out = pd.DataFrame(model_lr.predict(train_x),columns = ['predicted'],index = list(adata.obs.index))
        proba =  pd.DataFrame(model_lr.predict_proba(train_x),columns = model_lr.classes_,index = list(adata.obs.index))
        pred_out = pred_out.join(proba)
    
    else:
        print('{this partition modality is still under development!}')
    ## insert modules for low dim below

    # Simple dynamic confidence calling
    pred_out['confident_calls'] = pred_out['predicted']
    pred_out.loc[pred_out.max(axis=1)<(pred_out.mean(axis=1) + (1*pred_out.std(axis=1))),'confident_calls'] = pred_out.loc[pred_out.max(axis=1)<(pred_out.mean(axis=1) + (1*pred_out.std(axis=1))),'confident_calls'].astype(str) + '_uncertain'
    # means_ = self.model.scaler.mean_[lr_idx] if self.model.scaler.with_mean else 0
    return(pred_out,train_x,model_lr,adata_temp)

def freq_redist_68CI(adata,clusters_reassign):
    if freq_redist != False:
        print('Frequency redistribution commencing')
        cluster_prediction = "consensus_clus_prediction"
        lr_predicted_col = 'predicted'
        pred_out[clusters_reassign] = adata.obs[clusters_reassign].astype(str)
        reassign_classes = list(pred_out[clusters_reassign].unique())
        lm = 1 # lambda value
        pred_out[cluster_prediction] = pred_out[clusters_reassign]
        for z in pred_out[clusters_reassign][pred_out[clusters_reassign].isin(reassign_classes)].unique():
            df = pred_out
            df = df[(df[clusters_reassign].isin([z]))]
            df_count = pd.DataFrame(df[lr_predicted_col].value_counts())
            # Look for classificationds > 68CI
            if len(df_count) > 1:
                df_count_temp = df_count[df_count[lr_predicted_col]>int(int(df_count.mean()) + (df_count.std()*lm))]
                if len(df_count_temp >= 1):
                    df_count = df_count_temp
            #print(df_count)     
            freq_arranged = df_count.index
            cat = freq_arranged[0]
        #Make the cluster assignment first
            pred_out[cluster_prediction] = pred_out[cluster_prediction].astype(str)
            pred_out.loc[pred_out[clusters_reassign] == z, [cluster_prediction]] = cat
        # Create assignments for any classification >68CI
            for cats in freq_arranged:
                #print(cats)
                cats_assignment = cats#.replace(data1,'') + '_clus_prediction'
                pred_out.loc[(pred_out[clusters_reassign] == z) & (pred_out[lr_predicted_col] == cats),[cluster_prediction]] = cats_assignment
        min_counts = pd.DataFrame((pred_out[cluster_prediction].value_counts()))
        reassign = list(min_counts.index[min_counts[cluster_prediction]<=2])
        pred_out[cluster_prediction] = pred_out[cluster_prediction].str.replace(str(''.join(reassign)),str(''.join(pred_out.loc[pred_out[clusters_reassign].isin(list(pred_out.loc[(pred_out[cluster_prediction].isin(reassign)),clusters_reassign])),lr_predicted_col].value_counts().head(1).index.values)))
        return pred_out

### Feature importance notes
#- If we increase the x feature one unit, then the prediction will change e to the power of its weight. We can apply this rule to the all weights to find the feature importance.
#- We will calculate the Euler number to the power of its coefficient to find the importance.
#- To sum up an increase of x feature by one unit increases the odds of being versicolor class by a factor of x[importance] when all other features remain the same.

#- For low-dim, we look at the distribution of e^coef per class, we extract the 


# class coef_extract:
#     def __init__(self, model,features, pos):
# #         self.w = list(itertools.chain(*(model.coef_[pos]).tolist())) #model.coef_[pos]
#         self.w = model.coef_[class_pred_pos]
#         self.features = features 

def long_format_features(top_loadings):
    p = top_loadings.loc[:, top_loadings.columns.str.endswith("_e^coef")]
    p = pd.melt(p)
    n = top_loadings.loc[:, top_loadings.columns.str.endswith("_feature")]
    n = pd.melt(n)
    l = top_loadings.loc[:, top_loadings.columns.str.endswith("_coef")]
    l = pd.melt(l)
    n = n.replace(regex=r'_feature', value='')
    n = n.rename(columns={"variable": "class", "value": "feature"})
    p = (p.drop(["variable"],axis = 1)).rename(columns={ "value": "e^coef"})
    l = (l.drop(["variable"],axis = 1)).rename(columns={ "value": "coef"})
    concat = pd.concat([n,p,l],axis=1)
    return concat

def model_feature_sf(long_format_feature_importance, coef_use):
        long_format_feature_importance[str(coef_use) + '_pval'] = 'NaN'
        for class_lw in long_format_feature_importance['class'].unique():
            df_loadings = long_format_feature_importance[long_format_feature_importance['class'].isin([class_lw])]
            comps = coef_use #'e^coef'
            U = np.mean(df_loadings[comps])
            std = np.std(df_loadings[comps])
            med =  np.median(df_loadings[comps])
            mad = np.median(np.absolute(df_loadings[comps] - np.median(df_loadings[comps])))
            # Survival function scaled by 1.4826 of MAD (approx norm)
            pvals = scipy.stats.norm.sf(df_loadings[comps], loc=med, scale=1.4826*mad) # 95% CI of MAD <10,000 samples
            #pvals = scipy.stats.norm.sf(df_loadings[comps], loc=U, scale=1*std)
            df_loadings[str(comps) +'_pval'] = pvals
            long_format_feature_importance.loc[long_format_feature_importance.index.isin(df_loadings.index)] = df_loadings
        long_format_feature_importance['is_significant_sf'] = False
        long_format_feature_importance.loc[long_format_feature_importance[coef_use+ '_pval']<0.05,'is_significant_sf'] = True
        return long_format_feature_importance
# Apply SF to e^coeff mat data
#         pval_mat = pd.DataFrame(columns = mat.columns)
#         for class_lw in mat.index:
#             df_loadings = mat.loc[class_lw]
#             U = np.mean(df_loadings)
#             std = np.std(df_loadings)
#             med =  np.median(df_loadings)
#             mad = np.median(np.absolute(df_loadings - np.median(df_loadings)))
#             pvals = scipy.stats.norm.sf(df_loadings, loc=med, scale=1.96*U)

class estimate_important_features: # This calculates feature effect sizes of the model
    def __init__(self, model, top_n):
        print('Estimating feature importance')
        classes =  list(model.classes_)
         # get feature names
        try:
            model_features = list(itertools.chain(*list(model.features)))
        except:
            warnings.warn('no features recorded in data, naming features by position')
            print('if low-dim lr was submitted, run linear decoding function to obtain true feature set')
            model_features = list(range(0,model.coef_.shape[1]))
            model.features = model_features
        print('Calculating the Euler number to the power of coefficients')
        impt_ = pow(math.e,model.coef_)
        try:
            self.euler_pow_mat = pd.DataFrame(impt_,columns = list(itertools.chain(*list(model.features))),index = list(model.classes_))
        except:
            self.euler_pow_mat = pd.DataFrame(impt_,columns = list(model.features),index = list(model.classes_))
        self.top_n_features = pd.DataFrame(index = list(range(0,top_n)))
        # estimate per class feature importance
        
        print('Estimating feature importance for each class')
        mat = self.euler_pow_mat
        for class_pred_pos in list(range(0,len(mat.T.columns))):
            class_pred = list(mat.T.columns)[class_pred_pos]
            #     print(class_pred)
            temp_mat =  pd.DataFrame(mat.T[class_pred])
            temp_mat['coef'] = model.coef_[class_pred_pos]
            temp_mat = temp_mat.sort_values(by = [class_pred], ascending=False)
            temp_mat = temp_mat.reset_index()
            temp_mat.columns = ['feature','e^coef','coef']
            temp_mat = temp_mat[['feature','e^coef','coef']]
            temp_mat.columns =str(class_pred)+ "_" + temp_mat.columns
            self.top_n_features = pd.concat([self.top_n_features,temp_mat.head(top_n)], join="inner",ignore_index = False, axis=1)
            self.to_n_features_long = model_feature_sf(long_format_features(self.top_n_features),'e^coef')
            
    
    # plot class-wise features
def model_class_feature_plots(top_loadings, classes, comps):
    import matplotlib.pyplot as plt
    for class_temp in classes:
        class_lw = class_temp
        long_format = top_loadings
        df_loadings = long_format[long_format['class'].isin([class_lw])]
        plt.hist(df_loadings[comps])
        for i in ((df_loadings[comps][df_loadings[str(comps) +'_pval']<0.05]).unique()):
            plt.axvline(x=i,color='red')
        med = np.median(df_loadings[comps])
        plt.axvline(x=med,color='blue')
        plt.xlabel('feature_importance', fontsize=12)
        plt.title(class_lw)
        #plt.axvline(x=med,color='pink')
        df_loadings[comps][df_loadings[str(comps) +'_pval']<0.05]
        print(len(df_loadings[comps][df_loadings[str(comps) +'_pval']<0.05]))
        #Plot feature ranking
        plot_loading = pd.DataFrame(pd.DataFrame(df_loadings[comps][df_loadings[str(comps) +'_pval']<0.05]).iloc[:,0].sort_values(ascending=False))
        table = plt.table(cellText=plot_loading.values,colWidths = [1]*len(plot_loading.columns),
        rowLabels= list(df_loadings['feature'][df_loadings.index.isin(plot_loading.index)].reindex(plot_loading.index)), #plot_loading.index,
        colLabels=plot_loading.columns,
        cellLoc = 'center', rowLoc = 'center',
        loc='right', bbox=[1.4, -0.05, 0.5,1])
        table.scale(1, 2)
        table.set_fontsize(10)
        
def report_f1(model,train_x, train_label):
    ## Report accuracy score
    from sklearn.model_selection import cross_val_score
    from sklearn.model_selection import RepeatedStratifiedKFold
    from sklearn import metrics
    import seaborn as sn
    import pandas as pd
    from sklearn.metrics import confusion_matrix
    import matplotlib.pyplot as plt
    
    # cv = RepeatedStratifiedKFold(n_splits=2, n_repeats=2, random_state=1)
    # # evaluate the model and collect the scores
    # n_scores = cross_val_score(lr, train_x, train_label, scoring='accuracy', cv=cv, n_jobs=-1)
    # # report the model performance
    # print('Mean Accuracy: %.3f (%.3f)' % (np.mean(n_scores), np.std(n_scores)))

    # Report Precision score
    metric = pd.DataFrame((metrics.classification_report(train_label, model.predict(train_x), digits=2,output_dict=True))).T
    cm = confusion_matrix(train_label, model.predict(train_x))
    #cm = confusion_matrix(train_label, model.predict_proba(train_x))
    df_cm = pd.DataFrame(cm, index = model.classes_,columns = model.classes_)
    df_cm = (df_cm / df_cm.sum(axis=0))*100
    plt.figure(figsize = (20,15))
    sn.set(font_scale=1) # for label size
    pal = sns.diverging_palette(240, 10, n=10)
    #plt.suptitle(('Mean Accuracy 5 fold: %.3f std: %.3f' % (np.mean(n_scores),  np.std(n_scores))), y=1.05, fontsize=18)
    #Plot precision recall and recall
    table = plt.table(cellText=metric.values,colWidths = [1]*len(metric.columns),
    rowLabels=metric.index,
    colLabels=metric.columns,
    cellLoc = 'center', rowLoc = 'center',
    loc='bottom', bbox=[0.25, -0.6, 0.5, 0.3])
    table.scale(1, 2)
    table.set_fontsize(10)

    sn.heatmap(df_cm, annot=True, annot_kws={"size": 16},cmap=pal) # font size
    print(metrics.classification_report(train_label, model.predict(train_x), digits=2))

def subset_top_hvgs(adata_lognorm, n_top_genes):
    dispersion_norm = adata_lognorm.var['dispersions_norm'].values.astype('float32')

    dispersion_norm = dispersion_norm[~np.isnan(dispersion_norm)]
    dispersion_norm[
                ::-1
            ].sort()  # interestingly, np.argpartition is slightly slower

    disp_cut_off = dispersion_norm[n_top_genes - 1]
    gene_subset = adata_lognorm.var['dispersions_norm'].values >= disp_cut_off
    return(adata_lognorm[:,gene_subset])

def prep_scVI(adata, 
              n_hvgs = 5000,
              remove_cc_genes = True,
              remove_tcr_bcr_genes = False
             ):
    ## Remove cell cycle genes
    if remove_cc_genes:
        adata = panfetal_utils.remove_geneset(adata,genes.cc_genes)

    ## Remove TCR/BCR genes
    if remove_tcr_bcr_genes:
        adata = panfetal_utils.remove_geneset(adata, genes.IG_genes)
        adata = panfetal_utils.remove_geneset(adata, genes.TCR_genes)
        
    ## HVG selection
    adata = subset_top_hvgs(adata, n_top_genes=n_hvgs)
    return(adata)

# Modified LR train module, does not work with low-dim by default anymore, please use low-dim adapter
def LR_train(adata, train_x, train_label, penalty='elasticnet', sparcity=0.2,max_iter=200,l1_ratio =0.2,tune_hyper_params =False,n_splits=5, n_repeats=3,l1_grid = [0.2,0.5,0.8], c_grid = [0.2,0.4,0.6]):
    if tune_hyper_params == True:
        train_labels=train_label
        results = tune_lr_model(adata, train_x_partition = train_x, random_state = 42,  train_labels = train_labels, n_splits=n_splits, n_repeats=n_repeats,l1_grid = l1_grid, c_grid = c_grid)
        print('hyper_params tuned')
        sparcity = results.best_params_['C']
        l1_ratio = results.best_params_['l1_ratio']
    
    lr = LogisticRegression(penalty = penalty, C = sparcity, max_iter =  max_iter, n_jobs=thread_num)
    if (penalty == "l1"):
        lr = LogisticRegression(penalty = penalty, C = sparcity, max_iter =  max_iter, dual = True, solver = 'liblinear',multi_class = 'ovr', n_jobs=thread_num ) # one-vs-rest
    if (penalty == "elasticnet"):
        lr = LogisticRegression(penalty = penalty, C = sparcity, max_iter =  max_iter, dual=False,solver = 'saga',l1_ratio=l1_ratio,multi_class = 'ovr', n_jobs=thread_num)
    if train_x == 'X':
        subset_train = adata.obs.index
        # Define training parameters
        train_label = adata.obs[train_label].values
#        predict_label = train_label[subset_predict]
#        train_label = train_label[subset_train]
        train_x = adata.X#[adata.obs.index.isin(list(adata.obs[subset_train].index))]
#        predict_x = adata.X[adata.obs.index.isin(list(adata.obs[subset_predict].index))]
    elif train_x in adata.obsm.keys():
        # Define training parameters
        train_label = adata.obs[train_label].values
#        predict_label = train_label[subset_predict]
#         train_label = train_label[subset_train]
        train_x = adata.obsm[train_x]
#        predict_x = train_x
#        train_x = train_x[subset_train, :]
        # Define prediction parameters
#        predict_x = predict_x[subset_predict]
#        predict_x = pd.DataFrame(predict_x)
#        predict_x.index = adata.obs[subset_predict].index
    # Train predictive model using user defined partition labels (train_x ,train_label, predict_x)
    model = lr.fit(train_x, train_label)
    model.features = np.array(adata.var.index)
    return model

def tune_lr_model(adata, train_x_partition = 'X', random_state = 42,  train_labels = None, n_splits=5, n_repeats=3,l1_grid = [0.1,0.2,0.5,0.8], c_grid = [0.1,0.2,0.4,0.6]):
    import bless as bless
    from sklearn.gaussian_process.kernels import RBF
    from numpy import arange
    from sklearn.model_selection import RepeatedKFold
    from sklearn.datasets import make_classification
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import f1_score
    from sklearn.model_selection import GridSearchCV

    # If latent rep is provided, randomly sample data in spatially aware manner for initialisation
    r = np.random.RandomState(random_state)
    if train_x_partition in adata.obsm.keys():
        lvg = bless.bless(tune_train_x, RBF(length_scale=20), lam_final = 2, qbar = 2, random_state = r, H = 10, force_cpu=True)
    #     try:
    #         import cupy
    #         lvg_2 = bless(adata.obsm[train_x_partition], RBF(length_scale=10), 10, 10, r, 10, force_cpu=False)
    #     except ImportError:
    #         print("cupy not found, defaulting to numpy")
        adata_tuning = adata[lvg.idx]
        tune_train_x = adata_tuning.obsm[train_x_partition][:]
    else:
        print('no latent representation provided, random sampling instead')
        prop = 0.1
        random_vertices = []
        n_ixs = int(len(adata.obs) * prop)
        random_vertices = random.sample(list(range(len(adata.obs))), k=n_ixs)
        adata_tuning = adata[random_vertices]
        tune_train_x = adata_tuning.X
        
    if not train_labels == None:
        tune_train_label = adata_tuning.obs[train_labels]
    elif train_labels == None:
        try:
            print('no training labels provided, defaulting to unsuperived leiden clustering, updates will change this to voronoi greedy sampling')
            sc.tl.leiden(adata_tuning)
        except:
            print('no training labels provided, no neighbors, defaulting to unsuperived leiden clustering, updates will change this to voronoi greedy sampling')
            sc.pp.neighbors(adata_hm, n_neighbors=15, n_pcs=50)
            sc.tl.leiden(adata_tuning)
        tune_train_label = adata_tuning.obs['leiden']
    ## tune regularization for multinomial logistic regression
    print('starting tuning loops')
    X = tune_train_x
    y = tune_train_label
    grid = dict()
    # define model
    cv = RepeatedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=random_state)
    #model = LogisticRegression(penalty = penalty, max_iter =  200, dual=False,solver = 'saga', multi_class = 'multinomial',)
    model = LogisticRegression(penalty = penalty, C = sparcity, max_iter =  100, n_jobs=4)
    if (penalty == "l1"):
        model = LogisticRegression(penalty = penalty, C = sparcity, max_iter =  100, dual = True, solver = 'liblinear',multi_class = 'ovr', n_jobs=4 ) # one-vs-rest
    if (penalty == "elasticnet"):
        model = LogisticRegression(penalty = penalty, C = sparcity, max_iter =  100, dual=False,solver = 'saga',l1_ratio=l1_ratio,multi_class = 'ovr', n_jobs=4) # use multinomial class if probabilities are descrete
        grid['l1_ratio'] = l1_grid
    grid['C'] = c_grid
    # define search
    search = GridSearchCV(model, grid, scoring='neg_mean_absolute_error', cv=cv, n_jobs=-1)
    # perform the search
    results = search.fit(X, y)
    # summarize
    print('MAE: %.3f' % results.best_score_)
    print('Config: %s' % results.best_params_)
    return results

def prep_training_data(adata_temp,feat_use,batch_key, model_key, batch_correction=False, var_length = 7500,penalty='elasticnet',sparcity=0.2,max_iter = 200,l1_ratio = 0.1,partial_scale=True,train_x_partition ='X',theta = 3,tune_hyper_params=False ):
    model_name = model_key + '_lr_model'
    print('performing highly variable gene selection')
    sc.pp.highly_variable_genes(adata_temp, batch_key = batch_key, subset=False)
    adata_temp = subset_top_hvgs(adata_temp,var_length)
    #scale the input data
    if partial_scale == True:
        print('scaling input data, default option is to use incremental learning and fit in mini bulks!')
        # Partial scaling alg
        #adata_temp.X = (adata_temp.X)
        scaler = StandardScaler(with_mean=False)
        n = adata_temp.X.shape[0]  # number of rows
        # set dyn scale packet size
        x_len = len(adata_temp.var)
        y_len = len(adata_temp.obs)
        if y_len < 100000:
            dyn_pack = int(x_len/10)
            pack_size = dyn_pack
        else:
            # 10 pack for every 100,000
            dyn_pack = int((y_len/100000)*10)
            pack_size = int(x_len/dyn_pack)
        batch_size =  1000#pack_size#500  # number of rows in each call to partial_fit
        index = 0  # helper-var
        while index < n:
            partial_size = min(batch_size, n - index)  # needed because last loop is possibly incomplete
            partial_x = adata_temp.X[index:index+partial_size]
            scaler.partial_fit(partial_x)
            index += partial_size
        adata_temp.X = scaler.transform(adata_temp.X)
#     else:
#         sc.pp.scale(adata_temp, zero_center=True, max_value=None, copy=False, layer=None, obsm=None)
    if (train_x_partition != 'X') & (train_x_partition in adata_temp.obsm.keys()):
        print('train partition is not in OBSM, defaulting to PCA')
        # Now compute PCA
        sc.pp.pca(adata_temp, n_comps=100, use_highly_variable=True, svd_solver='arpack')
        sc.pl.pca_variance_ratio(adata_temp, log=True,n_pcs=100)
        
        # Batch correction options
        # The script will test later which Harmony values we should use 
        if(batch_correction == "Harmony"):
            print("Commencing harmony")
            adata_temp.obs['lr_batch'] = adata_temp.obs[batch_key]
            batch_var = "lr_batch"
            # Create hm subset
            adata_hm = adata_temp[:]
            # Set harmony variables
            data_mat = np.array(adata_hm.obsm["X_pca"])
            meta_data = adata_hm.obs
            vars_use = [batch_var]
            # Run Harmony
            ho = hm.run_harmony(data_mat, meta_data, vars_use,theta=theta)
            res = (pd.DataFrame(ho.Z_corr)).T
            res.columns = ['X{}'.format(i + 1) for i in range(res.shape[1])]
            # Insert coordinates back into object
            adata_hm.obsm["X_pca_back"]= adata_hm.obsm["X_pca"][:]
            adata_hm.obsm["X_pca"] = np.array(res)
            # Run neighbours
            #sc.pp.neighbors(adata_hm, n_neighbors=15, n_pcs=50)
            adata_temp = adata_hm[:]
            del adata_hm
        elif(batch_correction == "BBKNN"):
            print("Commencing BBKNN")
            sc.external.pp.bbknn(adata_temp, batch_key=batch_var, approx=True, metric='angular', copy=False, n_pcs=50, trim=None, n_trees=10, use_faiss=True, set_op_mix_ratio=1.0, local_connectivity=15) 
        print("adata1 and adata2 are now combined and preprocessed in 'adata' obj - success!")


    # train model
#    train_x = adata_temp.X
    #train_label = adata_temp.obs[feat_use]
    print('proceeding to train model')
    model = LR_train(adata_temp, train_x = train_x_partition, train_label=feat_use, penalty=penalty, sparcity=sparcity,max_iter=max_iter,l1_ratio = l1_ratio,tune_hyper_params = tune_hyper_params)
    model.features = list(adata_temp.var.index)
    return model

def regression_results(y_true, y_pred):
    # Regression metrics
    explained_variance=metrics.explained_variance_score(y_true, y_pred)
    mean_absolute_error=metrics.mean_absolute_error(y_true, y_pred) 
    mse=metrics.mean_squared_error(y_true, y_pred) 
    mean_squared_log_error=metrics.mean_squared_log_error(y_true, y_pred)
    median_absolute_error=metrics.median_absolute_error(y_true, y_pred)
    r2=metrics.r2_score(y_true, y_pred)
    print('explained_variance: ', round(explained_variance,4))    
    print('mean_squared_log_error: ', round(mean_squared_log_error,4))
    print('r2: ', round(r2,4))
    print('MAE: ', round(mean_absolute_error,4))
    print('MSE: ', round(mse,4))
    print('RMSE: ', round(np.sqrt(mse),4))

