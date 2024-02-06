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
    else:
        adata_temp = adata[:]
    # model projections
    print('Starting reference projection!')
    if train_x_partition == 'X':
        train_x = adata_temp.X
        pred_out = pd.DataFrame(model_lr.predict(train_x),columns = ['predicted'],index = list(adata.obs.index))
        proba =  pd.DataFrame(model_lr.predict_proba(train_x),columns = model_lr.classes_,index = list(adata.obs.index))
        pred_out = pred_out.join(proba)
        
    elif train_x_partition in list(adata_temp.obsm.keys()): 
        print('{low_dim: this partition modality is still under development!}')
        train_x = adata_temp.obsm[train_x_partition][:]
        pred_out = pd.DataFrame(model_lr.predict(train_x),columns = ['predicted'],index = list(adata.obs.index))
        proba =  pd.DataFrame(model_lr.predict_proba(train_x),columns = model_lr.classes_,index = list(adata.obs.index))
        pred_out = pred_out.join(proba)

    
    else:
        print('{this partition modality is still under development!}')
        print('Warning! No obsm partition detected! defaulting to PCA, if this is not self-projection, do not use the results!')
        if not 'X_pca' in adata.obsm.keys():
            print('performing highly variable gene selection')
            sc.pp.highly_variable_genes(adata_temp, batch_key = batch_key, subset=False)
            sc.pp.pca(adata_temp, n_comps=100, use_highly_variable=True, svd_solver='arpack')
            sc.pl.pca_variance_ratio(adata_temp, log=True,n_pcs=100)
            train_x = adata_temp.obsm['X_pca'][:]
            pred_out = pd.DataFrame(model_lr.predict(train_x),columns = ['predicted'],index = list(adata.obs.index))
            proba =  pd.DataFrame(model_lr.predict_proba(train_x),columns = model_lr.classes_,index = list(adata.obs.index))
            pred_out = pred_out.join(proba)

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
def model_class_feature_plots(top_loadings, classes, comps, p_lim, max_len,title):
    import matplotlib.pyplot as plt
    for class_temp in classes:
        class_lw = class_temp
        long_format = top_loadings
        df_loadings = long_format[long_format['class'].isin([class_lw])]
        plt.hist(df_loadings[comps])
        
        if len(((df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]).unique())) < max_len:
            for i in ((df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]).unique()):
                plt.axvline(x=i,color='red')
        elif len(((df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]).unique())) > max_len:
            for i in ((df_loadings[comps].nlargest(max_len)).unique()):
                plt.axvline(x=i,color='red')
        med = np.median(df_loadings[comps])
        plt.axvline(x=med,color='blue')
        plt.xlabel('feature_importance', fontsize=12)
        plt.title(title)
        #plt.axvline(x=med,color='pinkp_lim
        df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]
        print(len(df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]))
        
        
        #Plot feature ranking
        if len(((df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]).unique())) < max_len:
            plot_loading = pd.DataFrame(pd.DataFrame(df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]).iloc[:,0].sort_values(ascending=False))
        elif len(((df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]).unique())) > max_len:
            plot_loading = pd.DataFrame(pd.DataFrame(df_loadings[comps].nlargest(max_len)).iloc[:,0].sort_values(ascending=False))
        
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
    plt.figure(figsize = (20,20))
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
    table.set_fontsize(5)
    g = sn.heatmap(df_cm, annot=False, annot_kws={"size": 16},cmap=pal) # font size
    print(metrics.classification_report(train_label, model.predict(train_x), digits=2))
    plt.show()
    return metric
    
    
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
def LR_train(adata, train_x, train_label, penalty='elasticnet', sparcity=0.2,max_iter=200,l1_ratio =0.2,tune_hyper_params =False,n_splits=5, n_repeats=3,l1_grid = [0.05,0.2,0.5], c_grid = [0.05,0.2,0.4],sketch_obsm =None):
    if tune_hyper_params == True:
        train_labels = train_label
        results,adata_tuned = tune_lr_model(adata, train_x_partition = train_x, random_state = 42,  train_labels = train_labels, n_splits=n_splits, n_repeats=n_repeats,l1_grid = l1_grid, c_grid = c_grid,sketch_obsm = sketch_obsm)
        print('hyper_params tuned')
        sparcity = results.best_params_['C']
        l1_ratio = results.best_params_['l1_ratio']
        
    if not sketch_obsm == None:
        #sketch data
        try:
            adata = sketch_data(adata, train_x_partition = train_x, random_state = 42,  train_labels = train_label,sketch_obsm = sketch_obsm)
        except:
            print()

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
        model = lr.fit(train_x, train_label)
        model.features = np.array(adata.var.index)
    elif train_x in adata.obsm.keys():
        print('train with obsm')
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
        model.features = list(pd.DataFrame(train_x).columns)
        #model.features = list(range()) list(pd.DataFrame(adata.obsm[train_x]).columns)
    return model

def tune_lr_model(adata, train_x_partition = 'X', random_state = 42,  train_labels = None, n_splits=5, n_repeats=3,l1_grid = [0.05,0.2,0.5], c_grid = [0.05,0.2,0.4],sketch_obsm = None):
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
        tune_train_x = adata.obsm[train_x_partition][:]
        lvg = bless.bless(tune_train_x, RBF(length_scale=20), lam_final = 2, qbar = 10, random_state = r, H = 10, force_cpu=True)
        adata_tuning = adata[lvg.idx]
        if not train_labels == None:
            tune_train_label = adata_tuning.obs[train_labels]
            orig_labels = adata.obs[train_labels]
            # check if all labels are preserved
            non_int = len(list(set(tune_train_label)^set(orig_labels)))
            if ((non_int)) > 0 :
                print('{} labels lost in sketch, attempting re-sketch'.format(non_int))
                counter = 1
                while counter < 5:
                    q_bar = 10 + (5*counter)
                    lvg = bless.bless(tune_train_x, RBF(length_scale=20), lam_final = 2, qbar = q_bar, random_state = r, H = 10, force_cpu=True)
                    adata_tuning = adata[lvg.idx]
                    tune_train_label = adata_tuning.obs[train_labels]
                    non_int = len(list(set(tune_train_label)^set(orig_labels)))
                    if non_int > 0:
                        counter = counter + 1
                    else:
                        counter = 5
        tune_train_x = adata_tuning.obsm[sketch_obsm][:]
        print('Sketched data is {} long'.format(len(adata_tuning.obs)))

    elif sketch_obsm in adata.obsm.keys():
        sketch_obsm_id = adata.obsm[sketch_obsm][:]
        lvg = bless.bless(sketch_obsm_id, RBF(length_scale=20), lam_final = 2, qbar = 10, random_state = r, H = 10, force_cpu=True)
        adata_tuning = adata[lvg.idx]
        if not train_labels == None:
            tune_train_label = adata_tuning.obs[train_labels]
            orig_labels = adata.obs[train_labels]
            # check if all labels are preserved
            non_int = len(list(set(tune_train_label)^set(orig_labels)))
            if ((non_int)) > 0 :
                print('{} labels lost in sketch, attempting re-sketch'.format(non_int))
                counter = 1
                while counter < 5:
                    q_bar = 10 + (5*counter)
                    lvg = bless.bless(tune_train_x, RBF(length_scale=20), lam_final = 2, qbar = q_bar, random_state = r, H = 10, force_cpu=True)
                    adata_tuning = adata[lvg.idx]
                    tune_train_label = adata_tuning.obs[train_labels]
                    non_int = len(list(set(tune_train_label)^set(orig_labels)))
                    if non_int > 0:
                        counter = counter + 1
                    else:
                        counter = 5          
        print('Sketched data is {} long'.format(len(adata_tuning.obs)))
        tune_train_x = adata_tuning.obsm[sketch_obsm][:]
    #     try:
    #         import cupy
    #         lvg_2 = bless(adata.obsm[train_x_partition], RBF(length_scale=10), 10, 10, r, 10, force_cpu=False)
    #     except ImportError:
    #         print("cupy not found, defaulting to numpy")
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
    model = LogisticRegression(penalty = penalty, C = sparcity, max_iter =  500, n_jobs=4)
    if (penalty == "l1"):
        model = LogisticRegression(penalty = penalty, C = sparcity, max_iter =  500, dual = True, solver = 'liblinear',multi_class = 'ovr', n_jobs=4 ) # one-vs-rest
    if (penalty == "elasticnet"):
        model = LogisticRegression(penalty = penalty, C = sparcity, max_iter =  500, dual=False,solver = 'saga',l1_ratio=l1_ratio,multi_class = 'ovr', n_jobs=4) # use multinomial class if probabilities are descrete
        grid['l1_ratio'] = l1_grid
    grid['C'] = c_grid
    # define search
    search = GridSearchCV(model, grid, scoring='neg_mean_absolute_error', cv=cv, n_jobs=-1)
    # perform the search
    results = search.fit(X, y)
    # summarize
    print('MAE: %.3f' % results.best_score_)
    print('Config: %s' % results.best_params_)
    return results , adata_tuning

def prep_training_data(adata_temp,feat_use,batch_key, model_key, batch_correction=False, var_length = 7500,penalty='elasticnet',sparcity=0.2,max_iter = 200,l1_ratio = 0.1,partial_scale=True,train_x_partition ='X',theta = 3,tune_hyper_params=False,sketch_obsm = None ):
    model_name = model_key + '_lr_model'
    print('performing highly variable gene selection')
    #sc.pp.highly_variable_genes(adata_temp, batch_key = batch_key, subset=False)
#     #temp inclusion
#     sc.pp.pca(adata_temp, n_comps=100, use_highly_variable=True, svd_solver='arpack')
#     sc.pl.pca_variance_ratio(adata_temp, log=True,n_pcs=100)
        
    adata_temp = adata_temp#subset_top_hvgs(adata_temp,var_length)
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

    if (train_x_partition != 'X') & (train_x_partition not in adata_temp.obsm.keys()):
        print('train partition is not in OBSM, defaulting to PCA')
        # Now compute PCA
        sc.pp.pca(adata_temp, n_comps=100, use_highly_variable=True, svd_solver='arpack')
        sc.pl.pca_variance_ratio(adata_temp, log=True,n_pcs=100)
        
        # Batch correction options
        # The script will test later which Harmony values we should use 
        if(batch_correction == "Harmony"):
            print("Commencing harmony")
            if len(batch_key) == 1:
                adata_temp.obs['lr_batch'] = adata_temp.obs[batch_key]
                batch_var = "lr_batch"
            else:
                batch_var = batch_key
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

    # train model
#    train_x = adata_temp.X
    #train_label = adata_temp.obs[feat_use]
    print('proceeding to train model')
    model = LR_train(adata_temp, train_x = train_x_partition, train_label=feat_use, penalty=penalty, sparcity=sparcity,max_iter=max_iter,l1_ratio = l1_ratio,tune_hyper_params = tune_hyper_params,sketch_obsm = sketch_obsm)
#     model.features = list(adata_temp.var.index)
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
    
def sketch_data(adata, train_x_partition = 'X', sketch_obsm = None, random_state = 42,  train_labels = None):
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
        tune_train_x = adata.obsm[train_x_partition][:]
    elif sketch_obsm in adata.obsm.keys():
        tune_train_x = adata.obsm[sketch_obsm][:]
    else:
        print('No obsm partition detected! defaulting to PCA')
        if not 'X_pca' in adata.obsm.keys():
            print('performing highly variable gene selection')
            sc.pp.highly_variable_genes(adata_temp, batch_key = batch_key, subset=False)
            sc.pp.pca(adata_temp, n_comps=100, use_highly_variable=True, svd_solver='arpack')
            sc.pl.pca_variance_ratio(adata_temp, log=True,n_pcs=100)
        tune_train_x = adata.obsm['X_pca'][:]
    lvg = bless.bless(tune_train_x, RBF(length_scale=20), lam_final = 2, qbar = 2, random_state = r, H = 10, force_cpu=True)
    adata_tuning = adata[lvg.idx]
    print('sketched partition is {}, original is {}'.format(len(lvg.idx)),len(adata.obs))
    return adata_tuning


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

model = load_models(models,model_key)
model_lr =  model
# Estimate top model features for class descrimination
feature_importance = estimate_important_features(model_lr, 100)
mat = feature_importance.euler_pow_mat
top_loadings = feature_importance.to_n_features_long

org_key = 'SK'
model_key = 'overall_fsk_adt'
out_dir = './A2_V1_sk_model_pan_organ_heatmaps'
train_x_partition = 'X_scvi'
feat_use = 'age'
# if low_dim, else all below == False
use_varm = '/nfs/team205/ig7/projects/fetal_skin/3_160523_probabillistic_projection_organoid_adt_fetl/A2_V2_ldvae_models/v3_ldvae_obsm_weights.csv'


# Extra layer for low_dim models
if use_varm != False:
    print('Low_dim model detected, proceeding to translate linear components into feature sets')
model_varm = pd.read_csv(use_varm,index_col = 0)
if len(sorted(list((top_loadings['feature'].unique())))) == len(list(model_varm.columns)):
    print("length of model and varm features are equal, proceeding")
    feature_set = dict(zip(sorted(list((top_loadings['feature'].unique()))),list(model_varm.columns)))
else: 
    print('Lengths of model and varm features are unequal! proceeding but cannot be sure that ordering if preserved!')
    model_varm = model_varm.iloc[0:len(sorted(list((top_loadings['feature'].unique()))))]
    feature_set = dict(zip(sorted(list((top_loadings['feature'].unique()))),list(model_varm.columns)))
    
# Top predictive features for each component
varm_melt = pd.melt(model_varm.reset_index(), id_vars='index')
varm_melt['pvals']= np.NAN
for variable in varm_melt['variable'].unique():
    varm_loadings = varm_melt[varm_melt['variable'].isin([variable])]
    comps = 'value'
    U = np.mean(varm_loadings[comps])
    std = np.std(varm_loadings[comps])
    med =  np.median(varm_loadings[comps])
    mad = np.median(np.absolute(varm_loadings[comps] - np.median(varm_loadings[comps])))
    # Survival function scaled by 1.4826 of MAD (approx norm)
    pvals = scipy.stats.norm.sf(varm_loadings[comps], loc=med, scale=1.4826*mad)
    varm_melt.loc[varm_melt['variable'].isin([variable]),'pvals'] = pvals
varm_sig = varm_melt[varm_melt['pvals']<0.05]

# varm sig is the reference for weighting, report top 10 weighted features per class

top_loadings_lowdim = pd.DataFrame(columns =['class','feature', 'weighted_impact', 'e^coef_pval', 'e^coef', 'is_significant_sf'])
top_loadings_lw = top_loadings.groupby(['class']).head(10)
top_loadings_lw['feature'] = top_loadings_lw['feature'].map(feature_set)

for classes in top_loadings_lw['class'].unique():
    print(classes)
    # linearise in next update!
    for feature in top_loadings_lw.loc[top_loadings_lw['class'].isin([classes]) ,['feature','e^coef']].values:
        # multiply value by impact! 
        temp_varm_sig = varm_sig[varm_sig['variable'].isin([feature[0]])]
        temp_varm_sig['weighted_impact'] = temp_varm_sig['value'] * feature[1]
        temp_varm_sig = temp_varm_sig[['index','weighted_impact']]
        temp_varm_sig.columns = ['feature','weighted_impact']
        temp_varm_sig['class'] = classes
        temp_varm_sig['e^coef_pval'] = top_loadings_lw.loc[(top_loadings_lw['class'].isin([classes])) & (top_loadings_lw['feature'].isin([feature[0]])),'e^coef_pval'].values[0]
        temp_varm_sig['e^coef'] = top_loadings_lw.loc[(top_loadings_lw['class'].isin([classes])) & (top_loadings_lw['feature'].isin([feature[0]])),'e^coef'].values[0]
        temp_varm_sig['is_significant_sf'] = top_loadings_lw.loc[(top_loadings_lw['class'].isin([classes])) & (top_loadings_lw['feature'].isin([feature[0]])),'is_significant_sf'].values[0]
        top_loadings_lowdim = pd.concat([top_loadings_lowdim,temp_varm_sig],ignore_index=True)
top_loadings_lowdim = top_loadings_lowdim.sort_values(['weighted_impact'],ascending=False).groupby('class').head(100)

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

top_loadings_lowdim = compute_weighted_impact(varm_file = '/nfs/team205/ig7/projects/fetal_skin/3_160523_probabillistic_projection_organoid_adt_fetl/A2_V2_ldvae_models/v3_ldvae_obsm_weights.csv',top_loadings =  top_loadings, threshold=0.05)

def model_class_feature_plots(top_loadings, classes, comps, p_lim, max_len,title):
    import matplotlib.pyplot as plt
    for class_temp in classes:
        class_lw = class_temp
        long_format = top_loadings
        df_loadings = long_format[long_format['class'].isin([class_lw])]
        plt.hist(df_loadings[comps])
        
        if len(((df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]).unique())) < max_len:
            for i in ((df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]).unique()):
                plt.axvline(x=i,color='red')
        elif len(((df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]).unique())) > max_len:
            for i in ((df_loadings[comps].nlargest(max_len)).unique()):
                plt.axvline(x=i,color='red')
        med = np.median(df_loadings[comps])
        plt.axvline(x=med,color='blue')
        plt.xlabel('feature_importance', fontsize=12)
        plt.title(title)
        #plt.axvline(x=med,color='pinkp_lim
        df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]
        print(len(df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]))
        
        
        #Plot feature ranking
        if len(((df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]).unique())) < max_len:
            plot_loading = pd.DataFrame(pd.DataFrame(df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]).iloc[:,0].sort_values(ascending=False))
        elif len(((df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]).unique())) > max_len:
            plot_loading = pd.DataFrame(pd.DataFrame(df_loadings[comps].nlargest(max_len)).iloc[:,0].sort_values(ascending=False))
        
        table = plt.table(cellText=plot_loading.values,colWidths = [1]*len(plot_loading.columns),
        rowLabels= list(df_loadings['feature'][df_loadings.index.isin(plot_loading.index)].reindex(plot_loading.index)), #plot_loading.index,
        colLabels=plot_loading.columns,
        cellLoc = 'center', rowLoc = 'center',
        loc='right', bbox=[1.4, -0.05, 0.5,1])
        table.scale(1, 2)
        table.set_fontsize(10)

for class_lin in top_loadings_lowdim['class'].unique():
    model_class_feature_plots(top_loadings_lowdim, [class_lin], 'weighted_impact','e^coef',max_len= 20,title = lineage)
    plt.show()
    top_loadings_u =  top_loadings_lowdim_dic[lineage].head(400)
    top_loadings_u['gene'] = top_loadings_u['feature']
    #Convert to long-format
    # DE = convert_scanpy_degs(sc_de)
    Example_DE = top_loadings_u
    Example_DE.head(3)
    pre_ranked = True
    if pre_ranked == False:
        #Define gene list
        glist = as_gene_list(Example_DE,"gene")
        #Get possible enrichr queries
        names = gp.get_library_name('Human') # default: Human "e.g GO_Biological_Process_2018"
        #Define the databases to query
        database=['GO_Biological_Process_2018'] #'KEGG_2016',
        #No ranking Query databases
        enr = enrich_no_rank(glist,database)
        #Explore enr class obj
        #print(dir(enr))
        #print(enr.__dict__)
        # simple plotting function
        #from gseapy.plot import barplot, dotplot
        # to save your figure, make sure that ``ofname`` is not None
        #barplot(enr.res2d,title='KEGG_2013',)
        enr.results.head(5)
    else:
        #Get possible enrichr queries
        names = gp.get_library_name('Human') # default: Human "e.g GO_Biological_Process_2018"
        #Define the databases to query
        database='GO_Biological_Process_2021' #GO_Biological_Process_2023
        enr = pre_ranked_enr(Example_DE,["gene","weighted_impact"],permutation_num = 1000,database = database,cutoff=0.25,min_s=5)
        print("Normalised enrichment score ranges, for continous phenotype tests, a positive value indicates correlation with phenotype profile or top of the ranked list, negative values indicate inverse corelation with profile or corelation with the bottom of the list")
        # print(max(enr.res2d.nes))
        # print(min(enr.res2d.nes))
        print(enr.res2d.shape)
        enr.res2d.head(3)
        #For ranked sets of individual comparisons, it may be wise to determine an ES cut-off

        #enr.res2d.sort_values(by=['nes'], ascending=False)
    ## easy way
    terms = enr.res2d.Term
    # axs = enr.plot(terms=terms[1]) # v1.0.5
    axs = enr.plot(terms=terms[1:10],
                       legend_kws={'loc': (1.2, 0)}, # set the legend loc
                       show_ranking=True, # whether to show the second yaxis
                       figsize=(3,4))
    plt.show()

import gseapy as gp
from gseapy.plot import barplot, dotplot

def analyze_and_plot_feat_gsea(top_loadings_lowdim, class_name, max_len=20, pre_ranked=True, database='GO_Biological_Process_2021', cutoff=0.25, min_s=5):
    """
    Analyzes and plots the top loadings of a low-dimensional model.

    Parameters:
    top_loadings_lowdim (pd.DataFrame): A dataframe containing the top loadings for each class.
    class_name (str): The name of the class to analyze and plot.
    max_len (int): The maximum number of features to plot.
    pre_ranked (bool): Whether the data is pre-ranked.
    database (str): The name of the database to query for enrichment analysis.
    cutoff (float): The cutoff value for enrichment analysis.
    min_s (int): The minimum number of genes in a category to consider for enrichment analysis.

    Returns:
    None
    """
    # Filter the top loadings for the given class
    top_loadings_u = top_loadings_lowdim[top_loadings_lowdim['class'] == class_name].head(max_len)
    top_loadings_u['gene'] = top_loadings_u['feature']

    # Perform enrichment analysis
    if not pre_ranked == True:
        glist = as_gene_list(top_loadings_u, "gene")
        enr = enrich_no_rank(glist, [database])
    else:
        enr = pre_ranked_enr(top_loadings_u, ["gene", "weighted_impact"], permutation_num=1000, database=database, cutoff=cutoff, min_s=min_s)

    # Print the enrichment score range
    print("Normalised enrichment score ranges, for continuous phenotype tests, a positive value indicates correlation with phenotype profile or top of the ranked list, negative values indicate inverse correlation with profile or correlation with the bottom of the list")
    print(enr.res2d.shape)
    # Plot the enrichment results
    terms = enr.res2d.Term
    axs = enr.plot(terms=terms[1:10], legend_kws={'loc': (1.2, 0)}, show_ranking=True, figsize=(3, 4))
    plt.show()
    
for class_lin in top_loadings_lowdim['class'].unique():
    model_class_feature_plots(top_loadings_lowdim, [class_lin], 'weighted_impact','e^coef',max_len= 20,title = lineage)
    analyze_and_plot(top_loadings_lowdim,class_lin, max_len=20, pre_ranked=True, database='GO_Biological_Process_2021', cutoff=0.25, min_s=5)

top_loadings_lowdim_dic = {}
for lineage in lineages:
    model = load_models(models,lineage)
    model_lr =  model
    # Estimate top model features for class descrimination
    feature_importance = estimate_important_features(model_lr, 100)
    mat = feature_importance.euler_pow_mat
    top_loadings = feature_importance.to_n_features_long

    # Extra layer for low_dim models
    if use_varm != False:
        print('Low_dim model detected, proceeding to translate linear components into feature sets')
    model_varm = pd.read_csv(use_varm,index_col = 0)
    if len(sorted(list((top_loadings['feature'].unique())))) == len(list(model_varm.columns)):
        print("length of model and varm features are equal, proceeding")
        feature_set = dict(zip(sorted(list((top_loadings['feature'].unique()))),list(model_varm.columns)))
    else: 
        print('Lengths of model and varm features are unequal! proceeding but cannot be sure that ordering if preserved!')
        model_varm = model_varm.iloc[0:len(sorted(list((top_loadings['feature'].unique()))))]
        feature_set = dict(zip(sorted(list((top_loadings['feature'].unique()))),list(model_varm.columns)))

    # Top predictive features for each component
    varm_melt = pd.melt(model_varm.reset_index(), id_vars='index')
    varm_melt['pvals']= np.NAN
    for variable in varm_melt['variable'].unique():
        varm_loadings = varm_melt[varm_melt['variable'].isin([variable])]
        comps = 'value'
        U = np.mean(varm_loadings[comps])
        std = np.std(varm_loadings[comps])
        med =  np.median(varm_loadings[comps])
        mad = np.median(np.absolute(varm_loadings[comps] - np.median(varm_loadings[comps])))
        # Survival function scaled by 1.4826 of MAD (approx norm)
        pvals = scipy.stats.norm.sf(varm_loadings[comps], loc=med, scale=1.4826*mad)
        varm_melt.loc[varm_melt['variable'].isin([variable]),'pvals'] = pvals
    varm_sig = varm_melt[varm_melt['pvals']<0.05]

    # varm sig is the reference for weighting, report top 10 weighted features per class

    top_loadings_lowdim = pd.DataFrame(columns =['class','feature', 'weighted_impact', 'e^coef_pval', 'e^coef',
           'is_significant_sf'])
    top_loadings_lw = top_loadings.groupby(['class']).head(10)
    top_loadings_lw['feature'] = top_loadings_lw['feature'].map(feature_set)
    for classes in top_loadings_lw['class'].unique():
        print(classes)
        # linearise in next update!
        for feature in top_loadings_lw.loc[top_loadings_lw['class'].isin([classes]) ,['feature','e^coef']].values:
            # multiply value by impact! 
            temp_varm_sig = varm_sig[varm_sig['variable'].isin([feature[0]])]
            temp_varm_sig['weighted_impact'] = temp_varm_sig['value'] * feature[1]
            temp_varm_sig = temp_varm_sig[['index','weighted_impact']]
            temp_varm_sig.columns = ['feature','weighted_impact']
            temp_varm_sig['class'] = classes
            temp_varm_sig['e^coef_pval'] = top_loadings_lw.loc[(top_loadings_lw['class'].isin([classes])) & (top_loadings_lw['feature'].isin([feature[0]])),'e^coef_pval'].values[0]
            temp_varm_sig['e^coef'] = top_loadings_lw.loc[(top_loadings_lw['class'].isin([classes])) & (top_loadings_lw['feature'].isin([feature[0]])),'e^coef'].values[0]
            temp_varm_sig['is_significant_sf'] = top_loadings_lw.loc[(top_loadings_lw['class'].isin([classes])) & (top_loadings_lw['feature'].isin([feature[0]])),'is_significant_sf'].values[0]
            top_loadings_lowdim = pd.concat([top_loadings_lowdim,temp_varm_sig],ignore_index=True)
    top_loadings_lowdim = top_loadings_lowdim.sort_values(['weighted_impact'],ascending=False).groupby('class').head(100)

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
    top_loadings_lowdim_dic[lineage] = top_loadings_lowdim

top_loadings_lowdim_dic.keys()

top_loadings_lowdim_dic['Mural']['class'].unique()

top_loadings_lowdim_dic[lineage][top_loadings_lowdim_dic[lineage]['class'].str.contains('999')]

top_loadings_lowdim_dic[lineage]

top_loadings_lowdim_dic[lineage]['weighted_impact'].nlargest(10)

def model_class_feature_plots(top_loadings, classes, comps, p_lim, max_len,title):
    import matplotlib.pyplot as plt
    for class_temp in classes:
        class_lw = class_temp
        long_format = top_loadings
        df_loadings = long_format[long_format['class'].isin([class_lw])]
        plt.hist(df_loadings[comps])
        
        if len(((df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]).unique())) < max_len:
            for i in ((df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]).unique()):
                plt.axvline(x=i,color='red')
        elif len(((df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]).unique())) > max_len:
            for i in ((df_loadings[comps].nlargest(max_len)).unique()):
                plt.axvline(x=i,color='red')
        med = np.median(df_loadings[comps])
        plt.axvline(x=med,color='blue')
        plt.xlabel('feature_importance', fontsize=12)
        plt.title(title)
        #plt.axvline(x=med,color='pinkp_lim
        df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]
        print(len(df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]))
        
        
        #Plot feature ranking
        if len(((df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]).unique())) < max_len:
            plot_loading = pd.DataFrame(pd.DataFrame(df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]).iloc[:,0].sort_values(ascending=False))
        elif len(((df_loadings[comps][df_loadings[str(p_lim) +'_pval']<0.05]).unique())) > max_len:
            plot_loading = pd.DataFrame(pd.DataFrame(df_loadings[comps].nlargest(max_len)).iloc[:,0].sort_values(ascending=False))
        
        table = plt.table(cellText=plot_loading.values,colWidths = [1]*len(plot_loading.columns),
        rowLabels= list(df_loadings['feature'][df_loadings.index.isin(plot_loading.index)].reindex(plot_loading.index)), #plot_loading.index,
        colLabels=plot_loading.columns,
        cellLoc = 'center', rowLoc = 'center',
        loc='right', bbox=[1.4, -0.05, 0.5,1])
        table.scale(1, 2)
        table.set_fontsize(10)

%matplotlib inline
%config InlineBackend.figure_format='retina' # mac
%load_ext autoreload
%autoreload 2
importlibs()

for lineage in lineages:
    class_lin = top_loadings_lowdim_dic[lineage][top_loadings_lowdim_dic[lineage]['class'].str.contains('999')]['class'].unique()[0]
    model_class_feature_plots(top_loadings_lowdim_dic[lineage], [class_lin], 'weighted_impact','e^coef',max_len= 20,title = lineage)
    plt.show()
    
    top_loadings_u =  top_loadings_lowdim_dic[lineage].head(400)
    top_loadings_u['gene'] = top_loadings_u['feature']

    #sc_de = pd.read_csv("/nfs/team205/ig7/resources/scripts_dont_modify/gene_ontology_summarisation/resources/DEGS_Macrophage.csv")

    #Example convert ENSG to HGNC
    #convert_hgnc(['ENSG00000148795', 'ENSG00000165359', 'ENSG00000150676'])

    #Convert to long-format
    # DE = convert_scanpy_degs(sc_de)
    Example_DE = top_loadings_u
    Example_DE.head(3)
    pre_ranked = True
    if pre_ranked == False:
        #Define gene list
        glist = as_gene_list(Example_DE,"gene")
        #Get possible enrichr queries
        names = gp.get_library_name('Human') # default: Human "e.g GO_Biological_Process_2018"
        #Define the databases to query
        database=['GO_Biological_Process_2018'] #'KEGG_2016',
        #No ranking Query databases
        enr = enrich_no_rank(glist,database)
        #Explore enr class obj
        #print(dir(enr))
        #print(enr.__dict__)
        # simple plotting function
        #from gseapy.plot import barplot, dotplot
        # to save your figure, make sure that ``ofname`` is not None
        #barplot(enr.res2d,title='KEGG_2013',)
        enr.results.head(5)
    else:
        #Get possible enrichr queries
        names = gp.get_library_name('Human') # default: Human "e.g GO_Biological_Process_2018"
        #Define the databases to query
        database='GO_Biological_Process_2021' #GO_Biological_Process_2023

        enr = pre_ranked_enr(Example_DE,["gene","weighted_impact"],permutation_num = 1000,database = database,cutoff=0.25,min_s=5)

        print("Normalised enrichment score ranges, for continous phenotype tests, a positive value indicates correlation with phenotype profile or top of the ranked list, negative values indicate inverse corelation with profile or corelation with the bottom of the list")
        # print(max(enr.res2d.nes))
        # print(min(enr.res2d.nes))
        print(enr.res2d.shape)
        enr.res2d.head(3)
        #For ranked sets of individual comparisons, it may be wise to determine an ES cut-off

        #enr.res2d.sort_values(by=['nes'], ascending=False)

    ## easy way
    terms = enr.res2d.Term
    # axs = enr.plot(terms=terms[1]) # v1.0.5
    axs = enr.plot(terms=terms[1:10],
                       legend_kws={'loc': (1.2, 0)}, # set the legend loc
                       show_ranking=True, # whether to show the second yaxis
                       figsize=(3,4)
                      )
    # or use this to have more control on the plot
    # from gseapy import gseaplot2
    # terms = pre_res.res2d.Term[1:5]
    # hits = [pre_res.results[t]['hits'] for t in terms]
    # runes = [pre_res.results[t]['RES'] for t in terms]
    # fig = gseaplot2(terms=terms, ress=runes, hits=hits,
    #               rank_metric=gs_res.ranking,
    #               legend_kws={'loc': (1.2, 0)}, # set the legend loc
    #               figsize=(4,5)) # rank_metric=pre_res.ranking
    # from gseapy import dotplot
    # # to save your figure, make sure that ``ofname`` is not None
    # ax = dotplot(enr.res2d,
    #              column="FDR q-val",
    #              title='GO_Biological_Process_2018',
    #              cmap=plt.cm.viridis,
    #              size=3, # adjust dot size
    #              figsize=(4,5), cutoff=0.5, show_ring=False)
    plt.show()

def model_class_feature_plots_negative_impact(top_loadings, classes, comps, p_lim, max_len,title):
    import matplotlib.pyplot as plt
    for class_temp in classes:
        class_lw = class_temp
        long_format = top_loadings
        df_loadings = long_format[long_format['class'].isin([class_lw])]
        plt.hist(df_loadings[comps])
        
#         if len(((df_loadings[comps][df_loadings[str(p_lim) +'_pval']>0.05]).unique())) < max_len:
#             for i in ((df_loadings[comps][df_loadings[str(p_lim) +'_pval']>0.05]).unique()):
#                 plt.axvline(x=i,color='red')
#        elif len(((df_loadings[comps][df_loadings[str(p_lim) +'_pval']>0.05]).unique())) > max_len:
        for i in ((df_loadings[comps].nsmallest(max_len)).unique()):
            plt.axvline(x=i,color='red')
            
        med = np.median(df_loadings[comps])
        plt.axvline(x=med,color='blue')
        plt.xlabel('feature_importance', fontsize=12)
        plt.title(title)
        #plt.axvline(x=med,color='pinkp_lim
        df_loadings[comps][df_loadings[str(p_lim) +'_pval']>0.05]
        print(len(df_loadings[comps][df_loadings[str(p_lim) +'_pval']>0.05]))
        
        
        #Plot feature ranking
#         if len(((df_loadings[comps][df_loadings[str(p_lim) +'_pval']>0.05]).unique())) < max_len:
#             plot_loading = pd.DataFrame(pd.DataFrame(df_loadings[comps][df_loadings[str(p_lim) +'_pval']>0.05]).iloc[:,0].sort_values(ascending=False))
#         elif len(((df_loadings[comps][df_loadings[str(p_lim) +'_pval']>0.05]).unique())) > max_len:
        plot_loading = pd.DataFrame(pd.DataFrame(df_loadings[comps].nsmallest(max_len)).iloc[:,0].sort_values(ascending=False))
        
        table = plt.table(cellText=plot_loading.values,colWidths = [1]*len(plot_loading.columns),
        rowLabels= list(df_loadings['feature'][df_loadings.index.isin(plot_loading.index)].reindex(plot_loading.index)), #plot_loading.index,
        colLabels=plot_loading.columns,
        cellLoc = 'center', rowLoc = 'center',
        loc='right', bbox=[1.4, -0.05, 0.5,1])
        table.scale(1, 2)
        table.set_fontsize(10)

for lineage in lineages:
    class_lin = top_loadings_lowdim_dic[lineage][top_loadings_lowdim_dic[lineage]['class'].str.contains('999')]['class'].unique()[0]
    model_class_feature_plots_negative_impact(top_loadings_lowdim_dic[lineage], [class_lin], 'weighted_impact','e^coef',max_len= 20,title = lineage)
    plt.show()
    
    top_loadings_u =  top_loadings_lowdim_dic[lineage].sort_values('weighted_impact', axis=0, ascending=True, inplace=False, kind='quicksort', na_position='last', ignore_index=False, key=None).head(400)
    top_loadings_u['weighted_impact'] = 1 - top_loadings_u['weighted_impact']
    top_loadings_u['gene'] = top_loadings_u['feature']

    #sc_de = pd.read_csv("/nfs/team205/ig7/resources/scripts_dont_modify/gene_ontology_summarisation/resources/DEGS_Macrophage.csv")

    #Example convert ENSG to HGNC
    #convert_hgnc(['ENSG00000148795', 'ENSG00000165359', 'ENSG00000150676'])

    #Convert to long-format
    # DE = convert_scanpy_degs(sc_de)
    Example_DE = top_loadings_u
    Example_DE.head(3)
    pre_ranked = True
    if pre_ranked == False:
        #Define gene list
        glist = as_gene_list(Example_DE,"gene")
        #Get possible enrichr queries
        names = gp.get_library_name('Human') # default: Human "e.g GO_Biological_Process_2018"
        #Define the databases to query
        database=['GO_Biological_Process_2018'] #'KEGG_2016',
        #No ranking Query databases
        enr = enrich_no_rank(glist,database)
        #Explore enr class obj
        #print(dir(enr))
        #print(enr.__dict__)
        # simple plotting function
        #from gseapy.plot import barplot, dotplot
        # to save your figure, make sure that ``ofname`` is not None
        #barplot(enr.res2d,title='KEGG_2013',)
        enr.results.head(5)
    else:
        #Get possible enrichr queries
        names = gp.get_library_name('Human') # default: Human "e.g GO_Biological_Process_2018"
        #Define the databases to query
        database='GO_Biological_Process_2021'

        enr = pre_ranked_enr(Example_DE,["gene","weighted_impact"],permutation_num = 1000,database = database,cutoff=0.25,min_s=5)

        print("Normalised enrichment score ranges, for continous phenotype tests, a positive value indicates correlation with phenotype profile or top of the ranked list, negative values indicate inverse corelation with profile or corelation with the bottom of the list")
        # print(max(enr.res2d.nes))
        # print(min(enr.res2d.nes))
        print(enr.res2d.shape)
        enr.res2d.head(3)
        #For ranked sets of individual comparisons, it may be wise to determine an ES cut-off

        #enr.res2d.sort_values(by=['nes'], ascending=False)

    ## easy way
    terms = enr.res2d.Term
    # axs = enr.plot(terms=terms[1]) # v1.0.5
    axs = enr.plot(terms=terms[1:15],
                       legend_kws={'loc': (1.2, 0)}, # set the legend loc
                       show_ranking=True, # whether to show the second yaxis
                       figsize=(3,4)
                      )
    # or use this to have more control on the plot
    # from gseapy import gseaplot2
    # terms = pre_res.res2d.Term[1:5]
    # hits = [pre_res.results[t]['hits'] for t in terms]
    # runes = [pre_res.results[t]['RES'] for t in terms]
    # fig = gseaplot2(terms=terms, ress=runes, hits=hits,
    #               rank_metric=gs_res.ranking,
    #               legend_kws={'loc': (1.2, 0)}, # set the legend loc
    #               figsize=(4,5)) # rank_metric=pre_res.ranking
    # from gseapy import dotplot
    # # to save your figure, make sure that ``ofname`` is not None
    # ax = dotplot(enr.res2d,
    #              column="FDR q-val",
    #              title='GO_Biological_Process_2018',
    #              cmap=plt.cm.viridis,
    #              size=3, # adjust dot size
    #              figsize=(4,5), cutoff=0.5, show_ring=False)
    plt.show()
    

top_loadings_lowdim

model_class_feature_plots(top_loadings_lowdim, ['999'], 'e^coef')
plt.show()