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

models = {
'pan_fetal':'/nfs/team205/ig7/resources/scripts_dont_modify/logit_regression_models/adifa_lr/celltypist_model.Pan_Fetal_Human.pkl',
'pan_fetal_wget':'https://celltypist.cog.sanger.ac.uk/models/Pan_Fetal_Suo/v2/Pan_Fetal_Human.pkl',
'adata_scvi':'/nfs/team205/ig7/mount/gdrive/g_cloud/projects/amniontic_fluid/scvi_low_dim_model.sav',
'adata_ldvae':'/nfs/team205/ig7/mount/gdrive/g_cloud/projects/amniontic_fluid/ldvae_low_dim_model.sav',
'adata_harmony':'/nfs/team205/ig7/work_backups/backup_210306/projects/amiotic_fluid/train_low_dim_model/organ_low_dim_model.sav',
'test_low_dim_ipsc_ys':'/nfs/team205/ig7/work_backups/backup_210306/projects/YS/YS_030522_notebooks/Integrating_HM_data_030522/YS_logit/lr_model.sav',
'YS_X':'/nfs/team205/ig7/work_backups/backup_210306/projects/YS/resources/YS_X_model_080922.sav',
'YS_X_V3':'/nfs/team205/ig7/work_backups/backup_210306/projects/YS/rebuttal_figs_010922/train_YS_full_X_model/YS_X_A2_V12_lvl3_ELASTICNET_YS.sav',
'SK_model':'/nfs/team205/ig7/resources/scripts_dont_modify/logit_regression_models/LR_app_format/hudaa_skin/for_hudaa_A1_V2',
'Hudaa_model_trained':'/nfs/team298/hg6/Fetal_skin/LR_15012023/train-all_model.pkl',
'low_dim_sk_model':'/nfs/team205/ig7/resources/scripts_dont_modify/logit_regression_models/LR_app_format/hudaa_skin/SK_model_lowdim',
'low_dim_ldvae_model_YS_cross_organ':'/nfs/team205/ig7/work_backups/backup_210306/projects/YS/YS_data/YS_panf_gonads_brain_eliv_combined_060922/ldVAE_model_projections/YS_v3_macs_mod_ldvae_panf_model',
'overall_fsk_adt':'/nfs/team205/ig7/projects/fetal_skin/3_160523_probabillistic_projection_organoid_adt_fetl/A2_V2_ldvae_models/ldvae_linear_v5_adt_fsk_time_model_adt_subsamp',
}

lineages = ['Mural','Kera','Schwann','Fibro','Vascular','Melan']
model_key = 'ldvae_linear_v5_adt_fsk_time_model'
path = '/nfs/team205/ig7/projects/fetal_skin/3_160523_probabillistic_projection_organoid_adt_fetl/A2_V2_ldvae_models//'
for lineage in lineages:
    lineage_key = lineage
    models[lineage] = path + lineage + '_' + model_key

adatas_dict = {
'Fetal_skin_raw': '/nfs/team298/hg6/Fetal_skin/data/FS_raw_sub.h5ad',
'vascular_organoid': '/nfs/team298/hg6/Fetal_skin/data/vasc_org_raw.h5ad',
'YS':'/nfs/team205/ig7/work_backups/backup_210306/projects/YS/YS_data/Submission_2_data/A2_V5_scvi_YS_integrated/A2_V5_scvi_YS_integrated_raw_qc_scr_umap.h5ad',
'YS_test':'/nfs/team205/ig7/resources/scripts_dont_modify/logit_regression_models/LR_app_format/ys_test_data.h5ad',
'YS_A2_V10_X_raw':'/nfs/team205/ig7/work_backups/backup_210306/projects/YS/YS_data/Submission_2_data/A2_V10_scvi_YS_integrated/A2_V10_raw_counts_full_no_obs.h5ad',
'YS_A2_V10_X':'/nfs/team205/ig7/work_backups/backup_210306/projects/YS/YS_data/Submission_2_data/A2_V10_scvi_YS_integrated/A2_V10_qc_raw.h5ad',
'pan_f_YS_A1_V10':'/nfs/team205/ig7/work_backups/backup_210306/projects/YS/YS_data/YS_panf_gonads_brain_eliv_combined_060922/A1_Vx_pan_organ_integrations/A1_V10_raw_scvi_YS_updated_panf_gonads_brain_build_donor_organ_corrected_031022.h5ad',
'pan_f_YS_A1_V10_high_var_ldvae':'/nfs/team205/ig7/work_backups/backup_210306/projects/YS/YS_data/YS_panf_gonads_brain_eliv_combined_060922/ldVAE_model_projections/scvi_LDVAE_panf_pan_immune/A1_V10_raw_high_var_scvi_YS_updated_panf_gonads_brain_build_donor_organ_corrected_031022.h5ad',
'pan_f_YS_A1_V12_high_var_ldvae':'/nfs/team205/ig7/work_backups/backup_210306/projects/YS/YS_data/YS_panf_gonads_brain_eliv_combined_060922/ldVAE_model_projections/scvi_LDVAE_panf_pan_immune/A1_V12_raw_high_var_scvi_YS_updated_panf_gonads_brain_build_donor_organ_corrected_031022.h5ad',
'fsk_adt_org_int':'/nfs/team205/ig7/projects/fetal_skin/3_160523_probabillistic_projection_organoid_adt_fetl/v3_raw_concatenated_vae_xscvi_features.h5ad'
}

# Variable assignment
feat_use = 'age'
model_key = 'overall_fsk_aft'#'test_low_dim_ipsc_ys'# key for model of choice can be either url or local 
train_x_partition = 'X_scvi' # what partition was the data trained on? To keep simple, for now only accepts 'X'
dyn_std = 1.96 # Dynamic cutoffs using std of the mean for each celltype probability, gives a column notifying user of uncertain labels 1 == 68Ci, 1.96 = 95CI

list(models.keys())

org_key = 'SK'
model_key = 'overall_fsk_adt'
out_dir = './A2_V1_sk_model_pan_organ_heatmaps'
train_x_partition = 'X_scvi'
feat_use = 'age'
# if low_dim, else all below == False
use_varm = '/nfs/team205/ig7/projects/fetal_skin/3_160523_probabillistic_projection_organoid_adt_fetl/A2_V2_ldvae_models/v3_ldvae_obsm_weights.csv'