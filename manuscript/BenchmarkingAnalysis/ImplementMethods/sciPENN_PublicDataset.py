## This the following python
## /app/software/Python/3.8.2-GCCcore-9.3.0/bin/python
##
import numpy as np
from matplotlib import pyplot
import os
from copy import deepcopy

from time import time

from math import ceil
from scipy.stats import spearmanr, gamma, poisson

from anndata import AnnData, read_h5ad
import scanpy as sc
from scanpy import read
import pandas as pd

from torch.utils.data import DataLoader, TensorDataset
from torch import tensor
from torch.cuda import is_available

from sciPENN.sciPENN_API import sciPENN_API
import anndata
import json
import psutil

## original tutorial
# """Read in the data"""
# data_path = './manuscript/scripts/sciPENN_tutorial/data/'

# adata_gene = sc.read(data_path + "pbmc_gene.h5ad")
# adata_protein = sc.read(data_path + "pbmc_protein.h5ad")

# doublet_bool = (adata_gene.obs['celltype.l2'] != 'Doublet')

# adata_gene = adata_gene[doublet_bool].copy()
# adata_protein = adata_protein[doublet_bool].copy()


# """Create training and test"""

# train_bool = [x in ['P1', 'P3', 'P4', 'P7'] for x in adata_gene.obs['donor']]

# adata_gene_set1 = adata_gene[train_bool].copy()
# adata_protein_set1 = adata_protein[train_bool].copy()
# adata_gene_set2 = adata_gene[np.invert(train_bool)].copy()
# adata_protein_set2 = adata_protein[np.invert(train_bool)].copy()

# common_proteins = adata_protein_set1.var.index
# set1only_proteins = np.random.choice(common_proteins, len(common_proteins)//3, False)
# common_proteins = np.setdiff1d(common_proteins, set1only_proteins)
# set2only_proteins = np.random.choice(common_proteins, len(common_proteins)//2, False)

# set1only_proteins = set(set1only_proteins)
# set2only_proteins = set(set2only_proteins)

# keep_set1 = [x not in set2only_proteins for x in adata_protein_set1.var.index]
# keep_set2 = [x not in set1only_proteins for x in adata_protein_set1.var.index]

# adata_protein_set1 = adata_protein_set1[:, keep_set1].copy()
# adata_protein_set2 = adata_protein_set2[:, keep_set2].copy()


# sciPENN = sciPENN_API(gene_trainsets = [adata_gene_set1, adata_gene_set2], 
#                       protein_trainsets = [adata_protein_set1, adata_protein_set2], 
#                       train_batchkeys = ['donor', 'donor'])


# sciPENN.train(quantiles = [0.1, 0.25, 0.75, 0.9], n_epochs = 10000, ES_max = 12, decay_max = 6, 
#              decay_step = 0.1, lr = 10**(-3), weights_dir = "pbmc_to_pbmcINTEGRATE", load = True)


# imputed_test = sciPENN.impute()


# imputed_test.X

# imputed_test.obs['Dataset']
# imputed_test.obs['batch']



# proteins = imputed_test.var.index

# proteins1 = proteins[imputed_test.var['Dataset 1']] #get proteins sequenced in Dataset 1
# proteins2 = proteins[imputed_test.var['Dataset 2']] #get proteins sequenced in Dataset 1


# ds1_cells = imputed_test.obs['Dataset'] == 'Dataset 1'
# ds2_cells = imputed_test.obs['Dataset'] == 'Dataset 2'

# ds1_pred, ds1_seq = np.invert(imputed_test.var['Dataset 1']), imputed_test.var['Dataset 1']
# ds2_pred, ds2_seq = np.invert(imputed_test.var['Dataset 2']), imputed_test.var['Dataset 2']

# pred1 = imputed_test[ds1_cells, ds1_pred] #imputed protein array in dataset 1
# sequenced1 = imputed_test[ds1_cells, ds1_seq] #sequenced protein array in dataset 1
# pred2 = imputed_test[ds2_cells, ds2_pred] #imputed protein array in dataset 2
# sequenced2 = imputed_test[ds2_cells, ds2_seq] #sequenced protein array in dataset 2

# embedding = sciPENN.embed()
# embedding.obs['batch']

# q10_pred = imputed_test[ds1_cells, ds1_pred].layers['q10'] #get q10 for imputed proteins from reference 1
# q10_truth = imputed_test[ds1_cells, ds1_seq].layers['q10'] #get q10 for sequenced proteins from reference 1, not useful



# ## rewrite according to my scenario
# """Read in the data"""
# data_path = './manuscript/scripts/sciPENN_tutorial/data/'
# adata_gene = sc.read(data_path + "pbmc_gene.h5ad")
# adata_protein = sc.read(data_path + "pbmc_protein.h5ad")
# doublet_bool = (adata_gene.obs['celltype.l2'] != 'Doublet')
# adata_gene = adata_gene[doublet_bool].copy()
# adata_protein = adata_protein[doublet_bool].copy()

# sciPENN = sciPENN_API(gene_trainsets = [adata_gene], 
#                       protein_trainsets = [adata_protein], 
#                       train_batchkeys = ['donor'])


# sciPENN.train(quantiles = [0.1, 0.25, 0.75, 0.9], n_epochs = 10000, ES_max = 12, decay_max = 6, 
#              decay_step = 0.1, lr = 10**(-3), weights_dir = "pbmc_to_pbmcINTEGRATE", load = True)

# sciPENN.train(n_epochs = 10000, ES_max = 12, decay_max = 6, 
#              decay_step = 0.1, lr = 10**(-3), weights_dir = "pbmc_to_pbmcINTEGRATE")

# embedding = sciPENN.embed()

# embedding.obs['batch']

# imputed_test = sciPENN.impute()

## define memory usage
def log_memory_usage(file):
    memory = psutil.virtual_memory()
    file.write(f'Total memory: {memory.total / (1024.0 ** 3)} GB\n')
    file.write(f'Used memory: {memory.used / (1024.0 ** 3)} GB\n')
    file.write(f'Available memory: {memory.available / (1024.0 ** 3)} GB\n')


## Public Data 13datasets
out_path = './manuscript/results/public13Dataset_CITEseq/CSV/'


run_name = 'public13Dataset_CITEseq'
adts = pd.read_csv('./results/publicData_CITEseq/CSV/adt_data_common_RawCount_'+ run_name + '.csv', index_col=0)
gex = pd.read_csv('./results/publicData_CITEseq/CSV/rna_data_common_RawCount_'+ run_name + '.csv', index_col=0)
obs = pd.read_csv('./results/publicData_CITEseq/CSV/adt_feature_common_RawCount_'+ run_name + '.csv', index_col=0)
# obs.reset_index(drop=True, inplace=True)
# gex.reset_index(drop=True, inplace=True)
# adts.reset_index(drop=True, inplace=True)

adts.index = obs.index
gex.index = obs.index

adata_gene = anndata.AnnData(gex, obs)
adata_gene.index = adata_gene.obs_names
adata_protein = anndata.AnnData(adts, obs)
adata_protein.index = adata_protein.obs_names

batch_key = 'study_name'
with open(out_path + 'sciPENN_memory_usage_' + run_name + '_' + batch_key + '_GPU.csv', 'w') as f:
    log_memory_usage(f)
    start = time()
    
    sciPENN = sciPENN_API(gene_trainsets = [adata_gene], 
                        protein_trainsets = [adata_protein], 
                        train_batchkeys = [batch_key], use_gpu = True, min_cells = False, min_genes = False)
    log_memory_usage(f)
    sciPENN.train(quantiles = [0.1, 0.25, 0.75, 0.9], n_epochs = 10000, ES_max = 12, decay_max = 6, decay_step = 0.1, lr = 10**(-3), weights_dir = "sciPENN_public13Dataset_CITEseq_study_name_GPU", load = True)
    embedding = sciPENN.embed()
    imputed_test = sciPENN.impute()
    log_memory_usage(f)
    end = time()
    print(f'Seconds to complete: {round(end - start,0)}')
    train_time = round(end - start, 0)

with open(out_path + 'sciPENN_runtime_' + run_name + '_' + batch_key + '_GPU.csv', 'w') as f:
    json.dump(train_time, f)    

imputed_res = pd.DataFrame(imputed_test.X, index=imputed_test.obs_names, columns=imputed_test.var_names)
imputed_res.to_csv(out_path + 'sciPENN_Imputed_Protein_' + run_name + '_' + batch_key + '.csv')
embedding_res = pd.DataFrame(embedding.X, index=embedding.obs_names, columns=embedding.var_names)
embedding_res.to_csv(out_path + 'sciPENN_Embedding_' + run_name + '_' + batch_key + '.csv')
embedding_obs_res = pd.DataFrame(embedding.obs)
embedding_obs_res.to_csv(out_path + 'sciPENN_Embedding_DataFeature_' + run_name + '_' + batch_key + '.csv')

batch_key = 'sample'
with open(out_path + 'sciPENN_memory_usage_' + run_name + '_' + batch_key + '_CPU.csv', 'w') as f:
    log_memory_usage(f)
    start = time()
    sciPENN = sciPENN_API(gene_trainsets = [adata_gene], 
                        protein_trainsets = [adata_protein], 
                        train_batchkeys = [batch_key], use_gpu = False, min_cells = False, min_genes = False)
    log_memory_usage(f)
    sciPENN.train(quantiles = [0.1, 0.25, 0.75, 0.9], n_epochs = 10000, ES_max = 12, decay_max = 6, decay_step = 0.1, lr = 10**(-3), weights_dir = "sciPENN_public13Dataset_CITEseq_sample_CPU", load = True)
    embedding = sciPENN.embed()
    imputed_test = sciPENN.impute()
    log_memory_usage(f)
    end = time()
    print(f'Seconds to complete: {round(end - start,0)}')
    train_time = round(end - start, 0)

with open(out_path + 'sciPENN_runtime_' + run_name + '_' + batch_key + '_CPU.csv', 'w') as f:
    json.dump(train_time, f)  


batch_key = 'sample'
with open(out_path + 'sciPENN_memory_usage_' + run_name + '_' + batch_key + '_GPU.csv', 'w') as f:
    log_memory_usage(f)
    start = time()
    sciPENN = sciPENN_API(gene_trainsets = [adata_gene], 
                        protein_trainsets = [adata_protein], 
                        train_batchkeys = [batch_key], use_gpu = True, min_cells = False, min_genes = False)
    log_memory_usage(f)
    sciPENN.train(quantiles = [0.1, 0.25, 0.75, 0.9], n_epochs = 10000, ES_max = 12, decay_max = 6, decay_step = 0.1, lr = 10**(-3), weights_dir = "sciPENN_public13Dataset_CITEseq_sample_GPU", load = True)
    embedding = sciPENN.embed()
    imputed_test = sciPENN.impute()
    log_memory_usage(f)
    end = time()
    print(f'Seconds to complete: {round(end - start,0)}')
    train_time = round(end - start, 0)

with open(out_path + 'sciPENN_runtime_' + run_name + '_' + batch_key + '_GPU.csv', 'w') as f:
    json.dump(train_time, f)    

imputed_res = pd.DataFrame(imputed_test.X, index=imputed_test.obs_names, columns=imputed_test.var_names)
imputed_res.to_csv(out_path + 'sciPENN_Imputed_Protein_' + run_name + '_' + batch_key + '.csv')
embedding_res = pd.DataFrame(embedding.X, index=embedding.obs_names, columns=embedding.var_names)
embedding_res.to_csv(out_path + 'sciPENN_Embedding_' + run_name + '_' + batch_key + '.csv')
embedding_obs_res = pd.DataFrame(embedding.obs)
embedding_obs_res.to_csv(out_path + 'sciPENN_Embedding_DataFeature_' + run_name + '_' + batch_key + '.csv')


with open(out_path + 'sciPENN_memory_usage_' + run_name + '_' + batch_key + '_CPU.csv', 'w') as f:
    log_memory_usage(f)
    start = time.time()
    
    sciPENN = sciPENN_API(gene_trainsets = [adata_gene], 
                        protein_trainsets = [adata_protein], 
                        train_batchkeys = [batch_key], use_gpu = False, min_cells = False, min_genes = False)
    sciPENN.train(quantiles = [0.1, 0.25, 0.75, 0.9], n_epochs = 10000, ES_max = 12, decay_max = 6, decay_step = 0.1, lr = 10**(-3), weights_dir = "sciPENN_public13Dataset_CITEseq_sample_CPU", load = True)
    embedding = sciPENN.embed()
    imputed_test = sciPENN.impute()
    end = time.time()
    print(f'Seconds to complete: {round(end - start,0)}')
    train_time = round(end - start, 0)

with open(out_path + 'sciPENN_runtime_' + run_name + '_' + batch_key + '_CPU.csv', 'w') as f:
    json.dump(train_time, f)  



## T cell dataset
run_name = 'public13Dataset_CITEseq_Tcelllargeprop' #'public13Dataset_CITEseq_Tcelllargeprop_cd8T'
out_path = './manuscript/results/' + run_name + '/CSV/'

adts = pd.read_parquet('./manuscript/results/' + run_name + '/CSV/adt_data_common_RawCount_'+ run_name + '.parquet')
gex = pd.read_parquet('./manuscript/results/' + run_name + '/CSV/rna_data_common_RawCount_'+ run_name + '.parquet')
obs = pd.read_parquet('./manuscript/results/' + run_name + '/CSV/adt_feature_common_RawCount_'+ run_name + '.parquet')

adts.index = obs.index
gex.index = obs.index

adata_gene = anndata.AnnData(gex, obs)
adata_gene.index = adata_gene.obs_names
adata_protein = anndata.AnnData(adts, obs)
adata_protein.index = adata_protein.obs_names

batch_key = 'sample' # 'sample' #'study_name' #
with open(out_path + 'sciPENN_memory_usage_' + run_name + '_' + batch_key + '_GPU.csv', 'w') as f:
    log_memory_usage(f)
    start = time()
    sciPENN = sciPENN_API(gene_trainsets = [adata_gene], 
                        protein_trainsets = [adata_protein], 
                        train_batchkeys = [batch_key], use_gpu = True, min_cells = False, min_genes = False)
    log_memory_usage(f)
    sciPENN.train(quantiles = [0.1, 0.25, 0.75, 0.9], n_epochs = 10000, ES_max = 12, decay_max = 6, decay_step = 0.1, lr = 10**(-3), weights_dir = "sciPENN_public12Dataset_CITEseq_notriana_study_name", load = True)
    log_memory_usage(f)
    embedding = sciPENN.embed()
    imputed_test = sciPENN.impute()
    log_memory_usage(f)
    end = time()
    print(f'Seconds to complete: {round(end - start,0)}')
    train_time = round(end - start, 0)

with open(out_path + 'sciPENN_runtime_' + run_name + '_' + batch_key + '_GPU.csv', 'w') as f:
    json.dump(train_time, f)   


imputed_res = pd.DataFrame(imputed_test.X, index=imputed_test.obs_names, columns=imputed_test.var_names)
imputed_res.to_csv(out_path + 'sciPENN_Imputed_Protein_' + run_name + '_' + batch_key + '.csv')

embedding_res = pd.DataFrame(embedding.X, index=embedding.obs_names, columns=embedding.var_names)
embedding_res.to_csv(out_path + 'sciPENN_Embedding_' + run_name + '_' + batch_key + '.csv')

embedding_obs_res = pd.DataFrame(embedding.obs)
embedding_obs_res.to_csv(out_path + 'sciPENN_Embedding_DataFeature_' + run_name + '_' + batch_key + '.csv')
