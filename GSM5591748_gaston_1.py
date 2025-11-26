import sys
import os
from collections import defaultdict
import pandas as pd
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
from glmpca import glmpca
from itertools import combinations
import torch

import sys
from importlib import reload

import gaston
from gaston import neural_net,cluster_plotting, dp_related, segmented_fit, restrict_spots, model_selection
from gaston import binning_and_plotting, isodepth_scaling, run_slurm_scripts, parse_adata
from gaston import spatial_gene_classification, plot_cell_types, filter_genes, process_NN_output

import seaborn as sns
import math
######################################################################################################################################################################

data_folder='/home/qiaocy/bioinfo/GSM5591748' # folder containing the data
use_RGB=True # set to False if you do not want to use RGB as features
spot_umi_threshold=50 # filter out cells with low UMIs

counts_mat, coords_mat, gene_labels, rgb_mean=parse_adata.get_gaston_input_adata(data_folder, get_rgb=use_RGB, spot_umi_threshold=spot_umi_threshold) # get RGB mean if using RGB

# save matrices
# os.makedirs('GJB-LC', exist_ok=True)
np.save('/home/qiaocy/bioinfo/GSM5591748/counts_mat.npy', counts_mat) # save matrices
np.save('/home/qiaocy/bioinfo/GSM5591748/coords_mat.npy', coords_mat) # save matrices
np.save('/home/qiaocy/bioinfo/GSM5591748/gene_labels.npy', gene_labels) # save matrices
print('saved matrices')

######################################################################################################################################################################

# counts_mat = np.load('/home/qiaocy/bioinfo/GSM5591748/counts_mat.npy', allow_pickle = True)
# coords_mat = np.load('/home/qiaocy/bioinfo/GSM5591748/coords_mat.npy', allow_pickle = True)
# gene_labels = np.load('/home/qiaocy/bioinfo/GSM5591748/gene_labels_008um.npy', allow_pickle = True)

# counts_mat = counts_mat.item()
######################################################################################################################################################################

num_dims=5 # number of PCs to use
penalty=20 # may need to increase if this is too small

# CHANGE THESE PARAMETERS TO REDUCE RUNTIME
num_iters=1 # number of iterations for GLM-PCA
eps=1e-2 # tolerance for convergence
num_genes=5000 # number of genes to use for GLM-PCA

counts_mat_glmpca=counts_mat[:,np.argsort(np.sum(counts_mat, axis=0).A1)[-num_genes:]] 

counts_mat_glmpca = counts_mat_glmpca.toarray()
counts_mat = counts_mat.toarray()

######################################################################################################################################################################
# glmpca_res=glmpca.glmpca(counts_mat_glmpca.T, num_dims, fam="poi", penalty=penalty, verbose=True,
#                         ctl = {"maxIter":num_iters, "eps":eps, "optimizeTheta":True})
# A = glmpca_res['factors'] # should be of size N x num_dims, where each column is a PC

# if use_RGB:
#     A=np.hstack((A,rgb_mean)) # attach to RGB mean
# np.save('/home/qiaocy/bioinfo/GSM5591748/glmpca_008um_poi.npy', A) # save A

######################################################################################################################################################################
# visualize top GLM-PCs and RGB mean
# rotated_coords=dp_related.rotate_by_theta(coords_mat, -np.pi/2)
# R=2
# C=4
# fig,axs=plt.subplots(R,C,figsize=(20,10))
# for r in range(R):
#     for c in range(C):
#         i=r*C+c  
#         axs[r,c].scatter(rotated_coords[:,0], rotated_coords[:,1], c=A[:,i],cmap='Reds',s=3)
#         if i < num_dims:
#             axs[r,c].set_title(f'GLM-PC{i}')
#         else:
#             axs[r,c].set_title('RGB'[i-num_dims])

# plt.savefig('/home/qiaocy/bioinfo/GSM5591748/glmpca_008um_poi.png')

######################################################################################################################################################################
num_dims=5
clip=0.01 # have to clip values to be very small!

A = parse_adata.get_top_pearson_residuals(num_dims,counts_mat,coords_mat,clip=clip)
if use_RGB:
    A=np.hstack((A,rgb_mean)) # attach to RGB mean
np.save('/home/qiaocy/bioinfo/GSM5591748/analytic_pearson_poi.npy', A)
print('saved analytic pearson residuals')

######################################################################################################################################################################

# visualize top GLM-PCs
rotated_coords=dp_related.rotate_by_theta(coords_mat, -np.pi/2)
R=2
C=4
fig,axs=plt.subplots(R,C,figsize=(20,10))
for r in range(R):
    for c in range(C):
        i=r*C+c
        axs[r,c].scatter(rotated_coords[:,0], rotated_coords[:,1], c=A[:,i],cmap='Reds',s=3)
        axs[r,c].set_title(f'PC{i}')

plt.savefig('/home/qiaocy/bioinfo/GSM5591748/analytic_pearson_poi.png')
print('plotted analytic pearson residuals')