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

path_to_glmpca = '/home/qiaocy/bioinfo/GSM5591748/analytic_pearson_poi.npy'
path_to_coords = '/home/qiaocy/bioinfo/GSM5591748/coords_mat.npy'

A = np.load(path_to_glmpca)
S = np.load(path_to_coords)

S_torch, A_torch = neural_net.load_rescale_input_data(S, A)

conda_environment = "gaston_try"
path_to_conda_folder = "/home/qiaocy/miniconda3/envs/gaston_try/bin/activate"

isodepth_arch = [20, 20]
expression_arch = [20, 20]
epochs = 10000
checkpoint = 500
output_dir = '/home/qiaocy/bioinfo/GSM5591748/output'
optimizer = "adam"
num_restarts = 30

time_to_train="0-01:00:00"

run_slurm_scripts.train_NN_parallel(path_to_coords, path_to_glmpca, isodepth_arch, expression_arch, 
                      output_dir, conda_environment, path_to_conda_folder,
                      epochs=epochs, checkpoint=checkpoint, 
                      num_seeds=num_restarts,time=time_to_train)

seed_list = range(num_restarts)
for seed in seed_list:
    print(f'training neural network for seed {seed}')
    out_dir_seed = f"{output_dir}/rep{seed}"
    os.makedirs(out_dir_seed, exist_ok=True)
    mod, loss_list = neural_net.train(S_torch, A_torch,
                                      S_hidden_list=isodepth_arch, A_hidden_list=expression_arch,
                                      epochs=epochs, checkpoint=checkpoint,
                                      save_dir=out_dir_seed, optim=optimizer, seed=seed, save_final=True)

gaston_model, A, S = process_NN_output.process_files('/home/qiaocy/bioinfo/GSM5591748/output')

counts_mat = np.load('/home/qiaocy/bioinfo/GSM5591748/counts_mat.npy', allow_pickle=True)
coords_mat = np.load('/home/qiaocy/bioinfo/GSM5591748/coords_mat.npy', allow_pickle=True)
gene_labels = np.load('/home/qiaocy/bioinfo/GSM5591748/gene_labels.npy', allow_pickle=True)

# 获取最佳 num_layers
ll_values = model_selection.plot_ll_curve(gaston_model, A, S, max_domain_num=8, start_from=2)
save_path = '/home/qiaocy/bioinfo/GSM5591748/ll_curve.png'
plt.savefig(save_path)

# 提取 Kneedle 点作为 num_layers
# from kneed import KneeLocator
# domain_nums = range(2, 9)  # 从 2 到 8
# kneedle = KneeLocator(domain_nums, ll_values, curve='convex', direction='decreasing')
# num_layers = kneedle.knee
# print(f"Optimal number of layers (Kneedle): {num_layers}")

gaston_isodepth, gaston_labels = dp_related.get_isodepth_labels(gaston_model, A, S, num_layers)
gaston_isodepth = np.max(gaston_isodepth) - 1 * gaston_isodepth
gaston_labels = (num_layers - 1) - gaston_labels

counts_mat = counts_mat.item()

# 创建包含每个点信息的列表
spot_info = []
for i in range(len(gaston_labels)):
    spot = {
        'coords': coords_mat[i],          # 空间坐标 (x, y)
        'counts': counts_mat[i].toarray().flatten() if hasattr(counts_mat[i], 'toarray') else counts_mat[i],  # 基因表达
        'layer': gaston_labels[i],        # 层标签
        'isodepth': gaston_isodepth[i]    # 等深度值（可选）
    }
    spot_info.append(spot)

# 保存为 .npy 文件
output_file = '/home/qiaocy/bioinfo/GSM5591748/spot_info.npy'
np.save(output_file, spot_info, allow_pickle=True)
print(f"Saved spot information to {output_file}")

show_streamlines = True
rotate = np.radians(0)
arrowsize = 2
cluster_plotting.plot_isodepth(gaston_isodepth, S, gaston_model, figsize=(7, 6), streamlines=show_streamlines,
                               rotate=rotate, arrowsize=arrowsize, neg_gradient=True)
save_path = '/home/qiaocy/bioinfo/GSM5591748/isodepth.png'
plt.savefig(save_path)

domain_colors = ['plum', 'cadetblue', '#F3D9DC', 'dodgerblue', '#F44E3F']
cluster_plotting.plot_clusters(gaston_labels, S, figsize=(6, 6),
                               colors=domain_colors, s=20, lgd=False,
                               show_boundary=True, gaston_isodepth=gaston_isodepth, boundary_lw=5, rotate=rotate)
save_path = '/home/qiaocy/bioinfo/GSM5591748/clusters.png'
plt.savefig(save_path)

isodepth_min = (1/3)*gaston_isodepth.max()
isodepth_max = (2/3)*gaston_isodepth.max()

cluster_plotting.plot_clusters_restrict(gaston_labels, S, gaston_isodepth,
                                        isodepth_min=isodepth_min, isodepth_max=isodepth_max, figsize=(6, 6),
                                        colors=domain_colors, s=20, lgd=False, rotate=rotate)
save_path = '/home/qiaocy/bioinfo/GSM5591748/clusters_restricted.png'
plt.savefig(save_path)

if isinstance(counts_mat, np.ndarray) and counts_mat.dtype == object:
    counts_mat = counts_mat.item()
print(type(counts_mat))

adjust_physical = True
scale_factor = 100
plot_isodepth = True
show_streamlines = True
rotate = np.radians(0)
arrowsize = 1
n_neighbors = min(1000, len(S) - 1)

counts_mat_restrict, coords_mat_restrict, gaston_isodepth_restrict, gaston_labels_restrict, S_restrict = restrict_spots.restrict_spots(
    counts_mat, coords_mat, S, gaston_isodepth, gaston_labels,
    isodepth_min=isodepth_min, isodepth_max=isodepth_max,
    adjust_physical=adjust_physical, scale_factor=scale_factor,
    plot_isodepth=plot_isodepth, show_streamlines=show_streamlines,
    gaston_model=gaston_model, rotate=rotate, figsize=(6, 3),
    arrowsize=arrowsize, neg_gradient=True)

save_path = '/home/qiaocy/bioinfo/GSM5591748/clusters_restricted_streamlines.png'
plt.savefig(save_path)