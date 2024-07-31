'''Graph velocity results from scvelo-deterministic.

Sarah Ancheta, Ben Iovino   07/29/24    CZ-Biohub
'''

import argparse
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import os
import scvelo as scv


def main():
    """
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('--sample', type=str, help='File containing processed single-cell data')
    parser.add_argument('--data', type=str, help='Data type (e.g. rna, atac, joint)')
    args = parser.parse_args()

    # Read individual sample (single-cell data)
    sample = f'data/{args.sample}_processed_RNA.h5ad'
    sample = sc.read_h5ad(sample)

    # Add sample ID to cell ID's
    sample.obs['new_id'] = sample.obs.index.str.replace(r'-\d+.*$', '', regex=True) + f'-{args.sample}'
    sample.obs.set_index('new_id', inplace=True)

    # Read all samples (single-cell data)
    all = 'data/all_processed_RNA.h5ad'
    all = sc.read_h5ad(all)
    all.obs['dataset'] = all.obs['dataset'].astype(str)
    all.obs['new_id'] = all.obs.index.str.replace(r'-\d+.*$', '', regex=True) + '-' + all.obs['dataset']
    all.obs.set_index('new_id', inplace=True)

    # Remove cells from individual sample that are not in master dataset (all)
    sample = sample[sample.obs.index.isin(all.obs.index)]

    # Get annotation color from master dataset
    sample.obs['annotation_ML_coarse'] = all.obs.loc[sample.obs.index, 'annotation_ML_coarse']

    # Read velocity data
    velocity = f'data/velocity/{args.sample}_velocity.h5ad'
    velocity = sc.read_h5ad(velocity)

    # Make dictionary of UMAP coordinates for each cell in individual sample
    sample_umap = sample.obsm[f'X_umap.{args.data}']
    sample_umap_dict = dict(zip(sample.obs.index, sample_umap))

    # Assign UMAP coords to velocity data
    velocity.obs[f'X_umap.{args.data}'] = velocity.obs.index.map(sample_umap_dict)
    velocity.obsm[f'X_umap.{args.data}'] = np.array(list(velocity.obs[f'X_umap.{args.data}']))  # move to obsm

    # Plot velocity
    if not os.path.exists('data/graphs'):
        os.mkdir('data/graphs')
    scv.pl.velocity_embedding_stream(
        velocity,
        basis=f'X_umap.{args.data}',
        c='annotation_ML_coarse',
        legend_loc='right',
        title=f'RNA Velocity {args.sample}-{args.data}',
        save=f'data/graphs/{args.sample}-{args.data}.pdf')


if __name__ == '__main__':
    main()
