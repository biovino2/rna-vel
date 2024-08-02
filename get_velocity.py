'''Run scvelo-deterministic on your desired subset of samples.

Sarah Ancheta, Ben Iovino   07/29/24    CZ-Biohub
'''

import argparse
import scanpy as sc
import scvelo as scv
import os


def main():
    """
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('--subset', type=str, help='Path to sample subset file')
    args = parser.parse_args()

    if not os.path.exists('data/velocity'):
        os.mkdir('data/velocity')
    subset_fresh = sc.read_h5ad(args.subset)
    subset_fresh.var_names_make_unique()
    subset = subset_fresh.copy()

    # Preprocess data
    scv.pp.filter_genes(subset, min_shared_counts=20)
    scv.pp.normalize_per_cell(subset)
    scv.pp.filter_genes_dispersion(subset, n_top_genes=2000)
    scv.pp.log1p(subset)

    # Calculate moments
    scv.pp.filter_and_normalize(subset, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(subset, n_pcs=30, n_neighbors=30)
    sc.pp.neighbors(subset, n_neighbors=30, n_pcs=30)

    # Calculate velocities
    scv.tl.velocity(subset, mode='deterministic')
    scv.tl.velocity_graph(subset)

    # Write to file
    sample = os.path.basename(args.subset).split('.')[0]
    subset.write(f'data/velocity/{sample}_velocity.h5ad')


if __name__ == '__main__':
    main()
