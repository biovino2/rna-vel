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
    parser.add_argument('--splicecounts', type=str, help='Path to preprocessed splice counts file')
    args = parser.parse_args()

    if not os.path.exists('data/velocity'):
        os.mkdir('data/velocity')
    adata_fresh = sc.read_h5ad(args.splicecounts)
    adata_fresh.var_names_make_unique()
    adata = adata_fresh.copy()

    # Preprocess data
    scv.pp.filter_genes(adata, min_shared_counts=20)
    scv.pp.normalize_per_cell(adata)
    scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
    scv.pp.log1p(adata)

    # Calculate moments
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)

    # Calculate velocities
    scv.tl.velocity(adata, mode='deterministic')
    scv.tl.velocity_graph(adata)

    # Write to file
    sample = os.path.basename(args.splicecounts).split('_')[0]
    adata.write(f'data/velocity/{sample}_velocity.h5ad')


if __name__ == '__main__':
    main()
