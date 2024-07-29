'''This script combines the desired subset of samples and saves them as a combined adata object.

Sarah Ancheta, Ben Iovino   07/29/24    CZ-Biohub
'''

import argparse
import logging
import loompy
import os
import pandas as pd
import scanpy as sc
import scvelo as scv


def replace_sample_code(cell_id: str) -> str:
    """ Replace loom file ID with sample ID

    Args:
        cell_id (str): loom file ID

    Return:
        str: sample ID
    """

    sample_mapping = {
    'SK7J2': 'TDR118',
    '7TH98': 'TDR119',
    '47VDB': 'TDR124',
    'BOP08': 'TDR125',
    'H12VC': 'TDR126',
    'WX1JE': 'TDR127',
    '48LIH': 'TDR128'
    }

    parts = cell_id.split(':')
    sample_code = parts[0].split('_')[-1]  # Extract the sample code
    cell_barcode = parts[1].replace('x', '')  # Remove trailing 'x' if present

    return f"{cell_barcode}-{sample_mapping.get(sample_code, sample_code)}"


def combine_subset(subset: list[str], comb_file: str):
    """Saves a combined adata object of the desired subset of samples.

    Args:
        subset (list): List of loom files to combine
        filename (str): Name of the combined loom file
    """

    loom_files: list[str] = []
    for direc in os.listdir('data/loom_files'):
        if direc in subset:  # ignore samples not in subset
            for file in os.listdir(f'data/loom_files/{direc}'):
                logging.info('Adding %s to list of loom files...', file)
                loom_files.append(f'data/loom_files/{direc}/{file}')

    logging.info('Combining loom files...')
    loompy.combine(loom_files, comb_file, key="Accession")


def main():
    """
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('--subset', nargs='+', type=str, help='List of samples to subset')
    parser.add_argument('--filename', type=str, help='Name of the combined loom file')
    args = parser.parse_args()

    if not os.path.exists('data/subsets'):
        os.mkdir('data/subsets')
    if not os.path.exists('data/splice_counts'):
        os.mkdir('data/splice_counts')

    # Combine samples and set them to match the sample ID (rather than loom file ID)
    combine_subset(args.subset, f'data/subsets/{args.filename}.loom')
    adata  = scv.read_loom(f'data/subsets/{args.filename}.loom')
    adata.obs['new_id'] = adata.obs.index.map(replace_sample_code)
    adata.obs.set_index('new_id', inplace=True)

    # Match cell ID format for annotated "clean" adata
    clean_adata = sc.read_h5ad('data/proc_data.h5ad')
    clean_adata.obs['dataset'] = clean_adata.obs['dataset'].astype(str)
    clean_adata.obs['new_id'] = clean_adata.obs.index.str.replace(r'-\d+.*$', '', regex=True) + '-' + clean_adata.obs['dataset']
    clean_adata.obs.set_index('new_id', inplace=True)

    # Subset adata object to only include barcodes from the "clean" adata object
    barcodes = clean_adata.obs.index.tolist()
    adata = adata[adata.obs.index.isin(barcodes)].copy()

    # Copy over all obs columns from the clean adata object to the new object
    # First check that the indices are aligned
    additional_obs_df = clean_adata.obs.loc[adata.obs.index]
    adata.obs = pd.concat([adata.obs, additional_obs_df], axis=1)
    adata.write_h5ad(f'data/{args.filename}_splice_counts.h5ad')

if __name__ == '__main__':
    main()
