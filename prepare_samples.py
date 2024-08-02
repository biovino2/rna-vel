'''This script combines the desired subset of samples and saves them as a combined object.

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

    # This mapping has to be performed because the single cell files have different sample names
    sample_mapping = {
    'TDR118reseq': 'TDR118',
    'TDR119reseq': 'TDR119',
    'TDR124reseq': 'TDR124',
    'TDR125reseq': 'TDR125',
    'TDR126': 'TDR126',
    'TDR127': 'TDR127',
    'TDR128': 'TDR128'
    }

    parts = cell_id.split(':')
    cell_barcode = parts[1].replace('x', '')  # Remove trailing 'x' if present

    return f"{cell_barcode}-{sample_mapping[parts[0]]}"


def combine_subset(subset: 'list[str]', comb_file: str):
    """Saves a combined object of the desired subset of samples.

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

    # Combine samples and set them to match the sample ID (rather than loom file ID)
    combine_subset(args.subset, f'data/subsets/{args.filename}.loom')
    subset  = scv.read_loom(f'data/subsets/{args.filename}.loom')
    subset.obs['new_id'] = subset.obs.index.map(replace_sample_code)
    subset.obs.set_index('new_id', inplace=True)

    # Match cell ID format for annotated "clean" subset
    all = sc.read_h5ad('data/all_RNA.h5ad')
    all.obs['dataset'] = all.obs['dataset'].astype(str)
    all.obs['new_id'] = all.obs.index.str.replace(r'-\d+.*$', '', regex=True) + '-' + all.obs['dataset']
    all.obs.set_index('new_id', inplace=True)

    # Subset subset object to only include barcodes from the "clean" subset object
    barcodes = all.obs.index.tolist()
    subset = subset[subset.obs.index.isin(barcodes)].copy()

    # Copy over all obs columns from the clean subset object to the new object
    # First check that the indices are aligned
    additional_obs_df = all.obs.loc[subset.obs.index]
    subset.obs = pd.concat([subset.obs, additional_obs_df], axis=1)
    subset.write_h5ad(f'data/subsets/{args.filename}.h5ad')
    os.remove(f'data/subsets/{args.filename}.loom')  # not needed for downstream steps


if __name__ == '__main__':
    main()
