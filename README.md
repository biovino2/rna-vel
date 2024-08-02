# rna-vel
toolkit for working with RNA velocity

## Installation
You can install this repository with the following commands:

```
git clone https://github.com/biovino2/rna-vel
cd rna-vel/
conda env create --file env.yml
conda activate rna-vel
```

## Pipeline
There are several steps involved in calculating and visualizing velocity:

1) Preprocessing the raw BAM file from sequencing to get intron/exon counts
2) Copying over single-cell data
3) Preparing your sample(s) for velocity calculations
4) Calculating velocity
5) Plotting velocity alongside scRNA-seq/ATAC-seq data

## Preprocessing Data 
Copy the template configuration file:

```sh
cp config.template.cfg config.cfg
```

And edit 'config.cfg' to set your directory paths before preprocessing.

You can then run the preprocessing script as such:

```
bash preprocessing.sh
```

Which will run velocyto on any samples (bam files) that are located in the sample directory. If there is more than one sample located in this directory, they will be extracted and processed.

## Copying Data
Make sure to copy over your processed single-cell data to the data/ directory, which contain the desired cells and embedding coordinates. This pipeline assumes a certain naming schema in the .obsm frame for coordinates, but this can be changed in graph_velocity.py.

## Calculating Velocity
You can choose to use one or more samples when calculating velocity. Run the command below in order to combine your subset of samples:

```
python prepare_samples.py --subset <sample1> <sample2> <sample3> --filename combined_samples
```

NOTE: the sample names when calling --subset should match the folder name in loom_files/ or else they will not be recognized. This will save a .h5ad file ready for velocity calculations in subsets/. You can also run prepare_samples.py on a single sample.

You can now calculate velocity (sc-velo deterministic) by running the following command:

```
python get_velocity.py --subset data/subsets/combined_samples.h5ad
```

## Visualizing Velocity
You can then graph the velocity alongside your desired modality with the following command:

```
python graph_velocity.py --sample <combined_samples_velocity.h5ad> --data <rna/atac/joint>
```