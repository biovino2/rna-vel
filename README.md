# rna-vel
toolkit for working with RNA velocity

## Pipeline
There are several steps involved in calculating and visualizing velocity:

1) Preprocessing the raw BAM file from sequencing to get intron/exon counts
2) If you have multiple samples, preparing each one for velocity calculations
3) Calculating velocity
4) Plotting velocity alongside scRNA-seq/ATAC-seq data 

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

## Calculating Velocity
You can choose to use one or more samples when calculating velocity. Run the command below in order to combine your subset of samples:

```
python prepare_samples.py --subset <sample1> <sample2> <sample3> --filename combined_samples
```

This will save a .h5ad file containing all of the samples in 

You can now calculate velocity (sc-velo deterministic) by running the following command:

```
python get_velocity.py --splicecounts <data/splice_counts/combined_samples_splice_counts.h5ad>
```

## Visualizing Velocity
You can then graph the velocity alongside your desired modality with the following command:

```
python graph_velocity.py --sample <combined_samples_velocity.h5ad> --data <rna/atac/joint>
```