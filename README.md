# rna-vel
toolkit for working with RNA velocity

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