#!/bin/bash
: 'This file will look through a single folder, find any samples with BAM files nested, and run velocyto on them.
Folder structure: Sample path > Sample(s) > BAM files, filtered_feature_bc_matrix > barcodes
If you want to use your own barcodes, make a barcodes.tsv or barcodes.tsv.gz file, place it in the sample folder,
and edit the barcode_path variable accordingly.

You can edit the paths in the config.template.cfg file to match your own paths (explanation in README.md and config.template.cfg)

Leah Dorman, Sarah Ancheta, Ben Iovino  7/26/24 CZ-Biohub
'

source config.cfg

# Make data directory and copy processed data for next step
[ -d "data" ] || mkdir -p "data"

#1. Make a working directory (if it doesn't exist)
[ -d "$WDIR" ] || mkdir -p "$WDIR"
cd "$WDIR" || exit

#2. Set up directories for any scripts and slurm logs
[ -d "old_scripts" ] || mkdir -p "old_scripts"
[ -d "slurm_logs" ] || mkdir -p "slurm_logs"

#3. Grab the GTF file and copy into the working directory
if (file $GTF_FILE | grep -q compressed ) ; then
     cp $GTF_FILE annotation.gtf.gz
     gunzip annotation.gtf.gz
else
     cp $GTF_FILE annotation.gtf
fi

#4. Alow a repeat mask, or none
if test -f $REPEAT_MASK; then
  if (file $REPEAT_MASK | grep -q compressed ) ; then
     cp $REPEAT_MASK repeat_mask.gtf.gz
     gunzip repeat_mask.gtf.gz
  else
     cp $REPEAT_MASK repeat_mask.gtf
  fi
  MASK=Present
else
     echo "No Repeat Mask found, running unmasked. May take longer."
     MASK=None
fi

#5. Loop through directory containing samples to find bam files
# Assume directory contains multiple folders, each of which is one sample.
for FOLDER in $SAMPLE_PATH/*; do
     if [ ! -d "$FOLDER" ]; then
          echo "$FOLDER is not a directory, skipping"
          continue
     fi

     # Pull out the sample name from the full path
     SAMPLE=${FOLDER##*/}
     echo "Running velocity on sample $SAMPLE"

     # Move any old contents into a different folder (ignore if it doesn't exist)
     if [ -d "$DESTDIR/$SAMPLE" ]; then
          echo "$DESTDIR/$SAMPLE already exists, moving into a subfolder"
          mkdir -p "$DESTDIR/old_files/"
	     mv $DESTDIR/$SAMPLE/ $DESTDIR/old_files/
     fi

     # Make the output folder
     mkdir -p $DESTDIR/$SAMPLE/
     OUTPUT_PATH=$DESTDIR/$SAMPLE/

     # Copy the bam and barcodes files in here first, in case you don't have write permission in the source folder
     mkdir $SAMPLE
     cp -r $SAMPLE_PATH/$SAMPLE/outs/$BAM_FILE $SAMPLE/
     cp -r $SAMPLE_PATH/$SAMPLE/outs/$BARCODE_PATH $SAMPLE/
     BARCODE_FILE=$(basename $BARCODE_PATH)

# Make an sbatch script to run velocity - still in loop
     cat > $SAMPLE'_velocity.sbatch' <<EOF
#!/bin/bash
#SBATCH --job-name $SAMPLE-velocity # Name for your job
#SBATCH --ntasks 1              # Number of tasks
#SBATCH --cpus-per-task 32      # Number of cpus per task (64)
#SBATCH --time 24:00:00               # Runtime in hr:min:sec
#SBATCH --mem 256G            # Reserve RAM in gigabytes (G)
#SBATCH --partition gpu,cpu         # Partition to submit
#SBATCH --output slurm_logs/log-$SAMPLE-velocity-%j.txt       # Standard out goes to this file
#SBATCH --error slurm_logs/log-$SAMPLE-velocity-%j.txt        # Standard err goes to this file
#SBATCH --mail-user $USER@czbiohub.org     # this is the email you wish to be notified at
#SBATCH --mail-type ALL         # ALL will alert you of job beginning, completion, failure etc
#SBATCH --gpus 0            # Reserve 0 GPUs for usage


# RUN SCRIPT
module load data.science
source activate velocity
module load samtools

#make a sample path in the working directory
[ -d "$WDIR/$SAMPLE/" ] || mkdir -p "$WDIR/$SAMPLE/"


##write code here that identifies the folder that contains the bam file, and cd into it

cd $WDIR/$SAMPLE

##samtools sort -t CB -O BAM -o cellsorted_possorted_genome_bam.bam gex_possorted_bam.bam
##use the gex_possorted_bam.bam instead of the cellsorted, we remove cell sorting step

if [[ "$MASK" == "Present" ]]; then
     echo "Mask found, running velocity"
     velocyto run -b $BARCODE_FILE -o $OUTPUT_PATH -m ../repeat_mask.gtf $BAM_FILE $WDIR/annotation.gtf -e $SAMPLE
else
     echo "no mask file found, running unmasked"
     velocyto run -b $BARCODE_FILE -o $OUTPUT_PATH $BAM_FILE $WDIR/annotation.gtf -e $SAMPLE

fi


EOF

     # sbatch script and move it to old files
     sbatch $SAMPLE'_velocity.sbatch'
     mv $SAMPLE'_velocity.sbatch' old_scripts/
done
