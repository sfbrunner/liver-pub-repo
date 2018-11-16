#!/bin/bash
#BSUB -J pilepy[1-10]
#BSUB -q long
# -G team78-grp
#BSUB -e output/stderr_%I
#BSUB -o output/stdout_%I
#BSUB -n 1
#BSUB -M 4000
#BSUB -R "select[mem>4000] rusage[mem=4000] span[hosts=1]"

# Read list of all BAM files
echo $LSB_JOBINDEX
readarray bams < input/all_bams.txt
BAM="${bams[${LSB_JOBINDEX}-1]}"
BAM=$(printf $BAM)
echo $BAM

BAMNAME=$(basename $BAM)
echo $BAMNAME

# Define bed file
BED=input/example.bed

# Prepare output directory
OUTDIR=output
mkdir -p $OUTDIR

LOGDIR=$OUTDIR/logs
mkdir -p $LOGDIR

# Set Python env
unset PYTHONPATH
source /your/python/env/bin/activate

python ../../lib/pile.py \
	--bam_file $BAM \
    --bed_regions $BED \
	--output_file $OUTDIR/"$BAMNAME"_pilepy.csv \
	--log_file $LOGDIR/log.$BAMNAME
