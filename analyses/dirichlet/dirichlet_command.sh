#!/bin/bash
#BSUB -J dirichlet[13]
#BSUB -q long
#BSUB -G team78-grp
#BSUB -e /nfs/users/nfs_s/sb50/projects/main_repo/repo/projects/clonal_heterogeneity/analysis/pan_donor_analysis/dirichlet/logs/stderr_%I
#BSUB -o /nfs/users/nfs_s/sb50/projects/main_repo/repo/projects/clonal_heterogeneity/analysis/pan_donor_analysis/dirichlet/logs/stdout_%I
#BSUB -n 1
#BSUB -M 128000
#BSUB -R "select[mem>128000] rusage[mem=128000] span[hosts=1]"

# Read list of all donors
echo $LSB_JOBINDEX
readarray donors < ../pilepy/input/donor_list.txt
DONOR="${donors[${LSB_JOBINDEX}-1]}"
DONOR=$(printf $DONOR)
echo $DONOR

# Define input directory
MAIN_DIR=/nfs/users/nfs_s/sb50/mylustre/clonal_heterogeneity/donors/${DONOR}/dirichlet
INPUT_DIR=$MAIN_DIR/input

NUM_BURNIN=10000

# Run R code
R-3.3.0 < Run_Dirichlet_clustering_posthoc.R --vanilla --args \
	$INPUT_DIR/target_alt.csv \
	$INPUT_DIR/target_depth.csv \
	$INPUT_DIR/mut_contexts.txt \
	$MAIN_DIR \
	$NUM_BURNIN
