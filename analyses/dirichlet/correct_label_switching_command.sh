#!/bin/bash
#BSUB -J ls_correction[1-9]
#BSUB -q normal
#BSUB -e /nfs/users/nfs_s/sb50/projects/main_repo/repo/projects/clonal_heterogeneity/analysis/pan_donor_analysis/dirichlet/ls_logs/stderr_%I
#BSUB -o /nfs/users/nfs_s/sb50/projects/main_repo/repo/projects/clonal_heterogeneity/analysis/pan_donor_analysis/dirichlet/ls_logs/stdout_%I
#BSUB -n 1
#BSUB -M 16000
#BSUB -R "select[mem>16000] rusage[mem=16000] span[hosts=1]"

# Read list of all donors
echo $LSB_JOBINDEX
readarray LS_LIST < input/label_switching_sessions.txt
LS_PATH="${LS_LIST[${LSB_JOBINDEX}-1]}"
LS_PATH=$(printf $LS_PATH)
echo $LS_PATH

# Run R code
R-3.3.0 < Run_Dirichlet_label_switching_correction.R --vanilla --args $LS_PATH
