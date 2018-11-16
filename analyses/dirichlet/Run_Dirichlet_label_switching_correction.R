# Run Dirichlet clustering
# --
# Runs Dirichlet clustering and label switching on low input mutation data
#
# Args (command line):
#   1: Path to Rsession.dat which to pick up from
# 
# If command line args are not provided, default file locations are used.
# 
# --
# /// Author --- PETER CAMPBELL
# /// Edits --- Simon Brunner
# /// Edit date --- 11-MAY-2018

#############################################
## -- Receive command line args --
args = commandArgs(trailingOnly=TRUE)
args

# Check for command line args, else set defaults
if (length(args)==0) {
  fpath_rsession = ''
} else if (length(args)==1) {
  fpath_rsession = args[1]
} else {
  stop("Script takes either exactly 0 or 1 arguments.", call.=FALSE)
}

# For debugging
if(FALSE) {
  fpath_rsession = ''
}

# Report arguments
message(sprintf('Continuing R session stored at %s', fpath_rsession))

#############################################
## -- Library loading and sourcing --

load(fpath_rsession)
library(RColorBrewer)
library(label.switching)
library(philentropy)
library(ggplot2)
library(GGally)
library(dplyr)
library(data.table)

source("../../../../../lib/Mutation_Cluster_Dirichlet_Gibbs_sampler_functions.R")
source('../../../../../lib/Mutation_Cluster_Dirichlet_posthoc_merge_fx.R')

#############################################

# Merge clusters which differ only at level of sequencing errors

lymph.merge.ls <- merge.clusters.differ.by.seq.error(lymph.gs, lymph.ls, min.true.VAF = 0.025, overlap.frac=0.025, min.num.muts = 20, ls.type="ECR", number.burn.in, version=2)
ls_tbl = tibble(clust_assign = as.numeric(lymph.merge.ls$final.cluster)) %>%
  mutate(mut_id = row_number())
ls_tbl$chrom = lymph.depth$chrom
ls_tbl$pos = lymph.depth$pos
fwrite(ls_tbl, file.path(output_dir, 'clust_assign_posthoc.csv'))

#############################################

# Final spectrum and cluster location plots
pdf(file.path(output_dir, "Cluster_and_spectrum_plots.pdf"), paper="a4", width = 8, height=10)

par_min.threshold = 100
post.cluster.pos <- post.param.dist_posthoc(gs.out = lymph.gs, ls.merge.out = lymph.merge.ls, 
                                            centiles = c(0.025, 0.5, 0.975), ls.type = "ECR", 
                                            burn.in = number.burn.in, samp.names = short.names, 
                                            plot.heatmap = TRUE, min.threshold = par_min.threshold)  

mut.spec.by.clust <- context.extract_posthoc(mut.spec, lymph.merge.ls)

dev.off()

## Write heatmap data to disk
heatmap_tbl = tbl_df(post.cluster.pos$heat.dat*2)
#heatmap_tbl$cluster_id = rownames(post.cluster.pos$heat.dat)
heatmap_tbl$estimated.no.of.mutations = as.numeric(post.cluster.pos$num.mut[which(post.cluster.pos$num.muts>=par_min.threshold)])
heatmap_tbl$no.of.mutations.assigned = as.numeric(post.cluster.pos$num.mut[which(post.cluster.pos$num.muts>=par_min.threshold)])

# Add a root node 
heatmap_tbl = rbind(heatmap_tbl,
                    c(rep(1, dim(post.cluster.pos$heat.dat)[2]), sum(heatmap_tbl$no.of.mutations.assigned), sum(heatmap_tbl$no.of.mutations.assigned)))

fwrite(heatmap_tbl, file.path(output_dir, 'heatmap.tsv'), sep='\t')

## Export a table with mutations per cluster
mut_per_cluster = tibble(cluster_id = names(post.cluster.pos$num.mut), num_mut = post.cluster.pos$num.mut)
fwrite(mut_per_cluster, file.path(output_dir, 'muts_per_cluster.csv'))

#############################################

# Density plots of clusters with at least 100 mutations
pdf(file.path("Raw_mutation_VAF_by_cluster_plots.pdf"), width = 8, height = 10)

which.clusters <- as.double(names(which(table(lymph.merge.ls[["final.cluster"]]) >= 100)))

for (i in which.clusters) {
  temp.plot <- cluster.density.plot_posthoc(lymph.gs, lymph.merge.ls, post.cluster.pos, i, samp.names = short.names)
  print(temp.plot)
}

dev.off()


#save.image()
