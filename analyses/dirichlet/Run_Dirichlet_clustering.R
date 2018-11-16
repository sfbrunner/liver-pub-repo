# Run Dirichlet clustering
# --
# Runs Dirichlet clustering and label switching on low input mutation data
#
# Args (command line):
#   1: Path to CSV containing counts of alt base (rows are variants, any columns starting with 'PD' are samples)
#   2: Path to CSV containing depth (rows are variants, any columns starting with 'PD' are samples)
#   3: Path to table with mutation contexts (generated with pc8's "Context_pull_build37.pl" script)
#   4: Directory to output results
#   5: Number of burnin cycles, eg 10000
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
  fpath_alt_csv = '/nfs/users/nfs_s/sb50/mylustre/clonal_heterogeneity/donors/PD37113b/dirichlet/input/target_alt.csv'
  fpath_depth_csv = '/nfs/users/nfs_s/sb50/mylustre/clonal_heterogeneity/donors/PD37113b/dirichlet/input/target_depth.csv'
  fpath_mutant_locations = '/nfs/users/nfs_s/sb50/mylustre/clonal_heterogeneity/donors/PD37113b/dirichlet/input/'
  output_dir = '/nfs/users/nfs_s/sb50/mylustre/clonal_heterogeneity/donors/PD37113b/dirichlet'
  num_burnin = 10000
} else if (length(args)==5) {
  fpath_alt_csv = args[1]
  fpath_depth_csv = args[2]
  fpath_mutant_locations = args[3]
  output_dir = args[4]
  num_burnin = as.numeric(args[5])
} else {
  stop("Script takes either exactly 0 or 5 arguments.", call.=FALSE)
}

# For debugging
if(FALSE) {
  fpath_alt_csv = '/nfs/users/nfs_s/sb50/mylustre/clonal_heterogeneity/donors/PD37113b/dirichlet/input/target_alt.csv'
  fpath_depth_csv = '/nfs/users/nfs_s/sb50/mylustre/clonal_heterogeneity/donors/PD37113b/dirichlet/input/target_depth.csv'
  fpath_mutant_locations = '/nfs/users/nfs_s/sb50/mylustre/clonal_heterogeneity/donors/PD37113b/dirichlet/input/mut_contexts.txt'
  output_dir = '/nfs/users/nfs_s/sb50/mylustre/clonal_heterogeneity/donors/PD37113b/dirichlet'
}

# Report arguments
message(sprintf('Reading alt counts from %s', fpath_alt_csv))
message(sprintf('Reading depth from %s', fpath_depth_csv))
message(sprintf('Reading mutation contexts from %s', fpath_mutant_locations))
message(sprintf('Writing output to %s', output_dir))
message(sprintf('Number of burnin cycles is %s', num_burnin))

#############################################
## -- Library loading and sourcing --

library(RColorBrewer)
library(label.switching)
library(philentropy)
library(ggplot2)
library(GGally)
library(dplyr)
library(data.table)

source("../../../../../lib/Mutation_Cluster_Dirichlet_Gibbs_sampler_functions.R")

today = format(Sys.Date(), "%Y_%m_%d")
this_hash = paste(sample(letters,5,TRUE),collapse="")
output_dir = file.path(output_dir, paste('ndp', today, this_hash, sep="_"))

dir.create(output_dir, showWarnings = F)

# Define params
number.iterations <- 15000
number.burn.in <- 10000
#number.burn.in <- num_burnin

set.seed(28)

#############################################
## -- Run clustering --

lymph.var <- read.table(fpath_alt_csv, sep=",", header=TRUE, stringsAsFactors = FALSE)
lymph.var <- lymph.var[lymph.var$chrom != "X" & lymph.var$chrom != "Y",] # Remove X & Y chrom muts
lymph.y <- as.matrix(lymph.var[,which(grepl('PD', names(lymph.var)))])

lymph.depth <- read.table(fpath_depth_csv, sep=",", header=TRUE, stringsAsFactors = FALSE)
lymph.depth <- lymph.depth[lymph.depth$chrom != "X" & lymph.depth$chrom != "Y",]
lymph.N <- as.matrix(lymph.depth[,which(grepl('PD', names(lymph.var)))])

lymph.gs <- muts.cluster.dirichlet.gibbs(C = 100, y = lymph.y, N = lymph.N, iter=number.iterations,
                                         iter.chunk=100, outdir=output_dir)

mut.spec <- read.table(fpath_mutant_locations, sep="\t", header=TRUE, 
                       stringsAsFactors = FALSE)
mut.spec <- mut.spec[mut.spec$chrom != "X" & mut.spec$chrom != "Y",]

short.names <- names(lymph.var)[which(grepl('PD', names(lymph.var)))]
short.names[2:length(short.names)] <- substring(short.names[2:length(short.names)], 10)

save.image(file.path(output_dir, 'Rsession.dat'))

#############################################

# Generate convergence plots

pdf(file.path(output_dir, "Convergence_plots.pdf"), paper="a4", width = 8, height = 10)

for (i in c(10, 50, 100, 200, 500, seq(1000, number.iterations, 1000))) {
  heatmap.iter(lymph.gs, curr.iter = i, samp.names = short.names, min.threshold = 10)
}

plot(1:number.iterations, lymph.gs[["alpha"]], type="l", xlab="Iteration", ylab = "Alpha", main = )

par(mfrow = c(4,1))
for (i in 1:4) {
  label.switch.plot(lymph.gs, num.samples = 10)
}

dev.off()

#############################################

# Correct label switching phenomenon

lymph.ls <- correct.label.switch(lymph.gs, method.set = c("ECR"), burn.in = number.burn.in)
ls_tbl = tibble(clust_assign = as.numeric(lymph.ls$clusters)) %>%
  mutate(mut_id = row_number())
fwrite(ls_tbl, file.path(output_dir, 'clust_assign.csv'))

save.image(file.path(output_dir, 'Rsession_ls.dat'))

#############################################

# Final spectrum and cluster location plots
pdf(file.path(output_dir, "Cluster_and_spectrum_plots.pdf"), paper="a4", width = 8, height=10)

par_min.threshold = 100
post.cluster.pos <- post.param.dist(gs.out = lymph.gs, ls.out = lymph.ls, 
                                    centiles = c(0.025, 0.5, 0.975), ls.type = "ECR", 
                                    burn.in = number.burn.in, samp.names = short.names, plot.heatmap = TRUE,
                                    min.threshold = par_min.threshold)  

mut.spec.by.clust <- context.extract(mut.spec, lymph.ls)

dev.off()

## Write heatmap data to disk
heatmap_tbl = tbl_df(post.cluster.pos$heat.dat*2)
#heatmap_tbl$cluster_id = rownames(post.cluster.pos$heat.dat)
heatmap_tbl$estimated.no.of.mutations = as.numeric(post.cluster.pos$num.mut[which(post.cluster.pos$num.muts>=par_min.threshold)])
heatmap_tbl$no.of.mutations.assigned = as.numeric(post.cluster.pos$num.mut[which(post.cluster.pos$num.muts>=par_min.threshold)])

# Add a root node 
heatmap_tbl = rbind(heatmap_tbl,
                    c(rep(1, dim(post.cluster.pos$heat.dat)[2]), sum(heatmap_tbl$no.of.mutations.assigned), max(heatmap_tbl$no.of.mutations.assigned)))

fwrite(heatmap_tbl, file.path(output_dir, 'heatmap.tsv'), sep='\t')

## Export a table with mutations per cluster
mut_per_cluster = tibble(cluster_id = names(post.cluster.pos$num.mut), num_mut = post.cluster.pos$num.mut)
fwrite(mut_per_cluster, file.path(output_dir, 'muts_per_cluster.csv'))

#############################################

# Density plots of clusters with at least 100 mutations
pdf(file.path(output_dir, "Raw_mutation_VAF_by_cluster_plots.pdf"), width = 8, height = 10)

which.clusters <- as.double(names(which(table(lymph.ls$clusters) >= 100)))

for (i in which.clusters) {
  temp.plot <- cluster.density.plot(lymph.gs, lymph.ls, post.cluster.pos, i, samp.names = short.names)
  print(temp.plot)
}

dev.off()


#save.image()
