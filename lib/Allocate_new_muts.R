allocate.new.muts.to.cluster <- function(ls.merge.out, post.ls.dist, new.y, new.N, uniform_weights=FALSE) {
  # Function to allocate a new set of mutations to an existing DP clustering based on VAFs
  # Only allocates to non-zero clusters
  # Uses a naive Bayesian classifier - weighted by number of muts in cluster and VAF means
  # In this sense, represents a true posterior classification
  
  # ls.merge.out is the output from the DP clustering after merge.clusters.differ.by.seq.error() function
  # post.ls.dist is the posterior distribution of cluster positions from post.param.dist()
  # new.y and new.N are matrices of new mutations to be allocated to the best clusters. 
  # The columns are in the same order as the original y and N for the DP (ie samples are in same order)
  # new.y has number of reads reporting the variant; new.N has total depth
  
  # Returns vector of best allocations for the new mutations
  
  # Extract variables and data
  ls.orig <- ls.merge.out[["final.cluster"]]
  post.pi <- post.ls.dist[["centiles"]][,,"Centile.0.5"]
  post.pi <- t(post.pi)
  new.N.minus.y <- new.N - new.y
  num.new.muts <- nrow(new.y)
  
  # Generate weights for each cluster
  cluster.counts <- table(ls.orig)
  if(uniform_weights) {
    log.weights <- rep(1, length(cluster.counts))
  } else {
    log.weights <- log(cluster.counts)
  }
  num.clusters <- length(cluster.counts)
  
  # Extract mean VAFs for each non-zero cluster
  post.pi <- post.pi[, colnames(post.pi) %in% paste0("Cl.",names(cluster.counts))]
  log.post.pi <- log(post.pi)
  log.1.minus.pi <- log(1 - post.pi)
  
  # Allocate new mutations to clusters
  y.ij.log.post.pi <- new.y %*% log.post.pi
  N.minus.y.log.1.minus.pi <- new.N.minus.y %*% log.1.minus.pi
  log.Pr.S <- matrix(rep(log.weights, each=num.new.muts), nrow=num.new.muts, byrow = FALSE) +
    y.ij.log.post.pi +
    N.minus.y.log.1.minus.pi
  Pr.S <- exp(log.Pr.S) 
  Pr.S[rowSums(Pr.S) == 0,] <- 1/num.clusters
  
  S.i <- sapply(1:num.new.muts, function(k) {sum(rmultinom(1,1,Pr.S[k,]) * (1:num.clusters))})
  S.i <- as.double(names(cluster.counts)[S.i])
  return(S.i)
}

naive_bayes <- function(cluster_prob, muts_per_cluster, new.y, new.N) {
  # Assigns variants to highest probability cluster using Naive Bayes classifier relying on betabinomial model
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 26-JUN-2018
  #
  # Args
  #   cluster_prob: mxn matrix of m clusters and n samples, giving median VAF for each cluster-sample pair
  #   muts_per_cluster: vector of length m giving the number of variants previously assigned to each cluster
  #   new.y: table with counts per variant and sample. Each row is a variant. 
  #     Columns with names that overlap with those in cluster_prob are interpreted as sample columns.
  #   new.N: table with depth per variant and sample. Each row is a variant. 
  #     Columns with names that overlap with those in cluster_prob are interpreted as sample columns.
  # Returns:
  #   tibble where rows are variants and column "clust" gives the assigned cluster, "prob" the assignment probability
  
  # Identify overlapping samples
  inc_samp = intersect(colnames(new.y), colnames(cluster_prob))
  cluster_prob = cluster_prob[,inc_samp]
  mat_y = as.matrix(new.y[,inc_samp])
  mat_N = as.matrix(new.N[,inc_samp])
  
  # Calculate cluster probability
  p_cluster_log = log(muts_per_cluster)-log(sum(muts_per_cluster))
  names(p_cluster_log) = names(muts_per_cluster)
  # For each cluster, calculate probability for each variant using betabinom dist
  num_cluster = dim(cluster_prob)[1]
  num_samples = length(inc_samp)
  variant_prob = sapply(1:num_cluster, function(x) {
    rowSums(sapply(1:num_samples, function(s) {
      log(dbinom(mat_y[,inc_samp[s]], mat_N[,inc_samp[s]], cluster_prob[x,s]))
    }))
  })
  
  # Sum log of variant and cluster probability, take to exponent
  variant_prob_weighted = (matrix(rep(p_cluster_log, times=dim(variant_prob)[1]), ncol=length(p_cluster_log)) + variant_prob)
  
  # Find maximum probability cluster for each variant
  best_clust = sapply(1:dim(variant_prob_weighted)[1], function(x) {
    names(p_cluster_log)[which.max(variant_prob_weighted[x,])]
  })
  highest_prob = sapply(1:dim(variant_prob_weighted)[1], function(x) {
    variant_prob_weighted[x,which.max(variant_prob_weighted[x,])]
  })
  
  
  # Construct output table
  out_tbl = tibble(
    chrom = new.y$chrom,
    pos = new.y$pos,
    ref = new.y$ref,
    alt = new.y$alt,
    best_clust = best_clust,
    highest_prob = highest_prob
  )
  
  return(out_tbl)
  
}


