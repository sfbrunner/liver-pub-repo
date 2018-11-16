# Series of functions called to implement Dirichlet process Gibbs sampler for multi-sample VAFs

muts.cluster.dirichlet.gibbs <- function(C=30, y, N, iter=1000, A =0.01, B = 0.01, 
                                         include.split.merge = TRUE, iter.chunk = 1000, outdir = '.') {
  # y: matrix of number of reads reporting each variant: rows are mutations, columns samples
  # N: matrix of total depth each variant: rows are mutations, columns samples (same order as y)
  # C is the maximum number of clusters in the Dirichlet process
  # A and B are hyperparameters for the gamma prior on alpha 
  # iter is the number of iterations of the Gibbs sampler
  # include.split.merge is whether to include the SAMS split-merge step for clusters in the Gibbs sampler
  # iter.chunk is the number of iterations to run between writing iterations to file and resetting 
  
  num.muts <- nrow(y)
  num.samples <- ncol(y)
  N.minus.y <- N-y
  
  
  # Set up data formats for recording iterations
  pi.h.j <- array(NA, dim=c(iter, num.samples, C))
  V.h <- matrix(1, nrow=iter, ncol=C)
  Pr.S <- matrix(NA, nrow=num.muts, ncol=C)
  S.i <- matrix(NA, nrow=iter.chunk+1, ncol=num.muts)
  alpha <- rep(NA, iter)	
  merge.split <- data.frame(iter = 1:iter, first.cluster=c(1), second.cluster=c(1),
                            MH.ratio = c(0), decision = c(FALSE), 
                            new.cluster.loc = c(1), stringsAsFactors = FALSE)
  
  # Initialise
  pi.h.j[1,,] <- matrix(rep(seq(0.1, 0.9, 0.8/(C-1)), each = num.samples), byrow = TRUE)
  V.h[1,] <- c(rep(0.5,C-1), 1)
  alpha[1] <- 1
  V.h[1:iter, C] <- rep(1, iter)
  S.i[1,] <- S.i[2,] <- rep(1:C, num.muts / C + 1)[1:num.muts]
  S.i.filename <- file.path(outdir, paste0("S_i_", Sys.Date(), ".txt"))
  
  for (m in 2:iter) {
    m.S <- m - floor((m-1) / iter.chunk) * iter.chunk + 1
    
    # Step 1: Update cluster allocation for each individual mutation
    log.V.h <- log(V.h[m-1,])
    cumsum.log.1.minus.V.h <- c(0, cumsum(log(1 - V.h[m-1,1:(C-1)])))
    log.Vh.weights <- log.V.h + cumsum.log.1.minus.V.h
    log.pi.h.j <- log(pi.h.j[m-1,,])
    log.1.minus.pi <- log(1 - pi.h.j[m-1,,])
    
    y.ij.log.pi.h.j <- y %*% log.pi.h.j
    N.minus.y.log.1.minus.pi <- N.minus.y %*% log.1.minus.pi
    log.Pr.S <- matrix(rep(log.Vh.weights, each=num.muts), nrow=num.muts, byrow = FALSE) +
      y.ij.log.pi.h.j +
      N.minus.y.log.1.minus.pi
    Pr.S <- exp(log.Pr.S) 
    Pr.S[rowSums(Pr.S) == 0,] <- 1/C
    
    S.i[m.S,] <- sapply(1:num.muts, function(Pr, k) {sum(rmultinom(1,1,Pr[k,]) * (1:C))}, Pr=Pr.S)
    
    # Step 2: Update stick-breaking weights
    V.h[m,1:(C-1)]  <- sapply(1:(C-1), function(S, curr.m, curr.alpha, h) {
      rbeta(1, 1+sum(S[curr.m,] == h), curr.alpha+sum(S[curr.m,] > h))}, 
      S=S.i, curr.m=m.S, curr.alpha=alpha[m-1])
    V.h[m,c(V.h[m,1:(C-1)] == 1,FALSE)] <- 0.999 # Need to prevent one stick from taking all the remaining weight
    
    # Step 3: Update location of subclones / clusters
    curr.shape1 <- sapply(1:C, function(h) {colSums(matrix(y[S.i[m.S,] == h,], ncol=num.samples))})
    curr.shape2 <- sapply(1:C, function(h) {colSums(matrix(N[S.i[m.S,] == h,], ncol=num.samples))})
    pi.h.j[m,,] <- matrix(rbeta(n = num.samples * C,
                                shape1 = c(1 + curr.shape1),
                                shape2 = c(1 + curr.shape2 - curr.shape1)), 
                          nrow=num.samples, byrow=FALSE)
    
    # Step 4: If requested, test a split-merge step of clusters
    if (include.split.merge) {
      # Uniformly select pair of indices
      curr.pair <- sample.int(n=num.muts, size=2, replace=FALSE)
      
      # Check whether drawn pair of mutations are in same or different clusters
      if (S.i[m.S, curr.pair[1]] == S.i[m.S, curr.pair[2]]) {poss.split <- TRUE} else {poss.split <- FALSE}
      curr.cluster.1 <- S.i[m.S, curr.pair[1]]
      curr.cluster.2 <- S.i[m.S, curr.pair[2]]
      
      # Set up some baseline totals
      curr.sum.y.1 <- y[curr.pair[1],]
      curr.sum.N.1 <- N[curr.pair[1],]
      curr.sum.y.2 <- y[curr.pair[2],]
      curr.sum.N.2 <- N[curr.pair[2],]
      sum.wts.1 <- sum(log.p.beta.bin(y = curr.sum.y.1, n = curr.sum.N.1, alpha = 1, beta = 1))
      sum.wts.2 <- sum(log.p.beta.bin(y = curr.sum.y.2, n = curr.sum.N.2, alpha = 1, beta = 1))
      log.cluster.probs <- 0
      sum.all.y <- curr.sum.y.1 + curr.sum.y.2
      sum.all.N <- curr.sum.N.1 + curr.sum.N.2
      split.cl.1 <- curr.pair[1]
      split.cl.2 <- curr.pair[2]
      
      log.merge.given.y <- sum.wts.1 + sum(log.p.beta.bin(y = curr.sum.y.2, n = curr.sum.N.2, 
                                                          alpha = curr.sum.y.1 + 1, 
                                                          beta = curr.sum.N.1 - curr.sum.y.1 + 1))
      
      # Generate random permutation of other indices in the cluster(s)
      rand.indices <- which(S.i[m.S,] %in% c(curr.cluster.1, curr.cluster.2))[-curr.pair]
      if (length(rand.indices) > 0) {
        rand.indices <- sample(rand.indices, size = length(rand.indices), replace = FALSE)
        for (i in rand.indices) {
          # Work out weighting of the two clusters for allocating the next sample to 
          curr.weight.1 <- sum(log.p.beta.bin(y = y[i,], n = N[i,], 
                                              alpha = curr.sum.y.1 + 1, beta = curr.sum.N.1 - curr.sum.y.1 + 1))
          log.prob.1st.cluster <- log(length(split.cl.1)) + curr.weight.1
          
          curr.weight.2 <- sum(log.p.beta.bin(y = y[i,], n = N[i,], 
                                              alpha = curr.sum.y.2 + 1, beta = curr.sum.N.2 - curr.sum.y.2 + 1))
          log.prob.2nd.cluster <- log(length(split.cl.2)) + curr.weight.2
          
          prob.first.cluster <- 1/(1+exp(log.prob.2nd.cluster - log.prob.1st.cluster))
          
          # Determine which cluster to allocate a mutation to under a splitting scenario
          if (poss.split) {
            # Choose which cluster to allocate the next mutation to randomly
            draw <- runif(n = 1, min = 0, max = 1)
            if (draw <= prob.first.cluster) {alloc <- "first"} else {alloc <- "second"}
          } else {
            # Assign mutation to cluster from which it originated
            if (S.i[m.S, i] == curr.cluster.1) {alloc <- "first"} else {alloc <- "second"}
          }
          
          # Update the various running totals for the splitting scenario
          if (alloc == "first") {
            # Allocate to first cluster
            split.cl.1 <- c(split.cl.1, i)
            curr.sum.y.1 <- curr.sum.y.1 + y[i,]
            curr.sum.N.1 <- curr.sum.N.1 + N[i,]
            sum.wts.1 <- sum.wts.1 + curr.weight.1
            log.cluster.probs <- log.cluster.probs + log(prob.first.cluster)
          } else {
            # Allocate to second cluster
            split.cl.2 <- c(split.cl.2, i)
            curr.sum.y.2 <- curr.sum.y.2 + y[i,]
            curr.sum.N.2 <- curr.sum.N.2 + N[i,]
            sum.wts.2 <- sum.wts.2 + curr.weight.2
            log.cluster.probs <- log.cluster.probs + log(1 - prob.first.cluster)
          }
          
          # Update the various running totals for the merging scenario
          log.merge.given.y <- log.merge.given.y + 
            sum(log.p.beta.bin(y = y[i,], n = N[i,], alpha = sum.all.y + 1, 
                               beta = sum.all.N - sum.all.y + 1))
          sum.all.y <- sum.all.y + y[i,]
          sum.all.N <- sum.all.N + N[i,]
          
        }
      }
      
      # Now calculate the Metropolis-Hastings ratio
      log.split.given.y <- 2*log(alpha[m-1]) + lgamma(length(split.cl.1)) + 
        lgamma(length(split.cl.2)) + sum.wts.1 + sum.wts.2
      
      log.merge.given.y <- log.merge.given.y + log(alpha[m-1]) + lgamma(length(rand.indices)+2)
      
      if (poss.split) {
        MH.ratio <- exp(log.split.given.y - log.merge.given.y - log.cluster.probs)
      } else {
        MH.ratio <- exp(log.merge.given.y - log.split.given.y + log.cluster.probs)
      }
      
      # Check whether to split and, if so, do the necessary reassignments
      merge.split.decision <- FALSE
      if (runif(1,0,1) <= MH.ratio & poss.split) {
        # Keep the two proposed clusters
        # Find the cluster with the lowest membership
        new.cluster.loc <- as.double(names(
          which.min(table(factor(S.i[m.S,], levels = (1:C)[-curr.cluster.1])))))
        # Reassign cluster membership of one of the splits to this cluster
        S.i[m.S, split.cl.1] <- new.cluster.loc
        merge.split.decision <- TRUE
      }
      
      # Check whether to merge and, if so, do the necessary reassignments
      if (runif(1,0,1) <= MH.ratio & !poss.split) {
        # Merge the two proposed clusters
        # Find the cluster with the lowest index
        new.cluster.loc <- min(c(curr.cluster.1, curr.cluster.2))
        # Reassign cluster membership of the two original clusters to this single cluster
        S.i[m.S, split.cl.1] <- new.cluster.loc
        S.i[m.S, split.cl.2] <- new.cluster.loc
        merge.split.decision <- TRUE
      }
      
      if (merge.split.decision) {      
        # Re-draw stick-breaking weights
        V.h[m,1:(C-1)]  <- sapply(1:(C-1), function(S, curr.m, curr.alpha, h) {
          rbeta(1, 1+sum(S[curr.m,] == h), curr.alpha+sum(S[curr.m,] > h))}, 
          S=S.i, curr.m=m.S, curr.alpha=alpha[m-1])
        V.h[m,c(V.h[m,1:(C-1)] == 1,FALSE)] <- 0.999 # Prevent stick from taking remaining weight
        
        # Re-assign cluster locations
        curr.shape1 <- sapply(1:C, function(h) {colSums(matrix(y[S.i[m.S,] == h,], ncol=num.samples))})
        curr.shape2 <- sapply(1:C, function(h) {colSums(matrix(N[S.i[m.S,] == h,], ncol=num.samples))})
        pi.h.j[m,,] <- matrix(rbeta(n = num.samples * C,
                                    shape1 = c(1 + curr.shape1),
                                    shape2 = c(1 + curr.shape2 - curr.shape1)), 
                              nrow=num.samples, byrow=FALSE)
      }
      
      merge.split$iter[m] <- m
      merge.split$first.cluster[m] <- curr.cluster.1
      merge.split$second.cluster[m] <- curr.cluster.2
      merge.split$MH.ratio[m] <- MH.ratio
      merge.split$decision[m] <- merge.split.decision
      merge.split$new.cluster.loc[m] <- ifelse(test = merge.split.decision, 
                                               yes = new.cluster.loc, no = "NA")
    }
    
    # Step 5: Update alpha
    alpha[m] <- rgamma(1, shape=C+A-1, rate=B-sum(log(1-V.h[m,1:(C-1)]))) 
    
    # Check whether to write current state of S.i to file and reset
    if (m / iter.chunk == round(m/iter.chunk)) {
      print(m)
      if (m == iter.chunk) {
        write.table(S.i[2:(iter.chunk+1),], file=S.i.filename, sep="\t", quote=FALSE, 
                    append=FALSE, row.names=FALSE, col.names = TRUE)  
      } else {
        write.table(S.i[2:(iter.chunk+1),], file=S.i.filename, sep="\t", quote=FALSE, 
                    append=TRUE, row.names=FALSE, col.names = FALSE)
      }
      
      S.i[1,] <- S.i[iter.chunk+1,]
    }
    
  }
  
  return(list(V.h=V.h, pi.h.j=pi.h.j, alpha=alpha, y1=y, N1=N, merge.split = merge.split, 
              S.i.filename = S.i.filename))
} 

###########################################

log.p.beta.bin <- function(y, n, alpha, beta) {
  return(lchoose(n, y) + lbeta(y + alpha, n - y + beta) - lbeta(alpha, beta))
}

###########################################

# Some convergence plots

# Heatmap of cluster positions for a given iteration of the Gibbs sampler
heatmap.iter <- function(gs.out, curr.iter, min.threshold = 10, samp.names) {
  # curr.iter is the iteration number to plot the heatmap for
  # min.threshold is the minimum number of samples required in the cluster for inclusion in heatmap
  # samp.names is the vector of sample names to annotate the plot with
  S.i.file <- gs.out[["S.i.filename"]]
  pi.h <- gs.out[["pi.h.j"]][curr.iter,,]
  pi.h <- as.data.frame(pi.h)
  num.clusters <- ncol(pi.h)
  num.iter <- dim(gs.out[["pi.h.j"]])[1]
  names(pi.h) <- paste("Cl", 1:num.clusters, sep=".")
  row.names(pi.h) <- samp.names
  allocations <- t(read.table(file = S.i.file, header=TRUE, sep="\t", stringsAsFactors = FALSE,
                              skip = curr.iter-1, nrows=1))
  
  alloc.counts <- table(allocations)[table(allocations) >= min.threshold]
  
  temp.r <- as.matrix(pi.h[,as.double(names(alloc.counts))])
  temp.dendro <- as.dendrogram(hclust(as.dist(1 -
                                                distance(matrix(c(temp.r), ncol=ncol(temp.r)) / rowSums((temp.r)), method="cosine"))))
  heatmap(t(temp.r), Colv = temp.dendro, 
          scale="none", col=brewer.pal(9, "Oranges"),
          main = paste0("VAFs per cluster across samples: Iteration ", curr.iter),
          xlab = "Sample", ylab="Cluster", 
          RowSideColors = as.character(cut(alloc.counts, breaks = 9, 
                                           labels = brewer.pal(9,"Greens"))))
}

# Label switching
label.switch.plot <- function(gs.out, num.samples=11) {
  # Plot of cluster allocations for num.samples randomly chosen mutations
  
  S.i.file <- gs.out[["S.i.filename"]]
  merge.split <- gs.out[["merge.split"]]
  num.clusters <- dim(gs.out[["pi.h.j"]])[3]
  num.muts <- nrow(gs.out[["y1"]])
  samps <- sample.int(num.muts, size = num.samples)
  
  cols.to.read <- rep("NULL", num.muts)
  cols.to.read[samps] <- "integer"
  allocations <- read.table(S.i.file, colClasses = cols.to.read, header = TRUE, sep="\t")
  
  plot(1:nrow(allocations), allocations[,1] + runif(n = 1, min=-0.2, max=0.2), type="l", 
       col=brewer.pal(11, "BrBG")[1], ylim=c(-3,num.clusters), xlab="Iteration", 
       ylab="Cluster assignment")
  for (i in 2:num.samples) {
    lines(1:nrow(allocations), allocations[,i] + runif(n = 1, min=-0.2, max=0.2), 
          type="l", col=brewer.pal(11, "BrBG")[i])
  }
  points(merge.split$iter[merge.split$decision], 
         ifelse(merge.split$first.cluster[merge.split$decision] == merge.split$second.cluster[merge.split$decision], -1, -2), 
         pch=20, col=ifelse(merge.split$first.cluster[merge.split$decision] == merge.split$second.cluster[merge.split$decision], "plum4", "springgreen4"))
}

###########################################

# Correct label switching events

correct.label.switch <- function(gs.out, method.set, burn.in = 2000) {
  # gs.out is the output from the muts.cluster.dirichlet.gibbs function
  # burn.in is the number of iterations to drop
  # S.i.file is the file-name for file containing mutation allocations to clusters
  
  S.i.file <- gs.out[["S.i.filename"]]
  #S.i <- read.table(file = S.i.file, header=TRUE, sep="\t", stringsAsFactors = FALSE,
  #                  skip = burn.in)
  S.i <- fread(file = S.i.file, header=TRUE, skip = burn.in)
  S.i = as.data.frame(S.i)
  #randcols = sample(1:dim(S.i)[2], 1000)
  #S.i = S.i[,randcols]
  pi.h.j <- gs.out[["pi.h.j"]]
  V.h <- gs.out[["V.h"]]
  y <- gs.out[["y1"]]
  N <- gs.out[["N1"]]
  num.iter <- nrow(V.h)
  num.muts <- ncol(S.i)
  num.clusters <- ncol(V.h)
  num.samples <- dim(pi.h.j)[2]
  Pr.S <- array(NA, dim = c(num.iter-burn.in, num.muts, num.clusters))
  gs.params <- array(NA, dim = c(num.iter-burn.in, num.clusters, num.samples + 1))
  
  # Regenerate allocation probabilities
  for (m in (burn.in+1):num.iter) {
    log.V.h <- log(V.h[m-1,])
    cumsum.log.1.minus.V.h <- c(0, cumsum(log(1 - V.h[m-1,1:(num.clusters-1)])))
    log.Vh.weights <- log.V.h + cumsum.log.1.minus.V.h
    log.pi.h.j <- log(pi.h.j[m-1,,])
    log.1.minus.pi <- log(1 - pi.h.j[m-1,,])
    
    if (sum(method.set %in% c("STEPHENS", "ECR-ITERATIVE-2")) > 0) {
      y.ij.log.pi.h.j <- y %*% log.pi.h.j
      N.minus.y.log.1.minus.pi <- (N - y) %*% log.1.minus.pi
      log.Pr.S <- matrix(rep(log.Vh.weights, each=num.muts), nrow=num.muts, byrow = FALSE) +
        y.ij.log.pi.h.j +
        N.minus.y.log.1.minus.pi
      Pr.S[m - burn.in,,] <- exp(log.Pr.S) 
      Pr.S[m - burn.in,rowSums(Pr.S[m-burn.in,,]) == 0,] <- 1/num.clusters
    }
    gs.params[m - burn.in,,1:num.samples] <- t(pi.h.j[m,,])
    gs.params[m - burn.in,,num.samples + 1] <- log.Vh.weights
  }
  
  S.i <- as.matrix(S.i)
  ls <- label.switching(method = method.set, zpivot = S.i[nrow(S.i),], z = S.i, 
                        K = num.clusters, 
                        prapivot = gs.params[num.iter-burn.in,,], p = Pr.S, constraint = 1,
                        mcmc = gs.params)
  return(ls)
}

#############################################

# Merge clusters where differences are only related to sequencing errors

merge.clusters.differ.by.seq.error <- function(gs.out, ls.out, min.true.VAF = 0.025, overlap.frac=0.025, min.num.muts = 20, ls.type="ECR", burn.in, version=2) {
  # Function to merge clusters that are only separated from one another by seq error rates
  # Need for the function is that some clusters are falsely split by the level of seq error
  ## in samples where the mutations are absent (ie, some base positions noisier than others)
  
  # Two versions of this algorithm:
  # Version 1 = For 2 clusters, looks to see whether iterations of MCMC alternate sign for each non-zero sample
  # Version 2 = For 2 clusters, looks to see whether mean locations (pi.h.j) for each non-zero sample differ by less that a certain threshold 
  
  # gs.out is the output from the Gibbs sampler
  # ls.out is the output from the label.switching() function
  # min.true.VAF is the VAF threshold at which a sample will be counted as carrying the mutations
  # overlap.frac (VERSION 1) is the minimum fraction of MCMC estimates of pi that overlap in the two clusters for each non-zero sample required for merging - increase its value for less merging; decrease its value for more merging 
  # overlap.frac (VERSION 2) is the maximum difference between pi.h.j estimates for each non-zero sample
  # min.num.muts is the minimum number of mutations in the clusters to worry about merging them
  
  ls.final <- ls.out$clusters[1,]
  unique.clusters <- table(ls.final)
  kept_clusters = which(unique.clusters >= min.num.muts)
  unique.clusters <- unique.clusters[kept_clusters]
  unique.clusters <- names(unique.clusters)[order(unique.clusters, decreasing = TRUE)]
  
  y <- gs.out[["y1"]]
  N <- gs.out[["N1"]]
  merge.spl <- gs.out[["merge.split"]]
  
  ls.perm <- ls.out$permutations[ls.type][[1]]
  pi.h.j <- gs.out[["pi.h.j"]]
  num.clusters <- dim(pi.h.j)[3]
  num.samples <- dim(pi.h.j)[2]
  num.iter <- dim(pi.h.j)[1]
  
  pi.mcmc <- array(NA, dim = c(num.iter - burn.in, num.clusters, num.samples))
  for (i in (burn.in+1):num.iter) {pi.mcmc[i-burn.in,,] <- t(pi.h.j[i,,])}
  # Permute according to label-switching() defined optimal labelling 
  pi.mcmc <- permute.mcmc(mcmc = pi.mcmc, permutations = ls.perm)$output
  
  # For successful merge-split steps after burn.in, turn estimates of pi before the merge-split to NA
  merge.spl <- merge.spl[(burn.in+1):num.iter,]
  merge.spl <- merge.spl[merge.spl$decision,]
  merge.spl$new.cluster.loc <- as.numeric(merge.spl$new.cluster.loc)
  for (i in 1:nrow(merge.spl)) {
    final.cluster.1 <- ls.perm[merge.spl$iter[i] - burn.in, merge.spl$first.cluster[i]]
    final.cluster.2 <- ls.perm[merge.spl$iter[i] - burn.in, merge.spl$second.cluster[i]]
    final.cluster.3 <- ls.perm[merge.spl$iter[i] - burn.in, merge.spl$new.cluster.loc[i]]
    
    if (merge.spl$first.cluster[i] == merge.spl$second.cluster[i]) {
      pi.mcmc[1:(merge.spl$iter[i] - burn.in), final.cluster.1, ] <- NA
      pi.mcmc[1:(merge.spl$iter[i] - burn.in), final.cluster.3, ] <- NA
    } else {
      pi.mcmc[1:(merge.spl$iter[i] - burn.in), final.cluster.1, ] <- NA
      pi.mcmc[1:(merge.spl$iter[i] - burn.in), final.cluster.2, ] <- NA
    }
  }
  
  mean.vaf <- apply(pi.mcmc, MARGIN = 2, FUN = colMeans, na.rm=TRUE)
  if (version == 2) {
    S.i <- read.table(gs.out[["S.i.filename"]], header=TRUE, sep="\t", stringsAsFactors = FALSE,
                      skip = burn.in)
    size.clusters <- t(apply(S.i, MARGIN = 1, 
                             FUN = function(a) {table(factor(a, levels = 1:num.clusters))}))
  }
  
  # Check whether to merge clusters
  # Merge if in all samples with VAFs above the seq error threshold, MCMC pi overlaps > overlap.frac times 
  for (i in 1:(length(unique.clusters)-1)) {
    ind.i <- as.double(unique.clusters[i])
    for (j in (i+1):length(unique.clusters)) {
      ind.j <- as.double(unique.clusters[j])
      nonzero.samps <- which(mean.vaf[,ind.i] >= min.true.VAF | mean.vaf[,ind.j] >= min.true.VAF )
      num.overlaps <- sapply(1:dim(pi.mcmc)[3], function(a,k,x,y) {
        temp <- na.exclude(a[,x,k] > a[,y,k]);
        if (length(temp)>0) {sum(temp) / length(temp)} else {0}}, 
        a=pi.mcmc, x=ind.i, y=ind.j)
      if (version == 1 & 
          all(num.overlaps[nonzero.samps] > overlap.frac & 
              num.overlaps[nonzero.samps] < 1-overlap.frac)) {
        # Clusters should be merged
        # Merge them into the larger cluster (ind.i)
        ls.final[ls.final == ind.j] <- ind.i
      }
      if (version == 2 & 
          all(abs(mean.vaf[,ind.i] - mean.vaf[,ind.j]) < overlap.frac)) {
        # Clusters should be merged - weight the pi estimates by size of cluster
        ls.final[ls.final == ind.j] <- ind.i
        pi.mcmc[ , ind.i, ] <- (size.clusters[, ind.i] * pi.mcmc[ , ind.i, ] + 
                                  size.clusters[, ind.j] * pi.mcmc[ , ind.j, ]) /
          (size.clusters[, ind.i] + size.clusters[, ind.j])
      }
    }
  }
  return(list(final.cluster = ls.final, pi.mcmc = pi.mcmc))
}


#############################################

# Mutational spectra per cluster

context.extract <- function(muts, ls.out, plot.heatmap = TRUE, min.threshold = 100) {
  # Function to extract 96 mutation type and spectrum channels
  # muts is the data-frame of mutations
  # ls.out is the output from the label.switching() function
  # min.threshold is the minimum number of mutations in cluster for including in heatmap
  
  muts$TYPE <- rep("C>T", nrow(muts))
  muts$TYPE[(muts$ref=="C" & muts$alt == "G") | (muts$ref == "G" & muts$alt == "C")] <- "C>G"
  muts$TYPE[(muts$ref=="C" & muts$alt == "A") | (muts$ref == "G" & muts$alt == "T")] <- "C>A"
  muts$TYPE[(muts$ref=="T" & muts$alt == "A") | (muts$ref == "A" & muts$alt == "T")] <- "T>A"
  muts$TYPE[(muts$ref=="T" & muts$alt == "C") | (muts$ref == "A" & muts$alt == "G")] <- "T>C"
  muts$TYPE[(muts$ref=="T" & muts$alt == "G") | (muts$ref == "A" & muts$alt == "C")] <- "T>G"
  muts$CONTEXT[muts$ref == "G" | muts$ref == "A"] <- 
    sapply(muts$CONTEXT[muts$ref == "G" | muts$ref == "A"], rev.comp)
  muts$spec96 <- paste0(muts$TYPE, ".", substring(muts$CONTEXT, first = 10, last = 12))
  
  muts$cluster <- ls.out$clusters[1,]
  cluster.by.spec <- table(muts$cluster, muts$spec96)
  if (plot.heatmap) {
    temp.clust.spec <- cluster.by.spec[rowSums(cluster.by.spec) >= min.threshold,]
    row.names(temp.clust.spec) <- paste0("Cl.", row.names(temp.clust.spec))
    temp.dendro <- as.dendrogram(hclust(as.dist(1 -
                                                  distance(matrix(c(temp.clust.spec), ncol=96) / rowSums((temp.clust.spec)), method="cosine"))))
    
    heatmap(temp.clust.spec, scale="row",  
            col=brewer.pal(9, "PuBuGn"), Colv=NA, Rowv = temp.dendro, ylab="Cluster",
            xlab = "Mutation type", main = "96 channel mutation spectrum by cluster",
            ColSideColors = rep(brewer.pal(6, "BrBG"), each=16),
            RowSideColors = as.character(cut(rowSums(temp.clust.spec), breaks = 9, 
                                             labels = brewer.pal(9,"Greens"))))
  }
  
  return(muts)
  
}

rev.comp <- function(x) {
  x <- toupper(x)
  x <- rev(unlist(strsplit(x,"")))
  x <- paste(sapply(x, switch,  "A"="T", "T"="A","G"="C","C"="G"),collapse="")
  return(x)
}

#############################################

# Posterior distribution of cluster locations after resolving label switching

post.param.dist <- function(gs.out, ls.out, centiles = c(0.025, 0.5, 0.975), ls.type = "ECR", burn.in = 2000, samp.names, plot.heatmap = TRUE, min.threshold = 100) {
  # Function to return centiles of posterior distribution for parameters
  # gs.out is output from Gibbs sampler
  # ls.out is output from label switching correction
  # ls.type is the preferred method for label switching correction to be used
  # burn.in is the number of iterations of burn-in supplied to the label.switching() function
  # min.threshold is the minimum number of mutations in a cluster for plotting in heatmap
  
  ls.perm <- ls.out$permutations[ls.type][[1]]
  pi.h.j <- gs.out[["pi.h.j"]]
  num.clusters <- dim(pi.h.j)[3]
  num.samples <- dim(pi.h.j)[2]
  num.iter <- dim(pi.h.j)[1]
  
  pi.mcmc <- array(NA, dim = c(num.iter - burn.in, num.clusters, num.samples))
  for (i in (burn.in+1):num.iter) {pi.mcmc[i-burn.in,,] <- t(pi.h.j[i,,])}
  
  pi.mcmc <- permute.mcmc(mcmc = pi.mcmc, permutations = ls.perm)$output
  
  result <- array(NA, dim=c(num.clusters, num.samples, length(centiles)),
                  dimnames = list(clusters = paste0("Cl.", 1:num.clusters),
                                  samples = samp.names,
                                  centiles = paste0("Centile.", centiles)))
  
  for (i in 1:num.clusters) {
    for (j in 1:num.samples) {
      result[i,j,] <- quantile(pi.mcmc[,i,j], probs = centiles)
    }
  }
  
  num.muts <- table(ls.out$clusters)
  
  if (plot.heatmap) {
    temp.r <- as.matrix(result[as.double(names(num.muts))[num.muts >=min.threshold],,"Centile.0.5"])
    temp.dendro <- as.dendrogram(hclust(as.dist(1 -
                                                  distance(matrix(c(temp.r), ncol=num.samples) / rowSums((temp.r)), method="cosine"))))
    
    heatmap(temp.r, Rowv = temp.dendro, xlab="Sample", ylab="Cluster", 
            main="Posterior median VAFs by cluster and sample",
            scale="none", col=brewer.pal(9, "Oranges"),
            RowSideColors = as.character(cut(
              num.muts[num.muts >= min.threshold], breaks = 9, 
              labels = brewer.pal(9,"Greens"))))
  }
  
  return(list(centiles = result, num.muts = num.muts, pi.mcmc = pi.mcmc, heat.dat = temp.r))
}

#############################################

# Density plots for cluster locations

cluster.density.plot <- function(gs.out, ls.out, post.ls.dist, curr.cluster, vaf.threshold = 0.05, 
                                 samp.names) {
  y1 <- gs.out[["y1"]]
  N1 <- gs.out[["N1"]]
  vaf <- y1/N1
  colnames(vaf) <- samp.names
  which.muts <- which(ls.out$clusters == curr.cluster)
  
  pi.mcmc <- post.ls.dist[["pi.mcmc"]][,curr.cluster,]
  colnames(pi.mcmc) <- samp.names
  
  # Find which samples have median VAF greater than threshold; and add one random sample not part of cluster
  which.samples <- which(apply(pi.mcmc, MARGIN=2, FUN = median) >= vaf.threshold)
  if (length(which.samples) == 0) {which.samples <- sample(1:ncol(pi.mcmc), size=5)}
  if (length(which.samples) < ncol(pi.mcmc)) {
    which.samples <- c(which.samples, sample((1:ncol(pi.mcmc))[-which.samples], size = 1))
  }
  
  pi.df <- pi.mcmc[, which.samples]
  pi.df <- cbind(pi.df, data.frame(pi.est = factor(rep(20, nrow(pi.df))), raw = factor(rep(NA, nrow(pi.df)))))
  vaf.df <- vaf[which.muts, which.samples]
  vaf.df <- cbind(vaf.df, data.frame(pi.est = factor(rep(NA, nrow(vaf.df))), raw = factor(rep(20, nrow(vaf.df)))))
  combined.df <- rbind(pi.df, vaf.df)
  
  prs.plt <- ggpairs(vaf.df, columns = 1:(ncol(combined.df)-2), 
                     title = paste0("Cluster ", curr.cluster),
                     lower = list(continuous = wrap("points", alpha = 0.3, size=0.25, col="plum4")),
                     upper = list(continuous = wrap("density"))) +
    theme(panel.grid.minor = element_blank()) 
  pm2 <- prs.plt
  for(i in 2:prs.plt$nrow) {
    for(j in 1:(i-1)) {
      pm2[i,j] <- prs.plt[i,j] +
        scale_x_continuous(limits = c(0, 0.75)) +
        scale_y_continuous(limits = c(0, 0.75))
    }
  }
  
  return(pm2)
}





