# phylogeny_functions
# --
# This set of functions is useful for building and drawing phylogenetic trees.
# --
# /// Author --- SIMON FELIX BRUNNER
# /// Creation date --- 04-JUN-2018
library(ggrepel)
library(ape)
#library(sf)
#library(sp)

convert_edges_to_phylo <- function(tree_tbl, muts_per_cluster) {
  # Takes a table specifying tree edges and generates an object of type phylo (Ape tree)
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 04-JUN-2018
  #
  # Args
  #   tree_tbl: tibble with two columns: Edge.A and Edge.B. Each row specifies the connection
  #             from one node to another. The root node needs to be identified as a PD id (needs 'PD' in the name)
  #   muts_per_cluster: tibble with two columns: cluster_id and num_mut 
  #             (cluster IDs in muts_per_cluster and tree_tbl must be consistent)
  #   
  # Returns:
  #   tr: object of type phylo (Ape tree)
  # )
  # In the output of the function there is also a call to Ape's checkValidPhylo function. Please note that
  # if there are any FATAL errors then the tree is corrupt and cannot be plotted. In fact, plotting will
  # kill the current R session.
  
  # 1 --- Define nodes ---
  
  # Determine which nodes are tips (nodes that are only referenced once in edges)
  tip_labels = table(c(tree_tbl$Edge.A, tree_tbl$Edge.B))
  tip_labels = names(tip_labels)[tip_labels==1]
  
  # Construct table of tips
  tree_tips = tibble(old_label=tip_labels) %>%
    dplyr::mutate(new_label=dplyr::row_number(), node_type = 'tip')
  tip_labels_new = as.character(tree_tips$new_label) # For later use.
  num_tips = length(tip_labels)
  message(sprintf('Detected %s tree tips', num_tips))
  
  # Define root node
  donor_id = unique(tree_tbl$Edge.A[grepl('PD', tree_tbl$Edge.A)])
  root_node = tibble(old_label=donor_id, new_label=dim(tree_tips)[1]+1, node_type='root')
  message(sprintf('Detected root node %s', donor_id))
  
  # Define internal nodes
  all_nodes = unique(c(tree_tbl$Edge.A, tree_tbl$Edge.B))
  non_internal = unique(c(tip_labels, donor_id))
  internal_nodes = tibble(old_label=setdiff(all_nodes, non_internal)) %>%
    filter(old_label!=root_node$old_label[1]) %>%
    filter(old_label!=0) %>%
    dplyr::mutate(new_label = dim(tree_tips)[1]+1+row_number(), node_type='internal')
  message(sprintf('Detected %s internal nodes (non-root, non-tip)', dim(internal_nodes)[1]))
  
  # Build node table
  node_tbl = do.call('rbind', list(tree_tips, root_node, internal_nodes))
  
  # 2 --- Define edges ---
  
  # Build edge table
  edge_tbl = tree_tbl %>%
    dplyr::rename('old_child'=Edge.B, 'old_parent'=Edge.A) %>%
    mutate(old_child = as.character(old_child), old_parent = as.character(old_parent)) %>%
    left_join(node_tbl %>% dplyr::select(-node_type), by=c('old_child'='old_label')) %>%
    dplyr::rename('new_child'=new_label) %>%
    left_join(node_tbl %>% dplyr::select(-node_type), by=c('old_parent'='old_label')) %>%
    dplyr::rename('new_parent'=new_label)
  
  # 3 --- Generate full table containing all information about nodes, edges and labels  ---
  
  # Generate Ape table
  ape_tbl = edge_tbl %>%
    dplyr::rename('parent'=new_parent, 'child'=new_child) %>%
    dplyr::select(parent, child)
  # Correct cases where first edge column has tips instead of nodes
  ape_tbl = ape_tbl %>%
    mutate(new_parent = ifelse(parent<num_tips, child, parent)) %>%
    mutate(new_child = ifelse(parent<num_tips, parent, child)) %>%
    dplyr::select(-parent, -child) %>%
    dplyr::rename('parent'=new_parent, 'child'=new_child)
  # Make sure root is not in second edge column
  ape_tbl = ape_tbl %>%
    mutate(new_parent = ifelse(child==root_node$new_label, child, parent)) %>%
    mutate(new_child = ifelse(child==root_node$new_label, parent, child)) %>%
    dplyr::select(-parent, -child) %>%
    dplyr::rename('parent'=new_parent, 'child'=new_child) %>%
    dplyr::select(parent, child)
  
  # Merge cluster IDs and number of variants
  muts_per_new_cluster = left_join(
    node_tbl, muts_per_cluster %>% mutate(cluster_id = as.character(cluster_id)), by=c('old_label'='cluster_id')
  ) %>%
    mutate(old_label = paste('Cl', old_label, sep='.'))
  
  ape_tbl2 = ape_tbl %>%
    left_join(muts_per_new_cluster, by=c('child'='new_label'))
  
  num_nodes = as.integer(dim(internal_nodes)[1]+1)
  tr <- list(edge = (as.matrix(ape_tbl2[,1:2])), 
             tip.label = unlist(muts_per_new_cluster[1:length(tip_labels_new),1]), 
             Nnode = num_nodes,
             edge.length = ape_tbl2$num_mut,
             node.label = filter(muts_per_new_cluster, node_type!='tip')$old_label
  )
  class(tr) <- "phylo"
  checkValidPhylo(tr)
  
  return(tr)
}

draw_nice_tree <- function(tr, show_internal_nodes=FALSE, edge_strength=NA, ...) {
  if(is.na(edge_strength)) {
    plot(tr, ...)#, show.node.label=TRUE)#, direction='downwards')
  } else {
    col_weak = '#303030'# '#c1c1c1'
    col_strong = '#000000'
    edge_col = as.character(factor(edge_strength, levels=c('s', 'w'), labels=c(col_strong, col_weak)))
    edge_lty = as.numeric(as.character(factor(edge_strength, levels=c('s', 'w'), labels=c(1, 3))))
    plot(tr, edge.color=edge_col, edge.lty=edge_lty, ...)
  }
  axisPhylo(1, root.time=0, backward=FALSE)
  title(xlab = 'Number of mutations', cex=0.5, mgp=c(2.5,0,0))
  highlight_tree_nodes(tr, pch=23, bg='red', show_internal_nodes=show_internal_nodes)
  
  
}

draw_all_edge_sigs <- function(tr, sig_contrib_tbl, lwidth=0.2, sig_cols = NA, sig_order = NA) {
  # Draw signature contributions onto all branches in tree
  # (call after plot.phylo(tr))
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 04-JUN-2018
  #
  # Args
  #   tr: object of type phylo
  #   sig_contrib_tbl: table where each row represents a cluster (branch). Columns named
  #     'sig1', 'sig2', ... , 'sigN' store the fraction (between 0 and 1) that
  #     each signature contributes to the cluster (branch).
  #   lwidth: width of line to draw
  #   sig_cols: vector with colors for each signature. If left NA, then default colors are used.
  #   sig_order: vector containing signature names in a specified order
  #   
  # Returns: nothing
  
  par(xpd=TRUE, mar=c(12.1, 4.1, 4.1, 2.1))
  for(i in c(1:dim(tr$edge)[1])) {
    draw_edge_sigs(tr = tr, edge_id = i, sig_contrib_tbl = sig_contrib_tbl, lwidth=lwidth, sig_cols=sig_cols, sig_order=sig_order,
                   draw_legend = ifelse(i==1, TRUE, FALSE))
  }
}

draw_edge_sigs <- function(tr, edge_id, sig_contrib_tbl, lwidth=0.2, sig_cols = NA, sig_order = NA, draw_legend=FALSE) {
  # Draw signature contributions onto a specific branch in tree
  # (call after plot.phylo(tr))
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 04-JUN-2018
  #
  # Args
  #   tr: object of type phylo
  #   edge_id: id of edge to draw
  #   sig_contrib_tbl: table where each row represents a cluster (branch). Columns named
  #     'sig1', 'sig2', ... , 'sigN' store the fraction (between 0 and 1) that
  #     each signature contributes to the cluster (branch).
  #   lwidth: width of line to draw
  #   sig_cols: vector with colors for each signature. If left NA, then default colors are used.
  #   sig_order: vector containing signature names in a specified order
  #   
  # Returns: nothing
  
  # Get current edge
  this_edge = tr$edge[edge_id,]
  
  # Get parent coords
  node_depths = node.depth.edgelength(tr)
  num_tips = length(tr$tip.label)
  parent_node = this_edge[1]
  parent_x    = node_depths[parent_node]
  
  # Get child coords
  node_heights = node.height(tr)
  child_node = this_edge[2]
  child_x = node_depths[child_node]
  child_y = node_heights[child_node]
  
  ## Identify signatures to plot
  # Check if child is tip
  #child_type = filter(ape_tbl, child==this_edge[2])$node_type
  if(sum(tr$edge==this_edge[2]) == 1) {
    this_cluster = tr$tip.label[child_node]
  } else {
    this_cluster = tr$node.label[child_node-length(tr$tip.label)]
  }
  # Get sig contribs
  sig_contrib = unlist(filter(sig_contrib_tbl, cluster_id == this_cluster)[1,which(
    grepl('sig', names(sig_contrib_tbl)))])
  
  ## Plotting
  # Define colors
  if(is.na(sig_cols)) {
    sig_cols = c(brewer.pal(n=9, name="Set1"), brewer.pal(n=10, name="Set3")[-2])
  }
  total_len = child_x - parent_x
  
  if(is.na(sig_order)) {
    sig_order = names(sig_contrib)
  }
  # Reorder signature contributions
  sig_contrib = sig_contrib[sig_order]
  sig_cols = sig_cols[sig_order]
  
  for(i in c(1:length(sig_contrib))) {
    this_xleft = parent_x + (total_len * sum(sig_contrib[0:(i-1)]))
    this_xright = parent_x + (total_len * sum(sig_contrib[1:(i)]))
    rect(xleft = this_xleft, xright = this_xright, ybottom = child_y-(0.5*lwidth), ytop = child_y+(0.5*lwidth), col=sig_cols[i])
  }
  
  # Draw legend
  if(draw_legend) {
    # Create alphabetic order
    alpha_ord = order(rename_sigs(names(sig_cols)))
    leg_x = (par('usr')[2])/2
    leg_y = -0.25*par('usr')[4]
    legend(x = leg_x, y=leg_y, title=expression(bold("Signatures")), rename_sigs(names(sig_cols))[alpha_ord], fill=sig_cols[alpha_ord], cex=0.8, ncol=3, xjust=0.5)
  }
}

#' Title
#'
#' @param tr 
#' @param edge_id 
#' @param branch_strength_tbl 
#'
#' @return
draw_branch_strength <- function(tr, edge_id, branch_strength_tbl) {
  
}

highlight_tree_nodes <- function(tr, pch=15, show_internal_nodes=FALSE, ...) {
  # Highlight all internal nodes in a tree plot
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 04-JUN-2018
  #
  # Args
  #   tr: object of type phylo
  #   pch: symbol of point to plot
  #   ...: args will be forwarded to points() function
  #   
  # Returns: nothing
  
  # Identify internal node IDs
  node_idx = length(tr$tip.label) + 1:length(tr$node.label)
  
  # Get node coordinates
  node_depths = node.depth.edgelength(tr)
  node_heights = node.height(tr)
  x = node_depths[node_idx]
  y = node_heights[node_idx]
    
  # Plot the nodes
  points(x, y, pch=pch, ...)
  
  # Plot the labels
  if(show_internal_nodes) {
    #text(x,y,tr$node.label, cex=0.7)
    for(i in c(1:length(x))) {
      a=legend(x[i],y[i],str_replace(tr$node.label[i], 'Cl.', ''), box.col = "lightgrey", bg = "lightgrey", adj = 0, cex=0.5, xjust=0.5, yjust=0.5, plot=FALSE)
      # box size reduced by factor 0.75
      a=a$rect
      mid = a$top - 0.5*a$h
      reduction = 0.5
      
      #rect(xleft=a$left, ytop=mid+0.5*reduction*a$h, xright=a$left+a$w, ybottom=mid-0.5*reduction*a$h, bty='n')
      
      #legend(x[i],y[i],str_replace(tr$node.label[i], 'Cl.', ''), adj = 0, cex=0.5, xjust=0.5, yjust=0.5)
      rect(xleft=a$left, ytop=mid+0.5*reduction*a$h, xright=a$left+a$w, ybottom=mid-0.5*reduction*a$h, bty='n', col = "lightgrey", border = "lightgrey")
      text(x[i]-0.2*a$w,y[i],str_replace(tr$node.label[i], 'Cl.', ''), adj = 0, cex=0.5, xjust=0.5, yjust=0.5)
      #if(x[i]==0) { this_x=1 } else { this_x=x[i] }
      #textbox(x=this_x,y=y[i],str_replace(tr$node.label[i], 'Cl.', ''), justify='c', cex=0.5, border='lightgrey', fill='lightgrey', margin=0.1)
    }
  }
}

get_all_descendants <- function(tr, node_id, curr=NULL){
  # Returns all descendants of node with id node_id in tree tr
  if(is.null(curr)) curr<-vector()
  daughters<-tr$edge[which(tr$edge[,1]==node_id),2]
  curr<-c(curr,daughters)
  w<-which(daughters>=length(tr$tip))
  if(length(w)>0) for(i in 1:length(w)) 
    curr <- get_all_descendants(tr,daughters[w[i]],curr)
  return(curr)
}

get_terminal_descendants <- function(tr, node_id) {
  all_descendants = get_all_descendants(tr, node_id)
  # Filter out non-terminal descendants
  terminal_descendants = all_descendants[all_descendants <= Ntip(tree)]
  return(terminal_descendants)
}

get_root_node <- function(tr) {
  return(Ntip(tr) + 1)
}

get_root_descendants <- function(tr) {
  root_node = get_root_node(tr)
  return(tr$edge[which(tr$edge[,1] == root_node),2])
}

get_node_depth <- function(tr, node_id) {
  # Returns the depth of the node with ID node_id in tree tr
  return(node.depth.edgelength(tr)[node_id])
}

get_max_depth <- function(tr) {
  # Returns the maximum depth of the tree tr
  return(max(node.depth.edgelength(tr)))
}

get_max_nesting <- function(tr) {
  max_nesting = 0
  
  # Start from all root descendants
  root_des = get_root_descendants(tr)
  
  # Walk to furthest tip, record level of nesting
  for(i in c(1:length(root_des))) {
    has_descendant = TRUE
    this_nesting = 0
    this_ancestor = root_des[i]
    while(has_descendant) {
      new_anc = tr$edge[which(tr$edge[,1] %in% this_ancestor),2]
      if(length(new_anc)==0) { 
        has_descendant=FALSE
      } else { 
        this_ancestor = new_anc 
      }
      this_nesting = this_nesting + 1
    }
    if(this_nesting > max_nesting) {
      max_nesting = this_nesting
    }
  }
  
  return(max_nesting)
}

get_node_nesting <- function(tr, node_id) {
  # Returns the number of branchings that have occurred up to node with node_id
  node_nesting = -1
  has_ancestor = TRUE
  
  current_node = node_id
  while(has_ancestor) {
    current_node = tr$edge[which(tr$edge[,2] == current_node),1]
    if(length(current_node)==0) { has_ancestor = FALSE }
    node_nesting = node_nesting + 1
  }
  return(node_nesting)
}

get_root_of_node <- function(tr, node_id) {
  # Returns the node upstream in the tree that connects to the root node
  root_node = get_root_node(tr)
  this_node = node_id
  is_root = FALSE
  while(is_root==FALSE) {
    this_ancestor = tr$edge[which(tr$edge[,2]==this_node),1]
    if(this_ancestor == root_node) {
      is_root = TRUE
    } else {
      this_node = this_ancestor
    }
  }
  return(this_node)
}

load_sig_tables <- function(path_to_exposures, path_to_decode=NA, donor_id, max_contrib_thr = 0.1, sig_prefix='sigP', prefix_replace=NA) {
  # Wrapper function that loads tables that are required for signature plotting
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 04-JUN-2018
  #
  # Args
  #   path_to_exposures: path to a CSV containing the exposures of all clusters to all signatures
  #   donor_id: PD ID
  #   max_contrib_thr: see sig_exposure_per_cluster()
  #   sig_prefix: string to identify the name ofsignatures used as priors
  #   prefix_replace: string to replace the original prefix with, eg. COSMIC
  #   
  # Returns: named list
  # list(
  #  $exposures => table from path_to_exposures
  #  $sig_contrib => returned from sig_exposure_per_cluster
  #  $sig_colors => returned from load_sig_colors
  #  $sig_order => returned from get_sig_order
  # )

  exposures = tbl_df(fread(path_to_exposures))
  # Check if Signature column contains any prior signatures, if not rename Signatures.
  if(!grepl('P', paste(unique(exposures$Signature), collapse=''))) {
    exposures = exposures %>%
      mutate(Signature = ifelse(Signature==0, Signature, sprintf('P%s', Signature)))
  }
  
  if(!is.na(path_to_decode)) {
    hdp_decode = tbl_df(fread(path_to_decode))
    exposures = decode_hdp_sig_id(exposures, hdp_decode)
  }
  
  sig_contrib = sig_exposure_per_cluster(exposures, donor_id = donor_id, max_contrib_thr = max_contrib_thr)
  sig_colors = load_sig_colors(sig_contrib, sig_prefix = sig_prefix)
  sig_order = get_sig_order(sig_contrib)
  
  if(!is.na(prefix_replace)) {
    replace_pattern = str_replace(sig_prefix, 'sig', '')
    exposures = mutate(exposures, Signature = str_replace(Signature, replace_pattern, prefix_replace))
    names(sig_contrib) <- str_replace(names(sig_contrib), sig_prefix, prefix_replace)
    names(sig_colors) <- str_replace(names(sig_colors), sig_prefix, prefix_replace)
    sig_order <- str_replace(sig_order, sig_prefix, prefix_replace)
  }
  
  return(list('exposures'=exposures, 
              'sig_contrib'=sig_contrib, 
              'sig_colors'=sig_colors, 
              'sig_order'=sig_order))
}

sig_exposure_per_cluster <- function(exposures, donor_id, max_contrib_thr = 0.1) {
  # Loads a table with the exposure of each cluster to each signature
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 04-JUN-2018
  #
  # Args
  #   exposures: tibble with columns sample (PDID_clusterid), donor, Signature, Exposure
  #   donor_id: PD ID
  #   max_contrib_thr: signatures that have a maximal exposure lower than this value are removed
  #   
  # Returns: 
  # sig_contrib, a tibble with following columns:
  # - donor (PD ID)
  # - N signatures, name sigOther, sig0, sigNxx, sigPxx
  # - cluster_id, eg. Cl.1
  sig_contrib = exposures %>%
    filter(donor == donor_id) %>%
    mutate(Signature = paste('sig', Signature, sep='')) %>%
    group_by(Signature) %>%
    mutate(max_contrib = max(Exposure)) %>%
    ungroup() %>%
    filter(max_contrib > max_contrib_thr) %>%
    group_by(sample) %>%
    #filter(Signature != 'sig0') %>%
    mutate(sigOther = 1-sum(Exposure)) %>%
    ungroup() %>%
    dplyr::select(-max_contrib) %>%
    spread(key='Signature', value='Exposure') %>%
    mutate(sample = str_replace(sample, paste0(donor_id, '_'), '')) %>%
    mutate(cluster_id = paste('Cl', sample, sep='.')) %>%
    dplyr::select(-sample)
  sig_contrib
}

load_sig_colors <- function(sig_contrib, sig_prefix='sigP') {
  # Loads colors to draw signatures
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 04-JUN-2018
  #
  # Args
  #   sig_contrib: tibble, output from sig_exposure_per_cluster()
  #   
  # Returns: 
  # hdp_cols: named list of colors for each signature (signature ID being names of list)
  hdp_allcomps = names(sig_contrib)[grepl('sig', names(sig_contrib))]
  cosmic_cols = get_deconstructSigs_colors()
  cosmic_col_tbl = tibble(cosmic_sig = names(cosmic_cols), col_val = cosmic_cols) %>%
    mutate(cosmic_sig = paste(sig_prefix, cosmic_sig, sep=''))
  hdp_col_tbl = tibble(hdp_sig = hdp_allcomps) %>%
    left_join(cosmic_col_tbl, by=c('hdp_sig'='cosmic_sig')) %>%
    mutate(col_val = ifelse(hdp_sig == 'sigOther', 'grey', col_val))
  # Add new colors for new sigs
  num_new = sum(is.na(hdp_col_tbl$col_val))
  new_cols = setdiff(cosmic_col_tbl$cosmic_sig, hdp_col_tbl$hdp_sig)[2:(num_new+1)]
  hdp_col_tbl$col_val[which(is.na(hdp_col_tbl$col_val))] = cosmic_col_tbl$col_val[which(cosmic_col_tbl$cosmic_sig %in% new_cols)]
  hdp_cols = hdp_col_tbl$col_val
  names(hdp_cols) = hdp_col_tbl$hdp_sig
  hdp_cols
}

get_sig_order <- function(sig_contrib) {
  # Loads colors to draw signatures
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 04-JUN-2018
  #
  # Args
  #   sig_contrib: tibble, output from sig_exposure_per_cluster()
  #   
  # Returns: 
  # sig_order: list of signature IDs, in the order to plot
  
  sig_order_tbl = sig_contrib %>%
    gather(key='Signature', value='Exposure', names(sig_contrib)[which(grepl('sig', names(sig_contrib)))]) %>%
    group_by(Signature) %>% summarise(Exposure = median(Exposure)) %>%
    arrange(desc(Exposure))
  sig_order = sig_order_tbl$Signature
  return(sig_order)
}

rename_sigs <- function(sig_name) {
  if(grepl('PCAWG', sig_name)) {
    sig_name = str_replace(sig_name, 'sigPCAWG', 'PCAWG ')
  } else {
    sig_name = str_replace(sig_name, 'sigP', 'COSMIC ')
  }
  sig_name = str_replace(sig_name, 'sigN', 'New ')
  sig_name = str_replace(sig_name, 'sigOther', 'Other')
  sig_name = str_replace(sig_name, 'sig0', 'Noise')
  sig_name
}

assign_clusters_samples <- function(clusters, contrib_thr=0.1) {
  # Assigns samples to each cluster based on a cutoff
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 18-JUN-2018
  #
  # Args
  #   clusters: tibble where rows are clusters and columns are samples. Values are the contribution of a sample to a cluster
  #   contrib_thr: threshold value to assign a sample to a cluster.
  #   
  # Returns: 
  # clusters_bin: tibble with binary assignments of clusters and samples
  
  # Make sure all cluster names have same format (Cl.x)
  ch_name = which(grepl('PD', names(clusters)))
  if(length(ch_name)>0) { names(clusters)[ch_name] <- strsplit(names(clusters)[ch_name], '_')[[1]][2] }
  
  # Threshold each sample-cluster pair
  clusters_bin = clusters %>%
    gather(key='sample', value='contrib', names(clusters)[grepl('lo', names(clusters))]) %>%
    mutate(contrib = ifelse(contrib >= contrib_thr, 1, 0)) %>%
    group_by(cluster_id) %>% mutate(contrib_sum = sum(contrib)) %>% ungroup() %>%
    filter(contrib==1)
  
  return(clusters_bin)
}

get_poly_perimeter <- function(poly_mat) {
  num_coords = dim(poly_mat)[1]
  if(sum(poly_mat[1,]!=poly_mat[num_coords,])) {
    poly_mat = poly_mat[c(1:num_coords, 1),]
  }
  peri = 0
  for(i in c(2:num_coords)) {
    peri = peri + as.numeric(dist(poly_mat[c(i-1,i),]))
  }
  
  return(peri)
}

assign_phys_size <- function(clusters, phys_coord, contrib_thr=0.1, size_type='perimeter') {
  # Assigns physical coordinates to each cluster contained in tree.
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 16-JUN-2018
  #
  # Args
  #   clusters: tibble where rows are clusters and columns are samples. Values are the contribution of a sample to a cluster
  #   phys_coord: tibble giving for each sample (row) the x,y,z coordinates
  #   contrib_thr: threshold value to assign a sample to a cluster.
  #   
  # Returns: 
  # clusters_coords: tibble where each row is a cluster and contains x,y,z coordinates
  
  # Locate each cluster in physical space
  clusters_coords = locate_clusters(clusters, phys_coord, contrib_thr = contrib_thr)
  
  # Determine the size of each cluster
  uniq_clust = unique(clusters_coords$cluster_id)
  
  clust_dist = clusters_coords %>% dplyr::select(cluster_id, no.of.mutations.assigned) %>% distinct()
  clust_dist$med_dist = sapply(1:nrow(clust_dist), function(x) {
    #message(x)
    this_mat = as.matrix(clusters_coords[which(clusters_coords$cluster_id==clust_dist$cluster_id[x]), c('x', 'y', 'z')])
    polygon_points = chull(this_mat[,'x'], this_mat[,'y']) # the function identifies the points that form a convex hull polygon around all points included in the cluster
    
    if(size_type=='area') {
      if(length(polygon_points) == 1) {
        return(0)
      } else if(length(polygon_points) == 2) {
        return(as.numeric(dist(this_mat[polygon_points,])))
      } else {
        chullpoly = Polygon(this_mat[c(polygon_points, polygon_points[1]),c('x','y')])
        return(chullpoly@area)
      }
    }
    if(size_type=='perimeter') {
      if(length(polygon_points) == 1) {
        return(0)
      } else if(length(polygon_points) == 2) {
        return(2*as.numeric(dist(this_mat[polygon_points,])))
      } else {
        #sf::st_polygon(this_mat[c(polygon_points, polygon_points[1]),c('x','y')])
        #chullpoly = Polygon(this_mat[c(polygon_points, polygon_points[1]),c('x','y')])
        #geosphere::perimeter(chullpoly)
        #geosphere::perimeter(this_mat[c(polygon_points, polygon_points[1]),c('x','y')])
        return(get_poly_perimeter(this_mat[c(polygon_points, polygon_points[1]),c('x','y')]))
        #return(-1*pracma::polyarea(this_mat[c(polygon_points, polygon_points[1]),'x'], this_mat[c(polygon_points, polygon_points[1]),'y']))
      }
    }
    if(size_type=='max_radius') {
      if(length(polygon_points) == 1) {
        return(0)
      } else if(length(polygon_points) == 2) {
        return(as.numeric(dist(this_mat[polygon_points,]))/2)
      } else {
        poly_mat = this_mat[c(polygon_points),c('x','y')]
        # Get the centroid
        poly_centroid = colSums(poly_mat)/length(polygon_points)
        
        # Find the largest distance
        poly_dist = abs(poly_mat-matrix(rep(poly_centroid, each=length(polygon_points)), ncol=2))
        max_rad = max(sqrt(rowSums(poly_dist^2)))
        
        return(max_rad)
        #return(-1*pracma::polyarea(this_mat[c(polygon_points, polygon_points[1]),'x'], this_mat[c(polygon_points, polygon_points[1]),'y']))
      }
    }
    
    #all_dist = abs(this_mat[polygon_points,'x'] - this_mat[polygon_points,'y'])
    #median(sort(all_dist, decreasing = TRUE)[1:3], na.rm=TRUE)
    #polygon_points = polygon_points[1:3] # Take only top 3 points
    #pracma::polyarea(this_mat[polygon_points,'x'], this_mat[polygon_points,'y'])
    #median(as.numeric(dist(this_mat[polygon_points,])), na.rm=TRUE)
  })
  clust_dist = clust_dist %>%
    mutate(med_dist = ifelse(is.na(med_dist), 0, med_dist))
  
  return(clust_dist)
}

assign_divergence_time <- function(tr) {
  
  # For each leave node, divergence time is 1
  div_time_leaves = tibble(cluster_id = tr$tip.label, t=1)
  
  # Now handle internal nodes
  internal_nodes = unique(tr$edge[,1])
  internal_node_labels = tr$node.label
  
  # For each internal node:
  div_time_internal = list()
  for(i in c(1:length(internal_nodes))) {
    this_node = internal_nodes[i]
    this_label = internal_node_labels[i]
    if(grepl('PD', this_label)) {
      # Skip if this is the root node
      next
    }
    # - identify the length of the branch at the time where it splits (A)
    internal_depth = get_node_depth(tr, this_node)
    # - identify the median lengths of all branches diverging from it (B)
    descendants = get_terminal_descendants(tr = tr, node = this_node)
    tip_depth = median(get_node_depth(tr, descendants))
    # - calculate A/B - the relative divergence time.
    time_since_div = internal_depth / tip_depth
    # Append node infos to list
    div_time_internal[[i]] = tibble(cluster_id = this_label, t=time_since_div)
  }
  divergence_time = rbind(div_time_leaves, do.call('rbind', div_time_internal))
  
  return(divergence_time)
}

assign_divergence_time2 <- function(tr) {
  
  ## Define table of nodes and labels
  node_tbl = get_node_table(tr)
  node_tbl$ancestor = NA
  node_tbl$ancestor_depth = NA
  
  # Remove root node
  node_tbl = node_tbl # %>% filter(node_type != 'root')
  
  for(i in c(1:dim(node_tbl)[1])) {
    ## Identify position of each branch point
    # For each node, check what parental node is.
    this_node = node_tbl$node_id[i]
    this_ancestor = tr$edge[which(tr$edge[,2] == this_node), 1]
    
    if(node_tbl$node_type[i]!='root') {
      # Store position (number of muts) of parental node.
      node_tbl$ancestor[i] = this_ancestor
      node_tbl$ancestor_depth[i] = get_node_depth(tr, this_ancestor)
    } else {
      node_tbl$ancestor[i] = NA
      node_tbl$ancestor_depth[i] = 0
    }
  }
  
  ## Calculate median total branch length
  # For each leave node, determine total length of branch (num muts)
  tips = (filter(node_tbl, node_type=='tip'))$node_id
  tip_depth = sapply(tips, function(x) { get_node_depth(tr, x) })
  
  # Calculate median total branch length
  median_tot_len = median(tip_depth)
  
  ## Assign relative time: breakpoint time / median tot branch length
  node_tbl$median_tot_len = median_tot_len
  node_tbl = node_tbl %>%
    mutate(div_time = ancestor_depth / median_tot_len)
  
  return(node_tbl)
}

get_node_table <- function(tr) {
  # Returns a table containing for each node its ID in the edges table, as well as its label and type (tip, internal, root)
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 23-JUL-2018
  #
  # Args
  #   tr: ape object
  #   
  # Returns: 
  # node_tbl
  
  # Determine tips
  tip_id = 1:length(tr$tip.label)
  tip_tbl = tibble(node_id = tip_id, node_label = tr$tip.label, node_type = 'tip')
  
  # Determine internal nodes
  internal_tbl = tibble(node_id = sort(unique(tr$edge[, 1])), node_label = tr$node.label, node_type = 'internal')
  internal_tbl$node_type[1] = 'root'
  
  # Merge the two tables
  node_tbl = rbind(tip_tbl, internal_tbl)
  
  return(node_tbl)
}

get_divergence_spacetime <- function(tr, clusters_coords) {
  # Calculates the time of divergence of each node, and the distance covered since divergence
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 16-JUN-2018
  #
  # Args
  #   tr: ape object
  #   clusters_coords: table listing each cluster and its coordinates (generated with assign_phys_coord)
  #   
  # Returns: 
  # divergence_tbl: tibble with two columns: time_since_div (relative time since divergence) and dist_since_div (distance covered since divergence)
  #   contains minimum and maximum values.
  
  # Identify internal nodes
  internal_nodes = unique(tr$edge[,1])
  #internal_nodes = internal_nodes[2:length(internal_nodes)]
  
  # Loop through each node
  divergence_time = list()
  for(i in c(1:length(internal_nodes))) {
    this_node = internal_nodes[i]
    descendants = get_terminal_descendants(tr = tr, node = this_node)
    
    # Get node labels
    node_label = tr$node.label[i]
    descendant_labels = tr$tip.label[descendants]
    descendant_labels
    
    # -- Get the median distance between descendants
    descendant_coords = clusters_coords %>%
      filter(cluster_id %in% descendant_labels)
    # Skip this branch if there is only one descendant
    if(dim(descendant_coords)[1]==1) { next }
    median_dist = median(dist(as.matrix(descendant_coords[,c('x','y','z')])))
    
    # -- Calculate the time since divergence
    # First: internal node
    internal_depth = get_node_depth(tr, this_node)
    # Second: median branch length of all descendants
    tip_depth = median(get_node_depth(tr, descendants) - internal_depth)
    # Now calculate time since divergence
    time_since_div = tip_depth / (internal_depth + tip_depth)
    
    # 
    divergence_time[[i]] = tibble(time_since_div = time_since_div, dist_since_div = median_dist)
  }
  
  # Add a point for 0 divergence
  divergence_time[[length(divergence_time) + 1]] = tibble(time_since_div = 0, dist_since_div = 0)
  
  # Concatenate all node outputs
  divergence_tbl = tbl_df(do.call('rbind', divergence_time))
  
  return(divergence_tbl)
}

locate_clusters <- function(clusters, phys_coord, contrib_thr = 0.1) {
  #' Assign physical coordinates of samples to each cluster
  #' 
  #' Identifies contributing samples of each cluster and then assigns
  #' their physical coordinates to that cluster. Adds another 'synthetic cluster'
  #' for the root, which contains all samples.
  #' 
  #' @param clusters mxn matrix of m clusters and n samples. Values are average VAFs of each sample-cluster pair.
  #' @param phys_coord tibble of physical coordinates of each sample
  #' @param contrib_thr threshold above which a sample is considered as contributing to a cluster
  #' 
  #' @return tibble of clusters and physical coordinates
  #' 
  #' @author Simon F Brunner, \email{sb50@@sanger.ac.uk}
  
  # Make sure all cluster names have same format (Cl.x)
  ch_name = which(grepl('PD', names(clusters)))
  names(clusters)[ch_name] <- strsplit(names(clusters)[ch_name], '_')[[1]][2]
  
  # Threshold each sample-cluster pair
  clusters_bin = assign_clusters_samples(clusters = clusters, contrib_thr = contrib_thr)
  
  # Prepare physical coordinates table
  phys_coord_cl = phys_coord
  phys_coord_cl$sample = sapply(phys_coord_cl$sample, function(x) {
    strsplit(x, '_')[[1]][2]
  })
  phys_coord_cl = phys_coord_cl %>%
    dplyr::select(sample, x, y, z)
  
  # Merge coordinates with clusters
  clusters_coords = inner_join(clusters_bin, phys_coord_cl, by='sample') %>% filter(!is.na(x), !is.na(y), !is.na(z))
  
  # Add root
  root_coords = clusters_coords %>%
    dplyr::select(sample, x, y, z) %>% distinct() %>%
    mutate(cluster_id = 'root', no.of.mutations.assigned = 0, contrib=1, contrib_sum=NA)
  clusters_coords = rbind(clusters_coords, root_coords)
  
  return(clusters_coords)
  
}

get_polygon_from_cluster <- function(clusters_coords) {
  
  # Identify unique clusters
  uniq_clust = unique(clusters_coords$cluster_id)
  #clust_dist = clusters_coords %>% dplyr::select(cluster_id, no.of.mutations.assigned) %>% distinct()
  
  # Loop through each cluster
  polygon_lst = list()
  for(i in c(1:length(uniq_clust))) {
    this_mat = as.matrix(clusters_coords[which(clusters_coords$cluster_id==uniq_clust[i]), c('x', 'y', 'z')])
    
    # Determine convex hull of polygon made by cluster coords
    polygon_points = chull(this_mat[,'x'], this_mat[,'y']) # the function identifies the points that form a convex hull polygon around all points included in the cluster
    
    # Construct coordinates table of polygon
    this_polygon = tibble(x = this_mat[polygon_points,'x'], y = this_mat[polygon_points,'y'], z = this_mat[polygon_points,'z'])
    this_polygon$cluster_id = uniq_clust[i]
    
    polygon_lst[[i]] = this_polygon
  }
  polygon_tbl = do.call('rbind', polygon_lst)
  
  return(polygon_tbl)
  # clust_dist = clusters_coords %>% dplyr::select(cluster_id, no.of.mutations.assigned) %>% distinct()
  # clust_dist$med_dist = sapply(1:nrow(clust_dist), function(x) {
  #   #message(x)
  #   this_mat = as.matrix(clusters_coords[which(clusters_coords$cluster_id==clust_dist$cluster_id[x]), c('x', 'y', 'z')])
  #   polygon_points = chull(this_mat[,'x'], this_mat[,'y']) # the function identifies the points that form a convex hull polygon around all points included in the cluster
  #   
  #   if(length(polygon_points) == 1) {
  #     return(0)
  #   } else if(length(polygon_points) == 2) {
  #     return(as.numeric(dist(this_mat[polygon_points,])))
  #   } else {
  #     chullpoly = Polygon(this_mat[c(polygon_points, polygon_points[1]),c('x','y')])
  #     return(chullpoly@area)
  #   }
  # })
  # clust_dist = clust_dist %>%
  #   mutate(med_dist = ifelse(is.na(med_dist), 0, med_dist))
  # 
  # return(clust_dist)
}

draw_polygons <- function(polygon_tbl, outpath='', max_col=NA) {
  
  # Assign the number of coordinates of each cluster
  polygon_tbl = polygon_tbl %>% group_by(cluster_id) %>% mutate(num_coords = n()) %>% ungroup()
  
  # Determine the color range
  if(is.na(max_col)) {
    max_col = polygon_tbl$num_coords
  }
  col_range = rainbow(max_col)
  
  # Loop through each cluster, from biggest to smallest
  uniq_clust = (polygon_tbl %>% distinct(cluster_id, num_coords) %>% arrange(desc(num_coords)))$cluster_id
  label_tbl = polygon_tbl %>% group_by(x,y,z) %>% filter(n()==1) %>% ungroup()
  
  # Build table of coordinates where to place labels
  label_candidates = polygon_tbl %>% group_by(x,y,z) %>% mutate(num_clust = n()) %>% ungroup()
  label_tbl = label_candidates %>% filter(num_clust==1) %>% group_by(cluster_id) %>% filter(row_number()==1) %>% ungroup()
  missed_clust = uniq_clust[!(uniq_clust %in% label_tbl$cluster_id)]
  add_tbl = polygon_tbl %>% filter(cluster_id %in% missed_clust) %>% group_by(cluster_id) %>% filter(row_number()==1) %>% ungroup()
  label_tbl = rbind(label_tbl %>% dplyr::select(-num_clust), add_tbl)
  
  label_tbl = polygon_tbl %>% group_by(cluster_id) %>% summarise(x=mean(x), y=mean(y), num_coords=n()) %>% ungroup()
  
  col.groups = polygon_tbl %>% dplyr::select(cluster_id, col) %>% distinct()
  cols = col.groups$col
  names(cols)<-col.groups$cluster_id
  
  g = polygon_tbl %>%
    ggplot(aes(x=x, y=-y, group=cluster_id)) + 
    geom_polygon(size=3, aes(fill=cluster_id, color=cluster_id)) +
    geom_point(size=3, aes(fill=cluster_id, color=cluster_id)) +
    geom_point(color='black', shape='+', size=3) +
    theme_void() +
    geom_label_repel(data = label_tbl, aes(x=x, y=-y, label=cluster_id, fill=cluster_id), size=6, segment.color	='black') +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    coord_fixed(ratio=1) +
    theme(legend.position="none")
  # g = polygon_tbl %>%
  #   ggplot(aes(x=x, y=-y, group=cluster_id)) + 
  #   geom_polygon(size=3, aes(fill=as.factor(num_coords), color=as.factor(num_coords))) +
  #   geom_point(size=3, aes(fill=as.factor(num_coords), color=as.factor(num_coords))) +
  #   geom_point(color='black', shape='+', size=3) +
  #   theme_void() +
  #   geom_label_repel(data = label_tbl, aes(x=x, y=-y, label=cluster_id, fill=as.factor(num_coords)), size=3, segment.color	='black') +
  #   scale_fill_manual(values = cols) +
  #   coord_fixed(ratio=1) +
  #   theme(legend.position="none")
  
  return(g)
  # plot(x=NA, y=NA, type="n", xlab="dist", ylab="age", 
  #      xlim=c(min(polygon_tbl$x), max(polygon_tbl$x)), 
  #      ylim=c(min(polygon_tbl$y), max(polygon_tbl$y)))
  # 
  # for(i in c(1:length(uniq_clust))) {
  #   this_cluster = uniq_clust[i]
  #   this_coords = filter(polygon_tbl, cluster_id == this_cluster)
  #   
  #   num_coords = dim(this_coords)[1]
  #   
  #   if(num_coords>2) {
  #     polygon(this_coords$x, this_coords$y, col=col_range[num_coords], border=col_range[num_coords], lwd=5)
  #   } else if(num_coords==2) {
  #     lines(this_coords$x, this_coords$y, col=col_range[num_coords], lwd=5)
  #   } else if(num_coords==1) {
  #     points(this_coords$x, this_coords$y, pch=21, cex=1, col=col_range[num_coords], bg=col_range[num_coords])
  #   }
  #   
  #   # Draw cluster coordinate
  #   points(this_coords$x, this_coords$y, pch=3, cex=0.3)
  #   text(this_coords$x[1], this_coords$y[1], labels=this_cluster)
  #   #points(this_coords$x[1], this_coords$y[1], pch=3)
  # }
}

#' Extract branches of ape tree object
#' 
#' Starts at tip_node and walks back entire branch up to root.
#' 
#' @param tr ape tree object
#' @param tip_node ID of tip node to start walking backwards from
#' 
#' @return vector containing node IDs of nodes in branch (going from tip to root)
#' 
#' @author Simon F Brunner, \email{sb50@@sanger.ac.uk}
get_full_branch <- function(tr, tip_node) {
  
  # Walk branch backwards and store nodes
  edges = tr$edge
  
  has_ancestor = TRUE
  
  current_node = tip_node
  branch = c(current_node)
  counter = 2
  while(has_ancestor) {
    edge_id = which(edges[,2]==current_node)
    
    # Decide if there is an ancestor
    if(length(edge_id)>0) {
      # Ancestor present, thus make it the new current_node to continuing walking backwards
      current_node = edges[edge_id,1]
      branch[counter] = current_node
      counter = counter+1
    } else {
      # Ancestor absent, leave loop
      has_ancestor = FALSE
    }
  }
  
  return(branch)
}