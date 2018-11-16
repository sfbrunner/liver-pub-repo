# Pileup helper functions
# --
# This set of functions supports the handling of pileup tables. Examples: 
# - Extract alt counts at specified sites
# - Define the entropy at specified coordinates
# 
# --
# /// Author --- SIMON FELIX BRUNNER
# /// Creation date --- 10-MAY-2018

extract_mut_counts <- function(target_file_tbl, target_vars, input_dir, output_dir) {
  # Takes a table containing names and locations of pileup tables and extracts mutation calls at specified sites.
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 10-MAY-2018
  #
  # Args
  #   target_file_tbl: tibble with two columns, sample and fname
  #   target_vars: tibble with at least two columns: "coord_id" defining target sites in chr_pos format (eg. X_321535)
  #                   and "mut_id" defining target sites in chr_pos_ref_alt format (eg. X_321535_A_G)
  #   input_dir: directory containing pileup files
  #   output_dir: directory to store output files
  # Returns:
  #   list(
  #     target_alt_tbl, table containing the alt counts
  #     target_depth_tbl, table containing the depth at each site
  # )
  
  target_lst = list()
  for(i in c(1:dim(target_file_tbl)[1])) {
    # Append all samples into large table
    this_sample = target_file_tbl$sample[i]
    message(sprintf('%s Importing sample %s', Sys.time(), this_sample))
    this_tbl = tbl_df(fread(file.path(input_dir, target_file_tbl$fname[i])))
    this_tbl = this_tbl %>%
      mutate(coord_id = paste(chrom, pos, sep='_')) %>%
      filter(coord_id %in% target_vars$coord_id)
    this_tbl$sample = this_sample
    target_lst[[i]] = this_tbl
  }
  target_tbl = do.call('rbind', target_lst)
  target_tbl = remove_strand_info(target_tbl)
  
  # Calculate depth
  target_tbl = target_tbl %>% 
    mutate(depth = A+C+G+`T`+N+INS+DEL) %>%
    dplyr::select(-INS, -DEL, -N) %>%
    #distinct() %>%
    gather(key='alt', value='counts', c('A', 'C', 'G', 'T'))
  
  # Generate IDs and filter for target variants
  target_tbl = target_tbl %>%
    mutate(mut_id = paste(chrom, pos, ref, alt, sep='_')) %>%
    filter(mut_id %in% target_vars$mut_id)
  
  # Generate depth and alt tables
  target_alt_tbl = target_tbl %>%
    dplyr::select(-depth, -mapq) %>%
    spread(key='sample', value='counts', fill = 0)
  target_depth_tbl = target_tbl %>%
    dplyr::select(-counts, -mapq) %>%
    spread(key='sample', value='depth', fill = 0)
  
  # Export alt allele count and depth tables
  fpath_alt_tbl = file.path(output_dir, 'target_alt.csv')
  fwrite(target_alt_tbl, fpath_alt_tbl)
  message(sprintf('Wrote alt counts table to %s', fpath_alt_tbl))
  fpath_depth_tbl = file.path(output_dir, 'target_depth.csv')
  fwrite(target_depth_tbl, fpath_depth_tbl)
  message(sprintf('Wrote depth table to %s', fpath_depth_tbl))
  
  return(list(target_alt_tbl=target_alt_tbl, target_depth_tbl=target_depth_tbl))
}

generate_normal_panel <- function(target_file_tbl, target_vars, input_dir, output_dir) {
  # Takes a table containing names and locations of pileup tables and generates normal panel.
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 10-MAY-2018
  #
  # Args
  #   target_file_tbl: tibble with two columns, sample and fname
  #   target_vars: tibble with at least a column "coord_id" defining target sites in chr_pos format (eg. X_321535)
  #   input_dir: directory containing pileup files
  #   output_dir: directory to store output files
  # Returns:
  #   list(
  #     
  # )
  
  curr_tbl = tibble(coord_id = target_vars$coord_id, 
                    chrom = unlist(lapply(target_vars$mut_id, function(x) {
                      strsplit(x, '_')[[1]][1]
                    })), 
                    pos = unlist(lapply(target_vars$mut_id, function(x) {
                      strsplit(x, '_')[[1]][2]
                    })), 
                    ref = unlist(lapply(target_vars$mut_id, function(x) {
                      strsplit(x, '_')[[1]][3]
                    })), 
                    A=0, C=0, G=0, `T`=0, N=0, INS=0, DEL=0
                    ) %>% distinct()
  
  for(i in c(1:dim(target_file_tbl)[1])) {
    # Load table
    this_sample = target_file_tbl$sample[i]
    message(sprintf('%s Parsing normal table for sample %s', Sys.time(), this_sample))
    this_tbl = tbl_df(fread(file.path(input_dir, target_file_tbl$fname[i])))
    
    this_tbl = this_tbl %>%
      mutate(coord_id = paste(chrom, pos, sep="_")) %>%
      filter(coord_id %in% target_vars$coord_id) %>%
      dplyr::select(-mapq)
    this_tbl = remove_strand_info(this_tbl)
    
    this_matrix = as.matrix(this_tbl[,c('A', 'C', 'G', 'T', 'N', 'INS', 'DEL')])
    rownames(this_matrix)<-this_tbl$coord_id
    this_matrix[is.na(this_matrix)] <- 0
    
    match_idx = match(rownames(this_matrix), curr_tbl$coord_id, nomatch=0)
    
    curr_tbl[match_idx, c('A', 'C', 'G', 'T', 'N', 'INS', 'DEL')] = curr_tbl[match_idx, c('A', 'C', 'G', 'T', 'N', 'INS', 'DEL')] + this_matrix
  }
  
  # Postprocess
  curr_tbl = clean_pile_tbl(curr_tbl)
  
  # Write to disk
  fwrite(curr_tbl, file.path(output_dir, 'normal_panel.csv'))
  
  return(curr_tbl)
}



calculate_purity <- function(pile_tbl) {
  # Calculates two new columns: "gini" and "entropy" 
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 11-MAY-2018
  #
  # Args
  #   pile_tbl: tibble with at least following columns: coord_id (format chrom_pos), A, C, G, T, N, INS, DEL
  # Returns:
  #   list(
  #     purity_tbl
  # )
  purity_tbl = pile_tbl %>%
    gather(key='base', value='counts', c('A', 'C', 'G', 'T', 'N', 'INS', 'DEL')) %>%
    group_by(coord_id) %>%
    mutate(entropy = -sum(counts/sum(counts)*log2(counts/sum(counts)), na.rm=T)) %>%
    mutate(gini = 1-sum((counts/sum(counts))^2)) %>%
    ungroup() %>%
    spread(key='base', value='counts')
  return(purity_tbl)
}

remove_strand_info <- function(pile_tbl) {
  #pile_tbl$A = rowSums(pile_tbl[,c('A', 'a')], na.rm=T)
  pile_tbl$A = pile_tbl$A + pile_tbl$a
  pile_tbl$C = pile_tbl$C + pile_tbl$c
  pile_tbl$G = pile_tbl$G + pile_tbl$g
  pile_tbl$T = pile_tbl$T + pile_tbl$t
  pile_tbl$N = pile_tbl$N + pile_tbl$n
  pile_tbl$INS = pile_tbl$INS + pile_tbl$ins
  pile_tbl$DEL = pile_tbl$DEL + pile_tbl$del
  pile_tbl = pile_tbl %>% dplyr::select(-a, -c, -g, -t, -n, -ins, -del)
  pile_tbl
}

clean_pile_tbl <- function(pile_tbl) {
  pile_tbl = pile_tbl %>%
    group_by(coord_id) %>%
    mutate(A = sum(A, na.rm=T), C = sum(C, na.rm=T), G = sum(G, na.rm=T), `T` = sum(`T`, na.rm=T), 
           N = sum(N, na.rm=T), INS = sum(INS, na.rm=T), DEL = sum(DEL, na.rm=T)) %>%
    ungroup() %>%
    distinct() %>%
    mutate(A = ifelse(is.na(A), 0, A), 
           C = ifelse(is.na(C), 0, C), 
           G = ifelse(is.na(G), 0, G), 
           `T` = ifelse(is.na(`T`), 0, `T`), 
           N = ifelse(is.na(N), 0, N), 
           INS = ifelse(is.na(INS), 0, INS), 
           DEL = ifelse(is.na(DEL), 0, DEL))
  pile_tbl
}

extract_indel_counts <- function(target_file_tbl, target_vars, input_dir, output_dir) {
  # Takes a table containing names and locations of pileup tables and extracts mutation calls at specified sites.
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 10-MAY-2018
  #
  # Args
  #   target_file_tbl: tibble with two columns, sample and fname
  #   target_vars: tibble with at least two columns: "coord_id" defining target sites in chr_pos format (eg. X_321535)
  #                   and "mut_id" defining target sites in chr_pos_ref_alt format (eg. X_321535_A_G)
  #   input_dir: directory containing pileup files
  #   output_dir: directory to store output files
  # Returns:
  #   list(
  #     target_alt_tbl, table containing the alt counts
  #     target_depth_tbl, table containing the depth at each site
  # )
  
  target_lst = list()
  for(i in c(1:dim(target_file_tbl)[1])) {
    # Append all samples into large table
    this_sample = target_file_tbl$sample[i]
    message(sprintf('%s Importing sample %s', Sys.time(), this_sample))
    this_tbl = tbl_df(fread(file.path(input_dir, target_file_tbl$fname[i])))
    this_tbl = this_tbl %>%
      mutate(coord_id = paste(chrom, pos, sep='_')) %>%
      filter(coord_id %in% target_vars$coord_id)
    this_tbl$sample = this_sample
    target_lst[[i]] = this_tbl
  }
  target_tbl = do.call('rbind', target_lst)
  
  # Generate IDs and filter for target variants
  target_tbl = target_tbl %>%
    mutate(mut_id = paste(chrom, pos, ref, alt, sep='_')) %>%
    filter(mut_id %in% target_vars$mut_id)
  
  # Generate depth and alt tables
  target_alt_tbl = target_tbl %>%
    dplyr::select(-depth) %>%
    distinct() %>%
    spread(key='sample', value='count', fill = 0)
  target_depth_tbl = target_tbl %>%
    dplyr::select(-count) %>%
    distinct() %>%
    spread(key='sample', value='depth', fill = 0)
  
  # Export alt allele count and depth tables
  fpath_alt_tbl = file.path(output_dir, 'target_alt.csv')
  fwrite(target_alt_tbl, fpath_alt_tbl)
  message(sprintf('Wrote alt counts table to %s', fpath_alt_tbl))
  fpath_depth_tbl = file.path(output_dir, 'target_depth.csv')
  fwrite(target_depth_tbl, fpath_depth_tbl)
  message(sprintf('Wrote depth table to %s', fpath_depth_tbl))
  
  return(list(target_alt_tbl=target_alt_tbl, target_depth_tbl=target_depth_tbl))
}