# LCM_analysis_helper_functions
# --
# This set of functions is useful for multiple analysis steps of Laser Capture Microdissection (LCM) derived data.
# --
# /// Author --- SIMON FELIX BRUNNER
# /// Creation date --- 2-NOV-2017

library(useful)

get_spread_tbl <- function(caveman_var_calls_tbl, donor_name, minvaf_filter=0.1) {
  # Takes long-form variant calls and converts them into wide-form, indicating which sample contains which mutation.
  # Input table has a single column "sample", whereas output table has sample names as columns with each row being 
  # a distinct mutation. 
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 2-NOV-2017
  #
  # Args
  #   caveman_var_calls_tbl: table with variant calls. Required columns: donor, sample, chrom, ref, alt
  #   donor_name: name of donor to include in analysis, eg. 'PD37105b'
  #   minvaf_filter (default: 0.1): minimum VAF of variant across all samples
  # Returns:
  #   list(
  #     spread_tbl: table where each row is a different group of samples sharing mutations. 
  #       No individual mutations are contained.
  #     spread_tbl_full: table where each row is a different distinct mutation.
  #     spread_tbl_full_vaf: table like spread_tbl_full, but entries are VAFs of each mutation.
  #     spread_tbl_full_altreads: table like spread_tbl_full, but entries are VAFs*coverage of each mutation.
  #     spread_tbl_full_depth: table like spread_tbl_full, but entries are coverage at each mutation site.
  # )
  
  message(donor_name)
  spread_tbl_full = filter(caveman_var_calls_tbl, donor==donor_name) %>%
    mutate(sample=substr(sample,10,nchar(sample))) %>%
    group_by(chrom, pos, ref, alt, donor) %>%
    arrange(sample) %>%
    mutate(median_vaf=median(vaf), minvaf = min(vaf), shared_num_samples = n(), shared_id = paste(sample, collapse='_')) %>% ungroup() %>%
    dplyr::select(chrom, pos, ref, alt, donor, shared_id, shared_num_samples, sample, vaf, minvaf, median_vaf, coverage)
  spread_tbl_full = spread_tbl_full %>%
    mutate(shared_id = factor(shared_id, levels=unique(shared_id[order(nchar(shared_id), decreasing=T)]), ordered=T))
  spread_tbl_full_vaf = spread_tbl_full %>%
    dplyr::select(-coverage) %>%
    spread(key = sample, value=vaf, fill=0) %>%
    filter(alt!='INS', alt!='-', minvaf>minvaf_filter)
  spread_tbl_full_altreads = spread_tbl_full %>%
    mutate(altreads = round(vaf*coverage)) %>%
    dplyr::select(-vaf, -coverage) %>%
    spread(key = sample, value=altreads, fill=0) %>%
    filter(alt!='INS', alt!='-', minvaf>minvaf_filter)
  spread_tbl_full_depth = spread_tbl_full %>%
    dplyr::select(-vaf) %>%
    spread(key = sample, value=coverage, fill=0) %>%
    filter(alt!='INS', alt!='-', minvaf>minvaf_filter)
  spread_tbl_full = spread_tbl_full %>%
    mutate(is_var = 1, vaf=median_vaf) %>% dplyr::select(-median_vaf, -coverage) %>%
    spread(key = sample, value=is_var, fill=0) %>%
    filter(alt!='INS', alt!='-', minvaf>minvaf_filter)
  #spread_tbl = group_by(spread_tbl_full, shared_id) %>% mutate(shared_num_var=n()) %>% ungroup() 
  #spread_tbl = group_by(spread_tbl_full, donor, shared_id, shared_num_samples) %>% summarise(shared_num_var=n()) %>% ungroup() 
  spread_tbl = spread_tbl_full %>%
    group_by(shared_id) %>% mutate(shared_num_var = n()) %>% ungroup() %>%
    dplyr::select(-chrom, -pos, -ref, -alt, -donor, -vaf, -minvaf) %>% distinct()
  list(spread_tbl, spread_tbl_full, spread_tbl_full_vaf, spread_tbl_full_altreads, spread_tbl_full_depth)
}

plot_VAF_hist_per_donor <- function(var_calls_tbl, outpath=NA) {
  # Plots histogram of VAF distribution for each donor in var_call_tbl.
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 5-NOV-2017
  #
  # Args
  #   var_calls_tbl: table of variant calls. Required columns: donor, sample, chrom, ref, alt
  #   outpath (default: NA): path to save plot. If NA, then plot is not saved.
  # Returns:
  #   list(
  #     gg_handle: handle of ggplot object
  # )
  p = ggplot(var_calls_tbl, aes(x=vaf, fill=donor)) + geom_histogram(binwidth=0.01, alpha=0.7) +
    facet_grid(donor~., scales='free_y') +
    labs(x='VAF', y='Count', title='VAF distribution across all samples, binwidth = 0.01', fill='Donor') +
    theme(legend.position="top", legend.justification = 'right', legend.box = "horizontal")
  if(!is.na(outpath)) {
    p + ggsave(outpath, width=20, height=14, units='cm')
  }
  
  return(list(gg_handle=p))
}

plot_VAF_hist_per_sample <- function(var_calls_tbl, outpath=NA) {
  # Plots histogram of VAF distribution for each sample in var_call_tbl.
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 5-NOV-2017
  #
  # Args
  #   var_calls_tbl: table of variant calls. Required columns: donor, sample, chrom, ref, alt
  #   outpath (default: NA): path to save plot. If NA, then plot is not saved.
  # Returns:
  #   list(
  #     gg_handle: handle of ggplot object
  # )
  p = ggplot(var_calls_tbl, aes(x=vaf, fill=sample)) + geom_histogram(binwidth = 0.01, color='black', size=0.2) +
    facet_grid(sample~., scales='free_y') +
    theme(strip.text.y = element_text(angle = 0)) +
    theme(legend.position="none") +
    labs(x='VAF', y='Count', title='VAF distribution in each sample')
  if(!is.na(outpath)) {
    p + ggsave(outpath, width=14, height=20, units='cm')
  }
  
  return(list(gg_handle=p))
}

plot_coverage_hist_per_sample <- function(var_calls_tbl, outpath=NA) {
  # Plots histogram of coverage distribution for each sample in var_call_tbl.
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 5-NOV-2017
  #
  # Args
  #   var_calls_tbl: table of variant calls. Required columns: donor, sample, chrom, ref, alt
  #   outpath (default: NA): path to save plot. If NA, then plot is not saved.
  # Returns:
  #   list(
  #     gg_handle: handle of ggplot object
  # )
  p = ggplot(var_calls_tbl, aes(x=coverage, fill=sample)) + geom_histogram(binwidth = 1, color='black', size=0.2) +
    facet_grid(sample~., scales='free_y') +
    theme(strip.text.y = element_text(angle = 0)) +
    xlim(0,150) +
    theme(legend.position="none") +
    labs(x='Depth of coverage', y='Count', title='Coverage at variant call sites')
  if(!is.na(outpath)) {
    p + ggsave(outpath, width=14, height=20, units='cm')
  }
  
  return(list(gg_handle=p))
}

plot_vars_per_groupsize <- function(var_calls_tbl, outpath=NA) {
  # Bar plot of number of variants per groupsize (x=groupsize, y=num variants)
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 5-NOV-2017
  #
  # Args
  #   var_calls_tbl: table of variant calls. Required columns: donor, sample, chrom, ref, alt, groupsize
  #   outpath (default: NA): path to save plot. If NA, then plot is not saved.
  # Returns:
  #   list(
  #     gg_handle: handle of ggplot object
  # )
  p = ggplot(var_calls_tbl, aes(x=groupsize, fill=donor)) + geom_bar(alpha=0.7) +
    facet_grid(.~donor, scales='free', space = 'free') +
    scale_x_continuous(breaks = c(1:length(unique(var_calls_tbl$sample)))) +
    labs(x='Number of dissections sharing same mutation', 
         y='Number of shared mutations', title='Group sizes of shared mutations', fill='Donor')
  if(!is.na(outpath)) {
    p + ggsave(outpath, width=25, height=16, units='cm')
  }
  
  return(list(gg_handle=p))
}

plot_vars_per_mut_group <- function(var_calls_tbl, outpath=NA, filter_private=T) {
  # Bar plots of number of variants per sample and mutation group 
  # (x=number of variants, y=samples, fill=mutation group, eg. T>A, C>T etc)
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 5-NOV-2017
  #
  # Args
  #   var_calls_tbl: table of variant calls. Required columns: donor, sample, chrom, ref, alt, groupsize
  #   outpath (default: NA): path to save plots. If NA, then plot is not saved.
  #   filter_private (default: TRUE): if TRUE, only plots private variants
  # Returns:
  #   list(
  #     gg_handle: handle of ggplot object
  # )
  plt_dat = mutate(var_calls_tbl, mut_type = paste(ref, '>', alt, sep=''))
  if(filter_private) {
    plt_dat = plt_dat %>% filter(groupsize==1)
  }
  plt_dat = mutate(plt_dat, mut_type = ifelse(mut_type == 'G>T', 'C>A', 
                                              ifelse(mut_type == 'G>C', 'C>G', 
                                                     ifelse(mut_type == 'G>A', 'C>T', 
                                                            ifelse(mut_type == 'A>T', 'T>A', 
                                                                   ifelse(mut_type == 'A>G', 'T>C', 
                                                                          ifelse(mut_type == 'A>C', 'T>G', mut_type)))))))
  
  p <- ggplot(plt_dat, aes(x=sample, fill=mut_type)) + 
    coord_flip() +
    facet_grid(donor~., scales = 'free_y', space = 'free_y') +
    labs(fill='Donor', y='Number of mutations', x='Dissection', title='Mutations per dissection') +
    scale_fill_manual(values = c('blue', 'black', "red", 'gray', "green", "pink")) +
    theme(legend.position="top", legend.justification = 'right', legend.box = "horizontal") +
    guides(colour = guide_legend(nrow = 1)) +
    geom_bar(position = position_stack(), color='black', size=0.2)
  
  if(!is.na(outpath)) {
    p + ggsave(outpath, width=25, height=16, units='cm')
  }
  
  return(list(gg_handle=p))
}

plot_main_mut_sharing_groups <- function(spread_summary, outpath=NA, min_num_var=5) {
  # Plots groups of variants that share the same mutations.
  # (x: number of samples sharing same mutation, y: number of shared variants)
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 11-DEC-2017
  #
  # Args
  #   spread_summary: table of mutation sharing groups. Required columns: shared_id, shared_num_samples, shared_num_var
  #   outpath (default: NA): path to save plots. If NA, then plot is not saved.
  #   min_num_var (default: 5): minimum number of variants within a mutation sharing group to be included in plot
  # Returns:
  #   list(
  #     gg_handle: handle of ggplot object
  # )
  plt_dat = filter(spread_summary, shared_num_var>min_num_var, shared_num_samples>1)
  p = ggplot(plt_dat, aes(x=shared_num_samples, y=shared_num_var)) + geom_point(alpha=0.5, color='red') +
    scale_x_continuous(breaks = c(1:max(spread_summary$shared_num_samples)), minor_breaks=NULL) +
    labs(x="Number of samples sharing same mutation", y="Number of shared variants") +
    scale_y_log10() +
    xlim(0, max(plt_dat$shared_num_samples)+1) +
    geom_text(aes(label = shared_id, x=shared_num_samples), hjust=0, size=3, angle=45, nudge_x = 0.02, nudge_y = 0.02, check_overlap=T)
  
  if(!is.na(outpath)) {
    p + ggsave(outpath, width=40, height=30, units='cm')
  }
  
  return(list(gg_handle=p))
}

plot_shared_mut_boxplot <- function(spread_all_vaf, outpath=NA, min_num_var=50) {
  # Plots groups of variants that share the same mutations.
  # (x: number of samples sharing same mutation, y: number of shared variants)
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 11-DEC-2017
  #
  # Args
  #   spread_all: table of variant calls. Required columns: chrom, pos, ref, alt, donor, shared_id, shared_num_samples
  #   outpath (default: NA): path to save plots. If NA, then plot is not saved.
  #   min_num_var (default: 50): minimum number of variants within a mutation sharing group to be included in plot
  # Returns:
  #   list(
  #     gg_handle: handle of ggplot object
  # )
  plt_dat = group_by(spread_all_vaf, shared_id) %>% mutate(shared_num_var = n()) %>% ungroup() %>%
    filter(shared_num_var > min_num_var, shared_num_samples>1) %>%
    gather(key='sample', value='vaf', which(grepl('lo', names(spread_all_vaf)))) %>%
    filter(vaf>0)
  
  p = ggplot(plt_dat, aes(x=shared_id, y=vaf)) + 
    geom_boxplot() +
    geom_jitter(width = 0.2, size=0.2, alpha=0.2, color='red') +
    coord_flip() +
    facet_grid(factor(shared_num_samples)~., space='free', scales = 'free') +
    labs(title='VAF distribution of variants within dissection groups', y='VAF', x='')
  
  if(!is.na(outpath)) {
    p + ggsave(outpath, width=20, height=30, units='cm')
  }
  
  return(list(gg_handle=p))
}

plot_heatmap_tree <- function(spread_summary, outdir, min_num_var=100, max_iter=10) {
  # Plots heatmaps of shared mutation groups. 
  # Sample pairs that share mutations are the basis of a heatmap
  # of sample triplets that share mutations, which then in turn serve as a basis
  # for sample quartets etc, until no mutations are left.
  # (x: parental sample group, y: new sample group)
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 11-DEC-2017
  #
  # Args
  #   spread_summary: table of variant calls. Required columns: shared_id, shared_num_samples, shared_num_var
  #   outdir: directory in which to save plots. Is required for this function.
  #   min_num_var (default: 100): minimum number of variants within a mutation sharing group to be included in plot
  #   max_iter (default: 10): maximum number of loop iterations to plot heatmaps
  # Returns:
  #   nothing
  plt_tbl = NA
  for(k in c(1:max_iter)) {
    if(is.na(plt_tbl)) {
      start_tbl = filter(spread_summary, shared_num_samples==k, shared_num_var>min_num_var)
    } else {
      start_tbl = filter(spread_summary, shared_id %in% unique(plt_tbl$shared_id))
    }
    sample_cols = which(grepl('lo', names(start_tbl)))
    
    target_shared_id = list()
    for(i in c(1:dim(start_tbl)[1])) {
      print('Target table:')
      print(start_tbl$shared_id[i])
      this_group = start_tbl[i, sample_cols]
      this_mat = t(replicate(dim(spread_summary)[1], as.list(this_group)))
      this_tbl = as.matrix(spread_summary[,sample_cols])
      parental_idx = spread_summary$shared_id[which(rowSums(apply(apply(this_tbl, 2, as.numeric) - apply(this_mat, 2, as.numeric), 2, abs))==1)]
      this_target = filter(spread_summary, shared_id %in% parental_idx, shared_num_var>10, shared_num_samples == sum(this_group)+1)$shared_id
      if(length(this_target)>0) {
        this_out_tbl = filter(spread_summary, shared_id %in% this_target) %>% dplyr::select(shared_id, shared_num_var)
        this_out_tbl$daughter_id = start_tbl$shared_id[i]
        target_shared_id[[i]] = this_out_tbl
        print(this_target)
        plt_tbl = tbl_df(do.call('rbind', target_shared_id)) %>% distinct()
        
      }
    }
      if(dim(plt_tbl)[1]>0) {
        # Generate any combination of shared_id and daughter id
        all_combs = tbl_df(matrix(data = 0, nrow=length(unique(plt_tbl$shared_id)), ncol=length(unique(plt_tbl$daughter_id))))
        colnames(all_combs)<-unique(plt_tbl$daughter_id)
        all_combs$shared_id=unique(plt_tbl$shared_id)
        all_combs = all_combs %>%
          gather(key='daughter_id', value='is_shared', which(grepl('lo', names(all_combs)))) %>% 
          dplyr::select(-is_shared)
        all_combs = dplyr::setdiff(all_combs, plt_tbl %>% dplyr::select(shared_id, daughter_id)) %>%
          mutate(shared_num_var=NA)
        plt_tbl = do.call('rbind', list(plt_tbl, all_combs))
        
        p = ggplot(plt_tbl, aes(x=daughter_id, y=shared_id, fill=shared_num_var)) + geom_tile(color='white') +
          scale_fill_gradient2(low = 'white', mid='lightblue', high='red', na.value='grey90', midpoint=max(plt_tbl$shared_num_var, na.rm=T)/2) +
          theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
          geom_text(aes(label=shared_num_var), size=2, vjust='center') + 
          scale_x_discrete(expand=c(0,0)) +
          scale_y_discrete(expand=c(0,0)) +
          labs(x='Daugher ID', y='Parental ID', title='Dissection groups sharing mutations', fill='Mutations') +
          theme(legend.position="top", legend.justification = 'right', legend.box = "horizontal")
        p + ggsave(file.path(outdir, paste0('phylogenetic_heatmap_', sum(this_group)+1, '.pdf')), width=20, height=20, units='cm')
      }
  }
}

plot_pair_heatmap <- function(spread_summary, outpath=NA, min_num_var=10) {
  # Plots groups of variants that share the same mutations.
  # (x: number of samples sharing same mutation, y: number of shared variants)
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 11-DEC-2017
  #
  # Args
  #   spread_summary: table of variant calls. Required columns: shared_id, shared_num_samples, shared_num_var
  #   outpath (default: NA): path to save plots. If NA, then plot is not saved.
  #   min_num_var (default: 10): minimum number of variants in pair to display text label
  # Returns:
  #   list(
  #     gg_handle: handle of ggplot object
  # )
  # Generate every possible pair of samples
  all_pairs = tbl_df(as.data.frame(t(combn(x=names(spread_summary)[grepl('lo', names(spread_summary))], m=2))))
  all_pairs = all_pairs %>%
    mutate(shared_num_samples = 2, shared_num_var = 0, sample1 = V1, sample2 = V2, shared_id = paste(V1, V2, sep='_')) %>%
    dplyr::select(-V1, -V2)
  
  # Generate 
  pairs1 = spread_summary %>%
    filter(shared_num_samples==2) %>%
    gather(key='sample_id', value='is_shared', which(grepl('lo', names(spread_summary)))) %>%
    filter(is_shared==1) %>%
    group_by(shared_id, shared_num_samples) %>%
    summarise(sample1 = first(sample_id), sample2 = last(sample_id), shared_num_var = n()) %>% ungroup() %>% arrange(shared_id)
  
  plt_dat = do.call('rbind', list(pairs1, filter(all_pairs, !(shared_id %in% pairs1$shared_id)))) %>%
    mutate(shared_num_var = ifelse(shared_num_var==0, NA, shared_num_var))
  
  p=ggplot(plt_dat, aes(x=sample1, y=sample2, fill=shared_num_var)) + geom_tile(color='black') +
    theme_light() +
    scale_fill_gradient2(low = 'white', mid='lightblue', high='red', na.value='grey90', midpoint=max(plt_dat$shared_num_var, na.rm=T)/2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    geom_text(data = plt_dat %>% filter(shared_num_var>min_num_var), aes(label=shared_num_var), vjust='center', size=2) +
    labs(x='Sample 1', y='Sample 2', title='Pairs of dissections sharing mutations', fill='Mutations') +
    theme(legend.position="top", legend.justification = 'right', legend.box = "horizontal")
  
  if(!is.na(outpath)) {
    p + ggsave(outpath, width=20, height=20, units='cm')
  }
  
  return(list(gg_handle=p))
}

filter_muts_to_sharing_group <- function(spread_all, sample_groups, min_num_mut = 0, only_shared_by_all=F) {
  # Extract all those variants shared by the given samples
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 11-DEC-2017
  #
  # Args
  #   spread_all: table of variant calls. Required columns: ref, alt, donor, shared_id, shared_num_samples
  #   sample_groups: groups of samples for which to plot their shared mutations.
  #       example: list(c('lo001', 'lo002), c('lo001', lo005'))
  #   min_num_mut : minimal number of mutations required for group to be returned
  #   only_shared_by_all: only include mutations that are shared by ALL the samples in a given sample group.
  #       By default, any mutation affecting even just one of the samples will be included.
  # Returns:
  #   list(
  #     var_calls_lst => list of variant call tables
  # )
  
  var_calls_lst = list()
  for(i in c(1:length(sample_groups))) {
    this_group = sample_groups[[i]]
    group_name = paste(this_group, collapse='_')
    print(group_name)
    rm_samples = setdiff(names(spread_all)[which(grepl('lo', names(spread_all)))], this_group)
    context_tbl = spread_all
    context_tbl$sum_rm_samples = rowSums(context_tbl[,rm_samples])
    
    context_tbl = filter(context_tbl, sum_rm_samples==0) %>%
      group_by(shared_id) %>% mutate(shared_num_var = n()) %>% ungroup() %>%
      mutate(sample_group = group_name)
    if(only_shared_by_all == T) {
      context_tbl = filter(context_tbl, shared_num_samples==length(this_group))
    }
    if(dim(context_tbl)[1] < min_num_mut) { next }
    var_calls_lst[[i]] = context_tbl
  }
  return(var_calls_lst)
}

plot_signature_contribution_bars <- function(var_calls_lst, outpath=NA, min_num_mut = 0, par_width=14, par_height=20, 
                                             order_by_group_size=F) {
  # Plots contribution of signatures to each sample group
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 13-DEC-2017
  #
  # Args
  #   var_calls_lst: list of variant call tables. 
  #       The plot will contain one bar per distinct element of the 'sample_group' column. 
  #       Facets will be plotted per 'sample_collection' column.
  #   outpath (default: NA): path to save plots. If NA, then plot is not saved.
  #   min_num_mut : minimal number of mutations required for group to be plotted
  #   par_width (default: 14): width in centimeters
  #   par_height (default: 20): width in centimeters
  #   order_by_group_size (default: F): if TRUE, then the sample groups are ordered according to their length (ie. how many samples they contain)
  # Returns:
  #   list(
  #     gg_handle: handle of group weights ggplot
  #     weights_tbl: table with weights of each signature in each sample group
  # )
  weights_lst = list()
  
  # Build single table from list of calls
  all_calls = do.call('rbind', var_calls_lst) %>% group_by(sample_group) %>%
    mutate(is_rm = ifelse(n() < min_num_mut, T, F)) %>% ungroup() %>%
    filter(is_rm==F)
  if(!('sample_collection' %in% names(all_calls))) {
    all_calls$sample_collection = 'Sample collection'
  }
  sample_groups = unlist(lapply(var_calls_lst, function(x) {unique(x$sample_group)}))
  num_mut_lst = table(all_calls$sample_group)
  
  # Get mutation contexts
  mut_contexts = get_mutation_contexts(all_calls %>% mutate(sample = sample_group), generate_plot=F)
  
  # Get signature contributions
  for(i in c(1:length(sample_groups))) {
    group_name = sample_groups[[i]]
    sample_collection = unique(filter(all_calls, sample_group == group_name)$sample_collection)
    print(group_name)
    
    sigs = whichSignatures(tumor.ref = mut_contexts$trinuc, sample.id = group_name, signatures.ref = signatures.cosmic, contexts.needed = TRUE)
    weights = sigs$weights
    weights$sample_collection = sample_collection
    weights_lst[[group_name]] = weights
  }
  
  weights_df = do.call('rbind', weights_lst)
  weights_tbl = tbl_df(weights_df) %>% mutate(group_name=row.names(weights_df))
  weights_tbl = weights_tbl %>%  
    gather(key='Signature', value='Weight', which(grepl('Signature.', names(weights_tbl)))) %>%
    filter(Weight!=0) %>%
    mutate(Signature = gsub('Signature.', '', Signature)) %>%
    mutate(Signature = factor(Signature, levels=c(30:1), ordered=T)) %>%
    mutate(Sample_group = factor(group_name, levels=rev(sample_groups), ordered=T))
  if(order_by_group_size) {
    weights_tbl = weights_tbl %>% mutate(Sample_group = factor(group_name, levels=unique(sample_groups[order(nchar(sample_groups), decreasing=T)]), ordered=T)) 
  }
  weights_tbl$num_mut = unlist(as.list(num_mut_lst[weights_tbl$group_name]))
  #weights_tbl
  
  p = ggplot(weights_tbl, aes(x=Sample_group, y=Weight, fill=Signature)) + geom_bar(stat='identity', color='black') +
    scale_fill_manual(values = get_deconstructSigs_colors()) +
    coord_flip() +
    theme(legend.position="none") +
    geom_text(size = 3, position = position_stack(vjust = 0.5), aes(label=Signature)) +
    geom_label(aes(label=num_mut), y=1.1, fill='white') +
    ylim(0,1.2) +
    facet_grid(sample_collection ~ .)
  
  if(!is.na(outpath)) {
    p + ggsave(outpath, width=par_width, height=par_height, units='cm')
  }
  
  return(list(gg_handle=p, weights_tbl=weights_tbl))
}

get_var_calls_lst_all_combs <- function(spread_all_vaf, sample_ints, min_num_mut = 100) {
  # Calculates all possible combinations of sample_ints and extracts variants shared by each combination.
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 14-DEC-2017
  #
  # Args
  #   spread_all_vaf: table of variant calls. Required columns: ref, alt, donor, shared_id, shared_num_samples
  #   sample_ints: vector of samples to include in plot. Samples are given as integers. They will then
  #     be converted to actual sample names. ie. 1 => lo001
  #   min_num_mut (default: 100): minimal number of mutations required for group to be plotted
  # Returns:
  #   list(
  #     vaf_calls_lst: list of variant call tables
  # )
  samples = int_to_sample_id(sample_ints)
  all_combs = list()
  for(i in c(1:length(samples))) {
    this_combs = t(combn(samples, i))
    all_combs = c(all_combs, split(this_combs, seq(nrow(this_combs))))
  }
  
  print('Filtering variant calls.')
  var_calls_lst = filter_muts_to_sharing_group(spread_all_vaf, all_combs, only_shared_by_all = T, min_num_mut = min_num_mut)
  
  return(var_calls_lst = var_calls_lst)
}

plot_vaf_boxplot_from_tables <- function(var_calls_lst, outpath=NA, plot_jitter=T, par_width=14, par_height=20) {
  # Boxplot of VAFs in each table of var_calls_lst
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 14-DEC-2017
  #
  # Args
  #   var_calls_lst: list of variant call tables
  #   outpath (default: NA): path to save plots. If NA, then plot is not saved.
  #   plot_jitter (default: TRUE): if TRUE, then geom_jitter of each variant will be plotted on top of boxplot
  #   par_width (default: 14): width in centimeters
  #   par_height (default: 20): width in centimeters
  # Returns:
  #   list(
  #     gg_handle: handle of group weights ggplot
  # )
  
  all_calls = do.call('rbind', var_calls_lst) 
  sample_groups = unlist(lapply(var_calls_lst, function(x) {unique(x$sample_group)}))
  all_calls = all_calls %>%
    gather(key='sample', value='vaf', which(grepl('lo', names(all_calls)))) %>%
    mutate(shared_id = factor(shared_id, levels=rev(sample_groups), ordered=T)) %>%
    filter(vaf>0)
  p = ggplot(all_calls, aes(x=shared_id, y=vaf)) + geom_boxplot() +
    coord_flip() +
    labs(x='Sample group', y='VAF')
  if(plot_jitter) {
    p = p + geom_jitter(width=0.2, size=0.2, color='red', alpha=0.3)
  }
  
  if(!is.na(outpath)) {
    p + ggsave(outpath, width=par_width, height=par_height, units='cm')
  }
  
  return(list(gg_handle=p))
}

plot_mut_context_of_samples <- function(var_calls_tbl, outdir, plot_suffix='', context_plt_height=10) {
  # Plots mutation contexts and signature assignment of each sample.
  # (Note: if sample is set to donor column, then the function is performed on all samples of a donor at once)
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 14-DEC-2017
  #
  # Args
  #   var_calls_tbl: variant calls
  #   outdir: directory to save plots
  #   plot_suffix (default: ''): allows to add a suffix to the plot name.
  #   context_plt_height (default: 14): height of mutation context bar plot in centimeters
  # Returns:
  #   nothing
  mut_contexts = get_mutation_contexts(var_calls_tbl)
  mut_contexts$p + ggsave(file.path(outdir, paste0('mutation_contexts', plot_suffix, '.pdf')), width=20, height=context_plt_height, units='cm')
  
  samples = unique(var_calls_tbl$sample)
  for(sample in samples) {
    sigs = whichSignatures(tumor.ref = mut_contexts$trinuc, sample.id = sample, signatures.ref = signatures.cosmic, contexts.needed = TRUE)
    pdf(file.path(outdir, paste0('mutation_contexts_sigpie_', sample, '.pdf')), width = 2, height = 2)
    par(cex=0.5)
    par(col.main='white')
    makePie(sigs)
    par(col.main='black')
    title(main=paste0(sample, '\nNum mut: ', sum(mut_contexts$trinuc[sample,])))
    dev.off()
  }
  
  return(p = mut_contexts$p)
}

plot_mut_context_of_groups <- function(spread_all, sample_groups, outdir, min_num_mut = 0, only_shared_by_all=F, keep_order=F) {
  # Plots mutation contexts and signature assignment of the union of mutations within the given groups.
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 11-DEC-2017
  #
  # Args
  #   spread_all: table of variant calls. Required columns: ref, alt, donor, shared_id, shared_num_samples
  #   sample_groups: groups of samples for which to plot their shared mutations.
  #       example: list(c('lo001', 'lo002), c('lo001', lo005'))
  #   outdir : directory to save plots. Required.
  #   min_num_mut : minimal number of mutations required for group to be plotted
  #   only_shared_by_all: only include mutations that are shared by ALL the samples in a given sample group.
  #       By default, any mutation affecting even just one of the samples will be included.
  #   keep_order: keep the order of sample groups as it was specified in the sample_groups argument. Else will group by length of ID.
  # Returns:
  #   list(
  #     gg_handle: handle of group weights ggplot
  # )
  weights_lst = list()
  num_mut_lst = list()
  
  # Get mutation call tables
  context_tbl_lst = filter_muts_to_sharing_group(spread_all, sample_groups, min_num_mut, only_shared_by_all)
  
  for(i in c(1:length(sample_groups))) {
    context_tbl = context_tbl_lst[[i]]
    group_name = unique(context_tbl$sample_group)
    if(nchar(group_name)>30) {
      group_name = substr(group_name, 1, 30)
    }
    print(group_name)
    context_tbl$sample = group_name
    
    num_mut_lst[group_name] = dim(context_tbl)[1]
    mut_contexts = get_mutation_contexts(context_tbl)
    mut_contexts$p + ggsave(file.path(outdir, paste0('mutation_contexts_', group_name, '.pdf')), width=20, height=10, units='cm')
    
    sigs = whichSignatures(tumor.ref = mut_contexts$trinuc, sample.id = group_name, signatures.ref = signatures.cosmic, contexts.needed = TRUE)
    weights_lst[[group_name]] = sigs$weights
    
    pdf(file.path(outdir, paste0('mutation_contexts_', group_name, '_sigpie.pdf')), width = 2, height = 2)
    par(cex=0.5)
    par(col.main='white')
    makePie(sigs)
    par(col.main='black')
    title(main=paste0(group_name, '\nNum mut: ', dim(context_tbl)[1]))
    dev.off()
  }
  
  weights_df = do.call('rbind', weights_lst)
  weights_tbl = tbl_df(weights_df) %>% mutate(group_name=row.names(weights_df))
  weights_tbl = weights_tbl %>%  
    gather(key='Signature', value='Weight', which(grepl('Signature.', names(weights_tbl)))) %>%
    filter(Weight!=0) %>%
    mutate(Signature = gsub('Signature.', '', Signature)) %>%
    mutate(Signature = factor(Signature, levels=c(30:1), ordered=T))
  if(keep_order==F) {
    weights_tbl = weights_tbl %>% 
      mutate(Sample_group = factor(group_name, levels=unique(group_name[order(nchar(group_name), decreasing=T)]), ordered=T)) 
  } else {
    sample_group_levels = rev(unlist(lapply(sample_groups, paste, collapse='_')))
    weights_tbl = weights_tbl %>% 
      mutate(Sample_group = factor(group_name, levels=sample_group_levels, ordered=T)) 
  }
  weights_tbl$num_mut = unlist(num_mut_lst[weights_tbl$group_name])
  #weights_tbl
  
  p = ggplot(weights_tbl, aes(x=Sample_group, y=Weight, fill=Signature)) + geom_bar(stat='identity', color='black') +
    scale_fill_manual(values = get_deconstructSigs_colors()) +
    coord_flip() +
    theme(legend.position="none") +
    geom_text(size = 3, position = position_stack(vjust = 0.5), aes(label=Signature)) +
    geom_label(aes(label=num_mut), y=1.1, fill='white') +
    ylim(0,1.2) +
    ggsave(file.path(outdir, paste0('mutation_contexts_', paste(unique(unlist(sample_groups)), collapse='_'), '_groupweights.pdf')), width=20, height=14, units='cm')
  
  return(list(gg_handle=p, weights_lst=weights_lst))
}

plot_ApTpN_strand_bias <- function(context_df_stranded, outpath=NA) {
  # Bar plot of number of variants per groupsize (x=groupsize, y=num variants)
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 15-NOV-2017
  #
  # Args
  #   context_df_stranded: table with counts of variants at untranscribed and transcribed sites. 
  #       required columns: context, mut_type, strand, sample, counts
  #   outpath (default: NA): path to save plot. If NA, then plot is not saved.
  # Returns:
  #   list(
  #     gg_handle: handle of ggplot object
  # )
  strand_analysis = mutate(context_df_stranded, context=as.character(context)) %>%
    filter(grepl('A\\[T>', context)) %>% dplyr::select(-context) %>%
    group_by(sample, strand) %>% summarise(counts=sum(counts)) %>% ungroup()
  
  p_raw = ggplot(strand_analysis, aes(x=sample, y=counts, fill=strand)) + 
    coord_flip() +
    labs(x='Sample', y='Counts', fill='Strand', title='ApTpN strand bias')
  p = p_raw + geom_bar(stat='identity', position = "fill") +
    scale_y_continuous(labels = percent_format())
  if(!is.na(outpath)) {
    p + ggsave(outpath, width=25, height=16, units='cm')
  }
  
  return(list(gg_handle=p, strand_analysis=strand_analysis, p_raw=p_raw))
}

int_to_sample_id <- function(int_vec, len_id=3, prefix='lo') {
  # Convert integers to sample_id
  # Example: int_vec = c(1) returns c('lo001')
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 11-DEC-2017
  #
  # Args
  #   int_vec: vector of integers, eg. c(1,2,3)
  #   len_id: if sample id is of format 'lo001', then 3, if 'lo0001', then 4 etc
  # Returns:
  #   vector of samples, eg. c('lo001', 'lo002', 'lo003')
  paste0(prefix, sprintf(paste0("%0", len_id, "d"), int_vec))
}

load_phys_coord_tbl <- function(path_phys_coord, var_calls_tbl) {
  # Convert integers to sample_id
  # Example: int_vec = c(1) returns c('lo001')
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 11-DEC-2017
  #
  # Args
  #   path_phys_coord: path to CSV storing physical coordinates of each dissection
  #   var_calls_tbl: variant call table, necessary for removing physical coordinates of dissections that are not contained in variant calls
  # Returns:
  #   phys_coord_tbl
  phys_coord_tbl = tbl_df(fread(path_phys_coord)) %>%
    filter(!is.na(x+y+z))
  # Remove samples for which mutation data is not available
  phys_coord_tbl = filter(phys_coord_tbl, sample %in% unique(var_calls_tbl$sample))
  # Rename samples
  phys_coord_tbl$sample = unlist(lapply(phys_coord_tbl$sample, function(x) {
    strsplit(x, 'lo')[[1]][2]
  }))
  
  return(phys_coord_tbl)
}

plot_phys_coord_dendrogram <- function(phys_coord_tbl, outpath, width=6, height=5) {
  # Plot dendrogram of physical distances
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 14-DEC-2017
  #
  # Args
  #   phys_coord_tbl: tibble storing physical coordinates of dissections
  #   outpath: path to store plot
  # Returns:
  #   eucl_dist: euclidean distances
  
  dist_tbl = as.data.frame(dplyr::select(phys_coord_tbl, x, y, z))
  rownames(dist_tbl)<-phys_coord_tbl$sample
  
  # prepare hierarchical cluster
  eucl_dist = dist((dist_tbl), method='euclidean')
  hc = hclust(eucl_dist, method='average')
  # very simple dendrogram
  pdf(outpath, width=width, height=height)
  plot(hc, hang=-1, 
       main = 'Clustering of physical distances between dissections',
       xlab = 'Dissections',
       ylab = 'Distance [a.u.]')
  dev.off()
  
  return(eucl_dist)
}

plot_phys_coord <- function(phys_coord_tbl, outpath, width=14/2.54, height=7/2.54) {
  # Plot physical coordinates of dissections
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 14-DEC-2017
  #
  # Args
  #   phys_coord_tbl: tibble storing physical coordinates of dissections
  #   outpath: path to store plot
  # Returns:
  #   nothing
  
  pdf(outpath, width=width, height=height)
  par(mar=c(0,0,0,0))
  plot(phys_coord_tbl$x, phys_coord_tbl$y, xaxt='n', yaxt='n', ann=FALSE, type='p', pch=20, col='green', xlim=c(0,140), ylim=rev(range(seq(70,0))), cex=0.5)
  text(phys_coord$x+1, phys_coord$y, labels=phys_coord$sample, adj=0, cex=0.3, col='black')
  dev.off()
}

plot_phys_coord_num_mut <- function(phys_coord_tbl, var_calls_tbl, outpath, width=14/2.54, height=7/2.54) {
  # Plot physical coordinates of dissections and number of mutations at each coordinate
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 14-DEC-2017
  #
  # Args
  #   phys_coord_tbl: tibble storing physical coordinates of dissections
  #   var_calls_tbl: variant call table
  #   outpath: path to store plot
  # Returns:
  #   nothing
  
  num_muts = var_calls_tbl %>%
    group_by(sample) %>% summarise(num_mut = n()) %>% ungroup()
  num_muts$sample = unlist(lapply(num_muts$sample, function(x) { strsplit(x, '_lo')[[1]][2] }))
  plt_tbl = left_join(phys_coord_tbl, num_muts, by='sample')
  
  lwidth = sqrt(plt_tbl$num_mut)*0.5
  rainbow_cols = rainbow(max(plt_tbl$num_mut), start=0, end=1, alpha=0.6)
  
  pdf(outpath, width=width, height=height)
  layout(t(1:2),widths=c(1,1))
  par(mar=c(0,0,0,0))
  plot(plt_tbl$x, plt_tbl$y, xaxt='n', yaxt='n', ann=FALSE, type='p', pch=21, xlim=c(0,max(plt_tbl$x)*1.1), ylim=rev(range(seq(70,0))), cex=0.5, col=rainbow_cols[plt_tbl$num_mut], lwd=lwidth, bty='n')
  points(plt_tbl$x, plt_tbl$y, xaxt='n', yaxt='n', ann=FALSE, type='p', pch=19, xlim=c(0,max(plt_tbl$x)*1.1), ylim=rev(range(seq(70,0))), cex=0.5, col='black')
  text(plt_tbl$x+1, plt_tbl$y, labels=plt_tbl$sample, adj=0, cex=0.5, col='black')
  
  par(mar=c(5,3,5,10))
  image(y=1:max(plt_tbl$num_mut), z=t(1:max(plt_tbl$num_mut)), col=rainbow_cols[1:max(plt_tbl$num_mut)], axes=FALSE, cex.main=.8)
  axis(4,cex.axis=0.8,mgp=c(0,.5,0))
  
  dev.off()
}

plot_phys_coord_shared <- function(phys_coord, spread_all, num_samples, outpath, param_minvaf=NA, min_shared_var=20, param_max_shared=NA, label_coords=TRUE) {
  # Plot physical coordinates of dissections. Color by number of shared dissections.
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 15-NOV-2017
  #
  # Args
  #   phys_coord: tibble storing coordinate data
  #   spread_all: tibble storing shared mutation data
  #   num_samples: maximum number of sharedness to consider
  #   outpath: path to store plot
  #   param_minvaf (default: NA): minimum vaf to consider for a dissection / shared dissections to be plotted
  #   min_shared_vaf (default: 20): minimum number of shared variants
  #   param_max_shared (default: NA): maximum number of shared samples, important to scale colors
  #   label_coords (default: TRUE): whether or not to label the coordinates with numbers
  # Returns:
  #   list(
  #     polygon_tbl
  # )
  print(param_minvaf)
  print(min_shared_var)
  col_gradient_fx = colorRampPalette(c("yellow", "red"))
  col_gradient_fx = colorRampPalette(c("blue", "green", "yellow", "orange", "red"))
  col_gradient = col_gradient_fx(5)
  
  pdf(outpath, width=21/2.54, height=7/2.54)
  layout(t(1:2),widths=c(2,1))
  par(mar=c(0,0,0,0))
  plot(phys_coord$x, phys_coord$y, xaxt='n', yaxt='n', ann=FALSE, type='p', pch=20, col='black', xlim=c(0,max(phys_coord$x)*1.1), ylim=rev(range(seq(max(phys_coord$y)*1.1,0))), cex=0.5, bty="n")
  
  polygon_tbl = spread_all
  if(!is.na(param_minvaf)) {
    polygon_tbl = filter(polygon_tbl, vaf>param_minvaf)
  }
  polygon_tbl = polygon_tbl %>%
    group_by(shared_id) %>% mutate(shared_num_var = n()) %>% ungroup() %>%
    dplyr::select(-donor, -chrom, -pos, -ref, -alt, -vaf, -minvaf) %>% distinct() %>%
    filter(shared_num_samples>1, shared_num_var>min_shared_var)
  polygon_tbl = polygon_tbl %>%
    gather(key='sample', value='is_contained', names(polygon_tbl)[grepl('lo', names(polygon_tbl))]) %>% filter(is_contained==1) %>% arrange(shared_id) %>%
    mutate(sample = gsub('lo', '', sample))
  polygon_tbl = left_join(polygon_tbl, phys_coord, by ='sample')
  unique_groups = as.character(unique(polygon_tbl$shared_id))
  unique_groups = unique_groups[order(nchar(unique_groups), unique_groups)]
  head(polygon_tbl)
  max_shared_samples = ifelse(is.na(param_max_shared), max(round((nchar(unique_groups)+1)/6)), param_max_shared)
  rainbow_cols = rainbow(max_shared_samples, start=0, end=1, alpha=0.6)
  for(i in c(1:length(unique_groups))) {
    this_group = unique_groups[i]
    num_shared_samples = round((nchar(this_group)+1)/6)
    #if(nchar(this_group)!=41) { next }
    if(num_shared_samples>max_shared_samples) { next }
    this_tbl = filter(polygon_tbl, shared_id == this_group)
    lwidth = sqrt(this_tbl$shared_num_var[1])*0.5
    polygon(this_tbl$x, this_tbl$y, border=rainbow_cols[num_shared_samples], lwd=lwidth)
    #polygon(this_tbl$x, this_tbl$y, border=alpha(col_gradient[num_samples-1], 0.8), lwd=lwidth)
  }
  
  lines(phys_coord$x, phys_coord$y, xaxt='n', yaxt='n', ann=FALSE, type='p', pch=20, col='black', xlim=c(0,140), ylim=rev(range(seq(70,0))), cex=0.5)
  if(label_coords) {
    text(phys_coord$x+1, phys_coord$y, labels=phys_coord$sample, adj=0, cex=0.5, col='black')
  }
  # Colorbar
  par(mar=c(5,3,5,10))
  image(y=2:max_shared_samples, z=t(2:max_shared_samples), col=rainbow_cols[2:max_shared_samples], axes=FALSE, cex.main=.8)
  axis(4,cex.axis=0.8,mgp=c(0,.5,0))
  
  dev.off()
  
  return(polygon_tbl)
}

plot_phys_coord_group <- function(phys_coord, sample_groups, outpath, lcol='red', lwidth=3) {
  # Plot physical coordinates of a specific set of dissections.
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 29-MAY-2018
  #
  # Args
  #   phys_coord: tibble storing coordinate data
  #   sample_group: chosen groups of samples for which to plot coordinates (eg. list(c('001', '002', '005'), c('010', '011')))
  #   outpath: path to store plot
  #   lcol: polygon line color (can be vector of same length as sample_groups)
  #   lwidth: polygon line width
  # Returns:
  #   list(
  #     polygon_tbl
  # )
  
  # Draw plot
  pdf(outpath, width=14/2.54, height=7/2.54)
  #layout(t(1:2),widths=c(1,1))
  par(mar=c(0,0,0,0))
  
  # First draw coordinates to get figure set up
  plot(phys_coord$x, phys_coord$y, xaxt='n', yaxt='n', ann=FALSE, type='p', pch=20, col='black', xlim=c(0,max(phys_coord$x)*1.1), ylim=rev(range(seq(max(phys_coord$y)*1.1,0))), cex=0.5, bty="n")
  
  # Second draw the polygon(s)
  for(i in c(1:length(sample_groups))) {
    this_group = sample_groups[[i]]
    # Match the set of chosen samples with the physical coordinates
    this_polygon = phys_coord %>%
      mutate(draw_polygon = ifelse(sample %in% this_group, TRUE, FALSE)) %>%
      filter(draw_polygon == TRUE)
    
    # Get middle point of chosen coordinates
    med_x = median(this_polygon$x)
    med_y = median(this_polygon$y)
    
    this_polygon = this_polygon %>%
      mutate(cart_x = x-med_x, cart_y = y-med_y)
    pol_tbl = cart2pol(this_polygon$cart_x, this_polygon$cart_y, degrees = TRUE)
    this_polygon$theta = pol_tbl$theta
    this_polygon$r = pol_tbl$r
    this_polygon = this_polygon %>% 
      mutate(is_outer = ifelse(r > median(r), TRUE, FALSE)) %>%
      arrange(theta)
    
    if(typeof(lcol)=='list') {
      this_lcol = lcol[[i]]
    } else {
      this_lcol = lcol
    }
    polygon(this_polygon$x, this_polygon$y, col=rgb(this_lcol[1], this_lcol[2], this_lcol[3], alpha=0.3), border=rgb(this_lcol[1], this_lcol[2], this_lcol[3], alpha=0.3), lwd=lwidth)
    lines(this_polygon$x, this_polygon$y, xaxt='n', yaxt='n', ann=FALSE, type='p', pch=20, col=rgb(this_lcol[1], this_lcol[2], this_lcol[3]), xlim=c(0,max(phys_coord$x)*1.1), ylim=rev(range(seq(max(phys_coord$y)*1.1,0))), cex=0.5, bty="n")
    
  }
  
  # Then the coordinates of all dissections
  #lines(phys_coord$x, phys_coord$y, xaxt='n', yaxt='n', ann=FALSE, type='p', pch=20, col='black', xlim=c(0,max(phys_coord$x)*1.1), ylim=rev(range(seq(max(phys_coord$y)*1.1,0))), cex=0.5, bty="n")
  text(phys_coord$x+1, phys_coord$y, labels=phys_coord$sample, adj=0, cex=0.5, col='black')
  
  dev.off()
}

get_deconstructSigs_colors <- function(sig9_col = '#00FF7F') {
  # Returns color scheme of COSMIC signatures as used in deconstructSigs
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 12-DEC-2017
  #
  # Args
  #   sig9_col  if not NA, then this color is used to replace signature 9 color with.
  # Returns:
  #   named vector
  all.sigs = c("Signature.1", "Signature.1A", "Signature.1B", "Signature.2",  "Signature.3", "Signature.4",  
               "Signature.5", "Signature.6", "Signature.7", "Signature.8", "Signature.9",  
               "Signature.10", "Signature.11", "Signature.12", "Signature.13", "Signature.14", 
               "Signature.15", "Signature.16", "Signature.17", "Signature.18", "Signature.19", 
               "Signature.20", "Signature.21", "Signature.R1", "Signature.R2", "Signature.R3",
               "Signature.U1", "Signature.U2","Signature.22", "Signature.23", "Signature.24",
               "Signature.25", "Signature.26", "Signature.27", "Signature.28", "Signature.29", "Signature.30",
               "unknown")
  color_vector = c("#023FA5", "#023FA5","#7D87B9","#BEC1D4","#D6BCC0","#BB7784",
                   "gold","#4A6FE3","#8595E1","#B5BBE3","#E6AFB9",
                   "#E07B91","#D33F6A","#11C638","#8DD593","#C6DEC7",
                   "#EAD3C6","#F0B98D","#EF9708","#0FCFC0","#9CDED6",
                   "#D5EAE7","#F3E1EB","#F6C4E1","#F79CD4","#866097",
                   "#008941","#A30059","#F6C4E1","#F79CD4","#866097",
                   "#008941","#A30059","#008080","#8B0000","#F4A460","#663399",
                   "#706563")
  
  # Change sig9 => makes it easier to see in liver project
  if(!is.na(sig9_col)) {
    sig9_idx = which(all.sigs == 'Signature.9')
    color_vector[sig9_idx] = sig9_col
  }
  
  names(color_vector)<-gsub('Signature.', '', all.sigs)
  color_vector
}

decode_hdp_sig_id <- function(hdp_exposures, hdp_decodes) {
  # Decodes the IDs of signatures as returned from HDP using a decode table.
  #
  # /// Author --- SIMON FELIX BRUNNER
  # /// Creation date --- 21-JUN-2018
  #
  # Args
  #   hdp_exposures
  #   hdp_decodes
  # Returns:
  #   hdp_exposures: with decoded Signature column.
  sig_encode = hdp_decodes %>%
    mutate(hdp_sig_id = sprintf('P%s', sig_num))
  
  hdp_exposures = hdp_exposures %>%
    left_join(sig_encode %>% dplyr::select(hdp_sig_id, sig_name), by=c('Signature'='hdp_sig_id')) %>%
    mutate(Signature = ifelse(is.na(sig_name), Signature, sig_name)) %>% dplyr::select(-sig_name)
  
  return(hdp_exposures)
}

plot_signature <- function(exposures, sig_name='') {
  p = ggplot(exposures, aes(x=context, y=counts, fill=mut_type)) + 
    geom_bar(stat='identity', color='black', size=0.2, width=0.7) + 
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(axis.text.x=element_text(size=6, family="Courier")) +
    facet_grid(.~mut_type, scales='free') +
    theme(panel.margin.x = unit(0, "lines"), legend.position='none') +
    scale_fill_manual(values = c('blue', 'black', "red", 'gray', "green", "pink")) +
    labs(x='Context', y='Counts', fill='Mutation', title=sig_name)
  return(p)
}
