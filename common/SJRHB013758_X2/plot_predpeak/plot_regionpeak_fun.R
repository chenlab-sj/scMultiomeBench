library(dplyr)
library(zoo)
library(ggplot2)
library(seqminer)
library(EnsDb.Hsapiens.v75)
library(ggbio)
library(patchwork)
library(stringr)
library(Signac)
library(tidyverse)
plot_region_pileups = function(chrom, # chromosome
                               start, # start coordinate
                               end, # end coordinate
                               bam_bed_file, # tabix indexed BED file with chrom,start,end,cell
                               grouping, # dataframe with cell,group,total_reads
                               ensdb=EnsDb.Hsapiens.v75, # EnsDb to allow display of genes
                               smoothing_window=100, # window to use for averaging
                               x_steps=3000,
                               gene_annot = FALSE,
                               col_palette,
                               fixed_ylim = FALSE) { # limits number of actual to plot data for
  
  # Validate the grouping dataframe
  if (! all(c('cell', 'group', 'total_reads') %in% colnames(grouping))) {
    stop('Grouping DF must have cell, grouping, and total_reads columns.')
  }
  
  # Get reads for the requested region of BED file (using tabix)
  range_string = paste0(chrom, ':', start, '-', end)
  reads = seqminer::tabix.read.table(tabixRange=range_string, tabixFile=bam_bed_file)
  colnames(reads) = c('chrom', 'start', 'stop', 'cell',"count")
  
  # Get count of cells per group for normalization
  grouping = grouping %>%
    #dplyr::filter(cell %in% reads$cell)%>%
    group_by(group) %>%
    mutate(cells_in_group = n()) %>%
    ungroup()
  
  # Restrict to cell barcodes specified by user
  #  reads = reads[, -5]
  reads = dplyr::filter(reads, cell %in% grouping$cell)
  reads = dplyr::inner_join(reads, grouping, by='cell')
  

  ##################################################
  ## convert  count to coverage for each cell
  ## count.norm (by frag for bc) = count/ total_reads(for_bc) * 1e4
  
  #reads$count.norm <- reads$count/reads$total_reads * 1e4
  reads$count.norm <- reads$count 
  #reads$count.norm <- log2(reads$count.norm + 1 )
  

  
  # Expand reads so have one line per base
  expanded_ranges = dplyr::bind_rows(lapply(1:nrow(reads), function(i) {
    interval = reads[i, 'start']:reads[i, 'stop']
    
    return(data.frame(position=interval,
                      value=rep(reads[i, 'count.norm'],length(interval)),
                      cell=rep(reads[i, 'cell'],length(interval)),
                      group=rep(reads[i, 'group'], length(interval)),
                      #scaling_factor=rep(reads[i, 'total_reads'],length(interval) ),
                      cells_in_group=rep(reads[i, 'cells_in_group'], length(interval))))
  }))
  

  expanded_ranges %>%
#    left_join(frag_groupcount, by = "group")%>%
    group_by(position, group,cells_in_group) %>%
    summarize(total=sum(value)) %>%
    mutate(total_norm = total/mean (cells_in_group))%>%
#    mutate(total_norm = total_norm/cells_in_group)%>%
    group_by(group) %>%
    arrange(position)->peak_value
  


  # Finally, restrict to x_steps points for plotting to speed things up
  ## for plot
  
  total_range = end - start
  stepsize = ceiling(total_range / x_steps)
  positions_to_show = seq(start, end, by=stepsize)
  peak_value.plot <- peak_value %>%
  group_by(group)%>%
    mutate(total_smoothed=zoo::rollapply(total_norm, smoothing_window, sum, align='center',fill=NA)) %>%
     ungroup()%>%
    subset(position %in% positions_to_show)
  # 
  peak_value.vector <- peak_value.plot %>% 
     dplyr::select(position, group, total_smoothed)%>% 
     pivot_wider(names_from = group, values_from = total_smoothed)
    

  
  
  
  
  # # Make plots faceted by group
 if (isFALSE(fixed_ylim)){
   max_value = max(peak_value.plot$total_smoothed, na.rm = TRUE)
   max_ylim = max_value * 1.4
   }else{
     max_ylim = fixed_ylim
     }
  
  plot_object = ggplot() +
    geom_bar(data=peak_value.plot, aes(position, total_smoothed, color = group, fill = group), stat='identity') +
    facet_wrap(~group, ncol=1, strip.position="right") +
    ylim(0, max_ylim) +
    xlab('Position (bp)') +
    ylab(paste0('Log normalized Signal (', smoothing_window, ' bp window)')) +
    theme_classic()+
    scale_color_manual(values = col_palette)+
    scale_fill_manual(values = col_palette)+
    theme(legend.position="none")
  
  
  # Add gene track using ensembl ID
  if (gene_annot==TRUE){
    gr <- GRanges(seqnames = str_replace(chrom, 'chr', ''), IRanges(start, end), strand = "*")
    filters = AnnotationFilterList(GRangesFilter(gr), GenebiotypeFilter('protein_coding'))
    # 
    genes = autoplot(ensdb, filters, names.expr = "gene_name") + geom_text(size = 6)+
      theme(plot.margin = margin(t = 50, r = 10, b = 10, l = 10, unit = "pt"))
    # 
    genes_plot = genes@ggplot +
      xlim(start, end) +
      theme_classic() 
    # Now combine the two plots into one (using patchwork)
    return( list(peak_value.vector,genes_plot / plot_object + xlim(start, end) + plot_layout(heights = c(1, 9)) ))
  }else{
    return(list(peak_value.vector,plot_object))
  }
  
}

