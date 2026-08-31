# NOTE: paths below are placeholders. See config/config.R and the README.
library(dplyr)

## Fig7 revision: keep the 9 published methods' ORIGINAL results (read unchanged from the cluster
## sum_metrics_clean.csv / subset_stability.csv / reproducbility.csv / runtime_memory_sel.csv), and
## SPLICE IN the 3 new methods (MaxFuse/MIDAS/scButterfly). New-method grouped scores + BRCA rare-cell
## stability + reproducibility are pre-assembled in new_methods_grouped.csv; new-method runtime/memory
## (from the LSF run logs) are appended in runtime_memory_sel_plus.csv. Outputs -> this dir.
## paths auto-detect cluster vs Mac mount (so this runs both locally and as a bsub job)
ON_CLUSTER <- dir.exists("/path/to/multiomeBench")
SUM <- if (ON_CLUSTER) "/path/to/multiomeBench/sum_plot" else
                       "/path/to/multiomeBench/sum_plot"
LSA <- if (ON_CLUSTER) "/path/to/multiomeBench/common" else "/path/to/multiomeBench/common"

top10_methods <- c('scglue(multiome)','scVI','scglue','scJoint','Seurat(CCA)','Portal','simba','BindSC','scBridge',
                   'MaxFuse','MIDAS','scButterfly')
pbmc3k.file  <- file.path(LSA, 'pbmc3k/benchmark_matrix/sum_metrics_clean.csv')
pbmc3k <- read.csv(pbmc3k.file, row.names = 1)

RMS.file  <- file.path(LSA, 'Mast607A/benchmark_matrix/sum_metrics_clean.csv')
RMS <- read.csv(RMS.file,row.names = 1)
RMS <- RMS[which(RMS$method %in% top10_methods),]

BRCA_1samp.file  <- file.path(LSA, 'HT243B1-S1H4/benchmark_matrix/sum_metrics_clean.csv')
BRCA_1samp <- read.csv(BRCA_1samp.file ,row.names = 1)
BRCA_1samp <- BRCA_1samp[which(BRCA_1samp$method %in% top10_methods),]
BRCA_macrosub.file  <- file.path(LSA, 'HT243B1-S1H4/subset_stability.csv')
BRCA_macrosub <- read.csv(BRCA_macrosub.file  ,row.names = 2)
BRCA_1samp$Rare.cell.type.performance<-BRCA_macrosub[BRCA_1samp$method,'stability']

BRCA_rep.file  <- file.path(LSA, 'HT243B1-S1H4/reproducbility.csv')
BRCA_rep <- read.csv(BRCA_rep.file  ,row.names = 2)
BRCA_1samp$Computational.reproducibility<-BRCA_rep[BRCA_1samp$method,'reproducibility']


BMMC.file  <- file.path(LSA, 'BMMC_d1/sum_metrics_clean.csv')
BMMC <- read.csv(BMMC.file,row.names = 1)
BMMC <- BMMC[which(BMMC$method %in% top10_methods),]

# BRCA_2samp.file  <- '/path/to/multiomeBench/common/HT263_S1H1_HT243_S1H4/benchmark_matrix/benchmark_sum_clean.csv'
# BRCA_2samp <- read.csv(BRCA_2samp.file )

comp_res.file <- file.path(SUM, 'runtime_memory_sel_plus.csv')   ## cluster rows + 12 new-method rows
comp_res <- read.csv(comp_res.file)
comp_res_summary <-comp_res %>%
  group_by(methods) %>%
  summarise(
    avg_running_time = mean(running.time_raw, na.rm = TRUE),
    avg_memory_usage = mean(memory.usage_raw, na.rm = TRUE)
  )

comp_res_summary$avg_running_time<- ifelse(comp_res_summary$avg_running_time > 200000, 200000,comp_res_summary$avg_running_time) 

comp_res <-comp_res[,c('methods','GPU.requirment', 'Data.requirement')]

comp_res <- merge(comp_res, comp_res_summary, by ='methods')

comp_res$Data <- "Usability"


library(dplyr)
library(tidyr)
library(purrr)
library(RColorBrewer)


min_max_normalize <- function(x) {
  if (all(is.na(x))) {
    return(x)  # If all values are NA, return as is
  }
  return(ifelse(is.na(x), NA, (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))))
}

rev_min_max_normalize <- function(x) {
  if (all(is.na(x))) {
    return(x)  # If all values are NA, return as is
  }
  return(ifelse(is.na(x), NA, 1 - ((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))))
}

comp_res$`Running time`<-rev_min_max_normalize(comp_res$avg_running_time)
comp_res$`Memory usage`<-rev_min_max_normalize(comp_res$avg_memory_usage)

comp_res<- comp_res %>%
  mutate(`Overall usability` = (`GPU.requirment`+ `Running time` + `Memory usage`)/6 + (`Data.requirement`/2))

colnames(comp_res)<-c('method','GPU requirment', 'Data requirement',"avg_running_time","avg_memory_usage", 
                      "Data","Running time","Memory usage","Overall usability")

res_list <- list(pbmc3k,RMS, BRCA_1samp,BMMC,comp_res)
combined_res <- bind_rows(res_list)

## ---- SPLICE IN the 3 new methods (grouped scores + BRCA rare-cell stability + reproducibility) ----
## Data labels (pbmc3k / Mast607A / BRCA / BMMC_d1) match the cluster clean files, so the recode below
## treats them identically to the published rows. Runtime/memory (Usability) already come in via comp_res.
newg <- read.csv(file.path(SUM, "new_methods_grouped.csv"), check.names = TRUE)   # dotted col names
combined_res <- bind_rows(combined_res, newg)

# comp_res.file <- '/path/to/multiomeBench/common/runtime_memory.csv'
# comp_res <- read.csv(comp_res.file) 
# comp_res <- comp_res[comp_res$Dataset %in% c( "pbmc3k","BRCA","BMMC","SJRHB010958_X2" ),c('Dataset','methods', 'running.time', 'memory.usage')] ## use pbmc10k to represent
# colnames(comp_res)<-c('data',"Methods", "Running time","Memory usage")
# comp_res$data



sum_format <- function(data,dataname){
  data %>%
    pivot_longer(cols = -c(method,Data) ,
                 names_to = "metrics",
                 values_to = "value")->data.format

    #mutate(data = dataname)->data.format
  
  return(data.format)
  
}






combined_res$`Omics gap reduction` = min_max_normalize(combined_res$Omics.gap.reduction)
combined_res$`Biological diversity \n preservation` = min_max_normalize(combined_res$Bio.conservation)
combined_res$`Overall integration performance` = min_max_normalize(combined_res$Overall.Performance)
combined_res$`Optimal alignment \n between modalities` = min_max_normalize(combined_res$Alignment.accuracy)
combined_res$`Rare cell type performance` = min_max_normalize(combined_res$Rare.cell.type.performance)
combined_res$`Computational reproducibility` = min_max_normalize(combined_res$Computational.reproducibility)

combined_res$`Batch effects correction` = min_max_normalize(combined_res$Sample.batch.correction)







combined_res %>%
  select('GPU requirment', 'Data requirement', 'Running time','Memory usage','Overall usability' ,"Omics gap reduction","Biological diversity \n preservation",
         "Overall integration performance","Optimal alignment \n between modalities" ,'Computational reproducibility',"Rare cell type performance","Batch effects correction",
         "method","Data")%>%
  pivot_longer(cols =c('GPU requirment', 'Data requirement', 'Running time','Memory usage','Overall usability' ,"Omics gap reduction","Biological diversity \n preservation",
                       "Overall integration performance","Optimal alignment \n between modalities" ,'Computational reproducibility',"Rare cell type performance","Batch effects correction"), 
                                      names_to = "metrics",
                                      values_to = "value")%>%
  mutate(Data = case_when(
    Data == "pbmc3k"~"PBMC",
    Data== "BRCA"~"Cancer data",
    Data =="Mast607A"  ~ "Cell types with \n subtle differences",
    Data =="BMMC_d1"~"Additional sample \n batches",
    
    TRUE ~ Data
  ))->combined_res



############################################################
library(ggplot2)

## combined_res_order
# combined_res2 <- combined_res
#combined_res2$value <- ifelse(combined_res$metrics == 'Overall usability', (combined_res$value)* 0.2, combined_res$value)

combined_res%>%
  distinct()%>%
  drop_na() %>%
  filter( method %in% top10_methods)%>%
  filter(metrics %in% c('Overall integration performance', "Overall usability"))%>%
  #filter(metrics %in% c('Overall integration performance', "Rare cell type performance","Batch effects correction"))%>%
  
  group_by(method)%>%
  summarize(rank_score = sum(value, na.rm = TRUE))%>%
  arrange(desc(rank_score))->method_rank


combined_res$method <- factor(combined_res$method,levels=rev(c(method_rank$method) )) 
combined_res$Data <- factor(combined_res$Data, level = c("PBMC","Cancer data", "Cell types with \n subtle differences","Additional sample \n batches","Usability" ))  
combined_res$metrics <- factor(combined_res$metrics, levels = c("Biological diversity \n preservation","Omics gap reduction",
                                                                "Optimal alignment \n between modalities","Batch effects correction",'Computational reproducibility',"Rare cell type performance",
                                                                "GPU requirment", "Running time","Memory usage", "Data requirement",'Overall integration performance', "Overall usability"))  
combined_res %>% 
  drop_na()%>%
  filter(metrics %in% c("Biological diversity \n preservation","Omics gap reduction",
                        "Optimal alignment \n between modalities","Batch effects correction",'Computational reproducibility',"Rare cell type performance",
                        "GPU requirment", "Running time","Memory usage", "Data requirement"))%>%
  ggplot(aes(x = metrics, y = method) )+
  geom_point( aes(fill = value), alpha = 0.75, shape = 21, size = 8)+
  facet_grid(~Data,scales = "free_x",space = "free_x")+
  theme_classic()+
  theme(panel.spacing = unit(0.5, "lines"))+
  theme(text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 0.85, hjust=0.85))+
    scale_fill_gradientn(
      colors = brewer.pal(9, "YlOrRd") , 
      breaks = c(1,0), 
      labels = c("good", "poor"), 
      name = "Relative performance",
    )+
  xlab('')+
  ylab("")+
  theme(panel.spacing = unit(1, "lines"))+
  theme(
        axis.text.x = element_text(angle = 45, vjust = 0.85, hjust=0.85))+
  theme(legend.position = "bottom")+
  theme(text=element_text(size=14, color = "black"),
        plot.title = element_text(size = 14, color = "black", face = "bold"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        strip.text = element_text(
          size = 14),
        strip.background = element_rect(fill = NA, colour = NA))->p2

pdf(file.path(SUM,'sum_perform.pdf'), width=12.5, height=6.5)
plot(p2)
dev.off()

# tab20_colors <- c(
#   "#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78",
#            "#2ca02c", "#98df8a", "#d62728", "#ff9896",
#            "#9467bd", "#c5b0d5", "#8c564b", "#c49c94",
#            "#e377c2", "#f7b6d2", "#7f7f7f", "#c7c7c7",
#            "#bcbd22", "#dbdb8d", "#17becf", "#9edae5"
# )

method_colors <- c(
  'scglue' = "#1f77b4", 
  'scglue(multiome)' = "#aec7e8", 
  #'BindSC' = "#ff7f0e", 
  'Portal' = "#ffbb78", 
  'scBridge' = "#2ca02c", 
  'simba' = "#98df8a", 
  'Seurat(CCA)' = "#d62728", 
  'scJoint' = "#ff9896", 
  'scVI' = "#9467bd",
  'MinNet' = "#c5b0d5",
  'BindSC'="#e377c2",
  'MaxFuse' = "#8c564b",
  'MIDAS' = "#17becf",
  'scButterfly' = "#bcbd22"
)


combined_res %>%
  filter( method %in% top10_methods)%>%
  filter(metrics %in% c( "Overall integration performance"))%>%
  group_by(Data) %>%
  mutate(rank = rank(-value, ties.method = "first")) %>%
  ungroup()%>%
  as.data.frame()->combined_rank

# combined_res %>%
#   filter( method %in% top10_methods)%>%
#   filter(metrics == "	Overall usability")%>%
#   group_by(Data) %>%
#   mutate(rank = rank(-value, ties.method = "first")) %>%
#   ungroup()%>%
#   as.data.frame()

combined_rank%>%
  drop_na()%>%
  group_by(method)%>%
  summarise(value = mean(value, na.rm = TRUE))%>%
  mutate(ave_rank = rank(-value, ties.method = "first")) %>%
  ungroup()%>%
  as.data.frame()%>%
  arrange(ave_rank)->combined_rank_ave


combined_rank_overall <- data.frame(method = combined_rank_ave$method,
                                    Data = rep("Average rank", length(top10_methods) ),
                                    rank = combined_rank_ave$ave_rank,
                                    metrics = "Overall integration performance",
                                    value = combined_rank_ave$value)

combined_rank_all  <- rbind(combined_rank_overall,combined_rank)
combined_rank_all <- combined_rank_all %>%
  filter(Data %in% c("Average rank","PBMC","Cancer data","Cell types with \n subtle differences","Additional sample \n batches" ))
combined_rank_all$Data <- factor(combined_rank_all$Data, level = c("Average rank","PBMC","Cancer data","Cell types with \n subtle differences","Additional sample \n batches" ))  

## --- self-contained replacement for ggbump::geom_bump (no package install needed) ---
## Reproduces geom_bump(smooth = 8): between each pair of consecutive rank points, interpolate a
## logistic S-curve (same formula/defaults as ggbump), so the connecting lines are visually identical.
bump_sigmoid <- function(x, y, smooth = 8, n = 100) {
  do.call(rbind, lapply(seq_len(length(x) - 1L), function(i) {
    t <- seq(-smooth, smooth, length.out = n)
    s <- exp(t) / (exp(t) + 1)
    data.frame(x = (t + smooth) / (smooth * 2) * (x[i + 1L] - x[i]) + x[i],
               y = s * (y[i + 1L] - y[i]) + y[i])
  }))
}

## discrete Data axis -> numeric positions 1..5 (bump curves need a continuous x to interpolate over)
rank_x_levels <- c("Average rank", "PBMC", "Cancer data",
                   "Cell types with \n subtle differences", "Additional sample \n batches")
combined_rank_all$x_num <- as.integer(factor(combined_rank_all$Data, levels = rank_x_levels))

## one interpolated path per method (drawn with geom_path below, in place of geom_bump)
bump_paths <- do.call(rbind, lapply(split(combined_rank_all, combined_rank_all$method), function(d) {
  d <- d[order(d$x_num), ]
  if (nrow(d) < 2L) return(NULL)          # need >=2 points to interpolate; rbind() drops NULLs
  p <- bump_sigmoid(d$x_num, d$rank)
  p$method <- d$method[1]
  p
}))

ggplot() +
  geom_path(data = bump_paths, aes(x = x, y = y, group = method, color = method), size = 1) +
  geom_point(data = combined_rank_all, aes(x = x_num, y = rank, color = method), size = 5) +
  geom_text(data = combined_rank_all, aes(x = x_num, y = rank, label = rank),
            vjust = 0.6, size = 3.5, color = 'white') +
  scale_color_manual(values = method_colors) +
  scale_y_reverse(breaks = 1:nrow(combined_rank)) +
  labs(x = "", y = "") +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                legend.title = element_blank(), legend.position = 'bottom',    panel.border = element_blank()
  )  + scale_y_reverse(
    breaks = 1:n_distinct(combined_rank$method),
    labels = combined_rank_all %>% filter(Data == "Average rank") %>% arrange(rank) %>% pull(method)) +  
  # x expansion is set once on scale_x_continuous below
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    # axis.text.x = element_text(vjust = -1, hjust = 0.5), # Move x-axis text to the top
    axis.text.y = element_text(size = 10),
    legend.position = "bottom", # Remove the legend
    panel.background = element_blank(), # Remove background
    plot.background = element_blank()
  )+
  theme(text=element_text(size=14, color = "black"),
        plot.title = element_text(size = 14, color = "black", face = "bold"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        strip.text = element_text(
          size = 14),
        strip.background = element_rect(fill = NA, colour = NA))+
  scale_x_continuous(breaks = seq_along(rank_x_levels), labels = rank_x_levels,
                     position = "top", expand = c(0, 0.6)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) ->p3   # method color legend: 2 rows (was 3)


  pdf(file.path(SUM,'sum_perform_rank.pdf'), width=12, height=4.5)
  plot(p3)
  dev.off()
  
  
  ##################################################
  ## plot scability

comp_res.file <- file.path(SUM, 'runtime_memory_sel_plus.csv')
comp_res <- read.csv(comp_res.file)
comp_res <- comp_res %>%
  filter(methods %in% top10_methods) %>%
  filter(Dataset != "pbmc10k")   # dropped from the scalability figure (see legend)
library(reshape2)
library(ggplot2)
library(dplyr)
comp_res_sel <- melt(comp_res , id.vars = c("Dataset", "methods", "Test.cell.counts"), 
                measure.vars = c("running.time", "memory.usage"),
                variable.name = "Metric", value.name = "Value")

comp_res_sel$Metric <- recode(comp_res_sel$Metric, 
                              "running.time" = "Running time", 
                              "memory.usage" = "Memory usage")
ggplot(comp_res_sel, aes(x = Test.cell.counts, y = Value, color = methods)) +
  geom_point(size = 3) +  # Scatter plot
  geom_line(aes(group = methods), linetype = "dashed") +  # Connect points with dashed lines
  facet_wrap(~Metric, scales = "free") +  # Separate plots for running time & memory usage
  scale_y_log10() +  # Log scale for better visualization
  scale_x_log10() +
  theme_minimal() +
  labs(#title = "Performance Comparison by Test Cell Number",
       x = "Cell numbers",
       y = "Values") +
  theme(legend.position = "bottom", text = element_text(size = 12))+
  scale_color_manual(values = method_colors)+
  theme(text=element_text(size=14, color = "black"),
        plot.title = element_text(size = 14, color = "black", face = "bold"),
        axis.text.x = element_text(size = 14, color = "black"),
        axis.text.y = element_text(size = 14, color = "black"),
        strip.text = element_text(
          size = 14),
        strip.background = element_rect(fill = NA, colour = NA))->p3

pdf(file.path(SUM,'scalability.pdf'), width=9, height=4.5)
plot(p3)
dev.off()



