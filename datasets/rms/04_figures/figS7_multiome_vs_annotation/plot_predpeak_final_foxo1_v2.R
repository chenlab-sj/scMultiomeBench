## NOTE: paths below are placeholders. See config/config.R and the README
## for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
## Fig S7: Multiome-vs-Annotation pileup strip (one locus per script; the
## three strips were assembled manually into the published figure).
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(tidyr)
set.seed(123) 
mypalette <-c("#1F77B4", "#FF7F0E", "#2CA02C")

data.dir <- "/path/to/data/RMS/SJRHB013758_X2_scATAC/"
#data.dir <- "/path/to/data/RMS/SJRHB013758_X2_scATAC/"
fpath = paste0(data.dir,"GSM5390512_1951549_DYE3118_fragments.tsv.gz")
h5.file <- paste0(data.dir,"SJRHB013758_X2_scATAC_filtered_feature_bc_matrix.h5")

KNN.dir = "/path/to/multiomeBench/SJRHB013758_X2/knn_test/"
#KNN.dir = "/path/to/multiomeBench/SJRHB013758_X2/knn_test/"
KNN.label.file = paste0(KNN.dir,"knn_k10_pred_label.csv")
KNN.label <- read.csv(KNN.label.file, row.names = 1)
rownames(KNN.label) <- gsub("_atac", "", rownames(KNN.label))
bc <- rownames(KNN.label)


## add control and random
atac.annot.file = paste0(data.dir,'SJRHB013758_X2_clusters.csv') 
atac.annot <- read.csv(atac.annot.file, row.names =1)
colnames(atac.annot)[colnames(atac.annot) == "predicted.id"] <- "annotation"
atac.annot$random <-  sample(unique(atac.annot$annotation), nrow(atac.annot), replace = TRUE)
atac.annot <- atac.annot[bc,]
atac.meta <- merge(atac.annot, KNN.label,by = "row.names", all = TRUE )
row.names(atac.meta)<-atac.meta$Row.names
atac.meta$Row.names <- NULL


#3 random accuracy
type <- unique(atac.meta$annotation)
accuracy_list <- lapply(type, function(category) {
  actual <- atac.meta$annotation == category
  predicted <- atac.meta$random == category
  accuracy <- sum(actual & predicted) / sum(actual)
  count <- sum(predicted)
  data.frame(type= category, accuracy = accuracy, n = count)
})
accu_random <- do.call(rbind, accuracy_list)
accu_random$pipeline <- "Random"
###########################################
## create seurat obj
# library(EnsDb.Hsapiens.v75)
# annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
# seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
# genome(annotation) <- "hg19"
# 
# 
# counts <- Read10X_h5(filename = h5.file)
# counts <-counts[,rownames(atac.meta)]
# chrom_assay <- CreateChromatinAssay(
#   counts = counts,
#   sep = c(":", "-"),
#   genome = 'hg19',
#   fragments = fpath,
#   annotation = annotation
# )
# 
# 
# atac.label <- CreateSeuratObject(
#   counts = chrom_assay,
#   assay = "peaks",
#   meta.data = atac.meta
# )

#######################################
## get coverage from different labels

############################
# plot_list<-list()
# peak_list <- list()
# for (i in 1:ncol(atac.meta)){
#   pipeline <- colnames(atac.meta)[i]
#   print(pipeline)
#   cov_plot <- CoveragePlot(
#     object = atac.label,
#     region = "chr11-17725000-17745000",
#     annotation = FALSE,
#     peaks = FALSE,
#     group.by=c(pipeline),
#     ymax = 90, 
#     window = 10
#   )
#   peak<-cov_plot$data
#   peak$method <- pipeline
#   peak_list[[i]]<-peak
#   plot_list[[i]]<-cov_plot + ggtitle(pipeline)
#   
# }
# ## a quick check
# library(gridExtra)
# png(paste0(KNN.dir, "prdicted_ataclabel_myod1.png"), width=3000, height=1000, res= 100)
# grid.arrange(grobs = plot_list, ncol = 5)
# dev.off()
# 
# 
# peak_df <- do.call(rbind, peak_list)
# 
# peak_df %>%
#   dplyr::select(group, position, coverage, method)%>%
#   pivot_wider(names_from = group, values_from = coverage)->peak_out
# 
# write.csv(peak_out, file = paste0(KNN.dir, "preducted_ataclabel_cov.csv"), row.names = FALSE, quote = FALSE)

## module load macs2/2.1.1
## 

# pred_peaks<- list()
# for (i in 1:ncol(atac.meta)){
#   method = colnames(atac.meta)[i]
#   print(method)
#   peaks <- CallPeaks(
#     object = atac.label,
#     group.by = method,
#     macs2path="/path/to/tools/macs2/bin/macs2"
#   )
#   
#   pred_peaks[[method]] <- peaks
# }
# 
# save(pred_peaks,file = paste0(KNN.dir, "knn_k10_pred_macs2peak.RData"))
# 
# 
# typegr_sum <-function(grctrl, grpred, celltype){
#   grpred_sub <- subset(grpred, grepl(celltype, peak_called_in))
#   if(length(grpred_sub)!=0){
#     grctrl_sub <- subset(grctrl, grepl(celltype, peak_called_in))
#     overlap <- findOverlaps(grctrl_sub, grpred_sub, maxgap = 5000)
#     precision <- length(unique(subjectHits(overlap))) / length(grpred_sub)
#     recall <- length(unique(queryHits(overlap)))/ length(grctrl_sub)
#     f1 <- 2*(precision * recall) / (precision + recall)
#     return(c(celltype, precision, recall, f1))
#   }else{
#     return(c(celltype,0,0,0))
#   }
# 
# }
# grctrl <- pred_peaks[["annotation"]]
# result_dfs <- list()
# for(i in 2: ncol(atac.meta)){
#   method = colnames(atac.meta)[i]
#   print(method)
#   grpred <- pred_peaks[[method]]
#   result <- lapply(type, function(celltype) {
#     typegr_sum(grctrl, grpred, celltype)
#   })
#   result_df <- do.call(rbind, lapply(result, function(x) {
#     data.frame(t(x))
#   }))
#   result_df$method = method
#   result_dfs[[i -1]] <- result_df
#   
# }
# peakcall_sum <-  do.call(rbind, result_dfs)
# colnames(peakcall_sum)<-c('celltype', 'precision', 'recall', 'f1','method')
# 
# write.csv(peakcall_sum, paste0(KNN.dir, "knn_k10_peakcall_sum.csv"))
# 
# peakcall_sum %>%
#   dplyr::select(method, celltype, f1)%>%
#   pivot_wider(names_from = celltype,
#               values_from = f1)->peakcall_f1
# 
# write.csv(peakcall_f1, paste0(KNN.dir, "knn_k10_peakcall_f1.csv"))
# ## use customize function to generate peak_out
# source('/path/to/multiomeBench/SJRHB013758_X2/plot_predpeak/ataclabel2peak.R')

############################################
## finalized peak plot and statistics
###############################################
method.level <-c("Annotation","simba","Seurat.CCA.","scDART" ,"scglue" ,
                 "scBridge" ,"scglue.multiome.", "scJoint","BindSC" ,"scVI" ,
                 "Cobolt",  "Conos" ,"Portal" , "LIGER"  ,"Unioncom","scMoMaT"   ,
                 "MultiMAP" ,"Seurat.WNN.","Random"  )
method.level.sel <- c("Multiome","Annotation","Seurat \n (CCA)",'scDART','scglue \n (multiome)','scBridge','scVI',"Random")
method.level.sel <- c("Multiome","Annotation")

peak_out <- read.csv(paste0(KNN.dir, "predicted_ataclabel_foxo1value.csv"))
peak_out_train <- read.csv("/path/to/data/RMS/lca_celltype/Mast607A/predicted_ataclabel_foxo1value.csv")
peak_out <- rbind(peak_out_train, peak_out)
peak_out<-pivot_longer(peak_out, 
                       cols = c('Mesoderm', 'Myoblast','Myocyte'), 
                       names_to ="type")%>%
  mutate(pipeline =case_when(pipeline == "Seurat.CCA."~ "Seurat \n (CCA)",
                             pipeline == "Seurat.WNN."~ "Seurat(WNN)",
                             pipeline == "scglue.multiome."~ "scglue \n (multiome)",
                             pipeline == "annotation"~ "Annotation",
                             pipeline == "random"~ "Random",
                             pipeline == "lca_label"~ "Multiome",
                             
                             TRUE ~ pipeline))->peak_out
peak_out<- peak_out%>% filter( pipeline%in% method.level.sel)%>%
  mutate(pipeline= factor(pipeline, levels=c(method.level.sel)))%>%
  drop_na()




KNN.res <- read.csv(paste0(KNN.dir,"knn_k10_pred_label.csv"), row.name =1)
count_list <- list()
for(col in colnames(KNN.res)) {
  category_counts <- table(KNN.res[[col]])
  count_list[[col]] <- as.data.frame(category_counts)
}
count_df <- do.call(rbind, lapply(names(count_list), function(name) {
  df <- count_list[[name]]
  df$column <- name
  df
}))
colnames(count_df)<-c('type','n','pipeline')
count_df <- count_df%>%
  mutate(pipeline =case_when(pipeline == "Seurat.CCA."~ "Seurat \n (CCA)",
                             pipeline == "Seurat.WNN."~ "Seurat(WNN)",
                             pipeline == "scglue.multiome."~ "scglue \n (multiome)",
                             pipeline == "annotation"~ "Annotation",
                             pipeline == "random"~ "Random",
                             TRUE ~ pipeline))%>%
  dplyr::filter(pipeline %in% method.level.sel)


KNN.accu <- read.csv(paste0(KNN.dir,"knn_k10_pred_accu.csv"), row.names = 1)
KNN.accu <- KNN.accu %>%
  mutate(pipeline =case_when(pipeline == "Seurat(CCA)"~ "Seurat \n (CCA)",
                             pipeline == "scglue(multiome)"~ "scglue \n (multiome)",
                             TRUE ~ pipeline))%>%
  merge(count_df, by.x = c("pipeline", "type"), by.y = c("pipeline", "type"))%>%
  dplyr::filter(pipeline %in% method.level.sel)



#######################################################################
atac.annot %>%
  group_by(annotation)%>%
  summarise(n = n())%>% 
  mutate(pipeline = 'Annotation')%>%
  mutate(type = annotation)%>%
  mutate(accuracy = 1)%>%
  select('pipeline','type','accuracy','n')-> accu_annot

KNN.accu<-KNN.accu%>%
  rbind(accu_random)%>%
  rbind(accu_annot)

KNN.accu<-KNN.accu%>%
  mutate(label = paste0("accuracy:",round(accuracy,2),",n=",n))%>%
  mutate(pipeline = factor(pipeline, levels=c(method.level.sel)))%>%
  drop_na()


## add train cell type count
out.dir <- "/path/to/data/RMS/lca_celltype/Mast607A/"
lca_label <- read.csv(paste0(out.dir,"lca_label.csv"), row.names = 1)
train_annot <- as.data.frame(table(lca_label$label))
colnames(train_annot)<-c('type','n')
train_annot$pipeline <- 'Multiome'
train_annot$accuracy <- NA
train_annot$label <- paste0('n=',as.character(train_annot$n))
KNN.accu <- rbind(train_annot, KNN.accu)
KNN.accu$pipeline = factor(KNN.accu$pipeline, levels=c(method.level.sel))
# plot_list<-list()
# for (i in 1:length(method.level.sel)){
#   print(method.level.sel[i])
#   df_label <- KNN.accu %>%
#     dplyr::filter(method ==method.level.sel[i])%>%
#     dplyr::filter(type != "overall")
#   print(df_label)
#   peak_out %>% 
#     dplyr::filter(method ==method.level.sel[i])%>%
#     ggplot(aes(x =position, y = value, fill =type,color = type ))+
#     geom_bar( stat='identity')+
#     facet_wrap(~type, ncol =1,strip.position="right")+
#     geom_text(
#       data= df_label,
#       size = 3,
#       mapping = aes(x =41150000, y = 6, label = label),
# 
#     )+
#     ylim(0, 8) +
#     xlab('chr11 position (bp)') +
#     ylab(method.level.sel[i]) +
#     theme_classic()+
#     scale_color_manual(values = mypalette)+
#     scale_fill_manual(values = mypalette)+
#     theme(legend.position="none")+
#     #ggtitle(method.level.sel[i])+
#     theme(plot.title = element_text(hjust = 0.5, size = 10))+
#     theme(
#       strip.background = element_blank(),
#       strip.text.y = element_blank()
#     )->p
#   
#   plot_list[[i]]<-p
# }


# 
# library(gridExtra)
# pdf(paste0(KNN.dir, "predicted_ataclabel_myod1_sel.pdf"), width=20, height=1.8)
# grid.arrange(grobs = plot_list, ncol = 8)
# dev.off()
# 
# 
# 
# 
# df_label <- KNN.accu %>%
#   dplyr::filter(method %in% method.level.sel)%>%
#   dplyr::filter(type != "overall")

## CORRECTION (2026-08-31): the original script doubled `value` for the
## Annotation/Myocyte track within chr13:41,240,780-41,240,810 at this point.
## That transformation has been removed; this script plots the pileup values
## unmodified, so a re-render may differ from the published panel in that
## 30-bp window.
peak_out %>%
  
  ggplot(aes(x =position, y = value, fill =type,color = type ))+
  geom_bar( stat='identity')+
  
  #facet_grid(type~method)+
  facet_grid(pipeline~type)+
  geom_text(
    data= KNN.accu,
    size = 3.5,
    mapping = aes(x =41200000, y = 7.5, label = label),
    
  )+
  ylim(0, 8) +
  xlab('chr13 position (bp)') +
  theme_classic()+
  scale_color_manual(values = mypalette)+
  scale_fill_manual(values = mypalette)+
  theme(legend.position="none")+
  #ggtitle(method.level.sel[i])+
  theme(plot.title = element_text(hjust = 0.5, size = 11))+
  theme(panel.spacing.x = unit(0.05, "lines"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
  )+
  theme(text=element_text(size=11, color = "black"),
        plot.title = element_text(size = 11, color = "black", face = "bold"),
        axis.text.x = element_text(size = 11, color = "black"),
        
        strip.text = element_text(
          size = 10)
        #strip.text.y = element_blank()
  ) +
  theme(
    strip.background = element_blank(),   # Remove the facet box background
    panel.background = element_rect(fill = "transparent", color = NA),  # Make panel background transparent
    plot.background = element_rect(fill = "transparent", color = NA)    # Make plot background transparent
  )+
  #annotate("rect", xmin = 41129804, xmax = 41240734 , ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "blue")+
#   start = c(41145735, 41224650),  # List of SE region starts
# end = c(41171933, 41256947),    # List of SE region ends
  annotate("rect", xmin = 41145735, xmax = 41171933 , ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "blue")+
  annotate("rect", xmin = 41224650, xmax = 41256947 , ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "blue")+
  
  scale_x_continuous(breaks=seq(41100000,41300000,150000))->p


png(paste0(KNN.dir, "predicted_ataclabel_foxo1_sel.png"), width=2000, height=1600, res = 400)
pdf(paste0(KNN.dir, "predicted_ataclabel_foxo1_sel_v2.pdf"), width=7, height=6)
plot(p)
dev.off()

pdf(paste0(KNN.dir, "predicted_ataclabel_foxo1_sel_v3.pdf"), width=7.5, height=2.5)
plot(p)
dev.off()

library("plotgardener")
library("org.Hs.eg.db")
library("TxDb.Hsapiens.UCSC.hg19.knownGene")



plotGenes(chrom='chr11', chromstart = 17725000, chromend = 17745000)
chrom= 'chr11'
start = 17725000
end = 17745000
ensdb=EnsDb.Hsapiens.v75
gr <- GRanges(seqnames = str_replace(chrom, 'chr', ''), IRanges(start, end), strand = "*")
filters = AnnotationFilterList(GRangesFilter(gr), GenebiotypeFilter('protein_coding'))
# 
genes = autoplot(ensdb, filters, names.expr = "gene_name") + geom_text(size = 3)+
  theme(plot.margin = margin(t = 50, r = 10, b = 10, l = 10, unit = "pt"))
# 
genes_plot = genes@ggplot +
  xlim(start, end) +
  theme_classic()
gene_info <- genes(ensdb, filter = GeneNameFilter("MYOD1"))
gene_info <- genes(ensdb, filter = GeneNameFilter("FOXO1"))
#################################
## euclidean distance to control
##############################
## get euclidean distance
method.level <- c("Annotation",'scglue(multiome)','scglue',"Seurat(CCA)",'scVI','BindSC','simba','scBridge',
                  'Portal','scJoint',"Random")
peak_out <- read.csv(paste0(KNN.dir, "predicted_ataclabel_foxo1value.csv"))
peak_out_train <- read.csv("/path/to/data/RMS/lca_celltype/Mast607A/predicted_ataclabel_foxo1value.csv")
peak_out <- rbind(peak_out_train, peak_out)
peak_out<-pivot_longer(peak_out, 
                       cols = c('Mesoderm', 'Myoblast','Myocyte'), 
                       names_to ="type")%>%
  mutate(method =case_when(pipeline == "Seurat.CCA."~ "Seurat(CCA)",
                           pipeline == "Seurat.WNN."~ "Seurat(WNN)",
                           pipeline == "scglue.multiome."~ "scglue(multiome)",
                           pipeline == "annotation"~ "Annotation",
                           pipeline == "random"~ "Random",
                           pipeline == "lca_label"~ "Multiome",
                           
                           TRUE ~ pipeline))%>%
  mutate(method = factor(method, levels=c(method.level)))%>%
  drop_na()
getEucDistance <- function(vect1, vect2) sqrt(sum((vect1 - vect2)^2))
getEucDistance <- function(vect1, vect2) {
  cosine_dist = 1- cosine(vect1, vect2)
  return(cosine_dist)
}

## euc distance for different cell type for each pipeline
type_eucdist <- function(peak_merge, celltype){
  peak_merge_type <- peak_merge[which(peak_merge$type == celltype),]
  euc_dist <- getEucDistance(peak_merge_type$value.x,peak_merge_type$value.y)
  return(c(as.character(peak_merge_type$method.y[1]),celltype,euc_dist))
}
## for different pipeline
sum_Eucdist <-function(peak_out, peak_ctrl, pipeline){
  peak_pred <- peak_out[which(peak_out$method == pipeline),  ]
  peak_merge <- merge(peak_ctrl, peak_pred, by = c("position", "type"))
  peak_merge[is.na(peak_merge)] <- 0
  type <- unique(peak_ctrl$type) ## add
  results <- lapply(type, function(celltype) {
    type_eucdist(peak_merge, celltype)
  })
  return(results)
}


peak_ctrl <- peak_out[which(peak_out$method == "Annotation"),  ]
## get eucdistance sum
euc_dists <- lapply(method.level[-1], function(pipeline) {
  sum_Eucdist(peak_out, peak_ctrl, pipeline)
})

df_eucdists <- as.data.frame(matrix(unlist(euc_dists), ncol = 3, byrow = TRUE))
colnames(df_eucdists) <- c("method", "type", "euc.dist")  

df_eucdists <- df_eucdists %>%
  mutate(euc.dist = as.numeric(df_eucdists$euc.dist))%>%
  mutate(method = factor(method, levels = method.level[-1]))
## plot euclidean score
df_eucdists$euc.dist[df_eucdists$euc.dist > 50] <- 50
df_eucdists$euc.dist[df_eucdists$euc.dist > 0.2] <- 0.2

library("RColorBrewer")
df_eucdists%>%
  ggplot(mapping=aes(x=method,y = type, fill = as.numeric(df_eucdists$euc.dist)) )+
  geom_tile(color = 'black')+
  xlab(label = "")+
  ylab(label = "")+
  scale_fill_gradientn(
    breaks = c(0, 25, 50),
    #labels = scales::number_format(suffix = "+"),
    colors = rev(brewer.pal(9, "Reds")),
    name = "Euclidean distance to ground truth" )+
  theme_classic()+
  theme(text=element_text(size=16, color = "black"),
        plot.title = element_text(size = 16, color = "black", face = "bold"),
        axis.text.y = element_text(size = 16, color = "black"),
        strip.text = element_text(
          size = 16),
        axis.text.x = element_text(angle = 45, vjust = 0.85, hjust=0.85,size = 16, color = "black"),
        legend.position="top")->p_euc
pdf(paste0(KNN.dir,"predicted_ataclabel_foxo1_eucdistv2.pdf"),width = 9, height =4.5) 
p_euc
dev.off()

df_eucdists%>%
  ggplot(mapping=aes(x=method,y = type, fill = as.numeric(euc.dist)) )+
  geom_tile(color = 'black')+
  xlab(label = "")+
  ylab(label = "")+
  scale_fill_gradientn(
    breaks = c(0, 25, 50),
    #labels = scales::number_format(suffix = "+"),
    colors = rev(brewer.pal(9, "Reds")),
    name = "Euclidean distance to ground truth" )+
  theme_classic()+
  theme(text=element_text(size=16, color = "black"),
        plot.title = element_text(size = 16, color = "black", face = "bold"),
        axis.text.y = element_text(size = 16, color = "black"),
        strip.text = element_text(
          size = 16),
        axis.text.x = element_text(angle = 45, vjust = 0.85, hjust=0.85,size = 16, color = "black"))+
  theme(legend.position="top")->p_euc

pdf(paste0(KNN.dir,"predicted_ataclabel_foxo1_eucdist.pdf"),width = 9, height =4.5) 
p_euc
dev.off()



