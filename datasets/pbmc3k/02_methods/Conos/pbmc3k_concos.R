# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
#install.packages('devtools')
#devtools::install_github('kharchenkolab/conos')

## https://pklab.med.harvard.edu/peterk/conos/atac_rna/example.html
#install.packages('conos')
#install.packages('pagoda2')
library(conos)
library(pagoda2)
library(parallel)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
###############################################
bcmtx.h5 <- "/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix_test.h5"
fpath <-'/path/to/data/pbmc3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz'
n_lat <-30
project <- "pbmc3k"

outdir <- paste0("res_",project,"/")
if (!file.exists(outdir)) {
  dir.create(outdir)
}

###############################################
# load(url("http://pklab.med.harvard.edu/peterk/conos/atac_rna/data.RData"))
# p2l <- mclapply(data,basicP2proc,n.odgenes=3e3,min.cells.per.gene=-1,nPcs=30,make.geneknn=F,n.cores=30,mc.cores=1)
# l.con <- Conos$new(p2l,n.cores=30)
# l.con$buildGraph(k=15,k.self=5,k.self.weigh=0.01,ncomps=30,n.odgenes=5e3,space='PCA') 
# l.con$findCommunities(resolution=1.5)
# l.con$embedGraph(alpha=1/2)

t1 <- Sys.time()
inputdata.10x <- Read10X_h5(bcmtx.h5)
rna_counts  <- inputdata.10x$`Gene Expression`
atac_counts  <- inputdata.10x$Peaks

#bc<-colnames(rna_counts)
#atac_counts<-atac_counts[,bc]
# create atac object
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  fragments= fpath
)

all.atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# convert to UCSC style
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
Annotation(all.atac) <- annotations

# ## create RNA seurat object
# all.rna <- CreateSeuratObject(counts = rna_counts,assay = "RNA")
# all.rna <- NormalizeData(all.rna)
# all.rna <- FindVariableFeatures(all.rna)
#genes <- row.names(rna_counts)
gene.activities <- GeneActivity(all.atac)
# add gene activities as a new assay
all.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
# normalize gene activities
DefaultAssay(all.atac) <- "ACTIVITY"
activity.count <- all.atac@assays[["ACTIVITY"]]@counts
## renamed bc
colnames(rna_counts)<-paste0(colnames(rna_counts),"_rna")
colnames(activity.count)<-paste0(colnames(activity.count),"_atac")

pbmc3k.data<-list(rna=rna_counts, atac=activity.count)  
pbmc3k.p2l <- mclapply(pbmc3k.data,basicP2proc,n.odgenes=3e3,min.cells.per.gene=-1,nPcs=30,make.geneknn=F,n.cores=30,mc.cores=1)

pbmc3k.con <- Conos$new(pbmc3k.p2l,n.cores=30)
pbmc3k.con$buildGraph(k=15,k.self=5,k.self.weigh=0.01,ncomps=30,n.odgenes=5e3,space='PCA') 
pbmc3k.con$findCommunities(resolution=1.5)
pbmc3k.con$embedGraph(alpha=1/2,target.dims=10)
t2 <- Sys.time()
 
df_umap = as.data.frame(pbmc3k.con[["embedding"]])
colnames(df_umap) = paste0("latent_",1:ncol(df_umap))

write.csv(df_umap,paste0(outdir,"coembed_coor.csv"))
write.table(difftime(t2, t1, units = "secs")[[1]], 
            file = file.path(outdir,"runtime.txt"), 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
######################################################################
## clustering in conos

findres<-function(con.env, nclust){
  tmp <- con.env
  iteration =0
  resolutions = c(0,100)
  obtained_nclust = -1
  while ((obtained_nclust != nclust) & (iteration < 100)) {
    res = sum(resolutions)/2
    tmp$findCommunities(resolution=res)
    obtained_nclust<-length(levels(tmp[["clusters"]]$leiden$groups))
    #print(obtained_nclust)
    if (obtained_nclust < nclust){
      resolutions[1] = res
    }else{
      resolutions[2]=res
    }
    iteration = iteration + 1
    #print(iteration)
  }
  return(res)
}

res = findres(pbmc10k.con,10)
pbmc3k.con$findCommunities(resolution=res)

df_clust<- as.data.frame (pbmc3k.con[["clusters"]]$leiden$groups)
write.csv(df_clust,paste0(outdir,"conos_nclust.csv"))

