# NOTE: paths below are placeholders. See config/config.R and the README
# for the roots you must set (DATA_ROOT, PROJECT_ROOT, TOOLS_ROOT, REF_ROOT).
## from: 
## for GRCh38

########################################################
rm(list = ls())
gc()

library(Matrix)
library(BiocGenerics)
library(GenomicRanges)
library(IRanges)

#' Extend
#'
#' Resize GenomicRanges upstream and or downstream.
#' From \url{https://support.bioconductor.org/p/78652/}
#'
#' @param x A range
#' @param upstream Length to extend upstream
#' @param downstream Length to extend downstream
#' @param from.midpoint Count bases from region midpoint,
#' rather than the 5' or 3' end for upstream and downstream
#' respectively.
#'
#' @importFrom GenomicRanges trim
#' @importFrom BiocGenerics start strand end width
#' @importMethodsFrom GenomicRanges strand start end width
#' @importFrom IRanges ranges IRanges "ranges<-"
#' @export
#' @concept utilities
#' @return Returns a \code{\link[GenomicRanges]{GRanges}} object
#' @examples
#' Extend(x = blacklist_hg19, upstream = 100, downstream = 100)
Extend <- function(
    x,
    upstream = 0,
    downstream = 0,
    from.midpoint = FALSE
) {
  if (any(strand(x = x) == "*")) {
    warning("'*' ranges were treated as '+'")
  }
  on_plus <- strand(x = x) == "+" | strand(x = x) == "*"
  if (from.midpoint) {
    midpoints <- start(x = x) + (width(x = x) / 2)
    new_start <- midpoints - ifelse(
      test = on_plus, yes = upstream, no = downstream
    )
    new_end <- midpoints + ifelse(
      test = on_plus, yes = downstream, no = upstream
    )
  } else {
    new_start <- start(x = x) - ifelse(
      test = on_plus, yes = upstream, no = downstream
    )
    new_end <- end(x = x) + ifelse(
      test = on_plus, yes = downstream, no = upstream
    )
  }
  IRanges::ranges(x = x) <- IRanges::IRanges(start = new_start, end = new_end)
  x <- trim(x = x)
  return(x)
}

find_geneact <- function(peak.df, annotation.file, seq.levels, upstream = 2000, downstream = 0, verbose = FALSE){
  # peak.df is the regions
  peak = peak.df
  # reformualte peak.df of the form "chromosome", "start", "end"
  peak.df <- do.call(what = rbind, args = strsplit(x = peak.df, split = "_"))
  peak.df <- as.data.frame(x = peak.df)
  colnames(x = peak.df) <- c("chromosome", "start", "end")
  
  # peak.df -> peaks.gr
  peaks.gr <- GenomicRanges::makeGRangesFromDataFrame(df = peak.df)
  BiocGenerics::start(peaks.gr[BiocGenerics::start(peaks.gr) == 0, ]) <- 1
  
  # gtf stores the annotation (reference genome)
  gtf <- rtracklayer::import(con = annotation.file)
  gtf <- GenomeInfoDb::keepSeqlevels(x = gtf, value = seq.levels, pruning.mode = "coarse")
  if (!any(GenomeInfoDb::seqlevelsStyle(x = gtf) == GenomeInfoDb::seqlevelsStyle(x = peaks.gr))) {
    GenomeInfoDb::seqlevelsStyle(gtf) <- GenomeInfoDb::seqlevelsStyle(peaks.gr)
  }
  # gtf.genes stores the genes 
  gtf.genes <- gtf[gtf$type == "gene"]
  
  # update the regions correspond to each gtf.genes, gtf.body_prom
  gtf.body_prom <- Extend(x = gtf.genes, upstream = upstream, downstream = downstream)
  
  # assign peaks.gr to nearest gene region
  gene.distances <- GenomicRanges::distanceToNearest(x = peaks.gr, subject = gtf.body_prom)
  # only leave the ones(regions) overlap with the gene regions(distance = 0)
  keep.overlaps <- gene.distances[rtracklayer::mcols(x = gene.distances)$distance == 
                                    0]
  peak.ids <- peaks.gr[S4Vectors::queryHits(x = keep.overlaps)]
  gene.ids <- gtf.genes[S4Vectors::subjectHits(x = keep.overlaps)]
  gene.ids$gene_name[is.na(gene.ids$gene_name)] <- gene.ids$gene_id[is.na(gene.ids$gene_name)]
  peak.ids$gene.name <- gene.ids$gene_name
  peak.ids <- as.data.frame(x = peak.ids)
  peak.ids$peak <- peak[S4Vectors::queryHits(x = keep.overlaps)]
  # new annotations should include peaks and corresponding gene.name
  annotations <- peak.ids[, c("peak", "gene.name")]
  
  return(annotations)
}

################################################################################
## start running here
# hyper-parameters
species <- "Human"
# upstream region size (base-pair)
upstream <- 2000
# downstream region size (base-pair)
downstream <- 0


path <- '/path/to/data/HTAN/HT243B1-S1H4/scDART/'
# regions chrX_start_end

#####################
## add how to generate regions file


data.dir <- "/path/to/data/HTAN/HT243B1-S1H4/"
load(paste0(data.dir,"HT243_S1H4_seurat.RData"))
regions<-rownames(all.atac$ATAC$counts)
regions<-gsub("[:|-]", "_", regions)
writeLines(regions,  paste0(path,"regions.txt"))
#####################




regions <- read.table(file = paste0(path, "regions.txt"), sep = ",", header = FALSE)[[1]]

if(species == "Mouse"){
  A = find_geneact(peak.df = regions, annotation.file = "~/Dropbox (GaTech)/Research/Projects/pipeline_integration/reference_genome/Mus_musculus.GRCm38.84.gtf", 
                   seq.levels = c(1:19, "X", "Y"), upstream = upstream, downstream = downstream, verbose = TRUE)
} else if(species == "Human"){
  A = find_geneact(peak.df = regions, annotation.file = "/path/to/references/GRCh38-2020-A_build/gencode.v32.primary_assembly.annotation.filtered.gtf", 
                   seq.levels = paste0("chr",c(1:22, "X", "Y")), upstream = upstream, downstream = downstream, verbose = TRUE)
} else{
  stop("species can only be Human or Mouse")
}

# output gene activity matrix
write.table(A, file = paste0(path, "gact.csv"), sep = ",")

A <- read.csv(paste0(path, "gact.csv"))
## generate the region2gene csv in tutorial
#### remove the gact genes that not in the RNA matrix



genes<-rownames(all.rna$RNA$counts)
## filter top 10000 hvg due to memory
all.rna <- FindVariableFeatures(all.rna,nfeatures = 10000)
genes<-VariableFeatures(all.rna)
peak <- rownames(all.atac$ATAC$counts)
A_clean <- A[A$gene.name %in% genes,]
A_clean$peak <- gsub("_","-", A_clean$peak,fixed = TRUE)
t<-table(A_clean$peak,A_clean$gene.name)


write.csv(colnames(t), paste0(path,"gact_gene.csv"))
write.csv(rownames(t),paste0(path,"gact_peak.csv"))

# ## first check if t works?
t <- Matrix(t, sparse = TRUE)
writeMM(t, file = paste0( path, "region2gene.mtx"), colnames=colnames(t), rownames=rownames(t))
# 
# unmatch_genes <- genes[!genes %in% A_clean$gene.name]
# unmatch_genes.mtx <- matrix(0,
#                             nrow = length(A_clean$peak),
#                             ncol= length(unmatch_genes))
# rownames(unmatch_genes.mtx)<-A_clean$peak
# colnames(unmatch_genes.mtx)<-unmatch_genes
# unmatch_genes.mtx <- Matrix(unmatch_genes.mtx, sparse = TRUE)
# t2 <- cbind(t, unmatch_genes.mtx)
# 
# ## for features not present in coarse gam matrix
# unmatch_peak <- peak[!peak %in% A_clean$peak]
# unmatch_peak.mtx <- matrix(0,
#                            nrow = length(unmatch_peak),
#                            ncol= length(colnames(t2)))
# rownames(unmatch_peak.mtx)<-unmatch_peak
# colnames(unmatch_peak.mtx)<-colnames(t2)
# unmatch_peak.mtx <- Matrix(unmatch_peak.mtx, sparse = TRUE)
# t3 <- rbind(t2, unmatch_peak.mtx)
# t3<- t3[peak, genes]
# writeMM(t3, file = paste0(path, "pbmc3k_region2gene_t3.mtx"), colnames=colnames(t), rownames=rownames(t))

