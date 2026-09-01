library(Seurat)
library(Matrix)
library(rhdf5)
#https://rdrr.io/github/MarioniLab/DropletUtils/man/write10xCounts.html

## functions to prepare customize 10X h5
## subset barcode, for train and test data
sub_h5<-function(barcodes,feature_mtx, count_mtx,path,test=TRUE,module=NULL){
  ## create a new h5 file
  if (file.exists(path)){
    unlink(path, recursive = TRUE)
    h5createFile(path)
  } else {
    h5createFile(path)
  }
  
  group<-"matrix"
  h5createGroup(path, group)
  h5write(barcodes, file=path, name=paste0(group, "/barcodes")) ## input: barcodes
  
  # Saving feature information.
  h5createGroup(path, file.path(group, "features"))
  ## from feature matrix:
  features<-read.table(feature_mtx,sep='\t')
  if (is.null(module)){
    gene.id<-features$V1
    gene.symbol<-features$V2
    gene.type<-features$V3
  }else if ((module == 'Gene Expression') | (module == 'Peaks')){
    features <- features[features$V3==module,]
    gene.id<-features$V1
    gene.symbol<-features$V2
    gene.type<-features$V3
  }else{
    stop("stop: wrong module input")
  }
  
  h5write(gene.id, file=path, name=paste0(group, "/features/id"))
  h5write(gene.symbol, file=path, name=paste0(group, "/features/name"))
  h5write(gene.type,file=path, name=paste0(group, "/features/feature_type"))
  h5write("genome", file=path, name=paste0(group, "/features/_all_tag_keys")) ## prob alright
  genome ="GRCh38"
  h5write(rep(genome, length.out=length(gene.id)),
          file=path, name=paste0(group, "/features/genome"))
  
  # Writing attributes***okay to skip
  h5f <- H5Fopen(path)
  h5g <- H5Gopen(h5f, "/")
  chemistry="Single Cell 3' v3" ## v3 is default
  original.gem.groups=1L
  library.ids="custom"
  h5writeAttribute(chemistry, h5obj=h5g, name="chemistry_description", variableLengthString=TRUE, asScalar=TRUE, encoding ="UTF8")
  h5writeAttribute("matrix", h5obj=h5g, name="filetype", variableLengthString=TRUE, asScalar=TRUE, encoding="UTF8")
  h5writeAttribute(library.ids, h5obj=h5g, name="library_ids", variableLengthString=TRUE, asScalar=TRUE, encoding="UTF8")
  h5writeAttribute(original.gem.groups, h5obj=h5g, name="original_gem_groups")
  version = 3
  h5writeAttribute(as.integer(version) - 1L, h5obj=h5g, name="version") # this is probably correct.
  H5Gclose(h5g)
  H5Fclose(h5f)
  
  # Saving matrix information.
  x <- as(count_mtx, "CsparseMatrix")
  #print(gene.symbol)
  #print(barcodes)
  x<-x[gene.symbol,barcodes] ## re-order accord feature matrix
  h5write(x@x, file=path, name=paste0(group, "/data"))
  h5write(dim(x), file=path, name=paste0(group, "/shape"))
  h5write(x@i, file=path, name=paste0(group, "/indices")) # already zero-indexed.
  h5write(x@p, file=path, name=paste0(group, "/indptr"))
  
  ## test new10X_h5 using seurat function
  if(test==TRUE){
    new_10xh5<-Read10X_h5(path)
    return(new_10xh5)
  }else{
    print(paste0(path," has been created!"))
  }
}


###################################################
## function to create the custom_filter_feature_bc_mtx

sub_bc_matrix<-function(barcodes,feature.mtx, count.mtx,output_dir,module=NULL){
  if (!file.exists(output_dir)) {
    dir.create(output_dir)
  }
  ## barcodes.tsv
  writeLines(barcodes, paste0(output_dir,"barcodes.tsv"))
  
  ## features.tsv
  feature.mtx <- read.table(feature.mtx,header = FALSE, sep ="\t")
  
  if (is.null(module)){
    write.table(feature.mtx, file = paste0(output_dir,"features.tsv"),sep="\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
  }else if (module == 'Gene Expression'){
    feature_mtx <- feature.mtx[feature.mtx$V3 =='Gene Expression',]
    write.table(feature_mtx, file = paste0(output_dir,"features.tsv"),sep="\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
  }else if (module =='Peaks'){
    feature_mtx <- feature.mtx[feature.mtx$V3 =='Peaks',]
    write.table(feature_mtx, file = paste0(output_dir,"features.tsv"),sep="\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
  }
  
  ## matrix.mtx
  writeMM(count.mtx, paste0(output_dir,"matrix.mtx"))
  
  print(paste0(output_dir, " has been updated!"))
}

