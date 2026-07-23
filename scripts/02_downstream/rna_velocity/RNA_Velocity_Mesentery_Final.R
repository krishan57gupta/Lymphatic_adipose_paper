library(Seurat)
library(ggplot2)
library(ggpubr)
library(export)
library(Hmisc)
library(stringr)
library(MASS)
library(scales)
library(ggbreak) 

MTper=30
res=.1
res1=res*10^(nchar(as.character(res))-2)
dataType="Mesentery"
sampleNames=c("WT-HFD","LOF-HFD")
sn="/combined"
for(i in sampleNames){
  sn=paste(sn,i,sep="_")
}
sn=paste(sn,MTper,res1,sep="_")
sn=paste(sn,'/',sep="")
print(sn)
sampleList=list("WT-HFD"="WT-HFD-Mes",
                "LOF-HFD"="LOF-HFD-Mes")
mainDir <-  '/yourDataAndCodeFolder/'
dataDir <- paste(mainDir,'data/',sep="")

newDir=paste(mainDir,'loomFiles/RNA_Velocity/',dataType,sn,sep="")
plotDir <- paste(mainDir,'analysis/SCA_plots/',dataType,sn,sep="")
processedDataDir <- paste(mainDir,'analysis/processed_files/',dataType,sn,sep="")

if (!file.exists(plotDir)){
  dir.create(plotDir,recursive = TRUE)
}
if (!file.exists(paste(plotDir,'/otherPlots',sep=""))){
  dir.create(paste(plotDir,'/otherPlots',sep=""),recursive = TRUE)
}
if (!file.exists(processedDataDir))
  dir.create(processedDataDir,recursive = TRUE)

if (!file.exists(newDir))
  dir.create(newDir,recursive = TRUE)

Integrated="RNA"

seurat_obj<-readRDS(paste(processedDataDir,'combined_',dataType,'_ABC_F.rds',sep=''))


seurat_obj$barcode <- colnames(seurat_obj)
# seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
# seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
seurat_obj@meta.data$CMcor=NULL
seurat_obj@meta.data$CFcor=NULL
write.csv(seurat_obj@meta.data, file=paste0(newDir,'metadataF.csv'), quote=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')

writeMM(counts_matrix, file=paste0(newDir,'countsF.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seurat_obj@reductions$RNA_pca@cell.embeddings, file=paste0(newDir,'RNA_pcaF.csv'), quote=F)
write.csv(seurat_obj@reductions$ABC_pca@cell.embeddings, file=paste0(newDir,'ABC_pcaF.csv'), quote=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),paste0(newDir,file='gene_namesF.csv'),
  quote=F,row.names=F,col.names=F
)
write.table(
  data.frame('barcode'=colnames(counts_matrix)),paste0(newDir,file='cellBarcodesF.csv'),
  quote=F,row.names=F,col.names=F
)


