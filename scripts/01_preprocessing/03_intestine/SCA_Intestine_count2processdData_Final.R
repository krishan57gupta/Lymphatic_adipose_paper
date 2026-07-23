
source("/yourDataAndCodeFolder/functions.R")
library(Seurat)
library(ggplot2)
library(extrafont)
loadfonts(device = "win")
quartzFonts()
library(ggpubr)
library(export)
library(Hmisc)
library(stringr)
library(MASS)
library(scales)
library(ggbreak) 
library(ggridges)
library(Matrix)

######################### resolution of one problem
Csparse_validate = "CsparseMatrix_validate"

DoRun=0 # for running specific code if set to 1
options(future.globals.maxSize= 5147483648) #to increase max memory
MTper=30
res=.1
res1=res*10^(nchar(as.character(res))-2)
dataType="Intestine"
sampleNames=c("WT-HFD","LOF-HFD")
sn="/combined"
for(i in sampleNames){
  sn=paste(sn,i,sep="_")
}
sn=paste(sn,MTper,res1,sep="_")
sn=paste(sn,'/',sep="")
print(sn)

sampleList=list("WT-HFD"="WT-HFD-Gut",
                "LOF-HFD"="LOF-HFD-Gut")
mainDir <- "/yourDataAndCodeFolder/"
dataDir <- paste(mainDir,'data/',sep="")

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

redix=10
log_breaks = function(maj, radix=redix) {
  function(x) {
    minx         = floor(min(logb(x,radix), na.rm=T)) - 1
    maxx         = ceiling(max(logb(x,radix), na.rm=T)) + 1
    n_major      = maxx - minx + 1
    major_breaks = seq(minx, maxx, by=1)
    if (maj) {
      breaks = major_breaks
    } else {
      steps = logb(1:(radix-1),radix)
      breaks = rep(steps, times=n_major) +
        rep(major_breaks, each=radix-1)
    }
    radix^breaks
  }
}

Integrated="RNA"
width=3
height=3
pointSize=.5
family="sans" # for "Helvetica"
txtSize=10
plotTitleSize=10
axisSize=10
axisTitleSize=10
legendSize=10
legendTitleSize=10
dpi=600
ext='.pdf'
color1=c("blue","red","darkgreen","yellow","cyan","orange","magenta",'skyblue','purple','brown',
         'darkblue','green','grey','tan','pink','violet','black')
color2=c("green","red","magenta","orange","grey","pink","tan",'skyblue','lightcyan','violet',
         'darkblue','darkgreen','cyan','yellow','blue','brown','purple','black')
colorCode=c("#FFA500","#00AA00","#00FFAA")
colorCode=c("#FF6A00","#00AA00")
# c("#FF6A00","#FF7A00","#FF8A00")


plotList=list()
objectList=list()
objectList1=list()
for(condName in sampleNames){
  print(condName)
  combined <- Read10X(data.dir = paste(dataDir,sampleList[[condName]],'/outs/filtered_feature_bc_matrix/',sep=""))
  
  combined = CreateSeuratObject(counts = combined, project = paste(condName,"_0",sep=""),min.cells = 3, min.features = 200)
  combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")
  combined<-NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
  combined<-FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
  combined <- AddMetaData(object = combined,metadata = condName, col.name = 'conditions')
  combined <- AddMetaData(object = combined,metadata = paste(condName,"_0",sep=""), col.name = 'conditionsBatched')
  
  combined@meta.data[['CMcor']]<-paste('\n',condName,',_0\ncor = ',round(cor(combined@meta.data$nCount_RNA,combined@meta.data$percent.mt, method = 'pearson'),3),sep='') # CM means cor and MT
  combined@meta.data[['CFcor']]<-paste('\n',condName,',_0\ncor = ',round(cor(combined@meta.data$nCount_RNA,combined@meta.data$nFeature_RNA, method = 'pearson'),3),sep='') # CM means cor and Features
  
  objectList[[paste(condName,"_0",sep="")]]=combined
  
  ###########################################################################################################
  combined = CreateSeuratObject(counts = combined@assays$RNA, project = paste(condName,"_0",sep=""),min.cells = 3, min.features = 200)
  combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")
  combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < MTper)
  combined<-NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
  combined<-FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
  combined <- AddMetaData(object = combined,metadata = condName, col.name = 'conditions')
  combined <- AddMetaData(object = combined,metadata = paste(condName,"_0",sep=""), col.name = 'conditionsBatched')
  
  combined@meta.data[['CMcor']]<-paste('\n',condName,',_0\ncor = ',round(cor(combined@meta.data$nCount_RNA,combined@meta.data$percent.mt, method = 'pearson'),3),sep='') # CM means cor and MT
  combined@meta.data[['CFcor']]<-paste('\n',condName,',_0\ncor = ',round(cor(combined@meta.data$nCount_RNA,combined@meta.data$nFeature_RNA, method = 'pearson'),3),sep='') # CM means cor and Features
  
  objectList1[[paste(condName,"_0",sep="")]]=combined
  
  ###########################################################################################################
  #######################################               new samples.              ###########################
  ###########################################################################################################
  combined <- Read10X(data.dir = paste(dataDir,"06_01_2023/",sampleList[[condName]],'/outs/filtered_feature_bc_matrix/',sep=""))
  
  combined = CreateSeuratObject(counts = combined, project = paste(condName,"_1",sep=""),min.cells = 3, min.features = 200)
  combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")
  combined<-NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
  combined<-FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
  combined <- AddMetaData(object = combined,metadata = condName, col.name = 'conditions')
  combined <- AddMetaData(object = combined,metadata = paste(condName,"_1",sep=""), col.name = 'conditionsBatched')
  
  combined@meta.data[['CMcor']]<-paste('\n',condName,',_1\ncor = ',round(cor(combined@meta.data$nCount_RNA,combined@meta.data$percent.mt, method = 'pearson'),3),sep='') # CM means cor and MT
  combined@meta.data[['CFcor']]<-paste('\n',condName,',_1\ncor = ',round(cor(combined@meta.data$nCount_RNA,combined@meta.data$nFeature_RNA, method = 'pearson'),3),sep='') # CM means cor and Features
  
  objectList[[paste(condName,"_1",sep="")]]=combined
  
  ###########################################################################################################
  combined = CreateSeuratObject(counts = combined@assays$RNA, project = paste(condName,"_1",sep=""),min.cells = 3, min.features = 200)
  combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")
  combined <- subset(combined, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < MTper)
  combined<-NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
  combined<-FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
  combined <- AddMetaData(object = combined,metadata = condName, col.name = 'conditions')
  combined <- AddMetaData(object = combined,metadata = paste(condName,"_1",sep=""), col.name = 'conditionsBatched')
  
  combined@meta.data[['CMcor']]<-paste('\n',condName,',_1\ncor = ',round(cor(combined@meta.data$nCount_RNA,combined@meta.data$percent.mt, method = 'pearson'),3),sep='') # CM means cor and MT
  combined@meta.data[['CFcor']]<-paste('\n',condName,',_1\ncor = ',round(cor(combined@meta.data$nCount_RNA,combined@meta.data$nFeature_RNA, method = 'pearson'),3),sep='') # CM means cor and Features
  
  objectList1[[paste(condName,"_1",sep="")]]=combined
}

##########################################################################################  unfiltered graphs
tempObject<-objectList[[names(objectList)[1]]]
tempObjectList=objectList
tempObjectList[[names(objectList)[1]]]<-NULL
combined <- merge(tempObject, y = tempObjectList, add.cell.ids=names(objectList),project = dataType)
combined$conditions<-factor(x = combined$conditions, levels =sampleNames)

##########################################################################################

features <- SelectIntegrationFeatures(object.list = objectList1)
anchors <- FindIntegrationAnchors(object.list = objectList1, anchor.features = features)
combined <- IntegrateData(anchorset = anchors)
combined$conditions<-factor(x = combined$conditions, levels =sampleNames)

########################################################################################## Low dim embedding and clustering with ABC
DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)

combined$ABC_PCA_1 <- combined@reductions$pca@cell.embeddings[,1]
combined$ABC_PCA_2 <- combined@reductions$pca@cell.embeddings[,2]

combined$ABC_UMAP_1 <- combined@reductions$umap@cell.embeddings[,1]
combined$ABC_UMAP_2 <- combined@reductions$umap@cell.embeddings[,2]

combined@reductions$ABC_pca=combined@reductions$pca
combined@reductions$ABC_umap=combined@reductions$umap

combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = res)
combined <- AddMetaData(object = combined,metadata = combined@meta.data$seurat_clusters, col.name = 'ABC_clusters')

########################################################################################## 
cellTypesMarkersList =list()
cellTypesMarkersList[["ABC_cellTypes"]]=list(  # not performed
  '0'="Dividing cells",
  '1'="Enterocytes",
  '2'="Macrophages",
  '3'="B cells",
  '4'="Paneth cells",
  '5'='Pre-Collecting LECs',
  '6'="T cells",
  '7'="Capillary LECs",
  '8'="Goblet cells",
  '9'="Capillary LECs",
  '10'="Neurons",
  '11'='Valve LECs',
  '12'='Capillary LECs')
cellTypesMarkersList[["ABC_subCellTypes"]]=list(  # not performed
  '0'="DC",
  '1'="Ent",
  '2'="MP",
  '3'="BC",
  '4'="PC",
  '5'='Pre-Col LEC',
  '6'="TC",
  '7'="Cap LEC1",
  '8'="GC",
  '9'="Cap LEC2",
  '10'="Neu",
  '11'='Val LEC',
  '12'='Cap LEC3')
cellTypesMarkersList[["ABC_subCellTypes1"]]=list(  # not performed
  '0'="DC",
  '1'="Ent",
  '2'="MP",
  '3'="BC",
  '4'="PC",
  '5'='Pre-Col LEC',
  '6'="TC",
  '7'="Cap LEC",
  '8'="GC",
  '9'="Cap LEC",
  '10'="Neu",
  '11'='Val LEC',
  '12'='Cap LEC')

subCellTypeOrder=c('Cap LEC1','Cap LEC2','Cap LEC3',"Pre-Col LEC","Val LEC","DC","Ent","GC","PC","MP","BC","TC","Neu")
subCellTypeOrder1=c('Cap LEC',"Pre-Col LEC","Val LEC","DC","Ent","GC","PC","MP","BC","TC","Neu")
cellTypeOrder=c('Capillary LECs',"Pre-Collecting LECs","Valve LECs","Dividing cells","Enterocytes","Goblet cells","Paneth cells","Macrophages","B cells","T cells","Neurons")
subCellTypeOrderSelected=c('Cap LEC1','Cap LEC2','Cap LEC3',"Pre-Col LEC","Val LEC")

Idents(combined)<-combined@meta.data$ABC_clusters
combined<-RenameIdents(combined, cellTypesMarkersList$ABC_cellTypes)
combined@meta.data$ABC_cellTypes<-Idents(combined)
combined@meta.data$ABC_Cond_CT<-paste(combined@meta.data$ABC_cellTypes,combined@meta.data$conditions,sep="_")

Idents(combined)<-combined@meta.data$ABC_clusters
combined<-RenameIdents(combined, cellTypesMarkersList$ABC_subCellTypes)
combined@meta.data$ABC_subCellTypes<-Idents(combined)
combined@meta.data$ABC_Cond_SCT<-paste(combined@meta.data$ABC_subCellTypes,combined@meta.data$conditions,sep="_")

Idents(combined)<-combined@meta.data$ABC_clusters
combined<-RenameIdents(combined, cellTypesMarkersList$ABC_subCellTypes1)
combined@meta.data$ABC_subCellTypes1<-Idents(combined)
combined@meta.data$ABC_Cond_SCT1<-paste(combined@meta.data$ABC_subCellTypes1,combined@meta.data$conditions,sep="_")

combined@meta.data$ABC_subCellTypes=factor(combined@meta.data$ABC_subCellTypes,levels = subCellTypeOrder)
combined@meta.data$ABC_subCellTypes1=factor(combined@meta.data$ABC_subCellTypes1,levels = subCellTypeOrder1)
combined@meta.data$ABC_cellTypes=factor(combined@meta.data$ABC_cellTypes,levels = cellTypeOrder)
saveRDS(combined,paste(processedDataDir,'combined_',dataType,'_ABC_F.rds',sep=''))

