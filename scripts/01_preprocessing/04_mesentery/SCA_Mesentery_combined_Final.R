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

######################### resolution of one problem
dataType="Mesentery"
sampleNames=c("WT-HFD","LOF-HFD")
sampleList=list("WT-HFD"="WT-HFD-Mes",
                "LOF-HFD"="LOF-HFD-Mes")

Csparse_validate = "CsparseMatrix_validate"

DoRun=0 # for running specific code if set to 1
options(future.globals.maxSize= 5147483648) #to increase max memory
MTper=30
res=.1
res1=res*10^(nchar(as.character(res))-2)
sn="/combined"
for(i in sampleNames){
  sn=paste(sn,i,sep="_")
}
sn=paste(sn,MTper,res1,sep="_")
sn=paste(sn,'/',sep="")
print(sn)

mainDir <- '/yourDataAndCodeFolder/'
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
}

##########################################################################################  unfiltered graphs
tempObject<-objectList[[names(objectList)[1]]]
tempObjectList=objectList
tempObjectList[[names(objectList)[1]]]<-NULL
combined <- merge(tempObject, y = tempObjectList, add.cell.ids=names(objectList),project = dataType)
combined$conditions<-factor(x = combined$conditions, levels =sampleNames)
##########################################################################################

#########################################################################################  Integration and batch correction
features <- SelectIntegrationFeatures(object.list = objectList1)
anchors <- FindIntegrationAnchors(object.list = objectList1, anchor.features = features)
combined <- IntegrateData(anchorset = anchors)
combined$conditions<-factor(x = combined$conditions, levels =sampleNames)
##########################################################################################
########################################################################################## Low dim embedding and clustering with RNA
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
cellTypesMarkersList[["ABC_cellTypes"]]=list(  # not performed
  '0'="Enterocytes",
  '2'="Cap/Pre-Collecting LECs",
  '3'="Collecting LECs",
  '4'="Macrophages",
  '5'='Cap/Pre-Collecting LECs',
  '6'="Paneth cells",
  '7'="B cells",
  '8'="Collecting LECs", # Capillary LECs
  '9'="Valve LECs")
cellTypesMarkersList[["ABC_subCellTypes"]]=list(  # not performed
  '0'="Ent",
  '2'="Cap/Pre-Col LEC",
  '3'="Col LEC",
  '4'="MP",
  '5'='Cap/Pre-Col LEC',
  '6'="PC",
  '7'="BC",
  '8'="Col LEC", # Cap LEC
  '9'="Val LEC")
cellTypesMarkersList[["ABC_subCellTypes1"]]=list(  # not performed
  '0'="Ent",
  '2'="Cap/Pre-Col LEC",
  '3'="Col LEC",
  '4'="MP",
  '5'='Cap/Pre-Col LEC',
  '6'="PC",
  '7'="BC",
  '8'="Col LEC", # Cap LEC
  '9'="Val LEC")

subCellTypeOrder=c('Cap/Pre-Col LEC',"Col LEC","Val LEC","Ent","PC","MP","BC")
subCellTypeOrder1=c('Cap/Pre-Col LEC',"Col LEC","Val LEC","Ent","PC","MP","BC")
cellTypeOrder=c("Cap/Pre-Collecting LECs","Collecting LECs","Valve LECs","Enterocytes","Paneth cells","Macrophages","B cells")
subCellTypeOrderSelected=c('Cap/Pre-Col LEC',"Col LEC","Val LEC")



Idents(combined)<-combined@meta.data$ABC_clusters
combined=subset(x = combined, idents = c(0,2,3,4,5,6,7,8,9))
combined@meta.data$ABC_clusters=factor(combined@meta.data$ABC_clusters, levels=c(0,2,3,4,5,6,7,8,9))

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

Idents(combined)<-combined@meta.data$ABC_subCellTypes
combined=subset(x = combined, idents = c('Cap/Pre-Col LEC',"Col LEC","Val LEC","MP","BC"))
saveRDS(combined,paste(processedDataDir,'combined_',dataType,'_ABC_F.rds',sep=''))

##########################################################################
##########################################################################
####################### just for umap plot
combined=readRDS(paste(processedDataDir,'combined_',dataType,'_ABC_F.rds',sep=''))
DefaultAssay(combined) <- "RNA"
table(combined@meta.data$conditions)
minCells=min(table(combined@meta.data$conditions))
downSamples=c()
for(i in names(table(combined@meta.data$conditions))){
  downSamples=c(downSamples,sample(rownames(combined@meta.data[combined@meta.data$conditions==i,]),minCells))
}
combined=subset(combined,cells = downSamples)
table(combined@meta.data$conditions)

Idents(combined)<-combined@meta.data$ABC_subCellTypes
############################### Extended data Figure 8,a
DP<-Seurat::DimPlot(combined, reduction = "ABC_umap", pt.size=3,
                    # cols=color3,
                    # raster.dpi = c(512, 512),
                    raster = TRUE,
                    label = TRUE,
                    split.by = 'conditions')+
  theme(text = element_text(size=txtSize,family = family),
        plot.title = element_text(size = plotTitleSize, family = family),
        axis.text = element_text(size = axisSize,family = family),
        axis.title=element_text(size = axisTitleSize,family = family),
        legend.text=element_text(size=legendSize,family = family),
        legend.title=element_text(size=legendTitleSize,family = family))
print(DP)
# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*2, height=height*1, append=TRUE)

markersList=list()
markersList[["ABC_filtered1"]]=list('Cap/Pre-Collecting LECs'=c('Stmn2','Ccl21','Ano1','Fabp4','Plvap','Egfl7', 'Cdh5','Pdgfb'),
                                             'Collecting LECs'=c("Klf4","Eng",'Piezo2','Pdpn','Tppp3','Gja1'),
                                             'Valve LECs'=c('Prox1','Nfatc1','Itga9'),
                                             'MP'=c('Cxcl2', 'Lgals3','C1qb','Cd68','Lyz2'),
                                             'BC'=c('Iglv1','Mzb1','Iglc2','Cd79b'))

markersList1=markersList[["ABC_filtered1"]]
############################### Extended data Figure 8,b
HM1=DotPlot(combined, features = markersList1,
              col.min=-0.5,
              cols = c("white", "red"), 
              # split.by = 'conditions',
              dot.scale = 8) +RotatedAxis()+
  theme(text = element_text(size=txtSize,family = family),
        plot.title = element_text(size = plotTitleSize, family = family),
        axis.text = element_text(size = axisSize,family = family),
        axis.title=element_text(size = axisTitleSize,family = family),
        legend.text=element_text(size=legendSize,family = family),
        legend.title=element_text(size=legendTitleSize,family = family))
print(HM1)
# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*2.5, height=height*1, append=TRUE)

####################### Trajectory analysis (Monocle)###################################
combined5=combined
DefaultAssay(combined) <- "RNA"
cellTypes=c('Cap/Pre-Col LEC',"Col LEC","Val LEC")
plots=list()
cds_list=list()
for(num in c(1,2)){
  Idents(combined5)=combined5@meta.data$conditions
  combined6=subset(x = combined5, idents = sampleNames[num])
  Idents(combined6)<-combined6@meta.data$ABC_subCellTypes
  combined6=subset(x = combined6, idents = cellTypes)
  combined6@meta.data$ABC_subCellTypes=factor(x = combined6@meta.data$ABC_subCellTypes, levels =c('Cap/Pre-Col LEC',"Col LEC","Val LEC"))
  DefaultAssay(combined6)<-"RNA"
  # library(monocle)
  library(reticulate)
  library(ggplot2)
  library(SeuratWrappers)
  ################################################################################
  cds <- as.cell_data_set(combined6)
  cds<- monocle3::estimate_size_factors(cds)
  monocle3::fData(cds)$gene_short_name <- row.names(monocle3::fData(cds))
  dim(cds)
  
  cds <- monocle3::cluster_cells(cds)
  dim(cds)
  cds <- monocle3::learn_graph(cds,use_partition=FALSE)
  dim(cds)
  
  if(num==1){
    root_cell = colnames(cds)[cds@colData@listData$ABC_subCellTypes %in% c("Col LEC","Val LEC")]
  }else{
    root_cell = colnames(cds)[cds@colData@listData$ABC_subCellTypes %in% c('Cap/Pre-Col LEC')]
  }
  
  root_cell=na.omit(root_cell)
  cds <- monocle3::order_cells(cds,reduction_method = "UMAP",root_cells=root_cell)
  cds_list[[num]]=cds
  
  plot2=monocle3::plot_cells(cds, label_groups_by_cluster=FALSE, graph_label_size=0,label_leaves=F, color_cells_by = "pseudotime",
                             trajectory_graph_color = "red",label_branch_points = FALSE,label_roots = F, rasterize = TRUE,
                             alpha=1,cell_size = .25,trajectory_graph_segment_size = .5,group_label_size=5) # cell_size = 1,trajectory_graph_segment_size = 0.5)
  plots[[num]]=plot2
}
ggarrange(plots[[1]],plots[[2]],nrow=1,ncol=2,vjust=0.5,hjust=0.5,common.legend = TRUE)
# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*1.25, height=height*.75, append=TRUE)

##############################################
combined5=readRDS(paste(processedDataDir,'combined_',dataType,'_ABC_F.rds',sep=''))
combined5@meta.data$ABC_Cond_CT<-paste(combined5@meta.data$ABC_cellTypes,combined5@meta.data$conditions,sep="_")

Idents(combined5)<-combined5@meta.data$ABC_subCellTypes #ABC_subCellTypes

prop=table(combined5@meta.data$ABC_Cond_SCT) # ABC_Cond_SCT
cond=table(combined5@meta.data$conditions)
ct=sort(unique(combined5@meta.data$ABC_subCellTypes)) # ABC_subCellTypes

dff=matrix(0,length(cond),length(ct))
rownames(dff)=names(cond)
colnames(dff)=ct
for(i in 1:length(prop)){
  print(i)
  temp=str_split(names(prop)[i],'_',2)[[1]]
  abc=c(as.numeric(prop[i]))
  for(j in names(cond)){
    if(temp[2]==j){
      abc=as.numeric(abc)/cond[j]*100
    }
  }
  dff[temp[2],temp[1]]=abc
}

df=matrix(0,dim(dff)[1]*dim(dff)[2],3)
colnames(df)=c("cellType","condition","cellPercentage")
df=as.data.frame(df)
count=1
for(i in rownames(dff)){
  for(j in colnames(dff)){
    df[count,]=c(j,i,dff[i,j])
    count=count+1
  }
}
df$cellPercentage=as.numeric(df$cellPercentage)
df$cellPercentage=round(as.numeric(df$cellPercentage),2)

gg=list()
ct2=ct
for(j in ct2){
  print(j)
  df1=df[df$cellType==j,]
  df1$cellPercentage=as.numeric(df1$cellPercentage)
  maxP=max(df1$cellPercentage)*1.2
  df1$condition=factor(df1$condition,levels = sampleNames)
  gg[[j]]=ggplot(data=df1, aes(x=condition, y=cellPercentage)) +
    geom_bar(stat="identity", fill=colorCode[c(1,2)])+
    geom_text(aes(label=cellPercentage), position=position_dodge(width=1), size=3,vjust=0.5,hjust=0, color="black",angle=90)+theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+labs(y = "", x = j)+rremove("x.text")+ylim(0, maxP)+
    theme(text = element_text(size=txtSize,family = family),
          plot.title = element_text(size = plotTitleSize, family = family),
          axis.text = element_text(size = axisSize,family = family),
          axis.title=element_text(size = axisTitleSize,family = family),
          legend.text=element_text(size=legendSize,family = family),
          legend.title=element_text(size=legendTitleSize,family = family))
  # print(gg[[j]])
}

ggarrange(plotlist = gg, nrow = 1)
#graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*1.5, height=height*.7, append=TRUE)

for(i in 1:length(prop)){ # becouse DKO come first in prop and WT come in first in conditions
  if(i%%2==1){
    if(prop[i+1]/cond[1]>prop[i]/cond[2]){
      abc=prop.test(c(prop[i+1],prop[i]),c(cond[1],cond[2]),alternative="greater")$p.value
    }else{
      abc=prop.test(c(prop[i+1],prop[i]),c(cond[1],cond[2]),alternative="less")$p.value
    }
    print(prop[c(i,i+1)])
    print(abc)
  }
}