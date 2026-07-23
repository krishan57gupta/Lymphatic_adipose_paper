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
mainDir <- "yourDataAndCodeFolder"
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

###############################################################################
combined5=readRDS(paste(processedDataDir,'combined_',dataType,'_ABC_F.rds',sep=''))
Idents(combined5)<-combined5@meta.data$ABC_subCellTypes
############################### Extended data Figure 6,a
DP<-Seurat::DimPlot(combined5, reduction = "ABC_umap", pt.size=.1, label = TRUE, split.by = 'conditions')+
  theme(text = element_text(size=txtSize,family = family),
        plot.title = element_text(size = plotTitleSize, family = family),
        axis.text = element_text(size = axisSize,family = family),
        axis.title=element_text(size = axisTitleSize,family = family),
        legend.text=element_text(size=legendSize,family = family),
        legend.title=element_text(size=legendTitleSize,family = family))
print(DP)
# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*3, height=height*2, append=TRUE)

cellTypesMarkersList=list()
cellTypesMarkersList[["ABC_filtered11"]]=list('Capillary LECs'=c('Pdpn','Piezo2','Stmn2','Tppp3','Ano1','Lyve1'),
                                              'Pre-Collecting LECs'=c('Fabp4','Plvap','Egfl7', 'Cdh5','Pdgfb', 'Eng'),
                                              'Valve LECs'=c('Prox1','Nfatc1','Itga9'),
                                              'Dividing Cells'=c('Mki67',"Top2a","Bub1","Ccnb1",'Rb1'),
                                              'Enterocytes'=c('Anpep','Sult1b1','Gda','Btnl1','Ugt2b34'),
                                              'Goblet cells'=c('Muc2','Tff3','Agr2','Spink4'),
                                              'Paneth cells'=c('Reg4','Mptx2','Defa5','Defa26'),
                                              'Macrophages'=c('Cxcl2', 'Lgals3','C1qb','Cd68','Lyz2'),
                                              'B cells'=c('Iglv1','Mzb1','Iglc2','Cd79b'))
markersList=cellTypesMarkersList[["ABC_filtered11"]]
Idents(combined5)<-combined5@meta.data$ABC_subCellTypes
############################### Extended data Figure 6,b
HM1=DotPlot(combined5, features = markersList,
            cols = colorCode,
            dot.scale = 2,
            group.by = "ABC_subCellTypes") +RotatedAxis()+
  theme(text = element_text(size=txtSize,family = family),
        plot.title = element_text(size = plotTitleSize, family = family),
        axis.text = element_text(size = axisSize,family = family),
        axis.title=element_text(size = axisTitleSize,family = family),
        axis.text.x=element_text(angle=90),
        legend.text=element_text(size=legendSize,family = family),
        legend.title=element_text(size=legendTitleSize,family = family))
print(HM1)
# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*4, height=height*1.5, append=TRUE)


###############################################################
############################################################### only for LEC markers plot
###############################################################
Idents(combined5)<-combined5@meta.data$ABC_subCellTypes
combined6=subset(x = combined5, idents = c('Cap LEC1','Cap LEC2','Cap LEC3','Pre-Col LEC','Val LEC'))
cellTypesMarkersList[["ABC_filtered11"]]=list('Capillary LECs'=c('Pdpn','Piezo2','Stmn2','Tppp3','Ano1','Lyve1'),
                                              'Pre-Collecting LECs'=c('Fabp4','Plvap','Egfl7', 'Cdh5','Pdgfb', 'Eng'),
                                              'Valve LECs'=c('Prox1','Nfatc1','Itga9'))
markersList=cellTypesMarkersList[["ABC_filtered11"]]
############################### Figure 4,b
DP<-Seurat::DimPlot(combined6, reduction = "ABC_umap", pt.size=.1,
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

############################### Extended data Figure 6,c
HM1=DotPlot_1(combined6, features = markersList,
              cols =c("white","red"), # colorCode,
              group.by = "ABC_subCellTypes",
              dot.scale = 10) +RotatedAxis()+
  theme(text = element_text(size=txtSize,family = family),
        plot.title = element_text(size = plotTitleSize, family = family),
        axis.text = element_text(size = axisSize,family = family),
        axis.title=element_text(size = axisTitleSize,family = family),
        axis.text.x=element_text(angle=90),
        legend.text=element_text(size=legendSize,family = family),
        legend.title=element_text(size=legendTitleSize,family = family))
print(HM1)
# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*2, height=height*1, append=TRUE)
###############################################################
combined5@meta.data$ABC_Cond_SCT=factor(combined5@meta.data$ABC_Cond_SCT,levels=c('Cap LEC1_WT-HFD','Cap LEC1_LOF-HFD','Cap LEC2_WT-HFD','Cap LEC2_LOF-HFD',
                                                                                  'Cap LEC3_WT-HFD','Cap LEC3_LOF-HFD','Pre-Col LEC_WT-HFD','Pre-Col LEC_LOF-HFD',
                                                                                  'Val LEC_WT-HFD','Val LEC_LOF-HFD','DC_WT-HFD','DC_LOF-HFD',
                                                                                  'Ent_WT-HFD','Ent_LOF-HFD','GC_WT-HFD','GC_LOF-HFD',
                                                                                  'PC_WT-HFD','PC_LOF-HFD','MP_WT-HFD','MP_LOF-HFD',
                                                                                  'BC_WT-HFD','BC_LOF-HFD'))

capGenes1=c('Piezo2','Pdpn')
capGenes2=c('Ano1','Stmn2')
capGenes3=c('Stmn2','Krt18','Egfl7','Eng','Plvap','Fabp4','Tppp3','Pdgfb')
colGenes=c('Fabp4','Plvap')
valGenes=c('Sh2d6','Krt18')

Genes=list('Cap LEC1'=capGenes1,
           'Cap LEC2'=capGenes2,
           'Cap LEC3'=capGenes3,
           'Pre-Col LEC'=colGenes,
           'Val LEC'=valGenes)

Idents(combined5)<-combined5@meta.data$ABC_subCellTypes
point=.5
color1=colorCode
############################### Extended data Figure 6,c
VlnPlot(combined5, features = capGenes1,split.by = 'conditions',idents=c('Cap LEC1'),cols=color1,pt.size = point)# +labs(x = "")
# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*1, height=height*1, append=TRUE)
VlnPlot(combined5, features = capGenes2,split.by = 'conditions',idents=c('Cap LEC2'),cols=color1,pt.size = point)# +labs(x = "")
# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*1, height=height*1, append=TRUE)
VlnPlot(combined5, features = capGenes3,split.by = 'conditions',idents=c('Cap LEC3'),cols=color1,pt.size = point, ncol=8)# +labs(x = "")
# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*3, height=height*1, append=TRUE)
VlnPlot(combined5, features = colGenes,split.by = 'conditions',idents=c('Pre-Col LEC'),cols=color1,pt.size = point)# +labs(x = "")
# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*1, height=height*1, append=TRUE)
VlnPlot(combined5, features = valGenes,split.by = 'conditions',idents=c('Val LEC'),cols=color1,pt.size = point)# +labs(x = "")
# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*1, height=height*1, append=TRUE)

###############################################################################
combined5=readRDS(paste(processedDataDir,'combined_',dataType,'_ABC_F.rds',sep=''))
Idents(combined5)<-combined5@meta.data$ABC_subCellTypes #ABC_subCellTypes
cond=table(combined5@meta.data$conditions)
allCells=table(combined5@meta.data$conditions)
combined5=subset(x = combined5, idents = c('Cap LEC1','Cap LEC2','Cap LEC3','Pre-Col LEC','Val LEC'))
combined5@meta.data$ABC_Cond_CT<-paste(combined5@meta.data$ABC_cellTypes,combined5@meta.data$conditions,sep="_")



prop=table(combined5@meta.data$ABC_Cond_SCT) # ABC_Cond_SCT
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
############################### Extended data Figure 7,a
ggarrange(plotlist = gg, nrow = 1)
#graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*3, height=height*1, append=TRUE)

######################### proportion test
WT=allCells[1]
LOF=allCells[2]
# proportion test for each cell type
st=c()
for(i in unique(combined5@meta.data$ABC_subCellTypes)) # ABC_subCellTypes
{
  print(i)
  WTC=table(combined5@meta.data[combined5@meta.data$ABC_subCellTypes==i,]$conditions)[1]  # ABC_subCellTypes
  LOFC=table(combined5@meta.data[combined5@meta.data$ABC_subCellTypes==i,]$conditions)[2] # ABC_subCellTypes
  if((WTC/WT*100)>(LOFC/LOF*100)){
    alternative="greater"
  }else{
    alternative="less"
  }
  PT=prop.test(c(WTC,LOFC),c(WT,LOF),alternative=alternative)
  temp=c(i,PT$p.value,alternative,WTC,WT,LOFC,LOF,PT$estimate)
  st=rbind(st,temp)
}
colnames(st)=c("cellTypes","pValue",'alternative','WTC','WT','LOFC','LOF',"WTCprop","LOFprop")
st

####################### Trajectory analysis (Monocle)###################################
plots=list()
cds_list=list()
for(condN in c("WT-HFD","LOF-HFD")){
  combined5=readRDS(paste(processedDataDir,'combined_',dataType,'_ABC_F.rds',sep=''))
  Idents(combined5)=combined5@meta.data$ABC_subCellTypes
  combined6=subset(x = combined5, idents = c('Cap LEC1', 'Cap LEC2', 'Cap LEC3', 'Pre-Col LEC', 'Val LEC'))
  Idents(combined6)=combined6@meta.data$conditions
  combined6=subset(x = combined6, idents = c(condN))
  combined6@meta.data$ABC_subCellTypes=factor(x = combined6@meta.data$ABC_subCellTypes, levels =c('Cap LEC1', 'Cap LEC2', 'Cap LEC3', 'Pre-Col LEC', 'Val LEC'))
  DefaultAssay(combined6)<-"RNA"
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
  
  if(condN=="WT-HFD"){
    root_cell = colnames(cds)[cds@colData@listData$ABC_subCellTypes %in% c("Pre-Col LEC")]
  }
  if(condN=="LOF-HFD"){
    root_cell = colnames(cds)[cds@colData@listData$ABC_subCellTypes %in% c("Cap LEC1", "Cap LEC2", "Cap LEC3")]
  }
  
  root_cell=na.omit(root_cell)
  cds <- monocle3::order_cells(cds,reduction_method = "UMAP",root_cells=root_cell)
  
  cds_list[[condN]]=cds
  
  plot2=monocle3::plot_cells(cds, label_groups_by_cluster=FALSE, graph_label_size=0,label_leaves=F, color_cells_by = "pseudotime",
                             raster = TRUE, trajectory_graph_segment_size = .25, trajectory_graph_color = "red")
  plots[[condN]]=plot2
}
# saveRDS(cds_list, paste(processedDataDir,'combined_',dataType,'_cds_list.rds',sep=''))
############################### Extended data Figure 7,b
ggarrange(plotlist = plots ,ncol=2,nrow = 1, common.legend = TRUE)
# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*1.2, height=height*.7, append=TRUE)

####################################################################################
############################pathways analysis########################################
####################################################################################
combined5=readRDS(paste(processedDataDir,'combined_',dataType,'_ABC_F.rds',sep=''))
Idents(combined5)<-combined5@meta.data$ABC_Cond_SCT
nk.markers_separated=c()
for(i in unique(combined5@meta.data$ABC_subCellTypes)){
  print(i)
  abc1=FindMarkers(combined5, ident.1 = paste(i,"WT-HFD",sep="_"),ident.2=paste(i,"LOF-HFD",sep="_"),min.cells.group = 1)
  abc1['cellType']=i
  abc1['condition1']="WT-HFD"
  abc1['condition2']="LOF-HFD"
  abc1['genes']=rownames(abc1)
  nk.markers_separated=rbind(nk.markers_separated,abc1)
}
head(nk.markers_separated)
saveRDS(nk.markers_separated,paste(processedDataDir,'combined_',dataType,'_ABC_nk.markers_separated.rds',sep=''))
nk.markers_separated=readRDS(paste(processedDataDir,'combined_',dataType,'_ABC_nk.markers_separated.rds',sep=''))
# nk.markers_separated=nk.markers_separated[nk.markers_separated$p_val<0.05,]
nk.markers_separated=nk.markers_separated[order(abs(nk.markers_separated$avg_log2FC)),]
write.csv(nk.markers_separated,paste(plotDir,'combined_',dataType,'_cellSubPopulationConservedMarkersBetweenCondition.csv',sep=''))


library(clusterProfiler)
library(org.Mm.eg.db)
pathNum=10000 # 30
qvalue=1 # 0.05
ct =subCellTypeOrder
dataGO=c()
for(ii in ct){
  print(ii)
  
  for(pt in c("CC","BP","MF")){
    ################################################################################
    df=c()
    ct1='WT-HFD'
    ct2='LOF-HFD'
    nk.markers_0<-readRDS(paste(processedDataDir,'combined_',dataType,'_ABC_nk.markers_separated.rds',sep=''))
    nk.markers_0=nk.markers_0[nk.markers_0$condition1==ct1 & nk.markers_0$condition2==ct2,]
    nk.markers_0=nk.markers_0[nk.markers_0$avg_log2FC<0 & nk.markers_0$cellType==ii,]
    nk.markers_0=nk.markers_0[nk.markers_0$p_val<0.05,]
    if(dim(nk.markers_0)[1]>0){
      genes=nk.markers_0[,"genes"]
      gene.df <- bitr(genes, fromType = "SYMBOL",
                      toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                      OrgDb = org.Mm.eg.db)
      up <- enrichGO(gene         = gene.df$SYMBOL,
                     OrgDb         = org.Mm.eg.db,
                     pvalueCutoff = 1,
                     qvalueCutoff = 1,
                     minGSSize = 1,
                     maxGSSize = 10000,
                     keyType       = 'SYMBOL',
                     ont           = pt,
                     pAdjustMethod = "BH")
      # print(head(up,20))
      up=data.frame(up)
      up=up[up$qvalue<=qvalue,]
      print(dim(up))
      mm=paste("upregulated in ",ct2," as compared to ",ct1," ( ",ii," )",sep="")
      df=rbind(df,cbind("regu"="upregulated in LOF","regulation"=mm,up[1:min(pathNum,dim(up)[1]),]))
      df['pathwaysType']<- pt
      df['log2qvalue']<- ceil(-log2(df$qvalue))
      df['geneRatio']<-unlist(lapply(strsplit(df$GeneRatio,"/"),function(x) as.numeric(x[1])/as.numeric(x[2])))
      df['bgRatio']<-unlist(lapply(strsplit(df$BgRatio,"/"),function(x) as.numeric(x[1])/as.numeric(x[2])))
      df['foldChange']<- as.numeric(df$geneRatio)/as.numeric(df$bgRatio)
      df['log2FoldChange']<- log2(df['foldChange'])
      df=df[order(df$log2qvalue),]
      ggPlot=ggplot(data=df, aes(x=Description, y=log2FoldChange)) +
        geom_bar(stat="identity", position=position_dodge(1),fill="lightblue",color="white")+ theme_classic()+
        geom_text(aes(label=Description), position=position_dodge(width=1), size=3,vjust=0,hjust=1, color="Black",angle = 0)+
        geom_text(aes(label=Count), position=position_dodge(width=1), size=3,vjust=0,hjust=0, color="Black",angle = 0)+
        labs(y = "", x = "")+
        ggtitle(mm)+ coord_flip()+
        theme(text = element_text(size=txtSize,family = family),
              plot.title = element_text(size = plotTitleSize, family = family),
              axis.text = element_text(size = axisSize,family = family),
              # axis.text.x=element_blank(),
              # axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.title=element_text(size = axisTitleSize,family = family),
              legend.text=element_text(size=legendSize,family = family),
              legend.title=element_text(size=legendTitleSize,family = family))
      print(ggPlot)
      # graph2ppt(file=paste(plotDir,dataType,"_2",sep=''), width=width*2, height=.5+height*(dim(df)[1]*.15), append=TRUE)
      dataGO=rbind(dataGO,cbind(df,"cell type"=ii))
    }
    
    
    ##########
    df=c()
    ct1='WT-HFD'
    ct2='LOF-HFD'
    nk.markers_0<-readRDS(paste(processedDataDir,'combined_',dataType,'_ABC_nk.markers_separated.rds',sep=''))
    nk.markers_0=nk.markers_0[nk.markers_0$condition1==ct1 & nk.markers_0$condition2==ct2,]
    nk.markers_0=nk.markers_0[nk.markers_0$avg_log2FC>0 & nk.markers_0$cellType==ii,]
    nk.markers_0=nk.markers_0[nk.markers_0$p_val<0.05,]
    if(dim(nk.markers_0)[1]>0){
      genes=nk.markers_0[,"genes"]
      gene.df <- bitr(genes, fromType = "SYMBOL",
                      toType = c("ENSEMBL", "SYMBOL","ENTREZID"),
                      OrgDb = org.Mm.eg.db)
      down <- enrichGO(gene         = gene.df$SYMBOL,
                       OrgDb         = org.Mm.eg.db,
                       pvalueCutoff = 1,
                       qvalueCutoff = 1,
                       minGSSize = 1,
                       maxGSSize = 10000,
                       keyType       = 'SYMBOL',
                       ont           = pt,
                       pAdjustMethod = "BH")
      # print(head(down,20))
      down=data.frame(down)
      down=down[down$qvalue<=qvalue,]
      print(dim(down))
      mm=paste("downregulated in ",ct2," as compared to ",ct1," ( ",ii," )",sep="")
      df=rbind(df,cbind("regu"="downregulated in LOF","regulation"=mm,down[1:min(pathNum,dim(down)[1]),]))
      df['pathwaysType']<- pt
      df['log2qvalue']<- ceil(-log2(df$qvalue))
      df['geneRatio']<-unlist(lapply(strsplit(df$GeneRatio,"/"),function(x) as.numeric(x[1])/as.numeric(x[2])))
      df['bgRatio']<-unlist(lapply(strsplit(df$BgRatio,"/"),function(x) as.numeric(x[1])/as.numeric(x[2])))
      df['foldChange']<- as.numeric(df$geneRatio)/as.numeric(df$bgRatio)
      df['log2FoldChange']<- log(df['foldChange'])
      df=df[order(df$log2qvalue),]
      ggPlot=ggplot(data=df, aes(x=Description, y=log2FoldChange)) +
        geom_bar(stat="identity", position=position_dodge(1),fill="lightblue",color="white")+ theme_classic()+
        geom_text(aes(label=Description), position=position_dodge(width=1), size=3,vjust=0,hjust=1, color="Black",angle = 0)+
        geom_text(aes(label=Count), position=position_dodge(width=1), size=3,vjust=0,hjust=0, color="Black",angle = 0)+
        labs(y = "", x = "")+
        ggtitle(mm)+ coord_flip()+
        theme(text = element_text(size=txtSize,family = family),
              plot.title = element_text(size = plotTitleSize, family = family),
              axis.text = element_text(size = axisSize,family = family),
              # axis.text.x=element_blank(),
              # axis.ticks.x=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              axis.title=element_text(size = axisTitleSize,family = family),
              legend.text=element_text(size=legendSize,family = family),
              legend.title=element_text(size=legendTitleSize,family = family))
      print(ggPlot)
      # graph2ppt(file=paste(plotDir,dataType,"_2",sep=''), width=width*2, height=.5+height*(dim(df)[1]*.15), append=TRUE)
      dataGO=rbind(dataGO,cbind(df,"cell type"=ii))
    }
  }
}
dataGOTemp=dataGO
# write.csv(dataGO,paste(plotDir,dataType,"_pathwaysGO_separated.csv",sep=''))
dataGO<-read.csv(paste(plotDir,dataType,"_pathwaysGO_separated.csv",sep=''))
dataGO['Description1']=paste(dataGO$Description,dataGO$pathwaysType,sep="_")
dim(dataGO)

dataGO1<-dataGO
dim(dataGO1)
dataGO1=dataGO1[dataGO1$pvalue<0.05,]
dim(dataGO1)
dataGO1<-dataGO1[dataGO1$cell.type %in% c("Cap LEC1","Cap LEC2", "Cap LEC3","Pre-Col LEC", "Val LEC"),]
dim(dataGO1)

write.csv(dataGO1,paste(plotDir,dataType,"_pathwaysGO_separated_1.csv",sep=''))
dataGO1<-read.csv(paste(plotDir,dataType,"_pathwaysGO_separated_1.csv",sep=''))
# dataGO1=dataGO1[dataGO1$pathwaysType=="BP",]
dim(dataGO1)


dim(dataGO1)
dataGO1$log2FoldChange<-as.numeric(dataGO1$log2FoldChange)
dataGO2<-dataGO1[which(dataGO1$regu == "upregulated in LOF"),]
goIDS=c('tight junction_CC','adherens junction_CC','angiogenesis_BP','cell-cell junction assembly_BP','cholesterol efflux_BP',
        'reverse cholesterol transport_BP','cell-cell junction organization_BP','regulation of MAPK cascade_BP',
        'regulation of lipoprotein lipase activity_BP','triglyceride metabolic process_BP','fatty acid metabolic process_BP',
        'cholesterol metabolic process_BP')
goIDS=goIDS[length(goIDS):1]
print(length(goIDS))
goIDS=unique(goIDS)
print(length(goIDS))
goAll=dataGO2[dataGO2$Description1 %in% goIDS,]


length(goIDS)
goIDS=unique(goIDS)
length(goIDS)
dataGO2<-dataGO2[dataGO2$Description1 %in% goIDS,]
dataGO2$Description1<-factor(dataGO2$Description1, levels=goIDS)
dim(dataGO2)

############################### Extended data Figure 9,j
gg2=ggplot(dataGO2, aes(x=Description1, y=log2FoldChange, fill=cell.type))+
  geom_bar(stat='identity',position=position_dodge(.5),width = 0.1,aes(color=cell.type))+
  geom_point(stat='identity',position=position_dodge(.5),size=3,aes(color=cell.type))+
  theme_classic()+coord_flip()+
  ylab("GO enrichment ratio log2FoldChange")+ #NoLegend()+
  ggtitle("Upregulated")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))+
  theme(text = element_text(size=txtSize,family = family),
        plot.title = element_text(size = plotTitleSize, family = family),
        axis.text = element_text(size = axisSize,family = family),
        axis.title=element_text(size = axisTitleSize,family = family),
        legend.text=element_text(size=legendSize,family = family),
        legend.title=element_text(size=legendTitleSize,family = family))
gg2
# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*3, height=height*1.5, append=TRUE)





dataGO2<-dataGO1[which(dataGO1$regu == "downregulated in LOF"),]
goIDS=c('collagen-containing extracellular matrix_CC','banded collagen fibril_CC','fibrillar collagen trimer_CC',
        'response to lipid_BP','cell communication_BP','adipose tissue development_BP',
        'negative regulation of ERK1 and ERK2 cascade_BP','positive regulation of endocytosis_BP','collagen biosynthetic process_BP',
        'negative regulation of metabolic process_BP','intracellular transport_BP')
goIDS=goIDS[length(goIDS):1]

print(length(goIDS))
goIDS=unique(goIDS)
print(length(goIDS))

goAll=dataGO2[dataGO2$Description1 %in% goIDS,]


length(goIDS)
goIDS=unique(goIDS)
length(goIDS)
dataGO2<-dataGO2[dataGO2$Description1 %in% goIDS,]
dataGO2$Description1<-factor(dataGO2$Description1, levels=goIDS)
dim(dataGO2)

############################### Extended data Figure 9,k
gg2=ggplot(dataGO2, aes(x=Description1, y=log2FoldChange, fill=cell.type))+
  geom_bar(stat='identity',position=position_dodge(.5),width = 0.1,aes(color=cell.type))+
  geom_point(stat='identity',position=position_dodge(.5),size=3,aes(color=cell.type))+
  theme_classic()+coord_flip()+
  ylab("GO enrichment ratio log2FoldChange")+ #NoLegend()+
  ggtitle("Downregulated")+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))+
  theme(text = element_text(size=txtSize,family = family),
        plot.title = element_text(size = plotTitleSize, family = family),
        axis.text = element_text(size = axisSize,family = family),
        axis.title=element_text(size = axisTitleSize,family = family),
        legend.text=element_text(size=legendSize,family = family),
        legend.title=element_text(size=legendTitleSize,family = family))
gg2
# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*3, height=height*1.5, append=TRUE)


######################################################################################### Cell chat analysis
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
cellchatList=list()
cellchatPlot=list()
sampleNames=sampleNames[c(1,2)]
for(condName in sampleNames){
  print(condName)
  ########################### data preperation
  combined=readRDS(paste(processedDataDir,'combined_',dataType,'_ABC_F.rds',sep=''))
  Idents(combined)<-combined@meta.data$ABC_subCellTypes
  combined=subset(combined, idents = subCellTypeOrderSelected)
  combined@meta.data$ABC_subCellTypes=factor(combined@meta.data$ABC_subCellTypes,levels = c(subCellTypeOrderSelected)) #subCellTypeOrder for all)
  combined1=combined
  Idents(combined1)<-combined1@meta.data$conditions
  dim(combined1)
  data.input=combined1@assays$RNA@data
  dim(data.input)
  
  ##########################
  meta=combined1@meta.data
  meta$condition=meta$conditions
  meta$labels=meta$ABC_subCellTypes
  cell.use = rownames(meta)[meta$condition == condName]
  data.input = data.input[, cell.use]
  dim(data.input)
  meta = meta[cell.use, ]
  dim(meta)
  unique(meta$labels)
  
  ########################### cellchat object creation
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
  CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
  dplyr::glimpse(CellChatDB$interaction) # just show the data
  CellChatDB.use <- CellChatDB  # for all datasets
  cellchat@DB <- CellChatDB.use
  
  ########################### cellchat analysis
  cellchat <- subsetData(cellchat)
  future::plan("multisession", workers = 4)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- projectData(cellchat, PPI.mouse)
  
  cellchat <- computeCommunProb(cellchat,raw.use = FALSE)
  cellchat <- filterCommunication(cellchat, min.cells = 1)
  
  ########################### Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  ########################### Calculate the aggregated cell-cell communication network
  groupSize <- as.numeric(table(cellchat@idents))
  First=netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  Second=netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
  cellchatPlot[[paste(condName,"First",sep="_")]]=First
  cellchatPlot[[paste(condName,"Second",sep="_")]]=Second
  
  cellchatList[[condName]]=cellchat
}
ggarrange(cellchatPlot[[1]],cellchatPlot[[2]],cellchatPlot[[3]],cellchatPlot[[4]],labels = names(cellchatPlot))
# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*2, height=height*2, append=TRUE)
# saveRDS(cellchatList,paste(plotDir,"cellchatListAll.rds",sep=""))
cellchatPlot
########################## for comparing multiple datasets
library(CellChat)
library(patchwork)

cellchatList=readRDS(paste(plotDir,"cellchatListAll.rds",sep=""))

cellchat <- mergeCellChat(cellchatList, add.names = names(cellchatList))

cellchat1=cellchat

gg1_pathways=data.frame("interaction_name"=c('Lamc1_Dag1','Lamc1_Itga6_Itgb1','Sema3a_Nrp1_Plxna2','Kitl_Kit',
                                             'Hspg2_Dag1','Ncam1_Ncam1','Angptl2_Tlr4','Cdh1_Cdh1','Cdh2_Cdh2'))
gg1_pathways[,1]=toupper(gg1_pathways[,1])
size=8
############################### Extended data Figure 9,f
gg1 <- netVisual_bubble(cellchat1, sources.use =c(1:5) , targets.use = c(1:5),  comparison = c(1,2),max.dataset = 2, min.dataset = 1,
                        pairLR.use=gg1_pathways,
                        font.size=size,font.size.title=size,
                        color.text = colorCode[c(2,1)],color.grid = c("grey"),angle.x = 90,
                        title.name = paste("Increased signaling in ",sampleNames[2]," as compared to ",sampleNames[1],sep=""))
gg1

############################### Extended data Figure 9,h
ncg1=netVisual_chord_gene(cellchatList[[2]], sources.use =c(1:5) , targets.use = c(1:5),
                          pairLR.use=gg1_pathways,
                          scale=TRUE,link.target.prop=FALSE,annotationTrackHeight=0.01,
                          title.name = paste("Increased signaling in ",sampleNames[2]," as compared to ",sampleNames[1],sep=""))

# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*3, height=height*3, append=TRUE)
###########################################################
########################## for selected pathways and selected pairs of cell types and decreased signaling in LOF
############################################################
cellchat2=cellchat

gg2 <- netVisual_bubble(cellchat2, sources.use =c(1,2,3,4,5) , targets.use = c(1,2,3,4,5),  comparison = c(1,2),max.dataset = 1, min.dataset = 2, 
                        color.text = color1[c(1,2)],color.grid = c("grey"), remove.isolate = T)
gg2
unique(gg2$data$interaction_name)

gg2_pathways=data.frame("interaction_name"=c('TNXB_SDC4','FN1_SDC4',
                                             'COL6A3_SDC4','COL6A3_ITGA1_ITGB1','COL6A2_SDC4','COL6A2_ITGA1_ITGB1',
                                             'COL6A1_SDC4','COL6A1_ITGA1_ITGB1','COL4A1_SDC4','COL4A1_ITGA1_ITGB1',
                                             'COL1A2_SDC4','COL1A2_ITGA1_ITGB1','COL1A1_ITGA1_ITGB1'))
gg2_pathways[,1]=toupper(gg2_pathways[,1])
size=8
############################### Extended data Figure 9,g
gg2 <- netVisual_bubble(cellchat2, sources.use =c(1,2) , targets.use = c(2,3,4,5),  comparison = c(1,2),max.dataset = 1, min.dataset = 2, 
                        pairLR.use=gg2_pathways,
                        font.size=size,font.size.title=size,
                        color.text = colorCode[c(1,2)],color.grid = c("grey"),angle.x = 90,
                        title.name = paste("Decreased signaling in ",sampleNames[2]," as compared to ",sampleNames[1],sep=""))
gg2

############################### Extended data Figure 9,j
netVisual_chord_gene(cellchatList[[2]], sources.use =c(1:5) , targets.use = c(1:5),
                     pairLR.use=gg2_pathways,
                     scale=TRUE,link.target.prop=FALSE,annotationTrackHeight=0.01,
                     title.name = paste("Decreased signaling in ",sampleNames[2]," as compared to ",sampleNames[1],sep=""))
# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*3, height=height*3, append=TRUE)
############################################################

################################################################################# another analysis from mebocost reploting
ct=c('Cap LEC1','Cap LEC2','Cap LEC3',"Pre-Col LEC","Val LEC")
WT=data.frame("conditions"="WT","cellType"=c(rep(ct,each=2)),
              "comm"=c(rep(c("Sender","Recever"),length(ct))),count=c(24,4,27,13,0,0,24,12,11,29))
LOF=data.frame("conditions"="LOF","cellType"=c(rep(ct,each=2)),
               "comm"=c(rep(c("Sender","Recever"),length(ct))),count=c(8,4,18,2,3,2,19,2,8,2))

WT_LOF=rbind(WT,LOF)
WT_LOF$comm=factor(WT_LOF$comm,levels = c("Sender","Recever"))
WT_LOF$conditions=factor(WT_LOF$conditions,levels = c("WT","LOF"))
WT_LOF$cellType=factor(WT_LOF$cellType,levels = ct)

############################### Extended data Figure 9,a
ggplot(WT_LOF, aes(x=cellType, y=count,fill=conditions))+
  geom_bar(stat='identity',position=position_dodge(1))+
  geom_text(aes(label=count), position=position_dodge(.9), size=3,vjust=0,hjust=.5, color="black",angle=0)+
  scale_fill_manual(values = colorCode[c(1,2)])+
  facet_wrap(~comm, scales="free_y", nrow=2,strip.position="right")+
  theme_classic()+ # xlab("")+
  ylab("communication count")+ #NoLegend()+
  ylim(0,75)+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+# coord_flip()+
  theme(text = element_text(size=txtSize,family = family),
        plot.title = element_text(size = plotTitleSize, family = family),
        axis.text = element_text(size = axisSize,family = family),
        axis.title=element_text(size = axisTitleSize,family = family),
        legend.text=element_text(size=legendSize,family = family),
        legend.title=element_text(size=legendTitleSize,family = family))# +facet_grid(rows=vars(comm))


###################################################################################################
combined5=readRDS(paste(processedDataDir,'combined_',dataType,'_ABC_F.rds',sep=''))
genes=c(rev(c('Pdpn','Piezo2','Stmn2','Tppp3','Ano1','Fabp4','Plvap','Egfl7', 'Cdh5','Prox1','Nfatc1','Krt18','Sh2d6')),'Pdgfb', 'Eng','Lyve1')
plots1=list()
for(g in genes){
  print(g)
  combined10=combined5
  Idents(combined10)<-combined10@meta.data$conditions
  combined10=subset(combined10,idents=c("WT-HFD"))
  Idents(combined10)<-combined10@meta.data$ABC_subCellTypes
  cells=combined10@assays$RNA@data[g,]
  # View(cells)
  combinedA=subset(combined10,cells=names(cells)[cells<=0])
  combinedB=subset(combined10,cells=names(cells)[cells>0])
  combined10=merge(combinedA, y = c(combinedB), merge.dr=c('pca', 'umap','ABC_pca', 'ABC_umap'),
                   add.cell.ids = c("1", "2"), project = "reassign")
  # newCells=combined10@assays$RNA@data[g,]
  # View(newCells)
  combined10@meta.data$conditions=factor(combined10@meta.data$conditions,levels = c("WT-HFD","LOF-HFD"))
  plots1[[paste(g,1,sep="_")]]=FeaturePlot(combined10,g,cols=c("lightgray","red"),
                                           # split.by = 'conditions'
                                           # raster = TRUE, raster.dpi = c(256, 256),
                                           pt.size=2)+ labs(title="", x="",y="")+
    theme(text = element_text(size=txtSize,family = family),
          plot.title = element_text(size = plotTitleSize, family = family),
          axis.text = element_text(size = axisSize,family = family),
          axis.title=element_text(size = axisTitleSize,family = family),
          legend.text=element_text(size=legendSize,family = family),
          legend.title=element_text(size=legendTitleSize,family = family))
  
  print(g)
  combined10=combined5
  Idents(combined10)<-combined10@meta.data$conditions
  combined10=subset(combined10,idents=c("LOF-HFD"))
  Idents(combined10)<-combined10@meta.data$ABC_subCellTypes
  cells=combined10@assays$RNA@data[g,]
  # View(cells)
  combinedA=subset(combined10,cells=names(cells)[cells<=0])
  combinedB=subset(combined10,cells=names(cells)[cells>0])
  combined10=merge(combinedA, y = c(combinedB), merge.dr=c('pca', 'umap','ABC_pca', 'ABC_umap'),
                   add.cell.ids = c("1", "2"), project = "reassign")
  combined10@meta.data$conditions=factor(combined10@meta.data$conditions,levels = c("WT-HFD","LOF-HFD"))
  plots1[[paste(g,2,sep="_")]]=FeaturePlot(combined10,g,cols=c("lightgray","red"),
                                           pt.size=2)+ labs(title="", x="",y="")+
    theme(text = element_text(size=txtSize,family = family),
          plot.title = element_text(size = plotTitleSize, family = family),
          axis.text = element_text(size = axisSize,family = family),
          axis.title=element_text(size = axisTitleSize,family = family),
          legend.text=element_text(size=legendSize,family = family),
          legend.title=element_text(size=legendTitleSize,family = family))
}
############################### Extended data Figure 6,d
ggarrange(plots1[[1]],plots1[[2]],plots1[[3]],plots1[[4]],plots1[[5]],plots1[[6]],plots1[[7]],plots1[[8]],
          plots1[[9]],plots1[[10]],plots1[[11]],plots1[[12]],plots1[[13]],plots1[[14]],plots1[[15]],plots1[[16]],
          plots1[[17]],plots1[[18]],plots1[[19]],plots1[[20]],plots1[[21]],plots1[[22]],plots1[[23]],plots1[[24]],
          plots1[[25]],plots1[[26]],plots1[[27]],plots1[[28]],plots1[[29]],plots1[[30]],plots1[[31]],plots1[[32]],
          ncol = 4,nrow = 8, common.legend=TRUE)
# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*3, height=height*8, append=TRUE)

################################################################################
####################### MEBOCOST data updatetion
mebocost=read.csv("yourDataAndCodeFolder/MEBOCOST/ABC_WT_LOF/communicationPvalue_1.csv")
mebocost[,"X"]<-NULL
ct=lapply(mebocost[,'Cell_Pair'],function(x) strsplit(x,'_')[[1]][1])
cc=lapply(mebocost[,'Cell_Pair'],function(x) strsplit(x,'_')[[1]][2])
ct1=lapply(mebocost[,'Cell_Pair'],function(x) strsplit(strsplit(x,'_')[[1]][1],"→")[[1]][1])
ct2=lapply(mebocost[,'Cell_Pair'],function(x) strsplit(strsplit(x,'_')[[1]][1],"→")[[1]][2])
mebocost[,"cellPair"]<-as.vector(unlist(ct))
mebocost[,"cellCondition"]<-as.vector(unlist(cc))
mebocost[,"sender"]<-as.vector(unlist(ct1))
mebocost[,"receiver"]<-as.vector(unlist(ct2))
mebocost[,'signalCellPair']<-paste(mebocost[,"Signal_Pair"],mebocost[,"cellPair"],sep="_")

meboRowNames=sort(unique(mebocost[,'signalCellPair']))
meboColNames=c("metabolite","sensor","sender","receiver","sender_metabolite_sensor_receiver",
               "metanoliteProportion (WT)","metanoliteProportion (LOF)","metanoliteProportion (Log2FC)",
               "sensorProportion (WT)","sensorProportion (LOF)","sensorProportion (Log2FC)",
               "CommScore (WT)","CommScore (LOF)","CommScore (Log2FC)",
               "permutation_test_fdr (WT)","permutation_test_fdr (LOF)")
mebo=matrix(0,nrow=length(meboRowNames),ncol=length(meboColNames))
rownames(mebo)<-meboRowNames
colnames(mebo)<-meboColNames


for(i in mebocost[,"signalCellPair"]){
  WT=mebocost[mebocost[,'signalCellPair']==i & mebocost[,'cellCondition']=="WT",]
  LOF=mebocost[mebocost[,'signalCellPair']==i & mebocost[,'cellCondition']=="LOF",]
  if(length(WT[,'Metabolite_Name'])==0){
    mebo[i,"metabolite"]=LOF[,'Metabolite_Name']
  }else{
    mebo[i,"metabolite"]=WT[,'Metabolite_Name']
  }
  
  if(length(WT[,'Sensor'])==0){
    mebo[i,"sensor"]=LOF[,'Sensor']
  }else{
    mebo[i,"sensor"]=WT[,'Sensor']
  }
  
  if(length(WT[,"sender"])==0){
    mebo[i,"sender"]=LOF[,"sender"]
  }else{
    mebo[i,"sender"]=WT[,"sender"]
  }
  
  if(length(WT[,"receiver"])==0){
    mebo[i,"receiver"]=LOF[,"receiver"]
  }else{
    mebo[i,"receiver"]=WT[,"receiver"]
  }
  
  if(length(WT[,"signalCellPair"])==0){
    mebo[i,"sender_metabolite_sensor_receiver"]=paste('(',LOF[,"sender"],")_(",LOF[,'Metabolite_Name'],")_(",LOF[,'Sensor'],")_(",LOF[,"receiver"],")",sep="")
  }else{
    mebo[i,"sender_metabolite_sensor_receiver"]=paste("(",WT[,"sender"],")_(",WT[,'Metabolite_Name'],")_(",WT[,'Sensor'],")_(",WT[,"receiver"],")",sep="")
  }
  
  if(length(WT[,"metabolite_prop_in_sender"])==0){
    mebo[i,"metanoliteProportion (WT)"]=0
  }else{
    mebo[i,"metanoliteProportion (WT)"]=WT[,"metabolite_prop_in_sender"]
  }
  
  if(length(LOF[,"metabolite_prop_in_sender"])==0){
    mebo[i,"metanoliteProportion (LOF)"]=0
  }else{
    mebo[i,"metanoliteProportion (LOF)"]=LOF[,"metabolite_prop_in_sender"]
  }
  
  mebo[i,"metanoliteProportion (Log2FC)"]=log2(as.numeric(mebo[i,"metanoliteProportion (LOF)"])/as.numeric(mebo[i,"metanoliteProportion (WT)"]))
  
  if(length(WT[,"sensor_prop_in_receiver"])==0){
    mebo[i,"sensorProportion (WT)"]=0
  }else{
    mebo[i,"sensorProportion (WT)"]=WT[,"sensor_prop_in_receiver"]
  }
  
  if(length(LOF[,"sensor_prop_in_receiver"])==0){
    mebo[i,"sensorProportion (LOF)"]=0
  }else{
    mebo[i,"sensorProportion (LOF)"]=LOF[,"sensor_prop_in_receiver"]
  }
  
  mebo[i,"sensorProportion (Log2FC)"]=log2(as.numeric(mebo[i,"sensorProportion (LOF)"])/as.numeric(mebo[i,"sensorProportion (WT)"]))
  
  if(length(WT[,"Commu_Score"])==0){
    mebo[i,"CommScore (WT)"]=0
  }else{
    mebo[i,"CommScore (WT)"]=WT[,"Commu_Score"]
  }
  
  if(length(LOF[,"Commu_Score"])==0){
    mebo[i,"CommScore (LOF)"]=0
  }else{
    mebo[i,"CommScore (LOF)"]=LOF[,"Commu_Score"]
  }
  
  mebo[i,"CommScore (Log2FC)"]=log2(as.numeric(mebo[i,"CommScore (LOF)"])/as.numeric(mebo[i,"CommScore (WT)"]))
  
  if(length(WT[,"permutation_test_fdr"])==0){
    mebo[i,"permutation_test_fdr (WT)"]=0
  }else{
    mebo[i,"permutation_test_fdr (WT)"]=WT[,"permutation_test_fdr"]
  }
  
  if(length(LOF[,"permutation_test_fdr"])==0){
    mebo[i,"permutation_test_fdr (LOF)"]=0
  }else{
    mebo[i,"permutation_test_fdr (LOF)"]=LOF[,"permutation_test_fdr"]
  }
  
}
rownames(mebo)<-NULL
# write.csv(mebo,"/yourDataAndCodeFolder/MEBOCOST/ABC_WT_LOF/mebocostAllCommunicationFC2.csv")

################################################################################
metaCom=read.csv("/yourDataAndCodeFolder/MEBOCOST/ABC_WT_LOF/mebocostAllCommunicationFC2.csv")
dim(metaCom)
metaCom=metaCom[metaCom$CommScore..Log2FC.<0,]
dim(metaCom)
metaName=c("2-Hydroxyestradiol", "25-Hydroxycholesterol", "27-Hydroxycholesterol", "Docosahexaenoic acid", "Eicosapentaenoic acid", 
           "Sphingosine 1-phosphate","Stearic acid", "Glycerol","Cholesterol", "25-Hydroxycholesterol",
           "triglyceride", "phospholipids", "apolipoproteins", "lipid", "Chylomicrons", "Very low density lipoprotein VLDL", "LDL", "HDL", "IDL", 
           "Fatty acid", "cholesterol esters", "steriod", "triacylglycerols","saturated fat", "unsaturated fat", "sterol", "Palmitic acid", 
           "Margaric acid", "Stearic acid", "Nonadecylic acid", "docosahexaenoic acid", "alpha-linolenic acid", "eicosapentaenoic acid",
           "palmitoleic acid", "myristic acid", "palmitic acid", "linoleic acid")
cmm=intersect(unique(metaCom$metabolite),unique(metaName))
metaCom=metaCom[metaCom$metabolite %in% cmm,]
dim(metaCom)
metaCom['MS']<-paste(metaCom$metabolite,metaCom$sensor,sep="_")
metaCom['SR']<-paste(metaCom$sender,metaCom$receiver,sep="_")

LOF=data.frame("conditions"="LOF","comScore"=metaCom$CommScore..LOF.,
               "abundance"=metaCom$metanoliteProportion..WT.,"expression"=metaCom$sensorProportion..WT.,
               metaCom[,c("sender","receiver","metabolite","sensor")])
WT=data.frame("conditions"="WT","comScore"=metaCom$CommScore..WT.,
              "abundance"=metaCom$metanoliteProportion..LOF.,"expression"=metaCom$sensorProportion..LOF.,
              metaCom[,c("sender","receiver","metabolite","sensor")])
WT_LOF=rbind(WT,LOF)
WT_LOF['MS']<-paste(WT_LOF$metabolite,WT_LOF$sensor,sep="_")
WT_LOF['SR']<-paste(WT_LOF$sender,WT_LOF$receiver,sep="_")
WT_LOF['SMSR']<-paste(WT_LOF$sender,WT_LOF$metabolite,WT_LOF$sensor,WT_LOF$receiver,sep="_")
WT_LOF['CS']<-paste(WT_LOF$conditions,WT_LOF$sender,sep="_")
WT_LOF$conditions=factor(WT_LOF$conditions,levels = c("WT","LOF"))
WT_LOF['SM']<-paste(WT_LOF$sender,WT_LOF$metabolite,sep="_")
WT_LOF['RS']<-paste(WT_LOF$sensor,WT_LOF$receiver,sep="_")


interested_LEC=c("Cap LEC1","Cap LEC2","Cap LEC3","Pre-Col LEC","Val LEC")
# interested_NLEC=c("DC","Ent","GC","PC","MP","BC","TC","Neu")
interested_NLEC=c("Ent","GC")
interested_CT=c(interested_LEC,interested_NLEC)
WT_LOF=WT_LOF[WT_LOF$sender %in% interested_CT,]
WT_LOF=WT_LOF[WT_LOF$receiver %in% interested_CT,]
WT_LOF$sender=factor(WT_LOF$sender,levels=interested_CT)
WT_LOF$receiver=factor(WT_LOF$receiver,levels=interested_CT)

selected=c("Cap LEC1_Cholesterol_Rora_Cap LEC1","Cap LEC2_Cholesterol_Rora_Cap LEC1","Pre-Col LEC_Cholesterol_Rora_Cap LEC1","Val LEC_Cholesterol_Rora_Cap LEC1",
           "Cap LEC1_Cholesterol_Rorb_Cap LEC1","Cap LEC2_Cholesterol_Rorb_Cap LEC1","Pre-Col LEC_Cholesterol_Rorb_Cap LEC1","Val LEC_Cholesterol_Rorb_Cap LEC1",
           "Cap LEC1_Cholesterol_Ldlr_Cap LEC1","Cap LEC2_Cholesterol_Ldlr_Cap LEC1","Pre-Col LEC_Cholesterol_Ldlr_Cap LEC1","Val LEC_Cholesterol_Ldlr_Cap LEC1",
           "Cap LEC1_Cholesterol_Cd36_Cap LEC1","Cap LEC2_Cholesterol_Cd36_Cap LEC1","Pre-Col LEC_Cholesterol_Cd36_Cap LEC1","Val LEC_Cholesterol_Cd36_Cap LEC1",
           
           "Cap LEC1_Cholesterol_Rora_Cap LEC2","Cap LEC2_Cholesterol_Rora_Cap LEC2","Pre-Col LEC_Cholesterol_Rora_Cap LEC2","Val LEC_Cholesterol_Rora_Cap LEC2",
           "Cap LEC1_Cholesterol_Rorc_Cap LEC2","Cap LEC2_Cholesterol_Rorc_Cap LEC2","Pre-Col LEC_Cholesterol_Rorc_Cap LEC2","Val LEC_Cholesterol_Rorc_Cap LEC2",
           
           "Cap LEC1_Cholesterol_Scarb1_Pre-Col LEC","Cap LEC2_Cholesterol_Scarb1_Pre-Col LEC","Pre-Col LEC_Cholesterol_Scarb1_Pre-Col LEC","Val LEC_Cholesterol_Scarb1_Pre-Col LEC",
           "Cap LEC1_Cholesterol_Rora_Pre-Col LEC","Cap LEC2_Cholesterol_Rora_Pre-Col LEC","Pre-Col LEC_Cholesterol_Rora_Pre-Col LEC","Val LEC_Cholesterol_Rora_Pre-Col LEC",
           "Cap LEC1_Cholesterol_Ldlr_Pre-Col LEC","Cap LEC2_Cholesterol_Ldlr_Pre-Col LEC","Pre-Col LEC_Cholesterol_Ldlr_Pre-Col LEC","Val LEC_Cholesterol_Ldlr_Pre-Col LEC",
           "Cap LEC1_Cholesterol_Cd36_Pre-Col LEC","Cap LEC2_Cholesterol_Cd36_Pre-Col LEC","Pre-Col LEC_Cholesterol_Cd36_Pre-Col LEC","Val LEC_Cholesterol_Cd36_Pre-Col LEC",
           
           "Cap LEC1_Cholesterol_Ldlr_Val LEC","Cap LEC2_Cholesterol_Ldlr_Val LEC","Pre-Col LEC_Cholesterol_Ldlr_Val LEC","Val LEC_Cholesterol_Ldlr_Val LEC",
           "Cap LEC1_Cholesterol_Cd36_Val LEC","Cap LEC2_Cholesterol_Cd36_Val LEC","Pre-Col LEC_Cholesterol_Cd36_Val LEC","Val LEC_Cholesterol_Cd36_Val LEC",
           "Cap LEC1_Cholesterol_Npc1l1_Val LEC","Cap LEC2_Cholesterol_Npc1l1_Val LEC","Pre-Col LEC_Cholesterol_Npc1l1_Val LEC","Val LEC_Cholesterol_Npc1l1_Val LEC",
           "Cap LEC1_Cholesterol_Abcg3_Val LEC","Cap LEC2_Cholesterol_Abcg3_Val LEC","Pre-Col LEC_Cholesterol_Abcg3_Val LEC","Val LEC_Cholesterol_Abcg3_Val LEC",
           
           "Cap LEC1_Cholesterol_Abcg3_Ent","Cap LEC2_Cholesterol_Abcg3_Ent","Pre-Col LEC_Cholesterol_Abcg3_Ent","Val LEC_Cholesterol_Abcg3_Ent",
           "Cap LEC1_Cholesterol_Npc1l1_Ent","Cap LEC2_Cholesterol_Npc1l1_Ent","Pre-Col LEC_Cholesterol_Npc1l1_Ent","Val LEC_Cholesterol_Npc1l1_Ent",
           "Cap LEC1_Cholesterol_Ldlr_Ent","Cap LEC2_Cholesterol_Ldlr_Ent","Pre-Col LEC_Cholesterol_Ldlr_Ent","Val LEC_Cholesterol_Ldlr_Ent",
           
           "Cap LEC1_Cholesterol_Abcg3_GC","Cap LEC2_Cholesterol_Abcg3_GC","Pre-Col LEC_Cholesterol_Abcg3_GC","Val LEC_Cholesterol_Abcg3_GC",
           "Cap LEC1_Cholesterol_Npc1l1_GC","Cap LEC2_Cholesterol_Npc1l1_GC","Pre-Col LEC_Cholesterol_Npc1l1_GC","Val LEC_Cholesterol_Npc1l1_GC",
           "Cap LEC1_Cholesterol_Ldlr_GC","Cap LEC2_Cholesterol_Ldlr_GC","Pre-Col LEC_Cholesterol_Ldlr_GC","Val LEC_Cholesterol_Ldlr_GC",
           "Cap LEC1_Cholesterol_Rorc_GC","Cap LEC2_Cholesterol_Rorc_GC","Pre-Col LEC_Cholesterol_Rorc_GC","Val LEC_Cholesterol_Rorc_GC",
           "Cap LEC1_Cholesterol_Scarb1_GC","Cap LEC2_Cholesterol_Scarb1_GC","Pre-Col LEC_Cholesterol_Scarb1_GC","Val LEC_Cholesterol_Scarb1_GC")
WT_LOF=WT_LOF[WT_LOF$SMSR %in% selected,]
# WT_LOF$SMSR=factor(WT_LOF$SMSR,levels=selected)
WT_LOF$CS=factor(WT_LOF$CS,levels=unique(WT_LOF$CS)[c(1,4,2,5,3,6)])

# WT_LOF$normScore=WT_LOF$comScore
# for (i in WT_LOF$SMSR){
#   a=WT_LOF$normScore[WT_LOF$SMSR==i]
#   print(length(a))
#   if(length(a)!=0){
#     a=a-min(a)
#     if (max(a)!=0){
#       a=a/(max(a)-min(a))
#     }
#   }
#   WT_LOF$normScore[WT_LOF$SMSR==i]=a
#   print(paste(a[1],a[2],b[1],b[2],length(a)))
# }

WT_LOF$normScore=WT_LOF$comScore
for (i in WT_LOF$SMSR){
  a=WT_LOF$normScore[WT_LOF$SMSR==i]
  if (max(a)==0){
    b=a
  }else{
    b=a/max(a)
  }
  WT_LOF$normScore[WT_LOF$SMSR==i]=b
  print(paste(a[1],a[2],b[1],b[2],length(a)))
}
WT_LOF$normScore=as.numeric(WT_LOF$normScore)

# ggplot(WT_LOF, aes(x=MS, y=normScore,fill=conditions))+
#   geom_bar(stat='identity',position=position_dodge())+
#   scale_fill_manual(values = colorCode[c(1,2)])+
#   # facet_wrap(~MS, scales="free_y", nrow=1)+
#   theme_classic()+ # xlab("")+
#   ylab("communication score")+ #NoLegend()+
#   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+coord_flip()+
#   theme(text = element_text(size=txtSize,family = family),
#         plot.title = element_text(size = plotTitleSize, family = family),
#         axis.text = element_text(size = axisSize,family = family),
#         axis.title=element_text(size = axisTitleSize,family = family),
#         legend.text=element_text(size=legendSize,family = family),
#         legend.title=element_text(size=legendTitleSize,family = family))+facet_grid(rows=vars(receiver),cols=vars(sender))
color=c("darkred","red","orange","darkgreen","chartreuse3","aquamarine")
colorList=list()
for(i in 1:6){
  colorList[[unique(WT_LOF$CS)[c(1,4,2,5,3,6)][i]]]=color[i]
}
names(colorList)=unique(WT_LOF$CS)
WT_LOF$receiver=factor(WT_LOF$receiver, levels = c("Cap LEC1","Cap LEC2","Pre-Col LEC","Val LEC","Ent","GC"))

color=colorCode
colorList1=list()
for(i in 1:2){
  colorList1[[unique(WT_LOF$conditions)[i]]]=color[i]
}
WT_LOF$conditions=factor(WT_LOF$conditions,levels = c("LOF", "WT"))
ggplot(WT_LOF, aes(x=RS, y=normScore,fill=conditions))+
  # geom_segment(stat='identity',position=position_dodge(.5),aes(xend=Description1,yend=0,color=CS))+
  geom_bar(stat='identity',position=position_dodge(.5),width = 0.1,aes(color=conditions))+
  geom_point(stat='identity',position=position_dodge(.5),size=3,aes(color=conditions))+
  scale_fill_manual(values = colorList1)+
  scale_colour_manual(values = colorList1)+
  # facet_wrap(~MS, scales="free_y", nrow=1)+
  theme_classic()+ # xlab("")+
  ylab("communication score")+ #NoLegend()+
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+ coord_flip()+
  theme(text = element_text(size=txtSize,family = family),
        plot.title = element_text(size = plotTitleSize, family = family),
        axis.text = element_text(size = axisSize,family = family),
        axis.title=element_text(size = axisTitleSize,family = family),
        legend.text=element_text(size=legendSize,family = family),
        legend.title=element_text(size=legendTitleSize,family = family))+facet_grid(cols=vars(SM))
# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*3, height=height*2, append=TRUE)

################################################################################
####################### RNA velocity rank of transition probabilites
# new version cell rank result saved in the file Velo_Rank_Mesentery_Combined_HFD_LOF.ipynb
################################################
# for WT-HFD
mat <- matrix(
  c(
    0.339, 0.358, 0.000, 0.302, 0.000,
    0.234, 0.261, 0.000, 0.505, 0.000,
    0.362, 0.432, 0.000, 0.206, 0.000,
    0.163, 0.774, 0.000, 0.063, 0.000,
    0.604, 0.295, 0.000, 0.102, 0.000
  ),
  nrow = 5,
  byrow = TRUE
)
rownames(mat) <- c("Val LEC", "Pre-Col LEC", "Cap LEC3", "Cap LEC2", "Cap LEC1")
colnames(mat) <- c("Cap LEC1", "Cap LEC2", "Cap LEC3", "Pre-Col LEC", "Val LEC")
mat=mat[1:5,c(1,2,4)]
mat
pheatmap::pheatmap(mat,display_numbers=mat,cluster_rows = FALSE,cluster_cols = FALSE,breaks = seq(0, 1, length.out = 101))
# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*1, height=height*1, append=TRUE)

# for LOF-HFD
mat <- matrix(
  c(
    0.000, 0.000, 0.031, 0.025, 0.944,
    0.000, 0.000, 0.273, 0.724, 0.003,
    0.000, 0.000, 0.976, 0.024, 0.000,
    0.000, 0.000, 0.412, 0.583, 0.005,
    0.000, 0.000, 0.367, 0.628, 0.004
  ),
  nrow = 5,
  byrow = TRUE
)
rownames(mat) <- c("Val LEC", "Pre-Col LEC", "Cap LEC3", "Cap LEC2", "Cap LEC1")
colnames(mat) <- c("Cap LEC1", "Cap LEC2", "Cap LEC3", "Pre-Col LEC", "Val LEC")
mat=mat[1:5,c(3:5)]
pheatmap::pheatmap(mat,display_numbers=mat,cluster_rows = FALSE,cluster_cols = FALSE, breaks = seq(0, 1, length.out = 101))
# graph2ppt(file=paste(plotDir,dataType,sep=''), width=width*1, height=height*1, append=TRUE)

