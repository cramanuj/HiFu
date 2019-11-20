####################################################################
## Code to run multi-set CCA on mouse Organoid data (scRNA-seq)
##
## Author: Chaitanya Acharya
## Updated: Nov 7, 2019
####################################################################
## Set your working directory
setwd("~/Research/Takuya/") ## setwd("your directory")

## Install BiocManager if not yet installed
if(!requireNamespace("BiocManager")) install.packages("BiocManager")

## Load the required R libraries
## If the libraries are not installed, use BiocManager
lib.list=c("Seurat","plyr","ggpubr","dplyr","data.table","GSVA","biomaRt","monocle","genefu","survminer","survival","ggthemes","survutils","clustree","psych",
            "RColorBrewer","pheatmap")

for(i in 1:length(lib.list)){
	if(any(installed.packages()[,1]==lib.list[i])){
		library(lib.list[i],character.only=T)}else{
			BiocManager::install(lib.list[i])
			library(lib.list[i],character.only=T)
    };
}

## R scripts with utility functions
source("~/Research/scripts/r_scripts/useful_functions.R")
source("~/Research/scripts/r_scripts/plotfns.R")

## Cell cycle genes
cc_genes = scan("~/Research/Josh/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt",what="")

# We can segregate this list into markers of G2/M phase and markers of S phase and convert them from human to mouse markers
s_genes = convertHumanGeneList(cc_genes[1:43])
g2m_genes = convertHumanGeneList(cc_genes[44:97])

### Read all the possible .tsv files and identity select files
### Create a list object containing gene expression matrices of the select files
files = list.files(pattern=".tsv",path="scRNAseq_HiFU_TO",full.names=T)
data_list = llply(1:length(files),.progress="time",function(i){
  tmp = read.delim(files[i],header=T,row.names=1)
  # colnames(tmp)=paste(gsub(".tsv","",files[i]),colnames(tmp),sep="_")
  return(tmp)
})

### Create a SEURAT object from each gene expression matrix within the list
### Each SEURAT object is preprocessed separately and highly variable genes are computed
seu_list = gene_list = vector(mode="list",length=length(data_list))
for(i in 1:length(seu_list)){
  print(i)
  tmp = data_list[[i]]
  tmp = tmp[-c(grep("Rpl",rownames(tmp)),grep("Rps",rownames(tmp)),grep("^mt-",rownames(tmp))),]
  seu = CreateSeuratObject(raw.data=data_list[[i]],min.cells=5)
  seu = FilterCells(seu, subset.names = "nGene", low.thresholds = round(quantile(seu@meta.data$nGene,0.1)), high.thresholds = Inf)
  seu = NormalizeData(seu,display.progress = F)
  seu = CellCycleScoring(seu, s.genes = s_genes, g2m.genes = g2m_genes, set.ident = F)
  seu@meta.data$CC.Difference = seu@meta.data$S.Score - seu@meta.data$G2M.Score
  seu = ScaleData(seu, display.progress = F,vars.to.regress = c("nUMI"))
  seu = FindVariableGenes(seu,display.progress = F, do.plot = F, mean.function = ExpMean, dispersion.function = LogVMR)
  genes_use = head(rownames(seu@hvg.info), 2000)
  seu_list[[i]] = seu
  gene_list[[i]] = genes_use
}

#######################################################################################
## Gene selection for CCA
## We find genes that are highly variable in at least two datasets
#######################################################################################
genes_use <- names(which(table(unlist(gene_list)) > 1))
for (i in 1:length(seu_list)) {
  genes_use <- genes_use[genes_use %in% rownames(seu_list[[i]]@scale.data)]
}

#######################################################################################
### Run multi-set Canonical Correlation Analysis (CCA) on the common genes across all the genesets
### Calculate the ratio of total variance explained by PPCA vs total variance explained by CCA,
###       and filter cells based on these values
### CCA Align the matrices and compute the alignment metric score
### Map a t-SNE plot and find clusters
#######################################################################################

# pdf("CCA_plots.pdf",width=10,height=7)
cca_out = RunMultiCCA(seu_list,genes.use=genes_use,num.ccs = 20)
DimPlot(object = cca_out, reduction.use = "cca", group.by = "orig.ident", pt.size = 1, do.return = F)
VlnPlot(object = cca_out, features.plot = "CC1", group.by = "orig.ident", do.return = F)
MetageneBicorPlot(cca_out, grouping.var = "orig.ident", dims.eval = 1:20, display.progress = FALSE)
DimHeatmap(object = cca_out, reduction.type = "cca", cells.use = 500, dim.use = 1:9, do.balanced = TRUE)
cca_out = CalcVarExpRatio(cca_out,reduction.type = "pca", grouping.var = "orig.ident", dims.use = 1:20)
cca_out = SubsetData(cca_out, subset.name = "var.ratio.pca",accept.low = 0.40)
cca_out = AlignSubspace(cca_out, reduction.type = "cca", grouping.var = "orig.ident", dims.align = 1:10,num.genes=50)
metric = rep(0,10)
for(i in 1:10){
  metric[i]=CalcAlignmentMetric(cca_out,reduction.use = "cca.aligned",dims.use = 1:i, grouping.var =  "orig.ident")
  cat("Number of PCs: ",i," Alignment metric: ",metric[i],"\n")
}
cca_out@meta.data$orig.ident = factor(cca_out@meta.data$orig.ident,levels=c("UN","HF","MF","PD","PDL1"))

max_pc = 5
cca_out = FindClusters(cca_out, reduction.type = "cca.aligned", resolution = c(0.4,0.6,0.8,1), dims.use = 1:max_pc,print.output = F)
clustree(cca_out,prefix="res.",layout="sugiyama")
clustree(cca_out, prefix = "res.", node_colour = "sc3_stability")
cca_out = SetAllIdent(cca_out,id="res.0.8")
cca_out = BuildClusterTree(cca_out)
PlotClusterTree(cca_out)

cca_out@meta.data$res.0.8 = factor(cca_out@meta.data$res.0.8,levels=c(0:10))

cca_out = RunTSNE(cca_out, reduction.use = "cca.aligned", dims.use = 1:max_pc, do.fast = T,check_duplicates=F)
tsne_plot = TSNEPlot(object = cca_out, do.label =T,no.legend = FALSE,dark.theme = FALSE,label.size=5)
TSNEPlot(object = cca_out, do.label = F,group.by="orig.ident")

# cca_out = RunUMAP(cca_out,reduction.use = "cca.aligned",dims.use = 1:5,n_neighbors = 30L,min_dist=0.5,metric="cosine",umap.method = "umap-learn")
cca_out = RunUMAP(cca_out,reduction.use = "cca.aligned",dims.use = 1:4,n_neighbors = 10L,min_dist=0.1,metric="euclidean",umap.method = "umap-learn")
umap_plot = DimPlot(object = cca_out, reduction.use = 'umap',do.label=T,label.size=8) + theme_void() + theme(legend.position="none")
DimPlot(object = cca_out, reduction.use = 'umap',group.by="orig.ident") + theme_void() + theme(legend.title = element_blank())
DimPlot(object = cca_out, reduction.use = 'umap',group.by="orig.ident") + facet_wrap(~cca_out@meta.data$res.0.8)+ theme(legend.position="none")

# ggarrange(tsne_plot,umap_plot,nrow=1,ncol=2)
pdf("feature_plots.pdf",width=10,height=7)
DimPlot(object = cca_out, reduction.use = 'umap',do.label=T,label.size=5,no.legend=F)
## T cells
FeaturePlot(cca_out, no.legend=F,features.plot = c("Cd3d","Cd3g","Cd3e","Cd2"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
## CD8+ T cells
FeaturePlot(cca_out, no.legend=F,features.plot = c("Cd8a","Gzma"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
## CD4+ T cells and Treg
FeaturePlot(cca_out, no.legend=F,features.plot = c("Cd4","Foxp3"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
## NK cells
FeaturePlot(cca_out, no.legend=F,features.plot = c("Klrc1","Klrc3"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
## B cells
FeaturePlot(cca_out, no.legend=F,features.plot = c("Cd79a"), reduction.use="umap",nCol=1,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
## Plasma cells
FeaturePlot(cca_out, no.legend=F,features.plot = c("Slamf7","Igkc"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
## Dendritic cells
FeaturePlot(cca_out, no.legend=F,features.plot = c("Flt3"), reduction.use="umap",nCol=1,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
### Macrophages iNOS, CD80, MHCII, TLR2,4
FeaturePlot(cca_out, no.legend=F,features.plot = c("Apoe","H2-Aa","C1qa","C1qb","C1qb"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)

FeaturePlot(cca_out, no.legend=F,features.plot = c("Pdcd1","Ctla4","Havcr2","Cd274"), reduction.use="umap",nCol=2,min.cutoff = "q05", max.cutoff = "q95", cols.use = c("lightgrey", "blue"), pt.size = 1)
dev.off()

##########################################################
## Cluster-specific markers
##########################################################
markers = FindAllMarkers(object = cca_out,only.pos = FALSE, min.pct = 0.25, thresh.use = 0.25)
top10 = markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
top20 = markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
DoHeatmap(object = cca_out, genes.use = top10$gene, slim.col.label = TRUE, remove.key = F,cex.row = 5,group.cex = 7,col.low="blue",col.mid="white",col.high="red")

#################################################
## ssGSEA using genesets from Azizi et al
#################################################
azizi_gs = read.gmt.file("~/Research/pathways/mouse_pathways/azizi_genesets_mouse.gmt")
azizi_gsva = gsva(as.matrix(cca_out@data),azizi_gs$genesets,method="ssgsea",ssgsea.norm = F)
mmin = rowMins(azizi_gsva); mmax = rowMaxs(azizi_gsva)
azizi_out = azizi_gsva/(mmax - mmin)
rownames(azizi_out) = rownames(azizi_out)
azizi_out = data.frame(t(azizi_out))
# azizi_out = data.frame(t(Matrix::t(scale(Matrix::t(azizi_gsva)))))
cca_out@meta.data = cbind(cca_out@meta.data,t(azizi_out))
azizi_out$cluster = cca_out2@meta.data$res.0.8
azizi_out$Rx = as.factor(cca_out2@meta.data$Rx)
azizi_melt = melt(azizi_out)
names(azizi_melt)[3] = c("Genesets")
ES1_mean = azizi_melt[,-2] %>% dplyr::group_by(Genesets,cluster) %>% summarize_all(list(winsor.mean))
my_palette = rev(colorRampPalette(c("red", "white", "blue"))(n = 149))
col_breaks = c( seq(-1,0,length=50),               # for red
                seq(0.01,0.8,length=50),           # for yellow
                seq(0.81,1,length=50))
tmp = cca_out2@meta.data[,c("res.0.6","Rx")]
colnames(tmp)[1:2]=c("cluster","Rx")
tmp = tmp[order(tmp$Rx),]
tmp$cluster = as.factor(tmp$cluster)
tmp_dat = Matrix::t(scale(Matrix::t(azizi_gsva)))
tmp_dat = tmp_dat[,match(rownames(tmp),colnames(tmp_dat))]
mat_col = list(cluster=matlab.colors(nlevels(tmp$cluster)),Rx=genespring.colors(3))
names(mat_col[[2]])=levels(tmp$Rx)
names(mat_col[[1]])=levels(tmp$cluster)
pheatmap(tmp_dat,color=my_palette,breaks=col_breaks,cluster_row=T,cluster_col=F,annotation_col=tmp,annotation_colors =mat_col,show_colnames = F,cutree_rows=3,treeheight_row = 0,gaps_col=cumsum(table(tmp$Rx)))
