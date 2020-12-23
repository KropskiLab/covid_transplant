library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
options(future.globals.maxSize = 40000 * 1024^2)
library(RColorBrewer)
library(viridis)
library(hrbrthemes)

#read in covid sample

covid.data <- Read10X_h5("~/Desktop/COVID/filtered_feature_bc_matrix.h5")
adata <- CreateSeuratObject(counts = covid.data, min.features = 500, project = "COVID")

Idents(adata) <- 'orig.ident'
Idents(adata, cells = WhichCells(adata, idents = c('COVID'))) <- "COVID"
adata$Status <- Idents(adata)

adata <- PercentageFeatureSet(adata, pattern = "^MT-", col.name = "percent.mt")
adata <- subset(adata, subset = percent.mt < 15 & percent.mt >0.5 & nFeature_RNA >500)
adata <- SCTransform(adata, vars.to.regress = 'percent.mt')
adata <- RunPCA(adata, verbose = F)
adata <- RunUMAP(adata, dims = 1:25)
adata <- FindNeighbors(adata, dims = 1:25)
adata <- FindClusters(adata, resolution = 1.0)
DimPlot(adata, label = T)

marker_genes = c("PTPRC", "LYZ", "CD68", "CD14", "CD86", "CD3E", "CD4", "CD8A", "FOXP3", "GNLY", "CPA3", "CD19", "JCHAIN", "IRF7", "EPCAM", "FOXJ1", "SCGB1A1", "MUC5B", "SCGB3A2", "KRT5", "KRT17", "TP63", "SFTPC", "AGER", "PECAM1", "CA4", "CCL21", "WT1", "WNT2", "COL1A1", "HAS1", "DCN", "ACTA2", "CSPG4", "MKI67")
DotPlot(adata, features = marker_genes)+ theme(axis.text.x = element_text(angle = 45, hjust=1))

Idents(adata) <- 'seurat_clusters'
Idents(adata, cells = WhichCells(adata, idents = c(1,14,17,19,23))) <- "Epithelial"
Idents(adata, cells = WhichCells(adata, idents = c(0,22,4,12,15))) <- "Stromal"
Idents(adata, cells = WhichCells(adata, idents = c(2,3,5:11,13,16,18,20,21,24))) <- "Immune"
adata$population <- Idents(adata)

#Doublet exclusion epithelial cells
covid_epi <- subset(adata, idents = c('Epithelial'))
covid_epi <- SCTransform(covid_epi, vars.to.regress = 'percent.mt')
covid_epi <- RunPCA(covid_epi, verbose = F)
covid_epi <- RunUMAP(covid_epi, dims = 1:18)
covid_epi <- FindNeighbors(covid_epi, dims = 1:18)
covid_epi <- FindClusters(covid_epi, resolution = 0.6)
DimPlot(covid_epi, label = T)
DimPlot(covid_epi, group.by = 'orig.ident')
DotPlot(covid_epi, features = marker_genes)+ theme(axis.text.x = element_text(angle = 45, hjust=1))


#remove doublets
covid_epi <- subset(covid_epi, idents = c(0,1,3,6,7))
covid_epi <- SCTransform(covid_epi, vars.to.regress = 'percent.mt')
covid_epi <- RunPCA(covid_epi, verbose = F)
covid_epi <- RunUMAP(covid_epi, dims = 1:18)
covid_epi <- FindNeighbors(covid_epi, dims = 1:18)
covid_epi <- FindClusters(covid_epi, resolution = 1.0)
DimPlot(covid_epi, label = T)
DimPlot(covid_epi, group.by = 'orig.ident')
DotPlot(covid_epi, features = marker_genes)+ theme(axis.text.x = element_text(angle = 45, hjust=1))


#high-level-annotation
Idents(covid_epi) <- 'seurat_clusters'
Idents(covid_epi, cells = WhichCells(covid_epi, idents = c(7))) <- "Mesothelial"
Idents(covid_epi, cells = WhichCells(covid_epi, idents = c(6))) <- "Secretory - SCGB1A1+/SCGB3A2+"
Idents(covid_epi, cells = WhichCells(covid_epi, idents = c(5))) <- "KRT5-/KRT17+"
Idents(covid_epi, cells = WhichCells(covid_epi, idents = c(4))) <- "Ciliated"
Idents(covid_epi, cells = WhichCells(covid_epi, idents = c(2))) <- "AT1"
Idents(covid_epi, cells = WhichCells(covid_epi, idents = c(0,1,3))) <- "AT2"
covid_epi$celltype <- Idents(covid_epi)
DimPlot(covid_epi, label = T)



#Doublet exclusion stromal cells
covid_stromal <- subset(adata, idents = c('Stromal'))
covid_stromal <- SCTransform(covid_stromal, vars.to.regress = 'percent.mt')
covid_stromal <- RunPCA(covid_stromal, verbose = F)
covid_stromal <- RunUMAP(covid_stromal, dims = 1:20)
covid_stromal <- FindNeighbors(covid_stromal, dims = 1:20)
covid_stromal <- FindClusters(covid_stromal, resolution = 1.0)
DimPlot(covid_stromal, label = T)
DimPlot(covid_stromal, group.by = 'orig.ident')

DotPlot(covid_stromal, features = marker_genes)+ theme(axis.text.x = element_text(angle = 45, hjust=1))

#remove doublets
covid_stromal <- subset(covid_stromal, idents = c(0:2,4,5,6,8,9,10))
covid_stromal <- SCTransform(covid_stromal, vars.to.regress = 'percent.mt')
covid_stromal <- RunPCA(covid_stromal, verbose = F)
covid_stromal <- RunUMAP(covid_stromal, dims = 1:20)
covid_stromal <- FindNeighbors(covid_stromal, dims = 1:20)
covid_stromal <- FindClusters(covid_stromal, resolution = 1.0)
DimPlot(covid_stromal, label = T)
DimPlot(covid_stromal, group.by = 'orig.ident')

DotPlot(covid_stromal, features = marker_genes)+ theme(axis.text.x = element_text(angle = 45, hjust=1))

#remove doublets
covid_stromal <- subset(covid_stromal, idents = c(0:6,8,9))
covid_stromal <- SCTransform(covid_stromal, vars.to.regress = 'percent.mt')
covid_stromal <- RunPCA(covid_stromal, verbose = F)
covid_stromal <- RunUMAP(covid_stromal, dims = 1:20)
covid_stromal <- FindNeighbors(covid_stromal, dims = 1:20)
covid_stromal <- FindClusters(covid_stromal, resolution = 1.5)
DimPlot(covid_stromal, label = T)
DimPlot(covid_stromal, group.by = 'orig.ident')

stromal_genes <- c("PECAM1", 'CA4', 'HEY1', 'PLVAP', "ACKR1", 'VWF', 'APLN', 'APLNR', 'COL4A1', 'CCL21', 'ACTA2', 'CSPG4', 'COL1A1', 'COL3A1', 'DCN', 'PLIN2', 'HAS1', 'WNT5A', 'WNT2', 'TCF21', 'CTHRC1', 'ELN', 'PDGFRA', 'PDGFRB', 'WT1', "SFRP2", 'PI16')
DotPlot(covid_stromal, features = stromal_genes)+ theme(axis.text.x = element_text(angle = 45, hjust=1))

#high-level annotation
Idents(covid_stromal) <- 'seurat_clusters'
Idents(covid_stromal, cells = WhichCells(covid_stromal, idents = c(14))) <- "Venule"
Idents(covid_stromal, cells = WhichCells(covid_stromal, idents = c(13))) <- "Capillary - CA4+"
Idents(covid_stromal, cells = WhichCells(covid_stromal, idents = c(8,12))) <- "Pericyte"
Idents(covid_stromal, cells = WhichCells(covid_stromal, idents = c(11))) <- "Lymphatic"
Idents(covid_stromal, cells = WhichCells(covid_stromal, idents = c(10))) <- "SMC"
Idents(covid_stromal, cells = WhichCells(covid_stromal, idents = c(9))) <- "Endothelial - peribronchiolar"
Idents(covid_stromal, cells = WhichCells(covid_stromal, idents = c(7))) <- "Adventitial FB"
Idents(covid_stromal, cells = WhichCells(covid_stromal, idents = c(1,5,6))) <- "FB - WNT2+"
Idents(covid_stromal, cells = WhichCells(covid_stromal, idents = c(2,4))) <- "MyoFB"
Idents(covid_stromal, cells = WhichCells(covid_stromal, idents = c(3))) <- "Arteriole"
Idents(covid_stromal, cells = WhichCells(covid_stromal, idents = c(0))) <- "Capillary"
covid_stromal$celltype <- Idents(covid_stromal)


#Doublet exclusion immune cells
covid_immune <- subset(adata, idents = c('Immune'))
covid_immune <- SCTransform(covid_immune, vars.to.regress = 'percent.mt')
covid_immune <- RunPCA(covid_immune, verbose = F)
covid_immune <- RunUMAP(covid_immune, dims = 1:20)
covid_immune <- FindNeighbors(covid_immune, dims = 1:20)
covid_immune <- FindClusters(covid_immune, resolution = 1.0)
DimPlot(covid_immune, label = T)
DimPlot(covid_immune, group.by = 'orig.ident')

DotPlot(covid_immune, features = marker_genes)+ theme(axis.text.x = element_text(angle = 45, hjust=1))

#remove doublets
covid_immune <- subset(covid_immune, idents = c(0:11,14:18))
covid_immune <- SCTransform(covid_immune, vars.to.regress = 'percent.mt')
covid_immune <- RunPCA(covid_immune, verbose = F)
covid_immune <- RunUMAP(covid_immune, dims = 1:25)
covid_immune <- FindNeighbors(covid_immune, dims = 1:25)
covid_immune <- FindClusters(covid_immune, resolution = 1.0)
DimPlot(covid_immune, label = T)
DimPlot(covid_immune, group.by = 'orig.ident')

DotPlot(covid_immune, features = marker_genes)+ theme(axis.text.x = element_text(angle = 45, hjust=1))

immune_genes <- c("PTPRC", "LYZ", "CD68", "CD14", "FCGR3A", "CD86", 'HLA-DRA', "SPP1", "CD3E", "CD4", "CD8A", "PDCD1", "IL7R", "FOXP3", "GNLY", "CPA3", "CD19", "JCHAIN", "IRF7", "MKI67", 'MZB1', 'GZMA', 'GZMB','GZMK')
DotPlot(covid_immune, features = immune_genes)+ theme(axis.text.x = element_text(angle = 45, hjust=1))

#high-level annotation
Idents(covid_immune) <- 'seurat_clusters'
Idents(covid_immune, cells = WhichCells(covid_immune, idents = c(18))) <- "pDC"
Idents(covid_immune, cells = WhichCells(covid_immune, idents = c(17))) <- "Proliferating macrophage"
Idents(covid_immune, cells = WhichCells(covid_immune, idents = c(16))) <- "Mast"
Idents(covid_immune, cells = WhichCells(covid_immune, idents = c(11,15))) <- "Plasma"
Idents(covid_immune, cells = WhichCells(covid_immune, idents = c(14))) <- "Proliferating T cells"
Idents(covid_immune, cells = WhichCells(covid_immune, idents = c(13))) <- "B cells"
Idents(covid_immune, cells = WhichCells(covid_immune, idents = c(12))) <- "cDC1"
Idents(covid_immune, cells = WhichCells(covid_immune, idents = c(10))) <- "Tcell - stressed"
Idents(covid_immune, cells = WhichCells(covid_immune, idents = c(9))) <- "Monocyte"
Idents(covid_immune, cells = WhichCells(covid_immune, idents = c(5,7,8))) <- "Monocyte-derived macrophage"
Idents(covid_immune, cells = WhichCells(covid_immune, idents = c(6))) <- "Alveolar macrophage"
Idents(covid_immune, cells = WhichCells(covid_immune, idents = c(4))) <- "Inflammatory monocyte"
Idents(covid_immune, cells = WhichCells(covid_immune, idents = c(3))) <- "CD8"
Idents(covid_immune, cells = WhichCells(covid_immune, idents = c(2))) <- "CD4"
Idents(covid_immune, cells = WhichCells(covid_immune, idents = c(1))) <- "NKT"
Idents(covid_immune, cells = WhichCells(covid_immune, idents = c(0))) <- "NK"
covid_immune$celltype <- Idents(covid_immune)

#reallocate mesothelial
covid_merge <- merge(covid_epi, y = c(covid_stromal))
Idents(covid_merge, cells = WhichCells(covid_merge, idents = c('AT1', 'AT2', 'KRT5-/KRT17+', 'Secretory - SCGB1A1+/SCGB3A2+', 'Ciliated'))) <- "Epithelial"
Idents(covid_merge, cells = WhichCells(covid_merge, idents = c(
  "Venule",
  "Capillary - CA4+",
  "Pericyte",
  "Lymphatic",
  "SMC",
  "Endothelial - peribronchiolar",
  "Adventitial FB",
  "FB - WNT2+",
  "MyoFB",
  "Arteriole",
  "Capillary",
  'Mesothelial'))) <- 'Stromal'
covid_merge$population <- Idents(covid_merge)

Idents(covid_merge) <- 'celltype'

covid_epi2 <-subset(covid_merge, idents = c('AT1', 'AT2', 'KRT5-/KRT17+', 'Secretory - SCGB1A1+/SCGB3A2+', 'Ciliated'))
covid_stromal2 <-subset(covid_merge, idents = c(
"Venule",
"Capillary - CA4+",
"Pericyte",
"Lymphatic",
"SMC",
"Endothelial - peribronchiolar",
"Adventitial FB",
"FB - WNT2+",
"MyoFB",
"Arteriole",
"Capillary",
'Mesothelial'
))

#save subobjects
saveRDS(covid_epi2, file = '~/Desktop/COVID/scratch/covid_epi.rds')
saveRDS(covid_stromal2, file = '~/Desktop/COVID/scratch/covid_stromal.rds')
saveRDS(covid_immune, file = '~/Desktop/COVID/scratch/covid_immune.rds')



#####################################################################

#Integrate with Habermann Science Advances 2020 dataset

sci_adv <- readRDS(file = '~/Downloads/GSE135893_ILD_annotated_fullsize.rds')
DefaultAssay(sci_adv) <- 'RNA'
sci_adv[['SCT']] <- NULL

#assign disease status
Idents(sci_adv) <- 'Status'
Idents(sci_adv, cells = WhichCells(sci_adv, idents = c('Control'))) <- "Control"
Idents(sci_adv, cells = WhichCells(sci_adv, idents = c('ILD'))) <- "Chronic ILD"
sci_adv$Status <- Idents(sci_adv)


#create_population_subobjects
Idents(sci_adv) <- 'population'
sci_adv_epi <- subset(sci_adv, idents = c("Epithelial"))
sci_adv_immune <- subset(sci_adv, idents = c("Immune"))
sci_adv_stromal <- subset(sci_adv, idents = c("Endothelial", 'Mesenchymal'))

saveRDS(sci_adv_epi, file = '~/Desktop/COVID/scratch/sci_adv_epi.rds')
saveRDS(sci_adv_immune, file = '~/Desktop/COVID/scratch/sci_adv_immune.rds')
saveRDS(sci_adv_stromal, file = '~/Desktop/COVID/scratch/sci_adv_stromal.rds')





#integrate epithelial subset
epithelial <- merge(sci_adv_epi, y=c(covid_epi2))

#integrate
epithelial.list <- SplitObject(epithelial, split.by = "Status")
for (i in names(epithelial.list)) {
  epithelial.list[[i]] <- SCTransform(epithelial.list[[i]], verbose = FALSE)
}

epithelial.features <- SelectIntegrationFeatures(object.list = epithelial.list, nfeatures = 3000)

epithelial.list <- PrepSCTIntegration(object.list = epithelial.list, anchor.features = epithelial.features)

reference_dataset <- which(names(epithelial.list) == "Chronic ILD")

epithelial.anchors <- FindIntegrationAnchors(object.list = epithelial.list, normalization.method = "SCT", 
                                             anchor.features = epithelial.features, reference = reference_dataset)
epithelial.integrated <- IntegrateData(anchorset = epithelial.anchors, normalization.method = "SCT")
epithelial.integrated <- RunPCA(object = epithelial.integrated, verbose = FALSE)
epithelial.integrated <- RunUMAP(object = epithelial.integrated, dims = 1:20)

DimPlot(epithelial.integrated, group.by=c("Status"))

epithelial.integrated <- FindNeighbors(epithelial.integrated, dims=1:20)
epithelial.integrated <- FindClusters(epithelial.integrated, resolution=0.6)
DimPlot(epithelial.integrated, label = T)
DimPlot(epithelial.integrated, split.by = 'Status', label = T)

DefaultAssay(epithelial.integrated) <- 'SCT'

epithelial_markers <- c("EPCAM", "TP63", "KRT5", "KRT17", "FOXJ1", "TP73", "SCGB1A1", "SCGB3A1", "SCGB3A2", "MUC5B", "MUC5AC", 'AGER', "HOPX", "SFTPC", "SFTPA1", "ABCA3", "MKI67", "SFN", "COL1A1", "CALCA", "FOXI1")
DotPlot(epithelial.integrated, features = epithelial_markers)+ theme(axis.text.x = element_text(angle = 45, hjust=1))


Idents(epithelial.integrated) <- 'seurat_clusters'
Idents(epithelial.integrated, cells = WhichCells(epithelial.integrated, idents = c(19))) <- "Proliferating epithelial"
Idents(epithelial.integrated, cells = WhichCells(epithelial.integrated, idents = c(0,1,10,17,18))) <- "Ciliated"
Idents(epithelial.integrated, cells = WhichCells(epithelial.integrated, idents = c(7,16))) <- "Secretory - SCGB1A1+/SCGB3A2+"
Idents(epithelial.integrated, cells = WhichCells(epithelial.integrated, idents = c(9,14))) <- "Differentiating ciliated"
Idents(epithelial.integrated, cells = WhichCells(epithelial.integrated, idents = c(13))) <- "KRT5-/KRT17+"
Idents(epithelial.integrated, cells = WhichCells(epithelial.integrated, idents = c(12))) <- "AT1"
Idents(epithelial.integrated, cells = WhichCells(epithelial.integrated, idents = c(2,8,11,15))) <- "AT2"
Idents(epithelial.integrated, cells = WhichCells(epithelial.integrated, idents = c(6))) <- "Transitional"
Idents(epithelial.integrated, cells = WhichCells(epithelial.integrated, idents = c(5))) <- "Basal"
Idents(epithelial.integrated, cells = WhichCells(epithelial.integrated, idents = c(4))) <- "Secretory - MUC5B+"
Idents(epithelial.integrated, cells = WhichCells(epithelial.integrated, idents = c(3))) <- "Secretory - SCGB3A2+"
epithelial.integrated$types <- Idents(epithelial.integrated)

saveRDS(epithelial.integrated, file = '~/Desktop/COVID/scratch/epithelial_integrated.rds')

Idents(epithelial.integrated) <- 'Status'
covid_epi_integrated <- subset(epithelial.integrated, idents = c('COVID'))
saveRDS(covid_epi_integrated, file = '~/Desktop/COVID/scratch/covid_epi_integrated.rds')





#immune integration
immune <- merge(sci_adv_immune, y=c(covid_immune))

#integrate
immune.list <- SplitObject(immune, split.by = "Status")
for (i in names(immune.list)) {
  immune.list[[i]] <- SCTransform(immune.list[[i]], verbose = FALSE)
}

immune.features <- SelectIntegrationFeatures(object.list = immune.list, nfeatures = 3000)

immune.list <- PrepSCTIntegration(object.list = immune.list, anchor.features = immune.features)

reference_dataset <- which(names(immune.list) == "Chronic ILD")

immune.anchors <- FindIntegrationAnchors(object.list = immune.list, normalization.method = "SCT", 
                                         anchor.features = immune.features, reference = reference_dataset)
immune.integrated <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")
immune.integrated <- RunPCA(object = immune.integrated, verbose = FALSE)
immune.integrated <- RunUMAP(object = immune.integrated, dims = 1:25)

DimPlot(immune.integrated, group.by=c("Status"))

immune.integrated <- FindNeighbors(immune.integrated, dims=1:25)
immune.integrated <- FindClusters(immune.integrated, resolution=0.6)
DimPlot(immune.integrated, label = T)
DimPlot(immune.integrated, split.by = 'Status', label = T)
DimPlot(immune.integrated, group.by = 'Status')

DefaultAssay(immune.integrated) <- 'SCT'

immune_genes <- c("PTPRC", "LYZ", "CD68", 'PPARG', "CD14", "FCGR3A", "CD86", 'HLA-DRA', "SPP1", 'MRC1', 'IL1B', "CD3E", "CD4", "CD8A", "PDCD1", "IL7R", "FOXP3", "GNLY", "NKG7", "CPA3", "CD19", "JCHAIN", "IRF7", "MKI67")
DotPlot(immune.integrated, features = immune_genes)+ theme(axis.text.x = element_text(angle = 45, hjust=1))

Idents(immune.integrated) <- 'seurat_clusters'
Idents(immune.integrated, cells = WhichCells(immune.integrated, idents = c(23))) <- "pDC"
Idents(immune.integrated, cells = WhichCells(immune.integrated, idents = c(21))) <- "NKT"
Idents(immune.integrated, cells = WhichCells(immune.integrated, idents = c(20))) <- "Proliferating T-cells"
Idents(immune.integrated, cells = WhichCells(immune.integrated, idents = c(18))) <- "Mast"
Idents(immune.integrated, cells = WhichCells(immune.integrated, idents = c(17))) <- "Plasma"
Idents(immune.integrated, cells = WhichCells(immune.integrated, idents = c(16))) <- "B-cells"
Idents(immune.integrated, cells = WhichCells(immune.integrated, idents = c(5,15))) <- "CD4"
Idents(immune.integrated, cells = WhichCells(immune.integrated, idents = c(14))) <- "Proliferating macrophage"
Idents(immune.integrated, cells = WhichCells(immune.integrated, idents = c(10,12))) <- "Inflammatory monocyte"
Idents(immune.integrated, cells = WhichCells(immune.integrated, idents = c(11))) <- "NK"
Idents(immune.integrated, cells = WhichCells(immune.integrated, idents = c(8))) <- "CD8"
Idents(immune.integrated, cells = WhichCells(immune.integrated, idents = c(7))) <- "Monocyte"
Idents(immune.integrated, cells = WhichCells(immune.integrated, idents = c(6))) <- "Macrophage - CD206hi"
Idents(immune.integrated, cells = WhichCells(immune.integrated, idents = c(0,4,19))) <- "Alveolar macrophage"
Idents(immune.integrated, cells = WhichCells(immune.integrated, idents = c(3,9))) <- "cDC"
Idents(immune.integrated, cells = WhichCells(immune.integrated, idents = c(1))) <- "Macrophage - SPP1+"
Idents(immune.integrated, cells = WhichCells(immune.integrated, idents = c(2,13,22))) <- "Monocyte-derived macrophage"
immune.integrated$types <- Idents(immune.integrated)


saveRDS(immune.integrated, file = '~/Desktop/COVID/scratch/immune_integrated.rds')

immune.integrated <- UpdateSeuratObject(immune.integrated)
SaveH5Seurat(immune.integrated, filename = "immune_integrated.h5Seurat", overwrite = T)
Convert("immune_integrated.h5Seurat", dest = "h5ad")



#stromal integration

sci_adv_stromal <- readRDS(file = '~/Desktop/COVID/scratch/sci_adv_stromal.rds')
covid_stromal2 <- readRDS(file = '~/Desktop/COVID/scratch/covid_stromal.rds')

stromal <- merge(sci_adv_stromal, y=c(covid_stromal2))

#integrate
stromal.list <- SplitObject(stromal, split.by = "Status")
for (i in names(stromal.list)) {
  stromal.list[[i]] <- SCTransform(stromal.list[[i]], verbose = FALSE)
}

stromal.features <- SelectIntegrationFeatures(object.list = stromal.list, nfeatures = 3000)

stromal.list <- PrepSCTIntegration(object.list = stromal.list, anchor.features = stromal.features)

reference_dataset <- which(names(stromal.list) == "Chronic ILD")

stromal.anchors <- FindIntegrationAnchors(object.list = stromal.list, normalization.method = "SCT", 
                                          anchor.features = stromal.features, reference = reference_dataset)
stromal.integrated <- IntegrateData(anchorset = stromal.anchors, normalization.method = "SCT")
stromal.integrated <- RunPCA(object = stromal.integrated, verbose = FALSE)
stromal.integrated <- RunUMAP(object = stromal.integrated, dims = 1:25)

DimPlot(stromal.integrated, group.by=c("Status"))

stromal.integrated <- FindNeighbors(stromal.integrated, dims=1:25)
stromal.integrated <- FindClusters(stromal.integrated, resolution=1.0)
DimPlot(stromal.integrated, label = T)
DimPlot(stromal.integrated, split.by = 'Status', label = T)
DimPlot(stromal.integrated, group.by = 'Status')

stromal_genes <- c("PECAM1", 'HEY1', 'CA4', 'PLVAP', 'COL15A1', "ACKR1", 'VWF', 'APLNR', 'CCL21','ACTA2', 'CSPG4', 'COL1A1', 'COL3A1', 'DCN', 'PLIN2', 'HAS1', 'WNT5A', 'WNT2', 'TCF21', 'CTHRC1', 'PI16', 'PDGFRA', 'PDGFRB', 'WT1', "MKI67")
DotPlot(stromal.integrated, features = stromal_genes)+ theme(axis.text.x = element_text(angle = 45, hjust=1))

#remove doublets
stromal.integrated <- subset(stromal.integrated, idents = c(0:21,23,24))

#annotate cell clusters
Idents(stromal.integrated) <- 'seurat_clusters'
Idents(stromal.integrated, cells = WhichCells(stromal.integrated, idents = c(5,24))) <- "Capillary - CA4+"
Idents(stromal.integrated, cells = WhichCells(stromal.integrated, idents = c(0,15,17))) <- "Capillary"
Idents(stromal.integrated, cells = WhichCells(stromal.integrated, idents = c(3,21))) <- "FB - WNT2"
Idents(stromal.integrated, cells = WhichCells(stromal.integrated, idents = c(19))) <- "Mesothelial"
Idents(stromal.integrated, cells = WhichCells(stromal.integrated, idents = c(2,18))) <- "Arteriole"
Idents(stromal.integrated, cells = WhichCells(stromal.integrated, idents = c(14))) <- "Adventitial FB"
Idents(stromal.integrated, cells = WhichCells(stromal.integrated, idents = c(13,12))) <- "Lymphatic"
Idents(stromal.integrated, cells = WhichCells(stromal.integrated, idents = c(11))) <- "Pericyte"
Idents(stromal.integrated, cells = WhichCells(stromal.integrated, idents = c(10))) <- "SMC"
Idents(stromal.integrated, cells = WhichCells(stromal.integrated, idents = c(4,8,9,16,20,23))) <- "Peribronchiolar"
Idents(stromal.integrated, cells = WhichCells(stromal.integrated, idents = c(7))) <- "MyoFB"
Idents(stromal.integrated, cells = WhichCells(stromal.integrated, idents = c(6))) <- "Venule"
Idents(stromal.integrated, cells = WhichCells(stromal.integrated, idents = c(1))) <- "FB - HAS1+"
stromal.integrated$types <- Idents(stromal.integrated)


saveRDS(stromal.integrated, file = '~/Desktop/COVID/scratch/stromal_integrated.rds')

stromal.integrated <- UpdateSeuratObject(stromal.integrated)
SaveH5Seurat(stromal.integrated, filename = "stromal_integrated.h5Seurat", overwrite = T)
Convert("stromal_integrated.h5Seurat", dest = "h5ad")



#Extract and re-merge covid annotated dataset
Idents(epithelial.integrated) <- 'Status'
Idents(immune.integrated) <- 'Status'
Idents(stromal.integrated) <- 'Status'

covid_epi_annotated <- subset(epithelial.integrated, idents = c('COVID'))
covid_immune_annotated <- subset(immune.integrated, idents = c('COVID'))
covid_stromal_annotated <- subset(stromal.integrated, idents = c('COVID'))

covid <- merge(covid_epi_annotated, y=c(covid_immune_annotated, covid_stromal_annotated))
covid <- SCTransform(covid, vars.to.regress = 'percent.mt', return.only.var.genes = FALSE)
covid <- RunPCA(covid, verbose = F)
covid <- RunUMAP(covid, dims = 1:45)
Idents(covid) <- 'types'
DimPlot(covid, label = T)

setwd("~/Desktop/COVID/")

Idents(covid) <- 'types'

saveRDS(covid, file = "~/covid_annotated.rds")

#convert to h5ad
library(SeuratDisk)

#Set random seed to get reproducible results
set.seed(33)

#convert to h5ad files
covid <- UpdateSeuratObject(covid)
SaveH5Seurat(covid, filename = "covid.h5Seurat", overwrite = T)
Convert("covid.h5Seurat", dest = "h5ad")


#############
#cell proportion tables

Idents(stromal.integrated) <- 'types'

fb_subset <- subset(stromal.integrated, idents = c('FB - HAS1+', 'FB - WNT2', 'MyoFB', 'Adventitial FB'))
endo_subset <- subset(stromal.integrated, idents = c('Capillary', 'Capillary - CA4+', 'Arteriole', 'Venule', 'Peribronchiolar'))

endo_prop_table <- as.data.frame(prop.table(table(Idents(endo_subset), endo_subset$Status), margin = 2))
endo_prop_table <- endo_prop_table %>% 
  rename(
    Cell_type = Var1,
    Diagnosis = Var2
  )

ggplot(endo_prop_table, aes(fill=Cell_type, y=Freq, x=Diagnosis)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_brewer(palette = 'Paired') +
  ggtitle("Cell Types by Diagnosis") +
  theme_ipsum() +
  xlab("")


fb_prop_table <- as.data.frame(prop.table(table(Idents(fb_subset), fb_subset$Status), margin = 2))
fb_prop_table <- fb_prop_table %>% 
  rename(
    Cell_type = Var1,
    Diagnosis = Var2
  )

ggplot(fb_prop_table, aes(fill=Cell_type, y=Freq, x=Diagnosis)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_brewer(palette = 'Paired') +
  ggtitle("Cell Types by Diagnosis") +
  theme_ipsum() +
  xlab("")


Idents(immune.integrated) <- 'types'
myeloid <- subset(immune.integrated, idents = c('Alveolar macrophage', 'Macrophage - CD206hi', 'Monocyte-derived macrophage', 'cDC', 'Monocyte', 'Proliferating macrophage', 'Macrophage - SPP1+', 'Inflammatory monocyte'))
lymphoid <- subset(immune.integrated, idents = c('CD4', 'CD8', 'NKT', 'NK', 'B-cells', 'Plasma'))


myeloid_prop_table <- as.data.frame(prop.table(table(Idents(myeloid), myeloid$Status), margin = 2))
myeloid_prop_table <- myeloid_prop_table %>% 
  rename(
    Cell_type = Var1,
    Diagnosis = Var2
  )

ggplot(myeloid_prop_table, aes(fill=Cell_type, y=Freq, x=Diagnosis)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_brewer(palette = 'Paired') +
  ggtitle("Cell Types by Diagnosis") +
  theme_ipsum() +
  xlab("")


lymphoid_prop_table <- as.data.frame(prop.table(table(Idents(lymphoid), lymphoid$Status), margin = 2))
lymphoid_prop_table <- lymphoid_prop_table %>% 
  rename(
    Cell_type = Var1,
    Diagnosis = Var2
  )

ggplot(lymphoid_prop_table, aes(fill=Cell_type, y=Freq, x=Diagnosis)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_brewer(palette = 'Paired') +
  ggtitle("Cell Types by Diagnosis") +
  theme_ipsum() +
  xlab("")

Idents(epithelial.integrated) <- 'types'
epi_prop_table <- as.data.frame(prop.table(table(Idents(epithelial.integrated), epithelial.integrated$Status), margin = 2))
epi_prop_table <- epi_prop_table %>% 
  rename(
    Cell_type = Var1,
    Diagnosis = Var2
  )

ggplot(epi_prop_table, aes(fill=Cell_type, y=Freq, x=Diagnosis)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_brewer(palette = 'Paired') +
  ggtitle("Cell Types by Diagnosis") +
  theme_ipsum() +
  xlab("")



####cell number tables
Idents(epithelial.integrated) <- 'types'
Idents(stromal.integrated) <- 'types'
Idents(immune.integrated) <- 'types'

immune_table <- as.data.frame(table(Idents(immune.integrated), immune.integrated$Status), margin = 2)
immune_table <- immune_table %>% 
  rename(
    Cell_type = Var1,
    Status = Var2
  )
write.csv(immune_table, file = '~/Desktop/COVID/immune_cells_integrated.csv')


stromal_table <- as.data.frame(table(Idents(stromal.integrated), stromal.integrated$Status), margin = 2)
stromal_table <- stromal_table %>% 
  rename(
    Cell_type = Var1,
    Status = Var2
  )
write.csv(stromal_table, file = '~/Desktop/COVID/stromal_cells_integrated.csv')

epithelial_table <- as.data.frame(table(Idents(epithelial.integrated), epithelial.integrated$Status), margin = 2)
epithelial_table <- epithelial_table %>% 
  rename(
    Cell_type = Var1,
    Status = Var2
  )
write.csv(epithelial_table, file = '~/Desktop/COVID/epithelial_cells_integrated.csv')

#####
#output objects to h5ad files


#convert to h5ad
epithelial.integrated <- SCTransform(epithelial.integrated, batch_var='Status', return.only.var.genes = FALSE, vars.to.regress = 'percent.mt')

library(SeuratDisk)
setwd("~/Desktop/COVID/scratch/")

#Set random seed to get reproducible results
set.seed(33)

#convert to h5ad files
epithelial.integrated <- UpdateSeuratObject(epithelial.integrated)
SaveH5Seurat(epithelial.integrated, filename = "epithelial_integrated.h5Seurat", overwrite = T)
Convert("epithelial_integrated.h5Seurat", dest = "h5ad")

stromal.integrated <- SCTransform(stromal.integrated, batch_var='Status', return.only.var.genes = FALSE, vars.to.regress = 'percent.mt')
stromal.integrated <- UpdateSeuratObject(stromal.integrated)
SaveH5Seurat(stromal.integrated, filename = "stromal_integrated.h5Seurat", overwrite = T)
Convert("stromal_integrated.h5Seurat", dest = "h5ad")

immune.integrated <- SCTransform(immune.integrated, batch_var='Status', return.only.var.genes = FALSE)
immune.integrated <- UpdateSeuratObject(immune.integrated)
SaveH5Seurat(immune.integrated, filename = "immune_integrated.h5Seurat", overwrite = T)
Convert("immune_integrated.h5Seurat", dest = "h5ad")


VlnPlot(myeloid, features = c('FN1'), split.by = 'Status', pt.size = 0)
VlnPlot(myeloid, features = c('VCAN'), split.by = 'Status', pt.size = 0)
VlnPlot(myeloid, features = c('IL1B'), split.by = 'Status', pt.size = 0)



