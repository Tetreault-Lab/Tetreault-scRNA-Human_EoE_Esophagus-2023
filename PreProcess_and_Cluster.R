###########################################
######  Seurat Object Generation  #########
###########################################
  
  
# Load Libraries

library(ggplot2)
library(Seurat)
library(viridis)
library(tidyselect)
library(ggpubr)
library(tidyverse)
library(purrr)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(gtable)
library(reshape2)
library(tidyverse)
library(Biobase)
library(dplyr)
library(sctransform)
library(viridisLite)

# Global Functions -------------------------------------------------------------

# Loads an RData file, and returns it
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

# Negates %in% operator
`%!in%` = Negate(`%in%`)



# Create seurat object for each sample  ----------------------------------------
setwd("/projects/b1042/MPLab/mhc0155/EoE_Project/")

SeuratObjects <- list()
names <- c()
for(x in list.dirs(path = "./RawData",recursive = F, full.names = T)){
  y <- Read10X(paste0(x, "/filtered_feature_bc_matrix")) %>%
    CreateSeuratObject(project = "EoE_Esophageal_Biopsies", 
                       assay = "RNA", 
                       min.features = 150, 
                       min.cells =5)
  y[["percent.Mito"]] <- PercentageFeatureSet(y, pattern = "^MT-")
  y[["percent.Ribo"]] <- PercentageFeatureSet(y, pattern = "^RP")
  y <- y[,y@meta.data %>% 
           filter(nFeature_RNA > 100 & nCount_RNA > 400) %>% 
           row.names()] %>%
    Seurat::SCTransform(vst.flavor="v3") %>% 
    CellCycleScoring(s.features = cc.genes$s.genes, 
                     g2m.features = cc.genes$g2m.genes, 
                     set.ident = TRUE)
  y$Patient_Region <- gsub("./RawData/","",gsub("-", "_", x))
  SeuratObjects <- append(SeuratObjects,y)
  names <- c(names, gsub("./RawData/","",gsub("-", "_", x)))
}
names(SeuratObjects) <- names

# Merge initial samples for feature comparison
Initial_Merge <- Reduce(function(x,y) merge(x,y), SeuratObjects)

#Add metadata

PatientInfo <- read.csv("./Preprocess_Cluster/PatientInfo.csv")[,2:5]
Initial_Merge@meta.data <- Initial_Merge@meta.data %>%
  mutate(OriginalNames = rownames(.)) %>%
  mutate(Order = seq(1, nrow(.), by =1)) %>%
  merge(., PatientInfo, by = "Patient_Region") %>%
  arrange(Order) %>%
  column_to_rownames("OriginalNames") %>%
  dplyr::select(-Order)




# Visualize QC features  -------------------------------------------------------

palette <- c(colorRampPalette(c("darkgreen", "lightgreen"))(12), 
             colorRampPalette(c("darkblue", "lightblue"))(11))
features = c("nFeature_RNA", "nCount_RNA", "percent.Mito")

p1 <- VlnPlot(Initial_Merge,
              features = features[1],
              group.by = "Patient_Region",
              cols = palette,
              pt.size = 0) +
  theme(axis.text.x = element_blank(),
        legend.position = "none")+
  ylab(paste(features[1])) + xlab("") +
  ggtitle("")
p2 <- VlnPlot(Initial_Merge,
              features = features[2],
              group.by = "Patient_Region",
              cols = palette,
              pt.size = 0) +
  theme(axis.text.x = element_blank(),
        legend.position = "none")+
  ylab(paste(features[2])) + xlab("") +
  ggtitle("")
p3 <- VlnPlot(Initial_Merge,
              features = features[3],
              group.by = "Patient_Region",
              cols = palette,
              pt.size = 0) +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 9))+
  ylab(paste(features[3])) + xlab("") +
  ggtitle("") +
  geom_hline(yintercept = 20, linetype = "dashed", color = "grey")


p1/p2/p3



# Integrate samples ------------------------------------------------------------

features <- SelectIntegrationFeatures(object.list = SeuratObjects, 
                                      nfeatures = 3000)
SeuratObjects <- PrepSCTIntegration(object.list = SeuratObjects, 
                                    anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = SeuratObjects, 
                                  normalization.method = "SCT", 
                                  anchor.features = features, 
                                  reference = which(names(SeuratObjects) %in% c("P0597_P","HCE047_P")),
                                  reduction ="rpca")
Whole_Integrated_Obj <- IntegrateData(anchorset = anchors, 
                                      normalization.method = "SCT")
Whole_Integrated_Obj <- RunPCA(Whole_Integrated_Obj, npcs=150) %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(reduction = "pca", 
                dims = 1:20) %>%
  FindClusters(resolution = 0.5)

# Normalize RNA Assay
DefaultAssay(Whole_Integrated_Obj) <- "RNA"
Whole_Integrated_Obj <- NormalizeData(Whole_Integrated_Obj)

# Save obejct
save(Whole_Integrated_Obj, file = "./Integration_Steps/Objects/Whole_Integrated_Obj")





# Determine which seurat clusters are epithelial -------------------------------

## Calculate epithelial gene signature
Whole_Integrated_Obj <- AddModuleScore(object = Whole_Integrated_Obj,
                                       features = list("EpiGenes" = c(
                                         # Basal epithelial markers
                                         "KRT14",
                                         "DST",
                                         "COL17A1",
                                         "KRT15",
                                         "KRT5",
                                         # Suprabasal epithelial markers
                                         "KRT13",
                                         "KRT4",
                                         "KRT6A",
                                         # Superficial epithelial markers
                                         "KRT78",
                                         "FLG"
                                       )), 
                                       name = "Epithelial_gene_sig")

## Visualize Seurat clusters

UmapTheme <- theme(aspect.ratio = 1,
                   axis.ticks = element_blank(),
                   axis.text = element_blank(),
                   axis.title = element_blank(),
                   plot.title = element_text(size=8),
                   legend.text = element_text(size=7),
                   legend.title = element_text(size = 8))

p1 <- DimPlot(Whole_Integrated_Obj, 
              group.by = "seurat_clusters", 
              label = T,
              repel = T,
              raster = T,
              label.size = 3.5) + 
  UmapTheme +
  theme(legend.position = 'bottom', 
        legend.spacing.x = unit(.2, 'cm'),
        legend.text = element_text(margin = margin(t = -3)),
        legend.key.size = unit(.9,"line")) +
  guides(color = guide_legend(label.position = "right",
                              nrow =5,
                              override.aes = list(size=3))) +
  ggtitle("Seurat Clusters")

## Visualize epithelial gene signature
p2 <- FeaturePlot(Whole_Integrated_Obj, 
                  features = c("Epithelial_gene_sig1"),
                  min.cutoff = -2,
                  raster = T) + 
  UmapTheme +
  theme(legend.position = 'bottom', 
        legend.justification = "center") +
  scale_color_viridis(option = "turbo") + 
  ggtitle("Epithelial Gene Signature")

## Annotate clusters based on epithelial or non-epithelial status
Whole_Integrated_Obj$seurat_clusters_Annotated <- Whole_Integrated_Obj@meta.data %>%
  mutate(seurat_clusters_Annotated = ifelse(seurat_clusters %in% c(0:11), 
                                            paste0(as.character(comb_seurat_clusters), "_Epi"),
                                            paste0(as.character(comb_seurat_clusters), "_Non_Epi")),
         seurat_cluster_Annotated = factor(seurat_cluster_Annotated, levels = c(paste0(c(0:11), "_Epi"),
                                                                                paste0(c(12:17), "_Non_Epi")))) %>%
  select(seurat_cluster_Annotated)

## Visualize clusters based on epithelial or non-epithelial status
p3 <- DimPlot(Whole_Integrated_Obj, 
              group.by = "seurat_clusters_Annotated", 
              label = T,
              repel = T,
              raster = T,
              label.size = 3) + 
  UmapTheme +
  theme(legend.position = 'bottom', 
        legend.spacing.x = unit(.2, 'cm'),
        legend.text = element_text(margin = margin(t = -3)),
        legend.key.size = unit(.9,"line")) +
  guides(color = guide_legend(label.position = "right",
                              nrow =5,
                              override.aes = list(size=3))) +
  ggtitle("Seurat Clusters: Top level annotation")

## View cluster summary

p1|p2|p3




# Excise clusters with low depth  ----------------------------------------------


## Our dataset is predominantly epithelial, so epithelial and non-epithelial 
## cell types will have differential sequencing depth / nCount/nFeature thresholds

### Compare all epithelial cells to all non epithelial cells

Whole_Integrated_Obj$Epi_NonEpi <- Whole_Integrated_Obj@meta.data %>%
  mutate(Epi_NonEpi = ifelse(seurat_clusters %in% c(0:11), 
                             "Epi",
                             "Non_Epi"),
         Epi_NonEpi = factor(Epi_NonEpi, levels = c("Epi",
                                                    "Non_Epi"))) %>%
  select(Epi_NonEpi)



p1 <- DimPlot(Whole_Integrated_Obj, 
              group.by = "seurat_clusters_Annotated", 
              label = T,
              repel = T,
              raster = F,
              label.size = 2) + 
  UmapTheme +
  theme(legend.position = 'bottom', 
        legend.spacing.x = unit(.2, 'cm'),
        legend.text = element_text(margin = margin(t = -3)),
        legend.key.size = unit(.9,"line")) +
  guides(color = guide_legend(label.position = "right",
                              nrow =5,
                              override.aes = list(size=3))) +
  ggtitle("Epithelial cells vs Non-Epithelial cells")

p2 <- 
  # Proportions df
  as.data.frame(table(Whole_Integrated_Obj$Epi_NonEpi)) %>%
  mutate(percent = 100*(Freq / ncol(Whole_Integrated_Obj))) %>%
  arrange(desc(percent)) %>%
  filter(Var1 %in% levels(Whole_Integrated_Obj$Epi_NonEpi)) %>%
  mutate(caption = paste(round(percent,2), "% ", Var1, sep = "")) %>%
  mutate(caption = factor(caption, levels = .$caption)) %>%
  # Plot as pie chart
  ggplot(aes(x="", 
             y=Freq, 
             fill=caption)) +
  geom_bar(stat="identity", 
           width=1, 
           color="white", 
           size = 0.5) +
  coord_polar("y", 
              start=0) +
  theme_void() +
  theme(legend.title = element_text(size = 7),
        legend.text = element_text(size=7),
        legend.justification = c(.6,0.8),
        aspect.ratio =1,
        legend.key.size = unit(0.5, "cm")) +
  scale_fill_manual(name = "Cell Type",
                    values = c("red", 
                               "#2CB77B"))+
  ggtitle("") 



### Compare all epithelial cells to all non epithelial clusters

Whole_Integrated_Obj$Epi_OtherCellTypes <- Whole_Integrated_Obj@meta.data %>%
  mutate(Epi_OtherCellTypes = ifelse(Epi_NonEpi %in% "Epi", 
                                     as.character(Epi_NonEpi),
                                     as.character(seurat_clusters_Annotated)),
         Epi_OtherCellTypes = factor(Epi_OtherCellTypes, levels = c("Epi",
                                                                    levels(Whole_Integrated_Obj$seurat_clusters_Annotated)[13:18]))) %>%
  select(Epi_OtherCellTypes)

p3 <- DimPlot(Whole_Integrated_Obj, 
              group.by = "Epi_OtherCellTypes", 
              label = T,
              repel = T,
              raster = F,
              label.size = 2) + 
  UmapTheme +
  theme(legend.position = 'bottom', 
        legend.spacing.x = unit(.2, 'cm'),
        legend.text = element_text(margin = margin(t = -3)),
        legend.key.size = unit(.9,"line")) +
  guides(color = guide_legend(label.position = "right",
                              nrow =5,
                              override.aes = list(size=3))) +
  ggtitle("Epithelial cells vs Other cell type clusters")

p4 <-
  # Proportions df
  as.data.frame(table(Whole_Integrated_Obj$Epi_OtherCellTypes)) %>%
  mutate(percent = 100*(Freq / ncol(Whole_Integrated_Obj))) %>%
  arrange(desc(percent)) %>%
  filter(Var1 %in% levels(Whole_Integrated_Obj$Epi_OtherCellTypes)) %>%
  mutate(caption = paste(round(percent,2), "% ", Var1, sep = "")) %>%
  mutate(caption = factor(caption, levels = .$caption)) %>%
  # Pie chart
  ggplot(aes(x="", 
             y=Freq, 
             fill=caption)) +
  geom_bar(stat="identity", 
           width=1, 
           color="white", 
           size = 0.5) +
  coord_polar("y", start=0) +
  theme_void() +
  theme(legend.title = element_text(size = 7),
        legend.text = element_text(size=7),
        legend.justification = c(.6,0.8),
        aspect.ratio =1,
        legend.key.size = unit(0.5, "cm")) +
  scale_fill_manual(name = "Cell Type",
                    values = c("red",
                               rainbow(8)[2:8]))+
  ggtitle("") 


## View Proportion summary

(p1/p2) | (p3/p4) # epithelial cells are 85% of dataset, and thus will habe more reads, and higher nFeature / ncount then cell types with a <5% representation in the dataset




## Seperate object into epithelial & non-epithelial clusters


Whole_Integrated_Obj_Epi <- Whole_Integrated_Obj[,Whole_Integrated_Obj@meta.data %>% 
                                                   filter(seurat_clusters_Annotated %in% levels(Whole_Integrated_Obj$seurat_clusters_Annotated)[1:12]) %>%
                                                   row.names()]
Whole_Integrated_Obj_NonEpi <- Whole_Integrated_Obj[,Whole_Integrated_Obj@meta.data %>% 
                                                      filter(seurat_clusters_Annotated %in% levels(Whole_Integrated_Obj$seurat_clusters_Annotated)[13:18]) %>%
                                                      row.names()]

### Function to compare QC parameters across clusters

#' This function generates a stacked violin plot with QC features and additional thresholding annotations.
#'
#' @param object A Seurat object.
#' @param Title A title for the plot.
#' @param fHi Upper limit for the y-axis of nFeature_RNA.
#' @param cHi Upper limit for the y-axis of nCount_RNA.
#' @param mHi Upper limit for the y-axis of % mitochondrial transcripts.
#' @param fLine Y-coordinate for the thresholding dashed line for nFeature_RNA.
#' @param cLine Y-coordinate for the thresholding dashed line for nCount_RNA.
#' @param mLine Y-coordinate for the thresholding dashed line for % mitochondrial transcripts
#' @param AR    Aspect ratio of plot
#'
#' @return A stackeed violin plot.
ViolinParams <- function(object,
                         Title,
                         fHi,
                         cHi,
                         mHi,
                         fLine,
                         cLine,
                         mLine,
                         AR = 0.8){
  features = c("nFeature_RNA", "nCount_RNA", "percent.Mito")
  pa <- VlnPlot(object,
                features = features[1],
                group.by = "seurat_cluster_Annotated",
                pt.size = 0) +
    theme(axis.text.x = element_blank(),
          legend.position = "none",
          aspect.ratio = AR,
          axis.text.y = element_text(size = 6),
          plot.title = element_text(size = 9),
          axis.title.y = element_text(size = 7.5))+
    ylab(paste(features[1])) + xlab("") +
    ggtitle(paste(Title)) +
    geom_hline(yintercept = fLine, linetype = "dashed", color = "grey") + 
    ylim(c(0,fHi))+
    geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) 
  pb <- VlnPlot(object,
                features = features[2],
                group.by = "seurat_cluster_Annotated",
                pt.size = 0) +
    theme(axis.text.x = element_blank(),
          legend.position = "none",
          aspect.ratio = AR,
          axis.text.y = element_text(size = 6),
          axis.title.y = element_text(size = 7.5))+
    ylab(paste(features[2])) + xlab("") +
    ggtitle("") +
    ylim(c(0,cHi))  +
    geom_hline(yintercept = cLine, linetype = "dashed", color = "red")+
    geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) 
  pc <- VlnPlot(object,
                features = features[3],
                group.by = "seurat_cluster_Annotated",
                pt.size = 0) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 6),
          aspect.ratio = AR,
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(size = 7.5))+
    ylab(paste(features[3])) + xlab("") +
    ggtitle("") +
    geom_hline(yintercept = mLine, linetype = "dashed", color = "red") +
    ylim(0, mHi)+
    geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) 
  
  pfinal <- (pa/pb/pc)
  
  return(pfinal)
}


### Compare QC features across epithelial and non-epithelial objects


pa <- ViolinParams(Whole_Integrated_Obj_Epi,
                   Title = "Epithelial Cells",
                   fLine = quantile(Whole_Integrated_Obj_Epi$nFeature_RNA, 0.1),
                   cLine = quantile(Whole_Integrated_Obj_Epi$nCount_RNA, 0.1),
                   mLine = 15,
                   mHi = 20,
                   cHi = quantile(Whole_Integrated_Obj_Epi$nCount_RNA, 0.95),
                   fHi = quantile(Whole_Integrated_Obj_Epi$nFeature_RNA, 0.95),
                   AR = NULL)

pb <- ViolinParams(Whole_Integrated_Obj_NonEpi,
                   Title = "Non-Epithelial Cells",
                   fLine = quantile(Whole_Integrated_Obj_NonEpi$nFeature_RNA, 0.1),
                   cLine = quantile(Whole_Integrated_Obj_NonEpi$nCount_RNA, 0.1),
                   mLine = 15,
                   mHi = 20,
                   cHi = quantile(Whole_Integrated_Obj_NonEpi$nCount_RNA, 0.95),
                   fHi = quantile(Whole_Integrated_Obj_NonEpi$nFeature_RNA, 0.95))

pc <- FeaturePlot(Whole_Integrated_Obj_Epi, 
                  features = c("nFeature_RNA"),
                  max.cutoff = quantile(Whole_Integrated_Obj_Epi$nFeature_RNA, 0.95),
                  raster = T) + 
  UmapTheme +
  scale_color_viridis(option = "turbo") + 
  ggtitle("Epithelial Clusters: nFeature RNA")

pd <- FeaturePlot(Whole_Integrated_Obj_Epi, 
                  features = c("nCount_RNA"),
                  max.cutoff = quantile(Whole_Integrated_Obj_Epi$nCount_RNA, 0.95),
                  raster = T) + 
  UmapTheme +
  scale_color_viridis(option = "turbo") + 
  ggtitle("Epithelial Clusters: nCount RNA")


Whole_Integrated_Obj_Epi$lowSeqClusters <- Whole_Integrated_Obj_Epi@meta.data %>%
  mutate(lowSeqClusters = ifelse(seurat_clusters %in% c(8:9), 
                                 paste0(as.character(seurat_cluster_Annotated), "_Low_Seq_Depth"),
                                 "Normal_Seq_Depth"),
         lowSeqClusters = factor(lowSeqClusters, levels = c(paste0(c(8:9), "_Epi_Low_Seq_Depth"),
                                                            "Normal_Seq_Depth"))) %>%
  select(lowSeqClusters)

pe <- DimPlot(Whole_Integrated_Obj_Epi, 
              group.by = "lowSeqClusters", 
              label = F,
              raster = T) + 
  UmapTheme+
  ggtitle("Seurat Epi Clusters: by sequencing depth")

# View epithelial summary

pa|(pe/pc/pd)

# View non-epithelial summary

pb


# Confirm Low seq depth epithelial clusters do not exhibit maximum scaled intensity for any key epithelial marker 

gg_Legend <- function(df,
                      legend.textS = 9,
                      legend.titleS = 9,
                      title,
                      colors,
                      V,
                      H){
  p1 <- as_ggplot(cowplot::get_legend(df %>% 
                                        mutate(!!rlang::sym(colnames(.)[1]) := factor(!!rlang::sym(colnames(.)[1]), levels = unique(!!rlang::sym(colnames(.)[1])))) 
                                      %>%
                                        ggplot(aes(!!as.name(colnames(df)[1]), fill = !!as.name(colnames(df)[1]))) + 
                                        geom_bar() +
                                        theme(legend.text = element_text(size = legend.textS),
                                              legend.title = element_text(size = legend.titleS),
                                              legend.position = c(H, V),
                                              legend.justification = c(1, 1),
                                              legend.direction = "vertical") +
                                        scale_fill_manual(name = paste(title), 
                                                          values = colors)))
  return(p1)
}


## Epithelial gene markers
markers <- c(
  "KRT14", "COL17A1", "DST", # Progenitor / Basal
  "PCNA", "MKI67", "PTTG1", "CCND1", #Dividing
  "KRT13", "IVL", #Early Differentiation / Suprabasal
  "CNFN","KRT78","FLG" # Terminal Differentiation / Superficial
)

## Colors for axis labels 
gene_colors <- c(rep("orange", 3),
                 rep("blue", 4),
                 rep("#109E8A", 2),
                 rep("purple", 3))
cluster_colors <- c(rep("black", 10),
                    rep("red", 2))

## Construct dot plot
ggarrange(DotPlot(Whole_Integrated_Obj_Epi, 
                  group.by = "seurat_clusters_Annotated", 
                  features = markers ) +
            coord_flip(clip="off") + 
            theme(axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5, colour=cluster_colors, size =10),
                  panel.border = element_rect(color = "black", fill = NA, size = 1),
                  plot.margin = unit(c(2,0,2,0), "lines"),
                  aspect.ratio = 0.7,
                  legend.position = "top",
                  legend.text = element_text(size = 8),
                  legend.title = element_text(size = 9),
                  legend.justification = "center",
                  axis.text.y = element_text(colour=gene_colors, size =10)) +
            xlab("") + 
            ylab("") +
            
            ### Order clusters according to differentiation markers, with low seq depth clusters at the end
            scale_y_discrete(
              limits = levels(Whole_Integrated_Obj_Epi$seurat_clusters_Annotated)[c(7,1,8,6,12,4,5,2,11,3,9,10)],
              labels = c(levels(Whole_Integrated_Obj_Epi$seurat_clusters_Annotated)[c(7,1,8,6,12,4,5,2,11,3)],
                         paste0(levels(Whole_Integrated_Obj_Epi$seurat_clusters_Annotated)[9:10], "_low_seq_depth"))
            ) +
            scale_colour_gradient2(low="#382C84", mid="white", high="#EF0053", 
                                   breaks = c(-1.5,0,1.5)) + 
            
            ### add line annotation to indicate clusters with low seq depth
            annotation_custom(linesGrob(y=0:1, 
                                        x=1/12*10,  
                                        gp = gpar(col = "black", lty = 2, cex = 4)), 
                              xmin = -Inf, 
                              xmax = Inf, 
                              ymin = -Inf, 
                              ymax = Inf),
          
          ### add legend for gene category colors
          gg_Legend(df = data.frame("GeneCategory" = c(rep("Basal", 3),
                                                       rep("Dividing / Cell Cycle", 4),
                                                       rep("Suprabasal", 2),
                                                       rep("Superficial", 3)),
                                    "features.plot" = markers),
                    colors = c("orange", "blue", "#109E8A", "purple"),
                    title = "Gene Category",
                    H = 1,
                    V = 0.8),
          ncol = 2, 
          widths = c(4,1),
          heights = c(1,0.3))


# Confirm non-epi clusters cluster based on biological relevance and not low seq depth

## Gene markers
markers <- c(
  # T cells / NK
  "CD3D",
  "NKG7",
  #B Cells
  "CD79A",
  "IGHA1",
  # Macrophages / Dendritic Cells / Monocytes
  "CD68",
  "CD207",
  "CD14",
  # Mast Cells
  "KIT",
  "CPA3",
  # Endothelial Cells
  "VWF",
  "CDH5",
  # Fibroblasts
  "DCN", 
  "COL1A1", 
  "MYL9",
  # Smooth Muscle
  "MYH11",
  "CNN1")

## Colors for axis labels 
gene_colors <- c(rep(rainbow(10)[2], 4), # Lymphocytes
                 rep(rainbow(10)[4], 3), # MNPs
                 rep(rainbow(10)[6], 2), # Mast
                 rep(rainbow(10)[7], 2), # Endothelial
                 rep(rainbow(10)[9], 5)) # Mesenchymal
## Construct dot plot
ggarrange(DotPlot(Whole_Integrated_Obj_NonEpi, 
                  group.by = "seurat_clusters_Annotated", 
                  features = markers ) +
            coord_flip(clip="off") + 
            theme(axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5, size = 10),
                  panel.border = element_rect(color = "black", fill = NA, size = 1),
                  plot.margin = unit(c(2,0,2,0), "lines"),
                  aspect.ratio = 0.7,
                  legend.position = "top",
                  legend.text = element_text(size = 8),
                  legend.title = element_text(size = 9),
                  legend.justification = "center",
                  axis.text.y = element_text(colour=gene_colors, size = 10)) +
            xlab("") + 
            ylab("") +
            
            ### Order clusters according to differentiation markers, with low seq depth clusters at the end
            scale_y_discrete(limits = levels(Whole_Integrated_Obj_Epi$seurat_clusters_Annotated)[c(15,13,18,14, 17,16)]) +
            scale_colour_gradient2(low="#382C84", mid="white", high="#EF0053", 
                                   breaks = c(-1.5,0,1.5)),
          ### add legend for gene category colors
          gg_Legend(df = data.frame("GeneCategory" = c(rep("Lymphocytes", 4),
                                                       rep("MNPs", 3),
                                                       rep("Mast cells", 2),
                                                       rep("Endothelial", 2),
                                                       rep("Mesenchymal", 5)),
                                    "features.plot" = markers),
                    colors = c(
                      rainbow(10)[2],  # Lymphocytes
                      rainbow(10)[4],  # MNPs
                      rainbow(10)[6],  # Mast
                      rainbow(10)[7],  # Endothelial
                      rainbow(10)[9]   # Mesenchymal
                    ),
                    title = "Gene Category",
                    H = 1,
                    V = 0.8),
          ncol = 2, 
          widths = c(4,1),
          heights = c(1,0.3))


# Confirm non-epi clusters cluster based on biological relevance and not low seq depth

## Gene markers
markers <- c(
  # T cells / NK
  "CD3D",
  "NKG7",
  #B Cells
  "CD79A",
  "IGHA1",
  # Macrophages / Dendritic Cells / Monocytes
  "CD68",
  "CD207",
  "CD14",
  # Mast Cells
  "KIT",
  "CPA3",
  # Endothelial Cells
  "VWF",
  "CDH5",
  # Fibroblasts
  "DCN", 
  "COL1A1", 
  "MYL9",
  # Smooth Muscle
  "MYH11",
  "CNN1")

## Colors for axis labels 
gene_colors <- c(rep(rainbow(10)[2], 4), # Lymphocytes
                 rep(rainbow(10)[4], 3), # MNPs
                 rep(rainbow(10)[6], 2), # Mast
                 rep(rainbow(10)[7], 2), # Endothelial
                 rep(rainbow(10)[9], 5)) # Mesenchymal
## Construct dot plot
ggarrange(DotPlot(Whole_Integrated_Obj_NonEpi, 
                  group.by = "seurat_clusters_Annotated", 
                  features = markers ) +
            coord_flip(clip="off") + 
            theme(axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5, size = 10),
                  panel.border = element_rect(color = "black", fill = NA, size = 1),
                  plot.margin = unit(c(2,0,2,0), "lines"),
                  aspect.ratio = 0.7,
                  legend.position = "top",
                  legend.text = element_text(size = 8),
                  legend.title = element_text(size = 9),
                  legend.justification = "center",
                  axis.text.y = element_text(colour=gene_colors, size = 10)) +
            xlab("") + 
            ylab("") +
            
            ### Order clusters according to differentiation markers, with low seq depth clusters at the end
            scale_y_discrete(limits = levels(Whole_Integrated_Obj_Epi$seurat_clusters_Annotated)[c(15,13,18,14, 17,16)]) +
            scale_colour_gradient2(low="#382C84", mid="white", high="#EF0053", 
                                   breaks = c(-1.5,0,1.5)),
          ### add legend for gene category colors
          gg_Legend(df = data.frame("GeneCategory" = c(rep("Lymphocytes", 4),
                                                       rep("MNPs", 3),
                                                       rep("Mast cells", 2),
                                                       rep("Endothelial", 2),
                                                       rep("Mesenchymal", 5)),
                                    "features.plot" = markers),
                    colors = c(
                      rainbow(10)[2],  # Lymphocytes
                      rainbow(10)[4],  # MNPs
                      rainbow(10)[6],  # Mast
                      rainbow(10)[7],  # Endothelial
                      rainbow(10)[9]   # Mesenchymal
                    ),
                    title = "Gene Category",
                    H = 1,
                    V = 0.8),
          ncol = 2, 
          widths = c(4,1),
          heights = c(1,0.3))


# Epithelial clusters had two clusters exhibiting low sequencing depth and lack of biological significance

Whole_Integrated_Obj_Filtered <- Whole_Integrated_Obj[,Whole_Integrated_Obj@meta.data %>%
                                                   filter(seurat_clusters %!in% c(8:9)) %>%
                                                   row.names()]

#Assign Whole Object Clusters

DefaultAssay(Whole_Integrated_Obj_Filtered) <- "RNA"
Whole_Integrated_Obj_Filtered <- NormalizeData(Whole_Integrated_Obj_Filtered)


# Identify/confirm top markers per cluster
Idents(WholeObj_markers) <- "seurat_clusters"
WholeObj_markers <- FindAllMarkers(Whole_Integrated_Obj_Filtered,
                                   assay = "RNA",
                                   slot = "data",
                                   min.pct = 0.25,
                                   only.pos = T)

Markers_perCluster <- as.data.frame(WholeObj_markers) %>%
  mutate(dif_pct = pct.1 - pct.2) %>%
  group_by(cluster) %>%
  filter(pct.1 > 0.6) %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::slice(1:20) %>%
  ungroup() %>%
  arrange(desc(dif_pct)) %>%
  distinct(gene, .keep_all = T) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::slice(1:5) %>%
  ungroup() %>%
  as.data.frame()




markers <- c(# Epithelial Cells
  "KRT6A",
  "DSG3",
  # CD8 T cells
  "CD3D",
  #NK cells
  "NKG7",
  # Macrophages / Dendritic Cells / Monocytes
  "CD68",
  "CD207",
  "CD14",
  # Mast Cells
  "KIT",
  "CPA3",
  #B Cells
  "CD79A",
  "IGHA1",
  # Endothelial Cells
  "VWF",
  "CDH5",
  # Fibroblasts
  "DCN", 
  "COL1A1", 
  "MYL9",
  # Smooth Muscle
  "MYH11",
  "CNN1")


coloring <- rainbow(17)

p1 <- StackedVlnPlot(obj = Whole_Integrated_Obj_Filtered,
                     features = markers[1:9],
                     ident = "seurat_clusters_Annotated",
                     colors = coloring,
                     xlabs = levels(Whole_Integrated_Obj_Filtered$seurat_clusters_Annotated),
                     negNum = 12, 
                     font = 10)

p2 <- StackedVlnPlot(obj = Whole_Integrated_Obj_Filtered,
                     features = markers[10:18],
                     ident = "seurat_clusters_Annotated",
                     colors = coloring,
                     xlabs = levels(Whole_Integrated_Obj_Filtered$seurat_clusters_Annotated),
                     negNum = 12, 
                     font = 10)
p3 <- DimPlot(Whole_Integrated_Obj_Filtered, 
              group.by = "seurat_clusters_Annotated", 
              label = T,
              repel = T) + 
  theme(aspect.ratio = 1) +
  ggtitle("Seurat Clusters")

p1|p2|p3


# Reintegrate epithelial object & annotate -------------------------------------


#Integrate epithelial samples

Epi_Samples <- SplitObject(Whole_Integrated_Obj_Filtered[,Whole_Integrated_Obj_Filtered@meta.data %>% 
                                                  filter(CellType %in% "Epithelial cells") %>% 
                                                  row.names()],
                           split.by = "Patient_Region")

for(x in 1:length(Epi_Samples)){
  Epi_Samples[[x]] <- Seurat::SCTransform(Epi_Samples[[x]], 
                           vst.flavor="v2")
}

features <- SelectIntegrationFeatures(object.list = Epi_Samples, nfeatures = 3000)
Epi_Samples <- PrepSCTIntegration(object.list = Epi_Samples, anchor.features = features)
Anchors <- FindIntegrationAnchors(object.list = Epi_Samples, 
                                  normalization.method = "SCT", 
                                  anchor.features = features, 
                                  reference = which(names(Epi_Samples) %in% c("HCE048_P",
                                                                              "HCE048_D",
                                                                              "HCE049_P",
                                                                              "HCE049_D",
                                                                              "HCE051_P",
                                                                              "HCE051_D")))
Epis_integrated <- IntegrateData(anchorset = Anchors, normalization.method = "SCT")

Epis_integrated <-RunPCA(Epis_integrated, npcs=50) %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>%
  FindClusters(resolution = seq(0,1, by = 0.1))

# Normalize RNA assay

DefaultAssay(Epis_integrated) <- "RNA"
Epis_integrated <- NormalizeData(Epis_integrated)


# Identify top markers per epithelial cluster

Idents(Epis_integrated) <- "integrated_snn_res.0.5"
Epis_markers <- FindAllMarkers(Epis_integrated,
                               assay = "RNA",
                               slot = "data",
                               min.pct = 0.25,
                               only.pos = T)


Markers_perCluster <- as.data.frame(Epis_markers) %>%
  mutate(dif_pct = pct.1 - pct.2) %>%
  group_by(cluster) %>%
  filter(pct.1 > 0.6) %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::slice(1:20) %>%
  ungroup() %>%
  arrange(desc(dif_pct)) %>%
  distinct(gene, .keep_all = T) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  dplyr::slice(1:5) %>%
  ungroup() %>%
  as.data.frame()

# Visualize marker distribution in healthy patients

markers <- c("COL17A1",
             "KRT15", 
             "DST", 
             "CAV1",
             "SOX2",
             "TP63", 
             "KRT14", 
             "PCNA", 
             "MKI67", 
             "KRT5", 
             "KLF5", 
             "KRT13",
             "IVL",
             "SERPINB3",
             "CNFN",
             "KRT78",
             "SPRR2D",
             "FLG")

Idents(Epis_integrated) <- "DiseaseState"
HC_subset <- subset(Epis_integrated, idents = c("Healthy_Control"))

coloring <- c("#E16B63", #B1
              "#DB8E00", #B2
              "#4580F3", #B3
              "#76D100", #B4
              "#00A450", #B5
              "#00C19E", #SB1
              "#04D3FA", #SB2
              "#00A6FF", #SB3
              "#CC5DC9", #SF1
              "#FF63B6") #SF2

DefaultAssay(HC_subset) <- "RNA"
StackedVlnPlot_split.by(obj = HC_subset,
                        features = rev(markers),
                        ident = "integrated_snn_res.0.5",
                        colors =  coloring,
                        xlabs = unique(HC_subset$integrated_snn_res.0.5),
                        negNum = 12,
                        split.by = NULL,
                        ylabs = rev(markers),
                        font = 16) 

# Annotate epithelial clusters & compartments

Epi_Clusters <- c(paste0("Basal_", c(1:5)),
                  paste0("Suprabasal_", c(1:3)),
                  paste0("Superficial_", c(1:2)))

Epis_integrated$Clusters <- plyr::mapvalues(Epis_integrated$integrated_snn_res.0.5, from = c(0,6,7,5,11,3,4,1,10,2), to = Epi_Clusters)
Epis_integrated$Clusters <- factor(Epis_integrated$Clusters, levels = Epi_Clusters)

Epis_integrated$Compartments <- plyr::mapvalues(Epis_integrated$Clusters, from = Epi_Clusters, to = c(rep("Basal", 5),
                                                                                                      rep("Suprabasal", 3),
                                                                                                      rep("Superficial", 2)))
Epis_integrated$Compartments <- factor(Epis_integrated$Compartments, levels = c("Basal", 
                                                                                "Suprabasal", 
                                                                                "Superficial"))

save(Epis_integrated, file = "./Objects/Epis_integrated")
