############################################################################################################
# SingleCellAnalyst.org
# Phecda: Spatial
# Version: V1.0.0
# Creator: Lu Pan, Karolinska Institutet, lu.pan@ki.se
# Date: 2022-02-16
############################################################################################################
#' @importFrom Seurat
#' @importFrom ggplot2
#' @importFrom SingleR
#' @importFrom gridExtra
#' @importFrom
#' @importFrom
#' @importFrom
#' @importFrom
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Run a post-quantification 10X technology-format Spatial Transcriptomics pipeline
#'
#' This function will run a bioinformatics analysis of post-quantification Spatial
#' pipeline. Supports multiple samples analysis.
#'
#' @param project_name Project name. 'Polaris_Spatial' by default.
#' @param input_dir Directory to all input files. Current working directory by
#' default.
#' @param output_dir Output directory. Current working directory by default.
#' A new folder with the given project name with time stamp as suffix will be
#' created under the specified output directory.
#' @param pheno_file Meta data file directory. Accept only .csv/.txt format
#' files.
#' @export
#'
SpatialPip <- function(project_name = "Polaris_Spatial",
                     input_dir = "./",
                     output_dir = "./",
                     pheno_file){
  print("Initialising pipeline environment..")
  pheno_data <- pheno_ini(pheno_file, pipeline = "SPATIAL", isDir = T)
  color_conditions <- color_ini()
  ctime <- time_ini()
  hpca.se <- readRDS("DB/hpca.se.RDS")
  hs <- org.Hs.eg.db
  hgnc.table <- readRDS("DB/hgnc.table.RDS")
  p_val_adj <- 0.1
  project_name <- gsub("\\s+|\\(|\\)|-|\\/|\\?","",project_name)
  print(paste("Creating output folder ",project_name,"_",ctime," ..", sep = ""))
  cdir <- paste(output_dir,project_name,"_",ctime,"/", sep = "")
  dir.create(cdir)

  for(i in 1:nrow(pheno_data)){
    if(length(grep("\\.csv$|\\.txt$", pheno_data$FILE[i], ignore.case = T)) == 0){
      if(length(grep("\\.h5$", pheno_data$FILE[i], ignore.case = T)) == 0){
        pheno_data$FILE[i] <- paste(pheno_data$FILE[i],".h5", sep = "")
      }
    }
  }

  sample_files <- list.files(input_dir, recursive = T, full.names = T)
  sample_files <- sample_files[grep("_MACOSX",sample_files, ignore.case = T, invert = T)]
  sample_files <- sample_files[grep(gsub(".*\\/(.*)","\\1",pheno_file, ignore.case = T),sample_files, ignore.case = T, invert = T)]
  sample_files <- sample_files[grep("\\.h5$", sample_files, ignore.case = T)]

  data <- NULL
  sample_features <- NULL
  data_current <- NULL
  data_markers <- NULL
  plotx <- NULL
  results <- NULL
  annot_names <- NULL
  scrna_data <- NULL

  for(i in 1:length(sample_files)){
    current_sample <- gsub(".*\\/(.*\\.h5)$","\\1",sample_files[i], ignore.case = T)
    cpheno <- pheno_data[which(toupper(pheno_data$FILE) == toupper(current_sample)),]
    cname <- cpheno$SAMPLE_ID
    print(current_sample)
    print(paste("Running sample: ", current_sample, "..", sep = ""))
    current <- Load10X_Spatial(data.dir = gsub("(.*\\/).*","\\1",sample_files[i], ignore.case = T), filename = current_sample)
    current <- SCTransform(current, assay = "Spatial", verbose = FALSE)
    current@project.name <- project_name
    current$orig.ident <- cname
    names(current@images) <- cname
    # current@images[[1]]@key <- cname
    # current@images[[1]]@key <- paste(gsub("\\.","",cname),"_", sep = "")
    colnames(current@meta.data)[grep("nCount_Spatial", colnames(current@meta.data), ignore.case = T)] <- "Spatial_Counts"
    colnames(current@meta.data)[grep("nFeature_Spatial", colnames(current@meta.data), ignore.case = T)] <- "Feature_Counts"
    levels(current@active.ident) <- cname
    current <- add_names(current, cname, cname)
    current$CELL_ID <- row.names(current@meta.data)
    print(current)

    if(i==1){
      data <- current
      DefaultAssay(data) <- "SCT"
      VariableFeatures(data) <- VariableFeatures(data)
      sample_features[[i]] <- VariableFeatures(data)
      names(sample_features)[i] <- cname
    }else{
      data <- merge(data, current)
      DefaultAssay(data) <- "SCT"
      VariableFeatures(data) <- c(VariableFeatures(data),VariableFeatures(current))
      sample_features[[i]] <- VariableFeatures(current)
      names(sample_features)[i] <- cname
    }
    data_current[[i]] <- current
    names(data_current)[i] <- cname
    annot_names <- c(annot_names, cname)

  }

  print(paste("Running dimension reduction and clustering..", sep = ""))
  VariableFeatures(data) <- unlist(lapply(data_current, VariableFeatures))
  DefaultAssay(data) <- "SCT"
  data <- RunPCA(data, assay = "SCT", verbose = FALSE)
  data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
  data <- FindClusters(data, verbose = FALSE)
  data <- RunUMAP(data, reduction = "pca", dims = 1:30)

  for(i in 1:length(annot_names)){
    current <- subset(data, subset = orig.ident == annot_names[i])
    current@images <- current@images[names(current@images) == annot_names[i]]
    current_clusters <- as.numeric(as.character(unique(current@meta.data$seurat_clusters)))
    Idents(current) <- "seurat_clusters"
    DefaultAssay(current) <- "Spatial"
    current_de_markers <- NULL
    current_de_markers <- FindAllMarkers(current, min.pct = 0.25, logfc.threshold = 0.25)
    current_de_markers <- data.frame(SAMPLE = annot_names[i], current_de_markers)
    current_de_markers <- current_de_markers[which(current_de_markers$p_val_adj < 0.1),]
    current_de_markers <- current_de_markers[order(current_de_markers$avg_log2FC, decreasing = T),]
    data_markers[[i]] <- current_de_markers
    names(data_markers)[i] <- annot_names[i]

    Idents(current) <- "seurat_clusters"
    DefaultAssay(current) <- "Spatial"
    clu_ann <- SingleR(test = as.SingleCellExperiment(DietSeurat(current)),
                       clusters =  current$seurat_clusters,
                       ref = hpca.se, assay.type.test=1,
                       labels = hpca.se$label.main)
    current$Cell_Type <- clu_ann$labels[match(current$seurat_clusters,row.names(clu_ann))]
    current@meta.data[which(is.na(current$Cell_Type)),"Cell_Type"] <- "Unidentifiable"
    data@meta.data[which(data$orig.ident == annot_names[i]),"Cell_Type"] <- current$Cell_Type[match(row.names(data@meta.data[which(data$orig.ident == annot_names[i]),]), row.names(current@meta.data))]
    Idents(current) <- current$Cell_Type
    data_current[[which(names(data_current) %in% annot_names[i])]]$seurat_clusters <- current@meta.data[match(current$CELL_ID, data_current[[which(names(data_current) %in% annot_names[i])]]$CELL_ID),"seurat_clusters"]
  }

  cluster_colors <- gen_colors(color_conditions$tenx, length(unique(data$seurat_clusters)))
  names(cluster_colors) <- sort(unique(data$seurat_clusters), decreasing = F)
  ct_colors <- gen_colors(color_conditions$colorful, length(unique(data$Cell_Type)))
  names(ct_colors) <- sort(unique(data$Cell_Type))
  sample_colors <- gen_colors(color_conditions$bright, length(unique(data$orig.ident)))
  names(sample_colors) <- sort(unique(data$orig.ident), decreasing = F)

samples <- names(data_current)
for(i in 1:length(data_current)){
plotx <- data_current[[i]]
plotx$orig.ident <- ifelse(nchar(plotx$orig.ident) > 15, substr(plotx$orig.ident, 1, 15), plotx$orig.ident)
p1 <- NULL
p2 <- NULL
  p1 <- VlnPlot(plotx, group.by = "orig.ident", features = "Spatial_Counts", pt.size = 0.1, cols = color_conditions$general) + NoLegend() + xlab("SAMPLE_NAME")
  p2 <- SpatialFeaturePlot(plotx, features = "Spatial_Counts") + theme(legend.position = "right")

  somePNGPath <- paste(cdir,"1SCA_VIOLIN_IMAGE_FEATURES_",samples[i],"_",project_name, ".png", sep = "")
  png(somePNGPath, width = 4000, height =3000, units = "px", res = 300)
  print((p1|p2)+plot_annotation(title = gsub("\\."," - ",samples[i]), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
  dev.off()

  cfeatures <- sample_features[[i]]
  DefaultAssay(plotx) <- "SCT"
    p1 <- SpatialFeaturePlot(plotx, features = head(cfeatures, 5), ncol = 5)
    p2 <- SpatialFeaturePlot(plotx, features = head(cfeatures, 5), alpha = c(0.1, 1), ncol = 5)
    somePNGPath <- paste(cdir,"2SCA_SPATIALLY_VARIABLE_FEATURES_",samples[i],"_",project_name,".png", sep = "")
    png(somePNGPath, width = 4000, height =3000, units = "px", res = 300)
    print(p1/p2)
    dev.off()
}

plotx <- data.frame(UMAP_1 = data@reductions$umap@cell.embeddings[,"UMAP_1"],
                    UMAP_2 = data@reductions$umap@cell.embeddings[,"UMAP_2"],
                    CLUSTERS = data@meta.data$seurat_clusters,
                    SLIDES = data$orig.ident)
p1 <- NULL
p2 <- NULL
p1 <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CLUSTERS", point_size = 2, label_size = 6, numeric = T, annot = T,
                   plot_title = paste(project_name, ": UMAP CLUSTERS", sep = ""), col = cluster_colors)
p1 <- adjust_theme(p1)
p2 <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "SLIDES", point_size = 2,
                   plot_title = paste(project_name, ": SLIDES", sep = ""), col = sample_colors, annot = FALSE)
p2 <- p2+theme(legend.text = element_text(size = 10))
p2 <- adjust_theme(p2)

somePNGPath <- paste(cdir,"3SCA_UMAP_ALL_SAMPLES_",project_name, ".png", sep = "")
png(somePNGPath, width = 5000, height =2000, units = "px", res = 150)
print(p1|p2)
dev.off()

p1 <- own_facet_scatter(plotx, "UMAP_1", "UMAP_2", isfacet = T, title = paste(project_name, ": UMAP CLUSTERS", sep = ""),color_by = "CLUSTERS", group_by = "SLIDES", xlabel = "UMAP_1", ylabel = "UMAP_2", strip_size = 10, col = cluster_colors, ncol = 2)
p2 <- own_facet_scatter(plotx, "UMAP_1", "UMAP_2", isfacet = T, title = paste(project_name, ": SLIDES", sep = ""),color_by = "SLIDES", group_by = "SLIDES", xlabel = "UMAP_1", ylabel = "UMAP_2", strip_size = 10, col = sample_colors, ncol = 2)
p2 <- p2+theme(legend.text = element_text(size = 10))

somePNGPath <- paste(cdir,"4SCA_UMAP_ALL_SLIDES_",project_name, ".png", sep = "")
png(somePNGPath, width = 5000, height = ceiling(length(unique(plotx$SLIDES))/2)*1000, units = "px", res = 200)
grid.arrange(p1,p2,ncol=2,widths=c(0.45, 0.55))
dev.off()

for(i in 1:length(samples)){
p <- NULL
p <- SpatialDimPlot(data,ncol = 1,pt.size.factor = 1.6,
                    images = samples[i], cols = cluster_colors)+
  ggtitle(samples[i]) +
  labs(fill= "CLUSTERS")+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
        axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)))+
  guides(fill=guide_legend(override.aes = list(size = 5)))

somePNGPath <- paste(cdir,"5SCA_UMAP_SLIDE_CLUSTERS_",samples[i],"_",project_name,".png", sep = "")
png(somePNGPath, width = 4000, height =3000, units = "px", res = 200)
print(p)
dev.off()

current <- data_current[[i]]
Idents(current) <- "seurat_clusters"
DefaultAssay(current) <- "SCT"
col_left <- as.character(unique(current$seurat_clusters))[match(levels(data$seurat_clusters), as.character(unique(current$seurat_clusters)))]
current@images <- current@images[names(current@images) == samples[i]]
somePNGPath <- paste(cdir,"6SCA_SLIDE_IMAGE_BY_EACH_CLUSTER_",samples[i],"_",project_name, ".png", sep = "")
png(somePNGPath, width = 4000, height =2800, units = "px", res = 200)
print(SpatialDimPlot(current, cells.highlight = CellsByIdentities(object = current,
                                                                  idents = as.numeric(as.character(col_left[!is.na(col_left)]))), facet.highlight = TRUE, ncol = 5))
dev.off()

plotx <- subset(data, subset = orig.ident == samples[i])
Idents(plotx) <- "seurat_clusters"
plotx <- data.frame(UMAP_1 = plotx@reductions$umap@cell.embeddings[,"UMAP_1"],
                    UMAP_2 = plotx@reductions$umap@cell.embeddings[,"UMAP_2"],
                    CLUSTERS = plotx$seurat_clusters)

somePNGPath <- paste(cdir,"7SCA_UMAP_SAMPLE_CLUSTER_",samples[i],"_",project_name, ".png", sep = "")
png(somePNGPath, width = 4000, height = 3000, units = "px", res = 300)
print(plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CLUSTERS", plot_title = paste(samples[i], ": UMAP CLUSTERS", sep = ""), annot = T, label_size = 8, point_size = 2, numeric = T, col = data()[['cluster_colors']]))
dev.off()

plotx <- subset(data, subset = orig.ident == samples[i])
Idents(plotx) <- "seurat_clusters"
plotx <- data.frame(Cluster = plotx$seurat_clusters,
                    X = plotx@images[[samples[i]]]@coordinates$imagerow,
                    Y = plotx@images[[samples[i]]]@coordinates$imagecol)

col_left <- as.character(unique(plotx$Cluster))[match(levels(data@meta.data$seurat_clusters), as.character(unique(plotx$Cluster)))]

p1 <- NULL
p2 <- NULL
p1 <- plot_slice(plotx, pt_size = 3, plot_title = "SLICE CLUSTER VIEW", col = cluster_colors[which(!is.na(col_left))], is.facet = FALSE, annot = FALSE)
p1 <- p1+theme(legend.title = element_text(size = 15),
               legend.text = element_text(size = 15),
               legend.key.size = unit(0.2, "cm"))
p2 <- plot_slice(plotx, pt_size = 1.5, plot_title = "SLICE CLUSTER VIEW", col = cluster_colors[which(!is.na(col_left))], is.facet = TRUE, annot = FALSE, strip_size = 15)
somePNGPath <- paste(cdir,"8SCA_UMAP_SAMPLE_IMAGE_CLUSTER_",samples[i],"_",project_name, ".png", sep = "")
png(somePNGPath, width = 4000, height =2000, units = "px", res = 300)
print((p1|p2)+plot_annotation(title = samples[i], theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))))
dev.off()

current <- subset(data, subset = orig.ident == samples[i])
current_de <- data_markers[[samples[i]]]
cclusters <- as.character(unique(current_de$cluster))
for(j in 1:length(cclusters)){
ccurrent_de <- current_de[which(as.character(current_de$cluster) == as.character(cclusters[j])),]
if(nrow(ccurrent_de) > 0){
  ccurrent_de <- ccurrent_de[1:ifelse(nrow(ccurrent_de) < 5, nrow(ccurrent_de), 5),]
  p <- NULL
  p <- SpatialFeaturePlot(object = current, features = ccurrent_de$gene, images = samples[i], ncol = ifelse(nrow(ccurrent_de) < 5, nrow(ccurrent_de), 5))
  somePNGPath <- paste(cdir,"9SCA_SPATIALLY_TOP_FEATURES_",samples[i],"_CLUSTER_",cclusters[j],"_",project_name, ".png", sep = "")
  png(somePNGPath, width = 4000, height =2800, units = "px", res = 200)
  print(p)
  dev.off()
}
}

write.csv(data_markers[[i]], paste(cdir, "10SCA_TABLE_DATA_CLUSTER_DEG_",samples[i],".csv",sep = ""),row.names = F)

current_de <- current_de[current_de$p_val_adj < 0.1,]
current_de <- current_de[order(current_de$avg_log2FC, decreasing = T),]
current_de <- split(current_de, current_de$cluster)
top_markers <- NULL
for(j in 1:length(current_de)){
  if(nrow(current_de[[j]]) > 0){
    top_markers <- rbind(top_markers,current_de[[j]][1,])
  }
}

DefaultAssay(current) <- "Spatial"
p <- NULL
p <- VlnPlot(current, group.by = "seurat_clusters", features = unique(top_markers$gene),
             pt.size = 0, ncol = 2, cols = cluster_colors, log = T)
p <- p+plot_annotation(title = paste("TOP GENE IN EACH CLUSTER\nLog(Average Expression)", sep = ""), theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5)))

somePNGPath <- paste(cdir,"11SCA_VIOLINPLOT_TOP_MARKERS_",samples[i],"_",project_name, ".png", sep = "")
png(somePNGPath, width = 4000, height =400*ceiling(length(unique(top_markers$gene))/2), units = "px", res = 200)
print(p&xlab("CLUSTERS"))
dev.off()

Idents(current) <- current$Cell_Type
plotx <- data.frame(UMAP_1 = current@reductions$umap@cell.embeddings[,"UMAP_1"],
                    UMAP_2 = current@reductions$umap@cell.embeddings[,"UMAP_2"],
                    CELL_TYPE = current$Cell_Type)
p1 <- NULL
p2 <- NULL
p <- NULL
p1 <- plot_bygroup(plotx, x = "UMAP_1", y = "UMAP_2", group = "CELL_TYPE", plot_title = "SPATIAL UMAP - CELL TYPE",
                   col = ct_colors, annot = TRUE, legend_position = "right", point_size = 2)
p2 <- SpatialDimPlot(current,ncol = 1,pt.size.factor = 1.6,images = samples[i],
                     cols = ct_colors[which(names(ct_colors) %in% unique(current$Cell_Type))])+
  ggtitle("SPATIAL SLIDE - CELL TYPE") +
  labs(fill= "CELL TYPES")+
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, "cm"))+
  guides(fill=guide_legend(override.aes = list(size = 5)))
p <- (p1/p2)+plot_annotation(title = paste(project_name,": ",samples[i],sep = ""), theme = theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)))


somePNGPath <- paste(cdir,"12SCA_UMAP_AUTO_CELLTYPE_ANNOTATIONS_",samples[i],"_",project_name,".png", sep = "")
png(somePNGPath, width = 4000, height =4000, units = "px", res = 300)
print(p)
dev.off()

}






}
