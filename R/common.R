plot_ploidy <- function(results, output_dir){
  data <- results@data
  data_cell_stats <- results@cell_stats
  project_name <- results@project
  for(i in 1:length(data)){
  sample_id <- names(data)[i]
  current <- data[[i]]
  cell_stats <- data_cell_stats[which(data_cell_stats$Sample == sample_id),]
  cell_stats$id <- factor(cell_stats$id, levels = levels(current$id))
  print("Printing ploidy by chrom position plot ..")
  p1 <- ggplot(current, aes(range, id, color = Ploidy, fill = Ploidy))+ # [sample(1:nrow(data),size = 100, replace = F),]
    theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),axis.text.y=element_blank(),
          plot.title = element_text(size=25, face = "bold")) +
    # legend.text=element_text(size=18)) + # hjust = 0.5
    # geom_point(alpha = 0.6)+
    geom_tile(size = 0.8)+
    xlab("Chromosome Positions") + ylab("") +
    scale_color_manual(values = results@color_schemes$color_heatmap)+
    scale_fill_manual(values = results@color_schemes$color_heatmap)+
    facet_grid(.~chrom, scales = "free", switch = "x", space = "free_x")+
    ggtitle(paste(sample_id, ", ", signif(length(which(current$Ploidy == "2"))/nrow(current)*100,3), "% diploid"))

  p2 <- ggplot(cell_stats)+
    geom_bar(mapping = aes(x = 1, y = id, fill = Cell_Ploidy), stat = "identity", width = 1)+
    theme_void()+
    theme(panel.spacing.x = unit(1, "mm")) +
    scale_fill_manual(values = results@color_schemes$cell_ploidy_colors)

  legend <- plot_grid(get_legend(p2), get_legend(p1), ncol = 1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- p2 + theme(legend.position = "none")
  p3 <- plot_grid(p2, p1, align = "h", ncol = 2, axis = "tb", rel_widths = c(1,20), rel_heights = c(1,1))
  p <- plot_grid(p3, legend, nrow = 1, rel_widths = c(10,1),rel_heights = c(0.1, 1))

  somePNGPath <- paste(output_dir,"01SCA_PLOT_PLOIDY_INFO_CHROM_POSITION_",sample_id,"_",project_name,".png", sep = "")
  png(somePNGPath, width = 4000, height =3000, units = "px", res = 300)
  print(p)
  dev.off()

  }

}

plot_ploidy_binary <- function(results, output_dir){
  plotx <- results@binary_cnv
  p1 <- ggplot(plotx, aes(CNV_Event, cell_id, fill = Count))+
    theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),axis.text.y=element_blank(),
          plot.title = element_text(size=25, face = "bold"),
          axis.title = element_text(size=20, face = "bold"),
          legend.title = element_text(size =20, face = "bold"),
          legend.text = element_text(size = 15)) +
    geom_tile()+
    xlab("CNV Events") + ylab("") +
    scale_fill_manual(values = c("black","#eec643"))

  p2 <- ggplot(plotx)+
    geom_bar(mapping = aes(x = 1, y = cell_id, fill = Cluster), stat = "identity", width = 1)+
    theme_void()+
    theme(panel.spacing.x = unit(1, "mm"),
          legend.title = element_text(size =20, face = "bold"),
          legend.text = element_text(size = 15)) +
    scale_fill_manual(values = results@color_schemes$cluster_colors)

  legend <- plot_grid(get_legend(p2), get_legend(p1), ncol = 1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- p2 + theme(legend.position = "none")
  p3 <- plot_grid(p2, p1, align = "h", ncol = 2, axis = "tb", rel_widths = c(1,20), rel_heights = c(1,1))

  somePNGPath <- paste(output_dir,"5SCA_BINARY_CNV_EVENTS_BY_CLUSTERS_",results@project,".png", sep = "")
  png(somePNGPath, width = 3000, height =2000, units = "px", res = 300)
  print(plot_grid(p3, legend, nrow = 1, rel_widths = c(10,1),rel_heights = c(0.1, 1)))
  dev.off()

}

plot_selected_chroms <- function(results, output_dir){
  plotx <- results@selected_chrom
  plotx_bar <- results@chrom_bar

  p1 <- ggplot(plotx, aes(range, cell_id, color = Ploidy, fill = Ploidy))+
    theme(axis.ticks.x=element_blank(),axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),axis.text.y=element_blank(),
          plot.title = element_text(size=25, face = "bold"),
          axis.title = element_text(size=20, face = "bold"),
          legend.title = element_text(size =20, face = "bold"),
          legend.text = element_text(size = 15)) +
    geom_tile(size = 0.8)+
    scale_fill_manual(values = results@color_schemes$color_heatmap)+
    scale_color_manual(values = results@color_schemes$color_heatmap)+
    xlab("Chromosome Positions") + ylab("") +
    facet_grid(.~chrom, scales = "free", switch = "x", space = "free_x")

  p2 <- ggplot(plotx_bar)+
    geom_bar(mapping = aes(x = 1, y = cell_id, fill = Cluster), stat = "identity", width = 1)+
    theme_void()+
    theme(panel.spacing.x = unit(1, "mm"),
          legend.title = element_text(size =20, face = "bold"),
          legend.text = element_text(size = 15)) +
    scale_fill_manual(values = results@color_schemes$cluster_colors)

  legend <- plot_grid(get_legend(p2), get_legend(p1), ncol = 1)
  p1 <- p1 + theme(legend.position = "none")
  p2 <- p2 + theme(legend.position = "none")
  p3 <- plot_grid(p2, p1, align = "h", ncol = 2, axis = "tb", rel_widths = c(1,20), rel_heights = c(1,1))

  somePNGPath <- paste(output_dir,"6SCA_PLOIDY_INFO_TOP_50_CNVs_",results@project,".png", sep = "")
  png(somePNGPath, width = 4000, height =3000, units = "px", res = 300)
  print(plot_grid(p3, legend, nrow = 1, rel_widths = c(10,1),rel_heights = c(0.1, 1)))
  dev.off()
}

adjust_theme <- function(p, xangle = 0,legend = "right", title_size = 20, xsize=20, hejust = 0, vejust = 0, strip_size = 20, legend_title = 20, legend_text = 15){
  p <- p+ theme_classic()+
    # ggtitle(title) +
    # scale_fill_manual(values = col) +
    theme(axis.text.x = element_text(size = xsize, angle = xangle, hjust=hejust,vjust = vejust),
          axis.text.y = element_text(size = 20), legend.position = legend,
          axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)),
          legend.title = element_text(size =20, face = "bold"),
          legend.text = element_text(size = 15),
          strip.text.x = element_text(size = strip_size),
          strip.background = element_blank(),
          plot.title = element_text(size =title_size, face = "bold", hjust = 0.5))
  return(p)
}

plot_bygroup <- function(plotx, x, y, group, plot_title, col = NULL, annot = TRUE, legend_position = "right", numeric = FALSE,point_size = 4, label_size = 5, legendsize = 20){

  color_conditions <- color_ini()
  n <- length(unique(plotx[,group]))
  if(is.null(names(col)) | n > length(col)){
    # if(n > 1){
    col <- prepare_cols(col, n)
    # }
  }else{
    if(!is.null(names(col))){
      col <- col[which(names(col) %in% unique(plotx[,group]))]
    }else{
      col <- col[1:length(unique(plotx[,group]))]
    }
  }

  if(annot == TRUE){
    if(numeric == FALSE){
      plotx[,group] <- factor(plotx[,group], levels = sort(unique(as.character(plotx[,group]))))
    }else{
      plotx[,group] <- factor(plotx[,group], levels = sort(unique(as.numeric(as.character(plotx[,group])))))
    }
    centroids <- create_centroids(plotx, x, y, group)
    centroids$Col <- col

    plotx$Label <- ""
    for(i in 1:nrow(centroids)){
      plotx[nrow(plotx)+1,x] <- centroids[i,"Centroid_X"]
      plotx[nrow(plotx),y] <- centroids[i,"Centroid_Y"]
      plotx[nrow(plotx),group] <- centroids[i,"Cluster"]
      plotx[nrow(plotx),"Label"] <- centroids[i,"Cluster"]
    }

    p <- ggplot(plotx, aes(x = plotx[,x], y = plotx[,y], color = plotx[,group], label = Label)) +
      geom_point(alpha = ifelse(plotx$Label != "", 0, 1), size = point_size) +
      scale_color_manual(values = col) +
      theme_classic() +
      xlab(x) + ylab(y) +
      ggtitle(plot_title) +
      theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
            strip.text = element_text(size = 25),
            legend.position = legend_position,
            legend.title = element_blank(),
            legend.text = element_text(size = legendsize),
            legend.key.size = unit(0.8, "cm"),
            axis.text.x = element_text(size = 25),
            axis.text.y = element_text(size = 25),
            axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
            axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)))+
      guides(color=guide_legend(override.aes = list(size = 3), title=toupper(group), ncol=ifelse(length(unique(plotx[,group])) > 20, 2, 1)))

    if(numeric == FALSE){
      p <- p + geom_text_repel(max.overlaps = Inf,force=1,
                               point.padding = 0, # additional padding around each point
                               min.segment.length = 0, # draw all line segments
                               max.time = 1, max.iter = Inf, # stop after 1 second, or after 100,000 iterations
                               box.padding = 0.3, size = label_size, colour = "black")
      # p <- p + annotate("text", x=centroids$Centroid_X, y=centroids$Centroid_Y, label= centroids$Cluster, size = label_size, hjust = 0, fontface =1)
    }else{
      p <- p + annotate("text", x=centroids$Centroid_X, y=centroids$Centroid_Y, label= centroids$Cluster, size = label_size, hjust = 0, fontface =1)
    }
  }else{
    p <- ggplot(plotx, aes(x = plotx[,x], y = plotx[,y], color = plotx[,group])) +
      geom_point(alpha = 1, size = point_size) +
      scale_color_manual(values = col) +
      theme_classic() +
      xlab(x) + ylab(y) +
      ggtitle(plot_title) +
      theme(plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
            strip.text = element_text(size = 25),
            legend.position = legend_position,
            legend.title = element_blank(),
            legend.text = element_text(size = legendsize),
            legend.key.size = unit(1, "cm"),
            axis.text.x = element_text(size = 25),
            axis.text.y = element_text(size = 25),
            axis.title.x = element_text(size = 25, margin=margin(10,0,0,0)),
            axis.title.y = element_text(size = 25, margin=margin(0,10,0,0)))+
      guides(color=guide_legend(override.aes = list(size = 3), title=toupper(group), ncol=ifelse(length(unique(plotx[,group])) > 20, 2, 1)))
  }

  return(p)

}

create_centroids <- function(plotx, x, y, group){

  centroids <- split(plotx, plotx[,group])
  # centroid_groups <- names(centroids)
  centroids <- lapply(centroids, function(cent){
    cent <- data.frame(Cluster = ifelse(nrow(cent) > 0,as.character(unique(cent[,group])), NA),
                       Centroid_X = ifelse(nrow(cent) > 0, median(cent[,x]), NA),
                       Centroid_Y = ifelse(nrow(cent) > 0, median(cent[,y]), NA))})
  centroids <- do.call(rbind.data.frame, centroids)
  # centroids[,group] <- centroid_groups
  centroids <- centroids[!is.na(centroids[,"Cluster"]),]
  return(centroids)

}

complex_heatmap <- function(x, legend_title = NULL, col_title = NULL, row_title = NULL,
                            col = rev(rainbow(10)),nrow_clus = 1, row_annot = NULL, legendtitle = "EXPRESSION"){
  library("ComplexHeatmap")
  p <- Heatmap(x, name = legend_title, col = col, row_km = nrow_clus, right_annotation = row_annot,
               row_names_gp = gpar(fontsize = 8), row_names_side = "left",
               column_title = col_title,column_names_rot = 45,
               heatmap_legend_param = list(title = legendtitle))
  return(p)
}

violin_plot <- function(current, features, ncol = NULL, col, x_lab, log_status = TRUE){
  p <- VlnPlot(object = current, features = features, ncol = ncol, cols = col, pt.size = 0.05, log = log_status)  & xlab(x_lab)
  return(p)
}

own_facet_scatter <- function(plotx, feature1, feature2, isfacet = T,
                              title, col=color_conditions$bright, color_by,
                              group_by = NULL, xlabel = feature1, ylabel = feature2,
                              strip_size = 15, legend_pos = "right", ncol = 2){
  # set.seed(2022)
  p <- ggplot(plotx, aes(x = plotx[,feature1], y = plotx[,feature2], color = plotx[,color_by])) +
    geom_point(size = 0.9) +
    scale_color_manual(values = gen_colors(col, length(unique(plotx[,color_by])))) +
    # scale_color_manual(values = sample(gen_colors(col, length(unique(plotx[,color_by]))), size = length(unique(plotx[,color_by])), replace = F)) +
    theme_classic() +
    xlab(xlabel)+ylab(ylabel)+
    ggtitle(title) +
    guides(colour = guide_legend(title = color_by, override.aes = list(size=5)))+
    theme(legend.position = legend_pos,
          axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 20),
          axis.title.x = element_text(size = 20, margin=margin(10,0,0,0)),
          axis.title.y = element_text(size = 20, margin=margin(0,10,0,0)),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.8, "cm"),
          plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
          strip.text = element_text(size = strip_size, face = "bold"),
          strip.background = element_blank())
  if(isfacet == T & !is.null(group_by)){
    p <- p + facet_wrap(~plotx[,group_by], ncol = ncol)
  }
  return(p)
}

add_names <- function(current, sample_name,current_ident){
  sample_name <- gsub("_FILTERED_PEAK_BC|_FILTERED_FEATURE_BC","", sample_name, ignore.case = T)
  # sample_name <- ifelse(nchar(sample_name) > 8, gsub(".*\\.(.*)","\\1",sample_name), sample_name)
  current <- SetIdent(current, value = current_ident)
  levels(current@active.ident) <- sample_name
  current@project.name <- sample_name
  current$orig.ident <- sample_name
  return(current)
}

plot_slice <- function(plotx, plot_title, col = c(NULL,"default"), is.facet = FALSE, annot = TRUE, strip_size = 15, pt_size = 2){

  color_conditions <- color_ini()
  n <- length(unique(plotx$Cluster))
  col <- prepare_cols(col, n)

  # plotx$Cluster <- factor(plotx$Cluster, levels = sort(unique(as.character(plotx$Cluster))))

  if(annot == TRUE){
    centroids <- create_centroids(plotx, "X","Y", "Cluster")
    centroids$Col <- col
  }

  if(is.facet == TRUE){
    p <- ggplot(plotx, aes(X, Y, color = Cluster)) +
      # geom_point(shape = 15, aes(colour=Cluster)) +
      # geom_point(shape = 0,colour = "black", stroke = 0.05)+
      geom_point(shape = 16, size = pt_size) +
      coord_flip() + scale_x_reverse() + theme(legend.position = "right") +
      scale_color_manual(values = col) + facet_wrap(~Cluster) +
      theme_map() + guides(color = FALSE, size = FALSE) + ggtitle("") +
      theme(strip.text.x = element_text(face="bold", size = strip_size),
            strip.background = element_rect(colour="black", fill="white"))
  }else{
    p <- ggplot(plotx, aes(X, Y, color = Cluster)) +
      # geom_point(shape = 20, colour = "black", size = 1, stroke = 0.5) +
      geom_point(shape = 16, size = pt_size) +
      coord_flip() + scale_x_reverse() + theme(legend.position = "right") +
      scale_color_manual(values = col) +
      theme_map() + guides(color=guide_legend(title="Clusters", ncol=ifelse(length(unique(plotx$Cluster)) > 8, 2, 1))) +
      ggtitle(plot_title) + theme(plot.title = element_text(size=15, face = "bold", hjust = 0.5), legend.position = "right")
  }

  if(annot == TRUE){
    p <- p + guides(color = FALSE, size = FALSE) +
      annotate(geom="label", x=centroids$Centroid_X, y=centroids$Centroid_Y, label= centroids$Cluster, size = 3, fill=centroids$Col, col = "white")
  }
  return(p)

}

