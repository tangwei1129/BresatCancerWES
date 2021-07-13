rm(list=ls())


### Set up some variables
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(RColorBrewer)
library(scales)
#library(qualV)

## Input variables (default values for running in RStudio)
data_dir = "."
ancestry_data_file=file.path(data_dir, "admixture_table.tsv")

my_prefixes <- c("distance.snps_no_common","distance.snps","distance.snps_ancestry")

lapply(my_prefixes, function(dist_prefix) {
  # dist_prefix="distance.snps_no_common"

  dist_data_file=file.path(data_dir, paste0(dist_prefix,".mdist"))
  dist_ids_file=file.path(data_dir, paste0(dist_prefix,".mdist.id"))
  output_dir <- "figures"
  # output_file <- "mouse_ccbr981/sample_network.admix.pdf"
  mywd="/data/WES_breastcancer/analysis/admixture_plot/" ## working directory (only used if not running from terminal with 'Rscript')
  ## Or read arguments from command line
  # args <- commandArgs(trailingOnly = T)
  # if (length(args) > 0) {
  #   ancestry_data_file <- args[1]
  # } else {
  #   setwd(mywd)
  # }
  # 
  # if (length(args) > 1) {
  #   relatedness_data_file <- args[2]
  # }
  # 
  # if (length(args) > 2) {
  #   output_file <- args[3]
  # }
  # 
  
  if (!dir.exists(output_dir)){dir.create(output_dir, recursive = T)}
  
  ### Read admixture table and make it a matrix
  if (! file.exists(ancestry_data_file)) {
    stop(paste0("Can't find admixture data: ", ancestry_data_file))
  }
  admix_data <- read.table(ancestry_data_file, sep="\t", header=T)
  admix_plotdata <- as.matrix(t(admix_data[,-1]))
  colnames(admix_plotdata) <- admix_data[,1]
  
  library(ggplot2)
  full_admix_names <- c(
    AMR="[AMR] Ad Mixed American",
    AFR="[AFR] African",
    EAS="[EAS] East Asian",
    EUR="[EUR] European",
    SAS="[SAS] South Asian")
  # library(RColorBrewer)
  # admix_colors <- brewer_pal(palette = "Accent")(length(full_admix_names))
  library(wesanderson)
  admix_colors <- wes_palette("Darjeeling2", n=length(full_admix_names))
  names(admix_colors) <- names(full_admix_names)
  
  ### Read sample IDs for distances
  if (! file.exists(dist_ids_file)) {
    stop(paste0("Can't find sample IDs file: ", dist_ids_file))
  }
  sample_ids <- read.table(dist_ids_file, sep="\t", header=F)
  ### Read distance data
  if (! file.exists(dist_data_file)) {
    stop(paste0("Can't find distance data: ", dist_data_file))
  }
  sample_dists <- read.table(dist_data_file, sep=" ", header=F)
  allNAs <- apply(sample_dists, 2, function(x){sum(is.na(x))==length(x)})
  sample_dists <- sample_dists[, ! allNAs]
  
  id_col=which.max(apply(sample_ids, 2, function(x) {sum(x %in% colnames(admix_plotdata))}))
  
  rownames(sample_dists) <- sample_ids[,id_col]
  colnames(sample_dists) <- sample_ids[,id_col]
  
  ### Make dendogram from the distances
  hc <- hclust(as.dist(sample_dists))
  sample_labels = hc$labels[hc$order]  
  n = length(sample_labels)
  dend = as.dendrogram(hc)
  
  admix_legend <- Legend(labels = full_admix_names, type = "points", 
                         labels_gp = gpar(fontsize = 12),
                         pch=22,
                         size=unit(8,"mm"),
                         # legend_gp = gpar(fill = admix_colors, col = admix_colors,pch=22, size=20),
                         legend_gp = gpar(fill = admix_colors, col = "black"),
                         background = "white",
                         # ncol=ifelse(length(admix_colors) > 8, length(admix_colors) %% 8, 1),
                         ncol=1,
                         # title_position = "topleft", 
                         title = "1000 Genomes Super Population")
  
  
  display_labels <- colnames(admix_plotdata)
  ### Match admixture data order to the clustered order
  admix_plotdata <- admix_plotdata[,match(sample_labels, colnames(admix_plotdata))]
  
  
  admix_offset=0.5 ## Adjusts alignment between sample names and admixture plot
  admix_width=0.3 ## Sets width of admix bar
  
  ## Circilze
  ### To do: scale canvas size and/or track heights to get appropriate spacing between the admix barplot and sample labels
  # output_dimension <- length(sample_labels)/20
  output_dimension <- scales::rescale(length(sample_labels), from=c(2, 500), to=c(4, 16))
  output_dimension <- min(c(16, max(c(6, length(sample_labels)))))
  
  
  
  
  output_file <- file.path(output_dir, paste0(dist_prefix,".sample_similarity.pdf"))
  pdf(output_file, width=output_dimension, height=output_dimension)
  circos.par(cell.padding = c(0, 0, 0, 0),
             track.margin = c(0.01, 0.01),
             canvas.xlim = c(-1.5, 1)
             )
  circos.initialize(factors = "single_sector", xlim = c(0, n))
  
  ### Track 1: Admixture bar plot
  circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.2, 
               panel.fun = function(x, y) {
                 for(i in seq_len(n)) {
                   curr_ymin=0
                   for (admix_pop in rownames(admix_plotdata)) {
                     curr_ymax=curr_ymin+admix_plotdata[admix_pop,i]
                     # circos.segments(i-admix_offset,curr_ymin,i-admix_offset, curr_ymax, lwd=6,
                     #                 col=admix_colors[admix_pop])
                     circos.rect(xleft=i-admix_offset-admix_width, 
                                 xright=i-admix_offset+admix_width, 
                                 ybottom=curr_ymin, 
                                 ytop=curr_ymax,
                                 lwd=0.01,
                                 col=admix_colors[admix_pop])
                     curr_ymin=curr_ymax
                   }
                 }
               })
  
  ### Track 2: Sample name labels
  # space_per_char=0.012 * (length(display_labels)/150)
  # space_per_char=0.02
  space_per_char=rescale(max(nchar(display_labels)), from=c(5,40), to=c(0.02, 0.05))
  space_per_char <- min(c(0.05,max(c(0.01, space_per_char))))
  labels_height <- max(nchar(display_labels))*space_per_char
  # labels_height <- rescale(max(nchar(display_labels)), from=c(5,40), to=c(0.08, 0.5))
  labels_height <- min(c(0.4,max(c(0.1, labels_height))))
  # labels_height <- 0.3  ### Increase this to stop overlapping of labels and admixture; decrease to reduce whitespace between them
  
  circos.track(ylim = c(0, 1), bg.border = NA, track.height = labels_height, 
               panel.fun = function(x, y) {
                 for(i in seq_len(n)) {
                   circos.text(i-0.5, 0, display_labels[i], adj = c(0, 0.5),
                               facing = "clockwise", niceFacing = TRUE,
                               col = "grey30", cex = 1)
                 }
               })
  
  ### Track 3: Dendrogram
  dend_height = attr(dend, "height")
  circos.track(ylim = c(0, dend_height), bg.border = NA, 
               track.height = 0.3, panel.fun = function(x, y) {
                 circos.dendrogram(dend)
               })
  circos.clear()
  
  # draw(color_legend,x = unit(-1, "npc"), y = unit(1, "npc"), just = c("left","top"))
  draw(admix_legend,x = unit(0.1, "npc"), y = unit(0.5, "npc"))
  title("Admixture and Pairwise Similarity")
  
  dev.off()
  

})
