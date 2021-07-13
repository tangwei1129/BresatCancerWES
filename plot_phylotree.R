

### Set up some variables
library(ComplexHeatmap)
library(circlize)
library(dendextend)
library(RColorBrewer)
library(scales)
library(ape)



sample_dists <- read.table('plink.mdist', sep="\t", header=F)
sample_ids = read.table('plink.mdist.id' ,sep='\t')
graf=read.csv('GRAF_pop.result.csv',header=T,row.names = 1)
colnames(sample_dists) = sample_ids$V1
rownames(sample_dists) = sample_ids$V1


### construct cluster ###
hc <- hclust(as.dist(sample_dists))
labels = hc$labels[hc$order]  
n = length(labels)
dend = as.dendrogram(hc)

#### dend label text color ####
ct= graf$Self.reported.ancestry
ct= as.factor(ct)
levels(ct) = c("#046C9A","#D69C4E", "#ECCBAE")  # color for text match to bar plot
ct=as.character(ct)
names(ct) = rownames(graf)


### number of samples 
n = length(labels)  # number of samples


#### construct admixture plot from GRAF #####
### Read admixture table and make it a matrix
 
admix_data <- graf[,c('P_f....', 'P_e....', 'P_a....')]
admix_plotdata <- as.matrix(t(admix_data))/100
rownames(admix_plotdata) <- c('African','European','Asian')

library(ggplot2)
full_admix_names <- c(
  European = "European",
  African = "African", 
  Asian = "Asian")

# library(RColorBrewer)
# admix_colors <- brewer_pal(palette = "Accent")(length(full_admix_names))
library(wesanderson)
admix_colors <- wes_palette("Darjeeling2", n=length(full_admix_names))
names(admix_colors) <- names(full_admix_names)

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
                       title = "Population")


display_labels <- colnames(admix_plotdata)
### Match admixture data order to the clustered order
admix_plotdata <- admix_plotdata[,match(labels, colnames(admix_plotdata))]


admix_offset=0.5 ## Adjusts alignment between sample names and admixture plot
admix_width=0.3 ## Sets width of admix bar





## Circilze
### To do: scale canvas size and/or track heights to get appropriate spacing between the admix barplot and sample labels
# output_dimension <- length(sample_labels)/20
output_dimension <- scales::rescale(length(labels), from=c(2, 500), to=c(4, 16))
output_dimension <- min(c(16, max(c(6, length(labels)))))



pdf('adm_dis.pdf', width=output_dimension, height=output_dimension)
circos.par(cell.padding = c(0, 0, 0, 0),
           track.margin = c(0.01, 0.01),
           canvas.xlim = c(-1.5, 1)
          )


circos.initialize("a", xlim = c(0, n)) # only one sector


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

 
### track 2 dend sample name ####
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.15, 
             panel.fun = function(x, y) {
               for(i in seq_len(n)) {
                 circos.text(i-0.5, 0, labels[i], adj = c(0, 0.5), 
                             facing = "clockwise", niceFacing = TRUE,
                             col = ct[labels[i]], cex = 0.5)
               }
             })

 
### track 3 dend plot ####
suppressPackageStartupMessages(library(dendextend))
dend = color_branches(dend, k = 3, col = c("#046C9A","#D69C4E", "#ECCBAE"))
#dend = color_branches(dend, col = ct)
dend_height = attr(dend, "height")
circos.track(ylim = c(0, dend_height), bg.border = NA, 
             track.height = 0.3, panel.fun = function(x, y) {
               circos.dendrogram(dend)
             })
circos.clear()

# draw(color_legend,x = unit(-1, "npc"), y = unit(1, "npc"), just = c("left","top"))
draw(admix_legend,x = unit(0.3, "npc"), y = unit(0.8, "npc"))
title("Admixture and Pairwise Similarity")


dev.off()


