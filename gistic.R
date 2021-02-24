setwd("~/Desktop/WES/finalset")
rm(list=ls())
 

#######################################################################################
####  add ccnv data ############
################################

#gistic2 -b $basedir -seg $segfile -refgene $refgenefile -cnv $cnvfile -genegistic 1 -qvt 0.05 -js 20 -ta 0.3 -td 0.3 -smallmem 0 -broad 1 -brlen 0.98 -conf 0.99 -res 0.01 -armpeel 1 -savegene 1 -gcm extreme
# js100
file_dir = 'gistic2_results_0.3.js100.res.0.01.conf0.99'
all.lesions <- list.files(file_dir, pattern = 'all_lesions.conf_99.txt', full.names = T)
amp.genes <- list.files(file_dir, pattern = "amp_genes.conf_99.txt", full.names = T)
del.genes <- list.files(file_dir, pattern = "del_genes.conf_99.txt", full.names =T)
scores.gis <- list.files(file_dir, pattern = "scores.gistic", full.names = T)

gistic = readGistic(gisticAllLesionsFile = all.lesions, gisticAmpGenesFile = amp.genes, gisticDelGenesFile = del.genes, gisticScoresFile = scores.gis)
gistic
getSampleSummary(gistic)
getGeneSummary(gistic)
getCytobandSummary(gistic)
gisticChromPlot(gistic = gistic, markBands = "all", fdrCutOff = 0.05,ref.build = 'hg38',mutGenes = c("MYC","CCND1"))
gisticOncoPlot(gistic = gistic, sortByAnnotation = F, top = 50)
 








 








