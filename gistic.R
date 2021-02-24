setwd("~/Desktop/WES/finalset")
rm(list=ls())


library(maftools)
library(TCGAmutations)
TCGAmutations::tcga_available()
tcga_brca_mc3 = TCGAmutations::tcga_load(study = "BRCA",source="MC3") #tcga_brca_mc3
tcga_brca = TCGAmutations::tcga_load(study = "BRCA",source="Firehose")#tcga_brca
source('functions.R')


## check the titv for tcga dataset ##
#brca.titv = titv(maf = tcga_brca, plot = T, useSyn = TRUE)
#dev.new()
#brca.mc3.titv = titv(maf = tcga_brca_mc3, plot = T, useSyn = TRUE)


########################################################################################
############## data input
mutect2 = read.csv("../mutect2/keep.sets.cl.csv",header = T, row.names = 1)
mutect2.cgc = read.csv("../mutect2/CGC.mutation.csv", header = T, row.names = 1)
  

merge = read.csv("../merged/keep.sets.cl.csv",header = T, row.names = 1)
merge.cgc = read.csv("../merged/CGC.mutation.csv", header = T, row.names = 1)

  
library(VennDiagram)
venn.diagram(
  x = list(make.venn.list(mutect2), make.venn.list(mutect2.cgc), make.venn.list(merge),make.venn.list(merge.cgc)),
  category.names = c("mutect2" , "mutect2.cgc " , "merge","merge.cgc"),
  filename = 'SNVs.png',
  output=T
)

# missing important mutation from mutect2 calls but in multi call##
#GATA3.chr10.8074000.8074000.Sample_9520
#PIK3CA.chr3.179234286.179234286.Sample_24028
#PIK3CA.chr3.179234287.179234287.Sample_24028
#PTEN.chr10.87952155.87952155.Sample_20
#TP53.chr17.7674191.7674191.Sample_10669
###

### create two datasets, one mutect2, mutect2+merge+merge.cgc ####
mutect2.temp = rbind(mutect2[,1:45])
### remove FLAGS genes FLAGS - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267152/ ###
## remove flags 20 genes and Unknown ##
mutect2.final = subset(mutect2.temp, !Hugo_Symbol %in% c(flags(20),"Unknown"))
write.csv(mutect2.final,"mutect2.final.csv")

merge.temp = rbind(mutect2[,1:45],merge[,1:45],merge.cgc[,1:45])
merge.temp = merge.temp[!duplicated(merge.temp),]
### remove FLAGS genes FLAGS - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267152/ ###
## remove flags 20 genes and Unknown ##
merge.final = subset(merge.temp, !Hugo_Symbol %in% c(flags(20),"Unknown"))
write.csv(merge.final,"merge.final.csv")
###########

### read in maf ###
mutect2.maf = read.maf(mutect2.final, clinicalData = NULL, removeDuplicatedVariants = T)
merge.maf = read.maf(merge.final, clinicalData = NULL, removeDuplicatedVariants = T)
cgc = read.maf(merge.cgc)

#### write out for mutsig ####
#### go to /data/tangw3/MutSig
### ml MutSig MutSigCV **.maf $MUTSIG_REF/exome_full192.coverage.txt $MUTSIG_REF/gene.covariates.txt output $MUTSIG_REF/mutation_type_dictionary_file.txt $MUTSIG_REF/chr_files_hg38
###############

get_mutsig_file(mutect2)
get_mutsig_file(merge.temp)

prepareMutSig(maf = tcga_brca_mc3 , fn = "mc3")

### make oncoplot ##########################
### add mutsig results with filter 0.01  ###
### add titv results #######################

## mutect2
mutect2.mutsig = read.table('output.mutect2.mutSig.sig_genes.txt', header = T, sep="\t")[,c("gene","p")]
mutect2.mutsig$p = -log10(mutect2.mutsig$p + 10e-16)-2

oncoplot(maf = mutect2.maf, draw_titv = TRUE, top = 50, 
         leftBarData = mutect2.mutsig, leftBarLims = c(0, 16-2))
get.ti.tv(mutect2.maf)

## merge
merge.mutsig = read.table('output.merge.temp.mutSig.sig_genes.txt', header = T, sep="\t")[,c("gene","p")]
merge.mutsig$p = -log10(merge.mutsig$p + 10e-16)-2

oncoplot(maf = merge.maf, draw_titv = TRUE, top = 50, 
         leftBarData = merge.mutsig, leftBarLims = c(0, 16-2))
get.ti.tv(merge.maf)

## tcga brca mc3
mc3.mutsig = read.table('output.mc3.mutSig.sig_genes.txt', header = T, sep="\t")[,c("gene","q")]
mc3.mutsig$q = -log10(mc3.mutsig$q + 10e-16)

oncoplot(maf = tcga_brca_mc3, draw_titv = TRUE, top = 50, 
         leftBarData = mc3.mutsig, leftBarLims = c(0, 16))
get.ti.tv(tcga_brca_mc3)

###### plot mutation burden across TCGA cancer type  ######
par(mfrow=c(3,1))
tcgaCompare(maf = mutect2.maf, cohortName = 'LHC_mutect2', logscale = TRUE, capture_size = 50)
tcgaCompare(maf = merge.maf, cohortName = 'LHC_merge', logscale = TRUE, capture_size = 50)
tcgaCompare(maf = tcga_brca_mc3, cohortName = 'BRCA_MC3', logscale = TRUE, capture_size = 50)

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
 











#######################################################################################################
 
#test=subset(merge.only3caller)
final.set.maf=read.maf(final.set)
oncoplot(final.set.maf,top=50)
titv = titv(maf =merge.maf, plot = T, useSyn = TRUE)
plotTiTv(res =  titv)
 

 




 
library(maftools) 
### summary plots ###
laml = read.maf(maf ="mutect2_merged_cl.maf",clinicalData = "info.tsv")
plotmafSummary(maf = laml, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
write.mafSummary(maf = laml, basename = 'laml')

### oncoplot ###
#Changing colors for variant classifications   
col = RColorBrewer::brewer.pal(n = 8, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'In_Frame_Del')

#Color coding for RACE classification; try getAnnotations(x = laml) to see available annotations.
race.colors = c("black","yellow","grey")
names(race.colors) = c("African_American", "Asian", "Caucasian")
race.colors = list(RACE = race.colors)

TN.colors=c("black","grey","white")
names(TN.colors) = c("1", "0","unclassified")
TN.colors = list(Triple.negative = TN.colors)

study.colors = c("blue","red")
names(study.colors) = c("stress", "metabolon")
study.colors = list(study = study.colors)



oncoplot(maf = laml, top = 12, fontSize = 12,removeNonMutated = F,
		#colors=col,
		clinicalFeatures = c("RACE","Triple.negative"),annotationColor = c(race.colors,TN.colors),
		sortByAnnotation = T)


## compare EA and AA ##
AA.maf=read.maf(maf="AA.maf")
EA.maf=read.maf(maf="EA.maf")
aa.vs.ea <- mafCompare(m1 = AA.maf, m2 = EA.maf, m1Name = 'AA', m2Name = 'EA', minMut = 4)
print(aa.vs.ea)

lollipopPlot2(m1 = AA.maf, m2 = EA.maf, gene = "TP53", AACol1 = "HGVSp_Short", AACol2 = "HGVSp_Short", m1_name = "AA", m2_name = "EA")
lollipopPlot2(m1 = AA.maf, m2 = EA.maf, gene = "PIK3CA", AACol1 = "HGVSp_Short", AACol2 = "HGVSp_Short", m1_name = "AA", m2_name = "EA")
lollipopPlot2(m1 = AA.maf, m2 = EA.maf, gene = "TTN", AACol1 = "HGVSp_Short", AACol2 = "HGVSp_Short", m1_name = "AA", m2_name = "EA")
lollipopPlot2(m1 = AA.maf, m2 = EA.maf, gene = "MAP3K1", AACol1 = "HGVSp_Short", AACol2 = "HGVSp_Short", m1_name = "AA", m2_name = "EA")

lollipopPlot(maf = laml, gene = 'STARD8', AACol = 'HGVSp_Short', showMutationRate = TRUE)



























rainfallPlot(maf = laml, detectChangePoints = TRUE, fontSize = 12, pointSize = 0.6)






laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)
lollipopPlot(maf = laml, gene = 'TP53', AACol = 'HGVSp_Short', showMutationRate = TRUE)
lollipopPlot(maf = laml, gene = 'PIK3CA', AACol = 'HGVSp_Short', showMutationRate = T,labelPos = 511)
laml.mutload = tcgaCompare(maf = laml, cohortName = 'BreastCancer')



plotVaf(maf = laml, vafCol = 'tumor_freq',top=50)

beeswarm((subset(mutect2.cl,Hugo_Symbol=="ZIC1")$tumor_freq))
#PCLO,PLEC,muc17,4,16,CSMD3



laml.maf = system.file('extdata', 'tcga_laml.maf.gz', package = 'maftools')
laml.clin = system.file('extdata', 'tcga_laml_annot.tsv', package = 'maftools') # clinical inform
test = read.maf(maf = laml.maf, clinicalData = laml.clin)








