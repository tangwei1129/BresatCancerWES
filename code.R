
setwd("~/Desktop/WES/mutect2")
rm(list = ls())


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
mutect2 = read.csv("mutect2_variants_fixed_vep.maf",header=T,sep='\t', skip=1)
dat = mutect2
dim(dat)

### remove Tumor_Sample_Barcode == "Sample_6076"  ###
dat = subset(dat, !Tumor_Sample_Barcode == "Sample_6076")

#### fix the VEP custom database with vcf2maf problem #####
library(tidyverse)
CSQ = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|RefSeq|SOURCE|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|ExAC_AF|ExAC_AF_AFR|ExAC_AF_AMR|ExAC_AF_Adj|ExAC_AF_EAS|ExAC_AF_FIN|ExAC_AF_NFE|ExAC_AF_OTH|ExAC_AF_SAS|gnomADg|gnomADg_AF_afr|gnomADg_AF_amr|gnomADg_AF_asj|gnomADg_AF_eas|gnomADg_AF_fin|gnomADg_AF_nfe|gnomADg_AF_oth|gnomADg_AF_sas|gnomADg_AF_popmax|gnomADg_variant_type|gnomADg_segdup|gnomADg_lcr|gnomADe|gnomADe_non_cancer_AF_afr|gnomADe_non_cancer_AF_amr|gnomADe_non_cancer_AF_asj|gnomADe_non_cancer_AF_eas|gnomADe_non_cancer_AF_fin|gnomADe_non_cancer_AF_nfe|gnomADe_non_cancer_AF_oth|gnomADe_non_cancer_AF_sas|gnomADe_non_cancer_AF_popmax|gnomADe_AF_popmax|gnomADe_variant_type|gnomADe_segdup|gnomADe_lcr"
CSQ.header= paste0("vep.",strsplit(CSQ,"|",fixed = TRUE)[[1]])
length(CSQ.header)

dat.vep = separate(dat, CSQ, CSQ.header, "\\|", remove = T, convert = F)
#str(dat.vep, list.len=ncol(dat.vep))
#####

###
var.class = c('Frame_Shift_Del','Frame_Shift_Ins','Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'Translation_Start_Site','In_Frame_Del')
#
driver.genes = c(  "AKT1",	"ARID1A",	"ARID1B",	"BAP1"	,"BARD1"	,"BRCA1"	,"BRCA2",	"CASP8"	,"CCND1",	"CDH1"	,"CDKN1B"	,"CTCF",	"EP300"	,"ERBB2",	"ESR1",	"ETV6",	"FAT1",
                  "FBLN2"	,"FBXW7"	,"FLNA",	"FOXA1"	,"GATA3",	"IRS4",	"KEAP1"	,"KMT2C"	,"KMT2D"	,"LRP1B",	"MAP2K4",	"MAP3K1",	"MAP3K13"	,"MYH9",	"NCOR1",	"NF1",	"NOTCH1",	
                  "NTRK3",	"PBRM1",	"PCLO",	"PIK3CA",	"PPM1D"	,"PRKDC",	"PTEN",	"RB1"	,"RELN"	,"RUNX1",	"SALL4"	,"SMARCD1",	"SPEN"	,"TBX3"	,"TP53",	"ZMYM3")

################################################################################################# 

write.csv(subset(dat.vep,Hugo_Symbol %in% driver.genes & Variant_Classification %in% var.class),"drivergenes.csv")
topgene(dat.vep)






#############################################
### start to clean ##########################
# remove common_variant 0.001 from gnomADg and gnomADe and ExAC
af = 0.001
dat.vep.af = af.filter(dat.vep, af)
topgene(dat.vep.af)

# annotated duplication variants in cohort with dup times ##
# derived tumor_freq for mutations ##
#annotated duplication variants in cohort with tumor_freq ##
dat.vep.af = dat.vep.af %>% add_count(Chromosome,Start_Position,End_Position, Reference_Allele,Tumor_Seq_Allele2, name = "dup.n")
dat.vep.af$tumor_freq = dat.vep.af$t_alt_count/dat.vep.af$t_depth
dat.vep.af = dat.vep.af %>% 
            group_by(Chromosome,Start_Position,End_Position, Reference_Allele,Tumor_Seq_Allele2) %>% 
            mutate( dup.n.freq = mean(tumor_freq))
dat.vep.af = as.data.frame(dat.vep.af)

#### filter depth ###########
# >= 20 reads total mapped to the location of the mutant and 10 for normal
# >= 3 reads supporting the mutant allele, and 4 less for normal
depth=20
t_reads = 3
n_reads = 4
freq=0.01 ## need to create tumor_freq

dat.vep.af.depth=subset(dat.vep.af,
	(t_depth >= depth) & (n_depth >= 0.5*depth) & (t_alt_count >= t_reads) & (n_alt_count <= n_reads) & tumor_freq > freq)
topgene(dat.vep.af.depth)
##############################################



### add mc3 and genie function ###
match.column = c("Chromosome", "Start_Position")
#match.column2 = c("Chromosome", "Start_Position", "End_Position","Reference_Allele","Tumor_Seq_Allele2")
dat.vep.af.depth.mc3.genie = add.mc3.genie(dat.vep.af.depth, match.column)
table(is.na(dat.vep.af.depth.mc3.genie$mc3.n),is.na(dat.vep.af.depth.mc3.genie$center.number))

### keep mc3 and genie mutations ###
keep.mc3.genie = subset(dat.vep.af.depth.mc3.genie, (!is.na(mc3.n) | (!is.na(center.number))))


dim(keep.mc3.genie)
topgene(keep.mc3.genie)



### filter further for the uncertain mutations ####
filter.mc3.genie = subset(dat.vep.af.depth.mc3.genie, !((!is.na(mc3.n) | (!is.na(center.number)))))

#muts must freq 0.09 and t > 10 or, t > 5 and freq > 0.15 ###
keep.mc3.genie.freq0.1.t10 = subset(filter.mc3.genie,   (t_alt_count >= 10 & tumor_freq >= 0.09) | (t_alt_count > 5 & tumor_freq >= 0.15))
filter.mc3.genie.freq0.1.t10 = subset(filter.mc3.genie,   !((t_alt_count >= 10 & tumor_freq >= 0.09) | (t_alt_count > 5 & tumor_freq >= 0.15)))

#### COSV mutations, too noise ###
#keep.mc3.genie.freq0.1.t10.cosv = filter.mc3.genie.freq0.1.t10[grep('COSV', filter.mc3.genie.freq0.1.t10$Existing_variation),]
#filter.mc3.genie.freq0.1.t10.cosv = filter.mc3.genie.freq0.1.t10[grep('COSV', filter.mc3.genie.freq0.1.t10$Existing_variation, invert = T),]

### keep.mc3.genie   dup.n > 2 dup.n.freq > 0.1 ####
keep.mc3.genie.dup = subset(keep.mc3.genie, dup.n == 1 | (dup.n >= 2 & dup.n.freq > 0.1))
### keep.mc3.genie.freq.0.1.t10  dup.n >= 3 remove, and n >2 remove dup.n.freq > 0.1
keep.mc3.genie.freq0.1.t10.dup = subset(keep.mc3.genie.freq0.1.t10, dup.n == 1 |  (dup.n == 2 & dup.n.freq > 0.1) )
## combine keep sets ###
keep.sets = rbind(keep.mc3.genie.dup, keep.mc3.genie.freq0.1.t10.dup)
keep.sets.cl = subset(keep.sets, !(t_alt_count < 10 & tumor_freq < 0.05 & mc3.n == 2 & is.na(center.number)))
### remove FLAGS genes FLAGS - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4267152/ ###
## remove flags 20 genes ##
#keep.sets.cl = subset(keep.sets.cl, !Hugo_Symbol %in% flags(20))
### remove mt genes ###
keep.sets.cl = subset(keep.sets.cl, !Hugo_Symbol %in% grep("^MT-.*",keep.sets.cl$Hugo_Symbol,value = T))
write.csv(keep.sets.cl, "keep.sets.cl.csv")


#### extract CGC set from multi-callers ####
CGC.mutation = subset(filter.mc3.genie.freq0.1.t10, Hugo_Symbol %in% CGC)
dim(CGC.mutation)
write.csv(CGC.mutation,"CGC.mutation.csv")
#################################################
################  END ###########################
#################################################




topgene(keep.sets.cl)
a= keep.sets.cl
a=read.maf(a)
oncoplot(a,top=50)
titv(maf =a, plot = T, useSyn = TRUE)
 


test2=subset(merge.1caller,t_alt_count > 9 & tumor_freq > 0.15)
test=rbind(genie.set,test1,test2)
final.set=test[!duplicated(test),]
write.csv(final.set,"final.set.csv")
topgene(final.set)
dim(final.set)




#######################################################################################################









#test=subset(merge.only3caller)
final.set.maf=read.maf(final.set)
oncoplot(final.set.maf,top=50)
titv = titv(maf =final.set.maf, plot = T, useSyn = TRUE)
plotTiTv(res =  titv)
 

#plot titv summary
test.raw=merge.3caller.var
test.raw.maf=read.maf(test.raw)
oncoplot(test.raw.maf,top=25)





## get a gene matrix with tumor freq from 'merge.mutect2.strelka.af.depth' ##
gene.matrix.data=merge.mutect2.strelka.af.depth
gene.matrix.data=subset(gene.matrix.data,caller =="mutect2.strelka" & t_alt_count >= 9 )
gene.matrix.data=subset(gene.matrix.data, Variant_Classification %in% c('Frame_Shift_Del','Frame_Shift_Ins','Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Frame_Shift_Ins',
               'In_Frame_Ins', 'Splice_Site', 'Translation_Start_Site','In_Frame_Del'))
gene.matrix.merge=aggregate(gene.matrix.data$tumor_freq,by=list(gene.matrix.data$Hugo_Symbol),FUN=median)
colnames(gene.matrix.merge)=c("genename","median")
gene.matrix.merge$mean=aggregate(gene.matrix.data$tumor_freq,by=list(gene.matrix.data$Hugo_Symbol),FUN=mean)$x
gene.matrix.merge$mad=aggregate(gene.matrix.data$tumor_freq,by=list(gene.matrix.data$Hugo_Symbol),FUN=mad)$x
gene.matrix.merge$sd=aggregate(gene.matrix.data$tumor_freq,by=list(gene.matrix.data$Hugo_Symbol),FUN=sd)$x
gene.matrix.merge$n=aggregate(gene.matrix.data$tumor_freq,by=list(gene.matrix.data$Hugo_Symbol),FUN=length)$x
summary(gene.matrix.merge)
write.csv(gene.matrix.merge,"gene.matrix.merge.t10f0.1.csv")

## LABELing gene wise with mean < 0.1, with tumor_freq > 0.1 ##
tumor_freq.mean = 0.2
tumor_freq.median = 0.15
tumor_n = 6
freq2=0.05
freq3=0.1

gene.matrix.merge$cluster= ifelse(gene.matrix.merge$mean >= tumor_freq.mean & gene.matrix.merge$median >= tumor_freq.median & gene.matrix.merge$n >= tumor_n, freq2,freq3 )
## CDH1,ERBB2,GATA3,MAP3K1,PIK3CA,PPFIA3,TP53 survive ##
merge.mutect2.strelka.af.depth$freq.cutoff=gene.matrix.merge$cluster[match(merge.mutect2.strelka.af.depth$Hugo_Symbol,gene.matrix.merge$genename)]
merge.mutect2.strelka.af.depth$freq.filter=ifelse(merge.mutect2.strelka.af.depth$tumor_freq > merge.mutect2.strelka.af.depth$freq.cutoff,"PASS","Low.Tumor.Freq")

## cut tumor_freq 0.1 and t_alt_count >=5 ##
merge.mutect2.strelka.af.depth.0.1t5=subset(merge.mutect2.strelka.af.depth,tumor_freq > 0.1 & t_alt_count >= 5)




plotmafSummary(test.maf)
oncoplot(test)

#plot titv summary
test=read.maf(a)
titv = titv(maf =test, plot = FALSE, useSyn = TRUE)
plotTiTv(res =  titv)







test=read.maf(merge.mutect2.strelka.af.depth)
test0.1t5=subsetMaf(test,query='freq.filter =="PASS" & t_alt_count >= 5')

 
rescue=read.csv("rescue.maf.txt",header=T,sep="\t") 
rescue=read.maf(rescue)

merge.rescue.test0.1t5=merge_mafs(list(test0.1t5,rescue))
make.mutsig.maf(merge.rescue.test0.1t5,"merge.rescue.test0.1t5")
oncoplot(merge.rescue.test0.1t5,top=22)

#tumor.freq0.1=subsetMaf(test,query='tumor_freq > 0.1')
#tumor.freq0.1t5=subsetMaf(test,query='tumor_freq > 0.1 & t_alt_count >=5')
#merge.rescue=merge_mafs(list(tumor.freq0.1t5,rescue))
#tumor.freq0.15=subsetMaf(test,query='tumor_freq > 0.15')
#tumor.freq0.2=subsetMaf(test,query='tumor_freq > 0.2')


oncoplot(test1,top=24)# test0.1,test0.15,test0.2
oncoplot(tumor.freq0.1)
oncoplot(tumor.freq0.15)
oncoplot(tumor.freq0.2)

make.mutsig.maf=function(x,name){
x1=prepareMutSig(x)
x1=data.frame(x1)
x1$Variant_Classification=gsub("Splice_Region","Intron",x1$Variant_Classification)
write.table(x1,paste0(name,".maf"),sep="\t",quote=F,row.names=F)
}


make.mutsig.maf(tumor.freq0.1,"tumor.freq0.1")
make.mutsig.maf(tumor.freq0.15,"tumor.freq0.15")
make.mutsig.maf(tumor.freq0.2,"tumor.freq0.2")
make.mutsig.maf(test0.1,"test0.1")
make.mutsig.maf(test0.15,"test0.15")
make.mutsig.maf(test0.2,"test0.2")
##### ml MutSig  ####
##### MutSigCV xxxxxx.maf $MUTSIG_REF/exome_full192.coverage.txt $MUTSIG_REF/gene.covariates.txt xxxx.output $MUTSIG_REF/mutation_type_dictionary_file.txt $MUTSIG_REF/chr_files_hg19








write.csv(merge.mutect2.strelka.af.depth,"merge.mutect2.strelka.af.depth.csv")



## clean gene with mad < 0.1,with tumor freq > 0.1
cutoff=0.1
freq2=0.1 #0.05

gene.matrix.mutect2=aggregate(mutect2.cl$tumor_freq,by=list(mutect2.cl$Hugo_Symbol),FUN=median)
colnames(gene.matrix.mutect2)=c("genename","median")
gene.matrix.mutect2$mean=aggregate(mutect2.cl$tumor_freq,by=list(mutect2.cl$Hugo_Symbol),FUN=mean)$x
gene.matrix.mutect2$mad=aggregate(mutect2.cl$tumor_freq,by=list(mutect2.cl$Hugo_Symbol),FUN=mad)$x
gene.matrix.mutect2$sd=aggregate(mutect2.cl$tumor_freq,by=list(mutect2.cl$Hugo_Symbol),FUN=sd)$x
write.csv(gene.matrix.mutect2,"gene.matrix.mutect2.csv")

# choose mad to cut#
gene.matrix.mutect2$cluster= ifelse(gene.matrix.mutect2$mad >= cutoff, freq, freq2)

mutect2.cl.af.depth$freq.cutoff=gene.matrix.mutect2[match(mutect2.cl.af.depth$Hugo_Symbol,gene.matrix.mutect2$genename),]$cluster
mutect2.cl.af.depth$freq.filter=ifelse(mutect2.cl.af.depth$tumor_freq > mutect2.cl.af.depth$freq.cutoff,"PASS","Low.Tumor.Freq")
mutect2.cl.af.depth.freq=subset(mutect2.cl.af.depth,freq.filter =="PASS")

write.table(mutect2.cl.af.depth.freq,"mutect2_merged_cl.maf",sep="\t",quote=F,row.names=F)
AA.maf=subset(mutect2.cl.af.depth.freq,Tumor_Sample_Barcode %in% info[info$RACE=="African American",]$Tumor_Sample_Barcode)
EA.maf=subset(mutect2.cl.af.depth.freq,Tumor_Sample_Barcode %in% info[info$RACE=="Caucasian",]$Tumor_Sample_Barcode)
write.table(AA.maf,"AA.maf",sep="\t",quote=F,row.names=F)
write.table(EA.maf,"EA.maf",sep="\t",quote=F,row.names=F)




## finish clean up ##
#######################################################
info=read.csv("sampleinfo.222.csv",header=T)
info=subset(info,NTDESC =="tumor")
info=subset(info,!Tumor_Sample_Barcode %in% c("6075","6076"))
write.table(info, file='info.tsv', quote=FALSE, sep='\t', row.names = F)


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








