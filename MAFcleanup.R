
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


 







