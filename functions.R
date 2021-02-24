
### function of simple verison of dataset ###
simple=function(maf1,caller1,maf2,caller2,maf3,caller3){
maf1.s=subset(maf1,Variant_Classification %in% var.class & Hugo_Symbol %in% driver.genes)
maf1.s$caller=caller1
maf2.s=subset(maf2,Variant_Classification %in% var.class & Hugo_Symbol %in% driver.genes)
maf2.s$caller=caller2
maf3.s=subset(maf3,Variant_Classification %in% var.class & Hugo_Symbol %in% driver.genes)
maf3.s$caller=caller3
maf=rbind(maf1.s,maf2.s,maf3.s)
write.csv(maf,"simple.maf.csv")
return(maf)
}
#simple.maf=simple(mutect2.cl,"mutect2",mutect,"mutect",strelka.cl,"strelka")

################################################################################################
## merge strelka and mutect2 and mutect ##
##### function to merge two maf file with caller ####
merge.maf=function(maf1,caller1,maf2,caller2){
maf1$index=paste(maf1$Tumor_Sample_Barcode,maf1$Chromosome,maf1$Start_Position,maf1$End_Position,sep=".")
maf2$index=paste(maf2$Tumor_Sample_Barcode,maf2$Chromosome,maf2$Start_Position,maf2$End_Position,sep=".")
index.maf1.maf2=intersect(maf1$index,maf2$index)

maf1.maf2=subset(maf1,index %in% index.maf1.maf2)
maf1.maf2$caller=paste(caller1,caller2,sep=".")

maf1.diff=subset(maf1,!index %in% index.maf1.maf2)
maf2.diff=subset(maf2,!index %in% index.maf1.maf2)
maf1.diff$caller=caller1
maf2.diff$caller=caller2
maf=rbind(maf1.maf2,maf1.diff,maf2.diff)
return(maf)
}

##### function to merge THREE maf file with caller ####
merge.maf3=function(maf1,caller1,maf2,caller2,maf3,caller3){
maf1$index=paste(maf1$Tumor_Sample_Barcode,maf1$Chromosome,maf1$Start_Position,maf1$End_Position,sep=".")
maf2$index=paste(maf2$Tumor_Sample_Barcode,maf2$Chromosome,maf2$Start_Position,maf2$End_Position,sep=".")
maf3$index=paste(maf3$Tumor_Sample_Barcode,maf3$Chromosome,maf3$Start_Position,maf3$End_Position,sep=".")

index.maf1.maf2.maf3=intersect(intersect(maf1$index,maf2$index),maf3$index)
index.maf1.maf2=setdiff(intersect(maf1$index,maf2$index), index.maf1.maf2.maf3)
index.maf1.maf3=setdiff(intersect(maf1$index,maf3$index), index.maf1.maf2.maf3)
index.maf2.maf3=setdiff(intersect(maf2$index,maf3$index), index.maf1.maf2.maf3)

maf1.maf2.maf3=subset(maf1,index %in% index.maf1.maf2.maf3)
maf1.maf2.maf3$caller=paste(caller1,caller2,caller3,sep=".")

maf1.maf2=subset(maf1,index %in% index.maf1.maf2)
maf1.maf2$caller=paste(caller1,caller2,sep=".")

maf1.maf3=subset(maf1,index %in% index.maf1.maf3)
maf1.maf3$caller=paste(caller1,caller3,sep=".")

maf2.maf3=subset(maf2,index %in% index.maf2.maf3)
maf2.maf3$caller=paste(caller2,caller3,sep=".")

maf1.diff=subset(maf1,!index %in% c(index.maf1.maf2.maf3,index.maf1.maf2,index.maf1.maf3))
maf2.diff=subset(maf2,!index %in% c(index.maf1.maf2.maf3,index.maf1.maf2,index.maf2.maf3))
maf3.diff=subset(maf3,!index %in% c(index.maf1.maf2.maf3,index.maf1.maf3,index.maf2.maf3))

maf1.diff$caller=caller1
maf2.diff$caller=caller2
maf3.diff$caller=caller3

maf=rbind(maf1.maf2.maf3,maf1.maf2,maf1.maf3,maf2.maf3,maf1.diff,maf2.diff,maf3.diff)
return(maf)
}
##merge.3caller=merge.maf3(mutect2.cl,"mutect2",mutect,"mutect",strelka.cl,"strelka")
## write.csv(merge.3caller,"merge.3caller.csv")
#######################################################################################################


#################################################################################################
#### function to table key genes coding mutation ####
topgene=function(x){
  ###
  var.class = c('Frame_Shift_Del','Frame_Shift_Ins','Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Frame_Shift_Ins',
                'In_Frame_Ins', 'Splice_Site', 'Translation_Start_Site','In_Frame_Del')
  #
  driver.genes = c("TP53","PIK3CA","CDH1","GATA3","ERBB2","AKT1","FBXW7","MAP3K1","NCOR1","PTEN","MAP2K4","BRCA1","KMT2C","RUNX1")
  
a=subset(x,Hugo_Symbol %in% driver.genes)
a=subset(a,Variant_Classification %in% var.class)
a$Hugo_Symbol = as.factor(a$Hugo_Symbol)
b=table(droplevels(a$Hugo_Symbol))
print(b)
print(sum(b))
}
###################################################
#topgene(merge.3caller)
#CDH1  ERBB2  GATA3 MAP3K1 PIK3CA   TP53 
# 14     12     15     15     38     49
###########################################################


######################################################################
##### function to filter by af ##############
af.filter=function(x,af){
  #AF.name=c('AF','AFR_AF','AMR_AF','ASN_AF','EAS_AF','EUR_AF','SAS_AF','AA_AF','EA_AF',
	#   'ExAC_AF','ExAC_AF_AFR','ExAC_AF_AMR','ExAC_AF_EAS','ExAC_AF_FIN','ExAC_AF_NFE','ExAC_AF_OTH','ExAC_AF_SAS',
	#    'gnomAD_AF','gnomAD_AFR_AF','gnomAD_AMR_AF','gnomAD_ASJ_AF','gnomAD_EAS_AF','gnomAD_FIN_AF','gnomAD_NFE_AF','gnomAD_OTH_AF','gnomAD_SAS_AF'
	#   )
  
AF.name = grep('AF', colnames(x),value = T)  
af1=data.frame(x)[,colnames(x) %in% AF.name]
for (i in AF.name){
if(!is.numeric(af1[[i]])){af1[[i]]=as.numeric(as.character(af1[[i]]))}
}
af1[is.na(af1)]=0
af2= x[apply(af1, 1,function(x) !(any(x > af)) ),]
return(af2)
}


#########################################################################



############ funciton to merge cosm or genie database ########
merge.db = function(data,db,dbname){
          data$index = paste0(data$Chromosome,":", data$Start_Position)
	    db$index = paste0(db$Chromosome,":", db$Start_Position)
		db=db[,c(dbname,"index")]
	    data.db=merge(x=data,y=db,by.x="index",by.y="index",all.x=T)
          data.db=data.db[,-1]
		return(data.db)
}


## add mc3 and genie function ###

add.mc3.genie = function(x, match.column) {
  # match.column = c("Chromosome", "Start_Position", "End_Position","Reference_Allele","Tumor_Seq_Allele2")
  # match.column = c("Chromosome", "Start_Position") 
  library(tidyverse)  
  ### add MC3 label ###
  mc3=read.csv("../mc3.cl.csv",header = T)
  mc3$Chromosome = paste0("chr", mc3$Chromosome)
  mc3=mc3[, c(match.column,"mc3.n")]
  mc3 = mc3[!duplicated(mc3[,match.column]),]
  x.mc3 = x %>% left_join(mc3, by = match.column)
  
  ### add genie label ###
  genie=read.csv("../genie.cl.csv",header=T)
  genie$Chromosome = paste0("chr",genie$Chromosome)
  genie = genie[, c(match.column,"center.number")]
  genie = genie[!duplicated(genie[,match.column]),]
    x.mc3.genie = x.mc3 %>% left_join(genie, by = match.column)
  return(x.mc3.genie)
}

flags = function(top = NULL){
  top100flags = c("TTN", "MUC16", "OBSCN", "AHNAK2", "SYNE1", "FLG", "MUC5B",
                  "DNAH17", "PLEC", "DST", "SYNE2", "NEB", "HSPG2", "LAMA5", "AHNAK",
                  "HMCN1", "USH2A", "DNAH11", "MACF1", "MUC17", "DNAH5", "GPR98",
                  "FAT1", "PKD1", "MDN1", "RNF213", "RYR1", "DNAH2", "DNAH3", "DNAH8",
                  "DNAH1", "DNAH9", "ABCA13", "APOB", "SRRM2", "CUBN", "SPTBN5",
                  "PKHD1", "LRP2", "FBN3", "CDH23", "DNAH10", "FAT4", "RYR3", "PKHD1L1",
                  "FAT2", "CSMD1", "PCNT", "COL6A3", "FRAS1", "FCGBP", "DNAH7",
                  "RP1L1", "PCLO", "ZFHX3", "COL7A1", "LRP1B", "FAT3", "EPPK1",
                  "VPS13C", "HRNR", "MKI67", "MYO15A", "STAB1", "ZAN", "UBR4",
                  "VPS13B", "LAMA1", "XIRP2", "BSN", "KMT2C", "ALMS1", "CELSR1",
                  "TG", "LAMA3", "DYNC2H1", "KMT2D", "BRCA2", "CMYA5", "SACS",
                  "STAB2", "AKAP13", "UTRN", "VWF", "VPS13D", "ANK3", "FREM2",
                  "PKD1L1", "LAMA2", "ABCA7", "LRP1", "ASPM", "MYOM2", "PDE4DIP",
                  "TACC2", "MUC2", "TEP1", "HELZ2", "HERC2", "ABCA4")
  
  if(is.null(top)){
    top100flags
  }else{
    top100flags[1:top]
  }
}


CGC = read.csv("Census_allSat Dec 12 18_28_16 2020.csv",header =T)
CGC = CGC$Gene.Symbol

topgene.CGC =function(x){
  ###
  var.class = c('Frame_Shift_Del','Frame_Shift_Ins','Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Frame_Shift_Ins',
                'In_Frame_Ins', 'Splice_Site', 'Translation_Start_Site','In_Frame_Del')
  #
  driver.genes =  CGC
  a=subset(x,Hugo_Symbol %in% driver.genes)
  a=subset(a,Variant_Classification %in% var.class)
  a$Hugo_Symbol = as.factor(a$Hugo_Symbol)
  b=table(droplevels(a$Hugo_Symbol))
  print(b)
  print(sum(b))
}
























