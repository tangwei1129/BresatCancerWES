library(DESeq2)
library("vsn")
library("RColorBrewer")
library("gplots")


library(tidyverse)

###  merge QC files ####
files = list.files(path='./',recursive = T, pattern="RnaSeqMetrics.txt")
rna.metrix = do.call(rbind, lapply(files, function(x) read.table(x,header=T,sep="\t",nrow=1)))
rownames(rna.metrix) = files
write.csv(rna.metrix,'rna.matrix.csv')

files = list.files(path='./',recursive = T, pattern="sample.bam.metric")
bam.metrix = do.call(rbind, lapply(files, function(x) read.table(x,header=T,sep="\t",nrow=1)))
rownames(bam.metrix) = files
write.csv(bam.metrix,'bam.matrix.csv')



# Get the files names

setwd("/Volumes/data/Collaborationstudy/wangyongtao/RNAseq/")
files = list.files(path='./',recursive = T, pattern="sample.genes.results")
length(files)
# First apply read.csv, then rbind
raw.count = do.call(cbind, lapply(files, function(x) read.csv(x,row.names=1,header=T,sep="\t")[4]))
colnames(raw.count)=files
colnames(raw.count)=gsub("fastq/","",colnames(raw.count))
colnames(raw.count)=gsub("/sample.genes.results","",colnames(raw.count))

design=read.csv("design.csv",row.names=1,header=T)
design$name=rownames(design)

data=raw.count[, rownames(design)]
data=round(data)
### annotation ###
annotation=read.table("rat_annotation.txt",header=F,sep='\t')
colnames(annotation)= c('ENSRNOG','genename')
annotation$gene_id=annotation$ENSRNOG


######################################################### 
#######  cell1 ##########################################
#### clean data ####

# discard the genes with all samples less than 10 counts#
data_clean=data[apply(data, 1, function(x) !all(x <10)),]
dim(data_clean)

dds=DESeqDataSetFromMatrix(countData=data_clean,colData=design,design=~factor(treatment))

#### pca plot ####
rld <- rlogTransformation(dds, blind=TRUE)

pdf("PCA.cell1.pdf")
plotPCA(rld, intgroup=c("name"))
plotPCA(rld, intgroup=c("treatment"))
dev.off()

##### comparisons func #####
DESeq2test=function(dds,c1,c2){
dds_c2vsc1=dds[,dds$treatment %in% c(c2,c1)]
dds_c2vsc1$treatment <- droplevels(dds_c2vsc1$treatment)
dds_c2vsc1$treatment <- relevel(dds_c2vsc1$treatment, c1)
dds_c2vsc1 <- DESeq(dds_c2vsc1)

res_c2vsc1=results(dds_c2vsc1)
res_c2vsc1_anno=merge(x=data.frame(res_c2vsc1),y=annotation,by.x="row.names",by.y="gene_id",all.x=T)
write.csv(res_c2vsc1_anno,paste0("res_",c2,"vs",c1,"_anno.csv"))
return(dds_c2vsc1)
}



rnk=function(a){
        x=data.frame(results(a))
        x=x[!is.na(x$stat),]
        x=x[order(x$stat,decreasing=T),]
        x$gene=gsub(".*_","",rownames(x))
        x=x[,c("gene","stat")]
        name=deparse(substitute(a))
        
        write.table(x,paste0(name,".rnk"),row.names=F,sep="\t",quote=F)
        #return(x)
}

######   comparisons  ########
 
res_ZGN.1136.scvsVehicle.Control.sc=DESeq2test(dds,"Vehicle.Control.sc","ZGN.1136.sc")
res_ZGN.1136.scvsblank=DESeq2test(dds,"blank","ZGN.1136.sc")
res_Vehicle.Control.scvsblank=DESeq2test(dds,"blank","Vehicle.Control.sc")

res_ZGN1345.gavagevsWater.control.gavage=DESeq2test(dds,"Water.control.gavage","ZGN1345.gavage")
res_ZGN1345.gavagevsblank=DESeq2test(dds,"blank","ZGN1345.gavage")
res_Water.control.gavagecvsblank=DESeq2test(dds,"blank","Water.control.gavage")

rnk(res_ZGN.1136.scvsVehicle.Control.sc)
rnk(res_ZGN1345.gavagevsWater.control.gavage)




pdf("MA-plot.pdf")
plotMA(res_ZGN.1136.scvsVehicle.Control.sc,ylim=c(-2,2),main="ZGN.1136.scvsVehicle.Control.sc")
plotMA(res_ZGN.1136.scvsblank,ylim=c(-2,2),main="ZGN.1136.scvsblank")
plotMA(res_Vehicle.Control.scvsblank,ylim=c(-2,2),main="Vehicle.Control.scvsblank")

plotMA(res_ZGN1345.gavagevsWater.control.gavage,ylim=c(-2,2),main="ZGN1345.gavagevsWater.control.gavage")
plotMA(res_ZGN1345.gavagevsblank,ylim=c(-2,2),main="ZGN1345.gavagevsblank")
plotMA(res_Water.control.gavagecvsblank,ylim=c(-2,2),main="Water.control.gavagecvsblank")

dev.off()
 
####  plot genes ######
#### boxplot a gene #####
genelist=read.table('genelist.txt',header=F,sep='\t',row.names = 1)
colnames(genelist)='gene_id'
library(pheatmap)
library(beeswarm)
pheatmap(assay(rld)[rownames(rld) %in% rownames(genelist),],scale='row',annotation_col = design['treatment'])

genelist$id = rownames(genelist)
genelist.cl=genelist[rownames(genelist) %in% rownames(rld),]

pdf('geneplot.pdf')


for (i in seq(nrow(genelist.cl))) {
              
d <- plotCounts(dds, gene=rownames(genelist.cl)[i], intgroup="treatment",returnData=TRUE)
d$treatment=factor(d$treatment,levels=c('blank', 'Vehicle.Control.sc' , 'ZGN.1136.sc' ,'Water.control.gavage', 'ZGN1345.gavage'))

par(mar=c(13.1,4.1,4.1,2.1))
beeswarm(d$count~d$treatment,pch=16,cex=2,xlab="",ylab="count",col=c("lightblue","orange","red","orange","red"),cex.lab=1.5,cex.axis=1.5,las=3,main=genelist.cl$gene_id[i])
}
dev.off()















###
res=res_diffvsundiff
res=res[!is.na(res$padj),]


alpha <- 0.05 # Threshold on the adjusted p-value
cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
plot(res$log2FoldChange, -log10(res$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=1 )
with(subset(res, padj<0.05 & log2FoldChange > 0), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.9))
with(subset(res, padj<0.05 & log2FoldChange < 0), points(log2FoldChange, -log10(padj), pch=20, col="green", cex=0.9))

abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

gn.selected <- abs(res$log2FoldChange) > 2.5 & res$padj < alpha 
text(res$log2FoldChange[gn.selected],
     -log10(res$padj)[gn.selected],
     lab=rownames(res)[gn.selected ], cex=0.4)
 







### emt genes correlation in cell lines ###
gene11=read.table("genelist.txt",header=F)$V1
gene11=c("LINC00675",as.character(gene11))
gene11.id=annotation[annotation$gene_name %in% gene11,] 
rownames(rld)=gsub("\\..*","",rownames(rld))
gene11.exp=assay(rld)[rownames(rld) %in% gene11.id$gene_id,]
rownames(gene11.exp)=as.character(gene11.id$gene_name[match(rownames(gene11.exp),gene11.id$gene_id)])
library(pheatmap)
pheatmap(cor(t(gene11.exp[,-3]),method="pearson"),scale="none",treeheight_col=8,treeheight_row=8)
 


### emt genes correlation in TCGA ###
library(tidyverse)
TCGA=read_csv("COAD-TCGA-All data from Fathi-10-21-16.csv")
TCGA.id=grep(paste(paste0("\\|",gene11,"\\|"),collapse="|"),TCGA$X1, value=TRUE)
TCGA=data.frame(TCGA)
rownames(TCGA)=TCGA$X1
TCGA=TCGA[,-1]
TCGA.exp=TCGA[TCGA.id,]

TCGA.exp.N=TCGA.exp[,grep("_11$",colnames(TCGA.exp))]
TCGA.exp.T=TCGA.exp[,grep("_1$",colnames(TCGA.exp))]
pheatmap(cor(t(TCGA.exp.T),method="pearson"),scale="none",treeheight_col=8,treeheight_row=8,show_colnames=F)
pheatmap(cor(t(TCGA.exp.N),method="pearson"),scale="none",treeheight_col=8,treeheight_row=8,show_colnames=F)




















 

#### boxplot a gene #####
d <- plotCounts(dds, gene=7143, intgroup="treatment",returnData=TRUE)
d=subset(d,treatment %in% c("mock","DC20","DC50"))
d$treatment=factor(d$treatment,levels=c("mock","DC20","DC50"))
beeswarm(d$count~d$treatment,pch=16,cex=2,xlab="",ylab="count",col=c("lightblue","orange","red"),cex.lab=1.5,cex.axis=1.5,las=3,main="MKI67")

d <- plotCounts(dds, gene=14923, intgroup="treatment",returnData=TRUE)
d=subset(d,treatment %in% c("mock","DC20","DC50"))
d$treatment=factor(d$treatment,levels=c("mock","DC20","DC50"))
beeswarm(d$count~d$treatment,pch=16,cex=2,xlab="",ylab="count",col=c("lightblue","orange","red"),cex.lab=1.5,cex.axis=1.5,las=3,main="MALAT1")


### 11 gene index ###
gene11=c("BIRC5", "CCNB1", "CDC20", "CEP55", "MKI67", "NDC80", "NUF2", "PTTG1", "RRM2", "TYMS", "UBE2C")
gene11.id=annotation[annotation$gene_name %in% gene11,]$gene_id
gene11.index=cbind(design,matrix(colSums(assay(rld)[rownames(rld) %in% gene11.id,])/11))
colnames(gene11.index)[3]="Index"
d=subset(gene11.index,treatment %in% c("mock","DC20","DC50"))
d$treatment=factor(d$treatment,levels=c("mock","DC20","DC50"))
beeswarm(d$Index~d$treatment,pch=16,cex=2,xlab="",ylab="11-gene Index",col=c("lightblue","orange","red"),cex.lab=1.5,cex.axis=1.5,las=3,main="Proliferation")
text(3,10.2,"p=0.003",cex=1.2)




### significant heatmap ###
library(pheatmap)
select.dose=rownames(data.frame(subset(data.frame(results(ZR_DCdose)),padj<0.05 )))
pheatmap(assay(rld)[select.dose,c(1:12)],scale="row", cluster_rows=T, show_rownames=F,cluster_cols=T, annotation_col=design[,c(2,3)],cutree_rows=2)

select.dose=rownames(data.frame(subset(data.frame(results(M_DCdose)),padj<0.05 )))
mat=assay(rld)[select.dose,c(13:16,21:24,17:20)]
mat_cluster_cols=hclust(dist(t(mat)))
mat_cluster_cols$order=c(1,2,3,4,6,7,5,8,12,11,9,10)
pheatmap(assay(rld)[select.dose,c(13:16,21:24,17:20)],scale="row", cluster_rows=T, show_rownames=F,cluster_cols=mat_cluster_cols, annotation_col=design[3],cutree_rows=2)


res_DC20vsMock=data.frame(results(M.DC20.vs.mock))
res_DC50vsMock=data.frame(results(M.DC50.vs.mock))
res_DC50vsDC20=data.frame(results(M.DC50.vs.DC20))
res_DCdose=data.frame(results(M_DCdose))

select.union=Reduce(union, list(rownames(data.frame(subset(res_DC20vsMock,padj<0.05 & abs(log2FoldChange) >0.58 ))),
		rownames(data.frame(subset(res_DC50vsMock,padj<0.05 & abs(log2FoldChange) >0.58))),
		rownames(data.frame(subset(res_DC50vsDC20,padj<0.05 & abs(log2FoldChange) >0.58)))
		))
pheatmap(assay(rld)[select.union,c(13:24)],scale="row", cluster_rows=T, show_rownames=F,cluster_cols=T, annotation_col=design[,c(2,3)],cutree_rows=2)




select=Reduce(intersect, list(rownames(data.frame(subset(res_DC20vsMock,padj<0.05))),
		rownames(data.frame(subset(res_DC50vsMock,padj<0.05 ))),
		rownames(data.frame(subset(res_DC50vsDC20,padj<0.05 )))
		))

select=Reduce(setdiff, list(rownames(data.frame(subset(res_DC20vsMock,padj<0.05))),
		rownames(data.frame(subset(res_DC50vsMock,padj<0.05 ))),
		rownames(data.frame(subset(res_DC50vsDC20,padj<0.05 )))
		))

select.list=list(rownames(data.frame(subset(res_DC20vsMock,padj<0.05))),
		rownames(data.frame(subset(res_DC50vsMock,padj<0.05 ))),
		rownames(data.frame(subset(res_DC50vsDC20,padj<0.05 )))
		)
library(VennDiagram)
x=rownames(data.frame(subset(res_DC20vsMock,padj<0.05 & abs(log2FoldChange) >0.58)))
y=rownames(data.frame(subset(res_DC50vsMock,padj<0.05 & abs(log2FoldChange) >0.58)))
z=rownames(data.frame(subset(res_DC50vsDC20,padj<0.05 & abs(log2FoldChange) >0.58)))

venn.diagram(list(DC20=x,DC50=y,DC50vs20=z),
	main="Number of overlapped genes", filename="Number of overlapped genes",
	category.names = c("DC20","DC50","DC50vsDC20"),
	na="remove",lty=0,fill = c("violet", "red4","orangered"),scaled=T,euler.d=T,cex=2,cat.cex=2)

select.uniq=Reduce(union,list(Reduce(setdiff,list(x,y,z)),
		Reduce(setdiff,list(y,x,z)),
		Reduce(setdiff,list(z,x,y))
		))
pheatmap(assay(rld)[select.uniq,c(1:6,13,14,15)],scale="row", cluster_rows=T, show_rownames=F,cluster_cols=T, annotation_col=design)



#### boxplot a gene #####
d <- plotCounts(dds, gene="ENSG00000148773_MKI67", intgroup="treatment",returnData=TRUE)
d=subset(d,treatment %in% c("mock","DC20","DC50"))
d$treatment=factor(d$treatment,levels=c("mock","DC20","DC50"))
beeswarm(d$count~d$treatment,pch=16,cex=2,xlab="",ylab="count",col=c("lightblue","orange","red"),cex.lab=1.5,cex.axis=1.5,las=3,main="MKI67")

 


### 11 gene index ###
gene11=c("BIRC5", "CCNB1", "CDC20", "CEP55", "MKI67", "NDC80", "NUF2", "PTTG1", "RRM2", "TYMS", "UBE2C")
gene11.id=annotation[annotation$gene_name %in% gene11,]$gene_id
gene11.index=cbind(design,matrix(colSums(assay(rld)[rownames(rld) %in% gene11.id,])/11))
colnames(gene11.index)[5]="Index"

gene11.index$treatment=factor(gene11.index$treatment,levels=c("mock","DC20","DC50"))
d_ZR=subset(gene11.index, cell == "ZR7530")
d_M=subset(gene11.index, cell == "MDMBA175")

beeswarm(d_ZR$Index~d_ZR$treatment,pch=16,cex=2,xlab="",ylab="11-gene Index",col=c("lightblue","orange","red"),cex.lab=1.5,cex.axis=1.5,las=3,main="ZR7530 Proliferation")
text(3,9.80,"p = 0.24",cex=1.2)

beeswarm(d_M$Index~d_M$treatment,pch=16,cex=2,xlab="",ylab="11-gene Index",col=c("lightblue","orange","red"),cex.lab=1.5,cex.axis=1.5,las=3,main="MDMBA175 Proliferation")
text(2.5,9.95,"p = 3.8e-11",cex=1.2)




## pathway view ##
map="00100"
map="00140"

gt.data=50*(res_DCdose["log2FoldChange"])
rownames(gt.data)=substr(rownames(res_DCdose),1,15)
library(pathview)
pv.out <- pathview(
	 gene.data =gt.data,
	 gene.idtype = "ENSEMBL", limit=list(gene=1,cpd=1),
	 pathway.id = map, species = "hsa", out.suffix = map, keys.align = "y", 
       kegg.native = T, match.data = T, key.pos = "topright",same.layer=F)

plot.name <- paste(map, map, "eps", sep = ".") 
com_set=data.frame(cbind(res_DCdose,res_DC20vsMock[,c(2,5,6)],res_DC50vsMock[,c(2,5,6)]))
com_set$gene_name=gsub(".*_","",rownames(com_set))

dds_M <- estimateSizeFactors(dds_M)
counts_norm=data.frame(counts(dds_M,normalize=T))
counts_norm$gene_name=gsub(".*_","",rownames(counts_norm))

gene.00140=pv.out$plot.data.gene 
com_set00140=merge(x=gene.00140,y=com_set,by.x="labels",by.y="gene_name")
com_set00140=com_set00140[!duplicated(com_set00140$labels),]
com_set00140=merge(x=com_set00140,y=counts_norm,by.x="labels",by.y="gene_name",all.x=T)
write.csv(com_set00140,"com_set00140.csv")

gene.00100=pv.out$plot.data.gene 
com_set00100=merge(x=gene.00100,y=com_set,by.x="labels",by.y="gene_name")
com_set00100=com_set00100[!duplicated(com_set00100$labels),]
com_set00100=merge(x=com_set00100,y=counts_norm,by.x="labels",by.y="gene_name",all.x=T)
write.csv(com_set00100,"com_set00100.csv")
 



