library(DESeq2)
library("vsn")
library("RColorBrewer")
library("gplots")

 #
#### data read in ####
data=read.table("RawCountFile_rsemgenes.txt",header=T,row.names=1,sep= '\t')
design = read.csv('design.csv', header = T)
rownames(design) = paste0('X', 1:44, '_sample',1:44)
design



design.cl = subset(design, project == 'naringin project')
data.cl = data[,rownames(design.cl)]
data.cl = round(data.cl)


#### clean data ####
# discard the genes with all samples less than 10 counts#
data_clean=data.cl[apply(data.cl, 1, function(x) !all(x <10)),]
annotation=data.frame(do.call(rbind,strsplit(rownames(data_clean),"_")))
colnames(annotation)=c("ENSG","gene_name")
annotation$gene_id=rownames(data_clean)

dds=DESeqDataSetFromMatrix(countData=data_clean,colData=design.cl,design=~factor(treatment))


#### pca plot ####
vst <- vst(dds, blind=TRUE)
rld = rlogTransformation(dds, blind = T)
pdf("PCA.pdf")
plotPCA(vst[,1:6], intgroup=c("treatment"))
plotPCA(vst[,7:22], intgroup=c("treatment"))
plotPCA(vst[,c(7,8,9,10,15,16,17,18)], intgroup=c("treatment"))
plotPCA(vst[,c(11,12,13,14,19,20,21,22)], intgroup=c("treatment"))
dev.off()




hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), treatment)

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
library(pheatmap)
pheatmap(as.matrix(distsRL),
clustering_distance_rows=distsRL,
clustering_distance_cols=distsRL,
col=colors)



##### comparisons #####
##### comparisons func #####
DESeq2test=function(dds,c1,c2){
dds_c2vsc1=dds[,dds$treatment %in% c(c2,c1)]
dds_c2vsc1$treatment <- droplevels(dds_c2vsc1$treatment)
dds_c2vsc1$treatment <- relevel(dds_c2vsc1$treatment, c1)
dds_c2vsc1 <- DESeq(dds_c2vsc1)

res_c2vsc1=results(dds_c2vsc1)
#write.csv(data.frame(res_c2vsc1),paste0("res_",c2,"vs",c1,"_anno.fpkm.csv"))
return(dds_c2vsc1)
}

rnk=function(a){
x=data.frame(results(a))
x=x[!is.na(x$stat),]
x=x[order(x$stat,decreasing=T),]
x$gene=toupper(gsub(".*_","",rownames(x)))
x = x[!duplicated(x$gene),]
x=x[,c("gene","stat")]
name=deparse(substitute(a))

write.table(x,paste0(name,".rnk"),row.names=F,sep="\t",quote=F)
#return(x)
}


### comparison vehicle NT and Na NT, Vehicle T and Na T, Vehicle NT and Vehicle T, Na NT and Na T ###
# saline
dds_NaCTvsVehCT = dds[,1:6]
NaCTvsVehCT = DESeq2test(dds_NaCTvsVehCT,"Vehicle CT","Na CT")
write.csv(data.frame(results(NaCTvsVehCT)),"NaCTvsVehCT.csv")
rnk(NaCTvsVehCT)


#AOM/DSS
# vehicle NT and Na NT
dds_vehNTvsNaNT = dds[, colData(dds)$treatment %in% c("Vehicle NT","Na NT")]
vehNTvsNaNT = DESeq2test(dds_vehNTvsNaNT,"Vehicle NT","Na NT")

write.csv(data.frame(results(vehNTvsNaNT)),"vehNTvsNaNT.csv")
rnk(vehNTvsNaNT)

# Vehicle T and Na T
dds_vehTvsNaT = dds[, colData(dds)$treatment %in% c("Vehicle T","Na T")]
vehTvsNaT = DESeq2test(dds_vehTvsNaT,"Vehicle T","Na T")

write.csv(data.frame(results(vehTvsNaT)),"vehTvsNaT.csv")
rnk(vehTvsNaT)

# Vehicle NT and Vehicle T
dds_vehNTvsvehT = dds[, colData(dds)$treatment %in% c("Vehicle NT","Vehicle T")]
vehNTvsvehT = DESeq2test(dds_vehNTvsvehT,"Vehicle NT","Vehicle T")

write.csv(data.frame(results(vehNTvsvehT)),"vehNTvsvehT.csv")
rnk(vehNTvsvehT)

# Na NT and Na T
dds_NaNTvsNaT = dds[, colData(dds)$treatment %in% c("Na NT","Na T")]
NaNTvsNaT = DESeq2test(dds_NaNTvsNaT,"Na NT","Na T")

write.csv(data.frame(results(NaNTvsNaT)),"NaNTvsNaT.csv")
rnk(NaNTvsNaT)


# ma plot
pdf("MA-plot.pdf")
plotMA(NaCTvsVehCT,ylim=c(-2,2),main="NaCTvsVehCT") 
plotMA(vehNTvsNaNT,ylim=c(-2,2),main="vehNTvsNaNT") 
plotMA(vehTvsNaT,ylim=c(-2,2),main="vehTvsNaT") 
plotMA(vehNTvsvehT,ylim=c(-2,2),main="vehNTvsvehT") 
plotMA(NaNTvsNaT,ylim=c(-2,2),main="NaNTvsNaT") 
dev.off()


 



### significant heatmap ###
library(pheatmap)
select.dose=rownames(data.frame(subset(data.frame(results(shPRMT5.vs.Vector)),padj<0.05 & abs(log2FoldChange) > 1 )))
pheatmap(assay(rld)[select.dose,],scale="row", cluster_rows=T, show_rownames=F,cluster_cols=T, annotation_col=design[2],cutree_rows=2)
 
#### boxplot a gene #####
d <- plotCounts(dds, gene="ENSG00000148773_MKI67", intgroup="treatment",returnData=TRUE)
d=subset(d,treatment %in% c("mock","DC20","DC50"))
d$treatment=factor(d$treatment,levels=c("mock","DC20","DC50"))
beeswarm(d$count~d$treatment,pch=16,cex=2,xlab="",ylab="count",col=c("lightblue","orange","red"),cex.lab=1.5,cex.axis=1.5,las=3,main="MKI67")

### plot voc plot ###
voc=data.frame(results(shPRMT5.vs.Vector))
voc=voc[complete.cases(voc),]
## Sort by ordered padj
voc_ordered <- voc[order(voc$padj), ] 
## Create a column to indicate which genes to label
voc_ordered$genelabels <- ""
#voc_ordered$genelabels[1:20] <- gsub(".*_","",rownames(voc_ordered)[1:20])
voc_ordered$genelabels[1:30] <- gsub(".*_","",rownames(voc_ordered)[1:30])


ggplot(voc_ordered) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = padj<0.05)) +
  geom_text_repel(aes(x = log2FoldChange, y = -log10(padj), label = genelabels)) +
  ggtitle("shPRMT5.vs.Vector") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 

#### GSEA #####

gsea_folder = 'KEGG_vehNTvsNaNT.GseaPreranked.1622063853729'
gsea_folder = "KEGG_vehTvsNaT.GseaPreranked.1622063871909"
q = 0.05

# pos file
pos.file = list.files(path = gsea_folder, pattern = "^gsea_report_for_na_pos_.*tsv", all.files = FALSE,full.names = T )
pos =read.table(pos.file,sep='\t',header=T)
pos=pos[,c(1,5,8)]
rownames(pos)=gsub('KEGG_','', pos$NAME)
pos=subset(pos,FDR.q.val < q)
pos=pos[order(pos$FDR.q.val,decreasing = T),]

# neg file
neg.file = list.files(path = gsea_folder, pattern = "^gsea_report_for_na_neg_.*tsv", all.files = FALSE,full.names = T )
neg =read.table(neg.file,sep='\t',header=T)
neg=neg[,c(1,5,8)]
rownames(neg)=gsub('KEGG_','',neg$NAME)
neg=subset(neg,FDR.q.val < q)
neg=neg[order(neg$FDR.q.val,decreasing = T),]

# combine two
pos.neg=rbind(subset(neg,FDR.q.val < 0.05),subset(pos,FDR.q.val < 0.05))
pos.neg$FDR.q.val=-log10(pos.neg$FDR.q.val+0.00001)
pos.neg$FDR.q.val=ifelse(pos.neg$ES<0,-1*pos.neg$FDR.q.val,pos.neg$FDR.q.val)
pos.neg$NAME = rownames(pos.neg)

# plot #
ggplot(pos.neg,aes(x=reorder(NAME,ES), y= FDR.q.val,fill = ES)) +
  scale_y_continuous(breaks=c(-5,-2.5,0,2.5,5), labels=c("5", "2.5", "0","2.5","5")) +
  geom_col() + coord_flip() + ylab('-log10(FDR.q.val)') + xlab('') + 
  scale_fill_gradient2(low = "blue",mid="white", high = "red") +
  theme_classic() + ggtitle('KEGG Pathways')  +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),plot.margin = margin(0.1, 0.1, 0.1, 1, "pt") ) +
  geom_text(size = 3.5, aes(label= NAME, y =0),hjust= ifelse(pos.neg$ES > 0 ,1, 0))







##### clusterprofiler #####
df= results(shPRMT5.vs.Vector)
# we want the log2 fold change 
original_gene_list <- df$stat

# name the vector
names(original_gene_list) <- do.call(rbind,strsplit(rownames(df),"_"))[,1]

# omit any NA values 
gene_list <-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

library(clusterProfiler)
library(enrichplot)
organism = "org.Hs.eg.db"

gse <- gseGO(geneList=gene_list, 
             ont ="BP", 
             keyType = "ENSEMBL", 
             #nPerm = 100, 
             minGSSize = 20, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH")


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
 








