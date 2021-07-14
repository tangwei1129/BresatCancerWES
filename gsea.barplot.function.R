

#gsea.barplot(pos,1)

gsea.barplot<-function(d,cex){
library(RColorBrewer)
d$FDR.q.val=-log10(d$FDR.q.val+0.00001)
d$FDR.q.val=ifelse(d$ES<0,-1*d$FDR.q.val,d$FDR.q.val)

n1=nrow(d[d$ES>0,])
n2=nrow(d[d$ES<0,])
bp=barplot(d$FDR.q.val,col=ifelse(d$ES>0,colorRampPalette(brewer.pal(9,"Reds"))(n1+n2+1),colorRampPalette(brewer.pal(9,"Blues"))(n1+1)),main="",cex.lab=1.5,cex.axis=1.2,las=2, horiz=T,space=0,border="grey",axes = F,xlab='')

text(0.02,bp[1:n1+n2],rownames(d)[1:n1+n2],cex=cex,xpd=T,adj=1.05)
text(-0.02,bp[1:n2],rownames(d)[1:n2],cex=cex,xpd=T,adj=-0.05)

#par(new=T)
#plot(d$FDR.q.val,bp,type="b",axes = F, bty = "n", xlab = "", ylab = "",xlim=c(-15,15),ylim=c(0,nrow(d)),col="grey",pch=16)

axis(3,at=seq(-6,6,2),labels=c("","4","2","0","2","4",""))
mtext("-log10(FDR)",side=3,line=2)
#legend("topright",c("FDR","ES"),pch=c(16,0),cex=1.5,col=c("grey","black"),box.col=NA)
}
