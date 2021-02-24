d=read.table("M_DCdoseKEGG.GSEA.txt",header=T,row.names=1,sep="\t")
d=read.table("ZR_DCdoseKEGG.GSEA.txt",header=T,row.names=1,sep="\t")
d$FDR.q.val=-log10(d$FDR.q.val+0.000001)
d=d[rev(rownames(d)),]

gsea.barplot(d,1)

gsea.barplot<-function(d,cex){
library(RColorBrewer)
d$FDR.q.val=ifelse(d$ES<0,-1*d$FDR.q.val,d$FDR.q.val)

n1=nrow(d[d$ES>0,])
n2=nrow(d[d$ES<0,])
bp=barplot(d$ES,col=ifelse(d$ES>0,colorRampPalette(brewer.pal(9,"Reds"))(n1+1),colorRampPalette(brewer.pal(9,"Blues"))(n1+n2+1)),main="",xlab="Enrichment Score",xlim=c(-0.8,0.8),cex.lab=1.5,cex.axis=1.2,las=2, horiz=T,space=0,border="grey")

text(0.02,bp[n1+1:n2],rownames(d)[n1+1:n2],cex=cex,xpd=T,adj=0)
text(-0.02,bp[1:n1],rownames(d)[1:n1],cex=cex,xpd=T,adj=1)

par(new=T)
plot(d$FDR.q.val,bp,type="b",axes = F, bty = "n", xlab = "", ylab = "",xlim=c(-15,15),ylim=c(0,nrow(d)),col="grey",pch=16)

axis(3,at=seq(-6,6,2),labels=c("inf","4","2","0","2","4","inf"))
mtext("-log10(FDR)",side=3,line=2)
legend("topright",c("FDR","ES"),pch=c(16,0),cex=1.5,col=c("grey","black"),box.col=NA)
}
