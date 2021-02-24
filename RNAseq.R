library(DESeq2)
library("vsn")
library("RColorBrewer")
library("gplots")

 #
#### data read in ####
data=read.table("RawCountFile_rsemgenes.txt",header=T,sep="\t",row.names=1)
design=read.table("design.txt",row.names=1,header=T,sep="\t")
design$samplename=rownames(design)


data=data[, rownames(design)]
data=round(data)


#### clean data ####
# discard the genes with all samples less than 10 counts#
data_clean=data[apply(data, 1, function(x) !all(x <10)),]
annotation=data.frame(do.call(rbind,strsplit(rownames(data_clean),"_")))
colnames(annotation)=c("ENSG","gene","ENSG","gene_name")
annotation$gene_id=rownames(data_clean)

dds=DESeqDataSetFromMatrix(countData=data_clean,colData=design,design=~factor(treatment))


#### pca plot ####
rld <- rlogTransformation(dds, blind=TRUE)

pdf("PCA.pdf")
plotPCA(rld, intgroup=c("cell"))
plotPCA(rld[,c(1:12)], intgroup=c("cell","treatment"))
plotPCA(rld[,c(13:24)],intgroup=c("cell","treatment"))
dev.off()




hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), treatment)

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
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
x$gene=gsub(".*_","",rownames(x))
x=x[,c("gene","stat")]
name=deparse(substitute(a))

write.table(x,paste0(name,".rnk"),row.names=F,sep="\t",quote=F)
#return(x)
}


### comparison ###
dds_ZR=dds[,1:12]
dds_M=dds[,13:24]

## ZR ##
ZR.DC50.vs.mock=DESeq2test(dds_ZR,"mock","DC50")
write.csv(data.frame(results(ZR.DC50.vs.mock)),"ZR.DC50.vs.mock.csv")
ZR.DC20.vs.mock=DESeq2test(dds_ZR,"mock","DC20")
write.csv(data.frame(results(ZR.DC20.vs.mock)),"ZR.DC20.vs.mock.csv")
ZR.DC50.vs.DC20=DESeq2test(dds_ZR,"DC20","DC50")
write.csv(data.frame(results(ZR.DC50.vs.DC20)),"ZR.DC50.vs.DC20.csv")

rnk(ZR.DC50.vs.mock)
rnk(ZR.DC20.vs.mock)
rnk(ZR.DC50.vs.DC20)

ZR_DCdose=dds_ZR
design(ZR_DCdose)=~dose
ZR_DCdose <- DESeq(ZR_DCdose, betaPrior=F)
write.csv(data.frame(results(ZR_DCdose)),"ZR_DCdose.csv")
rnk(ZR_DCdose)

## M ##
M.DC50.vs.mock=DESeq2test(dds_M,"mock","DC50")
write.csv(data.frame(results(M.DC50.vs.mock)),"M.DC50.vs.mock.csv")
M.DC20.vs.mock=DESeq2test(dds_M,"mock","DC20")
write.csv(data.frame(results(M.DC20.vs.mock)),"M.DC20.vs.mock.csv")
M.DC50.vs.DC20=DESeq2test(dds_M,"DC20","DC50")
write.csv(data.frame(results(M.DC50.vs.DC20)),"M.DC50.vs.DC20.csv")

rnk(M.DC50.vs.mock)
rnk(M.DC20.vs.mock)
rnk(M.DC50.vs.DC20)

M_DCdose=dds_M
design(M_DCdose)=~dose
M_DCdose <- DESeq(M_DCdose, betaPrior=F)
write.csv(data.frame(results(M_DCdose)),"M_DCdose.csv")
rnk(M_DCdose)
 
 
 





pdf("MA-plot.pdf")
plotMA(ZR.DC50.vs.mock,ylim=c(-2,2),main="ZR7530 DC 50um vs mock")
plotMA(ZR.DC20.vs.mock,ylim=c(-2,2),main="ZR7530 DC 20um vs mock")
plotMA(M.DC50.vs.mock,ylim=c(-2,2),main="MDMBA175 DC 50um vs mock")
plotMA(M.DC20.vs.mock,ylim=c(-2,2),main="MDMBA175 DC 20um vs mock")
dev.off()



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
 








