### survival ###
mysurvival=function(survival_data,data){
surv_results=data.frame()
data=as.matrix(data)
library(survival)

for (i in 1:nrow(data)) {
qt=quantile(as.numeric(data[i,]),0.5)
surv_results[i,1]=as.numeric(qt)

group=ifelse(data[i,] > qt,"High","Low")
survival_data$group=group

surv=survfit(Surv(SurvivalTime,Death)~factor(group),data=survival_data)
surv_results[i,2]=surv$n[1]
surv_results[i,3]=surv$n[2]

coxph=coxph(Surv(SurvivalTime,Death)~factor(group),data=survival_data)
summary=summary(coxph)
surv_results[i,4]=summary$coefficients[2]
surv_results[i,5]=as.numeric(summary$logtest["pvalue"])
surv_results[i,6]=as.numeric(summary$sctest["pvalue"])
surv_results[i,7]=as.numeric(summary$waldtest["pvalue"])
surv_results[i,8]=cox.zph(coxph)$table[3]

coxph=coxph(Surv(SurvivalTime,Death)~data[i,],data=survival_data)
summary=summary(coxph)
surv_results[i,9]=as.numeric(summary$logtest["pvalue"])
surv_results[i,10]=as.numeric(summary$sctest["pvalue"])
surv_results[i,11]=as.numeric(summary$waldtest["pvalue"])
surv_results[i,12]=cox.zph(coxph)$table[3]
}

rownames(surv_results)=rownames(data)
colnames(surv_results)=c("qt","High","Low","HR","cox.LRT","cox.logrank","cox.wald","zph","cox.linear.LRT","cox.linear.logrank","cox.linear.wald","linear.zph")

surv_results.anno=merge(x=surv_results,y=metabolon_ID,by.x=0,by.y=0)
write.csv(surv_results.anno,"surv_results.csv")

print(head(surv_results.anno))
}


#### volcano plot ####

with(results.anno, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot"))

# Add colored points: red if padj<0.01, orange of log2FC>5, green if both)
with(subset(results.anno, adj.P.Val<.01 ), points(logFC, -log10(P.Value), pch=20, col="orange"))
with(subset(results.anno, logFC>5), points(logFC, -log10(P.Value), pch=20, col="red"))
with(subset(results.anno, logFC< -5), points(logFC, -log10(P.Value), pch=20, col="green"))

with(subset(results.anno, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
with(subset(res, padj<.05 & abs(log2FoldChange)>1), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=.8))



set.seed(123456)
library(ggplot2)
library(ggrepel)
#genes=surv_results.anno
genes <- read.csv("surv_results_scale.csv", header = T,row.names=2)
genes$log2HR=log2(1/genes$HR)
genes=genes[!grepl("X -", genes$BIOCHEMICAL),]


genes$Significant <- ifelse(genes$cox.logrank < 0.05, "p.value < 0.05", "N.S")

ggplot(genes, aes(x = log2HR, y = -log10(cox.logrank),size=2^abs(log2HR))) +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("grey", "red")) +
  scale_size(range = c(0, 6))+
  theme_bw(base_size = 16) +
  geom_text_repel(
    data = subset(genes, cox.logrank < 0.02),
    aes(label = BIOCHEMICAL,size = 2^abs(log2HR)),
    #size = 4,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )



































