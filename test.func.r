### survival modified ###
mysurvival.m=function(survival_data,data){
surv_results=data.frame()
data=as.matrix(data)
library(survival)

for (i in 1:nrow(data)) {
qt=quantile(as.numeric(data[i,]),0.5)
surv_results[i,1]=as.numeric(qt)

group=ifelse(data[i,] > qt,"High","Low")
survival_data$group=group

surv=survfit(Surv(SurvivalTime,Death)~group.GCDC + factor(group),data=survival_data)
surv_results[i,2]=surv$n[1]
surv_results[i,3]=surv$n[2]

coxph=coxph(Surv(SurvivalTime,Death)~group.GCDC+ factor(group),data=survival_data)
summary=summary(coxph)
surv_results[i,4]=summary$coefficients[1,2]
surv_results[i,5]=summary$coefficients[1,5]

surv_results[i,6]=as.numeric(summary$sctest["pvalue"])
surv_results[i,7]=as.numeric(summary$waldtest["pvalue"])
surv_results[i,8]=cox.zph(coxph)$table[3]

coxph=coxph(Surv(SurvivalTime,Death)~group.GCDC+ data[i,],data=survival_data)
summary=summary(coxph)
surv_results[i,9]=summary$coefficients[1,2]
surv_results[i,10]=summary$coefficients[1,5]

surv_results[i,11]=as.numeric(summary$waldtest["pvalue"])
surv_results[i,12]=cox.zph(coxph)$table[3]
}

rownames(surv_results)=rownames(data)
colnames(surv_results)=c("qt","High","Low","HR","cox.LRT","cox.logrank","cox.wald","zph","cox.linear.LRT","cox.linear.logrank","cox.linear.wald","linear.zph")

surv_results.anno=merge(x=surv_results,y=metabolon_ID,by.x=0,by.y=0)
write.csv(surv_results.anno,"surv_results2.csv")

print(head(surv_results.anno))
}

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


