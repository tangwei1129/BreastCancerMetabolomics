## data input ##
raw=read.csv("Metabolon data12.14.2015.csv",header=T,row.names=1)
scale=read.csv("Metabolon data12.14.2015_scale.csv",header=T,row.names=1)

metabolon_ID=read.csv("IDanno.csv",header=T,row.names=1)
phenotype=read.csv("phenotype.csv",header=T,row.names=1)
phenotype=phenotype[colnames(raw),]

## stats  of missing ##

## all samples 132 ##

numNAs <- apply(raw, 1, function(z) sum(is.na(z)))

# tumor 67 #
phenotype.T=phenotype[phenotype$Status=="Tumor",]
raw.T=raw[,colnames(raw) %in% rownames(phenotype.T)]
numNAs.T=apply(raw.T, 1, function(z) sum(is.na(z)))

missing=cbind(numNAs,numNAs.T)
missing=data.frame(missing)
missing$numNAs.p= missing$numNAs / ncol(raw)
missing$numNAs.T.p= missing$numNAs.T / ncol(raw.T)
write.csv(missing,"missing.csv")

###choose 75% missing in tumor as cut-off ###
cut=0.75
missing.cl=missing[missing$numNAs.T.p <=0.75,]

raw.cl=raw[rownames(raw) %in% rownames(missing.cl),]
scale.cl=scale[rownames(scale) %in% rownames(missing.cl),]


#### survival ####
survival_data=phenotype[phenotype$Status=="Tumor",]
survival_data$Death=survival_data$Death
survival_data$SurvivalTime=survival_data$Survival..Months


data=log2(scale.cl)

data=data[,colnames(data) %in% rownames(survival_data)]

table(colnames(data)==rownames(survival_data))




library("iDOS");
library("survival");
library("glmnet");
library("survcomp")

### survival specific datasets ####
cancer.types <- c("BRCA")
# create a map for stage variable
BRCA.stage.map <-  c("I" = "1", "II" = "2", "IIA" = "2","IIB"="2", "IIIa"="3","IIIA"="3","IIIB"="3")

# standardise BRCA TCGA age and stage
# age
#ann=survival_data[,c("AGE","Stage","Death","SurvivalTime")]
ann=survival_data[,c("AGE","Stage","Death","SurvivalTime","Tumor.grade","Tumor.size_recode","Menopause","NeoAdjuvant.Therapy","Hormone.therapy","Chemo.therapy","Node","BMI","RACE")]

ann$AGE[which(ann$AGE <= 50)] <- 0;
ann$AGE[which(ann$AGE > 50)] <- 1;
ann$age=ann$AGE
# stage
ann$stage <- ann$Stage;
ann$stage <- BRCA.stage.map[ann$stage];
ann$stage <- as.numeric(ann$stage);
ann$RACE=as.numeric(ann$RACE)
# rename #
ann=ann[,3:ncol(ann)]
#colnames(ann)=c("e.os","t.os","stage","age")
colnames(ann)=c("e.os","t.os","grade","size","meno","neo","horm","chemo","node","bmi","race","age","stage")
ann$bmi=ann$bmi+1
ann$size=ann$size+1
ann["LHC10249","t.os"]=0.1

######################################################################################
# data ready #
exp.data=list(BRCA=as.matrix(data))
ann.data=list(BRCA=as.data.frame(ann))
exp.data.T <- list();
ann.data.T <- list();
exp.data.V <- list();
ann.data.V <- list();

random.iterations <- 1000

# data-structures for multiple training and validation sets
CI.stats <- matrix(
	data = NA, 
	nrow = random.iterations,
	ncol = length(cancer.types),
	dimnames = list(
		1:random.iterations,
		cancer.types
		)
	);

CI.stats.T <- matrix(
	data = NA, 
	nrow = random.iterations,
	ncol = length(cancer.types),
	dimnames = list(
		1:random.iterations,
		cancer.types
		)
	);


CI.P <- CI.stats;

CI.inclusion <- matrix(
	data = 0, 
	#nrow = nrow(data) + 2,
	nrow = nrow(data) +11,
	ncol = length(cancer.types),
	dimnames = list(
		#c(rownames(data), c("age", "stage")),
		c(rownames(data), c("grade","size","meno","neo","horm","chemo","node","bmi","race","age","stage")),
		cancer.types
		)
	);

#########################################################################################
# split into random Training and Validation sets, and train/validate multivariate models
for (i in 1:random.iterations) {

	for(cancer.type in cancer.types) {
	# for(cancer.type in c("BRCA")) {

		cat("\niteration=", i, cancer.type);

		# create a random Training and Validation partition
		tmp.data <- create.training.validation.split(
			exp.data = exp.data[[cancer.type]], 
			ann.data = ann.data[[cancer.type]],
			seed.number = (512 + i - 1)
			);

		exp.data.T[[cancer.type]] <- tmp.data$exp.T;
		ann.data.T[[cancer.type]] <- tmp.data$ann.T;
		exp.data.V[[cancer.type]] <- tmp.data$exp.V;
		ann.data.V[[cancer.type]] <- tmp.data$ann.V;

		# select cancer specific features
		cox.features <- rownames(data);
		names(cox.features) <- cox.features;
		cox.features[names(cox.features)] <- 1


		# estimate UV value
		for (cox.feature in names(cox.features)) {

			# establish risk groups (median dichotomisation)
			risk.groups <- as.numeric(
				exp.data.T[[cancer.type]][cox.feature, ] >
				median(exp.data.T[[cancer.type]][cox.feature, ], na.rm = TRUE)
				);

			# fit Cox model
			coxmodel.fit <- coxph(
				Surv(ann.data.T[[cancer.type]]$t.os, ann.data.T[[cancer.type]]$e.os) ~ risk.groups
				);
			coxmodel.obj <- summary(coxmodel.fit);
			survdiff.obj <- survdiff(
				Surv(ann.data.T[[cancer.type]]$t.os, ann.data.T[[cancer.type]]$e.os) ~ risk.groups
				);

			# store P values for feature selection
			#cox.features[cox.feature] <- (1 - pchisq( survdiff.obj$chisq, df = (length(survdiff.obj$n) - 1) ));
			 cox.features[cox.feature] <- as.numeric(coxmodel.obj$logtest[3])
			}

		# select feature whereby P < 0.01
		cox.features.selected <- names(cox.features)[which(cox.features < 0.01)];

		if(length(cox.features.selected) > 1) {

			# make data structure for glmnet
			x.T <- t(exp.data.T[[cancer.type]][cox.features.selected, , drop=FALSE]);
			x.V <- t(exp.data.V[[cancer.type]][cox.features.selected, , drop=FALSE]);

			# add age, stage in BRCA
			if (cancer.type == "BRCA") {
				x.T <- cbind(x.T, ann.data.T[[cancer.type]][rownames(x.T), c("grade","size","meno","neo","horm","chemo","node","bmi","race","age","stage")]); 
				#x.T <- x.T
				x.V <- cbind(x.V, ann.data.V[[cancer.type]][rownames(x.V), c("grade","size","meno","neo","horm","chemo","node","bmi","race","age","stage")]); 
				#x.V <- x.V

				x.T <- x.T[complete.cases(x.T), ];
				x.V <- x.V[complete.cases(x.V), ];
				}

			y.T <- cbind(
				"time" = ann.data.T[[cancer.type]][rownames(x.T), "t.os"],  
				"status" = ann.data.T[[cancer.type]][rownames(x.T), "e.os"]
				);
			y.V <- cbind(
				"time" = ann.data.V[[cancer.type]][rownames(x.V), "t.os"],  
				"status" = ann.data.V[[cancer.type]][rownames(x.V), "e.os"]
				);

			# start multivariate models with embedded feature selection on pre-selected genes (UV coxph < 0.1)
			# FIT: fit cross validation based glmnet - 10 folds
			cv.glmnet.fit <- cv.glmnet(x = as.matrix(x.T), y = y.T, standardize = FALSE, family = "cox", nfolds = 10);

			# PREDICT: prediction is dependent on 's' parameter which is value of lambdas. The lamba improves
			# training set performance however, the best/robust performance is gained at lamba which has
			# the minimum mean cross validation error (fit$cvm)

			lambda.to.use <- cv.glmnet.fit$lambda.min;
			riskscores.T <- predict(cv.glmnet.fit, newx = as.matrix(x.T), type = "response", s = lambda.to.use);
			riskscores.V <- predict(cv.glmnet.fit, newx = as.matrix(x.V), type = "response", s = lambda.to.use);

			# selected features
			glm.features <- names(which(cv.glmnet.fit$glmnet.fit$beta[, 
				which(cv.glmnet.fit$lambda == lambda.to.use)
				] != 0));
			if (length(glm.features) > 0) {
				CI.inclusion[glm.features, cancer.type] <- CI.inclusion[glm.features, cancer.type] + 1;
				}
			}
		else if(length(cox.features.selected) == 1) {
			riskscores.T <- t(exp.data.T[[cancer.type]][cox.features.selected, , drop=FALSE]);
			riskscores.V <- t(exp.data.V[[cancer.type]][cox.features.selected, , drop=FALSE]);
			CI.inclusion[cox.features.selected, cancer.type] <- CI.inclusion[cox.features.selected, cancer.type] + 1;
			}
		else {
			cat("\nNot Training any model as no feature was selected after univariate cox modelling");
			next;
			}

		# C-Index
		c.ind.obj <- concordance.index(
			x = riskscores.V, 
			surv.time = ann.data.V[[cancer.type]][rownames(riskscores.V), "t.os"], 
			surv.event = ann.data.V[[cancer.type]][rownames(riskscores.V), "e.os"],
			method = "noether", 
			na.rm = TRUE
			);

		c.ind.obj.T <- concordance.index(
			x = riskscores.T, 
			surv.time = ann.data.T[[cancer.type]][rownames(riskscores.T), "t.os"], 
			surv.event = ann.data.T[[cancer.type]][rownames(riskscores.T), "e.os"],
			method = "noether", 
			na.rm = TRUE
			);

		CI.stats[i, cancer.type] <- c.ind.obj$c.index;
		CI.stats.T[i, cancer.type] <- c.ind.obj.T$c.index;

		CI.P[i, cancer.type] <- c.ind.obj$p.value;
		# cat("\nCox = ", summary(coxph.obj)$conf.int[1,1], "\tC.Index = ", c.ind.obj$c.index);
		}
	}



# store CI stats
write.table(
	CI.stats,
	file = paste("GLM_cox__CI__", random.iterations, ".txt", sep = ""),
	row.names = TRUE,
	col.names = NA,
	sep = "\t"
	);
write.table(
	CI.stats.T,
	file = paste("T.GLM_cox__CI__", random.iterations, ".txt", sep = ""),
	row.names = TRUE,
	col.names = NA,
	sep = "\t"
	);

write.table(
	CI.P,
	file = paste("GLM_cox__CI_P__", random.iterations, ".txt", sep = ""),
	row.names = TRUE,
	col.names = NA,
	sep = "\t"
	);

# store inclusion frequencies
write.table(
	CI.inclusion,
	file = paste("GLM_cox__inclusion_frequency_table__", random.iterations, ".txt", sep = ""),
	row.names = TRUE,
	col.names = NA,
	sep = "\t"
	);






























