random.iterations <- 10;
survival_data$AGE[which(survival_data$AGE <= 50)] <- 0;
survival_data$AGE[which(survival_data$AGE > 50)] <- 1;
data=as.matrix(data)

CI.stats <- matrix(
	data = NA, 
	nrow = random.iterations,
	ncol = 1);
CI.P <- CI.stats;
CI.inclusion <- matrix(
	data = 0, 
	nrow = nrow(data) + 1,
	ncol = 1);

# split into random Training and Validation sets, and train/validate multivariate models
for (i in 1:random.iterations) {	
	

		cat("\niteration=", i);

		# create a random Training and Validation partition
		tmp.data <- create.training.validation.split(
			exp.data = data, 
			ann.data = survival_data,
			seed.number = (512 + i - 1)
			);

		exp.data.T <- tmp.data$exp.T;
		ann.data.T <- tmp.data$ann.T;
		exp.data.V <- tmp.data$exp.V;
		ann.data.V <- tmp.data$ann.V;

		# select cancer specific features
		cox.features <- rownames(data);
		names(cox.features) <- cox.features;
		cox.features[names(cox.features)] <- 1;

		# estimate UV value
		for (cox.feature in names(cox.features)) {

			# establish risk groups (median dichotomisation)
			risk.groups <- as.numeric(
				exp.data.T[cox.feature, ] >
				median(as.numeric(exp.data.T[cox.feature, ]), na.rm = TRUE)
				);

			# fit Cox model
			coxmodel.fit <- coxph(
				Surv(ann.data.T$SurvivalTime, ann.data.T$Death) ~ risk.groups
				);
			coxmodel.obj <- summary(coxmodel.fit);
			survdiff.obj <- survdiff(
				Surv(ann.data.T$SurvivalTime, ann.data.T$Death) ~ risk.groups
				);

			# store P values for feature selection
			cox.features[cox.feature] <- (1 - pchisq( survdiff.obj$chisq, df = (length(survdiff.obj$n) - 1) ));
			}

		# select feature whereby P < 0.1
		cox.features.selected <- names(cox.features)[which(cox.features < 0.05)];

		if(length(cox.features.selected) > 1) {

			# make data structure for glmnet
			x.T <- t(exp.data.T[cox.features.selected, , drop=FALSE]);
			x.V <- t(exp.data.V[cox.features.selected, , drop=FALSE]);

			# add age, stage 
							
				x.T <- cbind(x.T, ann.data.T[rownames(x.T), c("AGE")]); 
				x.V <- cbind(x.V, ann.data.V[rownames(x.V), c("AGE")]); 
				x.T <- x.T[complete.cases(x.T), ];
				x.V <- x.V[complete.cases(x.V), ];
				

			y.T <- cbind(
				"time" = ann.data.T[rownames(x.T), "SurvivalTime"],  
				"status" = ann.data.T[rownames(x.T), "Death"]
				);
			y.V <- cbind(
				"time" = ann.data.V[rownames(x.V), "SurvivalTime"],  
				"status" = ann.data.V[rownames(x.V), "Death"]
				);

			# start multivariate models with embedded feature selection on pre-selected genes (UV coxph < 0.05)
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
				CI.inclusion[glm.features] <- CI.inclusion[glm.features] + 1;
				}
			}
		else if(length(cox.features.selected) == 1) {
			riskscores.T <- t(exp.data.T[cox.features.selected, , drop=FALSE]);
			riskscores.V <- t(exp.data.V[cox.features.selected, , drop=FALSE]);
			CI.inclusion[cox.features.selected] <- CI.inclusion[cox.features.selected] + 1;
			}
		else {
			cat("\nNot Training any model as no feature was selected after univariate cox modelling");
			next;
			}

		# C-Index
		c.ind.obj <- concordance.index(
			x = riskscores.V, 
			surv.time = ann.data.V[rownames(riskscores.V), "SurvivalTime"], 
			surv.event = ann.data.V[rownames(riskscores.V), "Death"],
			method = "noether", 
			na.rm = TRUE
			);
		CI.stats[i] <- c.ind.obj$c.index;
		CI.P[i] <- c.ind.obj$p.value;
		# cat("\nCox = ", summary(coxph.obj)$conf.int[1,1], "\tC.Index = ", c.ind.obj$c.index);
		}
	}






# store CI stats
write.table(
	CI.stats,
	file = paste(data.dir, "pan_cancer_GLM_cox__CI__", random.iterations, ".txt", sep = ""),
	row.names = TRUE,
	col.names = NA,
	sep = "\t"
	);

write.table(
	CI.P,
	file = paste(data.dir, "pan_cancer_GLM_cox__CI_P__", random.iterations, ".txt", sep = ""),
	row.names = TRUE,
	col.names = NA,
	sep = "\t"
	);

# store inclusion frequencies
write.table(
	CI.inclusion,
	file = paste(data.dir, "pan_cancer_GLM_cox__inclusion_frequency_table__", random.iterations, ".txt", sep = ""),
	row.names = TRUE,
	col.names = NA,
	sep = "\t"
	);

### SAVE SESSION ##############################################################
sink(file = paste(Sys.Date(), "-Session-Info-assess.prognostic.value.txt", sep = ""), type = c("output", "message"));
print(sessionInfo());
sink(NULL);
