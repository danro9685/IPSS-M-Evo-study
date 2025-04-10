# load the required libraries
library("glmnet")
library("openxlsx")
library("survival")
library("survminer")

# load the data
load("final_results/analysis_data_training.RData")
analysis_data = analysis_data_training[,-c(1,2)]

# perform model estimation via bootstrap
n_boot = 100
significant_features = matrix(NA, nrow = n_boot, ncol = length(colnames(analysis_data)[-c(1,2)]))
rownames(significant_features) = paste0("Iteration ",1:n_boot)
colnames(significant_features) = colnames(analysis_data)[-c(1,2)]
for(boot in 1:n_boot) {

    # set the seed
    set.seed(boot)

    # perform the analysis
    boot_analysis_data = analysis_data[sample(1:nrow(analysis_data),size=nrow(analysis_data),replace=TRUE),]
    analysis_cov = boot_analysis_data[,colnames(boot_analysis_data)[3:ncol(boot_analysis_data)]]
    string_test = paste0("~",paste0("analysis_cov$",colnames(analysis_cov),collapse="+"))
    x = model.matrix(as.formula(string_test))
    y = Surv(as.numeric(boot_analysis_data$LFS_MONTHS),as.numeric(boot_analysis_data$LFS_STATUS))
    cv.fit = cv.glmnet(x,y,family="cox",maxit=1000000000,alpha=1)
    coeff_values = coef(cv.fit,s=cv.fit$lambda.min)
    rownames(coeff_values) = gsub("analysis_cov\\$","",rownames(coeff_values))
    rownames(coeff_values) = gsub("\\."," ",rownames(coeff_values))
    invalid = which(rownames(coeff_values)=="(Intercept)")
    if(length(invalid)>0) {
        coeff_values = coeff_values[-invalid,,drop=FALSE]
    }
    coeff_names = rownames(coeff_values)
    coeff_values = as.numeric(coeff_values)
    names(coeff_values) = coeff_names

    # save the results
    significant_features[boot,names(coeff_values)] = as.numeric(coeff_values)
    
    cat(boot/n_boot,"\n")

}

# make some statistics
FEATURES = as.character(colnames(significant_features))
COEFFICIENTS = as.numeric((colSums(significant_features)/n_boot))
POSITIVE = apply(X=significant_features,MARGIN=2,FUN=function(x){length(which(x>0))})
POSITIVE = (POSITIVE/n_boot)
NEGATIVE = apply(X=significant_features,MARGIN=2,FUN=function(x){length(which(x<0))})
NEGATIVE = (NEGATIVE/n_boot)
NEUTRAL = apply(X=significant_features,MARGIN=2,FUN=function(x){length(which(x==0))})
NEUTRAL = (NEUTRAL/n_boot)
model_estimate = data.frame(FEATURES=FEATURES,COEFFICIENTS=COEFFICIENTS,POSITIVE=POSITIVE,NEGATIVE=NEGATIVE,NEUTRAL=NEUTRAL)
model_estimate = model_estimate[order(model_estimate$NEUTRAL),]
rownames(model_estimate) = 1:nrow(model_estimate)
bootstrap_model_estimate_training = model_estimate

# save the results
save(bootstrap_model_estimate_training,file="final_results/bootstrap_model_estimate_training.RData")
model_estimate = list()
model_estimate[["Training"]] = bootstrap_model_estimate_training
write.xlsx(x=model_estimate,file="final_results/model_estimates_bootstrap_training.xlsx",rowNames=FALSE,colNames=TRUE)
