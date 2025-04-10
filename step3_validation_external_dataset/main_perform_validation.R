# load the required libraries
library("openxlsx")
library("RColorBrewer")
library("survcomp")
library("survival")
library("survminer")

# load the data and the results
load("processed_data/clinical_data.RData")
load("processed_data/mutations.RData")
load("estimated_model/model_estimate.RData")

# process the survival data
survival_data = clinical_data[,c("SAMPLE_ID","AGE","OS_MONTHS","OS_STATUS","LFS_MONTHS","LFS_STATUS")]
colnames(survival_data) = c("SAMPLE_ID","AGE","OS_MONTHS","OS_STATUS","LFS_MONTHS","LFS_STATUS")
survival_data = survival_data[which(!is.na(survival_data$AGE)),]
survival_data[,"OS_STATUS"][which(as.numeric(survival_data[,"OS_MONTHS"])<1|as.numeric(survival_data[,"AGE"])>90)] = NA
survival_data[,"OS_MONTHS"][which(as.numeric(survival_data[,"OS_MONTHS"])<1|as.numeric(survival_data[,"AGE"])>90)] = NA
survival_data[,"OS_STATUS"][which(as.numeric(survival_data[,"OS_MONTHS"])>60)] = 0
survival_data[,"OS_MONTHS"][which(as.numeric(survival_data[,"OS_MONTHS"])>60)] = 60
survival_data[,"LFS_STATUS"][which(as.numeric(survival_data[,"LFS_MONTHS"])<1|as.numeric(survival_data[,"AGE"])>90)] = NA
survival_data[,"LFS_MONTHS"][which(as.numeric(survival_data[,"LFS_MONTHS"])<1|as.numeric(survival_data[,"AGE"])>90)] = NA
survival_data[,"LFS_STATUS"][which(as.numeric(survival_data[,"LFS_MONTHS"])>60)] = 0
survival_data[,"LFS_MONTHS"][which(as.numeric(survival_data[,"LFS_MONTHS"])>60)] = 60
rownames(survival_data) = survival_data$SAMPLE_ID
survival_data$SAMPLE_ID = NULL
survival_data$AGE = NULL
invalid_samples = unique(which(is.na(survival_data),arr.ind=TRUE)[,"row"])
if(length(invalid_samples)>0) {
    survival_data = survival_data[-invalid_samples,,drop=FALSE]
}

# process the analysis data
model_estimate = model_estimate[which(model_estimate$VARIABLE!="Root_to_ATRX"),]
rownames(model_estimate) = 1:nrow(model_estimate)
analysis_data = matrix(0, nrow = length(rownames(survival_data)), ncol = length(model_estimate$VARIABLE))
rownames(analysis_data) = rownames(survival_data)
colnames(analysis_data) = model_estimate$VARIABLE
IPSSM_SCORE = clinical_data$IPSSM_SCORE
names(IPSSM_SCORE) = clinical_data$SAMPLE_ID
AGE = clinical_data$AGE
names(AGE) = clinical_data$SAMPLE_ID
analysis_data[rownames(analysis_data),"IPSSM_SCORE"] = as.numeric(IPSSM_SCORE[rownames(analysis_data)])
analysis_data[rownames(analysis_data),"AGE"] = as.numeric(AGE[rownames(analysis_data)])
mutations = mutations[rownames(analysis_data),]
analysis_data[names(which(rowSums(mutations[,c("ASXL1","KRAS")])==2)),"ASXL1_to_KRAS"] = 1
analysis_data[names(which(rowSums(mutations[,c("NRAS","RUNX1")])==2)),"NRAS_and_RUNX1"] = 1
analysis_data[names(which(rowSums(mutations[,"JAK2",drop=FALSE])==1)),"Root_to_JAK2"] = 1
analysis_data[names(which(rowSums(mutations[,c("SRSF2","NRAS")])==2)),"SRSF2_to_NRAS"] = 1
analysis_data = data.frame(cbind(survival_data,analysis_data[rownames(survival_data),]))

# evaluate the model on the validation data
model_evaluation = matrix(NA, nrow = 1, ncol = 20)
rownames(model_evaluation) = "Validation"
colnames(model_evaluation) = c("OS ASCETIC c-index","OS ASCETIC c-index lower","OS ASCETIC c-index upper","OS IPSSM c-index","OS IPSSM c-index lower","OS IPSSM c-index upper","OS c-index pvalue","OS ASCETIC AIC","OS IPSSM AIC","OS AIC difference","LFS ASCETIC c-index","LFS ASCETIC c-index lower","LFS ASCETIC c-index upper","LFS IPSSM c-index","LFS IPSSM c-index lower","LFS IPSSM c-index upper","LFS c-index pvalue","LFS ASCETIC AIC","LFS IPSSM AIC","LFS AIC difference")

# compute ASCETIC risk score and clusters
set.seed(12345)
x = as.numeric(model_estimate$COEFFICIENT)
y = analysis_data[,model_estimate$VARIABLE]
beta = as.numeric(x%*%t(y))
risk_score = round(beta,digits=3)
names(risk_score) = rownames(analysis_data)
score_ascetic = as.numeric(risk_score)
clusters_ascetic = as.numeric(kmeans(x=score_ascetic,centers=6,iter.max=1000000000,nstart=1000)$cluster)
clusters_ascetic = paste0("C",clusters_ascetic)
score_ascetic_validation = score_ascetic
clusters_ascetic_validation = clusters_ascetic
score_ipssm = as.numeric(analysis_data$IPSSM_SCORE)

# compute the c-index comparing the two scores considering overall survival
set.seed(12345)
stime = as.numeric(analysis_data$OS_MONTHS)
sevent = as.numeric(analysis_data$OS_STATUS)
c1_score_ipssm = concordance.index(x=score_ipssm,surv.time=stime,surv.event=sevent,method="noether")
c2_score_ascetic = concordance.index(x=score_ascetic,surv.time=stime,surv.event=sevent,method="noether")
comparison_scores = cindex.comp(c2_score_ascetic,c1_score_ipssm)
comparison_scores_validation_os = comparison_scores
os_score_ipssm_lower = as.numeric(c1_score_ipssm$lower)
os_score_ipssm_upper = as.numeric(c1_score_ipssm$upper)
os_score_ascetic_lower = as.numeric(c2_score_ascetic$lower)
os_score_ascetic_upper = as.numeric(c2_score_ascetic$upper)

# compute the c-index comparing the two scores considering leukemia-free survival
set.seed(12345)
stime = as.numeric(analysis_data$LFS_MONTHS)
sevent = as.numeric(analysis_data$LFS_STATUS)
c1_score_ipssm = concordance.index(x=score_ipssm,surv.time=stime,surv.event=sevent,method="noether")
c2_score_ascetic = concordance.index(x=score_ascetic,surv.time=stime,surv.event=sevent,method="noether")
comparison_scores = cindex.comp(c2_score_ascetic,c1_score_ipssm)
comparison_scores_validation_pfs = comparison_scores
pfs_score_ipssm_lower = as.numeric(c1_score_ipssm$lower)
pfs_score_ipssm_upper = as.numeric(c1_score_ipssm$upper)
pfs_score_ascetic_lower = as.numeric(c2_score_ascetic$lower)
pfs_score_ascetic_upper = as.numeric(c2_score_ascetic$upper)

# compute the AIC score comparing the two scores considering overall survival
set.seed(12345)
stime = as.numeric(analysis_data$OS_MONTHS)
sevent = as.numeric(analysis_data$OS_STATUS)
fit_ipssm = coxph(Surv(stime, sevent) ~ score_ipssm)
res1 = extractAIC(fit_ipssm)
fit_ascetic = coxph(Surv(stime, sevent) ~ score_ascetic)
res2 = extractAIC(fit_ascetic)
diff_scores = (res1[2]-res2[2])
diff_scores_validation_os = diff_scores
diff_scores_validation_os_values = c(res1[2],res2[2])

# compute the AIC score comparing the two scores considering leukemia-free survival
set.seed(12345)
stime = as.numeric(analysis_data$LFS_MONTHS)
sevent = as.numeric(analysis_data$LFS_STATUS)
fit_ipssm = coxph(Surv(stime, sevent) ~ score_ipssm)
res1 = extractAIC(fit_ipssm)
fit_ascetic = coxph(Surv(stime, sevent) ~ score_ascetic)
res2 = extractAIC(fit_ascetic)
diff_scores = (res1[2]-res2[2])
diff_scores_validation_pfs = diff_scores
diff_scores_validation_pfs_values = c(res1[2],res2[2])

# save the results
model_evaluation["Validation","OS ASCETIC c-index"] = comparison_scores_validation_os$cindex1
model_evaluation["Validation","OS ASCETIC c-index lower"] = os_score_ascetic_lower
model_evaluation["Validation","OS ASCETIC c-index upper"] = os_score_ascetic_upper
model_evaluation["Validation","OS IPSSM c-index"] = comparison_scores_validation_os$cindex2
model_evaluation["Validation","OS IPSSM c-index lower"] = os_score_ipssm_lower
model_evaluation["Validation","OS IPSSM c-index upper"] = os_score_ipssm_upper
model_evaluation["Validation","OS c-index pvalue"] = comparison_scores_validation_os$p.value
model_evaluation["Validation","OS ASCETIC AIC"] = diff_scores_validation_os_values[2]
model_evaluation["Validation","OS IPSSM AIC"] = diff_scores_validation_os_values[1]
model_evaluation["Validation","OS AIC difference"] = diff_scores_validation_os
model_evaluation["Validation","LFS ASCETIC c-index"] = comparison_scores_validation_pfs$cindex1
model_evaluation["Validation","LFS ASCETIC c-index lower"] = pfs_score_ascetic_lower
model_evaluation["Validation","LFS ASCETIC c-index upper"] = pfs_score_ascetic_upper
model_evaluation["Validation","LFS IPSSM c-index"] = comparison_scores_validation_pfs$cindex2
model_evaluation["Validation","LFS IPSSM c-index lower"] = pfs_score_ipssm_lower
model_evaluation["Validation","LFS IPSSM c-index upper"] = pfs_score_ipssm_upper
model_evaluation["Validation","LFS c-index pvalue"] = comparison_scores_validation_pfs$p.value
model_evaluation["Validation","LFS ASCETIC AIC"] = diff_scores_validation_pfs_values[2]
model_evaluation["Validation","LFS IPSSM AIC"] = diff_scores_validation_pfs_values[1]
model_evaluation["Validation","LFS AIC difference"] = diff_scores_validation_pfs

# perform clustering on the validation dataset
dev.new(height=10,width=15)
clusters_ascetic_validation = rep(NA,length(score_ascetic_validation))
names(clusters_ascetic_validation) = rownames(analysis_data)
clusters_ascetic_validation[which(score_ascetic_validation<=0.50)] = "C1"
clusters_ascetic_validation[which(score_ascetic_validation>0.50&score_ascetic_validation<=1.00)] = "C2"
clusters_ascetic_validation[which(score_ascetic_validation>1.00&score_ascetic_validation<=1.50)] = "C3"
clusters_ascetic_validation[which(score_ascetic_validation>1.50&score_ascetic_validation<=2.00)] = "C4"
clusters_ascetic_validation[which(score_ascetic_validation>2.00&score_ascetic_validation<=3.00)] = "C5"
clusters_ascetic_validation[which(score_ascetic_validation>3.00)] = "C6"
Clusters = clusters_ascetic_validation
print(ggsurvplot(survfit(Surv(as.numeric(analysis_data$OS_MONTHS),as.numeric(analysis_data$OS_STATUS))~Clusters),
        data = analysis_data,
        xlab = "Months",
        ylab = "Overall survival probability",
        palette = "Dark2",
        mark.time = TRUE,
        pval = TRUE,
        risk.table = TRUE,
        ggtheme = theme_bw(),
        title = "Kaplan-Meier (Validation dataset)",
        font.main = 25,
        font.x = 25,
        font.y = 25,
        font.caption = 25,
        font.legend = 25,
        font.tickslab = 25))
print(ggsurvplot(survfit(Surv(as.numeric(analysis_data$LFS_MONTHS),as.numeric(analysis_data$LFS_STATUS))~Clusters),
        data = analysis_data,
        xlab = "Months",
        ylab = "Leukemia-free survival probability",
        palette = "Dark2",
        mark.time = TRUE,
        pval = TRUE,
        risk.table = TRUE,
        ggtheme = theme_bw(),
        title = "Kaplan-Meier (Validation dataset)",
        font.main = 25,
        font.x = 25,
        font.y = 25,
        font.caption = 25,
        font.legend = 25,
        font.tickslab = 25))

# save the results
res = list()
res[["Model Evaluation"]] = model_evaluation
write.xlsx(x=res,file="validation/model_evaluation.xlsx",rowNames=FALSE,colNames=TRUE)
save(model_evaluation,file="validation/model_evaluation.RData")
save(score_ascetic_validation,file="validation/score_ascetic_validation.RData")
save(clusters_ascetic_validation,file="validation/clusters_ascetic_validation.RData")
