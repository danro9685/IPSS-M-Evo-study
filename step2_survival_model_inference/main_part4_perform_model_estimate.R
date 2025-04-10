# load the required libraries
library("openxlsx")
library("RColorBrewer")
library("survcomp")
library("survival")
library("survminer")

# load the data and the results
load("final_data/IPSSM/clinical_data.RData")
ipssr_training = as.numeric(clinical_data$IPSSR_SCORE)
names(ipssr_training) = clinical_data$SAMPLE_ID
clinical_data = NULL
load("final_data/EUMDS/clinical_data.RData")
ipssr_testing = as.numeric(clinical_data$IPSSR_SCORE)
names(ipssr_testing) = clinical_data$SAMPLE_ID
clinical_data = NULL
load("final_results/analysis_data_training.RData")
load("final_results/analysis_data_testing.RData")
load("final_results/bootstrap_model_estimate_training.RData")
thr_sig = 0.30
is.sig = NULL
for(i in 1:nrow(bootstrap_model_estimate_training)) {
    is.sig = c(is.sig,max(bootstrap_model_estimate_training[i,c("POSITIVE","NEGATIVE")]))
}
valid_features = sort(unique(bootstrap_model_estimate_training$FEATURES[which(is.sig>thr_sig)]))
analysis_data_training = analysis_data_training[,which(colnames(analysis_data_training)%in%c("OS_MONTHS","OS_STATUS","LFS_MONTHS","LFS_STATUS",valid_features))]
analysis_data_testing = analysis_data_testing[,which(colnames(analysis_data_testing)%in%c("OS_MONTHS","OS_STATUS","LFS_MONTHS","LFS_STATUS",valid_features))]

# perform model estimation on the training dataset
analysis_data = analysis_data_training

# perform multivariate Cox regression to estimate the final model
set.seed(12345)
analysis_cov = analysis_data[,(5:ncol(analysis_data))]
string_test = paste0("analysis_cov$",colnames(analysis_cov),collapse="+")
string_test = gsub("analysis_cov\\$","",string_test)
time = as.numeric(analysis_data$LFS_MONTHS)
status = as.numeric(analysis_data$LFS_STATUS)
string_test = paste0("Surv(time, status) ~ ",string_test,collapse="")
string_test = as.formula(string_test)
res_cox = coxph(formula = string_test, data = analysis_data[,(5:ncol(analysis_data))])
res_cox = summary(res_cox)
res_cox = res_cox$coefficients[,c("coef","Pr(>|z|)")]
res_cox = res_cox[which(res_cox[,2]<0.05),]
model_estimate = data.frame(VARIABLE=rownames(res_cox),COEFFICIENT=as.numeric(res_cox[,"coef"]))
rownames(model_estimate) = 1:nrow(model_estimate)
analysis_data_training = analysis_data_training[,c("OS_MONTHS","OS_STATUS","LFS_MONTHS","LFS_STATUS",model_estimate$VARIABLE)]
analysis_data_testing = analysis_data_testing[,c("OS_MONTHS","OS_STATUS","LFS_MONTHS","LFS_STATUS",model_estimate$VARIABLE)]

# evaluate the model on the training data
analysis_data = analysis_data_training
model_evaluation = matrix(NA, nrow = 2, ncol = 28)
rownames(model_evaluation) = c("Training","Testing")
colnames(model_evaluation) = c("OS ASCETIC c-index","OS ASCETIC c-index lower","OS ASCETIC c-index upper","OS IPSSM c-index","OS IPSSM c-index lower","OS IPSSM c-index upper","OS c-index pvalue","OS ASCETIC AIC","OS IPSSM AIC","OS AIC difference","LFS ASCETIC c-index","LFS ASCETIC c-index lower","LFS ASCETIC c-index upper","LFS IPSSM c-index","LFS IPSSM c-index lower","LFS IPSSM c-index upper","LFS c-index pvalue","LFS ASCETIC AIC","LFS IPSSM AIC","LFS AIC difference","OS IPSSR c-index","OS IPSSR c-index lower","OS IPSSR c-index upper","OS IPSSR AIC","LFS IPSSR c-index","LFS IPSSR c-index lower","LFS IPSSR c-index upper","LFS IPSSR AIC")

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
score_ascetic_training = score_ascetic
clusters_ascetic_training = clusters_ascetic
score_ipssm = as.numeric(analysis_data$IPSSM_SCORE)

# compute the c-index comparing the two scores considering overall survival
set.seed(12345)
stime = as.numeric(analysis_data$OS_MONTHS)
sevent = as.numeric(analysis_data$OS_STATUS)
c1_score_ipssm = concordance.index(x=score_ipssm,surv.time=stime,surv.event=sevent,method="noether")
c2_score_ascetic = concordance.index(x=score_ascetic,surv.time=stime,surv.event=sevent,method="noether")
comparison_scores = cindex.comp(c2_score_ascetic,c1_score_ipssm)
comparison_scores_training_os = comparison_scores
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
comparison_scores_training_pfs = comparison_scores
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
diff_scores_training_os = diff_scores
diff_scores_training_os_values = c(res1[2],res2[2])

# compute the AIC score comparing the two scores considering leukemia-free survival
set.seed(12345)
stime = as.numeric(analysis_data$LFS_MONTHS)
sevent = as.numeric(analysis_data$LFS_STATUS)
fit_ipssm = coxph(Surv(stime, sevent) ~ score_ipssm)
res1 = extractAIC(fit_ipssm)
fit_ascetic = coxph(Surv(stime, sevent) ~ score_ascetic)
res2 = extractAIC(fit_ascetic)
diff_scores = (res1[2]-res2[2])
diff_scores_training_pfs = diff_scores
diff_scores_training_pfs_values = c(res1[2],res2[2])

# compute the c-index and AIC scores for IPSSR
set.seed(12345)
stime = as.numeric(analysis_data$OS_MONTHS)
sevent = as.numeric(analysis_data$OS_STATUS)
score_ipssr = ipssr_training[rownames(analysis_data)]
os_c_score_ipssr = concordance.index(x=score_ipssr,surv.time=stime,surv.event=sevent,method="noether")
set.seed(12345)
stime = as.numeric(analysis_data$LFS_MONTHS)
sevent = as.numeric(analysis_data$LFS_STATUS)
score_ipssr = ipssr_training[rownames(analysis_data)]
pfs_c_score_ipssr = concordance.index(x=score_ipssr,surv.time=stime,surv.event=sevent,method="noether")
set.seed(12345)
stime = as.numeric(analysis_data$OS_MONTHS)
sevent = as.numeric(analysis_data$OS_STATUS)
score_ipssr = ipssr_training[rownames(analysis_data)]
fit_ipssr = coxph(Surv(stime, sevent) ~ score_ipssr)
os_aic_score_ipssr = extractAIC(fit_ipssr)[2]
set.seed(12345)
stime = as.numeric(analysis_data$LFS_MONTHS)
sevent = as.numeric(analysis_data$LFS_STATUS)
score_ipssr = ipssr_training[rownames(analysis_data)]
fit_ipssr = coxph(Surv(stime, sevent) ~ score_ipssr)
pfs_aic_score_ipssr = extractAIC(fit_ipssr)[2]

# save the results
model_evaluation["Training","OS ASCETIC c-index"] = comparison_scores_training_os$cindex1
model_evaluation["Training","OS ASCETIC c-index lower"] = os_score_ascetic_lower
model_evaluation["Training","OS ASCETIC c-index upper"] = os_score_ascetic_upper
model_evaluation["Training","OS IPSSM c-index"] = comparison_scores_training_os$cindex2
model_evaluation["Training","OS IPSSM c-index lower"] = os_score_ipssm_lower
model_evaluation["Training","OS IPSSM c-index upper"] = os_score_ipssm_upper
model_evaluation["Training","OS c-index pvalue"] = comparison_scores_training_os$p.value
model_evaluation["Training","OS ASCETIC AIC"] = diff_scores_training_os_values[2]
model_evaluation["Training","OS IPSSM AIC"] = diff_scores_training_os_values[1]
model_evaluation["Training","OS AIC difference"] = diff_scores_training_os
model_evaluation["Training","LFS ASCETIC c-index"] = comparison_scores_training_pfs$cindex1
model_evaluation["Training","LFS ASCETIC c-index lower"] = pfs_score_ascetic_lower
model_evaluation["Training","LFS ASCETIC c-index upper"] = pfs_score_ascetic_upper
model_evaluation["Training","LFS IPSSM c-index"] = comparison_scores_training_pfs$cindex2
model_evaluation["Training","LFS IPSSM c-index lower"] = pfs_score_ipssm_lower
model_evaluation["Training","LFS IPSSM c-index upper"] = pfs_score_ipssm_upper
model_evaluation["Training","LFS c-index pvalue"] = comparison_scores_training_pfs$p.value
model_evaluation["Training","LFS ASCETIC AIC"] = diff_scores_training_pfs_values[2]
model_evaluation["Training","LFS IPSSM AIC"] = diff_scores_training_pfs_values[1]
model_evaluation["Training","LFS AIC difference"] = diff_scores_training_pfs
model_evaluation["Training","OS IPSSR c-index"] = os_c_score_ipssr$c.index
model_evaluation["Training","OS IPSSR c-index lower"] = os_c_score_ipssr$lower
model_evaluation["Training","OS IPSSR c-index upper"] = os_c_score_ipssr$upper
model_evaluation["Training","OS IPSSR AIC"] = os_aic_score_ipssr
model_evaluation["Training","LFS IPSSR c-index"] = pfs_c_score_ipssr$c.index
model_evaluation["Training","LFS IPSSR c-index lower"] = pfs_c_score_ipssr$lower
model_evaluation["Training","LFS IPSSR c-index upper"] = pfs_c_score_ipssr$upper
model_evaluation["Training","LFS IPSSR AIC"] = pfs_aic_score_ipssr

# evaluate the model on the testing data
analysis_data = analysis_data_testing

# compute ASCETIC risk score and clusters
set.seed(12345)
x = as.numeric(model_estimate$COEFFICIENT)
y = analysis_data[,model_estimate$VARIABLE]
beta = as.numeric(x%*%t(y))
risk_score = round(beta,digits=2)
names(risk_score) = rownames(analysis_data)
score_ascetic = as.numeric(risk_score)
score_ascetic_testing = score_ascetic
score_ipssm = as.numeric(analysis_data$IPSSM_SCORE)

# compute the c-index comparing the two scores considering overall survival
set.seed(12345)
stime = as.numeric(analysis_data$OS_MONTHS)
sevent = as.numeric(analysis_data$OS_STATUS)
c1_score_ipssm = concordance.index(x=score_ipssm,surv.time=stime,surv.event=sevent,method="noether")
c2_score_ascetic = concordance.index(x=score_ascetic,surv.time=stime,surv.event=sevent,method="noether")
comparison_scores = cindex.comp(c2_score_ascetic,c1_score_ipssm)
comparison_scores_testing_os = comparison_scores
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
comparison_scores_testing_pfs = comparison_scores
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
diff_scores_testing_os = diff_scores
diff_scores_testing_os_values = c(res1[2],res2[2])

# compute the AIC score comparing the two scores considering leukemia-free survival
set.seed(12345)
stime = as.numeric(analysis_data$LFS_MONTHS)
sevent = as.numeric(analysis_data$LFS_STATUS)
fit_ipssm = coxph(Surv(stime, sevent) ~ score_ipssm)
res1 = extractAIC(fit_ipssm)
fit_ascetic = coxph(Surv(stime, sevent) ~ score_ascetic)
res2 = extractAIC(fit_ascetic)
diff_scores = (res1[2]-res2[2])
diff_scores_testing_pfs = diff_scores
diff_scores_testing_pfs_values = c(res1[2],res2[2])

# compute the c-index and AIC scores for IPSSR
set.seed(12345)
stime = as.numeric(analysis_data$OS_MONTHS)
sevent = as.numeric(analysis_data$OS_STATUS)
score_ipssr = ipssr_testing[rownames(analysis_data)]
os_c_score_ipssr = concordance.index(x=score_ipssr,surv.time=stime,surv.event=sevent,method="noether")
set.seed(12345)
stime = as.numeric(analysis_data$LFS_MONTHS)
sevent = as.numeric(analysis_data$LFS_STATUS)
score_ipssr = ipssr_testing[rownames(analysis_data)]
pfs_c_score_ipssr = concordance.index(x=score_ipssr,surv.time=stime,surv.event=sevent,method="noether")
set.seed(12345)
stime = as.numeric(analysis_data$OS_MONTHS)
sevent = as.numeric(analysis_data$OS_STATUS)
score_ipssr = ipssr_testing[rownames(analysis_data)]
fit_ipssr = coxph(Surv(stime, sevent) ~ score_ipssr)
os_aic_score_ipssr = extractAIC(fit_ipssr)[2]
set.seed(12345)
stime = as.numeric(analysis_data$LFS_MONTHS)
sevent = as.numeric(analysis_data$LFS_STATUS)
score_ipssr = ipssr_testing[rownames(analysis_data)]
fit_ipssr = coxph(Surv(stime, sevent) ~ score_ipssr)
pfs_aic_score_ipssr = extractAIC(fit_ipssr)[2]

# save the results
model_evaluation["Testing","OS ASCETIC c-index"] = comparison_scores_testing_os$cindex1
model_evaluation["Testing","OS ASCETIC c-index lower"] = os_score_ascetic_lower
model_evaluation["Testing","OS ASCETIC c-index upper"] = os_score_ascetic_upper
model_evaluation["Testing","OS IPSSM c-index"] = comparison_scores_testing_os$cindex2
model_evaluation["Testing","OS IPSSM c-index lower"] = os_score_ipssm_lower
model_evaluation["Testing","OS IPSSM c-index upper"] = os_score_ipssm_upper
model_evaluation["Testing","OS c-index pvalue"] = comparison_scores_testing_os$p.value
model_evaluation["Testing","OS ASCETIC AIC"] = diff_scores_testing_os_values[2]
model_evaluation["Testing","OS IPSSM AIC"] = diff_scores_testing_os_values[1]
model_evaluation["Testing","OS AIC difference"] = diff_scores_testing_os
model_evaluation["Testing","LFS ASCETIC c-index"] = comparison_scores_testing_pfs$cindex1
model_evaluation["Testing","LFS ASCETIC c-index lower"] = pfs_score_ascetic_lower
model_evaluation["Testing","LFS ASCETIC c-index upper"] = pfs_score_ascetic_upper
model_evaluation["Testing","LFS IPSSM c-index"] = comparison_scores_testing_pfs$cindex2
model_evaluation["Testing","LFS IPSSM c-index lower"] = pfs_score_ipssm_lower
model_evaluation["Testing","LFS IPSSM c-index upper"] = pfs_score_ipssm_upper
model_evaluation["Testing","LFS c-index pvalue"] = comparison_scores_testing_pfs$p.value
model_evaluation["Testing","LFS ASCETIC AIC"] = diff_scores_testing_pfs_values[2]
model_evaluation["Testing","LFS IPSSM AIC"] = diff_scores_testing_pfs_values[1]
model_evaluation["Testing","LFS AIC difference"] = diff_scores_testing_pfs
model_evaluation["Testing","OS IPSSR c-index"] = os_c_score_ipssr$c.index
model_evaluation["Testing","OS IPSSR c-index lower"] = os_c_score_ipssr$lower
model_evaluation["Testing","OS IPSSR c-index upper"] = os_c_score_ipssr$upper
model_evaluation["Testing","OS IPSSR AIC"] = os_aic_score_ipssr
model_evaluation["Testing","LFS IPSSR c-index"] = pfs_c_score_ipssr$c.index
model_evaluation["Testing","LFS IPSSR c-index lower"] = pfs_c_score_ipssr$lower
model_evaluation["Testing","LFS IPSSR c-index upper"] = pfs_c_score_ipssr$upper
model_evaluation["Testing","LFS IPSSR AIC"] = pfs_aic_score_ipssr

# make the survival plots for the training dataset
dev.new(height=10,width=15)
old_clusters_ascetic = clusters_ascetic_training
clusters_ascetic_training[which(old_clusters_ascetic=="C1")] = "C5"
clusters_ascetic_training[which(old_clusters_ascetic=="C2")] = "C4"
clusters_ascetic_training[which(old_clusters_ascetic=="C3")] = "C2"
clusters_ascetic_training[which(old_clusters_ascetic=="C4")] = "C3"
clusters_ascetic_training[which(old_clusters_ascetic=="C5")] = "C1"
clusters_ascetic_training[which(old_clusters_ascetic=="C6")] = "C6"
print(fivenum(score_ascetic_training[which(clusters_ascetic_training=="C1")]))
print(fivenum(score_ascetic_training[which(clusters_ascetic_training=="C2")]))
print(fivenum(score_ascetic_training[which(clusters_ascetic_training=="C3")]))
print(fivenum(score_ascetic_training[which(clusters_ascetic_training=="C4")]))
print(fivenum(score_ascetic_training[which(clusters_ascetic_training=="C5")]))
print(fivenum(score_ascetic_training[which(clusters_ascetic_training=="C6")]))
clusters_ascetic_training_manual = rep(NA,length(score_ascetic_training))
names(clusters_ascetic_training_manual) = rownames(analysis_data_training)
clusters_ascetic_training_manual[which(score_ascetic_training<=0.50)] = "C1"
clusters_ascetic_training_manual[which(score_ascetic_training>0.50&score_ascetic_training<=1.00)] = "C2"
clusters_ascetic_training_manual[which(score_ascetic_training>1.00&score_ascetic_training<=1.50)] = "C3"
clusters_ascetic_training_manual[which(score_ascetic_training>1.50&score_ascetic_training<=2.00)] = "C4"
clusters_ascetic_training_manual[which(score_ascetic_training>2.00&score_ascetic_training<=3.00)] = "C5"
clusters_ascetic_training_manual[which(score_ascetic_training>3.00)] = "C6"
clusters_ascetic_training = clusters_ascetic_training_manual
Clusters = clusters_ascetic_training
print(ggsurvplot(survfit(Surv(as.numeric(analysis_data_training$OS_MONTHS),as.numeric(analysis_data_training$OS_STATUS))~Clusters),
        data = analysis_data_training,
        xlab = "Months",
        ylab = "Overall survival probability",
        palette = "Dark2",
        mark.time = TRUE,
        pval = TRUE,
        risk.table = TRUE,
        ggtheme = theme_bw(),
        title = "Kaplan-Meier (Training dataset)",
        font.main = 25,
        font.x = 25,
        font.y = 25,
        font.caption = 25,
        font.legend = 25,
        font.tickslab = 25))
print(ggsurvplot(survfit(Surv(as.numeric(analysis_data_training$LFS_MONTHS),as.numeric(analysis_data_training$LFS_STATUS))~Clusters),
        data = analysis_data_training,
        xlab = "Months",
        ylab = "Leukemia-free survival probability",
        palette = "Dark2",
        mark.time = TRUE,
        pval = TRUE,
        risk.table = TRUE,
        ggtheme = theme_bw(),
        title = "Kaplan-Meier (Training dataset)",
        font.main = 25,
        font.x = 25,
        font.y = 25,
        font.caption = 25,
        font.legend = 25,
        font.tickslab = 25))

# perform clustering on the testing dataset
dev.new(height=10,width=15)
print(fivenum(score_ascetic_training[which(clusters_ascetic_training=="C1")]))
print(fivenum(score_ascetic_training[which(clusters_ascetic_training=="C2")]))
print(fivenum(score_ascetic_training[which(clusters_ascetic_training=="C3")]))
print(fivenum(score_ascetic_training[which(clusters_ascetic_training=="C4")]))
print(fivenum(score_ascetic_training[which(clusters_ascetic_training=="C5")]))
print(fivenum(score_ascetic_training[which(clusters_ascetic_training=="C6")]))
clusters_ascetic_testing = rep(NA,length(score_ascetic_testing))
names(clusters_ascetic_testing) = rownames(analysis_data_testing)
clusters_ascetic_testing[which(score_ascetic_testing<=0.50)] = "C1"
clusters_ascetic_testing[which(score_ascetic_testing>0.50&score_ascetic_testing<=1.00)] = "C2"
clusters_ascetic_testing[which(score_ascetic_testing>1.00&score_ascetic_testing<=1.50)] = "C3"
clusters_ascetic_testing[which(score_ascetic_testing>1.50&score_ascetic_testing<=2.00)] = "C4"
clusters_ascetic_testing[which(score_ascetic_testing>2.00&score_ascetic_testing<=3.00)] = "C5"
clusters_ascetic_testing[which(score_ascetic_testing>3.00)] = "C6"
Clusters = clusters_ascetic_testing
print(ggsurvplot(survfit(Surv(as.numeric(analysis_data_testing$OS_MONTHS),as.numeric(analysis_data_testing$OS_STATUS))~Clusters),
        data = analysis_data_testing,
        xlab = "Months",
        ylab = "Overall survival probability",
        palette = "Dark2",
        mark.time = TRUE,
        pval = TRUE,
        risk.table = TRUE,
        ggtheme = theme_bw(),
        title = "Kaplan-Meier (Testing dataset)",
        font.main = 25,
        font.x = 25,
        font.y = 25,
        font.caption = 25,
        font.legend = 25,
        font.tickslab = 25))
print(ggsurvplot(survfit(Surv(as.numeric(analysis_data_testing$LFS_MONTHS),as.numeric(analysis_data_testing$LFS_STATUS))~Clusters),
        data = analysis_data_testing,
        xlab = "Months",
        ylab = "Leukemia-free survival probability",
        palette = "Dark2",
        mark.time = TRUE,
        pval = TRUE,
        risk.table = TRUE,
        ggtheme = theme_bw(),
        title = "Kaplan-Meier (Testing dataset)",
        font.main = 25,
        font.x = 25,
        font.y = 25,
        font.caption = 25,
        font.legend = 25,
        font.tickslab = 25))

# save the results
res = list()
res[["Model Estimate"]] = model_estimate
res[["Model Evaluation"]] = model_evaluation
write.xlsx(x=res,file="final_results/model_estimates.xlsx",rowNames=FALSE,colNames=TRUE)
save(model_estimate,file="final_results/model_estimate.RData")
save(model_evaluation,file="final_results/model_evaluation.RData")
save(score_ascetic_training,file="final_results/score_ascetic_training.RData")
save(clusters_ascetic_training,file="final_results/clusters_ascetic_training.RData")
save(score_ascetic_testing,file="final_results/score_ascetic_testing.RData")
save(clusters_ascetic_testing,file="final_results/clusters_ascetic_testing.RData")
