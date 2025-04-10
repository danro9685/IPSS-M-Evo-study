# load the required libraries
library("openxlsx")

# load the data
load(file="results/ascetic_results.RData")
load(file="results/ascetic_cv.RData")

# perform the analysis
all_ranks = list()
for(i in names(ascetic_results)) {
    
    # consider the current cancer type
    cancer_type = i
    cat(cancer_type,"\n")

    # get the ranking estimate
    rank_estimate = ascetic_results[[cancer_type]][["rankingEstimate"]]
    rank_cv = (rowSums(ascetic_cv[[cancer_type]][["rankingEstimate"]])/ncol(ascetic_cv[[cancer_type]][["rankingEstimate"]]))
    rank_estimate = cbind(rownames(rank_estimate),rank_estimate[,2],rank_cv[rownames(rank_estimate)])
    colnames(rank_estimate) = c("GENE","RANK","RANK_CV")
    rank_estimate = rank_estimate[order(rank_estimate[,"RANK_CV"]),]
    rownames(rank_estimate) = 1:nrow(rank_estimate)
    rank_estimate = data.frame(rank_estimate,check.names=FALSE,stringsAsFactors=FALSE)
    rank_estimate$RANK = as.numeric(rank_estimate$RANK)
    rank_estimate$RANK_CV = as.numeric(rank_estimate$RANK_CV)
    rownames(rank_estimate) = 1:nrow(rank_estimate)

    # save the results
    all_ranks[[cancer_type]] = rank_estimate

}
RankEstimate = all_ranks[[1]]

# save results
save(RankEstimate,file="results/RankEstimate.RData")
write.xlsx(RankEstimate,file="results/RankEstimate.xlsx",sheetName="Myelodysplastic Syndromes",rowNames=FALSE,colNames=TRUE)
