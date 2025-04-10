# load the required libraries
library("ASCETIC")

# set the seed
set.seed(12345)

# get the genes present both in the training and testing datasets
load(file="input_data_training/datasets.RData")
valid_genes = colnames(datasets[[1]])
datasets = NULL
load(file="input_data_testing/datasets.RData")
valid_genes = sort(unique(intersect(valid_genes,colnames(datasets[[1]]))))
datasets = NULL

# load the data for the analysis
load(file="input_data_training/clinical_data.RData")
load(file="input_data_training/datasets.RData")
datasets[[1]] = datasets[[1]][,valid_genes]
load(file="input_data_training/ccf.datasets.RData")
ccf.datasets[[1]] = ccf.datasets[[1]][,valid_genes]
load(file="input_data_training/vaf.datasets.RData")
vaf.datasets[[1]] = vaf.datasets[[1]][which(vaf.datasets[[1]]$GENE_ID%in%valid_genes),]
rownames(vaf.datasets[[1]]) = 1:nrow(vaf.datasets[[1]])

# perform the analysis
all_clinical_data = clinical_data
all_datasets = datasets
all_ccf.datasets = ccf.datasets
all_vaf.datasets = vaf.datasets
all_ascetic_results = list()
all_ascetic_cv = list()
for(i in sort(unique(names(all_clinical_data)))) {
    
    # consider the current cancer type
    dataset = NULL
    ccf.dataset = NULL
    vaf.dataset = NULL
    cancer_type = i
    cat(cancer_type,"\n")
    dataset = all_datasets[[cancer_type]]
    ccf.dataset = all_ccf.datasets[[cancer_type]]
    vaf.dataset = all_vaf.datasets[[cancer_type]]

    # run the ASCETIC inference
    sink("NUL")
    ascetic_results = asceticCCFResampling(dataset=dataset,ccfDataset=ccf.dataset,vafDataset=vaf.dataset,nsampling=100,regularization="aic",command="hc",restarts=50)
    sink()

    # perform a robust estimate via cross-validation
    valid_samples = sort(unique(rownames(dataset)))
    cv_size = round(length(valid_samples)*0.80)
    ranking_estimate = array(NA,c(nrow(ascetic_results$rankingEstimate),100))
    rownames(ranking_estimate) = rownames(ascetic_results$rankingEstimate)
    colnames(ranking_estimate) = paste0("Iteration ",1:ncol(ranking_estimate))
    poset = array(0,c(nrow(ranking_estimate),nrow(ranking_estimate)))
    rownames(poset) = rownames(ranking_estimate)
    colnames(poset) = rownames(ranking_estimate)
    aic = poset
    for(i in 1:100) {
        selected_samples = sort(sample(valid_samples,size=cv_size,replace=FALSE))
        curr_dataset = dataset[selected_samples,,drop=FALSE]
        curr_ccf.dataset = ccf.dataset[selected_samples,,drop=FALSE]
        curr_vaf.dataset = vaf.dataset[which(vaf.dataset$SAMPLE_ID%in%selected_samples),,drop=FALSE]
        rownames(curr_vaf.dataset) = 1:nrow(curr_vaf.dataset)
        sink("NUL")
        cv_ascetic_results = asceticCCFResampling(dataset=curr_dataset,ccfDataset=curr_ccf.dataset,vafDataset=curr_vaf.dataset,nsampling=100,regularization="aic",command="hc",restarts=50)
        sink()
        ranking_estimate[,i] = as.numeric(cv_ascetic_results$rankingEstimate[,"rank"])
        poset = poset + cv_ascetic_results$poset
        aic = aic + cv_ascetic_results$inference$aic
        cat(i/100,"\n")
    }
    ascetic_cv = list()
    ranking_estimate = data.frame(ranking_estimate,check.names=FALSE,stringsAsFactors=FALSE)
    poset = data.frame((poset/ncol(ranking_estimate)),check.names=FALSE,stringsAsFactors=FALSE)
    aic = data.frame((aic/ncol(ranking_estimate)),check.names=FALSE,stringsAsFactors=FALSE)
    ascetic_cv$rankingEstimate = ranking_estimate
    ascetic_cv$poset = poset
    ascetic_cv$aic = aic

    # save the results for the current cancer type
    all_ascetic_results[[cancer_type]] = ascetic_results
    all_ascetic_cv[[cancer_type]] = ascetic_cv

}
ascetic_results = all_ascetic_results
ascetic_cv = all_ascetic_cv

# save the results
save(ascetic_results,file="results/ascetic_results.RData")
save(ascetic_cv,file="results/ascetic_cv.RData")
