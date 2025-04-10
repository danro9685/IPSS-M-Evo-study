# load the data
load(file="final_data/IPSSM/clinical_data.RData")
load(file="final_data/IPSSM/mutations_data.RData")
load(file="final_data/IPSSM/mutations.RData")
load(file="final_data/common_genes.RData")
mutations = mutations[,common_genes]
mutations_data = mutations_data[which(mutations_data$GENE_NAME%in%common_genes),]
rownames(mutations_data) = 1:nrow(mutations_data)

# handle genes mutated multiple times in the same patient
invalid = NULL
cont = 0
for(i in sort(unique(mutations_data$SAMPLE_ID))) {
    curr_p_muts = mutations_data[which(mutations_data$SAMPLE_ID==i),,drop=FALSE]
    curr_double_muts = which(table(curr_p_muts$GENE_NAME)>1)
    if(length(curr_double_muts)>0) {
        curr_double_muts = names(curr_double_muts)
        for(j in curr_double_muts) {
            curr_double_muts_data = mutations_data[which(mutations_data$SAMPLE_ID==i&mutations_data$GENE_NAME==j),,drop=FALSE]
            keep = which.max(curr_double_muts_data$CCF)[1] # keep the mutation occurred first, i.e., with higher CCF
            invalid = c(invalid,which(mutations_data$SAMPLE_ID==i&mutations_data$GENE_NAME==j)[-keep])
        }
    }
    cont = cont + 1
    cat(cont/length(sort(unique(mutations_data$SAMPLE_ID))),"\n")
}
if(length(invalid)>0) {
    mutations_data = mutations_data[-invalid,,drop=FALSE]
    rownames(mutations_data) = 1:nrow(mutations_data)
}

# perform the analysis
all_clinical_data = list()
all_dataset = list()
all_ccf.dataset = list()
all_vaf.dataset = list()
for(i in "Myelodysplastic Syndromes") {
    
    # consider the current cancer type
    cancer_type = i
    cat(cancer_type,"\n")

    # get the data for the current cancer type
    ct_samples = sort(unique(clinical_data$SAMPLE_ID))
    ct_clinical_data = clinical_data[which(clinical_data$SAMPLE_ID%in%ct_samples),,drop=FALSE]
    rownames(ct_clinical_data) = 1:nrow(ct_clinical_data)
    ct_mutations = mutations_data[which(mutations_data$SAMPLE_ID%in%ct_clinical_data$SAMPLE_ID),,drop=FALSE]
    valid_genes_obs = sort(unique(names(which(table(ct_mutations$GENE_NAME)>=5))))
    ct_mutations = ct_mutations[which(ct_mutations$GENE_NAME%in%valid_genes_obs),,drop=FALSE]

    # process the mutations
    dataset = matrix(0,nrow=length(unique(ct_clinical_data$SAMPLE_ID)),ncol=length(unique(ct_mutations$GENE_NAME)))
    rownames(dataset) = sort(unique(ct_clinical_data$SAMPLE_ID))
    colnames(dataset) = sort(unique(ct_mutations$GENE_NAME))
    ccf.dataset = dataset
    for(j in 1:nrow(dataset)) {
        for(k in 1:ncol(dataset)) {
            curr_p = rownames(dataset)[j]
            curr_g = colnames(dataset)[k]
            curr_pos = which(ct_mutations$SAMPLE_ID==curr_p&ct_mutations$GENE_NAME==curr_g)
            if(length(curr_pos)>0) {
                dataset[j,k] = 1
                ccf.dataset[j,k] = mean(ct_mutations$CCF[curr_pos],na.rm=TRUE)
            }
        }
    }
    vaf.dataset = ct_mutations[,c("SAMPLE_ID","GENE_NAME","REF_COUNT","ALT_COUNT","COPY_NUMBER","NORMAL_PLOIDY","VAF","CCF")]
    colnames(vaf.dataset) = c("SAMPLE_ID","GENE_ID","REF_COUNT","ALT_COUNT","COPY_NUMBER","NORMAL_PLOIDY","VAF_ESTIMATE","CCF_ESTIMATE")
    ct_clinical_data = clinical_data[which(clinical_data$SAMPLE_ID%in%rownames(dataset)),]
    ct_clinical_data = ct_clinical_data[order(ct_clinical_data$SAMPLE_ID),]
    rownames(ct_clinical_data) = 1:nrow(ct_clinical_data)

    # save the results for the current cancer type
    all_clinical_data[[cancer_type]] = ct_clinical_data
    all_dataset[[cancer_type]] = dataset
    all_ccf.dataset[[cancer_type]] = ccf.dataset
    all_vaf.dataset[[cancer_type]] = vaf.dataset

}
clinical_data = all_clinical_data
datasets = all_dataset
ccf.datasets = all_ccf.dataset
vaf.datasets = all_vaf.dataset

# save the results
save(clinical_data,file="input_data_training/clinical_data.RData")
save(datasets,file="input_data_training/datasets.RData")
save(ccf.datasets,file="input_data_training/ccf.datasets.RData")
save(vaf.datasets,file="input_data_training/vaf.datasets.RData")
