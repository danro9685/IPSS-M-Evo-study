### PROCESS THE TRAINING DATASET ###

# load the data
load(file="input_data_training/datasets.RData")
load(file="input_data_training/ccf.datasets.RData")
load(file="results/RankEstimate.RData")
load(file="results/EvoSteps.RData")
load(file="results/EvoStepsValidation.RData")
EvoSteps = EvoSteps[which(EvoSteps$EVO_SCORE>0.50&EvoSteps$CO_OCCURRENCE>=10),]
EvoStepsValidation = EvoStepsValidation[which(EvoStepsValidation$EVO_SCORE>0.50&EvoStepsValidation$CO_OCCURRENCE>=10),]
id_evo_steps = paste(EvoSteps$PARENT,EvoSteps$CHILD,sep="_to_")
id_evo_steps_validation = paste(EvoStepsValidation$PARENT, EvoStepsValidation$CHILD,sep="_to_")
EvoSteps = EvoSteps[which(id_evo_steps%in%id_evo_steps_validation),]
rownames(EvoSteps) = 1:nrow(EvoSteps)
arcs = list()
arcs[["Myelodysplastic Syndromes"]] = EvoSteps
EvoSteps = NULL

# determine the valid co-occurences to be considered
co_occurences = NULL
for(i in sort(unique(RankEstimate$RANK))) {
    curr_genes = sort(unique(RankEstimate$GENE[which(RankEstimate$RANK==i)]))
    for(j in curr_genes) {
        start_g = which(curr_genes==j)
        if(start_g<length(curr_genes)) {
            for(k in (start_g+1):length(curr_genes)) {
                co_occurences = rbind(co_occurences,c(j,curr_genes[k]))
            }
        }
    }
}

# process the data
all_ascetic_features = list()
for(cancer_type in names(arcs)) {

    cat(cancer_type,"\n")

    # get ASCETIC evolutionary steps
    valid_arcs = arcs[[cancer_type]]
    PARENT = as.character(valid_arcs$PARENT)
    CHILD = as.character(valid_arcs$CHILD)
    is.valid = which(PARENT%in%colnames(datasets[[cancer_type]])&CHILD%in%colnames(datasets[[cancer_type]]))
    PARENT = PARENT[is.valid]
    CHILD = CHILD[is.valid]

    # add the arcs to Root
    to_root_child = colnames(datasets[[cancer_type]])
    to_root_parent = rep("Root",length(to_root_child))
    PARENT = c(PARENT,to_root_parent)
    CHILD = c(CHILD,to_root_child)
    curr_arcs = data.frame(PARENT=PARENT,CHILD=CHILD,check.names=FALSE,stringsAsFactors=FALSE)
    curr_datasets = cbind(rep(1,nrow(datasets[[cancer_type]])),datasets[[cancer_type]])
    colnames(curr_datasets)[1] = "Root"

    # build the input features for the survival analysis
    ascetic_features = array(0,c(nrow(curr_datasets),nrow(curr_arcs)))
    rownames(ascetic_features) = rownames(curr_datasets)
    colnames(ascetic_features) = paste0(curr_arcs$PARENT," to ",curr_arcs$CHILD)
    ascetic_features = ascetic_features[sort(unique(rownames(ascetic_features))),]
    ascetic_features = ascetic_features[,sort(unique(colnames(ascetic_features)))]
    for(i in rownames(ascetic_features)) {
        for(j in colnames(ascetic_features)) {
            curr_genes = strsplit(j,split=" to ")[[1]]
            if(sum(curr_datasets[i,curr_genes])==2) {
                ascetic_features[i,j] = 1
            }
        }
    }

    # add the co-occurences ad covariates
    for(cooc in 1:nrow(co_occurences)) {
        cooc_samples = which(rowSums(curr_datasets[,co_occurences[cooc,]])==2)
        if(length(cooc_samples)>0) {
            cooc_samples = sort(unique(intersect(names(cooc_samples),rownames(ascetic_features))))
            if(length(cooc_samples)>0) {
                curr_cooc_entry = rep(0,nrow(ascetic_features))
                names(curr_cooc_entry) = rownames(ascetic_features)
                curr_cooc_entry[cooc_samples] = 1
                ascetic_features = cbind(ascetic_features,curr_cooc_entry)
                colnames(ascetic_features)[ncol(ascetic_features)] = paste0(co_occurences[cooc,1]," and ",co_occurences[cooc,2])
            }
        }
    }

    # consider the features observed in at least 10 samples
    ascetic_features = ascetic_features[,sort(unique(names(which(colSums(ascetic_features)>=10))))]
    all_ascetic_features[[cancer_type]] = data.frame(ascetic_features,check.names=FALSE,stringsAsFactors=FALSE)

}
ascetic_features = all_ascetic_features[[cancer_type]]
ascetic_features_training = ascetic_features

### PROCESS THE TESTING DATASET ###

# load the data
load(file="input_data_testing/datasets.RData")
load(file="input_data_testing/ccf.datasets.RData")
load(file="results/EvoSteps.RData")
load(file="results/EvoStepsValidation.RData")
EvoSteps = EvoSteps[which(EvoSteps$EVO_SCORE>0.50&EvoSteps$CO_OCCURRENCE>=10),]
EvoStepsValidation = EvoStepsValidation[which(EvoStepsValidation$EVO_SCORE>0.50&EvoStepsValidation$CO_OCCURRENCE>=10),]
id_evo_steps = paste(EvoSteps$PARENT,EvoSteps$CHILD,sep="_to_")
id_evo_steps_validation = paste(EvoStepsValidation$PARENT, EvoStepsValidation$CHILD,sep="_to_")
EvoSteps = EvoSteps[which(id_evo_steps%in%id_evo_steps_validation),]
rownames(EvoSteps) = 1:nrow(EvoSteps)
arcs = list()
arcs[["Myelodysplastic Syndromes"]] = EvoSteps
EvoSteps = NULL

# consider each cancer type
all_ascetic_features = list()
for(cancer_type in names(arcs)) {

    cat(cancer_type,"\n")

    # get ASCETIC evolutionary steps
    valid_arcs = arcs[[cancer_type]]
    PARENT = as.character(valid_arcs$PARENT)
    CHILD = as.character(valid_arcs$CHILD)
    is.valid = which(PARENT%in%colnames(datasets[[cancer_type]])&CHILD%in%colnames(datasets[[cancer_type]]))
    PARENT = PARENT[is.valid]
    CHILD = CHILD[is.valid]

    # add the arcs to Root
    to_root_child = colnames(datasets[[cancer_type]])
    to_root_parent = rep("Root",length(to_root_child))
    PARENT = c(PARENT,to_root_parent)
    CHILD = c(CHILD,to_root_child)
    curr_arcs = data.frame(PARENT=PARENT,CHILD=CHILD,check.names=FALSE,stringsAsFactors=FALSE)
    curr_datasets = cbind(rep(1,nrow(datasets[[cancer_type]])),datasets[[cancer_type]])
    colnames(curr_datasets)[1] = "Root"

    # build the input features for the survival analysis
    ascetic_features = array(0,c(nrow(curr_datasets),nrow(curr_arcs)))
    rownames(ascetic_features) = rownames(curr_datasets)
    colnames(ascetic_features) = paste0(curr_arcs$PARENT," to ",curr_arcs$CHILD)
    ascetic_features = ascetic_features[sort(unique(rownames(ascetic_features))),]
    ascetic_features = ascetic_features[,sort(unique(colnames(ascetic_features)))]
    for(i in rownames(ascetic_features)) {
        for(j in colnames(ascetic_features)) {
            curr_genes = strsplit(j,split=" to ")[[1]]
            if(sum(curr_datasets[i,curr_genes])==2) {
                ascetic_features[i,j] = 1
            }
        }
    }

    # add the co-occurences ad covariates
    for(cooc in 1:nrow(co_occurences)) {
        cooc_samples = which(rowSums(curr_datasets[,co_occurences[cooc,]])==2)
        if(length(cooc_samples)>0) {
            cooc_samples = sort(unique(intersect(names(cooc_samples),rownames(ascetic_features))))
            if(length(cooc_samples)>0) {
                curr_cooc_entry = rep(0,nrow(ascetic_features))
                names(curr_cooc_entry) = rownames(ascetic_features)
                curr_cooc_entry[cooc_samples] = 1
                ascetic_features = cbind(ascetic_features,curr_cooc_entry)
                colnames(ascetic_features)[ncol(ascetic_features)] = paste0(co_occurences[cooc,1]," and ",co_occurences[cooc,2])
            }
        }
    }

    # consider the features observed in at least 10 samples
    ascetic_features = ascetic_features[,sort(unique(names(which(colSums(ascetic_features)>=10))))]
    all_ascetic_features[[cancer_type]] = data.frame(ascetic_features,check.names=FALSE,stringsAsFactors=FALSE)

}
ascetic_features = all_ascetic_features[[cancer_type]]
ascetic_features_testing = ascetic_features

# build the final datasets
valid_features = sort(unique(intersect(colnames(ascetic_features_training),colnames(ascetic_features_testing))))
ascetic_features_training = ascetic_features_training[,valid_features]
ascetic_features_testing = ascetic_features_testing[,valid_features]

# save the results
save(ascetic_features_training,file="final_results/ascetic_features_training.RData")
save(ascetic_features_testing,file="final_results/ascetic_features_testing.RData")
