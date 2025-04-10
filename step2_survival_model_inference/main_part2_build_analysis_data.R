### PROCESS THE TRAINING DATASET ###

# load and process the input data
load("final_data/IPSSM/clinical_data.RData")
load("final_data/IPSSM/ipss_mol_data.RData")
load("final_results/ascetic_features_training.RData")
ascetic_features = ascetic_features_training
colnames(ascetic_features) = gsub(" ","_",colnames(ascetic_features))
invalid_samples = unique(which(is.na(ipss_mol_data),arr.ind=TRUE)[,"row"])
ipss_mol_data = ipss_mol_data[-invalid_samples,]
ipss_mol_data = ipss_mol_data[,-which(colnames(ipss_mol_data)%in%gsub("Root_to_","",colnames(ascetic_features)))]
ipss_mol_data = ipss_mol_data[,-which(colnames(ipss_mol_data)=="FLT3_ITD")]

# process the survival data
survival_data = clinical_data[,c("SAMPLE_ID","AGE","OS_MONTHS","OS_STATUS","LFS_MONTHS","LFS_STATUS")]
colnames(survival_data) = c("SAMPLE_ID","AGE","OS_MONTHS","OS_STATUS","LFS_MONTHS","LFS_STATUS")
survival_data = survival_data[which(!is.na(survival_data$AGE)),]
survival_data[,"OS_STATUS"][which(survival_data[,"OS_STATUS"]=="0:LIVING")] = 0
survival_data[,"OS_STATUS"][which(survival_data[,"OS_STATUS"]=="1:DECEASED")] = 1
survival_data[,"OS_STATUS"][which(as.numeric(survival_data[,"OS_MONTHS"])<1|as.numeric(survival_data[,"AGE"])>90)] = NA
survival_data[,"OS_MONTHS"][which(as.numeric(survival_data[,"OS_MONTHS"])<1|as.numeric(survival_data[,"AGE"])>90)] = NA
survival_data[,"OS_STATUS"][which(as.numeric(survival_data[,"OS_MONTHS"])>60)] = 0
survival_data[,"OS_MONTHS"][which(as.numeric(survival_data[,"OS_MONTHS"])>60)] = 60
survival_data[,"LFS_STATUS"][which(survival_data[,"LFS_STATUS"]=="0:LeukemiaFree")] = 0
survival_data[,"LFS_STATUS"][which(survival_data[,"LFS_STATUS"]=="1:Transformed/Deceased")] = 1
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
valid_samples = sort(unique(intersect(rownames(ascetic_features),rownames(ipss_mol_data))))
valid_samples = sort(unique(intersect(valid_samples,rownames(survival_data))))
survival_data = survival_data[valid_samples,]
ascetic_features_curr = ascetic_features[rownames(survival_data),,drop=FALSE]
ascetic_features_curr = ascetic_features_curr[,sort(unique(names(which(colSums(ascetic_features_curr)>=10))))]
ipssmol_curr = ipss_mol_data[rownames(survival_data),,drop=FALSE]
analysis_data_training = cbind(survival_data,ipssmol_curr,ascetic_features_curr)

### PROCESS THE TESTING DATASET ###

# load the data
load("final_data/EUMDS/clinical_data.RData")
load("final_data/EUMDS/ipss_mol_data.RData")
load("final_results/ascetic_features_testing.RData")
ascetic_features = ascetic_features_testing
colnames(ascetic_features) = gsub(" ","_",colnames(ascetic_features))
ipss_mol_data$TP53_LOH = NULL
invalid_samples = unique(which(is.na(ipss_mol_data),arr.ind=TRUE)[,"row"])
ipss_mol_data = ipss_mol_data[-invalid_samples,]
ipss_mol_data = ipss_mol_data[,-which(colnames(ipss_mol_data)%in%gsub("Root_to_","",colnames(ascetic_features)))]
ipss_mol_data = ipss_mol_data[,-which(colnames(ipss_mol_data)%in%c("MLL_PTD","ETNK1","GNB1","PPM1D","PRPF8"))]

# process the survival data
survival_data = clinical_data[,c("SAMPLE_ID","AGE","OS_MONTHS","OS_STATUS","LFS_MONTHS","LFS_STATUS")]
colnames(survival_data) = c("SAMPLE_ID","AGE","OS_MONTHS","OS_STATUS","LFS_MONTHS","LFS_STATUS")
survival_data = survival_data[which(!is.na(survival_data$AGE)),]
survival_data[,"OS_STATUS"][which(survival_data[,"OS_STATUS"]=="0:LIVING")] = 0
survival_data[,"OS_STATUS"][which(survival_data[,"OS_STATUS"]=="1:DECEASED")] = 1
survival_data[,"OS_STATUS"][which(as.numeric(survival_data[,"OS_MONTHS"])<1|as.numeric(survival_data[,"AGE"])>90)] = NA
survival_data[,"OS_MONTHS"][which(as.numeric(survival_data[,"OS_MONTHS"])<1|as.numeric(survival_data[,"AGE"])>90)] = NA
survival_data[,"OS_STATUS"][which(as.numeric(survival_data[,"OS_MONTHS"])>60)] = 0
survival_data[,"OS_MONTHS"][which(as.numeric(survival_data[,"OS_MONTHS"])>60)] = 60
survival_data[,"LFS_STATUS"][which(survival_data[,"LFS_STATUS"]=="0:LeukemiaFree")] = 0
survival_data[,"LFS_STATUS"][which(survival_data[,"LFS_STATUS"]=="1:Transformed/Deceased")] = 1
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
valid_samples = sort(unique(intersect(rownames(ascetic_features),rownames(ipss_mol_data))))
valid_samples = sort(unique(intersect(valid_samples,rownames(survival_data))))
survival_data = survival_data[valid_samples,]
ascetic_features_curr = ascetic_features[rownames(survival_data),,drop=FALSE]
ascetic_features_curr = ascetic_features_curr[,sort(unique(names(which(colSums(ascetic_features_curr)>=10))))]
ipssmol_curr = ipss_mol_data[rownames(survival_data),,drop=FALSE]
analysis_data_testing = cbind(survival_data,ipssmol_curr,ascetic_features_curr)

# build the final datasets
valid_features = colnames(analysis_data_training)[which(colnames(analysis_data_training)%in%colnames(analysis_data_testing))]
analysis_data_training = analysis_data_training[,valid_features]
analysis_data_testing = analysis_data_testing[,valid_features]

# save the results
save(analysis_data_training,file="final_results/analysis_data_training.RData")
save(analysis_data_testing,file="final_results/analysis_data_testing.RData")
