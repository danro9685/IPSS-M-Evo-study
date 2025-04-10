# load the required libraries
library("openxlsx")

# load the data
load(file="final_data/EUMDS/mutations_data.RData")
load(file="results/ascetic_results.RData")
load(file="results/ascetic_cv.RData")

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
all_arcs = list()
set.seed(12345)
for(i in names(ascetic_results)) {
    
    # consider the current cancer type
    cancer_type = i
    cat(cancer_type,"\n")

    # get the arcs
    adj = (ascetic_cv[[cancer_type]]$aic * ascetic_results[[cancer_type]]$inference$aic)
    PARENT = rownames(adj)[which(adj>0,arr.ind=TRUE)[,"row"]]
    CHILD = colnames(adj)[which(adj>0,arr.ind=TRUE)[,"col"]]
    valid = which(PARENT%in%mutations_data$GENE_NAME&CHILD%in%mutations_data$GENE_NAME)
    PARENT = PARENT[valid]
    CHILD = CHILD[valid]
    CV_SCORE = adj[which(adj>0,arr.ind=TRUE)]
    CV_SCORE = CV_SCORE[valid]
    arcs = data.frame(PARENT=PARENT,CHILD=CHILD,CV_SCORE=CV_SCORE,check.names=FALSE,stringsAsFactors=FALSE)

    # compute and save the results for the current cancer type
    if(nrow(arcs)>0) {
        # compute confidence measures for each evolutionary step
        evo_steps = arcs
        CO_OCCURRENCE = rep(NA, nrow(evo_steps))
        VALID_MODELS = rep(NA, nrow(evo_steps))
        PVALUE = rep(NA, nrow(evo_steps))
        FIRST_PARENT = rep(NA, nrow(evo_steps))
        FIRST_CHILD = rep(NA, nrow(evo_steps))
        DISCONNECTED = rep(NA, nrow(evo_steps))
        for(j in 1:nrow(evo_steps)) {
            curr_p = evo_steps[j,"PARENT"]
            curr_c = evo_steps[j,"CHILD"]
            curr_co = sort(unique(names(which(table(mutations_data[which(mutations_data$GENE_NAME%in%c(curr_p,curr_c)),"SAMPLE_ID"])==2))))
            CO_OCCURRENCE[j] = length(curr_co)
            valid_models = sort(unique(names(which(table(mutations_data[which(mutations_data$GENE_NAME%in%c(curr_p,curr_c)),"SAMPLE_ID"])==2))))
            VALID_MODELS[j] = length(valid_models)
            if(length(valid_models)>=3) { # if I have at least 3 valid sample
                test_sig_p = mutations_data[which(mutations_data$SAMPLE_ID%in%valid_models&mutations_data$GENE_NAME==curr_p),"CCF"]
                test_sig_c = mutations_data[which(mutations_data$SAMPLE_ID%in%valid_models&mutations_data$GENE_NAME==curr_c),"CCF"]
                sig_pval = wilcox.test(x=test_sig_p,y=test_sig_c,alternative="greater")$p.value
                is_p = 0
                is_c = 0
                is_disc = 0
                for(k in valid_models) {
                    curr_mut_p = mutations_data[which(mutations_data$SAMPLE_ID==k&mutations_data$GENE_NAME==curr_p),,drop = FALSE]
                    curr_mut_c = mutations_data[which(mutations_data$SAMPLE_ID==k&mutations_data$GENE_NAME==curr_c),,drop = FALSE]
                    curr_support_p = length(which(curr_mut_p$CCF>curr_mut_c$CCF))
                    curr_support_c = length(which(curr_mut_c$CCF>curr_mut_p$CCF))
                    if(curr_support_p>0) {
                        is_p = is_p + 1
                    }
                    if(curr_support_c>0) {
                        is_c = is_c + 1
                    }
                    if(curr_support_p==0&&curr_support_c==0) {
                        is_disc = is_disc + 1
                    }
                }
                PVALUE[j] = sig_pval
                FIRST_PARENT[j] = is_p
                FIRST_CHILD[j] = is_c
                DISCONNECTED[j] = is_disc
            }
        }
        evo_steps$CO_OCCURRENCE = CO_OCCURRENCE
        evo_steps$VALID_MODELS = VALID_MODELS
        evo_steps$PVALUE = PVALUE
        evo_steps$FIRST_PARENT = FIRST_PARENT
        evo_steps$FIRST_CHILD = FIRST_CHILD
        evo_steps$UNCERTAIN = DISCONNECTED
        evo_steps$EVO_SCORE = (FIRST_PARENT/VALID_MODELS)
        evo_steps = evo_steps[order(evo_steps$EVO_SCORE,evo_steps$FIRST_PARENT,evo_steps$CV_SCORE,evo_steps$CO_OCCURRENCE,decreasing=TRUE),]
        rownames(evo_steps) = 1:nrow(evo_steps)
        evo_steps = evo_steps[,c(1,2,10,3,6,4,7,8,9)]
        all_arcs[[cancer_type]] = evo_steps
    }

}
EvoStepsValidation = all_arcs[[1]]
EvoStepsValidation = EvoStepsValidation[which(!is.na(EvoStepsValidation$EVO_SCORE)),]
rownames(EvoStepsValidation) = 1:nrow(EvoStepsValidation)

# save the results
save(EvoStepsValidation,file="results/EvoStepsValidation.RData")
write.xlsx(EvoStepsValidation,file="results/EvoStepsValidation.xlsx",sheetName="Myelodysplastic Syndromes",rowNames=FALSE,colNames=TRUE)
