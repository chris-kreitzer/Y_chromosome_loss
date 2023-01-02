##----------------+
## Subhi Version
##----------------+

library(tidyverse)
library(plyr)
######## Cleaned up Data Freeze####
data_freeze <- read.delim("~/Downloads/DF_CT_type.txt",header = TRUE)
unique(data_freeze$CANCER_TYPE)

data_freeze <- data_freeze %>% 
  dplyr::select(PATIENT_ID, SAMPLE_ID, CANCER_TYPE, AGE_AT_SEQ_REPORTED_YEARS, SAMPLE_TYPE, TUMOR_PURITY,  Y_call , ploidy  ,  purity,FGA,MSI_Type ) %>%
  group_by(CANCER_TYPE, Y_call) %>%
  dplyr::mutate(n_ycall = n()) %>%
  ungroup()



# Get cancer list
cancer_type_list <- unique(as.character(data_freeze$CANCER_TYPE))



# Load mutation status data
load("~/Desktop/work/IMPACT-40K-Working-Group/clinical/data/1.0.genomic_data_all.Rdata")
### For gene level oncogenic (CNA, MUT and fusion)###
GAM_input <-(data_gam_onc$gene)
### For pathway level)
GAM_input <-(data_gam_onc$Oncopath)

GAM_input[GAM_input>1]=1
GAM_input = cbind(rownames(GAM_input),GAM_input)
rownames(GAM_input) <- NULL
colnames(GAM_input)[1]<-"SAMPLE_ID"
GAM_input <- as.data.frame.matrix(GAM_input)

# Set up data frame with additional clinical variables
log_input_df <- data_freeze %>%
  left_join(GAM_input) %>%
  mutate(AGE_AT_SEQ_REPORTED_YEARS = as.numeric(AGE_AT_SEQ_REPORTED_YEARS))

# Change gene factors to numeric
log_input_df[,colnames(GAM_input)[2:ncol(GAM_input)]] <- log_input_df[,colnames(GAM_input)[2:ncol(GAM_input)]] %>%
  mutate_if(is.factor, as.character) %>%
  mutate_if(is.character, as.integer)

# Set up logistic regression function
# Filter for gene then for cancer type
# Run logistic regression if:
## 1. gene has at least 1 mutant and at least 1 wild type
## 2. gene has alteration frequency of >= 0.03
logistic_reg_fun <- function(data_frame, gene, cancer_type){
  data_frame <- data_frame[is.na(data_frame[,gene]) == F &
                             data_frame[,"CANCER_TYPE"] == cancer_type,]
  if (length(which(data_frame[,gene] == 1)) == 0 | 
      length(which(data_frame[,gene] == 0)) == 0){
    log_results_df <- data.frame(variable = "Not Tested",
                                 gene = gene, 
                                 cancer_type = cancer_type, 
                                 comments = "No Mutations in this gene") 
  } else if (length(which(data_frame[,gene] == 1))/nrow(data_frame) < 0.03) {
    log_results_df <- data.frame(variable = "Not Tested",
                                 gene = gene, 
                                 cancer_type = cancer_type, 
                                 comments = "Mutation frequency <3%") 
  } else { 
    
    formula <- as.formula(paste0(gene, "~ Y_call + SAMPLE_TYPE + FGA +MSI_Type"))
    log_results <- glm(formula, data = data_frame, family = binomial)
    log_results_df <- as.data.frame(summary(log_results)$coefficients)
    log_results_df$variable <- row.names(log_results_df)
    log_results_df$gene <- gene
    log_results_df$cancer_type <- cancer_type
    colnames(log_results_df)[1:4] <- c("estimate", "std_err", "z_value", "p_value")
    conf_df <- as.data.frame(confint.default(log_results))
    conf_df$variable <- row.names(conf_df)
    log_results_df <- left_join(log_results_df, conf_df, by = "variable")
    log_results_df <- log_results_df[c(5,7,6,1:4,8,9)]
    print(gene)
    print(cancer_type)
    return(log_results_df)
  }
}

# Set gene list
colnames(log_input_df) <- gsub("-", "_", colnames(log_input_df))
gene_list <- gsub("-", "_", colnames(GAM_input)[2:ncol(GAM_input)])

# Expand table for gene list and cancer list (get every combo)
gene_cancer_df <- expand.grid(gene_list, cancer_type_list)
gene_cancer_df <-  gene_cancer_df %>%
  mutate_if(is.factor, as.character)
gene_cancer_df

# Run logistic regression
log_results_df <- mapply(logistic_reg_fun,
                         gene = gene_cancer_df$Var1,
                         cancer_type = gene_cancer_df$Var2,
                         MoreArgs = list(data_frame = log_input_df),
                         SIMPLIFY = FALSE)
log_results_df <- do.call("rbind.fill", log_results_df)

# Adjust P value
log_results_df$p_adj <- p.adjust(log_results_df$p_value, method = "fdr")
#log_results_df1 <- filter(log_results_df, p_adj < 0.05)
# Save
write.table(log_results_df, "~/Desktop/Projects/Chr_Y/LGR_CT_Sample_Type_FGA_MSI_Type_Oncopath.txt", sep = "\t", row.names = F, quote = F)