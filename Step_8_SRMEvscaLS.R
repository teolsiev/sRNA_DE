# Introduction----

# Step 8
# This script is used to perform DE-analysis between sporadic rectal cancer patients (SRME) and path_MMR carriers with cancer.

# This script has been written by Tero Sievänen (Uni. of Jyväskylä) and Tia-Marje Korhonen (Uni. of Jyväskylä)
# Edited in 17.11.2021 by Tero Sievänen

# Data setup----

# Choose conditions of interest (path_MMR with cancer, n = 13 and SRME, n = 24)
select <- which(targets$Type=="LS" & targets$Healthy_now=="NO" | targets$Type=="SRME")

# Create a new filtered counts file
Counts <- counts[,select]

# New phenofile with only the variables of interest
Targets <- targets[select,]

# DE-analysis with DESeq2----

# Setup design matrix for DE-analysis
condition <- as.character (Targets$Type) # condition of interest to be tested
batch <- as.character (Targets$NGS) # batch effect
sex <- as.character(Targets$Sex) # Sex as a covariate
group_levels <- levels(as.factor(condition))
design <- data.frame(condition=as.factor(condition), batch=batch, sex=sex)
rownames(design) <- colnames(Counts)
dds <- DESeqDataSetFromMatrix(countData=Counts, colData=design, design = ~ batch + sex + condition)

# DESeq2 DE-analysis of the condition of interest, batch effect taken into account, sex as a covariate
dds <- DESeq(dds)

# Display results----
res <- results(dds, alpha=0.05)
resOrdered <- res[order(res$padj),]
summary(res)
resOrdered
mcols(res)$description

# Create a data frame of the ordered results (top 20)
subset<- head(resOrdered, 20) %>%
  as_tibble(rownames = "miR")

# Create a gene table of the subset
# Create a static gene table

#names(subset)<- c("miR","Mean", "log2FC", "SE", "Wald", "p-value", "FDR")
gt(subset)
