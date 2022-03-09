# Introduction----

# Step 3
# This script is used to perform DE-analysis between MLH1 and other variants (MSH2, MSH6 and PMS2) among LS subgroup
# Path_MMR variants differ from each other e.g. in their cancer predisposition risk, and therefore it is necessary to explore whether different path_MMR variant carriers have a detectable c-miR expression profile.

# This script has been written by Tero Sievänen (Uni. of Jyväskylä) and Tia-Marje Korhonen (Uni. of Jyväskylä)
# Edited in 17.11.2021 by Tero Sievänen

# Data setup----

# Choose only the healthy LS individuals (!= cancer | future_ca), n = 86
select <- which(targets$Type=="LS" & targets$Healthy_now=="YES")

# Create a new filtered counts file 
Counts <- counts[,select]

# New phenofile with only the variables of interest
Targets <- targets[select,]

# DE-analysis of path_MMR variants with DESeq2----

# Setup design matrix for DE-analysis
condition <- as.character(Targets$step_3) # condition of interest (path_MMR variant)
batch <- as.character (Targets$NGS) # Batch effect (in which NGS experiment the data was produced)
sex <- as.character(Targets$Sex) # Sex as a covariate due to sex difference in LS group
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
subset <- head(resOrdered, 20) %>%
  as_tibble(rownames = "miR")

# Create a gene table of the subset
gt(subset)
