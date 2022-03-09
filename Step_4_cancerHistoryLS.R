# Introduction----

# Step 4
# This script is used to perform DE-analysis between path_MMR carriers with or without previous cancer(s)

# This script has been written by Tero Sievänen (Uni. of Jyväskylä) and Tia-Marje Korhonen (Uni. of Jyväskylä)
# Edited in 17.11.2021 by Tero Sievänen

# # DE-analysis with DESeq2----

# Setup design matrix for DE-analysis
condition <- as.character (Targets$Previous_ca) # Condition of interest (cancer history, 0 = no, 1 = yes)
batch <- as.character (Targets$NGS) # Batch effect
sex <- as.character(Targets$Sex) # Sex used as covariate due to sex difference in LS group
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
