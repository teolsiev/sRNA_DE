# Introduction----

# Step 6
# This script is used to perform DE-analysis between healthy path_MMR carriers and non-LS control group.

# This script has been written by Tero Sievänen (Uni. of Jyväskylä) and Tia-Marje Korhonen (Uni. of Jyväskylä)
# Edited in 17.11.2021 by Tero Sievänen

# Data setup----

# Choose conditions of interest (healthy LS group, n = 81 and control group, n = 37)
select <- which(targets$Type=="LS" & targets$Healthy_now=="YES" | targets$Type=="CTRL")

# Create a new filtered counts file
Counts <- counts[,select]

# New phenofile with only the variables of interest
Targets <- targets[select,]

# DE-analysis with DESeq2----

# Setup design matrix for DE-analysis
condition <- as.character(Targets$Type) # condition of interest to be tested
batch <- as.character(Targets$NGS) # Batch effect
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

# Create a data frame of the ordered results (top 38)
padj.subset <- head(resOrdered, 40) %>%
 as_tibble(rownames = "miR")

# Create a static gene table
#names(padj.subset)<- c("miR","Mean", "log2FC", "SE", "Wald", "p-value", "FDR")
gt(padj.subset)

# Create a interactive gene table of the subset
datatable(padj.subset, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 2: Differentially expressed c-miRs in healthy path_MMR carriers vs non-LS control group',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 20, lengthMenu = c("5", "10", "15", "50"))) %>%
  formatRound(columns=c(2:7), digits=3)

# Create a data frame of all results for plotting
res.df <- as_tibble(res, rownames = "miR")

# Create a volcano plot of results
v2 <- ggplot(data=res.df,
             aes(y=-log10(res$padj), x=res$log2FoldChange, text = paste(rownames(Counts)))) +
  xlab("Log2FC") +
  ylab("-log10(Padj)") +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  labs(title="Volcano plot",
       subtitle = "LS healthy vs CTRL",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# Create an interactive volcano plot of results
ggplotly(v2)
