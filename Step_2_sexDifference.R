# Introduction ----

# Step 2
# This script is used to perform DE-analysis between sexes within all groups
# There are sex-based differences among men and women and therefore it should be confirmed whether this difference can be seen in c-miR expression.

# This script has been written by Tero Sievänen (Uni. of Jyväskylä) and Tia-Marje Korhonen (Uni. of Jyväskylä)
# Edited in 17.11.2021 by Tero Sievänen

# Data should have been imported already in Step 1 script

# Load packages----

library(tidyverse) # Tidyverse is used f or data wrangling
library(gt) # Tool for data visualization
library(DESeq2) # Tool for DE analysis
library(plotly) # Tool for interactive plots
library(DT) # Tool for interactive tables
library(ggplot2)

# Data import----

# Read in the filtered c-miR counts table from step 1
counts <- read.csv("FilteredCounts.txt",header=TRUE,sep="\t")

# Data setup----

# Add sample names (column names) to raw counts table
colnames(counts) <-targets$Filename

# Choose only the healthy LS individuals (!= cancer | future_ca), n = 86
select <- which(targets$Type=="LS" & targets$Healthy_now=="YES")

# AND/ OR 
#select <- which(targets$Type=="SRME") # No sex difference
#select <- which(targets$Type=="CTRL") # No sex difference

# Create a new filtered counts file 
Counts <- counts[,select]

# New phenofile with only the variables of interest
Targets <- targets[select,]

# DE-analysis of sex difference with DESeq2----

# Setup design matrix for DE-analysis
condition <- as.character(Targets$Sex) # condition of interest (sex)
batch <- as.character(Targets$NGS) # batch effect
group_levels <- levels(as.factor(condition))
design <- data.frame(condition=as.factor(condition), batch=batch) # DE-design for the analysis
rownames(design) <- colnames(Counts)
dds <- DESeqDataSetFromMatrix(countData=Counts, colData=design, design = ~ batch + condition)

# DESeq2 DE-analysis of the condition of interest, batch effect taken into account
dds <- DESeq(dds)

# Display results----

res <- results(dds, alpha=0.05) # Statistical significance at the level p=< 0.05
resOrdered <- res[order(res$padj),] # Order results based on adjusted p-value
summary(res)
resOrdered
mcols(res)$description

# Create a data frame of the ordered results (top 20)
padj.subset <- head(resOrdered, 20) %>%
  as_tibble(rownames = "miR")

# Create a gene table of the subset
#gt(subset)

# Create a interactive gene table of the subset
datatable(padj.subset, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Table 1: Differentially expressed c-miRs between sexes in LS',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 20, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:7), digits=3)

# Create a data frame of all results for plotting
res.df <- as_tibble(res, rownames = "miR")

# Create a volcano plot of results
v1 <- ggplot(data=res.df,
  aes(y=-log10(res$padj), x=res$log2FoldChange, text = paste(rownames(Counts)))) +
  xlab("Log2FC") +
ylab("-log10(Padj)") +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  labs(title="Volcano plot",
       subtitle = "Men vs women in LS group",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# Create an interactive volcano plot of results
ggplotly(v1)
