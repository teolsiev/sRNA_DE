# Introduction----

# The aim of the analysis is to investigate whether there are differences in circulating microRNA (c-miR) expression among (healthy) people with Lynch syndrome (LS) and between LS and healthy non-lS controls and sporadic rectal cancer patients.
# Since several traits (sex, path_MMR variant, disease history) generate heterogeneity in LS, and can affect on c-miR expression, the first step is to confirm if our study population can be treated as a single group in analyses based on the c-miR expression profile.

# Step 1
# This script is used to make a filtered c-miR raw count file for downstream analyses
# Genes with low counts must be removed before DE-analysis, since their biological significance is low or they might be false positives

# This script has been written by Tia-Marje Korhonen (Uni. of Jyväskylä)
# Edited in 17.11.2021 by Tero Sievänen (Uni. of Jyväskylä)

# Correspondence to: tero.o.sievanen@jyu.fi

# Load packages----
library(edgeR) # Tool for DE-analysis

# Data import----

# Import raw counts file
counts <- read.csv("rawCounts.tsv",header=TRUE,sep="\t")

# Import study design file
targets<- read.csv("phenodata_age.txt",sep="\t",header=TRUE)

# Add column names
colnames(counts) <-targets$Filename

# Filtering of raw c-miR counts----

# Filter phenodata based on sample type
Group1  <- which(targets$Type=="LS")
Group2  <- which(targets$Type=="CTRL")
Group3  <- which(targets$Type=="SRME")
ls <- targets[Group1,]
ctrl <- targets[Group2,]
srme <- targets[Group3,]

# Create a digital gene expression list (DGEList) using edgeR
myDGEList <- DGEList(counts = counts)

# Convert DGEList to counts per million (cpm)
cpm <- cpm(myDGEList)

# Keep only the c-miRs with >1 CPM in at least 70% of samples in a subgroup
# keepers <- rowSums(cpm>1)>=108  #70% of 155

# Filter the subgroups (LS,CTRL and SRME)
cpm1 <- cpm[,which(colnames(cpm)%in%ls$Sample)]
cpm2 <- cpm[,which(colnames(cpm)%in%ctrl$Sample)]
cpm3 <- cpm[,which(colnames(cpm)%in%srme$Sample)]

keep1 <- rowSums(cpm1>=1) >= 65   #  70% of 94
keep2 <- rowSums(cpm2>=1) >= 25   #  70% of 37
keep3<- rowSums(cpm3>=1) >= 16   #  70% of 24

# Make a new filtered raw c-miR counts file----

# Create a new filtered DGEList
myDGEList.filtered <- myDGEList[(keep1|keep2|keep3),] # OR myDGEList.filtered <- myDGEList[(keepers|keep1|keep2|keep3),]

# Create a new filtered raw counts matrix
FilteredCounts <- counts[which(rownames(counts)%in%rownames(myDGEList.filtered$counts)),]

# Create a new .txt file of filtered raw miR counts
write.table(FilteredCounts,"FilteredCounts.txt",sep="\t")
