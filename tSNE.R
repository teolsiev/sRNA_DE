# Introduction----
# This script is used to perform dimension reduction analysis with tSNE

# This script has been written by Tia-Marje Korhonen (Uni. of Jyväskylä)
# Edited by Tero Sievänen (Uni. of Jyväskylä)

# Correspondence to: tero.o.sievanen@jyu.fi

# Load libraries----
# Load libraries
library(Rtsne)
library(DESeq2)
library(cowplot)
library(ggplot2)
library(tidyverse)
library(gridGraphics)

# Data setup----
# Import count data
counts <- read.csv("FilteredCounts.txt",header=TRUE,sep="\t") # Filtered counts, this file must be generated with step_1 script

# Import phenotype data
targets<- read.csv("phenodata_age.txt",sep="\t",header=TRUE, dec = ",")

#Create a new filtered counts file
select <- which(targets$Type=="LS"| targets$Type=="SRME") # LS and SRME groups
Counts <- counts[,select]

# New phenofile with only the variables of interest
Targets <- targets[select,]

# Setup the design matrix, this is only to perform normalization on DESeq2 counts
condition <- as.character(Targets$Sex) # Irrelevant, can be any
batch <- as.character(Targets$NGS) # Batch effect
group_levels <- levels(as.factor(condition))
design <- data.frame(condition=as.factor(condition), batch=batch) # DE-design for the analysis
rownames(design) <- colnames(Counts)
dds <- DESeqDataSetFromMatrix(countData=Counts, colData=design, design = ~ batch + condition)

# DESeq2 DE-analysis of the condition of interest, batch effect taken into account
dds <- DESeq(dds)

# Export normalized counts from DESeq2
norm <- counts(dds, normalized=TRUE)

# Transpose the data for plotting
t.norm <- t(norm)

# Dimension reduction analysis----
# Run the t-SNE algorithm and store the results into an object called tsne_results
tsne_results <- Rtsne(t.norm, perplexity=35, check_duplicates = FALSE)

# Plotting----
# Generate a custom theme
theme_tero <-   
  theme_classic(base_size = 10)  +
  theme(legend.position = c(1,0.96), plot.margin = unit(c(1,4,1,1),"lines"),
        legend.box.spacing = unit(0.01, "cm")) +
  theme(axis.title = element_blank())  +
  theme(plot.tag = element_text(size = 30, face = "bold"))  +
  theme(axis.text = element_text(size = 15, color = "black")) 

# Generate a new data object for plotting
Data <- cbind(tsne_results$Y,Targets)

# Plot the images A-E
A <- ggplot(data=Data, aes(tsne_results$Y[,1],tsne_results$Y[,2])) + geom_point(size = 3) +  labs(tag = "A") + scale_color_discrete() + theme_tero
B <- ggplot(data=Data, aes(tsne_results$Y[,1],tsne_results$Y[,2], color= Targets$Health_2)) + geom_point(size = 3) + labs(tag = "B") +  scale_color_discrete(name = "Cancer status", labels =c("Healthy", "Path MMR cancer","Path SR cancer")) + theme_tero
C <- ggplot(data=Data, aes(tsne_results$Y[,1],tsne_results$Y[,2], color= Previous_ca_2)) + geom_point(size = 3) +  labs(tag = "C") + scale_color_discrete(name = "Cancer history",labels=c("Current cancer","Never","Previous cancer"))  + theme_tero
D <- ggplot(data=Data, aes(tsne_results$Y[,1],tsne_results$Y[,2], color= Targets$Variant_2)) + geom_point(size = 3) +  labs(tag = "D") + scale_color_discrete(name = "Mutation type")  + theme_tero
E <- ggplot(data=Data, aes(tsne_results$Y[,1],tsne_results$Y[,2], color= Targets$Age_2)) + geom_point(size = 3) +  labs(tag = "E") + scale_color_discrete(name = "Age")  + theme_tero

# Generate a new variable for age
AgeLevels <- factor(Targets$Age_3, levels = c("Over60", "fr50to60", "Under50"))  # Change order of legend items

# Plot the images F-H
F <- ggplot(data=Data, aes(tsne_results$Y[,1],tsne_results$Y[,2], color= AgeLevels)) + geom_point(size = 3) +  labs(tag = "F") + scale_color_discrete(name = "Age",labels = c( "Over 60", "Between 50 and 60", "Under 50"))  + theme_tero
G <- ggplot(data=Data, aes(tsne_results$Y[,1],tsne_results$Y[,2], color= Targets$BMI_2)) + geom_point(size = 3)  +  labs(tag = "G") + scale_color_discrete(name = "BMI") + theme_tero
H <- ggplot(data=Data, aes(tsne_results$Y[,1],tsne_results$Y[,2], color= Targets$Sex)) + geom_point(size = 3) +  labs(tag = "H") + scale_color_discrete(name = "Sex") + theme_tero

# Generate a new variable for batch effect
NGSLevels <- factor(Targets$NGS, levels = c("one", "two", "three"))  # Change order of legend items

# Plot the image I
I <- ggplot(data=Data, aes(tsne_results$Y[,1],tsne_results$Y[,2], color= NGSLevels)) + geom_point(size = 3) +   labs(tag = "I") +scale_color_discrete(name = "Batch") + theme_tero

# Generate the panel figure
plot1 <- plot_grid(A,B,C,D,E,F,G,H,I, ncol=3) # Panel using cowplot-package, three columns

# Save the figure
png(filename="tsne-panel3.png",height=700, width=1400)

# Print the figure
print(plot1)
dev.off()
