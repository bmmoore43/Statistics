#principal components analysis

setwd()

#install packages
library(devtools)
library(curl)
install_github("vqv/ggbiplot")
library(ggbiplot)

#try PCA with wine data and make plot
data(wine)
wine.pca <- prcomp(wine, scale. = TRUE)
wine.class
print(ggbiplot(wine.pca, obs.scale = 1, var.scale = 1, groups = wine.class, ellipse = TRUE, circle = TRUE))

# with RNA seq sample data
# Create a "data" directory
dir.create("data")

# Download the data provided by your collaborator
# using a for loop to automate this step
for(i in c("counts_raw.csv", "counts_transformed.csv", "sample_info.csv", "test_result.csv")){
  download.file(
    url = paste0("https://github.com/tavareshugo/data-carpentry-rnaseq/blob/master/data/", i, "?raw=true"),
    destfile = paste0("data/", i)
  )
}
raw_cts <- read.csv("./data/counts_raw.csv")
trans_cts <- read.csv("./data/counts_transformed.csv")
sample_info <- read.csv("./data/sample_info.csv")
test_result <- read.csv("./data/test_result.csv")

library(tidyverse)

#Check out raw data
#convert to "long form"
raw_cts_long <- raw_cts %>% 
  pivot_longer(wt_0_r1:mut_180_r3, names_to = "sample", values_to = "cts")
# Join with sample information table
raw_cts_long <- full_join(raw_cts_long, sample_info, by = ("sample"))
#plot boxplot
raw_cts_long %>%
  # make sure minute is specified as a factor
  ggplot(aes(factor(minute), log10(cts + 1), fill = strain)) + 
  geom_boxplot() + 
  facet_grid(cols = vars(replicate))

# Create a matrix from our table of counts
pca_matrix <- trans_cts %>% 
  # make the "gene" column become the rownames of the table
  column_to_rownames("gene") %>% 
  # coerce to a matrix
  as.matrix() %>% 
  # transpose the matrix so that rows = samples and columns = variables
  t()

# Perform the PCA
sample_pca <- prcomp(pca_matrix)

# elements of the PCA:
# PC scores: new coordinates of the data on the PC axis
# eigenvalues: variance explained by each PC
# variable loadings: “weight” that each original gene has on each PC axis
# check them out:
pc_scores <- sample_pca$x              # PC scores (a matrix)
pc_eigenvalues <- sample_pca$sdev^2    # eigenvalues (a vector) - notice we square the values
pc_loadings <- sample_pca$rotation     # variable loadings (a matrix)
# number of PCAs:
ncol(pc_scores)

# Annotate a PC plot
pca_plot <- sample_pca$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample") %>% 
  # join with "sample_info" table 
  full_join(sample_info, by = "sample") %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(minute), shape = strain)) +
  geom_point()

# print the result (in this case a ggplot)
pca_plot

# now try with real data
counts <- read.csv("Read_counts.txt_TPM.txt", sep="\t")
sample_info <- read.csv("sample_info.txt", sep="\t")

#plot replicates
#convert to "long form"
cts_long <- counts %>% 
  pivot_longer(flower.1_sam_sorted.sam:young.shoot.3_sam_sorted.sam, names_to = "sample", values_to = "cts")
# Join with sample information table
cts_long <- full_join(cts_long, sample_info, by = ("sample"))
#plot boxplot
cts_long %>%
  # make sure replicate is specified as a factor
  ggplot(aes(factor(replicate), cts, fill = strain)) + 
  geom_boxplot() + 
  facet_grid(cols = vars(strain))

#plot boxplot without outliers
cts_long %>%
  # make sure replicate is specified as a factor
  ggplot(aes(factor(replicate), cts, fill = strain)) + 
  geom_boxplot(outlier.shape = NA) + scale_y_continuous(limits = quantile(cts_long$cts, c(0.1, 0.9))) +
  facet_grid(cols = vars(strain))

# Create a matrix from our table of counts
pca_matrix <- counts %>% 
  # make the "gene" column become the rownames of the table
  column_to_rownames("gene") %>% 
  # coerce to a matrix
  as.matrix() %>% 
  # transpose the matrix so that rows = samples and columns = variables
  t()

# Perform the PCA
sample_pca <- prcomp(pca_matrix)

# elements of the PCA:
# PC scores: new coordinates of the data on the PC axis
# eigenvalues: variance explained by each PC
# variable loadings: “weight” that each original gene has on each PC axis
# check them out:
pc_scores <- sample_pca$x              # PC scores (a matrix)
pc_eigenvalues <- sample_pca$sdev^2    # eigenvalues (a vector) - notice we square the values
pc_loadings <- sample_pca$rotation     # variable loadings (a matrix)
# number of PCAs:
ncol(pc_scores)

# Annotate a PC plot
pca_plot <- sample_pca$x %>% # extract the loadings from prcomp
  # convert to a tibble retaining the sample names as a new column
  as_tibble(rownames = "sample") %>% 
  # join with "sample_info" table 
  full_join(sample_info, by = "sample") %>% 
  # create the plot
  ggplot(aes(x = PC1, y = PC2, colour = factor(strain), shape = replicate)) +
  geom_point()

# print the result (in this case a ggplot)
pca_plot
