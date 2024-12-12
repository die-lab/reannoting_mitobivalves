#setwd to the right dir

library(ggplot2)
library(car)
library(gridExtra)
library(tidyr)

df <- read.csv("class_tree_datasets.csv", sep=";")
table(as.factor(dataset$Class))

#fix the Order field, taking from taxonomy
extract_fourth_last_word <- function(taxonomy) {
  # Split the taxonomy string by spaces
  taxonomy_split <- unlist(strsplit(taxonomy, ";"))
  # Check if the taxonomy has at least 4 words
  if (length(taxonomy_split) >= 4) {
    # Return the fourth last word
    return(taxonomy_split[length(taxonomy_split) - 3])
  } else {
    # If there are less than 4 words, return NA or the full taxonomy as fallback
    return(NA)
  }
}
df$Order <- sapply(df$Taxonomy, extract_fourth_last_word)

#plotting mitochondrial length per class
ggplot(df, aes(x = Class, y = Mitochondrial_length)) +
  geom_boxplot() +           # Use boxplot to show the distribution
  theme_minimal() +          # Apply a clean theme
  labs(title = "Mitochondrial Length by Class",
       x = "Class",
       y = "Mitochondrial Length")

##ARE THEY SIGNIFICANTLY DIFFERENT?
anova_result <- aov(Mitochondrial_length ~ Class, data = df)
summary(anova_result)

#pair comparisons
tukey_result <- TukeyHSD(anova_result)
summary(tukey_result)

#gene_lenght plotting
gene_columns <- c("Mitochondrial_length", "ATP6_Gene_length", "ATP8_Gene_length", 
                  "COX1_Gene_length", "COX2_Gene_length", "COX3_Gene_length", 
                  "CYTB_gene_length", "ND1_gene_length", "ND2_gene_length", 
                  "ND3_gene_length", "ND4_gene_length", "ND4L_gene_length", 
                  "ND5_gene_length", "ND6_gene_length")
# Create a list to store the plots
plots <- list()
# Loop through each gene column and create a boxplot for each
for (gene in gene_columns) {
  p <- ggplot(df, aes(x = Class, y = df[[gene]])) +
    geom_boxplot() +
    theme_minimal() +
    labs(title = paste("Variance of", gene),
         x = "Class", y = gene)
  
  # Store each plot in the list
  plots[[gene]] <- p
}
# Combine all the plots into one image using grid.arrange()
do.call(grid.arrange, c(plots, ncol = 3))  # Adjust ncol for the number of columns in the grid


## find for each gene which classes has the highest weighted standard deviation
library(dplyr)

# Step 1: Calculate the weighted standard deviation for each gene and class (same as previous steps)
weighted_sd_per_class_gene <- data.frame(Class = character(), Gene = character(), Weighted_SD = numeric())

for (class in unique(df$Class)) {
  for (gene in gene_length_columns) {
    # Filter the data for the current class and gene length column
    class_data <- df[df$Class == class, ]
        # Extract the gene lengths for the current gene in the current class
    gene_lengths <- class_data[[gene]]
        # Get the number of observations in the current class (as weight for each observation)
    n_class <- length(gene_lengths)
        # Calculate the weighted mean for the gene length
    weighted_mean <- sum(gene_lengths, na.rm = TRUE) / n_class
        # Calculate the weighted variance for the gene lengths
    weighted_variance <- sum((gene_lengths - weighted_mean)^2, na.rm = TRUE) / n_class
        # Calculate the weighted standard deviation
    weighted_sd <- sqrt(weighted_variance)
        # Store the result in the data frame
    weighted_sd_per_class_gene <- rbind(weighted_sd_per_class_gene, 
                                        data.frame(Class = class, Gene = gene, Weighted_SD = weighted_sd))
  }
}
# Step 2: Calculate the difference in weighted standard deviation for each gene
sd_diff_per_gene <- weighted_sd_per_class_gene %>%
  group_by(Gene) %>%
  summarise(
    Max_SD = max(Weighted_SD, na.rm = TRUE),
    Min_SD = min(Weighted_SD, na.rm = TRUE),
    SD_Diff = Max_SD - Min_SD
  ) %>%
  ungroup()
# Step 3: Find the gene with the maximum difference in weighted standard deviation
max_diff_gene <- sd_diff_per_gene %>%
  filter(SD_Diff == max(SD_Diff, na.rm = TRUE))

# Output the full list of genes with max and min SD and the difference
sd_diff_per_gene
# Output the gene with the maximum difference in SD
max_diff_gene

##plot it
ggplot(weighted_sd_per_class_gene, aes(x = Class, y = Weighted_SD, fill = Class)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Gene, scales = "free_y") +  # One plot for each gene
  theme_minimal() +
  labs(title = "Weighted Standard Deviation of Gene Lengths per Gene and Class",
       x = "Class", y = "Weighted Standard Deviation") +
  scale_fill_brewer(palette = "Set3") +  # Use a color palette for distinction
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
        legend.title = element_text(size = 10),  # Adjust legend title size
        legend.text = element_text(size = 8))  # Adjust legend text size

####the same as before, but this time it gives me the Order along the most contributing to the standard deviation differences
# Step 1: Identify the gene with the maximum difference in weighted standard deviation (already found)
max_diff_gene <- sd_diff_per_gene %>%
  filter(SD_Diff == max(SD_Diff, na.rm = TRUE))
# Extract the gene with the maximum difference
max_diff_gene_name <- max_diff_gene$Gene
# Step 2: Calculate the residuals and squared residuals for this gene across its classes
residuals_per_gene_class <- data.frame(Class = character(), Gene = character(), 
                                       Observation = integer(), Residual = numeric(), 
                                       Squared_Residual = numeric(), Mitochondrial_genomes_codes = character(), 
                                       Order = character())
for (class in unique(df$Class)) {
  # Filter the data for the current gene and class
  class_data <- df[df$Class == class, ]
    # Extract the gene lengths for the current gene
  gene_lengths <- class_data[[max_diff_gene_name]]
    # Get the mitochondrial genome codes for the current class and gene
  mitochondrial_genomes <- class_data$Mitochondrial_genomes_codes
    # Get the Order values for the current class and gene
  order_values <- class_data$Order
    # Get the number of observations in the current class
  n_class <- length(gene_lengths)
    # Calculate the weighted mean for the gene length
  weighted_mean <- sum(gene_lengths, na.rm = TRUE) / n_class
    # Step 3: Calculate the residuals (differences between each observation and the weighted mean)
  residuals <- gene_lengths - weighted_mean
    # Step 4: Calculate the squared residuals
  squared_residuals <- residuals^2
    # Store the residuals, squared residuals, mitochondrial genome codes, and Order
  for (i in 1:n_class) {
    residuals_per_gene_class <- rbind(residuals_per_gene_class, 
                                      data.frame(Class = class, Gene = max_diff_gene_name, 
                                                 Observation = i, Residual = residuals[i], 
                                                 Squared_Residual = squared_residuals[i],
                                                 Mitochondrial_genomes_codes = mitochondrial_genomes[i],
                                                 Order = order_values[i]))
  }
}
# Step 5: Rank the observations based on squared residuals
ranked_contributors <- residuals_per_gene_class %>%
  arrange(desc(Squared_Residual)) %>%
  mutate(Rank = row_number())
# Step 6: Select only the columns that show the rank, mitochondrial genome codes, and Order
ranked_contributors_output <- ranked_contributors %>%
  select(Rank, Mitochondrial_genomes_codes, Order, Squared_Residual, Gene)
# Step 7: View the ranked contributors with the order, mitochondrial genome codes, and rank
head(ranked_contributors_output,20)
