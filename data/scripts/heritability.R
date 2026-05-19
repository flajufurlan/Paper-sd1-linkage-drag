
rm(list=ls())

#Packages

library(sommer)
library(dplyr)

# Load the genotype data and create relationship matrices
geno <- readRDS("geno.rds")

# List of datasets
datasets <- list(
  "21MP6-8_LA1_phenotypes.csv", 
  "21MP6-8_LA2_phenotypes.csv", 
  "22_MP6-8_LA1_phenotypes.csv", 
  "22_MP6-8_LA2_phenotypes.csv", 
  "23_MP6-8_RRS_phenotypes.csv", 
  "23_MP6-8_RRSL_phenotypes.csv"
)

# Traits to analyze
traits <- c("Plant.height.LSU_01.0000132", "Yield.LSU_01.0000138")

# Initialize a dataframe to store heritability results
h2_results <- data.frame(
  dataset = character(),
  trait = character(),
  h2_gb = numeric(),
  h2_gn = numeric(),
  stringsAsFactors = FALSE
)

# Loop over datasets and traits
for (dataset in datasets) {
  
  # Load phenotype data
  pheno.v0 <- read.csv(dataset)
  
  # Subset the genotype data to match the phenotypes
  geno.v <- geno[rownames(geno) %in% pheno.v0$germplasmName, ]
  
  # Create G-matrix and D-matrix
  G_mat <- A.mat(X = as.matrix(geno.v), min.MAF = 0.05)
  D_mat <- D.mat(X = as.matrix(geno.v), min.MAF = 0.05)
  
  # Subset phenotype data to match the G-matrix rows
  pheno.v1 <- pheno.v0[pheno.v0$germplasmName %in% rownames(G_mat), ]
  pheno.v1$germplasmName2 <- pheno.v1$germplasmName  # Clone ID for dominance
  
  # Loop over each trait to estimate heritability
  for (trait in traits) {
    
    # Check if the trait exists in the dataset
    if (trait %in% colnames(pheno.v1)) {
      
      # Fit the model using the mmer function
      mm.gp <- mmer(fixed = as.formula(paste(trait, "~ 1")),
                    random = ~vsr(germplasmName, Gu = G_mat) + vsr(germplasmName2, Gu = D_mat),
                    rcov = ~ vsr(units),
                    data = pheno.v1)
      
      # Extract variance components
      sig_g_mks <- as.numeric(mm.gp$sigma_scaled[[1]])  # Genotype (additive)
      sig_d_mks <- as.numeric(mm.gp$sigma_scaled[[2]])  # Dominance
      sig_r_mks <- as.numeric(mm.gp$sigma_scaled[[3]])  # Residual
      
      # Calculate broad-sense heritability (h2_gb) and narrow-sense heritability (h2_gn)
      h2_gb <- (sig_g_mks + sig_d_mks) / (sig_g_mks + sig_d_mks + sig_r_mks)  # Broad-sense heritability
      h2_gn <- sig_g_mks / (sig_g_mks + sig_r_mks)  # Narrow-sense heritability
      
      # Append the results to the dataframe
      h2_results <- h2_results %>%
        add_row(dataset = dataset, trait = trait, h2_gb = h2_gb, h2_gn = h2_gn)
    }
  }
}

# Save the results to a CSV or RDS file
write.csv(h2_results, "heritability_results.csv", row.names = FALSE)
saveRDS(h2_results, "heritability_results.rds")
