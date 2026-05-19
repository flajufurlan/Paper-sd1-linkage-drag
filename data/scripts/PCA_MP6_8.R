
#Packages


library(snpReady)
library(dplyr)
library(ggplot2)
library(factoextra)

## 2. Load and clean genotype data

# Load genotype
geno <- read.csv("geno.csv", stringsAsFactors = FALSE)

# Convert single-letter alleles to AA/CC/GG/TT
convert_alleles <- function(x) {
  x <- gsub("A", "AA", x)
  x <- gsub("C", "CC", x)
  x <- gsub("G", "GG", x)
  x <- gsub("T", "TT", x)
  x <- gsub("/", "", x)
  return(x)
}

geno <- as.data.frame(lapply(geno, convert_alleles), stringsAsFactors = FALSE)

# Load hapmap
map <- read.csv("hapmap.csv", stringsAsFactors = FALSE)

# Match marker order
matching_indices <- match(map[[1]], geno[[1]])
geno2 <- geno[matching_indices, ]

# Transpose to taxa x SNPs
M <- t(geno2[, -1])
colnames(M) <- geno2[, 1]
M[M == "NA"] <- NA
M <- as.data.frame(M)

## 3. SNP QC + Imputation (snpReady)

aux <- snpReady::raw.data(
  as.matrix(M),
  frame = "wide",
  hapmap = map,
  base = TRUE,
  maf = 0.05,
  imput.type = "wright",
  outfile = "012"
)

# Clean marker matrix
M.clean <- matrix(aux$M.clean, nrow = nrow(M), byrow = FALSE)
rownames(M.clean) <- rownames(M)

# Clean hapmap
hapmap.clean <- aux$Hapmap
colnames(M.clean) <- hapmap.clean$rs

# Save clean genotype
saveRDS(M.clean, "geno_clean.rds")

## 4. Load and align phenotype

pheno <- read.csv("BLUES_ALL_AR_LA.csv", stringsAsFactors = FALSE)

# Keep only needed columns
pheno <- pheno[, c(1, 3, 4)]
colnames(pheno)[1] <- "Taxa"

# Match genotype and phenotype
common <- intersect(rownames(M.clean), pheno$Taxa)

M.clean <- M.clean[common, ]
pheno <- pheno[match(common, pheno$Taxa), ]

stopifnot(all(rownames(M.clean) == pheno$Taxa))


## 5. PCA 

# Remove monomorphic SNPs (THIS FIXES NaN ISSUE)
snp_var <- apply(M.clean, 2, var, na.rm = TRUE)
M.var <- M.clean[, snp_var > 0]

# Check dimensions
dim(M.clean)   # original
dim(M.var)     # used for PCA

# PCA
pca <- prcomp(
  M.var,
  center = TRUE,
  scale. = FALSE
)

# Percent variance explained
pve <- (pca$sdev^2) / sum(pca$sdev^2) * 100
pve[1:5]

# PCA scores
pca_df <- data.frame(
  Taxa = rownames(M.var),
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  PC3 = pca$x[,3]
)


# PC1 vs PC2
p_pca_12 <- ggplot(pca_df, aes(PC1, PC2)) +
  geom_point(size = 2, alpha = 0.7) +
  theme_bw() +
  labs(
    x = paste0("PC1 (", round(pve[1], 2), "%)"),
    y = paste0("PC2 (", round(pve[2], 2), "%)")
  )


ggsave(
  filename = "PCA_MP6_8_PC1_PC2.tiff",
  plot = p_pca_12,
  device = "tiff",
  dpi = 300,
  width = 3.15,
  height = 3.15,
  units = "in",
  compression = "lzw"
)

# Scree plot
fviz_eig(
  pca,
  addlabels = TRUE,
  ylim = c(0, max(pve) + 5)
)

