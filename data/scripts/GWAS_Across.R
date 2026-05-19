rm(list=ls())

#Packages

library(snpReady)

# loading the marker file
geno <- read.csv("Raw_Genotype_DarT_MP6_8.csv")
# Function to convert single-letter alleles to double-letter and remove slashes
convert_alleles <- function(genotype) {
  genotype <- gsub("A", "AA", genotype)
  genotype <- gsub("C", "CC", genotype)
  genotype <- gsub("G", "GG", genotype)
  genotype <- gsub("T", "TT", genotype)
  genotype <- gsub("/", "", genotype) # Remove slashes
  return(genotype)
}

# Apply the function to each column in the data frame
geno <- as.data.frame(lapply(geno, convert_alleles))

#hapmap
map <- read.csv("hapmap.csv")

matching_indices <- match(map[[1]], geno[[1]])
geno2 <- geno[matching_indices, ]

# recoding the marker names
M <- t(geno2[, 2:ncol(geno2)])
colnames(M) <- geno2[,1]
M[M == "NA"] <- NA
M[1:6, 1:6]
dim(M)
M <- as.data.frame(M)

# cleaning the data 
aux <- snpReady::raw.data(as.matrix(M), frame = "wide", hapmap = map, base = T, maf = 0.05, imput.type = "wright", outfile = "012")
M.clean <- aux$M.clean

# reshaping the vector into a matrix
M.clean <- matrix(aux$M.clean, byrow = F, nrow = nrow(M))

# get the newest dataset of markers
dim(M.clean)
M.clean[1:5, 1:5]
rownames(M.clean) <- rownames(M)

hapmap.clean <- aux$Hapmap
dim(hapmap.clean)
head(hapmap.clean)

colnames(M.clean) <- hapmap.clean$rs
M.clean[1:5, 1:5]

saveRDS(M.clean, "geno_clean.rds")
################################### GWAS #############################


# loading required library
library(factoextra)
library(ggplot2)
library(dplyr)

# loading phenotype 
pheno <- read.csv("BLUES_ALL_AR_LA.csv", header = TRUE, sep = ",") 
head(pheno)
tail(pheno)
str(pheno)


#filtered_pheno <- pheno %>%
  #filter(pop == "MP_FOUNDERS" |pop == "MP6") %>% filter(!(germplasmName == "PSDO_MP" | germplasmName == "CTHL_MP"))

#pheno <- filtered_pheno

pheno <- pheno[,c(1,3,4)]


# selecting just the common genotypes
match(rownames(M.clean), pheno$germplasmName)
M.clean <- M.clean[(rownames(M.clean) %in% pheno$germplasmName),]
pheno <- pheno[pheno$germplasmName %in% rownames(M.clean),]
pheno <- pheno[match(pheno$germplasmName, rownames(M.clean)),]
M.clean <- M.clean[match(rownames(M.clean), pheno$germplasmName),]
dim(M.clean)
head(M.clean[,1:6])
tail(M.clean[,1:6])
all(pheno$germplasmName %in% rownames(M.clean)) # should the answer be true??

all(rownames(M.clean) %in% pheno$germplasmName)
M.clean <- M.clean[match(pheno$germplasmName, rownames(M.clean)), ]
all(rownames(M.clean) == pheno$germplasmName)

# PCA analysis
#res.pca <- prcomp(Ga, scale = F)
#fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 85))

#pc <- prcomp(Ga,
#center = F,
#scale. = F)

# add a new column (gid) 
M2 <- data.frame(taxa = rownames(M.clean), M.clean)
head(M2[,1:6])
dim(M2)

#myKI <- data.frame(Taxa = rownames(Ga), diag(ncol(Ga)))
myGD <- M2
myGM <- hapmap.clean 
colnames(myGM)[1:3] <- c("Name", "Chromosome", "Position")

#myCV <- pc$x
#myCV<- data.frame(Taxa = rownames(myCV), myCV)


myY <- pheno
str(myY)
colnames(myY)[1] <- "Taxa"



#library(bigmemory)
#library(biganalytics)
#library(compiler)
#source("http://zzlab.net/GAPIT/gapit_functions.txt")
#Source("http://zzlab.net/FarmCPU/FarmCPU_functions.txt")

#####Running GWAS###
library(GAPIT)


myGAPIT<- GAPIT((Y=myY),
                GD=myGD,
                GM=myGM,
                #CV=myCV[, c(1,2)],
                model= "FarmCPU")


