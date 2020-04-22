# Load environment
library(glmnet)
library(readr)
library(bigmemory)
library(doMC)
library(bigstatsr)

# Load arguments
arg = commandArgs(trailingOnly=TRUE)

# Arguments are
# 1) working directory - where files are located
# 2) file of vertex wise values
# 3) output folder (in working directory)
# 4) variable (to predict) name
# 5) lasso or ridge
# 6) ncpus

wd=arg[1]
vertexFile=arg[2]
outFolder=arg[3]
pheno=arg[4]
lassoOrRidge=arg[5]
ncpus=arg[6]

# Working directory
setwd(wd)

# Big statsr approach - attache previously created FBM
st=Sys.time()
#dat=big_attach(rdsfile = paste0("BodFiles/", vertexFile, ".rds"))
dat=read.big.matrix(paste0("BodFiles/", vertexFile), sep = " ", header = T )
et=Sys.time()
openingTime=et-st

# Write time
print(paste("Opening/attaching the table took"))
openingTime

# Check dimensions
print(dim(dat))
nfeatures=dim(dat)[2]

# Subset into discovery adn replication
disco=read.table(paste0("BodFiles/ID_discovery_UKBanalysis.txt"))
repli=read.table(paste0("BodFiles/ID_replication_UKBanalysis.txt"))
# Open phenotype and order according to the ids list
phenoVar=read.big.matrix(paste0("UKB_phenotypes_reg_Bsln/", pheno, "_reg.txt"), sep = "\t", header = F )

# exclude participants with no phenotype in training set
disco=disco[which(disco$V1 %in% phenoVar[,1]),]

datD=dat[which(dat[,1] %in% disco$V1 ),]
datR=dat[which(dat[,1] %in% repli$V1 ),]

# Check dimensions of subsets
dim(datD)
dim(datR)

# Format for bigstatsr
datD_FBM=as_FBM(datD)
datR_FBM=as_FBM(datR)

print(pheno)
print("Opening phenotype file")

# Re order and match IDs with phenotypes
phenoM=phenoVar[match(datD[,1], phenoVar[,1], nomatch = NULL),]

cor(datD[,1], phenoM[,1])
cbind(datD[1:10,1], phenoM[1:10,1])

print("Cross validation using bigstatsr (NB not real CV)")
st1=Sys.time()
cv.bstatsr <-  big_spLinReg(X = datD_FBM, y.train = phenoM[,3], ind.col = 2:nfeatures, alphas = 1, nlambda = 100, K = 5)
et1=Sys.time()
cvRunTime=et1-st1
# Write time
print(paste("Running the CV using bigstatsr"))
cvRunTime

# Print lasso details
summary(cv.bstatsr)

# Store predictor
res=cbind(datR[,1], as.vector(predict(cv.bstatsr, datR_FBM )) )
write.table(res, paste0( outFolder, "/", "lasso", "Predictor_", pheno, "_bigstatsr.txt"), col.names = F, row.names = F, quote=F )

# No Writing weights for scores calculation in OSCA
# Function not available (not found in package)
