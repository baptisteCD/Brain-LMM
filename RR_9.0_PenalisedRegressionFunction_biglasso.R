# Load environment
library(glmnet)
library(readr)
library(bigmemory)
library(doMC)
library(biglasso)

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

# Open using C++ big memory wrap up
st=Sys.time()
dat=read.big.matrix(paste0("BodFiles/", vertexFile), sep = " ", header = T )
et=Sys.time()
openingTime=et-st

# Write time
print(paste("Opening the table took"))
openingTime

# Check dimensions
print(dim(dat))

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

# Re order and match IDs with phenotypes
phenoM=phenoVar[match(datD[,1], phenoVar[,1]),]

cor(datD[,1], phenoM[,1])
cbind(datD[1:10,1], phenoM[1:10,1])

# Transform data for biglasso function
X2=as.big.matrix(datD[,-c(1,2)],  backingfile = "")
X2R=as.big.matrix(datR[,-c(1,2)],  backingfile = "")

# Check dimensions of subsets
dim(X2)
dim(X2R)

# Cross validation using biglasso
st1=Sys.time()
cv.blasso <-  cv.biglasso(y=phenoM[,3], X=X2, family = "gaussian", seed = 1234, ncores = ncpus, nfolds = 5, penalty =lassoOrRidge )
et1=Sys.time()
cvRunTime=et1-st1
# Write time
print(paste("Running the CV using biglasso took"))
cvRunTime

# Print lasso details
print(summary(cv.blasso))

# Store predictor
res=cbind(datD[,1], as.vector(predict(cv.blasso, X = X2R )) )
write.table(res, paste0( outFolder, "/", lassoOrRidge, "Predictor_", pheno, "_biglasso_noBin.txt"), col.names = F, row.names = F, quote=F )

# Write weights for scores calculation in OSCA
resW=as.data.frame(  as.matrix(coef(cv.blasso)))
resW$vertName=rownames(as.matrix(coef(cv.blasso)))
if(length(which(resW[,1]==0))>0){
resW=resW[-which(resW[,1]==0),] }
write.table(resW[-1,c(2,1)], paste0(outFolder, "/Weights_",lassoOrRidge,"_" , pheno, "_biglasso_noBin.txt"), col.names = F, quote=F, row.names = F)

