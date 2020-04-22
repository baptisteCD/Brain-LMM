# Load environment
library(glmnet)
library(readr)
library(bigmemory)
library(doMC)

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
# recruit cpus
registerDoMC(cores = ncpus)

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

# Cross validation for
st1=Sys.time()
cv.lasso <- cv.glmnet(y=phenoM[,3], x=datD[,-c(1,2)], alpha = ifelse(lassoOrRidge=="lasso",1,0), nfolds=5, parallel = TRUE)
et1=Sys.time()
cvRunTime=et1-st1
# Write time
print(paste("Running the CV took"))
cvRunTime

# Print lasso details
print(cv.lasso)

# Fit models with best lambda values
model.lasso <- glmnet(y=phenoM[,3], x=datD[,-c(1,2)], alpha = ifelse(lassoOrRidge=="lasso",1,0), lambda = cv.lasso$lambda.min)

# Store predictor
res=cbind(datR[,1], as.vector(predict(model.lasso, newx = datR[,-c(1,2)] )) )
write.table(res, paste0( outFolder, "/", lassoOrRidge, "Predictor_", pheno, ".txt"), col.names = F, row.names = F, quote=F )

# Write weights for scores calculation in OSCA
resW=as.data.frame(  as.matrix(coef(model.lasso)))
resW$vertName=rownames(as.matrix(coef(model.lasso)))
if(length(which(resW[,1]==0))>0){
resW=resW[-which(resW[,1]==0),] }
write.table(resW[-1,c(2,1)], paste0(outFolder, "/Weights_",lassoOrRidge,"_" , pheno, ".txt"), col.names = F, quote=F, row.names = F)

