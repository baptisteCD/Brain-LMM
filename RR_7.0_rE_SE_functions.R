
# Find the file
findVCOVmatFromGCTAOutput=function(var1, var2, path){
  flist=list.files(path)
   VCValues=matrix(NA, nrow = 6,ncol = 3)
    VCOVAR=matrix(NA, nrow = 6,ncol = 6)
  for (ffile in flist[grep(x = flist, pattern = ".log")]){

    rd=readLines(paste0(path, ffile))
    isV1=length(rd[grep(x = rd, pattern = var1)])
    isV2=length(rd[grep(x = rd, pattern = var2)])

    if (isV1 * isV2 > 0 & length(which(rd =="Sampling variance/covariance of the estimates of variance components:" ) ) >0){


    varNamesRowStart= which(rd =="Sampling variance/covariance of the estimates of variance components:" )+1
    varNamesRowEnd= which(rd =="Sampling variance/covariance of the estimates of variance components:" )+6
    matValues=as.numeric(c(unlist(strsplit(rd[varNamesRowStart:varNamesRowEnd], split="\t"))))
    VCOVAR=matrix(matValues,6,6)

    varNamesRowStart= which(rd =="Source\tVariance\tSE" )
    matValues=as.numeric(c(unlist(strsplit(rd[ (varNamesRowStart +1):(varNamesRowStart+6)], split="\t"))))
    VCValues=matrix(matValues,6,3, byrow = T)

    print(ffile)
    }
  }
  return( cbind(VCOVAR, VCValues) )
}

# Calculation
seRgFullApprox=function( VCX, VCY, CovXY , Var_CovXY, Var_VCX, Var_VCY, Cov_CovXY_VCX, Cov_CovXY_VCY, Cov_VCX_VCY , N){
  varRG= Var_CovXY / (VCX * VCY) - Cov_CovXY_VCX * (CovXY / (VCX**2)-VCY) - Cov_CovXY_VCY * (CovXY / (VCY**2)-VCX) + Var_VCX * ( CovXY **2 / (4* VCX **3 * VCY) ) + Var_VCY * ( CovXY **2 / (4* VCY **3 * VCX) ) + Cov_VCX_VCY * (CovXY **2 / (2* VCX **2 * VCY **2))

  seRG= sqrt(varRG)
  return(c( seRG ))
}
