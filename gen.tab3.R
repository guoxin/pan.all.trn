beta.set <- c(0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16)
lambda.set <- exp(seq(-17, -9))
dyn.load("GK3WN.so")
stringToK3 <- function(string.set, K1, beta){
  len <- length(string.set)
  K3 <- .C("GK3WN", 
    peptides = as.character(string.set), 
    Length = as.integer(len), 
    K3 = as.double(matrix(0.0, nrow = len, ncol = len)), 
    K1 = K1 ^ beta,
    weights = as.double(rep(1, max(nchar(string.set)))), 
    DUP = TRUE, #because of the characters
    PACKAGE = "GK3WN"
   )$K3
   dim(K3) <- c(len, len);gc()
   colnames(K3) <- rownames(K3) <- names(string.set)
   K3
}

getAUC.guoxin <- function(predicted, realValue, threshold){
    AN <- (realValue <= threshold)
    AP <- (!AN)
    compare <- matrix(predicted[AP], ncol = sum(AN), nrow = sum(AP))
    compare <- sweep(x = compare, MARGIN = 2, 
        STATS = predicted[AN], FUN = "-", check.margin = TRUE)
    auc0 <- mean(compare > 0)
    return(max(auc0, 1 - auc0))
}
getRMSE <- function(predicted, realValue){
    return(sqrt(mean( (predicted - realValue)^2 )))
}
getAccuracy.guoxin <- function(predicted, realValue, threshold){
    return( mean((predicted > threshold) == 
                 (realValue > threshold)) )
}

gpgx.time.stamp <- function(...){
    write(print(paste( system("date", intern = TRUE),":",  ...)), 
        file = "performanceRec.txt", append = TRUE)
}

LooPanReport <- function(peptideStrings, alleleStrings, 
                         peptideList, alleleList,
                         K1, rValues, cvInds,
                         beta.set, lambda.set, threshold){
    thisReport <- list()
    foldNames <- c("fold1", "fold2", "fold3", "fold4", "fold5")
    rmse.iniMatrix <- matrix(0, ncol = length(lambda.set), nrow = length(beta.set))
    rownames(rmse.iniMatrix) <- beta.set; colnames(rmse.iniMatrix) <- log(lambda.set)
    thisReport$predict <- double(length(peptideList))
    for(fold.id in 1:5){
        trn.ind <- (cvInds != fold.id); tst.ind <- !trn.ind
        thisReport[[ foldNames[fold.id] ]]$predict <- 
            array(dim = c(sum(tst.ind), length(beta.set), length(lambda.set)))
        thisReport[[ foldNames[fold.id] ]]$rmse.trn <- rmse.iniMatrix
    }
    peptide.K3Hat <- stringToK3(string.set = peptideStrings, K1 = K1, beta = 0.11387)
    for(beta.id in 1:length(beta.set)){
        allele.K3Hat <- stringToK3(string.set = alleleStrings, K1 = K1, beta = beta.set[beta.id])
        K3Hat <- peptide.K3Hat[peptideList, peptideList] * allele.K3Hat[alleleList, alleleList]
        for(fold.id in 1:5){
            trn.ind <- (cvInds != fold.id); tst.ind <- !trn.ind
            gpgx.time.stamp("NOW COMPUTING:: beta:", beta.set[beta.id], "fold:", fold.id)
            nowSolve <- K3Hat[trn.ind, trn.ind]
            for(lambda.id in 1:length(lambda.set)){
                diag(nowSolve) <- 1 + sum(trn.ind) * lambda.set[lambda.id]
                G.inv <- solve(nowSolve); gc();  gpgx.time.stamp("G.inv obtained")
                coeff <- G.inv %*% rValues[trn.ind]; G.inv <- diag(G.inv); gc(); gpgx.time.stamp("coeff obtained")
                thisReport[[ foldNames[fold.id] ]]$predict[ ,beta.id, lambda.id] <- 
                    K3Hat[tst.ind, trn.ind] %*% coeff; gpgx.time.stamp("report 1 done")
                thisReport[[ foldNames[fold.id] ]]$rmse.trn[beta.id, lambda.id] <- 
                    sqrt(mean((coeff / (G.inv))^2)); gpgx.time.stamp("rmse.trn found")
            }
        }
        save(thisReport, file = "thisReport.Rdata")
    }
    for(fold.id in 1:5){
        bestRMSE <- min(thisReport[[ foldNames[fold.id] ]]$rmse.trn)
        best.ind <- which(thisReport[[ foldNames[fold.id] ]]$rmse.trn == bestRMSE, arr.ind = T)
        if(dim(best.ind)[1]>1) gpgx.time.stamp("multiMinPt found:",fold.id, 
                               dim(best.ind), best.ind[1,1], best.ind[1,2])
        best.ind <- best.ind[1, ]
        thisReport[[ foldNames[fold.id] ]]$best.beta.id <- best.ind[1]
        thisReport[[ foldNames[fold.id] ]]$best.lambda.id <- best.ind[2]
        thisReport$predict[cvInds == fold.id] <- 
            thisReport[[ foldNames[fold.id] ]]$predict[ ,best.ind[1], best.ind[2]]
        thisReport[[ foldNames[fold.id] ]]$predict <- NULL
    }
    thisReport$auc <- getAUC.guoxin(predicted = thisReport$predict,
        realValue = rValues, threshold = threshold)
    thisReport$rmse <- getRMSE(predicted = thisReport$predict, realValue = rValues)
    gc(); return(thisReport)
}
load("readpan20120114 md5.cc4a88a33be3296333121e9c4d524a8f")
reportsPan <- LooPanReport(
    peptideStrings = readpan20120114$peptide.set, 
    alleleStrings = readpan20120114$allele.set,
    peptideList = readpan20120114$peptide.list,
    alleleList = readpan20120114$allele.list,
    K1 = readpan20120114$K1,
    rValues = readpan20120114$rValue,
    cvInds = readpan20120114$cvInd,
    beta.set = beta.set,
    lambda.set = lambda.set,
    threshold = 1-(log(500)/log(15000)) # 0.35370904587312
)
save(reportsPan, file = "reportsPan")
