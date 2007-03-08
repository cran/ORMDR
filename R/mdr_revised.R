##
## MDR
## 
## by SungGON Yi. <skon@kr.FreeBSD.org> and Eun-kyung Lee
##
##

require(combinat)

## x : a combination
errRate <- function(x, dataset, threshold) {
  t <- data.frame(table(dataset[, c(1, x)]))
  ratio <- t[t[, 1] == 1, "Freq"] / t[t[, 1] != 1, "Freq"]
  ratio[is.nan(ratio)] <- threshold
  weight <- ratio
  weight[ratio >= threshold] <- -0.5
  weight[ratio < threshold] <- 0.5
  weight <- rep(weight, rep(2, length(weight)))
  
  ## error rate of train set
  sum(t$Freq * abs(as.integer(t$resp) + weight - 1.5)) / nrow(dataset)
}

errRate.test <- function(x, dataset, test, threshold) {
  t <- data.frame(table(dataset[, c(1, x)]))
  ratio <- t[t[, 1] == 1, "Freq"] / t[t[, 1] != 1, "Freq"]
  ratio[is.nan(ratio)] <- threshold
  weight <- ratio
  weight[ratio >= threshold] <- -0.5
  weight[ratio < threshold] <- 0.5
  weight <- rep(weight, rep(2, length(weight)))

  ## error rate of test set
  t.test <- data.frame(table(test[, c(1, x)]))
  sum(t.test$Freq * abs(as.integer(t.test$resp) + weight - 1.5)) / nrow(test)
}

errRate.C <- function(comb, train, test, threshold) {
  z <- .C("err_rate", as.integer(comb), nrow(comb), ncol(comb),
          as.integer(unlist(train)), nrow(train), ncol(train),
          as.integer(unlist(test)), nrow(test),
          as.double(threshold), 
          err.train = double(ncol(comb)), 
          err.test = double(ncol(comb)), NAOK = FALSE,PACKAGE="ORMDR")
  
  cbind(train = z$err.train, test = z$err.test)
}  

##dataset <- read.csv(file1)

mdr.r <- function(dataset, colresp, cs, combi, cv.fold = 10) {
  resp <- dataset[, colresp]
  case <- which(resp == cs)
  ctl <- which(resp != cs)
  resp <- as.integer(resp == cs)
  snp <- dataset[, -colresp]

  cv.case <- matrix(0, nrow = length(case), ncol = 2)
  cv.case[, 1] <- 1:length(case) %% cv.fold + 1
  cv.case[, 2] <- sample(case)
  cv.ctl <- matrix(0, nrow = length(ctl), ncol = 2)
  cv.ctl[, 1] <- 1:length(ctl) %% cv.fold + 1
  cv.ctl[, 2] <- sample(ctl)
  cv <- rbind(cv.case, cv.ctl)
  d <- cbind(resp, snp)

  ## result
  z <- list()

  ## combinations
  comb <- list()
  for (k in 1:combi) {
    comb[[k]] <- combn(ncol(snp), k)
    comb[[k]] <- comb[[k]] + 1
    z[[k]] <- list()
    z[[k]]$min.comb <- matrix(0, nrow = k, ncol = cv.fold)
    z[[k]]$train.erate <- numeric(cv.fold)
    z[[k]]$test.erate <- numeric(cv.fold)
    z[[k]]$min.loc <- numeric(cv.fold)
  }

  for (i in 1:cv.fold) {
    train <- d[-cv[cv[, 1] == i, 2], ]
    test <- d[cv[cv[, 1] == i, 2], ]
    threshold <- length(which(train[, 1] == cs)) / 
      length(which(train[, 1] != cs)) 

    for (k in 1:combi) {
      ## train error ate for all combination of SNPs
      errs <- apply(comb[[k]], 2, errRate, train, threshold)
      z[[k]]$min.loc[i] <- which.min(errs)
      z[[k]]$min.comb[, i] <- comb[[k]][, z[[k]]$min.loc[i]]
      z[[k]]$train.erate[i] <- errs[z[[k]]$min.loc[i]]
      z[[k]]$test.erate[i] <- errRate.test(z[[k]]$min.comb[, i], train,
                                           test, threshold)
    }
  }

  ## maximum repeated selection
  for (k in 1:combi) {
    tmp <- table(z[[k]]$min.loc)
    max.cv.loc <- as.integer(names(tmp)[which.max(tmp)])
    z[[k]]$max.cv.comb <- comb[[k]][, max.cv.loc]
    z[[k]]$max.cv.num <- max(tmp)
  }
  
  z$data <- d
  z$comb <- comb

  z
}

## dataset object must be data.frame (or list).
## response data must be coded by 1 (case) or 0 (control).
## snp data must be coded by 0, 1, or 2.

#mdr.c <- function(dataset, colresp, cs, combi, cv.fold = 10) {
#  resp <- dataset[, colresp]
#  case <- which(resp == cs)
#  ctl <- which(resp != cs)
#  resp <- as.integer(resp == cs)
#  snp <- dataset[, -colresp]
#
#  cv.case <- matrix(0, nrow = length(case), ncol = 2)
#  cv.case[, 1] <- 1:length(case) %% cv.fold + 1
#  cv.case[, 2] <- sample(case)
#  cv.ctl <- matrix(0, nrow = length(ctl), ncol = 2)
#  cv.ctl[, 1] <- 1:length(ctl) %% cv.fold + 1
#  cv.ctl[, 2] <- sample(ctl)
#  cv <- rbind(cv.case, cv.ctl)
#  d <- cbind(resp, snp)
#
#  ## result
#  z <- list()
#
#  ## combinations
#  comb <- list()
#  for (k in 1:combi) {
#    comb[[k]] <- combn(ncol(snp), k)
#    comb[[k]] <- comb[[k]] + 1
#    z[[k]] <- list()
#    z[[k]]$min.comb <- matrix(0, nrow = k, ncol = cv.fold)
#    z[[k]]$train.erate <- numeric(cv.fold)
#    z[[k]]$test.erate <- numeric(cv.fold)
#   z[[k]]$min.loc <- numeric(cv.fold)
#  }
#
#  for (i in 1:cv.fold) {
#    train <- d[-cv[cv[, 1] == i, 2], ]
#    test <- d[cv[cv[, 1] == i, 2], ]
#    threshold <- length(which(train[, 1] == cs)) / 
#      length(which(train[, 1] != cs)) 
#
#    for (k in 1:combi) {
#      ## train and test error ate for all combination of SNPs
#      errs <- errRate.C(comb[[k]], train, test, threshold)
#      z[[k]]$min.loc[i] <- which.min(errs[, 1])
#      z[[k]]$min.comb[, i] <- comb[[k]][, z[[k]]$min.loc[i]]
#      z[[k]]$train.erate[i] <- errs[z[[k]]$min.loc[i], 1]
#      z[[k]]$test.erate[i] <- errs[z[[k]]$min.loc[i], 2]
#    }
#  }
#
#  ## maximum repeated selection
#  for (k in 1:combi) {
#    tmp <- table(z[[k]]$min.loc)
#    max.cv.loc <- as.integer(names(tmp)[which.max(tmp)])
#    z[[k]]$max.cv.comb <- comb[[k]][, max.cv.loc]
#    z[[k]]$max.cv.num <- max(tmp)
#  }
#  
#  z$data <- d
#  z$comb <- comb
#
#  z
#}

mdr.c <- function(dataset, colresp, cs, combi, cv.fold = 10) {
  resp <- dataset[, colresp]
  case <- which(resp == cs)
  ctl <- which(resp != cs)
  resp <- as.integer(resp == cs)
  snp <- dataset[, -colresp]

  cv.case <- matrix(0, nrow = length(case), ncol = 2)
  cv.case[, 1] <- 1:length(case) %% cv.fold + 1
  cv.case[, 2] <- sample(case)
  cv.ctl <- matrix(0, nrow = length(ctl), ncol = 2)
  cv.ctl[, 1] <- 1:length(ctl) %% cv.fold + 1
  cv.ctl[, 2] <- sample(ctl)
  cv <- rbind(cv.case, cv.ctl)
  d <- cbind(resp, snp)

  ## result
  z <- list()

  ## combinations
  comb <- list()
  k <-combi
    comb <- combn(ncol(snp), k)
    comb <- comb + 1
    z <- list()
    z$min.comb <- matrix(0, nrow = k, ncol = cv.fold)
    z$train.erate <- numeric(cv.fold)
    z$test.erate <- numeric(cv.fold)
    z$min.loc <- numeric(cv.fold)


  for (i in 1:cv.fold) {
    train <- d[-cv[cv[, 1] == i, 2], ]
    test <- d[cv[cv[, 1] == i, 2], ]
    threshold <- length(which(train[, 1] == cs)) / 
      length(which(train[, 1] != cs)) 

      ## train and test error ate for all combination of SNPs
      errs <- errRate.C(comb, train, test, threshold)
      z$min.loc[i] <- which.min(errs[, 1])
      z$min.comb[, i] <- comb[, z$min.loc[i]]
      z$train.erate[i] <- errs[z$min.loc[i], 1]
      z$test.erate[i] <- errs[z$min.loc[i], 2]
   }

  ## maximum repeated selection
    tmp <- table(z$min.loc)
    max.cv.loc <- as.integer(names(tmp)[which.max(tmp)])
    z$max.cv.comb <- comb[, max.cv.loc]
    z$max.cv.num <- max(tmp)
  
  z$data <- d
  z$comb <- comb

  name<-apply(t(z$min.comb),1,function(x){t<-NULL;
                                                                 for(i in 1:length(x)) 
                                                                       t<-paste(t,as.character(x[i]),sep="") ;return(t)})
  z$best.combi<-unlist(strsplit(names(sort(table(name),decreasing=T))[1],""))
  z
}

ormdr<-function(dataset,bestcombi,cs,colresp,CI.Asy=TRUE,CI.Boot=FALSE,B=5000)
{  
    case.id<-which(dataset[,colresp]==cs)
    control.id<-c(1:nrow(dataset))[-case.id]
  
     best.data.case<-list()
     best.data.cont<-list()
    for(i in 1:length(bestcombi))
    {     best.data.case[[i]]<-factor(dataset[case.id,bestcombi[i]],levels=c(0,1,2))
          best.data.cont[[i]]<-factor(dataset[control.id,bestcombi[i]],levels=c(0,1,2))
    }
   
    n.case<-length(case.id)
    n.cont<-length(control.id)
    t.case<-table(best.data.case)/n.case
    t.control<-table(best.data.cont)/n.cont
    Odds<-c(t.case/t.control)
   
    LU.Asy<-LU.Boot<-matrix(NA,ncol=2,nrow=3**length(bestcombi))
    if(CI.Asy)
    {         L<-c(exp(log(Odds)-1.96*sqrt((t.case)/(t.case*n.case)+(1-t.control)/(t.control*n.cont))))
            U<-c(exp(log(Odds)+1.96*sqrt((1-t.case)/(t.case*n.case)+(1-t.control)/(t.control*n.cont))))
            LU.Asy<-cbind(L,U)
    }
    if(CI.Boot)
    {     Odds.keep<-NULL
           for(i in 1:B)
          {     case.id<-sample(1:n.case,n.case,replace=T)
                 cont.id<-sample(1:n.cont,n.cont,replace=T)
                 best.data.case<-list()
                 best.data.cont<-list()
                 for(i in 1:length(bestcombi))
                {     best.data.case[[i]]<-factor(dataset[case.id,bestcombi[i]],levels=c(0,1,2))
                       best.data.cont[[i]]<-factor(dataset[control.id,bestcombi[i]],levels=c(0,1,2))
                }
                 P1.t<-c(table(best.data.case)/n.case)
                 P2.t<-c(table(best.data.cont)/n.cont)
                 Odds.t<-P1.t/P2.t
                 Odds.keep<-rbind(Odds.keep,c(Odds.t))
          }
          LU.Boot<-NULL
          for(id in 1:ncol(Odds.keep))
         {     L.t<-sort(Odds.keep[,id])[round(B*0.025)]
               U.t<-sort(Odds.keep[,id])[round(B*0.975)]
               LU.Boot<-rbind(LU.Boot,c(L.t,U.t))
          }

     }
     classID<-cbind(rep(0:2,3),rep(0:2,each=3))
     if(length(bestcombi)>2)
     {  for(i in 3:length(bestcombi))
             classID<-cbind(rbind(classID,classID,classID),rep(0:2,each=nrow(classID)))
     }
     cell.freq<-cbind(c(table(best.data.case)),c(table(best.data.cont)))
     Hi.Low<-ifelse(Odds>=1,"High","Low")
     ORMDR.table<-cbind(classID,cell.freq,Hi.Low,round(Odds,3),rank(Odds),round(LU.Asy,3),round(LU.Boot,3))
     colnames(ORMDR.table)<-c(colnames(RAW.data)[SNPlist.id][bestcombi],                                        "case.freq","cont.freq","Hi.Low","Odds.ratio","Rank","Asy.L","Asy.U","Boot.L","Boot.U")
     return(ORMDR.table)
}

