#############################################################################
#### This code can be used to find the designs used in the manuscript.
#### Depending on computing power, this may take days. As an alternative,
#### one may simply load the file "data.RData".
#### Another alternative is to decrease the value of max.combns from 
#### 1e6 to 1e5 or even 1e4, though this will result in poorer designs.
#############################################################################
 

source("1-reproduce-results-functions-find-designs.R")
library(doParallel)
library(data.table)
library(rlist)
library(ph2rand) # For Jung results

########### Section 3.1: Comparing design approaches using multiple criteria ############
# (also used for Section 3.2: Misspecification of response rates)

#### Find designs for block size 2 and block size 8:

pc <- c(0.1, 0.2, 0.3, 0.4, 0.5)
pt <- pc + 0.2

block2.n <- seq(from=8, to=60, by=1)
block8.n <- seq(from=8, to=60, by=4)

# NOTE: 30 cores were used for original timings.
cores <- detectCores()-2

max.combns <- 1e6
# NOTE: The simulations take a long time using 1e6 max combinations. Recommend max.combns <- 1e4 or 1e5 for checking, or load attached .RData file.

j <- 1
registerDoParallel(cores)
timerow1.to60.1e6.noskip <- system.time({
  row1.to60.1e6.noskip <- foreach(i=block2.n, .combine = rbind) %dopar% {
    find2armDesigns(nmin=i, nmax=i, block.size=2, pc=pc[j], pt=pt[j], alpha=0.15, power=0.8, maxtheta0=pt[j], mintheta1=0.7, bounds="ahern", max.combns=max.combns)
  }
})[3] # 45mins for 8 to 50, 1e5 combns; 12 hours for 1e6
stopImplicitCluster()
registerDoParallel(cores)
timerow1.to60.1e6.noskip.block8 <- system.time({
  row1.to60.1e6.noskip.block8 <- foreach(i=block8.n, .combine = rbind) %dopar% {
    find2armDesigns(nmin=i, nmax=i, block.size=8, pc=pc[j], pt=pt[j], alpha=0.15, power=0.8, maxtheta0=pt[j], mintheta1=0.7, bounds="ahern", max.combns=max.combns)
  }
})[3] 
stopImplicitCluster()


j <- 2
registerDoParallel(cores)
timerow2.to60.1e6.noskip <- system.time({
  row2.to60.1e6.noskip <- foreach(i=block2.n, .combine = rbind) %dopar% {
    find2armDesigns(nmin=i, nmax=i, block.size=2, pc=pc[j], pt=pt[j], alpha=0.15, power=0.8, maxtheta0=pt[j], mintheta1=0.7, bounds="ahern", max.combns=max.combns)
  }
})[3] 
stopImplicitCluster()

registerDoParallel(cores)
timerow2.to60.1e6.noskip.block8 <- system.time({
  row2.to60.1e6.noskip.block8 <- foreach(i=block8.n, .combine = rbind) %dopar% {
    find2armDesigns(nmin=i, nmax=i, block.size=8, pc=pc[j], pt=pt[j], alpha=0.15, power=0.8, maxtheta0=pt[j], mintheta1=0.7, bounds="ahern", max.combns=max.combns)
  }
})[3] 
stopImplicitCluster()


j <- 3
registerDoParallel(cores)
timerow3.to60.1e6.noskip <- system.time({
  row3.to60.1e6.noskip <- foreach(i=block2.n, .combine = rbind) %dopar% {
    find2armDesigns(nmin=i, nmax=i, block.size=2, pc=pc[j], pt=pt[j], alpha=0.15, power=0.8, maxtheta0=pt[j], mintheta1=0.7, bounds="ahern", max.combns=max.combns)
  }
})[3] 
stopImplicitCluster()
registerDoParallel(cores)
timerow3.to60.1e6.noskip.block8 <- system.time({
  row3.to60.1e6.noskip.block8 <- foreach(i=block8.n, .combine = rbind) %dopar% {
    find2armDesigns(nmin=i, nmax=i, block.size=8, pc=pc[j], pt=pt[j], alpha=0.15, power=0.8, maxtheta0=pt[j], mintheta1=0.7, bounds="ahern", max.combns=max.combns)
  }
})[3] 
stopImplicitCluster()


j <- 4
registerDoParallel(cores)
timerow4.to60.1e6.noskip <- system.time({
  row4.to60.1e6.noskip <- foreach(i=block2.n, .combine = rbind) %dopar% {
    find2armDesigns(nmin=i, nmax=i, block.size=2, pc=pc[j], pt=pt[j], alpha=0.15, power=0.8, maxtheta0=pt[j], mintheta1=0.7, bounds="ahern", max.combns=max.combns)
  }
})[3] 
stopImplicitCluster()
registerDoParallel(cores)
timerow4.to60.1e6.noskip.block8 <- system.time({
  row4.to60.1e6.noskip.block8 <- foreach(i=block8.n, .combine = rbind) %dopar% {
    find2armDesigns(nmin=i, nmax=i, block.size=8, pc=pc[j], pt=pt[j], alpha=0.15, power=0.8, maxtheta0=pt[j], mintheta1=0.7, bounds="ahern", max.combns=max.combns)
  }
})[3] 
stopImplicitCluster()


j <- 5
registerDoParallel(cores)
timerow5.to60.1e6.noskip <- system.time({
  row5.to60.1e6.noskip <- foreach(i=block2.n, .combine = rbind) %dopar% {
    find2armDesigns(nmin=i, nmax=i, block.size=2, pc=pc[j], pt=pt[j], alpha=0.15, power=0.8, maxtheta0=pt[j], mintheta1=0.7, bounds="ahern", max.combns=max.combns)
  }
})[3] 
stopImplicitCluster()
registerDoParallel(cores)
timerow5.to60.1e6.noskip.block8 <- system.time({
  row5.to60.1e6.noskip.block8 <- foreach(i=block8.n, .combine = rbind) %dopar% {
    find2armDesigns(nmin=i, nmax=i, block.size=8, pc=pc[j], pt=pt[j], alpha=0.15, power=0.8, maxtheta0=pt[j], mintheta1=0.7, bounds="ahern", max.combns=max.combns)
  }
})[3] 
stopImplicitCluster()



block2.1e6 <- list(row1.to60.1e6.noskip, row2.to60.1e6.noskip, row3.to60.1e6.noskip, row4.to60.1e6.noskip, row5.to60.1e6.noskip)
block8.1e6 <- list(row1.to60.1e6.noskip.block8, row2.to60.1e6.noskip.block8, row3.to60.1e6.noskip.block8, row5.to60.1e6.noskip.block8, row4.to60.1e6.noskip.block8)

for(i in 1:5){
  block2.1e6[[i]] <- rmDominatedDesigns(df=block2.1e6[[i]], n="eff.n")
  block8.1e6[[i]] <- rmDominatedDesigns(df=block8.1e6[[i]], n="eff.n")  
}






################### Find Carsten designs ##################
# Initial (parallel) run:

df1 <- findN1N2R1R2twoarm(nmin=5, nmax=100)
nrow(df1)
df1 <- df1[df1$r1>=0 & df1$r2>=0 & df1$r1<=df1$n1 & df1$r1<=df1$r2,]
nrow(df1)
# Split into lists, by r2:
df2 <- split(df1, list(df1$n1, df1$n, df1$r1))
# remove empty lists
# library(rlist)
nr.combn.list <- list.clean(df2, function(x) nrow(x) == 0L, recursive = F)
rm(df1, df2)

pc <- c(0.1, 0.2, 0.3, 0.4, 0.5)
pt <- pc + 0.2

init.runs <- 100
init.alpha <- 0.25
init.power <- 0.7

j <- 1
registerDoParallel(cores)
row1results.carsten <- foreach(i=1:length(nr.combn.list)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=nr.combn.list[[i]], pc=pc[j], pt=pt[j], runs=init.runs, alpha=init.alpha, power=init.power, seed=2686, method = "carsten")
  } # 100 runs +/- 0.1 alpha and power -- 

j <- 2
row2results.carsten <- foreach(i=1:length(nr.combn.list)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=nr.combn.list[[i]], pc=pc[j], pt=pt[j], runs=init.runs, alpha=init.alpha, power=init.power, seed=2686, method = "carsten")
  } #  100 runs +/- 0.1 alpha and power -- 

j <- 3
row3results.carsten <- foreach(i=1:length(nr.combn.list)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=nr.combn.list[[i]], pc=pc[j], pt=pt[j], runs=init.runs, alpha=init.alpha, power=init.power, seed=2686, method = "carsten")
  } #  100 runs +/- 0.1 alpha and power --

j <- 4
  row4results.carsten <- foreach(i=1:length(nr.combn.list)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=nr.combn.list[[i]], pc=pc[j], pt=pt[j], runs=init.runs, alpha=init.alpha, power=init.power, seed=2686, method = "carsten")
  } #  100 runs +/- 0.1 alpha and power -- 

j <- 5
row5results.carsten <- foreach(i=1:length(nr.combn.list)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=nr.combn.list[[i]], pc=pc[j], pt=pt[j], runs=init.runs, alpha=init.alpha, power=init.power, seed=2686, method = "carsten")
  }  #  100 runs +/- 0.1 alpha and power -- 

########################### Finer (parallel) run:

final.runs <- 1000
final.alpha <- 0.15
final.power <- 0.8

row1results.carsten <- list.clean(row1results.carsten, is.null, recursive = F)
row2results.carsten <- list.clean(row2results.carsten, is.null, recursive = F)
row3results.carsten <- list.clean(row3results.carsten, is.null, recursive = F)
row4results.carsten <- list.clean(row4results.carsten, is.null, recursive = F)
row5results.carsten <- list.clean(row5results.carsten, is.null, recursive = F)

j <- 1
row1results.carsten.final <- foreach(i=1:length(row1results.carsten)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=row1results.carsten[[i]], pc=pc[j], pt=pt[j], runs=final.runs, alpha=final.alpha, power=final.power, seed=2686, method = "carsten")
  } # 2 hours

j <- 2
  row2results.carsten.final <- foreach(i=1:length(row2results.carsten)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=row2results.carsten[[i]], pc=pc[j], pt=pt[j], runs=final.runs, alpha=final.alpha, power=final.power, seed=2686, method = "carsten")
  } # 2 hours

j <- 3
  row3results.carsten.final <- foreach(i=1:length(row3results.carsten)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=row3results.carsten[[i]], pc=pc[j], pt=pt[j], runs=final.runs, alpha=final.alpha, power=final.power, seed=2686, method = "carsten")
  } # 1 hour

j <- 4
  row4results.carsten.final <- foreach(i=1:length(row4results.carsten)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=row4results.carsten[[i]], pc=pc[j], pt=pt[j], runs=final.runs, alpha=final.alpha, power=final.power, seed=2686, method = "carsten")
  } # 1 hour

j <- 5

  row5results.carsten.final <- foreach(i=1:length(row5results.carsten)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=row5results.carsten[[i]], pc=pc[j], pt=pt[j], runs=final.runs, alpha=final.alpha, power=final.power, seed=2686, method = "carsten")
  } # 15 mins


##################### Collate, and remove empty lists:

carsten.final <- list(row1results.carsten.final, row2results.carsten.final, row3results.carsten.final, row4results.carsten.final, row5results.carsten.final)

for(i in 1:5){
  carsten.final[[i]] <- list.clean(carsten.final[[i]], is.null, recursive = F)
  carsten.final[[i]] <- do.call(rbind, carsten.final[[i]])
}

####################### Remove dominated designs:
  carsten.admissible <- foreach(i=1:5) %dopar% {
    rmDominatedDesigns(carsten.final[[i]], "EssH0", "Ess", "n")
  }




########################## Find Chen designs: ###################
################# Initial (parallel) run

init.runs <- 100
init.alpha <- 0.25
init.power <- 0.7

j <- 1
row1results.chen <- foreach(i=1:length(nr.combn.list)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=nr.combn.list[[i]], pc=pc[j], pt=pt[j], runs=init.runs, alpha=init.alpha, power=init.power, seed=2686, method = "chen")
  } # 100 runs +/- 0.1 alpha and power -- 15mins

j <- 2
row2results.chen <- foreach(i=1:length(nr.combn.list)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=nr.combn.list[[i]], pc=pc[j], pt=pt[j], runs=init.runs, alpha=init.alpha, power=init.power, seed=2686, method = "chen")
  }  #  100 runs +/- 0.1 alpha and power -- 15mins

j <- 3
row3results.chen <- foreach(i=1:length(nr.combn.list)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=nr.combn.list[[i]], pc=pc[j], pt=pt[j], runs=init.runs, alpha=init.alpha, power=init.power, seed=2686, method = "chen")
  } #  100 runs +/- 0.1 alpha and power -- 15mins

j <- 4
row4results.chen <- foreach(i=1:length(nr.combn.list)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=nr.combn.list[[i]], pc=pc[j], pt=pt[j], runs=init.runs, alpha=init.alpha, power=init.power, seed=2686, method = "chen")
  } #  100 runs +/- 0.1 alpha and power -- 15mins

j <- 5
  row5results.chen <- foreach(i=1:length(nr.combn.list)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=nr.combn.list[[i]], pc=pc[j], pt=pt[j], runs=init.runs, alpha=init.alpha, power=init.power, seed=2686, method = "chen")
  } #  100 runs +/- 0.1 alpha and power -- 15mins

  
########################### Finer (parallel) run:

final.runs <- 10000
final.alpha <- 0.15
final.power <- 0.8

row1results.chen <- list.clean(row1results.chen, is.null, recursive = F)
row2results.chen <- list.clean(row2results.chen, is.null, recursive = F)
row3results.chen <- list.clean(row3results.chen, is.null, recursive = F)
row4results.chen <- list.clean(row4results.chen, is.null, recursive = F)
row5results.chen <- list.clean(row5results.chen, is.null, recursive = F)

j <- 1
  row1results.chen.final <- foreach(i=1:length(row1results.chen)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=row1results.chen[[i]], pc=pc[j], pt=pt[j], runs=final.runs, alpha=final.alpha, power=final.power, seed=2686, method = "chen")
  } # 3 hours

j <- 2
  row2results.chen.final <- foreach(i=1:length(row2results.chen)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=row2results.chen[[i]], pc=pc[j], pt=pt[j], runs=final.runs, alpha=final.alpha, power=final.power, seed=2686, method = "chen")
  } # 2 hours

j <- 3
  row3results.chen.final <- foreach(i=1:length(row3results.chen)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=row3results.chen[[i]], pc=pc[j], pt=pt[j], runs=final.runs, alpha=final.alpha, power=final.power, seed=2686, method = "chen")
  } # 2 hours

j <- 4
  row4results.chen.final <- foreach(i=1:length(row4results.chen)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=row4results.chen[[i]], pc=pc[j], pt=pt[j], runs=final.runs, alpha=final.alpha, power=final.power, seed=2686, method = "chen")
  } # 2 hours

j <- 5
  row5results.chen.final <- foreach(i=1:length(row5results.chen)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=row5results.chen[[i]], pc=pc[j], pt=pt[j], runs=final.runs, alpha=final.alpha, power=final.power, seed=2686, method = "chen")
  } # 2 hours

##################### Collate, and remove empty lists:

chen.final <- list(row1results.chen.final, row2results.chen.final, row3results.chen.final, row4results.chen.final, row5results.chen.final)

for(i in 1:5){
  chen.final[[i]] <- list.clean(chen.final[[i]], is.null, recursive = F)
  chen.final[[i]] <- do.call(rbind, chen.final[[i]])
}

####################### Remove dominated designs:
chen.admissible <- foreach(i=1:5) %dopar% {
    rmDominatedDesigns(chen.final[[i]], "EssH0", "Ess", "n")
  }






##################### Find Jung designs ##################### 

pc <- c(0.1, 0.2, 0.3, 0.4, 0.5)
pt <- pc+0.2

jung <- vector("list", 5)
  for(i in 1:5){
  jung[[i]] <- des_comp(J = 2,
                        type = "binomial",
                        alpha = 0.15,
                        beta = 0.2,
                        delta = 0.2,
                        pi_null = pc[i],
                        pi_alt = pc[i],
                        two_stage = list(equal = F,
                                         w=c(1, 0, 0, 0, 0),
                                         pi_ess = pc[i],
                                         efficacy = F,
                                         futility = T))$feasible
  } # 5 mins


########## Subset Jung's feasible results to just the admissible results:

jung.df <- vector("list", 5)
for(i in 1:5){
  jung.df[[i]] <- as.data.frame(jung[[i]])
}

jung.admissible <- foreach(i=1:5) %dopar% {
    rmDominatedDesigns(jung.df[[i]], "ESS0", "ESS1", "max(n)")
  } # 12 hours / 1.5m rows


############ COMBINE ALL DESIGNS INTO SINGLE LIST: #############


for(i in 1:5){
  jung.admissible[[i]]$narm <- jung.admissible[[i]]$`max(n)`/2
  jung.admissible[[i]] <- jung.admissible[[i]][, c("n01", "narm", "max(n)", "f1", "e2", "max alpha", "min power", "ESS0", "ESS1")]
  names(jung.admissible[[i]]) <- c("n1", "narm", "n", "r1", "r2", "typeIerr", "pwr", "EssH0", "Ess")
  jung.admissible[[i]]$theta0 <- NA
  jung.admissible[[i]]$theta1 <- NA
  jung.admissible[[i]]$design <- "jung"
  
  names(carsten.admissible[[i]])[2] <- "narm"
  names(chen.admissible[[i]])[2] <- "narm"
  carsten.admissible[[i]]$n <- 2*carsten.admissible[[i]]$narm
  chen.admissible[[i]]$n <- 2*chen.admissible[[i]]$narm
  carsten.admissible[[i]] <- carsten.admissible[[i]][, c("n1", "narm", "n", "r1", "r2", "typeIerr", "pwr", "EssH0", "Ess")]
  chen.admissible[[i]] <- chen.admissible[[i]][, c("n1", "narm", "n", "r1", "r2", "typeIerr", "pwr", "EssH0", "Ess")]
  carsten.admissible[[i]]$theta0 <- NA
  carsten.admissible[[i]]$theta1 <- NA
  chen.admissible[[i]]$theta0 <- NA
  chen.admissible[[i]]$theta1 <- NA
  carsten.admissible[[i]]$design <- "carsten"
  chen.admissible[[i]]$design <- "chen"
  
  block2.1e6[[i]]$n1 <- NA
  block2.1e6[[i]]$r1 <- NA
  block2.1e6[[i]] <- block2.1e6[[i]][, c("n1", "n", "eff.n", "r1", "r", "alpha", "power", "EssH0", "Ess", "theta0", "theta1")]
  block2.1e6[[i]]$design <- "block2"
    names(block2.1e6[[i]]) <- c("n1", "narm", "n", "r1", "r2", "typeIerr", "pwr", "EssH0", "Ess", "theta0", "theta1", "design")
  
  block8.1e6[[i]]$n1 <- NA
  block8.1e6[[i]]$r1 <- NA
  block8.1e6[[i]] <- block8.1e6[[i]][, c("n1", "n", "eff.n", "r1", "r", "alpha", "power", "EssH0", "Ess", "theta0", "theta1")]
  block8.1e6[[i]]$design <- "block8"
  names(block8.1e6[[i]]) <- c("n1", "narm", "n", "r1", "r2", "typeIerr", "pwr", "EssH0", "Ess", "theta0", "theta1", "design")
}


all.list <- vector("list", 5)
for(i in 1:5){
  all.list[[i]] <- rbind(jung.admissible[[i]], carsten.admissible[[i]], chen.admissible[[i]],  block2.1e6[[i]], block8.1e6[[i]])
}













################ Section 3.3: Based on example 1 from Jung (2008) #############
# Designs for block size 2 and block size 8:

block2.n <- seq(from=8, to=100, by=1)
block8.n <- seq(from=8, to=100, by=4)
max.combns <- 1e5

  ex1.100.1e5 <- foreach(i=block2.n, .combine = rbind) %dopar% {
    find2armDesigns(nmin=i, nmax=i, block.size=2, pc=0.7, pt=0.85, alpha=0.15, power=0.8, maxtheta0=0.85, mintheta1=0.7, bounds="ahern", max.combns=max.combns)
  }
  ex1.100.1e5.block8 <- foreach(i=block8.n, .combine = rbind) %dopar% {
    find2armDesigns(nmin=i, nmax=i, block.size=8, pc=0.7, pt=0.85, alpha=0.15, power=0.8, maxtheta0=0.85, mintheta1=0.7, bounds="ahern", max.combns=max.combns)
  }

ex1.100.1e5 <- rmDominatedDesigns(ex1.100.1e5, n = "eff.n")
ex1.100.1e5.block8 <- rmDominatedDesigns(ex1.100.1e5.block8, n = "eff.n")



# Designs for Carsten:

df1 <- findN1N2R1R2twoarm(nmin=5, nmax=200)
df1 <- df1[df1$r1>=0 & df1$r2>=0 & df1$r1<=df1$n1 & df1$r1<=df1$r2,]
# Split into lists, by r2:
df2 <- split(df1, list(df1$n1, df1$n, df1$r1))
# remove empty lists
nr.combn.list <- list.clean(df2, function(x) nrow(x) == 0L, recursive = F)
rm(df1, df2)

init.runs <- 100
init.alpha <- 0.25
init.power <- 0.7

  ex1.carsten.init <- foreach(i=1:length(nr.combn.list)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=nr.combn.list[[i]], pc=0.7, pt=0.85, runs=init.runs, alpha=init.alpha, power=init.power, seed=2686, method = "carsten")
  }# 100 runs +/- 0.1 alpha and power -- 

final.runs <- 10000
final.alpha <- 0.15
final.power <- 0.8

ex1.carsten.init <- list.clean(ex1.carsten.init, is.null, recursive = F)

j <- 1
ex1.carsten.final <- foreach(i=1:length(ex1.carsten.init)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=ex1.carsten.init[[i]], pc=0.7, pt=0.85, runs=final.runs, alpha=final.alpha, power=final.power, seed=2686, method = "carsten")
  }

do.call(rbind, ex1.carsten.final)


# Designs for Chen:
ex1.chen.init <- foreach(i=1:length(nr.combn.list)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=nr.combn.list[[i]], pc=0.7, pt=0.85, runs=init.runs, alpha=init.alpha, power=init.power, seed=2686, method = "chen")
  } # 100 runs +/- 0.1 alpha and power -- 

ex1.chen.init <- list.clean(ex1.chen.init, is.null, recursive = F)

# list of 139,000 -- far too much. Some lists have 20 rows.
# Create subset of trials with max N per arm <=100:

index <- lapply(ex1.chen.init, function(x) x$n[1]<=100)
sum(unlist(index))
ex1.chen.init100  <- subset(ex1.chen.init, n[1]<=100)
length(ex1.chen.init100)

final.runs <- 10000

j <- 1
  ex1.chen.final <- foreach(i=1:length(ex1.chen.init100)) %dopar% {
    findCarstenChenTypeITypeIIRmRows(nr.list=ex1.chen.init100[[i]], pc=0.7, pt=0.85, runs=final.runs, alpha=final.alpha, power=final.power, seed=2686, method = "chen")
  }

ex1.chen.final <- do.call(rbind, ex1.chen.final)
ex1.chen.final <- rmDominatedDesigns(ex1.chen.final)
nrow(ex1.chen.final)




# What is the expected number of successes on Carsten when pt=pc=p0
p0 <- seq(0.1, 0.5, by=0.1)
q0 <- 1-p0
# P(success) = success on trt AND failure on control, ie pt*qc, ie p0*q0
# P(any one pair is a success)=p0*q0
p0*q0
# For p0=0.1, h0-opt design has narm=16,
# p0=0.2 -> 29
# p0=0.3 -> 44
# p0=0.4 -> 48
# p0=0.5 -> 83
n <- c(16, 29, 44, 48, 83)
n*p0*q0

#### What about the E(successes) using our approach
# P(one success on treatment) = p0
# P(one success on control) = q0
n.block2 <- c(52, 48, 58, 60, 60)
p0*n.block2/2 + q0*n.block2/2
n.block2 * (p0+q0)/2
2*p0*q0 + p0^2 + q0^2


# load("example1_14jun.RData")
# Block designs:
# Minimax, H0-opt and H1-opt:
block2.mini <- ex1.100.1e5[ex1.100.1e5$n==min(ex1.100.1e5$n),]
block2.h0.opt <- ex1.100.1e5[ex1.100.1e5$EssH0==min(ex1.100.1e5$EssH0),]
block2.h1.opt <- ex1.100.1e5[ex1.100.1e5$Ess==min(ex1.100.1e5$Ess),]
block2 <- rbind(block2.mini, block2.h0.opt, block2.h1.opt)
block2$opt <- c("mini", "h0", "h1")
block2$design <- "Block 2"
names(block2)[4:5] <- c("typeIerr", "pwr")

block8.mini <- ex1.100.1e5.block8[ex1.100.1e5.block8$n==min(ex1.100.1e5.block8$n),]
block8.h0.opt <- ex1.100.1e5.block8[ex1.100.1e5.block8$EssH0==min(ex1.100.1e5.block8$EssH0),]
block8.h1.opt <- ex1.100.1e5.block8[ex1.100.1e5.block8$Ess==min(ex1.100.1e5.block8$Ess),]
block8 <- rbind(block8.mini, block8.h0.opt, block8.h1.opt)
block8$opt <- c("mini", "h0", "h1")
block8$design <- "Block 8"
names(block8)[4:5] <- c("typeIerr", "pwr")

chen.mini <- ex1.chen.final[ex1.chen.final$n==min(ex1.chen.final$n),][1,]
chen.h0.opt <- ex1.chen.final[ex1.chen.final$EssH0==min(ex1.chen.final$EssH0),]
chen.h1.opt <- ex1.chen.final[ex1.chen.final$Ess==min(ex1.chen.final$Ess),]
chen <- rbind(chen.mini, chen.h0.opt, chen.h1.opt)
chen$opt <- c("mini", "h0", "h1")
names(chen)[names(chen)=="r2"] <- "r"
chen$design <- "Chen"


# Designs for Jung:
# From Jung's paper:
# Minimax: (n1, n, a1, a)=31, 63, -1, 6) EN=52.16 ie 104.32
# H0-optimal: (n1, n, a1, a)=(27, 73, 1, 6) EN=47.28 ie 94.56
# load("ex1_jung.RData")
jung.mini <- findSingle2arm2stageJungDesignFast(n1=31, n2=63-31, n=63, a1=-1, r2=6, p0=0.7, p1=0.85, alpha=0.15, power=0.8)
jung.h0.opt <- findSingle2arm2stageJungDesignFast(n1=27, n2=73-27, n=73, a1=1, r2=6, p0=0.7, p1=0.85, alpha=0.15, power=0.8)
# From Michael Grayling's code:
jung.h1.opt <- findSingle2arm2stageJungDesignFast(n1=56, n2=62-56, n=62, a1=5, r2=5, p0=0.7, p1=0.85, alpha=0.15, power=0.8)

jung <- data.frame(rbind(jung.mini, jung.h0.opt, jung.h1.opt))
jung$n <- jung$max.n./2  
jung$r1 <- jung$f1-1
jung <- jung[, c("n01", "n", "r1", "e2", "max.alpha", "min.power", "ESS0", "ESS1")]
jung$opt <- c("mini", "h0", "h1")

names(jung) <- c("n1", "n", "r1", "r", "typeIerr", "pwr", "EssH0", "Ess", "opt")
jung$design <- "Jung"





## Put results together:
ex1.df <- merge(jung, chen, all=T)
ex1.df <- merge(ex1.df, block2, all=T)
ex1.df <- merge(ex1.df, block8, all=T)

ex1.df$design <- factor(ex1.df$design, levels = c("Jung", "Chen", "Block 2", "Block 8"))
ex1.df$opt <- factor(ex1.df$opt, levels = c("h0", "h1", "mini"))

ex1.df <- ex1.df[order(ex1.df$opt, ex1.df$design), ]
ex1.df$n.total <- 2*ex1.df$n
rownames(ex1.df) <- paste(ex1.df$design, ex1.df$opt)






######## Double check example 1: ###########
# 
# registerDoParallel(cores)
# maxr <- ceiling(64*0.85)
# time.ex1.100.1e5 <- system.time({
#   check <- foreach(i=1:maxr, .combine = rbind) %dopar% {
#     find2armDesigns(nmin=62, nmax=65, block.size=2, pc=0.7, pt=0.85, alpha=0.15, power=0.8, maxtheta0=0.85, mintheta1=0.7, bounds=NULL, fixed.r=i, max.combns=1e5)
#   }
# })[3]
# stopImplicitCluster()
# 
# check
# 
# registerDoParallel(cores)
# time.ex1.100.1e5 <- system.time({
#   check.block8 <- foreach(i=1:maxr, .combine = rbind) %dopar% {
#     find2armDesigns(nmin=60, nmax=64, block.size=8, pc=0.7, pt=0.85, alpha=0.15, power=0.8, maxtheta0=0.85, mintheta1=0.7, bounds=NULL, fixed.r=i, max.combns=1e5)
#   }
# })[3]
# stopImplicitCluster()








