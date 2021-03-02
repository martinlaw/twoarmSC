source("1-reproduce-results-functions-find-designs.R")


pc <- c(0.1, 0.2, 0.3, 0.4, 0.5)
pt <- pc + 0.2

# NOTE: 30 cores were used for original timings.
cores <- detectCores()-2



registerDoParallel(cores)
timerow1.to60.1e6.noskip <- system.time({
  row1.to60.1e6.noskip <- foreach(i=block2.n, .combine = rbind) %dopar% {
    findSCdesSlow(nmin=16, nmax=40, block.size=8, pc=0.1, pt=0.3, alpha=0.15, power=0.8, maxtheta0=0.3, mintheta1=0.7, bounds="wald", max.combns=1e3)
  }
})[3] 

set.seed(1)
system.time({
  x <- findSCdes(nmin=32, nmax=56, block.size=8, pc=0.1, pt=0.3, alpha=0.15, power=0.8, maxtheta0=0.3, mintheta1=0.7, bounds="wald", max.combns=1e3)
})
set.seed(1)
system.time({
  xx <- findSCdesXX(nmin=32, nmax=56, block.size=8, pc=0.1, pt=0.3, alpha=0.15, power=0.8, maxtheta0=0.3, mintheta1=0.7, bounds="wald", max.combns=1e3)
})
set.seed(1)
system.time({
  xxx <- findSCdesXXX(nmin=32, nmax=56, block.size=8, pc=0.1, pt=0.3, alpha=0.15, power=0.8, maxtheta0=0.3, mintheta1=0.7, bounds="wald", max.combns=1e3)
})

for(i in 1:nrow(x)) x[i,]==xx[i,]

set.seed(1)
system.time({
  x2 <- find2armBlockDesigns(nmin=32, nmax=56, block.size=8, pc=0.1, pt=0.3, alpha=0.15, power=0.8, maxtheta0=0.3, mintheta1=0.7, bounds="ahern", max.combns=1e4)
})
set.seed(1)
system.time({
  x2x <- find2armBlockDesignsX(nmin=32, nmax=56, block.size=8, pc=0.1, pt=0.3, alpha=0.15, power=0.8, maxtheta0=0.3, mintheta1=0.7, bounds="ahern", max.combns=1e4)
})
set.seed(1)
system.time({
  xxxfast <- findSCdesXXX(nmin=32, nmax=56, block.size=8, pc=0.1, pt=0.3, alpha=0.15, power=0.8, maxtheta0=0.3, mintheta1=0.7, bounds="ahern", max.combns=1e4, fast.method = TRUE)
})
x2x
xxxfast
