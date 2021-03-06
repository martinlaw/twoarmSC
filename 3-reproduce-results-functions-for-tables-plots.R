####
#### These functions are used to run the code that analyses the designs used.
####

library(gridExtra)
library(grid)
library(gdata)
library(ggplot2)
library(ggthemes)
library(xtable)
library(rpact)


#### FIGURE 2: ####

findDifferences <- function(loss.df, designs.df){
  min.loss <-  aggregate(loss.df, list(designs.df$design), min)
  name <- min.loss[,1]
  min.loss <- t(min.loss[,-1])
  colnames(min.loss) <- name
  min.loss <- as.data.frame(min.loss)
  head(min.loss)
  
  # find differences:
  min.loss$block2_block8 <- min.loss$block2 - min.loss$block8
  min.loss$block2_carsten <- min.loss$block2 - min.loss$carsten
  min.loss$block2_chen <- min.loss$block2 - min.loss$chen
  min.loss$block2_jung <- min.loss$block2 - min.loss$jung
  
  min.loss$block8_carsten <- min.loss$block8 - min.loss$carsten
  min.loss$block8_chen <- min.loss$block8 - min.loss$chen
  min.loss$block8_jung <- min.loss$block8 - min.loss$jung
  
  min.loss$carsten_chen <- min.loss$carsten - min.loss$chen
  min.loss$carsten_jung <- min.loss$carsten - min.loss$jung
  
  min.loss$chen_jung <- min.loss$chen - min.loss$jung
  min.loss
}

#################### SECTION 3.2: MISSPECIFICATION OF RESPONSE RATES -- FIGS 3 AND 4 ################

plotPower <- function(x, power.ss="power"){
  theplot <- ggplot(x, aes(true.pt, true.pc)) +
    theme_tufte() +
    geom_raster(aes(fill = switch(power.ss, "power" = rejectH0, "ss" = ss))) +
    scale_fill_gradient2(low="yellow", high="red", guide="colorbar") +
    geom_text(aes(true.pt, true.pc, label = switch(power.ss, "power" = round(rejectH0, 2), "ss" = round(ss,1))), color = "black", size = 5) +
    labs(fill = switch(power.ss, "power" = "P(reject H0)", "ss" = "ESS"),
         y=expression("True p"[C]),
         x=expression("True p"[T])) +
    theme(axis.text=element_text(size=rel(1.25)),
          axis.title=element_text(size=rel(1.5)))
  theplot
}

plotPowerBoth <- function(x, power.ss="power"){
  theplot <- ggplot(x, aes(true.pt, true.pc)) +
    theme_tufte() +
    ggtitle(bquote(paste("P(reject H"[0], "), H"[0], "-optimal design for p"[0], "=",.(x$pc[1]), ", p"[1], "=", .(x$pt[1])))) +
    geom_raster(aes(fill = switch(power.ss, "power" = rejectH0, "ss" = ss))) +
    scale_fill_gradient2(low="yellow", high="red", guide="colorbar") +
    geom_text(aes(true.pt, true.pc,
                  label = switch(power.ss, "power" = round(rejectH0, 2), "ss" = round(ss,1))),
              color = "black",
              size = rel(4)) +
    labs(fill = switch(power.ss, "power" = expression("P(reject H"[0],")"), "ss" = "ESS"),
         y=expression("True p"[C]),
         x=expression("True p"[T])) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.text=element_text(size=rel(1.25)),
          axis.title=element_text(size=rel(1.25)),
          strip.text = element_text(size=rel(1.25))
    )+
    facet_grid(. ~ design)
  theplot
}








# Write a program that will find the sample size using our design and Carsten's design, for a given set of data #

# for p0=0.1, p1=0.3, alpha=0.15, power=0.8, h0-optimal designs are:
# Carsten: 
# n1=7; n=16; r1=1; r2=3
# This design should have the following OCs:
# Type I error: 0.148, Power: 0.809, EssH0: 21.27400, EssH1: 19.44200
#
# Our design (not strictly H0-optimal, but is within 0.5 of optimal wrt EssH0 and has N=31, vs N=71 for the actual H0-optimal):
# n=31; r2=3; theta0=0.1277766; theta1=0.9300000

# This function, varypPlot, requires others: 
########## Function to find CP matrix for our design: ##########
findBlockCP <- function(n, r, pc, pt, theta0, theta1){
  
  Bsize <- 2
  pat.cols <- seq(from=2*n, to=2, by=-2)[-1]
  qc <- 1-pc
  qt <- 1-pt
  prob0 <- qt*pc
  prob1 <- pt*pc + qt*qc
  prob2 <- pt*qc
  prob.vec <- c(prob0, prob1, prob2)
  
  # CREATE UNCURTAILED MATRIX
  mat <- matrix(3, ncol=2*n, nrow=min(n+r+Bsize+2, 2*n+1))
  rownames(mat) <- 0:(nrow(mat)-1)
  mat[(n+r+2):nrow(mat),] <- 1 
  mat[1:(n+r+1),2*n] <- 0 # Fail at end
  
  for(i in (n+r+1):1){  
    for(j in pat.cols){  # Only look every C patients (no need to look at final col)
      if(i-1<=j){ # Condition: Sm<=m
        mat[i,j] <- ifelse(test=j-(i-1) > n-r+1, yes=0, no=sum(prob.vec*mat[i:(i+Bsize), j+Bsize])) 
        # IF success is not possible (i.e. [total no. of pats-Xa+Ya-Xb] > n-r+1), THEN set CP to zero. Otherwise, calculate it based on "future" CPs.
      } 
    }
  }
  
  for(i in 3:nrow(mat)){
    mat[i, 1:(i-2)] <- NA
  }
  uncurt <- mat
  
  ### CREATE CURTAILED MATRIX
  for(i in (n+r+1):1){  
    for(j in pat.cols){  # Only look every Bsize patients (no need to look at final col)
      if(i-1<=j){ # Condition: Sm<=m
        newcp <- sum(prob.vec*mat[i:(i+Bsize), j+Bsize])
        if(newcp > theta1) mat[i,j] <- 1
        if(newcp < theta0) mat[i,j] <- 0
        if(newcp <= theta1 & newcp >= theta0) mat[i,j] <- newcp
      } 
    }
  }
  return(mat)
} # END OF FUNCTION


######################## FUNCTION TO FIND BASIC CP MATRIX (CP=0/1/neither) FOR CARSTEN ##################
carstenCPonly <- function(n1, n, r1, r, pair=FALSE){
  cpmat <- matrix(0.5, ncol=2*n, nrow=r+1)
  rownames(cpmat) <- 0:(nrow(cpmat)-1)
  cpmat[r+1, (2*r):(2*n)] <- 1 # Success is a PAIR of results s.t. Xt-Xc=1
  cpmat[1:r,2*n] <- 0 # Fail at end
  pat.cols <- seq(from=(2*n)-2, to=2, by=-2)
  for(i in nrow(cpmat):1){  
    for(j in pat.cols){  # Only look every C patients (no need to look at final col)
      if(i-1<=j/2){ # Condition: number of successes must be <= pairs of patients so far
        # IF success is not possible (i.e. [total no. of failures] > n1-r1+1 at stage 1 or [total no. of failures] > n-r+1), THEN set CP to zero:
        if((i<=r1 & j/2 - (i-1) >= n1-r1+1) | (j/2-(i-1) >= n-r+1)){
          cpmat[i,j] <- 0
        }
      } else{
        cpmat[i,j] <- NA # impossible to have more successes than pairs
      }
    }
  }
  
  if(pair==TRUE){
    cpmat <- cpmat[, seq(2, ncol(cpmat), by=2)]
  }
  return(cpmat)
} # END OF FUNCTION



varypPlot <- function(h0=TRUE, n, r, theta0=NULL, theta1=NULL, pc, pt, true.pc, true.pt,
                      runs, seed=sample(1:1e5, 1), method){
    set.seed(seed)  

  
  #### Simulate data:
  max.n <- max(unlist(n))
  h1.dat <- matrix(runif(max.n*runs), ncol=max.n, byrow=TRUE)
  # Matrix of successes (0, 1 or 2 successes in each pair):
  if(h0==TRUE){
    true.pt <- true.pc
  }
  
  true.qc <- 1-true.pc
  true.qt <- 1-true.pt

  # Find CPs for block and Carsten designs:
  if(method=="block"){
    block.mat <- findBlockCP(n=n, r=r, pc=pc, pt=pt, theta0=theta0, theta1=theta1)
    
  # Find lower and upper stopping boundaries for block and Carsten designs: 
    lower <- rep(-Inf, ncol(block.mat)/2)
    upper <- rep(Inf, ncol(block.mat)/2)
    looks <- seq(2, ncol(block.mat), by=2)
    for(j in 1:length(looks)){
      if(any(block.mat[,looks[j]]==0, na.rm = T)){
        lower[j] <- max(which(block.mat[,looks[j]]==0))-1 # Minus 1 to account for row i == [number of responses-1]
      } 
      if(any(block.mat[,looks[j]]==1, na.rm = T)){
        upper[j] <- which.max(block.mat[,looks[j]])-1 # Minus 1 to account for row i == [number of responses-1]
      }
    }
    
    block.success.mat <- matrix(rep(0, max.n*runs), ncol=max.n)
    block.success.mat[h1.dat < true.pt*true.qc] <- 2
    block.success.mat[h1.dat >= true.pt*true.qc & h1.dat < true.pt*true.qc + true.pt*true.pc + true.qt*true.qc] <- 1
    # Cumulative successes:
    block.success.mat.cum <- t(apply(block.success.mat, 1, cumsum))
    
    #### Note where a boundary is reached for each simulation, for block and Carsten designs -- this is the number of pairs observed before stopping, i.e. the sample size: 
    block.ss.arm <- apply(block.success.mat.cum, 1, function(x) match(TRUE, x <= lower | x >= upper))
    block.ss.mean <- 2*sum(block.ss.arm)/runs
    
    block.lower.match <- apply(block.success.mat.cum, 1, function(x) match(TRUE, x <= lower, nomatch=1e5))
    block.upper.match <- apply(block.success.mat.cum, 1, function(x) match(TRUE, x >= upper, nomatch=1e5))
    block.prob.reject.h0 <- sum(block.upper.match < block.lower.match)/runs
    block.se <- sqrt(block.prob.reject.h0*(1-block.prob.reject.h0))/runs
    
    ss <- block.ss.mean
    rejectH0 <- block.prob.reject.h0
  }
  
  if(method=="carsten"){
    carsten.mat <- carstenCPonly(n1=n[[1]], n=n[[2]], r1=r[[1]], r=r[[2]], pair=F)
    lower.carsten <- rep(-Inf, ncol(carsten.mat)/2)
    upper.carsten <- rep(Inf, ncol(carsten.mat)/2)
    looks.carsten <- seq(2, ncol(carsten.mat), by=2)
    for(j in 1:length(looks.carsten)){
      if(any(carsten.mat[,looks.carsten[j]]==0, na.rm = T)){
        lower.carsten[j] <- max(which(carsten.mat[,looks.carsten[j]]==0))-1 # Minus 1 to account for row i == [number of responses-1]
      } 
      if(any(carsten.mat[,looks.carsten[j]]==1, na.rm = T)){
        upper.carsten[j] <- which.max(carsten.mat[,looks.carsten[j]])-1 # Minus 1 to account for row i == [number of responses-1]
      }
    }
    # Matrix of successes for Carsten:
    carsten.success.mat <- matrix(rep(0, max.n*runs), ncol=max.n)
    carsten.success.mat[h1.dat < true.pt*true.qc] <- 1
    
    # Cumulative successes:
    carsten.success.mat.cum <- t(apply(carsten.success.mat, 1, cumsum))
    
    #### Note where a boundary is reached for each simulation, for block and Carsten designs -- this is the number of pairs observed before stopping, i.e. the sample size: 
    carsten.ss.arm <- apply(carsten.success.mat.cum, 1, function(x) match(TRUE, x <= lower.carsten | x >= upper.carsten))
    
    carsten.ss.mean <- 2*sum(carsten.ss.arm)/runs
    
    carsten.lower.match <- apply(carsten.success.mat.cum, 1, function(x) match(TRUE, x <= lower.carsten, nomatch=1e5))
    carsten.upper.match <- apply(carsten.success.mat.cum, 1, function(x) match(TRUE, x >= upper.carsten, nomatch=1e5))
    carsten.prob.reject.h0 <- sum(carsten.upper.match < carsten.lower.match)/runs
    carsten.se <- sqrt(carsten.prob.reject.h0*(1-carsten.prob.reject.h0))/runs
    
    ss <- carsten.ss.mean
    rejectH0 <- carsten.prob.reject.h0
  }
  
  output <- data.frame(ss=ss, rejectH0=rejectH0, pc=pc, pt=pt, true.pc=true.pc, true.pt=true.pt)
  return(output)
} # END OF FUNCTION varypPlot



varypPowerSS <- function(method, n, r, pc, pt, theta0=NULL, theta1=NULL, p.range, runs, lower.triangle=TRUE,
                         seed=sample(1:1e5, 1), power.ss="power", title.text=FALSE, plot=TRUE){
  df <- data.frame(expand.grid(seq(p.range[1], p.range[2], by=0.1), seq(p.range[1], p.range[2], by=0.1)))
  names(df) <- c("true.pt", "true.pc")
  if(lower.triangle==TRUE){
    df <- df[df$true.pt >= df$true.pc, ]
  }
  x <- apply(df, 1, function(x) varypPlot(h0=FALSE, n=n, r=r, pc=pc, pt=pt, true.pc=x["true.pc"], true.pt=x["true.pt"], theta0=theta0, theta1=theta1, runs=runs, seed=seed, method=method))
  x <- do.call(rbind, x)
  
  theplot <- ggplot(x, aes(true.pt, true.pc)) +
    geom_raster(aes(fill = switch(power.ss, "power" = rejectH0, "ss" = ss))) +
    scale_fill_gradient2(low="yellow", high="red", guide="colorbar") +
    geom_text(aes(true.pt, true.pc, label = switch(power.ss, "power" = round(rejectH0, 2), "ss" = round(ss,1))), color = "black", size = 5) +
    labs(fill = switch(power.ss, "power" = "P(reject H0)", "ss" = "ESS"),
         y=expression("True p"[C]),
         x=expression("True p"[T])) +
    theme(legend.position = c(0.2, 0.7),
          axis.text=element_text(size=rel(1.25)),
          axis.title=element_text(size=rel(1.5)))
  
  if(title.text==TRUE){
    if(method=="carsten"){
      theplot <- theplot + 
        labs(title = paste("Carsten design: n1=", n[[1]], ", n=", n[[2]], ", r1=", r[[1]], ", r=", r[[2]], ". pc=", pc, ", pt=", pt, sep=""))
    }
    
    if(method=="block"){
      theplot <- theplot + 
        labs(title = paste("Block design: n=", n, ", r=", r, ", theta0=", round(theta0,3), ", theta1=", round(theta1,3),  ". pc=", pc, ", pt=", pt, sep=""))
    }
  }
  
  if(plot==TRUE){
    return(theplot)
  } else {
    return(x)
  }
} # END OF FUNCTION varypPowerSS

# Write a program that will find the sample size using our design and Carsten's design, for a given set of data #

# for p0=0.1, p1=0.3, alpha=0.15, power=0.8, h0-optimal designs are:
# Carsten: 
# n1=7; n=16; r1=1; r2=3
# This design should have the following OCs:
# Type I error: 0.148, Power: 0.809, EssH0: 21.27400, EssH1: 19.44200
#
# Our design (not strictly H0-optimal, but is within 0.5 of optimal wrt EssH0 and has N=31, vs N=71 for the actual H0-optimal):
# n=31; r2=3; theta0=0.1277766; theta1=0.9300000

rejectionRegions <- function(n, r, theta0=NULL, theta1=NULL, pc, pt, method=NULL){
  
  library(ggplot2)

  #### Find CPs for block and Carsten designs:
  if(method=="block"){
    block.mat <- findBlockCP(n=n, r=r, pc=pc, pt=pt, theta0=theta0, theta1=theta1)
    #### Find lower and upper stopping boundaries for block and Carsten designs:
    lower <- rep(-Inf, ncol(block.mat)/2)
    upper <- rep(Inf, ncol(block.mat)/2)
    looks <- seq(2, ncol(block.mat), by=2)
    for(j in 1:length(looks)){
      if(any(block.mat[,looks[j]]==0, na.rm = T)){
        lower[j] <- max(which(block.mat[,looks[j]]==0))-1 # Minus 1 to account for row i == [number of responses-1]
      } 
      if(any(block.mat[,looks[j]]==1, na.rm = T)){
        upper[j] <- which.max(block.mat[,looks[j]])-1 # Minus 1 to account for row i == [number of responses-1]
      }
    }
    
    m.list <- vector("list", max(looks))
    for(i in 1:length(looks)){
      m.list[[i]] <- data.frame(m=looks[i], successes=0:max(looks), outcome=NA)
      m.list[[i]]$outcome[m.list[[i]]$successes <= lower[i]] <- "NOGO"
      m.list[[i]]$outcome[m.list[[i]]$successes >= lower[i]+1 & m.list[[i]]$successes <= upper[i]-1] <- "Continue"
      m.list[[i]]$outcome[m.list[[i]]$successes >= upper[i] & m.list[[i]]$successes <= looks[i]] <- "GO"
      m.list[[i]]$outcome[m.list[[i]]$successes > looks[i]] <- NA
    }
  }
  
  
  
  if(method=="carsten"){
    carsten.mat <- carstenCPonly(n1=n[[1]], n=n[[2]], r1=r[[1]], r=r[[2]], pair=F)
    #### Find lower and upper stopping boundaries for block and Carsten designs:
    lower.carsten <- rep(-Inf, ncol(carsten.mat)/2)
    upper.carsten <- rep(Inf, ncol(carsten.mat)/2)
    looks.carsten <- seq(2, ncol(carsten.mat), by=2)
    for(j in 1:length(looks.carsten)){
      if(any(carsten.mat[,looks.carsten[j]]==0, na.rm = T)){
        lower.carsten[j] <- max(which(carsten.mat[,looks.carsten[j]]==0))-1 # Minus 1 to account for row i == [number of responses-1]
      } 
      if(any(carsten.mat[,looks.carsten[j]]==1, na.rm = T)){
        upper.carsten[j] <- which.max(carsten.mat[,looks.carsten[j]])-1 # Minus 1 to account for row i == [number of responses-1]
      }
    }
    
    looks <- looks.carsten
    lower <- lower.carsten
    upper <- upper.carsten
    
    #browser()
    
    m.list <- vector("list", max(looks))
    for(i in 1:length(looks)){
      m.list[[i]] <- data.frame(m=looks[i], successes=0:(max(looks)/2), outcome=NA)
      m.list[[i]]$outcome[m.list[[i]]$successes <= lower[i]] <- "NOGO"
      m.list[[i]]$outcome[m.list[[i]]$successes >= lower[i]+1 & m.list[[i]]$successes <= upper[i]-1] <- "Continue"
      m.list[[i]]$outcome[m.list[[i]]$successes >= upper[i] & m.list[[i]]$successes <= looks[i]/2] <- "GO"
      m.list[[i]]$outcome[m.list[[i]]$successes > looks[i]/2] <- NA
      
    }
  }
  
  m.df <- do.call(rbind, m.list)
  
  ggplot(m.df, aes(m, successes)) +
    geom_raster(aes(fill = outcome)) +
    theme(legend.position = c(0.2, 0.7))
}

##################### END OF SECTION 3.2 ###################


# Plot thetaE, thetaF of all feasible designs ####
plotFeasible <- function(df, criterion=0){
  library(ggplot2)
  theme_set(theme_bw(base_size = 20)) # Increase font size and set theme for plots.
  plt <- ggplot(data=df,
                aes(x = theta0, y = theta1))
  if(criterion==0){
    main.text <- expression(paste(theta[F], ", ", theta[E], " and ESS(", p[0], ", ", p[0], ") for all feasible designs"))
    legend.text <- expression(paste("ESS(", p[0], ", ", p[0], ")"))
    plt <- plt + geom_point(aes(color=EssH0), size=3)
  }
  if(criterion==1){
    main.text <- expression(paste(theta[F], ", ", theta[E], " and ESS(", p[0], ", ", p[1], ") for all feasible designs"))
    legend.text <- expression(paste("ESS(", p[0], ", ", p[1], ")")) 
    plt <- plt + geom_point(aes(color=Ess), size=3)
  }
  plt <- plt + 
    labs(title=main.text,
         x=expression(theta[F]),
         y=expression(theta[E]),
         col=legend.text)
  return(plt)
}


# Plot rejection regions for proposed design and Carsten design: ####
plotRejectionRegions <- function(data.f, method){
  library(ggplot2)
  theme_set(theme_bw(base_size = 20))# Increase font size and set theme for plots.
  rejection.plot <- ggplot(data.f, aes(m, successes)) +
    geom_raster(aes(fill = outcome)) +
    theme(legend.position = c(0.2, 0.8))+
    labs(fill="Decision")
  if(method=="carsten"){
    rejection.plot <- rejection.plot +
      labs(title="Decision space for Carsten design",
           x="Participants (total)",
           y="Successes (pairs)")
  }
  if(method=="block"){
    rejection.plot <- rejection.plot +
      labs(title="Decision space for block size two",
           x="Participants (total)",
           y="Successes (total)")
  }
  return(rejection.plot)
}
