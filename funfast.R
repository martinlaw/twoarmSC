find2armBlockDesignsX <- function(nmin,
                                 nmax,
                                 block.size,
                                 pc,
                                 pt,
                                 alpha,
                                 power,
                                 maxtheta0=NULL,
                                 mintheta1=0.7,
                                 bounds=NULL,
                                 fixed.r=NULL,
                                 max.combns=1e6)
{
  require(tcltk)
  require(data.table)
  
  Bsize <- block.size 
  
  if(Bsize%%2!=0) stop("Block size must be an even number")
  
  if((2*nmin)%%Bsize!=0) stop("2*nmin must be a multiple of block size")
  if((2*nmax)%%Bsize!=0) stop("2*nmax must be a multiple of block size")
  
  nposs <- seq(from=nmin, to=nmax, by=Bsize/2)
  
  qc <- 1-pc
  qt <- 1-pt
  
  ################ Function for finding the Prob(reponses on treatment + non-responses on control)=0, 1, 2,... Bsize:
  findProbVec <- function(Bsize, pt=pt, qt=qt, pc=pc, qc=qc){
    prob.vec <- rep(NA, Bsize+1)
    for(i in 1:(Bsize+1)){
      positives <- i-1
      full.vec <- expand.grid(rep(list(0:1), Bsize))
      positive.mat <- full.vec[rowSums(full.vec) == positives,]
      negative.mat <- -1*(positive.mat-1)
      
      positive.vec <- rep(c(pt,qc), each=Bsize/2)
      negative.vec <- rep(c(qt,pc), each=Bsize/2)
      
      posneg.mat <- t(t(positive.mat)*positive.vec) + t(t(negative.mat)*negative.vec)
      prob.vec[i] <- sum(apply(posneg.mat, 1, prod))
    }
    if(sum(prob.vec)-1 > 1e-8) stop("Probabilities do not sum to 1.")
    prob.vec
  }
  
  ################ Function for finding the uncurtailed CP matrix:
  findBlock2armUncurtailedMatrix <- function(n, r, Bsize, pat.cols, prob.vec){
    
    cpmat <- matrix(3, ncol=2*n, nrow=min(n+r+Bsize+2, 2*n+1))
    rownames(cpmat) <- 0:(nrow(cpmat)-1)
    cpmat[(n+r+2):nrow(cpmat),] <- 1 
    cpmat[1:(n+r+1),2*n] <- 0 # Fail at end
    
    for(i in (n+r+1):1){  
      for(j in pat.cols){  # Only look every C patients (no need to look at final col)
        if(i-1<=j){ # Condition: Sm<=m
          cpmat[i,j] <- ifelse(test=j-(i-1) >= n-r+1, yes=0, no=sum(prob.vec*cpmat[i:(i+Bsize), j+Bsize])) 
          # IF success is not possible (i.e. [total no. of pats-Xa+Ya-Xb] >= n-r+1), THEN set CP to zero. Otherwise, calculate it based on "future" CPs.
        } 
      }
    }
    
    for(i in 3:nrow(cpmat)){
      cpmat[i, 1:(i-2)] <- NA
    }
    cpmat
  }
  
  
  
  prob.vec <- findProbVec(Bsize=Bsize, pt=pt, qt=qt, pc=pc, qc=qc)
  prob.vec.p0 <- findProbVec(Bsize=Bsize, pt=pc, qt=qc, pc=pc, qc=qc)
  
  pat.cols.list <- lapply(nposs, function(x) seq(from=2*x, to=Bsize, by=-Bsize)[-1])
  names(pat.cols.list) <- nposs
  
  if(is.null(maxtheta0)){
    maxtheta0 <- pt
  }
  
  r.list <- list()
  for(i in 1:length(nposs))
  {
    r.list[[i]] <- 0:(nposs[i]-2) # r values: 0 to nposs[i]-2
  }
  
  ns <- NULL
  
  for(i in 1:length(nposs)){
    ns <- c(ns, rep(nposs[i], length(r.list[[i]])))
  }
  
  sc.subset <- data.frame(n=ns, r=unlist(r.list))
  
  if(!is.null(bounds)){
    # Incorporate A'Hern's bounds:
    if(bounds=="ahern")  { 
      #sc.subset <- sc.subset[sc.subset$r >= pc*sc.subset$n & sc.subset$r <= pt*sc.subset$n, ] # One-arm case
      sc.subset <- sc.subset[sc.subset$r >= 1 & sc.subset$r <= pt*sc.subset$n, ] # Try this for two-arm case -- interval [1, pt*Narm]
    }
    
    if(bounds=="wald"){
      # Even better to incorporate Wald's bounds:
      denom <- log(pt/pc) - log((1-pt)/(1-pc))
      accept.null <- log((1-power)/(1-alpha)) / denom  + nposs * log((1-pc)/(1-pt))/denom
      accept.null <- floor(accept.null)
      
      reject.null <- log((power)/alpha) / denom  + nposs * log((1-pc)/(1-pt))/denom
      reject.null <- ceiling(reject.null)
      
      r.wald <- NULL
      ns.wald <- NULL
      for(i in 1:length(nposs)){
        r.wald <- c(r.wald, accept.null[i]:reject.null[i])
        ns.wald <- c(ns.wald, rep(nposs[i], length(accept.null[i]:reject.null[i])))
      }
      sc.subset <- data.frame(n=ns.wald, r=r.wald)
      sc.subset <- sc.subset[sc.subset$n - sc.subset$r >=2, ]
    }
  }
  
  # In case you want to specify values for r:
  if(!is.null(fixed.r))  {
    sc.subset <- sc.subset[sc.subset$r %in% fixed.r,]
  }
  
  ###### Find thetas for each possible {r, N} combn:
  
  
  mat.list <- vector("list", nrow(sc.subset))
  
  for(i in 1:nrow(sc.subset)){
    mat.list[[i]] <- findBlock2armUncurtailedMatrix(n=sc.subset[i,"n"], r=sc.subset[i,"r"], Bsize=Bsize, pat.cols=pat.cols.list[[paste(sc.subset$n[i])]], prob.vec=prob.vec)
  }
  
  store.all.thetas <- lapply(mat.list, function(x) {sort(unique(c(x))[unique(c(x)) <= 1])})
  
  
  ##### To cut down on computation, try cutting down the number of thetas used:
  ##### max.combns:=max. number of (theta0, theta1) combinations.
  ##### n.thetas*(n.thetas-1)/2 = n.combns, so if n.thetas > sqrt(2*max.combns), take out every other value, excluding 0 and 1.
  ##### Note: further below, more combns are removed if constraints on maxtheta0 and mintheta1 are specified.
  # check ####
  if(max.combns!=Inf){
    maxthetas <- sqrt(2*max.combns)
    for(i in 1:nrow(sc.subset))
    {
      while(length(store.all.thetas[[i]]) > maxthetas)
      {
        every.other.element <- rep(c(FALSE, TRUE), 0.5*(length(store.all.thetas[[i]])-2))
        store.all.thetas[[i]] <-  store.all.thetas[[i]][c(TRUE, every.other.element, TRUE)]
      }
    }
  }

  if(!is.null(exact.theta0) & !is.null(exact.theta1)){ # if exact thetas are given (to speed up a result check):
    for(i in 1:length(store.all.thetas)){
      keep <- abs(store.all.thetas[[i]]-exact.theta0)<1e-3 | abs(store.all.thetas[[i]]-exact.theta1)<1e-3
      store.all.thetas[[i]] <- store.all.thetas[[i]][keep]
    }
  }
  
  h.results.list <- vector("list", nrow(sc.subset)) # 
  
  pb <- txtProgressBar(min = 0, max = nrow(sc.subset), style = 3)
  
  # Now, find the designs, looping over each possible {r, N} combination, and within each {r, N} combination, loop over all combns of {theta0, theta1}:
  for(h in 1:nrow(sc.subset)){
    k <- 1
    blank.mat <- matrix(NA, nrow=nrow(mat.list[[h]]), ncol=ncol(mat.list[[h]]))
    rownames(blank.mat) <- 0:(nrow(blank.mat)-1)
    zero.mat <- matrix(0, nrow=nrow(mat.list[[h]]), ncol=ncol(mat.list[[h]]))
    rownames(zero.mat) <- rownames(blank.mat)
    pat.cols.single <- pat.cols.list[[paste(sc.subset$n[h])]]
    
    ########### START 2D BISECTION
    theta0.vec <- store.all.thetas[[h]][store.all.thetas[[h]]<=maxtheta0]
    theta1.vec <- store.all.thetas[[h]][store.all.thetas[[h]]>=mintheta1]
    h.results <- vector("list", length(theta0.vec)*length(theta1.vec)) 
    # Bounds for theta0:
    a0 <- 1
    b0 <- length(theta0.vec)
    d0 <- ceiling((b0-a0)/2)
    # Bounds for theta1:
    a1 <- 1
    b1 <- length(theta1.vec)
    d1 <- ceiling((b1-a1)/2)
    mintheta0 <- NA
    maxtheta1 <- NA
    while(min((b0-a0),(b1-a1))>1 & is.na(mintheta0)){ # Break/move on when bisection method fails to find anything OR when final feasible design is found.
      output <- find2armBlockOCs(n=sc.subset$n[h], r=sc.subset$r[h], Bsize=Bsize, theta0=theta0.vec[d0], theta1=theta1.vec[d1], mat=mat.list[[h]],  power=power, alpha=alpha,
                                 pat.cols=pat.cols.single, prob.vec=prob.vec, prob.vec.p0=prob.vec.p0, blank.mat=blank.mat, zero.mat=zero.mat)
      if(!is.na(output[6])){ # If ESS is not NA, then design IS feasible, and do:
        feasible <- TRUE
        maxtheta1 <- theta1.vec[d1]
        while((feasible==TRUE) & d0<length(theta0.vec)){
          d0 <- d0+1
          h.results[[k]] <- find2armBlockOCs(n=sc.subset$n[h], r=sc.subset$r[h], Bsize=Bsize, theta0=theta0.vec[d0], theta1=maxtheta1, mat=mat.list[[h]],  power=power, alpha=alpha,
                                             pat.cols=pat.cols.single, prob.vec=prob.vec, prob.vec.p0=prob.vec.p0, blank.mat=blank.mat, zero.mat=zero.mat)
          feasible <- !is.na(h.results[[k]][6])
          k <- k+1
        } # Once the final feasible design for the given theta0/1 is found (or we reach the largest theta0), record theta0 and make it a limit:
        mintheta0 <- theta0.vec[d0-1]
      } else { # If design isn't feasible, decrease theta0, increase theta1 and test again: 
        b0 <- d0
        a1 <- d1
        d0 <- a0 + floor((b0-a0)/2)
        d1 <- a1 + floor((b1-a1)/2)
      }
    }
    if(!is.na(mintheta0)){ # If at least one feasible design was found, then mintheta0 exists, and we search over all theta0/1 combinations subject to the new limits we have just created:
      theta0.vec <- theta0.vec[theta0.vec>=mintheta0]
      theta1.vec <- theta1.vec[theta1.vec<=maxtheta1]
      for(i in 1:length(theta1.vec)){
        for(j in 1:length(theta0.vec)){
          #  print(paste(theta0.vec[i], theta1.vec[j]))
          h.results[[k]] <- find2armBlockOCs(n=sc.subset$n[h], r=sc.subset$r[h], Bsize=Bsize, theta0=theta0.vec[j], theta1=theta1.vec[i], mat=mat.list[[h]],
                                             power=power, alpha=alpha, pat.cols=pat.cols.single, prob.vec=prob.vec, prob.vec.p0=prob.vec.p0, blank.mat=blank.mat, zero.mat=zero.mat)
          k <- k+1
        }
      }
    } # if no feasible designs found, do nothing and let loop end.
    
    setTxtProgressBar(pb, h)
    h.results.df <- do.call(rbind, h.results)
    
    if(!is.null(h.results.df)){
      # Remove all "skipped" results:
      colnames(h.results.df) <- c("n", "r", "block", "alpha", "power", "EssH0", "Ess", "theta0", "theta1", "eff.n")
      h.results.df <- h.results.df[!is.na(h.results.df[, "Ess"]),]
      if(nrow(h.results.df)>0){
        # Remove dominated and duplicated designs:
        discard <- rep(NA, nrow(h.results.df))
        for(i in 1:nrow(h.results.df)){
          discard[i] <- sum(h.results.df[i, "EssH0"] > h.results.df[, "EssH0"] & h.results.df[i, "Ess"] > h.results.df[, "Ess"] & h.results.df[i, "n"] >= h.results.df[, "n"])
        }
        h.results.df <- h.results.df[discard==0,, drop=FALSE]
        
        
        #if(is.matrix(h.results.df)){ # i.e. if there is more than one design (if not, h.results.df is a vector)
        # duplicates <- duplicated(h.results.df[, c("n", "Ess", "EssH0"), drop=FALSE])
        # h.results.df <- h.results.df[!duplicates,, drop=FALSE]
        # }
        h.results.list[[h]] <- h.results.df
      }
    }
  } # End of "h" loop
  
  
  full.results <- do.call(rbind, h.results.list)
  #if(length(full.results)==0) stop("There are no feasible designs for this combination of design parameters" , call. = FALSE)
  if(length(full.results)>0){
    # Discard all "inferior" designs:
    discard <- rep(NA, nrow(full.results))
    for(i in 1:nrow(full.results)){
      discard[i] <- sum(full.results[i, "EssH0"] > full.results[, "EssH0"] & full.results[i, "Ess"] > full.results[, "Ess"] & full.results[i, "n"] >= full.results[, "n"])
      #print(i)
    }
    subset.results <- full.results[discard==0,,drop=FALSE]
    
    
    # Remove duplicates:
    duplicates <- duplicated(subset.results[, c("n", "EssH0", "Ess"), drop=FALSE])
    admissible.ds <- subset.results[!duplicates,,drop=FALSE]
    admissible.ds$looks <- admissible.ds[,"eff.n"]/admissible.ds[,"block"]
    admissible.ds$pc <- rep(pc, nrow(admissible.ds))
    admissible.ds$pt <- rep(pt, nrow(admissible.ds))
    return(admissible.ds)
  }
}
