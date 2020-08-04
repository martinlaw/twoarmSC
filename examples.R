source("1-reproduce-results-functions-find-designs.R")

####
#
# function: findSCdes
#
# findSCdes is the function used to find designs. Arguments:
#
# nmin, nmax:   min and max sample size per arm to search over.
# pc, pt:       anticipated response rates on control and treatment arms.
# alpha, power: required levels for type I error and complement of type II error.
# maxtheta0:    maximum permitted lower CP threshold (for making a no go decision),
#               denoted theta_F_max in manuscript.
# mintheta1:    minimum permitted upper CP threshold  (for making a go decision), 
#               denoted theta_E_min in manuscript.
# bounds:       restricts the possible values of final stopping boundary r based
#               on the work of A'Hern ("ahern"), ("wald") or neither (NULL).
#               Using either "ahern" or "wald" will considerably speed up the code
#               with little (if any) loss of admissible designs. Further, although
#               "ahern" is more restrictive (and thus faster) than "wald", this also
#               seems to have no negative effect on the number of admissible designs found.
# fixed.r:      restricts the possible values of final stopping boundary r to any value(s) 
#               chosen. Mostly useful to rerun a particular analysis.
# max.combns:   maximum number of CP combinations searched over per (r,N). Reducing
#               this value greatly improves run time.
#
# The output is a dataframe of admissible designs, one per row. 
# Each design has the following output:
#
# n: max sample size per arm
# r: difference in number of responses at final analysis required to reject H0
# block: Block size -- the total number of patients between interim analyses
# alpha, power: type I error and complement of type II error
# EssH0: ESS(p0, p0) as described in manuscript
# Ess: ESS(p0, p1) as described in manuscript
# theta0: thetaF as described in manuscript
# theta1: thetaE as described in manuscript
# eff.n: max sample size total
# looks: max number of interim analyses
# pc, pt: anticipated response rates on control and treatment arms respectively
#
####

system.time({
des <- findSCdes(nmin=40,
          nmax=44,
          block.size=8,
          pc=0.3,
          pt=0.5,
          alpha=0.15,
          power=0.8,
          maxtheta0=0.5,
          mintheta1=0.9,
          bounds="ahern",
          fixed.r=NULL,
          max.combns=1e5)
})


####
#
# function: findBounds
# 
# findBounds is used to obtain the stopping boundaries of a design found using findSCdes.
# It takes a single argument: a single design, i.e. single row from the output of findSCdes.
# The output is a dataframe containing the lower and upper stopping boundaries in terms of
# number of successes so far, i.e. responses on treatment arm plus non-responses on control arm.
# 
# Note: "NA" means that there is no stopping boundary at that point. If there is no lower or
# upper stopping boundary, no interim analysis need take place.
#
####

boundaries <- findBounds(des[1,])





# Parallelisation is strongly recommended. Below is an example of parallelised code:

block4.n <- seq(from=16, to=48, by=2)

library(doParallel)
cores <- detectCores()-2
registerDoParallel(cores)

results.par <- foreach(i=block4.n, .combine = rbind) %dopar% {
  findSCdes(nmin=i,
            nmax=i,
            block.size=4,
            pc=0.3,
            pt=0.5,
            alpha=0.15,
            power=0.8,
            maxtheta0=0.5,
            mintheta1=0.9,
            bounds="ahern",
            max.combns=1e5)
} # Approx. 15 minutes using 20 cores.