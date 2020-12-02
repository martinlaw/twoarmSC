####
#### This code can be used to analyse the designs found from running the file
#### "2-reproduce-results-find-designs.R", or alternatively, the designs saved
#### in the file "data.RData". 
####

source("1-reproduce-results-functions-find-designs.R")
source("3-reproduce-results-functions-for-tables-plots.R")
# load("data.RData")

############### Create weights for loss function ###############
 
grid.size <- 0.01
qs <- expand.grid( seq(0, 1, grid.size),  seq(0, 1, grid.size))
qs <- cbind(qs, round(1-qs[, 1]-qs[, 2], 10))

colnames(qs) <- c("w0", "w1", "w2")
hold.wts <- qs[, "w2"] >= 0
qs <- qs[hold.wts, ]
qs <- qs[-1,] # Removed the very first row (all weight on "n") -- because it's possible to get multiple minimums
rownames(qs) <- 1:nrow(qs)

# All loss functions for every design and set of weights:
all.loss <- vector("list", 5)
tf.all.loss <- vector("list", 5)
type.loss <- vector("list", 5)

for(i in 1:5){
# all.loss[[1]] contains all loss scores. It has one row for each admissible design, and one column for each weight combination {w0, w1, w2}
# tf.all.loss[[1]] also has one row for each admissible design, and one column for each weight combination {w0, w1, w2};
# With the exception of ties, each column contains one TRUE value, pertaining to the admissible design with the lowest loss score
  all.loss[[i]] <- apply(qs, 1, function(x) {x["w0"]*all.list[[i]][,"EssH0"] + x["w1"]*all.list[[i]][,"Ess"] + x["w2"]*all.list[[i]][,"n"]})
  tf.all.loss[[i]] <- apply(all.loss[[i]], 2, function(x) x==min(x))
# For each set of weights (i.e. columns), what is the design type of the design with the lowest loss score?
  type.loss[[i]] <- apply(tf.all.loss[[i]], 2, function(x) all.list[[i]]$design[x])
}

# all.loss[[1]] contains all loss scores. It has one row for each admissible design, and one column for each weight combination {w0, w1, w2}
# tf.all.loss[[1]] also has one row for each admissible design, and one column for each weight combination {w0, w1, w2};
# With the exception of ties, each column contains one TRUE value, pertaining to the admissible design with the lowest loss score

table(type.loss[[1]])
table(type.loss[[2]])
table(type.loss[[3]])
table(type.loss[[4]])
table(type.loss[[5]])

######## FIGURE 1 ###########
type.loss.all.df <- data.frame(loss=unlist(type.loss),
                               w0=c(rep(qs[,"w0"], times=5)),
                               w1=c(rep(qs[,"w1"], times=5)),
                               #p0=c(rep(seq(0.1, 0.5, by=0.1), each=5150))
                               p0=rep(c("p0=0.1", "p0=0.2","p0=0.3","p0=0.4","p0=0.5"), each=5150)
                               )

plot1 <- ggplot(type.loss.all.df, aes(w1, w0)) +
  ggtitle('Omni-admissible design type') +
  theme_tufte() +
  xlab('w1') +
  ylab('w0') +
  geom_raster(aes(fill = loss)) +
  scale_fill_discrete(name="",
                      breaks=levels(type.loss.all.df$loss),
                      labels=levels(type.loss.all.df$loss))+
  theme(axis.text.x=element_text(angle=0),
        axis.title=element_text(size=rel(1)),
        axis.text=element_text(size=rel(1)),
       # legend.title = element_blank(),
       legend.position = c(0.75, 0.15),
        plot.title = element_text(hjust = 0.5, size = rel(1.4)))+
  facet_wrap(p0 ~ ., ncol=2)
plot1
# dev.print(pdf, "omni-admissible_design.pdf")


################### Difference in loss score between each design ####################

# Want to compare the best design of each type at each weight combination.
# 5 design types == 10 comparisons. Consider a heatmap, or boxplots.
# To begin, we want a matrix with 8 columns: w0, w1, 1-w0-w1, and the loss score for each of the five design types:


loss.diff <- vector("list", 5)
for(i in 1:5){
  loss.diff[[i]] <- findDifferences(all.loss[[i]], all.list[[i]])
}



loss.diff.block2.carsten <- NULL
loss.diff.block2.chen <- NULL
loss.diff.block2.jung <- NULL
loss.diff.block8.carsten <- NULL
loss.diff.block8.chen <- NULL
loss.diff.block8.jung <- NULL
for(i in 1:5){
  loss.diff.block2.jung <- c(loss.diff.block2.jung, loss.diff[[i]]$block2_jung)
  loss.diff.block2.carsten <- c(loss.diff.block2.carsten, loss.diff[[i]]$block2_carsten)
  loss.diff.block2.chen <- c(loss.diff.block2.chen, loss.diff[[i]]$block2_chen)
  loss.diff.block2.block8 <- c(loss.diff.block2.chen, loss.diff[[i]]$block2_chen)
  
  
  loss.diff.block8.jung <- c(loss.diff.block8.jung, loss.diff[[i]]$block8_jung)
  loss.diff.block8.carsten <- c(loss.diff.block8.carsten, loss.diff[[i]]$block8_carsten)
  loss.diff.block8.chen <- c(loss.diff.block8.chen, loss.diff[[i]]$block8_chen)
}

################ Put the above all together into one data frame. It will all look better using facets ################

loss.diff.full.df <- data.frame(diff=c(loss.diff.block2.jung, loss.diff.block2.carsten, loss.diff.block2.chen),
                                w0=c(rep(qs[,"w0"], times=5), rep(qs[,"w0"], times=5), rep(qs[,"w0"], times=5)),
                                w1=c(rep(qs[,"w1"], times=5), rep(qs[,"w1"], times=5), rep(qs[,"w1"], times=5)),
                                p0=c(rep(seq(0.1, 0.5, by=0.1), each=5150), rep(seq(0.1, 0.5, by=0.1), each=5150), rep(seq(0.1, 0.5, by=0.1), each=5150)),
                                diffdesign=rep(c("Jung", "Carsten", "Chen"), each=25750)
                                )

loss.diff.full.block8.df <- data.frame(diff=c(loss.diff.block8.jung, loss.diff.block8.carsten, loss.diff.block8.chen),
                                w0=c(rep(qs[,"w0"], times=5), rep(qs[,"w0"], times=5), rep(qs[,"w0"], times=5)),
                                w1=c(rep(qs[,"w1"], times=5), rep(qs[,"w1"], times=5), rep(qs[,"w1"], times=5)),
                                p0=c(rep(seq(0.1, 0.5, by=0.1), each=5150), rep(seq(0.1, 0.5, by=0.1), each=5150), rep(seq(0.1, 0.5, by=0.1), each=5150)),
                                diffdesign=rep(c("Jung", "Carsten", "Chen"), each=25750)
)

########## FIGURE 2 ############
allplot <- ggplot(loss.diff.full.df, aes(w1, w0, p0)) +
  theme_tufte() +
  labs(title = "Difference in loss scores: \n  Our approach with block size 2 compared to others, for p0=0.1,...,0.5") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('q1') +
  ylab('q0') +
  geom_tile(aes(fill = diff)) +
  scale_fill_gradient2(midpoint=0, low = 'darkred', mid="white", high = 'darkblue') +
  facet_grid(p0 ~ diffdesign)
allplot
#dev.print(pdf, "loss_score_difference_block2.pdf")

########### SUPPLEMENTARY FIGURE ##########
allplot.block8 <- ggplot(loss.diff.full.block8.df, aes(w1, w0, p0)) +
  theme_tufte() +
  labs(title = "Difference in loss scores: \n  Our approach with block size 8 compared to others, for p0=0.1,...,0.5") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab('q1') +
  ylab('q0') +
  geom_tile(aes(fill = diff)) +
  scale_fill_gradient2(midpoint=0, low = 'darkred', mid="white", high = 'darkblue') +
  facet_grid(p0 ~ diffdesign)
allplot.block8
#dev.print(pdf, "loss_score_difference_block8.pdf")



#################### SECTION 3.2: MISSPECIFICATION OF RESPONSE RATES -- FIGS 3 AND 4 ################
# Take the H0-optimal design for Carsten and block2:
h0.opt.row1 <- do.call(rbind, by(all.list[[1]], all.list[[1]]$design, function(x) x[which.min(x$EssH0), ] ))
h0.opt.row1.carsten <- h0.opt.row1[h0.opt.row1$design=="carsten",]
tt <- varypPowerSS(n=list(7,16), r=list(1,3), pc=0.1, pt=0.3, theta0=NULL, theta1=NULL, method="carsten", lower.triangle = FALSE,
                   p.range=c(0.1, 0.7), runs=10000, seed=1, power.ss="power", title.text=FALSE, plot=F)
tt2 <- varypPowerSS(n=31, r=3, pc=0.1, pt=0.3, theta0=0.1277766, theta1=0.9316456, method="block", lower.triangle = FALSE,
                   p.range=c(0.1, 0.7), runs=10000, seed=1, power.ss="power", title.text=FALSE, plot=F)
 
tt$design <- "Carsten"
tt2$design <- "Block size 2"
ttboth <- rbind(tt, tt2)

plotPowerBoth(ttboth)
#dev.print(pdf, "h0-optimal_p01.pdf")


all.list[[2]]
# Take the H0-optimal design for Carsten and block2:
h0.opt.row2 <- do.call(rbind, by(all.list[[2]], all.list[[2]]$design, function(x) x[which.min(x$EssH0), ] ))
h0.opt.row2.carsten <- h0.opt.row2[h0.opt.row2$design=="carsten",]
h0.carsten.row2 <- varypPowerSS(n=list(10,29), r=list(2,7), pc=0.2, pt=0.4, theta0=NULL, theta1=NULL, method="carsten", lower.triangle = FALSE,
                   p.range=c(0.1, 0.7), runs=1000, seed=2, power.ss="power", title.text=FALSE, plot=F)
h0.block.row2 <- varypPowerSS(n=48, r=5, pc=0.2, pt=0.4, theta0=0.1150793, theta1=0.9635876, method="block", lower.triangle = FALSE,
                    p.range=c(0.1, 0.7), runs=1000, seed=2, power.ss="power", title.text=FALSE, plot=F)
h0.carsten.row2$design <- "Carsten"
h0.block.row2$design <- "Block 2"
h0.row <- rbind(h0.carsten.row2, h0.block.row2)
plotPowerBoth(h0.row)
#dev.print(pdf, "h0-optimal_p02.pdf")







################# SECTION 3.1 -- TABLE 2 #####################
# Create table for H0/1-optimal, H0/1-minimax, for p0=0.3

# Take the H0-optimal design for Carsten and block2:
df3 <- all.list[[3]]
h0.opt.row3 <- do.call(rbind, by(df3, df3$design, function(x) x[which.min(x$EssH0), ] ))
h0.opt.row3 <- h0.opt.row3[order(factor(h0.opt.row3$design, levels = c("jung", "carsten","chen", "block2", "block8"))),]
rownames(h0.opt.row3) <- c("Jung", "Carsten", "Chen", "Block 2", "Block 8")

h1.opt.row3 <- do.call(rbind, by(df3, df3$design, function(x) x[which.min(x$Ess), ] ))
h1.opt.row3 <- h1.opt.row3[order(factor(h1.opt.row3$design, levels = c("jung", "carsten","chen", "block2", "block8"))),]
rownames(h1.opt.row3) <- c("Jung", "Carsten", "Chen", "Block 2", "Block 8")

h0.mini.row3 <- do.call(rbind, by(df3, df3$design, function(x) x[which.min(x$narm), ] ))
h0.mini.row3 <- h0.mini.row3[order(factor(h0.mini.row3$design, levels = c("jung", "carsten","chen", "block2", "block8"))),]
rownames(h0.mini.row3) <- c("Jung", "Carsten", "Chen", "Block 2", "Block 8")

p3 <- rbind(h0.opt.row3, h1.opt.row3, h0.mini.row3)
p3 <- p3[, c("r1", "n1", "r2", "narm", "n", "EssH0", "Ess", "theta0", "theta1")]
head(p3)

colnames(p3) <- c("$r_1$", "$n_1$", "$r$", "$N_{arm}$", "N", "$E(N)_{H_0}$", "$E(N)_{H_1}$", "$\\theta_F$", "$\\theta_E$")
library(xtable)

addtorow <- list(pos=list(0, 5, 10),
                 command=c("$\\mathbf{H_0}$\\textbf{-optimal} \\\\", "\\midrule $\\mathbf{H_1}$\\textbf{-optimal} \\\\",
                           "\\midrule $\\mathbf{H_{0/1}}$\\textbf{-minimax} \\\\"))
print(xtable(p3, caption="$H_0$-optimal, $H_1$-optimal and $H_{0/1}$-minimax designs, for $p_0=0.3, p_1=0.5, \\alpha=0.15, \\text{ power}=0.8$.",
             align=c("r", rep("c", 9)), digits=c(0, 0, 0, 0, 0, 0, 1, 1, 4, 4), label="tab:optim_p03"),
      booktabs=TRUE,
      file="optim_p03.tex",
      sanitize.text.function=function(x){x},
      add.to.row = addtorow,
      hline.after = c(-1, 15))

################# COMPARISON TO GROUP SEQUENTIAL DESIGN USING RPACT ################# 
library(rpact)
design <- getDesignGroupSequential(sided = 1, alpha = 0.15, beta = 0.2,
                                   informationRates = 1:10 / 10,
                                   typeOfDesign = "asOF",
                                   typeBetaSpending = "bsOF")
ss <- getSampleSizeRates(design, pi2 = 0.3, pi1 = 0.5)
# ESS0=59.1, ESS1=62.8, max N=107.2. Compare to Block 8 in Table 2.
################# END OF SECTION 3.1 -- TABLE 2 #####################


############## FIGURES 3 AND 4 ##########
# Take the H0-optimal design for Carsten and block2:
h0.opt.row1 <- do.call(rbind, by(all.list[[1]], all.list[[1]]$design, function(x) x[which.min(x$EssH0), ] ))
h0.opt.row1.carsten <- h0.opt.row1[h0.opt.row1$design=="carsten",]
tt <- varypPowerSS(n=list(7,16), r=list(1,3), pc=0.1, pt=0.3, theta0=NULL, theta1=NULL, method="carsten", lower.triangle = FALSE,
                   p.range=c(0.1, 0.7), runs=10000, seed=1, power.ss="power", title.text=FALSE, plot=F)
tt2 <- varypPowerSS(n=31, r=3, pc=0.1, pt=0.3, theta0=0.1277766, theta1=0.9316456, method="block", lower.triangle = FALSE,
                    p.range=c(0.1, 0.7), runs=10000, seed=1, power.ss="power", title.text=FALSE, plot=F)


tt$design <- "Carsten"
tt2$design <- "Block size 2"
ttboth <- rbind(tt, tt2)

plotPowerBoth(ttboth)
#dev.print(pdf, "h0-optimal_p01.pdf")


# Take the H0-optimal design for Carsten and block2:
h0.opt.row2 <- do.call(rbind, by(all.list[[2]], all.list[[2]]$design, function(x) x[which.min(x$EssH0), ] ))
h0.opt.row2.carsten <- h0.opt.row2[h0.opt.row2$design=="carsten",]
h0.carsten.row2 <- varypPowerSS(n=list(10,29), r=list(2,7), pc=0.2, pt=0.4, theta0=NULL, theta1=NULL, method="carsten", lower.triangle = FALSE,
                                p.range=c(0.1, 0.7), runs=1000, seed=2, power.ss="power", title.text=FALSE, plot=F)
h0.block.row2 <- varypPowerSS(n=48, r=5, pc=0.2, pt=0.4, theta0=0.1150793, theta1=0.9635876, method="block", lower.triangle = FALSE,
                              p.range=c(0.1, 0.7), runs=1000, seed=2, power.ss="power", title.text=FALSE, plot=F)
h0.carsten.row2$design <- "Carsten"
h0.block.row2$design <- "Block 2"
h0.row <- rbind(h0.carsten.row2, h0.block.row2)
plotPowerBoth(h0.row)
#dev.print(pdf, "h0-optimal_p02.pdf")
########## END OF FIGURES 3 AND 4 ###########






######## Section 3.3 ##########
ex1.sub <- ex1.df[, c("r1", "n1", "r", "n", "n.total", "EssH0", "Ess", "theta0", "theta1")]
  
colnames(ex1.sub) <- c("$r_1$", "$n_1$", "$r$", "$N_{arm}$", "N", "$E(N)_{H_0}$", "$E(N)_{H_1}$", "$\\theta_F$", "$\\theta_E$")

addtorow <- list(pos=list(0, 4, 8),
                 command=c("$\\mathbf{H_0}$\\textbf{-optimal} \\\\", "\\midrule $\\mathbf{H_1}$\\textbf{-optimal} \\\\",
                           "\\midrule $\\mathbf{H_{0/1}}$\\textbf{-minimax} \\\\"))
print(xtable(ex1.sub, caption="$H_0$-optimal, $H_1$-optimal and $H_{0/1}$-minimax designs, for $p_0=0.7, p_1=0.85, \\alpha=0.15, \\text{ power}=0.8$.",
             align=c("r", rep("c", 9)), digits=c(0, 0, 0, 0, 0, 0, 1, 1, 4, 4), label="tab:ex1"),
      booktabs=TRUE,
      file="ex1.tex",
      sanitize.text.function=function(x){x},
      add.to.row = addtorow,
      hline.after = c(-1, 12),
      include.rownames = TRUE)
################# END OF SECTION 3.3 -- TABLE 3 #####################

##### Comparison to group sequential design #####
design <- getDesignGroupSequential(sided = 1, alpha = 0.15, beta = 0.2,
                                   informationRates = 1:10 / 10,
                                   typeOfDesign = "asOF",
                                   typeBetaSpending = "bsOF")
# asOF : alpha spending O'Brien & Fleming
# bsOF : beta spending O'Brien & Fleming
ss <- getSampleSizeRates(design, pi2 = 0.3, pi1 = 0.5)



### Comparison of sample size for Simon's design, single-stage single-arm design and single-stage two-arm dseign:
library(clinfun)

# Table 2:
ph2simon(pu = 0.3, pa = 0.5, ep1 = 0.15, ep2 = 0.2)
# H0-optimal: E(N)=17
power.prop.test(p1 = 0.3, p2 = 0.5, sig.level = 0.15, power = .80, alternative = "two.sided")
# Two-arm TWO-SIDED: E(N)=N=124
power.prop.test(p1 = 0.3, p2 = 0.5, sig.level = 0.15, power = .80, alternative = "one.sided")
# Two-arm ONE-SIDED: E(N)=N=84
mean(c(17, 84))
# 50.5
mean(c(25.1, 84))
# 54.55
mean(c(17.6, 84))
# 50.8
mean(c(20.6, 84))
# 52.3

# Single arm:
ph2single(pu = 0.3, pa = 0.5, ep1 = 0.15, ep2 = 0.2) # n=21
mean(c(21, 84))
# 52.5

source("singlearmSC/findSimonN1N2R1R2.R")
source("singlearmSC/findSimonDesigns.R")

ph2simon(pu = 0.7, pa = 0.85, ep1 = 0.15, ep2 = 0.2)
# H0-optimal: E(N)=20.7
power.prop.test(p1 = 0.7, p2 = 0.85, sig.level = 0.15, power = .80, alternative = "two.sided")
# Two-arm TWO-SIDED: E(N)=N=160
power.prop.test(p1 = 0.7, p2 = 0.85, sig.level = 0.15, power = .80, alternative = "one.sided")
# Two-arm ONE-SIDED: E(N)=N=108
mean(c(21, 108))
# 64.5
mean(c(26.5, 108))
# 67.25
mean(c(28.4, 108))
# 68.2

# Single arm:
ph2single(pu = 0.7, pa = 0.85, ep1 = 0.15, ep2 = 0.2) # n=31
mean(c(31, 108))
# 69.5



# Find minimum sample size for obtained designs: #####
source("returnCPmat.R")
#load("H:/PhD/twoarm/tidy_block_results.RData")
#load("tidy_block_results.RData")
findMinN <- function(df, opt, p0, p1){
  if(opt=="h0opt"){    des <- df[which.min(df$EssH0),]  }
  if(opt=="h1opt"){    des <- df[which.min(df$Ess),]  }
  if(opt=="h0mini"){   des <- df[order(df$n, df$EssH0),][1,]  }
  if(opt=="h1mini"){   des <- df[order(df$n, df$Ess),][1,]  }
  cpmat <- returnCPmat(n = des$n, r = des$r, Bsize = des$block, pc = p0, pt = p1, theta0 = des$theta0, theta1 = des$theta1, alpha = 0.15, power = 0.8) 
  first.row <- cpmat[1,]
  minN <- which.min(first.row)
  minN
}

p0 <- c(0.1, 0.2, 0.3, 0.4, 0.5)
p1 <- p0+0.2
minN.block2 <- matrix(NA, nrow=5, ncol=4)
colnames(minN.block2) <- c("h0opt", "h1opt", "h0mini", "h1mini")
for(i in 1:5){
  minN.block2[i,] <- c(findMinN(block2.1e6[[i]], "h0opt", p0=p0[i], p1=p1[i]),
                       findMinN(block2.1e6[[i]], "h1opt", p0=p0[i], p1=p1[i]),
                       findMinN(block2.1e6[[i]], "h0mini", p0=p0[i], p1=p1[i]),
                       findMinN(block2.1e6[[i]], "h1mini", p0=p0[i], p1=p1[i])
                       )
}
median(minN.block2)
quantile(minN.block2, 0.25)
quantile(minN.block2, 0.75)

minN.block8 <- matrix(NA, nrow=5, ncol=4)
colnames(minN.block8) <- c("h0opt", "h1opt", "h0mini", "h1mini")
for(i in 1:5){
  minN.block8[i,] <- c(findMinN(block8.1e6[[i]], "h0opt", p0=p0[i], p1=p1[i]),
                       findMinN(block8.1e6[[i]], "h1opt", p0=p0[i], p1=p1[i]),
                       findMinN(block8.1e6[[i]], "h0mini",p0=p0[i], p1=p1[i]),
                       findMinN(block8.1e6[[i]], "h1mini",p0=p0[i], p1=p1[i])
                       )
}
median(minN.block8)
quantile(minN.block8, 0.25)
quantile(minN.block8, 0.75)



# Plot thetaF vs thetaE for all feasible designs: ####
system.time({
  n40 <- findSCdes(nmin=20,
                   nmax=20,
                   pc=0.1,
                   pt=0.4,
                   alpha=0.15,
                   power=0.8,
                   maxtheta0=0.4,
                   mintheta1=0.7,
                   bounds=NULL,
                   fixed.r=3,
                   block.size = 2,
                   max.combns=1e6,
                   rm.dominated.designs = FALSE)
}) # 1 min for 1e6
n40.best <- rmDominatedDesigns(df=n40)

system.time({
  n60 <- findSCdes(nmin=30,
                   nmax=30,
                   pc=0.1,
                   pt=0.4,
                   alpha=0.15,
                   power=0.8,
                   maxtheta0=0.4,
                   mintheta1=0.7,
                   bounds=NULL,
                   fixed.r=5,
                   block.size = 2,
                   max.combns=1e6,
                   rm.dominated.designs = FALSE)
}) # 30secs for 1e4, 2mins for 1e5, 6mins/500secs for 1e6
n60.best <- rmDominatedDesigns(df=n60)

# pdf(file="figs/thetaF_vs_thetaE_ESS0_n40_2.pdf", width = 8, height = 7)
# plotFeasible(n40, criterion = 0)
# dev.off()
# pdf(file="figs/thetaF_vs_thetaE_ESS1_n40_2.pdf", width = 8, height = 7)
# plotFeasible(n40, criterion = 1)
# dev.off()
# pdf(file="figs/thetaF_vs_thetaE_ESS0_n60_2.pdf", width = 8, height = 7)
# plotFeasible(n60, criterion = 0)
# dev.off()
# pdf(file="figs/thetaF_vs_thetaE_ESS1_n60_2.pdf", width = 8, height = 7)
# plotFeasible(n60, criterion = 1)
# dev.off()


# Plot rejection regions for Carsten and for proposed design: ####
system.time({
  n40.fast <- find2armBlockDesigns(nmin=20,
                                   nmax=20,
                                   pc=0.1,
                                   pt=0.4,
                                   alpha=0.15,
                                   power=0.8,
                                   maxtheta0=0.4,
                                   mintheta1=0.7,
                                   bounds=NULL,
                                   fixed.r=3,
                                   block.size = 2,
                                   max.combns=1e6)
}) # 6 secs for 1e6

system.time({
  n60.fast <- find2armBlockDesigns(nmin=30,
                                   nmax=30,
                                   pc=0.1,
                                   pt=0.4,
                                   alpha=0.15,
                                   power=0.8,
                                   maxtheta0=0.4,
                                   mintheta1=0.7,
                                   bounds=NULL,
                                   fixed.r=5,
                                   block.size = 2,
                                   max.combns=1e6)
}) # 10 secs for 1e5, 35 secs for 1e6

n40.best
n40.fast
n60.best
n60.fast

rej.region.carsten <- findRejectionRegions(n=list(20, 34), r=list(4, 10), pc=0.3, pt=0.5, method="carsten")
rej.region.block2 <- findRejectionRegions(n=40, r=4, pc=0.3, pt=0.5, method="block", theta0 = 0.04276781, theta1 = 0.9841904)

# pdf("figs/rejection_region_carsten.pdf", width = 7, height = 7)
# plotRejectionRegions(rej.region.carsten, "carsten")
# dev.off()
# pdf("figs/rejection_region_block3.pdf", width = 7, height = 7)
# plotRejectionRegions(rej.region.block2, "block")
# dev.off()
