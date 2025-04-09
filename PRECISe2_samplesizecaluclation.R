###############################################################################|
### PRECISe 2: Sample size calculation 
### September 2024
###
### This script is for the sample size calculation
###
###############################################################################|


mb <- 0.78
m30 <- 0.33
m90 <- 0.38

sdb <- 0.25
sd30 <- 0.33
sd90 <- 0.38

mu    <- c(mb, m30, m90) # c(0.78, 0.33, 0.38)
sd    <- c(sdb, sd30, sd90) # c(0.25, 0.33, 0.38)
es    <- 0.06

#     b      30     90     
# b   0.25^2 0.22   0.23   
# 30  0.22   0.33^2 0.47   
# 90  0.23   0.47   0.38^2 
# Let op! covariance = correlation*sd1*sd2
sig   <- matrix(c(sd[1]^2, 0.22*0.25*0.33, 0.23*0.25*0.38,
                  0.22*0.33*0.25, sd[2]^2, 0.47*0.33*0.38, 
                  0.23*0.38*0.25, 0.47*0.38*0.33, sd[3]^2),
                nrow = 3, byrow = TRUE)
is.positive.definite(sig)


set.seed(7181)
# Simulation

samplesize <- function(n, mu, sig, iters){
  set.seed(7181)
  beta    <- rep(NA, iters)
  pval    <- rep(NA, iters)
  for(i in 1:iters){
    
    Xi <- mvrnorm(n, mu = mu, Sigma = sig, empirical = TRUE)
    Xi[, 1] <- Xi[, 1]
    Xi[, 2] <- Xi[, 2]
    Xi[, 3] <- Xi[, 3]
    
    d  <- data.frame(ID = seq(1:n), Xi)
    dl <- melt(d, id.vars = c("ID", "X1"), measure.vars = c("X2", "X3"),
               variable.name = "fu_moment", value.name = "qol")
    names(dl)[2] <- "baseline"
    dl <- dl[order(dl$ID), ]
    dl$fu_t <- as.numeric(rep(c(30, 90), n))
    dl$rand <- c(rep(0, 0.5*length(dl$ID)), rep(1, 0.5*length(dl$ID)))
    dl$qol  <- ifelse(dl$rand == 1, dl$qol + 0.5*es, dl$qol - 0.5*es)
    
    fit <- lme(qol ~ rand + baseline, data = dl, rand = ~ 1|ID,
               correlation = corAR1(form = ~ fu_t | ID),
               control = list(opt = "optim"))
    beta[i] <- fixef(fit)[2]
    pval[i] <- summary(fit)$tTable[2, 5]
    
  }
  my_list <- list(beta, pval, "power" = (sum(pval <= 0.05)/iters)*100)
  return(my_list) 
}

   

# Results
set.seed(7181)
power <- c()
n <- 466
for (i in n) {
  p <- samplesize(n = 2*i, mu = mu, sig = sig, iters = 1000)$power
  power <- cbind(power, p)
}

# Adjustment for loss to follow-up

N <- (n*2)/(1 - 0.092)
