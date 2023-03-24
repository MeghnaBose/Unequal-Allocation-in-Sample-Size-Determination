library(dplyr)
library(pbapply)
library(resample)
library(data.table)

p_A <- 0.1
VE <- 0.4
p_B <- (1-VE)*p_A
p <- c(p_A,p_B)
W <- 0.5

sample_size <- function(p_A, VE, type=c("Equal", "Double", "Neyman", "RSIHR"), W, alpha = 0.05)
{
  p_B <- (1-VE)*p_A
  d <- log((W/(1-VE) + sqrt((W/(1-VE))^2+4))/2)
  rho <- case_when(type == "Equal" ~ 0.5,
                   type == "Double" ~ 2/3,
                   type == "Neyman" ~ sqrt((1-p_A)/p_A)/(sqrt((1-p_A)/p_A)+sqrt((1-p_B)/p_B)),
                   type == "RSIHR" ~ (sqrt(1-p_A)*p_B)/(sqrt(1-p_A)*p_B + sqrt(1-p_B)*p_A))
  z <- qnorm(alpha/2, lower.tail = FALSE)
  n <- (z/d)^2*((1-p_A)/(p_A*rho)+(1-p_B)/(p_B*(1-rho)))
  return(c(round(n*rho),round(n*(1-rho))))
}

n <- sample_size(p_A, VE, type =c("Equal", "Double", "Neyman", "RSIHR"), W);n

n_Neyman <- sample_size(p_A, VE, type = "RSIHR", W);n_Neyman

power <- function(p, n, alpha = 0.05)
{
  beta <- log(p[2]/p[1])
  sigma <- (1-p[1])/(p[1]*n[1])+(1-p[2])/(p[2]*n[2])
  p_0 <- (n[1]*p[1]+n[2]*p[2])/(n[1]+n[2])
  sigma_0 <- sqrt((1-p_0)/p_0*(1/n[1]+1/n[2]))
  z <- qnorm(alpha, lower.tail = FALSE)
  power <- 1 - pnorm((z*sigma_0+beta)/sigma)
  return(power)
}
power(p, n_Neyman)
g <- function(x,y,gamma = 0){
  z <- case_when(x == 1 ~ 0,
                 x == 0 ~ 1, 
                 TRUE ~ (y*(y/x)^gamma)/(y*(y/x)^gamma+(1-y)*((1-y)/(1-x))^gamma))
  return(z)
}

ERADE_g <- function(x,y, delta=0.5)
{
  z <- case_when(x > y ~ delta*y,
                 x == y ~ y,
                 x < y ~ 1-delta*(1-y))
  return(z)
}

rho <- function(p, target = c("Equal", "Neyman", "RSIHR")){
  r <- case_when(target == "Equal" ~ 0.5,
                 target == "Neyman" ~ sqrt((1-p[1])*p[2])/(sqrt((1-p[1])*p[2])+sqrt((1-p[2])*p[1])),
                 target == "RSIHR" ~ (sqrt(1-p[1])*p[2])/(p[1]*sqrt(1-p[2])+sqrt(1-p[1])*p[2]))
  return(r)
} 
# Simulation
sample_size_sim <- function(seed = 5, p, alpha = 0.05, desired_W, target = c("Equal", "Neyman", "RSIHR"), gamma = 0){
  
alloc <- c()
response <- c()
alloc[1:(2*seed)] <- 1
seed_vac <- sample(c(1:(2*seed)),seed, replace = FALSE)
alloc[seed_vac] <- 2
response <- rbinom(length(alloc),1, prob = p[alloc])

#est_p <- c(mean(response[which(alloc == 1)]), mean(response[which(alloc == 2)]))
est_p <- c((sum(response[which(alloc == 1)]))/(length(response[which(alloc == 1)])), (sum(response[which(alloc == 2)]))/(length(response[which(alloc == 2)])))
n <- c(length(response[which(alloc == 1)]), length(response[which(alloc == 2)]))

i = 2*seed
while(all(est_p == 0) || all(est_p == 1))
{
  alloc[i+1] <- sample(c(1,2), 1)
  response[i+1] <- rbinom(1, 1, p[alloc[i+1]])

est_p <- c((sum(response[which(alloc == 1)]))/(length(response[which(alloc == 1)])), (sum(response[which(alloc == 2)]))/(length(response[which(alloc == 2)])))
n <- c(length(response[which(alloc == 1)]), length(response[which(alloc == 2)]))
i = i+1
}

est_VE <- 1 - est_p[2]/est_p[1]
est_sigma <- sqrt((1-est_p[1])/(est_p[1]*n[1])+(1-est_p[2])/(est_p[2]*n[2]))
z <- qnorm(alpha/2, lower.tail = FALSE)
d <- z*est_sigma

W <- (1-est_VE)*(exp(d)-exp(-d))
while(W > desired_W || is.na(W))
{
  alloc[i+1] <- sample(c(1,2), 1, prob = c(g(mean(alloc == 1),rho(est_p, target = target), gamma = gamma), 1 - g(mean(alloc == 1),rho(est_p, target = target), gamma = gamma)))
  response[i+1] <- rbinom(1, 1, p[alloc[i+1]])
  
  est_p <- c((sum(response[which(alloc == 1)]))/(length(response[which(alloc == 1)])), (sum(response[which(alloc == 2)]))/(length(response[which(alloc == 2)])))
  n <- c(length(response[which(alloc == 1)]), length(response[which(alloc == 2)]))
  
  est_VE <- 1 - est_p[2]/est_p[1]
  est_sigma <- sqrt((1-est_p[1])/(est_p[1]*n[1])+(1-est_p[2])/(est_p[2]*n[2]))
  print(est_sigma)
  d <- z*est_sigma
  
  W <- (1-est_VE)*(exp(d)-exp(-d))
  #print(W)
  i = i + 1
}
est_beta <- log(est_p[2]/est_p[1])
est_p_0 <- (n[1]*est_p[1]+n[2]*est_p[2])/(n[1]+n[2])
est_sigma_0 <- sqrt((1-est_p_0)/est_p_0*(1/n[1]+1/n[2]))
Z <- est_beta/est_sigma_0
z_c <- -qnorm(0.05, lower.tail = FALSE)
rej_ind <- if_else(Z < z_c, 1, 0)
ret <- c(sum(n), n[2], rej_ind)
return(ret)
}

#s <- sample_size_sim(p = p, desired_W = 0.5, target = "Equal", gamma = 0);s
#pboptions(type="txt", char = "=")

# ERADE Simulation
sample_size_sim_ERADE <- function(seed = 5, p, alpha = 0.05, desired_W, target = c("Equal", "Neyman", "RSIHR"), delta = 0.5){
  
  alloc <- c()
  response <- c()
  alloc[1:(2*seed)] <- 1
  seed_vac <- sample(c(1:(2*seed)),seed, replace = FALSE)
  alloc[seed_vac] <- 2
  response <- rbinom(length(alloc),1, prob = p[alloc])
  
  #est_p <- c(mean(response[which(alloc == 1)]), mean(response[which(alloc == 2)]))
  est_p <- c((sum(response[which(alloc == 1)]))/(length(response[which(alloc == 1)])), (sum(response[which(alloc == 2)]))/(length(response[which(alloc == 2)])))
  n <- c(length(response[which(alloc == 1)]), length(response[which(alloc == 2)]))
  
  i = 2*seed
  while(all(est_p == 0) || all(est_p == 1))
  {
    alloc[i+1] <- sample(c(1,2), 1)
    response[i+1] <- rbinom(1, 1, p[alloc[i+1]])
    
    est_p <- c((sum(response[which(alloc == 1)]))/(length(response[which(alloc == 1)])), (sum(response[which(alloc == 2)]))/(length(response[which(alloc == 2)])))
    n <- c(length(response[which(alloc == 1)]), length(response[which(alloc == 2)]))
    i = i+1
  }
  
  est_VE <- 1 - est_p[2]/est_p[1]
  est_sigma <- sqrt((1-est_p[1])/(est_p[1]*n[1])+(1-est_p[2])/(est_p[2]*n[2]))
  z <- qnorm(alpha/2, lower.tail = FALSE)
  d <- z*est_sigma
  
  W <- (1-est_VE)*(exp(d)-exp(-d))
  while(W > desired_W || is.na(W))
  {
    #print(W)
    alloc[i+1] <- sample(c(1,2), 1, prob = c(ERADE_g(mean(alloc == 1),rho(est_p, target = target), delta = delta), 1 - ERADE_g(mean(alloc == 1),rho(est_p, target = target), delta = delta)))
    response[i+1] <- rbinom(1, 1, p[alloc[i+1]])
    
    est_p <- c((sum(response[which(alloc == 1)]))/(length(response[which(alloc == 1)])), (sum(response[which(alloc == 2)]))/(length(response[which(alloc == 2)])))
    n <- c(length(response[which(alloc == 1)]), length(response[which(alloc == 2)]))
    
    est_VE <- 1 - est_p[2]/est_p[1]
    est_sigma <- sqrt((1-est_p[1])/(est_p[1]*n[1])+(1-est_p[2])/(est_p[2]*n[2]))
    d <- z*est_sigma
    
    W <- (1-est_VE)*(exp(d)-exp(-d))
    #print(W)
    i = i + 1
  }
  est_beta <- log(est_p[2]/est_p[1])
  est_p_0 <- (n[1]*est_p[1]+n[2]*est_p[2])/(n[1]+n[2])
  est_sigma_0 <- sqrt((1-est_p_0)/est_p_0*(1/n[1]+1/n[2]))
  Z <- est_beta/est_sigma_0
  z_c <- -qnorm(0.05, lower.tail = FALSE)
  rej_ind <- if_else(Z < z_c, 1, 0)
  ret <- c(sum(n), n[2], rej_ind)
  return(ret)
}

sample_size_sim_Double <- function(seed = 5, p, alpha = 0.05, desired_W){
  
  alloc <- c()
  response <- c()
  alloc[1:(2*seed)] <- 1
  seed_vac <- sample(c(1:(2*seed)),seed, replace = FALSE)
  alloc[seed_vac] <- 2
  response <- rbinom(length(alloc),1, prob = p[alloc])
  
  #est_p <- c(mean(response[which(alloc == 1)]), mean(response[which(alloc == 2)]))
  est_p <- c((sum(response[which(alloc == 1)]))/(length(response[which(alloc == 1)])), (sum(response[which(alloc == 2)]))/(length(response[which(alloc == 2)])))
  n <- c(length(response[which(alloc == 1)]), length(response[which(alloc == 2)]))
  
  est_VE <- 1 - est_p[2]/est_p[1]
  est_sigma <- sqrt((1-est_p[1])/(est_p[1]*n[1])+(1-est_p[2])/(est_p[2]*n[2]))
  z <- qnorm(alpha/2, lower.tail = FALSE)
  d <- z*est_sigma
  
  W <- (1-est_VE)*(exp(d)-exp(-d))
  i = 2*seed
  while(W > desired_W || is.na(W))
  {
    #print(W)
    alloc[i+1] <- sample(c(1,2), 1, prob = c(1/3,2/3))
    response[i+1] <- rbinom(1, 1, p[alloc[i+1]])
    
    est_p <- c((sum(response[which(alloc == 1)]))/(length(response[which(alloc == 1)])), (sum(response[which(alloc == 2)]))/(length(response[which(alloc == 2)])))
    n <- c(length(response[which(alloc == 1)]), length(response[which(alloc == 2)]))
    
    est_VE <- 1 - est_p[2]/est_p[1]
    est_sigma <- sqrt((1-est_p[1])/(est_p[1]*n[1])+(1-est_p[2])/(est_p[2]*n[2]))
    z <- qnorm(alpha/2, lower.tail = FALSE)
    d <- z*est_sigma
    
    W <- (1-est_VE)*(exp(d)-exp(-d))
    i = i+1
  }
  
  est_beta <- log(est_p[2]/est_p[1])
  est_p_0 <- (n[1]*est_p[1]+n[2]*est_p[2])/(n[1]+n[2])
  est_sigma_0 <- sqrt((1-est_p_0)/est_p_0*(1/n[1]+1/n[2]))
  Z <- est_beta/est_sigma_0
  z_c <- -qnorm(0.05, lower.tail = FALSE)
  rej_ind <- if_else(Z < z_c, 1, 0)
  ret <- c(sum(n), n[2], rej_ind)
  return(ret)
}

nsim <- 10
#Sample_Sizes_sim <- pbreplicate(nsim, {sample_size_sim(p = p, desired_W = 0.24)})

#SMLE
Sample_Sizes_sim_Equal <- matrix(NA, nrow = nsim, ncol = 3)
Sample_Sizes_sim_Double <- matrix(NA, nrow = nsim, ncol = 3)
Sample_Sizes_sim_Neyman_SMLE <- matrix(NA, nrow = nsim, ncol = 3)
Sample_Sizes_sim_RSIHR_SMLE <- matrix(NA, nrow = nsim, ncol = 3)

for(i in c(1:nsim))
{
  print(i)
  Sample_Sizes_sim_Equal[i,] <- sample_size_sim(p=p,desired_W = W,target = "Equal", gamma = 0)
  #Sample_Sizes_sim_Double[i,] <- sample_size_sim_Double(p = p, desired_W = W)
  #Sample_Sizes_sim_Neyman_SMLE[i,] <- sample_size_sim(p=p,desired_W = W,target = "Neyman", gamma = 0)
  #Sample_Sizes_sim_RSIHR_SMLE[i,] <- sample_size_sim(p=p,desired_W = W,target = "RSIHR", gamma = 0)
}

Table_SMLE <- data.frame(cbind(Sample_Sizes_sim_Equal, Sample_Sizes_sim_Double, Sample_Sizes_sim_Neyman_SMLE, Sample_Sizes_sim_RSIHR_SMLE))

colnames(Table_SMLE) <- c(rep("Equal",3),rep("Double",3),rep("Neyman",3),rep("RSIHR",3))

write.csv(Table_SMLE, "SMLE_Asy_Sim_Results.csv", row.names = FALSE)


# Table_SMLE <- data.frame(rbind(colMeans(Sample_Sizes_sim_Equal), colMeans(Sample_Sizes_sim_Double), colMeans(Sample_Sizes_sim_Neyman_SMLE), colMeans(Sample_Sizes_sim_RSIHR_SMLE)))
# 
# colnames(Table_SMLE) <- c("Total", "Vaccine", "Simulated Power")
# 
# Table_SMLE <- cbind(Allocation_Type = c("Equal", "Double", "Neyman", "RSIHR"), Table_SMLE)

#DBCD
Sample_Sizes_sim_Neyman_DBCD <- matrix(NA, nrow = nsim, ncol = 3)
Sample_Sizes_sim_RSIHR_DBCD <- matrix(NA, nrow = nsim, ncol = 3)    

for(i in c(1:nsim))
{
  print(i)
  Sample_Sizes_sim_Neyman_DBCD[i,] <- sample_size_sim(p=p,desired_W = W,target = "Neyman", gamma = 2)
  Sample_Sizes_sim_RSIHR_DBCD[i,] <- sample_size_sim(p=p,desired_W = W,target = "RSIHR", gamma = 2)
}

Table_DBCD <- data.frame(cbind(Sample_Sizes_sim_Neyman_DBCD, Sample_Sizes_sim_RSIHR_DBCD))

colnames(Table_DBCD) <- c(rep("Neyman",3),rep("RSIHR",3))

write.csv(Table_DBCD, "DBCD_Asy_Sim_Results.csv", row.names = FALSE)

# ERADE
Sample_Sizes_sim_Neyman_ERADE <- matrix(NA, nrow = nsim, ncol = 3)
Sample_Sizes_sim_RSIHR_ERADE <- matrix(NA, nrow = nsim, ncol = 3)

for(i in c(1:nsim))
{
  print(i)
  Sample_Sizes_sim_Neyman_ERADE[i,] <- sample_size_sim_ERADE(p=p,desired_W = W,target = "Neyman")
  Sample_Sizes_sim_RSIHR_ERADE[i,] <- sample_size_sim(p=p,desired_W = W,target = "RSIHR")
}

Table_ERADE <- data.frame(cbind(Sample_Sizes_sim_Neyman_ERADE, Sample_Sizes_sim_RSIHR_ERADE))

colnames(Table_ERADE) <- c(rep("Neyman",3),rep("RSIHR",3))

write.csv(Table_ERADE, "ERADE_Asy_Sim_Results.csv", row.names = FALSE)
#===================================================
# Test simulation
# W <- c()
# for (i in c(1:10000))
# {
# vac_response <- rbinom(8563, 1, p_B)
# unvac_response <- rbinom(3814, 1, p_A)
# est_p <- c(mean(unvac_response), mean(vac_response))
# n <- c(length(response[which(alloc == 1)]), length(response[which(alloc == 2)]))
# est_VE <- 1 - est_p[2]/est_p[1]
# est_sigma <- sqrt((1-est_p[1])/(est_p[1]*n[1])+(1-est_p[2])/(est_p[2]*n[2]))
# z <- qnorm(alpha/2, lower.tail = FALSE)
# d <- z*est_sigma
# 
# W[i] <- (1-est_VE)*(exp(d)-exp(-d))
# }
#===================================================
# Simulation considering drift

sample_size_sim_drift <- function(seed = 5, p, alpha = 0.05, desired_W, drift = 1.001, target = c("Equal", "Neyman", "RSIHR"), gamma = 0){
  
  alloc <- c()
  response <- c()
  alloc[1:(2*seed)] <- 1
  
  seed_vac <- sample(c(1:(2*seed)),seed, replace = FALSE)
  alloc[seed_vac] <- 2
  response <- c()
  for(i in c(1:(2*seed)))
  {
    p <- p*drift
    response[i] <- rbinom(1,1, prob = p[alloc[i]])
  }
  
  
  #est_p <- c(mean(response[which(alloc == 1)]), mean(response[which(alloc == 2)]))
  est_p <- c((sum(response[which(alloc == 1)]))/(length(response[which(alloc == 1)])), (sum(response[which(alloc == 2)]))/(length(response[which(alloc == 2)])))
  n <- c(length(response[which(alloc == 1)]), length(response[which(alloc == 2)]))
  
  i = 2*seed
  while(all(est_p == 0) || all(est_p == 1))
  {
    alloc[i+1] <- sample(c(1,2), 1)
    p <- p*drift
    #print(p)
    response[i+1] <- rbinom(1, 1, p[alloc[i+1]])
    
    est_p <- c((sum(response[which(alloc == 1)]))/(length(response[which(alloc == 1)])), (sum(response[which(alloc == 2)]))/(length(response[which(alloc == 2)])))
    n <- c(length(response[which(alloc == 1)]), length(response[which(alloc == 2)]))
    i = i+1
  }
  
  est_VE <- 1 - est_p[2]/est_p[1]
  est_sigma <- sqrt((1-est_p[1])/(est_p[1]*n[1])+(1-est_p[2])/(est_p[2]*n[2]))
  z <- qnorm(alpha/2, lower.tail = FALSE)
  d <- z*est_sigma
  
  W <- (1-est_VE)*(exp(d)-exp(-d))

  while(W > desired_W || is.na(W))
  {
    p <- p*drift
    #print(p)
    alloc[i+1] <- sample(c(1,2), 1, prob = c(g(mean(alloc == 1),rho(est_p, target = target), gamma = gamma), 1 - g(mean(alloc == 1),rho(est_p, target =  target), gamma = gamma)))
    response[i+1] <- rbinom(1, 1, p[alloc[i+1]])
    
    est_p <- c((sum(response[which(alloc == 1)]))/(length(response[which(alloc == 1)])), (sum(response[which(alloc == 2)]))/(length(response[which(alloc == 2)])))
    n <- c(length(response[which(alloc == 1)]), length(response[which(alloc == 2)]))
    
    est_VE <- 1 - est_p[2]/est_p[1]
    est_sigma <- sqrt((1-est_p[1])/(est_p[1]*n[1])+(1-est_p[2])/(est_p[2]*n[2]))
    d <- z*est_sigma
    
    W <- (1-est_VE)*(exp(d)-exp(-d))
    print(W)
    #print(W)
    #print(W)
    i = i + 1
  }
  est_beta <- log(est_p[2]/est_p[1])
  est_p_0 <- (n[1]*est_p[1]+n[2]*est_p[2])/(n[1]+n[2])
  est_sigma_0 <- sqrt((1-est_p_0)/est_p_0*(1/n[1]+1/n[2]))
  Z <- est_beta/est_sigma_0
  z_c <- -qnorm(0.05, lower.tail = FALSE)
  rej_ind <- if_else(Z < z_c, 1, 0)
  ret <- c(sum(n), n[2],rej_ind)
  return(ret)
}

#ERADE

sample_size_sim_drift_ERADE <- function(seed = 5, p, alpha = 0.05, desired_W, drift = 1.001, target = c("Equal", "Neyman", "RSIHR"), delta = 0.5){
  
  alloc <- c()
  response <- c()
  alloc[1:(2*seed)] <- 1
  
  seed_vac <- sample(c(1:(2*seed)),seed, replace = FALSE)
  alloc[seed_vac] <- 2
  response <- c()
  for(i in c(1:(2*seed)))
  {
    p <- p*drift
    response[i] <- rbinom(1,1, prob = p[alloc[i]])
  }
  
  
  #est_p <- c(mean(response[which(alloc == 1)]), mean(response[which(alloc == 2)]))
  est_p <- c((sum(response[which(alloc == 1)]))/(length(response[which(alloc == 1)])), (sum(response[which(alloc == 2)]))/(length(response[which(alloc == 2)])))
  n <- c(length(response[which(alloc == 1)]), length(response[which(alloc == 2)]))
  
  est_VE <- 1 - est_p[2]/est_p[1]
  est_sigma <- sqrt((1-est_p[1])/(est_p[1]*n[1])+(1-est_p[2])/(est_p[2]*n[2]))
  z <- qnorm(alpha/2, lower.tail = FALSE)
  d <- z*est_sigma
  
  W <- (1-est_VE)*(exp(d)-exp(-d))
  
  i = 2*seed
  while(all(est_p == 0) || all(est_p == 1))
  {
    alloc[i+1] <- sample(c(1,2), 1)
    p <- p*drift
    #print(p)
    response[i+1] <- rbinom(1, 1, p[alloc[i+1]])
    
    est_p <- c((sum(response[which(alloc == 1)]))/(length(response[which(alloc == 1)])), (sum(response[which(alloc == 2)]))/(length(response[which(alloc == 2)])))
    n <- c(length(response[which(alloc == 1)]), length(response[which(alloc == 2)]))
    i = i+1
  }
  
  est_VE <- 1 - est_p[2]/est_p[1]
  est_sigma <- sqrt((1-est_p[1])/(est_p[1]*n[1])+(1-est_p[2])/(est_p[2]*n[2]))
  z <- qnorm(alpha/2, lower.tail = FALSE)
  d <- z*est_sigma
  
  W <- (1-est_VE)*(exp(d)-exp(-d))
  
  while(W > desired_W || is.na(W))
  {
    p <- p*drift
    #print(p)
    alloc[i+1] <- sample(c(1,2), 1, prob = c(ERADE_g(mean(alloc == 1),rho(est_p, target = target), delta = delta), 1 - ERADE_g(mean(alloc == 1),rho(est_p, target =  target), delta = delta)))
    response[i+1] <- rbinom(1, 1, p[alloc[i+1]])
    
    est_p <- c((sum(response[which(alloc == 1)]))/(length(response[which(alloc == 1)])), (sum(response[which(alloc == 2)]))/(length(response[which(alloc == 2)])))
    n <- c(length(response[which(alloc == 1)]), length(response[which(alloc == 2)]))
    
    est_VE <- 1 - est_p[2]/est_p[1]
    est_sigma <- sqrt((1-est_p[1])/(est_p[1]*n[1])+(1-est_p[2])/(est_p[2]*n[2]))
    d <- z*est_sigma
    
    W <- (1-est_VE)*(exp(d)-exp(-d))
    #print(W)
    #print(W)
    i = i + 1
  }
  est_beta <- log(est_p[2]/est_p[1])
  est_p_0 <- (n[1]*est_p[1]+n[2]*est_p[2])/(n[1]+n[2])
  est_sigma_0 <- sqrt((1-est_p_0)/est_p_0*(1/n[1]+1/n[2]))
  Z <- est_beta/est_sigma_0
  z_c <- -qnorm(0.05, lower.tail = FALSE)
  rej_ind <- if_else(Z < z_c, 1, 0)
  ret <- c(sum(n), n[2],rej_ind)
  return(ret)
}

# Double Allocation

sample_size_sim_drift_Double <- function(seed = 5, p, alpha = 0.05, desired_W, drift = 1.001){
  
  alloc <- c()
  response <- c()
  alloc[1:(2*seed)] <- 1
  
  seed_vac <- sample(c(1:(2*seed)),seed, replace = FALSE)
  alloc[seed_vac] <- 2
  response <- c()
  for(i in c(1:(2*seed)))
  {
    p <- p*drift
    response[i] <- rbinom(1,1, prob = p[alloc[i]])
  }
  
  
  #est_p <- c(mean(response[which(alloc == 1)]), mean(response[which(alloc == 2)]))
  est_p <- c((sum(response[which(alloc == 1)]))/(length(response[which(alloc == 1)])), (sum(response[which(alloc == 2)]))/(length(response[which(alloc == 2)])))
  n <- c(length(response[which(alloc == 1)]), length(response[which(alloc == 2)]))
  
  est_VE <- 1 - est_p[2]/est_p[1]
  est_sigma <- sqrt((1-est_p[1])/(est_p[1]*n[1])+(1-est_p[2])/(est_p[2]*n[2]))
  z <- qnorm(alpha/2, lower.tail = FALSE)
  d <- z*est_sigma
  
  W <- (1-est_VE)*(exp(d)-exp(-d))
  
  i = 2*seed
  while(W > desired_W || is.na(W))
  {
    alloc[i+1] <- sample(c(1,2), 1, prob = c(1/3,2/3))
    p <- p*drift
    #print(p)
    response[i+1] <- rbinom(1, 1, p[alloc[i+1]])
    
    est_p <- c((sum(response[which(alloc == 1)]))/(length(response[which(alloc == 1)])), (sum(response[which(alloc == 2)]))/(length(response[which(alloc == 2)])))
    n <- c(length(response[which(alloc == 1)]), length(response[which(alloc == 2)]))
    
    est_VE <- 1 - est_p[2]/est_p[1]
    est_sigma <- sqrt((1-est_p[1])/(est_p[1]*n[1])+(1-est_p[2])/(est_p[2]*n[2]))
    z <- qnorm(alpha/2, lower.tail = FALSE)
    d <- z*est_sigma
    
    W <- (1-est_VE)*(exp(d)-exp(-d))
    
    i = i+1
  }
  
  est_beta <- log(est_p[2]/est_p[1])
  est_p_0 <- (n[1]*est_p[1]+n[2]*est_p[2])/(n[1]+n[2])
  est_sigma_0 <- sqrt((1-est_p_0)/est_p_0*(1/n[1]+1/n[2]))
  Z <- est_beta/est_sigma_0
  z_c <- -qnorm(0.05, lower.tail = FALSE)
  rej_ind <- if_else(Z < z_c, 1, 0)
  ret <- c(sum(n), n[2], rej_ind)
  return(ret)
}
#s <- sample_size_sim_drift(p = p, desired_W = 0.5, target = "Equal", gamma = 0)
#pboptions(type="txt", char = "=")
#Sample_Sizes_sim <- pbreplicate(1000, {sample_size_sim(p = p, desired_W = 0.24)})
nsim <- 1000

# SMLE with drift
Sample_Sizes_sim_drift_Equal <- matrix(NA, nrow = nsim, ncol = 3)
Sample_Sizes_sim_drift_Double <- matrix(NA, nrow = nsim, ncol = 3)
Sample_Sizes_sim_drift_Neyman_SMLE <- matrix(NA, nrow = nsim, ncol = 3)
Sample_Sizes_sim_drift_RSIHR_SMLE <- matrix(NA, nrow = nsim, ncol = 3)

for(i in c(1:nsim))
{
  print(i)
  Sample_Sizes_sim_drift_Equal[i,] <- sample_size_sim_drift(p=p,desired_W = W,drift = 0.999, target = "Equal", gamma = 0)
  
  #Sample_Sizes_sim_drift_Double[i,] <- sample_size_sim_drift_Double(p = p, desired_W = W)
  #Sample_Sizes_sim_drift_Neyman_SMLE[i,] <- sample_size_sim_drift(p=p,desired_W = W,target = "Neyman", gamma = 0)
  #Sample_Sizes_sim_drift_RSIHR_SMLE[i,] <- sample_size_sim_drift(p=p,desired_W = W,target = "RSIHR", gamma = 0)
}

Table_drift_SMLE <- data.frame(cbind(Sample_Sizes_sim_drift_Equal, Sample_Sizes_sim_drift_Double, Sample_Sizes_sim_drift_Neyman_SMLE, Sample_Sizes_sim_drift_RSIHR_SMLE))

colnames(Table_drift_SMLE) <- c(rep("Equal",3),rep("Double",3),rep("Neyman",3),rep("RSIHR",3))

write.csv(Table_drift_SMLE, "Drift_SMLE_Asy_Sim_Results.csv", row.names = FALSE)

#DBCD with drift
Sample_Sizes_sim_drift_Neyman_DBCD <- matrix(NA, nrow = nsim, ncol = 3)
Sample_Sizes_sim_drift_RSIHR_DBCD <- matrix(NA, nrow = nsim, ncol = 3)    

for(i in c(1:nsim))
{
  print(i)
  Sample_Sizes_sim_drift_Neyman_DBCD[i,] <- sample_size_sim_drift(p=p,desired_W = W,target = "Neyman", gamma = 2)
  Sample_Sizes_sim_drift_RSIHR_DBCD[i,] <- sample_size_sim_drift(p=p,desired_W = W,target = "RSIHR", gamma = 2)
}

Table_drift_DBCD <- data.frame(cbind(Sample_Sizes_sim_drift_Neyman_DBCD, Sample_Sizes_sim_drift_RSIHR_DBCD))

colnames(Table_drift_DBCD) <- c(rep("Neyman",3),rep("RSIHR",3))

write.csv(Table_drift_DBCD, "Drift_DBCD_Asy_Sim_Results.csv", row.names = FALSE)

# ERADE with drift
Sample_Sizes_sim_drift_Neyman_ERADE <- matrix(NA, nrow = nsim, ncol = 3)
Sample_Sizes_sim_drift_RSIHR_ERADE <- matrix(NA, nrow = nsim, ncol = 3)

for(i in c(1:nsim))
{
  print(i)
  Sample_Sizes_sim_drift_Neyman_ERADE[i,] <- sample_size_sim_drift_ERADE(p=p,desired_W = W,target = "Neyman")
  Sample_Sizes_sim_drift_RSIHR_ERADE[i,] <- sample_size_sim_drift_ERADE(p=p,desired_W = W,target = "RSIHR")
}

Table_drift_ERADE <- data.frame(cbind(Sample_Sizes_sim_drift_Neyman_ERADE, Sample_Sizes_sim_drift_RSIHR_ERADE))

colnames(Table_drift_ERADE) <- c(rep("Neyman",3),rep("RSIHR",3))

write.csv(Table_drift_ERADE, "Drift_ERADE_Asy_Sim_Results.csv", row.names = FALSE)

Table_SMLE_sum <- Table_SMLE %>% 
  summarise(Equal_mean_tot = mean(Equal),
            Equal_mean_vac = mean(Equal.1),
            Equal_sim_power = mean(Equal.2),
            Double_mean_tot = mean(Double),
            Double_mean_vac = mean(Double.1),
            Double_sim_power = mean(Double.2),
            Neyman_mean_tot = mean(Neyman),
            Neyman_mean_vac = mean(Neyman.1),
            Neyman_sim_power = mean(Neyman.2),
            RSIHR_mean_tot = mean(RSIHR),
            RSIHR_mean_vac = mean(RSIHR.1),
            RSIHR_sim_power = mean(RSIHR.2))

# for(i in c(1:100)){
#   p <- p*drift
#   print((1-p)/p)
# }
