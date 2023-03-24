library(dplyr)
library(data.table)
#Observed proportions
p_A <- c(101/5829,71/4455,30/1374,112/5829,153/5829)
p_B <- c(30/5807,27/4440,3/1367,37/5807,68/5807)
#Estimate of vaccine efficacy
VE <- 1 - p_B/p_A
#Number of individuals in the vaccinated and unvaccinated group
n_A <- c(5829,4455,1374,5829,5829)
n_B <- c(5807,4440,1367,5807,5807)
sigma2 <- (1-p_A)/(p_A*n_A) + (1-p_B)/(p_B*n_B)
z <- qnorm(0.025, lower.tail = FALSE)
d <- z*sqrt(sigma2)
#Width and Relative width of the 95% CI
W <- (1-VE)*(exp(d)-exp(-d))
RW <- W/VE;RW

sample_size <- function(p_A, VE, RW, type = c("Equal", "Double", "Neyman", "RSIHR"), alpha = 0.05)
{
  p_B <- (1-VE)*p_A
  d <- log((RW*VE/(1-VE) + sqrt((RW*VE/(1-VE)) ^2+4))/2)
  rho <- case_when(type == "Equal" ~ 0.5,
                   type == "Double" ~ (1/3),
                   type == "Neyman" ~ sqrt((1-p_A)/p_A)/(sqrt((1-p_A)/p_A)+sqrt((1-p_B)/p_B)),
                   type == "RSIHR" ~ (sqrt(1-p_A)*p_B)/(sqrt(1-p_A)*p_B + sqrt(1-p_B)*p_A))
  z <- qnorm(alpha/2, lower.tail = FALSE)
  n <- (z/d)^2*((1-p_A)/(p_A*rho)+(1-p_B)/(p_B*(1-rho)))
  Data <- cbind(round(n), round(n*(1-rho)))
  colnames(Data) <- c(paste("Total",type), paste("Vaccine",type))
  return(Data)
}

# Total Allocation vs allocation to the vaccinated group
Equal_n <- sample_size(p_A, VE, type = "Equal", RW);Equal_n
Double_n <- sample_size(p_A, VE, type = "Double", RW);Double_n
Neyman_n <- sample_size(p_A, VE, type = "Neyman", RW); Neyman_n
RSIHR_n <- sample_size(p_A, VE, type = "RSIHR", RW); RSIHR_n

VE <- round(VE,3)
p_A <- round(p_A,4)
RW <- round(RW,4)

All_Data <- as.data.frame(cbind(VE, ARU = p_A, RW,Actual_Total = n_A+n_B, Actual_Vaccine = n_B, Equal_n, Double_n, Neyman_n, RSIHR_n))
#write.csv(All_Data, "Alloc_Data_Sim_results.csv", row.names = FALSE)
