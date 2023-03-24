library(dplyr)
library(data.table)

#Risk Ratio
p_A <- c(0.1, 0.01, 0.001)
VE <- 0.8
#si <- 0.6
#p_B <- (1-VE)*p_A
#rho_Neyman <- sqrt((1-p_A)/p_A)/(sqrt((1-p_A)/p_A)+sqrt((1-p_B)/p_B))
#rho_equal <- 0.5
#rho_RSIHR <- sqrt(p_B)/(sqrt(p_B)+sqrt(p_A))
#rho_RPW <- (1-p_B)/(2-p_A-p_B)
RW <- seq(1, 0.4, by = -0.2)
#d <- log((RW*VE/(1-VE) + sqrt((RW*VE/(1-VE))^2+4))/2)
# z <- qnorm(0.025, lower.tail = FALSE)
# n <- (z/d)^2*((1-p_A)/(p_A*rho)+(1-p_B)/(p_B*(1-rho)))

sample_size <- function(p_A, VE, type=c("Equal", "Double", "Neyman", "RSIHR"), RW, sam = c("Total", "Vaccine"), alpha = 0.05)
{
  p_B <- (1-VE)*p_A
  d <- log((RW*VE/(1-VE) + sqrt((RW*VE/(1-VE))^2+4))/2)
  rho <- case_when(type == "Equal" ~ 0.5,
                   type == "Double" ~ 1/3,
                   type == "Neyman" ~ sqrt((1-p_A)/p_A)/(sqrt((1-p_A)/p_A)+sqrt((1-p_B)/p_B)),
                   type == "RSIHR" ~ (sqrt(1-p_A)*p_B)/(sqrt(1-p_A)*p_B + sqrt(1-p_B)*p_A))
  z <- qnorm(alpha/2, lower.tail = FALSE)
  n <- (z/d)^2*((1-p_A)/(p_A*rho)+(1-p_B)/(p_B*(1-rho)))
  if(sam == "Total")
  {
    return(round(n))
  }else{
    return(round(n*(1-rho)))
  }
}
#n <- sample_size(p_A = 0.01, VE = 0.4, type = "RSIHR", RW = 1, sam = "Total")
Equal_n_total <- sapply(p_A, sample_size, VE=VE, type = "Equal", RW = RW, sam = "Total")
Equal_n_vac <- sapply(p_A, sample_size, VE=VE, type = "Equal", RW = RW, sam = "Vaccine")
Equal_n <- matrix(paste(Equal_n_total, " (", Equal_n_vac, ")", sep = ""), nrow = length(RW));Equal_n
colnames(Equal_n) <- p_A
Equal_n_Data <- cbind(data.frame(RW = RW, Type = "Equal"), Equal_n)

Double_n_total <- sapply(p_A, sample_size, VE=VE, type = "Double", RW = RW, sam = "Total")
Double_n_vac <- sapply(p_A, sample_size, VE=VE, type = "Double", RW = RW, sam = "Vaccine")
Double_n <- matrix(paste(Double_n_total, " (", Double_n_vac, ")", sep = ""), nrow = length(RW));Double_n
colnames(Double_n) <- p_A
Double_n_Data <- cbind(data.frame(RW = RW, Type = "Double"), Double_n)

Neyman_n_total <- sapply(p_A, sample_size, VE=VE, type = "Neyman", RW = RW, sam = "Total")
Neyman_n_vac <- sapply(p_A, sample_size, VE=VE, type = "Neyman", RW = RW, sam = "Vaccine")
Neyman_n <- matrix(paste(Neyman_n_total, " (", Neyman_n_vac, ")", sep = ""), nrow = length(RW));Neyman_n
colnames(Neyman_n) <- p_A
Neyman_n_Data <- cbind(data.frame(RW = RW, Type = "Neyman"), Neyman_n)

RSIHR_n_total <- sapply(p_A, sample_size, VE=VE, type = "RSIHR", RW = RW, sam = "Total")
RSIHR_n_vac <- sapply(p_A, sample_size, VE=VE, type = "RSIHR", RW = RW, sam = "Vaccine")
RSIHR_n <- matrix(paste(RSIHR_n_total, " (", RSIHR_n_vac, ")", sep = ""), nrow = length(RW));RSIHR_n
colnames(RSIHR_n) <- p_A
RSIHR_n_Data <- cbind(data.frame(RW = RW, Type = "RSIHR"), RSIHR_n)


All_n <- rbind(Equal_n_Data, Double_n_Data, Neyman_n_Data, RSIHR_n_Data)
All_n_Data <- All_n %>% 
  arrange(desc(RW))

#write.csv(All_n_Data, "VE0.8.csv", row.names = FALSE)

VE <- c(0.8,0.6,0.4,0.3)
p_A <- 0.01
p_B <- (1-VE)*p_A
W <- 0.24
d <- log((W/(1-VE) + sqrt((W/(1-VE))^2+4))/2)
z <- qnorm(0.025, lower.tail = FALSE)
LCL <- 1-(1-VE)*exp(d)
UCL <- 1-(1-VE)*exp(-d)
rho_Neyman <- sqrt((1-p_A)/p_A)/(sqrt((1-p_A)/p_A)+sqrt((1-p_B)/p_B))
rho_equal <- 0.5
rho_RSIHR <- (sqrt(1-p_A)*p_B)/(sqrt(1-p_A)*p_B + sqrt(1-p_B)*p_A)
rho_Double <- 1/3
#n <- (z/d)^2*((1-p_A)/(p_A*rho)+(1-p_B)/(p_B*(1-rho)))
Table <- data.frame(VE) %>% 
  mutate("RW" = W/VE,
         "LCL" = round(1-(1-VE)*exp(d),2),
         "UCL" = round(1-(1-VE)*exp(-d),2),
         "Equal_Total" = round((z/d)^2*((1-p_A)/(p_A*rho_equal)+(1-p_B)/(p_B*(1-rho_equal)))),
         "Week_E" = floor(Equal_Total/1000)+1,
         "Vaccine_prop_Equal" = round(1-rho_equal,2),
         "Double_Total" = round((z/d)^2*((1-p_A)/(p_A*rho_Double)+(1-p_B)/(p_B*(1-rho_Double)))),
         "Week_D" = floor(Double_Total/1000)+1,
         "Vaccine_prop_Double" = round(1-rho_Double, 2),
         "Percent_reduction_Double" = round((Equal_Total-Double_Total)/Equal_Total * 100,2),
         "Neyman_Total" = round((z/d)^2*((1-p_A)/(p_A*rho_Neyman)+(1-p_B)/(p_B*(1-rho_Neyman)))),
         "Week_N" = floor(Neyman_Total/1000)+1,
         "Vaccine_prop_Neyman" = round(1-rho_Neyman, 2),
         "Percent_reduction_Neyman" = round((Equal_Total-Neyman_Total)/Equal_Total * 100,2),
         "RSIHR_Total" = round((z/d)^2*((1-p_A)/(p_A*rho_RSIHR)+(1-p_B)/ (p_B*(1-rho_RSIHR)))),
         "Week_R" = floor(RSIHR_Total/1000 )+1,
         "Vaccine_prop_RSIHR" = round(1-rho_RSIHR, 2),
         "Percent_reduction_RSIHR" = round((Equal_Total-RSIHR_Total)/Equal_Total * 100,2))   
#write.csv(Table, "CI_comp.csv", row.names = FALSE)

#Odds ratio

p_A <- seq(0.1,0.5, by = 0.1)
VE <- 0.8
# p_B <- p_A*(1-VE)/(1-p_A*VE)
# c <- 1
RW <- 1
# d <- log((RW*VE/(1-VE) + sqrt((RW*VE/(1-VE))^2+4))/2)
# z <- qnorm(0.025, lower.tail = FALSE)
#
# n_B <- (z/d)^2*(1/(p_B*(1-p_B))+1/(c*p_A*(1-p_A)))

case_sample_size <- function(p_A, VE, type = c("Equal", "Double", "Quadruple", "Neyman"), RW, alpha = 0.05){
  
  p_B <- p_A*(1-VE)/(1-p_A*VE)
  d <- log((RW*VE/(1-VE) + sqrt((RW*VE/(1-VE))^2+4))/2)
  c <- case_when(type == "Equal" ~ 1,
                 type == "Double" ~ 2,
                 type == "Quadruple" ~ 4,
                 type == "Neyman" ~ sqrt((1-p_B)*p_B)/sqrt((1-p_A)*p_A))
  z <- qnorm(alpha/2, lower.tail = FALSE)
  
  n_B <- (z/d)^2*(1/(p_B*(1-p_B))+1/(c*p_A*(1-p_A)))
  return(round(c(n_B*(1+c),n_B)))
}

type =  c("Equal", "Double", "Quadruple", "Neyman")

#n_B_Neyman <- case_sample_size(p_A = 0.1, VE, type = "Neyman", RW); n_B_Neyman
n_B_1 <- sapply(type, case_sample_size, p_A = p_A, VE = VE, RW = RW);n_B_1

n_B_2 <- sapply(type, case_sample_size, p_A = p_A, VE = VE, RW = 0.5);n_B_2

# par(mfrow=c(1,2))
# matplot(p_A,n_B_1, type = "l")
# legend(0.4,100, legend = colnames(n_B_All), fill = c(1:4))
# 
# matplot(p_A,n_B_2, type = "l")
# legend(0.4,300, legend = colnames(n_B_All), fill = c(1:4))

Plot_Table_1 <- cbind(RW = 1, p_A=rep(p_A,2),Prop = c(rep("Total",length(p_A)), rep("Cases",length(p_A))), n_B_1)
Plot_Table_1 <- dcast(as.data.table(Plot_Table_1), RW+p_A ~ Prop, value.var = c("Equal", "Double", "Quadruple", "Neyman"))

Plot_Table_2 <- cbind(RW = 0.5, p_A=rep(p_A,2),Prop = c(rep("Total",length(p_A)), rep("Cases",length(p_A))), n_B_2)
Plot_Table_2 <- dcast(as.data.table(Plot_Table_2), RW+p_A ~ Prop, value.var = c("Equal", "Double", "Quadruple", "Neyman"))

#write.csv(rbind(Plot_Table_1, Plot_Table_2), "VE0.8OR.csv", row.names = FALSE) 
