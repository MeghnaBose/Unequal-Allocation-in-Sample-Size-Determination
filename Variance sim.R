library(dplyr)
library(ggplot2)
# Asymptotic variance plot

Asymp_var_DBCD <- function(p_A, VE, alpha, target = c("Neyman", "RSIHR")){
  q_A <- 1-p_A
  p_B <- p_A*(1-VE)
  q_B <- 1-p_B
  var <- case_when(target == "Neyman" ~  (sqrt(p_A*q_A^3)*(p_B*q_B+0.5*(1+alpha))+sqrt(p_B*q_B^3)*(p_A*q_A+0.5*(1+alpha)))/((1+2*alpha)*(sqrt(p_A*q_B)+sqrt(p_B*q_A))^3*sqrt(q_A*q_B)),
                   target == "RSIHR" ~ (p_A*p_B*q_B*sqrt(q_B)*(p_A*q_A + 0.5*(1+alpha)*(1+q_A)^2)+p_A*p_B*q_A*sqrt(q_A)*(p_B*q_B + 0.5*(1+alpha)*(1+q_B)^2))/((1+2*alpha)*(p_A*sqrt(q_B)+sqrt(q_A)*p_B)^3 * sqrt(q_A*q_B)))
    return(var)
}

Asymp_var_ERADE <- function(p_A, VE, target = c("Neyman", "RSIHR")){
  q_A <- 1-p_A
  p_B <- p_A*(1-VE)
  q_B <- 1-p_B
  var <- case_when(target == "Neyman" ~  (sqrt(q_B^3*p_B)+sqrt(q_A^3*p_A))/(4*sqrt(q_A*q_B)*(sqrt(q_A*p_B)+sqrt(p_A*q_B))^3),
                   target == "RSIHR" ~ (p_A*p_B*(sqrt(q_B^3)*(1+q_A)^2+sqrt(q_A^3)*(1+q_B)^2))/(4*sqrt(q_A*q_B)*(sqrt(q_A)*p_B+sqrt(q_B)*p_A)^3))
  return(var)
}

p_A <- seq(0.1, 0.9, by = 0.1)
VE <- seq(0.1, 0.9, by = 0.1)

#SMLE_Neyman <- Asymp_var_DBCD(  p_A=0.9, VE, alpha = 0, target = "Neyman")
#DBCD_Neyman <- Asymp_var_DBCD(p_A, VE, alpha = 2, target = "Neyman")

SMLE_Neyman <- cbind(VE, sapply(p_A, Asymp_var_DBCD, VE = VE, alpha = 0, target = "Neyman"))
colnames(SMLE_Neyman) <- c("VE", as.character(p_A))
DBCD_Neyman <- cbind(VE, sapply(p_A, Asymp_var_DBCD, VE = VE, alpha = 2, target = "Neyman"))
colnames(DBCD_Neyman) <- c("VE", as.character(p_A))
ERADE_Neyman <- cbind(VE, sapply(p_A, Asymp_var_ERADE, VE = VE, target = "Neyman"))
colnames(ERADE_Neyman) <- c("VE", as.character(p_A))

Final_Data <- bind_rows(mutate(data.frame(SMLE_Neyman), type = "SMLE"), mutate(data.frame(DBCD_Neyman), type = "DBCD"), mutate(data.frame(ERADE_Neyman), type = "ERADE"))

Final_Data_Neyman <- arrange(Final_Data, VE)


SMLE_RSIHR <- cbind(VE, sapply(p_A, Asymp_var_DBCD, VE = VE, alpha = 0, target = "RSIHR"))
colnames(SMLE_RSIHR) <- c("VE", p_A)
DBCD_RSIHR <- cbind(VE, sapply(p_A, Asymp_var_DBCD, VE = VE, alpha = 2, target = "RSIHR"))
colnames(DBCD_RSIHR) <- c("VE", p_A)
ERADE_RSIHR <- cbind(VE, sapply(p_A, Asymp_var_ERADE, VE = VE, target = "RSIHR"))
colnames(ERADE_RSIHR) <- c("VE", p_A)

Final_Data <- bind_rows(mutate(data.frame(SMLE_RSIHR), type = "SMLE"), mutate(data.frame(DBCD_RSIHR), type = "DBCD"), mutate(data.frame(ERADE_RSIHR), type = "ERADE"))

Final_Data_RSIHR <- arrange(Final_Data, VE)

#write.csv(Final_Data_Neyman, "Neyman Var.csv", row.names = FALSE)
#write.csv(Final_Data_RSIHR, "RSIHR Var.csv", row.names = FALSE)

# df <- data.frame(p_A, SMLE_Neyman, DBCD_Neyman)
# g <- ggplot(df, aes(p_A))
# g <- g + geom_line(aes(y=SMLE_Neyman), colour="red")
# g <- g + geom_line(aes(y=DBCD_Neyman), colour="green") + labs(x = "p_A", y = "Asymptotic Variance", color = "Legend") + scale_color_manual(values = colors) 
# g

Asymp_var_DBCD <- function(p_B, p_A, alpha, target = c("Neyman", "RSIHR")){
  q_A <- 1-p_A
  #p_B <- x
  q_B <- 1-p_B
  var <- case_when(target == "Neyman" ~  (sqrt(p_A*q_A^3)*(p_B*q_B+0.5*(1+alpha))+sqrt(p_B*q_B^3)*(p_A*q_A+0.5*(1+alpha)))/((1+2*alpha)*(sqrt(p_A*q_B)+sqrt(p_B*q_A))^3*sqrt(q_A*q_B)),
                   target == "RSIHR" ~ (p_A*p_B*q_B*sqrt(q_B)*(p_A*q_A + 0.5*(1+alpha)*(1+q_A)^2)+p_A*p_B*q_A*sqrt(q_A)*(p_B*q_B + 0.5*(1+alpha)*(1+q_B)^2))/((1+2*alpha)*(p_A*sqrt(q_B)+sqrt(q_A)*p_B)^3 * sqrt(q_A*q_B)))
  return(var)
} 

Asymp_var_ERADE <- function(p_B, p_A, target = c("Neyman", "RSIHR")){
  q_A <- 1-p_A
  q_B <- 1-p_B
  var <- case_when(target == "Neyman" ~  (sqrt(q_B^3*p_B)+sqrt(q_A^3*p_A))/(4*sqrt(q_A*q_B)*(sqrt(q_A*p_B)+sqrt(p_A*q_B))^3),
                   target == "RSIHR" ~ (p_A*p_B*(sqrt(q_B^3)*(1+q_A)^2+sqrt(q_A^3)*(1+q_B)^2))/(4*sqrt(q_A*q_B)*(sqrt(q_A)*p_B+sqrt(q_B)*p_A)^3))
  return(var)
}

ggplot(data.frame(p_B=c(0, 0.01)), aes(x=p_B)) + 
  stat_function(fun=Asymp_var_DBCD, args = list(p_A=0.01, alpha = 0, target = "Neyman"))+stat_function(fun=Asymp_var_DBCD, args = list(p_A=0.01, alpha = 2, target = "Neyman"))+ stat_function(fun=Asymp_var_ERADE, args = list(p_A=0.01, target = "Neyman"))

ggplot(data.frame(p_A=c(0.01, 1)), aes(x=p_A)) + 
  stat_function(fun=Asymp_var_DBCD, args = list(p_B=0.01, alpha = 0, target = "Neyman"))+stat_function(fun=Asymp_var_DBCD, args = list(p_B=0.01, alpha = 2, target = "Neyman"))+ stat_function(fun=Asymp_var_ERADE, args = list(p_B=0.01, target = "Neyman"))



#eq = function(x){x*x}
ggplot(data.frame(x=c(1, 50)), aes(x=x)) + 
  stat_function(fun=eq)
