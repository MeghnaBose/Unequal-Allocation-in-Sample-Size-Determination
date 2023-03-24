library(dplyr)

sample_size_Z <- function(rho,p_A,VE,alpha = 0.05, gamma = 0.05){
  p_B <- (1-VE)*p_A
  beta <- log(p_B/p_A)
  n <- ((qnorm(alpha, lower.tail = FALSE)*sqrt((rho*(1-p_A)+(1-rho)*(1-p_B))/(rho*p_A+(1-rho)*p_B))*sqrt(1/(rho*(1-rho)))+qnorm(gamma, lower.tail = FALSE)*sqrt((1-p_A)/(rho*p_A)+(1-p_B)/((1-rho)*p_B)))/beta)^2
  return(n)
}

Exp_case_Z <- function(rho,p_A,VE,alpha= 0.05,gamma = 0.05){
  p_B <- (1-VE)*p_A
  beta <- log(p_B/p_A)
  n <- ((qnorm(alpha, lower.tail = FALSE)*sqrt((rho*(1-p_A)+(1-rho)*(1-p_B))/(rho*p_A+(1-rho)*p_B))*sqrt(1/(rho*(1-rho)))+qnorm(gamma, lower.tail = FALSE)*sqrt((1-p_A)/(rho*p_A)+(1-p_B)/((1-rho)*p_B)))/beta)^2
  E <- n*(rho*p_A+(1-rho)*p_B)
  return(E)
}

power <- function(T, theta_0, theta_1, alpha = 0.05){
  Y_q <- qbinom(alpha, size = T, prob = theta_0)-1
  Y_c <- if_else(pbinom(Y_q, size = T, prob = theta_0)<=0.05, Y_q, Y_q-1)
  power <- pbinom(Y_c, size = T, prob = theta_1)
  return(power)
}

sample_size_C <- function(rho, p_A, VE, alpha = 0.05, gamma = 0.05, seed = 1)
{
  p_B <- (1-VE)*p_A
  c <- (1-rho)/rho
  
  theta_0 <- 1/(1+c)
  theta_1 <- (1-VE)/(1+c-VE)
  
  T_range <- seed
  repeat{
    if(power(T_range, theta_0 , theta_1) < (1-gamma) & power(T_range+10, theta_0, theta_1)> (1-gamma)){
      break
    }
    T_range <- T_range+10
  }
  
  power_range <- c(T_range:(T_range+10))
  power_seq <- sapply(c(T_range:(T_range+10)), power, theta_0, theta_1)
  T <- power_range[abs(power_seq-(1-gamma))==min(abs(power_seq-(1-gamma)))]
  
  n_B <- T/((c+1-VE)*p_A)
  n <- (1+c)*n_B
  return(n)
}

VE <- 0.5
p_A <- c(0.1,0.05,0.01)

#optim_Neyman <- optimize(sample_size_Z, interval = c(0,1), p_A = 0.001, VE = VE);optim_Neyman
optim_Neyman <- lapply(p_A, function(x) {optimize(sample_size_Z, interval = c(0,1), p_A = x, VE = VE)})
#optim_RSIHR <- optimize(Exp_case_Z, interval = c(0,1), p_A = p_A, VE = VE);optim_RSIHR
optim_RSIHR <- lapply(p_A, function(x) {optimize(Exp_case_Z, interval = c(0,1), p_A = x, VE = VE)})

size_Z_Equal <- sample_size_Z(0.5,p_A,VE)
size_Z_Double <- sample_size_Z(1/3,p_A,VE)
size_Z_Neyman <- sapply(optim_Neyman,"[[",2)
size_Z_RSIHR <- sample_size_Z(sapply(optim_RSIHR,"[[",1),p_A,VE)

rho_Neyman <- sapply(optim_Neyman,"[[",1)
rho_RSIHR <- sapply(optim_RSIHR,"[[",1)

size_C_Equal <- sample_size_C(0.5, p_A, VE)
size_C_Double <- sample_size_C(1/3,p_A,VE)

Table <- data.frame(p_A=p_A) %>% 
  mutate(Equal_Z = paste(round(size_Z_Equal),"(",0.5,")", sep =""),
         Double_Z = paste(round(size_Z_Double),"(",0.67,")", sep =""),
         Neyman_Z = paste(round(size_Z_Neyman),"(",round(1-rho_Neyman,2),")", sep =""),
         RSIHR_Z = paste(round(size_Z_RSIHR),"(",round(1-rho_RSIHR,2),")", sep =""),
         Equal_C = paste(round(size_C_Equal),"(",0.5,")", sep =""),
         Double_C = paste(round(size_C_Double),"(",0.67,")", sep =""))

#write.csv(Table,"Cond_VE0.9.csv", row.names = FALSE)

