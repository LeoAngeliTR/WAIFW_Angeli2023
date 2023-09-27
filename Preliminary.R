#LIBRARIES----
library(rstudioapi)
library(matrixcalc)
library(popdemo)
library(demography)
library(smoothAPC)
library('plot.matrix')
library(plot3D)
library(blockmatrix)
library(lubridate)
library(readxl)
library(geometry)
library(dplyr)
library(ggplot2)
library(viridisLite)
library(viridis)
library(hrbrthemes)
library(reshape2)
library(ie2misc)
library(ggpubr)
library(expm)
library(jsonlite)
library(xtable)
library(BSDA)
library (readr)
library(scales)
library(data.table)
library(tidyverse)
library(qs)
#PRELIMINARY OPERATIONS-------
##### The following code sets the working directory to the one where the file is

rm(list=ls(all=TRUE))
if(require(rstudioapi) && isAvailable()){
  current_path <- getActiveDocumentContext()$path 
  setwd(dirname(current_path ))
}

#Load data on susceptibility and dates of the CoMix waves
load("./Preliminary_dates_susceptibles.RData")
#Call the function needed
source('./add_functions.R')

##Epidemiological rates## -Using posterior mean values of Abrams et al.------
params <- set_epi_params(3) # 3 corresponds to the age breaks mimicking the school system in Belgium
age_breaks<-params$Age    #age intervals
Asusc_vec <- params$A     #age-specific q-susceptibility
Hinf_vec <- params$H      #age-specific q-infectiousness
p <- params$p_asym        #age-specific proportion of asymptomatic (assuming 50% of the overall population is asymptomatic)
phi_0 <- params$phi_0     #age-specific proportion of mildly symptomatic
omega_1<-params$omega_1   #severe symptoms removal rate
p_mask <- params$p_mask   #age specific mask-usage
prop_ratio <- 0.51        #Infectiousness ratio (Symptomatic/Asymptomatic)
gamma= 0.729              #exposed removal rate 
theta= 0.475              #presym. removal rate
sigma1=0.24               #asympt. removal rate
sigma2=0.756              #recovery rate of mildly symptomatic infections

#Age specific vectors & matrices
Sigma <- sigma2*phi_0
Psi <- sigma2*(1-phi_0)

##Q-SUSCEPTIBILITY
A <- diag(Asusc_vec)

##Q-INFECTIVITY
H <- diag(Hinf_vec)

nbreaks <- length(age_breaks)
I <- diag(rep(1,nbreaks))#identity matrix

#PRE-PANDEMIC CONTACT MATRICES#####
load("./Preliminary_contact_demography.RData")

#Adapt numerically simulated susceptible numbers to the age_group of case=3.
#We define proportion of susceptible in the two different age groups set
#By imposing that the age-distribution over the non-coinciding age intervals is the same
#as the demographic one.
p1 <- N_vec[1]/N_0_90[1]
p2 <- 1-p1
n1 <- N_0_90[1]-N_vec[1]
p3 <- (N_vec[2]-n1)/N_0_90[2]
p4 <- N_vec[3]/N_0_90[2]
n2 <- N_vec[4]-N_0_90[3]
p5 <- n2/N_0_90[2]

props <- c(p1,p2,p3,p4,p5)

SW <- adapt_susc(M=SW_70,age_groups=age_breaks,dates=w_dates,proportions = props)
SW_low <- adapt_susc(M=SW_70_low,age_groups=age_breaks,dates=w_dates,proportions = props)
SW_up <- adapt_susc(M=SW_70_up,age_groups=age_breaks,dates=w_dates,proportions = props)

#NGM COMPOSITION-----


###Calibrate q-factor (residual) on initial R_0=3.4####
##Here evaluate a proportionality factor, non age-specific, such that 
# the correspondent evaluated R_0=3.4 (as in Coletti et al. 2021)
R_0 <- 0
q_sym <- 0.1  
q_asym <- prop_ratio*q_sym
Beta_asym <- matrix()
Beta_sym <- matrix()
K_1 <- matrix()
G=diag(gamma, nbreaks,nbreaks)
O=diag(theta, nbreaks,nbreaks)
S1=diag(sigma1, nbreaks,nbreaks)
S2=diag(Sigma)
R=diag(Psi)
Omega=diag(omega_1)
factor1 <- (theta*sigma1)**(-1)
factor2 <- solve(Omega*(R+S2))
D_1_asym <- (S1+(O*p))%*%diag(factor1,nbreaks,nbreaks)
D_2_sym<- factor2*((Omega+R)*(1-p))

#composing transmission terms
Beta_sym <-A%*%B_sym%*%(q_sym*H)
Beta_asym <- A%*%B_asym%*%(q_asym*H)
Beta_asym <- diag(S_vec)%*%Beta_asym
Beta_sym <- diag(S_vec)%*%Beta_sym
K_1 <- Beta_asym%*%D_1_asym +Beta_sym%*%D_2_sym
R_0 <-max(Re(eigen(K_1)$values))
q=3.4/R_0 
q_asym <- q_asym*q
q_sym <- q_sym*q
Beta_asym <- q*Beta_asym
Beta_sym <- q*Beta_sym
K_1 <- K_1*q

#Resulting dominant eigenstructure
R_0_1 <- max(Re(eigen(K_1)$values))
w_1 <- abs(eigen(K_1)$vectors[,1])
v_1 <- abs(eigen(t(K_1))$vector[,1])
w_1 <- w_1/sum(w_1)  
v_1 <- v_1/sum(v_1*w_1)


###Damping ratio ####
damp_ratio <- R_0_1/abs(eigen(K_1)$values)[2]

#SENSITIVITY ANALYSIS----

#initialize matrices
age_names <- matrixCoMix_start$demography$age.group
di_names <- list(age_names,age_names)
S_built_1  <- matrix(calc.sens(K_1),nbreaks,nbreaks,dimnames=di_names)
E_built_1<- matrix(hadamard.prod(S_built_1,K_1)/R_0_1,nbreaks,nbreaks,dimnames=di_names)

#Overall sensitivity index
S_tot_index<-sqrt(sum(hadamard.prod(S_built_1,S_built_1))) 


#Damping ratio sensitivity
Sens_dr <- damp_ratio_sens(K_1)$Sens
El_dr <- damp_ratio_sens(K_1)$El
rownames(Sens_dr) <- age_names
colnames(Sens_dr)<- age_names
rownames(El_dr) <- age_names
colnames(El_dr)<- age_names






##Compute measures of sensitvity####
#Generates a list of matrices embedding the Rt sensitivity w.r.t. model params
mat_sens <- lower_sens_matrices(K=K_1,               
                                M_asym = B_mean_asym,
                                M_sym = B_mean_sym,
                                Susc = S_vec,
                                q_factor = q_sym)


# #Figure 2a-----
# tiff(filename = "Figure2a.tiff", width = 8, height = 6, units = 'in', res = 500)
# plot_matrix(round(t(K_1), 3), paste("NGM"),age_names,xlab = 'Age infected',ylab = 'Age susceptible')
# dev.off()
# #Figure 2b----
# tiff(filename = "Figure2b.tiff", width = 8, height = 6, units = 'in', res = 500)
# plot_matrix(round(t(mat_sens$Sens_ratio),3), paste("Sensitivity ratios"),age_names, ylab = 'Age (numerator)',xlab = 'Age (denominator)')
# dev.off()
# 
#Table A.3-----
# 
# Sum_S       <- c(Sum=round(sum(S_built_1),2),Norm=round(norm(S_built_1,type = "F"),2) , M= 'Sensitivity', Age_spec= T ,AVG=round(mean(K_1),2), Term="Overall")
# Sum_S_A     <-c(Sum=round(sum(mat_sens$S_A),2),Norm=round(norm(mat_sens$S_A,type = "F"),2), M= 'Sensitivity', Age_spec= T , AVG=round(mean(Asusc_vec),2),Term="q-susceptibility")
# Sum_S_H     <-c(Sum=round(sum(mat_sens$S_H),2),Norm=round(norm(mat_sens$S_H,type = "F"),2), M= 'Sensitivity', Age_spec= T ,AVG= round(mean(Hinf_vec),2), Term="q-infectiousness")
# Sum_S_C     <-c(Sum=round(sum(mat_sens$S_C_asym+mat_sens$S_C_sym),2),Norm=round(norm(mat_sens$S_C_asym+mat_sens$S_C_sym,type = "F"),2), M= 'Sensitivity', Age_spec= T ,AVG=round(mean(avg_cont_matrix),2), Term="Reported contacts")
# Sum_S_C_asym<-c(Sum=round(sum(mat_sens$S_C_asym),2),Norm=round(norm(mat_sens$S_C_asym,type = "F"),2), M= 'Sensitivity', Age_spec= T , AVG=round(mean(B_mean_asym),2) ,Term="Contacts (asym)")
# Sum_S_C_sym <-c(Sum=round(sum(mat_sens$S_C_sym),2),Norm=round(norm(mat_sens$S_C_sym,type = "F"),2), M= 'Sensitivity', Age_spec= T , AVG=round(mean(B_mean_sym),2) ,Term="Contacts (sym)")
# 
# Sum_Theta_S <- c(Sum=round(sum( mat_sens$S_T),2),Norm=round(norm( mat_sens$S_T,type = "F"),2), M= 'Sensitivity', Age_spec= F ,AVG=theta ,Term="Pre_sym. removal rate")
# Sum_D1_S      <-c(Sum=round(sum(mat_sens$S_S1),2),Norm=round(norm(mat_sens$S_S1,type = "F"),2), M= 'Sensitivity', Age_spec= F ,AVG=sigma1 ,Term="Asym. recovery rate")
# Sum_Psi_S     <-c(Sum=round(sum(mat_sens$S_Psi),2),Norm=round(norm(mat_sens$S_Psi,type = "F"),2), M= 'Sensitivity', Age_spec= T ,AVG=round(mean(Psi),2) ,Term="Severe birth rate")
# Sum_D2_S      <-c(Sum=round(sum(mat_sens$S_D2),2),Norm=round(norm(mat_sens$S_D2,type = "F"),2), M= 'Sensitivity', Age_spec= T ,AVG=round(mean(Sigma),2) ,Term="Mild. recovery rate")
# Sum_O_S       <- c(Sum=round(sum(mat_sens$S_O),2),Norm=round(norm(mat_sens$S_O,type = "F"),2), M= 'Sensitivity', Age_spec= T ,AVG=round(mean(omega_1),2) ,Term="Sev. removal rate")
# Sum_Prob_S    <- c(Sum=round(sum(mat_sens$S_P),2),Norm=round(norm(mat_sens$S_P,type = "F"),2), M= 'Sensitivity', Age_spec= T ,AVG=mean(p) ,Term="Prob. remaining asym")
# data_sens<- as.data.frame(rbind(Sum_S,Sum_S_C_asym,Sum_S_C_sym, Sum_S_A,Sum_S_H,Sum_S_C ,Sum_S_C_asym,Sum_S_C_sym,  
#                                 Sum_Theta_S,Sum_D1_S,Sum_Psi_S,Sum_D2_S,Sum_O_S,Sum_Prob_S)) 
# xtable(data_sens, caption = "Parameter-specific $R_0$ sensitivity", label = "tab:tot_sens")
# 
# #Table A.4-----
# Sum_E <- c(Sum=round(sum(E_built_1),2), Norm=round(norm(E_built_1,type = "F"),2), M = 'Elasticity', Age_spec= T ,Term="Overall")
# Sum_E_asym <- c(Sum=round(sum(mat_sens$E_asym),2),Norm=round(norm(mat_sens$E_asym,type = "F"),2), M= 'Elasticity', Age_spec= T ,Term="Overall asym")
# Sum_E_sym <- c(Sum=round(sum(mat_sens$E_sym),2),Norm=round(norm(mat_sens$E_sym,type = "F"),2), M= 'Elasticity', Age_spec= T ,Term="Overall sym")
# Sum_E_Theta <- c(Sum=round(sum( mat_sens$E_T),2),Norm=round(norm( mat_sens$E_T,type = "F"),2), M= 'Elasticity', Age_spec= F ,Term="Pre_sym. removal rate")
# Sum_D1 <-c(Sum=round(sum(mat_sens$E_S1),2), Norm=round(norm(mat_sens$E_S1,type = "F"),2), M= 'Elasticity', Age_spec= F ,Term="Asym. recovery rate")
# Sum_Psi <-c(Sum=round(sum(mat_sens$E_Psi),5),Norm=round(norm(mat_sens$E_Psi,type = "F"),5), M= 'Elasticity', Age_spec= T ,Term="Severe birth rate")
# Sum_D2<-c(Sum=round(sum(mat_sens$E_D2),2), Norm=round(norm(mat_sens$E_D2,type = "F"),2), M= 'Elasticity', Age_spec= T ,Term="Mild. recovery rate")
# Sum_O <- c(Sum=round(sum(mat_sens$E_O),2), Norm=round(norm(mat_sens$E_O,type = "F"),3), M= 'Elasticity', Age_spec= T ,Term="Sev. removal rate")
# Sum_Prob <- c(Sum=round(sum(mat_sens$E_P),2),Norm=round(norm(mat_sens$E_P,type = "F"),2), M= 'Elasticity', Age_spec= T ,Term="Prob. remaining asym")
# 
# data_sens<- as.data.frame(rbind(Sum_E,Sum_E_asym,Sum_E_sym,
#                                 Sum_E_Theta,Sum_D1,Sum_Psi,Sum_D2,Sum_O,Sum_Prob))
# 
# 
# xtable(data_sens, caption = "Parameter-specific $R_0$  elasticity.", label = "tab:tot_el")
# 
# 





#Reproduction numbers####
url <- 'https://gist.githubusercontent.com/gjbex/b495f98d58ea50bb8ad4632f1504972b/raw/r_data_cases.json'
# read url and convert to data.frame
R_est <- fromJSON(txt=url)
#continuous time R_t
trasl <-0
R_cont <- as.data.frame(R_est)
R_cont <- R_cont[which(R_cont$Region=='Belgique'),]
R_cont <- R_cont[which(R_cont$Date>=w_dates[1]+trasl & R_cont$Date<=w_dates[length(w_dates)]+trasl),]

##Discrete time-series----

##R_0
R_series_mean <-c()
R_series_up <-c() 
R_series_low <-c()
R_series_mean_mid <-c()
R_series_up_mid <-c() 
R_series_low_mid <-c()
S_series <- c()
W_series <- c()
W_series_perc <- c()
Flag_series <- c()
days_range <- 3
for(i in 1:length(w_dates)){
  confr_date_Rt <- w_dates[i]+trasl
  R_series_mean[i] <- mean(R_cont$Restim[which(R_cont$Date>confr_date_Rt-days_range & R_cont$Date<confr_date_Rt+days_range)])
  R_series_up[i] <- mean(R_cont$RCIup[which(R_cont$Date>confr_date_Rt-days_range & R_cont$Date<confr_date_Rt+days_range)])
  R_series_low[i] <- mean(R_cont$RCIlow[which(R_cont$Date>confr_date_Rt-days_range & R_cont$Date<confr_date_Rt+days_range)])
  R_series_mean_mid[i] <- mean(R_cont$Restim[which(R_cont$Date>confr_date_Rt-(days_range*2) & R_cont$Date<confr_date_Rt)])
  R_series_up_mid [i]<- mean(R_cont$RCIup[which(R_cont$Date>confr_date_Rt-(days_range*2) & R_cont$Date<confr_date_Rt)])
  R_series_low_mid[i]<- mean(R_cont$RCIlow[which(R_cont$Date>confr_date_Rt-(days_range*2) & R_cont$Date<confr_date_Rt)])
}




