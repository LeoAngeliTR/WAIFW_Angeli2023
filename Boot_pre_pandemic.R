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

#LOAD DATA---- 
#Load baseline data workspaces
load('./Preliminary_dates_susceptibles.RData')
load('./Preliminary_contact_demography.RData')
load('./Pre_pandemic_naive_susc.RData')

##Load additional functions
source('./add_functions.R')


nboot <-1000
# #we use the social contact data from 2011 (Van Hoang et al.)
Bootstrapped <- Boot_baseline(nboot)
matrices_boot <- Bootstrapped$Matrices
eigen_boot<- Bootstrapped$Eigen_str


######Bootstrapped CIs limits####
summed_c_up <- vector(mode = "double", length = nbreaks)
summed_c_mean <-vector(mode = "double", length = nbreaks)
summed_c_median <-vector(mode = "double", length = nbreaks)
summed_c_low <- vector(mode = "double", length = nbreaks)
store_c <-matrix(nrow = nbreaks,ncol = nboot)

# 
summed_csym_up <- vector(mode = "double", length = nbreaks)
summed_csym_mean <-vector(mode = "double", length = nbreaks)
summed_csym_low <- vector(mode = "double", length = nbreaks)
store_csym <-matrix(nrow = nbreaks,ncol = nboot)

# 
summed_el_up <- vector(mode = "double", length = nbreaks)
summed_el_mean <- vector(mode = "double", length = nbreaks)
summed_el_median <- vector(mode = "double", length = nbreaks)
summed_el_low <- vector(mode = "double", length = nbreaks)
store_e <-matrix(nrow = nbreaks,ncol = nboot)

summed_sens_up <- vector(mode = "double", length = nbreaks)
summed_sens_mean <- vector(mode = "double", length = nbreaks)
summed_sens_median <- vector(mode = "double", length = nbreaks)
summed_sens_low <- vector(mode = "double", length = nbreaks)
store_s <-matrix(nrow = nbreaks,ncol = nboot)

summed_sens_asym_up <- vector(mode = "double", length = nbreaks)
summed_sens_asym_mean <- vector(mode = "double", length = nbreaks)
summed_sens_asym_low <- vector(mode = "double", length = nbreaks)
summed_sens_sym_up <- vector(mode = "double", length = nbreaks)
summed_sens_sym_mean <- vector(mode = "double", length = nbreaks)
summed_sens_sym_low <- vector(mode = "double", length = nbreaks)
store_s_asym <-matrix(nrow = nbreaks,ncol = nboot)
store_s_sym <-matrix(nrow = nbreaks,ncol = nboot)


summed_e_asym_up <- vector(mode = "double", length = nbreaks)
summed_e_asym_mean <- vector(mode = "double", length = nbreaks)
summed_e_asym_median <- vector(mode = "double", length = nbreaks)
summed_e_asym_low <- vector(mode = "double", length = nbreaks)

summed_e_sym_up <- vector(mode = "double", length = nbreaks)
summed_e_sym_mean <- vector(mode = "double", length = nbreaks)
summed_e_sym_median <- vector(mode = "double", length = nbreaks)
summed_e_sym_low <- vector(mode = "double", length = nbreaks)

store_e_asym <-matrix(nrow = nbreaks,ncol = nboot)
store_e_sym <-matrix(nrow = nbreaks,ncol = nboot)

#consider susceptible side
summed_e_asym_susc_up<- vector(mode = "double", length = nbreaks)
summed_e_asym_susc_mean<- vector(mode = "double", length = nbreaks)
summed_e_asym_susc_low<- vector(mode = "double", length = nbreaks)

summed_e_sym_susc_up<- vector(mode = "double", length = nbreaks)
summed_e_sym_susc_mean<- vector(mode = "double", length = nbreaks)
summed_e_sym_susc_low<- vector(mode = "double", length = nbreaks)

store_e_asym_susc <-matrix(nrow = nbreaks,ncol = nboot)
store_e_sym_susc <-matrix(nrow = nbreaks,ncol = nboot)

# ##
boot_e_asym_18_60 <- vector(mode = "double", length = nboot)
boot_e_asym_18_50 <- vector(mode = "double", length = nboot)
boot_e_asym_30_50 <- vector(mode = "double", length = nboot)
boot_e_asym_tot <- vector(mode = "double", length = nboot)
boot_e_sym_tot <- vector(mode = "double", length = nboot)

#Eigenvectors
list_boot_w <- vector(mode = "list", length = nboot)
list_boot_v <- vector(mode = "list", length = nboot)



for (i in seq_len(nboot)) {
  store_c[,i] <- apply(matrices_boot[[i]]$C_mean_asym,1,sum)
  store_csym[,i] <- apply(matrices_boot[[i]]$C_mean_sym,1,sum)
  store_e[,i] <- apply(matrices_boot[[i]]$E,2,sum)
  store_s[,i] <- apply(matrices_boot[[i]]$S,2,sum)
  store_e_asym[,i] <- apply(matrices_boot[[i]]$E_asym,2,sum)#(apply(matrices_boot[[i]]$E_asym,2,sum)+apply(matrices_boot[[i]]$E_asym,1,sum))#/2
  store_e_sym[,i] <- apply(matrices_boot[[i]]$E_sym,2,sum)#(apply(matrices_boot[[i]]$E_sym,2,sum)+apply(matrices_boot[[i]]$E_sym,1,sum))#/2
  store_s_asym[,i]  <- apply(matrices_boot[[i]]$S_asym,2,sum)
  store_s_sym[,i]  <- apply(matrices_boot[[i]]$S_sym,2,sum)
  
  store_e_asym_susc[,i] <- apply(matrices_boot[[i]]$E_asym,1,sum)
  store_e_sym_susc[,i] <- apply(matrices_boot[[i]]$E_sym,1,sum)
  
  #Populate eigenvectors lists
  list_boot_w[[i]] <- eigen_boot[[i]]$wb
  list_boot_v[[i]] <- eigen_boot[[i]]$vb
  
  
}
for (j in seq_len(nbreaks)) {
  summed_c_up[j]   <- quantile(sort(store_c[j,]), probs=c(0.975))
  summed_c_mean[j] <- mean(store_c[j,])
  summed_c_median[j] <- median(store_c[j,])
  summed_c_low[j]  <- quantile(sort(store_c[j,]), probs=c(0.025))
  
  summed_csym_up[j]   <- quantile(sort(store_csym[j,]), probs=c(0.975))
  summed_csym_mean[j] <- mean(store_csym[j,])
  summed_csym_low[j]  <- quantile(sort(store_csym[j,]), probs=c(0.025))
  
  summed_el_up[j]  <- quantile(sort(store_e[j,]), probs=c(0.975))
  summed_el_mean[j] <- mean(store_e[j,])
  summed_el_median[j] <-median(store_e[j,])
  summed_el_low[j] <- quantile(sort(store_e[j,]), probs=c(0.025))
  
  
  summed_sens_up[j]   <- quantile(sort(store_s[j,]), probs=c(0.975))
  summed_sens_mean[j] <- mean(store_s[j,])
  summed_sens_median[j] <-median(store_s[j,])
  summed_sens_low[j]  <- quantile(sort(store_s[j,]), probs=c(0.025))

  summed_e_asym_up[j]  <- quantile(sort(store_e_asym[j,]), probs=c(0.975))
  summed_e_asym_mean[j] <- mean(store_e_asym[j,])
  summed_e_asym_median[j] <- median(store_e_asym[j,])
  summed_e_asym_low[j] <- quantile(sort(store_e_asym[j,]), probs=c(0.025))
  
  
 
  summed_e_sym_up[j]   <- quantile(sort(store_e_sym[j,]), probs=c(0.975))
  summed_e_sym_mean[j] <- mean(store_e_sym[j,])
  summed_e_sym_median[j] <- median(store_e_sym[j,])
  summed_e_sym_low[j]  <- quantile(sort(store_e_sym[j,]), probs=c(0.025))
  
 
  
  
  
}

# 

##Eigenstructure 

right_eigs_w <- do.call(rbind,list_boot_w)
left_eigs_v <- do.call(rbind,list_boot_v)

eig_w_mean <- apply(right_eigs_w, 2, mean)
eig_w_low <- apply(right_eigs_w, 2, quantile,prob=0.025)
eig_w_up <- apply(right_eigs_w, 2, quantile,prob=0.975)

eig_v_mean <- apply(left_eigs_v, 2, mean)
eig_v_low <- apply(left_eigs_v, 2, quantile,prob=0.025)
eig_v_up <- apply(left_eigs_v, 2, quantile,prob=0.975)


#Preparing data to plot#####

measures_pre_pandemic <-data.frame(Age=character(),
                                   Contacts=double(),
                                   Cont_up=double(),
                                   Cont_mean=double(),
                                   Cont_median=double(),
                                   Cont_low=double(),
                                   
                                   Contacts_sym=double(),
                                   Cont_sym_up=double(),
                                   Cont_sym_mean=double(),
                                   Cont_sym_low=double(),
                                   
                                   Sens=double(),
                                   El=double(),
                                   
                                   S_up=double(),
                                   S_mean=double(),
                                   S_median=double(),
                                   S_low=double(),
                                   E_up=double(),
                                   E_mean=double(),
                                   E_median=double(),
                                   E_low=double(),
                                   
                                   E_asym_up=double(),
                                   E_asym_mean=double(),
                                   E_asym_median=double(),
                                   E_asym_low=double(),
                                  
                                   E_sym_up=double(),
                                   E_sym_mean=double(),
                                   E_sym_median=double(),
                                   E_sym_low=double(),
                                   stringsAsFactors=FALSE
                                   )

for (i in 1:length(age_names)){
  new_record <- list(age_names[i],
                     summed_contacts[i],
                     summed_c_up[i],
                     summed_c_mean[i],
                     summed_c_median[i],
                     summed_c_low[i],
                     
                     summed_contacts_sym[i],
                     summed_csym_up[i],
                     summed_csym_mean[i],
                     summed_csym_low[i],
                     
                     summed_sens[i],
                     summed_el[i],
                     
                     summed_sens_up[i],
                     summed_sens_mean[i],
                     summed_sens_median[i],
                     summed_sens_low[i],
                     summed_el_up[i],
                     summed_el_mean[i],
                     summed_el_median[i],
                     summed_el_low[i],
                   
                     summed_e_asym_up[i],
                     summed_e_asym_mean[i],
                     summed_e_asym_median[i],
                     summed_e_asym_low[i],
                    
                     summed_e_sym_up[i],
                     summed_e_sym_mean[i],
                     summed_e_sym_median[i],
                     summed_e_sym_low[i]
                     
  )
  
  measures_pre_pandemic[nrow(measures_pre_pandemic) + 1,] = new_record
}

measures_pre_pandemic$Age<-factor(measures_pre_pandemic$Age, levels =age_names)


#PLOTS----

#Sens2
graph <- ggplot(measures_pre_pandemic, aes(x=Age, y=S_mean,fill=S_mean)) +geom_bar(stat = "identity")
graph + geom_text(aes(label=round(S_mean,2)), vjust=-0.65, hjust=0.001, size=7.5)+
  geom_errorbar(aes(x=Age, ymin=S_low, ymax=S_up), width=0.3, colour="darkred", alpha=0.9, size=1.2)+
  labs(title="", x ="Age classes", y = bquote('Cumulative sensitivities '~tilde(s[i])))+
  scale_fill_gradient(low = "gold1", high = 'red2')+
  theme_minimal()+coord_flip()+geom_point(aes(x=Age, y=S_mean))+
  theme(panel.grid.minor = element_blank(),axis.text = element_text(size = 20),
        axis.title.y =element_text(face = 'bold', size = 20),axis.title.x = element_text(size=20),legend.position="none")


#Elast # Figure 2c
graph <- ggplot(measures_pre_pandemic, aes(x=Age, y=E_mean,fill=E_mean)) +geom_bar(stat = "identity")
graph + geom_text(aes(label=paste0(round(E_mean,2)*100,'%')), vjust=-0.65, hjust=0.001, size=7.5)+
  geom_errorbar(aes(x=Age, ymin=E_low, ymax=E_up), width=0.3, colour="darkred", alpha=0.9, size=1.2)+
  labs(title="", x ="Age classes", y = bquote('Cumulative elasticities '~tilde(e[i])))+
  scale_fill_gradient(low = "gold1", high = 'red2')+
  theme_minimal()+coord_flip()+
  theme(panel.grid.minor = element_blank(),axis.text = element_text(size = 20),
        axis.title.y =element_text(face = 'bold', size = 20),axis.title.x = element_text(size=20),legend.position="none")

#Contacts # Figure 2d
graph <- ggplot(measures_pre_pandemic, aes(x=Age, y=Cont_mean/sum(Cont_mean),fill=Cont_mean/sum(Cont_mean))) + geom_bar(stat = "identity")
graph + geom_text(aes(label=paste0(round(Cont_mean/sum(Cont_mean),2)*100,'%')), vjust=-0.65, hjust=1, size=7.5)+
#  geom_errorbar(aes(x=Age, ymin=Cont_low/sum(Cont_up), ymax=Cont_up/sum(Cont_low)), width=0.3, colour="darkred", alpha=0.9, size=1.2)+ #worst case
  geom_errorbar(aes(x=Age, ymin=Cont_low/sum(Cont_mean), ymax=Cont_up/sum(Cont_mean)), width=0.3, colour="darkred", alpha=0.9, size=1.2)+
      labs(title="", x ="Age classes", y = 'Cumulative contacts %')+
  scale_fill_gradient(low = "gold1", high = 'red2')+scale_y_reverse()+
  theme_minimal()+coord_flip()+geom_point(aes(x=Age,y=Cont_mean/sum(Cont_mean)))+
  scale_x_discrete(name = "Age classes", position = "top") +
  theme(panel.grid.minor = element_blank(),axis.text = element_text(size = 20),
        axis.title.y =  element_text(face = 'bold', size = 20),axis.title.x = element_text(size=20),legend.position="none")

# 
# # #ratio elasticity to contact rate (asym)/elasticity to contact rate (sym)
# # 
# graph <- ggplot(measures_pre_pandemic, aes(x=Age, y=E_asym_mean/E_sym_mean,fill=E_asym_mean/E_sym_mean)) +geom_bar(stat = "identity")
# graph + geom_text(aes(label=round(E_asym_mean/E_sym_mean,2)), vjust=-0.65, hjust=0.001, size=5.5)+
#   geom_errorbar(aes(x=Age, ymin=E_asym_low/E_sym_up, ymax=E_asym_up/E_sym_low), width=0.3, colour="darkred", alpha=0.9, size=1.2)+
#   labs(title="", x ="Age classes", y = expression(paste(tilde(e[i])^'asym','/',tilde(e[i])^'sym')))+
#   scale_fill_gradient(low = "gold1", high = 'red2')+
#   theme_minimal()+coord_flip()+
#   theme(panel.grid.minor = element_blank(),axis.text = element_text(size = 13),
#         axis.title.y =element_text(face = 'bold', size = 13),axis.title.x = element_text(size=13),legend.position="none")
# 
# # # 
# # # #ratios contacts asym vs sym 
# graph <- ggplot(measures_pre_pandemic, aes(x=Age, y=Cont_mean/Cont_sym_mean,fill=Cont_mean/Cont_sym_mean)) +geom_bar(stat = "identity")
# graph + geom_text(aes(label=round(Cont_mean/Cont_sym_mean,2)), vjust=-.65, hjust=1, size=5.5)+
#   geom_errorbar(aes(x=Age, ymin=Cont_low/Cont_sym_up, ymax=Cont_up/Cont_sym_low), width=0.3, colour="darkred", alpha=0.9, size=1.2)+
#   labs(title="", x ="Age classes", y = expression(paste(tilde(c[i])^'asym','/',tilde(c[i])^'sym')))+
#   scale_fill_gradient(low = "gold1", high = 'red2')+
#   theme_minimal()+coord_flip()+scale_y_reverse()+
#   scale_x_discrete(name = "Age classes", position = "top") +
#   theme(panel.grid.minor = element_blank(),axis.text = element_text(size = 13),
#         axis.title.y =element_text(face = 'bold', size = 13),axis.title.x = element_text(size=13),legend.position="none")
# 
# 
