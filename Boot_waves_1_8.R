
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

R_t=T   #set TRUE to evaluate effective reproduction number
#CoMix waves
wrange <-1:8

###Initialize the list containing wave-specific measures of sensitivities, contacts and matrices(NGM)
list_boot_matrices_wave <- vector(mode = "list", length = length(wrange))
list_boot_eigen_wave<- vector(mode = "list", length = length(wrange))


##Bootsraps all locations matrices
nboot <-1000
##Bootsraps location specific matrices
nboot_loc <-1000

cicle <- 1

##Store wild-type q factor----
wild_q_sym <- q_sym

for(wave in wrange){
  
  set.seed(2345)
  
  ##Update Susceptibles----
  S_vec <- S_vec*(!R_t) + SW[,wave]*R_t
  
  list_up_matrices<- list()
  list_low_matrices<- list()
  list_up_eigen<- list()
  list_low_eigen<- list()
  
  
  survey_object_CoMix<-get_survey_object(country=country,
                                         daytype=daytype,
                                         touch=touch,
                                         duration=duration,
                                         gender=gender,
                                         cnt_location=cnt_location,
                                         bool_reciprocal=bool_reciprocal,              ### should the data generate symmetrical matrix?
                                         bool_suppl_professional_cnt=bool_suppl_professional_cnt,   ### should supplementary working contacts be included?
                                         bool_hhmatrix_selection=bool_hhmatrix_selection,      ###select just household members contact
                                         wave=id_waves[wave],                      ### wave to be included (optional)
                                         quiet = TRUE)
  
  
  matrixCoMix<-contact_matrix(survey=survey_object_CoMix, ## First and most important parameter: the survey object
                              age.limits = age_breaks,
                              n=1,                             ## number of bootstraps
                              missing.contact.age="sample",    ## What to do with contacts with missing contact age (possibilities: sample (sample from population), remove (remove them))
                              estimated.contact.age="sample",
                              weigh.age=T,                  ## TRUE if weight on the age is performed
                              weigh.dayofweek=T,            ## TRUE if weight on the day of the week is performed
                              return.demography   = TRUE,
                              weight.threshold = 3             ## maximum weight for one participant
  )
  
  # output_folder <- getwd()
  output_folder="./output_matrices_imputed/"
  age_names <- matrixCoMix$demography$age.group
  di_names <- list(age_names,age_names)
  nbreaks <- length(age_breaks)
  pop_matrix <- matrix(rep(N_vec,nbreaks),ncol=nbreaks,byrow = T)
  check_location <- all(cnt_location==opt_location)
  
  
  
  
  
  list_boot <- NGM_boot(nboot=nboot, popmatrix=pop_matrix,
                        S_vec = S_vec,
                        N_vec=N_vec,
                        reciprocal= bool_reciprocal,
                        check_loc=check_location, 
                        output_folder=output_folder,
                        wave=wave, q_susc=A, q_inf=H, 
                        prop_ratio=prop_ratio,
                        D1=D_1_asym,D2=D_2_sym,
                        age_breaks=age_breaks, R_pcr=T)
  
  
  
  
  #NGM wave specific----
  
  list_matrices <- vector(mode = "list", length =length(list_boot$K))
  list_eigen <- vector(mode = "list", length =length(list_boot$K))
  
  for (iboot in seq_len(length(list_boot$K))) {
    
    K_1 <- list_boot$K[[iboot]]
    R_0_1 <- max(Re(eigen(K_1)$values))
    w_1 <- abs(eigen(K_1)$vectors[,1])
    v_1 <- abs(eigen(t(K_1))$vector[,1])
    w_1 <- w_1/sum(w_1)  
    v_1 <- v_1/sum(v_1*w_1)
    
    
    ####Damping Ratio####
    damping_factor <-  R_0_1/Re(eigen(K_1)$values)[2]
    
    
    # #SENSITIVITY&ELASTICITY MATRICES####
    age_names <- colnames(matrixCoMix$matrix)
    di_names  <- list( age_names, age_names)
    S_built<- matrix(calc.sens(K_1),nbreaks,nbreaks,dimnames=di_names)
    E_built <- matrix(hadamard.prod(S_built,K_1)/R_0_1,nbreaks,nbreaks,dimnames=di_names)
    
    
     list_matrices[[iboot]] <- list(K=list_boot$K[[iboot]],
                                   C_mean_asym=list_boot$C_asym[[iboot]],
                                   C_mean_sym=list_boot$C_sym[[iboot]],
                                   S1=list_boot$Sens[[iboot]],
                                   E1=list_boot$El[[iboot]],
                                   S=S_built,
                                   E=E_built
                                   )
    
    list_eigen[[iboot]] <-  list(R_0=R_0_1,
                                 w=w_1,
                                 v=v_1,
                                 damp_ratio=damping_factor)
    
    
    
  }
  
  
  
  list_boot_matrices_wave[[cicle]] <- list_matrices  
  list_boot_eigen_wave[[cicle]] <- list_eigen
  
  
  cicle <- cicle+1
  
  print(wave)
  
}

##BUILDING PLOT LISTS
plot_scatter_sens <- list()
plot_scatter_el <- list()
plot_bars_el_sens <- list()
plot_bars_contacts <- list()



cicle <- 1
for (surv in wrange[1:8]) {
  #REMEMBER: the positions of the contact matrices do not 
  #correspond to the wave number; we need to readapt as follows:
  this_wave <- which(wrange==surv)
  
  list_contacts <- list()
  list_contacts_f <- list()
  list_sens <- list()
  list_el <- list()
  list_Rt <- list()
  list_ngms <- list()
  
  for (i in seq_len(nboot)) {
    list_contacts[[i]] <- list_boot_matrices_wave[[this_wave]][[i]]$C_mean_asym
    if(this_wave!=8 & this_wave!=length(wrange)){
      list_contacts_f[[i]] <- list_boot_matrices_wave[[this_wave+1]][[i]]$C_mean_asym   
      
    }else{
      #This will generate an identical list to list_contacts, so that for the last wave of the observation 
      #we dont have information on the difference in contacts
      list_contacts_f[[i]] <- list_contacts[[i]]
      
    }
    list_sens[[i]] <- list_boot_matrices_wave[[this_wave]][[i]]$S1
    list_el[[i]] <- list_boot_matrices_wave[[this_wave]][[i]]$E1
    list_ngms[[i]] <- list_boot_matrices_wave[[this_wave]][[i]]$K
    list_Rt[[i]] <- list_boot_eigen_wave[[this_wave]][[i]]$R_0
  }
  
 
  S<- apply(simplify2array(list_sens),1:2,mean)
  E<- apply(simplify2array(list_el),1:2,mean)
  NGM <- apply(simplify2array(list_ngms),1:2,mean)

  Rt <-R_series_mean[wrange[this_wave]]
  Rt_up <- R_series_up[wrange[this_wave]]
  Rt_low <- R_series_low[wrange[this_wave]]
  
  #CONTACTS
  list_c1_sum <- lapply(list_contacts,colSums)
  list_c2_sum <- lapply(list_contacts_f, colSums)
  
  list_diff <-   lapply(seq_along(list_contacts),
                        function(i) list_c2_sum[[i]]-list_c1_sum[[i]])
  list_diff_perc <- lapply(seq_along(list_contacts),
                           function(i) list_diff[[i]]/list_c1_sum[[i]])
  list_ratio_cont <- lapply(seq_along(list_contacts),
                            function(i) list_c2_sum[[i]]/list_c1_sum[[i]])
  
  mean_C<-apply(do.call(rbind, list_c1_sum),2,mean)
  Q_C <-apply(do.call(rbind, list_c1_sum), 2, quantile, probs = c(0.025, 0.975))
  
  mean_diff<-apply(do.call(rbind, list_diff),2,mean)
  mean_diff_perc <-apply(do.call(rbind, list_diff_perc),2,mean)
  mean_ratio_cont <- apply(do.call(rbind, list_ratio_cont),2,mean)
  
  Q_diff <-apply(do.call(rbind, list_diff), 2, quantile, probs = c(0.025, 0.975))
  Q_diff_perc<-apply(do.call(rbind, list_diff_perc), 2, quantile, probs = c(0.025, 0.975))
  Q_ratio <-apply(do.call(rbind, list_ratio_cont), 2, quantile, probs = c(0.025, 0.975))
  
  #SENSITIVITY
  list_ngm_sum<-lapply(list_ngms, colSums)
  list_sens_sum<-lapply(list_sens, colSums)
  list_el_sum<-  lapply(list_el, colSums)
  
  mean_ngm_sum<-apply(do.call(rbind, list_ngm_sum),2,mean)
  Q_K<-apply(do.call(rbind, list_ngm_sum), 2, quantile, probs = c(0.025, 0.975))
  
  mean_sens_sum<-apply(do.call(rbind, list_sens_sum),2,mean)
  Q_S<-apply(do.call(rbind, list_sens_sum), 2, quantile, probs = c(0.025, 0.975))
  mean_el_sum<-apply(do.call(rbind, list_el_sum),2,mean)
  Q_E<-apply(do.call(rbind, list_el_sum), 2, quantile, probs = c(0.025, 0.975))
  
  
  
  
  
  graph_df<-  data.frame(Wave=integer(),
                         Age=character(),
                         K_mean=double(),
                         Cont=double(),
                         Diff_perc=double(),
                         Diff_abs=double(),
                         Sens=double(),
                         El=double(),
                         R_t=double(),
                         
                         K_up=double(),
                         K_low=double(),
                         C_up=double(),
                         C_low=double(),
                         Diff_perc_up=double(),
                         Diff_perc_low=double(),
                         Diff_abs_up=double(),
                         Diff_abs_low=double(),
                         Sens_up=double(),
                         Sens_low=double(),
                         El_up=double(),
                         El_low=double(),
                         R_t_up=double(),
                         R_t_low=double(),
                         # 
                         stringsAsFactors=F)
  
  for (k in 1:length(age_breaks)) {
    new_record <- list(as.integer(wrange[this_wave]),
                       age_names[k],
                       #mean values
                       mean_ngm_sum[k],
                       mean_C[k],
                       mean_diff_perc[k],
                       mean_diff[k],
                       mean_sens_sum[k],
                       mean_el_sum[k],
                       round(Rt,2),
                       #quantiles 
                       Q_K[2,][k],
                       Q_K[2,][k],
                       
                       Q_C[2,][k],
                       Q_C[1,][k],
                       Q_diff_perc[2,][k],
                       Q_diff_perc[1,][k],
                       Q_diff[2,][k],
                       Q_diff[1,][k],
                       Q_S[2,][k],
                       Q_S[1,][k],
                       Q_E[2,][k],
                       Q_E[1,][k],
                       round(Rt_up,2),
                       round(Rt_low,2)
                       #  
                       
    )
    
    graph_df[nrow(graph_df) + 1,] = new_record
  }
  graph_df$Age <- factor(graph_df$Age, levels=age_names)
  
  my_colors <- c( "Elasticity" = "#000033","Sensitivity" = "#F9FF00")
  col_blind_pal <- c( "#F0E442","#D55E00","#E69F00", "#117733","#44AA99","#0072B2", "#56B4E9", "#CC6677", "#882255") 
  
  plot_bars_contacts[[cicle]] <- ggplot(graph_df, aes(x=Age, y=Cont)) +geom_bar(stat = "identity", fill="azure3")+
    geom_errorbar(aes(x=Age, ymin=C_low, ymax=C_up), width=0.4, colour="black", alpha=0.9, linewidth=0.8)+
    labs(title=paste0("Wave ",graph_df[1,]$Wave,' - Rt= ',graph_df[1,]$R_t), x ="", y = "")+
    geom_point(aes(x=Age, y=Cont))+
    theme(axis.text.x=element_text(size=28, angle=90, vjust=0.3),
          axis.text.y=element_text(size=28),
          legend.key.size =  unit(1,'cm'),
          plot.title=element_text(size=30),
          axis.title.x = element_text(size=17,face ='plain',vjust=-3),
          legend.title = element_text(size=30, ),  # Title
          legend.text  = element_text(size=28, family ='sans'), 
          plot.margin = unit(c(1,1,1,1), "cm"))
  
  
  df1 <- data.frame(Age=age_names,Value=mean_sens_sum/sum(mean_sens_sum),UP=Q_S[2,]/sum(mean_sens_sum),LOW=Q_S[1,]/sum(mean_sens_sum),Measure='Sensitivity')
  df2 <- data.frame(Age=age_names,Value=mean_el_sum,UP=Q_E[2,],LOW=Q_E[1,],Measure='Elasticity')
  df_graph <- rbind(df1,df2)
  df_graph$Age <- factor(df_graph$Age, levels=age_names)
  
  
  # Set a fixed size for the plots
  plot_size <- theme(plot.title = element_text(size = 20),
                     axis.title = element_text(size = 20),
                     axis.text = element_text(size = 16),
                     legend.title = element_text(size=20),
                     legend.text = element_text(size=20)
  )
  
  #Set the margins for the axis labels and titles
  plot_margin <- theme(plot.margin = unit(c(1.2, 0.7, 0.7, 0.7), "lines"),
                       axis.title.y = element_text(margin = margin(r = 10)),
                       axis.title.x = element_text(margin = margin(t = 10)))
  
  plot_title <- paste0("Wave ",graph_df[1,]$Wave,' (',w_dates[wrange[this_wave]],
                       ')    Rt=',graph_df$R_t)
  
  plot_bars_el_sens[[cicle]] <- ggplot(df_graph, aes(x=as.factor(Age), y=Value, fill=Measure)) +
    scale_size_continuous(range=c(1, 5)) +
    scale_fill_manual(values = c("#767676","#C8C8C8")) +
    theme(legend.text = element_text(size=20))+
    geom_bar(position=position_dodge(), stat="identity", colour='black') +
    geom_errorbar(aes(ymin=LOW, ymax=UP), width=.2,position=position_dodge(.9))+
    labs(title=plot_title, fill='', x='Age', y='')+
    plot_size 
  
  
  plot_scatter_el[[cicle]] <- ggplot(graph_df, aes(x=Cont, y=El, color=Age)) +
    geom_point(size=8)+ 
    plot_size +
    plot_margin+
    scale_size_continuous(range=c(1, 5)) +
    scale_color_manual(values = col_blind_pal) +
    labs(y=bquote('Elasticity'~tilde(e)[j]),x='Contacts', 
         title = plot_title )
  
  
  plot_scatter_sens[[cicle]]<- ggplot(graph_df, aes(x=Diff_abs, y=Sens, color=Age)) +
    geom_point(size=8)+
    plot_size +
    plot_margin+
    scale_size_continuous(range=c(1, 5)) +
    scale_color_manual(values = col_blind_pal) +
    labs(y=bquote('Sensitivity'~tilde(s)[j]), x='Contacts variation', title=plot_title )
  
  
  cicle <- cicle+1  
}

##Figure 3a----
tiff(filename = "Figure3a.tiff", width = 16, height = 11, units = 'in', res = 500)
ggarrange(plotlist=plot_bars_el_sens[1:4],nrow=2,ncol=2,common.legend = TRUE, legend="bottom")
dev.off()

p3a <-ggarrange(plotlist=plot_bars_el_sens[1:4],nrow=2,ncol=2,common.legend = TRUE, legend="bottom")
ggsave("ElSens_error1to4.png", plot = p3a, width = 15, height = 11, dpi = 500)


##Figure 3b----
tiff(filename = "Figure3b.tiff", width = 16, height = 11, units = 'in', res = 500)
ggarrange(plotlist=plot_scatter_el[1:4],nrow=2,ncol=2, common.legend = TRUE, legend="bottom")
dev.off()

p3b <-ggarrange(plotlist=plot_scatter_el[1:4],nrow=2,ncol=2, common.legend = TRUE, legend="bottom")
ggsave("El_vs_Cont1to4.png", plot = p3b, width = 15, height = 11, dpi = 500)

