##CLASSES NEEDED
setClass("MatrixProp", slots=list(
  determinant="numeric", 
  primitive="logical", 
  irreducible="logical",
  ev="eigen",
  ev_left="eigen",
  R_0='numeric'))



###### functions needed#####
#to move in another file to call at the beginning of the script
get_percapita_matrix <- function (cnt_matrix, pop_vec, byrow=F){
  res <-matrix(0,nrow(cnt_matrix),ncol(cnt_matrix))
  #Under the hypothesis of NGM=SABC set:
  #byrow=False when the generic element of the contact matrix m_ij 
  #represents the mean number of contacts per unit time of an individual from group i
  #with individuals of group j (ComiX Case)
  #byrow=True when m_ij is intended other way around (the contact is made infected side)
  for(i in 1:nrow(cnt_matrix)){
    for (j in 1:ncol(cnt_matrix)) {
      
      if(byrow==TRUE){
        res[i,j] <- cnt_matrix[i,j]/pop_vec[i]
      }
      else{
        res[i,j] <- cnt_matrix[i,j]/pop_vec[j] 
      }
    }
    
  }
  return(res)
}

built_ngm <- function(S,qsusc,qinf,C_asym,C_sym,D1,D2,q){
  q_a <- prop_ratio*q
  Asym <- S%*%qsusc%*%C_asym%*%(q_a*qinf)
  Sym <- S%*%qsusc%*%C_sym%*%(q*qinf)
  Kbuilt <- Asym%*%D1+Sym%*%D2
  return(list(Asymp=Asym, Symp=Sym, NGM=Kbuilt))
  
}
symmetrize <- function(A,N=NULL){
  res <- matrix(0,nrow(A),ncol(A))
  if(is.null(N)){
    s <- is.symmetric.matrix(A)
    res <- A*s + (!s)*(A+t(A))/2
  }else{
    
    for (i in 1:nrow(A)) {
      for (j in 1:ncol(A)) {
        res[i,j] <- (A[i,j]*N[i]+A[j,i]*N[j])/(2*N[i])
      }
      
    }
    
  }
  return(res)
}



mat_obj <- function(A){
  obj <- new("MatrixProp")
  obj@determinant=det(A)
  #obj@primitive=isPrimitive(A) 
  #obj@irreducible=isIrreducible(A)
  obj@ev= eigen(A)
  obj@ev_left<- eigen(t(A))
  return(obj)
}

set_epi_params <- function(case) {
  
  age_breaks <- switch(case,
                       c(0,10,20,30,40,50,60,70),
                       c(0,10,20,30,40,50,60,70,80,90),
                       c(0,6,12,18,30,40,50,60,70)
  )
  Asusc_vec <- switch(case, #age-specific q-susceptibility
                      c(0.4,0.38,0.79,0.86,0.8,0.82,0.88,0.74),
                      c(0.4,0.38,0.79,0.86,0.8,0.82,0.88,0.74,0.74,0.74),
                      c(0.4,0.39,0.38,0.79,0.86,0.8,0.82,0.88,0.74)
  )
  
  Hinf_vec <- switch(case, #age-specific q-infectiousness
                     c(0.55,0.55,0.59,0.7,0.76,0.9,0.9,0.9),
                     c(0.55,0.55,0.59,0.7,0.76,0.9,0.9,0.9,0.9,0.9),
                     c(0.54,0.55,0.56,0.59,0.7,0.76,0.9,0.99,0.99)
  )
  
  p <- switch(case,   #age specific probabilty of asymptomatic infections
              c(0.94,0.90,0.84,0.61,0.49,0.21,0.02,0.02),
              c(0.94,0.90,0.84,0.61,0.49,0.21,0.02,0.02,0.2,0.2),
              #c(0.9,0.9,0.84,0.61,0.49,0.21,0.02,0.02,0.02)
              c(0.94,0.92,0.90,0.84,0.61,0.49,0.21,0.02,0.02))
  
  phi_0 <- switch(case, #Age-specific probability of deveeloping only mild symptoms
                  c(0.972,0.992,0.984,0.987,0.977,0.971,0.958,0.936),
                  c(0.972,0.992,0.984,0.987,0.977,0.971,0.958,0.926,0.956,0.926),
                  c(0.972,0.982,0.992,0.984,0.987,0.977,0.971,0.958,0.936))
  
  omega_1 <- switch(case, #Severe symptomatic removal rate
                    c(0.167,0.095,0.099,0.162,0.338,0.275,0.343,0.338),
                    c(0.167,0.095,0.099,0.162,0.338,0.275,0.343,0.378,0.334,0.302),
                    c(0.167,0.131,0.095,0.099,0.162,0.338,0.275,0.343,0.338))
  p_mask <- switch(case,
                   c(0.37,0.37,0.37,0.41,0.41,0.41,0.57,0.57),
                   c(0.37,0.37,0.37,0.41,0.41,0.41,0.57,0.57,0.57,0.57),
                   c(0.37,0.37,0.37,0.37,0.41,0.41,0.41,0.57,0.57))
  
  # SW <- switch(case, #per wave susceptible number
  #              SW_70,
  #              SW_70,
  #              round(SW_school))
  # 
  params <- list(Age=age_breaks, A=Asusc_vec, H=Hinf_vec, p_asym=p,phi_0=phi_0, omega_1=omega_1, p_mask=p_mask)
  
  ifelse(any(sapply(params, is.null)) || (anyNA(params)),stop("An error occurred: check the value of 'case'"), print("Epi-parameter inizialized")) 
  
  return(params)
}


adapt_susc <- function(M,age_groups, dates, proportions){
  
  R <- matrix(NA, length(age_groups),length(dates))
  for(i in 1:ncol(M)){
    #we just assume the susceptibles are distributed across the
    #age classes of the youngest following the demographic distribution
    R[1,i] <- M[1,i]*proportions[1]
    R[2,i] <- M[1,i]*proportions[2]+M[2,i]*proportions[3]
    R[3,i] <- M[2,i]*proportions[4]
    R[4,i] <- M[2,i]*proportions[5]+M[3,i]
    R[5:9,i] <- M[4:nrow(M),i]
  }
  return(round(R))
}

plot_matrix <- function(mij,title= '',axisnames=NULL, xlab=NULL,ylab=NULL){
  if(is.null(axisnames)==FALSE){
    colnames(mij) <- axisnames
  }
  
  if(all(is.na(mij))){
    return(NA)
  }
  redc <- rev(heat.colors(100))
  par(mar=c(5, 6, 2, 2),mgp=c(3,0.5,0))
  p <- simage(s = mij, 
              xlab=ifelse(is.null(xlab),"Age of participant (year)",xlab),
              ylab=ifelse(is.null(xlab),"Age of contact (year)",ylab),
              legend.width=1,
              slim=c(min(mij,na.rm=T), max(mij,na.rm=T)), 
              cex.lab=1.2,
              cex.main=1.2, 
              las=0.1,
              col=redc, 
              main=title, 
              xaxt="n", 
              yaxt="n")
  # set axis 
  plt_ticks <- seq(0,1,length=nrow(mij))
  axis(2, at=plt_ticks, labels = c(colnames(mij)),cex.axis=0.9,tick = FALSE,las=1)
  axis(1, at=plt_ticks, labels = c(colnames(mij)),cex.axis=0.9,tick = FALSE)
  
  # format results (rounding/scientific)
  if(any(max(mij,na.rm=T)>1)){
    mij <- round(mij,digits=format_num_digits)
  } else{
    mij <- format(mij,digits = format_num_digits)
  }
  # get grid centers and add value
  e_grid <- expand.grid(plt_ticks,plt_ticks)
  text(e_grid, labels = mij)
}


non_neg <- function(M=matrix()){
  
    if(any(M<0)){
    M[which(M<0,arr.ind = T)] <- 0.1
    }
  return(M)
}


deal_nas <- function(A,B,perc=1){
  #A and B are two objects 'contact_matrix' (lists returned by contact_matrix function)
  #A contact matrix with missing partecipants in certain age groups
  #B pre_pandemic contact matrix from which one can get proxies of contact for missing ages
  #perc a coefficient [0,1] to tune the imputed contact rates
  missing_pop <- NULL
  for(j in 1:nrow(A$matrix)){
    if (all(is.na(A$matrix[j,]))){
      A$matrix[j,] <- B$matrix[j,]*perc
      missing_pop <- c(missing_pop,B$demography$population[j])
    }
  }
  ret_list <- list(A$matrix,c(missing_pop,A$demography$population))
  names(ret_list) <- c('contact_matrix','population')
  return(ret_list)
}

plot_bars <- function(mij,title= '', extratitle=NULL, lab_y='', col="lightcyan", adj_ylim=1,
                      adj_labels=0,setylim=F,yaxis,barnames=NULL,change_xlab=F,x_lab=NULL){
  if(all(is.na(mij))){
    return(NA)
  }
  if(all(mij>=0)){
    lim_y=c(0,max(mij))*adj_ylim
    
  }else{
    tunen <- min(mij)< -max(abs(mij))/2
    tunep <- max(mij)>  max(abs(mij))/2
    lim_y=c(min(mij)*(tunen)-max(abs(mij))/2*(!tunen),max(mij)*(tunep)+max(abs(mij))/2*(!tunep))*adj_ylim
    
  }
  y_lim <- NULL
  ifelse(setylim==T,y_lim <- yaxis, y_lim <-lim_y) 
  if(!is.null(barnames)){
    names(mij) <- barnames
  }
  if(change_xlab==F){
    x <- barplot(mij,main =  paste(title,as.character(extratitle),sep = ' '),
                 ylim =y_lim,
                 xlab = paste('Age groups-partecipants'),
                 ylab =lab_y, 
                 col  = col)
    lab_pos <- mij
  }
  else{
    x <- barplot(mij,main =  paste(title,as.character(extratitle),sep = ' '),
                 ylim =y_lim,
                 xlab =x_lab,
                 ylab =lab_y, 
                 col  = col)
    lab_pos <- mij
  }
 
  for(i in 1:length(mij)){
    lab_pos[i] <- mij[i]+sign(mij[i])*adj_labels
  }
  text(x,lab_pos, labels = as.character(signif(mij, digits = 2)),pos = 3)
  
}
build_perturbation_matrices <- function(nrows,ncols=NULL,w,v,dimnames,KL,R_0){
  #assemble matrices
  if(is.null(ncols)){ncols <- nrows}
  S <- matrix(0,ncol = ncol(KL),nrow=nrow(KL))
  E <- matrix(0, ncol = ncol(KL),nrow=nrow(KL))
  for (i in 1:nrows){
    for(j in 1:ncols){
      S[i,j] <-v[i]*w[j]
      E[i,j] <-v[i]*w[j]*KL[i,j]/R_0
      # S[i,j] <-v[i+(position*nrows)]*w[j]
      # E[i,j] <-S[i,j]*KL[i,j+(position*nrows)]/R_0
    }
  }
  return(list(S=S,E=E))
}
calc.sens <- function(G, normalize=TRUE){
  ev <- eigen(G)
  lmax <- which(Re(ev$values)==max(Re(ev$values)))
  U <- ev$vectors
  u <- abs(Re(U[,lmax]))
  ev_left <- eigen(t(G))
  lmax_l <- which(Re(ev_left$values)==max(Re(ev_left$values)))
  V <- ev_left$vectors
  v <- abs(Re(V[,lmax_l]))
  if(normalize){
    u <- u/sum(u)  
    v <- v/(as.numeric(v%*%u))
  }
  
  s <- as.matrix(v%o%u)
  return(s)
}
calc.el <- function(G, normalize=TRUE){
  M <- calc.sens(G,normalize)
  l <- Re(eigen(G)$values[1])
  E <- M*G
  E <- E/l
  return(E)
}

#the function below takes values from the global environment and it is assumend 
# to be called within the execution of the script "wave_specific_an.R"
lower_sens_matrices <- function(K, #NGM
                                M_asym, #Contact Matrix Asym
                                M_sym,  #Contact Matrix Sympt
                                Susc,   #Susceptibles
                                q_factor #q factor for symptomatic cases
                                ){     
  
  #Initialize matrices
  R_t <-Re(eigen(K)$values[1]) 
  Sens <- calc.sens(K)
  w <- top_eigenvectors(K)$right
  v <- top_eigenvectors(K)$left
  
  #Asymp and Sympt dependent parts of K (next gen matrix)
  K_asym<- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  K_sym<- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  
  #Lower Level Seinsitivities: we consider just 
  A_pert<- matrix(0,nbreaks,nbreaks,dimnames=di_names)       #
  Sens_A<- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  H_pert<- matrix(0,nbreaks,nbreaks,dimnames=di_names)       
  Sens_H<- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  
  #Asympt vs Sympt sensitivities & elasticities (linear terms)
  C_pert_asym<- matrix(0,nbreaks,nbreaks,dimnames=di_names)       # initialize matrix of perturbations of K entries w.r.t. contact rates variations
  Sens_C_asym<- matrix(0,nbreaks,nbreaks,dimnames=di_names)      # initialize sensitivity matrix to contact rates variations
  C_pert_sym<- matrix(0,nbreaks,nbreaks,dimnames=di_names)       # initialize matrix of perturbations of K entries w.r.t. contact rates variations
  Sens_C_sym<- matrix(0,nbreaks,nbreaks,dimnames=di_names)      # initialize sensitivity matrix to contact rates variations
  El_C_asym <- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  El_C_sym <- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  
  
  #Sensitivity to transition rates
  #Theta
  Theta_pert <- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  Sen_Theta  <- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  El_Theta   <- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  
  
  #Delta1
  S1_pert <- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  Sen_S1 <- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  El_S1  <- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  
  #Omega
  O_pert <- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  Sen_O  <- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  El_O   <- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  
  #Psi
  Psi_pert <- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  Sen_Psi  <- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  El_Psi   <- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  
  #Delta2
  D2_pert <- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  Sen_D2  <- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  El_D2  <- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  
  #Pob being asymptomatic
  Prob_pert <- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  Sen_Prob  <- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  El_Prob   <- matrix(0,nbreaks,nbreaks,dimnames=di_names)
  
  #Sensitivity ratio
  Sens_ratio <-matrix(0,nbreaks,nbreaks,dimnames=di_names)  
  
  ##Compute measures of sensitvity##
  for (i in 1:nbreaks){
    for(j in 1:nbreaks){
      
      #Store components of pre& non symptomatic cases and that of symptomatic cases
      K_asym[i,j] <- Susc[i]*q_factor*prop_ratio*M_asym[i,j]/N_vec[j]*A[i,i]*H[j,j]*(sigma1+p[j]*theta)/(sigma1*theta)
      K_sym[i,j]  <- Susc[i]*q_factor*M_sym[i,j]/N_vec[j]*A[i,i]*H[j,j]*(Psi[j]+Omega[j,j])*(1-p[j])/(Omega[j,j]*(Psi[j]+Sigma[j]))
      #sensitivity to q-infectiousness H[i]
      H_pert[i,j] <- K[i,j]/H[j,j]
      Sens_H[i,j] <- Sens[i,j]*H_pert[i,j]
      
      #sensitivity to q-susceptibility A[i]
      A_pert[i,j] <- K[i,j]/A[i,i]
      Sens_A[i,j] <- Sens[i,j]*A_pert[i,j]
      
      #sensitivity to cont rate (asym) - generic perturbation in transmission terms of the asymptomatic
      C_pert_asym[i,j]<-K_asym[i,j]/M_asym[i,j]
      Sens_C_asym[i,j]<-Sens[i,j]*C_pert_asym[i,j] 
      El_C_asym[i,j] <-Sens_C_asym[i,j]*M_asym[i,j]/R_t
      
      #sensitivity to cont rate (sym) - generic perturbation in transmission terms of the symptomatic
      C_pert_sym[i,j]<- K_sym[i,j]/M_sym[i,j]
      Sens_C_sym[i,j]<- Sens[i,j]*C_pert_sym[i,j]    
      El_C_sym[i,j] <-Sens_C_sym[i,j]*M_sym[i,j]/R_t
      
      #store asym and sym transmission terms
      T_ij_asym <- Susc[i]*q_factor*prop_ratio*M_asym[i,j]/N_vec[j]*A[i,i]*H[j,j]
      T_ij_sym  <- Susc[i]*q_factor*M_sym[i,j]/N_vec[j]*A[i,i]*H[j,j]
      
      #Sensitivity to removal rate of pre-symptomatic (Theta)
      Theta_pert[i,j] <- -1/(theta^2)*T_ij_asym
      Sen_Theta[i,j] <- Sens[i,j]*Theta_pert[i,j]
      El_Theta[i,j] <- Sen_Theta[i,j]*theta/R_t
      
      #To recovery rate of asymptomatic (Sigma1)
      S1_pert[i,j] <- -p[j]/(sigma1^2)*T_ij_asym
      Sen_S1[i,j] <- Sens[i,j]*S1_pert[i,j]
      El_S1[i,j] <- Sen_S1[i,j]*sigma1/R_t
      
      #To recovery of mild cases (delta2)
      D2_pert[i,j] <- -T_ij_sym*(Psi[j]+Omega[j,j])*(1-p[j])/(Omega[j,j]*(Psi[j]+Sigma[j])^2)
      Sen_D2[i,j]  <- Sens[i,j]*D2_pert[i,j]
      El_D2[i,j]  <- Sen_D2[i,j]*Sigma[j]/R_t
      
      #To transition of mild cases to severe (Psi)
      Psi_pert[i,j] <- T_ij_sym*(Sigma[j]-Omega[j,j])*(1-p[j])/(Omega[j,j]*(Psi[j]+Sigma[j])^2)
      Sen_Psi[i,j]  <- Sens[i,j]*Psi_pert[i,j]
      El_Psi[i,j]  <- Sen_Psi[i,j]*Psi[j]/R_t
      
      #To removal severe (omega)
      O_pert[i,j] <- -T_ij_sym*Psi[j]*(1-p[j])/(Omega[j,j]^2*(Psi[j]+Sigma[j]))
      Sen_O[i,j]  <- Sens[i,j]*O_pert[i,j]
      El_O[i,j]  <- Sen_O[i,j]*Omega[j,j]/R_t
      
      #To prob staying asymp
      Prob_pert[i,j] <- T_ij_asym/sigma1 - T_ij_sym*(Psi[j]+Omega[j,j])/(Omega[j,j]*(Psi[j]+Sigma[j]))
      Sen_Prob[i,j] <- Sens[i,j]*Prob_pert[i,j]
      El_Prob[i,j] <- Sen_Prob[i,j]*p[j]/R_t
      #sens. ratio
      Sens_ratio[i,j] <- v[i]/v[j]
      
    }
  }
  ret <- list(S_H=Sens_H,
              S_A=Sens_A,
              S_C_asym=Sens_C_asym,
              S_C_sym= Sens_C_sym,
              S_T=Sen_Theta,
              S_S1=Sen_S1,
              S_D2=Sen_D2,
              S_Psi=Sen_Psi,
              S_O=Sen_O,
              S_P=Sen_Prob,
              Sens_ratio=Sens_ratio,
              E_asym=El_C_asym,
              E_sym=El_C_sym,
              E_T=El_Theta,
              E_S1=El_S1,
              E_D2=El_D2,
              E_Psi=El_Psi,
              E_O=El_O,
              E_P=El_Prob
              )
  
  return(ret)
  
}




top_eigenvectors <- function(G,scale=TRUE){
  ev <- eigen(G)
  lmax <- which(Re(ev$values)==max(Re(ev$values)))
  U <- ev$vectors
  u <- abs(Re(U[,lmax]))
  ev_left <- eigen(t(G))
  lmax_l <- which(Re(ev_left$values)==max(Re(ev_left$values)))
  V <- ev_left$vectors
  v <- abs(Re(V[,lmax_l]))
  if(scale==T){
    u <- u/sum(u)  
    v <- v/sum(v*u)
  }
  
  return(list(right=u,left=v))
}

gen_w_dates <- function(){
  
  waves <- data.frame(DATE=character(0),
                      id=integer(0),
                      stringsAsFactors=FALSE)
  
  
  survey_object_CoMix<-get_survey_object(country="Belgium 2020 CoMix (Coletti 2020)",
                                         daytype="All contacts",
                                         touch="All contacts",
                                         duration="All contacts",
                                         gender="All",
                                         cnt_location=opt_location,
                                         bool_reciprocal=FALSE,              ### should the data generate symmetrical matrix?
                                         bool_suppl_professional_cnt=TRUE,   ### should supplementary working contacts be included?
                                         bool_hhmatrix_selection=FALSE,
                                         wave='All',                      ### wave to be included  
                                         quiet = TRUE)
  
  id_waves<-gtools::mixedsort(unique(survey_object_CoMix$participants$wave))
  
  
  
  for (w in id_waves) {
    
    survey_object_CoMix<-get_survey_object(country="Belgium 2020 CoMix (Coletti 2020)",
                                           daytype="All contacts",
                                           touch="All contacts",
                                           duration="All contacts",
                                           gender="All",
                                           cnt_location=opt_location,
                                           bool_reciprocal=FALSE,              ### should the data generate symmetrical matrix?
                                           bool_suppl_professional_cnt=TRUE,   ### should supplementary working contacts be included?
                                           bool_hhmatrix_selection=FALSE,
                                           wave=w,                      ### wave to be included  
                                           quiet = TRUE)
    
    newrows<-data.frame(table(survey_object_CoMix[["participants"]][["sday_id"]])) %>%
      filter(Freq != 0) %>%
      mutate(DATE=format(as.Date(as.character(Var1),format="%Y.%m.%d"),"%d/%m/%Y"))%>%
      select(DATE,Freq)
    newrows$Id <-as.integer(strsplit(w,':')[[1]][1])
    
    waves <- rbind(waves,newrows)
  }
  return(waves)
}


tol_sign <- function(a,b,tol){
  ret <- 0
  if((a-b)>= tol){
    ret <- 1
  }else if((a-b)<=(-tol)){
    ret <- -1
  }
  return(ret)
  
}

w_row_str <- function(M,tolerance){
  ret <- vector(mode = "list", length = nrow(M))
  for (i in seq_len(nrow(M))) {
    ret[[i]] <- matrix(0,ncol(M),ncol(M))
    for (j in seq_len(ncol(M))) {
      for (k in seq_len(ncol(M))) {
        if(j!=k & k!=i){
          ret[[i]][j,k] <- tol_sign(a=M[i,k],b=M[i,j],tol=tolerance)
          
        }
      }
    }
    
  }
  return(ret)
}

w_row_gradient <- function(list_a,list_b){
  ret <-0 
  num <- 0
  den <- 0
  for (i in seq_len(length(list_a))){
    num <- num+sum(hadamard.prod(list_a[[i]],list_b[[i]]))
    den <- den+sum(abs(hadamard.prod(list_a[[i]],list_b[[i]])))
    
  }
  ret <- ret+(num/den)
  return(ret)
}


#takes as input the projection matrix and other parameters and 
#returns a otput the vector with proportions of infected in each age-class
#The initial condition is randomly generated
project_gen <- function(Op,niter,low_bound,up_bound,replace=T){
  I_0 <-  sample(low_bound:up_bound, size = nrow(Op), replace = replace)
  I_n <- rep(0,nrow(Op))
  n_iter <- niter
  I_series <- vector(mode='list',length = n_iter+1)
  I_series[[1]] <-as.array(I_0/sum(I_0))
  
  for (i  in 2:(n_iter+1)) {
    I_n <- Op%*%I_0
    I_0 <- I_n
    prop_n <- I_n/sum(I_n)
    I_series[[i]] <- prop_n
  }
  return(I_series)
 
}

##
err_conv <- function(M,w,niter,N,replace=T,bounds){
  errors_iter <- vector(mode='list',length = N)
  for (j in 1:N) {
    errors_iter[[j]] <- rep(0,niter+1)
    I_series <- project_gen(M, niter,
                            low_bound = bounds[1], 
                            up_bound =bounds[2],replace)
    for (i in 1:(n_iter+1)) {
      errors_iter[[j]][i] <- sum(abs(I_series[[i]]-w))/2
    }
    
  }
  # convert the list to a matrix
  error_matrix <- do.call(rbind, errors_iter)
  # calculate the matrix columns mean
  err_mean <- apply(error_matrix, 2, mean)
  err_sd <- apply(error_matrix, 2, sd)
  err_low<- apply(error_matrix, 2, quantile,prob=0.025)
  err_up<- apply(error_matrix, 2, quantile,prob=0.975)
  
  return(list(e_mean=err_mean,e_sd=err_sd, e_low=err_low,e_up=err_up) )
}

err_conv_D <-function(M,niter,N,replace=T,bounds){
  # S_diff <- vector(mode='list',length = N)
  # R_diff <- vector(mode='list',length = N)
  
  l <- Re(eigen(M)$values[1])
  w <- abs(Re(eigen(M)$vectors[,1]))
  v <- abs(Re(eigen(t(M))$vectors[,1]))
  w <- w/sum(w)  
  v <- v/(as.numeric(v%*%w))

  errors_D1 <- vector(mode='list',length = N)
  errors_D2 <- vector(mode='list',length = N)
  
  for (j in 1:N) {
    errors_D1[[j]] <- rep(0,niter+1)
    errors_D2[[j]] <- rep(0,niter+1)
    I_series <- vector(mode='list',length = niter+1)
    S_diff<- vector(mode='list',length = niter+1)
    R_diff<- vector(mode='list',length = niter+1)
    
    I_0 <-  sample(bounds[1]:bounds[2], size = nrow(M), replace = replace)
    #I_0 <- I_0/sum(I_0)
    c1 <- sum(v*I_0)
    I_n <- rep(0,nrow(M))
   
    I_series[[1]] <-as.array(I_0)
    S_diff[[1]] <- I_0-(c1*w)
    R_diff[[1]] <- abs(I_0-(c1*w))
    errors_D1[[j]][1] <- sum(abs(S_diff[[1]]))
    errors_D2[[j]][1] <- sum(R_diff[[1]])
   
      
    for (i  in 2:(niter+1)) {
      I_n <- M%*%I_0
      I_0 <- I_n
      #prop_n <- I_n/sum(I_n)
      I_series[[i]] <-I_n
      S_diff[[i]] <- S_diff[[i-1]]+(I_series[[i]]/(l^i))-(c1*w)
      R_diff[[i]] <- R_diff[[i-1]]+abs((I_series[[i]]/(l^i))-(c1*w))
      
      errors_D1[[j]][i] <- sum(abs(S_diff[[i]]))
      errors_D2[[j]][i] <- sum(R_diff[[i]])
    }
  }
  # convert the list to a matrix
  error_matrix_D1 <- do.call(rbind, errors_D1)
  error_matrix_D2 <- do.call(rbind, errors_D2)
  # calculate the matrix columns mean
  err_mean1 <- apply(error_matrix_D1, 2, mean)
  err_sd1 <- apply(error_matrix_D1, 2, sd)
  err_low1<- apply(error_matrix_D1, 2, quantile,prob=0.025)
  err_up1<- apply(error_matrix_D1, 2, quantile,prob=0.975)
  
  err_mean2 <- apply(error_matrix_D2, 2, mean)
  err_sd2 <- apply(error_matrix_D2, 2, sd)
  err_low2<- apply(error_matrix_D2, 2, quantile,prob=0.025)
  err_up2<- apply(error_matrix_D2, 2, quantile,prob=0.975)
  
  return(list(e_mean1=err_mean1, e_sd1=err_sd1, e_low1=err_low1, e_up1=err_up1,
              e_mean2=err_mean2, e_sd2=err_sd2, e_low2=err_low2, e_up2=err_up2) )
}





Boot_baseline <- function(nboot){
  
  
  
  
  set.seed(2345)
  survey_object_CoMix<-get_survey_object(country="Belgium 2010* (Van Hoang 2020)",
                                         daytype=daytype,
                                         touch=touch,
                                         duration=duration,
                                         gender=gender,
                                         cnt_location=cnt_location,
                                         bool_reciprocal=bool_reciprocal,              ### should the data generate symmetrical matrix?
                                         bool_suppl_professional_cnt=bool_suppl_professional_cnt,   ### should supplementary working contacts be included?
                                         bool_hhmatrix_selection=bool_hhmatrix_selection,      ###select just household members contact
                                         wave="All",                      ### wave to be included (optional)
                                         quiet =T)
  #### The following code generates the matrix
  matrixCoMix_start<-contact_matrix(survey=survey_object_CoMix, ## First and most important parameter: the survey object
                                    age.limits = age_breaks,
                                    n=nboot,                             ## number of bootstraps
                                    missing.contact.age="sample",    ## What to do with contacts with missing contact age (possibilities: sample (sample from population), remove (remove them))
                                    estimated.contact.age="sample",
                                    weigh.age=TRUE,                  ## TRUE if weight on the age is performed
                                    weigh.dayofweek=TRUE,            ## TRUE if weight on the day of the week is performed
                                    return.demography   = TRUE,
                                    weight.threshold = 3             ## maximum weight for one participant
  )
  
  age_names <- matrixCoMix_start$demography$age.group
  di_names <- list(age_names,age_names)
  nbreaks <- length(age_breaks)
  S_vec <-matrixCoMix_start$demography[['population']]
  ###Initialize the list containing wave-specific measures of sensitivities
  matrices_boot <- vector(mode = "list", length = nboot)
  eigen_boot <- vector(mode = "list", length = nboot)
  
  
  for (boot in seq_len(nboot)) {
    
    M_mean_asym <- matrix(0,nbreaks,nbreaks)
    M_asym<-matrix(0,nbreaks,nbreaks)
    
    if('demography' %in% names(matrixCoMix_start) && !any(is.na(matrixCoMix_start$matrices[[boot]]))){
      num_age_groups <- nrow(matrixCoMix_start$demography)
      pop_matrix     <- matrix(rep(S_vec,num_age_groups),ncol=num_age_groups,byrow = T)
      
    }
    if(bool_reciprocal==TRUE && all(cnt_location==opt_location)){
      M_mean_asym <- symmetrize(matrixCoMix_start$matrices[[boot]]$matrix,N=S_vec)
      M_asym <-  M_mean_asym/pop_matrix
    }else if(bool_reciprocal==TRUE && !all(cnt_location==opt_location)){
      M_mean_asym <- symmetrize(matrixCoMix_start$matrices[[boot]]$matrix,N=S_vec)
      M_asym <-  M_mean_asym/pop_matrix
      warning('The contact matrices have been symmetrize, but some of the reported contacts may not be reciprocal: e.g. costumer(work) - client(leisure)')
    }else{
      M_mean_asym <- matrixCoMix_start$matrices[[boot]]$matrix
      M_asym <-  M_mean_asym/pop_matrix
      
    }
    
    
    
    M_mean_sym <- matrix(0,nbreaks,nbreaks)
    M_sym<-matrix(0,nbreaks,nbreaks)
    for(loc in opt_location){
      survey_object_CoMix_loc<-get_survey_object(country="Belgium 2010* (Van Hoang 2020)",
                                                 daytype=daytype,
                                                 touch=touch,
                                                 duration=duration,
                                                 gender=gender,
                                                 cnt_location=loc,
                                                 bool_reciprocal=F,              ### should the data generate symmetrical matrix?
                                                 bool_suppl_professional_cnt=bool_suppl_professional_cnt,   ### should supplementary working contacts be included?
                                                 bool_hhmatrix_selection=bool_hhmatrix_selection,      ###select just household members contact
                                                 wave="All",                      ### wave to be included (optional)
                                                 quiet = TRUE)
      
      
      matrixCoMix_loc<-contact_matrix(survey=survey_object_CoMix_loc, ## First and most important parameter: the survey object
                                      age.limits = age_breaks,
                                      n=1,                             ## number of bootstraps
                                      missing.contact.age="sample",    ## What to do with contacts with missing contact age (possibilities: sample (sample from population), remove (remove them))
                                      estimated.contact.age="sample",
                                      weigh.age=TRUE,                  ## TRUE if weight on the age is performed
                                      weigh.dayofweek=TRUE,            ## TRUE if weight on the day of the week is performed
                                      return.demography   = TRUE,
                                      weight.threshold = 3             ## maximum weight for one participant
      )
      
      
      
      
      if('demography' %in% names(matrixCoMix_loc) && !any(is.na(matrixCoMix_loc$matrix))){
        num_age_groups <- nrow(matrixCoMix_loc$demography)
        pop_matrix     <- matrix(rep( S_vec ,num_age_groups),ncol=num_age_groups,byrow = T)
        #matrixCoMix_start$matrix_per_capita <- matrixCoMix_start$matrix/ pop_matrix
      }
      
      C_mean <- matrixCoMix_loc$matrix
      M_mean_sym <-M_mean_sym + (1*(loc=="Home")+ 0.09*(loc=="Work")+0.09*(loc=="School")+0.13*(loc=="Transport")+0.06*(loc=="Leisure")+0.25*(loc=="Otherplace"))*C_mean
      
    }
    if(bool_reciprocal==TRUE && all(cnt_location==opt_location)){
      M_mean_sym <- symmetrize(M_mean_sym,N=S_vec)
      M_sym <-  M_mean_sym/pop_matrix
    }else if(bool_reciprocal==TRUE && !all(cnt_location==opt_location)){
      M_mean_sym <- symmetrize(M_mean_sym,N=S_vec)
      M_sym <-  M_mean_sym/pop_matrix
      warning('The contact matrices have been symmetrize, but some of the reported contacts may not be reciprocal: e.g. costumer(work) - client(leisure)')
    }else{
      M_sym <-  M_mean_sym/pop_matrix
    }
    
    
    
    Tboot_sym <-A%*%M_sym%*%(q_sym*H)
    Tboot_asym <- A%*%M_asym%*%(q_asym*H)
    
    
    for(i in 1:ncol(Tboot_sym)){
      for(j in 1:ncol(Tboot_sym)){
        Tboot_asym[i,j] <-  Tboot_asym[i,j]*S_vec[i]
        Tboot_sym[i,j] <- Tboot_sym[i,j]*S_vec[i]
      }
    }
    
    
    K_1_b <- Tboot_asym%*%D_1_asym +Tboot_sym%*%D_2_sym
    R_0b <- max(Re(eigen(K_1_b)$values))
    damp_ratio <- R_0b/Re(eigen(K_1_b)$values)[2]
    w_b <- abs(eigen(K_1_b)$vectors[,1])
    v_b <- abs(eigen(t(K_1_b))$vectors[,1])
    w_b <- w_b/sum(w_b)
    v_b <- v_b/sum(v_b*w_b)
    
    #initialize matrices
    S_built_b<- matrix(calc.sens(K_1_b),nbreaks,nbreaks,dimnames=di_names)
    E_built_b<- matrix(hadamard.prod(S_built_b,K_1_b)/R_0b,nbreaks,nbreaks,dimnames=di_names)
    
   
    analysis_mat <- lower_sens_matrices(K= K_1_b,               
                                        M_mean_asym,
                                        M_mean_sym,
                                        Susc = S_vec,
                                        q_factor = q_sym)
    
    
    
    matrices_boot[[boot]] <- list(C_mean_asym=M_mean_asym,
                                  C_mean_sym=M_mean_sym,
                                  C_asym=M_asym,
                                  C_sym=M_sym,
                                  NGM=K_1_b,
                                  S=S_built_b,
                                  E=E_built_b,
                                  S_asym=analysis_mat$S_C_asym,
                                  S_sym=analysis_mat$S_C_sym,
                                  E_asym= analysis_mat$E_asym,
                                  E_sym= analysis_mat$E_sym,
                                  Sens_A=analysis_mat$S_A,
                                  Sens_H=analysis_mat$S_H,
                                  El_Theta=analysis_mat$E_T,
                                  El_S1=analysis_mat$E_S1 ,
                                  El_D2=analysis_mat$E_D2 ,
                                  El_O=analysis_mat$E_O,
                                  El_Psi=analysis_mat$E_Psi,
                                  El_Prob=analysis_mat$E_P,
                                  Sens_Ratio=analysis_mat$Sens_ratio
    )
    
    
    
    eigen_boot[[boot]] <- list(R_0=R_0b,
                               wb=w_b,
                               vb=v_b,
                               damp_ratio_b=damp_ratio)
    
    print(boot)
  }
  return(list(Matrices=matrices_boot,Eigen_str=eigen_boot))
  
  
  
  
}


mult_mat <- function(x=list(),M) {
  # Perform the multiplication using the row times column method
  for (i in seq_len(x)) {
    x[[i]] <- x[[i]] %*% M
  }
  
}



NGM_boot <- function(nboot,popmatrix,N_vec,
                          S_vec,
                          reciprocal,check_loc, 
                          output_folder, wave,
                          prop_ratio,
                          q_susc, q_inf, D1, D2, age_breaks,R_pcr){
  
  
  ##Asymptomatic contact matrices
  
  if(wave<9){
    
    fname=paste0(output_folder,"boot_all_wave_",wave,".rds")
    mat_list <- readRDS(fname)
    #nboot <- length(mat_list)
    
  }else{
    set.seed(2345)
    
    
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
    
    
    matrixCoMix_all<-contact_matrix(survey=survey_object_CoMix, ## First and most important parameter: the survey object
                                    age.limits = age_breaks,
                                    n=nboot,                             ## number of bootstraps
                                    missing.contact.age="sample",    ## What to do with contacts with missing contact age (possibilities: sample (sample from population), remove (remove them))
                                    estimated.contact.age="sample",
                                    weigh.age=T,                  ## TRUE if weight on the age is performed
                                    weigh.dayofweek=T,            ## TRUE if weight on the day of the week is performed
                                    return.demography   = TRUE,
                                    weight.threshold = 3             ## maximum weight for one participant
    )
    
    mat_list <- vector(mode = "list", length =nboot)
    
    
    for(i_boot in 1:nboot){
      
        mat_list[[i_boot]] <- matrixCoMix_all$matrices[[i_boot]]$matrix
      
    }
  }
  
  
  
  ##Symptomatic contact matrices
  list_cont_matrices_Home <- list()
  list_cont_matrices_Work<- list()
  list_cont_matrices_School<- list()
  list_cont_matrices_Transport<- list()
  list_cont_matrices_Leisure<- list()
  list_cont_matrices_Otherplace<- list()
  
  for(loc in opt_location){
    
    set.seed(2345)
    survey_object_CoMix_loc<-get_survey_object(country=country,
                                               daytype=daytype,
                                               touch=touch,
                                               duration=duration,
                                               gender=gender,
                                               cnt_location=loc,
                                               bool_reciprocal=F,              ### should the data generate symmetrical matrix?
                                               bool_suppl_professional_cnt=bool_suppl_professional_cnt,   ### should supplementary working contacts be included?
                                               bool_hhmatrix_selection=bool_hhmatrix_selection,      ###select just household members contact
                                               wave=id_waves[wave],                      ### wave to be included (optional)
                                               quiet = TRUE)
    
    matrixCoMix_loc<-contact_matrix(survey=survey_object_CoMix_loc, ## First and most important parameter: the survey object
                                    age.limits = age_breaks,
                                    n=nboot,                             ## number of bootstraps
                                    missing.contact.age="sample",    ## What to do with contacts with missing contact age (possibilities: sample (sample from population), remove (remove them))
                                    estimated.contact.age="sample",
                                    weigh.age=TRUE,                  ## TRUE if weight on the age is performed
                                    weigh.dayofweek=TRUE,            ## TRUE if weight on the day of the week is performed
                                    return.demography   = TRUE,
                                    weight.threshold = 3             ## maximum weight for one participant
    )
    
    
    matrices_boot_loc <- vector(mode = "list", length = nboot_loc)
    
    for(i_boot in 1:nboot){
      
      matrices_boot_loc[[i_boot]] <- matrixCoMix_loc$matrices[[i_boot]]$matrix
      
    }
    
    
    switch (loc,
            "Home" = {fname=paste0(output_folder,"boot_Home_wave_",wave,".rds")
            #list_cont_matrices_Home <- ifelse(wave<9,readRDS(fname), matrices_boot_loc) 
            if (wave<9) {
              list_cont_matrices_Home <- readRDS(fname) 
              
            }else{
              list_cont_matrices_Home <- matrices_boot_loc
            }
            
            },
            "Work" = {fname=paste0(output_folder,"boot_Work_wave_",wave,".rds")
            if (wave<9) {
              list_cont_matrices_Work <- readRDS(fname) 
              
            }else{
              list_cont_matrices_Work <- matrices_boot_loc
            }
            },   
            "School" = {fname=paste0(output_folder,"boot_School_wave_",wave,".rds")
            if (wave<9) {
              list_cont_matrices_School <- readRDS(fname) 
              
            }else{
              list_cont_matrices_School <- matrices_boot_loc
            }
            },
            "Transport"  = {fname=paste0(output_folder,"boot_Transport_wave_",wave,".rds")
            if (wave<9) {
              list_cont_matrices_Transport <- readRDS(fname) 
              
            }else{
              list_cont_matrices_Transport <- matrices_boot_loc
            }
            },
            "Leisure" = {fname=paste0(output_folder,"boot_Leisure_wave_",wave,".rds")
            if (wave<9) {
              list_cont_matrices_Leisure <- readRDS(fname) 
              
            }else{
              list_cont_matrices_Leisure <- matrices_boot_loc
            }
            },
            "Otherplace" = {fname=paste0(output_folder,"boot_Otherplace_wave_",wave,".rds")
            if (wave<9) {
              list_cont_matrices_Otherplace <- readRDS(fname) 
              
            }else{
              list_cont_matrices_Otherplace <- matrices_boot_loc
            }
            }
    )
    
    
    
  }
  
  
   mat_asym_list <-  list()
  mat_sym_list <- list()
  mat_asym_list_pc <-  list()
  mat_sym_list_pc <- list()
  
   if(reciprocal==TRUE && check_loc){
      mat_asym_list <- lapply(mat_list, symmetrize, N=N_vec)
      mat_asym_list <- lapply(mat_list, non_neg)
      mat_asym_list_pc <- lapply(mat_asym_list, function(x) x / pop_matrix)
        
      
    }else if(reciprocal==TRUE && !check_loc){
      mat_asym_list <- lapply(mat_list, symmetrize, N=N_vec)
      mat_asym_list <- lapply(mat_list, non_neg)
      mat_asym_list_pc <- lapply(mat_asym_list, function(x) x/pop_matrix)
      warning('The contact matrices have been symmetrize, but some of the reported contacts may not be reciprocal: e.g. costumer(work) - client(leisure)')
    }else{
      mat_asym_list <- lapply(mat_list, non_neg)
      mat_asym_list_pc <- lapply(mat_asym_list, function(x) x/pop_matrix)
    }
   
    
    
  mat_sym_list <-lapply(seq_along(list_cont_matrices_Home), 
                        function(i) list_cont_matrices_Home[[i]]+ 
                          0.09*list_cont_matrices_Work[[i]]+
                          0.09*list_cont_matrices_School[[i]]+
                          0.13*list_cont_matrices_Transport[[i]]+
                          0.06*list_cont_matrices_Leisure[[i]]+
                          0.25*list_cont_matrices_Otherplace[[i]])
    
    
  mat_sym_list_pc <- list()
    
    
    if(reciprocal==TRUE && check_loc){
      mat_sym_list <- lapply(mat_sym_list, symmetrize, N=N_vec)
      mat_sym_list <- lapply(mat_sym_list, non_neg)
      mat_sym_list_pc <- lapply(mat_sym_list, function(x) x/pop_matrix)
      
    }else if(reciprocal==TRUE && !check_loc){
      mat_sym_list <- lapply(mat_sym_list, symmetrize, N=N_vec)
      mat_sym_list <- lapply(mat_sym_list, non_neg)
      mat_sym_list_pc <- lapply(mat_sym_list, function(x) x/pop_matrix)
      warning('The contact matrices have been symmetrize, but some of the reported contacts may not be reciprocal: e.g. costumer(work) - client(leisure)')
    }else{
      mat_sym_list <- lapply(mat_sym_list, non_neg)
      mat_sym_list_pc <- lapply(mat_sym_list, function(x) x/pop_matrix)
    }
    
    q_asym <- prop_ratio*q_sym
    
    
    Tboot_left <- diag(S_vec)%*%q_susc
    Tboot_asym_right <- q_asym*q_inf
    Tboot_sym_right <- q_sym*q_inf
    
    list_asym <- lapply(mat_asym_list_pc, function(x) Tboot_left%*%x%*%Tboot_asym_right%*%D1)
    list_sym <- lapply(mat_sym_list_pc, function(x) Tboot_left%*%x%*%Tboot_sym_right%*%D2)
    
    K_matrices_boot<- lapply(seq_along(list_asym), 
                             function(i) list_asym[[i]]+list_sym[[i]])
      
      
    if(R_pcr){
      for (k in seq_len(nboot)) {
        R_boot <-max(Re(eigen(K_matrices_boot[[k]])$values)) 
        q_boot <- R_series_mean[wave]/R_boot 
        K_matrices_boot[[k]] <- K_matrices_boot[[k]]*q_boot
      }
      
    }  
    
    
      mapply(function(x,y) x+y,list_asym,list_sym)
    C_matrices_boot_asym<- mat_asym_list
    C_matrices_boot_sym<- mat_sym_list
    Sens_boot <- lapply(K_matrices_boot, calc.sens)
    El_boot <- lapply(K_matrices_boot, calc.el)
    
     
  return(list(K=K_matrices_boot,          #matrix evaluated on each bootstapped contact matrices 
              C_asym=C_matrices_boot_asym, #after having rendered reciprocal every bootstrapped contact matrices
              C_sym=C_matrices_boot_sym,
              C_row=mat_list,
              Sens=Sens_boot,
              El=El_boot))
              
              
}
  
  


