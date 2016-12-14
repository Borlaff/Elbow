library("doParallel")
library("foreach")
library("tools")
library("plotrix")
library("modeest")

###########################################################################################################
# NAME:
#   Elbow
#
# PURPOSE:
#   Fit double exponential functions to surface brightness profiles and
#   calculate confidence values. 
#
# The program requires the 'doParallel', 'foreach', 'tools', 'plotrix' and 'modeest' libraries to run.  
# 
#
# :Categories:
#    Statistics, Surface Brightness profiles, Galaxy structure. 
#    
# INPUTS (MANDATORY):
#    r: Radius - X: An n-element vector containing the independent variable values.
#                X may be of type integer, floating point, or double-precision floating-point.
# 
#         An n-element integer, single-, or double-precision floating-point vector.
#    r_down,r_up: 1sigma confidence limits for the r (radius) value. Their values should be r_down < r < r_up.
#                 A common error is to introduce relative uncertainities to the central error, not absolute values.  
# 
#    mu: Magnitude - Y: An n-element integer, single-, or double-precision floating-point vector.
# 
#    mu_down,mu_up: 1sigma confidence limits for the mu (magnitude) value. Their values should be mu_down > mu > mu_up.  
#                   A common error is to mistake down with lower (numeric) magnitudes, which are brighter intensities. 
#                   The _down stands for the brightness, not the numeric (inverse logarithmic) scale. 
# 
#   zeropoint: An arbitrary zp to transform between intensity and magnitudes. The output magnitude values will be consistent 
#              with this zp. 
#   pix_size:  A pixel size to transform intensity and magnitudes.The output magnitude values will be consistent 
#              with this pix_size.
# 
# INPUTS (OPTIONAL): 
#   min_lim (Default: NA): If supplied, Elbow will internally exclude those points with r < min_lim from the fit. 
#   max_lim (Default: NA): If supplied, Elbow will internally exclude those points with r > max_lim from the fit.  
#   nsimul_break (Default=1000): The number of simulations that will be used for the Bootstrapping + Monte Carlo fit. 
#                                A reasonable number would be at least 10^4 simulations.   
#          
#       
# :Author:
#       
#           Alejandro S. Borlaff 
#           IAC Researcher
#           Instituto de Astrofísica de Canarias
#           C/ Vía Láctea, s/n
#           E38205 - San Cristobal de La Laguna (Santa Cruz de Tenerife). Spain
#           E-mail: asborlaff@ucm.es - asborlaff@iac.es 
# :History:
#     Change History::
#        Written, January - December 2016. First release. 
###############################################################################################################

unregister <- function(){
  # Unregister is a function to clean the multithread pool cores used in parallel processing. 
  # 
  env<-foreach:::.foreachGlobals
  rm(list=ls(name=env),pos=env)
}




elbow<-function(r,r_down,r_up,mu,mu_down,mu_up,zeropoint, pix_size=1, min_lim=NA, max_lim=NA,nsimul_break=1000){

  print("Elbow v.1.1")
  
  filter_params <- function(muoi,muoo,hi,ho){
    # Add filter conditions to remove extreme values in the distributions: CAUTION, too restrictive filters will affect the p-values. USE AT YOUR OWN RISK! 
    # Examples 
    filter1 <- (muoo > 999)
    filter2 <- (muoi > 999)
    filter3 <- (muoo < 1)
    filter4 <- (muoi < 1)
    filter5 <- (ho < -999)
    filter6 <- (ho > 999)
    filter7 <- (hi < -999)
    filter8 <- (hi > 999)
    return(c(filter1 | filter2 | filter3 | filter4 | filter5 | filter6 | filter7 | filter8))
  }
  
  
  
  fit_broken_lm <- function(x,y, plot_mode) {
    f <- function (Cx){
      lhs <- function(x) ifelse(x < Cx,Cx-x,0) # Given a certain break, this function returns the distance to that point of the points with values below rbreak. The rest is equal to 0.  
      rhs <- function(x) ifelse(x < Cx,0,x-Cx) # Given a certain break, this function returns the distance to that point of the points with values over rbreak. The rest is equal to 0.  
      fit <- lm(y ~ lhs(x) + rhs(x)) # Here we apply a double linear profile 
      c(summary(fit)$r.squared, 
        summary(fit)$coef[1], summary(fit)$coef[2],
        summary(fit)$coef[3]) 
    }
    
    r2 <- function(x) -(f(x)[1]) # This function return the r2 value for optimize. 
    
    res <- optimize(r2,interval=c(min(x),max(x))) # Perform optimize to find best rbreak 
    res <- c(res$minimum,f(res$minimum)) # We calculate the coeficients for the best rbreak. 
    
    best_Cx <- res[1]
    coef1 <- res[3]
    coef2 <- res[4]
    coef3 <- res[5]
    if(plot_mode!="none"){
      if(plot_mode=="plot") plot(x,y,pch=20,col="black")
      if(plot_mode=="points") points(x,y,pch=20,col="black")
      abline(coef1+best_Cx*coef2,-coef2) # lhs = Inner profile
      abline(coef1-best_Cx*coef3,coef3)  # rhs = Outer profile
    }
    return(c(coef1+best_Cx*coef2,-coef2,coef1-best_Cx*coef3,coef3,best_Cx))
  }
  
  sigma_error <- 0.682689492137086 # 1sigma  
  
  # We identify the section of the profile to study.
  if(is.na(min_lim)) min_lim = min(r)
  if(is.na(max_lim)) max_lim = max(r)
  valid_interval = which(r>=min_lim)[1]:which(r>=max_lim)[1]
  
  radius<-r[valid_interval]
  radius_err_up<-abs(radius-r_up[valid_interval])
  radius_err_down<-abs(radius-r_down[valid_interval])
  magnitude<-mu[valid_interval]
  magnitude_up<-mu_up[valid_interval]
  magnitude_down<-mu_down[valid_interval]
  intens<-pix_size^2*10^((zeropoint-magnitude)/2.5)
  intens_down<-pix_size^2*10^((zeropoint-magnitude_down)/2.5)
  intens_up<-pix_size^2*10^((zeropoint-magnitude_up)/2.5)
  
  # Preparing the storage arrays. 
  magnitude_bootMC<-matrix(nrow=nsimul_break, ncol=length(radius))
  radius_bootMC<-matrix(nrow=nsimul_break, ncol=length(radius))
  fit_broken_MC<-matrix(nrow=nsimul_break, ncol=4)


  # Multithreading pool - Only with parallel processing. 
    cat("\nHow many cores do you have?: ",detectCores())
    cl<-makeCluster(detectCores()-2)
    registerDoParallel(cl)
    cat("\nHow many cores joined the cluster?",getDoParWorkers(),"\n")
  
      
  # Now performing the real fitting procedure. 
  print("Stand by one, this may take a while...")
    fit_broken_MC<-foreach(i=1:nsimul_break, .combine="cbind") %dopar% {
      # Bootstrapping: We select the points chosen a sample of size N with replacement. 
      index<-sample(seq(1:length(radius)), size=length(radius), replace=TRUE) 
      
      radius_MC<-radius[index]
      radius_err_MC_down<-radius_err_down[index]
      radius_err_MC_up<-radius_err_up[index]
      intens_MC<-intens[index]
      intens_down_MC<-intens_down[index]
      intens_up_MC<-intens_up[index]
      
      # Monte Carlo: We move the data points within their respective errors. 
      sigma_intens<-abs(intens_up_MC-intens_down_MC)/2
      radius_MC<-radius_MC + rnorm(length(index),0,(radius_err_MC_up+radius_err_MC_down)/2) 
      intens_MC<-intens_MC + rnorm(length(index),0,sigma_intens)

      # We pass the profile to fit_broken_lm and fit this simulation. 
      fit_broken_MC_vector<-fit_broken_lm(radius_MC,-2.5*log10(intens_MC/pix_size^2)+zeropoint,plot_mode="none")
      fit_broken_MC_vector[2]<-2.5/(log(10)*fit_broken_MC_vector[2])
      fit_broken_MC_vector[4]<-2.5/(log(10)*fit_broken_MC_vector[4])
      fit_broken_MC <- fit_broken_MC_vector
    }
   print("Simulations done!")
    unregister() # We remove the core pool. 
  
    # Extracting values from the fit. 
    muoi=fit_broken_MC[1,]
    hi=fit_broken_MC[2,]
    muoo=fit_broken_MC[3,]
    ho=fit_broken_MC[4,]

  
  # We calculate the distribution of Rbrk by assuming mu,i(Rbrk) == mu,o(Rbrk)
  rbrk <- (ho*log(10)/2.5)*((muoi - ho*muoo/hi)/(1-ho/hi)-muoo)
  sd_rbrk<-sd(rbrk,na.rm=TRUE)
  rbrk_down<-quantile(rbrk, (1-sigma_error)/2, na.rm=TRUE)
  rbrk_up<-quantile(rbrk, (1-(1-sigma_error)/2), na.rm=TRUE)
  
  # The distribution of mubrk is calculated by interpolation over the profile.
  mubrk<-(approx(r, mu, xout=rbrk))$y
  sd_mubrk<-sd(mubrk,na.rm=TRUE)
  mubrk_up<-quantile(mubrk, (1-sigma_error)/2 ,na.rm=TRUE)
  mubrk_down<-quantile(mubrk, (1-(1-sigma_error)/2),na.rm=TRUE)
  
  filter <- filter_params(muoi=muoi, muoo=muoo, hi=hi, ho=ho)
  
  if(length(which(filter))!=0){
    muoo_clean <- muoo[-which(filter)]
    muoi_clean <- muoi[-which(filter)]
    ho_clean <- ho[-which(filter)]
    hi_clean <- hi[-which(filter)]
  } else {
    muoo_clean <- muoo
    muoi_clean <- muoi
    ho_clean <- ho
    hi_clean <- hi
  }
  # We calculate the best muoi and muoo values with a mode (most probable solution)
  median_muoi<-mlv(muoi, method = "venter", type = "shorth", na.rm=TRUE)$M
  median_muoo<-mlv(muoo,  method = "venter", type = "shorth",na.rm=TRUE)$M
  #median_muoi<-median(muoi, na.rm=TRUE)
  #median_muoo<-median(muoo, na.rm=TRUE)
  
  # The p-values of muo are calculated as the fraction of simulations that give a different results than the "medium" one. 
  if(median_muoo > median_muoi){
    p_muo<-max(1/nsimul_break,length(which((muoo_clean-muoi_clean)<0))/(length(muoo_clean)))  
  } else {
    p_muo<-max(1/nsimul_break,length(which((muoi_clean-muoo_clean)<0))/(length(muoo_clean)))  
  }
  muoi_up<-quantile(muoi, (1-sigma_error)/2,na.rm=TRUE)
  muoi_down<-quantile(muoi, (1-(1-sigma_error)/2),na.rm=TRUE)
  muoo_up<-quantile(muoo, (1-sigma_error)/2,na.rm=TRUE)
  muoo_down<-quantile(muoo, (1-(1-sigma_error)/2),na.rm=TRUE)
  
  # We calculate the best hi and ho values with a mode (most probable solution)
  median_hi<-mlv(hi,  method = "venter", type = "shorth", na.rm=TRUE)$M
  median_ho<-mlv(ho,  method = "venter", type = "shorth", na.rm=TRUE)$M
  #median_hi<-median(hi, na.rm=TRUE)
  #median_ho<-median(ho, na.rm=TRUE)
  
  # The p-values of muo are calculated as the fraction of simulations that give a different results than the "medium" one. 
  if(median_ho > median_hi){
    p_ho<-max(1/nsimul_break,length(which((ho_clean-hi_clean)<0))/(length(ho_clean))) 
  } else {
    p_ho<-max(1/nsimul_break,length(which((hi_clean-ho_clean)<0))/(length(ho_clean))) 
  }
  hi_down<-quantile(hi, (1-sigma_error)/2,na.rm=TRUE)
  hi_up<-quantile(hi, (1-(1-sigma_error)/2),na.rm=TRUE)
  ho_down<-quantile(ho, (1-sigma_error)/2,na.rm=TRUE)
  ho_up<-quantile(ho, (1-(1-sigma_error)/2),na.rm=TRUE)
  
  # We calculate the distribution of Rbrk by assuming mu,i(Rbrk) == mu,o(Rbrk)
  median_rbrk <- (median_ho/2.5*log(10))*((median_muoi - median_ho*median_muoo/median_hi)/(1- median_ho/median_hi)-median_muoo)
  median_mubrk <- (approx(r, mu, median_rbrk))$y
  

  results<-c(median_muoi,muoi_up,muoi_down,median_muoo,muoo_up,muoo_down,median_hi,hi_up,hi_down,median_ho,ho_up,ho_down, p_muo, p_ho, median_rbrk, rbrk_up, rbrk_down,median_mubrk, mubrk_up, mubrk_down)
  names(results)<-c("median_muoi","muoi_up","muoi_down","median_muoo","muoo_up","muoo_down","median_hi","hi_up","hi_down","median_ho","ho_up","ho_down","p_muo","p_ho","median_rbrk","rbrk_up","rbrk_down","median_mubrk","mubrk_up","mubrk_down")
  print(results)
  write.table(data.frame(hi,muoi,ho,muoo,rbrk,mubrk), file="PDD.dat", col.names = TRUE, row.names = FALSE, quote=FALSE, sep="\t")
  
  return(as.data.frame(t(results)))
}
