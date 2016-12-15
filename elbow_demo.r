rm(list=ls())
source('elbow_1.1.R') 

#######################################################################################
# Elbow 1.1 - Demo 
# This is a demonstrator/tutorial for the use of Elbow: automated break analysis of disc surface brightness profiles 
# The following program generates a double exponential profile with variable parameters
# and tries to fit it with Elbow. The output of the program will be a named vector with the parameters 
# and will create a table in the same directory with the simulation results. 
# This demo program will load the table to create several plots with the probability distributions of the breaks. 
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
# mu1(R)= mui + 2.5/log(10) * R / hi 
# mu2(R)= muo + 2.5/log(10) * R / ho
# mu = -2.5*log10(int) + 25

# Play with these parameters to generate an antitruncated profile that Elbow will try to fit. 

mui = 20
muo = 22.5
hi = 1
ho = 3 
maglim = 26
nsimul = 1000

########################################################################################
r <- seq(0,30, by=0.1)
mu1 = mui + 2.5*r/log(10)/hi
mu2 = muo + 2.5*r/log(10)/ho
mu_model <- c()
for(i in 1:length(r)){
  mu_model[i]<-min(mu1[i],mu2[i])
}

int = 10^((25 - mu_model)/2.5)
int_lim = 10^((25 - maglim)/2.5)

int_noise = int + rnorm(length(int),0,int_lim)
int_up <- int_noise + int_lim
int_down <- int_noise - int_lim


mu =  -2.5*log10(int_noise) + 25
mu_up =  -2.5*log10(int_up) + 25
mu_down =  -2.5*log10(int_down) + 25
r_up <- r + 0.1
r_down <- r - 0.1 

max_lim = r[which(mu>maglim)][1]

plot(r,mu, pch=20, ylim=c(maglim+2,18), xlim=c(0,15), ylab=expression(paste(mu,"(r) (mag arcsec"^"-2"*")")), xlab=expression('R (arcsec)'))
abline(v=max_lim, lwd=2)

segments(y0=mu_down, y1=mu_up, x0=r, lwd=1, col="black")
segments(x0=r_up, x1=r_down, y0=mu, lwd=1, col="black")
abline(h=maglim, lwd=3, col="red")

elbow_fit <- elbow(r = r, r_down = r_down, r_up = r_up, 
                   mu = mu, mu_down = mu_down, mu_up,
                   zeropoint = 25, pix_size = 1, max_lim=max_lim, nsimul_break = nsimul)

PDD <- read.delim("PDD.dat")


abline(a=elbow_fit$median_muoi, b=2.5/log(10)/elbow_fit$median_hi, lwd=2, col="blue")
abline(a=elbow_fit$median_muoo, b=2.5/log(10)/elbow_fit$median_ho, lwd=2, col="red")

hist(PDD$hi,breaks=100,col="grey", main=expression(paste("h",""["i"]*" (arcsec)")))
abline(v=elbow_fit$median_hi,lwd=2)
box(lwd=2)
rug(PDD$hi)

hist(PDD$muoi,breaks=100,col="grey", main=expression(paste(mu,""["o"]*" (mag arcsec"^"-2"*")")))
abline(v=elbow_fit$median_muoi,lwd=2)
box(lwd=2)
rug(PDD$muoi)


hist(PDD$ho,breaks=100,col="grey", main=expression(paste("h",""["o"]*" (arcsec)")))
abline(v=elbow_fit$median_ho,lwd=2)
box(lwd=2)
rug(PDD$ho)

hist(PDD$muoo,breaks=100,col="grey", main=expression(paste(mu,""["o"]*" (mag arcsec"^"-2"*")")))
abline(v=elbow_fit$median_muoo,lwd=2)
box(lwd=2)
rug(PDD$muoo)

xlab=expression('h (arcsec)')
ylab=expression(paste(mu,"(r) (mag arcsec"^"-2"*")"))

colramp = colorRampPalette(c('white', 'blue', 'green', 'yellow', 'orange', 'red', 'magenta','purple'))
smoothScatter(c(PDD$hi,PDD$ho), c(PDD$muoi,PDD$muoo), pch=".",cex=1, xlim=c(0,15),
              ylim=c(24,20), ylab=ylab, xlab=xlab, colramp=colramp, nbin=512, bandwidth=c(0.1,0.1))

points(elbow_fit$median_hi, elbow_fit$median_muoi, pch=19)
points(hi, mui, pch=1, cex=2, col="black")
segments(x0=elbow_fit$hi_down, x1=elbow_fit$hi_up, y0=elbow_fit$median_muoi, lwd=2, col="black")
segments(x0=elbow_fit$median_hi, y0=elbow_fit$muoi_down, y1=elbow_fit$muoi_up, lwd=2, col="black")

points(elbow_fit$median_ho, elbow_fit$median_muoo, pch=19)
points(ho, muo, pch=1, cex=2, col="black")
segments(x0=elbow_fit$ho_down, x1=elbow_fit$ho_up, y0=elbow_fit$median_muoo, lwd=2, col="black")
segments(x0=elbow_fit$median_ho, y0=elbow_fit$muoo_down, y1=elbow_fit$muoo_up, lwd=2, col="black")

#########################################################################
signif_mubrk_up <- floor(-log10(elbow_fit$median_mubrk-elbow_fit$mubrk_up) + 2)
signif_mubrk_down <- floor(-log10(elbow_fit$mubrk_down-elbow_fit$median_mubrk) + 2)
signif_mubrk <- max(signif_mubrk_up, signif_mubrk_down)

mubrk_l <- formatC(elbow_fit$median_mubrk,format="f",digits=signif_mubrk)
mubrk_up_l <- formatC((elbow_fit$mubrk_up-elbow_fit$median_mubrk),format="f",digits=signif_mubrk)
mubrk_down_l <-  formatC((elbow_fit$mubrk_down-elbow_fit$median_mubrk), format="f",digits=signif_mubrk)
#########################################################################  
signif_rbrk_up <- floor(-log10(elbow_fit$rbrk_up-elbow_fit$median_rbrk) + 2)
signif_rbrk_down <- floor(-log10(elbow_fit$median_rbrk-elbow_fit$rbrk_down) + 2)
signif_rbrk <- max(signif_rbrk_up, signif_rbrk_down)

rbrk_l <- formatC(elbow_fit$median_rbrk,format="f",digits=signif_rbrk)
rbrk_up_l <- formatC((elbow_fit$rbrk_up-elbow_fit$median_rbrk),format="f",digits=signif_rbrk)
rbrk_down_l <-  formatC((elbow_fit$median_rbrk-elbow_fit$rbrk_down), format="f",digits=signif_rbrk)
#########################################################################  
signif_muoi_up <- floor(-log10(elbow_fit$median_muoi-elbow_fit$muoi_up) + 2)
signif_muoi_down <- floor(-log10(elbow_fit$muoi_down-elbow_fit$median_muoi) + 2)
signif_muoi <- max(signif_muoi_up, signif_muoi_down)

muoi_l <- formatC(elbow_fit$median_muoi,format="f",digits=signif_muoi)
muoi_up_l <- formatC((elbow_fit$median_muoi-elbow_fit$muoi_up),format="f",digits=signif_muoi)
muoi_down_l <-  formatC((elbow_fit$muoi_down-elbow_fit$median_muoi), format="f",digits=signif_muoi)
#########################################################################  
signif_hi_up <- floor(-log10(elbow_fit$hi_up-elbow_fit$median_hi) + 2)
signif_hi_down <- floor(-log10(elbow_fit$median_hi-elbow_fit$hi_down) + 2)
signif_hi <- max(signif_hi_up, signif_hi_down)

hi_l <- formatC(elbow_fit$median_hi,format="f",digits=signif_hi)
hi_up_l <- formatC((elbow_fit$hi_up-elbow_fit$median_hi),format="f",digits=signif_hi)
hi_down_l <-  formatC((elbow_fit$median_hi-elbow_fit$hi_down), format="f",digits=signif_hi)
#########################################################################  
signif_muoo_up <- floor(-log10(elbow_fit$median_muoo-elbow_fit$muoo_up) + 2)
signif_muoo_down <- floor(-log10(elbow_fit$muoo_down-elbow_fit$median_muoo) + 2)
signif_muoo <- max(signif_muoo_up, signif_muoo_down)

muoo_l <- formatC(elbow_fit$median_muoo,format="f",digits=signif_muoo)
muoo_up_l <- formatC((elbow_fit$median_muoo-elbow_fit$muoo_up),format="f",digits=signif_muoo)
muoo_down_l <-  formatC((elbow_fit$muoo_down-elbow_fit$median_muoo), format="f",digits=signif_muoo)
#########################################################################  
signif_ho_up <- floor(-log10(elbow_fit$ho_up-elbow_fit$median_ho) + 2)
signif_ho_down <- floor(-log10(elbow_fit$median_ho-elbow_fit$ho_down) + 2)
signif_ho <- max(signif_ho_up, signif_ho_down)

ho_l <- formatC(elbow_fit$median_ho,format="f",digits=signif_ho)
ho_up_l <- formatC((elbow_fit$ho_up-elbow_fit$median_ho),format="f",digits=signif_ho)
ho_down_l <-  formatC((elbow_fit$median_ho-elbow_fit$ho_down), format="f",digits=signif_ho)
#########################################################################  

legend("bottomright", legend=c("Model's inner profile","Model's outer profile","Elbow's inner profile fit","Elbow's outer profile fit"), pch=c(1,19))
legend(legend=c(paste("Results:"), 
                as.expression(bquote(paste("R",""["break"]*" = ",.(rbrk_l)[paste("-",.(rbrk_down_l),sep="")]^paste("+",.(rbrk_up_l),sep="")," arcsec"))),
                as.expression(bquote(paste(mu,""["break"]*" = ",.(mubrk_l)[paste("-",.(mubrk_down_l),sep="")]^paste("+",.(mubrk_up_l),sep="")," mag arcsec"^"-2"*""))),
                as.expression(bquote(paste("h",""["0,i"]*" = ",.(hi_l)[paste("-",.(hi_down_l),sep="")]^paste("+",.(hi_up_l),sep="")," arcsec"))),
                as.expression(bquote(paste(mu,""["0,i"]*" = ",.(muoi_l)[paste("-",.(muoi_down_l),sep="")]^paste("+",.(muoi_up_l),sep="")," mag arcsec"^"-2"*""))),
                as.expression(bquote(paste("h",""["0,o"]*" = ",.(ho_l)[paste("-",.(ho_down_l),sep="")]^paste("+",.(ho_up_l),sep="")," arcsec"))),
                as.expression(bquote(paste(mu,""["0,o"]*" = ",.(muoo_l)[paste("-",.(muoo_down_l),sep="")]^paste("+",.(muoo_up_l),sep="")," mag arcsec"^"-2"*"")))),
       "topright",y.intersp=1.5, bg="white")

