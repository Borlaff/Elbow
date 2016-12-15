# Elbow
Elbow: a statistically robust and automated method to fit and classify the surface brightness profiles

Note to the users: 

Elbow runs in R. 

  R Core Team (2016). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna,
  Austria. URL https://www.R-project.org/.

In order to install R, follow the instructions at https://cloud.r-project.org/

Once to installed it, just type in your terminal at the Elbow directory:

> R
>> source("elbow_demo.r")

for a working example. Follow the notes in both .r files (Elbow_demo.r & Elbow_X.X.r") 

Don't hesitate to contact me at: 
Alejandro S. Borlaff: E-mail: asborlaff@ucm.es - asborlaff@iac.es 

###########################################################################################################
 NAME:
   Elbow

 PURPOSE:
   Fit double exponential functions to surface brightness profiles and
   calculate confidence values. 

 The program requires the 'doParallel', 'foreach', 'tools', 'plotrix' and 'modeest' libraries to run.  
 

 :Categories:
    Statistics, Surface Brightness profiles, Galaxy structure. 
    
 INPUTS (MANDATORY):
 
 r: Radius - X: The vector containing the radius of your galaxy. You can use either kpc, arcsec or pix, but be 
aware that the inner and outer scalelengths, (hi, ho), and the break radius (rbreak) will be in the same units. 
                 An n-element vector containing the independent variable values.
                X may be of type integer, floating point, or double-precision floating-point.
                
 r_down,r_up: 1sigma confidence limits for the r (radius) value. Their values should be r_down < r < r_up.
                 A common error is to introduce relative uncertainities to the central error, not absolute values.  
 
 mu: Magnitude - Y: The vector containing the surface brightness profile of your object, for each radius element. The r, r_up, r_down and mu, mu_up, mu_down must be of the same length. The units are (mag/arcsec^2).  An n-element integer, single-, or double-precision floating-point vector.
 
 mu_down,mu_up: 1sigma confidence limits for the mu (magnitude) value. Their values should be mu_down > mu > mu_up.  
                   A common error is to mistake down with lower (numeric) magnitudes, which are brighter intensities. 
                   The _down stands for the brightness, not the numeric (inverse logarithmic) scale. 
 
 zeropoint: An arbitrary zp to transform between intensity and magnitudes. The output magnitude values will be consistent 
              with this zp. 
 pix_size:  A pixel size to transform intensity and magnitudes.The output magnitude values will be consistent 
              with this pix_size.
 
 INPUTS (OPTIONAL): 
   min_lim (Default: NA): If supplied, Elbow will internally exclude those points with r < min_lim from the fit.
   
   max_lim (Default: NA): If supplied, Elbow will internally exclude those points with r > max_lim from the fit.  
   
   nsimul_break (Default=1000): The number of simulations that will be used for the Bootstrapping + Monte Carlo fit. 
                                A reasonable number would be at least 10^4 simulations.   
          
          
 OUTPUTS: 
 The output will be a numeric named vector with the following parameters:
 
names(results)<-c("median_muoi","muoi_up","muoi_down","median_muoo","muoo_up","muoo_down","median_hi","hi_up","hi_down","median_ho","ho_up","ho_down","p_muo","p_ho","median_rbrk","rbrk_up","rbrk_down","median_mubrk","mubrk_up","mubrk_down")

median_muoi: Median central surface brightness (\muoi) for the inner profile. 
muoi_up / muoi_down: 1sigma confidence interval for the central surface brightness (\muoi) for the inner profile. 
 
median_muoo: Median central surface brightness (\muoo) for the outer profile. 
muoo_up / muoo_down: 1sigma confidence interval for the central surface brightness (\muoo) for the outer profile. 

median_hi: Median scale-length (\hi) for the inner profile. 
hi_up / hi_down: 1sigma confidence interval for the scale-length (\hi) for the inner profile. 
 
median_ho: Median scale-length (\ho) for the outer profile. 
ho_up / ho_down: 1sigma confidence interval for the scale-length (\ho) for the outer profile. 

p_muo: Likelihood that the central surface brightness of the inner and outer profiles are coming from the same parent distribution (i.e, if p_muo < 0.05, you have >95% likelihood that the break is real. A Type-I, pure exponential profile will return p_muo >> 0.05). 
p_ho: Likelihood that the scale-lenghts of the inner and outer profiles are coming from the same parent distribution.

median_rbrk: Median break radius for the profile. 
rbrk_up / rbrk_down: 1sigma confidence interval for the break radius for the profile.

median_mubrk: Median surface brightness at the break radius for the profile. 
mubrk_up / mubrk_down: 1sigma confidence interval for the surface brightness at the break radius for the profile.

       
 :Author:
       
           Alejandro S. Borlaff 
           IAC Researcher
           Instituto de Astrofísica de Canarias
           C/ Vía Láctea, s/n
           E38205 - San Cristobal de La Laguna (Santa Cruz de Tenerife). Spain
           E-mail: asborlaff@ucm.es - asborlaff@iac.es 
 :History:
     Change History::
        Written, January - December 2016. First release. 
#############################################################################################################
