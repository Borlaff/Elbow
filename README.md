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
 r: Radius - X: An n-element vector containing the independent variable values.
                X may be of type integer, floating point, or double-precision floating-point.
 r_down,r_up: 1sigma confidence limits for the r (radius) value. Their values should be r_down < r < r_up.
                 A common error is to introduce relative uncertainities to the central error, not absolute values.  
 
 mu: Magnitude - Y: An n-element integer, single-, or double-precision floating-point vector.
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
