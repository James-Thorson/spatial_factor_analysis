Description
=======================

Spatial factor analysis: a tool for decomposing species density into a small number of latent maps


Instructions
=============
First, install the "devtools" package from CRAN

    # Install and load devtools package
    install.packages("devtools")
    library("devtools")

Second, please install the following:
* TMB (Template Model Builder): https://github.com/kaskr/adcomp
* INLA (integrated nested Laplace approximations): http://www.r-inla.org/download

Note: at the moment, TMB and INLA can be installed using the commands 

    # devtools command to get TMB from GitHub
    install_github("kaskr/adcomp/TMB") 
    # source script to get INLA from the web
    source("http://www.math.ntnu.no/inla/givemeINLA.R")  
    
Next, please install the geostatistical_delta-GLMM package from this GitHub repository using a function in the "devtools" package:

    # Install package
    install_github("James-Thorson/spatial_factor_analysis") 
    # Load package
    library(SpatialFA)

Please see the examples folder for an example of how to run the model:
https://github.com/James-Thorson/spatial_factor_analysis/tree/master/examples

Known installation/usage issues
=============
none

Further reading
=============

For more details regarding development and testing of this software please see:

* Thorson, J.T., Scheuerell, M., Shelton, A.O., See, K., Skaug, H.J., and Kristensen, K. In press. Spatial factor analysis: a new tool for estimating joint species distributions and correlations in species range. Methods Ecol. Evol.

For more details regarding estimating spatial variation via Gaussian random fields please see:

* Thorson, J.T., Shelton, A.O., Ward, E.J., and Skaug, H. In press. Geostatistical delta-generalized linear mixed models improve precision for estimated abundance indices for West Coast groundfishes. ICES J. Mar. Sci.
* Thorson, J.T., Skaug, H., Kristensen, K., Shelton, A.O., Ward, E.J., Harms, J., and Benante, J. In press. The importance of spatial models for estimating the strength of density dependence. Ecology. doi: http://dx.doi.org/10.1890/14-0739.1.

For more details regarding Template Model Builder please see:

* Kristensen, K. 2014. TMB: General random effect model builder tool inspired by ADMB. Available from https://github.com/kaskr/adcomp.
* Kristensen, K., Nielsen, A., Berg, C.W., and Skaug, H. In press. Template Model Builder TMB. J. Stat. Softw.

