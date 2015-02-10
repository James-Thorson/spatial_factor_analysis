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

