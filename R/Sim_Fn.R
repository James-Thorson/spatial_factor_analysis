Sim_Fn <-
function( n_species, n_stations=200, n_factors=2, SpatialScale=0.1, SD_O=1.0, logMeanDens=3.0, Psi=NULL, Loc=NULL ){
  # Parameters
  if( is.null(Psi) ){
    Psi = matrix( rnorm(n_factors*n_species), nrow=n_factors, ncol=n_species)
    for(i in 1:nrow(Psi)) Psi[i,seq(from=1,to=i-1,length=i-1)] = 0
  }
  Beta = rep(logMeanDens, n_species)

  # Spatial model
  if( is.null(Loc) ) Loc = cbind( "x"=runif(n_stations, min=0,max=1), "y"=runif(n_stations, min=0,max=1) )
  model_O <- RMgauss(var=SD_O^2, scale=SpatialScale)

  # Simulate fields
  Omega = matrix(NA, ncol=n_factors, nrow=n_stations)
  for(i in 1:n_factors){
    Omega[,i] = RFsimulate(model = model_O, x=Loc[,'x'], y=Loc[,'y'])@data[,1]
  }
  ln_Y_exp = Omega%*%Psi + outer(rep(1,n_stations),Beta)

  # Simulate data
  Y = matrix(rpois( n_stations*n_species, lambda=exp(ln_Y_exp) ), ncol=n_species, byrow=FALSE)
  X = cbind( rep(1,nrow(Y)) )

  # Return stuff
  Sim_List = list("X"=X, "Y"=Y, "Psi"=Psi, "Loc"=Loc, "Omega"=Omega, "ln_Y_exp"=ln_Y_exp)
  return(Sim_List)
}
