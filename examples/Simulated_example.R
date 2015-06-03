
# Step 1: Install TMB and INLA packages
# Step 2: Install SpatialFA package
devtools::install_github("james-thorson/spatial_factor_analysis")

# Load libraries
library(TMB)
library(INLA)
library(SpatialFA)

# Set working directory to location of TMB code
  TmbFile = paste0(system.file("executables", package="SpatialFA"),"/")
  setwd(TmbFile)

# Simulation settings
  n_factors_true = 2
  n_species = 10
  n_stations = 100

# Estimation settings
  Predictive = 1 # 0: None;  1: Predictive mean;  2: Predictive probability
  VaryingKappa = 0 # 0: No; 1: Yes
  Aniso = 0 # 0: No; 1: Yes
  Version = "spatial_factor_analysis_v24"
  n_factors = 2
  ObsModel = c("Poisson", "NegBin0", "NegBin1", "NegBin12", "Lognorm_Pois")[5]

# Compile if necessary
  compile( paste0(Version,".cpp") )

# Simulate data
  Data_List = Sim_Fn( n_species=n_species, n_stations=n_stations, n_factors=n_factors_true )
  Y = Data_List[["Y"]]
  X = Data_List[["X"]]
  Loc = Data_List[["Loc"]]

# Create SPDE mesh
  mesh = inla.mesh.create( Loc[,c('x','y')], plot.delay=NULL, refine=FALSE)
  spde = inla.spde2.matern(mesh,alpha=2)
  spde$param.inla$M0_inv = as(diag(1/diag(spde$param.inla$M0)),"dgTMatrix")

# Generate inputs
  InputList = MakeInput_Fn(Y=Y, X=X, n_factors=n_factors, n_stations=n_stations, Loc=Loc, Version=Version, isPred=rep(0,nrow(Y)), ObsModel=ObsModel, VaryingKappa=VaryingKappa, Use_REML=FALSE, Aniso=Aniso, mesh=mesh, spde=spde)

# Link DLL
  dyn.load( dynlib(Version) )                                                         # log_tau=0.0,

# Initialization
  obj <- MakeADFun(data=InputList[['Data']], parameters=InputList[['Params']], random=InputList[['Random']], map=InputList[['Map']], hessian=FALSE, inner.control=list(maxit=1000) )
  obj$control <- c( obj$control, list(trace=1, parscale=1, REPORT=1, reltol=1e-12, maxit=100) )
  obj$env$inner.control <- c(obj$env$inner.control, list("step.tol"=1e-8, "tol10"=1e-6, "grad.tol"=1e-8) )

# First run
  Init = obj$fn( obj$par )
  Initial_gradient = obj$gr( obj$par )

# Run model
  opt = nlminb(start=obj$env$last.par.best[-c(obj$env$random)], objective=obj$fn, gradient=obj$gr, control=list(eval.max=1e4, iter.max=1e4, trace=1, rel.tol=1e-14) )

# Summarize
  SD = try( sdreport(obj) )
  Report = obj$report()

# Compare with simulated values
  # Psi (loadings matrix)
  Report$Psi; Data_List$Psi
  # Omega (prediction for each species)
  par(mfrow=c(1,n_factors))
  for(i in 1:n_factors) plot(y=Report$Omega[1:n_stations,i], x=Data_List$Omega[,i], xlab="True", ylab="Estimated", main=paste("Factor",i))
