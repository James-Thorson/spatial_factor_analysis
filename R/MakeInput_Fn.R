MakeInput_Fn <-
function(Version, Y, X, Y_Report=NULL, LastRun=NULL, Loc, isPred, ObsModel, VaryingKappa, n_factors, n_stations, Use_REML, Aniso, mesh, spde){

  # Pre-processing in R:  Assume 2-Dimensional
  Dset = 1:2
  # Triangle info
  TV = mesh$graph$tv       # Triangle to vertex indexing
  V0 = mesh$loc[TV[,1],Dset]   # V = vertices for each triangle
  V1 = mesh$loc[TV[,2],Dset]
  V2 = mesh$loc[TV[,3],Dset]
  E0 = V2 - V1                      # E = edge for each triangle
  E1 = V0 - V2
  E2 = V1 - V0
  # Calculate Areas
  TmpFn = function(Vec1,Vec2) abs(det( rbind(Vec1,Vec2) ))
  Tri_Area = rep(NA, nrow(E0))
  for(i in 1:length(Tri_Area)) Tri_Area[i] = TmpFn( E0[i,],E1[i,] )/2   # T = area of each triangle

  # Other pre-processing
  NAind = ifelse( is.na(Y), 1, 0)
  if(is.null(Y_Report)) Y_Report = array(0, dim=dim(Y))
  n_species = ncol(Y)
  n_fit = nrow(Y)

  # Data
  ErrorDist = as.integer(switch(ObsModel, "Poisson"=0, "NegBin0"=1, "NegBin1"=1, "NegBin2"=1, "NegBin12"=1, "Lognorm_Pois"=0))
  if(Version=="spatial_factor_analysis_v18") Data = list(Y=Y, Y_Report=Y_Report, NAind=NAind, isPred=isPred, ErrorDist=ErrorDist, VaryingKappa=as.integer(VaryingKappa), n_species=n_species, n_stations=n_fit, n_factors=n_factors, meshidxloc=mesh$idx$loc-1, n_p=ncol(X), X=X, G0=spde$param.inla$M0, G1=spde$param.inla$M1, G2=spde$param.inla$M2)
  if(Version %in% c("spatial_factor_analysis_v19","spatial_factor_analysis_v20") ) Data = list(Pen_Vec=c(0,0), Y=Y, Y_Report=Y_Report, NAind=NAind, isPred=isPred, ErrorDist=ErrorDist, VaryingKappa=as.integer(VaryingKappa), n_species=n_species, n_stations=n_fit, n_factors=n_factors, meshidxloc=mesh$idx$loc-1, n_p=ncol(X), X=X, G0=spde$param.inla$M0, G1=spde$param.inla$M1, G2=spde$param.inla$M2)
  if(Version %in% c("spatial_factor_analysis_v21","spatial_factor_analysis_v22") ) Data = list(Aniso=as.integer(Aniso), Pen_Vec=c(0,0), Y=Y, Y_Report=Y_Report, NAind=NAind, isPred=isPred, ErrorDist=ErrorDist, VaryingKappa=as.integer(VaryingKappa), n_species=n_species, n_stations=n_fit, n_factors=n_factors, n_i=mesh$n, meshidxloc=mesh$idx$loc-1, n_p=ncol(X), X=X, n_tri=nrow(mesh$graph$tv), Tri_Area=Tri_Area, E0=E0, E1=E1, E2=E2, TV=TV-1, G0=spde$param.inla$M0, G0_inv=spde$param.inla$M0_inv, G1=spde$param.inla$M1, G2=spde$param.inla$M2)
  if(Version %in% c("spatial_factor_analysis_v23","spatial_factor_analysis_v24") ) Data = list(Aniso=as.integer(Aniso), Options_Vec=NA, Pen_Vec=c(0,0), Y=Y, Y_Report=Y_Report, NAind=NAind, isPred=isPred, ErrorDist=ErrorDist, VaryingKappa=as.integer(VaryingKappa), n_species=n_species, n_stations=n_fit, n_factors=n_factors, n_i=mesh$n, meshidxloc=mesh$idx$loc-1, n_p=ncol(X), X=X, n_tri=nrow(mesh$graph$tv), Tri_Area=Tri_Area, E0=E0, E1=E1, E2=E2, TV=TV-1, G0=spde$param.inla$M0, G0_inv=spde$param.inla$M0_inv, G1=spde$param.inla$M1, G2=spde$param.inla$M2)
  if(Version %in% c("spatial_factor_analysis_v23","spatial_factor_analysis_v24") ){
    if(ObsModel!="Lognorm_Pois") Data[["Options_Vec"]] = as.integer(0)
    if(ObsModel=="Lognorm_Pois" & Use_REML==TRUE) Data[["Options_Vec"]] = as.integer(2)
    if(ObsModel=="Lognorm_Pois" & Use_REML==FALSE) Data[["Options_Vec"]] = as.integer(1)
  }
  # Parameters
  if(is.null(LastRun) || StartValue=="Default" ){
    if(Version %in% c("spatial_factor_analysis_v18","spatial_factor_analysis_v19") ) Params = list( beta=matrix(0,nrow=ncol(X),ncol=n_species,byrow=TRUE), Psi_val=rep(1,n_factors*n_species-n_factors*(n_factors-1)/2), log_kappa_input=rep(0,ifelse(VaryingKappa==0,1,n_factors)), ln_VarInfl=matrix(-2,nrow=3,ncol=n_species), Omega_input=matrix(rep(0,spde$n.spde*n_factors),ncol=n_factors))
    if(Version %in% "spatial_factor_analysis_v20") Params = list( beta=matrix(0,nrow=ncol(X),ncol=n_species,byrow=TRUE), Psi_val=rep(1,n_factors*n_species-n_factors*(n_factors-1)/2), log_kappa_input=rep(0,ifelse(VaryingKappa==0,1,n_factors)), ln_VarInfl_NB1=NA, ln_VarInfl_NB2=NA, ln_VarInfl_ZI=NA, Omega_input=matrix(rep(0,spde$n.spde*n_factors),ncol=n_factors))
    if(Version %in% "spatial_factor_analysis_v21") Params = list( ln_H_input=c(0,-10,0), beta=matrix(0,nrow=ncol(X),ncol=n_species,byrow=TRUE), Psi_val=rep(1,n_factors*n_species-n_factors*(n_factors-1)/2), log_kappa_input=rep(0,ifelse(VaryingKappa==0,1,n_factors)), ln_VarInfl_NB1=NA, ln_VarInfl_NB2=NA, ln_VarInfl_ZI=NA, Omega_input=matrix(rep(0,spde$n.spde*n_factors),ncol=n_factors))
    if(Version %in% "spatial_factor_analysis_v22") Params = list( ln_H_input=c(0,0), beta=matrix(0,nrow=ncol(X),ncol=n_species,byrow=TRUE), Psi_val=rep(1,n_factors*n_species-n_factors*(n_factors-1)/2), log_kappa_input=rep(0,ifelse(VaryingKappa==0,1,n_factors)), ln_VarInfl_NB1=NA, ln_VarInfl_NB2=NA, ln_VarInfl_ZI=NA, Omega_input=matrix(rep(0,spde$n.spde*n_factors),ncol=n_factors))
    if(Version %in% "spatial_factor_analysis_v23") Params = list( ln_H_input=c(0,0), beta=matrix(0,nrow=ncol(X),ncol=n_species,byrow=TRUE), Psi_val=rep(1,n_factors*n_species-n_factors*(n_factors-1)/2), log_kappa_input=rep(0,ifelse(VaryingKappa==0,1,n_factors)), ln_VarInfl_NB1=NA, ln_VarInfl_NB2=NA, ln_VarInfl_ZI=NA, ln_VarInfl_Lognorm=NA, Omega_input=matrix(rep(0,spde$n.spde*n_factors),ncol=n_factors), Lognorm_input=array(0,dim=dim(Data[["Y"]])) )
    if(Version %in% "spatial_factor_analysis_v24") Params = list( ln_H_input=c(0,0), beta=matrix(0,nrow=ncol(X),ncol=n_species,byrow=TRUE), Psi_val=rep(1,n_factors*n_species-n_factors*(n_factors-1)/2), log_kappa_input=rep(0,ifelse(VaryingKappa==0,1,n_factors)), ln_VarInfl_NB0=NA, ln_VarInfl_NB1=NA, ln_VarInfl_NB2=NA, ln_VarInfl_ZI=NA, ln_VarInfl_Lognorm=NA, Omega_input=matrix(rep(0,spde$n.spde*n_factors),ncol=n_factors), Lognorm_input=array(0,dim=dim(Data[["Y"]])) )
    if(ObsModel=="Poisson"){
      if('ln_VarInfl_NB0'%in%names(Params)) Params[['ln_VarInfl_NB0']] = rep(-20,n_species) 
      Params[['ln_VarInfl_NB1']] = rep(-20,n_species)
      Params[['ln_VarInfl_NB2']] = rep(-20,n_species)
      Params[['ln_VarInfl_ZI']] = rep(-20,n_species)
      if('ln_VarInfl_Lognorm'%in%names(Params)) Params[['ln_VarInfl_Lognorm']] = rep(-20,n_species)
    }
    if(ObsModel=="NegBin0"){
      if('ln_VarInfl_NB0'%in%names(Params)) Params[['ln_VarInfl_NB0']] = rep(-2,n_species)
      Params[['ln_VarInfl_NB1']] = rep(-20,n_species)
      Params[['ln_VarInfl_NB2']] = rep(-20,n_species)
      Params[['ln_VarInfl_ZI']] = rep(-20,n_species)
      if('ln_VarInfl_Lognorm'%in%names(Params)) Params[['ln_VarInfl_Lognorm']] = rep(-20,n_species)
    }
    if(ObsModel=="NegBin1"){
      if('ln_VarInfl_NB0'%in%names(Params)) Params[['ln_VarInfl_NB0']] = rep(-20,n_species) 
      Params[['ln_VarInfl_NB1']] = rep(-2,n_species)
      Params[['ln_VarInfl_NB2']] = rep(-20,n_species)
      Params[['ln_VarInfl_ZI']] = rep(-20,n_species)
      if('ln_VarInfl_Lognorm'%in%names(Params)) Params[['ln_VarInfl_Lognorm']] = rep(-20,n_species)
    }
    if(ObsModel=="NegBin2"){
      if('ln_VarInfl_NB0'%in%names(Params)) Params[['ln_VarInfl_NB0']] = rep(-20,n_species) 
      Params[['ln_VarInfl_NB1']] = rep(-20,n_species)
      Params[['ln_VarInfl_NB2']] = rep(-2,n_species)
      Params[['ln_VarInfl_ZI']] = rep(-20,n_species)
      if('ln_VarInfl_Lognorm'%in%names(Params)) Params[['ln_VarInfl_Lognorm']] = rep(-20,n_species)
    }
    if(ObsModel=="NegBin12"){
      if('ln_VarInfl_NB0'%in%names(Params)) Params[['ln_VarInfl_NB0']] = rep(-20,n_species) 
      Params[['ln_VarInfl_NB1']] = rep(-2,n_species)
      Params[['ln_VarInfl_NB2']] = rep(-2,n_species)
      Params[['ln_VarInfl_ZI']] = rep(-20,n_species)
      if('ln_VarInfl_Lognorm'%in%names(Params)) Params[['ln_VarInfl_Lognorm']] = rep(-20,n_species)
    }
    if(ObsModel=="Lognorm_Pois"){
      if('ln_VarInfl_NB0'%in%names(Params)) Params[['ln_VarInfl_NB0']] = rep(-20,n_species) 
      Params[['ln_VarInfl_NB1']] = rep(-20,n_species)
      Params[['ln_VarInfl_NB2']] = rep(-20,n_species)
      Params[['ln_VarInfl_ZI']] = rep(-20,n_species)
      if('ln_VarInfl_Lognorm'%in%names(Params)) Params[['ln_VarInfl_Lognorm']] = rep(-5,n_species)
    }
  }
  if(!is.null(LastRun) && StartValue=="Last_Run"){
    Par_last = LastRun$opt$par
    Par_best = LastRun$ParBest
    # names(Save)
    if(Version %in% c("spatial_factor_analysis_v18","spatial_factor_analysis_v19") ) Params = list( beta=matrix(Par_last[which(names(Par_last)=="beta")],nrow=ncol(X),ncol=n_species,byrow=TRUE), Psi_val=c(Par_last[which(names(Par_last)=="Psi_val")],rep(1,n_species-n_factors+1)), log_kappa_input=c(Par_last[which(names(Par_last)=="log_kappa_input")],list(NULL,1)[[ifelse(VaryingKappa==0,1,2)]]), ln_VarInfl=Par_last[which(names(Par_last)=="ln_VarInfl")], Omega_input=matrix(c(Par_best[which(names(Par_best)=="Omega_input")],rep(0,spde$n.spde)),ncol=n_factors))
    if(Version %in% "spatial_factor_analysis_v20") Params = list( beta=matrix(Par_last[which(names(Par_last)=="beta")],nrow=ncol(X),ncol=n_species,byrow=TRUE), Psi_val=c(Par_last[which(names(Par_last)=="Psi_val")],rep(1,n_species-n_factors+1)), log_kappa_input=c(Par_last[which(names(Par_last)=="log_kappa_input")],list(NULL,1)[[ifelse(VaryingKappa==0,1,2)]]), ln_VarInfl_NB1=NA, ln_VarInfl_NB2=NA, ln_VarInfl_ZI=NA, Omega_input=matrix(c(Par_best[which(names(Par_best)=="Omega_input")],rep(0,spde$n.spde)),ncol=n_factors))
    if(Version %in% "spatial_factor_analysis_v21") Params = list( ln_H_input=Par_last[which(names(Par_last)=="ln_H_input")], beta=matrix(Par_last[which(names(Par_last)=="beta")],nrow=ncol(X),ncol=n_species,byrow=TRUE), Psi_val=c(Par_last[which(names(Par_last)=="Psi_val")],rep(1,n_species-n_factors+1)), log_kappa_input=c(Par_last[which(names(Par_last)=="log_kappa_input")],list(NULL,1)[[ifelse(VaryingKappa==0,1,2)]]), ln_VarInfl_NB1=NA, ln_VarInfl_NB2=NA, ln_VarInfl_ZI=NA, Omega_input=matrix(c(Par_best[which(names(Par_best)=="Omega_input")],rep(0,spde$n.spde)),ncol=n_factors))
    if(Version %in% "spatial_factor_analysis_v22") Params = list( ln_H_input=Par_best[which(names(Par_best)=="ln_H_input")], beta=matrix(Par_best[which(names(Par_best)=="beta")],nrow=ncol(X),ncol=n_species,byrow=TRUE), Psi_val=c(Par_last[which(names(Par_last)=="Psi_val")],rep(1,n_species-n_factors+1)), log_kappa_input=c(Par_last[which(names(Par_last)=="log_kappa_input")],list(NULL,1)[[ifelse(VaryingKappa==0,1,2)]]), ln_VarInfl_NB1=NA, ln_VarInfl_NB2=NA, ln_VarInfl_ZI=NA, Omega_input=matrix(c(Par_best[which(names(Par_best)=="Omega_input")],rep(0,spde$n.spde)),ncol=n_factors))
    if(Version %in% "spatial_factor_analysis_v23") Params = list( ln_H_input=Par_best[which(names(Par_best)=="ln_H_input")], beta=matrix(Par_best[which(names(Par_best)=="beta")],nrow=ncol(X),ncol=n_species,byrow=TRUE), Psi_val=c(Par_last[which(names(Par_last)=="Psi_val")],rep(1,n_species-n_factors+1)), log_kappa_input=c(Par_last[which(names(Par_last)=="log_kappa_input")],list(NULL,1)[[ifelse(VaryingKappa==0,1,2)]]), ln_VarInfl_NB1=NA, ln_VarInfl_NB2=NA, ln_VarInfl_ZI=NA, ln_VarInfl_Lognorm=NA, Omega_input=matrix(c(Par_best[which(names(Par_best)=="Omega_input")],rep(0,spde$n.spde)),ncol=n_factors), Lognorm_input=matrix(Par_best[which(names(Par_best)=="Lognorm_input")],ncol=n_species))
    if(Version %in% c("spatial_factor_analysis_v24")) Params = list( ln_H_input=Par_best[which(names(Par_best)=="ln_H_input")], beta=matrix(Par_best[which(names(Par_best)=="beta")],nrow=ncol(X),ncol=n_species,byrow=TRUE), Psi_val=c(Par_last[which(names(Par_last)=="Psi_val")],rep(1,n_species-n_factors+1)), log_kappa_input=c(Par_last[which(names(Par_last)=="log_kappa_input")],list(NULL,1)[[ifelse(VaryingKappa==0,1,2)]]), ln_VarInfl_NB0=NA, ln_VarInfl_NB1=NA, ln_VarInfl_NB2=NA, ln_VarInfl_ZI=NA, ln_VarInfl_Lognorm=NA, Omega_input=matrix(c(Par_best[which(names(Par_best)=="Omega_input")],rep(0,spde$n.spde)),ncol=n_factors), Lognorm_input=matrix(Par_best[which(names(Par_best)=="Lognorm_input")],ncol=n_species))
    if(ObsModel=="Poisson"){
      if('ln_VarInfl_NB0'%in%names(Params)) Params[['ln_VarInfl_NB0']] = rep(-20,n_species) 
      Params[['ln_VarInfl_NB1']] = rep(-20,n_species)
      Params[['ln_VarInfl_NB2']] = rep(-20,n_species)
      Params[['ln_VarInfl_ZI']] = rep(-20,n_species)
      if('ln_VarInfl_Lognorm'%in%names(Params)) Params[['ln_VarInfl_Lognorm']] = rep(-20,n_species)
    }
    if(ObsModel=="NegBin0"){
      if('ln_VarInfl_NB0'%in%names(Params)) Params[['ln_VarInfl_NB0']] = Par_best[which(names(Par_best)=="ln_VarInfl_NB0")]
      Params[['ln_VarInfl_NB1']] = rep(-20,n_species)
      Params[['ln_VarInfl_NB2']] = rep(-20,n_species)
      Params[['ln_VarInfl_ZI']] = rep(-20,n_species)
      if('ln_VarInfl_Lognorm'%in%names(Params)) Params[['ln_VarInfl_Lognorm']] = rep(-20,n_species)
    }
    if(ObsModel=="NegBin1"){
      if('ln_VarInfl_NB0'%in%names(Params)) Params[['ln_VarInfl_NB0']] = rep(-20,n_species) 
      Params[['ln_VarInfl_NB1']] = Par_best[which(names(Par_best)=="ln_VarInfl_NB1")]
      Params[['ln_VarInfl_NB2']] = rep(-20,n_species)
      Params[['ln_VarInfl_ZI']] = rep(-20,n_species)
      if('ln_VarInfl_Lognorm'%in%names(Params)) Params[['ln_VarInfl_Lognorm']] = rep(-20,n_species)
    }
    if(ObsModel=="NegBin2"){
      if('ln_VarInfl_NB0'%in%names(Params)) Params[['ln_VarInfl_NB0']] = rep(-20,n_species) 
      Params[['ln_VarInfl_NB1']] = rep(-20,n_species)
      Params[['ln_VarInfl_NB2']] = Par_best[which(names(Par_best)=="ln_VarInfl_NB2")]
      Params[['ln_VarInfl_ZI']] = rep(-20,n_species)
      if('ln_VarInfl_Lognorm'%in%names(Params)) Params[['ln_VarInfl_Lognorm']] = rep(-20,n_species)
    }
    if(ObsModel=="NegBin12"){
      if('ln_VarInfl_NB0'%in%names(Params)) Params[['ln_VarInfl_NB0']] = rep(-20,n_species) 
      Params[['ln_VarInfl_NB1']] = Par_best[which(names(Par_best)=="ln_VarInfl_NB1")]
      Params[['ln_VarInfl_NB2']] = Par_best[which(names(Par_best)=="ln_VarInfl_NB2")]
      Params[['ln_VarInfl_ZI']] = rep(-20,n_species)
      if('ln_VarInfl_Lognorm'%in%names(Params)) Params[['ln_VarInfl_Lognorm']] = rep(-20,n_species)
    }
    if(ObsModel=="Lognorm_Pois"){
      if('ln_VarInfl_NB0'%in%names(Params)) Params[['ln_VarInfl_NB0']] = rep(-20,n_species) 
      Params[['ln_VarInfl_NB1']] = rep(-20,n_species)
      Params[['ln_VarInfl_NB2']] = rep(-20,n_species)
      Params[['ln_VarInfl_ZI']] = rep(-20,n_species)
      if('ln_VarInfl_Lognorm'%in%names(Params)) Params[['ln_VarInfl_Lognorm']] = Par_best[which(names(Par_best)=="ln_VarInfl_Lognorm")]
    }
    if(Version %in% c("spatial_factor_analysis_v20","spatial_factor_analysis_v21") ){
      if(Aniso==0) Params[['ln_H_input']] = c(0,-10,0)
    }
    if(Version %in% c("spatial_factor_analysis_v22","spatial_factor_analysis_v23","spatial_factor_analysis_v24") ){
      if(Aniso==0) Params[['ln_H_input']] = c(0,0)
    }
  }
  # Fix parameters
  Map = list()
  if(ObsModel=="Poisson"){
    Map[['ln_VarInfl_NB0']] = factor( rep(NA,n_species) )
    Map[['ln_VarInfl_NB1']] = factor( rep(NA,n_species) )
    Map[['ln_VarInfl_NB2']] = factor( rep(NA,n_species) )
    Map[['ln_VarInfl_ZI']] = factor( rep(NA,n_species) )
    Map[['ln_VarInfl_Lognorm']] = factor( rep(NA,n_species) )
    Map[['Lognorm_input']] = factor( array(NA,dim=dim(Params[['Lognorm_input']])) )
  }
  if(ObsModel=="NegBin0"){
    Map[['ln_VarInfl_ZI']] = factor( rep(NA,n_species) )
    Map[['ln_VarInfl_NB1']] = factor( rep(NA,n_species) )
    Map[['ln_VarInfl_NB2']] = factor( rep(NA,n_species) )
    Map[['ln_VarInfl_Lognorm']] = factor( rep(NA,n_species) )
    Map[['Lognorm_input']] = factor( array(NA,dim=dim(Params[['Lognorm_input']])) )
  }
  if(ObsModel=="NegBin1"){
    Map[['ln_VarInfl_NB0']] = factor( rep(NA,n_species) )
    Map[['ln_VarInfl_ZI']] = factor( rep(NA,n_species) )
    Map[['ln_VarInfl_NB2']] = factor( rep(NA,n_species) )
    Map[['ln_VarInfl_Lognorm']] = factor( rep(NA,n_species) )
    Map[['Lognorm_input']] = factor( array(NA,dim=dim(Params[['Lognorm_input']])) )
  }
  if(ObsModel=="NegBin2"){
    Map[['ln_VarInfl_NB0']] = factor( rep(NA,n_species) )
    Map[['ln_VarInfl_NB1']] = factor( rep(NA,n_species) )
    Map[['ln_VarInfl_ZI']] = factor( rep(NA,n_species) )
    Map[['ln_VarInfl_Lognorm']] = factor( rep(NA,n_species) )
    Map[['Lognorm_input']] = factor( array(NA,dim=dim(Params[['Lognorm_input']])) )
  }
  if(ObsModel=="NegBin12"){
    Map[['ln_VarInfl_NB0']] = factor( rep(NA,n_species) )
    Map[['ln_VarInfl_ZI']] = factor( rep(NA,n_species) )
    Map[['ln_VarInfl_Lognorm']] = factor( rep(NA,n_species) )
    Map[['Lognorm_input']] = factor( array(NA,dim=dim(Params[['Lognorm_input']])) )
  }
  if(ObsModel=="Lognorm_Pois"){
    Map[['ln_VarInfl_NB0']] = factor( rep(NA,n_species) )
    Map[['ln_VarInfl_NB1']] = factor( rep(NA,n_species) )
    Map[['ln_VarInfl_NB2']] = factor( rep(NA,n_species) )
    Map[['ln_VarInfl_ZI']] = factor( rep(NA,n_species) )
  }
  if(Version %in% c("spatial_factor_analysis_v18","spatial_factor_analysis_v19","spatial_factor_analysis_v20","spatial_factor_analysis_v21") ){
    if(Aniso==0) Map[['ln_H_input']] = factor( rep(NA,3) )
    if(Aniso==1 & VaryingKappa==0) Map[['log_kappa_input']] = factor( NA )
  }
  if(Version %in% c("spatial_factor_analysis_v22","spatial_factor_analysis_v23","spatial_factor_analysis_v24") ){
    if(Aniso==0) Map[['ln_H_input']] = factor( rep(NA,2) )
  }
  # Declare random
  Random = c("Omega_input")
  if(Version %in% c("spatial_factor_analysis_v23","spatial_factor_analysis_v24") ) Random = c(Random, "Lognorm_input")
  if(Use_REML==TRUE) Random = c(Random, "beta", "ln_VarInfl_NB0", "ln_VarInfl_NB1", "ln_VarInfl_NB2", "ln_VarInfl_ZI", "ln_VarInfl_Lognorm")
    
  # Return
  Return = list("Map"=Map, "Random"=Random, "Params"=Params, "Data"=Data, "spde"=spde)
}
