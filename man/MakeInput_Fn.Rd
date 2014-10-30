\name{MakeInput_Fn}
\alias{MakeInput_Fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
This function generates the stochastic partial differential equation mesh and and other inputs necessary for parameter estimation via Template Model Builder
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MakeInput_Fn(Version, Y, X, Y_Report = NULL, LastRun = NULL, Loc, isPred, ErrorDist, ObsModel, VaryingKappa, n_factors, n_stations, Use_REML, Aniso)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Version}{
%%     ~~Describe \code{Version} here~~
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Y_Report}{
%%     ~~Describe \code{Y_Report} here~~
}
  \item{LastRun}{
%%     ~~Describe \code{LastRun} here~~
}
  \item{Loc}{
%%     ~~Describe \code{Loc} here~~
}
  \item{isPred}{
%%     ~~Describe \code{isPred} here~~
}
  \item{ErrorDist}{
%%     ~~Describe \code{ErrorDist} here~~
}
  \item{ObsModel}{
%%     ~~Describe \code{ObsModel} here~~
}
  \item{VaryingKappa}{
%%     ~~Describe \code{VaryingKappa} here~~
}
  \item{n_factors}{
%%     ~~Describe \code{n_factors} here~~
}
  \item{n_stations}{
%%     ~~Describe \code{n_stations} here~~
}
  \item{Use_REML}{
%%     ~~Describe \code{Use_REML} here~~
}
  \item{Aniso}{
%%     ~~Describe \code{Aniso} here~~
}
  \item{mesh}{
%%     ~~Describe \code{Aniso} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Version, Y, X, Y_Report = NULL, LastRun = NULL, Loc, 
    isPred, ErrorDist, ObsModel, VaryingKappa, n_factors, n_stations, 
    Use_REML, Aniso) 
{
    mesh = inla.mesh.create(Loc[, c("x", "y")], plot.delay = NULL, 
        refine = FALSE)
    spde = inla.spde2.matern(mesh, alpha = 2)
    Dset = 1:2
    TV = mesh$graph$tv
    V0 = mesh$loc[TV[, 1], Dset]
    V1 = mesh$loc[TV[, 2], Dset]
    V2 = mesh$loc[TV[, 3], Dset]
    E0 = V2 - V1
    E1 = V0 - V2
    E2 = V1 - V0
    TmpFn = function(Vec1, Vec2) abs(det(rbind(Vec1, Vec2)))
    Tri_Area = rep(NA, nrow(E0))
    for (i in 1:length(Tri_Area)) Tri_Area[i] = TmpFn(E0[i, ], 
        E1[i, ])/2
    NAind = ifelse(is.na(Y), 1, 0)
    if (is.null(Y_Report)) 
        Y_Report = array(0, dim = dim(Y))
    n_species = ncol(Y)
    n_fit = nrow(Y)
    if (Version == "spatial_factor_analysis_v18") 
        Data = list(Y = Y_fit, Y_Report = Y_Report, NAind = NAind_fit, 
            isPred = isPred, ErrorDist = as.integer(2), VaryingKappa = as.integer(VaryingKappa), 
            n_species = n_species, n_stations = n_fit, n_factors = n_factors, 
            meshidxloc = mesh$idx$loc - 1, n_p = ncol(X_fit), 
            X = X_fit, G0 = spde$param.inla$M0, G1 = spde$param.inla$M1, 
            G2 = spde$param.inla$M2)
    if (Version \%in\% c("spatial_factor_analysis_v19", "spatial_factor_analysis_v20")) 
        Data = list(Pen_Vec = c(0, 0), Y = Y_fit, Y_Report = Y_Report, 
            NAind = NAind_fit, isPred = isPred, ErrorDist = as.integer(2), 
            VaryingKappa = as.integer(VaryingKappa), n_species = n_species, 
            n_stations = n_fit, n_factors = n_factors, meshidxloc = mesh$idx$loc - 
                1, n_p = ncol(X_fit), X = X_fit, G0 = spde$param.inla$M0, 
            G1 = spde$param.inla$M1, G2 = spde$param.inla$M2)
    if (Version \%in\% c("spatial_factor_analysis_v21", "spatial_factor_analysis_v22")) 
        Data = list(Aniso = as.integer(Aniso), Pen_Vec = c(0, 
            0), Y = Y_fit, Y_Report = Y_Report, NAind = NAind_fit, 
            isPred = isPred, ErrorDist = as.integer(2), VaryingKappa = as.integer(VaryingKappa), 
            n_species = n_species, n_stations = n_fit, n_factors = n_factors, 
            n_i = mesh$n, meshidxloc = mesh$idx$loc - 1, n_p = ncol(X_fit), 
            X = X_fit, n_tri = nrow(mesh$graph$tv), Tri_Area = Tri_Area, 
            E0 = E0, E1 = E1, E2 = E2, TV = TV - 1, G0 = spde$param.inla$M0, 
            G0_inv = as(diag(1/diag(spde$param.inla$M0)), "dgTMatrix"), 
            G1 = spde$param.inla$M1, G2 = spde$param.inla$M2)
    if (StartValue == "Default" | is.null(LastRun)) {
        if (Version \%in\% c("spatial_factor_analysis_v18", "spatial_factor_analysis_v19")) 
            Params = list(beta = matrix(0, nrow = ncol(X), ncol = n_species, 
                byrow = TRUE), Psi_val = rep(1, n_factors * n_species - 
                n_factors * (n_factors - 1)/2), log_kappa_input = rep(0, 
                ifelse(VaryingKappa == 0, 1, n_factors)), ln_VarInfl = matrix(-2, 
                nrow = 3, ncol = n_species), Omega_input = matrix(rep(0, 
                spde$n.spde * n_factors), ncol = n_factors))
        if (Version \%in\% "spatial_factor_analysis_v20") 
            Params = list(beta = matrix(0, nrow = ncol(X), ncol = n_species, 
                byrow = TRUE), Psi_val = rep(1, n_factors * n_species - 
                n_factors * (n_factors - 1)/2), log_kappa_input = rep(0, 
                ifelse(VaryingKappa == 0, 1, n_factors)), ln_VarInfl_NB1 = NA, 
                ln_VarInfl_NB2 = NA, ln_VarInfl_ZI = NA, Omega_input = matrix(rep(0, 
                  spde$n.spde * n_factors), ncol = n_factors))
        if (Version \%in\% "spatial_factor_analysis_v21") 
            Params = list(ln_H_input = c(0, -10, 0), beta = matrix(0, 
                nrow = ncol(X), ncol = n_species, byrow = TRUE), 
                Psi_val = rep(1, n_factors * n_species - n_factors * 
                  (n_factors - 1)/2), log_kappa_input = rep(0, 
                  ifelse(VaryingKappa == 0, 1, n_factors)), ln_VarInfl_NB1 = NA, 
                ln_VarInfl_NB2 = NA, ln_VarInfl_ZI = NA, Omega_input = matrix(rep(0, 
                  spde$n.spde * n_factors), ncol = n_factors))
        if (Version \%in\% "spatial_factor_analysis_v22") 
            Params = list(ln_H_input = c(0, 0), beta = matrix(0, 
                nrow = ncol(X), ncol = n_species, byrow = TRUE), 
                Psi_val = rep(1, n_factors * n_species - n_factors * 
                  (n_factors - 1)/2), log_kappa_input = rep(0, 
                  ifelse(VaryingKappa == 0, 1, n_factors)), ln_VarInfl_NB1 = NA, 
                ln_VarInfl_NB2 = NA, ln_VarInfl_ZI = NA, Omega_input = matrix(rep(0, 
                  spde$n.spde * n_factors), ncol = n_factors))
        if (ObsModel == "Poisson") {
            Params[["ln_VarInfl_NB1"]] = rep(-20, n_species)
            Params[["ln_VarInfl_NB2"]] = rep(-20, n_species)
            Params[["ln_VarInfl_ZI"]] = rep(-20, n_species)
        }
        if (ObsModel == "NegBin1") {
            Params[["ln_VarInfl_NB1"]] = rep(-2, n_species)
            Params[["ln_VarInfl_NB2"]] = rep(-20, n_species)
            Params[["ln_VarInfl_ZI"]] = rep(-20, n_species)
        }
        if (ObsModel == "NegBin12") {
            Params[["ln_VarInfl_NB1"]] = rep(-2, n_species)
            Params[["ln_VarInfl_NB2"]] = rep(-2, n_species)
            Params[["ln_VarInfl_ZI"]] = rep(-20, n_species)
        }
    }
    if (StartValue == "Last_Run" & !is.null(LastRun)) {
        Par_last = LastRun$opt$par
        Par_best = LastRun$obj$env$last.par.best
        if (Version \%in\% c("spatial_factor_analysis_v18", "spatial_factor_analysis_v19")) 
            Params = list(beta = matrix(Par_last[which(names(Par_last) == 
                "beta")], nrow = ncol(X), ncol = n_species, byrow = TRUE), 
                Psi_val = c(Par_last[which(names(Par_last) == 
                  "Psi_val")], rep(1, n_species - n_factors + 
                  1)), log_kappa_input = c(Par_last[which(names(Par_last) == 
                  "log_kappa_input")], list(NULL, 1)[[ifelse(VaryingKappa == 
                  0, 1, 2)]]), ln_VarInfl = Par_last[which(names(Par_last) == 
                  "ln_VarInfl")], Omega_input = matrix(c(Par_best[which(names(Par_best) == 
                  "Omega_input")], rep(0, spde$n.spde)), ncol = n_factors))
        if (Version \%in\% "spatial_factor_analysis_v20") 
            Params = list(beta = matrix(Par_last[which(names(Par_last) == 
                "beta")], nrow = ncol(X), ncol = n_species, byrow = TRUE), 
                Psi_val = c(Par_last[which(names(Par_last) == 
                  "Psi_val")], rep(1, n_species - n_factors + 
                  1)), log_kappa_input = c(Par_last[which(names(Par_last) == 
                  "log_kappa_input")], list(NULL, 1)[[ifelse(VaryingKappa == 
                  0, 1, 2)]]), ln_VarInfl_NB1 = NA, ln_VarInfl_NB2 = NA, 
                ln_VarInfl_ZI = NA, Omega_input = matrix(c(Par_best[which(names(Par_best) == 
                  "Omega_input")], rep(0, spde$n.spde)), ncol = n_factors))
        if (Version \%in\% "spatial_factor_analysis_v21") 
            Params = list(ln_H_input = Par_last[which(names(Par_last) == 
                "ln_H_input")], beta = matrix(Par_last[which(names(Par_last) == 
                "beta")], nrow = ncol(X), ncol = n_species, byrow = TRUE), 
                Psi_val = c(Par_last[which(names(Par_last) == 
                  "Psi_val")], rep(1, n_species - n_factors + 
                  1)), log_kappa_input = c(Par_last[which(names(Par_last) == 
                  "log_kappa_input")], list(NULL, 1)[[ifelse(VaryingKappa == 
                  0, 1, 2)]]), ln_VarInfl_NB1 = NA, ln_VarInfl_NB2 = NA, 
                ln_VarInfl_ZI = NA, Omega_input = matrix(c(Par_best[which(names(Par_best) == 
                  "Omega_input")], rep(0, spde$n.spde)), ncol = n_factors))
        if (Version \%in\% "spatial_factor_analysis_v22") 
            Params = list(ln_H_input = Par_best[which(names(Par_best) == 
                "ln_H_input")], beta = matrix(Par_best[which(names(Par_best) == 
                "beta")], nrow = ncol(X), ncol = n_species, byrow = TRUE), 
                Psi_val = c(Par_last[which(names(Par_last) == 
                  "Psi_val")], rep(1, n_species - n_factors + 
                  1)), log_kappa_input = c(Par_last[which(names(Par_last) == 
                  "log_kappa_input")], list(NULL, 1)[[ifelse(VaryingKappa == 
                  0, 1, 2)]]), ln_VarInfl_NB1 = NA, ln_VarInfl_NB2 = NA, 
                ln_VarInfl_ZI = NA, Omega_input = matrix(c(Par_best[which(names(Par_best) == 
                  "Omega_input")], rep(0, spde$n.spde)), ncol = n_factors))
        if (ObsModel == "Poisson") {
            Params[["ln_VarInfl_NB1"]] = rep(-20, n_species)
            Params[["ln_VarInfl_NB2"]] = rep(-20, n_species)
            Params[["ln_VarInfl_ZI"]] = rep(-20, n_species)
        }
        if (ObsModel == "NegBin1") {
            Params[["ln_VarInfl_NB1"]] = Par_best[which(names(Par_best) == 
                "ln_VarInfl_NB1")]
            Params[["ln_VarInfl_NB2"]] = rep(-20, n_species)
            Params[["ln_VarInfl_ZI"]] = rep(-20, n_species)
        }
        if (ObsModel == "NegBin12") {
            Params[["ln_VarInfl_NB1"]] = Par_best[which(names(Par_best) == 
                "ln_VarInfl_NB1")]
            Params[["ln_VarInfl_NB2"]] = Par_best[which(names(Par_best) == 
                "ln_VarInfl_NB2")]
            Params[["ln_VarInfl_ZI"]] = rep(-20, n_species)
        }
        if (Version \%in\% c("spatial_factor_analysis_v20", "spatial_factor_analysis_v21")) {
            if (Aniso == 0) 
                Params[["ln_H_input"]] = c(0, -10, 0)
        }
        if (Version \%in\% c("spatial_factor_analysis_v22")) {
            if (Aniso == 0) 
                Params[["ln_H_input"]] = c(0, 0)
        }
    }
    Map = list()
    if (ObsModel == "Poisson") {
        Map[["ln_VarInfl_NB1"]] = factor(rep(NA, n_species))
        Map[["ln_VarInfl_NB2"]] = factor(rep(NA, n_species))
        Map[["ln_VarInfl_ZI"]] = factor(rep(NA, n_species))
    }
    if (ObsModel == "NegBin1") {
        Map[["ln_VarInfl_ZI"]] = factor(rep(NA, n_species))
        Map[["ln_VarInfl_NB2"]] = factor(rep(NA, n_species))
    }
    if (ObsModel == "NegBin12") {
        Map[["ln_VarInfl_ZI"]] = factor(rep(NA, n_species))
    }
    if (Version \%in\% c("spatial_factor_analysis_v18", "spatial_factor_analysis_v19", 
        "spatial_factor_analysis_v20", "spatial_factor_analysis_v21")) {
        if (Aniso == 0) 
            Map[["ln_H_input"]] = factor(rep(NA, 3))
        if (Aniso == 1 & VaryingKappa == 0) 
            Map[["log_kappa_input"]] = factor(NA)
    }
    if (Version \%in\% c("spatial_factor_analysis_v22")) {
        if (Aniso == 0) 
            Map[["ln_H_input"]] = factor(rep(NA, 2))
    }
    Random = c("Omega_input")
    if (Use_REML == TRUE) 
        Random = c(Random, "beta", "ln_VarInfl_NB1", "ln_VarInfl_NB2", 
            "ln_VarInfl_ZI")
    Return = list(Map = Map, Random = Random, Params = Params, 
        Data = Data, spde = spde)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
