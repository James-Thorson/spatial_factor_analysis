// Space time 
#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_FACTOR(Options_Vec);
  DATA_VECTOR(Pen_Vec);
  DATA_MATRIX(Y);       	// Responses
  DATA_MATRIX(Y_Report);
  DATA_MATRIX(NAind);		// 1 = Y is NA, 0 = is not NA
  DATA_VECTOR(isPred);
  DATA_INTEGER(ErrorDist);
  DATA_INTEGER(VaryingKappa); // 0: No; 1:Yes
  DATA_INTEGER(Aniso); // 0: No; 1:Yes
  DATA_INTEGER(n_species)	// Number of stations = 24
  DATA_INTEGER(n_stations)	// Number of stations = 24
  DATA_INTEGER(n_factors)	// Number of stations = 24
  DATA_INTEGER(n_i) // Number of vertices in mesh
  DATA_FACTOR(meshidxloc);	// Pointers into random effects vector x
  DATA_INTEGER(n_p)          	// # columns X
  DATA_MATRIX(X);		// Covariate desing matrix

  // Aniso objects
  DATA_INTEGER(n_tri);      //  Number of triangles
  DATA_VECTOR(Tri_Area);
  DATA_MATRIX(E0);
  DATA_MATRIX(E1);
  DATA_MATRIX(E2);
  DATA_FACTOR(TV);        //  This already includes the -1 for indexing in C

  // SPDE objects
  DATA_SPARSE_MATRIX(G0);
  DATA_SPARSE_MATRIX(G0_inv);
  DATA_SPARSE_MATRIX(G1);
  DATA_SPARSE_MATRIX(G2);

  PARAMETER_VECTOR(ln_H_input);
  PARAMETER_MATRIX(beta);
  PARAMETER_VECTOR(Psi_val);
  //PARAMETER(log_tau);
  PARAMETER_VECTOR(log_kappa_input);
  PARAMETER_VECTOR(ln_VarInfl_NB0);
  PARAMETER_VECTOR(ln_VarInfl_NB1);
  PARAMETER_VECTOR(ln_VarInfl_NB2);
  PARAMETER_VECTOR(ln_VarInfl_ZI);
  PARAMETER_VECTOR(ln_VarInfl_Lognorm);
  PARAMETER_MATRIX(Omega_input);
  PARAMETER_MATRIX(Lognorm_input);

  using namespace density;
  int i,j;
  Type g=0;
  
  matrix<Type> H(2,2);
  H(0,0) = exp(ln_H_input(0));
  H(1,0) = ln_H_input(1);
  H(0,1) = ln_H_input(1);
  H(1,1) = (1+ln_H_input(1)*ln_H_input(1)) / exp(ln_H_input(0));
  vector<Type> Range_raw(n_factors);
  Type H_trace = H(0,0)+H(1,1);
  Type H_det = H(0,0)*H(1,1)-H(0,1)*H(1,0);
  vector<Type> log_kappa(n_factors);
  vector<Type> kappa2(n_factors);
  vector<Type> kappa4(n_factors);
  for( int k=0;k<n_factors;k++){
    if(VaryingKappa==0) log_kappa(k) = log_kappa_input(0);
    if(VaryingKappa==1) log_kappa(k) = log_kappa_input(k);
    kappa2(k) = exp(2.0*log_kappa(k));
    kappa4(k) = kappa2(k)*kappa2(k);
    Range_raw(k) = sqrt(8) / exp( log_kappa(k) );
    REPORT( Range_raw(k)*(H_trace/2+sqrt(H_trace*H_trace/4-H_det)) );
    REPORT( Range_raw(k)*(H_trace/2-sqrt(H_trace*H_trace/4-H_det)) );
  }

  // Calculate adjugate of H
  matrix<Type> adj_H(2,2);
  adj_H(0,0) = H(1,1);
  adj_H(0,1) = -1 * H(0,1);
  adj_H(1,0) = -1 * H(1,0);
  adj_H(1,1) = H(0,0);

  Eigen::SparseMatrix<Type> G1_aniso(n_i,n_i); 
  Eigen::SparseMatrix<Type> G2_aniso(n_i,n_i); 
  if(Aniso==1){
    // Calculate G1 - pt. 1
    array<Type> Gtmp(n_tri,3,3);
    for(i=0; i<n_tri; i++){    
      // 1st line: E0(i,) %*% adjH %*% t(E0(i,)), etc.    
      Gtmp(i,0,0) = (E0(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E0(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));  
      Gtmp(i,0,1) = (E1(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E1(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));  
      Gtmp(i,0,2) = (E2(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E2(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
      Gtmp(i,1,1) = (E1(i,0)*(E1(i,0)*adj_H(0,0)+E1(i,1)*adj_H(1,0)) + E1(i,1)*(E1(i,0)*adj_H(0,1)+E1(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
      Gtmp(i,1,2) = (E2(i,0)*(E1(i,0)*adj_H(0,0)+E1(i,1)*adj_H(1,0)) + E2(i,1)*(E1(i,0)*adj_H(0,1)+E1(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
      Gtmp(i,2,2) = (E2(i,0)*(E2(i,0)*adj_H(0,0)+E2(i,1)*adj_H(1,0)) + E2(i,1)*(E2(i,0)*adj_H(0,1)+E2(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    }
    // Calculate G1 - pt. 2
    for(i=0; i<n_tri; i++){
      int i0 = i;
      int i1 = i + n_tri; 
      int i2 = i + 2*n_tri; 
      G1_aniso.coeffRef(TV(i1),TV(i0)) = G1_aniso.coeffRef(TV(i1),TV(i0)) + (Gtmp(i,0,1));  
      G1_aniso.coeffRef(TV(i0),TV(i1)) = G1_aniso.coeffRef(TV(i0),TV(i1)) + (Gtmp(i,0,1));  
      G1_aniso.coeffRef(TV(i2),TV(i1)) = G1_aniso.coeffRef(TV(i2),TV(i1)) + (Gtmp(i,1,2));  
      G1_aniso.coeffRef(TV(i1),TV(i2)) = G1_aniso.coeffRef(TV(i1),TV(i2)) + (Gtmp(i,1,2));  
      G1_aniso.coeffRef(TV(i2),TV(i0)) = G1_aniso.coeffRef(TV(i2),TV(i0)) + (Gtmp(i,0,2));  
      G1_aniso.coeffRef(TV(i0),TV(i2)) = G1_aniso.coeffRef(TV(i0),TV(i2)) + (Gtmp(i,0,2));  
      G1_aniso.coeffRef(TV(i0),TV(i0)) = G1_aniso.coeffRef(TV(i0),TV(i0)) + (Gtmp(i,0,0));  
      G1_aniso.coeffRef(TV(i1),TV(i1)) = G1_aniso.coeffRef(TV(i1),TV(i1)) + (Gtmp(i,1,1));  
      G1_aniso.coeffRef(TV(i2),TV(i2)) = G1_aniso.coeffRef(TV(i2),TV(i2)) + (Gtmp(i,2,2));  
    }
    
    G2_aniso = G1_aniso * G0_inv * G1_aniso; 
  }
  
  vector<Type> h(n_species);
  vector<Type> log_tau(n_factors);
  matrix<Type> Omega(n_i,n_factors);
  using namespace density;
  Eigen::SparseMatrix<Type> Q;
  for (int k=0;k<n_factors;k++){
    if(Aniso==0) Q = kappa4(k)*G0 + Type(2.0)*kappa2(k)*G1 + G2;
    if(Aniso==1) Q = kappa4(k)*G0 + Type(2.0)*kappa2(k)*G1_aniso + G2_aniso;
    g += GMRF(Q)(Omega_input.col(k));
    if(Aniso==0) log_tau(k) = log(1 / (exp(log_kappa(k)) * sqrt(4*3.141592)) );  // Ensures that MargSD = 1
    if(Aniso==1) log_tau(k) = log(1 / (exp(log_kappa(k)) * sqrt(4*3.141592*sqrt(H(0,0)*H(1,1)-H(0,1)*H(1,0)))) );  // Ensures that MargSD = 1
    for(int k2=0;k2<n_i;k2++) Omega(k2,k) = Omega_input(k2,k) / exp(log_tau(k));
  }
  
  matrix<Type> Psi(n_factors,n_species);
  int Count = 0;
  for(int k=0;k<n_factors;k++){
    for(int i=0;i<n_species;i++){
      if(i>=k){
        Psi(k,i) = Psi_val(Count);
        Count++;
      }else{
        Psi(k,i) = 0.0;
      }
    }
  }
  
  // Hyper-distribution for lognormal variance inflation
  if( Options_Vec(0)==1 | Options_Vec(0)==2 ){
    for(int i=0;i<n_species;i++){
    for(int j=0;j<n_stations;j++){
      g -= dnorm( Lognorm_input(j,i), Type(0.0), exp(ln_VarInfl_Lognorm(i)), true );
    }}
    // Jacobian for log-transformation (only necessary if using REML to integrate across ln_VarInfl_Lognorm)
    if( Options_Vec(0)==2 ){ 
      for(int i=0;i<n_species;i++){
        g -= ln_VarInfl_Lognorm(i);
      }
    }
  }
  
  // Likelihood contribution from observations
  Type Tmp;
  Type pred_tmp;
  matrix<Type> ln_mean_y( n_stations, n_species );
  matrix<Type> lin_pred( n_stations, n_species );
  lin_pred = X * beta; 
  for(int i=0;i<n_species;i++){
    //lin_pred = X * beta.col(i).matrix();
    h(i) = 0;
    for(int j=0;j<n_stations;j++){
      if(NAind(j,i) < 0.5){
        // Predictor
        ln_mean_y(j,i) = lin_pred(j,i);
        for(int k=0;k<n_factors;k++){
          ln_mean_y(j,i) = ln_mean_y(j,i) + Omega_input(meshidxloc(j),k) / exp(log_tau(k)) * Psi(k,i);
        }
        // Additional variance
        pred_tmp = ln_mean_y(j,i);
        if( Options_Vec(0)==1 ) pred_tmp += Lognorm_input(j,i);
        // Likelihood
        Type var_y = exp(ln_VarInfl_NB0(i)) + exp(pred_tmp)*(1.0+exp(ln_VarInfl_NB1(i))) + exp(2.0*pred_tmp)*exp(ln_VarInfl_NB2(i));
        if(ErrorDist==0){
          g -= dpois( Y(j,i), exp(pred_tmp), true ) * (1-isPred(j)); 
          h(i) -= dpois( Y(j,i), exp(pred_tmp), true ) * isPred(j); 
        }
        if(ErrorDist==1){
          g -= dnbinom2( Y(j,i), exp(pred_tmp), var_y, true ) * (1-isPred(j)); 
          h(i) -= dnbinom2( Y(j,i), exp(pred_tmp), var_y, true ) * isPred(j); 
        }
        if(ErrorDist==2){ 
          Tmp = -log( exp(ln_VarInfl_ZI(i)) + 1) + dnbinom2( Y(j,i), exp(pred_tmp), var_y, true );
          if(Y(j,i)==0){ Tmp = log( 1/(1+exp(-ln_VarInfl_ZI(i))) + exp(Tmp) ); }
          g -= Tmp * (1-isPred(j)); 
          h(i) -= Tmp * isPred(j); 
        }
        if(Y_Report(j,i) > 0.5){
          REPORT( lin_pred(j,i) );
          REPORT( exp(ln_mean_y(j,i)) );
        }
      }
    }
  }
  
  // Add penalties to increase stiffness -- Mean = 0
  vector<Type> Mean(n_factors); 
  if( Pen_Vec(0)>0.001 | Pen_Vec(0)<(-0.001) ){
    for(int k=0;k<n_factors;k++){
      Mean(k) = 0;
      for(int l=0;l<Omega.col(k).size();l++){
        Mean(k) += Omega(l,k);
      }
      Mean(k) = Mean(k) / Omega.col(k).size();
      if( Pen_Vec(0)>0.001 ) g += Pen_Vec(0) * pow(Mean(k), 2);
    }
    REPORT( Mean );
  }
  // Add penalties to increase stiffness -- Cov = 0
  matrix<Type> Cov(n_factors,n_factors);
  if( (Pen_Vec(1)>0.001 | Pen_Vec(1)<(-0.001)) & n_factors>=2 ){
    for(int k1=0;k1<n_factors;k1++){
    for(int k2=k1+1;k2<n_factors;k2++){
      Cov(k1,k2) = 0; 
      for(int l=0;l<Omega.col(k1).size();l++){
        Cov(k1,k2) += (Omega(l,k1)-Mean(k1)) * (Omega(l,k2)-Mean(k2));
      }
      Cov(k1,k2) = Cov(k1,k2) / Omega.col(k1).size();
      if( Pen_Vec(1)>0.001 ) g += Pen_Vec(1) * pow(Cov(k1,k2), 2);
    }}
    REPORT( Cov );
  }
  
  // Reporting
  for(int k=0;k<n_factors;k++){
    for(int i=0;i<n_species;i++){
      Type Marg_SD = 1/sqrt(4*3.141592*exp(2*log_tau(k))*exp(2*log_kappa(k)));
      Type Eff_SD = 1/sqrt(4*3.141592*exp(2*log_tau(k))*exp(2*log_kappa(k))) * Psi(k,i);
      ADREPORT( Marg_SD ); 
      ADREPORT( Eff_SD ); 
    }
  }
  for(int i1=0;i1<n_species;i1++){
  for(int i2=0;i2<n_species;i2++){
    Type Num = 0;
    Type Denom1 = 0;
    Type Denom2 = 0;
    for(int k=0;k<n_factors;k++){
      Num += Psi(k,i1) * Psi(k,i2);
      Denom1 += Psi(k,i1) * Psi(k,i1);
      Denom2 += Psi(k,i2) * Psi(k,i2);
    }
    Type Corr = Num / pow(Denom1,0.5) / pow(Denom2,0.5);
    ADREPORT( Corr ); 
  }}
  ADREPORT( log_kappa );
  ADREPORT( log_tau );
  REPORT( Range_raw );
  ADREPORT( Omega );
  REPORT( Omega );
  for(int k=0;k<n_factors;k++){
    ADREPORT( Psi.row(k) );
  }
  ADREPORT( h );
  REPORT( h );
  REPORT( ln_mean_y );
  REPORT( Psi );
  REPORT( H );
  REPORT( adj_H );
  if( Options_Vec(0)==1 | Options_Vec(0)==2) REPORT(Lognorm_input);

  return g;
}
