#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace Rcpp;
using namespace arma;

static double lfactorial_cpp(int i);
//static int max_cpp(int a, int b);
//static int min_cpp(int a, int b);
static NumericVector sort_cpp(NumericVector v);
static IntegerVector sort_int(IntegerVector v);
static double mean_cpp(NumericVector x);

static double ODE_fun(int K, double gamma, double p, double alpha, int C);
static IntegerVector delta_randomize( int T, int M );
static IntegerVector Vector_equal(IntegerVector x);
static IntegerVector Locator(IntegerVector x);
static IntegerVector mcmc_moves(IntegerVector delta_vec, NumericVector w);
static CharacterVector csample_char( CharacterVector x, int size, bool replace, NumericVector prob );



// [[Rcpp::export]]
Rcpp::List growth_cp(IntegerVector C, int M, double is_p, double is_alpha, int POP, int T_fin, NumericVector w, bool store = 0, bool ppm = 0, bool cred = 0, bool predict = 0) {
  // Read data information
  int T = C.length();
  //int T = DeltaC.length();
  //IntegerVector C (T);
  bool p_unknown = (is_p == -1.0);
  bool alpha_unknown = (is_alpha == -1.0);
  
  // Set algorithm settings
  int iter = 100000; 
  int burn = iter*0.5;
  int m;
  //int K_start = (max(C) + 1)*2;
  //IntegerVector a_K (M);
  
  //for(m = 0; m < M; m++){
  //  a_K(m) = a_K(m) + (max(C)/M);
  //  K_start(m) = 2*a_K(m);
  //}
  
  double phi_start = 10, lambda_start = 0.5, p_start = 0.5, alpha_start = 1.0;
  IntegerVector delta_start = delta_randomize(T, M);
  IntegerVector m_start (M);
  IntegerVector K_start (M);
  int t;
  
  //C(0) = 100; // (only for simulation study)Change to something else later
  //for(t = 1; t < T; t++){
  //  C(t) = C(t - 1) + DeltaC(t);
  //}
  //Rcout<<C;
  if(M == 1){
    K_start(0) = (max(C) + 1)*2;
  }
  else{
    m_start = Locator(delta_start);
    m_start.push_back(T - 1);
    for(m = 0; m < M; m++){
      if(m == M - 1){
        K_start(m) = (max(C) + 1)*2;
      }
      else{
        K_start(m) = (C(m_start(m + 1)) + 1)*2;
      }
    }
  }
  
  
  double tau_logphi = 1.0, tau_logK = 1.0, tau_loglambda = 0.1, tau_logp = 0.1, tau_logalpha = 0.1;
  
  // Set hyperparameters
  double a_phi = 0.001, b_phi = 0.001; 
  int b_K = POP;
  
  //int a_K = max(C) + 1;
  double a_lambda, b_lambda;
  if (is_p == 1.0 && is_alpha == 1.0)
  {
    a_lambda = 1.0;
    b_lambda = 1.0;
  }
  else
  {
    a_lambda = 0.001;
    b_lambda = 0.001;
  }
  double a_p = 1, b_p = 1;
  double a_alpha = 0.001, b_alpha = 0.001;
  
  // Non informative omega
  NumericVector omega (T);
  omega(0) = 1;
  int ii;
  for(ii = 1; ii < T; ii++){  // Time points 2 - 6 and T - 5 to T cannot be changepoints 
    if(ii < 7){
      omega(ii) = 0.0;
    }
    else if(ii < T - 6 && ii >= 7){
      omega(ii) = 0.001;
    }
    else{
      omega(ii) = 0.0;
    }
  }
  
  // Set temporary variables
  int it, i, count = 0, count_2 = 0;
  NumericVector phi_temp (M);
  NumericVector lambda_temp (M);
  NumericVector p_temp (M);
  NumericVector alpha_temp (M);
  IntegerVector K_temp (M);
  IntegerVector diffK_temp (M);
  IntegerVector delta_temp (T);
  
  NumericVector phi_map (M);
  NumericVector lambda_map (M);
  NumericVector p_map (M);
  NumericVector alpha_map (M);
  if(is_p == 1.0){
    p_map.fill(1.0);
  }
  if(is_alpha == 1.0){
    alpha_map.fill(1.0);
  }
  IntegerVector K_map (M);
  IntegerVector delta_map (T);
  
  double hastings, logposterior = 0, logposterior_map, logomega = 0, lognum, logden, logomega_temp;
  NumericVector logposterior_m (M); // logposterior of sub-segment m
  //NumericVector hastings (M);
  NumericVector ODE (T - 1); // g function
  NumericVector ODE_temp (T - 1);
  IntegerVector accept_phi (M);
  IntegerVector accept_K (M);
  IntegerVector accept_alpha (M);
  IntegerVector accept_p (M);
  IntegerVector accept_lambda (M);
  double accept_delta = 0;
  NumericVector logposterior_store(iter);
  NumericMatrix phi_store(M, iter - burn);
  NumericMatrix K_store(M, iter - burn);
  NumericMatrix lambda_store(M, iter - burn);
  NumericMatrix p_store(M, iter - burn);
  NumericMatrix alpha_store(M, iter - burn);
  IntegerMatrix delta_store(iter - burn, T);
  IntegerVector z_store(T);
  NumericVector delta_prop(T);
  IntegerMatrix C_predict(200, T_fin);
  IntegerMatrix N_predict(200, T_fin);
  Mat<int> CI_mat(iter - burn, T, fill::zeros);
  
  // Initialization
  
  IntegerVector K = Vector_equal(K_start);
  IntegerVector diffK (M);
  //IntegerVector K = IntegerVector::create(10000, 8000, 12000);
  NumericVector p (M, p_start);
  NumericVector alpha (M, alpha_start);
  if(is_p == 1.0){
    p.fill(1.0);
  }
  if(is_alpha == 1.0){
    alpha.fill(1.0);
  }
  NumericVector phi (M, phi_start);
  //NumericVector phi = NumericVector::create(100.0, 100.0);
  NumericVector lambda (M, lambda_start);
  //NumericVector lambda = NumericVector::create(0.15, 0.2);
  IntegerVector delta = Vector_equal(delta_start);
  //IntegerVector DeltaC (T - 1);
  IntegerVector m_loc (T);
  IntegerVector m_loc_temp (T);
  arma::mat ppm_store(T, T, arma::fill::zeros);

  for(m = 0; m < M; m++){
    phi_map[m] = phi[m];
    K_map[m] = K[m];
    lambda_map[m] = lambda[m];
    if(p_unknown){
      p_map[m] = p[m];
    }
    if(alpha_unknown){
      alpha_map[m] = alpha[m];
    }
  }
  
  delta_map = Vector_equal(delta);
  
  
  m_loc = Locator(delta);  // To extract the location of changepoints from delta vector
  m_loc.push_back(T - 1);  // To bound the last sub-segment by adding 1 after time T (required for inner loop)
  
  //Rcout<<m_loc;
  
  logposterior = 0;
  
   //DeltaC(0) = C(0);
   IntegerVector DeltaC(T - 1);
   for(t = 0; t < T - 1; t++){
   DeltaC(t) = C(t + 1) - C(t);
   }
   
   //Rcout<< DeltaC.length();
   //Rcout<< C.length();
     
  
  for( m = 0; m < M; m++){ // loop for outer sum
    for(t = m_loc(m); t < m_loc(m + 1); t++){ // loop for sub-segment specific inner sum 
      ODE(t) = ODE_fun(K(m), lambda(m), p(m), alpha(m), C(t));
      if(ODE(t) <= 0.0){
        //ODE(t) = std::numeric_limits<double>::epsilon(); // machine epsilon
        ODE(t) = 2.220446e-16; 
      }
      logposterior = logposterior + lgamma(DeltaC(t) + phi(m)) - lgamma(phi(m)) - lfactorial_cpp(DeltaC(t)) + phi(m)*(log(phi(m)) - log(ODE(t) + phi(m))) + DeltaC(t)*(log(ODE(t)) - log(ODE(t) + phi(m)));
    }
    logposterior = logposterior + (a_phi - 1)*log(phi(m)) - b_phi*phi(m);
    
    if (alpha(m) == 1.0 && p(m) == 1.0){
      logposterior = logposterior + (a_lambda - 1)*log(lambda(m)) + (b_lambda - 1)*log(1 - lambda(m));
    }
    else{
      logposterior = logposterior + (a_lambda - 1)*log(lambda(m)) - b_lambda*lambda(m);
    }
    
    if(alpha_unknown){
      logposterior = logposterior + (a_alpha - 1)*log(alpha(m)) - b_alpha*alpha(m);
    }
    
    if(p_unknown){
      logposterior = logposterior + (a_p - 1)*log(p(m)) + (b_p - 1)*log(1 - p(m));
    }
    //logposterior = logposterior + logposterior_m(m);
  }
  
  if(M > 1){
    for(t = 7; t < T - 6; t++){ // bernoulli product prior for delta
      logomega = logomega + delta(t)*log(omega(t)) + (1 - delta(t))*log(1 - omega(t));
    }
    logposterior = logposterior + logomega;
  }
  //Rcout<<logposterior;
  
  // MCMC
  
  for (it = 0; it < iter; it++){
    
    if(M > 1){
      // Update delta (cp vector)
      
      delta_temp = mcmc_moves(delta, w); // update delta
      
      m_loc_temp = Locator(delta_temp);
      m_loc_temp.push_back(T - 1);
      
      hastings = 0;
      logden = 0;
      lognum = 0;
      
      for( m = 0; m < M; m++){ // Outer sum
        for(t = m_loc(m); t < m_loc(m + 1); t++){ // Inner sum for old delta
          logden = logden + lgamma(DeltaC(t) + phi(m)) - lgamma(phi(m)) - lfactorial_cpp(DeltaC(t)) + phi(m)*(log(phi(m)) - log(ODE(t) + phi(m))) + DeltaC(t)*(log(ODE(t)) - log(ODE(t) + phi(m)));
        }
        
        for(t = m_loc_temp(m); t < m_loc_temp(m + 1); t++){ // Inner sum for updated delta
          ODE_temp(t) = ODE_fun(K(m), lambda(m), p(m), alpha(m), C(t));
          if(ODE_temp(t) <= 0.0){
            //ODE(t) = std::numeric_limits<double>::epsilon(); // machine epsilon
            ODE_temp(t) = 2.220446e-16; 
          }
          lognum = lognum + lgamma(DeltaC(t) + phi(m)) - lgamma(phi(m)) - lfactorial_cpp(DeltaC(t)) + phi(m)*(log(phi(m)) - log(ODE_temp(t) + phi(m))) + DeltaC(t)*(log(ODE_temp(t)) - log(ODE_temp(t) + phi(m)));
        }
      }
      
      //logomega_temp = 0;
      //logomega = 0;
      
      /*for(t = 7; t < T - 6; t++){  // Not needed if all omegas are equal
       logomega_temp = logomega_temp + delta_temp(t)*log(omega(t)) + (1 - delta_temp(t))*log(1 - omega(t));
       logomega = logomega + delta(t)*log(omega(t)) + (1 - delta(t))*log(1 - omega(t));
      }*/
      
      hastings = hastings + lognum - logden;
      //hastings = hastings + logomega_temp - logomega; // Not needed if all omegas are equal
      
      if(hastings >= log(double(rand()%10001)/10000)){
        delta = Vector_equal(delta_temp);
        m_loc = Locator(delta);
        m_loc.push_back(T - 1);
        for( t = 0; t < T - 1; t++){
          ODE(t) = ODE_temp(t);
        }
        logposterior = logposterior + hastings;
        if (it > burn) {
          accept_delta++;
        }
      }
      
      // For change point proportions plot
      
      if(store){
        if(it > burn){
          for(t = 0; t < T; t++){
            if(delta(t) == 1.0){
              delta_prop(t) = delta_prop(t) + 1.0;
            }
          }
        }
      }
    }
    
    
    hastings = 0;
    
    for( m = 0; m < M; m++){
      // Update phi_m
      
      phi_temp(m) = exp(r_truncnorm(log(phi(m)), tau_logphi, log(1), log(100)));
      hastings = 0;
      for(t = m_loc(m); t < m_loc(m + 1); t++){
        
        hastings = hastings + phi_temp(m)*log(phi_temp(m)) - lgamma(phi_temp(m)) + lgamma(phi_temp(m) + DeltaC(t)) - (phi_temp(m) + DeltaC(t))*log(phi_temp(m) + ODE(t));
        hastings = hastings - (phi(m)*log(phi(m)) - lgamma(phi(m)) + lgamma(phi(m) + DeltaC(t)) - (phi(m) + DeltaC(t))*log(phi(m) + ODE(t)));
        
      }
      hastings = hastings + (a_phi - 1)*log(phi_temp(m)) - b_phi*phi_temp(m);
      hastings = hastings - ((a_phi - 1)*log(phi(m)) - b_phi*phi(m));
      
      if(hastings >= log(double(rand()%10001)/10000)){
        phi(m) = phi_temp(m);
        logposterior = logposterior + hastings;
        if (it > burn) {
          accept_phi(m) = accept_phi(m) + 1;
        }
      }
      
      /*
       // Update K_m
       
       K_temp(m) = exp(r_truncnorm(log(K(m)), tau_logK, log(a_K(m)), log(b_K)));
       hastings = 0;
       
       for(t = m_loc(m); t < m_loc(m + 1); t++){
       
       ODE_temp(t) = ODE_fun(K_temp(m), lambda(m), p(m), alpha(m), C(t));
       if(ODE_temp(t) <= 0.0){
       //ODE(t) = std::numeric_limits<double>::epsilon(); // machine epsilon
       ODE_temp(t) = 2.220446e-16; 
       }
       //hastings = hastings + DeltaC(t)*log(1 - 1.0*pow(1.0*C(t)/K_temp(m), alpha(m))) - (phi(m) + DeltaC(t))*log(ODE_temp(t) + phi(m));
       //hastings = hastings - (DeltaC(t)*log(1 - 1.0*pow(1.0*C(t)/K(m), alpha(m))) - (phi(m) + DeltaC(t))*log(ODE(t) + phi(m)));
       
       hastings = hastings + DeltaC(t)*log(ODE_temp(t)) - (phi(m) + DeltaC(t))*log(ODE_temp(t) + phi(m));
       hastings = hastings - (DeltaC(t)*log(ODE(t)) - (phi(m) + DeltaC(t))*log(ODE(t) + phi(m)));
       
       }
       if(hastings >= log(double(rand()%10001)/10000)){
       K(m) = K_temp(m);
       logposterior = logposterior + hastings;
       for(t = m_loc(m); t < m_loc(m + 1); t++){
       ODE(t) = ODE_temp(t);
       }
       if (it > burn) {
       accept_K(m) = accept_K(m) + 1;
       }
       }
       
       */
      // Update K_m
      if(m < M - 1){
        K_temp(m) = exp(r_truncnorm(log(K(m)), tau_logK, log(C(m_loc(m + 1))), log(b_K)));
      }
      else if(m == M - 1){
        K_temp(m) = exp(r_truncnorm(log(K(m)), tau_logK, log(1.0*max(C)), log(b_K)));
      }
      hastings = 0;
      
      for(t = m_loc(m); t < m_loc(m + 1); t++){
        
        ODE_temp(t) = ODE_fun(K_temp(m), lambda(m), p(m), alpha(m), C(t));
        if(ODE_temp(t) <= 0.0){
          //ODE(t) = std::numeric_limits<double>::epsilon(); // machine epsilon
          ODE_temp(t) = 2.220446e-16; 
        }
        //hastings = hastings + DeltaC(t)*log(1 - 1.0*pow(1.0*C(t)/K_temp(m), alpha(m))) - (phi(m) + DeltaC(t))*log(ODE_temp(t) + phi(m));
        //hastings = hastings - (DeltaC(t)*log(1 - 1.0*pow(1.0*C(t)/K(m), alpha(m))) - (phi(m) + DeltaC(t))*log(ODE(t) + phi(m)));
        
        hastings = hastings + DeltaC(t)*log(ODE_temp(t)) - (phi(m) + DeltaC(t))*log(ODE_temp(t) + phi(m));
        hastings = hastings - (DeltaC(t)*log(ODE(t)) - (phi(m) + DeltaC(t))*log(ODE(t) + phi(m)));
        
      }
      if(hastings >= log(double(rand()%10001)/10000)){
        K(m) = K_temp(m);
        logposterior = logposterior + hastings;
        for(t = m_loc(m); t < m_loc(m + 1); t++){
          ODE(t) = ODE_temp(t);
        }
        if (it > burn) {
          accept_K(m) = accept_K(m) + 1;
        }
      }
      
      // Update lambda_m
      if (p(m) == 1.0 && alpha(m) == 1.0) 
      {
        lambda_temp(m) = exp(r_truncnorm(log(lambda(m)), tau_loglambda, log(10e-9), log(1)));
      }
      else{
        lambda_temp(m) = exp(rnorm(1, log(lambda(m)), tau_loglambda)(0));
      }
      hastings = 0;
      
      for(t = m_loc(m); t < m_loc(m + 1); t++){
        
        ODE_temp(t) = ODE_fun(K(m), lambda_temp(m), p(m), alpha(m), C(t));
        if(ODE_temp(t) <= 0.0){
          //ODE(t) = std::numeric_limits<double>::epsilon(); // machine epsilon
          ODE_temp(t) = 2.220446e-16; 
        }
        
        hastings = hastings + DeltaC(t)*log(lambda_temp(m)) - (phi(m) + DeltaC(t))*log(ODE_temp(t) + phi(m));
        hastings = hastings - (DeltaC(t)*log(lambda(m)) - (phi(m) + DeltaC(t))*log(ODE(t) + phi(m)));
        
      }
      if(p(m) == 1.0 && alpha(m) == 1.0){
        hastings = hastings + (a_lambda - 1)*log(lambda_temp(m)) + (b_lambda - 1)*log(1 - lambda_temp(m));
        hastings = hastings - ((a_lambda - 1)*log(lambda(m)) + (b_lambda - 1)*log(1 - lambda(m)));
      }
      else{
        hastings = hastings + (a_lambda - 1)*log(lambda_temp(m)) - b_lambda*lambda_temp(m);
        hastings = hastings - ((a_lambda - 1)*log(lambda(m)) - b_lambda*lambda(m));
      }
      
      if(hastings >= log(double(rand()%10001)/10000)){
        lambda(m) = lambda_temp(m);
        logposterior = logposterior + hastings;
        for(t = m_loc(m); t < m_loc(m + 1); t++){
          ODE(t) = ODE_temp(t);
        }
        if (it > burn) {
          accept_lambda(m) = accept_lambda(m) + 1;
        }
      }
      
      
      // Update p_m
      if(p_unknown){
        p_temp(m) = exp(r_truncnorm(log(p(m)), tau_logp, log(10e-9), log(1)));
        hastings = 0;
        
        for(t = m_loc(m); t < m_loc(m + 1); t++){
          
          ODE_temp(t) = ODE_fun(K(m), lambda(m), p_temp(m), alpha(m), C(t));
          if(ODE_temp(t) <= 0.0){
            //ODE(t) = std::numeric_limits<double>::epsilon(); // machine epsilon
            ODE_temp(t) = 2.220446e-16; 
          }
          //hastings = hastings + p_temp(m)*DeltaC(t)*log(C(t)) - (phi(m) + DeltaC(t))*log(ODE_temp(t) + phi(m));
          //hastings = hastings - (p(m)*DeltaC(t)*log(C(t)) - (phi(m) + DeltaC(t))*log(ODE(t) + phi(m)));
          
          hastings = hastings + DeltaC(t)*log(ODE_temp(t)) - (phi(m) + DeltaC(t))*log(ODE_temp(t) + phi(m));
          hastings = hastings - (DeltaC(t)*log(ODE(t)) - (phi(m) + DeltaC(t))*log(ODE(t) + phi(m)));
          
        }
        hastings = hastings + (a_p - 1)*log(p_temp(m)) + (b_p - 1)*log(1 - p_temp(m));
        hastings = hastings - ((a_p - 1)*log(p(m)) + (b_p - 1)*log(1 - p(m)));
        
        if(hastings >= log(double(rand()%10001)/10000)){
          p(m) = p_temp(m);
          logposterior = logposterior + hastings;
          for(t = m_loc(m); t < m_loc(m + 1); t++){
            ODE(t) = ODE_temp(t);
          }
          if (it > burn) {
            accept_p(m) = accept_p(m) + 1;
          }
        }
      }
      
      
      // Update alpha_m
      if(alpha_unknown){
        alpha_temp(m) = exp(rnorm(1, log(alpha(m)), tau_logalpha)(0));
        hastings = 0;
        
        for(t = m_loc(m); t < m_loc(m + 1); t++){
          
          ODE_temp(t) = ODE_fun(K(m), lambda(m), p(m), alpha_temp(m), C(t));
          if(ODE_temp(t) <= 0.0){
            //ODE(t) = std::numeric_limits<double>::epsilon(); // machine epsilon
            ODE_temp(t) = 2.220446e-16; 
          }
          //hastings = hastings + DeltaC(t)*log(1 - 1.0*pow(1.0*C(t)/K(m), alpha_temp(m))) - (phi(m) + DeltaC(t))*log(ODE_temp(t) + phi(m));
          //hastings = hastings - (DeltaC(t)*log(1 - 1.0*pow(1.0*C(t)/K(m), alpha(m))) - (phi(m) + DeltaC(t))*log(ODE(t) + phi(m)));
          
          hastings = hastings + DeltaC(t)*log(ODE_temp(t)) - (phi(m) + DeltaC(t))*log(ODE_temp(t) + phi(m));
          hastings = hastings - (DeltaC(t)*log(ODE(t)) - (phi(m) + DeltaC(t))*log(ODE(t) + phi(m)));
          
        }
        hastings = hastings + (a_alpha - 1)*log(alpha_temp(m)) - b_alpha*alpha_temp(m);
        hastings = hastings - ((a_alpha - 1)*log(alpha(m)) - b_alpha*alpha(m));
        
        if(hastings >= log(double(rand()%10001)/10000)){
          alpha(m) = alpha_temp(m);
          logposterior = logposterior + hastings;
          for(t = m_loc(m); t < m_loc(m + 1); t++){
            ODE(t) = ODE_temp(t);
          }
          if (it > burn) {
            accept_alpha(m) = accept_alpha(m) + 1;
          }
        }
      }
      // Storing logposterior
      
      if(store){
        logposterior_store(it) = logposterior;
      }
      
      // Monitor the process
      if(it*100/iter == count)
      {
        Rcout<<count<< "% has been done\n";
        count = count + 10;
      }
    }
    
    // Obtaining MAP estimates
    if( it >= burn){
      if( it == burn){
        logposterior_map = logposterior;
        delta_map = Vector_equal(delta);
        for( m = 0; m < M; m++){
          phi_map(m) = phi(m);
          K_map(m) = K(m);
          lambda_map(m) = lambda(m);
          if(p_unknown){
            p_map(m) = p(m);
          }
          if(alpha_unknown){
            alpha_map(m) = alpha(m);
          }
        }
      }
      else{
        if(logposterior > logposterior_map){
          logposterior_map = logposterior;
          delta_map = Vector_equal(delta);
          for( m = 0; m < M; m++){
            phi_map(m) = phi(m);
            K_map(m) = K(m);
            lambda_map(m) = lambda(m);
            if(p_unknown){
              p_map(m) = p(m);
            }
            if(alpha_unknown){
              alpha_map(m) = alpha(m);
            }
          }
        }
      }
      if(store){
        //delta_store.row(it - burn) = delta;
        for( m = 0; m < M; m++){
          phi_store(m, it - burn) = phi(m);
          K_store(m, it - burn) = K(m);
          lambda_store(m, it - burn) = lambda(m);
          if(p_unknown){
            p_store(m, it - burn) = p(m);
          }
          if(alpha_unknown){
            alpha_store(m, it - burn) = alpha(m);
          }
        }
        for(int t = 0; t < T; t++){
          if(delta(t) == 1){
            delta_prop(t) = delta_prop(t) + 1.0;
            CI_mat(it - burn, t) = 1;
          }
        }
      }
      // Prediction
      if(predict){
        if (it % 250 == 0)
        {
          for (t = 0; t < T_fin; t++) 
          {
            if (t == 0)
            {
              N_predict(count_2, t) = rnbinom_mu(1, phi(M - 1), ODE_fun(K(M - 1), lambda(M - 1), p(M - 1), alpha(M - 1), C(T - 1)))(0);
              C_predict(count_2, t) = C(T - 1) + N_predict(count_2, t);
            }
            else
            {
              N_predict(count_2, t) = rnbinom_mu(1, phi(M - 1), ODE_fun(K(M - 1), lambda(M - 1), p(M - 1), alpha(M - 1), C_predict(count_2, t - 1)))(0);
              C_predict(count_2, t) = C_predict(count_2, t - 1) + N_predict(count_2, t);
              
            }
          }
          count_2 = count_2 + 1;
        }
      }
    }
    // Obtaining delta ppm
    
    if(ppm){
      if(it >= burn){
        z_store(0) = delta(0);
        for(t = 0; t < T - 1; t++){
          z_store(t + 1) = z_store(t) + delta(t + 1);
        }
        for(t = 0; t < T - 1; t++){
          for(int tt = 1; tt < T; tt++ ){
            if(z_store[t] == z_store[tt]){
              ppm_store(t, tt) = ppm_store(t, tt) + 1;
              ppm_store(tt, t) = ppm_store(tt, t) + 1;
            }
          }
        }
      }
    }
    
  }
  
  
  // Credible Intervals
  NumericVector phi_upp(M);
  NumericVector phi_lwr(M);
  NumericVector K_upp(M);
  NumericVector K_lwr(M);
  NumericVector lambda_upp(M);
  NumericVector lambda_lwr(M);
  NumericVector p_upp(M);
  NumericVector p_lwr(M);
  NumericVector alpha_upp(M);
  NumericVector alpha_lwr(M);
  NumericVector vtemp;
  if(cred){
    for(m = 0; m < M; m++){
      vtemp = sort_cpp(phi_store.row(m));
      phi_lwr(m) = vtemp(round(0.025*(iter - burn)) - 1);
      phi_upp(m) = vtemp(round(0.975*(iter - burn)) - 1);
      vtemp = sort_cpp(K_store.row(m));
      K_lwr(m) = vtemp(round(0.025*(iter - burn)) - 1);
      K_upp(m) = vtemp(round(0.975*(iter - burn)) - 1);
      vtemp = sort_cpp(lambda_store.row(m));
      lambda_lwr(m) = vtemp(round(0.025*(iter - burn)) - 1);
      lambda_upp(m) = vtemp(round(0.975*(iter - burn)) - 1);
      vtemp = sort_cpp(p_store.row(m));
      p_lwr(m) = vtemp(round(0.025*(iter - burn)) - 1);
      p_upp(m) = vtemp(round(0.975*(iter - burn)) - 1);
      vtemp = sort_cpp(alpha_store.row(m));
      alpha_lwr(m) = vtemp(round(0.025*(iter - burn)) - 1);
      alpha_upp(m) = vtemp(round(0.975*(iter - burn)) - 1);
    }
  } 
  
  // Update C0 
  
  int c, C0 = 1;
  double lambda_0, logprob_max, logprob, ODE_0;
  for (c = 1; c < C(0); c++)
  {
    ODE_0 = ODE_fun(K_map[0], lambda_map[0], p_map[0], alpha_map[0], c);
    if (c == 1)
    {  // Would it be the parameters of first subsegment or last??
      logprob_max = lgamma(C(0) - C0 + phi_map[0]) - lfactorial_cpp(C(0) - C0) - phi_map[0]*log(ODE_0 + phi_map[0]) + (C(0) - C0)*(log(ODE_0) - log(ODE_0 + phi_map[0]));
    }
    else
    {
      logprob = lgamma(C(0) - C0 + phi_map[0]) - lfactorial_cpp(C(0) - C0) - phi_map[0]*log(ODE_0 + phi_map[0]) + (C(0) - C0)*(log(ODE_0) - log(ODE_0 + phi_map[0]));
      if (logprob > logprob_max)
      {
        logprob_max = logprob;
        C0 = c;
      }
    }
  }
  
  // Fit the MAP model
  // int peak; //= K_map*pow(1.0*p_map/(alpha_map + p_map), 1/alpha_map);
  NumericVector C_fit(T + T_fin, 1.0*POP);
  NumericVector N_fit(T + T_fin, 0.0);
  m = 0;
  for (t = 0; t < T + T_fin; t++)
  {
    if (t == 0)
    {
      N_fit(t) = ODE_fun(K_map[0], lambda_map[0], p_map[0], alpha_map[0], C0);
      C_fit(t) = C0 + N_fit(t);
    }
    else
    {
      if( t < T){
        if( delta_map(t) == 1){
          m++;
        }
      }
     
      N_fit(t) = ODE_fun(K_map[m], lambda_map[m], p_map[m], alpha_map[m], C_fit(t - 1));
      if (C_fit(t - 1) + N_fit(t) >= POP)
      {
        N_fit(t) = POP - C_fit(t - 1);
        break;
      }
      C_fit(t) = C_fit(t - 1) + N_fit(t);
    }
  }

  // Prediction

  NumericVector C_mean(T_fin);
  IntegerVector C_upp(T_fin);
  IntegerVector C_lwr(T_fin);
  NumericVector N_mean(T_fin);
  IntegerVector N_upp(T_fin);
  IntegerVector N_lwr(T_fin);
  

  if(predict){
    for (t = 0; t < T_fin; t++) 
    {
      vtemp = sort_int(C_predict.column(t));
      C_lwr(t) = vtemp(4);
      C_upp(t) = vtemp(194);
      C_mean(t) = mean_cpp(vtemp);
      vtemp = sort_int(N_predict.column(t));
      N_lwr(t) = vtemp(4);
      N_upp(t) = vtemp(194);
      N_mean(t) = mean_cpp(vtemp);
    }
  }
  
  // Result wrap-up
  /*for ( m = 0; m < M; m++){
   accept_phi(m) = accept_phi(m)/(iter - burn);
   accept_K(m) = accept_K(m)/(iter - burn);
   accept_lambda(m) = accept_lambda(m)/(iter - burn);
   accept_p(m) = accept_p(m)/(iter - burn);
   accept_alpha(m) = accept_alpha(m)/(iter - burn);
  }
   accept_delta = accept_delta/(iter - burn);*/
  
  Rcpp::List accept = Rcpp::List::create(
    Rcpp::Named("accept_delta") = accept_delta,
    Rcpp::Named("accept_phi") = accept_phi,
    Rcpp::Named("accept_K") = accept_K,
    Rcpp::Named("accept_lambda") = accept_lambda,
    Rcpp::Named("accept_p") = accept_p,
    Rcpp::Named("accept_alpha") = accept_alpha
  );
  
  Rcpp::List map;
  Rcpp::List store_list;
  if(cred){
    map = Rcpp::List::create(
      Rcpp::Named("logposterior_map") = logposterior_map, 
      Rcpp::Named("phi_map") = phi_map,
      Rcpp::Named("phi_upp") = phi_upp,
      Rcpp::Named("phi_lwr") = phi_lwr,
      Rcpp::Named("K_map") = K_map, 
      Rcpp::Named("K_upp") = K_upp, 
      Rcpp::Named("K_lwr") = K_lwr, 
      Rcpp::Named("lambda_map") = lambda_map, 
      Rcpp::Named("lambda_upp") = lambda_upp, 
      Rcpp::Named("lambda_lwr") = lambda_lwr, 
      Rcpp::Named("p_map") = p_map, 
      Rcpp::Named("p_upp") = p_upp, 
      Rcpp::Named("p_lwr") = p_lwr,
      Rcpp::Named("alpha_map") = alpha_map,
      Rcpp::Named("alpha_upp") = alpha_upp, 
      Rcpp::Named("alpha_lwr") = alpha_lwr,
      Rcpp::Named("delta_map") = delta_map
    );
  }
  
  else{
    map = Rcpp::List::create(
      Rcpp::Named("logposterior_map") = logposterior_map, 
      Rcpp::Named("phi_map") = phi_map, 
      Rcpp::Named("K_map") = K_map, 
      Rcpp::Named("lambda_map") = lambda_map, 
      Rcpp::Named("p_map") = p_map, 
      Rcpp::Named("alpha_map") = alpha_map,
      Rcpp::Named("delta_map") = delta_map
    );
  }
  
  if(store){
    store_list = Rcpp::List::create(
      Rcpp::Named("phi_store") = phi_store, 
      Rcpp::Named("K_store") = K_store, 
      Rcpp::Named("lambda_store") = lambda_store, 
      Rcpp::Named("p_store") = p_store, 
      Rcpp::Named("alpha_store") = alpha_store,
      Rcpp::Named("delta_store") = delta_store
    );
  }
  
  
  return Rcpp::List::create(Rcpp::Named("C_pred_mean") = C_mean, 
    Rcpp::Named("C_pred_upp") = C_upp, 
    Rcpp::Named("C_pred_lwr") = C_lwr, 
    Rcpp::Named("N_pred_mean") = N_mean, 
    Rcpp::Named("N_pred_upp") = N_upp, 
    Rcpp::Named("N_pred_lwr") = N_lwr, 
    Rcpp::Named("C_fit") = C_fit, 
    Rcpp::Named("N_fit") = N_fit, 
    Rcpp::Named("C0") = C0, 
    Rcpp::Named("accept") = accept,
    Rcpp::Named("map") = map,
    Rcpp::Named("iter") = iter,
    Rcpp::Named("store_list") = store_list,
    Rcpp::Named("logposterior_store") = logposterior_store,
    Rcpp::Named("CI_mat") = CI_mat,
    Rcpp::Named("ppm_store") = ppm_store,
    Rcpp::Named("delta_prop") = delta_prop);
}





// [[Rcpp::export]]
IntegerVector mcmc_moves(IntegerVector delta_vec, NumericVector w){
  int t, iter = 0;
  int T = delta_vec.length();
  CharacterVector method = {"g", "l", "s"}; // global swap, local swap, shift
  IntegerVector loc = Locator(delta_vec); // location of change points
  int M = loc.length();
  IntegerVector rownames = seq(0, T - 7); // all time points excluding the last 6
  IntegerVector loc_0 = setdiff(rownames, loc); // location of zeros
  IntegerVector loc_minus0 = Vector_equal(loc); // location of change points excluding the first one
  loc_minus0.erase(0); 
  CharacterVector sam_method = csample_char(method, 1, -1, w); // sample methods with weight w
  IntegerVector sam = sample(loc_minus0, 1, -1); // sample change point excluding the first one
  IntegerVector delta_vec_temp = Vector_equal(delta_vec);
  
  // Global Swap
  if( sam_method[0] == "g" ){
    delta_vec_temp(sam(0)) = 0;
    IntegerVector sam_2 = sample(loc_0, 1, -1);
    delta_vec_temp(sam_2(0)) = 1;
    IntegerVector loc_temp = Locator( delta_vec_temp );
    IntegerVector loc_diff = diff( loc_temp ); // distance between change points in terms of time t
    LogicalVector loc_logic = loc_diff < 7; 
    while(is_true(any(loc_logic))){ // continuous global swap until the distance between all change points is greater than 7
      IntegerVector delta_vec_temp2 = Vector_equal(delta_vec_temp);
      delta_vec_temp2(sam_2(0)) = 0;
      IntegerVector loc_temp = Locator(delta_vec_temp2);
      IntegerVector loc_0_temp = setdiff(rownames, loc_temp);
      sam_2 = sample(loc_0_temp, 1, -1);
      delta_vec_temp2(sam_2(0)) = 1;
      loc_temp = Locator(delta_vec_temp2);
      loc_diff = diff( loc_temp );
      loc_logic = loc_diff < 7;
      delta_vec_temp = Vector_equal(delta_vec_temp2);
      iter++;
      if ( iter > 20 ){ // break the loop if it gets stuck
        Rcout<< "Not enough available sample points \n";
        break;
      }
    }
  }
  // Local Swap
  else if( sam_method(0) == "l" ){
    NumericVector r = rbinom(1, 1, 0.5);
    if( r(0) == 1 ){
      delta_vec_temp[sam(0)] = 0;
      delta_vec_temp[sam(0) + 1] = 1;
      IntegerVector loc_temp = Locator( delta_vec_temp );
      IntegerVector loc_diff = diff( loc_temp );
      LogicalVector loc_logic = loc_diff < 7;
      int sam_plus = sam(0) + 1;
      while(is_true(any(loc_logic)) || sam_plus > T - 7){ // continuous global swap until the distance between all change points is greater than 7 
        IntegerVector delta_vec_temp2 = Vector_equal(delta_vec_temp);
        delta_vec_temp2(sam_plus) = 0;
        IntegerVector loc_temp = Locator(delta_vec_temp2);
        IntegerVector loc_0_temp = setdiff(rownames, loc_temp);
        IntegerVector sam_temp = sample(loc_0_temp, 1, -1);
        delta_vec_temp2(sam_temp(0)) = 1;
        loc_temp = Locator(delta_vec_temp2);
        sam_plus = sam_temp(0);
        loc_diff = diff( loc_temp );
        loc_logic = loc_diff < 7;
        delta_vec_temp = Vector_equal(delta_vec_temp2);
      } 
    }
    else if( r(0) == 0 ){
      delta_vec_temp[sam(0)] = 0;
      delta_vec_temp[sam(0) - 1] = 1;
      IntegerVector loc_temp = Locator( delta_vec_temp );
      IntegerVector loc_diff = diff( loc_temp );
      LogicalVector loc_logic = loc_diff < 7;
      int sam_minus = sam(0) - 1;
      while(is_true(any(loc_logic))){
        IntegerVector delta_vec_temp2 = Vector_equal(delta_vec_temp);
        delta_vec_temp2(sam_minus) = 0;
        IntegerVector loc_temp = Locator(delta_vec_temp2);
        IntegerVector loc_0_temp = setdiff(rownames, loc_temp);
        IntegerVector sam_temp = sample(loc_0_temp, 1, -1);
        delta_vec_temp2(sam_temp(0)) = 1;
        loc_temp = Locator(delta_vec_temp2);
        sam_minus = sam_temp(0);
        loc_diff = diff( loc_temp );
        loc_logic = loc_diff < 7;
        delta_vec_temp = Vector_equal(delta_vec_temp2);
      } 
    }
  }
  // Shift
  else if( sam_method(0) == "s" ){
    NumericVector sam_3 = rbinom(1, 1, 0.5);
    if( sam_3(0) == 1 || delta_vec(7) == 1 ){
      int j;
      for( j = 1; j < M; j++){
        delta_vec_temp(loc(j)) = 0;
        delta_vec_temp(loc(j) + 1) = 1;
      }
      IntegerVector loc_temp = Locator( delta_vec_temp );
      IntegerVector loc_diff = diff( loc_temp );
      LogicalVector loc_logic = loc_diff < 7;
      int sam_plus = loc_temp(M - 1);
      while(is_true(any(loc_logic)) || sam_plus > T - 7){
        IntegerVector delta_vec_temp2 = Vector_equal(delta_vec_temp);
        delta_vec_temp2(sam_plus) = 0;
        IntegerVector loc_temp = Locator(delta_vec_temp2);
        IntegerVector loc_0_temp = setdiff(rownames, loc_temp);
        IntegerVector sam_temp = sample(loc_0_temp, 1, -1);
        delta_vec_temp2(sam_temp(0)) = 1;
        loc_temp = Locator(delta_vec_temp2);
        sam_plus = sam_temp(0);
        loc_diff = diff( loc_temp );
        loc_logic = loc_diff < 7;
        delta_vec_temp = Vector_equal(delta_vec_temp2);
      } 
    }
    else if( sam_3(0) == 0 ){
      int j;
      for( j = 1; j < M; j++){
        delta_vec_temp(loc(j)) = 0;
        delta_vec_temp(loc(j) - 1) = 1;
      }
      IntegerVector loc_temp = Locator( delta_vec_temp );
      IntegerVector loc_diff = diff( loc_temp );
      LogicalVector loc_logic = loc_diff < 7;
      int sam_minus = loc_temp(1);
      while(is_true(any(loc_logic)) || sam_minus > T - 7 ){
        IntegerVector delta_vec_temp2 = Vector_equal(delta_vec_temp);
        delta_vec_temp2(sam_minus) = 0;
        IntegerVector loc_temp = Locator(delta_vec_temp2);
        IntegerVector loc_0_temp = setdiff(rownames, loc_temp);
        IntegerVector sam_temp = sample(loc_0_temp, 1, -1);
        delta_vec_temp2(sam_temp(0)) = 1;
        loc_temp = Locator(delta_vec_temp2);
        sam_minus = sam_temp(0);
        loc_diff = diff( loc_temp );
        loc_logic = loc_diff < 7;
        delta_vec_temp = Vector_equal(delta_vec_temp2);
      } 
    }
  }
  //delta_vec = Vector_equal(delta_vec_temp);
  return delta_vec_temp;
}


// [[Rcpp::export]]
CharacterVector csample_char( CharacterVector x, int size, bool replace, NumericVector prob ) {
  CharacterVector ret = sample(x, size, replace, prob) ;
  return ret ;
}  


// [[Rcpp::export]]
IntegerVector Locator(IntegerVector x) {
  IntegerVector loc (1);
  int t;
  for( t = 1; t < x.length(); t++ ){
    if( x[t] == 1 ){
      loc.push_back(t);
    }
  } 
  return loc;
}

// [[Rcpp::export]]
NumericVector sort_cpp(NumericVector v) {
  std::sort(v.begin(), v.end());
  return(v);
}

// [[Rcpp::export]]
IntegerVector sort_int(IntegerVector v) {
  std::sort(v.begin(), v.end());
  return(v);
}

// [[Rcpp::export]]
IntegerVector delta_randomize( int T, int M ) {
  int m;
  IntegerVector delta (T);
  IntegerVector t_seq = seq(7, T - 7);
  IntegerVector cp_sam = sample(t_seq, M - 1, -1);
  cp_sam = sort_int(cp_sam);
  delta(0) = 1;
  cp_sam.push_front(0);
  LogicalVector cp_logic = diff(cp_sam) < 7;
  while(is_true(any(cp_logic))){
    cp_sam = sample(t_seq, M - 1, -1);
    cp_sam = sort_int(cp_sam);
    cp_sam.push_front(0);
    cp_logic = diff(cp_sam) < 7;
  }
  for(m = 0; m < M; m++){
    delta(cp_sam(m)) = 1;
  }
  return(delta);
}
// [[Rcpp::export]]
IntegerVector Vector_equal(IntegerVector x) { // equate 2 vectors
  int i;
  IntegerVector y (x.length());
  for(i = 0; i < x.length(); i++ ){
    y[i] = x[i];
  }
  return y;
}  

// [[Rcpp::export]]
double lfactorial_cpp(int i) { //Stirling Approximation for log(n!)
  if( i <= 0){
    i = 1;
  }
  double temp;
  double pi = 3.141592653589793238462643383280;
  temp = i*log(i) - i + log(sqrt(2*pi*i));
  return temp;
}

// [[Rcpp::export]]
double ODE_fun(int K, double lambda, double p, double alpha, int C) { // g function
  
  lambda = lambda*pow(C, p)*(1.0 - 1.0*pow(1.0*C/K, alpha));
  if (lambda < 0.0) {
    lambda = 0.0;
  }
  return lambda;
}

double mean_cpp(NumericVector x) {
  int n = x.size(); // Size of vector
  double sum = 0; // Sum value
  // For loop, note cpp index shift to 0
  for(int i = 0; i < n; i++){
    // Shorthand for sum = sum + x[i]
    sum = sum + x(i);
  }
  return sum/n; // Obtain and return the Mean
}

