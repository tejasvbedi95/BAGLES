#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]

using namespace Rcpp;
using namespace arma;

static double lfactorial_cpp(int i);
static double ldtrun_pois(int M, double lambda, int min, int max);
//static int max_cpp(int a, int b);
//static int min_cpp(int a, int b);
static NumericVector sort_cpp(NumericVector v);
static IntegerVector sort_int(IntegerVector v);
static double mean_cpp(NumericVector x);

static double ODE_fun(int K, double gamma, double p, int C);
static IntegerVector delta_randomize( int T, int M );
static CharacterVector csample_char( CharacterVector x, int size, bool replace, NumericVector prob );

// [[Rcpp::export]]
Rcpp::List growth_cp_rj(Col<int> DeltaC, int M_max, int POP, int T_fin, double alpha, bool store) {
  
  // Read data information
  //int T = C.length();
  int T = DeltaC.size();
  Col<int> C(T, fill::zeros);
  int M = 2, M_temp = 0;
  
  C(0) = 100; // (only for simulation study) Change to something else later
  for(int t = 1; t < T; t++){
    C(t) = C(t - 1) + DeltaC(t);
  }
  
  // Set algorithm settings
  int iter = 100000; 
  int burn = iter*0.5;
  
  // Starting Parameter Values & Temporary vars
  Col<int> delta_start(T);
  
  if(M == 1){
    delta_start.fill(0);
    delta_start(0) = 1;
  }else{
    delta_start = as<arma::Col<int>>(delta_randomize(T, M));
  }
  
  NumericVector w (5);
  
  double phi_start = 10.0;
  
  Col<double> lambda_start(M);
  lambda_start.fill(0.5);
  
  Col<double> p_start(M);
  p_start.fill(0.5);
  
  Col<double> K_start(M);
  NumericVector m_start (M);
  if(M == 1){
    K_start(0) = (max(C) + 1)*2;
  }
  else{
    m_start = wrap(arma::find(delta_start == 1));
    m_start.push_back(T);
    for(int m = 0; m < M; m++){
      if(m == M - 1){
        K_start(m) = (max(C) + 1)*2;
      }
      else{
        K_start(m) = (C(m_start(m + 1)) + 1)*2;
      }
    }
  }
  
  
  // Set hyperparameters
  double tau_logphi = 1.0, tau_logK = 1.0, tau_loglambda = 0.1, tau_logp = 0.1;
  double a_phi = 0.001, b_phi = 0.001; 
  double a_lambda = 0.001, b_lambda = 0.001;
  double a_p = 1.0, b_p = 1.0;
  int b_K = POP;
  
  Col<double> omega(T);
  omega(0) = 1.0;
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
  int it, i, m, m_temp, m_loc = 0, count = 0, count_2 = 0, count_3 = 0;
  double phi_temp = 0.0, p_h = 0.0;
  Col<double> lambda_temp = lambda_start;
  Col<double> p_temp = p_start;
  Col<double> K_temp = K_start;
  Col<int> delta_temp(T, fill::zeros);
  Col<int> delta_prop(T, fill::zeros);
  Col<int> z(T, fill::zeros);
  Col<double> K_birth(1, fill::zeros);
  Col<double> lambda_birth(1, fill::zeros);
  Col<double> p_birth(1, fill::zeros);
  double logpost_birth = 0, logpost_death = 0, log_dK = 0, log_dp = 0, log_dlambda = 0;
  
  
  double hastings, logposterior = 0, logposterior_map, logomega = 0, lognum, logden, logomega_temp;
  Col<double> ODE(T, fill::zeros); // g function
  Col<double> ODE_temp(T, fill::zeros);
  
  // Store
  Col<double> logposterior_store(iter, fill::zeros);
  Col<double> M_store(iter - burn, fill::zeros);
  //Col<double> M_store(iter, fill::zeros);
  Col<double> K_store(iter - burn, fill::zeros);
  // Initialization
  double phi = phi_start;
  Col<double> lambda = lambda_start;
  Col<double> K = K_start;
  Col<double> p = p_start;
  Col<int> delta = delta_start;
  
  double phi_map = phi;
  Col<double> lambda_map = lambda;
  Col<double> p_map = p;
  Col<double> K_map = K;
  int M_map = M;
  Col<int> delta_map = delta;
  
  double accept_phi = 0, accept_delta = 0, accept_b = 0, accept_d = 0, accept_g = 0, accept_l = 0;
  
  CharacterVector method  = {"b", "d", "g", "l", "s"}; // birth, death, global swap, local swap, stay
  IntegerVector loc_1 (M); // location of change points
  IntegerVector loc_1_temp (M); // location of change points
  IntegerVector loc_0 (T - M - 12); // location of legal non cps
  CharacterVector sam_method; // sample methods with weight w
  IntegerVector sam_1 (1); // sample change point 
  IntegerVector sam_0 (1); // sample non change point
  IntegerVector loc_diff (M - 2); // distance between change points in terms of time t
  LogicalVector loc_logic (M - 2); 
  Mat<int> K_expand(iter - burn - 1, T, fill::zeros);
  Mat<double> p_expand(iter - burn - 1, T, fill::zeros);
  Mat<double> lambda_expand(iter - burn - 1, T, fill::zeros);
  Mat<int> CI_mat(iter - burn, T, fill::zeros);
  IntegerMatrix C_predict(200, T_fin);
  IntegerMatrix N_predict(200, T_fin);
  //arma::mat ppm_store(T, T, arma::fill::zeros);
  //IntegerVector z_store(T);
  
  logposterior = 0;
  m = -1;
  // correct ODE fun as it should depend on C(t - 1) not C(t)
  
  for(int t = 0; t < T; t++){  
    if(delta(t) == 1){
      m++;
      logposterior = logposterior + (a_lambda - 1)*log(lambda(m)) - b_lambda*lambda(m);
      logposterior = logposterior + (a_p - 1)*log(p(m)) + (b_p - 1)*log(1 - p(m));
    }
    ODE(t) = ODE_fun(K(m), lambda(m), p(m), C(t));
    logposterior = logposterior + lgamma(DeltaC(t) + phi) - lgamma(phi) - lfactorial_cpp(DeltaC(t)) + phi*(log(phi) - log(ODE(t) + phi)) + DeltaC(t)*(log(ODE(t)) - log(ODE(t) + phi));
  }
  logposterior = logposterior + (a_phi - 1)*log(phi) - b_phi*phi;
  
  if(M > 1){
    for(int t = 7; t < T - 6; t++){ // bernoulli product prior for delta
      logomega = logomega + delta(t)*log(omega(t)) + (1 - delta(t))*log(1 - omega(t));
    }
    logposterior = logposterior + logomega;
  }
  logposterior = logposterior + ldtrun_pois(M, alpha, 2, M_max);
  
  // MCMC
  
  for (it = 0; it < iter; it++){
    
    // Update phi
    phi_temp = exp(r_truncnorm(log(phi), tau_logphi, log(1), log(100)));
    hastings = 0;
    
    for(int t = 0; t < T; t++){
      hastings = hastings + phi_temp*log(phi_temp) - lgamma(phi_temp) + lgamma(phi_temp + DeltaC(t)) - (phi_temp + DeltaC(t))*log(phi_temp + ODE(t));
      hastings = hastings - (phi*log(phi) - lgamma(phi) + lgamma(phi + DeltaC(t)) - (phi + DeltaC(t))*log(phi + ODE(t)));
    }
    
    hastings = hastings + (a_phi - 1)*log(phi_temp) - b_phi*phi_temp;
    hastings = hastings - ((a_phi - 1)*log(phi) - b_phi*phi);
    
    if(hastings >= log(double(rand()%10001)/10000)){
      phi = phi_temp;
      logposterior = logposterior + hastings;
      if (it > burn) {
        accept_phi = accept_phi + 1.0;
      }
    } 
    
    if( M == 2 ){
      //w = {1, 0, 0, 0, 3};
      //w = {1, 0, 0, 0, 0};
      w = {3, 0, 1, 1, 1};
    }
    else if(M == M_max){
      w = {0, 3, 1, 1, 1};
      //w = {0, 1, 0, 0, 0};
      
    }
    else{
      w = {3/2, 3/2, 1, 1, 1};
      //w = {0, 3, 1, 1, 1};
    }
    sam_method = csample_char(method, 1, 0, w)(0);
    
    loc_1 = wrap(find(delta == 1));
    loc_1.erase(0); // Excluding the first cp
    loc_0 = wrap(find(delta == 0));
    loc_0.erase(1, 6); // Excluding zeros at t = 2,...,7
    loc_0.erase(loc_0.length() - 6, loc_0.length()); // Excluding last 6 zeros
    
    delta_temp = delta;
    
    //Rcout<<sam_method<<"\n";
    // Birth Move
    if(sam_method(0) == "b"){
      M_temp = M + 1;
      
      sam_0 = sample(loc_0, 1, -1);
      delta_temp(sam_0(0)) = 1;
      
      
      loc_1_temp = wrap(find(delta_temp == 1));
      loc_diff = diff(loc_1_temp);
      loc_logic = loc_diff < 7;
      
      if(is_true(any(loc_logic))){
        delta_temp(sam_0(0)) = 0;
      }
      else{
        z = cumsum(delta);
        m_loc = z(sam_0(0)) - 1;
        
        K_birth = exp(r_truncnorm(log(K(m_loc)), tau_logK, log(C(sam_0(0) + 1)), log(b_K)));
        lambda_birth = exp(rnorm(1, log(lambda(m_loc)), tau_loglambda)(0));
        p_birth = exp(r_truncnorm(log(p(m_loc)), tau_logp, log(10e-9), log(1)));
        
        
        // Insert / delete rows isnt working as expected
        
        K_temp.insert_rows(m_loc + 1, K_birth);
        lambda_temp.insert_rows(m_loc + 1, lambda_birth);
        p_temp.insert_rows(m_loc + 1, p_birth);
        
        //K_temp(m_loc + 1) =  K_birth(0);
        //lambda_temp(m_loc + 1) =  lambda_birth(0);
        //p_temp(m_loc + 1) =  p_birth(0);
        
        
        //Rcout<<m_loc + 1<<"\n";
        //Rcout<<K<<"\n";
        //Rcout<<K_temp<<"\n";
        
        
        log_dK = dtruncnorm(wrap(log(K_birth)), log(K(m_loc)), tau_logK, log(C(sam_0(0) + 1)), log(b_K), 1)(0);
        log_dlambda = R::dnorm(as_scalar(log(lambda_birth)), log(lambda(m_loc)), tau_loglambda, 1);
        log_dp = dtruncnorm(wrap(log(p_birth)), log(p(m_loc)), tau_logp, log(10e-9), log(1), 1)(0);
        
        //log_dK = log_dK + log(10.0);
        //log_dlambda = log_dlambda + log(10.0);
        //log_dp = log_dp + log(10.0);
        
        logpost_birth = 0;
        hastings = 0;
        m = -1;
        m_temp = -1;
        
        for(int t = 0; t < T; t++){
          if(delta_temp(t) == 1){
            m_temp++;
            logpost_birth = logpost_birth + (a_lambda - 1)*log(lambda_temp(m_temp)) + b_lambda*lambda_temp(m_temp) + (a_p - 1)*log(p_temp(m_temp)) + (b_p - 1)*log(1 - p_temp(m_temp));
          }
          ODE_temp(t) = ODE_fun(K_temp(m_temp), lambda_temp(m_temp), p_temp(m_temp), C(t));
          if(delta(t) == 1){
            m++;
            logpost_birth = logpost_birth - ((a_lambda - 1)*log(lambda(m)) + b_lambda*lambda(m) + (a_p - 1)*log(p(m)) + (b_p - 1)*log(1 - p(m)));
          }
          logpost_birth = logpost_birth + DeltaC(t)*log(ODE_temp(t)) - (phi + DeltaC(t))*log(ODE_temp(t) + phi);
          logpost_birth = logpost_birth - (DeltaC(t)*log(ODE(t)) - (phi + DeltaC(t))*log(ODE(t) + phi));
        }
        logpost_birth = logpost_birth + ldtrun_pois(M_temp, alpha, 2, M_max) - ldtrun_pois(M, alpha, 2, M_max);
        
        hastings = log(T - M - 12) - log(M_temp - 1) + logpost_birth - (log_dK + log_dlambda + log_dp);
        
        if(hastings >= log(double(rand()%10001)/10000)){
          M = M_temp;
          lambda.set_size(M);
          K.set_size(M);
          p.set_size(M);
          
          lambda = lambda_temp;
          K = K_temp;
          p = p_temp;
          delta = delta_temp;
          logposterior = logposterior + logpost_birth;
          ODE = ODE_temp;
          if (it > burn) {
            accept_b = accept_b + 1;
          }
        }
        else{
          lambda_temp = lambda;
          K_temp = K;
          p_temp = p;
        }
        //Rcout<<M<<"\n";
        //Rcout<<p<<"\n";
      }
    }
    
    
    // Death Move
    if(sam_method(0) == "d"){
      M_temp = M - 1;
      
      sam_1 = sample(loc_1, 1, -1);
      delta_temp(sam_1(0)) = 0;
      
      //loc_1_temp = wrap(find(delta_temp == 1));
      
      z = cumsum(delta_temp);
      m_loc = z(sam_1(0)) - 1;
      
      
      K_birth = exp(r_truncnorm(log(K(m_loc)), tau_logK, log(C(sam_1(0) + 1)), log(b_K)));
      lambda_birth = exp(rnorm(1, log(lambda(m_loc)), tau_loglambda)(0));
      p_birth = exp(r_truncnorm(log(p(m_loc)), tau_logp, log(10e-9), log(1)));
      
      K_temp.shed_row(m_loc + 1);
      lambda_temp.shed_row(m_loc + 1);
      p_temp.shed_row(m_loc + 1);
      
      // Density estimate of K_birth is really small, figure something out
      // May have to try split/merge moves
      
      log_dK = dtruncnorm(  as_scalar(log(K_birth)),as_scalar(log(K(m_loc))),tau_logK, log(C(sam_0(0) + 1)), log(b_K))(0);
      /*
       if(log_dK == 0){
       log_dK = log(0.001);
       }*/
      log_dlambda = R::dnorm( as_scalar(log(lambda(m_loc))), as_scalar(log(lambda_birth)), tau_loglambda, 1);
      log_dp = d_truncnorm( log(p(m_loc)), as_scalar(log(p_birth)), tau_logp, log(10e-9), log(1), 1);
      
      //log_dK = log_dK + log(10.0);
      //log_dlambda = log_dlambda + log(10.0);
      //log_dp = log_dp + log(10.0);
      
      logpost_death = 0;
      hastings = 0;
      m = -1;
      m_temp = -1;
      for(int t = 0; t < T; t++){
        if(delta_temp(t) == 1){
          m_temp++;
          logpost_death = logpost_death + (a_lambda - 1)*log(lambda_temp(m_temp)) + b_lambda*lambda_temp(m_temp) + (a_p - 1)*log(p_temp(m_temp)) + (b_p - 1)*log(1 - p_temp(m_temp));
        }
        ODE_temp(t) = ODE_fun(K_temp(m_temp), lambda_temp(m_temp), p_temp(m_temp), C(t));
        if(delta(t) == 1){
          m++;
          logpost_death = logpost_death - ((a_lambda - 1)*log(lambda(m)) + b_lambda*lambda(m) + (a_p - 1)*log(p(m)) + (b_p - 1)*log(1 - p(m)));
        }
        logpost_death = logpost_death + DeltaC(t)*log(ODE_temp(t)) - (phi + DeltaC(t))*log(ODE_temp(t) + phi);
        logpost_death = logpost_death - (DeltaC(t)*log(ODE(t)) - (phi + DeltaC(t))*log(ODE(t) + phi));
      }
      logpost_death = logpost_death + ldtrun_pois(M_temp, alpha, 2, M_max) - ldtrun_pois(M, alpha, 2, M_max);
      
      hastings = log(M - 1) - log(T - M_temp - 12) + logpost_death + log_dK + log_dlambda + log_dp;
      
      if(hastings >= log(double(rand()%10001)/10000)){
        M = M_temp;
        lambda.set_size(M);
        K.set_size(M);
        p.set_size(M);
        
        lambda = lambda_temp;
        K = K_temp;
        p = p_temp;
        delta = delta_temp;
        logposterior = logposterior + logpost_death;
        ODE = ODE_temp;
        if (it > burn) {
          accept_d = accept_d + 1;
        }
      }
      else{
        lambda_temp = lambda;
        K_temp = K;
        p_temp = p;
      }
      //Rcout<<M<<"\n";
      //Rcout<<p<<"\n";
    }
    
    
    // Update delta for fixed M
    if(sam_method(0) == "g"){
      // Global swap
      sam_1 = sample(loc_1, 1, -1);
      sam_0 = sample(loc_0, 1, -1);
      delta_temp(sam_1(0)) = 0;
      delta_temp(sam_0(0)) = 1;
      
      loc_1_temp = wrap(find(delta_temp == 1));
      loc_diff = diff(loc_1_temp);
      loc_logic = loc_diff < 7;
      if(is_true(any(loc_logic))){
        delta_temp(sam_1(0)) = 1;
        delta_temp(sam_0(0)) = 0;
      }
      else{
        hastings = 0;
        m_temp = -1;
        for(int t = 0; t < T; t++){
          if(delta_temp(t) == 1){
            m_temp++;
          }
          ODE_temp(t) = ODE_fun(K(m_temp), lambda(m_temp), p(m_temp), C(t));
          hastings = hastings + lgamma(DeltaC(t) + phi) - lgamma(phi) - lfactorial_cpp(DeltaC(t)) + phi*(log(phi) - log(ODE_temp(t) + phi)) + DeltaC(t)*(log(ODE_temp(t)) - log(ODE_temp(t) + phi));
          hastings = hastings - (lgamma(DeltaC(t) + phi) - lgamma(phi) - lfactorial_cpp(DeltaC(t)) + phi*(log(phi) - log(ODE(t) + phi)) + DeltaC(t)*(log(ODE(t)) - log(ODE(t) + phi)));
        }
        if(hastings >= log(double(rand()%10001)/10000)){
          delta = delta_temp;
          ODE = ODE_temp;
          logposterior = logposterior + hastings;
          if (it > burn) {
            accept_delta++;
            accept_g++;
          }
        }
      }
      
    }
    
    if(sam_method(0) == "l"){
      // Local Swap
      if(runif(1)(0) > 0.5){
        sam_1 = sample(loc_1, 1, -1);
        delta_temp(sam_1(0)) = 0;
        delta_temp(sam_1(0) + 1) = 1;
        
        loc_1_temp = wrap(find(delta_temp == 1));
        loc_diff = diff(loc_1_temp);
        loc_logic = loc_diff < 7;
        
        if(is_true(any(loc_logic)) || sam_1(0) > T - 8){
          delta_temp(sam_1(0)) = 1;
          delta_temp(sam_1(0) + 1) = 0;
        }
        else{
          hastings = 0;
          m_temp = -1;
          for(int t = 0; t < T; t++){
            if(delta_temp(t) == 1){
              m_temp++;
            }
            ODE_temp(t) = ODE_fun(K(m_temp), lambda(m_temp), p(m_temp), C(t));
            hastings = hastings + lgamma(DeltaC(t) + phi) - lgamma(phi) - lfactorial_cpp(DeltaC(t)) + phi*(log(phi) - log(ODE_temp(t) + phi)) + DeltaC(t)*(log(ODE_temp(t)) - log(ODE_temp(t) + phi));
            hastings = hastings - (lgamma(DeltaC(t) + phi) - lgamma(phi) - lfactorial_cpp(DeltaC(t)) + phi*(log(phi) - log(ODE(t) + phi)) + DeltaC(t)*(log(ODE(t)) - log(ODE(t) + phi)));
          }
          if(hastings >= log(double(rand()%10001)/10000)){
            delta = delta_temp;
            ODE = ODE_temp;
            logposterior = logposterior + hastings;
            if (it > burn) {
              accept_delta++;
              accept_l++;
            }
          }
        }
      }
      else{
        sam_1 = sample(loc_1, 1, -1);
        delta_temp(sam_1(0)) = 0;
        delta_temp(sam_1(0) - 1) = 1;
        
        loc_1_temp = wrap(find(delta_temp == 1));
        loc_diff = diff(loc_1_temp);
        loc_logic = loc_diff < 7;
        
        if(is_true(any(loc_logic)) || sam_1(0) < 8){
          delta_temp(sam_1(0)) = 1;
          delta_temp(sam_1(0) - 1) = 0;
        }
        else{
          hastings = 0;
          m_temp = -1;
          for(int t = 0; t < T; t++){
            if(delta_temp(t) == 1){
              m_temp++;
            }
            ODE_temp(t) = ODE_fun(K(m_temp), lambda(m_temp), p(m_temp), C(t));
            hastings = hastings + lgamma(DeltaC(t) + phi) - lgamma(phi) - lfactorial_cpp(DeltaC(t)) + phi*(log(phi) - log(ODE_temp(t) + phi)) + DeltaC(t)*(log(ODE_temp(t)) - log(ODE_temp(t) + phi));
            hastings = hastings - (lgamma(DeltaC(t) + phi) - lgamma(phi) - lfactorial_cpp(DeltaC(t)) + phi*(log(phi) - log(ODE(t) + phi)) + DeltaC(t)*(log(ODE(t)) - log(ODE(t) + phi)));
          }
          if(hastings >= log(double(rand()%10001)/10000)){
            delta = delta_temp;
            ODE = ODE_temp;
            logposterior = logposterior + hastings;
            if (it > burn) {
              accept_delta++;
              accept_l++;
            }
          }
        }
      }
    }
    
    
    if(sam_method(0) == "g" || sam_method(0) == "l" || sam_method(0) == "s"){
      loc_1 = wrap(find(delta == 1));
      loc_1.push_back(T);
      
      for(m = 0; m < M; m++){
        
        // Update K
        
        if(m < M - 1 && M > 1){
          K_temp(m) = exp(r_truncnorm(log(K(m)), tau_logK, log(C(loc_1(m + 1))), log(b_K)));
        }
        else if(m == M - 1){
          K_temp(m) = exp(r_truncnorm(log(K(m)), tau_logK, log(1.0*max(C)), log(b_K)));
        }
        
        hastings = 0;
        
        for(int t = loc_1(m); t < loc_1(m + 1); t++){
          ODE_temp(t) = ODE_fun(K_temp(m), lambda(m), p(m), C(t));
          hastings = hastings + DeltaC(t)*log(ODE_temp(t)) - (phi + DeltaC(t))*log(ODE_temp(t) + phi);
          hastings = hastings - (DeltaC(t)*log(ODE(t)) - (phi + DeltaC(t))*log(ODE(t) + phi));
        }
        if(hastings >= log(double(rand()%10001)/10000)){
          K(m) = K_temp(m);
          logposterior = logposterior + hastings;
          for(int t = loc_1(m); t < loc_1(m + 1); t++){
            ODE(t) = ODE_temp(t);
          }
          //    if (it > burn) {
          //     accept_K(m) = accept_K(m) + 1;
          //   }
        }
        
        // Update lambda
        
        lambda_temp(m) = exp(rnorm(1, log(lambda(m)), tau_loglambda)(0));
        
        hastings = 0;
        
        for(int t = loc_1(m); t < loc_1(m + 1); t++){
          ODE_temp(t) = ODE_fun(K(m), lambda_temp(m), p(m), C(t));
          hastings = hastings + DeltaC(t)*log(ODE_temp(t)) - (phi + DeltaC(t))*log(ODE_temp(t) + phi);
          hastings = hastings - (DeltaC(t)*log(ODE(t)) - (phi + DeltaC(t))*log(ODE(t) + phi));
        }
        if(hastings >= log(double(rand()%10001)/10000)){
          lambda(m) = lambda_temp(m);
          logposterior = logposterior + hastings;
          for(int t = loc_1(m); t < loc_1(m + 1); t++){
            ODE(t) = ODE_temp(t);
          }
          // if (it > burn) {
          //    accept_lambda(m) = accept_lambda(m) + 1;
          //  }
        }
        
        // Update p
        
        p_temp(m) = exp(r_truncnorm(log(p(m)), tau_logp, log(10e-9), log(1)));
        
        hastings = 0;
        
        for(int t = loc_1(m); t < loc_1(m + 1); t++){
          ODE_temp(t) = ODE_fun(K(m), lambda(m), p_temp(m), C(t));
          hastings = hastings + DeltaC(t)*log(ODE_temp(t)) - (phi + DeltaC(t))*log(ODE_temp(t) + phi);
          hastings = hastings - (DeltaC(t)*log(ODE(t)) - (phi + DeltaC(t))*log(ODE(t) + phi));
        }
        if(hastings >= log(double(rand()%10001)/10000)){
          p(m) = p_temp(m);
          logposterior = logposterior + hastings;
          for(int t = loc_1(m); t < loc_1(m + 1); t++){
            ODE(t) = ODE_temp(t);
          }
          //  if (it > burn) {
          //    accept_p(m) = accept_p(m) + 1;
          //  }
        }
      }
    }
    
    
    // Monitor the process
    
    if(it*100/iter == count){
      Rcout<<count<< "% has been done\n";
      count = count + 10;
    }
    if(store){
      logposterior_store(it) = logposterior;
      if(it >= burn){
        M_store(it - burn) = M;
        K_store(it - burn) = K(M - 1);
      }
      //if(it >= burn){
      //  M_store(it - burn) = M;
      //}
    }
    
    
    // Obtaining MAP estimates
    if( it >= burn){
      if( it == burn){
        logposterior_map = logposterior;
        delta_map = delta;
        phi_map = phi;
        K_map = K;
        lambda_map = lambda;
        p_map = p;
        M_map = M;
      }
      else{
        if(logposterior > logposterior_map){
          logposterior_map = logposterior;
          delta_map = delta;
          phi_map = phi;
          K_map = K;
          lambda_map = lambda;
          p_map = p;
          M_map = M;
        }
      }
      
      // Obtaining PPI
      if(store){
        for(int t = 0; t < T; t++){
          if(delta(t) == 1){
            delta_prop(t) = delta_prop(t) + 1.0;
            CI_mat(it - burn, t) = 1;
          }
          // if(delta_map(t) == 1){
          //   m++;
          // }
          // K_expand(it - burn - 1, t) = K_map(m);
          // p_expand(it - burn - 1, t) = p_map(m);
          // lambda_expand(it - burn - 1, t) = lambda_map(m);
        }
        
      }
      
      // Obtaining PPM 
      // if(store){
      //   if(it >= burn){
      //     z_store(0) = delta(0);
      //     for(int t = 0; t < T - 1; t++){
      //       z_store(t + 1) = z_store(t) + delta(t + 1);
      //     }
      //     for(int t = 0; t < T - 1; t++){
      //       for(int tt = 1; tt < T; tt++ ){
      //         if(z_store[t] == z_store[tt]){
      //           ppm_store(t, tt) = ppm_store(t, tt) + 1;
      //           ppm_store(tt, t) = ppm_store(tt, t) + 1;
      //         }
      //       }
      //     }
      //   }
      // }
      
      // Prediction
      int t;
      if (it % 250 == 0)
      {
        for (t = 0; t < T_fin; t++) 
        {
          if (t == 0)
          {
            N_predict(count_2, t) = rnbinom_mu(1, phi, ODE_fun(K(M - 1), lambda(M - 1), p(M - 1), C(T - 1)))(0);
            C_predict(count_2, t) = C(T - 1) + N_predict(count_2, t);
          }
          else
          {
            N_predict(count_2, t) = rnbinom_mu(1, phi, ODE_fun(K(M - 1), lambda(M - 1), p(M - 1), C_predict(count_2, t - 1)))(0);
            C_predict(count_2, t) = C_predict(count_2, t - 1) + N_predict(count_2, t);
            
          }
        }
        count_2 = count_2 + 1;
      }
      
    }
    
    
  }
  
  // Update C0 
  
  int c, C0 = 1;
  double lambda_0, logprob_max, logprob, ODE_0;
  for (c = 1; c < C(0); c++)
  {
    ODE_0 = ODE_fun(K_map[0], lambda_map[0], p_map[0], c);
    if (c == 1)
    {  // Would it be the parameters of first subsegment or last??
      logprob_max = lgamma(C(0) - C0 + phi_map) - lfactorial_cpp(C(0) - C0) - phi_map*log(ODE_0 + phi_map) + (C(0) - C0)*(log(ODE_0) - log(ODE_0 + phi_map));
    }
    else
    {
      logprob = lgamma(C(0) - C0 + phi_map) - lfactorial_cpp(C(0) - C0) - phi_map*log(ODE_0 + phi_map) + (C(0) - C0)*(log(ODE_0) - log(ODE_0 + phi_map));
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
  for (int t = 0; t < T + T_fin; t++)
  {
    if (t == 0)
    {
      N_fit(t) = ODE_fun(K_map[0], lambda_map[0], p_map[0], C0);
      C_fit(t) = C0 + N_fit(t);
    }
    else
    {
      if(t < T){
        if( delta_map(t) == 1){
          m++;
        }
      }

      N_fit(t) = ODE_fun(K_map[m], lambda_map[m], p_map[m], C_fit(t - 1));
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
  NumericVector vtemp;
  
  for (int t = 0; t < T_fin; t++) 
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
  
  // Result wrap-up
  
  Rcpp::List accept = Rcpp::List::create(
    Rcpp::Named("accept_delta") = accept_delta,
    Rcpp::Named("accept_phi") = accept_phi,
    Rcpp::Named("accept_b") = accept_b,
    Rcpp::Named("accept_d") = accept_d,
    Rcpp::Named("accept_g") = accept_g, 
    Rcpp::Named("accept_l") = accept_l  
  );
  
  Rcpp::List map;
  Rcpp::List store_list;
  
  map = Rcpp::List::create(
    Rcpp::Named("logposterior_map") = logposterior_map, 
    Rcpp::Named("M_map") = M_map,
    Rcpp::Named("phi_map") = phi_map, 
    Rcpp::Named("K_map") = K_map, 
    Rcpp::Named("lambda_map") = lambda_map, 
    Rcpp::Named("p_map") = p_map, 
    Rcpp::Named("delta_map") = trans(delta_map)
  );
  
  
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
    //Rcpp::Named("ppm_store") = ppm_store,
    Rcpp::Named("delta_prop") = delta_prop,
    Rcpp::Named("M_store") = M_store,
    Rcpp::Named("K_store") = K_store,
    // Rcpp::Named("K_expand") = K_expand,
    // Rcpp::Named("p_expand") = p_expand,
    // Rcpp::Named("lambda_expand") = lambda_expand,
    Rcpp::Named("CI_mat") = CI_mat,
    Rcpp::Named("logposterior_store") = logposterior_store
  );
  
  //Rcout<<hastings<<"\n";
  //Rcout<<log_dK<<"\n";
  //Rcout<<log_dp<<"\n";
  //Rcout<<log_dlambda<<"\n";
  //Rcout<<log(lambda(m_loc))<<"\n";
  //Rcout<<log(lambda_birth)<<"\n";
  //Rcout<<log(lambda_birth)<<"\n";
  //Rcout<<sam_method<<"\n";
}

// [[Rcpp::export]]
CharacterVector csample_char( CharacterVector x, int size, bool replace, NumericVector prob ) {
  CharacterVector ret = sample(x, size, replace, prob) ;
  return ret ;
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
double ODE_fun(int K, double lambda, double p, int C) { // g function
  double ODE;
  ODE = lambda*pow(C, p)*(1.0 - 1.0*(1.0*C/K));
  if (ODE <= 0.0) {
    ODE = 0.001;
  }
  return ODE;
}

// [[Rcpp::export]]
double ldtrun_pois(int M, double lambda, int min, int max) { // g function
  double d;
  d = R::dpois(M, lambda, true) - log(R::ppois(max, lambda, 1, false) - R::ppois(min - 1, lambda, 1, false));
  return d;
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













