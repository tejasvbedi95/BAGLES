## Function to extract real-time data from covid-tracking-project api

# data.extract<-function(){
#   raw<-httr::GET("https://covidtracking.com/api/v1/states/daily.json")
#   data<-jsonlite::fromJSON(rawToChar(raw$content))
#   data$date<-as.Date(as.character(data[,1]),format="%Y%m%d")
#   data<-subset(data,select=c(date,state,positive,death,recovered,hospitalizedCurrently,hospitalizedCumulative,inIcuCurrently,inIcuCumulative,onVentilatorCurrently,onVentilatorCumulative))
#   data$state<-usdata::abbr2state(data$state)
#   data<-data[order(data$state),]
#   data[is.na(data)]<-0
#   return(data)
# }

set.data<-function(data,st){
  data<-subset(data,state==st,select=c(date,positive))
  #N<-pop.data["Texas",]
  date<-as.Date(data$date)
  data.dfr<-data[order(data$positive),]
  data<-data.dfr[,-1]
  #browser()
  names(data)<-data.dfr$date[order(date, decreasing = T)]
  return(data[which(data > 100)])
}

## SMEG simulator
g_fun <- function( C, lambda, p, alpha, K ){
  g <- lambda*(C^p)*(1 - (C/K)^alpha)
  return(g)
}

cp_simulator <- function( T, M, phi, lambda, p, alpha, K, delta, C0 = 100, seed = 7){
  set.seed(seed)
  C.vec <- DeltaC.vec <- I <- R <- c()
  gam <- 0.1
  DeltaC.vec[1] <- rnbinom(1, size = phi[1], mu = g_fun( C0, lambda[1], p[1], alpha[1], K[1] ))
  C.vec[1] <- DeltaC.vec[1] + C0
  I[1] <- C.vec[1]
  R[1] <- 0
  
  m = 1
  for( t in 2:T ){
    if( delta[t] == 1 ){
      m <- m + 1
    }
    DeltaC.vec[t] = rnbinom(1, size = phi[m], mu = g_fun( C.vec[t - 1], lambda[m], p[m], alpha[m], K[m] ))
    C.vec[t] = DeltaC.vec[t] + C.vec[t - 1]
    R[t] = R[t - 1] + ceiling(gam*I[t - 1])
    I[t] = I[t - 1] + (C.vec[t] - C.vec[t - 1]) - (R[t] - R[t - 1])
  }
  return(list(DeltaC = DeltaC.vec, C = C.vec, I = I, R = R))
}

## SMILES sim
sim_smiles <- function(T, N, delta, beta, gamma, phiS, phiR, I0, R0, seed = 1 ){
  set.seed(seed)
  S0 <- N - I0 - R0
  C0 <- I0
  S <- I <- R <- C <- DeltaC <- c()
  S[1] <- S0 - rnbinom(1, mu = beta[1]*(1/N)*S0*I0, size = phiS)
  R[1] <- R0 + rnbinom(1, mu = gamma*I0, size = phiR)
  I[1] <- N - S[1] - R[1]
  DeltaC[1] <- S0 - S[1]
  C[1] <- C0 + DeltaC[1]
  m <- 1
  for( t in 2:T){
    if(delta[t] == 1){
      m <- m + 1
    }
    S[t] <- S[t - 1] - rnbinom(1, mu = beta[m]*(1/N)*S[t - 1]*I[t - 1], size = phiS)
    R[t] <- R[t - 1] + rnbinom(1, mu = gamma*I[t - 1], size = phiR)
    I[t] <- N - S[t] - R[t]
    DeltaC[t] <- S[t - 1] - S[t]
    C[t] <- C[t - 1] + DeltaC[t]
  }
  return(list(S = S, I = I, R = R, C = C, DeltaC = DeltaC))
}

# Generate credible intervals of fits
cp_intervals <- function(C, T, M, phi, lambda, p, alpha, K, delta, C0 = 100, seed = 1){
  set.seed(seed)
  #DeltaC <- c(C0, diff(C))
  C.ll <- DeltaC.ll <- c()
  C.ul <- DeltaC.ul <- c()
  
  DeltaC.ll[1] <- quantile(rnbinom(1000, size = phi[1], mu = g_fun( C0, lambda[1], p[1], alpha[1], K[1] )), probs = 0.025, na.rm = T)
  DeltaC.ul[1] <- quantile(rnbinom(1000, size = phi[1], mu = g_fun( C0, lambda[1], p[1], alpha[1], K[1] )), probs = 0.975, na.rm = T)
  
  C.ll[1] <- DeltaC.ll[1] + C0
  C.ul[1] <- DeltaC.ul[1] + C0
  
  m = 1
  for( t in 2:T ){
    if( delta[t] == 1 ){
      m <- m + 1
    }
    DeltaC.ll[t] = quantile(rnbinom(1000, size = phi[m], mu = g_fun( C[t - 1], lambda[m], p[m], alpha[m], K[m] )), probs = 0.025, na.rm = T)
    DeltaC.ul[t] = quantile(rnbinom(1000, size = phi[m], mu = g_fun( C[t - 1], lambda[m], p[m], alpha[m], K[m] )), probs = 0.975, na.rm = T)
    
    C.ll[t] = DeltaC.ll[t] + C.ll[t - 1]
    C.ul[t] = DeltaC.ul[t] + C.ul[t - 1]
  }
  Intervals <- cbind.data.frame(DeltaC.ll = DeltaC.ll, DeltaC.ul = DeltaC.ul, C.ll = C.ll, C.ul = C.ul)
  return(Intervals)
}

## Clustering Measures

F.measure <- function(clust_true, clust_est){
  t <- length(clust_true)
  M1 <- length(unique(clust_true))
  M2 <- length(unique(clust_est)) ## May need different M for rjmcmc
  f.measure <- 0
  for(m1 in 1:M1){
    lm1 <- sum(clust_true == m1)
    f_temp <- 0
    f.vec <- c()
    for(m2 in 1:M2){
      lm2 <- sum(clust_est == m2)
      f_new <- sum((clust_true == m1)*(clust_est == m2))/(lm1 + lm2)
      if(f_new > f_temp){
        f_temp <- f_new
      }
    }
    f.vec[m1] <- f_temp
    f.measure <- sum(f.measure, (2/t)*(lm1*f.vec[m1]))
  }
  return(f.measure)
}

NVI <- function(clust_true, clust_est){
  t <- length(clust_true)
  M1 <- length(unique(clust_true))
  M2 <- length(unique(clust_est)) ## May need different M for rjmcmc
  e1 <- e2 <- mi <- 0
  for(m1 in 1:M1){
    lm1 <- sum(clust_true == m1)
    e1 <- sum(e1, lm1*log(lm1/t))
  }
  for(m2 in 1:M2){
    lm2 <- sum(clust_est == m2)
    e2 <- sum(e2, lm2*log(lm2/t))
  }
  
  for(m1 in 1:M1){
    lm1 <- sum(clust_true == m1)
    for(m2 in 1:M2){
      lm2 <- sum(clust_est == m2)
      cp <- sum((clust_true == m1)*(clust_est == m2))
      if(cp != 0){
        mi <- sum(mi, cp*log((cp*t)/(lm1*lm2)))
      }
    }
  }
  #browser()
  nvi <- -(1/(t*log(t)))*(e1 + e2 + 2*mi)
  return(nvi)
}

MI <- function(clust_true, clust_est){
  t <- length(clust_true)
  M1 <- length(unique(clust_true))
  M2 <- length(unique(clust_est)) ## May need different M for rjmcmc
  mi <- 0
  for(m1 in 1:M1){
    lm1 <- sum(clust_true == m1)
    for(m2 in 1:M2){
      lm2 <- sum(clust_est == m2)
      cp <- sum((clust_true == m1)*(clust_est == m2))
      if(cp != 0){
        mi <- sum(mi, cp*log((cp*t)/(lm1*lm2)))
      }
    }
  }
  mi <- mi/t
  return(mi)
}

NMI <- function(clust_true, clust_est){
  t <- length(clust_true)
  M1 <- length(unique(clust_true))
  M2 <- length(unique(clust_est)) ## May need different M for rjmcmc
  e1 <- e2 <- mi <- 0
  
  for(m1 in 1:M1){
    lm1 <- sum(clust_true == m1)
    e1 <- sum(e1, lm1*log(lm1/t))
  }
  for(m2 in 1:M2){
    lm2 <- sum(clust_est == m2)
    e2 <- sum(e2, lm2*log(lm2/t))
  }
  
  for(m1 in 1:M1){
    lm1 <- sum(clust_true == m1)
    for(m2 in 1:M2){
      lm2 <- sum(clust_est == m2)
      cp <- sum((clust_true == m1)*(clust_est == m2))
      if(cp != 0){
        mi <- sum(mi, cp*log((cp*t)/(lm1*lm2)))
      }
    }
  }
  nmi <- mi/(sqrt(e1 * e2))
  return(nmi)
}

# NVI(c(rep(1, 4), rep(3, 6), rep(2, 2)), c(rep(1, 4), rep(2, 4), rep(3, 4)))
# NMI(c(rep(1, 40), rep(2, 60), rep(3, 20)), c(rep(1, 10), rep(3, 90), rep(2, 20)))
# MI(c(rep(1, 40), rep(2, 60), rep(3, 20)), c(rep(1, 10), rep(3, 90), rep(2, 20)))


ARI <- function(clust_true, clust_est){
  T1 <- length(clust_true)
  T2 <- length(clust_est) ## May need different M for rjmcmc
  a <- b <- c <- d <- 0
  for(t1 in 2:T1){
    t2 <- 1
    while(t2 < t1){
      a <- sum(a, (clust_true[t1] == clust_true[t2])*(clust_est[t1] == clust_est[t2]))
      b <- sum(b, (clust_true[t1] == clust_true[t2])*(clust_est[t1] != clust_est[t2]))
      c <- sum(c, (clust_true[t1] != clust_true[t2])*(clust_est[t1] == clust_est[t2]))
      d <- sum(d, (clust_true[t1] != clust_true[t2])*(clust_est[t1] != clust_est[t2]))
      t2 <- t2 + 1
    }
  }
  ARI <- (choose(T1, 2)*(a + d) - ((a + b)*(a + c) + (c + d)*(b + d)))/(choose(T1, 2)^2 - ((a + b)*(a + c) + (c + d)*(b + d)))
  return(ARI)
}


segment_builder <- function(T, cp_vec){
  if(cp_vec[1] != 1){
    cp_vec <- c(1, cp_vec)
  }
  z <- c()
  for(m in 1:length(cp_vec)){
    if(m + 1 > length(cp_vec)){
      z_m <- rep(m, T - cp_vec[m] + 1)
    }else{
      z_m <- rep(m, cp_vec[m + 1] - cp_vec[m])
    }
    z <- c(z, z_m)
  }
  return(z)
}

# To compute delta CIs
CI_cp <- function(res, cp_vec){
  
  
  cp_vec <- cp_vec[-1]
  cp_mat <- NULL
  for(i in 1:length(cp_vec)){
    #browser()
    j <- k <- 1
    cp_id <- cp_vec[i]
    p_val_ll <- suppressWarnings(cor.test(res$CI_mat[,cp_id], res$CI_mat[,cp_id - j], alternative = "less")$p.value)
    if(p_val_ll < 0.05 & !is.na(p_val_ll)){
      ll <- cp_id - j
    }else{
      ll <- cp_id
    }
    while(p_val_ll < 0.05 & !is.na(p_val_ll)){
      j <- j + 1
      p_val_ll <- suppressWarnings(cor.test(res$CI_mat[,cp_id], res$CI_mat[,cp_id - j], alternative = "less")$p.value)
      ll <- cp_id - j
      if(ll == 1){
        break
      }
      #browser()
    }
    p_val_ul <- suppressWarnings(cor.test(res$CI_mat[,cp_id], res$CI_mat[,cp_id + k], alternative = "less")$p.value)
    while(p_val_ul < 0.05 & !is.na(p_val_ul)){
      #browser()
      k <- k + 1
      p_val_ul <- suppressWarnings(cor.test(res$CI_mat[,cp_id], res$CI_mat[,cp_id + k], alternative = "less")$p.value)
      ul <- cp_id + k
    }
    cp_mat <- rbind(cp_mat, cbind(cp_id, ll, ul))
  }
  return(cp_mat)
}


