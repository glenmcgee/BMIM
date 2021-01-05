#include "RcppArmadillo.h"
using namespace Rcpp;
#include <bsmimheader.h> // to import GIG


// draw interval-truncated normal (truncated such that it cannot be between lo and hi)
double rinttruncnorm(const double&     mu,               // mean of normal
                     const double&     sig,              // sd of normal
                     const double&     lo,               // lower boundary of truncated interval
                     const double&     hi) {             // upper boundary of truncated interval

  double x;
  
  x = rnorm(1,mu,sig)(0);
  
  while(x>lo & x<hi){
    x = rnorm(1,mu,sig)(0);
  }
  
  return x; 
}

// log density of interval-truncated normal
double logdinttruncnorm(const double&     x,                // value
                        const double&     mu,               // mean of normal
                        const double&     sig,              // sd of normal
                        const double&     lo,               // lower boundary of truncated interval
                        const double&     hi) {             // upper boundary of truncated interval
  double logdens = 0;
  
  logdens = R::dnorm(x,mu,sig,TRUE) 
    -log( 1.0 - (R::pnorm(hi,mu,sig,TRUE,FALSE)-R::pnorm(lo,mu,sig,TRUE,FALSE)) );
  
  return logdens; 
}



// Construct Sigma, A1, A2, C1, C2
arma::field<arma::mat> get_Vcomps(const bool&       poly,               // 0=gaussian kernel / 1=polynomial kernel
                                  const int&        d,                  // degree of polynomial kernel
                                  const int&        N,                  // no. of observations
                                  const int&        P_z,                // no. of covariates
                                  const arma::mat&  yz,                 // [y,z]
                                  const double&     logLambdaInverse,   // log(lamda^{-1})
                                  const arma::mat&  XthetaStar,         // X times theta*
                                  const bool&       randint,            // indicator for random intercepts model
                                  const double&     logLambdaBInverse,  // log(lamda^{-1}_B)
                                  const arma::mat&  Bmat) {             // block diagonal matrix of cluster membership
  
  // Construct \lambda^{-1}*K where K is the kernel matrix
  arma::mat laminvK(N,N,arma::fill::zeros);
  
  if(poly==TRUE){ // polynomial kernel
    for(int i=0; i<N; i++){
      for(int j=i; j<N; j++){
        laminvK(i,j) = exp(logLambdaInverse)*pow(1 + as_scalar(XthetaStar.row(i)*XthetaStar.row(j).t()), d); 
        laminvK(j,i) = laminvK(i,j); // symmetry
      }  
      laminvK(i,i) = laminvK(i,i); // 
    }
    
  }else{ // gaussian kernel
    for(int i=0; i<N-1; i++){
      for(int j=i+1; j<N; j++){
        laminvK(i,j) = exp(logLambdaInverse - arma::sum(square(vectorise(XthetaStar.row(i)-XthetaStar.row(j)))));
        laminvK(j,i) = laminvK(i,j); // symmetry
      }  
      laminvK(i,i) = exp(logLambdaInverse); // 
    }
    laminvK(N-1,N-1) = exp(logLambdaInverse); // 
  }
  
  // Construct Sigma=I+\lambda^{-1}*K+\lambda_B^{-1}B where K is the kernel matrix
  arma::mat Sigma(N,N,arma::fill::eye);           // initialize Sigma = I
  Sigma = Sigma + laminvK;                        // Sigma = I + \lambda^{-1}K
  if(randint==TRUE){                              // Add \lambda_B^{-1}B for random intercepts model
    Sigma = Sigma + exp(logLambdaBInverse)*Bmat;  // Sigma = I + \lambda^{-1}K + \lambda_B^{-1}B
  } 
  
  // Computing components for posterior (ll)
  arma::mat C1 = arma::chol(Sigma, "lower");                                  // sqrt of Sigma (taking cholesky decomposition loosely to mean sqrt)
  arma::mat Cyz = arma::solve(trimatl(C1),yz);                                // sqrt of [y,Z]^T * Sigma^{-1} * [y,Z]
  arma::vec B = Cyz.cols(1,P_z).t() * Cyz.col(0);                             // B = z^T * Sigma^{-1} * y
  arma::mat C2 = arma::chol(Cyz.cols(1,P_z).t() * Cyz.cols(1,P_z), "lower");  // C2 = sqrt of z^T * Sigma^{-1} * z
  arma::vec C2B = vectorise(arma::solve(trimatl(C2),B));                      // C2B = C2 * B = sqrt of A2
  
  arma::mat A(2,1);                                                           // A=(A1,A2)^T
  A(0,0) = arma::sum(square(vectorise(Cyz.col(0))));                          // A1 = y^T * Sigma^{-1} * y
  A(1,0) = arma::sum(square(C2B));                                            // A2 = y^T Sigma^{-1} z (z^T Sigma^{-1} z)^{-1} z^T Sigma^{-1} y
  // double A1 = arma::sum(square(vectorise(Cyz.col(0))));                       // A1 = y^T * Sigma^{-1} * y
  // double A2 = arma::sum(square(C2B));                                         // A2 = y^T Sigma^{-1} z (z^T Sigma^{-1} z)^{-1} z^T Sigma^{-1} y
  
  arma::field<arma::mat> comps(6);
  comps[0] = Sigma;                           // Sigma=I+\lambda^{-1}*K where K is the kernel matrix
  comps[1] = A;                               // A=(A1,A2)^T // A1 = y^T * Sigma^{-1} * y // A2 = y^T Sigma^{-1} z (z^T Sigma^{-1} z)^{-1} z^T Sigma^{-1} y
  comps[2] = B;                               // B = z^T * Sigma^{-1} * y
  comps[3] = C1;                              // sqrt of Sigma (taking cholesky decomposition loosely to mean sqrt)
  comps[4] = C2;                              // C2 = sqrt of z^T * Sigma^{-1} * z
  comps[5] = laminvK;                         // \lambda^{-1}*K where K is the kernel matrix
  
  return comps; 
}


// [[Rcpp::depends(RcppArmadillo)]]
//' MCMC for BSMIM
//' 
//' Returns posterior draws for \eqn{\theta},\eqn{\lambda^{-1}},\eqn{\rho},\eqn{\sigma^2},\eqn{\gamma}, as well as the mean(\eqn{h}) and cov(\eqn{h})
//' Adapted from "bkmrdlm_multi.cpp" in the "regimes" package by Ander Wilson.
//' Note that this implementation with horseshoe2 prior returns posterior draws of \eqn{\rho} as a list of M matrices (unlike main MCMC code, which returns one matrix of M columns)
//' 
//' @param yz matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
//' @param Xlist list of (N by L_m) matrices representing X_m
//' @param b_lambda hyperparameter for \eqn{\lambda^{-1}}
//' @param b_lambdaB hyperparameter for \eqn{\lambda^{-1}_B}
//' @param a_sig first hyperparameter for \eqn{\sigma^{-2}}
//' @param b_sig second hyperparameter for \eqn{\sigma^{-2}}
//' @param tau02 hyperparameter for tau2; classic horseshoe is 1
//' @param kappa vector of hyperparameters \eqn{\kappa_m} for \eqn{\theta^*_m}
//' @param poly 0 = gaussian kernel, 1 = polynomial kernel
//' @param d degree of polynomial kernel
//' @param horseshoe 0 = no selection, 1 = componentwise horseshoe priors
//' @param randint 0 = no random intercepts, 1 = random intercepts model
//' @param Bmat N xN block diagonal matrix indicating cluster membership for random intercepts model
//' @param draw_h 0 = dont draw h, 1 = draw h
//' @param n_inner no. of MCMC iterations to run in the inner loop. n_outer*n_inner iteraction will be run.
//' @param n_outer no. of MCMC iterations to run in the outer loop. n_outer iterations will be saved.
//' @param n_burn no. of MCMC iterations to discard as burn-in
//' @author Glen McGee and Ander Wilson (adapted from the "regimes" package by Ander Wilson).
//' @export
// [[Rcpp::export]]
List bsmim_mcmc2(const arma::mat&    yz,         // matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
                 const Rcpp::List&   Xlist,      // list of (N by L_m) matrices representing X_m
                 const double&       b_lambda,   // hyperparameter for \lambda^{-1}
                 const double&       b_lambdaB,  // hyperparameter for \lambda^{-1}_B
                 const double&       a_sig,      // first hyperparameter for \sigma^{-2}
                 const double&       b_sig,      // second hyperparameter for \sigma^{-2}
                 const double&       tau02,      // hyperparameter for tau2
                 const arma::vec&    kappa,      // vector of hyperparameters kappa_m for \theta^*_m
                 const bool&         poly,       // 0=gaussian kernel / 1=polynomial kernel
                 const int&          d,          // degree of polynomial kernel
                 const bool&         horseshoe,  // 0=no variable selection / 1=horseshoe prior
                 const bool&         randint,    // 0=no random intercepts / 1=random intercepts
                 const arma::mat&    Bmat,       // block diagonal matrix B indicating cluster membership for random intercepts model
                 const bool&         draw_h,     // 0=dont sample h, 1= sample h
                 const int&          n_inner,    // no. of MCMC iterations to run in the inner loop. n_outer*n_inner iteraction will be run.
                 const int&          n_outer,    // no. of MCMC iterations to run in the outer loop. n_outer iterations will be saved.
                 const int&          n_burn) {   // no. of MCMC iterations to discard as burn-in
  
  // step size
  double step_tau = 1;
  double step_nu = 1;
  
  // Dimensions
  int N = yz.n_rows;              // no. of observations
  int M = Xlist.size();           // no. of indices (dimension of h(.))
  int P_z = yz.n_cols-1;          // no. of covariates
  IntegerVector Lvec = rep(0,M);  // vector of L_m representing the no. of columns of X_m
  
  // Initialize
  arma::field<arma::mat> X(M);                    // field of (N by L_m) matrices representing X_m
  arma::field<arma::vec> thetaStar(M);            // field of L_m vectors representing \thetaStar_m
  arma::field<arma::vec> lognu(M);                // field of L_m vectors representing log \nu_{ml}
  arma::mat  XthetaStar(N,M,arma::fill::zeros);   // (N by M) matrix representing [X_1\thetaStar_1,...,X_M\thetaStar_M]
  for(int m=0; m<M; m++){
    X[m] = as<arma::mat>(Xlist[m]);                 // mth matrix in X is X_m
    Lvec(m) = X[m].n_cols;                          // mth element of Lvec is L_m (the length of \thetaStar_m)
    lognu[m] = as<arma::vec>(rnorm(Lvec(m)));       // set starting value for lognu_m
    thetaStar[m] = as<arma::vec>(rnorm(Lvec(m)));   // set starting value for \thetaStar_m
    XthetaStar.col(m) = X[m] * thetaStar[m];
  }
  
  double logLambdaInverse = rnorm(1)(0)*sqrt(b_lambda)/10;    // set starting value for log(\lambda^{-1})
  double logLambdaBInverse = 0;                               // set log(\lambda^{-1}_B) to 0 if no random intercepts
  if(randint==TRUE){
    logLambdaBInverse = rnorm(1)(0)*sqrt(b_lambdaB)/10;       // set starting value for log(\lambda^{-1}_B) if using random intercepts
  }
  double logtau = 0;                                          // set log(tau) to 0 if no random intercepts
  if(horseshoe==TRUE){
    logtau = rcauchy(1)(0)*sqrt(tau02)/10;                    // set starting value for log(tau) if using horseshoe
  }
  
  
  // Construct Sigma=I+\lambda^{-1}*K, A1, A2, C1, C2
  arma::field<arma::mat> Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat);
  
  // First evaluation of the log-liklihood (ll)
  double ll = -logLambdaInverse*logLambdaInverse/(2*b_lambda) 
    - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
    - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
      
   
  // Initialize other values
  double logLambdaInverse_PROP, angle, amax, amin, llu, ll_PROP, logratio;
  double lltau, logtau_PROP, lltau_PROP, temp_sumthet;
  double llnu, lognu_PROP, llnu_PROP;
  double llB, logLambdaBInverse_PROP, llB_PROP;
  arma::vec ons(N, arma::fill::ones);           // vector of length N of ones
  arma::vec vec_PROP;                           // proposal vector to define ellipse
  arma::vec thetaStar_PROP;                     // proposal vector for thetaStar

  
  // Store posterior samples
  arma::vec logLambdaInversekeep(n_outer,arma::fill::zeros);  // vector to store posterior draws of lambda^{-1}
  arma::vec logLambdaBInversekeep(n_outer,arma::fill::zeros); // vector to store posterior draws of lambda_B^{-1}
  arma::vec logtaukeep(n_outer,arma::fill::zeros);            // vector to store posterior draws of tau
  arma::mat rhokeep(n_outer,M,arma::fill::zeros);             // matrix to store posterior draws of rho // each row is rho_1,...,rho_m
  // arma::mat nukeep(n_outer,M,arma::fill::zeros);           // matrix to store posterior draws of nu // each row is rho_1,...,rho_m
  arma::vec sig2invkeep(n_outer,arma::fill::zeros);           // vector to store posterior draws of sigma^{-2}
  arma::mat gammakeep(n_outer,P_z,arma::fill::zeros);         // matrix to store posterior draws of gamma // each row is gamma_1,...,gamma_{P_z}
  arma::field<arma::mat>  thetakeepfield(M);                  // field of matrices to store posterior draws of theta (for computation)
  arma::field<arma::mat>  nukeepfield(M);                     // field of matrices to store posterior draws of nu (for computation)
  Rcpp::List  thetakeep(M);                                   // List of matrices to store posteror draws of theta (same as above but for exporting to R)
  Rcpp::List  nukeep(M);                                      // List of matrices to store posteror draws of nu (same as above but for exporting to R)
  for(int m=0; m<M; m++){
    arma::mat thetakeeptemp(n_outer,Lvec(m),arma::fill::zeros); // temporary matrix with L_m columns to store posteror draws of theta_m
    thetakeepfield[m] = thetakeeptemp;                          // initialize elements of thetakeepfield to be the right size
    thetakeep[m] = thetakeeptemp;                               // initialize elements of thetakeepfield to be the right size
    
    nukeepfield[m] = thetakeeptemp;                             // initialize elements of nukeepfield to be the right size
    nukeep[m] = thetakeeptemp;                                  // initialize elements of nukeepfield to be the right size
  } 
  
  // Initialize values for posterior of h(.)
  arma::mat laminvK(N,N,arma::fill::zeros);           // lambda^{-1}K where K is the kernel matrix
  arma::mat laminvKSiginv(N,N,arma::fill::zeros);     // \lambda^{-1}KSigma^{-1}=laminvK*Sigma^{-1} or t(solve(Sigma)*laminvK)
  arma::vec hmean(N,arma::fill::zeros);               // marginal mean vector for h(.)
  arma::mat hcov(N,N,arma::fill::zeros);              // marginal variance-covariance matrix of h(.)
  arma::vec hmeantemp(N,arma::fill::zeros);           // temporary vector for computing marginal covariance of h(.)
  arma::mat hcovtemp(N,N,arma::fill::zeros);          // marginal variance-covariance matrix of h(.)
  arma::mat SIG(P_z,P_z,arma::fill::zeros);           // temporary matrix
  arma::mat hsampkeep(n_outer,N,arma::fill::zeros);   // matrix to store posterior draws of h: rows are samples of h
  arma::mat svdU;                                     // for SVD decomposition
  arma::vec svdD;                                     // for SVD decomposition
  arma::mat svdV;                                     // for SVD decomposition
  
  // -------------------------------------
  // MCMC
  // -------------------------------------
  for(int s_outer=0; s_outer<n_outer; s_outer++){                       // outer loop
    
    if(s_outer % 10 == 0) Rcout << "#" << s_outer << std::endl;         // output progress every 10 iterations
    
    for(int s_inner=0; s_inner<n_inner; s_inner++){                     // inner loop
      
      
      // -------------------------------------------------------------------------- //
      // UPDATE LAMBDA^{-1} via random walk metropolis
      logLambdaInverse_PROP = logLambdaInverse + rnorm(1)(0); // random walk step
      
      // update Sigma=I+\lambda^{-1}*K, A1, A2, C1, C2
      Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse_PROP, XthetaStar, randint, logLambdaBInverse, Bmat);   
      
      // proposed ll 
      ll_PROP = -logLambdaInverse_PROP*logLambdaInverse_PROP/(2*b_lambda)
        - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
        - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
        
      // log ratio for MH acceptance
      logratio = ll_PROP-ll;
      if(logratio>0){
        logratio = 0;
      }
      
      // accept or reject    
      if(logratio >  log(runif(1)(0)) ){                
        logLambdaInverse = logLambdaInverse_PROP; 
        ll = ll_PROP;
      }

      // Reset matrix Sigma and ll
      Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat); 

      ll = -logLambdaInverse*logLambdaInverse/(2*b_lambda) 
        - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
        - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
        
        
        
        
      // -------------------------------------------------------------------------- //
      // UPDATE LAMBDA_B^{-1} via random walk metropolis
      if(randint==TRUE){ // only update lambda_B^{-1} if using random intercepts model
        
        logLambdaBInverse_PROP = logLambdaBInverse + rnorm(1)(0); // random walk step
        
        // current llB (uses the prior for logLambdaBInverse, hence the first line below is different than ll elsewhere)
        llB = -logLambdaBInverse*logLambdaBInverse/(2*b_lambdaB) 
          - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
          - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 

        // update Sigma=I+\lambda^{-1}*K, A1, A2, C1, C2
        Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse_PROP, Bmat);   

        // proposed llB (uses the prior for logLambdaBInverse, hence the first line below is different than ll elsewhere)
        llB_PROP = -logLambdaBInverse_PROP*logLambdaBInverse_PROP/(2*b_lambdaB)
          - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
          - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
          
        // log ratio for MH acceptance
        logratio = llB_PROP-llB;
        if(logratio>0){
          logratio = 0;
        }
        
        // accept or reject    
        if(logratio >  log(runif(1)(0)) ){                
          logLambdaBInverse = logLambdaBInverse_PROP; 
        }
          
        // Reset matrix Sigma and ll
        Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat); 
        
        // resetting ll (not llB)
        ll = -logLambdaInverse*logLambdaInverse/(2*b_lambda) 
          - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
          - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
          
      }
      
          
        
      // -------------------------------------------------------------------------- //
      // UPDATE TAU via random walk metropolis
      if(horseshoe==TRUE){  // only update tau if using horseshoe priors
        
        logtau_PROP = logtau + step_tau*rnorm(1)(0); // random walk step
        
        // temporary sum of (theta*/nu)^2 for lltau
        temp_sumthet = 0;                
        for(int m=0; m<M; m++){    // THERE WAS A BUG HERE: M-1
          temp_sumthet += sum(square(thetaStar[m]/exp(lognu[m])));
        }
        
        // current lltau (updating to include correct thetaStar)
        lltau =  (1-sum(Lvec))*logtau
          - 0.5*exp(-2*logtau)*temp_sumthet
          - log(tau02) - log(1+(exp(2*logtau)/(tau02*tau02))); 
          
        // proposed lltau
        lltau_PROP =  (1-sum(Lvec))*logtau_PROP
          - 0.5*exp(-2*logtau_PROP)*temp_sumthet
          - log(tau02) - log(1+(exp(2*logtau_PROP)/(tau02*tau02))); 
          
        // log ratio for MH acceptance
        logratio = lltau_PROP-lltau;
        if(logratio>0){
          logratio = 0;
        }
        
        // accept or reject    
        if(logratio >  log(runif(1)(0)) ){                
          logtau = logtau_PROP; 
        }
        
      }
      
       
      
      // -------------------------------------------------------------------------- //
      // UPDATE NU 
      if(horseshoe==TRUE){ // if using horseshoe priors, update via random walk metropolis
        
        for(int m=0; m<M; m++){
          
          for(int l=0; l<Lvec(m); l++){
            lognu_PROP = lognu[m](l) + step_nu*rnorm(1)(0); // random walk step
            
            // current llnu (updating to include correct thetastar)
            llnu =  - 0.5*exp(-2*logtau-2*lognu[m](l))*pow(thetaStar[m](l),2.0)
              - log( 1 + exp(2*lognu[m](l)) );
            
            // proposed llnu
            llnu_PROP =  - 0.5*exp(-2*logtau-2*lognu_PROP)*pow(thetaStar[m](l),2.0)
              - log( 1 + exp(2*lognu_PROP) );
            
            // log ratio for MH acceptance
            logratio = llnu_PROP-llnu;
            if(logratio>0){
              logratio = 0;
            }
            
            // accept or reject    
            if(logratio >  log(runif(1)(0)) ){                
              lognu[m](l) = lognu_PROP; 
            }
            
          }       // l loop
          
        }       // m loop
        
      }else{ // if not using horseshoe, draw directly from generalized inverse gaussian
        
        for(int m=0; m<M; m++){ // note there are only m nus under the non-horseshoe setup (excess columns are deleted in the R wrapper function)
          
          lognu[m](0) = log( do_rgiglocal(-(Lvec(m)-1)/2, sum(square(thetaStar[m]))/kappa[m], 1) );    // draw directly from generalized inverse gaussian
          
          while(lognu[m](0)==R_PosInf){ // if its being sent to Inf because sum(square(thetaStar[m])) is too close to zero
            lognu[m](0) = log( do_rgiglocal(-(Lvec(m)-1)/2, pow(10,-14.0)/kappa[m], 1) ); // set it just close enough to zero
          }
          // nu[m] = do_rgiglocal(-(Lvec(m)-1)/2, sum(square(thetaStar[m]))/kappa[m], 1);    // draw directly from generalized inverse gaussian
          
        }       // m loop
        
      } // end IF statement for NU
      
          
          
             
              
              
      // -------------------------------------------------------------------------- //
      // UPDATE thetaStar via elliptical slice sampler 
      for(int m=0; m<M; m++){   // loop over M indices
        
        if(horseshoe==TRUE){      // under horseshoe prior
          
          vec_PROP = as<arma::vec>(rnorm(Lvec(m)))%exp(lognu[m])*exp(logtau);   // draw proposed V to define ellipse
        
        }else{  // under prior with no variable selection
          
          vec_PROP = rnorm(Lvec(m))*sqrt(exp(lognu[m](0)))*sqrt(kappa[m]);   // draw proposed V to define ellipse
          
        }
        
        // draw random angle phi from (0,pi)
        angle = runif(1)(0)*M_PI;                                             
        amax = angle;                           // set max angle phi_min=phi
        amin = angle - M_PI;                    // set min angle phi_max=phi-pi

        llu = log(runif(1)(0)) + ll;            // threshold = ll+log(U)
        ll_PROP = llu-1;                        // initialize ll_PROP
        
        // shrink slice
        while(ll_PROP < llu){   // accept when ll_PROP>ll+log(U)  
          
          // propose thetaStar
          thetaStar_PROP = thetaStar[m]*cos(angle) + vec_PROP*sin(angle); 
          
          // update components that depend on thetaStar
          XthetaStar.col(m) = X[m] * thetaStar_PROP;                          // update XthetaStar with proposed thetaStar

          // update Sigma with proposed thetaStar
          Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat); 

          // proposed ll 
          ll_PROP = -logLambdaInverse*logLambdaInverse/(2*b_lambda) 
            - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
            - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
            
          if(angle > 0){
            amax = angle;
          }else{
            amin = angle;
          }
          
          angle = runif(1)(0)*(amax-amin)+amin;   // shrink slice
          
        } // end loop shrinking elliptical slice
        
        // save results
        ll =  ll_PROP;
        thetaStar[m] = thetaStar_PROP;
        XthetaStar.col(m) = X[m] * thetaStar_PROP;                      
        
      } // end loop through update of thetaStar
      
      
      
      
    } // end inner loop
    
    
    
    
    // save posterior samples
    for(int m=0; m<M; m++){
      // partition thetaStar* into rho and theta
      rhokeep.submat(s_outer,m,s_outer,m) = sum(square(vectorise(thetaStar[m])));                             // rho=||theta*||^2
      thetakeepfield[m].row(s_outer) =  thetaStar[m].t()/sqrt(as_scalar(rhokeep.submat(s_outer,m,s_outer,m)));// theta=||theta||^{-1}theta*
      
      nukeepfield[m].row(s_outer) = exp(lognu[m].t());            // save nu
      // nukeep.submat(s_outer,m,s_outer,m) = exp(lognu[m]); // save nu
    }
    logLambdaInversekeep(s_outer) = logLambdaInverse;       // save lambda^{-1}  
    logLambdaBInversekeep(s_outer) = logLambdaBInverse;     // save lambda^{-1}  
    logtaukeep(s_outer) = logtau;                           // save tau
    
    // draw sigma^{-2} and gamma directly and save them
    sig2invkeep(s_outer) = rgamma( 1 , a_sig + .5*(N-P_z) , 1/(b_sig + 0.5*(Vcomps[1](0,0)-Vcomps[1](1,0))) )(0);  // function generates a gamma distribution parameterized with mean alpha*beta so we just inverted the second component;
    // gammakeep.row(s_outer) = arma::solve(trimatu(C2) , C2B +  as<arma::vec>(rnorm(P_z))/sqrt(sig2invkeep(s_outer) )).t(); // draw gamma by drawing MVN(0,I) then scaling by sqrt of variance and shifting by mean
    SIG = inv(Vcomps[4]*Vcomps[4].t());

    gammakeep.row(s_outer) = (SIG*Vcomps[2]+chol(SIG,"lower")*as<arma::vec>(rnorm(P_z))/sqrt(sig2invkeep(s_outer) )).t();

    // compute fitted h after burn-in
    if(s_outer>=n_burn){
      laminvK = Vcomps[5];      
      laminvKSiginv = solve(Vcomps[0] , laminvK).t();                                     // lambda^{-1} K Sigma^{-1} which is equal to t(Sigma^{-1} lambda^{-1} K) 
      hmeantemp = laminvKSiginv * (yz.col(0)-yz.cols(1,P_z)*gammakeep.row(s_outer).t());  // temporary vector to be used in mean and covariance
      // hcovtemp = laminvKSiginv*(Vcomps[0]-laminvK)/sig2invkeep(s_outer);               // temporary covariance matrix to be used in marginal covariance
      hcovtemp = (laminvK-laminvKSiginv*laminvK.t())/sig2invkeep(s_outer);                // temporary covariance matrix to be used in marginal covariance
      // compute marginal mean and covariance for h
      hmean += hmeantemp/(n_outer-n_burn);                                                                    // law of total expecation, E[E[h|s]]
      hcov += hcovtemp/(n_outer-n_burn) + hmeantemp*hmeantemp.t()/(n_outer-n_burn); // first two components of law of total variance: E(Var(h|s))+E[E[h|s]^2]-E[E[h|s]]^2
          // old (no B) : hcov += laminvKSiginv/(n_outer-n_burn)/sig2invkeep(s_outer) + hmeantemp*hmeantemp.t()/(n_outer-n_burn); // first two components of law of total variance: E(Var(h|s))+E[E[h|s]^2]-E[E[h|s]]^2
          // hcov += laminvKSiginv*(Vcomps[0]-laminvK)/(n_outer-n_burn)/sig2invkeep(s_outer) + hmeantemp*hmeantemp.t()/(n_outer-n_burn); // first two components of law of total variance: E(Var(h|s))+E[E[h|s]^2]-E[E[h|s]]^2
      if(draw_h==TRUE){
        arma::svd(svdU,svdD,svdV,hcovtemp);                                                 // SVD decomposition for covariance of h (since it's more stable than cholesky)
        // draw h
        hsampkeep.row(s_outer) = ( hmeantemp+(svdU*diagmat(svdD))*arma::randn(N) ).t();     // hsampkeep.row(s_outer) = ( hmeantemp+arma::chol(hcovtemp)*arma::randn(N) ).t();
      }
    }

  } // end outer loop
  
  
  
  
  // Convert field to list
  for(int m=0; m<M; m++){
    thetakeep[m] = thetakeepfield[m];
    nukeep[m] = nukeepfield[m];
  }
  
  // return a list
  return List::create(Named("theta")          = thetakeep,
                      Named("lambdaInverse")  = exp(logLambdaInversekeep),
                      Named("lambdaBInverse") = exp(logLambdaBInversekeep),
                      Named("rho")            = rhokeep,
                      Named("nu")             = nukeep,
                      Named("tau")            = exp(logtaukeep),
                      Named("sigma2")         = 1/sig2invkeep,
                      Named("gamma")          = gammakeep,
                      Named("hsamp")          = hsampkeep,
                      Named("hmean")          = hmean,
                      Named("hcov")           = hcov - ( hmean*hmean.t())); // subtracting third component of law of total variance (i.e. E[E[h|s]]^2)
}









// [[Rcpp::depends(RcppArmadillo)]]
//' MCMC for BSMIM with spike and slab prior (inv-uniform slab for rho version)
//' 
//' Returns posterior draws for \eqn{\theta},\eqn{\lambda^{-1}},\eqn{\rho},\eqn{\sigma^2},\eqn{\gamma}, as well as the mean(\eqn{h}) and cov(\eqn{h})
//' Adapted from "bkmrdlm_multi.cpp" in the "regimes" package by Ander Wilson.
//' Note that this implementation with horseshoe2 prior returns posterior draws of \eqn{\rho} as a list of M matrices (unlike main MCMC code, which returns one matrix of M columns)
//' 
//' @param yz matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
//' @param Xlist list of (N by L_m) matrices representing X_m
//' @param a_lam shape hyperparameter for \eqn{\lambda^{-1}}
//' @param b_lam rate hyperparameter for \eqn{\lambda^{-1}}
//' @param b_lambdaB hyperparameter for \eqn{\lambda^{-1}_B}
//' @param a_sig first hyperparameter for \eqn{\sigma^{-2}}
//' @param b_sig second hyperparameter for \eqn{\sigma^{-2}}
//' @param a_theta lower bound of unif for invunif prior on \eqn{\theta^{*2}_m}
//' @param b_theta upper bound of unif for invunif prior on \eqn{\theta^{*2}_m}
//' @param poly 0 = gaussian kernel, 1 = polynomial kernel
//' @param d degree of polynomial kernel
//' @param horseshoe 0 = no selection, 1 = componentwise horseshoe priors
//' @param randint 0 = no random intercepts, 1 = random intercepts model
//' @param Bmat N xN block diagonal matrix indicating cluster membership for random intercepts model
//' @param draw_h 0 = dont draw h, 1 = draw h
//' @param n_inner no. of MCMC iterations to run in the inner loop. n_outer*n_inner iteraction will be run.
//' @param n_outer no. of MCMC iterations to run in the outer loop. n_outer iterations will be saved.
//' @param n_burn no. of MCMC iterations to discard as burn-in
//' @author Glen McGee and Ander Wilson (adapted from the "regimes" package by Ander Wilson).
//' @export
// [[Rcpp::export]]
List bsmim_spikeslab_mcmc2(const arma::mat&    yz,         // matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
                           const Rcpp::List&   Xlist,      // list of (N by L_m) matrices representing X_m
                           const double&       a_lam,   // shape hyperparameter for \lambda^{-1}
                           const double&       b_lam,    // rate hyperparameter for \lambda^{-1}
                           const double&       b_lambdaB,  // hyperparameter for \lambda^{-1}_B
                           const double&       a_sig,      // first hyperparameter for \sigma^{-2}
                           const double&       b_sig,      // second hyperparameter for \sigma^{-2}
                           const double&       a_theta,    // hyperparameter for unif of invunif prior for thetastar^2 (slab component)
                           const double&       b_theta,    // hyperparameter for unif of invunif prior for thetastar^2 (slab component)
                           const double&       step_theta, // step size for theta* random walk (move 2)
                           const double&       a_pi,       // first hyperparameter for beta distribution of pi 
                           const double&       b_pi,       // second hyperparameter for beta distribution of pi 
                           const bool&         poly,       // 0=gaussian kernel / 1=polynomial kernel
                           const int&          d,          // degree of polynomial kernel
                           const bool&         randint,    // 0=no random intercepts / 1=random intercepts
                           const arma::mat&    Bmat,       // block diagonal matrix B indicating cluster membership for random intercepts model
                           const bool&         draw_h,     // 0=dont sample h, 1= sample h
                           const int&          n_inner,    // no. of MCMC iterations to run in the inner loop. n_outer*n_inner iteraction will be run.
                           const int&          n_outer,    // no. of MCMC iterations to run in the outer loop. n_outer iterations will be saved.
                           const int&          n_burn) {   // no. of MCMC iterations to discard as burn-in
  
  // step size for random walk in move 2 of MH step (spike/slab)
  // double step_theta = 0.2;
  // // TRYING DIFFERENT LAMBDA PRIOR
  // double a_lam = 1.0;
  // double b_lam = 0.1;
  
  // inf-unif prior for rho
  double ar = a_theta;
  double br = b_theta;
  
  
  
  // Dimensions
  int N = yz.n_rows;              // no. of observations
  int M = Xlist.size();           // no. of indices (dimension of h(.))
  int P_z = yz.n_cols-1;          // no. of covariates
  IntegerVector Lvec = rep(0,M);  // vector of L_m representing the no. of columns of X_m
  
  // Initialize
  arma::field<arma::mat> X(M);                    // field of (N by L_m) matrices representing X_m
  arma::field<arma::vec> thetaStar(M);            // field of L_m vectors representing \thetaStar_m
  arma::field<arma::vec> delta(M);                // field of L_m vectors representing \delta_m
  arma::mat  XthetaStar(N,M,arma::fill::zeros);   // (N by M) matrix representing [X_1\thetaStar_1,...,X_M\thetaStar_M]
  for(int m=0; m<M; m++){
    X[m] = as<arma::mat>(Xlist[m]);                 // mth matrix in X is X_m
    Lvec(m) = X[m].n_cols;                          // mth element of Lvec is L_m (the length of \thetaStar_m)
    thetaStar[m] = as<arma::vec>(rnorm(Lvec(m)));   // set starting value for \thetaStar_m
    XthetaStar.col(m) = X[m] * thetaStar[m];
    
    arma::vec deltatemp(Lvec(m),arma::fill::ones); // temporary vec with L_m elements to store posteror draws of delta
    delta[m] = deltatemp;                          // set starting value for \delta_m
  }
  // starting vals must respect bounds of prior
  for(int m=0; m<M; m++){
    for(int l=0; l<Lvec(m); l++){
      // Rcout << thetaStar[m](l) << std::endl;
      // Rcout << pow(br,-0.5) << std::endl;
      if(thetaStar[m](l)>-pow(br,-0.5) & thetaStar[m](l)<pow(br,-0.5)){
        // Rcout << "m " << m << "l " << l << std::endl;      
        thetaStar[m](l) = rinttruncnorm(0,step_theta,-pow(br,-0.5),pow(br,-0.5));
        // Rcout << thetaStar[m](l) << std::endl;
      }
    }
  }
  
  
  double logLambdaInverse = rnorm(1)(0)*sqrt(1)/10;    // set starting value for log(\lambda^{-1})
  double logLambdaBInverse = 0;                               // set log(\lambda^{-1}_B) to 0 if no random intercepts
  if(randint==TRUE){
    logLambdaBInverse = rnorm(1)(0)*sqrt(b_lambdaB)/10;       // set starting value for log(\lambda^{-1}_B) if using random intercepts
  }
  
  
  
  // Construct Sigma=I+\lambda^{-1}*K, A1, A2, C1, C2
  arma::field<arma::mat> Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat);
  
  // First evaluation of the log-liklihood (ll)
  double ll = logLambdaInverse+R::dgamma(exp(logLambdaInverse),a_lam,1/b_lam,1) //-logLambdaInverse*logLambdaInverse/(2*b_lambda) 
    - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
    - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
    
    
    // Initialize other values
    double logLambdaInverse_PROP, llu, ll_PROP, logratio;
    double llB, logLambdaBInverse_PROP, llB_PROP;
    arma::vec ons(N, arma::fill::ones);           // vector of length N of ones
    arma::vec vec_PROP;                           // proposal vector to define ellipse
    
    int MHmove, sumdelta, sumdelta_PROP, move_which, move_id;    // move type for joint MH draw of delta and theta* (1 or 2); sum of delta=1; which component MHmove applies to; id to sum
    double lltheta, lltheta_PROP, logdiff_PROP;
    arma::field<arma::vec> thetaStar_PROP = thetaStar;
    arma::field<arma::vec> delta_PROP = delta;
    
    
    // Store posterior samples
    arma::vec logLambdaInversekeep(n_outer,arma::fill::zeros);  // vector to store posterior draws of lambda^{-1}
    arma::vec logLambdaBInversekeep(n_outer,arma::fill::zeros); // vector to store posterior draws of lambda_B^{-1}
    arma::mat rhokeep(n_outer,M,arma::fill::zeros);             // matrix to store posterior draws of rho // each row is rho_1,...,rho_m
    arma::vec sig2invkeep(n_outer,arma::fill::zeros);           // vector to store posterior draws of sigma^{-2}
    arma::mat gammakeep(n_outer,P_z,arma::fill::zeros);         // matrix to store posterior draws of gamma // each row is gamma_1,...,gamma_{P_z}
    arma::field<arma::mat>  thetakeepfield(M);                  // field of matrices to store posterior draws of theta (for computation)
    Rcpp::List  thetakeep(M);                                   // List of matrices to store posteror draws of theta (same as above but for exporting to R)
    for(int m=0; m<M; m++){
      arma::mat thetakeeptemp(n_outer,Lvec(m),arma::fill::zeros); // temporary matrix with L_m columns to store posteror draws of theta_m
      thetakeepfield[m] = thetakeeptemp;                          // initialize elements of thetakeepfield to be the right size
      thetakeep[m] = thetakeeptemp;                               // initialize elements of thetakeepfield to be the right size
      
    } 
    
    // Initialize values for posterior of h(.)
    arma::mat laminvK(N,N,arma::fill::zeros);           // lambda^{-1}K where K is the kernel matrix
    arma::mat laminvKSiginv(N,N,arma::fill::zeros);     // \lambda^{-1}KSigma^{-1}=laminvK*Sigma^{-1} or t(solve(Sigma)*laminvK)
    arma::vec hmean(N,arma::fill::zeros);               // marginal mean vector for h(.)
    arma::mat hcov(N,N,arma::fill::zeros);              // marginal variance-covariance matrix of h(.)
    arma::vec hmeantemp(N,arma::fill::zeros);           // temporary vector for computing marginal covariance of h(.)
    arma::mat hcovtemp(N,N,arma::fill::zeros);          // marginal variance-covariance matrix of h(.)
    arma::mat SIG(P_z,P_z,arma::fill::zeros);           // temporary matrix
    arma::mat hsampkeep(n_outer,N,arma::fill::zeros);   // matrix to store posterior draws of h: rows are samples of h
    arma::mat svdU;                                     // for SVD decomposition
    arma::vec svdD;                                     // for SVD decomposition
    arma::mat svdV;                                     // for SVD decomposition
    
    // -------------------------------------
    // MCMC
    // -------------------------------------
    for(int s_outer=0; s_outer<n_outer; s_outer++){                       // outer loop
      
      if(s_outer % 10 == 0) Rcout << "#" << s_outer << std::endl;         // output progress every 10 iterations
      
      for(int s_inner=0; s_inner<n_inner; s_inner++){                     // inner loop
        
        
        // -------------------------------------------------------------------------- //
        // UPDATE LAMBDA^{-1} via random walk metropolis
        logLambdaInverse_PROP = logLambdaInverse + rnorm(1)(0); // random walk step
        
        // update Sigma=I+\lambda^{-1}*K, A1, A2, C1, C2
        Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse_PROP, XthetaStar, randint, logLambdaBInverse, Bmat);   
        
        // proposed ll  (inverting b_lam, since this takes scale instead of rate)
        ll_PROP = logLambdaInverse_PROP+R::dgamma(exp(logLambdaInverse_PROP),a_lam,1/b_lam,1) //-logLambdaInverse*logLambdaInverse/(2*b_lambda) -logLambdaInverse_PROP*logLambdaInverse_PROP/(2*b_lambda)
          - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
          - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
          
          // log ratio for MH acceptance
          logratio = ll_PROP-ll;
          if(logratio>0){
            logratio = 0;
          }
          
          // accept or reject    
          if(logratio >  log(runif(1)(0)) ){                
            logLambdaInverse = logLambdaInverse_PROP; 
            ll = ll_PROP;
          }
          
          // Reset matrix Sigma and ll
          Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat); 
          
          ll = logLambdaInverse+R::dgamma(exp(logLambdaInverse),a_lam,1/b_lam,1) //-logLambdaInverse*logLambdaInverse/(2*b_lambda)
            - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
            - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
            
            
            
            
            // -------------------------------------------------------------------------- //
            // UPDATE LAMBDA_B^{-1} via random walk metropolis
            if(randint==TRUE){ // only update lambda_B^{-1} if using random intercepts model
              
              logLambdaBInverse_PROP = logLambdaBInverse + rnorm(1)(0); // random walk step
              
              // current llB (uses the prior for logLambdaBInverse, hence the first line below is different than ll elsewhere)
              llB = logLambdaBInverse+R::dgamma(exp(logLambdaBInverse),a_lam,1/b_lam,1) //-logLambdaBInverse*logLambdaBInverse/(2*b_lambdaB) 
                - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
                
                // update Sigma=I+\lambda^{-1}*K, A1, A2, C1, C2
                Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse_PROP, Bmat);   
                
                // proposed llB (uses the prior for logLambdaBInverse, hence the first line below is different than ll elsewhere)
                llB_PROP = logLambdaBInverse_PROP+R::dgamma(exp(logLambdaBInverse_PROP),a_lam,1/b_lam,1) //-logLambdaBInverse_PROP*logLambdaBInverse_PROP/(2*b_lambdaB)
                  - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                  - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
                  
                  // log ratio for MH acceptance
                  logratio = llB_PROP-llB;
                  if(logratio>0){
                    logratio = 0;
                  }
                  
                  // accept or reject    
                  if(logratio >  log(runif(1)(0)) ){                
                    logLambdaBInverse = logLambdaBInverse_PROP; 
                  }
                  
                  // Reset matrix Sigma and ll
                  Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat); 
                  
                  // resetting ll (not llB)
                  ll = logLambdaInverse+R::dgamma(exp(logLambdaInverse),a_lam,1/b_lam,1) //-logLambdaInverse*logLambdaInverse/(2*b_lambda) 
                    - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                    - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
                    
            }
            
            
            
            // -------------------------------------------------------------------------- //
            // UPDATE delta and thetaStar via Metropolis Hastings
            
            // count number of delta==1
            sumdelta = 0;
            for(int m=0; m<M; m++){
              sumdelta += sum(delta[m]);  
            }
            
            // select move type
            if(sumdelta==0){ // must choose move 1 if all delta==0 (move 2 samples from the deltas equal to 1)
              MHmove = 1;
            }else{ // randomly select move type
              MHmove = Rcpp::rbinom(1,1,0.5)(0);
            }
            
            // Make move 1 or 2
            if(MHmove==1){  // MHmove version 1
              
              // reset proposals
              delta_PROP = delta;
              thetaStar_PROP = thetaStar;
              
              // randomly select a component
              move_which = sample(sum(Lvec), 1)(0);
              
              move_id = 0;                                // set index
              for(int m=0; m<M; m++){
                for(int l=0; l<Lvec(m); l++){
                  
                  move_id += 1;
                  if(move_id==move_which){   // make move for randomly selected component
                    
                    
                    // Rcout << "MOVE 1: m-" << m << " l-" << l << std::endl;
                    
                    // compute lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
                    lltheta = log(tgamma(sumdelta+a_pi))+log(tgamma(sum(Lvec)-sumdelta+b_pi))
                    // + delta[m](l)*(-3*log( fabs(thetaStar[m](l))) - log( 1*( (pow(br,-0.5)<fabs(thetaStar[m](l))) & (fabs(thetaStar[m](l))<pow(ar,-0.5))  ) ))  //R::dnorm(thetaStar_PROP[m](l),0.0,s_theta,TRUE) //-logLambdaInverse*logLambdaInverse/(2*b_lambda)R::dnorm(thetaStar[m](l),0.0,s_theta,TRUE)//+(1-delta[m][l])*log(1)//-logLambdaInverse*logLambdaInverse/(2*b_lambda) 
                      - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                      - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
                    if(delta[m](l)==1){//adding component to lltheta (avoiding NAs)
                      lltheta += delta[m](l)*(-3*log( fabs(thetaStar[m](l))) - log( 1*( (pow(br,-0.5)<fabs(thetaStar[m](l))) & (fabs(thetaStar[m](l))<pow(ar,-0.5))  ) ));
                    }  
                      // propose new delta: 1 to 0 or 0 to 1
                      delta_PROP[m](l) = 1-delta[m](l);
                      if(delta_PROP[m](l)==1){ // adjust sum delta for switch (either increase by 1 or decrease by 1)
                        sumdelta_PROP = sumdelta+1;
                      }else{
                        sumdelta_PROP = sumdelta-1;
                      }
                      
                      //propose new theta* (based on new delta)
                      if(delta_PROP[m](l)==0){    // if delta_PROP, then theta*=0
                        thetaStar_PROP[m](l) = 0;
                      }else{                      // otherwise draw from Q1 (which is the prior)
                        thetaStar_PROP[m](l) =  1/sqrt(runif(1,ar,br)(0));// thetaStar_PROP[m](l) = s_theta*rnorm(1)(0);
                        if(runif(1)(0)<0.5){ // flip coin for pos/neg
                          thetaStar_PROP[m](l) = -thetaStar_PROP[m](l);
                        }
                      }
                      
                      // update components that depend on thetaStar
                      XthetaStar.col(m) = X[m] * thetaStar_PROP[m];                          // update XthetaStar with proposed thetaStar
                      
                      // update Sigma with proposed thetaStar
                      Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat); 
                      
                      // proposed lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
                      lltheta_PROP = log(tgamma(sumdelta_PROP+a_pi))+log(tgamma(sum(Lvec)-sumdelta_PROP+b_pi))
                        // + delta_PROP[m](l)*(-3*log( fabs(thetaStar_PROP[m](l))) - log( 1*( (pow(br,-0.5)<fabs(thetaStar_PROP[m](l))) & (fabs(thetaStar_PROP[m](l))<pow(ar,-0.5))  ) ))  //R::dnorm(thetaStar_PROP[m](l),0.0,s_theta,TRUE) //-logLambdaInverse*logLambdaInverse/(2*b_lambda)//+(1-delta[m][l])*log(1)//-logLambdaInverse*logLambdaInverse/(2*b_lambda) 
                          - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                          - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
                      if(delta_PROP[m](l)==1){//adding component to lltheta_PROP (avoiding NAs)
                        lltheta_PROP += delta_PROP[m](l)*(-3*log( fabs(thetaStar_PROP[m](l))) - log( 1*( (pow(br,-0.5)<fabs(thetaStar_PROP[m](l))) & (fabs(thetaStar_PROP[m](l))<pow(ar,-0.5))  ) ));
                      }    
                          
                          
                      // difference between proposal distributions: logP(theta_prop,delta_prop|theta,delta)-logP(theta,delta|theta_prop,delta_prop)
                      logdiff_PROP = //delta_PROP[m](l)*(-3*log( fabs(thetaStar_PROP[m](l))) - log( 1*( (pow(br,-0.5)<fabs(thetaStar_PROP[m](l))) & (fabs(thetaStar_PROP[m](l))<pow(ar,-0.5))  ) ))  //R::dnorm(thetaStar_PROP[m](l),0.0,s_theta,TRUE) //   logP(theta_prop,delta_prop|theta,delta)
                        //- delta[m](l)*(-3*log( fabs(thetaStar[m](l))) - log( 1*( (pow(br,-0.5)<fabs(thetaStar[m](l))) & (fabs(thetaStar[m](l))<pow(ar,-0.5))  ) ))  //R::dnorm(thetaStar[m](l),0.0,s_theta,TRUE) // logP(theta,delta|theta_prop,delta_prop)
                          +((sumdelta_PROP==0)*log(1)+(sumdelta_PROP!=0)*log(0.5)) // extra difference component because sometimes we must choose move 1.
                          -((sumdelta==0)*log(1)+(sumdelta!=0)*log(0.5));
                    
                      if(delta_PROP[m](l)==1){ //adding component to logdiff (avoiding NAs)
                        logdiff_PROP += delta_PROP[m](l)*(-3*log( fabs(thetaStar_PROP[m](l))) - log( 1*( (pow(br,-0.5)<fabs(thetaStar_PROP[m](l))) & (fabs(thetaStar_PROP[m](l))<pow(ar,-0.5))  ) ));
                      }else{
                        logdiff_PROP -= delta[m](l)*(-3*log( fabs(thetaStar[m](l))) - log( 1*( (pow(br,-0.5)<fabs(thetaStar[m](l))) & (fabs(thetaStar[m](l))<pow(ar,-0.5))  ) ));
                      }  
                          
                      // log ratio for MH acceptance
                      logratio = (lltheta_PROP-lltheta)-logdiff_PROP; // the extra (negative!) difference in log proposals since this is MH not metropolis
                      if(logratio>0){
                        logratio = 0;
                      }
                              
                              
                      // accept or reject    
                      if(logratio >  log(runif(1)(0)) ){    
                        delta[m](l) = delta_PROP[m](l);
                        thetaStar[m](l) = thetaStar_PROP[m](l);
                        // ll = ll_PROP;
                      }else{
                        XthetaStar.col(m) = X[m] * thetaStar[m];
                      }
                              
                              
                              
                  } // end (if move_id)
                  
                  
                  
                } // l loop
              } // m loop
              
              
              
              
              
            }else{ // MHmove version 2
              
              move_which = sample(sumdelta, 1)(0);    // randomly select component for move (such that delta==1)
              
              move_id = 0;                                // set index
              for(int m=0; m<M; m++){
                for(int l=0; l<Lvec(m); l++){
                  
                  if(delta[m](l)==1){
                    move_id += 1;
                    if(move_id==move_which){   //RANDOM WALK STEP
                      
                      // Rcout << "MOVE 2: m-" << m << " l-" << l << std::endl;
                      
                      // compute lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
                      lltheta = -3*log( fabs(thetaStar[m](l))) - log( 1*( (pow(br,-0.5)<fabs(thetaStar[m](l))) & (fabs(thetaStar[m](l))<pow(ar,-0.5))  ) )  //R::dnorm(thetaStar_PROP[m](l),0.0,s_theta,TRUE) //-logLambdaInverse*logLambdaInverse/(2*b_lambda)
                      - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                      - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
                      
                      // propose new theta
                      thetaStar_PROP[m](l) = rinttruncnorm(thetaStar[m](l),step_theta,-pow(br,-0.5),pow(br,-0.5));  // thetaStar_PROP[m](l) = thetaStar[m](l) + step_theta*rnorm(1)(0);
                      
                      // update components that depend on thetaStar
                      XthetaStar.col(m) = X[m] * thetaStar_PROP[m];                          // update XthetaStar with proposed thetaStar
                      
                      // update Sigma with proposed thetaStar
                      Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat); 
                      
                      // proposed lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
                      lltheta_PROP = -3*log( fabs(thetaStar_PROP[m](l))) - log( 1*( (pow(br,-0.5)<fabs(thetaStar_PROP[m](l))) & (fabs(thetaStar_PROP[m](l))<pow(ar,-0.5))  ) )  //R::dnorm(thetaStar_PROP[m](l),0.0,s_theta,TRUE) //-logLambdaInverse*logLambdaInverse/(2*b_lambda) 
                        - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                        - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
                        
                      logdiff_PROP = logdinttruncnorm(thetaStar_PROP[m](l),thetaStar[m](l),step_theta,-pow(br,-0.5),pow(br,-0.5))
                        -logdinttruncnorm(thetaStar[m](l),thetaStar_PROP[m](l),step_theta,-pow(br,-0.5),pow(br,-0.5));
                        
                      // log ratio for MH acceptance
                      logratio = (lltheta_PROP-lltheta)-logdiff_PROP;
                      if(logratio>0){
                        logratio = 0;
                      }
                        
                      // accept or reject    
                      if(logratio >  log(runif(1)(0)) ){       
                        thetaStar[m](l) = thetaStar_PROP[m](l);
                        // ll = ll_PROP;
                      }else{
                        XthetaStar.col(m) = X[m] * thetaStar[m];
                      }
                        
                        
                        
                    } // end (if move_id)
                  } // end (if delta=1)
                  
                  
                  
                } // l loop
              } // m loop
            } // end move 2
            
            // Reset matrix Sigma and ll
            Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat); 
            
            ll = logLambdaInverse+R::dgamma(exp(logLambdaInverse),a_lam,1/b_lam,1) //-logLambdaInverse*logLambdaInverse/(2*b_lambda) 
              - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
              - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
              
              
              
              
              
      } // end inner loop
      
      
      
      
      // save posterior samples
      for(int m=0; m<M; m++){
        // partition thetaStar* into rho and theta
        rhokeep.submat(s_outer,m,s_outer,m) = sum(square(vectorise(thetaStar[m])));                             // rho=||theta*||^2
        thetakeepfield[m].row(s_outer) =  thetaStar[m].t()/sqrt(as_scalar(rhokeep.submat(s_outer,m,s_outer,m)));// theta=||theta||^{-1}theta*
        
        
      }
      logLambdaInversekeep(s_outer) = logLambdaInverse;       // save lambda^{-1}  
      logLambdaBInversekeep(s_outer) = logLambdaBInverse;     // save lambda^{-1}  
      
      
      // draw sigma^{-2} and gamma directly and save them
      sig2invkeep(s_outer) = rgamma( 1 , a_sig + .5*(N-P_z) , 1/(b_sig + 0.5*(Vcomps[1](0,0)-Vcomps[1](1,0))) )(0);  // function generates a gamma distribution parameterized with mean alpha*beta so we just inverted the second component;
      // gammakeep.row(s_outer) = arma::solve(trimatu(C2) , C2B +  as<arma::vec>(rnorm(P_z))/sqrt(sig2invkeep(s_outer) )).t(); // draw gamma by drawing MVN(0,I) then scaling by sqrt of variance and shifting by mean
      SIG = inv(Vcomps[4]*Vcomps[4].t());
      
      gammakeep.row(s_outer) = (SIG*Vcomps[2]+chol(SIG,"lower")*as<arma::vec>(rnorm(P_z))/sqrt(sig2invkeep(s_outer) )).t();
      
      // compute fitted h after burn-in
      if(s_outer>=n_burn){
        laminvK = Vcomps[5];      
        laminvKSiginv = solve(Vcomps[0] , laminvK).t();                                     // lambda^{-1} K Sigma^{-1} which is equal to t(Sigma^{-1} lambda^{-1} K) 
        hmeantemp = laminvKSiginv * (yz.col(0)-yz.cols(1,P_z)*gammakeep.row(s_outer).t());  // temporary vector to be used in mean and covariance
        // hcovtemp = laminvKSiginv*(Vcomps[0]-laminvK)/sig2invkeep(s_outer);               // temporary covariance matrix to be used in marginal covariance
        hcovtemp = (laminvK-laminvKSiginv*laminvK.t())/sig2invkeep(s_outer);                // temporary covariance matrix to be used in marginal covariance
        // compute marginal mean and covariance for h
        hmean += hmeantemp/(n_outer-n_burn);                                                                    // law of total expecation, E[E[h|s]]
        hcov += hcovtemp/(n_outer-n_burn) + hmeantemp*hmeantemp.t()/(n_outer-n_burn); // first two components of law of total variance: E(Var(h|s))+E[E[h|s]^2]-E[E[h|s]]^2
        // old (no B) : hcov += laminvKSiginv/(n_outer-n_burn)/sig2invkeep(s_outer) + hmeantemp*hmeantemp.t()/(n_outer-n_burn); // first two components of law of total variance: E(Var(h|s))+E[E[h|s]^2]-E[E[h|s]]^2
        // hcov += laminvKSiginv*(Vcomps[0]-laminvK)/(n_outer-n_burn)/sig2invkeep(s_outer) + hmeantemp*hmeantemp.t()/(n_outer-n_burn); // first two components of law of total variance: E(Var(h|s))+E[E[h|s]^2]-E[E[h|s]]^2
        if(draw_h==TRUE){
          // note: "std" method used to avoid errors---to return to default, remove "std"
          arma::svd(svdU,svdD,svdV,hcovtemp,"std"); //arma::svd(svdU,svdD,svdV,hcovtemp);                                                 // SVD decomposition for covariance of h (since it's more stable than cholesky)
          // draw h
          hsampkeep.row(s_outer) = ( hmeantemp+(svdU*diagmat(svdD))*arma::randn(N) ).t();     // hsampkeep.row(s_outer) = ( hmeantemp+arma::chol(hcovtemp)*arma::randn(N) ).t();
        }
      }
      
    } // end outer loop
    
    
    
    
    // Convert field to list
    for(int m=0; m<M; m++){
      thetakeep[m] = thetakeepfield[m];
    }
    
    // return a list
    return List::create(Named("theta")          = thetakeep,
                        Named("lambdaInverse")  = exp(logLambdaInversekeep),
                        Named("lambdaBInverse") = exp(logLambdaBInversekeep),
                        Named("rho")            = rhokeep,
                        Named("sigma2")         = 1/sig2invkeep,
                        Named("gamma")          = gammakeep,
                        Named("hsamp")          = hsampkeep,
                        Named("hmean")          = hmean,
                        Named("hcov")           = hcov - ( hmean*hmean.t())); // subtracting third component of law of total variance (i.e. E[E[h|s]]^2)
}



// [[Rcpp::depends(RcppArmadillo)]]
//' MCMC for BSMIM with spike and slab prior (inv-uniform slab for rho version)
//' 
//' Returns posterior draws for \eqn{\theta},\eqn{\lambda^{-1}},\eqn{\rho},\eqn{\sigma^2},\eqn{\gamma}, as well as the mean(\eqn{h}) and cov(\eqn{h})
//' Adapted from "bkmrdlm_multi.cpp" in the "regimes" package by Ander Wilson.
//' Note that this implementation with horseshoe2 prior returns posterior draws of \eqn{\rho} as a list of M matrices (unlike main MCMC code, which returns one matrix of M columns)
//' 
//' @param yz matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
//' @param Xlist list of (N by L_m) matrices representing X_m
//' @param a_lam shape hyperparameter for \eqn{\lambda^{-1}}
//' @param b_lam rate hyperparameter for \eqn{\lambda^{-1}}
//' @param b_lambdaB hyperparameter for \eqn{\lambda^{-1}_B}
//' @param a_sig first hyperparameter for \eqn{\sigma^{-2}}
//' @param b_sig second hyperparameter for \eqn{\sigma^{-2}}
//' @param tau02 hyperparameter for tau2; classic horseshoe is 1
//' @param kappa vector of hyperparameters \eqn{\kappa_m} for \eqn{\theta^*_m}
//' @param poly 0 = gaussian kernel, 1 = polynomial kernel
//' @param d degree of polynomial kernel
//' @param horseshoe 0 = no selection, 1 = componentwise horseshoe priors
//' @param randint 0 = no random intercepts, 1 = random intercepts model
//' @param Bmat N xN block diagonal matrix indicating cluster membership for random intercepts model
//' @param draw_h 0 = dont draw h, 1 = draw h
//' @param n_inner no. of MCMC iterations to run in the inner loop. n_outer*n_inner iteraction will be run.
//' @param n_outer no. of MCMC iterations to run in the outer loop. n_outer iterations will be saved.
//' @param n_burn no. of MCMC iterations to discard as burn-in
//' @author Glen McGee and Ander Wilson (adapted from the "regimes" package by Ander Wilson).
//' @export
// [[Rcpp::export]]
List bsmim_spikeslab_gaussprior_mcmc2(const arma::mat&    yz,         // matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
                                      const Rcpp::List&   Xlist,      // list of (N by L_m) matrices representing X_m
                                      const double&       a_lam,   // shape hyperparameter for \lambda^{-1}
                                      const double&       b_lam,    // rate hyperparameter for \lambda^{-1}
                                      const double&       b_lambdaB,  // hyperparameter for \lambda^{-1}_B
                                      const double&       a_sig,      // first hyperparameter for \sigma^{-2}
                                      const double&       b_sig,      // second hyperparameter for \sigma^{-2}
                                      const double&       s_theta,    // hyperparameter for sd of gaussian prior for thetastar (slab component)
                                      const double&       step_theta, // step size for theta* random walk (move 2)
                                      const double&       a_pi,       // first hyperparameter for beta distribution of pi 
                                      const double&       b_pi,       // second hyperparameter for beta distribution of pi 
                                      const bool&         poly,       // 0=gaussian kernel / 1=polynomial kernel
                                      const int&          d,          // degree of polynomial kernel
                                      const bool&         randint,    // 0=no random intercepts / 1=random intercepts
                                      const arma::mat&    Bmat,       // block diagonal matrix B indicating cluster membership for random intercepts model
                                      const bool&         draw_h,     // 0=dont sample h, 1= sample h
                                      const int&          n_inner,    // no. of MCMC iterations to run in the inner loop. n_outer*n_inner iteraction will be run.
                                      const int&          n_outer,    // no. of MCMC iterations to run in the outer loop. n_outer iterations will be saved.
                                      const int&          n_burn) {   // no. of MCMC iterations to discard as burn-in
  
  // step size for random walk in move 2 of MH step (spike/slab)
  // double step_theta = 0.2;
  // // TRYING DIFFERENT LAMBDA PRIOR
  // double a_lam = 1.0;
  // double b_lam = 0.1;
  
  // // inf-unif prior for rho
  // double ar = 0;
  // double br = s_theta;
  
  
  
  // Dimensions
  int N = yz.n_rows;              // no. of observations
  int M = Xlist.size();           // no. of indices (dimension of h(.))
  int P_z = yz.n_cols-1;          // no. of covariates
  IntegerVector Lvec = rep(0,M);  // vector of L_m representing the no. of columns of X_m
  
  // Initialize
  arma::field<arma::mat> X(M);                    // field of (N by L_m) matrices representing X_m
  arma::field<arma::vec> thetaStar(M);            // field of L_m vectors representing \thetaStar_m
  arma::field<arma::vec> delta(M);                // field of L_m vectors representing \delta_m
  arma::mat  XthetaStar(N,M,arma::fill::zeros);   // (N by M) matrix representing [X_1\thetaStar_1,...,X_M\thetaStar_M]
  for(int m=0; m<M; m++){
    X[m] = as<arma::mat>(Xlist[m]);                 // mth matrix in X is X_m
    Lvec(m) = X[m].n_cols;                          // mth element of Lvec is L_m (the length of \thetaStar_m)
    thetaStar[m] = as<arma::vec>(rnorm(Lvec(m)));   // set starting value for \thetaStar_m
    XthetaStar.col(m) = X[m] * thetaStar[m];
    
    arma::vec deltatemp(Lvec(m),arma::fill::ones); // temporary vec with L_m elements to store posteror draws of delta
    delta[m] = deltatemp;                          // set starting value for \delta_m
  }
  // // starting vals must respect bounds of prior
  // for(int m=0; m<M; m++){
  //   for(int l=0; l<Lvec(m); l++){
  //     // Rcout << thetaStar[m](l) << std::endl;
  //     // Rcout << pow(br,-0.5) << std::endl;
  //     if(thetaStar[m](l)>-pow(br,-0.5) & thetaStar[m](l)<pow(br,-0.5)){
  //       // Rcout << "m " << m << "l " << l << std::endl;      
  //       thetaStar[m](l) = rinttruncnorm(0,step_theta,-pow(br,-0.5),pow(br,-0.5));
  //       // Rcout << thetaStar[m](l) << std::endl;
  //     }
  //   }
  // }
  
  
  double logLambdaInverse = rnorm(1)(0)*sqrt(1)/10;    // set starting value for log(\lambda^{-1})
  double logLambdaBInverse = 0;                               // set log(\lambda^{-1}_B) to 0 if no random intercepts
  if(randint==TRUE){
    logLambdaBInverse = rnorm(1)(0)*sqrt(b_lambdaB)/10;       // set starting value for log(\lambda^{-1}_B) if using random intercepts
  }
  
  
  
  // Construct Sigma=I+\lambda^{-1}*K, A1, A2, C1, C2
  arma::field<arma::mat> Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat);
  
  // First evaluation of the log-liklihood (ll)
  double ll = logLambdaInverse+R::dgamma(exp(logLambdaInverse),a_lam,1/b_lam,1) //-logLambdaInverse*logLambdaInverse/(2*b_lambda) 
    - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
    - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
    
    
    // Initialize other values
    double logLambdaInverse_PROP, llu, ll_PROP, logratio;
    double llB, logLambdaBInverse_PROP, llB_PROP;
    arma::vec ons(N, arma::fill::ones);           // vector of length N of ones
    arma::vec vec_PROP;                           // proposal vector to define ellipse
    
    int MHmove, sumdelta, sumdelta_PROP, move_which, move_id;    // move type for joint MH draw of delta and theta* (1 or 2); sum of delta=1; which component MHmove applies to; id to sum
    double lltheta, lltheta_PROP, logdiff_PROP;
    arma::field<arma::vec> thetaStar_PROP = thetaStar;
    arma::field<arma::vec> delta_PROP = delta;
    
    
    // Store posterior samples
    arma::vec logLambdaInversekeep(n_outer,arma::fill::zeros);  // vector to store posterior draws of lambda^{-1}
    arma::vec logLambdaBInversekeep(n_outer,arma::fill::zeros); // vector to store posterior draws of lambda_B^{-1}
    arma::mat rhokeep(n_outer,M,arma::fill::zeros);             // matrix to store posterior draws of rho // each row is rho_1,...,rho_m
    arma::vec sig2invkeep(n_outer,arma::fill::zeros);           // vector to store posterior draws of sigma^{-2}
    arma::mat gammakeep(n_outer,P_z,arma::fill::zeros);         // matrix to store posterior draws of gamma // each row is gamma_1,...,gamma_{P_z}
    arma::field<arma::mat>  thetakeepfield(M);                  // field of matrices to store posterior draws of theta (for computation)
    Rcpp::List  thetakeep(M);                                   // List of matrices to store posteror draws of theta (same as above but for exporting to R)
    for(int m=0; m<M; m++){
      arma::mat thetakeeptemp(n_outer,Lvec(m),arma::fill::zeros); // temporary matrix with L_m columns to store posteror draws of theta_m
      thetakeepfield[m] = thetakeeptemp;                          // initialize elements of thetakeepfield to be the right size
      thetakeep[m] = thetakeeptemp;                               // initialize elements of thetakeepfield to be the right size
      
    } 
    
    // Initialize values for posterior of h(.)
    arma::mat laminvK(N,N,arma::fill::zeros);           // lambda^{-1}K where K is the kernel matrix
    arma::mat laminvKSiginv(N,N,arma::fill::zeros);     // \lambda^{-1}KSigma^{-1}=laminvK*Sigma^{-1} or t(solve(Sigma)*laminvK)
    arma::vec hmean(N,arma::fill::zeros);               // marginal mean vector for h(.)
    arma::mat hcov(N,N,arma::fill::zeros);              // marginal variance-covariance matrix of h(.)
    arma::vec hmeantemp(N,arma::fill::zeros);           // temporary vector for computing marginal covariance of h(.)
    arma::mat hcovtemp(N,N,arma::fill::zeros);          // marginal variance-covariance matrix of h(.)
    arma::mat SIG(P_z,P_z,arma::fill::zeros);           // temporary matrix
    arma::mat hsampkeep(n_outer,N,arma::fill::zeros);   // matrix to store posterior draws of h: rows are samples of h
    arma::mat svdU;                                     // for SVD decomposition
    arma::vec svdD;                                     // for SVD decomposition
    arma::mat svdV;                                     // for SVD decomposition
    
    // -------------------------------------
    // MCMC
    // -------------------------------------
    for(int s_outer=0; s_outer<n_outer; s_outer++){                       // outer loop
      
      if(s_outer % 10 == 0) Rcout << "#" << s_outer << std::endl;         // output progress every 10 iterations
      
      for(int s_inner=0; s_inner<n_inner; s_inner++){                     // inner loop
        
        
        // -------------------------------------------------------------------------- //
        // UPDATE LAMBDA^{-1} via random walk metropolis
        logLambdaInverse_PROP = logLambdaInverse + rnorm(1)(0); // random walk step
        
        // update Sigma=I+\lambda^{-1}*K, A1, A2, C1, C2
        Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse_PROP, XthetaStar, randint, logLambdaBInverse, Bmat);   
        
        // proposed ll  (inverting b_lam, since this takes scale instead of rate)
        ll_PROP = logLambdaInverse_PROP+R::dgamma(exp(logLambdaInverse_PROP),a_lam,1/b_lam,1) //-logLambdaInverse*logLambdaInverse/(2*b_lambda) -logLambdaInverse_PROP*logLambdaInverse_PROP/(2*b_lambda)
          - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
          - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
          
          // log ratio for MH acceptance
          logratio = ll_PROP-ll;
          if(logratio>0){
            logratio = 0;
          }
          
          // accept or reject    
          if(logratio >  log(runif(1)(0)) ){                
            logLambdaInverse = logLambdaInverse_PROP; 
            ll = ll_PROP;
          }
          
          // Reset matrix Sigma and ll
          Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat); 
          
          ll = logLambdaInverse+R::dgamma(exp(logLambdaInverse),a_lam,1/b_lam,1) //-logLambdaInverse*logLambdaInverse/(2*b_lambda)
            - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
            - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
            
            
            
            
            // -------------------------------------------------------------------------- //
            // UPDATE LAMBDA_B^{-1} via random walk metropolis
            if(randint==TRUE){ // only update lambda_B^{-1} if using random intercepts model
              
              logLambdaBInverse_PROP = logLambdaBInverse + rnorm(1)(0); // random walk step
              
              // current llB (uses the prior for logLambdaBInverse, hence the first line below is different than ll elsewhere)
              llB = logLambdaBInverse+R::dgamma(exp(logLambdaBInverse),a_lam,1/b_lam,1) //-logLambdaBInverse*logLambdaBInverse/(2*b_lambdaB) 
                - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
                
                // update Sigma=I+\lambda^{-1}*K, A1, A2, C1, C2
                Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse_PROP, Bmat);   
                
                // proposed llB (uses the prior for logLambdaBInverse, hence the first line below is different than ll elsewhere)
                llB_PROP = logLambdaBInverse_PROP+R::dgamma(exp(logLambdaBInverse_PROP),a_lam,1/b_lam,1) //-logLambdaBInverse_PROP*logLambdaBInverse_PROP/(2*b_lambdaB)
                  - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                  - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
                  
                  // log ratio for MH acceptance
                  logratio = llB_PROP-llB;
                  if(logratio>0){
                    logratio = 0;
                  }
                  
                  // accept or reject    
                  if(logratio >  log(runif(1)(0)) ){                
                    logLambdaBInverse = logLambdaBInverse_PROP; 
                  }
                  
                  // Reset matrix Sigma and ll
                  Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat); 
                  
                  // resetting ll (not llB)
                  ll = logLambdaInverse+R::dgamma(exp(logLambdaInverse),a_lam,1/b_lam,1) //-logLambdaInverse*logLambdaInverse/(2*b_lambda) 
                    - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                    - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
                    
            }
            
            
            
            // -------------------------------------------------------------------------- //
            // UPDATE delta and thetaStar via Metropolis Hastings
            
            // count number of delta==1
            sumdelta = 0;
            for(int m=0; m<M; m++){
              sumdelta += sum(delta[m]);
            }
            
            // select move type
            if(sumdelta==0){ // must choose move 1 if all delta==0 (move 2 samples from the deltas equal to 1)
              MHmove = 1;
            }else{ // randomly select move type
              MHmove = Rcpp::rbinom(1,1,0.5)(0);
            }
            
            // Make move 1 or 2
            if(MHmove==1){  // MHmove version 1
              
              // reset proposals
              delta_PROP = delta;
              thetaStar_PROP = thetaStar;
              
              // randomly select a component
              move_which = sample(sum(Lvec), 1)(0);
              
              move_id = 0;                                // set index
              for(int m=0; m<M; m++){
                for(int l=0; l<Lvec(m); l++){
                  
                  move_id += 1;
                  if(move_id==move_which){   // make move for randomly selected component
                    
                    // compute lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
                    lltheta = log(tgamma(sumdelta+a_pi))+log(tgamma(sum(Lvec)-sumdelta+b_pi))
                    + delta[m](l)*R::dnorm(thetaStar[m](l),0.0,s_theta,TRUE)//+(1-delta[m][l])*log(1)//-logLambdaInverse*logLambdaInverse/(2*b_lambda)
                      - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                      - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
                      
                      // propose new delta: 1 to 0 or 0 to 1
                      delta_PROP[m](l) = 1-delta[m](l);
                      if(delta_PROP[m](l)==1){ // adjust sum delta for switch (either increase by 1 or decrease by 1)
                        sumdelta_PROP = sumdelta+1;
                      }else{
                        sumdelta_PROP = sumdelta-1;
                      }
                      
                      //propose new theta* (based on new delta)
                      if(delta_PROP[m](l)==0){    // if delta_PROP, then theta*=0
                        thetaStar_PROP[m](l) = 0;
                      }else{                      // otherwise draw from Q1 (which is the prior)
                        thetaStar_PROP[m](l) = s_theta*rnorm(1)(0);
                      }
                      
                      // update components that depend on thetaStar
                      XthetaStar.col(m) = X[m] * thetaStar_PROP[m];                          // update XthetaStar with proposed thetaStar
                      
                      // update Sigma with proposed thetaStar
                      Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat);
                      
                      // proposed lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
                      lltheta_PROP = log(tgamma(sumdelta_PROP+a_pi))+log(tgamma(sum(Lvec)-sumdelta_PROP+b_pi))
                        + delta_PROP[m](l)*R::dnorm(thetaStar_PROP[m](l),0.0,s_theta,TRUE)//+(1-delta[m][l])*log(1)//-logLambdaInverse*logLambdaInverse/(2*b_lambda)
                          - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                          - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
                          
                          // difference between proposal distributions: logP(theta_prop,delta_prop|theta,delta)-logP(theta,delta|theta_prop,delta_prop)
                          logdiff_PROP = delta_PROP[m](l)*R::dnorm(thetaStar_PROP[m](l),0.0,s_theta,TRUE) //   logP(theta_prop,delta_prop|theta,delta)
                            - delta[m](l)*R::dnorm(thetaStar[m](l),0.0,s_theta,TRUE) // logP(theta,delta|theta_prop,delta_prop)
                              +((sumdelta_PROP==0)*log(1)+(sumdelta_PROP!=0)*log(0.5)) // extra difference component because sometimes we must choose move 1.
                              -((sumdelta==0)*log(1)+(sumdelta!=0)*log(0.5));
                              
                              // log ratio for MH acceptance
                              logratio = (lltheta_PROP-lltheta)-logdiff_PROP; // the extra (negative!) difference in log proposals since this is MH not metropolis
                              if(logratio>0){
                                logratio = 0;
                              }
                              
                              // accept or reject
                              if(logratio >  log(runif(1)(0)) ){
                                delta[m](l) = delta_PROP[m](l);
                                thetaStar[m](l) = thetaStar_PROP[m](l);
                                // ll = ll_PROP;
                              }else{
                                XthetaStar.col(m) = X[m] * thetaStar[m];
                              }
                              
                              
                  } // end (if move_id)
                  
                  
                  
                } // l loop
              } // m loop
              
              
              
              
              
            }else{ // MHmove version 2
              
              
              move_which = sample(sumdelta, 1)(0);    // randomly select component for move (such that delta==1)
              
              move_id = 0;                                // set index
              for(int m=0; m<M; m++){
                for(int l=0; l<Lvec(m); l++){
                  
                  if(delta[m](l)==1){
                    move_id += 1;
                    if(move_id==move_which){   //RANDOM WALK STEP
                      
                      
                      // compute lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
                      lltheta = R::dnorm(thetaStar[m](l),0.0,s_theta,TRUE) //-logLambdaInverse*logLambdaInverse/(2*b_lambda)
                      - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                      - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
                      
                      // propose new theta
                      thetaStar_PROP[m](l) = thetaStar[m](l) + step_theta*rnorm(1)(0);
                      
                      // update components that depend on thetaStar
                      XthetaStar.col(m) = X[m] * thetaStar_PROP[m];                          // update XthetaStar with proposed thetaStar
                      
                      // update Sigma with proposed thetaStar
                      Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat);
                      
                      // proposed lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
                      lltheta_PROP = R::dnorm(thetaStar_PROP[m](l),0.0,s_theta,TRUE) //-logLambdaInverse*logLambdaInverse/(2*b_lambda)
                        - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                        - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
                        
                        
                        // log ratio for MH acceptance
                        logratio = lltheta_PROP-lltheta;
                        if(logratio>0){
                          logratio = 0;
                        }
                        
                        // accept or reject
                        if(logratio >  log(runif(1)(0)) ){
                          thetaStar[m](l) = thetaStar_PROP[m](l);
                          // ll = ll_PROP;
                        }else{
                          XthetaStar.col(m) = X[m] * thetaStar[m];
                        }
                        
                        
                        
                    } // end (if move_id)
                  } // end (if delta=1)
                  
                  
                  
                } // l loop
              } // m loop
            } // end move 2
            
            // Reset matrix Sigma and ll
            Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat); 
            
            ll = logLambdaInverse+R::dgamma(exp(logLambdaInverse),a_lam,1/b_lam,1) //-logLambdaInverse*logLambdaInverse/(2*b_lambda) 
              - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
              - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
              
              
              
              
              
      } // end inner loop
      
      
      
      
      // save posterior samples
      for(int m=0; m<M; m++){
        // partition thetaStar* into rho and theta
        rhokeep.submat(s_outer,m,s_outer,m) = sum(square(vectorise(thetaStar[m])));                             // rho=||theta*||^2
        thetakeepfield[m].row(s_outer) =  thetaStar[m].t()/sqrt(as_scalar(rhokeep.submat(s_outer,m,s_outer,m)));// theta=||theta||^{-1}theta*
        
        
      }
      logLambdaInversekeep(s_outer) = logLambdaInverse;       // save lambda^{-1}  
      logLambdaBInversekeep(s_outer) = logLambdaBInverse;     // save lambda^{-1}  
      
      
      // draw sigma^{-2} and gamma directly and save them
      sig2invkeep(s_outer) = rgamma( 1 , a_sig + .5*(N-P_z) , 1/(b_sig + 0.5*(Vcomps[1](0,0)-Vcomps[1](1,0))) )(0);  // function generates a gamma distribution parameterized with mean alpha*beta so we just inverted the second component;
      // gammakeep.row(s_outer) = arma::solve(trimatu(C2) , C2B +  as<arma::vec>(rnorm(P_z))/sqrt(sig2invkeep(s_outer) )).t(); // draw gamma by drawing MVN(0,I) then scaling by sqrt of variance and shifting by mean
      SIG = inv(Vcomps[4]*Vcomps[4].t());
      
      gammakeep.row(s_outer) = (SIG*Vcomps[2]+chol(SIG,"lower")*as<arma::vec>(rnorm(P_z))/sqrt(sig2invkeep(s_outer) )).t();
      
      // compute fitted h after burn-in
      if(s_outer>=n_burn){
        laminvK = Vcomps[5];      
        laminvKSiginv = solve(Vcomps[0] , laminvK).t();                                     // lambda^{-1} K Sigma^{-1} which is equal to t(Sigma^{-1} lambda^{-1} K) 
        hmeantemp = laminvKSiginv * (yz.col(0)-yz.cols(1,P_z)*gammakeep.row(s_outer).t());  // temporary vector to be used in mean and covariance
        // hcovtemp = laminvKSiginv*(Vcomps[0]-laminvK)/sig2invkeep(s_outer);               // temporary covariance matrix to be used in marginal covariance
        hcovtemp = (laminvK-laminvKSiginv*laminvK.t())/sig2invkeep(s_outer);                // temporary covariance matrix to be used in marginal covariance
        // compute marginal mean and covariance for h
        hmean += hmeantemp/(n_outer-n_burn);                                                                    // law of total expecation, E[E[h|s]]
        hcov += hcovtemp/(n_outer-n_burn) + hmeantemp*hmeantemp.t()/(n_outer-n_burn); // first two components of law of total variance: E(Var(h|s))+E[E[h|s]^2]-E[E[h|s]]^2
        // old (no B) : hcov += laminvKSiginv/(n_outer-n_burn)/sig2invkeep(s_outer) + hmeantemp*hmeantemp.t()/(n_outer-n_burn); // first two components of law of total variance: E(Var(h|s))+E[E[h|s]^2]-E[E[h|s]]^2
        // hcov += laminvKSiginv*(Vcomps[0]-laminvK)/(n_outer-n_burn)/sig2invkeep(s_outer) + hmeantemp*hmeantemp.t()/(n_outer-n_burn); // first two components of law of total variance: E(Var(h|s))+E[E[h|s]^2]-E[E[h|s]]^2
        if(draw_h==TRUE){
          // note: "std" method used to avoid errors---to return to default, remove "std"
          arma::svd(svdU,svdD,svdV,hcovtemp,"std"); //arma::svd(svdU,svdD,svdV,hcovtemp);                                                 // SVD decomposition for covariance of h (since it's more stable than cholesky)
          // draw h
          hsampkeep.row(s_outer) = ( hmeantemp+(svdU*diagmat(svdD))*arma::randn(N) ).t();     // hsampkeep.row(s_outer) = ( hmeantemp+arma::chol(hcovtemp)*arma::randn(N) ).t();
        }
      }
      
    } // end outer loop
    
    
    
    
    // Convert field to list
    for(int m=0; m<M; m++){
      thetakeep[m] = thetakeepfield[m];
    }
    
    // return a list
    return List::create(Named("theta")          = thetakeep,
                        Named("lambdaInverse")  = exp(logLambdaInversekeep),
                        Named("lambdaBInverse") = exp(logLambdaBInversekeep),
                        Named("rho")            = rhokeep,
                        Named("sigma2")         = 1/sig2invkeep,
                        Named("gamma")          = gammakeep,
                        Named("hsamp")          = hsampkeep,
                        Named("hmean")          = hmean,
                        Named("hcov")           = hcov - ( hmean*hmean.t())); // subtracting third component of law of total variance (i.e. E[E[h|s]]^2)
}



// [[Rcpp::depends(RcppArmadillo)]]
//' MCMC for BSMIM with spike and slab prior (inv-uniform slab for rho version)
//' 
//' Returns posterior draws for \eqn{\theta},\eqn{\lambda^{-1}},\eqn{\rho},\eqn{\sigma^2},\eqn{\gamma}, as well as the mean(\eqn{h}) and cov(\eqn{h})
//' Adapted from "bkmrdlm_multi.cpp" in the "regimes" package by Ander Wilson.
//' Note that this implementation with horseshoe2 prior returns posterior draws of \eqn{\rho} as a list of M matrices (unlike main MCMC code, which returns one matrix of M columns)
//' 
//' @param yz matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
//' @param Xlist list of (N by L_m) matrices representing X_m
//' @param a_lam shape hyperparameter for \eqn{\lambda^{-1}}
//' @param b_lam rate hyperparameter for \eqn{\lambda^{-1}}
//' @param b_lambdaB hyperparameter for \eqn{\lambda^{-1}_B}
//' @param a_sig first hyperparameter for \eqn{\sigma^{-2}}
//' @param b_sig second hyperparameter for \eqn{\sigma^{-2}}
//' @param tau02 hyperparameter for tau2; classic horseshoe is 1
//' @param kappa vector of hyperparameters \eqn{\kappa_m} for \eqn{\theta^*_m}
//' @param poly 0 = gaussian kernel, 1 = polynomial kernel
//' @param d degree of polynomial kernel
//' @param horseshoe 0 = no selection, 1 = componentwise horseshoe priors
//' @param randint 0 = no random intercepts, 1 = random intercepts model
//' @param Bmat N xN block diagonal matrix indicating cluster membership for random intercepts model
//' @param draw_h 0 = dont draw h, 1 = draw h
//' @param n_inner no. of MCMC iterations to run in the inner loop. n_outer*n_inner iteraction will be run.
//' @param n_outer no. of MCMC iterations to run in the outer loop. n_outer iterations will be saved.
//' @param n_burn no. of MCMC iterations to discard as burn-in
//' @author Glen McGee and Ander Wilson (adapted from the "regimes" package by Ander Wilson).
//' @export
// [[Rcpp::export]]
List bsmim_informative_mcmc2(const arma::mat&    yz,         // matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
                                      const Rcpp::List&   Xlist,      // list of (N by L_m) matrices representing X_m
                                      const double&       a_lam,   // shape hyperparameter for \lambda^{-1}
                                      const double&       b_lam,    // rate hyperparameter for \lambda^{-1}
                                      const double&       b_lambdaB,  // hyperparameter for \lambda^{-1}_B
                                      const double&       a_sig,      // first hyperparameter for \sigma^{-2}
                                      const double&       b_sig,      // second hyperparameter for \sigma^{-2}
                                      const double&       s_theta,    // hyperparameter for sd of gaussian prior for thetastar (slab component)
                                      const double&       step_theta, // step size for theta* random walk (move 2)
                                      const double&       a_pi,       // first hyperparameter for beta distribution of pi 
                                      const double&       b_pi,       // second hyperparameter for beta distribution of pi 
                                      const bool&         poly,       // 0=gaussian kernel / 1=polynomial kernel
                                      const int&          d,          // degree of polynomial kernel
                                      const bool&         randint,    // 0=no random intercepts / 1=random intercepts
                                      const arma::mat&    Bmat,       // block diagonal matrix B indicating cluster membership for random intercepts model
                                      const bool&         draw_h,     // 0=dont sample h, 1= sample h
                                      const arma::vec&    thetaconstraint, // M-vector for type of constraints (0 is none, 1 is positive, 2 is dirichlet)
                                      const double        a_slabpos,  // shape for gamma distribution (for gamma slab of positivity constraint); default=4
                                      const double        b_slabpos,  // rate for gamma distribution(for gamma slab of positivity constraint); default=2
                                      const Rcpp::List&   alphas,     // list of Lm vectors representing alpha hyperparameters for dirichlet prior
                                      const double        a_rho,      // shape for gamma distribution (rho^{1/2}) ## default=?
                                      const double        b_rho,      // rate for gamma distribution (rho^{1/2}) ## default=?
                                      const int&          n_inner,    // no. of MCMC iterations to run in the inner loop. n_outer*n_inner iteraction will be run.
                                      const int&          n_outer,    // no. of MCMC iterations to run in the outer loop. n_outer iterations will be saved.
                                      const int&          n_burn) {   // no. of MCMC iterations to discard as burn-in
  
  
  // Dimensions
  int N = yz.n_rows;              // no. of observations
  int M = Xlist.size();           // no. of indices (dimension of h(.))
  int P_z = yz.n_cols-1;          // no. of covariates
  IntegerVector Lvec = rep(0,M);  // vector of L_m representing the no. of columns of X_m
  
  
  // Initialize
  arma::field<arma::mat> X(M);                    // field of (N by L_m) matrices representing X_m
  arma::field<arma::vec> thetaStar(M);            // field of L_m vectors representing \thetaStar_m
  arma::field<arma::vec> delta(M);                // field of L_m vectors representing \delta_m
  arma::mat  XthetaStar(N,M,arma::fill::zeros);   // (N by M) matrix representing [X_1\thetaStar_1,...,X_M\thetaStar_M]
  arma::field<arma::vec> alpha(M);                // field of Lm vectors representing alpha hyperparameters for dirichlet
  for(int m=0; m<M; m++){
    X[m] = as<arma::mat>(Xlist[m]);                 // mth matrix in X is X_m
    Lvec(m) = X[m].n_cols;                          // mth element of Lvec is L_m (the length of \thetaStar_m)
    thetaStar[m] = as<arma::vec>(rnorm(Lvec(m)));   // set starting value for \thetaStar_m
    if(thetaconstraint[m]>0){
      thetaStar[m] = pow(thetaStar[m],2);   // set starting value for \thetaStar_m
    }
    XthetaStar.col(m) = X[m] * thetaStar[m];
    alpha[m] = as<arma::vec>(alphas[m]);   // set starting value for \thetaStar_m
    
    arma::vec deltatemp(Lvec(m),arma::fill::ones); // temporary vec with L_m elements to store posteror draws of delta
    delta[m] = deltatemp;                          // set starting value for \delta_m
  }
  
  // count constraints
  int sumM0 = 0;
  int sumM1 = 0;
  int sumM2 = 0;
  for(int m=0; m<M; m++){
    // only loop through components with positivity constraints
    if(thetaconstraint[m]==0){
      sumM0 += Lvec(m);
    }
    if(thetaconstraint[m]==1){
      sumM1 += Lvec(m);
    }
    if(thetaconstraint[m]==2){
      sumM2 += Lvec(m);
    }
  }
  
  // // starting vals must respect bounds of prior
  // for(int m=0; m<M; m++){
  //   for(int l=0; l<Lvec(m); l++){
  //     // Rcout << thetaStar[m](l) << std::endl;
  //     // Rcout << pow(br,-0.5) << std::endl;
  //     if(thetaStar[m](l)>-pow(br,-0.5) & thetaStar[m](l)<pow(br,-0.5)){
  //       // Rcout << "m " << m << "l " << l << std::endl;      
  //       thetaStar[m](l) = rinttruncnorm(0,step_theta,-pow(br,-0.5),pow(br,-0.5));
  //       // Rcout << thetaStar[m](l) << std::endl;
  //     }
  //   }
  // }
  
  
  double logLambdaInverse = rnorm(1)(0)*sqrt(1)/10;    // set starting value for log(\lambda^{-1})
  double logLambdaBInverse = 0;                               // set log(\lambda^{-1}_B) to 0 if no random intercepts
  if(randint==TRUE){
    logLambdaBInverse = rnorm(1)(0)*sqrt(b_lambdaB)/10;       // set starting value for log(\lambda^{-1}_B) if using random intercepts
  }
  
  
  
  // Construct Sigma=I+\lambda^{-1}*K, A1, A2, C1, C2
  arma::field<arma::mat> Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat);
  
  // First evaluation of the log-liklihood (ll)
  double ll = logLambdaInverse+R::dgamma(exp(logLambdaInverse),a_lam,1/b_lam,1) //-logLambdaInverse*logLambdaInverse/(2*b_lambda) 
    - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
    - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
    
    
    // Initialize other values
    double logLambdaInverse_PROP, llu, ll_PROP, logratio;
    double llB, logLambdaBInverse_PROP, llB_PROP;
    arma::vec ons(N, arma::fill::ones);           // vector of length N of ones
    arma::vec vec_PROP;                           // proposal vector to define ellipse
    
    int MHmove, sumdelta0, sumdelta1, sumdelta_PROP, move_which, move_id;    // move type for joint MH draw of delta and theta* (1 or 2); sum of delta=1; which component MHmove applies to; id to sum
    double lltheta, lltheta_PROP, logdiff_PROP;
    arma::field<arma::vec> thetaStar_PROP = thetaStar;
    arma::field<arma::vec> delta_PROP = delta;  //
    
    
    // Store posterior samples
    arma::vec logLambdaInversekeep(n_outer,arma::fill::zeros);  // vector to store posterior draws of lambda^{-1}
    arma::vec logLambdaBInversekeep(n_outer,arma::fill::zeros); // vector to store posterior draws of lambda_B^{-1}
    arma::mat rhokeep(n_outer,M,arma::fill::zeros);             // matrix to store posterior draws of rho // each row is rho_1,...,rho_m
    arma::vec sig2invkeep(n_outer,arma::fill::zeros);           // vector to store posterior draws of sigma^{-2}
    arma::mat gammakeep(n_outer,P_z,arma::fill::zeros);         // matrix to store posterior draws of gamma // each row is gamma_1,...,gamma_{P_z}
    arma::field<arma::mat>  thetakeepfield(M);                  // field of matrices to store posterior draws of theta (for computation)
    Rcpp::List  thetakeep(M);                                   // List of matrices to store posteror draws of theta (same as above but for exporting to R)
    arma::field<arma::mat>  thetaPOSkeepfield(M);                  // field of matrices to store posterior draws of theta (for computation)
    Rcpp::List  thetaPOSkeep(M);                                   // List of matrices to store posteror draws of theta (same as above but for exporting to R)
    for(int m=0; m<M; m++){
      arma::mat thetakeeptemp(n_outer,Lvec(m),arma::fill::zeros); // temporary matrix with L_m columns to store posteror draws of theta_m
      thetakeepfield[m] = thetakeeptemp;                          // initialize elements of thetakeepfield to be the right size
      thetakeep[m] = thetakeeptemp;                               // initialize elements of thetakeep to be the right size 
      if(thetaconstraint(m)>0){
        thetaPOSkeepfield[m] = thetakeeptemp;                          // initialize elements of thetaPOSkeepfield to be the right size
        thetaPOSkeep[m] = thetakeeptemp;                               // initialize elements of thetaPOSkeep to be the right size
      }
    } 
    
    // Initialize values for posterior of h(.)
    arma::mat laminvK(N,N,arma::fill::zeros);           // lambda^{-1}K where K is the kernel matrix
    arma::mat laminvKSiginv(N,N,arma::fill::zeros);     // \lambda^{-1}KSigma^{-1}=laminvK*Sigma^{-1} or t(solve(Sigma)*laminvK)
    arma::vec hmean(N,arma::fill::zeros);               // marginal mean vector for h(.)
    arma::mat hcov(N,N,arma::fill::zeros);              // marginal variance-covariance matrix of h(.)
    arma::vec hmeantemp(N,arma::fill::zeros);           // temporary vector for computing marginal covariance of h(.)
    arma::mat hcovtemp(N,N,arma::fill::zeros);          // marginal variance-covariance matrix of h(.)
    arma::mat SIG(P_z,P_z,arma::fill::zeros);           // temporary matrix
    arma::mat hsampkeep(n_outer,N,arma::fill::zeros);   // matrix to store posterior draws of h: rows are samples of h
    arma::mat svdU;                                     // for SVD decomposition
    arma::vec svdD;                                     // for SVD decomposition
    arma::mat svdV;                                     // for SVD decomposition
    
    // -------------------------------------
    // MCMC
    // -------------------------------------
    for(int s_outer=0; s_outer<n_outer; s_outer++){                       // outer loop
      
      if(s_outer % 10 == 0) Rcout << "#" << s_outer << std::endl;         // output progress every 10 iterations
      
      for(int s_inner=0; s_inner<n_inner; s_inner++){                     // inner loop
        
        
        // -------------------------------------------------------------------------- //
        // UPDATE LAMBDA^{-1} via random walk metropolis
        logLambdaInverse_PROP = logLambdaInverse + rnorm(1)(0); // random walk step
        
        // update Sigma=I+\lambda^{-1}*K, A1, A2, C1, C2
        Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse_PROP, XthetaStar, randint, logLambdaBInverse, Bmat);   
        
        // proposed ll  (inverting b_lam, since this takes scale instead of rate)
        ll_PROP = logLambdaInverse_PROP+R::dgamma(exp(logLambdaInverse_PROP),a_lam,1/b_lam,1) //-logLambdaInverse*logLambdaInverse/(2*b_lambda) -logLambdaInverse_PROP*logLambdaInverse_PROP/(2*b_lambda)
          - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
          - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
          
          // log ratio for MH acceptance
          logratio = ll_PROP-ll;
          if(logratio>0){
            logratio = 0;
          }
          
          // accept or reject    
          if(logratio >  log(runif(1)(0)) ){                
            logLambdaInverse = logLambdaInverse_PROP; 
            ll = ll_PROP;
          }
          
          // Reset matrix Sigma and ll
          Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat); 
          
          ll = logLambdaInverse+R::dgamma(exp(logLambdaInverse),a_lam,1/b_lam,1) //-logLambdaInverse*logLambdaInverse/(2*b_lambda)
            - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
            - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
            
            
            
            
            // -------------------------------------------------------------------------- //
            // UPDATE LAMBDA_B^{-1} via random walk metropolis
            if(randint==TRUE){ // only update lambda_B^{-1} if using random intercepts model
              
              logLambdaBInverse_PROP = logLambdaBInverse + rnorm(1)(0); // random walk step
              
              // current llB (uses the prior for logLambdaBInverse, hence the first line below is different than ll elsewhere)
              llB = logLambdaBInverse+R::dgamma(exp(logLambdaBInverse),a_lam,1/b_lam,1) //-logLambdaBInverse*logLambdaBInverse/(2*b_lambdaB) 
                - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
                
                // update Sigma=I+\lambda^{-1}*K, A1, A2, C1, C2
                Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse_PROP, Bmat);   
                
                // proposed llB (uses the prior for logLambdaBInverse, hence the first line below is different than ll elsewhere)
                llB_PROP = logLambdaBInverse_PROP+R::dgamma(exp(logLambdaBInverse_PROP),a_lam,1/b_lam,1) //-logLambdaBInverse_PROP*logLambdaBInverse_PROP/(2*b_lambdaB)
                  - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                  - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
                  
                  // log ratio for MH acceptance
                  logratio = llB_PROP-llB;
                  if(logratio>0){
                    logratio = 0;
                  }
                  
                  // accept or reject    
                  if(logratio >  log(runif(1)(0)) ){                
                    logLambdaBInverse = logLambdaBInverse_PROP; 
                  }
                  
                  // Reset matrix Sigma and ll
                  Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat); 
                  
                  // resetting ll (not llB)
                  ll = logLambdaInverse+R::dgamma(exp(logLambdaInverse),a_lam,1/b_lam,1) //-logLambdaInverse*logLambdaInverse/(2*b_lambda) 
                    - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                    - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
                    
            }
            
            
            
            
            // -------------------------------------------------------------------------- //
            // UPDATE thetaStar for DIRICHLET CONSTRAINTS (thetaconstraints=2)
            if(sumM2>0){ 
              
              // reset proposals
              thetaStar_PROP = thetaStar;
              
              for(int m=0; m<M; m++){
                if(thetaconstraint(m)==2){// only loop through components with dirichlet constraints
                  for(int l=0; l<Lvec(m); l++){
                    
                    // first term is from change of vars from thetastar to logthetastar (for random walk)
                    lltheta = log(thetaStar[m](l)) //for change of variable from thetastar to logthetastar
                    +log(tgamma(sum(alpha[m])))-log(tgamma(alpha[m](l)))+(alpha[m](l)-1)*log(thetaStar[m](l))+(1-sum(alpha[m]))*log(sum(thetaStar[m]))  //dirichlet component for prior thetas (and jacobian from change of vars from original thetas)
                    +R::dgamma(sum(thetaStar[m]),a_rho,1/b_rho,TRUE) // prior component for sum thetastar (i.e. rho^{1/2})
                    - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                    - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
                    
                    // propose new log theta
                    thetaStar_PROP[m](l) = exp(log(thetaStar[m](l)) + step_theta*rnorm(1)(0));
                    
                    // update components that depend on thetaStar
                    XthetaStar.col(m) = X[m] * thetaStar_PROP[m];                          // update XthetaStar with proposed thetaStar
                    
                    // update Sigma with proposed thetaStar
                    Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat);
                    
                    // proposed lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
                    lltheta_PROP = log(thetaStar_PROP[m](l)) //for change of variable from thetastar to logthetastar
                      +log(tgamma(sum(alpha[m])))-log(tgamma(alpha[m](l)))+(alpha[m](l)-1)*log(thetaStar_PROP[m](l))+(1-sum(alpha[m]))*log(sum(thetaStar_PROP[m]))  //dirichlet component for prior thetas (and jacobian from change of vars from original thetas)
                      +R::dgamma(sum(thetaStar_PROP[m]),a_rho,1/b_rho,TRUE) // prior component for sum thetastar (i.e. rho^{1/2})
                      - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                      - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
                      
                      
                    // log ratio for MH acceptance
                    logratio = lltheta_PROP-lltheta;
                    if(logratio>0){
                      logratio = 0;
                    }
                    
                    // accept or reject
                    if(logratio >  log(runif(1)(0)) ){
                      thetaStar[m](l) = thetaStar_PROP[m](l);
                      // ll = ll_PROP;
                    }else{
                      XthetaStar.col(m) = X[m] * thetaStar[m];
                    }
                    
                  } // end l loop
                } // end if (thetaconstraints)
              } // end m loop
              
            } // end update (dirichlet constraints)
            
            
            
            
            
            
            // -------------------------------------------------------------------------- //
            // UPDATE delta and thetaStarFOR POSITIVE CONSTRAINTS

            if(sumM1>0){
              // count number of delta==1
              sumdelta1 = 0;
              for(int m=0; m<M; m++){
                // only loop through components with positivity constraints
                if(thetaconstraint(m)==1){
                  sumdelta1 += sum(delta[m]);
                }
              }

              // select move type
              if(sumdelta1==0){ // must choose move 1 if all delta==0 (move 2 samples from the deltas equal to 1)
                MHmove = 1;
              }else{ // randomly select move type
                MHmove = Rcpp::rbinom(1,1,0.5)(0);
              }
              
              // Make move 1 or 2
              if(MHmove==1){  // MHmove version 1
                
                // reset proposals
                delta_PROP = delta;
                thetaStar_PROP = thetaStar;
                
                // randomly select a component
                move_which = sample(sumM1, 1)(0);
                
                move_id = 0;                                // set index
                for(int m=0; m<M; m++){
                  if(thetaconstraint(m)==1){// only loop through components with positivity constraints
                    for(int l=0; l<Lvec(m); l++){
                      
                      
                      move_id += 1;
                      
                      if(move_id==move_which){   // make move for randomly selected component
                         
                        // compute lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
                        lltheta = log(tgamma(sumdelta1+a_pi))+log(tgamma(sumM1-sumdelta1+b_pi))
                        //+ delta[m](l)*R::dgamma(thetaStar[m](l),a_slabpos,1/b_slabpos,TRUE)//+(1-delta[m][l])*log(1)//-logLambdaInverse*logLambdaInverse/(2*b_lambda)
                          - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                          - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
                        if(delta[m](l)==1){ //adding component (avoiding NAs, because dgamma at 0 can give Inf)
                          lltheta += delta[m](l)*R::dgamma(thetaStar[m](l),a_slabpos,1/b_slabpos,TRUE);
                        }  
                        
                        // propose new delta: 1 to 0 or 0 to 1
                        delta_PROP[m](l) = 1-delta[m](l);
                        if(delta_PROP[m](l)==1){ // adjust sum delta for switch (either increase by 1 or decrease by 1)
                          sumdelta_PROP = sumdelta1+1;
                        }else{
                          sumdelta_PROP = sumdelta1-1;
                        }
                        
                        //propose new theta* (based on new delta)
                        if(delta_PROP[m](l)==0){    // if delta_PROP, then theta*=0
                          thetaStar_PROP[m](l) = 0;
                        }else{                      // otherwise draw from Q1 (which is the prior)
                          thetaStar_PROP[m](l) = rgamma(1,a_slabpos,1/b_slabpos)(0); //takes scale not rate
                        }
                        // update components that depend on thetaStar
                        XthetaStar.col(m) = X[m] * thetaStar_PROP[m];                          // update XthetaStar with proposed thetaStar
                        
                        // update Sigma with proposed thetaStar
                        Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat);
                        
                        // proposed lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
                        lltheta_PROP = log(tgamma(sumdelta_PROP+a_pi))+log(tgamma(sumM1-sumdelta_PROP+b_pi))
                          //+ delta_PROP[m](l)*R::dgamma(thetaStar_PROP[m](l),a_slabpos,1/b_slabpos,TRUE)//+(1-delta[m][l])*log(1)//-logLambdaInverse*logLambdaInverse/(2*b_lambda)
                            - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                            - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
                        if(delta_PROP[m](l)==1){ //adding component (avoiding NAs, because dgamma at 0 can give Inf)
                          lltheta_PROP += delta_PROP[m](l)*R::dgamma(thetaStar_PROP[m](l),a_slabpos,1/b_slabpos,TRUE);
                        }   
                        
                        // difference between proposal distributions: logP(theta_prop,delta_prop|theta,delta)-logP(theta,delta|theta_prop,delta_prop)
                        logdiff_PROP = //delta_PROP[m](l)*R::dgamma(thetaStar_PROP[m](l),a_slabpos,1/b_slabpos,TRUE) //   logP(theta_prop,delta_prop|theta,delta)
                          //- delta[m](l)*R::dgamma(thetaStar[m](l),a_slabpos,1/b_slabpos,TRUE) // logP(theta,delta|theta_prop,delta_prop)
                            +((sumdelta_PROP==0)*log(1)+(sumdelta_PROP!=0)*log(0.5)) // extra difference component because sometimes we must choose move 1.
                            -((sumdelta1==0)*log(1)+(sumdelta1!=0)*log(0.5));
                      
                        if(delta_PROP[m](l)==1){ //adding component to logdiff (avoiding NAs, because dgamma at 0 can give Inf)
                          logdiff_PROP += delta_PROP[m](l)*R::dgamma(thetaStar_PROP[m](l),a_slabpos,1/b_slabpos,TRUE);
                        }else{
                          logdiff_PROP -= delta[m](l)*R::dgamma(thetaStar[m](l),a_slabpos,1/b_slabpos,TRUE);
                        }  
                        
                        
                        // log ratio for MH acceptance
                        logratio = (lltheta_PROP-lltheta)-logdiff_PROP; // the extra (negative!) difference in log proposals since this is MH not metropolis
                        if(logratio>0){
                          logratio = 0;
                        }
                          
                        // accept or reject
                        if(logratio >  log(runif(1)(0)) ){
                          delta[m](l) = delta_PROP[m](l);
                          thetaStar[m](l) = thetaStar_PROP[m](l);
                          // ll = ll_PROP;
                        }else{
                          XthetaStar.col(m) = X[m] * thetaStar[m];
                        }
                                  
                                  
                      } // end (if move_id)
                      
                      
                      
                    } // l loop
                  } // end if (thetaconstraint)
                } // m loop
                
                
                
                
                
              }else{ // MHmove version 2
                
                
                move_which = sample(sumdelta1, 1)(0);    // randomly select component for move (such that delta==1)
                
                move_id = 0;                                // set index
                for(int m=0; m<M; m++){
                  if(thetaconstraint(m)==1){     // only loop through components with positivity constraints
                    for(int l=0; l<Lvec(m); l++){
                      
                      
                      if(delta[m](l)==1){
                        move_id += 1;
                        if(move_id==move_which){   //RANDOM WALK STEP
                          
                          
                          // compute lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
                          // first term is from change of vars  
                          lltheta = log(thetaStar[m](l))+R::dgamma(thetaStar[m](l),a_slabpos,1/b_slabpos,TRUE) //-logLambdaInverse*logLambdaInverse/(2*b_lambda)
                          - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                          - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
                          
                          // propose new theta
                          thetaStar_PROP[m](l) = exp(log(thetaStar[m](l)) + step_theta*rnorm(1)(0));
                          
                          // update components that depend on thetaStar
                          XthetaStar.col(m) = X[m] * thetaStar_PROP[m];                          // update XthetaStar with proposed thetaStar
                          
                          // update Sigma with proposed thetaStar
                          Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat);
                          
                          // proposed lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
                          lltheta_PROP = log(thetaStar_PROP[m](l))+R::dgamma(thetaStar_PROP[m](l),a_slabpos,1/b_slabpos,TRUE) //-logLambdaInverse*logLambdaInverse/(2*b_lambda)
                            - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                            - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
                            
                            
                            // log ratio for MH acceptance
                            logratio = lltheta_PROP-lltheta;
                            if(logratio>0){
                              logratio = 0;
                            }
                            
                            // accept or reject
                            if(logratio >  log(runif(1)(0)) ){
                              thetaStar[m](l) = thetaStar_PROP[m](l);
                              // ll = ll_PROP;
                            }else{
                              XthetaStar.col(m) = X[m] * thetaStar[m];
                            }
                            
                            
                            
                        } // end (if move_id)
                      } // end (if delta=1)
                    } // l loop
                    
                    
                  } // end (if thetaconstraint) 
                } // m loop
              } // end move 2
              
            } // end positive constraint update
            
            
            
            
            
            // -------------------------------------------------------------------------- //
            // UPDATE delta and thetaStar via Metropolis Hastings NO CONSTRAINTS
            
            if(sumM0>1){
              // count number of delta==1
              sumdelta0 = 0;
              for(int m=0; m<M; m++){
                // only loop through components with no constraints
                if(thetaconstraint(m)==0){
                  sumdelta0 += sum(delta[m]);
                }
              }
              
              // select move type
              if(sumdelta0==0){ // must choose move 1 if all delta==0 (move 2 samples from the deltas equal to 1)
                MHmove = 1;
              }else{ // randomly select move type
                MHmove = Rcpp::rbinom(1,1,0.5)(0);
              }
              
              // Make move 1 or 2
              if(MHmove==1){  // MHmove version 1
                
                // reset proposals
                delta_PROP = delta;
                thetaStar_PROP = thetaStar;
                
                // randomly select a component
                move_which = sample(sumM0, 1)(0);
                
                move_id = 0;                                // set index
                for(int m=0; m<M; m++){
                  if(thetaconstraint(m)==0){// only loop through components with no constraints
                    for(int l=0; l<Lvec(m); l++){
                      
                      move_id += 1;
                      if(move_id==move_which){   // make move for randomly selected component
                        
                        // compute lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
                        lltheta = log(tgamma(sumdelta0+a_pi))+log(tgamma(sumM0-sumdelta0+b_pi))
                        + delta[m](l)*R::dnorm(thetaStar[m](l),0.0,s_theta,TRUE)//+(1-delta[m][l])*log(1)//-logLambdaInverse*logLambdaInverse/(2*b_lambda)
                          - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                          - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
                          
                        // propose new delta: 1 to 0 or 0 to 1
                        delta_PROP[m](l) = 1-delta[m](l);
                        if(delta_PROP[m](l)==1){ // adjust sum delta for switch (either increase by 1 or decrease by 1)
                          sumdelta_PROP = sumdelta0+1;
                        }else{
                          sumdelta_PROP = sumdelta0-1;
                        }
                        
                        //propose new theta* (based on new delta)
                        if(delta_PROP[m](l)==0){    // if delta_PROP, then theta*=0
                          thetaStar_PROP[m](l) = 0;
                        }else{                      // otherwise draw from Q1 (which is the prior)
                          thetaStar_PROP[m](l) = s_theta*rnorm(1)(0);
                        }
                        
                        // update components that depend on thetaStar
                        XthetaStar.col(m) = X[m] * thetaStar_PROP[m];                          // update XthetaStar with proposed thetaStar
                        
                        // update Sigma with proposed thetaStar
                        Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat);
                        
                        // proposed lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
                        lltheta_PROP = log(tgamma(sumdelta_PROP+a_pi))+log(tgamma(sumM0-sumdelta_PROP+b_pi))
                          + delta_PROP[m](l)*R::dnorm(thetaStar_PROP[m](l),0.0,s_theta,TRUE)//+(1-delta[m][l])*log(1)//-logLambdaInverse*logLambdaInverse/(2*b_lambda)
                            - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                            - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
                            
                        // difference between proposal distributions: logP(theta_prop,delta_prop|theta,delta)-logP(theta,delta|theta_prop,delta_prop)
                        logdiff_PROP = delta_PROP[m](l)*R::dnorm(thetaStar_PROP[m](l),0.0,s_theta,TRUE) //   logP(theta_prop,delta_prop|theta,delta)
                          - delta[m](l)*R::dnorm(thetaStar[m](l),0.0,s_theta,TRUE) // logP(theta,delta|theta_prop,delta_prop)
                            +((sumdelta_PROP==0)*log(1)+(sumdelta_PROP!=0)*log(0.5)) // extra difference component because sometimes we must choose move 1.
                            -((sumdelta0==0)*log(1)+(sumdelta0!=0)*log(0.5));
                            
                        // log ratio for MH acceptance
                        logratio = (lltheta_PROP-lltheta)-logdiff_PROP; // the extra (negative!) difference in log proposals since this is MH not metropolis
                        if(logratio>0){
                          logratio = 0;
                        }
                        
                        // accept or reject
                        if(logratio >  log(runif(1)(0)) ){
                          delta[m](l) = delta_PROP[m](l);
                          thetaStar[m](l) = thetaStar_PROP[m](l);
                          // ll = ll_PROP;
                        }else{
                          XthetaStar.col(m) = X[m] * thetaStar[m];
                        }
                                
                                
                      } // end (if move_id)
                      
                      
                      
                    } // l loop
                  } // end if (thetaconstraint)
                } // m loop
                
                
                
                
                
              }else{ // MHmove version 2
                
                
                move_which = sample(sumdelta0, 1)(0);    // randomly select component for move (such that delta==1)
                
                move_id = 0;                                // set index
                for(int m=0; m<M; m++){
                  if(thetaconstraint(m)==0){// only loop through components with no constraints
                    for(int l=0; l<Lvec(m); l++){
                      
                      if(delta[m](l)==1){
                        move_id += 1;
                        if(move_id==move_which){   //RANDOM WALK STEP
                          
                          
                          // compute lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
                          lltheta = R::dnorm(thetaStar[m](l),0.0,s_theta,TRUE) //-logLambdaInverse*logLambdaInverse/(2*b_lambda)
                          - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                          - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
                          
                          // propose new theta
                          thetaStar_PROP[m](l) = thetaStar[m](l) + step_theta*rnorm(1)(0);
                          
                          // update components that depend on thetaStar
                          XthetaStar.col(m) = X[m] * thetaStar_PROP[m];                          // update XthetaStar with proposed thetaStar
                          
                          // update Sigma with proposed thetaStar
                          Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat);
                          
                          // proposed lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
                          lltheta_PROP = R::dnorm(thetaStar_PROP[m](l),0.0,s_theta,TRUE) //-logLambdaInverse*logLambdaInverse/(2*b_lambda)
                            - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
                            - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
                            
                            
                          // log ratio for MH acceptance
                          logratio = lltheta_PROP-lltheta;
                          if(logratio>0){
                            logratio = 0;
                          }
                          
                          // accept or reject
                          if(logratio >  log(runif(1)(0)) ){
                            thetaStar[m](l) = thetaStar_PROP[m](l);
                            // ll = ll_PROP;
                          }else{
                            XthetaStar.col(m) = X[m] * thetaStar[m];
                          }
                            
                            
                            
                        } // end (if move_id)
                      } // end (if delta=1)
                      
                      
                      
                    } // l loop
                  } // end if (thetaconstraint)
                } // m loop
              } // end move 2
            } // end delta/theta update (no constraints)
            
            // Reset matrix Sigma and ll
            Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat); 
            
            ll = logLambdaInverse+R::dgamma(exp(logLambdaInverse),a_lam,1/b_lam,1) //-logLambdaInverse*logLambdaInverse/(2*b_lambda) 
              - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
              - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
              
              
              
              
              
      } // end inner loop
      
      
      
      
      // save posterior samples
      for(int m=0; m<M; m++){
        // partition thetaStar* into rho and theta
        rhokeep.submat(s_outer,m,s_outer,m) = sum(square(vectorise(thetaStar[m])));                             // rho=||theta*||^2
        thetakeepfield[m].row(s_outer) =  thetaStar[m].t()/sqrt(as_scalar(rhokeep.submat(s_outer,m,s_outer,m)));// theta=||theta||^{-1}theta*
        
        // get dirichlet weights (i.e. use a different identifiability constraint when the thetastar are constrained to be positive; thetaPOS sum to 1) 
        if(thetaconstraint(m)>0){
          thetaPOSkeepfield[m].row(s_outer) =  thetaStar[m].t()/sum(vectorise(thetaStar[m]));// theta=||theta||^{-1}theta*
        }
        
      }
      logLambdaInversekeep(s_outer) = logLambdaInverse;       // save lambda^{-1}  
      logLambdaBInversekeep(s_outer) = logLambdaBInverse;     // save lambda^{-1}  
      
      
      // draw sigma^{-2} and gamma directly and save them
      sig2invkeep(s_outer) = rgamma( 1 , a_sig + .5*(N-P_z) , 1/(b_sig + 0.5*(Vcomps[1](0,0)-Vcomps[1](1,0))) )(0);  // function generates a gamma distribution parameterized with mean alpha*beta so we just inverted the second component;
      // gammakeep.row(s_outer) = arma::solve(trimatu(C2) , C2B +  as<arma::vec>(rnorm(P_z))/sqrt(sig2invkeep(s_outer) )).t(); // draw gamma by drawing MVN(0,I) then scaling by sqrt of variance and shifting by mean
      SIG = inv(Vcomps[4]*Vcomps[4].t());
      
      gammakeep.row(s_outer) = (SIG*Vcomps[2]+chol(SIG,"lower")*as<arma::vec>(rnorm(P_z))/sqrt(sig2invkeep(s_outer) )).t();
      
      // compute fitted h after burn-in
      if(s_outer>=n_burn){
        laminvK = Vcomps[5];      
        laminvKSiginv = solve(Vcomps[0] , laminvK).t();                                     // lambda^{-1} K Sigma^{-1} which is equal to t(Sigma^{-1} lambda^{-1} K) 
        hmeantemp = laminvKSiginv * (yz.col(0)-yz.cols(1,P_z)*gammakeep.row(s_outer).t());  // temporary vector to be used in mean and covariance
        // hcovtemp = laminvKSiginv*(Vcomps[0]-laminvK)/sig2invkeep(s_outer);               // temporary covariance matrix to be used in marginal covariance
        hcovtemp = (laminvK-laminvKSiginv*laminvK.t())/sig2invkeep(s_outer);                // temporary covariance matrix to be used in marginal covariance
        // compute marginal mean and covariance for h
        hmean += hmeantemp/(n_outer-n_burn);                                                                    // law of total expecation, E[E[h|s]]
        hcov += hcovtemp/(n_outer-n_burn) + hmeantemp*hmeantemp.t()/(n_outer-n_burn); // first two components of law of total variance: E(Var(h|s))+E[E[h|s]^2]-E[E[h|s]]^2
        // old (no B) : hcov += laminvKSiginv/(n_outer-n_burn)/sig2invkeep(s_outer) + hmeantemp*hmeantemp.t()/(n_outer-n_burn); // first two components of law of total variance: E(Var(h|s))+E[E[h|s]^2]-E[E[h|s]]^2
        // hcov += laminvKSiginv*(Vcomps[0]-laminvK)/(n_outer-n_burn)/sig2invkeep(s_outer) + hmeantemp*hmeantemp.t()/(n_outer-n_burn); // first two components of law of total variance: E(Var(h|s))+E[E[h|s]^2]-E[E[h|s]]^2
        if(draw_h==TRUE){
          // note: "std" method used to avoid errors---to return to default, remove "std"
          arma::svd(svdU,svdD,svdV,hcovtemp,"std"); //arma::svd(svdU,svdD,svdV,hcovtemp);                                                 // SVD decomposition for covariance of h (since it's more stable than cholesky)
          // draw h
          hsampkeep.row(s_outer) = ( hmeantemp+(svdU*diagmat(svdD))*arma::randn(N) ).t();     // hsampkeep.row(s_outer) = ( hmeantemp+arma::chol(hcovtemp)*arma::randn(N) ).t();
        }
      }
      
    } // end outer loop
    
    
    
    
    // Convert field to list
    for(int m=0; m<M; m++){
      thetakeep[m] = thetakeepfield[m];
      thetaPOSkeep[m] = thetaPOSkeepfield[m];
    }
    
    // return a list
    return List::create(Named("theta")          = thetakeep,
                        Named("thetaPOS")       = thetaPOSkeep,                 // weights defined differently when theta constrained to be positive
                        Named("lambdaInverse")  = exp(logLambdaInversekeep),
                        Named("lambdaBInverse") = exp(logLambdaBInversekeep),
                        Named("rho")            = rhokeep,
                        Named("sigma2")         = 1/sig2invkeep,
                        Named("gamma")          = gammakeep,
                        Named("hsamp")          = hsampkeep,
                        Named("hmean")          = hmean,
                        Named("hcov")           = hcov - ( hmean*hmean.t())); // subtracting third component of law of total variance (i.e. E[E[h|s]]^2)
}



// Normal prior version
// // [[Rcpp::depends(RcppArmadillo)]]
// //' MCMC for BSMIM with spike and slab prior
// //' 
// //' Returns posterior draws for \eqn{\theta},\eqn{\lambda^{-1}},\eqn{\rho},\eqn{\sigma^2},\eqn{\gamma}, as well as the mean(\eqn{h}) and cov(\eqn{h})
// //' Adapted from "bkmrdlm_multi.cpp" in the "regimes" package by Ander Wilson.
// //' Note that this implementation with horseshoe2 prior returns posterior draws of \eqn{\rho} as a list of M matrices (unlike main MCMC code, which returns one matrix of M columns)
// //' 
// //' @param yz matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
// //' @param Xlist list of (N by L_m) matrices representing X_m
// //' @param b_lambda hyperparameter for \eqn{\lambda^{-1}}
// //' @param b_lambdaB hyperparameter for \eqn{\lambda^{-1}_B}
// //' @param a_sig first hyperparameter for \eqn{\sigma^{-2}}
// //' @param b_sig second hyperparameter for \eqn{\sigma^{-2}}
// //' @param tau02 hyperparameter for tau2; classic horseshoe is 1
// //' @param kappa vector of hyperparameters \eqn{\kappa_m} for \eqn{\theta^*_m}
// //' @param poly 0 = gaussian kernel, 1 = polynomial kernel
// //' @param d degree of polynomial kernel
// //' @param horseshoe 0 = no selection, 1 = componentwise horseshoe priors
// //' @param randint 0 = no random intercepts, 1 = random intercepts model
// //' @param Bmat N xN block diagonal matrix indicating cluster membership for random intercepts model
// //' @param draw_h 0 = dont draw h, 1 = draw h
// //' @param n_inner no. of MCMC iterations to run in the inner loop. n_outer*n_inner iteraction will be run.
// //' @param n_outer no. of MCMC iterations to run in the outer loop. n_outer iterations will be saved.
// //' @param n_burn no. of MCMC iterations to discard as burn-in
// //' @author Glen McGee and Ander Wilson (adapted from the "regimes" package by Ander Wilson).
// //' @export
// // [[Rcpp::export]]
// List bsmim_spikeslab_mcmc2(const arma::mat&    yz,         // matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
//                            const Rcpp::List&   Xlist,      // list of (N by L_m) matrices representing X_m
//                            const double&       b_lambda,   // hyperparameter for \lambda^{-1}
//                            const double&       b_lambdaB,  // hyperparameter for \lambda^{-1}_B
//                            const double&       a_sig,      // first hyperparameter for \sigma^{-2}
//                            const double&       b_sig,      // second hyperparameter for \sigma^{-2}
//                            const double&       s_theta,    // hyperparameter for sd of slab component
//                            const double&       step_theta, // step size for theta* random walk (move 2)
//                            const double&       a_pi,       // first hyperparameter for beta distribution of pi 
//                            const double&       b_pi,       // second hyperparameter for beta distribution of pi 
//                            const bool&         poly,       // 0=gaussian kernel / 1=polynomial kernel
//                            const int&          d,          // degree of polynomial kernel
//                            const bool&         randint,    // 0=no random intercepts / 1=random intercepts
//                            const arma::mat&    Bmat,       // block diagonal matrix B indicating cluster membership for random intercepts model
//                            const bool&         draw_h,     // 0=dont sample h, 1= sample h
//                            const int&          n_inner,    // no. of MCMC iterations to run in the inner loop. n_outer*n_inner iteraction will be run.
//                            const int&          n_outer,    // no. of MCMC iterations to run in the outer loop. n_outer iterations will be saved.
//                            const int&          n_burn) {   // no. of MCMC iterations to discard as burn-in
//   
//   // step size for random walk in move 2 of MH step (spike/slab)
//   // double step_theta = 0.2;
//   // TRYING DIFFERENT LAMBDA PRIOR
//   double a_lam = 1.0;
//   double b_lam = 0.1;
//   // Dimensions
//   int N = yz.n_rows;              // no. of observations
//   int M = Xlist.size();           // no. of indices (dimension of h(.))
//   int P_z = yz.n_cols-1;          // no. of covariates
//   IntegerVector Lvec = rep(0,M);  // vector of L_m representing the no. of columns of X_m
//   
//   // Initialize
//   arma::field<arma::mat> X(M);                    // field of (N by L_m) matrices representing X_m
//   arma::field<arma::vec> thetaStar(M);            // field of L_m vectors representing \thetaStar_m
//   arma::field<arma::vec> delta(M);                // field of L_m vectors representing \delta_m
//   arma::mat  XthetaStar(N,M,arma::fill::zeros);   // (N by M) matrix representing [X_1\thetaStar_1,...,X_M\thetaStar_M]
//   for(int m=0; m<M; m++){
//     X[m] = as<arma::mat>(Xlist[m]);                 // mth matrix in X is X_m
//     Lvec(m) = X[m].n_cols;                          // mth element of Lvec is L_m (the length of \thetaStar_m)
//     thetaStar[m] = as<arma::vec>(rnorm(Lvec(m)));   // set starting value for \thetaStar_m
//     XthetaStar.col(m) = X[m] * thetaStar[m];
//     
//     arma::vec deltatemp(Lvec(m),arma::fill::ones); // temporary vec with L_m elements to store posteror draws of delta
//     delta[m] = deltatemp;                          // set starting value for \delta_m
//   }
// 
//   
//   double logLambdaInverse = rnorm(1)(0)*sqrt(b_lambda)/10;    // set starting value for log(\lambda^{-1})
//   double logLambdaBInverse = 0;                               // set log(\lambda^{-1}_B) to 0 if no random intercepts
//   if(randint==TRUE){
//     logLambdaBInverse = rnorm(1)(0)*sqrt(b_lambdaB)/10;       // set starting value for log(\lambda^{-1}_B) if using random intercepts
//   }
// 
//   
//   
//   // Construct Sigma=I+\lambda^{-1}*K, A1, A2, C1, C2
//   arma::field<arma::mat> Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat);
//   
//   // First evaluation of the log-liklihood (ll)
//   double ll = logLambdaInverse+R::dgamma(exp(logLambdaInverse),a_lam,1/b_lam,1) //-logLambdaInverse*logLambdaInverse/(2*b_lambda) 
//     - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
//     - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
//     
//     
//     // Initialize other values
//     double logLambdaInverse_PROP, llu, ll_PROP, logratio;
//     double llB, logLambdaBInverse_PROP, llB_PROP;
//     arma::vec ons(N, arma::fill::ones);           // vector of length N of ones
//     arma::vec vec_PROP;                           // proposal vector to define ellipse
//     
//     int MHmove, sumdelta, sumdelta_PROP, move_which, move_id;    // move type for joint MH draw of delta and theta* (1 or 2); sum of delta=1; which component MHmove applies to; id to sum
//     double lltheta, lltheta_PROP, logdiff_PROP;
//     arma::field<arma::vec> thetaStar_PROP = thetaStar;
//     arma::field<arma::vec> delta_PROP = delta;
// 
//     
//     // Store posterior samples
//     arma::vec logLambdaInversekeep(n_outer,arma::fill::zeros);  // vector to store posterior draws of lambda^{-1}
//     arma::vec logLambdaBInversekeep(n_outer,arma::fill::zeros); // vector to store posterior draws of lambda_B^{-1}
//     arma::mat rhokeep(n_outer,M,arma::fill::zeros);             // matrix to store posterior draws of rho // each row is rho_1,...,rho_m
//     arma::vec sig2invkeep(n_outer,arma::fill::zeros);           // vector to store posterior draws of sigma^{-2}
//     arma::mat gammakeep(n_outer,P_z,arma::fill::zeros);         // matrix to store posterior draws of gamma // each row is gamma_1,...,gamma_{P_z}
//     arma::field<arma::mat>  thetakeepfield(M);                  // field of matrices to store posterior draws of theta (for computation)
//     Rcpp::List  thetakeep(M);                                   // List of matrices to store posteror draws of theta (same as above but for exporting to R)
//     for(int m=0; m<M; m++){
//       arma::mat thetakeeptemp(n_outer,Lvec(m),arma::fill::zeros); // temporary matrix with L_m columns to store posteror draws of theta_m
//       thetakeepfield[m] = thetakeeptemp;                          // initialize elements of thetakeepfield to be the right size
//       thetakeep[m] = thetakeeptemp;                               // initialize elements of thetakeepfield to be the right size
// 
//     } 
//     
//     // Initialize values for posterior of h(.)
//     arma::mat laminvK(N,N,arma::fill::zeros);           // lambda^{-1}K where K is the kernel matrix
//     arma::mat laminvKSiginv(N,N,arma::fill::zeros);     // \lambda^{-1}KSigma^{-1}=laminvK*Sigma^{-1} or t(solve(Sigma)*laminvK)
//     arma::vec hmean(N,arma::fill::zeros);               // marginal mean vector for h(.)
//     arma::mat hcov(N,N,arma::fill::zeros);              // marginal variance-covariance matrix of h(.)
//     arma::vec hmeantemp(N,arma::fill::zeros);           // temporary vector for computing marginal covariance of h(.)
//     arma::mat hcovtemp(N,N,arma::fill::zeros);          // marginal variance-covariance matrix of h(.)
//     arma::mat SIG(P_z,P_z,arma::fill::zeros);           // temporary matrix
//     arma::mat hsampkeep(n_outer,N,arma::fill::zeros);   // matrix to store posterior draws of h: rows are samples of h
//     arma::mat svdU;                                     // for SVD decomposition
//     arma::vec svdD;                                     // for SVD decomposition
//     arma::mat svdV;                                     // for SVD decomposition
//     
//     // -------------------------------------
//     // MCMC
//     // -------------------------------------
//     for(int s_outer=0; s_outer<n_outer; s_outer++){                       // outer loop
//       
//       if(s_outer % 10 == 0) Rcout << "#" << s_outer << std::endl;         // output progress every 10 iterations
//       
//       for(int s_inner=0; s_inner<n_inner; s_inner++){                     // inner loop
//         
//         
//         // -------------------------------------------------------------------------- //
//         // UPDATE LAMBDA^{-1} via random walk metropolis
//         logLambdaInverse_PROP = logLambdaInverse + rnorm(1)(0); // random walk step
//         
//         // update Sigma=I+\lambda^{-1}*K, A1, A2, C1, C2
//         Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse_PROP, XthetaStar, randint, logLambdaBInverse, Bmat);   
//         
//         // proposed ll 
//         ll_PROP = logLambdaInverse_PROP+R::dgamma(exp(logLambdaInverse_PROP),a_lam,1/b_lam,1) //-logLambdaInverse*logLambdaInverse/(2*b_lambda) -logLambdaInverse_PROP*logLambdaInverse_PROP/(2*b_lambda)
//           - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
//           - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
//           
//           // log ratio for MH acceptance
//           logratio = ll_PROP-ll;
//           if(logratio>0){
//             logratio = 0;
//           }
//           
//           // accept or reject    
//           if(logratio >  log(runif(1)(0)) ){                
//             logLambdaInverse = logLambdaInverse_PROP; 
//             ll = ll_PROP;
//           }
//           
//           // Reset matrix Sigma and ll
//           Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat); 
//           
//           ll = logLambdaInverse+R::dgamma(exp(logLambdaInverse),a_lam,1/b_lam,1) //-logLambdaInverse*logLambdaInverse/(2*b_lambda)
//             - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
//             - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
//             
//             
//             
//             
//             // -------------------------------------------------------------------------- //
//             // UPDATE LAMBDA_B^{-1} via random walk metropolis
//             if(randint==TRUE){ // only update lambda_B^{-1} if using random intercepts model
//               
//               logLambdaBInverse_PROP = logLambdaBInverse + rnorm(1)(0); // random walk step
//               
//               // current llB (uses the prior for logLambdaBInverse, hence the first line below is different than ll elsewhere)
//               llB = logLambdaBInverse+R::dgamma(exp(logLambdaBInverse),a_lam,1/b_lam,1) //-logLambdaBInverse*logLambdaBInverse/(2*b_lambdaB) 
//                 - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
//                 - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
//                 
//                 // update Sigma=I+\lambda^{-1}*K, A1, A2, C1, C2
//                 Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse_PROP, Bmat);   
//                 
//                 // proposed llB (uses the prior for logLambdaBInverse, hence the first line below is different than ll elsewhere)
//                 llB_PROP = logLambdaBInverse_PROP+R::dgamma(exp(logLambdaBInverse_PROP),a_lam,1/b_lam,1) //-logLambdaBInverse_PROP*logLambdaBInverse_PROP/(2*b_lambdaB)
//                   - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
//                   - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) );
//                   
//                   // log ratio for MH acceptance
//                   logratio = llB_PROP-llB;
//                   if(logratio>0){
//                     logratio = 0;
//                   }
//                   
//                   // accept or reject    
//                   if(logratio >  log(runif(1)(0)) ){                
//                     logLambdaBInverse = logLambdaBInverse_PROP; 
//                   }
//                   
//                   // Reset matrix Sigma and ll
//                   Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat); 
//                   
//                   // resetting ll (not llB)
//                   ll = logLambdaInverse+R::dgamma(exp(logLambdaInverse),a_lam,1/b_lam,1) //-logLambdaInverse*logLambdaInverse/(2*b_lambda) 
//                     - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
//                     - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
//                     
//             }
//             
//             
//             
//             
//             // -------------------------------------------------------------------------- //
//             // UPDATE delta and thetaStar via Metropolis Hastings
//             
//             // count number of delta==1
//             sumdelta = 0;
//             for(int m=0; m<M; m++){
//               sumdelta += sum(delta[m]);  
//             }
//             
//             // select move type
//             if(sumdelta==0){ // must choose move 1 if all delta==0 (move 2 samples from the deltas equal to 1)
//               MHmove = 1;
//             }else{ // randomly select move type
//               MHmove = Rcpp::rbinom(1,1,0.5)(0);
//             }
//             
//             // Make move 1 or 2
//             if(MHmove==1){  // MHmove version 1
//               
//               // reset proposals
//               delta_PROP = delta;
//               thetaStar_PROP = thetaStar;
//               
//               // randomly select a component
//               move_which = sample(sum(Lvec), 1)(0);
// 
//               move_id = 0;                                // set index
//               for(int m=0; m<M; m++){
//                 for(int l=0; l<Lvec(m); l++){
//                   
//                     move_id += 1;
//                     if(move_id==move_which){   // make move for randomly selected component
//                       
//                       // compute lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
//                       lltheta = log(tgamma(sumdelta+a_pi))+log(tgamma(sum(Lvec)-sumdelta+b_pi))
//                       + delta[m](l)*R::dnorm(thetaStar[m](l),0.0,s_theta,TRUE)//+(1-delta[m][l])*log(1)//-logLambdaInverse*logLambdaInverse/(2*b_lambda) 
//                       - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
//                       - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
//                       
//                       // propose new delta: 1 to 0 or 0 to 1
//                       delta_PROP[m](l) = 1-delta[m](l);
//                       if(delta_PROP[m](l)==1){ // adjust sum delta for switch (either increase by 1 or decrease by 1)
//                         sumdelta_PROP = sumdelta+1;
//                       }else{
//                         sumdelta_PROP = sumdelta-1;
//                       }
//                       
//                       //propose new theta* (based on new delta)
//                       if(delta_PROP[m](l)==0){    // if delta_PROP, then theta*=0
//                         thetaStar_PROP[m](l) = 0;
//                       }else{                      // otherwise draw from Q1 (which is the prior)
//                         thetaStar_PROP[m](l) = s_theta*rnorm(1)(0);
//                       }
//                       
//                       // update components that depend on thetaStar
//                       XthetaStar.col(m) = X[m] * thetaStar_PROP[m];                          // update XthetaStar with proposed thetaStar
//                       
//                       // update Sigma with proposed thetaStar
//                       Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat); 
//                       
//                       // proposed lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
//                       lltheta_PROP = log(tgamma(sumdelta_PROP+a_pi))+log(tgamma(sum(Lvec)-sumdelta_PROP+b_pi))
//                         + delta_PROP[m](l)*R::dnorm(thetaStar_PROP[m](l),0.0,s_theta,TRUE)//+(1-delta[m][l])*log(1)//-logLambdaInverse*logLambdaInverse/(2*b_lambda) 
//                         - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
//                         - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
//                       
//                       // difference between proposal distributions: logP(theta_prop,delta_prop|theta,delta)-logP(theta,delta|theta_prop,delta_prop)
//                       logdiff_PROP = delta_PROP[m](l)*R::dnorm(thetaStar_PROP[m](l),0.0,s_theta,TRUE) //   logP(theta_prop,delta_prop|theta,delta)
//                         - delta[m](l)*R::dnorm(thetaStar[m](l),0.0,s_theta,TRUE) // logP(theta,delta|theta_prop,delta_prop)
//                         +((sumdelta_PROP==0)*log(1)+(sumdelta_PROP!=0)*log(0.5)) // extra difference component because sometimes we must choose move 1.
//                         -((sumdelta==0)*log(1)+(sumdelta!=0)*log(0.5));
//                         
//                       // log ratio for MH acceptance
//                       logratio = (lltheta_PROP-lltheta)-logdiff_PROP; // the extra (negative!) difference in log proposals since this is MH not metropolis
//                       if(logratio>0){
//                         logratio = 0;
//                       }
//                         
//                       // accept or reject    
//                       if(logratio >  log(runif(1)(0)) ){    
//                         delta[m](l) = delta_PROP[m](l);
//                         thetaStar[m](l) = thetaStar_PROP[m](l);
//                         // ll = ll_PROP;
//                       }else{
//                         XthetaStar.col(m) = X[m] * thetaStar[m];
//                       }
//                         
//                         
//                     } // end (if move_id)
// 
//                 
//                   
//                 } // l loop
//               } // m loop
//               
//               
//               
//  
//               
//             }else{ // MHmove version 2
//               
// 
//               move_which = sample(sumdelta, 1)(0);    // randomly select component for move (such that delta==1)
//               
//               move_id = 0;                                // set index
//               for(int m=0; m<M; m++){
//                 for(int l=0; l<Lvec(m); l++){
//                   
//                   if(delta[m](l)==1){
//                     move_id += 1;
//                     if(move_id==move_which){   //RANDOM WALK STEP
//                       
//                       
//                       // compute lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
//                       lltheta = R::dnorm(thetaStar[m](l),0.0,s_theta,TRUE) //-logLambdaInverse*logLambdaInverse/(2*b_lambda) 
//                         - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
//                         - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
//                       
//                       // propose new theta
//                       thetaStar_PROP[m](l) = thetaStar[m](l) + step_theta*rnorm(1)(0);
//                       
//                       // update components that depend on thetaStar
//                       XthetaStar.col(m) = X[m] * thetaStar_PROP[m];                          // update XthetaStar with proposed thetaStar
//                       
//                       // update Sigma with proposed thetaStar
//                       Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat); 
//                       
//                       // proposed lltheta (excludes the lambda prior, and the other thetas/deltas since they are unchanged)
//                       lltheta_PROP = R::dnorm(thetaStar_PROP[m](l),0.0,s_theta,TRUE) //-logLambdaInverse*logLambdaInverse/(2*b_lambda) 
//                         - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
//                         - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
//                         
//                       
//                       // log ratio for MH acceptance
//                       logratio = lltheta_PROP-lltheta;
//                       if(logratio>0){
//                         logratio = 0;
//                       }
//                       
//                       // accept or reject    
//                       if(logratio >  log(runif(1)(0)) ){       
//                         thetaStar[m](l) = thetaStar_PROP[m](l);
//                         // ll = ll_PROP;
//                       }else{
//                         XthetaStar.col(m) = X[m] * thetaStar[m];
//                       }
//                       
//                         
//   
//                     } // end (if move_id)
//                   } // end (if delta=1)
//                   
//               
//                   
//                 } // l loop
//               } // m loop
//             } // end move 2
//             
//             // Reset matrix Sigma and ll
//             Vcomps = get_Vcomps(poly, d, N, P_z, yz, logLambdaInverse, XthetaStar, randint, logLambdaBInverse, Bmat); 
//             
//             ll = logLambdaInverse+R::dgamma(exp(logLambdaInverse),a_lam,1/b_lam,1) //-logLambdaInverse*logLambdaInverse/(2*b_lambda) 
//               - arma::accu(log(Vcomps[3].diag())) - arma::accu(log(Vcomps[4].diag()))       // computing determinants via cholesky decompositions
//               - ( 0.5 * (N - P_z) + a_sig ) * log( b_sig + 0.5*Vcomps[1](0,0)-0.5*Vcomps[1](1,0) ); 
//               
//             
//             
//           
//             
//       } // end inner loop
//       
//       
//       
//       
//       // save posterior samples
//       for(int m=0; m<M; m++){
//         // partition thetaStar* into rho and theta
//         rhokeep.submat(s_outer,m,s_outer,m) = sum(square(vectorise(thetaStar[m])));                             // rho=||theta*||^2
//         thetakeepfield[m].row(s_outer) =  thetaStar[m].t()/sqrt(as_scalar(rhokeep.submat(s_outer,m,s_outer,m)));// theta=||theta||^{-1}theta*
//         
// 
//       }
//       logLambdaInversekeep(s_outer) = logLambdaInverse;       // save lambda^{-1}  
//       logLambdaBInversekeep(s_outer) = logLambdaBInverse;     // save lambda^{-1}  
// 
//       
//       // draw sigma^{-2} and gamma directly and save them
//       sig2invkeep(s_outer) = rgamma( 1 , a_sig + .5*(N-P_z) , 1/(b_sig + 0.5*(Vcomps[1](0,0)-Vcomps[1](1,0))) )(0);  // function generates a gamma distribution parameterized with mean alpha*beta so we just inverted the second component;
//       // gammakeep.row(s_outer) = arma::solve(trimatu(C2) , C2B +  as<arma::vec>(rnorm(P_z))/sqrt(sig2invkeep(s_outer) )).t(); // draw gamma by drawing MVN(0,I) then scaling by sqrt of variance and shifting by mean
//       SIG = inv(Vcomps[4]*Vcomps[4].t());
//       
//       gammakeep.row(s_outer) = (SIG*Vcomps[2]+chol(SIG,"lower")*as<arma::vec>(rnorm(P_z))/sqrt(sig2invkeep(s_outer) )).t();
//       
//       // compute fitted h after burn-in
//       if(s_outer>=n_burn){
//         laminvK = Vcomps[5];      
//         laminvKSiginv = solve(Vcomps[0] , laminvK).t();                                     // lambda^{-1} K Sigma^{-1} which is equal to t(Sigma^{-1} lambda^{-1} K) 
//         hmeantemp = laminvKSiginv * (yz.col(0)-yz.cols(1,P_z)*gammakeep.row(s_outer).t());  // temporary vector to be used in mean and covariance
//         // hcovtemp = laminvKSiginv*(Vcomps[0]-laminvK)/sig2invkeep(s_outer);               // temporary covariance matrix to be used in marginal covariance
//         hcovtemp = (laminvK-laminvKSiginv*laminvK.t())/sig2invkeep(s_outer);                // temporary covariance matrix to be used in marginal covariance
//         // compute marginal mean and covariance for h
//         hmean += hmeantemp/(n_outer-n_burn);                                                                    // law of total expecation, E[E[h|s]]
//         hcov += hcovtemp/(n_outer-n_burn) + hmeantemp*hmeantemp.t()/(n_outer-n_burn); // first two components of law of total variance: E(Var(h|s))+E[E[h|s]^2]-E[E[h|s]]^2
//         // old (no B) : hcov += laminvKSiginv/(n_outer-n_burn)/sig2invkeep(s_outer) + hmeantemp*hmeantemp.t()/(n_outer-n_burn); // first two components of law of total variance: E(Var(h|s))+E[E[h|s]^2]-E[E[h|s]]^2
//         // hcov += laminvKSiginv*(Vcomps[0]-laminvK)/(n_outer-n_burn)/sig2invkeep(s_outer) + hmeantemp*hmeantemp.t()/(n_outer-n_burn); // first two components of law of total variance: E(Var(h|s))+E[E[h|s]^2]-E[E[h|s]]^2
//         if(draw_h==TRUE){
//           arma::svd(svdU,svdD,svdV,hcovtemp);                                                 // SVD decomposition for covariance of h (since it's more stable than cholesky)
//           // draw h
//           hsampkeep.row(s_outer) = ( hmeantemp+(svdU*diagmat(svdD))*arma::randn(N) ).t();     // hsampkeep.row(s_outer) = ( hmeantemp+arma::chol(hcovtemp)*arma::randn(N) ).t();
//         }
//       }
//       
//     } // end outer loop
//     
//     
//     
//     
//     // Convert field to list
//     for(int m=0; m<M; m++){
//       thetakeep[m] = thetakeepfield[m];
//     }
//     
//     // return a list
//     return List::create(Named("theta")          = thetakeep,
//                         Named("lambdaInverse")  = exp(logLambdaInversekeep),
//                         Named("lambdaBInverse") = exp(logLambdaBInversekeep),
//                         Named("rho")            = rhokeep,
//                         Named("sigma2")         = 1/sig2invkeep,
//                         Named("gamma")          = gammakeep,
//                         Named("hsamp")          = hsampkeep,
//                         Named("hmean")          = hmean,
//                         Named("hcov")           = hcov - ( hmean*hmean.t())); // subtracting third component of law of total variance (i.e. E[E[h|s]]^2)
// }

