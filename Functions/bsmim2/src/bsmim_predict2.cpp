#include "RcppArmadillo.h"
using namespace Rcpp;


// Construct \lambda^{-1}*K where K is the kernel matrix (do not add )
arma::field<arma::mat> get_Kmat(const bool&       poly,               // 0=gaussian kernel / 1=polynomial kernel
                                const int&        d,                  // degree of polynomial kernel
                                const int&        N,                  // no. of observations
                                const int&        points,             // no. of new grid points
                                const double&     logLambdaInverse,   // log(lam^{-1})
                                const arma::mat&  XthetaStar,
                                const arma::mat&  XthetaStar_new,               
                                const bool&       randint,            // indicator for random intercepts model
                                const double&     logLambdaBInverse,  // log(lamda^{-1}_B)
                                const arma::mat&  Bmat) {             // block diagonal matrix of cluster membership
  
  // Initialize matrices
  arma::mat lamInv_K(N,N,arma::fill::zeros);                // lambda^{-1} times kernel matrix for old data
  arma::mat lamInv_Knew(points,points,arma::fill::zeros);   // lambda^{-1} times kernel matrix for new data
  arma::mat lamInv_Knewold(points,N,arma::fill::zeros);     // lambda^{-1} times kernel matrix for new AND old data (off-diagonal)
  
  if(poly==TRUE){ // polynomial kernel
    
    // Construct lambda^{-1} times Kernel matrix for old values
    for(int i=0; i<N; i++){
      for(int j=i; j<N; j++){
        lamInv_K(i,j) = exp(logLambdaInverse)*pow(1 + as_scalar(XthetaStar.row(i)*XthetaStar.row(j).t()), d); 
        lamInv_K(j,i) = lamInv_K(i,j); // symmetry
      }  
      lamInv_K(i,i) = lamInv_K(i,i) ; // do not add I
    }
    
    // Construct lambda^{-1} times Kernel matrix for new values
    for(int p=0; p<points; p++){
      for(int q=p; q<points; q++){
        lamInv_Knew(p,q) = exp(logLambdaInverse)*pow(1 + as_scalar(XthetaStar_new.row(p)*XthetaStar_new.row(q).t()), d);
        lamInv_Knew(q,p) = lamInv_Knew(p,q);
      }
    }
    
    // Construct lambda^{-1} times Kernel matrix for new values and old values
    for(int p=0; p<points; p++){
      for(int j=0; j<N; j++){
        lamInv_Knewold(p,j) = exp(logLambdaInverse)*pow(1 + as_scalar(XthetaStar_new.row(p)*XthetaStar.row(j).t()), d);
      }
    }
    
    
  }else{ // gaussian kernel
    
    // Construct lambda^{-1} times Kernel matrix for old values
    for(int i=0; i<N-1; i++){
      for(int j=i+1; j<N; j++){
        lamInv_K(i,j) = exp(logLambdaInverse - arma::sum(square(vectorise(XthetaStar.row(i)-XthetaStar.row(j)))));
        lamInv_K(j,i) = lamInv_K(i,j); // symmetry
      }  
      lamInv_K(i,i) = exp(logLambdaInverse); // do not add I
    }
    lamInv_K(N-1,N-1) = exp(logLambdaInverse); // do not add I
    
    // Construct lambda^{-1} times Kernel matrix for new values
    for(int p=0; p<points-1; p++){
      for(int q=p+1; q<points; q++){
        lamInv_Knew(p,q) = exp(logLambdaInverse - arma::sum(square(vectorise(XthetaStar_new.row(p)-XthetaStar_new.row(q)))));
        lamInv_Knew(q,p) = lamInv_Knew(p,q);
      }
      lamInv_Knew(p,p) = exp(logLambdaInverse);
    }
    lamInv_Knew(points-1,points-1) = exp(logLambdaInverse);
    
    // Construct lambda^{-1} times Kernel matrix for new values and old values
    for(int p=0; p<points; p++){
      for(int j=0; j<N; j++){
        lamInv_Knewold(p,j) = exp(logLambdaInverse - arma::sum(square(vectorise(XthetaStar_new.row(p)-XthetaStar.row(j)))));
      } 
    }
    
  }

  // Construct Sigma=I+\lambda^{-1}*K+\lambda_B^{-1}B where K is the kernel matrix
  arma::mat Sigma(N,N,arma::fill::eye);           // initialize Sigma = I
  Sigma = Sigma + lamInv_K;                        // Sigma = I + \lambda^{-1}K
  if(randint==TRUE){                              // Add \lambda_B^{-1}B for random intercepts model
    Sigma = Sigma + exp(logLambdaBInverse)*Bmat;  // Sigma = I + \lambda^{-1}K + \lambda_B^{-1}B
  } 


  // output 
  arma::field<arma::mat> comps(3);
  comps[0] = Sigma;                         // Sigma=I+\lambda^{-1}*K+\lambda_B^{-1}B
  comps[1] = lamInv_Knew;                   // lambda^{-1} times kernel matrix for new data
  comps[2] = lamInv_Knewold;                // lambda^{-1} times kernel matrix for new AND old data
 
  return comps; 
}



// [[Rcpp::depends(RcppArmadillo)]]
//' Predicting hnew for BSMIM by Component
//' 
//' Returns a list containing (points by (sum_m Lm) ) matrices of hmean and hvar
//' 
//' @param yz matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
//' @param Xlist list of (N by L_m) matrices representing X_m
//' @param thetalist list of (N by L_m) matrices representing X_m
//' @param psilist list of basis matrices
//' @param rho (S by M) matrix of rho_m draws
//' @param gamma list of (S by P_z) matrices with gamma draws
//' @param lambdaInverse S-vector of lambda^{-1} draws
//' @param lambdaBInverse S-vector of lambda^{-1}_B draws
//' @param sigma2 S-vector of sigma^2 draws
//' @param weightslist list of (S by L_m) weight matrices to apply to exposures X_m (basically A theta^*, so that X*theta*=x(Atheta*)) 
//' @param gridpoints (points by M) matrix containing grid of new index levels for each index m
//' @param Xqlist list of L_m-vectors containing exposure quantiles
//' @param poly 0 = gaussian kernel, 1 = polynomial kernel
//' @param d degree of polynomial kernel
//' @param randint 0 = no random intercepts, 1 = random intercepts model
//' @param Bmat N xN block diagonal matrix indicating cluster membership for random intercepts model
//' @author Glen McGee and Ander Wilson (adapted from the "regimes" package by Ander Wilson)
//' @export
// [[Rcpp::export]]
List bsmim_predict_old_cpp2(const arma::mat&  yz,               // matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
                        const Rcpp::List& Xlist,            // list of (N by L_m) matrices representing X_m
                        const Rcpp::List& thetalist,        // list of (S by L_m) matrices with theta draws
                        const Rcpp::List& psilist,          // list of basis matrices
                        const arma::mat&  rho,              // (S by M) matrix of rho_m draws
                        const arma::mat&  gamma,            // list of (S by P_z) matrices with gamma draws
                        const arma::vec&  lambdaInverse,    // S-vector of lambda^{-1} draws
                        const arma::vec&  lambdaBInverse, // S-vector of lambda^{-1}_B draws
                        const arma::vec&  sigma2,           // S-vector of sigma^2 draws
                        const Rcpp::List& weightslist,      // list of (S by L_m) weight matrices to apply to exposures X_m (basically A theta^*, so that X*theta*=x(Atheta*)) 
                        const Rcpp::List& gridpointslist,   // (points by M) matrix containing grid of new index levels for each index m
                        const Rcpp::List& Xqlist,           // list of L_m-vectors containing componentwise exposure quantiles    
                        const bool&       poly,             // 0=gaussian kernel / 1=polynomial kernel
                        const int&        d,              // degree of polynomial kernel
                        const bool&       randint,        // 0=no random intercepts / 1=random intercepts
                        const arma::mat&  Bmat) {         // block diagonal matrix B indicating cluster membership for random intercepts model 
  // Dimensions
  int N = yz.n_rows;              // no. of observations
  int M = Xlist.size();           // no. of indices (dimension of h(.))
  int P_z = yz.n_cols-1;          // no. of covariates
  int S_iter = gamma.n_rows;      // no. of posterior draws 
  IntegerVector Lvec = rep(0,M);  // vector of L_m representing the no. of columns of X_m
  int col_lm;                     // iterator to index the column of output
  
  // Change lists into fields
  arma::field<arma::mat> X(M);                          // field of (N by L_m) matrices representing X_m
  arma::field<arma::mat> theta(M);                      // field of (S by L_m) vectors representing \thetaStar_m draws
  arma::field<arma::mat> psi(M);                        // field of basis matrices
  arma::field<arma::mat> weights(M);                    // field of weights matrices
  arma::field<arma::mat> gridpoints(M);                 // field of (points by Lm) matrices representing grid of new values
  arma::field<arma::vec> Xq(M);                         // field of exposure median vectors
  for(int m=0; m<M; m++){
    X[m] = as<arma::mat>(Xlist[m]);                     // mth matrix in X is X_m
    Lvec(m) = X[m].n_cols;                              // mth element of Lvec is L_m (the length of \thetaStar_m)
    theta[m] = as<arma::mat>(thetalist[m]);             // mth matrix in theta is (S by L_m) matrix of theta_m draws
    psi[m] = as<arma::mat>(psilist[m]);                 // mth matrix in psi is basis matrix for index m
    weights[m] = as<arma::mat>(weightslist[m]);         // mth matrix in weights is (S by L_m) matrix of A theta^*, so that X*theta*=x(Atheta*)
    gridpoints[m] = as<arma::mat>(gridpointslist[m]);   // mth matrix in gridpoints is (points by L_m) matrix of new values to predict at
    Xq[m] = as<arma::vec>(Xqlist[m]);                   // mth vector in Xmedian is L_m vector of exposure medians
  }
  int points = gridpoints[0].n_rows;                    // no. of grid points (new values)
  
  // Initialize other structures
  arma::vec y_zgam(N,arma::fill::zeros);                    // N-vector of residuals (Y-Z\gamma)
  arma::mat XthetaStar(N,M,arma::fill::zeros);              // (N by M) matrix representing [X_1\thetaStar_1,...,X_M\thetaStar_M]
  arma::mat XthetaStar_new(points,M,arma::fill::zeros);     // (points by M) matrix of Xnew theta* based on gridpoints (NEW values)
  arma::mat Sigma_inv(N,N,arma::fill::zeros);               // Construct matrix Sigma=I+\lambda^{-1}*K where K is the kernel matrix
  arma::mat lamInv_Knew(points,points,arma::fill::zeros);   // lambda^{-1} times kernel matrix for new data
  arma::mat lamInv_Knewold(points,N,arma::fill::zeros);     // lambda^{-1} times kernel matrix for new AND old data (off-diagonal)
  arma::vec hnew(points,arma::fill::zeros);                 // temporary vector for computing marginal covariance of h(.)
  arma::mat hmean(points,sum(Lvec),arma::fill::zeros);      // marginal mean vector for h(.)
  arma::mat hvar(points,sum(Lvec),arma::fill::zeros);       // marginal variance-covariance matrix of h(.)
  arma::vec ons(N, arma::fill::ones);                       // vector of length N of ones
  arma::field<arma::mat> Kmats;                             // field to contain output of get_Kmat()

  for(int s=0; s<S_iter; s++){                            // loop over posterior draws
  
    // residual vector (Y-Z\gamma)
    y_zgam = yz.col(0) - (yz.cols(1,P_z) * gamma.row(s).t()); // N-vector for Y-Z\gamma at the sth iteration
    
    // X\theta*
    for(int m=0; m<M; m++){    
      XthetaStar.col(m) =   X[m] * weights[m].row(s).t(); // set columns XthetaStar to Xm times theta_m* at the sth iteration
    }
    
    
    // Construct Xthetastar_new with gridpoints 
    for(int m=0; m<M; m++){             // loop over indices 
      for(int l=0; l<Lvec(m); l++){     // loop over components    
        
        if(m==0 & l==0){
          col_lm = 0;                   // reset iterator for column index
        }
        
        // Construct Xthetastar_new with gridpoints 
        for(int mm=0; mm<M; mm++){                      // loop over indices
          // for(int ll=0; ll<Lvec(mm); ll++){             // loop over components
          
          // if(mm==m & ll==l){                          // set to grid of new values for exposure of interest  
          if(mm==m){                          // set to grid of new values for exposure of interest 
            for(int p=0; p<points; p++){
              // XthetaStar_new(p,mm) = sqrt(rho(s,mm))*( arma::sum(vectorise(Xq[mm])%vectorise(weights[mm].row(s))) + weights[mm](s,ll)*(gridpoints[mm](p,ll) - Xq[mm][ll]) ); //rho^{1/2}*(Xm** thetam_s) where Xm** is equal to Xm at all median values except new values for exposure of interest
              // XthetaStar_new(p,mm) = sqrt(rho(s,mm))*( arma::sum(vectorise(Xq[mm])%vectorise(weights[mm].row(s))) + weights[mm](s,l)*(gridpoints[mm](p,l) - Xq[mm][l]) ); //rho^{1/2}*(Xm** thetam_s) where Xm** is equal to Xm at all median values except new values for exposure of interest
              // sqrt(rho) already included in weights
              XthetaStar_new(p,mm) =  arma::sum(vectorise(Xq[mm])%vectorise(weights[mm].row(s))) + weights[mm](s,l)*(gridpoints[mm](p,l) - Xq[mm][l]) ; //rho^{1/2}*(Xm** thetam_s) where Xm** is equal to Xm at all median values except new values for exposure of interest
            }
          }else{                                      // otherwise set to quantiles for comparison
            for(int p=0; p<points; p++){
              XthetaStar_new(p,mm) = arma::sum(vectorise(Xq[mm])%vectorise(weights[mm].row(s)));//Xmedian[mm].t() * weights[mm].row(s).t();
            }
          }
          
          // } // ll
        } // mm
        
        
        
        // Construct kernel matrices
        Kmats = get_Kmat(poly, d, N, points, log(lambdaInverse[s]), XthetaStar, XthetaStar_new, randint, log(lambdaBInverse[s]), Bmat);
        Sigma_inv = inv(Kmats[0]);      // Construct Sigma^{-1}
        lamInv_Knew = Kmats[1];         // Construct lambda^{-1} times Kernel matrix for new values
        lamInv_Knewold = Kmats[2];      // Construct lambda^{-1} times Kernel matrix for new and old values
        
        
        // Predict hnew
        hnew = lamInv_Knewold * Sigma_inv * y_zgam;                                                                 // mean of hnew at each draw is   lambda^{-1} Koldnew Sigma^{-1} (Y-Zgamma) 
        hmean.col(col_lm) += hnew/S_iter;                                                                           // law of total expectation over all draws
        hvar.col(col_lm) +=                                                                                         // law of total variance except for hmean^2 which will be subtracted post hoc (see below)
          sigma2[s] * (diagvec(lamInv_Knew) - diagvec(lamInv_Knewold * Sigma_inv * lamInv_Knewold.t()))/S_iter  +   // mean(var(hnew)), where var(hnew) is sigma^2(lambda^{-1}Knew-\lambda^{-2}KoldnewSigma^{-1}Koldnew^T)
          square(hnew)/S_iter;                                                                                      // mean of square (and subtract square of mean below)
        
        
        
        // index the columns of results
        col_lm+=1;     
        
      } // end l loop over components
      
    } // end m loop over indices
    
  } // end s loop over posterior draws
  
  
  
  
  
  // return a list
  return List::create(Named("hmean") = hmean,                    // mean of hnew 
                      Named("hsd") =  sqrt(hvar-square(hmean))); // subtracting third component of law of total variance (i.e. E[E[h|s]]^2)
  
  
}


// [[Rcpp::depends(RcppArmadillo)]]
//' Predicting hnew for BSMIM by Component
//' 
//' Returns a list containing (points by (sum_m Lm) ) matrices of hmean and hvar
//' 
//' @param yz matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
//' @param Xlist list of (N by L_m) matrices representing X_m
//' @param thetalist list of (N by L_m) matrices representing X_m
//' @param psilist list of basis matrices
//' @param rho (S by M) matrix of rho_m draws
//' @param gamma list of (S by P_z) matrices with gamma draws
//' @param lambdaInverse S-vector of lambda^{-1} draws
//' @param lambdaBInverse S-vector of lambda^{-1}_B draws
//' @param sigma2 S-vector of sigma^2 draws
//' @param weightslist list of (S by L_m) weight matrices to apply to exposures X_m (basically A theta^*, so that X*theta*=x(Atheta*)) 
//' @param gridpoints (points by M) matrix containing grid of new index levels for each index m
//' @param poly 0 = gaussian kernel, 1 = polynomial kernel
//' @param d degree of polynomial kernel
//' @param randint 0 = no random intercepts, 1 = random intercepts model
//' @param Bmat N xN block diagonal matrix indicating cluster membership for random intercepts model
//' @author Glen McGee and Ander Wilson (adapted from the "regimes" package by Ander Wilson)
//' @export
// [[Rcpp::export]]
List bsmim_predict_cpp2(const arma::mat&  yz,               // matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
                        const Rcpp::List& Xlist,            // list of (N by L_m) matrices representing X_m
                        const Rcpp::List& thetalist,        // list of (S by L_m) matrices with theta draws
                        const Rcpp::List& psilist,          // list of basis matrices
                        const arma::mat&  rho,              // (S by M) matrix of rho_m draws
                        const arma::mat&  gamma,            // list of (S by P_z) matrices with gamma draws
                        const arma::vec&  lambdaInverse,    // S-vector of lambda^{-1} draws
                        const arma::vec&  lambdaBInverse, // S-vector of lambda^{-1}_B draws
                        const arma::vec&  sigma2,           // S-vector of sigma^2 draws
                        const Rcpp::List& weightslist,      // list of (S by L_m) weight matrices to apply to exposures X_m (basically A theta^*, so that X*theta*=x(Atheta*)) 
                        const Rcpp::List& gridpointslist,   // (points by M) matrix containing grid of new index levels for each index m
                        const bool&       poly,             // 0=gaussian kernel / 1=polynomial kernel
                        const int&        d,              // degree of polynomial kernel
                        const bool&       randint,        // 0=no random intercepts / 1=random intercepts
                        const arma::mat&  Bmat) {         // block diagonal matrix B indicating cluster membership for random intercepts model 
  // Dimensions
  int N = yz.n_rows;              // no. of observations
  int M = Xlist.size();           // no. of indices (dimension of h(.))
  int P_z = yz.n_cols-1;          // no. of covariates
  int S_iter = gamma.n_rows;      // no. of posterior draws 
  IntegerVector Lvec = rep(0,M);  // vector of L_m representing the no. of columns of X_m
  // int col_lm;                     // iterator to index the column of output
  
  // Change lists into fields
  arma::field<arma::mat> X(M);                          // field of (N by L_m) matrices representing X_m
  arma::field<arma::mat> theta(M);                      // field of (S by L_m) vectors representing \thetaStar_m draws
  arma::field<arma::mat> psi(M);                        // field of basis matrices
  arma::field<arma::mat> weights(M);                    // field of weights matrices
  arma::field<arma::mat> gridpoints(M);                 // field of (points by Lm) matrices representing grid of new values
  arma::field<arma::vec> Xq(M);                         // field of exposure median vectors
  for(int m=0; m<M; m++){
    X[m] = as<arma::mat>(Xlist[m]);                     // mth matrix in X is X_m
    Lvec(m) = X[m].n_cols;                              // mth element of Lvec is L_m (the length of \thetaStar_m)
    theta[m] = as<arma::mat>(thetalist[m]);             // mth matrix in theta is (S by L_m) matrix of theta_m draws
    psi[m] = as<arma::mat>(psilist[m]);                 // mth matrix in psi is basis matrix for index m
    weights[m] = as<arma::mat>(weightslist[m]);         // mth matrix in weights is (S by L_m) matrix of A theta^*, so that X*theta*=x(Atheta*)
    gridpoints[m] = as<arma::mat>(gridpointslist[m]);   // mth matrix in gridpoints is (points by L_m) matrix of new values to predict at
  }
  int points = gridpoints[0].n_rows;                    // no. of grid points (new values)
  
  // Initialize other structures
  arma::vec y_zgam(N,arma::fill::zeros);                    // N-vector of residuals (Y-Z\gamma)
  arma::mat XthetaStar(N,M,arma::fill::zeros);              // (N by M) matrix representing [X_1\thetaStar_1,...,X_M\thetaStar_M]
  arma::mat XthetaStar_new(points,M,arma::fill::zeros);     // (points by M) matrix of Xnew theta* based on gridpoints (NEW values)
  arma::mat Sigma_inv(N,N,arma::fill::zeros);               // Construct matrix Sigma=I+\lambda^{-1}*K where K is the kernel matrix
  arma::mat lamInv_Knew(points,points,arma::fill::zeros);   // lambda^{-1} times kernel matrix for new data
  arma::mat lamInv_Knewold(points,N,arma::fill::zeros);     // lambda^{-1} times kernel matrix for new AND old data (off-diagonal)
  // arma::vec hnew(points,arma::fill::zeros);                 // temporary vector for computing marginal covariance of h(.)
  // arma::mat hmean(points,sum(Lvec),arma::fill::zeros);      // marginal mean vector for h(.)
  // arma::mat hvar(points,sum(Lvec),arma::fill::zeros);       // marginal variance-covariance matrix of h(.)
  arma::vec hmeantemp(points,arma::fill::zeros);            // temporary vector for computing marginal covariance of h(.)
  arma::vec hmean(points,arma::fill::zeros);                // marginal mean vector for h(.)
  arma::mat hcovtemp(points,points,arma::fill::zeros);      // marginal variance-covariance matrix of h(.)
  arma::mat hcov(points,points,arma::fill::zeros);          // marginal variance-covariance matrix of h(.)
  arma::vec ons(N, arma::fill::ones);                       // vector of length N of ones
  arma::field<arma::mat> Kmats;                             // field to contain output of get_Kmat()
  

  for(int s=0; s<S_iter; s++){                            // loop over posterior draws
    
    // residual vector (Y-Z\gamma)
    y_zgam = yz.col(0) - (yz.cols(1,P_z) * gamma.row(s).t()); // N-vector for Y-Z\gamma at the sth iteration
    
    // X\theta*
    for(int m=0; m<M; m++){    
      XthetaStar.col(m) =   X[m] * weights[m].row(s).t(); // set columns XthetaStar to Xm times theta_m* at the sth iteration
    }
    
    // Construct Xthetastar_new with gridpoints 
    for(int m=0; m<M; m++){    
      XthetaStar_new.col(m) =   gridpoints[m] * weights[m].row(s).t(); // set columns XthetaStar_new to Xm times theta_m* at the sth iteration
    }
    
    // for(int m=0; m<M; m++){             // loop over indices 
    //   for(int l=0; l<Lvec(m); l++){     // loop over components    
    
    
    
    // Construct kernel matrices
    Kmats = get_Kmat(poly, d, N, points, log(lambdaInverse[s]), XthetaStar, XthetaStar_new, randint, log(lambdaBInverse[s]), Bmat);
    Sigma_inv = inv(Kmats[0]);      // Construct Sigma^{-1}
    lamInv_Knew = Kmats[1];         // Construct lambda^{-1} times Kernel matrix for new values
    lamInv_Knewold = Kmats[2];      // Construct lambda^{-1} times Kernel matrix for new and old values
    
    
    // Predict hnew
    hmeantemp = lamInv_Knewold * Sigma_inv * y_zgam;                                          // temporary vector to be used in mean and covariance
    hcovtemp = sigma2[s] * (lamInv_Knew - lamInv_Knewold * Sigma_inv * lamInv_Knewold.t());   // temporary covariance matrix to be used in marginal covariance
    // compute marginal mean and covariance for h                                                                            
    hmean += hmeantemp/(S_iter);                                                              // law of total expecation, E[E[h|s]]
    hcov += hcovtemp/(S_iter) + hmeantemp*hmeantemp.t()/(S_iter);                             // first two components of law of total variance: E(Var(h|s))+E[E[h|s]^2]-E[E[h|s]]^2

    
    
  } // end s loop over posterior draws
  

  
  // return a list
  return List::create(Named("hmean") = hmean,                             // mean of hnew 
                      Named("hcov") = hcov - ( hmean*hmean.t()));          // subtracting third component of law of total variance (i.e. E[E[h|s]]^2)

  
  
}


// col_Means
arma::rowvec col_Means(const arma::mat&  X) {             
  
  arma::rowvec meanvec(X.n_cols,arma::fill::zeros);             
  for(int j=0; j<X.n_cols; j++){
    for(int i=0; i<X.n_rows; i++){
      meanvec(j) += X(i,j)/(X.n_rows);
    }  
  }
 
  return meanvec; 
}

// col_Medians
arma::rowvec col_Medians(const arma::mat&  X) {             
  
  arma::rowvec meanvec(X.n_cols,arma::fill::zeros);             
  for(int j=0; j<X.n_cols; j++){
      meanvec(j) = median(X.col(j));
  }
  
  return meanvec; 
}



// [[Rcpp::depends(RcppArmadillo)]]
//' Predicting hnew for BSMIM by Component using APPROXIMATE method, computing the conditional mean and variance of hnew given posterior means of everything else
//' 
//' Returns a list containing (points by (sum_m Lm) ) matrices of hmean and hvar
//' 
//' @param yz matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
//' @param Xlist list of (N by L_m) matrices representing X_m
//' @param thetalist list of (N by L_m) matrices representing X_m
//' @param psilist list of basis matrices
//' @param rho (S by M) matrix of rho_m draws
//' @param gamma list of (S by P_z) matrices with gamma draws
//' @param lambdaInverse S-vector of lambda^{-1} draws
//' @param lambdaBInverse S-vector of lambda^{-1}_B draws
//' @param sigma2 S-vector of sigma^2 draws
//' @param weightslist list of (S by L_m) weight matrices to apply to exposures X_m (basically A theta^*, so that X*theta*=x(Atheta*)) 
//' @param gridpoints (points by M) matrix containing grid of new index levels for each index m
//' @param Xqlist list of L_m-vectors containing exposure quantiles
//' @param poly 0 = gaussian kernel, 1 = polynomial kernel
//' @param d degree of polynomial kernel
//' @param randint 0 = no random intercepts, 1 = random intercepts model
//' @param Bmat N xN block diagonal matrix indicating cluster membership for random intercepts model
//' @author Glen McGee (adapted from the "regimes" package by Ander Wilson)
//' @export
// [[Rcpp::export]]
List bsmim_predict_approx_old_cpp2(const arma::mat&  yz,               // matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
                               const Rcpp::List& Xlist,            // list of (N by L_m) matrices representing X_m
                               const Rcpp::List& thetalist,        // list of (S by L_m) matrices with theta draws
                               const Rcpp::List& psilist,          // list of basis matrices
                               const arma::mat&  rho,              // (S by M) matrix of rho_m draws
                               const arma::mat&  gamma,            // list of (S by P_z) matrices with gamma draws
                               const arma::vec&  lambdaInverse,    // S-vector of lambda^{-1} draws
                               const arma::vec&  lambdaBInverse, // S-vector of lambda^{-1}_B draws
                               const arma::vec&  sigma2,           // S-vector of sigma^2 draws
                               const Rcpp::List& weightslist,      // list of (S by L_m) weight matrices to apply to exposures X_m (basically A theta^*, so that X*theta*=x(Atheta*)) 
                               const Rcpp::List& gridpointslist,   // (points by M) matrix containing grid of new index levels for each index m
                               const Rcpp::List& Xqlist,           // list of L_m-vectors containing componentwise exposure quantiles    
                               const bool&       poly,             // 0=gaussian kernel / 1=polynomial kernel
                               const int&        d,              // degree of polynomial kernel
                               const bool&       randint,        // 0=no random intercepts / 1=random intercepts
                               const arma::mat&  Bmat) {         // block diagonal matrix B indicating cluster membership for random intercepts model 
  // Dimensions
  int N = yz.n_rows;              // no. of observations
  int M = Xlist.size();           // no. of indices (dimension of h(.))
  int P_z = yz.n_cols-1;          // no. of covariates
  // int S_iter = gamma.n_rows;      // no. of posterior draws 
  IntegerVector Lvec = rep(0,M);  // vector of L_m representing the no. of columns of X_m
  int col_lm;                     // iterator to index the column of output
  
  // Change lists into fields
  arma::field<arma::mat> X(M);                          // field of (N by L_m) matrices representing X_m
  arma::field<arma::mat> theta(M);                      // field of (S by L_m) vectors representing \thetaStar_m draws
  arma::field<arma::mat> psi(M);                        // field of basis matrices
  arma::field<arma::mat> weights(M);                    // field of weights matrices
  arma::field<arma::mat> gridpoints(M);                 // field of (points by Lm) matrices representing grid of new values
  arma::field<arma::vec> Xq(M);                         // field of exposure median vectors
  for(int m=0; m<M; m++){
    X[m] = as<arma::mat>(Xlist[m]);                     // mth matrix in X is X_m
    Lvec(m) = X[m].n_cols;                              // mth element of Lvec is L_m (the length of \thetaStar_m)
    theta[m] = as<arma::mat>(thetalist[m]);             // mth matrix in theta is (S by L_m) matrix of theta_m draws
    psi[m] = as<arma::mat>(psilist[m]);                 // mth matrix in psi is basis matrix for index m
    weights[m] = as<arma::mat>(weightslist[m]);         // mth matrix in weights is (S by L_m) matrix of A theta^*, so that X*theta*=x(Atheta*)
    gridpoints[m] = as<arma::mat>(gridpointslist[m]);   // mth matrix in gridpoints is (points by L_m) matrix of new values to predict at
    Xq[m] = as<arma::vec>(Xqlist[m]);                   // mth vector in Xmedian is L_m vector of exposure medians
  }
  int points = gridpoints[0].n_rows;                    // no. of grid points (new values)
  
  // Initialize other structures
  arma::vec y_zgam(N,arma::fill::zeros);                    // N-vector of residuals (Y-Z\gamma)
  arma::mat XthetaStar(N,M,arma::fill::zeros);              // (N by M) matrix representing [X_1\thetaStar_1,...,X_M\thetaStar_M]
  arma::mat XthetaStar_new(points,M,arma::fill::zeros);     // (points by M) matrix of Xnew theta* based on gridpoints (NEW values)
  arma::mat Sigma_inv(N,N,arma::fill::zeros);               // Construct matrix Sigma=I+\lambda^{-1}*K where K is the kernel matrix
  arma::mat lamInv_Knew(points,points,arma::fill::zeros);   // lambda^{-1} times kernel matrix for new data
  arma::mat lamInv_Knewold(points,N,arma::fill::zeros);     // lambda^{-1} times kernel matrix for new AND old data (off-diagonal)
  arma::vec hnew(points,arma::fill::zeros);                 // temporary vector for computing marginal covariance of h(.)
  arma::mat hmean(points,sum(Lvec),arma::fill::zeros);      // marginal mean vector for h(.)
  arma::mat hvar(points,sum(Lvec),arma::fill::zeros);       // marginal variance-covariance matrix of h(.)
  arma::vec ons(N, arma::fill::ones);                       // vector of length N of ones
  arma::field<arma::mat> Kmats;                             // field to contain output of get_Kmat()
  

    
    // residual vector (Y-Z\gamma)
    y_zgam = yz.col(0) - (yz.cols(1,P_z) * col_Means(gamma).t()); // N-vector for Y-Z\gamma at the sth iteration
    
    // X\theta*
    for(int m=0; m<M; m++){    
      XthetaStar.col(m) =   X[m] * col_Means(weights[m]).t(); // set columns XthetaStar to Xm times theta_m* at the sth iteration
    }
    
    
    // Construct Xthetastar_new with gridpoints 
    for(int m=0; m<M; m++){             // loop over indices 
      for(int l=0; l<Lvec(m); l++){     // loop over components    
        
        if(m==0 & l==0){
          col_lm = 0;                   // reset iterator for column index
        }
        
        // Construct Xthetastar_new with gridpoints 
        for(int mm=0; mm<M; mm++){                      // loop over indices
          // for(int ll=0; ll<Lvec(mm); ll++){             // loop over components
          
          // if(mm==m & ll==l){                          // set to grid of new values for exposure of interest  
          if(mm==m){                          // set to grid of new values for exposure of interest 
            for(int p=0; p<points; p++){
              // XthetaStar_new(p,mm) = sqrt(rho(s,mm))*( arma::sum(vectorise(Xq[mm])%vectorise(weights[mm].row(s))) + weights[mm](s,ll)*(gridpoints[mm](p,ll) - Xq[mm][ll]) ); //rho^{1/2}*(Xm** thetam_s) where Xm** is equal to Xm at all median values except new values for exposure of interest
              // XthetaStar_new(p,mm) = sqrt(rho(s,mm))*( arma::sum(vectorise(Xq[mm])%vectorise(weights[mm].row(s))) + weights[mm](s,l)*(gridpoints[mm](p,l) - Xq[mm][l]) ); //rho^{1/2}*(Xm** thetam_s) where Xm** is equal to Xm at all median values except new values for exposure of interest
              // sqrt(rho) already included in weights
              XthetaStar_new(p,mm) =  arma::sum(vectorise(Xq[mm])%vectorise(col_Means(weights[mm]))) + col_Means(weights[mm])[l]*(gridpoints[mm](p,l) - Xq[mm][l]) ; //rho^{1/2}*(Xm** thetam_s) where Xm** is equal to Xm at all median values except new values for exposure of interest
            }
          }else{                                      // otherwise set to quantiles for comparison
            for(int p=0; p<points; p++){
              XthetaStar_new(p,mm) = arma::sum(vectorise(Xq[mm])%vectorise(col_Means(weights[mm])));//Xmedian[mm].t() * weights[mm].row(s).t();
            }
          }
          
          // } // ll
        } // mm
        
        
        
        // Construct kernel matrices
        Kmats = get_Kmat(poly, d, N, points, log(mean(lambdaInverse)), XthetaStar, XthetaStar_new, randint, log(mean(lambdaBInverse)), Bmat);
        Sigma_inv = inv(Kmats[0]);      // Construct Sigma^{-1}
        lamInv_Knew = Kmats[1];         // Construct lambda^{-1} times Kernel matrix for new values
        lamInv_Knewold = Kmats[2];      // Construct lambda^{-1} times Kernel matrix for new and old values
        
        
        // Predict hnew
        hmean.col(col_lm) = lamInv_Knewold * Sigma_inv * y_zgam;           // mean of hnew at each draw is   lambda^{-1} Koldnew Sigma^{-1} (Y-Zgamma) 
        hvar.col(col_lm) =  mean(sigma2) * (diagvec(lamInv_Knew) - diagvec(lamInv_Knewold * Sigma_inv * lamInv_Knewold.t()));   // var(hnew) is sigma^2(lambda^{-1}Knew-\lambda^{-2}KoldnewSigma^{-1}Koldnew^T)

        
        // index the columns of results
        col_lm+=1;     
        
      } // end l loop over components
      
    } // end m loop over indices
    

  

  
  // return a list
  return List::create(Named("hmean") = hmean,      // conditional mean of hnew 
                      Named("hsd") =  sqrt(hvar)); // conditional sd of hnew (no iterated expectation needed)
  
  
}

// [[Rcpp::depends(RcppArmadillo)]]
//' Predicting hnew for BSMIM by Component using APPROXIMATE method, computing the conditional mean and variance of hnew given posterior means of everything else
//' 
//' Returns a list containing (points by (sum_m Lm) ) matrices of hmean and hvar
//' 
//' @param yz matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
//' @param Xlist list of (N by L_m) matrices representing X_m
//' @param thetalist list of (N by L_m) matrices representing X_m
//' @param psilist list of basis matrices
//' @param rho (S by M) matrix of rho_m draws
//' @param gamma list of (S by P_z) matrices with gamma draws
//' @param lambdaInverse S-vector of lambda^{-1} draws
//' @param lambdaBInverse S-vector of lambda^{-1}_B draws
//' @param sigma2 S-vector of sigma^2 draws
//' @param weightslist list of (S by L_m) weight matrices to apply to exposures X_m (basically A theta^*, so that X*theta*=x(Atheta*)) 
//' @param gridpoints (points by M) matrix containing grid of new index levels for each index m
//' @param Xqlist list of L_m-vectors containing exposure quantiles
//' @param poly 0 = gaussian kernel, 1 = polynomial kernel
//' @param d degree of polynomial kernel
//' @param randint 0 = no random intercepts, 1 = random intercepts model
//' @param Bmat N xN block diagonal matrix indicating cluster membership for random intercepts model
//' @author Glen McGee (adapted from the "regimes" package by Ander Wilson)
//' @export
// [[Rcpp::export]]
List bsmim_predict_approx_cpp2(const arma::mat&  yz,               // matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
                               const Rcpp::List& Xlist,            // list of (N by L_m) matrices representing X_m
                               const Rcpp::List& thetalist,        // list of (S by L_m) matrices with theta draws
                               const Rcpp::List& psilist,          // list of basis matrices
                               const arma::mat&  rho,              // (S by M) matrix of rho_m draws
                               const arma::mat&  gamma,            // list of (S by P_z) matrices with gamma draws
                               const arma::vec&  lambdaInverse,    // S-vector of lambda^{-1} draws
                               const arma::vec&  lambdaBInverse, // S-vector of lambda^{-1}_B draws
                               const arma::vec&  sigma2,           // S-vector of sigma^2 draws
                               const Rcpp::List& weightslist,      // list of (S by L_m) weight matrices to apply to exposures X_m (basically A theta^*, so that X*theta*=x(Atheta*)) 
                               const Rcpp::List& gridpointslist,   // (points by M) matrix containing grid of new index levels for each index m
                               const bool&       poly,             // 0=gaussian kernel / 1=polynomial kernel
                               const int&        d,              // degree of polynomial kernel
                               const bool&       randint,        // 0=no random intercepts / 1=random intercepts
                               const arma::mat&  Bmat) {         // block diagonal matrix B indicating cluster membership for random intercepts model 

  // Dimensions
  int N = yz.n_rows;              // no. of observations
  int M = Xlist.size();           // no. of indices (dimension of h(.))
  int P_z = yz.n_cols-1;          // no. of covariates
  // int S_iter = gamma.n_rows;      // no. of posterior draws 
  IntegerVector Lvec = rep(0,M);  // vector of L_m representing the no. of columns of X_m
  // int col_lm;                     // iterator to index the column of output
  
  // Change lists into fields
  arma::field<arma::mat> X(M);                          // field of (N by L_m) matrices representing X_m
  arma::field<arma::mat> theta(M);                      // field of (S by L_m) vectors representing \thetaStar_m draws
  arma::field<arma::mat> psi(M);                        // field of basis matrices
  arma::field<arma::mat> weights(M);                    // field of weights matrices
  arma::field<arma::mat> gridpoints(M);                 // field of (points by Lm) matrices representing grid of new values
  arma::field<arma::vec> Xq(M);                         // field of exposure median vectors
  for(int m=0; m<M; m++){
    X[m] = as<arma::mat>(Xlist[m]);                     // mth matrix in X is X_m
    Lvec(m) = X[m].n_cols;                              // mth element of Lvec is L_m (the length of \thetaStar_m)
    theta[m] = as<arma::mat>(thetalist[m]);             // mth matrix in theta is (S by L_m) matrix of theta_m draws
    psi[m] = as<arma::mat>(psilist[m]);                 // mth matrix in psi is basis matrix for index m
    weights[m] = as<arma::mat>(weightslist[m]);         // mth matrix in weights is (S by L_m) matrix of A theta^*, so that X*theta*=x(Atheta*)
    gridpoints[m] = as<arma::mat>(gridpointslist[m]);   // mth matrix in gridpoints is (points by L_m) matrix of new values to predict at
  }
  int points = gridpoints[0].n_rows;                    // no. of grid points (new values)
  
  // Initialize other structures
  arma::vec y_zgam(N,arma::fill::zeros);                    // N-vector of residuals (Y-Z\gamma)
  arma::mat XthetaStar(N,M,arma::fill::zeros);              // (N by M) matrix representing [X_1\thetaStar_1,...,X_M\thetaStar_M]
  arma::mat XthetaStar_new(points,M,arma::fill::zeros);     // (points by M) matrix of Xnew theta* based on gridpoints (NEW values)
  arma::mat Sigma_inv(N,N,arma::fill::zeros);               // Construct matrix Sigma=I+\lambda^{-1}*K where K is the kernel matrix
  arma::mat lamInv_Knew(points,points,arma::fill::zeros);   // lambda^{-1} times kernel matrix for new data
  arma::mat lamInv_Knewold(points,N,arma::fill::zeros);     // lambda^{-1} times kernel matrix for new AND old data (off-diagonal)
  // arma::vec hnew(points,arma::fill::zeros);                 // temporary vector for computing marginal covariance of h(.)
  // arma::mat hmean(points,sum(Lvec),arma::fill::zeros);      // marginal mean vector for h(.)
  // arma::mat hvar(points,sum(Lvec),arma::fill::zeros);       // marginal variance-covariance matrix of h(.)
  // arma::vec hmeantemp(points,arma::fill::zeros);            // temporary vector for computing marginal covariance of h(.)
  arma::vec hmean(points,arma::fill::zeros);                // marginal mean vector for h(.)
  // arma::mat hcovtemp(points,points,arma::fill::zeros);      // marginal variance-covariance matrix of h(.)
  arma::mat hcov(points,points,arma::fill::zeros);          // marginal variance-covariance matrix of h(.)
  arma::vec ons(N, arma::fill::ones);                       // vector of length N of ones
  arma::field<arma::mat> Kmats;                             // field to contain output of get_Kmat()
  
  
  
  // residual vector (Y-Z\gamma)
  y_zgam = yz.col(0) - (yz.cols(1,P_z) * col_Medians(gamma).t()); // N-vector for Y-Z\gamma at the sth iteration
  
  // X\theta*
  for(int m=0; m<M; m++){    
    XthetaStar.col(m) =   X[m] * col_Medians(weights[m]).t(); // set columns XthetaStar to Xm times theta_m* at the sth iteration
  }
  
  // Construct Xthetastar_new with gridpoints 
  for(int m=0; m<M; m++){    
    XthetaStar_new.col(m) =   gridpoints[m] * col_Medians(weights[m]).t(); // set columns XthetaStar_new to Xm times theta_m* at the sth iteration
  }
  
  // for(int m=0; m<M; m++){             // loop over indices 
  //   for(int l=0; l<Lvec(m); l++){     // loop over components    
  
  
  
  // Construct kernel matrices
  Kmats = get_Kmat(poly, d, N, points, log(median(lambdaInverse)), XthetaStar, XthetaStar_new, randint, log(median(lambdaBInverse)), Bmat);
  Sigma_inv = inv(Kmats[0]);      // Construct Sigma^{-1}
  lamInv_Knew = Kmats[1];         // Construct lambda^{-1} times Kernel matrix for new values
  lamInv_Knewold = Kmats[2];      // Construct lambda^{-1} times Kernel matrix for new and old values
  
  
  // Predict hnew
  hmean = lamInv_Knewold * Sigma_inv * y_zgam;                                                    // mean of hnew at each draw is   lambda^{-1} Koldnew Sigma^{-1} (Y-Zgamma)           
  hcov = median(sigma2) * (lamInv_Knew - lamInv_Knewold * Sigma_inv * lamInv_Knewold.t());        // var(hnew) is sigma^2(lambda^{-1}Knew-\lambda^{-2}KoldnewSigma^{-1}Koldnew^T)           
  
  
  


  
  // return a list
  return List::create(Named("hmean") = hmean,         // mean of hnew 
                      Named("hcov") = hcov );         // var of hnew

  
}







// [[Rcpp::depends(RcppArmadillo)]]
//' Predicting hnew for BSMIM by Index
//' 
//' 
//' 
//' @param yz matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
//' @param Xlist list of (N by L_m) matrices representing X_m
//' @param thetalist list of (N by L_m) matrices representing X_m
//' @param psilist list of basis matrices
//' @param rho (S by M) matrix of rho_m draws
//' @param gamma list of (S by P_z) matrices with gamma draws
//' @param lambdaInverse S-vector of lambda^{-1} draws
//' @param lambdaBInverse S-vector of lambda^{-1}_B draws
//' @param sigma2 S-vector of sigma^2 draws
//' @param weightslist list of (S by L_m) weight matrices to apply to exposures X_m (basically A theta^*, so that X*theta*=x(Atheta*)) 
//' @param gridpoints (points by M) matrix containing grid of new index levels for each index m
//' @param Xmedianlist list of L_m-vectors containing exposure medians   
//' @param Eq M-vector of index quantiles
//' @param crossM index to set cross section
//' @param poly 0 = gaussian kernel, 1 = polynomial kernel
//' @param d degree of polynomial kernel
//' @param randint 0 = no random intercepts, 1 = random intercepts model
//' @param Bmat N xN block diagonal matrix indicating cluster membership for random intercepts model
//' @author Glen McGee and Ander Wilson (adapted from the "regimes" package by Ander Wilson)
//' @export
// [[Rcpp::export]]
List bsmim_predict_indexwise_cpp2(const arma::mat&  yz,             // matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
                                  const Rcpp::List& Xlist,          // list of (N by L_m) matrices representing X_m
                                  const Rcpp::List& thetalist,      // list of (S by L_m) matrices with theta draws
                                  const Rcpp::List& psilist,        // list of basis matrices
                                  const arma::mat&  rho,            // (S by M) matrix of rho_m draws
                                  const arma::mat&  gamma,          // list of (S by P_z) matrices with gamma draws
                                  const arma::vec&  lambdaInverse,  // S-vector of lambda^{-1} draws
                                  const arma::vec&  lambdaBInverse, // S-vector of lambda^{-1}_B draws
                                  const arma::vec&  sigma2,         // S-vector of sigma^2 draws
                                  const Rcpp::List& weightslist,    // list of (S by L_m) weight matrices to apply to exposures X_m (basically A theta^*, so that X*theta*=x(Atheta*)) 
                                  const arma::mat&  gridpoints,     // (points by M) matrix containing grid of new index levels for each index m
                                  const Rcpp::List& Xmedianlist,    // list of L_m-vectors containing exposure medians    
                                  const arma::vec&  Eq,             // M-vector of index quantiles
                                  const int&        crossM,         // index to set cross section
                                  const bool&       poly,           // 0=gaussian kernel / 1=polynomial kernel
                                  const int&        d,              // degree of polynomial kernel
                                  const bool&       randint,        // 0=no random intercepts / 1=random intercepts
                                  const arma::mat&  Bmat) {         // block diagonal matrix B indicating cluster membership for random intercepts model         
  
  // Dimensions
  int N = yz.n_rows;              // no. of observations
  int M = Xlist.size();           // no. of indices (dimension of h(.))
  int P_z = yz.n_cols-1;          // no. of covariates
  int S_iter = gamma.n_rows;      // no. of posterior draws 
  int points = gridpoints.n_rows; // no. of grid points (new values)
  
  // Change lists into fields
  arma::field<arma::mat> X(M);            // field of (N by L_m) matrices representing X_m
  arma::field<arma::mat> theta(M);        // field of (S by L_m) vectors representing \thetaStar_m draws
  arma::field<arma::mat> psi(M);          // field of basis matrices
  arma::field<arma::mat> weights(M);      // field of weights matrices
  arma::field<arma::vec> Xmedian(M);      // field of exposure median vectors
  for(int m=0; m<M; m++){
    X[m] = as<arma::mat>(Xlist[m]);                 // mth matrix in X is X_m
    theta[m] = as<arma::mat>(thetalist[m]);         // mth matrix in theta is (S by L_m) matrix of theta_m draws
    psi[m] = as<arma::mat>(psilist[m]);             // mth matrix in psi is basis matrix for index m
    weights[m] = as<arma::mat>(weightslist[m]);     // mth matrix in weights is (S by L_m) matrix of A theta^*, so that X*theta*=x(Atheta*)
    Xmedian[m] = as<arma::vec>(Xmedianlist[m]);     // mth vector in Xmedian is L_m vector of exposure medians
  }
  
  // Initialize other structures
  arma::vec y_zgam(N,arma::fill::zeros);                    // N-vector of residuals (Y-Z\gamma)
  arma::mat XthetaStar(N,M,arma::fill::zeros);              // (N by M) matrix representing [X_1\thetaStar_1,...,X_M\thetaStar_M]
  arma::mat XthetaStar_new(points,M,arma::fill::zeros);     // (points by M) matrix of Xnew theta* based on gridpoints (NEW values)
  arma::mat Sigma_inv(N,N,arma::fill::zeros);               // Construct matrix Sigma=I+\lambda^{-1}*K where K is the kernel matrix
  arma::mat lamInv_Knew(points,points,arma::fill::zeros);   // lambda^{-1} times kernel matrix for new data
  arma::mat lamInv_Knewold(points,N,arma::fill::zeros);     // lambda^{-1} times kernel matrix for new AND old data (off-diagonal)
  arma::vec hnew(points,arma::fill::zeros);                 // temporary vector for computing marginal covariance of h(.)
  arma::mat hmean(points,M,arma::fill::zeros);              // marginal mean vector for h(.)
  arma::mat hvar(points,M,arma::fill::zeros);               // marginal variance-covariance matrix of h(.)
  arma::vec ons(N, arma::fill::ones);                       // vector of length N of ones
  arma::field<arma::mat> Kmats;                             // field to contain output of get_Kmat()
    
  
  for(int s=0; s<S_iter; s++){                            // loop over posterior draws
    // if(s % 10 == 0) Rcout << "#" << s << std::endl;           // output progress every 100 iterations
    
    // residual vector (Y-Z\gamma)
    y_zgam = yz.col(0) - (yz.cols(1,P_z) * gamma.row(s).t()); // N-vector for Y-Z\gamma at the sth iteration
    
    // X\theta*
    for(int m=0; m<M; m++){    
      XthetaStar.col(m) =   X[m] * weights[m].row(s).t(); // set columns XthetaStar to Xm times theta_m* at the sth iteration
    }
    

    // Construct Xthetastar_new with gridpoints 
    for(int m=0; m<M; m++){     // loop over indices 
      // Construct Xthetastar_new with gridpoints 
      for(int mm=0; mm<M; mm++){                                        // gridpoints for exposure of interest
        if(mm == m) {
          XthetaStar_new.col(mm) = sqrt(rho(s,mm))*gridpoints.col(mm);
        } else if(mm == (crossM-1)){                                    // and with cross sectional values for interaction
          for(int p=0; p<points; p++){     
            XthetaStar_new(p,mm) = sqrt(rho(s,mm))*Eq[mm];
          }
        } else {                                                        // and with median values elsewhere 
          for(int p=0; p<points; p++){     
            XthetaStar_new(p,mm) = arma::sum(vectorise(Xmedian[mm])%vectorise(weights[mm].row(s)));//Xmedian[mm].t() * weights[mm].row(s).t();
          }
        }
      }
      
      
      // Construct kernel matrices
      Kmats = get_Kmat(poly, d, N, points, log(lambdaInverse[s]), XthetaStar, XthetaStar_new, randint, log(lambdaBInverse[s]), Bmat);
      Sigma_inv = inv(Kmats[0]);      // Construct Sigma^{-1}
      lamInv_Knew = Kmats[1];         // Construct lambda^{-1} times Kernel matrix for new values
      lamInv_Knewold = Kmats[2];      // Construct lambda^{-1} times Kernel matrix for new and old values
      
    
      // Predict hnew
      hnew = lamInv_Knewold * Sigma_inv * y_zgam;                                                                 // mean of hnew at each draw is   lambda^{-1} Koldnew Sigma^{-1} (Y-Zgamma) 
      hmean.col(m) += hnew/S_iter;                                                                                // law of total expectation over all draws
      hvar.col(m) +=                                                                                              // law of total variance except for hmean^2 which will be subtracted post hoc (see below)
        sigma2[s] * (diagvec(lamInv_Knew) - diagvec(lamInv_Knewold * Sigma_inv * lamInv_Knewold.t()))/S_iter  +   // mean(var(hnew)), where var(hnew) is sigma^2(lambda^{-1}Knew-\lambda^{-2}KoldnewSigma^{-1}Koldnew^T)
        square(hnew)/S_iter;                                                                                      // mean of square (and subtract square of mean below)
      
      
      
    } // end m loop over indices
    
  } // end s loop over posterior draws
  
  
  
  
  
  // return a list
  return List::create(Named("hmean") = hmean,                    // mean of hnew 
                      Named("hsd") =  sqrt(hvar-square(hmean))); // subtracting third component of law of total variance (i.e. E[E[h|s]]^2)
  
  
}


// [[Rcpp::depends(RcppArmadillo)]]
//' Generic function for predicting hnew for BSMIM given new exposure levels X
//' 
//' Returns a matrix containing hmean and hvar
//' 
//' @param yz matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
//' @param Xlist list of (N by L_m) matrices representing X_m
//' @param thetalist list of (N by L_m) matrices representing X_m
//' @param psilist list of basis matrices
//' @param rho (S by M) matrix of rho_m draws
//' @param gamma list of (S by P_z) matrices with gamma draws
//' @param lambdaInverse S-vector of lambda^{-1} draws
//' @param lambdaBInverse S-vector of lambda^{-1}_B draws
//' @param sigma2 S-vector of sigma^2 draws
//' @param weightslist list of (S by L_m) weight matrices to apply to exposures X_m (basically A theta^*, so that X*theta*=x(Atheta*)) 
//' @param gridpoints (points by M) matrix containing grid of new index levels for each index m
//' @param poly 0 = gaussian kernel, 1 = polynomial kernel
//' @param d degree of polynomial kernel
//' @param randint 0 = no random intercepts, 1 = random intercepts model
//' @param Bmat N xN block diagonal matrix indicating cluster membership for random intercepts model
//' @author Glen McGee and Ander Wilson (adapted from the "regimes" package by Ander Wilson)
//' @export
// [[Rcpp::export]]
List bsmim_predict_X_cpp2(const arma::mat&  yz,               // matrix [Y,Z], Z does not include a vector of 1s (i.e. no intercept)
                        const Rcpp::List& Xlist,            // list of (N by L_m) matrices representing X_m
                        const Rcpp::List& thetalist,        // list of (S by L_m) matrices with theta draws
                        const Rcpp::List& psilist,          // list of basis matrices
                        const arma::mat&  rho,              // (S by M) matrix of rho_m draws
                        const arma::mat&  gamma,            // list of (S by P_z) matrices with gamma draws
                        const arma::vec&  lambdaInverse,    // S-vector of lambda^{-1} draws
                        const arma::vec&  lambdaBInverse, // S-vector of lambda^{-1}_B draws
                        const arma::vec&  sigma2,           // S-vector of sigma^2 draws
                        const Rcpp::List& weightslist,      // list of (S by L_m) weight matrices to apply to exposures X_m (basically A theta^*, so that X*theta*=x(Atheta*)) 
                        const Rcpp::List& gridpointslist,   // (points by M) matrix containing grid of new index levels for each index m
                        const bool&       poly,             // 0=gaussian kernel / 1=polynomial kernel
                        const int&        d,              // degree of polynomial kernel
                        const bool&       randint,        // 0=no random intercepts / 1=random intercepts
                        const arma::mat&  Bmat) {         // block diagonal matrix B indicating cluster membership for random intercepts model 
  // Dimensions
  int N = yz.n_rows;              // no. of observations
  int M = Xlist.size();           // no. of indices (dimension of h(.))
  int P_z = yz.n_cols-1;          // no. of covariates
  int S_iter = gamma.n_rows;      // no. of posterior draws 
  IntegerVector Lvec = rep(0,M);  // vector of L_m representing the no. of columns of X_m
  
  // Change lists into fields
  arma::field<arma::mat> X(M);                          // field of (N by L_m) matrices representing X_m
  arma::field<arma::mat> theta(M);                      // field of (S by L_m) vectors representing \thetaStar_m draws
  arma::field<arma::mat> psi(M);                        // field of basis matrices
  arma::field<arma::mat> weights(M);                    // field of weights matrices
  arma::field<arma::mat> gridpoints(M);                 // field of (points by Lm) matrices representing grid of new values
  for(int m=0; m<M; m++){
    X[m] = as<arma::mat>(Xlist[m]);                     // mth matrix in X is X_m
    Lvec(m) = X[m].n_cols;                              // mth element of Lvec is L_m (the length of \thetaStar_m)
    theta[m] = as<arma::mat>(thetalist[m]);             // mth matrix in theta is (S by L_m) matrix of theta_m draws
    psi[m] = as<arma::mat>(psilist[m]);                 // mth matrix in psi is basis matrix for index m
    weights[m] = as<arma::mat>(weightslist[m]);         // mth matrix in weights is (S by L_m) matrix of A theta^*, so that X*theta*=x(Atheta*)
    gridpoints[m] = as<arma::mat>(gridpointslist[m]);   // mth matrix in gridpoints is (points by L_m) matrix of new values to predict at
  }
  int points = gridpoints[0].n_rows;                    // no. of grid points (new values)
  
  // Initialize other structures
  arma::vec y_zgam(N,arma::fill::zeros);                    // N-vector of residuals (Y-Z\gamma)
  arma::mat XthetaStar(N,M,arma::fill::zeros);              // (N by M) matrix representing [X_1\thetaStar_1,...,X_M\thetaStar_M]
  arma::mat XthetaStar_new(points,M,arma::fill::zeros);     // (points by M) matrix of Xnew theta* based on gridpoints (NEW values)
  arma::mat Sigma_inv(N,N,arma::fill::zeros);               // Construct matrix Sigma=I+\lambda^{-1}*K where K is the kernel matrix
  arma::mat lamInv_Knew(points,points,arma::fill::zeros);   // lambda^{-1} times kernel matrix for new data
  arma::mat lamInv_Knewold(points,N,arma::fill::zeros);     // lambda^{-1} times kernel matrix for new AND old data (off-diagonal)
  arma::vec hmeantemp(points,arma::fill::zeros);            // temporary vector for computing marginal covariance of h(.)
  arma::vec hmean(points,arma::fill::zeros);                // marginal mean vector for h(.)
  arma::mat hcovtemp(points,points,arma::fill::zeros);      // marginal variance-covariance matrix of h(.)
  arma::mat hcov(points,points,arma::fill::zeros);          // marginal variance-covariance matrix of h(.)
  arma::mat hsamp(S_iter,points,arma::fill::zeros);         // matrix to store posterior draws of h: rows are samples of h
  arma::mat svdU;                                           // for SVD decomposition
  arma::vec svdD;                                           // for SVD decomposition
  arma::mat svdV;                                           // for SVD decomposition
  arma::vec ons(N, arma::fill::ones);                       // vector of length N of ones
  arma::field<arma::mat> Kmats;                             // field to contain output of get_Kmat()
  
  for(int s=0; s<S_iter; s++){                            // loop over posterior draws
    
    // residual vector (Y-Z\gamma)
    y_zgam = yz.col(0) - (yz.cols(1,P_z) * gamma.row(s).t()); // N-vector for Y-Z\gamma at the sth iteration
    
    // X\theta*
    for(int m=0; m<M; m++){    
      XthetaStar.col(m) =   X[m] * weights[m].row(s).t(); // set columns XthetaStar to Xm times theta_m* at the sth iteration
    }
    
    // Construct Xthetastar_new with gridpoints 
    for(int m=0; m<M; m++){    
      XthetaStar_new.col(m) =   gridpoints[m] * weights[m].row(s).t(); // set columns XthetaStar_new to Xm times theta_m* at the sth iteration
    }
    
    // for(int m=0; m<M; m++){             // loop over indices 
    //   for(int l=0; l<Lvec(m); l++){     // loop over components    
      
      
      
    // Construct kernel matrices
    Kmats = get_Kmat(poly, d, N, points, log(lambdaInverse[s]), XthetaStar, XthetaStar_new, randint, log(lambdaBInverse[s]), Bmat);
    Sigma_inv = inv(Kmats[0]);      // Construct Sigma^{-1}
    lamInv_Knew = Kmats[1];         // Construct lambda^{-1} times Kernel matrix for new values
    lamInv_Knewold = Kmats[2];      // Construct lambda^{-1} times Kernel matrix for new and old values
    
    
    // Predict hnew
    hmeantemp = lamInv_Knewold * Sigma_inv * y_zgam;                                          // temporary vector to be used in mean and covariance
    hcovtemp = sigma2[s] * (lamInv_Knew - lamInv_Knewold * Sigma_inv * lamInv_Knewold.t());   // temporary covariance matrix to be used in marginal covariance
    // compute marginal mean and covariance for h                                                                            
    hmean += hmeantemp/(S_iter);                                                              // law of total expecation, E[E[h|s]]
    hcov += hcovtemp/(S_iter) + hmeantemp*hmeantemp.t()/(S_iter);                             // first two components of law of total variance: E(Var(h|s))+E[E[h|s]^2]-E[E[h|s]]^2
    // draw h
    arma::svd(svdU,svdD,svdV,hcovtemp,"std");         // ADDED STD TO AVOID ERRORS            // SVD decomposition for covariance of h (since it's more stable than cholesky)                       
    hsamp.row(s) = ( hmeantemp+(svdU*diagmat(svdD))*arma::randn(points) ).t();                // hsampkeep.row(s_outer) = ( hmeantemp+arma::chol(hcovtemp)*arma::randn(N) ).t();
    
   
    
  } // end s loop over posterior draws
  
  
  
  
  // return a list
  return List::create(Named("hmean") = hmean,                             // mean of hnew 
                      Named("hcov") = hcov - ( hmean*hmean.t()),          // subtracting third component of law of total variance (i.e. E[E[h|s]]^2)
                      Named("hsamp") = hsamp );                            // predicted hnew
  
  
}
