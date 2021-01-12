#' Run Single-Trait with Annotation Model
#'
#' @param data A dataframe with columns: Gene name, mutability, de novo mutation count for traits
#' @param dn_col A character indicating the column name of de novo mutation counts for the trait
#' @param Anno Annotation file for the trait
#' @param N_1 Cohort size for the trait
#' @param pi_init Initial value for risk gene proportion. Default=0.1
#' @param threshold_1 Threshold for EM algorithm. Default=1e-3
#' @param threshold_2 Threshold for Newton's method. Default=1e-3
#' @param max_iter Maximum iteration for Newton's method. Default=200
#' @return The estimated model parameters and the posterior probabilities of genes under different assumptions
#' \item{result}{A dataframe that includes estimated posterior probabilities of risk genes for each trait and estimated posterior probability for shared risk gene}
#' \item{pi}{Estimated proportion vector, the second value represents risk gene proportion}
#' \item{beta}{Estimated beta vector}
#' \item{Z_mat}{Estimated posterior probabilities of genes under different assumptions}
#' @import stats
#' @export
#'
Single_Anno<-function(data,dn_col,N_1,Anno,
                    pi_init=0.1,beta0_init,
                    threshold_1=1e-3,threshold_2=1e-3,
                    max_iter=200){
  dnm<-data

  #number of gene
  P <- dim(dnm)[1]
  rownames(dnm) <- 1:P

  #mutability of genes
  mu  <- dnm$mut

  #annotation
  X_cand<-as.matrix(cbind(1,Anno))
  X_scale<-apply(X_cand[,-1], 2, scale)
  X<-as.matrix(cbind(1,X_scale))

  Q<-dim(X)[2]

  #count for the trait
  Y_1 <- dnm[,dn_col]

  Z <- matrix(NA, P, 2)

  #prior probability
  pi_old <- rep(0, 2)
  pi<-c(1-pi_init,pi_init)

  #effect sizes of annotation
  beta_1 <- beta_1_old<-rep(99999, Q)
  beta_1_new <- c(beta0_init,rep(0,Q-1))

  k=0
  while(sum(abs(pi-pi_old))+sum(abs(beta_1_new-beta_1_old)) > threshold_1)
  {
    print(paste("round", k+1, "start!"))

    prob_1 <- matrix(NA, P, 2)
    prob_1[,1] <- dpois(Y_1, 2*N_1*mu)
    prob_1[,2] <- dpois(Y_1, 2*N_1*mu*exp(X %*% beta_1_new))

    prob <- prob_1

    #update z (gene i)
    prob_weighted <- prob %*% diag(pi)
    Z <- prob_weighted / apply(prob_weighted, 1, sum)

    #update pi
    pi_old <- pi
    pi <- apply(Z, 2, sum) / P

    #update beta
    m <- 0
    beta_1_old<-beta_1_new
    while(sum(abs(beta_1_new - beta_1)) > threshold_2 & m < max_iter)
    {
      # print(paste("beta1:", beta_1_new))
      d_beta_1 <- t(Y_1 * X - as.vector(2*N_1*mu*exp(X %*% beta_1_new)) * X) %*% Z[, 2]
      dd_beta_1 <- -t(X) %*% diag(Z[, 2] * as.vector(2*N_1*mu*exp(X %*% beta_1_new))) %*% X

      beta_1 <- beta_1_new
      beta_1_new <- beta_1 - solve(dd_beta_1) %*% d_beta_1
      m <- m + 1
    }
    k <- k + 1
  }
  print("Complete!")
  result<-data.frame(Gene=dnm$Gene,
                     dn_count=Y_1,
                     Z=Z[,2])
  return(list(result=result,pi=pi,beta=beta_1,Z_mat=Z))
}
