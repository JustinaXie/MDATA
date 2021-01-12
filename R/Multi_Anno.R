#' Run Multi-Trait with Annotation Model
#'
#' @param data A dataframe with columns: Gene name, mutability, de novo mutation count for trait1, de novo mutation count for trait2
#' @param Anno1 Annotation file for trait 1
#' @param Anno2 Annotation file for trait 2
#' @param N_1 Cohort size for trait 1
#' @param N_2 Cohort size for trait 2
#' @param pi_init Initial value for probabilities of genes under different assumptions. Default=c(0.85,0.05,0.05,0.05)
#' @param beta0_trait1_init Initial value for log gamma for trait 1
#' @param beta0_trait2_init Initial value for log gamma for trait 2
#' @param threshold_1 Threshold for EM algorithm. Default=1e-3
#' @param threshold_2 Threshold for Newton's method. Default=1e-3
#' @param max_iter Maximum iteration for Newton's method. Default=200
#' @return The estimated model parameters and the posterior probabilities of genes under different assumptions
#' \item{result}{A dataframe that includes estimated posterior probabilities of risk genes for each trait and estimated posterior probability for shared risk gene}
#' \item{pi}{Estimated proportion of genes under different assumptions}
#' \item{beta_trait1}{Estimated beta vector for trait 1}
#' \item{beta_trait2}{Estimated beta vector for trait 2}
#' \item{Z_mat}{Estimated posterior probabilities of genes under different assumptions}
#' @import stats
#' @export
#'
Multi_Anno<-function(data,Anno1,Anno2,N_1,N_2,
                     pi_init=c(0.85,0.05,0.05,0.05),
                     beta0_trait1_init,beta0_trait2_init,
                     threshold_1=1e-3,threshold_2=1e-3,
                     max_iter=200
                     ){
  #number of gene
  P <- dim(dnm)[1]
  rownames(dnm) <- 1:P

  #mutability of genes
  mu  <- dnm$mut

  #count for trait 1
  Y_1 <- dnm$dn_trait1
  #count for trait 2
  Y_2 <- dnm$dn_trait2

  X1_cand<-as.matrix(cbind(1,Anno1))
  X2_cand<-as.matrix(cbind(1,Anno2))

  X1_scale<-apply(X1_cand[,-1], 2, scale)
  X1_scale<-as.matrix(cbind(1,X1_scale))
  X2_scale<-apply(X2_cand[,-1], 2, scale)
  X2_scale<-as.matrix(cbind(1,X2_scale))

  X_1<-X1_scale
  X_2<-X2_scale
  Q_1<-dim(X_1)[2]
  Q_2<-dim(X_2)[2]

  #indicator
  Z <- matrix(NA, P, 4)

  #prior probability
  pi <- pi_init
  pi_old <- rep(0, 4)

  #effect sizes of annotation
  beta_1 <- beta_1_old<-rep(99999, Q_1)
  beta_1_new <- c(beta0_trait1_init,rep(0,Q_1-1))
  beta_2 <- beta_2_old<-rep(99999, Q_2)
  beta_2_new <- c(beta0_trait2_init,rep(0,Q_2-1))

  k=0
  while(sum(abs(pi-pi_old))+sum(abs(beta_1-beta_1_old))+sum(abs(beta_2_new-beta_2_old)) > threshold_1)
  {
    print(paste("round", k+1, "start!"))
    prob_1 <- matrix(NA, P, 4)
    prob_1[,1] <- dpois(Y_1, 2*N_1*mu)*dpois(Y_2, 2*N_2*mu)
    prob_1[,2] <- dpois(Y_1, 2*N_1*mu*exp(X_1 %*% beta_1_new))*dpois(Y_2, 2*N_2*mu)
    prob_1[,3] <- dpois(Y_1, 2*N_1*mu)*dpois(Y_2, 2*N_2*mu*exp(X_2 %*% beta_2_new))
    prob_1[,4] <- dpois(Y_1, 2*N_1*mu*exp(X_1 %*% beta_1_new))*dpois(Y_2, 2*N_2*mu*exp(X_2 %*% beta_2_new))
    prob <- prob_1

    #update z (gene i)
    prob_weighted <- prob %*% diag(pi)
    Z <- prob_weighted / apply(prob_weighted, 1, sum)

    #update pi_no
    pi_old <- pi
    pi <- apply(Z, 2, sum) / P

    #update beta_1
    m1 <- 0
    beta_1_old<-beta_1_new
    while(sum(abs(beta_1_new - beta_1))> threshold_2 & m1 < max_iter)
    {
      d_beta_1 <- t(Y_1 * X_1 - as.vector(2*N_1*mu*exp(X_1%*%beta_1_new)) * X_1) %*% (Z[, 2]+Z[,4])
      dd_beta_1 <- -t(X_1) %*% diag((Z[, 2]+Z[, 4]) * as.vector(2*N_1*mu*exp(X_1%*%beta_1_new))) %*% X_1

      beta_1 <- beta_1_new
      beta_1_new <- beta_1 - solve(dd_beta_1) %*% d_beta_1
      m1 <- m1 + 1
      # print(paste("m1:",m1))
    }

    #update beta_2
    m2 <- 0
    beta_2_old<-beta_2_new
    while(sum(abs(beta_2_new - beta_2))> threshold_2 & m2 < max_iter)
    {
      d_beta_2 <- t(Y_2 * X_2 - as.vector(2*N_2*mu*exp(X_2%*%beta_2_new)) * X_2) %*% (Z[, 3]+Z[,4])
      dd_beta_2 <- -t(X_2) %*% diag((Z[, 3]+Z[, 4]) * as.vector(2*N_2*mu*exp(X_2%*%beta_2_new))) %*% X_2

      beta_2 <- beta_2_new
      beta_2_new <- beta_2 - solve(dd_beta_2) %*% d_beta_2
      m2 <- m2 + 1
      # print(paste("m2:",m2))
    }
    k <- k + 1
  }
  print("Complete!")
  result<-data.frame(Gene=dnm$Gene,
                     dn_trait1=Y_1,
                     dn_trait2=Y_2,
                     dn_both=pmin(Y_1,Y_2),
                     Z_trait1=Z[,2]+Z[,4],
                     Z_trait2=Z[,3]+Z[,4],
                     Z_both=Z[,4])
  return(list(result=result,pi=pi,beta_trait1=beta_1,beta_trait2=beta_2,Z_mat=Z))
}
