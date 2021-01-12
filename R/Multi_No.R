#' Run Multi-Trait without Annotation Model
#'
#' @param data A dataframe with columns: Gene name, mutability, de novo mutation count for trait1, de novo mutation count for trait2
#' @param N_1 Cohort size for trait 1
#' @param N_2 Cohort size for trait 2
#' @param pi_init Initial value for probabilities of genes under different assumptions. Default=c(0.85,0.05,0.05,0.05)
#' @param threshold Threshold for EM algorithm. Default=1e-6
#' @return The estimated model parameters and the posterior probabilities of genes under different assumptions
#' \item{result}{A dataframe that includes estimated posterior probabilities of risk genes for each trait and estimated posterior probability for shared risk gene}
#' \item{pi}{Estimated proportion of genes under different assumptions}
#' \item{beta_trait1}{Estimated log of gamma for trait 1}
#' \item{beta_trait2}{Estimated log of gamma for trait 2}
#' \item{Z_mat}{Estimated posterior probabilities of genes under different assumptions}
#' @import stats
#' @export
#'

Multi_No<-function(data,N_1,N_2,
                   pi_init=c(0.85,0.05,0.05,0.05),
                   threshold=1e-6){
  dnm<-data

  #number of gene
  P <- dim(dnm)[1]
  rownames(dnm) <- 1:P

  #mutability of genes
  mu  <- dnm$mut

  #count for trait 1
  Y_1 <- dnm$dn_trait1
  #count for trait 2
  Y_2 <- dnm$dn_trait2

  #indicator
  Z <- matrix(NA, P, 4)

  #prior probability
  pi <- pi_init
  pi_old <- rep(0, 4)

  #effect sizes of annotation
  beta0_1 <- 99999
  beta0_1_new <- 0.1
  beta0_2 <- 99999
  beta0_2_new <- 0.1

  k=0
  while(sum(abs(pi-pi_old)) > threshold)
  {
    print(paste("round", k+1, "start!"))
    prob_1 <- matrix(NA, P, 4)
    prob_1[,1] <- dpois(Y_1, 2*N_1*mu)*dpois(Y_2, 2*N_2*mu)
    prob_1[,2] <- dpois(Y_1, 2*N_1*mu*exp(beta0_1_new))*dpois(Y_2, 2*N_2*mu)
    prob_1[,3] <- dpois(Y_1, 2*N_1*mu)*dpois(Y_2, 2*N_2*mu*exp(beta0_2_new))
    prob_1[,4] <- dpois(Y_1, 2*N_1*mu*exp(beta0_1_new))*dpois(Y_2, 2*N_2*mu*exp(beta0_2_new))
    prob <- prob_1

    #update z (gene i)
    prob_weighted <- prob %*% diag(pi)
    Z <- prob_weighted / apply(prob_weighted, 1, sum)

    #update pi
    pi_old <- pi
    pi <- apply(Z, 2, sum) / P

    #update beta0_1
    beta0_1 <- beta0_1_new
    beta0_1_new <- log(sum((Z[,2]+Z[,4]) * Y_1) / sum((Z[,2]+Z[,4]) * 2*N_1*mu))

    #update beta0_2
    beta0_2 <- beta0_2_new
    beta0_2_new <- log(sum((Z[,3]+Z[,4]) * Y_2) / sum((Z[,3]+Z[,4]) * 2*N_2*mu))

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
  return(list(result=result,pi=pi,beta_trait1=beta0_1,beta_trait2=beta0_2,Z_mat=Z))
}
