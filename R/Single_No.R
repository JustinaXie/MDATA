#' Run Single-Trait without Annotation Model
#'
#' @param data A dataframe with columns: Gene name, mutability, de novo mutation count for traits
#' @param dn_col A character indicating the column name of de novo mutation counts for the trait
#' @param N_1 Cohort size for the trait
#' @param pi_init Initial value for risk gene proportion. Default=0.1
#' @param threshold Threshold for EM algorithm. Default=1e-6
#' @param maxiter Maximum number of iterations. Default=5e4
#' @return The estimated model parameters and the posterior probabilities of genes under different assumptions
#' \item{result}{A dataframe that includes estimated posterior probabilities of risk genes for each trait and estimated posterior probability for shared risk gene}
#' \item{pi}{Estimated proportion vector, the second value represents risk gene proportion}
#' \item{beta}{Estimated log of gamma for risk genes}
#' \item{Z_mat}{Estimated posterior probabilities of genes under different assumptions}
#' @import stats
#' @export
#'
Single_No<-function(data,dn_col,N_1,
                    pi_init=0.1,
                    threshold=1e-6,
                    maxiter=50000){
  dnm<-data

  #number of gene
  P <- dim(dnm)[1]
  rownames(dnm) <- 1:P

  #mutability of genes
  mu  <- dnm$mut

  #count for the trait
  Y_1 <- dnm[,dn_col]

  #Initialization scheme: fixed approach
  pi <- c(1-pi_init,pi_init)
  pi_old <- rep(0, 2)

  # Initialization of gamma
  beta_0 <- 99999
  beta_0_new <- 0.1

  #first round estimation
  #####loop start######
  k=0
  while(sum(abs(pi-pi_old)) > threshold)
  {
    print(paste("round", k+1, "start!"))

    prob_1 <- matrix(NA, P, 2)
    prob_1[,1] <- dpois(Y_1, 2*N_1*mu)
    prob_1[,2] <- dpois(Y_1, 2*N_1*mu*exp(beta_0_new))

    prob <- prob_1

    #update z (gene i)
    prob_weighted <- prob %*% diag(pi)
    Z <- prob_weighted / apply(prob_weighted, 1, sum)

    #update pi
    pi_old <- pi
    pi<- apply(Z, 2, sum) / P

    #update beta
    beta_0 <- beta_0_new
    beta_0_new <- log(sum(Z[,2] * Y_1) / sum(Z[,2] * 2*N_1*mu))

    k <- k + 1
    if(k>maxiter) {
      print(paste0("Over ",maxiter," iterations"))
      return(0)
    }
  }
  print("Complete!")
  result<-data.frame(Gene=dnm$Gene,
                     dn_count=Y_1,
                     Z=Z[,2])
  return(list(result=result,pi=pi,beta=beta_0,Z_mat=Z))
}
