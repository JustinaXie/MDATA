#' Get significant gene list from single-trait models
#'
#' @param output A list object from the output of single-trait models
#' @param threshold Threshold for FDR. Default=0.05
#' @return The significant gene vector for the trait under FDR<threshold
#' \item{Gene_trait}{Significant gene list for the trait}
#' @import stats
#' @export
#'
Get_Single_Gene<-function(output,threshold=0.05){
  result<-output[["result"]]
  Z<-output[["Z_mat"]]
  P<-length(result$Gene)

  lfdr<-data.frame(index=1:P,Z=Z[,1])
  lfdr_order<-lfdr[order(lfdr$Z,decreasing = FALSE),]
  lfdr_order$FDR<-cumsum(lfdr_order$Z)/(1:P)
  tmp<-lfdr_order[order(lfdr_order$index,decreasing = FALSE),]
  tmp_FDR<-tmp$FDR
  result$FDR<-tmp_FDR
  return(list(Gene=result$Gene[which(tmp_FDR<0.05)],Output=result))
}

