#' Get significant gene list from Multi_Trait Models
#'
#' @param output A list object from the output of multi-trait models
#' @param threshold Threshold for FDR. Default=0.05
#' @return The significant gene list for two traits under FDR<threshold
#' \item{Gene_trait1}{Significant gene list for trait 1}
#' \item{Gene_trait2}{Significant gene list for trait 2}
#' @import stats
#' @export
#'
Get_Multi_Gene<-function(output,threshold=0.05){
  result<-output[["result"]]
  Z<-output[["Z_mat"]]
  P<-length(result$Gene)
  lfdr<-data.frame(index=1:P,Z_1=Z[,1]+Z[,3],Z_2=Z[,1]+Z[,2])
  lfdr_order_trait1<-lfdr[order(lfdr$Z_1,decreasing = FALSE),]
  lfdr_order_trait1$FDR_trait1<-cumsum(lfdr_order_trait1$Z_1)/(1:P)
  tmp_trait1<-lfdr_order_trait1[order(lfdr_order_trait1$index,decreasing = FALSE),]
  tmp_FDR_trait1<-tmp_trait1$FDR_trait1
  Gene_trait1=result$Gene[which(tmp_FDR_trait1<threshold)]

  lfdr_order_trait2<-lfdr[order(lfdr$Z_2,decreasing = FALSE),]
  lfdr_order_trait2$FDR_trait2<-cumsum(lfdr_order_trait2$Z_2)/(1:P)
  tmp_trait2<-lfdr_order_trait2[order(lfdr_order_trait2$index,decreasing = FALSE),]
  tmp_FDR_trait2<-tmp_trait2$FDR_trait2
  Gene_trait2=result$Gene[which(tmp_FDR_trait2<threshold)]

  return(list(Gene_trait1=Gene_trait1,Gene_trait2=Gene_trait2))
}

