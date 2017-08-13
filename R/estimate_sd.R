#' Estimation of the standard deviation
#'
#' The function is able to estimate the standard deviation based on selected columns from the input table in two different ways. The first one is based on the MAD and the second one on the standard deviation of the single dataset values. Each estimation uses specific correction factors.
#' @param dmat Dataframe columns from that the standard deviation should be computed.
#' @param use.mad Logical value. TRUE returns the estimation of the standard deviation based on the MAD. FALSE returns the estimation of the standard deviation based on the standard deviation of the single dataset values.


estimate.sd <- function(dmat,use.mad){
  mads <- c()
  
  if (use.mad){
    # mads <- c()
    for (i in 1:nrow(dmat)){
      mads <- c(mads,mad(t(dmat[i,]),constant=1, na.rm = T))
    }
    return(mean(mads,na.rm=T)*correction.factor.mad(ncol(dmat)))
  }

  # mads <- c()
  for (i in 1:nrow(dmat)){
    mads <- c(mads,sd(t(dmat[i,]), na.rm = T))
  }
  return(mean(mads,na.rm=T)*correction.factor.sd(ncol(dmat)))
}




