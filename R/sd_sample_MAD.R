sd_sample_MAD <- function (input_table, columns, use.mad){
  

source("R/estimate_sd.R")
source("R/correction_factor_mad.R")
source("R/correction_factor_sd.R")


dmat <- input_table[columns]

estimateSd <- estimate.sd(dmat, use.mad)

return (estimateSd)

}




