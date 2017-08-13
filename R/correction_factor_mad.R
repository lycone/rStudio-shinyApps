correction.factor.mad <- function(n){
  if (n<11){
    cf <- c(NA,1.7744,2.2053,2.0164,1.8056,1.7628,1.6876,1.6717,1.6326,1.6254)
    return(cf[n])
  }
  if (n<30){return(1.5454)}
  if (n<40){return(1.5228)}
  if (n<60){return(1.5129)}
  if (n<80){return(1.5017)}
  if (n<100){return(1.4970)}
  if (n<120){return(1.4942)}
  if (n<150){return(1.4922)}
  if (n<200){return(1.4903)}
  if (n<300){return(1.4884)}
  if (n<500){return(1.4864)}
  if (n<1000){return(1.4849)}
  return(1.4836)
}
