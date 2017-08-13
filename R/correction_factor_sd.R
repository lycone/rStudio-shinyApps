correction.factor.sd <- function(n){
  k <- floor(n/2)
  if ((n%%2)==0){
    return(1/(sqrt(2/(pi*(2*k-1)))*((2^(2*k-2))*(factorial(k-1))^2)/factorial(2*k-2)))
  }
  return(1/(sqrt(pi/k)*(factorial(2*k-1))/((2^(2*k-1))*(factorial(k-1))^2)))
}
