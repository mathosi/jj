#' Prime factorization
#'
#' Returns the prime numbers as vector
#'
#' @param x single number, which is to be split into prime numbers
#' @keywords prime
#' @export
#' @examples
#' jj_prime_factorization(45)


jj_prime_factorization=function(x){
  if(length(x) > 1){
    warning('Only first entry of x will be used')
    x = x[1]
  }
  stopifnot(is.numeric(x))
  n=c()
  i=2
  r=x
  while(prod(n)!=x){
    if(!r%%i) {n=c(n,i);r=r/i;i=1}
    i=i+1
  }
  n
}
