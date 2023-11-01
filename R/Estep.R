Estep <- function(x, omega, mu, sigma)
{
  n <- length(x); k <- length(omega)
  gamma <- matrix(0,n,k)
  for(j in 1:n)
  {
    den=0
    for(i in 1:k)
    {
      gamma[j,i] <- omega[i]*dnorm(x[j],mu[i],sigma[i])
      den        <- den+omega[i]*dnorm(x[j],mu[i],sigma[i])
    }
    gamma[j,]  <- gamma[j,]/den
  }
  return(gamma)
}