EM <- function(x, k, er, q)
{
  omega <- c(rep(1/k,k)); mu <- numeric(k); sigma <- numeric(k)
  for(i in 1:k)
  {
    mu[i] <- quantile(x,i/(k+1))
  }
  sigma <- c(rep(sd(x),k)); dist <- 100
  while(dist>=er)
  {
    gamma <- Estep(x, omega, mu, sigma)
    fit   <- Mstep(x, gamma, q)
    omega <- fit$omega; mu1 <- fit$mu; sigma <- fit$sigma
    dist <- sum((mu-mu1)^2); mu <- mu1
  }
  fit <- list("omega"=omega,"mu"=mu, "sigma"=sigma, "gamma"=gamma)
  return(fit)
}