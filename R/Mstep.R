Mstep <- function(x, gamma, q)
{
  n <- dim(gamma)[1]; k <- dim(gamma)[2]
  omega <- apply(gamma,2,sum)/n
  mu <- numeric(k); sigma <- numeric(k)
  #for(i in 1:k)
  #{
  #  mu[i] <- sum(gamma[,i]*x)/sum(gamma[,i])
  #  sigma[i]=(sum(gamma[,i]*(x-mu[i])^2)/sum(gamma[,i]))^.5
  #}
  mu[1]    <- mean(x); mu[2]    <- mean(x)/q
  sigma[1] <- sd(x);   sigma[2] <- sd(x)/q
  return(list("omega"=omega,"mu"=mu, "sigma"=sigma))
}