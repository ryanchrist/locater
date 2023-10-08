
# Variance Normalization does not fix Gaussian outliers

n_reps <- 1e3
y <- matrix(0,n_reps,2)

for(i in 1:n_reps){
  xx <- c(rnorm(1e3),6)
  y[i,1] <- -pnorm(6,sd = sd(xx),lower.tail = FALSE,log.p = TRUE)/log(10)
  y[i,2] <- -pnorm(max(locater:::rank2gauss(rank(xx))),lower.tail = FALSE,log.p = TRUE)/log(10)
print(i)
}



plot(ecdf(y[,1]))
abline(v=min(y[,2]),col="red")
abline(v=max(y[,2]),col="red")

# logistic <- function(x){1/(1+exp(-x))}
# logit <- function(p){log(p/(1-p))}
# xx <- seq(-2,2,len=1e3)
# plot(xx, logistic(log(5) * sqrt(5)* xx))
#
# mean(logistic(logit(0.0015) + log(5) * sqrt(5)*rnorm(1e5)))
