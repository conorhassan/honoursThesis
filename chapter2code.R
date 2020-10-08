# loading required packages
library(rjags)
library(MCMCpack)
library(RColorBrewer)
library(latex2exp)
library(GISTools)

# Chapter 2 code

# function regarding mu when sigma is known
sigma_known <- function(mean_prior, var_prior, var_known, sample_mean, n_obs){
  dist_mean = ((mean_prior/var_prior) + ((n_obs*sample_mean)/var_known)) / 
              ((1/var_prior) + (n_obs/var_known))
  dist_var = 1 / ((1/var_prior) +(n_obs/var_known))
  prob_great_0 = pnorm(q=0, mean=dist_mean, sd=sqrt(dist_var), lower.tail=FALSE)
  cred_int = qnorm(p=c(0.025, 0.5, 0.975), mean=dist_mean, sd=sqrt(dist_var))
  return(list(posterior_prob=prob_great_0,
              credible_interval=cred_int))
}

# helper function for when sigma is unknown
find_mu_probability <- function(x, sample_mean, sample_sd, n){
  sample_mean + (sample_sd / sqrt(n)) * qt(x, n-1, lower.tail=FALSE)
}

# function regarding mu when sigma unknown
sigma_unknown <- function(data_vector){
  sample_mean = mean(data_vector)
  sample_sd = sd(data_vector)
  n_obs = length(data_vector)
  posterior_prob = uniroot(find_mu_probability, 
                           interval = c(0, 1),
                           sample_mean = sample_mean,
                           sample_sd = sample_sd,
                           n = n_obs)$root
  credible_interval <- find_mu_probability(c(0.025, 0.5, 0.975), 
                                           sample_mean,
                                           sample_sd, 
                                           n_obs)
  return(list(posterior_probability=posterior_prob,
              credible_interval=credible_interval))
} 

# setting mu and the two sigmas
set.seed(7)
mu1 <- 0.5
sd1 <- 1
sd2 <- 3

# plots for the strong signal
#xi, mu posterior with var known and unknown
plot(seq(0, 30), mu1*seq(0, 30), 
     type="l",
     col=coul[1],
     bty="n",
     xlab="Time",
     ylab=expression(xi[t]),
     ylim=c(-30, 55),
     xaxt = "n",
     yaxt= "n")
axis(2, cex.axis=0.7)
axis(1, cex.axis=0.8)
for(i in 1:30){
  observed_sequence <- c(0, cumsum(rnorm(30, mu1, sd1)))
  lines(seq(0, 30), observed_sequence, col=coul.t[1])
}
expected_xi1 <- mu1*seq(1, 30)
post_expected_unknown1 <- numeric(30)
set.seed(8)
for(t in 1:30){
  post_expected_unknown1[t] <- sigma_known(0, 1000, sd1^2, mu1, t)$posterior_prob
}
plot(seq(1, 30), post_expected_unknown1, 
     ylim=c(0, 1), type="l", col=coul[2],
     bty="n",
     xlab="Time",
     ylab=expression(P(mu > 0)),
     xaxt = "n",
     yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  observed_info <- rnorm(30, mu1, sd1)
  post <- numeric(30)
  for(j in 1:30){
    post[j] <- sigma_known(0, 1000, sd1^2, mean(observed_info[1:j]), j)$posterior_prob
  }
  lines(seq(1, 30), post, col=coul.t[2])
}
set.seed(8)
plot(1, type="n", xlim=c(0, 30),
     ylim=c(0, 1),
     bty="n",
     xlab="Time",
     ylab=expression(P(mu > 0)),
     xaxt = "n",
     yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(m in 1:30){
  observed_info <- rnorm(30, mu1, sd1)
  prob_unknown <- numeric(29)
  for(k in 2:30){
    prob_unknown[(k-1)] <- sigma_unknown(observed_info[1:k])$posterior_prob
  }
  lines(seq(2, 30), prob_unknown, col=coul.t[3])
}

# plots for the noisy signal
#xi, mu posterior with var known and unknown
plot(seq(0, 30), mu1*seq(0, 30), 
     type="l",
     col=coul[1],
     bty="n",
     xlab="Time",
     ylab=expression(xi[t]),
     ylim=c(-30, 55),
     xaxt = "n",
     yaxt= "n")
axis(2, cex.axis=0.7)
axis(1, cex.axis=0.8)
for(i in 1:30){
  observed_sequence <- c(0, cumsum(rnorm(30, mu1, sd2)))
  lines(seq(0, 30), observed_sequence, col=coul.t[1])
}
expected_xi1 <- mu1*seq(1, 30)
post_expected_unknown1 <- numeric(30)
set.seed(8)
for(t in 1:30){
  post_expected_unknown1[t] <- sigma_known(0, 1000, sd2^2, mu1, t)$posterior_prob
}
plot(seq(1, 30), post_expected_unknown1, 
     ylim=c(0, 1), type="l", col=coul[2],
     bty="n",
     xlab="Time",
     ylab=expression(P(mu > 0)),
     xaxt = "n",
     yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  observed_info <- rnorm(30, mu1, sd2)
  post <- numeric(30)
  for(j in 1:30){
    post[j] <- sigma_known(0, 1000, sd2^2, mean(observed_info[1:j]), j)$posterior_prob
  }
  lines(seq(1, 30), post, col=coul.t[2])
}
set.seed(8)
plot(1, type="n", xlim=c(0, 30),
     ylim=c(0, 1),
     bty="n",
     xlab="Time",
     ylab=expression(P(mu > 0)),
     xaxt = "n",
     yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(m in 1:30){
  observed_info <- rnorm(30, mu1, sd2)
  prob_unknown <- numeric(29)
  for(k in 2:30){
    prob_unknown[(k-1)] <- sigma_unknown(observed_info[1:k])$posterior_prob
  }
  lines(seq(2, 30), prob_unknown, col=coul.t[3])
}

# functions to estimate sigma
rInvGamma <- function (n, shape, scale = 1) {
  return(1/rgamma(n = n, shape = shape, rate = scale))
}
rScaledInvChiSq <- function (n, nu, tau2) {
  return(rInvGamma(n, shape = nu / 2, scale = (nu * tau2) / 2))
}
# estimate of first sd
sd_estimate <- matrix(nrow=29, ncol=30)
for(j in 1:30){
  current_obs <- rnorm(30, mu1, sd1)
  for(i in 2:30){
    inv_sample <- rScaledInvChiSq(1000000, i-1, sd(current_obs[1:i]))
    sd_estimate[(i-1), j] <- quantile(inv_sample, probs = c(0.5))
  }
}
# estimate of second sd
sd_estimate2 <- matrix(nrow=29, ncol=30)
for(j in 1:30){
  current_obs <- rnorm(30, mu1, sd2)
  for(i in 2:30){
    inv_sample <- rScaledInvChiSq(1000000, i-1, sd(current_obs[1:i]))
    sd_estimate2[(i-1), j] <- quantile(inv_sample, probs = c(0.5))
  }
}
# plotting the two estimates of sigma
plot(1, type="n", xlim=c(0, 30), ylim=c(0, 10),
     bty="n", xlab="Time", ylab=TeX("$\\hat{\\sigma^2}$"),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), sd_estimate[, s], type="l", col=coul.t[2])
}
plot(1, type="n", xlim=c(0, 30), ylim=c(0, 10),
     bty="n", xlab="Time", ylab=TeX("$\\hat{\\sigma^2}$"),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), sd_estimate2[, s], type="l", col=coul.t[3])
}