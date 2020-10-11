# chapter 1 
setwd("~/Documents/sem2rcode/JAGSmodeltxtfiles")

# set seed
set.seed(7)

# generating simulations for the Brody model
time <- seq(0, 30, by =1)

BrodySim1 <- 0.2 * time + c(0, cumsum(rnorm(n=30)))
BrodySim2 <- 0.2 * time + c(0, cumsum(rnorm(n=30)))
BrodySim3 <- 0.2 * time + c(0, cumsum(rnorm(n=30)))

# function that returns posterior for Brody model
posterior <- function(xi, t, X, sigma, p){
  numerator = p * exp(sigma*X*xi - 0.5*sigma^2*X*t)
  denominator = (1-p) + p * exp(sigma*xi - 0.5*sigma^2*t)
  probability = numerator / denominator
  return(probability)
}

BrodyPost1 <- posterior(BrodySim1, time, 1, 0.2, 0.5)
BrodyPost2 <- posterior(BrodySim2, time, 1, 0.2, 0.5)
BrodyPost3 <- posterior(BrodySim3, time, 1, 0.2, 0.5)

# function that adds the fake noise term proposed by Brody
Brody_fake_news <- function(time, alpha, beta, tau){
  N <- length(time)
  fake_news_adj <- numeric(N)
  for(i in 1:N){
    if(time[i] >= tau){
      fake_news_adj[i] <- alpha*(time[i]-tau)*
                          exp(-beta*(time[i]-tau))
    }
  }
  return(fake_news_adj)
}

BrodyFake1 <- Brody_fake_news(time, 1, 0.05, 1)
BrodyFake2 <- Brody_fake_news(time, 2, 1, 25)
BrodyFake3 <- Brody_fake_news(time, 5, 0.5, 15)

obs1 <- BrodySim1 + BrodyFake1
obs2 <- BrodySim2 + BrodyFake2
obs3 <- BrodySim3 + BrodyFake3

BrodyAdjPost1 <- posterior(obs1, time, 1, 0.2, 0.5)
BrodyAdjPost2 <- posterior(obs2, time, 1, 0.2, 0.5)
BrodyAdjPost3 <- posterior(obs3, time, 1, 0.2, 0.5)
