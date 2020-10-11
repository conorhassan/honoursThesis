second_fake_effect = function(time_seq, A, beta, tau){
  effect <- numeric(length(time_seq))
  for(t in 1:length(time_seq)){
    if(time_seq[t] < tau){
      effect[t] <- 0
    }
    if(time_seq[t] >= tau){
      effect[t] <- A*exp(-beta*(time_seq[t]-tau))
    }
  }
  return(effect)
}
init_robust <- function(){
  return(list(.RNG.name="base::Super-Duper",
              .RNG.seed=7))
}

# effect of the fake news
fake_effect_1 <- second_fake_effect(1:30, -10, 0.5, 10)

# t=1 to t=30.... simulating from mu=0.5, sd=2
set.seed(7)

scen1.s1.m1.mu.post <- matrix(nrow=29, ncol=30)
scen1.s1.m2.mu.post <- matrix(nrow=29, ncol=30)
scen1.s1.m3.mu.post <- matrix(nrow=29, ncol=30)
for(s in 1:30){
  current_sim = rnorm(30, mean=0.5, sd=2) + fake_effect_1
  for(t in 2:30){
    # model 1 
    moment <- sigma_known(0, 1000, 4, mean(current_sim[1:t]), t)
    scen1.s1.m1.mu.post[t-1, s] <- moment$posterior_prob
    # model 2
    scen1.s1.m2.mu.post[t-1, s] <- sigma_unknown(current_sim[1:t])$posterior_prob
    # model 3
    robust.model = jags.model(file = "bias_model.txt",
                            data = list(y = current_sim[1:t],
                                        N = t), 
                            n.chains = 3,
                            inits = init_robust,
                            n.adapt = 10000)  
    robust.output = coda.samples(model=robust.model, 
                               variable.names=c("mu"),
                               n.iter=10000)
    robust.output.m <- as.matrix(robust.output)
    scen1.s1.m3.mu.post[t-1, s] = sum(robust.output.m[, 1] > 0) / 30000
  }
}

set.seed(7)
scen1.s1.z <- matrix(nrow=30, ncol=30)
scen1.s1.m4.mu.post <- matrix(nrow=29, ncol=30)
scen1.s1.m4.A <- matrix(nrow=29, ncol=90)
scen1.s1.m4.beta <- matrix(nrow=29, ncol=90)
scen1.s1.m5.mu.post <- matrix(nrow=29, ncol=30)
scen1.s1.m5.A <- matrix(nrow=29, ncol=90)
scen1.s1.m5.beta <- matrix(nrow=29, ncol=90)
scen1.s1.m5.fake <- matrix(nrow=29, ncol=90)
for(s in 1:30){
  current_sim = rnorm(30, mean=0.5, sd=2) + fake_effect_1
  scen1.s1.z[, s] <- current_sim
  for(t in 2:30){
      estimation.model.tau.k <- jags.model(file = "normal_tau_known.txt",
                                           data = list(y=current_sim[1:t],
                                                       T=t,
                                                       fake=10),
                                           n.chains = 3,
                                           inits = init_robust,
                                           n.adapt = 2500)
      estimation.output.tau.k <- coda.samples(model = estimation.model.tau.k,
                                              variable.names=c("mu",
                                                               "A",
                                                               "beta"),
                                              n.iter=10000)
      column_index = 3*(s-1) + c(1, 2, 3)
      estimation.output.tau.m <- as.matrix(estimation.output.tau.k)
      scen1.s1.m4.mu.post[(t-1), s] = sum(estimation.output.tau.m[, 3] > 0) / 30000
      scen1.s1.m4.A[(t-1), column_index] <- summary(estimation.output.tau.k)$quantiles[1, c(1, 3, 5)]
      scen1.s1.m4.beta[(t-1), column_index] <- summary(estimation.output.tau.k)$quantiles[2, c(1, 3, 5)]
      estimation.model = jags.model(file = "normal_estimation_model.txt",
                                    data = list(y=current_sim[1:t], 
                                                T=t), 
                                    n.chains = 3,
                                    inits = init_robust,
                                    n.adapt = 2500)  
      estimation.output = coda.samples(model=estimation.model, 
                                       variable.names=c("mu",
                                                        "A",
                                                        "beta", 
                                                        "fake"),
                                       n.iter=10000)
      estimation.output.m <- as.matrix(estimation.output)
      scen1.s1.m5.mu.post[(t-1), s] = sum(estimation.output.m[, 4] > 0) / 30000
      scen1.s1.m5.A[(t-1), column_index] <- summary(estimation.output)$quantiles[1, c(1, 3, 5)]
      scen1.s1.m5.beta[(t-1), column_index] <- summary(estimation.output)$quantiles[2, c(1, 3, 5)]
      scen1.s1.m5.fake[(t-1), column_index] <- summary(estimation.output)$quantiles[3, c(1, 3, 5)]
    }
  }
}

##### saving down data.frame from the last simulation above at t=30 to plot
## A vs. beta, beta vs fake, fake vs a
set.seed(7)
current_sim = rnorm(30, 0.5, 2) + fake_effect_1
current_sim = rnorm(30, 0.5, 2) + fake_effect_1
estimation.model = jags.model(file = "normal_estimation_model.txt",
                              data = list(y=current_sim, 
                                          T=30), 
                              n.chains = 3,
                              inits = init_robust,
                              n.adapt = 2500)  
estimation.output = coda.samples(model=estimation.model, 
                                 variable.names=c("mu",
                                                  "A",
                                                  "beta", 
                                                  "fake"),
                                 n.iter=10000)
estimation.output.S1.m <- as.matrix(estimation.output)
plot(estimation.output.S1.m[, 1], estimation.output.S1.m[, 2], xlab="A", ylab="beta")
plot(estimation.output.S1.m[, 1], estimation.output.S1.m[, 3], xlab="A", ylab="fake")
plot(estimation.output.S1.m[, 2], estimation.output.S1.m[, 3], xlab="beta", ylab="fake")

plot(estimation.output.S1.m[, 1], estimation.output.S1.m[, 2],
     bty="n", 
     xlab=TeX('$\\hat{\\alpha}$'), ylab=TeX('$\\hat{\\beta}$'),
     col=coul[4])
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)

####### 2nd Fake news effect with Dorje like error in CH
fake_effect_2 <- second_fake_effect(1:30, -20, 1, 10)

# t=1 to t=30.... simulating from mu=1, sd=2
set.seed(8)

# specifying data frames
scen1.s2.z <- matrix(nrow=30, ncol=30)
scen1.s2.m1.mu.post <- matrix(nrow=29, ncol=30)
scen1.s2.m2.mu.post <- matrix(nrow=29, ncol=30)
scen1.s2.m3.mu.post <- matrix(nrow=29, ncol=30)
scen1.s2.m4.mu.post <- matrix(nrow=29, ncol=30)
scen1.s2.m4.A <- matrix(nrow=29, ncol=90)
scen1.s2.m4.beta <- matrix(nrow=29, ncol=90)
scen1.s2.m5.mu.post <- matrix(nrow=29, ncol=30)
scen1.s2.m5.A <- matrix(nrow=29, ncol=90)
scen1.s2.m5.beta <- matrix(nrow=29, ncol=90)
scen1.s2.m5.fake <- matrix(nrow=29, ncol=90)


for(s in 1:30){
  current_sim = rnorm(30, mean=1, sd=2) + fake_effect_2
  scen1.s2.z[, s] <- current_sim
  for(t in 2:30){
    # model 1 
    moment <- sigma_known(0, 1000, 4, mean(current_sim[1:t]), t)
    scen1.s2.m1.mu.post[t-1, s] <- moment$posterior_prob
    # model 2
    scen1.s2.m2.mu.post[t-1, s] <- sigma_unknown(current_sim[1:t])$posterior_prob
    # model 3
    robust.model = jags.model(file = "bias_model.txt",
                              data = list(y = current_sim[1:t],
                                          N = t), 
                              n.chains = 3,
                              inits = init_robust,
                              n.adapt = 10000)  
    robust.output = coda.samples(model=robust.model, 
                                 variable.names=c("mu"),
                                 n.iter=10000)
    robust.output.m <- as.matrix(robust.output)
    scen1.s2.m3.mu.post[t-1, s] = sum(robust.output.m[, 1] > 0) / 30000
    estimation.model.tau.k <- jags.model(file = "normal_tau_known.txt",
                                         data = list(y=current_sim[1:t],
                                                     T=t,
                                                     fake=10),
                                         n.chains = 3,
                                         inits = init_robust,
                                         n.adapt = 2500)
    estimation.output.tau.k <- coda.samples(model = estimation.model.tau.k,
                                            variable.names=c("mu",
                                                             "A",
                                                             "beta"),
                                            n.iter=10000)
    column_index = 3*(s-1) + c(1, 2, 3)
    estimation.output.tau.m <- as.matrix(estimation.output.tau.k)
    scen1.s2.m4.mu.post[(t-1), s] = sum(estimation.output.tau.m[, 3] > 0) / 30000
    scen1.s2.m4.A[(t-1), column_index] <- summary(estimation.output.tau.k)$quantiles[1, c(1, 3, 5)]
    scen1.s2.m4.beta[(t-1), column_index] <- summary(estimation.output.tau.k)$quantiles[2, c(1, 3, 5)]
    estimation.model = jags.model(file = "normal_estimation_model.txt",
                                  data = list(y=current_sim[1:t], 
                                              T=t), 
                                  n.chains = 3,
                                  inits = init_robust,
                                  n.adapt = 2500)  
    estimation.output = coda.samples(model=estimation.model, 
                                     variable.names=c("mu",
                                                      "A",
                                                      "beta", 
                                                      "fake"),
                                     n.iter=10000)
    estimation.output.m <- as.matrix(estimation.output)
    scen1.s2.m5.mu.post[(t-1), s] = sum(estimation.output.m[, 4] > 0) / 30000
    scen1.s2.m5.A[(t-1), column_index] <- summary(estimation.output)$quantiles[1, c(1, 3, 5)]
    scen1.s2.m5.beta[(t-1), column_index] <- summary(estimation.output)$quantiles[2, c(1, 3, 5)]
    scen1.s2.m5.fake[(t-1), column_index] <- summary(estimation.output)$quantiles[3, c(1, 3, 5)]
  }
}

##### saving down data.frame from the last simulation above at t=30 to plot
## A vs. beta, beta vs fake, fake vs a
set.seed(8)
current_sim = rnorm(30, 1, 2) + fake_effect_2
#current_sim = rnorm(30, 0.5, 2) + fake_effect_1
estimation.model = jags.model(file = "normal_estimation_model.txt",
                              data = list(y=current_sim, 
                                          T=30), 
                              n.chains = 3,
                              inits = init_robust,
                              n.adapt = 2500)  
estimation.output = coda.samples(model=estimation.model, 
                                 variable.names=c("mu",
                                                  "A",
                                                  "beta", 
                                                  "fake"),
                                 n.iter=10000)
estimation.output.S2.m <- as.matrix(estimation.output)
plot(estimation.output.S2.m[, 1], estimation.output.S2.m[, 2], xlab="A", ylab="beta")
plot(estimation.output.S2.m[, 1], estimation.output.S2.m[, 3], xlab="A", ylab="fake")
plot(estimation.output.S2.m[, 2], estimation.output.S2.m[, 3], xlab="beta", ylab="fake")
plot(estimation.output.S2.m[, 1], estimation.output.S2.m[, 2],
     bty="n", 
     xlab=TeX('$\\hat{\\alpha}$'), ylab=TeX('$\\hat{\\beta}$'),
     col=coul[5])


##### Third Fake News effect from Dorje like fake news #####


fake_effect_3 <- second_fake_effect(1:30, -1, 0.01, 10)

# t=1 to t=30.... simulating from mu=0.5, sd=2
set.seed(9)

# specifying data frames
scen1.s3.z <- matrix(nrow=30, ncol=30)
scen1.s3.m1.mu.post <- matrix(nrow=29, ncol=30)
scen1.s3.m2.mu.post <- matrix(nrow=29, ncol=30)
scen1.s3.m3.mu.post <- matrix(nrow=29, ncol=30)
scen1.s3.m4.mu.post <- matrix(nrow=29, ncol=30)
scen1.s3.m4.A <- matrix(nrow=29, ncol=90)
scen1.s3.m4.beta <- matrix(nrow=29, ncol=90)
scen1.s3.m5.mu.post <- matrix(nrow=29, ncol=30)
scen1.s3.m5.A <- matrix(nrow=29, ncol=90)
scen1.s3.m5.beta <- matrix(nrow=29, ncol=90)
scen1.s3.m5.fake <- matrix(nrow=29, ncol=90)


for(s in 1:30){
  current_sim = rnorm(30, mean=0.5, sd=2) + fake_effect_3
  scen1.s3.z[, s] <- current_sim
  for(t in 2:30){
    # model 1 
    moment <- sigma_known(0, 1000, 4, mean(current_sim[1:t]), t)
    scen1.s3.m1.mu.post[t-1, s] <- moment$posterior_prob
    # model 2
    scen1.s3.m2.mu.post[t-1, s] <- sigma_unknown(current_sim[1:t])$posterior_prob
    # model 3
    robust.model = jags.model(file = "bias_model.txt",
                              data = list(y = current_sim[1:t],
                                          N = t), 
                              n.chains = 3,
                              inits = init_robust,
                              n.adapt = 10000)  
    robust.output = coda.samples(model=robust.model, 
                                 variable.names=c("mu"),
                                 n.iter=10000)
    robust.output.m <- as.matrix(robust.output)
    scen1.s3.m3.mu.post[t-1, s] = sum(robust.output.m[, 1] > 0) / 30000
    estimation.model.tau.k <- jags.model(file = "normal_tau_known.txt",
                                         data = list(y=current_sim[1:t],
                                                     T=t,
                                                     fake=10),
                                         n.chains = 3,
                                         inits = init_robust,
                                         n.adapt = 2500)
    estimation.output.tau.k <- coda.samples(model = estimation.model.tau.k,
                                            variable.names=c("mu",
                                                             "A",
                                                             "beta"),
                                            n.iter=10000)
    column_index = 3*(s-1) + c(1, 2, 3)
    estimation.output.tau.m <- as.matrix(estimation.output.tau.k)
    scen1.s3.m4.mu.post[(t-1), s] = sum(estimation.output.tau.m[, 3] > 0) / 30000
    scen1.s3.m4.A[(t-1), column_index] <- summary(estimation.output.tau.k)$quantiles[1, c(1, 3, 5)]
    scen1.s3.m4.beta[(t-1), column_index] <- summary(estimation.output.tau.k)$quantiles[2, c(1, 3, 5)]
    estimation.model = jags.model(file = "normal_estimation_model.txt",
                                  data = list(y=current_sim[1:t], 
                                              T=t), 
                                  n.chains = 3,
                                  inits = init_robust,
                                  n.adapt = 2500)  
    estimation.output = coda.samples(model=estimation.model, 
                                     variable.names=c("mu",
                                                      "A",
                                                      "beta", 
                                                      "fake"),
                                     n.iter=10000)
    estimation.output.m <- as.matrix(estimation.output)
    scen1.s3.m5.mu.post[(t-1), s] = sum(estimation.output.m[, 4] > 0) / 30000
    scen1.s3.m5.A[(t-1), column_index] <- summary(estimation.output)$quantiles[1, c(1, 3, 5)]
    scen1.s3.m5.beta[(t-1), column_index] <- summary(estimation.output)$quantiles[2, c(1, 3, 5)]
    scen1.s3.m5.fake[(t-1), column_index] <- summary(estimation.output)$quantiles[3, c(1, 3, 5)]
  }
}

# getting alpha versus beta comparsion plot for parameter set 3
set.seed(9)
current_sim = rnorm(30, 0.5, 2) + fake_effect_3 
#current_sim = rnorm(30, 0.5, 2) + fake_effect_1
estimation.model = jags.model(file = "normal_estimation_model.txt",
                              data = list(y=current_sim[1:15], 
                                          T=15), 
                              n.chains = 3,
                              inits = init_robust,
                              n.adapt = 2500)  
estimation.output = coda.samples(model=estimation.model, 
                                 variable.names=c("mu",
                                                  "A",
                                                  "beta", 
                                                  "fake"),
                                 n.iter=10000)
estimation.output.S3.m <- as.matrix(estimation.output)
plot(estimation.output.S2.m[, 1], estimation.output.S2.m[, 2], xlab="A", ylab="beta")
plot(estimation.output.S2.m[, 1], estimation.output.S2.m[, 3], xlab="A", ylab="fake")
plot(estimation.output.S2.m[, 2], estimation.output.S2.m[, 3], xlab="beta", ylab="fake")
plot(estimation.output.S3.m[, 1], estimation.output.S3.m[, 2],
     bty="n", 
     xlab=TeX('$\\hat{\\alpha}$'), ylab=TeX('$\\hat{\\beta}$'),
     col=coul[6])

quantile(estimation.output.S3.m[, "A"],
         probs = c(0.025, 0.5, 0.975))










fade <- add.alpha(brewer.pal(8, "Dark2"), 0.1)

apply(scen1.s1.z, 2, cumsum)
## this shows we want columns 1, 2, 8 highlighted
plot(1, type="n", xlim=c(0, 30), ylim=c(-40, 30),
     bty="n", xlab="Time", ylab=expression(xi[t]),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  lines(seq(0, 30), c(0, cumsum(scen1.s1.z[, i])), "l", col=fade[8])
}
lines(seq(0, 30), c(0, cumsum(scen1.s1.z[, 1])), col=coul[1])
lines(seq(0, 30), c(0, cumsum(scen1.s1.z[, 2])), col=coul[2])
lines(seq(0, 30), c(0, cumsum(scen1.s1.z[, 8])), col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu > 0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(i in 1:30){
  lines(seq(2, 30), scen1.s1.m1.mu.post[, i], "l", col=fade[8])
}
lines(seq(2, 30), scen1.s1.m1.mu.post[, 1], col=coul[1])
lines(seq(2, 30), scen1.s1.m1.mu.post[, 2], col=coul[2])
lines(seq(2, 30), scen1.s1.m1.mu.post[, 8], col=coul[3])



plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu > 0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(i in 1:30){
  lines(seq(2, 30), scen1.s1.m2.mu.post[, i], "l", col=fade[8])
}
lines(seq(2, 30), scen1.s1.m2.mu.post[, 1], col=coul[1])
lines(seq(2, 30), scen1.s1.m2.mu.post[, 2], col=coul[2])
lines(seq(2, 30), scen1.s1.m2.mu.post[, 8], col=coul[3])


plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu > 0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(i in 1:30){
  lines(seq(2, 30), scen1.s1.m3.mu.post[, i], "l", col=fade[8])
}
lines(seq(2, 30), scen1.s1.m3.mu.post[, 1], col=coul[1])
lines(seq(2, 30), scen1.s1.m3.mu.post[, 2], col=coul[2])
lines(seq(2, 30), scen1.s1.m3.mu.post[, 8], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu > 0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(i in 1:30){
  lines(seq(2, 30), scen1.s1.m4.mu.post[, i], "l", col=fade[8])
}
lines(seq(2, 30), scen1.s1.m4.mu.post[, 1], col=coul[1])
lines(seq(2, 30), scen1.s1.m4.mu.post[, 2], col=coul[2])
lines(seq(2, 30), scen1.s1.m4.mu.post[, 8], col=coul[3])


plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu > 0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(i in 1:30){
  lines(seq(2, 30), scen1.s1.m5.mu.post[, i], "l", col=fade[8])
}
lines(seq(2, 30), scen1.s1.m5.mu.post[, 1], col=coul[1])
lines(seq(2, 30), scen1.s1.m5.mu.post[, 2], col=coul[2])
lines(seq(2, 30), scen1.s1.m5.mu.post[, 8], col=coul[3])

### plotting A for m4 and m5 simulation 1
plot(1, type="n", xlim=c(10, 30), ylim=c(-20, 15),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\alpha}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  #lines(seq(10, 30), scen1.s1.m4.A[9:29, ind[1]], col="red")
  lines(seq(10, 30), scen1.s1.m4.A[9:29, ind[2]], col=fade[8])
}
lines(seq(10, 30), scen1.s1.m4.A[9:29, 2], col=coul[1])
lines(seq(10, 30), scen1.s1.m4.A[9:29, 5], col=coul[2])
lines(seq(10, 30), scen1.s1.m4.A[9:29, 23], col=coul[3])

plot(1, type="n", xlim=c(2, 30), ylim=c(-20, 15),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\alpha}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  lines(seq(2, 30), scen1.s1.m5.A[1:29, ind[2]], col=fade[8])
}
lines(seq(2, 30), scen1.s1.m5.A[1:29, 2], col=coul[1])
lines(seq(2, 30), scen1.s1.m5.A[1:29, 5], col=coul[2])
lines(seq(2, 30), scen1.s1.m5.A[1:29, 23], col=coul[3])

# plotting beta for m4 and m5 

plot(1, type="n", xlim=c(10, 30), ylim=c(0, 150),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\beta}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  #lines(seq(10, 30), scen1.s1.m4.A[9:29, ind[1]], col="red")
  lines(seq(10, 30), scen1.s1.m4.beta[9:29, ind[2]], col=fade[8])
}
lines(seq(10, 30), scen1.s1.m4.beta[9:29, 2], col=coul[1])
lines(seq(10, 30), scen1.s1.m4.beta[9:29, 5], col=coul[2])
lines(seq(10, 30), scen1.s1.m4.beta[9:29, 23], col=coul[3])

plot(1, type="n", xlim=c(10, 30), ylim=c(0, 5),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\beta}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  #lines(seq(10, 30), scen1.s1.m4.A[9:29, ind[1]], col="red")
  lines(seq(10, 30), scen1.s1.m4.beta[9:29, ind[2]], col=fade[8])
}
lines(seq(10, 30), scen1.s1.m4.beta[9:29, 2], col=coul[1])
lines(seq(10, 30), scen1.s1.m4.beta[9:29, 5], col=coul[2])
lines(seq(10, 30), scen1.s1.m4.beta[9:29, 23], col=coul[3])


plot(1, type="n", xlim=c(2, 30), ylim=c(0, 5),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\beta}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  lines(seq(2, 30), scen1.s1.m5.beta[1:29, ind[2]], col=fade[8])
}
lines(seq(2, 30), scen1.s1.m5.beta[1:29, 2], col=coul[1])
lines(seq(2, 30), scen1.s1.m5.beta[1:29, 5], col=coul[2])
lines(seq(2, 30), scen1.s1.m5.beta[1:29, 23], col=coul[3])

# m5 predicting the fake news
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 30),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\tau}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  lines(seq(2, 30), scen1.s1.m5.fake[1:29, ind[2]], col=fade[8])
}
lines(seq(2, 30), scen1.s1.m5.fake[1:29, 2], col=coul[1])
lines(seq(2, 30), scen1.s1.m5.fake[1:29, 5], col=coul[2])
lines(seq(2, 30), scen1.s1.m5.fake[1:29, 23], col=coul[3])














####### doing the plots for the second simulation
apply(scen1.s2.z, 2, cumsum)[30, ]
## this shows we want columns 1, 2, 8 highlighted
plot(1, type="n", xlim=c(0, 30), ylim=c(-40, 30),
     bty="n", xlab="Time", ylab=expression(xi[t]),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  lines(seq(0, 30), c(0, cumsum(scen1.s2.z[, i])), "l", col=fade[8])
}
lines(seq(0, 30), c(0, cumsum(scen1.s1.z[, 20])), col=coul[1])
lines(seq(0, 30), c(0, cumsum(scen1.s1.z[, 21])), col=coul[2])
lines(seq(0, 30), c(0, cumsum(scen1.s1.z[, 22])), col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu > 0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(i in 1:30){
  lines(seq(2, 30), scen1.s2.m1.mu.post[, i], "l", col=fade[8])
}
lines(seq(2, 30), scen1.s2.m1.mu.post[, 20], col=coul[1])
lines(seq(2, 30), scen1.s2.m1.mu.post[, 21], col=coul[2])
lines(seq(2, 30), scen1.s2.m1.mu.post[, 22], col=coul[3])



plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu > 0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(i in 1:30){
  lines(seq(2, 30), scen1.s2.m2.mu.post[, i], "l", col=fade[8])
}
lines(seq(2, 30), scen1.s2.m2.mu.post[, 20], col=coul[1])
lines(seq(2, 30), scen1.s2.m2.mu.post[, 21], col=coul[2])
lines(seq(2, 30), scen1.s2.m2.mu.post[, 22], col=coul[3])


plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu > 0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(i in 1:30){
  lines(seq(2, 30), scen1.s2.m3.mu.post[, i], "l", col=fade[8])
}
lines(seq(2, 30), scen1.s2.m3.mu.post[, 20], col=coul[1])
lines(seq(2, 30), scen1.s2.m3.mu.post[, 21], col=coul[2])
lines(seq(2, 30), scen1.s2.m3.mu.post[, 22], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu > 0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(i in 1:30){
  lines(seq(2, 30), scen1.s2.m4.mu.post[, i], "l", col=fade[8])
}
lines(seq(2, 30), scen1.s2.m4.mu.post[, 20], col=coul[1])
lines(seq(2, 30), scen1.s2.m4.mu.post[, 21], col=coul[2])
lines(seq(2, 30), scen1.s2.m4.mu.post[, 22], col=coul[3])


plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu > 0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(i in 1:30){
  lines(seq(2, 30), scen1.s2.m5.mu.post[, i], "l", col=fade[8])
}
lines(seq(2, 30), scen1.s2.m5.mu.post[, 20], col=coul[1])
lines(seq(2, 30), scen1.s2.m5.mu.post[, 21], col=coul[2])
lines(seq(2, 30), scen1.s2.m5.mu.post[, 22], col=coul[3])

### plotting A for m4 and m5 simulation 2
plot(1, type="n", xlim=c(10, 30), ylim=c(-50, 0),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\alpha}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  #lines(seq(10, 30), scen1.s1.m4.A[9:29, ind[1]], col="red")
  lines(seq(10, 30), scen1.s2.m4.A[9:29, ind[2]], col=fade[8])
}
lines(seq(10, 30), scen1.s2.m4.A[9:29, 59], col=coul[1])
lines(seq(10, 30), scen1.s2.m4.A[9:29, 62], col=coul[2])
lines(seq(10, 30), scen1.s2.m4.A[9:29, 65], col=coul[3])

plot(1, type="n", xlim=c(2, 30), ylim=c(-50, 0),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\alpha}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  lines(seq(2, 30), scen1.s2.m5.A[1:29, ind[2]], col=fade[8])
}
lines(seq(2, 30), scen1.s2.m5.A[1:29, 59], col=coul[1])
lines(seq(2, 30), scen1.s2.m5.A[1:29, 62], col=coul[2])
lines(seq(2, 30), scen1.s2.m5.A[1:29, 65], col=coul[3])

# plotting beta for m4 and m5 

plot(1, type="n", xlim=c(10, 30), ylim=c(0, 75),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\beta}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  #lines(seq(10, 30), scen1.s1.m4.A[9:29, ind[1]], col="red")
  lines(seq(10, 30), scen1.s2.m4.beta[9:29, ind[2]], col=fade[8])
}
lines(seq(10, 30), scen1.s2.m4.beta[9:29, 59], col=coul[1])
lines(seq(10, 30), scen1.s2.m4.beta[9:29, 62], col=coul[2])
lines(seq(10, 30), scen1.s2.m4.beta[9:29, 65], col=coul[3])

plot(1, type="n", xlim=c(10, 30), ylim=c(0, 5),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\beta}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  #lines(seq(10, 30), scen1.s1.m4.A[9:29, ind[1]], col="red")
  lines(seq(10, 30), scen1.s2.m4.beta[9:29, ind[2]], col=fade[8])
}
lines(seq(10, 30), scen1.s2.m4.beta[9:29, 59], col=coul[1])
lines(seq(10, 30), scen1.s2.m4.beta[9:29, 62], col=coul[2])
lines(seq(10, 30), scen1.s2.m4.beta[9:29, 65], col=coul[3])


plot(1, type="n", xlim=c(2, 30), ylim=c(0, 5),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\beta}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  lines(seq(2, 30), scen1.s2.m5.beta[1:29, ind[2]], col=fade[8])
}
lines(seq(2, 30), scen1.s2.m5.beta[1:29, 59], col=coul[1])
lines(seq(2, 30), scen1.s2.m5.beta[1:29, 62], col=coul[2])
lines(seq(2, 30), scen1.s2.m5.beta[1:29, 65], col=coul[3])

# m5 predicting the fake news
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 30),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\tau}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  lines(seq(2, 30), scen1.s2.m5.fake[1:29, ind[2]], col=fade[8])
}
lines(seq(2, 30), scen1.s2.m5.fake[1:29, 59], col=coul[1])
lines(seq(2, 30), scen1.s2.m5.fake[1:29, 62], col=coul[2])
lines(seq(2, 30), scen1.s2.m5.fake[1:29, 65], col=coul[3])














####### doing the plots for the third simulation
apply(scen1.s3.z, 2, cumsum)[30, ]
## this shows we want columns 1, 2, 8 highlighted
plot(1, type="n", xlim=c(0, 30), ylim=c(-30, 30),
     bty="n", xlab="Time", ylab=expression(xi[t]),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  lines(seq(0, 30), c(0, cumsum(scen1.s3.z[, i])), "l", col=fade[8])
}
lines(seq(0, 30), c(0, cumsum(scen1.s3.z[, 1])), col=coul[1])
lines(seq(0, 30), c(0, cumsum(scen1.s3.z[, 29])), col=coul[2])
lines(seq(0, 30), c(0, cumsum(scen1.s3.z[, 30])), col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu > 0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(i in 1:30){
  lines(seq(2, 30), scen1.s3.m1.mu.post[, i], "l", col=fade[8])
}
lines(seq(2, 30), scen1.s3.m1.mu.post[, 1], col=coul[1])
lines(seq(2, 30), scen1.s3.m1.mu.post[, 29], col=coul[2])
lines(seq(2, 30), scen1.s3.m1.mu.post[, 30], col=coul[3])



plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu > 0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(i in 1:30){
  lines(seq(2, 30), scen1.s3.m2.mu.post[, i], "l", col=fade[8])
}
lines(seq(2, 30), scen1.s3.m2.mu.post[, 1], col=coul[1])
lines(seq(2, 30), scen1.s3.m2.mu.post[, 29], col=coul[2])
lines(seq(2, 30), scen1.s3.m2.mu.post[, 30], col=coul[3])


plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu > 0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(i in 1:30){
  lines(seq(2, 30), scen1.s3.m3.mu.post[, i], "l", col=fade[8])
}
lines(seq(2, 30), scen1.s3.m3.mu.post[, 1], col=coul[1])
lines(seq(2, 30), scen1.s3.m3.mu.post[, 29], col=coul[2])
lines(seq(2, 30), scen1.s3.m3.mu.post[, 30], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu > 0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(i in 1:30){
  lines(seq(2, 30), scen1.s3.m4.mu.post[, i], "l", col=fade[8])
}
lines(seq(2, 30), scen1.s3.m4.mu.post[, 1], col=coul[1])
lines(seq(2, 30), scen1.s3.m4.mu.post[, 29], col=coul[2])
lines(seq(2, 30), scen1.s3.m4.mu.post[, 30], col=coul[3])


plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu > 0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(i in 1:30){
  lines(seq(2, 30), scen1.s3.m5.mu.post[, i], "l", col=fade[8])
}
lines(seq(2, 30), scen1.s3.m5.mu.post[, 1], col=coul[1])
lines(seq(2, 30), scen1.s3.m5.mu.post[, 29], col=coul[2])
lines(seq(2, 30), scen1.s3.m5.mu.post[, 30], col=coul[3])

### plotting A for m4 and m5 simulation 2
plot(1, type="n", xlim=c(10, 30), ylim=c(-10, 10),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\alpha}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  #lines(seq(10, 30), scen1.s1.m4.A[9:29, ind[1]], col="red")
  lines(seq(10, 30), scen1.s3.m4.A[9:29, ind[2]], col=fade[8])
}
lines(seq(10, 30), scen1.s3.m4.A[9:29, 2], col=coul[1])
lines(seq(10, 30), scen1.s3.m4.A[9:29, 86], col=coul[2])
lines(seq(10, 30), scen1.s3.m4.A[9:29, 89], col=coul[3])

plot(1, type="n", xlim=c(2, 30), ylim=c(-10, 10),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\alpha}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  lines(seq(2, 30), scen1.s3.m5.A[1:29, ind[2]], col=fade[8])
}
lines(seq(2, 30), scen1.s3.m5.A[1:29, 2], col=coul[1])
lines(seq(2, 30), scen1.s3.m5.A[1:29, 86], col=coul[2])
lines(seq(2, 30), scen1.s3.m5.A[1:29, 89], col=coul[3])

# plotting beta for m4 and m5 

plot(1, type="n", xlim=c(10, 30), ylim=c(0, 60),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\beta}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  #lines(seq(10, 30), scen1.s1.m4.A[9:29, ind[1]], col="red")
  lines(seq(10, 30), scen1.s3.m4.beta[9:29, ind[2]], col=fade[8])
}
lines(seq(10, 30), scen1.s3.m4.beta[9:29, 2], col=coul[1])
lines(seq(10, 30), scen1.s3.m4.beta[9:29, 86], col=coul[2])
lines(seq(10, 30), scen1.s3.m4.beta[9:29, 89], col=coul[3])

plot(1, type="n", xlim=c(10, 30), ylim=c(0, 0.10),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\beta}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  #lines(seq(10, 30), scen1.s1.m4.A[9:29, ind[1]], col="red")
  lines(seq(10, 30), scen1.s3.m4.beta[9:29, ind[2]], col=fade[8])
}
lines(seq(10, 30), scen1.s3.m4.beta[9:29, 2], col=coul[1])
lines(seq(10, 30), scen1.s3.m4.beta[9:29, 86], col=coul[2])
lines(seq(10, 30), scen1.s3.m4.beta[9:29, 89], col=coul[3])


plot(1, type="n", xlim=c(2, 30), ylim=c(0, 0.1),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\beta}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  lines(seq(2, 30), scen1.s3.m5.beta[1:29, ind[2]], col=fade[8])
}
lines(seq(2, 30), scen1.s3.m5.beta[1:29, 59], col=coul[1])
lines(seq(2, 30), scen1.s3.m5.beta[1:29, 62], col=coul[2])
lines(seq(2, 30), scen1.s3.m5.beta[1:29, 65], col=coul[3])

# m5 predicting the fake news
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 30),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\tau}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  lines(seq(2, 30), scen1.s3.m5.fake[1:29, ind[2]], col=fade[8])
}
lines(seq(2, 30), scen1.s3.m5.fake[1:29, 2], col=coul[1])
lines(seq(2, 30), scen1.s3.m5.fake[1:29, 86], col=coul[2])
lines(seq(2, 30), scen1.s3.m5.fake[1:29, 89], col=coul[3])


















































plot(1, type="n", xlab="", ylab="", xlim=c(2, 30), ylim=c(0, 30))
for(i in c(2, 5, 8, 11, 14)){
  lines(seq(2, 30), scen1.s1.m5.fake[, i], "l")
}
plot(1, type="n", xlab="", ylab="", xlim=c(2, 30), ylim=c(-15, 0))
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  lines(seq(2, 30), scen1.s1.m4.A[, ind[2]])
}
plot(1, type="n", xlab="", ylab="", xlim=c(2, 30), ylim=c(0, 30))
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  lines(seq(2, 30), scen1.s3.m5.fake[, ind[2]])
}


plot(1, type="n", xlab="", ylab="", xlim=c(2, 30), ylim=c(0, 30))
for(i in 1:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  lines(seq(2, 30), scen1.s3.m5.beta[, ind[2]])
}


######## Simulation 2
plot(1, type="n", xlab="", ylab="", xlim=c(0, 30), ylim=c(-40, 30))
for(i in 1:30){
  lines(seq(0, 30), c(0, cumsum(scen1.s2.z[, i])), "l")
}


plot(1, type="n", xlab="", ylab="", xlim=c(2, 30), ylim=c(0, 1))
for(i in 1:30){
  lines(seq(2, 30), scen1.s2.m1.mu.post[, i], "l")
}
plot(1, type="n", xlab="", ylab="", xlim=c(2, 30), ylim=c(0, 1))
for(i in 1:30){
  lines(seq(2, 30), scen1.s2.m2.mu.post[, i], "l")
}
plot(1, type="n", xlab="", ylab="", xlim=c(2, 30), ylim=c(0, 1))
for(i in 1:30){
  lines(seq(2, 30), scen1.s2.m3.mu.post[, i], "l")
}
plot(1, type="n", xlab="", ylab="", xlim=c(2, 30), ylim=c(0, 1))
for(i in 1:30){
  lines(seq(2, 30), scen1.s2.m4.mu.post[, i], "l")
}
plot(1, type="n", xlab="", ylab="", xlim=c(2, 30), ylim=c(0, 1))
for(i in 1:30){
  lines(seq(2, 30), scen1.s2.m5.mu.post[, i], "l")
}

### plotting A for m4 and m5 simulation 1
plot(1, type="n", xlab="", ylab="", xlim=c(10, 30), ylim=c(-40, 0))
for(i in 10:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  #lines(seq(10, 30), scen1.s1.m4.A[9:29, ind[1]], col="red")
  lines(seq(10, 30), scen1.s2.m4.A[9:29, ind[2]])
  #lines(seq(10, 30), scen1.s1.m4.A[9:29, ind[3]], col="blue")
}
plot(1, type="n", xlab="", ylab="", xlim=c(2, 30), ylim=c(-60, 40))
for(i in 10:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  #lines(seq(2, 30), scen1.s1.m5.A[1:29, ind[1]], col="red")
  lines(seq(2, 30), scen1.s2.m5.A[1:29, ind[2]])
  #lines(seq(2, 30), scen1.s1.m5.A[1:29, ind[3]], col="blue")
  abline(a=-10, b =0, col="red")
}

## plotting beta for m4 and m5
plot(1, type="n", xlab="", ylab="", xlim=c(10, 30), ylim=c(0, 15))
for(i in 10:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  #lines(seq(10, 30), scen1.s1.m4.beta[9:29, ind[1]], col="red")
  lines(seq(10, 30), scen1.s2.m4.beta[9:29, ind[2]])
  #lines(seq(10, 30), scen1.s1.m4.beta[9:29, ind[3]], col="blue")
}

plot(1, type="n", xlab="", ylab="", xlim=c(2, 30), ylim=c(0, 5))
for(i in 10:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  #lines(seq(2, 30), scen1.s1.m5.A[1:29, ind[1]], col="red")
  lines(seq(2, 30), scen1.s2.m5.beta[1:29, ind[2]])
  #lines(seq(2, 30), scen1.s1.m5.A[1:29, ind[3]], col="blue")
  abline(a=-10, b =0, col="red")
}

# m5 predicting the fake news
plot(1, type="n", xlab="", ylab="", xlim=c(2, 30), ylim=c(0, 30))
for(i in 10:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  #lines(seq(2, 30), scen1.s1.m5.A[1:29, ind[1]], col="red")
  lines(seq(2, 30), scen1.s2.m5.fake[1:29, ind[2]])
  #lines(seq(2, 30), scen1.s1.m5.A[1:29, ind[3]], col="blue")
  abline(a=-10, b =0, col="red")
}

#### Simulation 3
######## Simulation 2
plot(1, type="n", xlab="", ylab="", xlim=c(0, 30), ylim=c(-40, 30))
for(i in 1:30){
  lines(seq(0, 30), c(0, cumsum(scen1.s3.z[, i])), "l")
}


plot(1, type="n", xlab="", ylab="", xlim=c(2, 30), ylim=c(0, 1))
for(i in 1:30){
  lines(seq(2, 30), scen1.s3.m1.mu.post[, i], "l")
}
plot(1, type="n", xlab="", ylab="", xlim=c(2, 30), ylim=c(0, 1))
for(i in 1:30){
  lines(seq(2, 30), scen1.S3.m2.mu.post[, i], "l")
}
plot(1, type="n", xlab="", ylab="", xlim=c(2, 30), ylim=c(0, 1))
for(i in 1:30){
  lines(seq(2, 30), scen1.s3.m3.mu.post[, i], "l")
}
plot(1, type="n", xlab="", ylab="", xlim=c(2, 30), ylim=c(0, 1))
for(i in 1:30){
  lines(seq(2, 30), scen1.s3.m4.mu.post[, i], "l")
}
plot(1, type="n", xlab="", ylab="", xlim=c(2, 30), ylim=c(0, 1))
for(i in 1:30){
  lines(seq(2, 30), scen1.s3.m5.mu.post[, i], "l")
}

### plotting A for m4 and m5 simulation 1
plot(1, type="n", xlab="", ylab="", xlim=c(10, 30), ylim=c(-40, 0))
for(i in 10:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  #lines(seq(10, 30), scen1.s1.m4.A[9:29, ind[1]], col="red")
  lines(seq(10, 30), scen1.s3.m4.A[9:29, ind[2]])
  #lines(seq(10, 30), scen1.s1.m4.A[9:29, ind[3]], col="blue")
}
plot(1, type="n", xlab="", ylab="", xlim=c(2, 30), ylim=c(-60, 40))
for(i in 10:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  #lines(seq(2, 30), scen1.s1.m5.A[1:29, ind[1]], col="red")
  lines(seq(2, 30), scen1.s3.m5.A[1:29, ind[2]])
  #lines(seq(2, 30), scen1.s1.m5.A[1:29, ind[3]], col="blue")
  abline(a=-10, b =0, col="red")
}

## plotting beta for m4 and m5
plot(1, type="n", xlab="", ylab="", xlim=c(10, 30), ylim=c(0, 15))
for(i in 10:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  #lines(seq(10, 30), scen1.s1.m4.beta[9:29, ind[1]], col="red")
  lines(seq(10, 30), scen1.s2.m4.beta[9:29, ind[2]])
  #lines(seq(10, 30), scen1.s1.m4.beta[9:29, ind[3]], col="blue")
}

plot(1, type="n", xlab="", ylab="", xlim=c(2, 30), ylim=c(0, 5))
for(i in 10:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  #lines(seq(2, 30), scen1.s1.m5.A[1:29, ind[1]], col="red")
  lines(seq(2, 30), scen1.s3.m5.beta[1:29, ind[2]])
  #lines(seq(2, 30), scen1.s1.m5.A[1:29, ind[3]], col="blue")
  abline(a=-10, b =0, col="red")
}

# m5 predicting the fake news
plot(1, type="n", xlab="", ylab="", xlim=c(2, 30), ylim=c(0, 30))
for(i in 10:30){
  ind <- (i-1)*3 + c(1, 2, 3)
  #lines(seq(2, 30), scen1.s1.m5.A[1:29, ind[1]], col="red")
  lines(seq(2, 30), scen1.s3.m5.fake[1:29, ind[2]])
  #lines(seq(2, 30), scen1.s1.m5.A[1:29, ind[3]], col="blue")
  abline(a=-10, b =0, col="red")
}


set.seed(7)
scen1.s1.m2.mu.post <- 
  
  
  
  
  
  
  
normal_estimation_model <- "  
model{
  for(t in 1:T){
    bias[t] <- ifelse(t < fake, 0, A / (exp(beta * (t - fake))))
    y[t] ~ dnorm(mu + bias[t], tau)
  }
  mu ~ dnorm(0, 0.001)
  A ~ dunif(-100, 100)
  beta ~ dgamma(0.001, 0.001)
  tau ~ dexp(0.001)
  sigma <- 1 / sqrt(tau)
  fake ~ dunif(0, 30)
}
"
cat(normal_estimation_model, file="normal_estimation_model.txt")
  
normal_tau_known <- "    
model{
  for(t in 1:T){
    bias[t] <- ifelse(t < fake, 0, A / (exp(beta * (t - fake))))
    y[t] ~ dnorm(mu + bias[t], tau)
  }
  mu ~ dnorm(0, 0.001)
  A ~ dunif(-100, 100)
  beta ~ dgamma(0.001, 0.001)
  tau ~ dexp(0.001)
  sigma <- 1 / sqrt(tau)
}
"
cat(normal_tau_known, file="normal_tau_known.txt")  
