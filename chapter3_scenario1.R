###### Chapter 3 Autoregressive fake news simulation


set.seed(2020)
###### ARsim1 <- numeric(30)

# first simulation parameter set
prior_time <- rep(1/30, 30)
ARmu1 <- 0.5
ARsig1 <- 2
ARphi1 <- 0.2
ARfake1 <- c(rep(0, 8), -15, rep(0, 20))

# specifying AR simulations for the first parameter set
ARsim1 <- matrix(nrow=30, ncol=30)
for(s in 1:30){
  ARsim1[1, s] <- rnorm(1, ARmu1, ARsig1)
  noise <- rnorm(29, 0, sqrt(ARsig1*(1-ARphi1^2)))
  for(t in 2:30){
    ARsim1[t, s] <- (1-ARphi1)*ARmu1 + ARphi1*ARsim1[t-1, s] + noise[t-1] + ARfake1[t-1]
  }
}

AR.S1.M1.post <- matrix(nrow=29, ncol=30)
AR.S1.M1.mu <- matrix(nrow=29, ncol=90)
AR.S1.M2.post <- matrix(nrow=29, ncol=30)
AR.S1.M2.mu <- matrix(nrow=29, ncol=90)
AR.S1.M2.sigma <- matrix(nrow=29, ncol=90)
AR.S1.M3.post <- matrix(nrow=29, ncol=30)
AR.S1.M3.mu <- matrix(nrow=29, ncol=90)
AR.S1.M3.sigma <- matrix(nrow=29, ncol=90)
AR.S1.M4.post <- matrix(nrow=29, ncol=30) 
AR.S1.M4.mu <- matrix(nrow=29, ncol=90)
AR.S1.M4.sigma <- matrix(nrow=29, ncol=90)
AR.S1.M4.A <- matrix(nrow=29, ncol=90)
AR.S1.M5.post <- matrix(nrow=29, ncol=30) 
AR.S1.M5.mu <- matrix(nrow=29, ncol=90)
AR.S1.M5.sigma <- matrix(nrow=29, ncol=90)
AR.S1.M5.A <- matrix(nrow=29, ncol=90)
AR.S1.M5.tau <- matrix(nrow=29, ncol=90)

for(s in 1:30){
  column_index <- (s-1)*3 + c(1, 2, 3)
  for(t in 2:30){
  # model 1
      M1.model = jags.model(file = "ar_known_sigma.txt",
                            data = list(y=ARsim1[1:t, s], 
                                        T=t,
                                        tau=1/3), 
                            inits = init_robust,  
                            n.chains = 3,
                            n.adapt = 2500)  
      
      # creating output object
      M1.output = coda.samples(model=M1.model, 
                               variable.names=c("mu"),
                               n.iter=10000)
      M1.output.m <- as.matrix(M1.output)
      AR.S1.M1.post[(t-1), s] <- sum(M1.output.m[, 1] > 0) / 30000
      AR.S1.M1.mu[(t-1), column_index] <- summary(M1.output)$quantiles[c(1, 3, 5)]
      
      
      # model 2
      M2.model = jags.model(file = "ar_model.txt",
                            data = list(y=ARsim1[1:t, s], 
                                        T=t), 
                            inits = init_robust,
                            n.chains = 3,
                            n.adapt = 2500)  
      M2.output = coda.samples(model=M2.model, 
                               variable.names=c("mu",
                                                "sigma"),
                               n.iter=10000)
      M2.output.m <- as.matrix(M2.output)
      AR.S1.M2.post[(t-1), s] <- sum(M2.output.m[, 1] > 0) / 30000
      AR.S1.M2.mu[(t-1), column_index] <- summary(M2.output)$quantiles[1, c(1, 3, 5)]
      AR.S1.M2.sigma[(t-1), column_index] <- summary(M2.output)$quantiles[2, c(1, 3, 5)]
      
      
      # model 3
      M3.model = jags.model(file = "var_ar_model.txt",
                            data = list(y=ARsim1[1:t, s], 
                                        T=t), 
                            inits = init_robust,
                            n.chains = 3,
                            n.adapt = 2500)  
      M3.output = coda.samples(model=M3.model, 
                               variable.names=c("mu",
                                                "sigma"),
                               n.iter=10000)
      M3.output.m <- as.matrix(M3.output)
      AR.S1.M3.post[(t-1), s] = sum(M3.output.m[, 1] > 0) / 30000
      AR.S1.M3.mu[(t-1), column_index] <- summary(M3.output)$quantiles[1, c(1, 3, 5)]
      AR.S1.M3.sigma[(t-1), column_index] <- summary(M3.output)$quantiles[2, c(1, 3, 5)]
      
      # model 4
      M4.model = jags.model(file = "AR_fake_known.txt",
                            data = list(y=ARsim1[1:t, s], 
                                        T=t,
                                        fake=10), 
                            inits = init_robust, 
                            n.chains = 3,
                            n.adapt = 2500)  
      M4.output = coda.samples(model=M4.model, 
                               variable.names=c("mu",
                                                "sigma",
                                                "A"),
                               n.iter=10000)
      M4.output.m <- as.matrix(M4.output)
      AR.S1.M4.post[(t-1), s] = sum(M4.output.m[, 2] > 0) / 30000
      AR.S1.M4.A[(t-1), column_index] <- summary(M4.output)$quantiles[1, c(1, 3, 5)]
      AR.S1.M4.mu[(t-1), column_index] <- summary(M4.output)$quantiles[2, c(1, 3, 5)]
      AR.S1.M4.sigma[(t-1), column_index] <- summary(M4.output)$quantiles[3, c(1, 3, 5)]
      
      # model 5
      M5.model = jags.model(file = "AR_estimation.txt",
                            data = list(y=ARsim1[1:t, s], 
                                        T=t,
                                        prior_time=prior_time), 
                            inits = init_robust,
                            n.chains = 3,
                            n.adapt = 2500)  
      M5.output = coda.samples(model=M5.model, 
                               variable.names=c("mu",
                                                "sigma",
                                                "A",
                                                "fake"),
                               n.iter=10000)
      M5.output.m <- as.matrix(M5.output)
      AR.S1.M5.post[(t-1), s] = sum(M5.output.m[, 3] > 0) / 30000
      AR.S1.M5.A[(t-1), column_index] <- summary(M5.output)$quantiles[1, c(1, 3, 5)]
      AR.S1.M5.tau[(t-1), column_index] <- summary(M5.output)$quantiles[2, c(1, 3, 5)]
      AR.S1.M5.mu[(t-1), column_index] <- summary(M5.output)$quantiles[3, c(1, 3, 5)]
      AR.S1.M5.sigma[(t-1), column_index] <- summary(M5.output)$quantiles[4, c(1, 3, 5)]
  }
}

########### Simulation 2 ###############
# second simulation parameter set
prior_time <- rep(1/30, 30)
ARmu2 <- 0.7
ARsig2 <- 2
ARphi2 <- 0.2
ARfake2 <- c(rep(0, 8), -5, rep(0, 20))

# specifying AR simulations for the second parameter set
ARsim2 <- matrix(nrow=30, ncol=30)
for(s in 1:30){
  ARsim2[1, s] <- rnorm(1, ARmu2, ARsig2)
  noise <- rnorm(29, 0, sqrt(ARsig2*(1-ARphi2^2)))
  for(t in 2:30){
    ARsim2[t, s] <- (1-ARphi2)*ARmu2 + ARphi2*ARsim2[t-1, s] + noise[t-1] + ARfake2[t-1]
  }
}

AR.S2.M1.post <- matrix(nrow=29, ncol=30)
AR.S2.M1.mu <- matrix(nrow=29, ncol=90)
AR.S2.M2.post <- matrix(nrow=29, ncol=30)
AR.S2.M2.mu <- matrix(nrow=29, ncol=90)
AR.S2.M2.sigma <- matrix(nrow=29, ncol=90)
AR.S2.M3.post <- matrix(nrow=29, ncol=30)
AR.S2.M3.mu <- matrix(nrow=29, ncol=90)
AR.S2.M3.sigma <- matrix(nrow=29, ncol=90)
AR.S2.M4.post <- matrix(nrow=29, ncol=30) 
AR.S2.M4.mu <- matrix(nrow=29, ncol=90)
AR.S2.M4.sigma <- matrix(nrow=29, ncol=90)
AR.S2.M4.A <- matrix(nrow=29, ncol=90)
AR.S2.M5.post <- matrix(nrow=29, ncol=30) 
AR.S2.M5.mu <- matrix(nrow=29, ncol=90)
AR.S2.M5.sigma <- matrix(nrow=29, ncol=90)
AR.S2.M5.A <- matrix(nrow=29, ncol=90)
AR.S2.M5.tau <- matrix(nrow=29, ncol=90)

for(s in 1:30){
  column_index <- (s-1)*3 + c(1, 2, 3)
  for(t in 2:30){
    # model 1
    M1.model = jags.model(file = "ar_known_sigma.txt",
                          data = list(y=ARsim2[1:t, s], 
                                      T=t,
                                      tau=1/4), 
                          inits = init_robust,  
                          n.chains = 3,
                          n.adapt = 2500)  
    
    # creating output object
    M1.output = coda.samples(model=M1.model, 
                             variable.names=c("mu"),
                             n.iter=10000)
    M1.output.m <- as.matrix(M1.output)
    AR.S2.M1.post[(t-1), s] <- sum(M1.output.m[, 1] > 0) / 30000
    AR.S2.M1.mu[(t-1), column_index] <- summary(M1.output)$quantiles[c(1, 3, 5)]
    
    
    # model 2
    M2.model = jags.model(file = "ar_model.txt",
                          data = list(y=ARsim2[1:t, s], 
                                      T=t), 
                          inits = init_robust,
                          n.chains = 3,
                          n.adapt = 2500)  
    M2.output = coda.samples(model=M2.model, 
                             variable.names=c("mu",
                                              "sigma"),
                             n.iter=10000)
    M2.output.m <- as.matrix(M2.output)
    AR.S2.M2.post[(t-1), s] <- sum(M2.output.m[, 1] > 0) / 30000
    AR.S2.M2.mu[(t-1), column_index] <- summary(M2.output)$quantiles[1, c(1, 3, 5)]
    AR.S2.M2.sigma[(t-1), column_index] <- summary(M2.output)$quantiles[2, c(1, 3, 5)]
    
    
    # model 3
    M3.model = jags.model(file = "var_ar_model.txt",
                          data = list(y=ARsim2[1:t, s], 
                                      T=t), 
                          inits = init_robust,
                          n.chains = 3,
                          n.adapt = 2500)  
    M3.output = coda.samples(model=M3.model, 
                             variable.names=c("mu",
                                              "sigma"),
                             n.iter=10000)
    M3.output.m <- as.matrix(M3.output)
    AR.S2.M3.post[(t-1), s] = sum(M3.output.m[, 1] > 0) / 30000
    AR.S2.M3.mu[(t-1), column_index] <- summary(M3.output)$quantiles[1, c(1, 3, 5)]
    AR.S2.M3.sigma[(t-1), column_index] <- summary(M3.output)$quantiles[2, c(1, 3, 5)]
    
    # model 4
    M4.model = jags.model(file = "AR_fake_known.txt",
                          data = list(y=ARsim2[1:t, s], 
                                      T=t,
                                      fake=10), 
                          inits = init_robust, 
                          n.chains = 3,
                          n.adapt = 2500)  
    M4.output = coda.samples(model=M4.model, 
                             variable.names=c("mu",
                                              "sigma",
                                              "A"),
                             n.iter=10000)
    M4.output.m <- as.matrix(M4.output)
    AR.S2.M4.post[(t-1), s] = sum(M4.output.m[, 2] > 0) / 30000
    AR.S2.M4.A[(t-1), column_index] <- summary(M4.output)$quantiles[1, c(1, 3, 5)]
    AR.S2.M4.mu[(t-1), column_index] <- summary(M4.output)$quantiles[2, c(1, 3, 5)]
    AR.S2.M4.sigma[(t-1), column_index] <- summary(M4.output)$quantiles[3, c(1, 3, 5)]
    
    # model 5
    M5.model = jags.model(file = "AR_estimation.txt",
                          data = list(y=ARsim2[1:t, s], 
                                      T=t,
                                      prior_time=prior_time), 
                          inits = init_robust,
                          n.chains = 3,
                          n.adapt = 2500)  
    M5.output = coda.samples(model=M5.model, 
                             variable.names=c("mu",
                                              "sigma",
                                              "A",
                                              "fake"),
                             n.iter=10000)
    M5.output.m <- as.matrix(M5.output)
    AR.S2.M5.post[(t-1), s] = sum(M5.output.m[, 3] > 0) / 30000
    AR.S2.M5.A[(t-1), column_index] <- summary(M5.output)$quantiles[1, c(1, 3, 5)]
    AR.S2.M5.tau[(t-1), column_index] <- summary(M5.output)$quantiles[2, c(1, 3, 5)]
    AR.S2.M5.mu[(t-1), column_index] <- summary(M5.output)$quantiles[3, c(1, 3, 5)]
    AR.S2.M5.sigma[(t-1), column_index] <- summary(M5.output)$quantiles[4, c(1, 3, 5)]
  }
}

###### third simulation parameter set #######3
prior_time <- rep(1/30, 30)
ARsim3 <- numeric(30)
ARmu3 <- 0.7
ARsig3 <- 2.5
ARphi3 <- 0.4
ARfake3 <- c(rep(0, 8), -5, rep(0, 20))

# specifying AR simulations for the second parameter set
ARsim3 <- matrix(nrow=30, ncol=30)
for(s in 1:30){
  ARsim3[1, s] <- rnorm(1, ARmu3, ARsig3)
  noise <- rnorm(29, 0, sqrt(ARsig3*(1-ARphi3^2)))
  for(t in 2:30){
    ARsim3[t, s] <- (1-ARphi3)*ARmu3 + ARphi3*ARsim3[t-1, s] + noise[t-1] + ARfake3[t-1]
  }
}

AR.S3.M1.post <- matrix(nrow=29, ncol=30)
AR.S3.M1.mu <- matrix(nrow=29, ncol=90)
AR.S3.M2.post <- matrix(nrow=29, ncol=30)
AR.S3.M2.mu <- matrix(nrow=29, ncol=90)
AR.S3.M2.sigma <- matrix(nrow=29, ncol=90)
AR.S3.M3.post <- matrix(nrow=29, ncol=30)
AR.S3.M3.mu <- matrix(nrow=29, ncol=90)
AR.S3.M3.sigma <- matrix(nrow=29, ncol=90)
AR.S3.M4.post <- matrix(nrow=29, ncol=30) 
AR.S3.M4.mu <- matrix(nrow=29, ncol=90)
AR.S3.M4.sigma <- matrix(nrow=29, ncol=90)
AR.S3.M4.A <- matrix(nrow=29, ncol=90)
AR.S3.M5.post <- matrix(nrow=29, ncol=30) 
AR.S3.M5.mu <- matrix(nrow=29, ncol=90)
AR.S3.M5.sigma <- matrix(nrow=29, ncol=90)
AR.S3.M5.A <- matrix(nrow=29, ncol=90)
AR.S3.M5.tau <- matrix(nrow=29, ncol=90)

for(s in 1:30){
  column_index <- (s-1)*3 + c(1, 2, 3)
  for(t in 2:30){
    # model 1
    M1.model = jags.model(file = "ar_known_sigma.txt",
                          data = list(y=ARsim3[1:t, s], 
                                      T=t,
                                      tau=1), 
                          inits = init_robust,  
                          n.chains = 3,
                          n.adapt = 2500)  
    
    # creating output object
    M1.output = coda.samples(model=M1.model, 
                             variable.names=c("mu"),
                             n.iter=10000)
    M1.output.m <- as.matrix(M1.output)
    AR.S3.M1.post[(t-1), s] <- sum(M1.output.m[, 1] > 0) / 30000
    AR.S3.M1.mu[(t-1), column_index] <- summary(M1.output)$quantiles[c(1, 3, 5)]
    
    
    # model 2
    M2.model = jags.model(file = "ar_model.txt",
                          data = list(y=ARsim3[1:t, s], 
                                      T=t), 
                          inits = init_robust,
                          n.chains = 3,
                          n.adapt = 2500)  
    M2.output = coda.samples(model=M2.model, 
                             variable.names=c("mu",
                                              "sigma"),
                             n.iter=10000)
    M2.output.m <- as.matrix(M2.output)
    AR.S3.M2.post[(t-1), s] <- sum(M2.output.m[, 1] > 0) / 30000
    AR.S3.M2.mu[(t-1), column_index] <- summary(M2.output)$quantiles[1, c(1, 3, 5)]
    AR.S3.M2.sigma[(t-1), column_index] <- summary(M2.output)$quantiles[2, c(1, 3, 5)]
    
    
    # model 3
    M3.model = jags.model(file = "var_ar_model.txt",
                          data = list(y=ARsim3[1:t, s], 
                                      T=t), 
                          inits = init_robust,
                          n.chains = 3,
                          n.adapt = 2500)  
    M3.output = coda.samples(model=M3.model, 
                             variable.names=c("mu",
                                              "sigma"),
                             n.iter=10000)
    M3.output.m <- as.matrix(M3.output)
    AR.S3.M3.post[(t-1), s] = sum(M3.output.m[, 1] > 0) / 30000
    AR.S3.M3.mu[(t-1), column_index] <- summary(M3.output)$quantiles[1, c(1, 3, 5)]
    AR.S3.M3.sigma[(t-1), column_index] <- summary(M3.output)$quantiles[2, c(1, 3, 5)]
    
    # model 4
    M4.model = jags.model(file = "AR_fake_known.txt",
                          data = list(y=ARsim3[1:t, s], 
                                      T=t,
                                      fake=10), 
                          inits = init_robust, 
                          n.chains = 3,
                          n.adapt = 2500)  
    M4.output = coda.samples(model=M4.model, 
                             variable.names=c("mu",
                                              "sigma",
                                              "A"),
                             n.iter=10000)
    M4.output.m <- as.matrix(M4.output)
    AR.S3.M4.post[(t-1), s] = sum(M4.output.m[, 2] > 0) / 30000
    AR.S3.M4.A[(t-1), column_index] <- summary(M4.output)$quantiles[1, c(1, 3, 5)]
    AR.S3.M4.mu[(t-1), column_index] <- summary(M4.output)$quantiles[2, c(1, 3, 5)]
    AR.S3.M4.sigma[(t-1), column_index] <- summary(M4.output)$quantiles[3, c(1, 3, 5)]
    
    # model 5
    M5.model = jags.model(file = "AR_estimation.txt",
                          data = list(y=ARsim3[1:t, s], 
                                      T=t,
                                      prior_time=prior_time), 
                          inits = init_robust,
                          n.chains = 3,
                          n.adapt = 2500)  
    M5.output = coda.samples(model=M5.model, 
                             variable.names=c("mu",
                                              "sigma",
                                              "A",
                                              "fake"),
                             n.iter=10000)
    M5.output.m <- as.matrix(M5.output)
    AR.S3.M5.post[(t-1), s] = sum(M5.output.m[, 3] > 0) / 30000
    AR.S3.M5.A[(t-1), column_index] <- summary(M5.output)$quantiles[1, c(1, 3, 5)]
    AR.S3.M5.tau[(t-1), column_index] <- summary(M5.output)$quantiles[2, c(1, 3, 5)]
    AR.S3.M5.mu[(t-1), column_index] <- summary(M5.output)$quantiles[3, c(1, 3, 5)]
    AR.S3.M5.sigma[(t-1), column_index] <- summary(M5.output)$quantiles[4, c(1, 3, 5)]
  }
}

plot(1, type="n", xlim=c(0, 30), ylim=c(-12, 35),
     bty="n", xlab="Time", ylab=expression(xi[t]),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(0, 30), c(0, cumsum(ARsim3[, s])), col=coul.t[1])
}

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu>0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), AR.S1.M1.post[, s], col=coul.t[1])
}

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu>0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), AR.S1.M4.post[, s], col=coul.t[1])
}

plot(1, type="n", xlim=c(0, 30), ylim=c(-12, 35),
     bty="n", xlab="Time", ylab=expression(xi[t]),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(0, 30), c(0, cumsum(ARsim2[, s])), col=coul.t[1])
}

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu>0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), AR.S1.M1.post[, s], col=coul.t[1])
}

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu>0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), AR.S1.M4.post[, s], col=coul.t[1])
}





# Plotting plots for AR Simulation 1 
plot(1, type="n", xlim=c(0, 30), ylim=c(-40, 30),
     bty="n", xlab="Time", ylab=expression(xi[t]),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(0, 30), c(0, cumsum(ARsim3[, s])), col=fade[8])
}
lines(seq(0, 30), c(0, cumsum(ARsim1[, 3])), col=coul[1])
lines(seq(0, 30), c(0, cumsum(ARsim1[, 5])), col=coul[2])
lines(seq(0, 30), c(0, cumsum(ARsim1[, 10])), col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu>0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), AR.S1.M1.post[, s], col=fade[8])
}
lines(seq(2, 30), AR.S1.M1.post[, 3], col=coul[1])
lines(seq(2, 30), AR.S1.M1.post[, 5], col=coul[2])
lines(seq(2, 30),  AR.S1.M1.post[, 10], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu>0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), AR.S1.M2.post[, s], col=fade[8])
}
lines(seq(2, 30), AR.S1.M2.post[, 3], col=coul[1])
lines(seq(2, 30), AR.S1.M2.post[, 5], col=coul[2])
lines(seq(2, 30),  AR.S1.M2.post[, 10], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu>0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), AR.S1.M3.post[, s], col=fade[8])
}
lines(seq(2, 30), AR.S1.M3.post[, 3], col=coul[1])
lines(seq(2, 30), AR.S1.M3.post[, 5], col=coul[2])
lines(seq(2, 30), AR.S1.M3.post[, 10], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu>0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), AR.S1.M4.post[, s], col=fade[8])
}
lines(seq(2, 30), AR.S1.M4.post[, 3], col=coul[1])
lines(seq(2, 30), AR.S1.M4.post[, 5], col=coul[2])
lines(seq(2, 30),  AR.S1.M4.post[, 10], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu>0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), AR.S1.M5.post[, s], col=fade[8])
}
lines(seq(2, 30), AR.S1.M5.post[, 3], col=coul[1])
lines(seq(2, 30), AR.S1.M5.post[, 5], col=coul[2])
lines(seq(2, 30),  AR.S1.M5.post[, 10], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(-25, 0),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\alpha}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lol <- 3*(s-1) + c(1, 2, 3)
  lines(seq(2, 30), AR.S1.M4.A[, lol[2]], col=fade[8])
}
lines(seq(2, 30), AR.S1.M4.A[, 8], col=coul[1])
lines(seq(2, 30), AR.S1.M4.A[, 14], col=coul[2])
lines(seq(2, 30), AR.S1.M4.A[, 29], col=coul[3])


plot(1, type="n", xlim=c(0, 30), ylim=c(-25, 0),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\alpha}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lol <- 3*(s-1) + c(1, 2, 3)
  lines(seq(2, 30), AR.S1.M5.A[, lol[2]], col=fade[8])
}
lines(seq(2, 30), AR.S1.M5.A[, 8], col=coul[1])
lines(seq(2, 30), AR.S1.M5.A[, 14], col=coul[2])
lines(seq(2, 30), AR.S1.M5.A[, 29], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 30),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\tau}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lol <- 3*(s-1) + c(1, 2, 3)
  lines(seq(2, 30), AR.S1.M5.tau[, lol[2]], col=fade[8])
}
lines(seq(2, 30), AR.S1.M5.tau[, 8], col=coul[1])
lines(seq(2, 30), AR.S1.M5.tau[, 14], col=coul[2])
lines(seq(2, 30), AR.S1.M5.tau[, 29], col=coul[3])



# plotting AR sim 2
plot(1, type="n", xlim=c(0, 30), ylim=c(-20, 30),
     bty="n", xlab="Time", ylab=expression(xi[t]),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(0, 30), c(0, cumsum(ARsim2[, s])), col=fade[8])
}
lines(seq(0, 30), c(0, cumsum(ARsim2[, 3])), col=coul[1])
lines(seq(0, 30), c(0, cumsum(ARsim2[, 5])), col=coul[2])
lines(seq(0, 30), c(0, cumsum(ARsim2[, 10])), col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu>0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), AR.S2.M1.post[, s], col=fade[8])
}
lines(seq(2, 30), AR.S2.M1.post[, 3], col=coul[1])
lines(seq(2, 30), AR.S2.M1.post[, 5], col=coul[2])
lines(seq(2, 30),  AR.S2.M1.post[, 10], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu>0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), AR.S2.M2.post[, s], col=fade[8])
}
lines(seq(2, 30), AR.S2.M2.post[, 3], col=coul[1])
lines(seq(2, 30), AR.S2.M2.post[, 5], col=coul[2])
lines(seq(2, 30),  AR.S2.M2.post[, 10], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu>0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), AR.S2.M3.post[, s], col=fade[8])
}
lines(seq(2, 30), AR.S2.M3.post[, 3], col=coul[1])
lines(seq(2, 30), AR.S2.M3.post[, 5], col=coul[2])
lines(seq(2, 30), AR.S2.M3.post[, 10], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu>0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), AR.S2.M4.post[, s], col=fade[8])
}
lines(seq(2, 30), AR.S2.M4.post[, 3], col=coul[1])
lines(seq(2, 30), AR.S2.M4.post[, 5], col=coul[2])
lines(seq(2, 30),  AR.S2.M4.post[, 10], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu>0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), AR.S2.M5.post[, s], col=fade[8])
}
lines(seq(2, 30), AR.S2.M5.post[, 3], col=coul[1])
lines(seq(2, 30), AR.S2.M5.post[, 5], col=coul[2])
lines(seq(2, 30),  AR.S2.M5.post[, 10], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(-10, 0),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\alpha}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lol <- 3*(s-1) + c(1, 2, 3)
  lines(seq(2, 30), AR.S2.M4.A[, lol[2]], col=fade[8])
}
lines(seq(2, 30), AR.S2.M4.A[, 8], col=coul[1])
lines(seq(2, 30), AR.S2.M4.A[, 14], col=coul[2])
lines(seq(2, 30), AR.S2.M4.A[, 29], col=coul[3])


plot(1, type="n", xlim=c(0, 30), ylim=c(-10, 5),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\alpha}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lol <- 3*(s-1) + c(1, 2, 3)
  lines(seq(2, 30), AR.S2.M5.A[, lol[2]], col=fade[8])
}
lines(seq(2, 30), AR.S2.M5.A[, 8], col=coul[1])
lines(seq(2, 30), AR.S2.M5.A[, 14], col=coul[2])
lines(seq(2, 30), AR.S2.M5.A[, 29], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 30),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\tau}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lol <- 3*(s-1) + c(1, 2, 3)
  lines(seq(2, 30), AR.S2.M5.tau[, lol[2]], col=fade[8])
}
lines(seq(2, 30), AR.S2.M5.tau[, 8], col=coul[1])
lines(seq(2, 30), AR.S2.M5.tau[, 14], col=coul[2])
lines(seq(2, 30), AR.S2.M5.tau[, 29], col=coul[3])


### doing simulation 3 of the ARsim
plot(1, type="n", xlim=c(0, 30), ylim=c(-25, 30),
     bty="n", xlab="Time", ylab=expression(xi[t]),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(0, 30), c(0, cumsum(ARsim3[, s])), col=fade[8])
}
lines(seq(0, 30), c(0, cumsum(ARsim3[, 1])), col=coul[1])
lines(seq(0, 30), c(0, cumsum(ARsim3[, 6])), col=coul[2])
lines(seq(0, 30), c(0, cumsum(ARsim3[, 26])), col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu>0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), AR.S3.M1.post[, s], col=fade[8])
}
lines(seq(2, 30), AR.S3.M1.post[, 1], col=coul[1])
lines(seq(2, 30), AR.S3.M1.post[, 6], col=coul[2])
lines(seq(2, 30),  AR.S3.M1.post[, 26], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu>0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), AR.S3.M2.post[, s], col=fade[8])
}
lines(seq(2, 30), AR.S3.M2.post[, 1], col=coul[1])
lines(seq(2, 30), AR.S3.M2.post[, 6], col=coul[2])
lines(seq(2, 30),  AR.S3.M2.post[, 26], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu>0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), AR.S3.M3.post[, s], col=fade[8])
}
lines(seq(2, 30), AR.S3.M3.post[, 1], col=coul[1])
lines(seq(2, 30), AR.S3.M3.post[, 6], col=coul[2])
lines(seq(2, 30), AR.S3.M3.post[, 26], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu>0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), AR.S3.M4.post[, s], col=fade[8])
}
lines(seq(2, 30), AR.S3.M4.post[, 1], col=coul[1])
lines(seq(2, 30), AR.S3.M4.post[, 6], col=coul[2])
lines(seq(2, 30),  AR.S3.M4.post[, 26], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=expression(P(mu>0)),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lines(seq(2, 30), AR.S3.M5.post[, s], col=fade[8])
}
lines(seq(2, 30), AR.S3.M5.post[, 1], col=coul[1])
lines(seq(2, 30), AR.S3.M5.post[, 6], col=coul[2])
lines(seq(2, 30),  AR.S3.M5.post[, 26], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(-10, 0),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\alpha}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lol <- 3*(s-1) + c(1, 2, 3)
  lines(seq(2, 30), AR.S3.M4.A[, lol[2]], col=fade[8])
}
lines(seq(2, 30), AR.S3.M4.A[, 2], col=coul[1])
lines(seq(2, 30), AR.S3.M4.A[, 17], col=coul[2])
lines(seq(2, 30), AR.S3.M4.A[, 77], col=coul[3])


plot(1, type="n", xlim=c(0, 30), ylim=c(-10, 5),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\alpha}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lol <- 3*(s-1) + c(1, 2, 3)
  lines(seq(2, 30), AR.S3.M5.A[, lol[2]], col=fade[8])
}
lines(seq(2, 30), AR.S3.M5.A[, 2], col=coul[1])
lines(seq(2, 30), AR.S3.M5.A[, 17], col=coul[2])
lines(seq(2, 30), AR.S3.M5.A[, 77], col=coul[3])

plot(1, type="n", xlim=c(0, 30), ylim=c(0, 30),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\tau}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(s in 1:30){
  lol <- 3*(s-1) + c(1, 2, 3)
  lines(seq(2, 30), AR.S3.M5.tau[, lol[2]], col=fade[8])
}
lines(seq(2, 30), AR.S3.M5.tau[, 2], col=coul[1])
lines(seq(2, 30), AR.S3.M5.tau[, 17], col=coul[2])
lines(seq(2, 30), AR.S3.M5.tau[, 77], col=coul[3])
