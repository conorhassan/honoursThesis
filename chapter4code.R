# loading libraries
library(rjags)
library(MCMCpack)
library(RColorBrewer)
library(latex2exp)

# creating colour scheme(s) 
coul <- brewer.pal(8, "Dark2")
coul.t <- add.alpha(coul, 0.1)

# Chaper 4 code

# triangle simplex function
simplex <- function(x, y, label = expression(z[t]^1, z[t]^2, z[t]^3),...){
  op <- par(mar=c(1,0,0,0) + 0.1, pty="s")
  plot(x=c(-0.2,1.2), y=c(-0.2,1.2), type="n", axes=FALSE, xlab="", ylab="")
  points(x=c(0,0.5,1,0), y=c(0,0.5*sqrt(3),0,0), type="l")
  points(x=c(0.5, 0.5), y=c(0, sqrt(3)/6), type="l", lty="dashed")
  points(x=c(0.25, 0.5), y=c(0.25*sqrt(3), sqrt(3)/6), type="l", lty="dashed")
  points(x=c(0.75, 0.5), y=c(0.25*sqrt(3), sqrt(3)/6), type="l", lty="dashed")
  xx <- 0.5*(1-x+y)
  yy <- 0.5*sqrt(3)*(1-x-y)
  points(x=xx, y=yy, ..., pch=19)
  if (!is.null(label)) {
    text(x=0.5, y=0.5*sqrt(3), pos=3, labels=label[3])
    text(x=0.0, y=0.0, pos=2, labels=label[1])
    text(x=1.0, y=0.0, pos=4, labels=label[2])
  }
  
  # restore plotting parameters
  par(op)
}

# setting RNG seed
set.seed(7)

# plotting Figure 4.1
instances = rdirichlet(30, 20*c(0.4, 0.35, 0.25))
simplex(x=instances[,1], y=instances[,2], type="p", lwd=0.7, col=coul[1])
points(x=c(0.5, 0.5), y=c(0, sqrt(3)/6), type="l", lty="dashed")
points(x=c(0.25, 0.5), y=c(0.25*sqrt(3), sqrt(3)/6), type="l", lty="dashed")
points(x=c(0.75, 0.5), y=c(0.25*sqrt(3), sqrt(3)/6), type="l", lty="dashed")
points(x=c(0.475), y=c(0.125*sqrt(3)), pch=18, cex=2)
instances = rdirichlet(30, 20*c(0.4, 0.35, 0.25))
simplex(x=instances[,1], y=instances[,2], type="p", lwd=0.7, col=coul[2])
points(x=c(0.5, 0.5), y=c(0, sqrt(3)/6), type="l", lty="dashed")
points(x=c(0.25, 0.5), y=c(0.25*sqrt(3), sqrt(3)/6), type="l", lty="dashed")
points(x=c(0.75, 0.5), y=c(0.25*sqrt(3), sqrt(3)/6), type="l", lty="dashed")
points(x=c(0.475), y=c(0.125*sqrt(3)), pch=18, cex=2)
instances = rdirichlet(30, 20*c(0.4, 0.35, 0.25))
simplex(x=instances[,1], y=instances[,2], type="p", lwd=0.7, col=coul[3])
points(x=c(0.5, 0.5), y=c(0, sqrt(3)/6), type="l", lty="dashed")
points(x=c(0.25, 0.5), y=c(0.25*sqrt(3), sqrt(3)/6), type="l", lty="dashed")
points(x=c(0.75, 0.5), y=c(0.25*sqrt(3), sqrt(3)/6), type="l", lty="dashed")
points(x=c(0.475), y=c(0.125*sqrt(3)), pch=18, cex=2)
instances = rdirichlet(30, 20*c(0.4, 0.35, 0.25))
simplex(x=instances[,1], y=instances[,2], type="p", lwd=0.7, col=coul[4])
points(x=c(0.5, 0.5), y=c(0, sqrt(3)/6), type="l", lty="dashed")
points(x=c(0.25, 0.5), y=c(0.25*sqrt(3), sqrt(3)/6), type="l", lty="dashed")
points(x=c(0.75, 0.5), y=c(0.25*sqrt(3), sqrt(3)/6), type="l", lty="dashed")
points(x=c(0.475), y=c(0.125*sqrt(3)), pch=18, cex=2)
instances = rdirichlet(30, 20*c(0.4, 0.35, 0.25))
simplex(x=instances[,1], y=instances[,2], type="p", lwd=0.7, col=coul[5])
points(x=c(0.5, 0.5), y=c(0, sqrt(3)/6), type="l", lty="dashed")
points(x=c(0.25, 0.5), y=c(0.25*sqrt(3), sqrt(3)/6), type="l", lty="dashed")
points(x=c(0.75, 0.5), y=c(0.25*sqrt(3), sqrt(3)/6), type="l", lty="dashed")
points(x=c(0.475), y=c(0.125*sqrt(3)), pch=18, cex=2)
instances = rdirichlet(30, 20*c(0.4, 0.35, 0.25))
simplex(x=instances[,1], y=instances[,2], type="p", lwd=0.7, col=coul[6])
points(x=c(0.5, 0.5), y=c(0, sqrt(3)/6), type="l", lty="dashed")
points(x=c(0.25, 0.5), y=c(0.25*sqrt(3), sqrt(3)/6), type="l", lty="dashed")
points(x=c(0.75, 0.5), y=c(0.25*sqrt(3), sqrt(3)/6), type="l", lty="dashed")
points(x=c(0.475), y=c(0.125*sqrt(3)), pch=18, cex=2)

# creating Figure 4.2
set.seed(7)
dirpost <- matrix(NA, nrow=29, ncol=90)
pi1est <- matrix(NA, nrow=29, ncol=30)
pi2est <- matrix(NA, nrow=29, ncol=30)
pi3est <- matrix(NA, nrow=29, ncol=30)
lambdaest <- matrix(NA, nrow=29, ncol=30)
for(s in 1:30){
  current_dir <- rdirichlet(30, 20*c(0.40, 0.35, 0.25))
  column_index = 3*(s-1) + c(1, 2, 3)
  for(t in 2:30){
    dir.model = jags.model(file = "dir_lambda_unkn.txt",
                           data = list(y=current_dir[1:t, 1:3], 
                                       n=t,
                                       K=3), 
                           n.chains = 3,
                           n.adapt = 1000)  
    dir.output = coda.samples(model=dir.model, 
                              variable.names=c("lambda",
                                               "pi"),
                              n.iter=10000)
    dir.output.m <- as.matrix(dir.output[, c(2, 3, 4)])
    logical.matrix <- matrix(NA, nrow=dim(dir.output.m)[1], ncol=dim(dir.output.m)[2])
    row_max <- apply(dir.output.m, c(1), max)
    for(i in 1:nrow(logical.matrix)){
      for(j in 1:ncol(logical.matrix)){
        if(dir.output.m[i, j]==row_max[i]){
          logical.matrix[i, j] <- TRUE
        } else {
          logical.matrix[i, j] <- FALSE
        }
      }
    }
    dirpost[(t-1), column_index] <- colSums(logical.matrix) / 
      dim(logical.matrix)[1]
    lambdaest[(t-1), s] <- summary(dir.output)$quantiles[1, 3]
    pi1est[(t-1), s] <- summary(dir.output)$quantiles[2, 3]
    pi2est[(t-1), s] <- summary(dir.output)$quantiles[3, 3]
    pi3est[(t-1), s] <- summary(dir.output)$quantiles[4, 3]
  }
}
set.seed(7)
plot(1, type="n", xlim=c(0, 30), ylim=c(0, 15),
     bty="n", xlab="Time", ylab=expression(xi[t]),
     xaxt = "n", yaxt= "n", cex.lab=0.7)
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(i in 1:30){
  current_dir <- rdirichlet(30, 20*c(0.40, 0.35, 0.25))
  cumsum_current <- apply(current_dir, 2, cumsum)
  lines(seq(1, 30), cumsum_current[, 1], "l", col=coul.t[1])
  lines(seq(1, 30), cumsum_current[, 2], "l", col=coul.t[2])
  lines(seq(1, 30), cumsum_current[, 3], "l", col=coul.t[3])
}
plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab="Predicted probability \ncandidate 1 wins",
     xaxt = "n", yaxt= "n", cex.lab=0.7)
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(i in 1:30){
  column_index = 3*(i-1) + c(1, 2, 3)
  lines(seq(2, 30), dirpost[, column_index[1]], "l", col=coul.t[1])
}
plot(1, type="n", xlim=c(0, 30), ylim=c(0, 60),
     bty="n", xlab="Time", ylab=expression(lambda),
     xaxt = "n", yaxt= "n", cex.lab=0.7)
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(i in 1:30){
  lines(seq(2, 30), lambdaest[, i], "l", col=coul.t[4])
}
plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab="Predicted probability \ncandidate 2 wins",
     xaxt = "n", yaxt= "n", cex.lab=0.7)
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(i in 1:30){
  column_index = 3*(i-1) + c(1, 2, 3)
  lines(seq(2, 30), dirpost[, column_index[2]], "l", col=coul.t[2])
}
plot(1, type="n", xlim=c(0, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab="Predicted probability \ncandidate 3 wins",
     xaxt = "n", yaxt= "n", cex.lab=0.7)
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
for(i in 1:30){
  column_index = 3*(i-1) + c(1, 2, 3)
  lines(seq(2, 30), dirpost[, column_index[3]], "l", col=coul.t[3])
}

######
# creating dataframes for Figures 4.2-4.8
######
set.seed(2020)
# strength of fake news effect nu function
scenario_5_effect <- function(time_seq, tau, A, beta){
  T <- length(time_seq)
  alpha <- numeric(T)
  for(t in 1:T){
    if(t < tau){
      alpha[t] <- 0
    }
    if(t >= tau){
      alpha[t] <- A / exp(beta * (t - tau))
    }
  }
  return(alpha)
}
# JAGS initial value function
init_function <- function(){
  return(list(lambda=runif(0, 50),
              .RNG.name="base::Super-Duper",
              .RNG.seed=7))
}
# Parameter set 1
true_pi <- c(0.6, 0.2, 0.2)
true_lambda <- 15
rho = c(0.1, 0.7, 0.2)
tau = 15
A = 1
beta = 0.1
# simulation 1 fake news effect
scen1_effect <- scenario_5_effect(seq(1, 30), tau=tau, A=A, beta=beta)
# initialising simulations
scen1 <- matrix(NA, nrow = 30, ncol = 9)
for(i in 1:3){
  column_ind <- 3*(i-1) + c(1, 2, 3)
  for(t in 1:30){
    scen1[t, column_ind] <- rdirichlet(1, true_lambda*(scen1_effect[t]*rho+(1-scen1_effect[t])*true_pi))[1, ]
  }
}
# initialising data.frames
scen1_m1_pi <- matrix(NA, nrow=30, ncol=9)
scen1_m1_lambda <- matrix(NA, nrow=30, ncol=3)
scen1_m1_posterior <- matrix(NA, nrow=30, ncol=9)
scen1_m2_pi <- matrix(NA, nrow=30, ncol=9)
scen1_m2_lambda <- matrix(NA, nrow=30, ncol=3)
scen1_m2_posterior <- matrix(NA, nrow=30, ncol=9)
scen1_m3_pi <- matrix(NA, nrow=30, ncol=9)
scen1_m3_lambda <- matrix(NA, nrow=30, ncol=3)
scen1_m3_rho <- matrix(NA, nrow=30, ncol=9)
scen1_m3_A <- matrix(NA, nrow=30, ncol=3)
scen1_m3_beta = matrix(NA, nrow=30, ncol=3)
scen1_m3_posterior <- matrix(NA, nrow=30, ncol=9)
scen1_m4_pi <- matrix(NA, nrow=30, ncol=9)
scen1_m4_rho <- matrix(NA, nrow=30, ncol=9)
scen1_m4_lambda <- matrix(NA, nrow=30, ncol=3)
scen1_m4_A <- matrix(NA, nrow=30, ncol=3)
scen1_m4_beta = matrix(NA, nrow=30, ncol=3)
scen1_m4_tau <- matrix(NA, nrow=30, ncol=3)
scen1_m4_posterior <- matrix(NA, nrow=30, ncol=9)
scen1_m5_pi <- matrix(NA, nrow=30, ncol=9)
scen1_m5_rho <- matrix(NA, nrow=30, ncol=9)
scen1_m5_lambda <- matrix(NA, nrow=30, ncol=3)
scen1_m5_beta = matrix(NA, nrow=30, ncol=3)
scen1_m5_posterior <- matrix(NA, nrow=30, ncol=9)
for(i in 1:3){
  column_ind = 3*(i-1) + c(1, 2, 3)
  for(k in 2:30){
    # M(3, 1)
    model = jags.model(file = "dir_lambda_unkn.txt",
                       data = list(y=scen1[1:k, column_ind], 
                                   n=k, 
                                   K=3),
                       n.chains = 3,
                       inits = init_function,
                       n.adapt = 1000)  
    # creating output object
    output = coda.samples(model=model, 
                          variable.names=c("pi",
                                           "lambda"),
                          n.iter=10000)
    scen1_m1_lambda[k, i] <- summary(output)$quantiles[1, 3]
    scen1_m1_pi[k, column_ind[1]] <- summary(output)$quantiles[2, 3]
    scen1_m1_pi[k, column_ind[2]] <- summary(output)$quantiles[3, 3]
    scen1_m1_pi[k, column_ind[3]] <- summary(output)$quantiles[4, 3]
    output_matrix <- as.matrix(output[, c(2, 3, 4)])
    logical_matrix <- matrix(NA, nrow=dim(output_matrix)[1], ncol=dim(output_matrix)[2])
    row_max <- apply(output_matrix, c(1), max)
    for(m in 1:nrow(logical_matrix)){
      for(j in 1:ncol(logical_matrix)){
        if(output_matrix[m, j]==row_max[m]){
          logical_matrix[m, j] <- TRUE
        } else {
          logical_matrix[m, j] <- FALSE
        }
      }
    }
    scen1_m1_posterior[k, column_ind] <- colSums(logical_matrix) / 
      dim(logical_matrix)[1]
    
    # M(3, 2)
    model = jags.model(file = "dirichlet_robust_hack.txt",
                       data = list(y=scen1[1:k, column_ind], 
                                   n=k, 
                                   K=3), 
                       n.chains = 3,
                       inits = init_function,
                       n.adapt = 1000)  
    # creating output object
    output = coda.samples(model=model, 
                          variable.names=c("pi",
                                           "lambda"),
                          n.iter=10000)
    scen1_m2_lambda[k, i] <- summary(output)$quantiles[1, 3]
    scen1_m2_pi[k, column_ind[1]] <- summary(output)$quantiles[2, 3]
    scen1_m2_pi[k, column_ind[2]] <- summary(output)$quantiles[3, 3]
    scen1_m2_pi[k, column_ind[3]] <- summary(output)$quantiles[4, 3]
    output_matrix <- as.matrix(output[, c(2, 3, 4)])
    logical_matrix <- matrix(NA, nrow=dim(output_matrix)[1], ncol=dim(output_matrix)[2])
    row_max <- apply(output_matrix, c(1), max)
    for(m in 1:nrow(logical_matrix)){
      for(j in 1:ncol(logical_matrix)){
        if(output_matrix[m, j]==row_max[m]){
          logical_matrix[m, j] <- TRUE
        } else {
          logical_matrix[m, j] <- FALSE
        }
      }
    }
    scen1_m2_posterior[k, column_ind] <- colSums(logical_matrix) / 
      dim(logical_matrix)[1]
    
    # M(3, 4)
    if(k >= 15){
      model = jags.model(file = "fake_param_est_dir.txt",
                         data = list(y=scen1[1:k, column_ind], 
                                     T=k, 
                                     K=3,
                                     tau = tau,
                                     tau_adj = tau-1), 
                         n.chains = 3,
                         inits = init_function,
                         n.adapt = 10000)  
      # creating output object
      output = coda.samples(model=model, 
                            variable.names=c("pi",
                                             "rho", 
                                             "A",
                                             "beta", 
                                             "lambda"),
                            n.iter=30000)
      scen1_m3_A[k, i] <- summary(output)$quantiles[1, 3]
      scen1_m3_beta[k, i] <- summary(output)$quantiles[2, 3]
      scen1_m3_lambda[k, i] <- summary(output)$quantiles[3, 3]
      scen1_m3_pi[k, column_ind[1]] <- summary(output)$quantiles[4, 3]
      scen1_m3_pi[k, column_ind[2]] <- summary(output)$quantiles[5, 3]
      scen1_m3_pi[k, column_ind[3]] <- summary(output)$quantiles[6, 3]
      scen1_m3_rho[k, column_ind[1]] <- summary(output)$quantiles[7, 3]
      scen1_m3_rho[k, column_ind[2]] <- summary(output)$quantiles[8, 3]
      scen1_m3_rho[k, column_ind[3]] <- summary(output)$quantiles[9, 3]
      output_matrix <- as.matrix(output[, c(4, 5, 6)])
      logical_matrix <- matrix(NA, nrow=dim(output_matrix)[1], ncol=dim(output_matrix)[2])
      row_max <- apply(output_matrix, c(1), max)
      for(m in 1:nrow(logical_matrix)){
        for(j in 1:ncol(logical_matrix)){
          if(output_matrix[m, j]==row_max[m]){
            logical_matrix[m, j] <- TRUE
          } else {
            logical_matrix[m, j] <- FALSE
          }
        }
      }
      scen1_m3_posterior[k, column_ind] <- colSums(logical_matrix) / 
        dim(logical_matrix)[1]
    }
    # M(3, 5)
    model = jags.model(file = "estimating_tau.txt",
                       data = list(y=scen1[1:k, column_ind], 
                                   T=k, 
                                   K=3), 
                       n.chains = 3,
                       inits = init_function,
                       n.adapt = 1000)  
    # creating output object
    output = coda.samples(model=model, 
                          variable.names=c("pi",
                                           "rho", 
                                           "A",
                                           "beta", 
                                           "lambda",
                                           "tau"),
                          n.iter=10000)
    scen1_m4_A[k, i] <- summary(output)$quantiles[1, 3]
    scen1_m4_beta[k, i] <- summary(output)$quantiles[2, 3]
    scen1_m4_lambda[k, i] <- summary(output)$quantiles[3, 3]
    scen1_m4_pi[k, column_ind[1]] <- summary(output)$quantiles[4, 3]
    scen1_m4_pi[k, column_ind[2]] <- summary(output)$quantiles[5, 3]
    scen1_m4_pi[k, column_ind[3]] <- summary(output)$quantiles[6, 3]
    scen1_m4_rho[k, column_ind[1]] <- summary(output)$quantiles[7, 3]
    scen1_m4_rho[k, column_ind[2]] <- summary(output)$quantiles[8, 3]
    scen1_m4_rho[k, column_ind[3]] <- summary(output)$quantiles[9, 3]
    scen1_m4_tau[k, i] <- summary(output)$quantiles[10, 3]
    output_matrix <- as.matrix(output[, c(4, 5, 6)])
    logical_matrix <- matrix(NA, nrow=dim(output_matrix)[1], ncol=dim(output_matrix)[2])
    row_max <- apply(output_matrix, c(1), max)
    for(m in 1:nrow(logical_matrix)){
      for(j in 1:ncol(logical_matrix)){
        if(output_matrix[m, j]==row_max[m]){
          logical_matrix[m, j] <- TRUE
        } else {
          logical_matrix[m, j] <- FALSE
        }
      }
    }
    scen1_m4_posterior[k, column_ind] <- colSums(logical_matrix) / 
      dim(logical_matrix)[1]
    # M(3, 3)
    if(k >= 15){
      model = jags.model(file = "alpha_known.txt",
                         data = list(y=scen1[1:k, column_ind], 
                                     T=k, 
                                     K=3,
                                     tau = tau,
                                     tau_adj = tau-1,
                                     A=1), 
                         n.chains = 3,
                         inits = init_function,
                         n.adapt = 10000)  
      # creating output object
      output = coda.samples(model=model, 
                            variable.names=c("pi",
                                             "rho", 
                                             "beta", 
                                             "lambda"),
                            n.iter=30000)
      scen1_m5_beta[k, i] <- summary(output)$quantiles[1, 3]
      scen1_m5_lambda[k, i] <- summary(output)$quantiles[2, 3]
      scen1_m5_pi[k, column_ind[1]] <- summary(output)$quantiles[3, 3]
      scen1_m5_pi[k, column_ind[2]] <- summary(output)$quantiles[4, 3]
      scen1_m5_pi[k, column_ind[3]] <- summary(output)$quantiles[5, 3]
      scen1_m5_rho[k, column_ind[1]] <- summary(output)$quantiles[6, 3]
      scen1_m5_rho[k, column_ind[2]] <- summary(output)$quantiles[7, 3]
      scen1_m5_rho[k, column_ind[3]] <- summary(output)$quantiles[8, 3]
      output_matrix <- as.matrix(output[, c(3, 4, 5)])
      logical_matrix <- matrix(NA, nrow=dim(output_matrix)[1], ncol=dim(output_matrix)[2])
      row_max <- apply(output_matrix, c(1), max)
      for(m in 1:nrow(logical_matrix)){
        for(j in 1:ncol(logical_matrix)){
          if(output_matrix[m, j]==row_max[m]){
            logical_matrix[m, j] <- TRUE
          } else {
            logical_matrix[m, j] <- FALSE
          }
        }
      }
      scen1_m5_posterior[k, column_ind] <- colSums(logical_matrix) / 
        dim(logical_matrix)[1]
    }
  }
}
# Updating beta for simulation set 2
beta = 0.3
# simulation 2 fake news effect
scen2_effect <- scenario_5_effect(seq(1, 30), tau=tau, A=A, beta=beta)
# initialising data
scen2 <- matrix(NA, nrow = 30, ncol = 9)
for(i in 1:3){
  column_ind <- 3*(i-1) + c(1, 2, 3)
  for(t in 1:30){
    scen2[t, column_ind] <- rdirichlet(1, true_lambda*(scen2_effect[t]*rho+(1-scen2_effect[t])*true_pi))[1, ]
  }
}
# initialising data.frames
scen2_m1_pi <- matrix(NA, nrow=30, ncol=9)
scen2_m1_lambda <- matrix(NA, nrow=30, ncol=3)
scen2_m1_posterior <- matrix(NA, nrow=30, ncol=9)
scen2_m2_pi <- matrix(NA, nrow=30, ncol=9)
scen2_m2_lambda <- matrix(NA, nrow=30, ncol=3)
scen2_m2_posterior <- matrix(NA, nrow=30, ncol=9)
scen2_m3_pi <- matrix(NA, nrow=30, ncol=9)
scen2_m3_lambda <- matrix(NA, nrow=30, ncol=3)
scen2_m3_rho <- matrix(NA, nrow=30, ncol=9)
scen2_m3_A <- matrix(NA, nrow=30, ncol=3)
scen2_m3_beta = matrix(NA, nrow=30, ncol=3)
scen2_m3_posterior <- matrix(NA, nrow=30, ncol=9)
scen2_m4_pi <- matrix(NA, nrow=30, ncol=9)
scen2_m4_rho <- matrix(NA, nrow=30, ncol=9)
scen2_m4_lambda <- matrix(NA, nrow=30, ncol=3)
scen2_m4_A <- matrix(NA, nrow=30, ncol=3)
scen2_m4_beta = matrix(NA, nrow=30, ncol=3)
scen2_m4_tau <- matrix(NA, nrow=30, ncol=3)
scen2_m4_posterior <- matrix(NA, nrow=30, ncol=9)
scen2_m5_pi <- matrix(NA, nrow=30, ncol=9)
scen2_m5_rho <- matrix(NA, nrow=30, ncol=9)
scen2_m5_lambda <- matrix(NA, nrow=30, ncol=3)
scen2_m5_beta = matrix(NA, nrow=30, ncol=3)
scen2_m5_posterior <- matrix(NA, nrow=30, ncol=9)
for(i in 1:3){
  column_ind = 3*(i-1) + c(1, 2, 3)
  for(k in 2:30){
    # M(3, 1)
    model = jags.model(file = "dir_lambda_unkn.txt",
                       data = list(y=scen2[1:k, column_ind], 
                                   n=k, 
                                   K=3),
                       n.chains = 3,
                       inits = init_function,
                       n.adapt = 1000)  
    # creating output object
    output = coda.samples(model=model, 
                          variable.names=c("pi",
                                           "lambda"),
                          n.iter=10000)
    scen2_m1_lambda[k, i] <- summary(output)$quantiles[1, 3]
    scen2_m1_pi[k, column_ind[1]] <- summary(output)$quantiles[2, 3]
    scen2_m1_pi[k, column_ind[2]] <- summary(output)$quantiles[3, 3]
    scen2_m1_pi[k, column_ind[3]] <- summary(output)$quantiles[4, 3]
    output_matrix <- as.matrix(output[, c(2, 3, 4)])
    logical_matrix <- matrix(NA, nrow=dim(output_matrix)[1], ncol=dim(output_matrix)[2])
    row_max <- apply(output_matrix, c(1), max)
    for(m in 1:nrow(logical_matrix)){
      for(j in 1:ncol(logical_matrix)){
        if(output_matrix[m, j]==row_max[m]){
          logical_matrix[m, j] <- TRUE
        } else {
          logical_matrix[m, j] <- FALSE
        }
      }
    }
    scen2_m1_posterior[k, column_ind] <- colSums(logical_matrix) / 
      dim(logical_matrix)[1]
    # M(3, 2)
    model = jags.model(file = "dirichlet_robust_hack.txt",
                       data = list(y=scen2[1:k, column_ind], 
                                   n=k, 
                                   K=3), 
                       n.chains = 3,
                       inits = init_function,
                       n.adapt = 1000)  
    # creating output object
    output = coda.samples(model=model, 
                          variable.names=c("pi",
                                           "lambda"),
                          n.iter=10000)
    scen2_m2_lambda[k, i] <- summary(output)$quantiles[1, 3]
    scen2_m2_pi[k, column_ind[1]] <- summary(output)$quantiles[2, 3]
    scen2_m2_pi[k, column_ind[2]] <- summary(output)$quantiles[3, 3]
    scen2_m2_pi[k, column_ind[3]] <- summary(output)$quantiles[4, 3]
    output_matrix <- as.matrix(output[, c(2, 3, 4)])
    logical_matrix <- matrix(NA, nrow=dim(output_matrix)[1], ncol=dim(output_matrix)[2])
    row_max <- apply(output_matrix, c(1), max)
    for(m in 1:nrow(logical_matrix)){
      for(j in 1:ncol(logical_matrix)){
        if(output_matrix[m, j]==row_max[m]){
          logical_matrix[m, j] <- TRUE
        } else {
          logical_matrix[m, j] <- FALSE
        }
      }
    }
    scen2_m2_posterior[k, column_ind] <- colSums(logical_matrix) / 
      dim(logical_matrix)[1]
    # M(3, 4)
    if(k >= 15){
      model = jags.model(file = "fake_param_est_dir.txt",
                         data = list(y=scen2[1:k, column_ind], 
                                     T=k, 
                                     K=3,
                                     tau = 15,
                                     tau_adj = 14), 
                         n.chains = 3,
                         inits = init_function,
                         n.adapt = 10000)  
      # creating output object
      output = coda.samples(model=model, 
                            variable.names=c("pi",
                                             "rho", 
                                             "A",
                                             "beta", 
                                             "lambda"),
                            n.iter=30000)
      scen2_m3_A[k, i] <- summary(output)$quantiles[1, 3]
      scen2_m3_beta[k, i] <- summary(output)$quantiles[2, 3]
      scen2_m3_lambda[k, i] <- summary(output)$quantiles[3, 3]
      scen2_m3_pi[k, column_ind[1]] <- summary(output)$quantiles[4, 3]
      scen2_m3_pi[k, column_ind[2]] <- summary(output)$quantiles[5, 3]
      scen2_m3_pi[k, column_ind[3]] <- summary(output)$quantiles[6, 3]
      scen2_m3_rho[k, column_ind[1]] <- summary(output)$quantiles[7, 3]
      scen2_m3_rho[k, column_ind[2]] <- summary(output)$quantiles[8, 3]
      scen2_m3_rho[k, column_ind[3]] <- summary(output)$quantiles[9, 3]
      output_matrix <- as.matrix(output[, c(4, 5, 6)])
      logical_matrix <- matrix(NA, nrow=dim(output_matrix)[1], ncol=dim(output_matrix)[2])
      row_max <- apply(output_matrix, c(1), max)
      for(m in 1:nrow(logical_matrix)){
        for(j in 1:ncol(logical_matrix)){
          if(output_matrix[m, j]==row_max[m]){
            logical_matrix[m, j] <- TRUE
          } else {
            logical_matrix[m, j] <- FALSE
          }
        }
      }
      scen2_m3_posterior[k, column_ind] <- colSums(logical_matrix) / 
        dim(logical_matrix)[1]
    }
    # M(3, 5)
    model = jags.model(file = "estimating_tau.txt",
                       data = list(y=scen2[1:k, column_ind], 
                                   T=k, 
                                   K=3), 
                       n.chains = 3,
                       inits = init_function,
                       n.adapt = 1000)  
    # creating output object
    output = coda.samples(model=model, 
                          variable.names=c("pi",
                                           "rho", 
                                           "A",
                                           "beta", 
                                           "lambda",
                                           "tau"),
                          n.iter=10000)
    scen2_m4_A[k, i] <- summary(output)$quantiles[1, 3]
    scen2_m4_beta[k, i] <- summary(output)$quantiles[2, 3]
    scen2_m4_lambda[k, i] <- summary(output)$quantiles[3, 3]
    scen2_m4_pi[k, column_ind[1]] <- summary(output)$quantiles[4, 3]
    scen2_m4_pi[k, column_ind[2]] <- summary(output)$quantiles[5, 3]
    scen2_m4_pi[k, column_ind[3]] <- summary(output)$quantiles[6, 3]
    scen2_m4_rho[k, column_ind[1]] <- summary(output)$quantiles[7, 3]
    scen2_m4_rho[k, column_ind[2]] <- summary(output)$quantiles[8, 3]
    scen2_m4_rho[k, column_ind[3]] <- summary(output)$quantiles[9, 3]
    scen2_m4_tau[k, i] <- summary(output)$quantiles[10, 3]
    output_matrix <- as.matrix(output[, c(4, 5, 6)])
    logical_matrix <- matrix(NA, nrow=dim(output_matrix)[1], ncol=dim(output_matrix)[2])
    row_max <- apply(output_matrix, c(1), max)
    for(m in 1:nrow(logical_matrix)){
      for(j in 1:ncol(logical_matrix)){
        if(output_matrix[m, j]==row_max[m]){
          logical_matrix[m, j] <- TRUE
        } else {
          logical_matrix[m, j] <- FALSE
        }
      }
    }
    scen2_m4_posterior[k, column_ind] <- colSums(logical_matrix) / 
      dim(logical_matrix)[1]
    # M(3, 3)
    if(k >= 15){
      model = jags.model(file = "alpha_known.txt",
                         data = list(y=scen2[1:k, column_ind], 
                                     T=k, 
                                     K=3,
                                     tau = tau,
                                     tau_adj = tau-1,
                                     A=1), 
                         n.chains = 3,
                         inits = init_function,
                         n.adapt = 10000)  
      # creating output object
      output = coda.samples(model=model, 
                            variable.names=c("pi",
                                             "rho", 
                                             "beta", 
                                             "lambda"),
                            n.iter=30000)
      scen2_m5_beta[k, i] <- summary(output)$quantiles[1, 3]
      scen2_m5_lambda[k, i] <- summary(output)$quantiles[2, 3]
      scen2_m5_pi[k, column_ind[1]] <- summary(output)$quantiles[3, 3]
      scen2_m5_pi[k, column_ind[2]] <- summary(output)$quantiles[4, 3]
      scen2_m5_pi[k, column_ind[3]] <- summary(output)$quantiles[5, 3]
      scen2_m5_rho[k, column_ind[1]] <- summary(output)$quantiles[6, 3]
      scen2_m5_rho[k, column_ind[2]] <- summary(output)$quantiles[7, 3]
      scen2_m5_rho[k, column_ind[3]] <- summary(output)$quantiles[8, 3]
      output_matrix <- as.matrix(output[, c(3, 4, 5)])
      logical_matrix <- matrix(NA, nrow=dim(output_matrix)[1], ncol=dim(output_matrix)[2])
      row_max <- apply(output_matrix, c(1), max)
      for(m in 1:nrow(logical_matrix)){
        for(j in 1:ncol(logical_matrix)){
          if(output_matrix[m, j]==row_max[m]){
            logical_matrix[m, j] <- TRUE
          } else {
            logical_matrix[m, j] <- FALSE
          }
        }
      }
      scen2_m5_posterior[k, column_ind] <- colSums(logical_matrix) / 
        dim(logical_matrix)[1]
    }
  }
}
# updating beta for the third simulation
beta = 1
# Simulation 3 fake news effect
scen3_effect <- scenario_5_effect(seq(1, 30), tau=tau, A=A, beta=beta)
# initialising data
scen3 <- matrix(NA, nrow = 30, ncol = 9)
for(i in 1:3){
  column_ind <- 3*(i-1) + c(1, 2, 3)
  for(t in 1:30){
    scen3[t, column_ind] <- rdirichlet(1, true_lambda*(scen3_effect[t]*rho+(1-scen3_effect[t])*true_pi))[1, ]
  }
}
# initialising dataframes
scen3_m1_pi <- matrix(NA, nrow=30, ncol=9)
scen3_m1_lambda <- matrix(NA, nrow=30, ncol=3)
scen3_m1_posterior <- matrix(NA, nrow=30, ncol=9)
scen3_m2_pi <- matrix(NA, nrow=30, ncol=9)
scen3_m2_lambda <- matrix(NA, nrow=30, ncol=3)
scen3_m2_posterior <- matrix(NA, nrow=30, ncol=9)
scen3_m3_pi <- matrix(NA, nrow=30, ncol=9)
scen3_m3_lambda <- matrix(NA, nrow=30, ncol=3)
scen3_m3_rho <- matrix(NA, nrow=30, ncol=9)
scen3_m3_A <- matrix(NA, nrow=30, ncol=3)
scen3_m3_beta = matrix(NA, nrow=30, ncol=3)
scen3_m3_posterior <- matrix(NA, nrow=30, ncol=9)
scen3_m4_pi <- matrix(NA, nrow=30, ncol=9)
scen3_m4_rho <- matrix(NA, nrow=30, ncol=9)
scen3_m4_lambda <- matrix(NA, nrow=30, ncol=3)
scen3_m4_A <- matrix(NA, nrow=30, ncol=3)
scen3_m4_beta = matrix(NA, nrow=30, ncol=3)
scen3_m4_tau <- matrix(NA, nrow=30, ncol=3)
scen3_m4_posterior <- matrix(NA, nrow=30, ncol=9)
scen3_m5_pi <- matrix(NA, nrow=30, ncol=9)
scen3_m5_rho <- matrix(NA, nrow=30, ncol=9)
scen3_m5_lambda <- matrix(NA, nrow=30, ncol=3)
scen3_m5_beta = matrix(NA, nrow=30, ncol=3)
scen3_m5_posterior <- matrix(NA, nrow=30, ncol=9)
for(i in 1:3){
  column_ind = 3*(i-1) + c(1, 2, 3)
  for(k in 2:30){
    # M(3, 1)
    model = jags.model(file = "dir_lambda_unkn.txt",
                       data = list(y=scen3[1:k, column_ind], 
                                   n=k, 
                                   K=3),
                       n.chains = 3,
                       inits = init_function,
                       n.adapt = 1000)  
    # creating output object
    output = coda.samples(model=model, 
                          variable.names=c("pi",
                                           "lambda"),
                          n.iter=10000)
    scen3_m1_lambda[k, i] <- summary(output)$quantiles[1, 3]
    scen3_m1_pi[k, column_ind[1]] <- summary(output)$quantiles[2, 3]
    scen3_m1_pi[k, column_ind[2]] <- summary(output)$quantiles[3, 3]
    scen3_m1_pi[k, column_ind[3]] <- summary(output)$quantiles[4, 3]
    output_matrix <- as.matrix(output[, c(2, 3, 4)])
    logical_matrix <- matrix(NA, nrow=dim(output_matrix)[1], ncol=dim(output_matrix)[2])
    row_max <- apply(output_matrix, c(1), max)
    for(m in 1:nrow(logical_matrix)){
      for(j in 1:ncol(logical_matrix)){
        if(output_matrix[m, j]==row_max[m]){
          logical_matrix[m, j] <- TRUE
        } else {
          logical_matrix[m, j] <- FALSE
        }
      }
    }
    scen3_m1_posterior[k, column_ind] <- colSums(logical_matrix) / 
      dim(logical_matrix)[1]
    # M(3, 2)
    model = jags.model(file = "dirichlet_robust_hack.txt",
                       data = list(y=scen3[1:k, column_ind], 
                                   n=k, 
                                   K=3), 
                       n.chains = 3,
                       inits = init_function,
                       n.adapt = 1000)  
    # creating output object
    output = coda.samples(model=model, 
                          variable.names=c("pi",
                                           "lambda"),
                          n.iter=10000)
    scen3_m2_lambda[k, i] <- summary(output)$quantiles[1, 3]
    scen3_m2_pi[k, column_ind[1]] <- summary(output)$quantiles[2, 3]
    scen3_m2_pi[k, column_ind[2]] <- summary(output)$quantiles[3, 3]
    scen3_m2_pi[k, column_ind[3]] <- summary(output)$quantiles[4, 3]
    output_matrix <- as.matrix(output[, c(2, 3, 4)])
    logical_matrix <- matrix(NA, nrow=dim(output_matrix)[1], ncol=dim(output_matrix)[2])
    row_max <- apply(output_matrix, c(1), max)
    for(m in 1:nrow(logical_matrix)){
      for(j in 1:ncol(logical_matrix)){
        if(output_matrix[m, j]==row_max[m]){
          logical_matrix[m, j] <- TRUE
        } else {
          logical_matrix[m, j] <- FALSE
        }
      }
    }
    scen3_m2_posterior[k, column_ind] <- colSums(logical_matrix) / 
      dim(logical_matrix)[1]
    # M(3, 4)
    if(k >= 15){
      model = jags.model(file = "fake_param_est_dir.txt",
                         data = list(y=scen3[1:k, column_ind], 
                                     T=k, 
                                     K=3,
                                     tau = 15,
                                     tau_adj = 14), 
                         n.chains = 3,
                         inits = init_function,
                         n.adapt = 10000)  
      # creating output object
      output = coda.samples(model=model, 
                            variable.names=c("pi",
                                             "rho", 
                                             "A",
                                             "beta", 
                                             "lambda"),
                            n.iter=30000)
      scen3_m3_A[k, i] <- summary(output)$quantiles[1, 3]
      scen3_m3_beta[k, i] <- summary(output)$quantiles[2, 3]
      scen3_m3_lambda[k, i] <- summary(output)$quantiles[3, 3]
      scen3_m3_pi[k, column_ind[1]] <- summary(output)$quantiles[4, 3]
      scen3_m3_pi[k, column_ind[2]] <- summary(output)$quantiles[5, 3]
      scen3_m3_pi[k, column_ind[3]] <- summary(output)$quantiles[6, 3]
      scen3_m3_rho[k, column_ind[1]] <- summary(output)$quantiles[7, 3]
      scen3_m3_rho[k, column_ind[2]] <- summary(output)$quantiles[8, 3]
      scen3_m3_rho[k, column_ind[3]] <- summary(output)$quantiles[9, 3]
      output_matrix <- as.matrix(output[, c(4, 5, 6)])
      logical_matrix <- matrix(NA, nrow=dim(output_matrix)[1], ncol=dim(output_matrix)[2])
      row_max <- apply(output_matrix, c(1), max)
      for(m in 1:nrow(logical_matrix)){
        for(j in 1:ncol(logical_matrix)){
          if(output_matrix[m, j]==row_max[m]){
            logical_matrix[m, j] <- TRUE
          } else {
            logical_matrix[m, j] <- FALSE
          }
        }
      }
      scen3_m3_posterior[k, column_ind] <- colSums(logical_matrix) / 
        dim(logical_matrix)[1]
    }
    # M(3, 5)
    model = jags.model(file = "estimating_tau.txt",
                       data = list(y=scen3[1:k, column_ind], 
                                   T=k, 
                                   K=3), 
                       n.chains = 3,
                       inits = init_function,
                       n.adapt = 1000)  
    # creating output object
    output = coda.samples(model=model, 
                          variable.names=c("pi",
                                           "rho", 
                                           "A",
                                           "beta", 
                                           "lambda",
                                           "tau"),
                          n.iter=10000)
    scen3_m4_A[k, i] <- summary(output)$quantiles[1, 3]
    scen3_m4_beta[k, i] <- summary(output)$quantiles[2, 3]
    scen3_m4_lambda[k, i] <- summary(output)$quantiles[3, 3]
    scen3_m4_pi[k, column_ind[1]] <- summary(output)$quantiles[4, 3]
    scen3_m4_pi[k, column_ind[2]] <- summary(output)$quantiles[5, 3]
    scen3_m4_pi[k, column_ind[3]] <- summary(output)$quantiles[6, 3]
    scen3_m4_rho[k, column_ind[1]] <- summary(output)$quantiles[7, 3]
    scen3_m4_rho[k, column_ind[2]] <- summary(output)$quantiles[8, 3]
    scen3_m4_rho[k, column_ind[3]] <- summary(output)$quantiles[9, 3]
    scen3_m4_tau[k, i] <- summary(output)$quantiles[10, 3]
    output_matrix <- as.matrix(output[, c(4, 5, 6)])
    logical_matrix <- matrix(NA, nrow=dim(output_matrix)[1], ncol=dim(output_matrix)[2])
    row_max <- apply(output_matrix, c(1), max)
    for(m in 1:nrow(logical_matrix)){
      for(j in 1:ncol(logical_matrix)){
        if(output_matrix[m, j]==row_max[m]){
          logical_matrix[m, j] <- TRUE
        } else {
          logical_matrix[m, j] <- FALSE
        }
      }
    }
    scen3_m4_posterior[k, column_ind] <- colSums(logical_matrix) / 
      dim(logical_matrix)[1]
    # M(3, 3)
    if(k >= 15){
      model = jags.model(file = "alpha_known.txt",
                         data = list(y=scen3[1:k, column_ind], 
                                     T=k, 
                                     K=3,
                                     tau = tau,
                                     tau_adj = tau-1,
                                     A=1), 
                         n.chains = 3,
                         inits = init_function,
                         n.adapt = 10000)  
      # creating output object
      output = coda.samples(model=model, 
                            variable.names=c("pi",
                                             "rho", 
                                             "beta", 
                                             "lambda"),
                            n.iter=30000)
      scen3_m5_beta[k, i] <- summary(output)$quantiles[1, 3]
      scen3_m5_lambda[k, i] <- summary(output)$quantiles[2, 3]
      scen3_m5_pi[k, column_ind[1]] <- summary(output)$quantiles[3, 3]
      scen3_m5_pi[k, column_ind[2]] <- summary(output)$quantiles[4, 3]
      scen3_m5_pi[k, column_ind[3]] <- summary(output)$quantiles[5, 3]
      scen3_m5_rho[k, column_ind[1]] <- summary(output)$quantiles[6, 3]
      scen3_m5_rho[k, column_ind[2]] <- summary(output)$quantiles[7, 3]
      scen3_m5_rho[k, column_ind[3]] <- summary(output)$quantiles[8, 3]
      output_matrix <- as.matrix(output[, c(3, 4, 5)])
      logical_matrix <- matrix(NA, nrow=dim(output_matrix)[1], ncol=dim(output_matrix)[2])
      row_max <- apply(output_matrix, c(1), max)
      for(m in 1:nrow(logical_matrix)){
        for(j in 1:ncol(logical_matrix)){
          if(output_matrix[m, j]==row_max[m]){
            logical_matrix[m, j] <- TRUE
          } else {
            logical_matrix[m, j] <- FALSE
          }
        }
      }
      scen3_m5_posterior[k, column_ind] <- colSums(logical_matrix) / 
        dim(logical_matrix)[1]
    }
  }
}
# plotting figures 4.3 and 4.4
# plotting cumulative sentiment
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(0, 30), ylim=c(0, 15),
     bty="n", xlab="Time", ylab=expression(xi[t]),
     xaxt = "n", yaxt= "n", cex.lab=0.7)
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(0, 30), c(0, cumsum(scen1[, 1])), col=coul[1])
lines(seq(0, 30), c(0, cumsum(scen1[, 2])), col=coul[2])
lines(seq(0, 30), c(0, cumsum(scen1[, 3])), col=coul[3])
lines(seq(0, 30), c(0, cumsum(scen1[, 4])), col=coul[1], lty=3)
lines(seq(0, 30), c(0, cumsum(scen1[, 5])), col=coul[2], lty=3)
lines(seq(0, 30), c(0, cumsum(scen1[, 6])), col=coul[3], lty=3)
lines(seq(0, 30), c(0, cumsum(scen1[, 7])), col=coul[1], lty=5)
lines(seq(0, 30), c(0, cumsum(scen1[, 8])), col=coul[2], lty=5)
lines(seq(0, 30), c(0, cumsum(scen1[, 9])), col=coul[3], lty=5)
# posterior of each model
load("dirichlet_data.Rdata")
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab="Predicted probability \n candidate wins",
     xaxt = "n", yaxt= "n", cex.lab=0.7)
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen1_m1_posterior[2:30, 1], col=coul[1])
lines(seq(2, 30), scen1_m1_posterior[2:30, 2], col=coul[2])
lines(seq(2, 30), scen1_m1_posterior[2:30, 3], col=coul[3])
lines(seq(2, 30), scen1_m1_posterior[2:30, 4], col=coul[1], lty=3)
lines(seq(2, 30), scen1_m1_posterior[2:30, 5], col=coul[2], lty=3)
lines(seq(2, 30), scen1_m1_posterior[2:30, 6], col=coul[3], lty=3)
lines(seq(2, 30), scen1_m1_posterior[2:30, 7], col=coul[1], lty=5)
lines(seq(2, 30), scen1_m1_posterior[2:30, 8], col=coul[2], lty=5)
lines(seq(2, 30), scen1_m1_posterior[2:30, 9], col=coul[3], lty=5)
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab="Predicted probability \n candidate wins",
     xaxt = "n", yaxt= "n", cex.lab=0.7)
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen1_m2_posterior[2:30, 1], col=coul[1])
lines(seq(2, 30), scen1_m2_posterior[2:30, 2], col=coul[2])
lines(seq(2, 30), scen1_m2_posterior[2:30, 3], col=coul[3])
lines(seq(2, 30), scen1_m2_posterior[2:30, 4], col=coul[1], lty=3)
lines(seq(2, 30), scen1_m2_posterior[2:30, 5], col=coul[2], lty=3)
lines(seq(2, 30), scen1_m2_posterior[2:30, 6], col=coul[3], lty=3)
lines(seq(2, 30), scen1_m2_posterior[2:30, 7], col=coul[1], lty=5)
lines(seq(2, 30), scen1_m2_posterior[2:30, 8], col=coul[2], lty=5)
lines(seq(2, 30), scen1_m2_posterior[2:30, 9], col=coul[3], lty=5)
# alpha estimates
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(15, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\alpha}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(15, 30), scen1_m3_A[15:30, 1], col=coul[4])
lines(seq(15, 30), scen1_m3_A[15:30, 2], col=coul[4], lty=3)
lines(seq(15, 30), scen1_m3_A[15:30, 3], col=coul[4], lty=5)
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\alpha}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen1_m4_A[2:30, 1], col=coul[4])
lines(seq(2, 30), scen1_m4_A[2:30, 2], col=coul[4], lty=3)
lines(seq(2, 30), scen1_m4_A[2:30, 3], col=coul[4], lty=5)
# beta estimates for M(3, 3), M(3, 4), M(3, 5)
plot(1, type="n", xlim=c(15, 30), ylim=c(0, 0.3),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\beta}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(15, 30), scen1_m3_beta[15:30, 1], col=coul[4])
lines(seq(15, 30), scen1_m3_beta[15:30, 2], col=coul[4], lty=3)
lines(seq(15, 30), scen1_m3_beta[15:30, 3], col=coul[4], lty=5)
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 0.3),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\beta}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen1_m4_beta[2:30, 1], col=coul[4])
lines(seq(2, 30), scen1_m4_beta[2:30, 2], col=coul[4], lty=3)
lines(seq(2, 30), scen1_m4_beta[2:30, 3], col=coul[4], lty=5)
plot(1, type="n", xlim=c(15, 30), ylim=c(0, 0.3),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\beta}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen1_m5_beta[2:30, 1], col=coul[4])
lines(seq(2, 30), scen1_m5_beta[2:30, 2], col=coul[4], lty=3)
lines(seq(2, 30), scen1_m5_beta[2:30, 3], col=coul[4], lty=5)
# estimating tau for model 4
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 30),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\tau}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen1_m4_tau[2:30, 1], col=coul[4])
lines(seq(2, 30), scen1_m4_tau[2:30, 2], col=coul[4], lty=3)
lines(seq(2, 30), scen1_m4_tau[2:30, 3], col=coul[4], lty=5)
# estimating rho
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(15, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\rho}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(15, 30), scen1_m3_rho[15:30, 1], col=coul[1])
lines(seq(15, 30), scen1_m3_rho[15:30, 2], col=coul[2])
lines(seq(15, 30), scen1_m3_rho[15:30, 3], col=coul[3])
lines(seq(15, 30), scen1_m3_rho[15:30, 4], col=coul[1], lty=3)
lines(seq(15, 30), scen1_m3_rho[15:30, 5], col=coul[2], lty=3)
lines(seq(15, 30), scen1_m3_rho[15:30, 6], col=coul[3], lty=3)
lines(seq(15, 30), scen1_m3_rho[15:30, 7], col=coul[1], lty=5)
lines(seq(15, 30), scen1_m3_rho[15:30, 8], col=coul[2], lty=5)
lines(seq(15, 30), scen1_m3_rho[15:30, 9], col=coul[3], lty=5)
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\rho}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen1_m4_rho[2:30, 1], col=coul[1])
lines(seq(2, 30), scen1_m4_rho[2:30, 2], col=coul[2])
lines(seq(2, 30), scen1_m4_rho[2:30, 3], col=coul[3])
lines(seq(2, 30), scen1_m4_rho[2:30, 4], col=coul[1], lty=3)
lines(seq(2, 30), scen1_m4_rho[2:30, 5], col=coul[2], lty=3)
lines(seq(2, 30), scen1_m4_rho[2:30, 6], col=coul[3], lty=3)
lines(seq(2, 30), scen1_m4_rho[2:30, 7], col=coul[1], lty=5)
lines(seq(2, 30), scen1_m4_rho[2:30, 8], col=coul[2], lty=5)
lines(seq(2, 30), scen1_m4_rho[2:30, 9], col=coul[3], lty=5)
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(15, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\rho}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(15, 30), scen1_m5_rho[15:30, 1], col=coul[1])
lines(seq(15, 30), scen1_m5_rho[15:30, 2], col=coul[2])
lines(seq(15, 30), scen1_m5_rho[15:30, 3], col=coul[3])
lines(seq(15, 30), scen1_m5_rho[15:30, 4], col=coul[1], lty=3)
lines(seq(15, 30), scen1_m5_rho[15:30, 5], col=coul[2], lty=3)
lines(seq(15, 30), scen1_m5_rho[15:30, 6], col=coul[3], lty=3)
lines(seq(15, 30), scen1_m5_rho[15:30, 7], col=coul[1], lty=5)
lines(seq(15, 30), scen1_m5_rho[15:30, 8], col=coul[2], lty=5)
lines(seq(15, 30), scen1_m5_rho[15:30, 9], col=coul[3], lty=5)
# plotting figures 4.5 and 4.6
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(0, 30), ylim=c(0, 20),
     bty="n", xlab="Time", ylab=expression(xi[t]),
     xaxt = "n", yaxt= "n", cex.lab=0.7)
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(0, 30), c(0, cumsum(scen2[, 1])), col=coul[1])
lines(seq(0, 30), c(0, cumsum(scen2[, 2])), col=coul[2])
lines(seq(0, 30), c(0, cumsum(scen2[, 3])), col=coul[3])
lines(seq(0, 30), c(0, cumsum(scen2[, 4])), col=coul[1], lty=3)
lines(seq(0, 30), c(0, cumsum(scen2[, 5])), col=coul[2], lty=3)
lines(seq(0, 30), c(0, cumsum(scen2[, 6])), col=coul[3], lty=3)
lines(seq(0, 30), c(0, cumsum(scen2[, 7])), col=coul[1], lty=5)
lines(seq(0, 30), c(0, cumsum(scen2[, 8])), col=coul[2], lty=5)
lines(seq(0, 30), c(0, cumsum(scen2[, 9])), col=coul[3], lty=5)
# posterior of each model
#load("dirichlet_data.Rdata")
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab="Predicted probability \n candidate wins",
     xaxt = "n", yaxt= "n", cex.lab=0.7)
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen2_m1_posterior[2:30, 1], col=coul[1])
lines(seq(2, 30), scen2_m1_posterior[2:30, 2], col=coul[2])
lines(seq(2, 30), scen2_m1_posterior[2:30, 3], col=coul[3])
lines(seq(2, 30), scen2_m1_posterior[2:30, 4], col=coul[1], lty=3)
lines(seq(2, 30), scen2_m1_posterior[2:30, 5], col=coul[2], lty=3)
lines(seq(2, 30), scen2_m1_posterior[2:30, 6], col=coul[3], lty=3)
lines(seq(2, 30), scen2_m1_posterior[2:30, 7], col=coul[1], lty=5)
lines(seq(2, 30), scen2_m1_posterior[2:30, 8], col=coul[2], lty=5)
lines(seq(2, 30), scen2_m1_posterior[2:30, 9], col=coul[3], lty=5)
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab="Predicted probability \n candidate wins",
     xaxt = "n", yaxt= "n", cex.lab=0.7)
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen2_m2_posterior[2:30, 1], col=coul[1])
lines(seq(2, 30), scen2_m2_posterior[2:30, 2], col=coul[2])
lines(seq(2, 30), scen2_m2_posterior[2:30, 3], col=coul[3])
lines(seq(2, 30), scen2_m2_posterior[2:30, 4], col=coul[1], lty=3)
lines(seq(2, 30), scen2_m2_posterior[2:30, 5], col=coul[2], lty=3)
lines(seq(2, 30), scen2_m2_posterior[2:30, 6], col=coul[3], lty=3)
lines(seq(2, 30), scen2_m2_posterior[2:30, 7], col=coul[1], lty=5)
lines(seq(2, 30), scen2_m2_posterior[2:30, 8], col=coul[2], lty=5)
lines(seq(2, 30), scen2_m2_posterior[2:30, 9], col=coul[3], lty=5)
# alpha estimates
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(15, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\alpha}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(15, 30), scen2_m3_A[15:30, 1], col=coul[4])
lines(seq(15, 30), scen2_m3_A[15:30, 2], col=coul[4], lty=3)
lines(seq(15, 30), scen2_m3_A[15:30, 3], col=coul[4], lty=5)
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\alpha}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen2_m4_A[2:30, 1], col=coul[4])
lines(seq(2, 30), scen2_m4_A[2:30, 2], col=coul[4], lty=3)
lines(seq(2, 30), scen2_m4_A[2:30, 3], col=coul[4], lty=5)
# beta estimates for m3, m4, m5
plot(1, type="n", xlim=c(15, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\beta}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(15, 30), scen2_m3_beta[15:30, 1], col=coul[4])
lines(seq(15, 30), scen2_m3_beta[15:30, 2], col=coul[4], lty=3)
lines(seq(15, 30), scen2_m3_beta[15:30, 3], col=coul[4], lty=5)
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\beta}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen2_m4_beta[2:30, 1], col=coul[4])
lines(seq(2, 30), scen2_m4_beta[2:30, 2], col=coul[4], lty=3)
lines(seq(2, 30), scen2_m4_beta[2:30, 3], col=coul[4], lty=5)
plot(1, type="n", xlim=c(15, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\beta}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen2_m5_beta[2:30, 1], col=coul[4])
lines(seq(2, 30), scen2_m5_beta[2:30, 2], col=coul[4], lty=3)
lines(seq(2, 30), scen2_m5_beta[2:30, 3], col=coul[4], lty=5)
# estimating tau for model 4
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 30),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\tau}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen2_m4_tau[2:30, 1], col=coul[4])
lines(seq(2, 30), scen2_m4_tau[2:30, 2], col=coul[4], lty=3)
lines(seq(2, 30), scen2_m4_tau[2:30, 3], col=coul[4], lty=5)
# estimating rho
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(15, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\rho}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(15, 30), scen2_m3_rho[15:30, 1], col=coul[1])
lines(seq(15, 30), scen2_m3_rho[15:30, 2], col=coul[2])
lines(seq(15, 30), scen2_m3_rho[15:30, 3], col=coul[3])
lines(seq(15, 30), scen2_m3_rho[15:30, 4], col=coul[1], lty=3)
lines(seq(15, 30), scen2_m3_rho[15:30, 5], col=coul[2], lty=3)
lines(seq(15, 30), scen2_m3_rho[15:30, 6], col=coul[3], lty=3)
lines(seq(15, 30), scen2_m3_rho[15:30, 7], col=coul[1], lty=5)
lines(seq(15, 30), scen2_m3_rho[15:30, 8], col=coul[2], lty=5)
lines(seq(15, 30), scen2_m3_rho[15:30, 9], col=coul[3], lty=5)
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\rho}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen2_m4_rho[2:30, 1], col=coul[1])
lines(seq(2, 30), scen2_m4_rho[2:30, 2], col=coul[2])
lines(seq(2, 30), scen2_m4_rho[2:30, 3], col=coul[3])
lines(seq(2, 30), scen2_m4_rho[2:30, 4], col=coul[1], lty=3)
lines(seq(2, 30), scen2_m4_rho[2:30, 5], col=coul[2], lty=3)
lines(seq(2, 30), scen2_m4_rho[2:30, 6], col=coul[3], lty=3)
lines(seq(2, 30), scen2_m4_rho[2:30, 7], col=coul[1], lty=5)
lines(seq(2, 30), scen2_m4_rho[2:30, 8], col=coul[2], lty=5)
lines(seq(2, 30), scen2_m4_rho[2:30, 9], col=coul[3], lty=5)
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(15, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\rho}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(15, 30), scen2_m5_rho[15:30, 1], col=coul[1])
lines(seq(15, 30), scen2_m5_rho[15:30, 2], col=coul[2])
lines(seq(15, 30), scen2_m5_rho[15:30, 3], col=coul[3])
lines(seq(15, 30), scen2_m5_rho[15:30, 4], col=coul[1], lty=3)
lines(seq(15, 30), scen2_m5_rho[15:30, 5], col=coul[2], lty=3)
lines(seq(15, 30), scen2_m5_rho[15:30, 6], col=coul[3], lty=3)
lines(seq(15, 30), scen2_m5_rho[15:30, 7], col=coul[1], lty=5)
lines(seq(15, 30), scen2_m5_rho[15:30, 8], col=coul[2], lty=5)
lines(seq(15, 30), scen2_m5_rho[15:30, 9], col=coul[3], lty=5)
# plotting figures 4.7 and 4.8
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(0, 30), ylim=c(0, 20),
     bty="n", xlab="Time", ylab=expression(xi[t]),
     xaxt = "n", yaxt= "n", cex.lab=0.7)
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(0, 30), c(0, cumsum(scen3[, 1])), col=coul[1])
lines(seq(0, 30), c(0, cumsum(scen3[, 2])), col=coul[2])
lines(seq(0, 30), c(0, cumsum(scen3[, 3])), col=coul[3])
lines(seq(0, 30), c(0, cumsum(scen3[, 4])), col=coul[1], lty=3)
lines(seq(0, 30), c(0, cumsum(scen3[, 5])), col=coul[2], lty=3)
lines(seq(0, 30), c(0, cumsum(scen3[, 6])), col=coul[3], lty=3)
lines(seq(0, 30), c(0, cumsum(scen3[, 7])), col=coul[1], lty=5)
lines(seq(0, 30), c(0, cumsum(scen3[, 8])), col=coul[2], lty=5)
lines(seq(0, 30), c(0, cumsum(scen3[, 9])), col=coul[3], lty=5)
# posterior of each model
#load("dirichlet_data.Rdata")
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab="Predicted probability \n candidate wins",
     xaxt = "n", yaxt= "n", cex.lab=0.7)
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen3_m1_posterior[2:30, 1], col=coul[1])
lines(seq(2, 30), scen3_m1_posterior[2:30, 2], col=coul[2])
lines(seq(2, 30), scen3_m1_posterior[2:30, 3], col=coul[3])
lines(seq(2, 30), scen3_m1_posterior[2:30, 4], col=coul[1], lty=3)
lines(seq(2, 30), scen3_m1_posterior[2:30, 5], col=coul[2], lty=3)
lines(seq(2, 30), scen3_m1_posterior[2:30, 6], col=coul[3], lty=3)
lines(seq(2, 30), scen3_m1_posterior[2:30, 7], col=coul[1], lty=5)
lines(seq(2, 30), scen3_m1_posterior[2:30, 8], col=coul[2], lty=5)
lines(seq(2, 30), scen3_m1_posterior[2:30, 9], col=coul[3], lty=5)
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab="Predicted probability \n candidate wins",
     xaxt = "n", yaxt= "n", cex.lab=0.7)
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen3_m2_posterior[2:30, 1], col=coul[1])
lines(seq(2, 30), scen3_m2_posterior[2:30, 2], col=coul[2])
lines(seq(2, 30), scen3_m2_posterior[2:30, 3], col=coul[3])
lines(seq(2, 30), scen3_m2_posterior[2:30, 4], col=coul[1], lty=3)
lines(seq(2, 30), scen3_m2_posterior[2:30, 5], col=coul[2], lty=3)
lines(seq(2, 30), scen3_m2_posterior[2:30, 6], col=coul[3], lty=3)
lines(seq(2, 30), scen3_m2_posterior[2:30, 7], col=coul[1], lty=5)
lines(seq(2, 30), scen3_m2_posterior[2:30, 8], col=coul[2], lty=5)
lines(seq(2, 30), scen3_m2_posterior[2:30, 9], col=coul[3], lty=5)
# alpha estimates
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(15, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\alpha}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(15, 30), scen3_m3_A[15:30, 1], col=coul[4])
lines(seq(15, 30), scen3_m3_A[15:30, 2], col=coul[4], lty=3)
lines(seq(15, 30), scen3_m3_A[15:30, 3], col=coul[4], lty=5)
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\alpha}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen3_m4_A[2:30, 1], col=coul[4])
lines(seq(2, 30), scen3_m4_A[2:30, 2], col=coul[4], lty=3)
lines(seq(2, 30), scen3_m4_A[2:30, 3], col=coul[4], lty=5)
# beta estimates for m3, m4, m5
plot(1, type="n", xlim=c(15, 30), ylim=c(0, 2),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\beta}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(15, 30), scen3_m3_beta[15:30, 1], col=coul[4])
lines(seq(15, 30), scen3_m3_beta[15:30, 2], col=coul[4], lty=3)
lines(seq(15, 30), scen3_m3_beta[15:30, 3], col=coul[4], lty=5)
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 2),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\beta}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen3_m4_beta[2:30, 1], col=coul[4])
lines(seq(2, 30), scen3_m4_beta[2:30, 2], col=coul[4], lty=3)
lines(seq(2, 30), scen3_m4_beta[2:30, 3], col=coul[4], lty=5)
plot(1, type="n", xlim=c(15, 30), ylim=c(0, 2),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\beta}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen3_m5_beta[2:30, 1], col=coul[4])
lines(seq(2, 30), scen3_m5_beta[2:30, 2], col=coul[4], lty=3)
lines(seq(2, 30), scen3_m5_beta[2:30, 3], col=coul[4], lty=5)
# estimating tau for model 4
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 30),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\tau}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen3_m4_tau[2:30, 1], col=coul[4])
lines(seq(2, 30), scen3_m4_tau[2:30, 2], col=coul[4], lty=3)
lines(seq(2, 30), scen3_m4_tau[2:30, 3], col=coul[4], lty=5)
# estimating rho
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(15, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\rho}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(15, 30), scen3_m3_rho[15:30, 1], col=coul[1])
lines(seq(15, 30), scen3_m3_rho[15:30, 2], col=coul[2])
lines(seq(15, 30), scen3_m3_rho[15:30, 3], col=coul[3])
lines(seq(15, 30), scen3_m3_rho[15:30, 4], col=coul[1], lty=3)
lines(seq(15, 30), scen3_m3_rho[15:30, 5], col=coul[2], lty=3)
lines(seq(15, 30), scen3_m3_rho[15:30, 6], col=coul[3], lty=3)
lines(seq(15, 30), scen3_m3_rho[15:30, 7], col=coul[1], lty=5)
lines(seq(15, 30), scen3_m3_rho[15:30, 8], col=coul[2], lty=5)
lines(seq(15, 30), scen3_m3_rho[15:30, 9], col=coul[3], lty=5)
plot(1, type="n", xlim=c(2, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\rho}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(2, 30), scen3_m4_rho[2:30, 1], col=coul[1])
lines(seq(2, 30), scen3_m4_rho[2:30, 2], col=coul[2])
lines(seq(2, 30), scen3_m4_rho[2:30, 3], col=coul[3])
lines(seq(2, 30), scen3_m4_rho[2:30, 4], col=coul[1], lty=3)
lines(seq(2, 30), scen3_m4_rho[2:30, 5], col=coul[2], lty=3)
lines(seq(2, 30), scen3_m4_rho[2:30, 6], col=coul[3], lty=3)
lines(seq(2, 30), scen3_m4_rho[2:30, 7], col=coul[1], lty=5)
lines(seq(2, 30), scen3_m4_rho[2:30, 8], col=coul[2], lty=5)
lines(seq(2, 30), scen3_m4_rho[2:30, 9], col=coul[3], lty=5)
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
plot(1, type="n", xlim=c(15, 30), ylim=c(0, 1),
     bty="n", xlab="Time", ylab=TeX('$\\hat{\\rho}$'),
     xaxt = "n", yaxt= "n")
axis(2, cex.axis=0.65)
axis(1, cex.axis=0.8)
lines(seq(15, 30), scen3_m5_rho[15:30, 1], col=coul[1])
lines(seq(15, 30), scen3_m5_rho[15:30, 2], col=coul[2])
lines(seq(15, 30), scen3_m5_rho[15:30, 3], col=coul[3])
lines(seq(15, 30), scen3_m5_rho[15:30, 4], col=coul[1], lty=3)
lines(seq(15, 30), scen3_m5_rho[15:30, 5], col=coul[2], lty=3)
lines(seq(15, 30), scen3_m5_rho[15:30, 6], col=coul[3], lty=3)
lines(seq(15, 30), scen3_m5_rho[15:30, 7], col=coul[1], lty=5)
lines(seq(15, 30), scen3_m5_rho[15:30, 8], col=coul[2], lty=5)
lines(seq(15, 30), scen3_m5_rho[15:30, 9], col=coul[3], lty=5)