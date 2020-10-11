library(RColorBrewer)
library(sciencesPo)
par(mgp=c(2,0.5,0), mar=c(5,4,4,2)+0.1)
coul <- brewer.pal(8, "Dark2") 

plot(time, 
     BrodySim1, 
     type = "l",
     col=coul[1],
     bty="n",
     xlab="Time",
     ylab=expression(xi[t]),
     ylim=c(-1, 20))
lines(time, BrodySim2, col=coul[2])
lines(time, BrodySim3, col=coul[3])
legend(0.5, 19.5, legend=c("Simulation 1", "Simulation 2", "Simulation 3"),
       col=c(coul[1], coul[2], coul[3]), pch=19, cex=0.6, bty = "n")

plot(time, 
     BrodyPost1, 
     type = "l",
     col=coul[1],
     bty="n",
     xlab="Time",
     ylab="P(X=1)",
     ylim=c(0, 1))
lines(time, BrodyPost2, col=coul[2])
lines(time, BrodyPost3, col=coul[3])
legend(0.5, 0.98, legend=c("Simulation 1", "Simulation 2", "Simulation 3"),
       col=c(coul[1], coul[2], coul[3]), pch=19, cex=0.6, bty = "n")

plot(time, 
     obs1, 
     type = "l",
     col=coul[1],
     bty="n",
     xlab="Time",
     ylab=expression(xi[t]))
lines(time, obs2, col=coul[2])
lines(time, obs3, col=coul[3])
lines(time, BrodySim1, lty="dashed", col=coul[1])
lines(time, BrodySim2, lty="dashed", col=coul[2])
lines(time, BrodySim3, lty="dashed", col=coul[3])
legend(0.5, 24, legend=c("Simulation 1", "Simulation 2", "Simulation 3"),
       col=c(coul[1], coul[2], coul[3]), pch=19, cex=0.6, bty = "n")

plot(time, 
     BrodyAdjPost1, 
     type = "l",
     col=coul[1],
     bty="n",
     xlab="Time",
     ylab="P(X=1)",
     ylim=c(0, 1))
lines(time, BrodyAdjPost2, col=coul[2])
lines(time, BrodyAdjPost3, col=coul[3])
lines(time, BrodyPost1, col=coul[1], lty="dashed")
lines(time, BrodyPost2, col=coul[2], lty ="dashed")
lines(time, BrodyPost3, col=coul[3], lty="dashed")
legend(0.5, 0.98, legend=c("Simulation 1", "Simulation 2", "Simulation 3"),
       col=c(coul[1], coul[2], coul[3]), pch=19, cex=0.6, bty = "n")




# idea
library(GISTools)
time <- seq(0, 30)
expected_xi <- seq(0, 30) * 0.2
set.seed(7)
coul.t = add.alpha(coul,0.25)

plot(time, 
     expected_xi, 
     type = "l",
     col=coul[1],
     bty="n",
     xlab="Time",
     ylab=expression(xi[t]),
     ylim=c(-10, 20),
     xaxt = "n",
     yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:50){
observed_sequence <- 0.2 * time + c(0, cumsum(rnorm(n=30)))
lines(time, observed_sequence, col=coul.t[1])
}

expected_obs <- posterior(expected_xi, time, 1, 0.2, 0.5)

set.seed(7)

plot(time, 
     expected_obs, 
     type = "l",
     col=coul[2],
     bty="n",
     xlab="Time",
     ylab="P(X=1)",
     ylim=c(0, 1),
     xaxt = "n",
     yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:50){
        observed_sequence <- 0.2 * time + c(0, cumsum(rnorm(n=30)))
        post_obs <- posterior(observed_sequence, time, 1, 0.2, 0.5)
        lines(time, post_obs, col=coul.t[2])
}

fake_effect_1 <- Brody_fake_news(time, 1, 0.05, 1)
set.seed(7)
expected_xi_1 <- expected_xi + fake_effect_1
plot(time, 
     expected_xi_1, 
     type = "l",
     col=coul[1],
     bty="n",
     xlab="Time",
     ylab=expression(xi[t]),
     ylim=c(-5, 27),
     xaxt = "n",
     yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:50){
        observed_sequence <- 0.2 * time + c(0, cumsum(rnorm(n=30))) + fake_effect_1
        lines(time, observed_sequence, col=coul.t[1])
}

expected_obs_1 <- posterior(expected_xi_1, time, 1, 0.2, 0.5)

set.seed(7)

plot(time, 
     expected_obs_1, 
     type = "l",
     col=coul[2],
     bty="n",
     xlab="Time",
     ylab="P(X=1)",
     ylim=c(0, 1),
     xaxt = "n",
     yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:50){
        observed_sequence <- 0.2 * time + c(0, cumsum(rnorm(n=30))) + fake_effect_1
        post_obs <- posterior(observed_sequence, time, 1, 0.2, 0.5)
        lines(time, post_obs, col=coul.t[2])
}

fake_effect_2 <- Brody_fake_news(time, 10, 1, 25)
set.seed(7)
expected_xi_2 <- expected_xi + fake_effect_2
plot(time, 
     expected_xi_2, 
     type = "l",
     col=coul[1],
     bty="n",
     xlab="Time",
     ylab=expression(xi[t]),
     ylim=c(-10, 27),
     xaxt = "n",
     yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:50){
        observed_sequence <- 0.2 * time + c(0, cumsum(rnorm(n=30))) + fake_effect_2
        lines(time, observed_sequence, col=coul.t[1])
}

expected_obs_2 <- posterior(expected_xi_2, time, 1, 0.2, 0.5)

set.seed(7)

plot(time, 
     expected_obs_2, 
     type = "l",
     col=coul[2],
     bty="n",
     xlab="Time",
     ylab="P(X=1)",
     ylim=c(0, 1),
     xaxt = "n",
     yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:50){
        observed_sequence <- 0.2 * time + c(0, cumsum(rnorm(n=30))) + fake_effect_2
        post_obs <- posterior(observed_sequence, time, 1, 0.2, 0.5)
        lines(time, post_obs, col=coul.t[2])
}

fake_effect_3 <- Brody_fake_news(time, 20, 0.5, 15)
set.seed(7)
expected_xi_3 <- expected_xi + fake_effect_3
plot(time, 
     expected_xi_3, 
     type = "l",
     col=coul[1],
     bty="n",
     xlab="Time",
     ylab=expression(xi[t]),
     ylim=c(-10, 30),
     xaxt = "n",
     yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:50){
        observed_sequence <- 0.2 * time + c(0, cumsum(rnorm(n=30))) + fake_effect_3
        lines(time, observed_sequence, col=coul.t[1])
}

expected_obs_3 <- posterior(expected_xi_3, time, 1, 0.2, 0.5)
set.seed(7)
plot(time, 
     expected_obs_3, 
     type = "l",
     col=coul[2],
     bty="n",
     xlab="Time",
     ylab="P(X=1)",
     ylim=c(0, 1),
     xaxt = "n",
     yaxt= "n")
axis(2, cex.axis=0.8)
axis(1, cex.axis=0.8)
for(i in 1:50){
        observed_sequence <- 0.2 * time + c(0, cumsum(rnorm(n=30))) + fake_effect_3
        post_obs <- posterior(observed_sequence, time, 1, 0.2, 0.5)
        lines(time, post_obs, col=coul.t[2])
}
