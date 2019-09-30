## Source code to analyze cav example
## Under multi-censored observations
## Made every year while participant is alive

# Source Functions and observation data
source("cav-functions.R")

library(msm)
person.data = data.frame(cbind(cav$PTNUM, cav$years, cav$state))
names(person.data) = c("user", "time", "state")
unique.users = unique(person.data$user)

## Remove agreement of failure and state obs times
## Add a small tolerance term epsilon = 10^-6
epsilon = 10^-6
death.times = person.data$time[person.data$state == 4]
obs.times = person.data$time[person.data$state != 4]
bad.death.times = unique(death.times[(is.element(death.times, obs.times))])
for(times in bad.death.times) {
  bad.obs = which(person.data$time == times & person.data$state == 4)
  person.data$time[bad.obs] = person.data$time[bad.obs] + epsilon
}

## Figure out who died and who is censored times
death.data = c()
for(i in 1:length(unique.users)) {
  person.data$user[person.data$user == unique.users[i]] = i
  temp.data = person.data[person.data$user == i,]
  final.obs.temp = temp.data[which(temp.data$time == max(temp.data$time)),]
  death.data = rbind(death.data, c(final.obs.temp$user, final.obs.temp$time, final.obs.temp$state == 4))
}

death.data = data.frame(death.data)
names(death.data) = c("user", "time", "status")
censoring.data = death.data

global.data = convert.to.global(person.data, death.data)

rho = 10 # Assume rho is fixed throughout

## Generate the survival curve under non-parametric
## estimator where observations = transition times.
result.distribution.nocav = nonparametric_estimator(c(1,0,0,0), global.data)
result.distribution.severe = nonparametric_estimator(c(0,0,1,0), global.data)


### FULLY PARAMETRIC MODEL
library(expm)

llik.parametric = logLik.iid(person.data)
parameters = rep(0.1, 7)
lower.bound = rep(0.00001, 9)
upper.bound = rep(1,9)
# mle.parameter = optim(par = parameters, fn = llik.parametric, lower = lower.bound, upper = upper.bound)

mle.parameters = c(0.04116865, 0.01627384, 0.12346032,
                  0.35751591, 0.04557379, 0.13354901, 0.30483947)

## CONSTRUCT Q MATRIX
Q = matrix(0, nrow = 4, ncol = 4)
Q[1,c(2,4)] = mle.parameters[1:2]; Q[1,1] = - sum(Q[1,2:4])
Q[2,c(1,3:4)] = mle.parameters[3:5]; Q[2,2] = - sum(Q[2,c(1,3:4)])
Q[3,c(2,4)] = mle.parameters[6:7]; Q[3,3] = - sum(Q[3, c(1:2,4)])
Q[4,1:4] = 0

survival.list.nocav = survival.list.severe = vector(length = length(global.data$time))

for(iter in 1:length(global.data$time)) {
  current.time = global.data$time[iter]
  P = expm(Q*current.time)
  survival.list.nocav[iter] = sum(P[1,1:3])
  survival.list.severe[iter] = sum(P[3,1:3])
}


## My results
Survival.list = readRDS("C:/Users/wdem/Dropbox/Dissertation Work/Markov State Process/multi-state/cav/survival-fits.RDS")
max.iter = 1000; min.iter = 50
survival.grid = seq(0,20,0.1)
nocav.survival.plot.matrix = mild.survival.plot.matrix = 
  severe.survival.plot.matrix = matrix(nrow = length(survival.grid), ncol = max.iter - min.iter + 1)


for (iter in min.iter:max.iter) {
  temp.surv = Survival.list[[iter]]
  
  for (i in 1:length(survival.grid)) {
    temp.abs.diff = abs(temp.surv[,1] - survival.grid[i])
    temp.entry.to.use = min(which( temp.abs.diff == min(temp.abs.diff) ))
    nocav.survival.plot.matrix[i, iter - min.iter + 1] = temp.surv[temp.entry.to.use,2]
    mild.survival.plot.matrix[i, iter - min.iter + 1] = temp.surv[temp.entry.to.use,3]
    severe.survival.plot.matrix[i, iter - min.iter + 1] = temp.surv[temp.entry.to.use,4]
  }
}

quantiles.survival.mild = quantiles.survival.nocav = 
  quantiles.survival.severe = matrix(nrow = length(survival.grid), ncol = 3)

for (i in 1:nrow(quantiles.survival.nocav)) {
  quantiles.survival.nocav[i,] = as.vector(quantile(nocav.survival.plot.matrix[i,], prob = c(0.05, 0.50, 0.95)))
  quantiles.survival.mild[i,] = as.vector(quantile(mild.survival.plot.matrix[i,], prob = c(0.05, 0.50, 0.95)))
  quantiles.survival.severe[i,] = as.vector(quantile(severe.survival.plot.matrix[i,], prob = c(0.05, 0.50, 0.95)))
}


### FULL PLOT 

png(filename = "C:/Users/wdem/Dropbox/Dissertation Work/Markov State Process/multi-state/cav/cav-comparison-plot.png", 
    width = 480, height = 480, units = "px", pointsize = 12, bg = "white")

par(mar = c(5,4,2,1)+ 0.1)
plot(global.data$time, 1 - result.distribution.nocav[,4], 
     xlim = c(0,20), ylim = c(0,1), type = "l", bty = 'n',
     xlab = "Time in study", ylab = "Survival function")
lines(global.data$time, 1 - result.distribution.severe[,4], lty = 2)

lines(global.data$time, survival.list.nocav, col = "red")
lines(global.data$time, survival.list.severe, col = "red", lty = 2)

lines(survival.grid, 1 - quantiles.survival.nocav[,2], col = "blue")
lines(survival.grid, 1 - quantiles.survival.severe[,2], col = "blue", lty = 2)
dev.off()