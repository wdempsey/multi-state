## Build pretty trace plots
par(mfrow = c(2,2), mar = c(4,4,1,1) + 0.1)
plot(posterior.results[1,], ylab = "Values", xlab = "Iteration", axes = FALSE); axis(side = 1); axis(side = 2)
plot(posterior.results[2,], ylab = "Values", xlab = "Iteration", axes = FALSE); axis(side = 1); axis(side = 2)
plot(posterior.results[3,], ylab = "Values", xlab = "Iteration", axes = FALSE); axis(side = 1); axis(side = 2)
plot(posterior.results[4,], ylab = "Values", xlab = "Iteration", axes = FALSE); axis(side = 1); axis(side = 2)

## Trace plots -> Toss first 100
min.iter = 100; max.iter = ncol(posterior.results)

## Build 5% - 95% quantiles
## And plot histograms
quantile(posterior.results[1,min.iter:max.iter], probs = c(0.05, 0.50, 0.95)); mean(posterior.results[1,min.iter:max.iter])
quantile(posterior.results[2,min.iter:max.iter], probs = c(0.05, 0.50, 0.95)); mean(posterior.results[2,min.iter:max.iter])
quantile(posterior.results[3,min.iter:max.iter], probs = c(0.05, 0.50, 0.95)); mean(posterior.results[3,min.iter:max.iter])
quantile(posterior.results[4,min.iter:max.iter], probs = c(0.05, 0.50, 0.95)); mean(posterior.results[4,min.iter:max.iter])

plot(density(posterior.results[1,min.iter:max.iter]), axes = FALSE, main = "", xlab = ""); axis(side = 1); axis(side = 2)
plot(density(posterior.results[2,min.iter:max.iter]), axes = FALSE, main = "", xlab = ""); axis(side = 1); axis(side = 2)
plot(density(posterior.results[3,min.iter:max.iter]), axes = FALSE, main = "", xlab = ""); axis(side = 1); axis(side = 2)
plot(density(posterior.results[4,min.iter:max.iter]), axes = FALSE, main = "", xlab = ""); axis(side = 1); axis(side = 2)

## Survival curves from baseline
require(Matrix); Survival.list = list()

if(!file.exists("survival-fits.RDS")) {
  for (iter in min.iter:max.iter) {
    current.global = person.global.data.list[[iter]]
    current.params = posterior.results[,iter]
    Survival.list[[iter]] = Survival.fit(current.global, current.params, grid.length = 0.1)
  }
  saveRDS(Survival.list, file = "survival-fits.RDS")
} else {
  Survival.list = loadRDS("survival-fits.RDS")
}

### Plot truth
simdata = readRDS("simdata.rds"); simdata = data.frame(simdata)
names(simdata) = c("time", "num.healthy", "num.ill", "num.dead", "switch1", "switch2", "failure.time")
truth = c(1/2, 1/5, 0.65, 1.71)
true.Survival = Survival.fit(simdata, truth, grid.length = 0.1)

survival.grid = seq(0,10,0.1)
healthy.survival.plot.matrix = ill.survival.plot.matrix = matrix(nrow = length(survival.grid), ncol = max.iter - min.iter + 1)

for (iter in min.iter:max.iter) {
  temp.surv = Survival.list[[iter]]
  
  for (i in 1:length(survival.grid)) {
    temp.abs.diff = abs(temp.surv[,1] - survival.grid[i])
    temp.entry.to.use = min(which( temp.abs.diff == min(temp.abs.diff) ))
    healthy.survival.plot.matrix[i, iter - min.iter + 1] = temp.surv[temp.entry.to.use,2]
    ill.survival.plot.matrix[i, iter - min.iter + 1] = temp.surv[temp.entry.to.use,3]
  }
}

quantiles.survival.ill = quantiles.survival.healthy = matrix(nrow = length(survival.grid), ncol = 3)

for (i in 1:nrow(quantiles.survival.ill)) {
  quantiles.survival.ill[i,] = as.vector(quantile(ill.survival.plot.matrix[i,], prob = c(0.05, 0.50, 0.95)))
  quantiles.survival.healthy[i,] = as.vector(quantile(healthy.survival.plot.matrix[i,], prob = c(0.05, 0.50, 0.95)))
}

obs.Survival = true.Survival[,1] < 10

par(mfrow = c(2,1), mar = c(4,4,1,1) + 0.1)
plot(survival.grid, 1 - quantiles.survival.healthy[,2], type = "l", axes = FALSE, ylab = "Survival", xlab = "Time")
axis(side = 1); axis(side = 2)
lines(survival.grid, 1 - quantiles.survival.healthy[,1], lty = 2)
lines(survival.grid, 1 - quantiles.survival.healthy[,3], lty = 2)
lines(true.Survival[obs.Survival,1], 1 - true.Survival[obs.Survival,2], col = "red")

plot(survival.grid, 1 - quantiles.survival.ill[,2], type = "l", axes = FALSE, ylab = "Survival", xlab = "Time")
axis(side = 1); axis(side = 2)
lines(survival.grid, 1 - quantiles.survival.ill[,1], lty = 2)
lines(survival.grid, 1 - quantiles.survival.ill[,3], lty = 2)
lines(true.Survival[obs.Survival,1], 1 - true.Survival[obs.Survival,3], col = "red")

