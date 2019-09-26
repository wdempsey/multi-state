posterior.results = readRDS(file = "posterior-params-results.RDS")
person.global.data.list = readRDS(file = "posterior-transdata-results.RDS")
acceptance.rate.results = readRDS(file = "acceptance-results.RDS")
source("functions.R")

special.posterior.results = readRDS(file = "special-case/special-posterior-params-results.RDS")
special.person.global.data.list = readRDS(file = "special-case/special-posterior-transdata-results.RDS")
special.acceptance.rate.results = readRDS(file = "special-case/special-acceptance-results.RDS")

## Build pretty trace plots
png(filename = "cav-traceplots.png", width = 480, height = 480, units = "px", pointsize = 12,
    bg = "white")
par(mfrow = c(2,2), mar = c(4,4,1,1) + 0.1)
plot(posterior.results[,1], ylab = "Values", 
     xlab = "", main = "", axes = FALSE); axis(side = 1); axis(side = 2)
mtext(expression(paste("Trace plot for ", nu[11])), side = 3, line = 0, cex = 0.75)
plot(posterior.results[,2], ylab = "Values", xlab = "", axes = FALSE); axis(side = 1); axis(side = 2)
mtext(expression(paste("Trace plot for ", nu[12])), side = 3, line = 0, cex = 0.75)
plot(posterior.results[,3], ylab = "Values", xlab = "Iteration", axes = FALSE); axis(side = 1); axis(side = 2)
mtext(expression(paste("Trace plot for ", gamma[21])), side = 3, line = 0, cex = 0.75)
plot(posterior.results[,4], ylab = "Values", xlab = "Iteration", axes = FALSE); axis(side = 1); axis(side = 2)
mtext(expression(paste("Trace plot for ", gamma[22])), side = 3, line = 0, cex = 0.75)
dev.off()

## Trace plots -> Toss first 100
min.iter = 100; max.iter = nrow(posterior.results)

## Build 5% - 95% quantiles
## And plot histograms
quantile(posterior.results[min.iter:max.iter,1], probs = c(0.05, 0.50, 0.95)); mean(posterior.results[min.iter:max.iter,1])
quantile(special.posterior.results[min.iter:max.iter,1], probs = c(0.05, 0.50, 0.95)); mean(special.posterior.results[min.iter:max.iter,1])

quantile(posterior.results[min.iter:max.iter,2], probs = c(0.05, 0.50, 0.95)); mean(posterior.results[min.iter:max.iter,2])
quantile(special.posterior.results[min.iter:max.iter,2], probs = c(0.05, 0.50, 0.95)); mean(special.posterior.results[min.iter:max.iter,2])

quantile(posterior.results[min.iter:max.iter,3], probs = c(0.05, 0.50, 0.95)); mean(posterior.results[min.iter:max.iter,3])
quantile(special.posterior.results[min.iter:max.iter,3], probs = c(0.05, 0.50, 0.95)); mean(special.posterior.results[min.iter:max.iter,3])

quantile(posterior.results[min.iter:max.iter,4], probs = c(0.05, 0.50, 0.95)); mean(posterior.results[min.iter:max.iter,4])
quantile(special.posterior.results[min.iter:max.iter,4], probs = c(0.05, 0.50, 0.95)); mean(special.posterior.results[min.iter:max.iter,4])


png(filename = "cav-quantile-plots.png", width = 480, height = 480, units = "px", pointsize = 12,
    bg = "white")
par(mfrow = c(2,2), mar = c(4,4,1,1) + 0.1)
plot(density(posterior.results[min.iter:max.iter,1]), axes = FALSE, main = "", xlab = ""); axis(side = 1); axis(side = 2)
mtext(expression(paste("Density for ", nu[11])), side=1, line = 3, cex=0.75)
lines(density(special.posterior.results[min.iter:max.iter,1]), col = "red")

plot(density(posterior.results[min.iter:max.iter,2]), axes = FALSE, main = "", xlab = ""); axis(side = 1); axis(side = 2)
mtext(expression(paste("Density for ", nu[12])), side=1, line = 3, cex=0.75)
lines(density(special.posterior.results[min.iter:max.iter,2]), col = "red")

plot(density(posterior.results[min.iter:max.iter,3]), axes = FALSE, main = "", xlab = ""); axis(side = 1); axis(side = 2)
mtext(expression(paste("Density for ", gamma[21])), side=1, line = 3, cex=0.75)
lines(density(special.posterior.results[min.iter:max.iter,3]), col = "red")

plot(density(posterior.results[min.iter:max.iter,4]), axes = FALSE, main = "", xlab = ""); axis(side = 1); axis(side = 2)
mtext(expression(paste("Density for ", gamma[22])), side=1, line = 3, cex=0.75)
lines(density(special.posterior.results[min.iter:max.iter,4]), col = "red")
dev.off()

## Survival curves from baseline
require(Matrix); Survival.list = list(); rho = 10

if(!file.exists("survival-fits.RDS")) {
  for (iter in min.iter:max.iter) {
    current.global = person.global.data.list[[iter]]
    current.params = posterior.results[,iter]
    Survival.list[[iter]] = Survival.fit(current.global, current.params, grid.length = 0.1)
  }
  saveRDS(Survival.list, file = "survival-fits.RDS")
} else {
  Survival.list = readRDS("survival-fits.RDS")
}

### Plot truth
simdata = readRDS("simdata.rds"); simdata = data.frame(simdata)
names(simdata) = c("time", "num.healthy", "num.ill", "num.dead", "switch1", "switch2", "failure.time")
truth = c(1/2, 1/5, 0.65, 1.71)
source("functions.R"); rho = 10
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

png(filename = "survival-plots.png", width = 480, height = 240, units = "px", pointsize = 12,
    bg = "white")
par(mfrow = c(1,1), mar = c(3,4,1,1) + 0.1)
plot(survival.grid, 1 - quantiles.survival.healthy[,2], type = "l", axes = FALSE, ylab = "Survival", xlab = "Time", cex.axis = 0.75)
axis(side = 1); axis(side = 2)
lines(survival.grid, 1 - quantiles.survival.healthy[,1], lty = 2)
lines(survival.grid, 1 - quantiles.survival.healthy[,3], lty = 2)
lines(true.Survival[obs.Survival,1], 1 - true.Survival[obs.Survival,2], col = "red")

lines(survival.grid, 1 - quantiles.survival.ill[,2], lty = 1)
lines(survival.grid, 1 - quantiles.survival.ill[,1], lty = 2)
lines(survival.grid, 1 - quantiles.survival.ill[,3], lty = 2)
lines(true.Survival[obs.Survival,1], 1 - true.Survival[obs.Survival,3], col = "red")
text(2.4,0.85, "Healthy at baseline", cex = 0.75)
text(4,0.16, "Ill at baseline", cex = 0.75)

dev.off()
