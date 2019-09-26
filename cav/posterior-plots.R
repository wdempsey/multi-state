posterior.results = readRDS(file = "posterior-params-results.RDS")
person.global.data.list = readRDS(file = "posterior-transdata-results.RDS")
acceptance.rate.results = readRDS(file = "acceptance-results.RDS")
source("cav-functions.R")

max.iter = 1000

## Build pretty trace plots
png(filename = "cav-trace-plots.png", width = 480, height = 480, units = "px", pointsize = 12,
    bg = "white")
par(mfrow = c(4,2), mar = c(5,4,1,1) + 0.1)
plot(posterior.results[1,1:max.iter], ylab = "Values", xlab = "Iteration", axes = FALSE); 
axis(side = 1); axis(side = 2)
mtext(expression(paste("Trace plot for ", nu[11])), side = 3, line = -1, cex = 0.75)
plot(posterior.results[2,1:max.iter], ylab = "Values", xlab = "Iteration", axes = FALSE); 
axis(side = 1); axis(side = 2)
mtext(expression(paste("Trace plot for ", nu[12])), side = 3, line = -1, cex = 0.75)
plot(posterior.results[3,1:max.iter], ylab = "Values", xlab = "Iteration", axes = FALSE); 
axis(side = 1); axis(side = 2)
mtext(expression(paste("Trace plot for ", gamma[21])), side = 3, line = -1, cex = 0.75)
plot(posterior.results[4,1:max.iter], ylab = "Values", xlab = "Iteration", axes = FALSE); 
axis(side = 1); axis(side = 2)
mtext(expression(paste("Trace plot for ", gamma[31])), side = 3, line = -1, cex = 0.75)
plot(posterior.results[5,1:max.iter], ylab = "Values", xlab = "Iteration", axes = FALSE); 
axis(side = 1); axis(side = 2)
mtext(expression(paste("Trace plot for ", gamma[22])), side = 3, line = -1, cex = 0.75)
plot(posterior.results[6,1:max.iter], ylab = "Values", xlab = "Iteration", axes = FALSE); 
axis(side = 1); axis(side = 2)
mtext(expression(paste("Trace plot for ", gamma[32])), side = 3, line = -1, cex = 0.75)
plot(posterior.results[7,1:max.iter], ylab = "Values", xlab = "Iteration", axes = FALSE); 
axis(side = 1); axis(side = 2)
mtext(expression(paste("Trace plot for ", alpha)), side = 3, line = -1, cex = 0.75)
dev.off()
## Trace plots -> Toss first 100
min.iter = 50 #; max.iter = ncol(posterior.results)

## Build 5% - 95% quantiles
## And plot histograms
for (i in 1:7) {
  print( c(quantile(posterior.results[i,min.iter:max.iter], probs = c(0.05, 0.50, 0.95)), mean(posterior.results[i,min.iter:max.iter]) ) )
}

png(filename = "cav-quantile-plots.png", width = 480, height = 480, units = "px", pointsize = 12,
    bg = "white")
par(mfrow = c(4,2), mar = c(4,4,1,1) + 0.1)
plot(density(posterior.results[1,min.iter:max.iter]), axes = FALSE, main = "", xlab = expression(paste("Trace plot for ", nu[11]))); 
axis(side = 1); axis(side = 2)
plot(density(posterior.results[2,min.iter:max.iter]), axes = FALSE, main = "", xlab = expression(paste("Trace plot for ", nu[12]))); 
axis(side = 1); axis(side = 2)
plot(density(posterior.results[3,min.iter:max.iter]), axes = FALSE, main = "", xlab = expression(paste("Trace plot for ", gamma[21]))); 
axis(side = 1); axis(side = 2)
plot(density(posterior.results[4,min.iter:max.iter]), axes = FALSE, main = "", xlab = expression(paste("Trace plot for ", gamma[31]))); 
axis(side = 1); axis(side = 2)
plot(density(posterior.results[5,min.iter:max.iter]), axes = FALSE, main = "", xlab = expression(paste("Trace plot for ", gamma[22]))); 
axis(side = 1); axis(side = 2)
plot(density(posterior.results[6,min.iter:max.iter]), axes = FALSE, main = "", xlab = expression(paste("Trace plot for ", gamma[32])));
axis(side = 1); axis(side = 2)
plot(density(posterior.results[7,min.iter:max.iter]), axes = FALSE, main = "", xlab = expression(paste("Trace plot for ", alpha))); 
axis(side = 1); axis(side = 2)
dev.off()

## Survival curves from baseline
require(Matrix); Survival.list = Survival.list.at.five = list()
rho = 10

if(!file.exists("survival-fits.RDS")) {
  for (iter in min.iter:max.iter) {
    init.time = 0
    current.global = person.global.data.list[[iter]]
    current.params = posterior.results[,iter]
    Survival.list[[iter]] = Survival.fit(current.global, current.params, rho, init.time, grid.length = 0.1)
    
    init.time = 5
    global.at.five = current.global[current.global$time > init.time,]
    temp = current.global[max(which(current.global$time <= init.time)),]
    temp$time = init.time; temp[6:10] = 0
    global.at.five = rbind(temp, global.at.five)
    Survival.list.at.five[[iter]] = Survival.fit(global.at.five, current.params, rho, init.time, grid.length = 0.1)
  }
  saveRDS(Survival.list, file = "survival-fits.RDS")
  saveRDS(Survival.list.at.five, file = "survival-fits-at-five.RDS")
} else {
  Survival.list = readRDS("survival-fits.RDS")
  Survival.list.at.five = readRDS("survival-fits-at-five.RDS")
}

### Compute quantiles
survival.grid = seq(5,20,0.1)
nocav.survival.plot.matrix = mild.survival.plot.matrix = 
  severe.survival.plot.matrix = matrix(nrow = length(survival.grid), ncol = max.iter - min.iter + 1)


for (iter in min.iter:max.iter) {
  temp.surv = Survival.list.at.five[[iter]]
  
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

png(filename = "cav-survival-plots-at-five.png", width = 480, height = 320, units = "px", pointsize = 12,
    bg = "white")
par(mfrow = c(1,1), mar = c(4,4,1,1) + 0.1)
plot(survival.grid, 1 - quantiles.survival.nocav[,2], type = "l", 
     axes = FALSE, ylab = "Survival", xlab = "Time", ylim = c(0.0,1.0))
axis(side = 1); axis(side = 2)
# lines(survival.grid, 1 - quantiles.survival.nocav[,1], lty = 2)
# lines(survival.grid, 1 - quantiles.survival.nocav[,3], lty = 2)

lines(survival.grid, 1 - quantiles.survival.mild[,2], col = "red")
# lines(survival.grid, 1 - quantiles.survival.mild[,1], col = "red", lty = 2)
# lines(survival.grid, 1 - quantiles.survival.mild[,3], col = "red", lty = 2)

lines(survival.grid, 1 - quantiles.survival.severe[,2], col = "blue")
lines(survival.grid, 1 - quantiles.survival.severe[,1], col = "blue", lty = 2)
lines(survival.grid, 1 - quantiles.survival.severe[,3], col = "blue", lty = 2)
dev.off()

sum( (1-quantiles.survival.nocav[,2]) * 0.1)

sum( (1-quantiles.survival.mild[,2]) * 0.1)

sum( (1-quantiles.survival.severe[,2]) * 0.1)


### Compute quantiles
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

## ADD KM
library(msm)
person.data = data.frame(cbind(cav$PTNUM, cav$years, cav$state))
names(person.data) = c("user", "time", "state")
unique.users = unique(person.data$user)

death.data = c()

for(i in 1:length(unique.users)) {
  person.data$user[person.data$user == unique.users[i]] = i
  temp.data = person.data[person.data$user == i,]
  final.obs.temp = temp.data[which(temp.data$time == max(temp.data$time)),]
  death.data = rbind(death.data, c(final.obs.temp$user, final.obs.temp$time, final.obs.temp$state == 4))
}

death.data = data.frame(death.data)
names(death.data) = c("user", "time", "status")

library(survival)

km.fit = survfit(Surv(time,status) ~ 1, death.data)

summary.km.fit = summary(km.fit)

sum(diff(summary.km.fit$time)*summary.km.fit$surv[-242])

plot(summary.km.fit$time, summary.km.fit$surv, type = "l")
lines(summary.km.fit$time, summary.km.fit$lower, lty = 2)
lines(summary.km.fit$time, summary.km.fit$upper, lty = 2)

png(filename = "cav-survival-plots.png", width = 480, height = 320, units = "px", pointsize = 12,
    bg = "white")
par(mfrow = c(1,1), mar = c(4,4,1,1) + 0.1)
plot(survival.grid, 1 - quantiles.survival.nocav[,2], type = "l", lty = 1,
     axes = FALSE, ylab = "Survival", xlab = "Time", ylim = c(0.0,1.0))
axis(side = 1); axis(side = 2)
# lines(survival.grid, 1 - quantiles.survival.nocav[,1], lty = 2)
# lines(survival.grid, 1 - quantiles.survival.nocav[,3], lty = 2)

lines(summary.km.fit$time, summary.km.fit$surv, type = "l", lwd = 2, lty = 2)

lines(survival.grid, 1 - quantiles.survival.mild[,2], col = "red")
# lines(survival.grid, 1 - quantiles.survival.mild[,1], col = "red", lty = 2)
# lines(survival.grid, 1 - quantiles.survival.mild[,3], col = "red", lty = 2)

lines(survival.grid, 1 - quantiles.survival.severe[,2], col = "blue")
lines(survival.grid, 1 - quantiles.survival.severe[,1], col = "blue", lty = 2)
lines(survival.grid, 1 - quantiles.survival.severe[,3], col = "blue", lty = 2)
dev.off()

sum( (1-quantiles.survival.nocav[,2]) * 0.1)

sum( (1-quantiles.survival.mild[,2]) * 0.1)

sum( (1-quantiles.survival.severe[,2]) * 0.1)
