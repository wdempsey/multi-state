source("functions.R")
person.simdata.LOCF = data.frame(readRDS("person-simdata-LOCF.rds"))
names(person.simdata.LOCF) = c("user", "time", "state")

simdata = convert.to.global(person.simdata.LOCF)

simdata = simdata[-which(simdata$switch1 == 0 & simdata$switch2 == 0)[-1],] # Remove non-switch times

rho = 10 # Assume rho is fixed throughout

nu.HI = nu.ID = beta.HI = beta.ID = 1
params = c(nu.HI, nu.ID, beta.HI, beta.ID)

temp = optim(params, fn = total.llik(simdata,rho), lower = rep(0.01,4), upper = rep(5,4))
print(temp$par)
truth = c(1/2, 1/5, 0.65, 1.71)


library(numDeriv)
Hess.hat = hessian(func = total.llik(simdata,rho), x = temp$par)

std.err = sqrt(diag(solve(Hess.hat)))

temp$par + 1.96*std.err
temp$par - 1.96*std.err
truth

