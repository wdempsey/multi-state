## Source code to analyze toy example
## Under multi-censored observations
## Made every year while participant is alive

# Source Functions and observation data
source("functions.R")
person.simdata.LOCF = data.frame(readRDS("person-simdata-LOCF.rds"))
names(person.simdata.LOCF) = c("user", "time", "state")

global.simdata.LOCF = convert.to.global(person.simdata.LOCF)

## Remove components of global that are non-switch times
## Except for the initial observatoin
global.simdata.LOCF = global.simdata.LOCF[-which(global.simdata.LOCF$switch1 == 0 & 
                                                   global.simdata.LOCF$switch2 == 0)[-1],] 

rho = 10 # Assume rho is fixed throughout

## Initialize the true observations
## equal to the transition data
all.global.trans.data = global.simdata.LOCF
all.person.trans.data = all.person.obs.data = person.simdata.LOCF

## Initialize using trans.data (i.e., under LOCF assumption)
## Gives biased but OK starting values
nu.HI = nu.ID = beta.HI = beta.ID = 1
params = c(nu.HI, nu.ID, beta.HI, beta.ID)

temp = optim(params, fn = total.llik(all.global.trans.data,rho), lower = rep(0.01,4), upper = rep(5,4))
init.par = temp$par
truth = c(1/2, 1/5, 0.65, 1.71)

### CHECK: RUN UNIFORMIZATION WITH TRUTH 
### SHOULD CONVERGE TO SENSIBLE SAMPLED TRAJECTORIES
### RUN 500 times, and then estimate the parameters 
### using total.llik -> should be close to truth
iter = 1; max.iter = 1000; smooth.to.iter = 100
init.multiple = 5; final.multiple = 2
current.params = truth
prior.params = data.frame(alpha = 1, beta = 1, log.mu = 0, log.sd = 1, step.size = 0.15)
all.users = unique(all.person.obs.data$user)
posterior.results = matrix(nrow = length(truth), ncol = max.iter)
acceptance.rate.results = vector(length = max.iter)
person.global.data.list = list()
set.seed("1283714")

for (iter in 1:max.iter) {
  current.multiple = max(final.multiple, init.multiple + (final.multiple - init.multiple) / (smooth.to.iter-1) * (iter-1))
  ## Cycle through all users once for now
  ## GENERATE INITIAL TRANSITIONS
  is.it.initial.iter = as.logical(iter == 1)
  new.all.person.trans.data = sample.all.new.users(all.person.obs.data, all.person.trans.data, 
                                                   current.params, current.multiple, rho, 
                                                   initial.iter = is.it.initial.iter) 
  new.all.person.global.data = convert.to.global(new.all.person.trans.data)
  new.nu.params = nu.gibbs.updates(new.all.person.global.data, current.params, rho, prior.params)
  new.beta.params = beta.gibbs.updates(new.all.person.global.data, new.nu.params, current.params, rho, prior.params)
  print(paste("For iteration = ", iter,". We have parameters", sep = ""))
  print(new.beta.params$new.params)
  posterior.results[,iter] = new.beta.params$new.params
  acceptance.rate.results[iter] = new.beta.params$accept.rate
  current.params = new.beta.params$new.params
  all.person.trans.data = new.all.person.trans.data
  person.global.data.list[[iter]] = new.all.person.global.data
}

saveRDS(posterior.results, file = "posterior-params-results.RDS")
saveRDS(person.global.data.list, file = "posterior-transdata-results.RDS")
saveRDS(acceptance.rate.results, file = "acceptance-results.RDS")

