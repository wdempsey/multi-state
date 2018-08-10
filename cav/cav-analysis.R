## Source code to analyze cav example
## Under multi-censored observations
## Made every year while participant is alive

# Source Functions and observation data
source("cav-functions.R")

library(msm)
person.data = data.frame(cbind(cav$PTNUM, cav$years, cav$state))
names(person.data) = c("user", "time", "state")
unique.users = unique(person.data$user)

for(i in 1:length(unique.users)) {
  person.data$user[person.data$user == unique.users[i]] = i
}

global.data = convert.to.global(person.data)

rho = 10 # Assume rho is fixed throughout

## Initialize the true observations
## equal to the transition data
all.global.trans.data = global.data
all.person.trans.data = all.person.obs.data = person.data

## Initialize using trans.data (i.e., under LOCF assumption)
## Gives biased but OK starting values
nu.11 = nu.12 = 0.2; beta.21 = beta.22 = beta.31 = beta.32 = 1
p.21 = 1/2; psi = 1/2
params = c(nu.11, nu.12, beta.21, beta.22, beta.31, beta.32, p.21)

bad.obs = rowSums(all.global.trans.data[,6:9]) == 0
all.global.trans.data = all.global.trans.data[!bad.obs, ]

### CHECK: RUN UNIFORMIZATION WITH TRUTH 
### SHOULD CONVERGE TO SENSIBLE SAMPLED TRAJECTORIES
### RUN 500 times, and then estimate the parameters 
### using total.llik -> should be close to truth
iter = 1; max.iter = 1000; smooth.to.iter = 100
init.multiple = 5; final.multiple = 2
current.params = params
prior.params = data.frame(alpha = 1, beta = 1, log.mu = 0, log.sd = 1, 
                          step.size = 0.1, alpha.step.size = 10)
all.users = unique(all.person.obs.data$user)
posterior.results = matrix(nrow = length(current.params), ncol = max.iter)
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
  new.nu.params = nu.gibbs.updates(new.all.person.global.data, current.params, rho, psi, prior.params)
  new.beta.params = beta.gibbs.updates(new.all.person.global.data, new.nu.params, 
                                       current.params, rho, psi, prior.params)
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

