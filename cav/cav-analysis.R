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

## Initialize the true observations
## equal to the transition data
all.global.trans.data = global.data
all.person.trans.data = all.person.obs.data = person.data

## Initialize using trans.data (i.e., under LOCF assumption)
## Gives biased but OK starting values
nu.11 = nu.12 = 0.2; beta.21 = 1.0; beta.22 = 1; beta.31 = 0.5; beta.32 = 1
p.21 = 0.2
params = c(nu.11, nu.12, beta.21, beta.31, beta.22, beta.32, p.21)

### CHECK: RUN UNIFORMIZATION WITH TRUTH 
### SHOULD CONVERGE TO SENSIBLE SAMPLED TRAJECTORIES
### RUN 500 times, and then estimate the parameters 
### using total.llik -> should be close to truth
iter = 1; max.iter = 1000; smooth.to.iter = 500
init.multiple = 5; final.multiple = 2
current.params = params
prior.params = data.frame(alpha = 1/2, beta = 1/2, log.mu = 0, log.sd = 1, 
                          step.size = 0.2, alpha.alive = 0.2*10, beta.alive = 0.8*10)
all.users = unique(all.person.obs.data$user)
posterior.results = matrix(nrow = length(current.params), ncol = max.iter)
acceptance.rate.results = vector(length = max.iter)
person.global.data.list = list()
set.seed("1283714")
save.every.num.iters = 50 
update.forwback.every.num.iters = 25

for (iter in 1:max.iter) {
  print(iter)
  current.multiple = max(final.multiple, init.multiple + (final.multiple - init.multiple) / (smooth.to.iter-1) * (iter-1))
  ## Cycle through all users once for now
  ## GENERATE INITIAL TRANSITIONS
  is.it.initial.iter = as.logical(iter == 1)
  if(iter %% update.forwback.every.num.iters == 1) {
    new.all.person.trans.data = sample.all.new.users(all.person.obs.data, all.person.trans.data, 
                                                     censoring.data, current.params, current.multiple, rho, 
                                                     initial.iter = is.it.initial.iter) 
    new.all.person.global.data = convert.to.global(new.all.person.trans.data, censoring.data)
  }
  new.nu.params = nu.gibbs.updates(new.all.person.global.data, current.params, rho, prior.params) # also updates.alive.alive P
  new.alpha.alive = alpha.alive.gibbs.updates(new.all.person.global.data, prior.params)
  new.beta.params = beta.gibbs.updates(new.all.person.global.data, new.nu.params, 
                                       new.alpha.alive, current.params, rho, prior.params)
  print(paste("For iteration = ", iter,". We have parameters", sep = ""))
  print(new.beta.params$new.params)
  posterior.results[,iter] = new.beta.params$new.params
  acceptance.rate.results[iter] = mean(new.beta.params$accept.rate)
  current.params = new.beta.params$new.params
  all.person.trans.data = new.all.person.trans.data
  person.global.data.list[[iter]] = new.all.person.global.data
  
  if ((iter %% save.every.num.iters) == 0) {
    saveRDS(posterior.results, file = "posterior-params-results.RDS")
    saveRDS(person.global.data.list, file = "posterior-transdata-results.RDS")
    saveRDS(acceptance.rate.results, file = "acceptance-results.RDS")
  }
}


