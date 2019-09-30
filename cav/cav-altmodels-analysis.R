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
init.distribution = c(1,0,0,0)
current.distribution = init.distribution
result.distribution = matrix(nrow = nrow(global.data), ncol = 4)
result.distribution[1,] = current.distribution

for (iter in 2:nrow(global.data)) {
  
  # IS IT A FAILURE TIME
  is_fail_time = global.data$failure.time[iter] == 1
  ## Calc transition to failure from each state
  if(is_fail_time) {
    stay_nocav = global.data$num.nocav[iter]/global.data$num.nocav[iter-1]
    stay_mild = global.data$num.mild[iter]/global.data$num.mild[iter-1]; if(is.nan(stay_mild)) {stay_mild=1}
    stay_severe = global.data$num.severe[iter]/global.data$num.severe[iter-1]; if(is.nan(stay_severe)) {stay_severe=1}
    
    current.distribution[4] = current.distribution[4] + sum(current.distribution[1:3] * (1 - c(stay_nocav, stay_mild, stay_severe)))
    current.distribution[1:3] = current.distribution[1:3] * c(stay_nocav, stay_mild, stay_severe)
    current.distribution
  } else {
    nocav_to_mild = global.data$switch1[iter]/global.data$num.nocav[iter-1]; if(is.nan(nocav_to_mild)){nocav_to_mild = 0}
    mild_to_nocav = global.data$switch2a[iter]/global.data$num.mild[iter-1]; if(is.nan(mild_to_nocav)){mild_to_nocav = 0}
    mild_to_severe = global.data$switch2b[iter]/global.data$num.mild[iter-1]; if(is.nan(mild_to_severe)){mild_to_severe = 0}
    severe_to_mild = global.data$switch3[iter]/global.data$num.severe[iter-1]; if(is.nan(severe_to_mild)){severe_to_mild = 0}
    
    if(mild_to_severe == Inf) {
      current.distribution[1] = current.distribution[1] * (1-nocav_to_mild) 
      current.distribution[3] = current.distribution[3] + current.distribution[1] * nocav_to_mild
    } else {
      current.distribution[1] = current.distribution[1] * (1-nocav_to_mild) + current.distribution[2] * mild_to_nocav
      current.distribution[2] = current.distribution[1] * nocav_to_mild + current.distribution[2] * (1- mild_to_nocav - mild_to_severe) + current.distribution[3] * severe_to_mild
      current.distribution[3] = current.distribution[3] * (1-severe_to_mild) + current.distribution[2] * mild_to_severe
    }
    
  }
  result.distribution[iter,] = current.distribution
}


plot(global.data$time, 1-result.distribution[,4], xlim = c(0,20), ylim = c(0,1), type = "l")


### FULLY PARAMETRIC MODEL
library(expm)

llik.parametric = logLik.iid(person.data)
parameters = rep(1/3, 9)
lower.bound = rep(0.00001, 9)
upper.bound = rep(10,9)
mle.parameter = optim(par = parameters, fn = llik.parametric, lower = lower.bound, upper = upper.bound)

mle.parameter
