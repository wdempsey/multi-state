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

plot(summary.km.fit$time, summary.km.fit$surv, type = "l")
lines(summary.km.fit$time, summary.km.fit$lower, lty = 2)
lines(summary.km.fit$time, summary.km.fit$upper, lty = 2)

print(km.fit, print.rmean=TRUE)


### Death data at five analysis using Cox PH
death.data.after.cutoff = c(); cutoff = 5

for(i in 1:length(unique.users)) {
  temp.data = person.data[person.data$user == i,]
  if( any(temp.data$time > cutoff) ) {
    prior.to.cutoff = temp.data[max( which(temp.data$time <= cutoff)) ,]
    final.obs.temp = temp.data[which(temp.data$time == max(temp.data$time)),]
    death.data.after.cutoff = rbind(death.data.after.cutoff, c(final.obs.temp$user, final.obs.temp$time, prior.to.cutoff$state, final.obs.temp$state == 4))
  }
}

death.data.after.cutoff = data.frame(death.data.after.cutoff)
names(death.data.after.cutoff) = c("user", "time", "prior.state", "status")

death.data.after.cutoff$time = death.data.after.cutoff$time - cutoff

library(survival)

cox.fit = coxph(Surv(time,status) ~ as.factor(prior.state), death.data.after.cutoff)

summary(cox.fit)

plot(survfit(cox.fit, newdata=data.frame(prior.state=1)),
     xlab = "Years", ylab="Survival") 

lines(survfit(cox.fit, newdata=data.frame(prior.state=2)),
      xlab = "Years", ylab="Survival", col = "red") 

lines(survfit(cox.fit, newdata=data.frame(prior.state=3)),
      xlab = "Years", ylab="Survival", col = "blue") 

### Let's look at global and see if we can figure out the issue
posterior.results = readRDS(file = "posterior-params-results.RDS")
person.global.data.list = readRDS(file = "posterior-transdata-results.RDS")
