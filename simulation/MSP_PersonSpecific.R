## Source stuff
source("functions.R")
simdata = readRDS("simdata.rds")
simdata = data.frame(simdata)
names(simdata) = c("time", "num.healthy", "num.ill", "num.dead", "switch1", "switch2", "failure.time")
rho = 10 # Assume rho is fixed throughout

## LABELLING IS 1 = HEALTHY, 2 = ILL, 3 = DEATH
set.seed(140710)
N = 250
all.persondata = c(); simdata.new = simdata
for (i in 1:N) {
  temp = construct.person(simdata.new)
  current.total = rowSums(temp[[2]][,2:4])[1]
  all.persondata = rbind(all.persondata, cbind(N-current.total,temp[[1]]) )
  simdata.new = temp[[2]]
}

all.persondata = data.frame(all.persondata)
names(all.persondata) = c("person", "time", "state")

saveRDS(all.persondata, "person-simdata.rds")

all.persondata.LOCF = c()

for (user in 1:N){ 
  user.temp = all.persondata[all.persondata$person == user,2:3]
  time.temp = c(0:floor(max(user.temp$time)), max(user.temp$time))
  obs.temp = sapply(time.temp, previous.obs, user.temp)
  all.persondata.LOCF = rbind(all.persondata.LOCF, cbind(user, time.temp, obs.temp))
}

saveRDS(all.persondata.LOCF, "person-simdata-LOCF.rds")

