## Data Generating Mechanism
## Source code to generate toy example

source("functions.R")

set.seed("1813")
# N = 200
# 140 healthy; 60 ill
n.H = 150
n.I = 100
n.D = 0
rho = 10 # Set for this example
num.iters = 1000
time = 0
record = c(time, n.H, n.I, n.D, 0, 0, 0)

# Marginal holding rate healthy to ill is 2 years
# Set beta.HI so holding rate from ill to healthy after 3 years 
nu.HI = 1/2 * rho; beta.HI = 0.65
nu.HI * zeta(1,0,beta.HI,rho) # Exp rate (1/mean) for switching from healthy to ill
nu.HI * zeta(0,1,beta.HI,rho) # Exp rate (1/mean) for switching from ill to healthy

# Marginal survival given healthy is 5 years
# Marginal survival given ill is only 3 years
nu.HD = 1/5 * rho; beta.ID = 1.71
nu.HD * zeta(1,0,beta.ID,rho) # Exp rate (1/mean) for switching from healthy to ill
nu.HD * zeta(0,1,beta.ID,rho) # Exp rate (1/mean) for switching from ill to healthy

for(iter in 1:num.iters) {
  # Generate holding times
  rate.HI = nu.HI * zeta(n.H, n.I, beta.HI, rho)
  rate.HD = nu.HD * zeta(n.H, n.I, beta.ID, rho)
  
  randomhold.HI = rexp(1, rate.HI)
  randomhold.HD = rexp(1, rate.HD)
  minhold = min(randomhold.HI, randomhold.HD)
  
  if(randomhold.HI == minhold) {
    switch = sample.split(n.H,n.I,beta.HI,rho)
    n.H = n.H - switch[1] + switch[2]; n.I = n.I - switch[2] + switch[1]
  } else {
    switch = sample.split(n.H,n.I,beta.ID,rho)
    n.D = n.D + sum(switch)
    n.H = n.H - switch[1]; n.I = n.I - switch[2]
  }
  time = time + minhold
  state = c(n.H, n.I, n.D, switch)
  death.time = as.numeric(minhold == randomhold.HD)
  record = rbind(record, c(time, state, death.time))
  
  if(n.H == 0 & n.I == 0 ) {break}
  print(iter)
  
}

record = data.frame(as.matrix(record))
names(record) = c("time", "num.healthy", "num.ill", "num.dead", "switch1", "switch2", "failure.time")

saveRDS(record, "simdata.RDS")
# 
# death.record = record[record[,5] == 1,1:4]
# 
# death.record = rbind(c(0,140,60,0), death.record)
# 
# #  MLE estimate of nu
# N = 200
# num.blocks = total = 0
# for(row in 2:nrow(death.record)) {
#   gap = death.record[row,1] - death.record[row-1,1]
#   value = digamma(N-death.record[row-1,4] + rho) - digamma(rho)
#   total = total + value*gap
#   num.blocks = num.blocks + 1
# }
# 
# nu.hat = (total/num.blocks)^(-1)
# 
# print(nu.hat)
