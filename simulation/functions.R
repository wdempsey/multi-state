compute.pdf <- function(r1, d1, r2,d2,beta, rho) {
  f2 <- function(x) {
    # if(d1 > 0) {
      return( x^(r1) * (1-x)^(d1) * x^(r2*beta) * (1-x^(beta))^(d2) * x^(rho-1)  * (1-x)^(-1) )
    # } else {
      # return( x^(r1) * (1-x)^(d1) * x^(r2*beta) * (1-x^(beta))^(d2-1) * x^(rho-1) )
    # }
  }
  return(f2)
}

q.split <- function(r1,d1,r2,d2, beta, rho, tol = 0.000001) {
  temp = compute.pdf(r1,d1,r2,d2, beta, rho)
  return(integrate(temp,0.00,tol)$value + integrate(temp,tol,1)$value)
}

compute.normalizing <- function(n1, n2, beta, rho) {
  f2 <- function(x) {
    return( x^(rho-1) * ( 1-x^(n1) * x^(beta * n2) ) * (1-x)^(-1) )
  }
  return(f2)
}

zeta <- function(n1, n2, beta, rho) {
  temp = compute.normalizing(n1, n2, beta, rho)
  return(integrate(temp,0,1)$value)
}

sample.split <- function(n1,n2,beta,rho) {
  results = c()
  
  for(d1 in 0:n1) {
    for(d2 in 0:n2) {
      if( !(d1 == 0 & d2 == 0) ) {
        temp = q.split(n1-d1,d1,n2-d2,d2, beta, rho) * choose(n1, d1) * choose(n2,d2)
        results = rbind(results, c(d1, d2, temp))
      }
    }
  }
  # sum(results[,3]) 
  # zeta(n1, n2, beta, rho)
  results[,3] = results[,3]/ sum(results[,3])
  row.sample = sample(1:nrow(results), size = 1, prob = results[,3])
  return(results[row.sample,1:2])
}

lik.component <- function(params, n1, n2, d1, d2, isdeathtime, window.time, rho) {
  r1 = n1-d1; r2 = n2-d2;
  rate1 = params[1] * zeta(n1, n2, params[3], rho) * rho
  rate2 = params[2] * zeta(n1, n2, params[4], rho) * rho
  total.rate =  rate1 + rate2
  if(isdeathtime) {
    nu = params[2]*rho; beta = params[4]
  } else {
    nu = params[1]*rho; beta = params[3]
  }
  
  return(
    log(nu) - total.rate*window.time + log(q.split(r1,d1,r2,d2, beta, rho, tol = 0.001))
  )
}

total.llik <- function(data, rho) {
  subfunction <- function(params) {
    total = 0
    for(t in 2:nrow(data)) {
      isdeathtime = data$failure.time[t]
      n1 = data$num.healthy[t-1]; n2 = data$num.ill[t-1]
      d1 = data$switch1[t]; d2 = data$switch2[t]
      window.length = data$time[t] - data$time[t-1]
      total = total + lik.component(params, n1, n2, d1, d2, isdeathtime, window.length, rho)
    }
    return(-total)
  }
  return(subfunction)
}

construct.person <- function(simdata) {
  simdata.new = simdata
  
  if(all(simdata[1,2:3] > 0)) {
    init = sample(1:2,size = 1)
  } else {
    init = which(simdata[1,2:3] > 0)
  }
  
  persondata = c(simdata$time[1], init)
  current.state = init; alive = TRUE; row = 1
  simdata.new[row, 1+current.state] = simdata.new[row, 1+current.state] - 1
  row = row+1
  
  
  while(alive) {
    prob = simdata[row, 4+current.state]/simdata[row-1, 1+current.state]
    did.switch.happen = runif(1) < prob
    if (did.switch.happen) {
      ## Record User Data
      simdata.new[row, 4+current.state] = simdata.new[row, 4+current.state] - 1
      if(simdata$failure.time[row]) {
        current.state = 3; alive = FALSE
      } else{
        if(current.state == 2) {current.state = 1} else {current.state = 2}
      }
      temp = c(simdata$time[row],current.state)
      persondata = rbind(persondata, temp)
    } 
    simdata.new[row,1 + current.state] = simdata.new[row,1 + current.state] - 1
    row = row + 1
  }
  if (row < nrow(simdata.new)) {
    simdata.new[row:nrow(simdata.new),1 + current.state] = simdata.new[row:nrow(simdata.new),1 + current.state] - 1  
  }
  
  return( list("persondata"=persondata,
               "simdata"=simdata.new)
  )
  
}

previous.obs <- function(time, user.temp) {
  temp = user.temp$state[user.temp$time <= time]
  return(temp[length(temp)])
}

convert.to.global <- function(person.data) {
  ## Convert data to global configuration
  ## Assume each participant starts with initial obs
  ## At t = 0
  current.time = 0.0
  states = person.data$state[person.data$time == current.time]
  num.healthy = sum(states == 1); num.ill = sum(states == 2); num.dead = sum(states == 3)
  time = 0.0; switch1 = 0; switch2 = 0; failure.time = 0
  init = matrix(c(time, num.healthy, num.ill, num.dead, switch1, switch2, failure.time),
                nrow = 1)
  simdata = data.frame(init)
  names(simdata) = c("time", "num.healthy", "num.ill", "num.dead", "switch1", "switch2", "failure.time")
  if(nrow(person.data) == 0) {max.time = 0} else {max.time = max(person.data[,2])}
  
  
  while (current.time < max.time) {
    next.time = min(person.data$time[person.data$time > current.time])
    users.list = person.data$user[person.data$time == next.time]
    next.states = person.data$state[person.data$time == next.time]
    
    temp = person.data[is.element(person.data$user,users.list) & person.data$time < next.time,]
    switch1 = 0; switch2 = 0; failure.time = as.numeric(any(next.states==3))
    
    for (user in users.list) {
      temp.time = temp$time[temp$user == user]
      temp.state = temp$state[temp$user == user]
      previous.state = temp.state[temp.time == max(temp.time)]
      next.state = next.states[users.list == user]
      if(previous.state == 1 & next.state != 1) {
        switch1 = switch1 + 1
      } else if (previous.state == 2 & next.state != 2) {
        switch2 = switch2 + 1
      }
    }
    if (failure.time == 1) {
      num.healthy = num.healthy - switch1; num.ill = num.ill - switch2; num.dead = num.dead + switch1 + switch2
    } else {
      num.healthy = num.healthy - switch1 + switch2; num.ill = num.ill - switch2 + switch1; num.dead = num.dead
    }
    simdata = rbind(simdata, c(next.time, num.healthy, num.ill, 
                               num.dead, switch1, switch2, 
                               failure.time)
    )
    current.time = next.time
  }
  
  return(simdata)
}


### Gibbs functions
### Should move to separate file 
### Will help with readability

## Functions for computing a single iteration 
## of the MCMC sampler

extract.patient <- function(all.person.obs.data, all.person.trans.data, 
                            patient.id) {
  ## Extract information for chosen patient id 
  ## Construct other.trans.data of non-patient id
  ## records
  
  chosen.person.obs.data = all.person.obs.data[all.person.obs.data$user == patient.id,]
  chosen.person.trans.data = all.person.trans.data[all.person.trans.data$user == patient.id,]
  other.person.trans.data = all.person.trans.data[all.person.trans.data$user != patient.id,]
  return(list("chosen.person.obs.data" = chosen.person.obs.data,
              "chosen.person.trans.data" = chosen.person.trans.data,
              "other.person.trans.data" = other.person.trans.data)
  )
}

extract.Omega.and.A <- function(chosen.person.obs.data, chosen.person.trans.data, other.person.trans.data, current.params, current.multiple = 5, rho = 10) {
  ## Generate the vector of Omega's for each window of time
  ## This changes at each transition time
  other.global.trans.data = convert.to.global(other.person.trans.data)
  person.failure.time = chosen.person.obs.data$time[chosen.person.obs.data$state==3]
  other.global.trans.data.shortened = other.global.trans.data[other.global.trans.data$time < person.failure.time, ]
  
  ## Construct rate functions 
  rate1.addhealthy = current.params[1] * mapply(zeta, n1 = other.global.trans.data.shortened$num.healthy+1, 
                                                n2 = other.global.trans.data.shortened$num.ill, 
                                                beta = current.params[3], rho = rho) * rho
  
  rate2.addhealthy = current.params[2] * mapply(zeta, n1 = other.global.trans.data.shortened$num.healthy+1,
                                                n2 = other.global.trans.data.shortened$num.ill,
                                                beta = current.params[4], rho = rho) * rho
  
  rate1.healthy = current.params[1] * mapply(zeta, n1 = other.global.trans.data.shortened$num.healthy, 
                                             n2 = other.global.trans.data.shortened$num.ill, 
                                             beta = current.params[3], rho = rho) * rho
  
  rate2.healthy = current.params[2] * mapply(zeta, n1 = other.global.trans.data.shortened$num.healthy,
                                             n2 = other.global.trans.data.shortened$num.ill,
                                             beta = current.params[4], rho = rho) * rho
  
  total.rate.addhealthy = rate1.addhealthy + rate2.addhealthy - (rate1.healthy + rate2.healthy)
  
  
  rate1.addill = current.params[1] * mapply(zeta, n1 = other.global.trans.data.shortened$num.healthy, 
                                            n2 = other.global.trans.data.shortened$num.ill+1, 
                                            beta = current.params[3], rho = rho) * rho
  
  rate2.addill = current.params[2] * mapply(zeta, n1 = other.global.trans.data.shortened$num.healthy,
                                            n2 = other.global.trans.data.shortened$num.ill+1,
                                            beta = current.params[4], rho = rho) * rho
  
  rate1.ill = current.params[1] * mapply(zeta, n1 = other.global.trans.data.shortened$num.healthy, 
                                         n2 = other.global.trans.data.shortened$num.ill, 
                                         beta = current.params[3], rho = rho) * rho
  
  rate2.ill = current.params[2] * mapply(zeta, n1 = other.global.trans.data.shortened$num.healthy,
                                         n2 = other.global.trans.data.shortened$num.ill,
                                         beta = current.params[4], rho = rho) * rho
  
  total.rate.addill = rate1.addill + rate2.addill - (rate1.ill + rate2.ill)
  
  both.rates = cbind(total.rate.addhealthy, total.rate.addill)
  bar.Omega = current.multiple * apply(both.rates,1,max)
  bar.Omega.time = other.global.trans.data.shortened$time
  bar.Omega.time = c(bar.Omega.time, person.failure.time)
  
  ## NOW DO FOR ALL PEOPLE TO CONSTRUCT A.t
  all.global.trans.data = convert.to.global(rbind(other.person.trans.data, chosen.person.trans.data))
  all.global.trans.data.shortened = all.global.trans.data[all.global.trans.data$time < person.failure.time, ]
  
  ## Construct rate functions 
  rate1 = current.params[1] * mapply(zeta, n1 = all.global.trans.data.shortened$num.healthy, 
                                     n2 = all.global.trans.data.shortened$num.ill, 
                                     beta = current.params[3], rho = rho) * rho
  
  rate2 = current.params[2] * mapply(zeta, n1 = all.global.trans.data.shortened$num.healthy, 
                                     n2 = all.global.trans.data.shortened$num.ill, 
                                     beta = current.params[4], rho = rho) * rho
  
  ## Remove chosen person from all.global.trans.data.shortened
  temp.trans.data.shortened = all.global.trans.data.shortened
  for (i in 1:nrow(temp.trans.data.shortened)) {
    temp.which.obs = max(which(chosen.person.trans.data$time <= temp.trans.data.shortened$time[i]))
    temp.is.healthy = as.numeric(chosen.person.trans.data$state[temp.which.obs] == 1)
    if( temp.is.healthy == 1 ) {
      temp.trans.data.shortened$num.healthy[i] = temp.trans.data.shortened$num.healthy[i] - 1
    } else {
      temp.trans.data.shortened$num.ill[i] = temp.trans.data.shortened$num.ill[i] - 1
    }
  }
  
  rate1.At.baseline = current.params[1] * mapply(zeta, n1 = temp.trans.data.shortened$num.healthy, 
                                     n2 = temp.trans.data.shortened$num.ill, 
                                     beta = current.params[3], rho = rho) * rho
  
  rate2.At.baseline = current.params[2] * mapply(zeta, n1 = temp.trans.data.shortened$num.healthy, 
                                     n2 = temp.trans.data.shortened$num.ill, 
                                     beta = current.params[4], rho = rho) * rho
  
  bar.A = rate1 + rate2 - rate1.At.baseline - rate2.At.baseline
  bar.A.time = all.global.trans.data.shortened$time
  bar.A.time = c(bar.A.time, person.failure.time)

  return(list("bar.Omega" = bar.Omega, "bar.Omega.time" = bar.Omega.time, 
              "bar.A" = bar.A, "bar.A.time" = bar.A.time))
  
}

generate.points <- function(bar.Omega, bar.Omega.time, bar.A, bar.A.time, initial.iter) {
  ## Generate the vector of potential "new" transition times
  all.times = unique(c(bar.Omega.time, bar.A.time))
  all.times = all.times[order(all.times)]
  current.time = all.times[1]
  new.potential.trans.times = c()
  for (t in 1:(length(all.times)-1)) {
    current.Omega = bar.Omega[max(which(bar.Omega.time <= all.times[t]))]
    current.A = bar.A[max(which(bar.A.time <= all.times[t]))]
    if (initial.iter) {
      current.holding.rate = current.Omega
    }
    else{
      current.holding.rate = current.Omega - current.A
    }
    if(current.holding.rate < 0 ) {print ("Negative rate, by accident!")}
    while(current.time < all.times[t+1]) {
      new.potential.trans.times = c(new.potential.trans.times, current.time)
      holding.time = rexp(1,current.holding.rate)
      current.time = current.time + holding.time
    }
  }
  new.potential.trans.times = new.potential.trans.times[new.potential.trans.times!=0]
  new.potential.trans.times
  return(new.potential.trans.times)
}

construct.info <- function(new.potential.trans.times, chosen.person.trans.data, other.person.trans.data, chosen.person.obs.data, current.params) {
  ## Combine all the information necessary to perform 
  ## the filter forward, backward sampling algorithm.
  ## These are all trans times, obs times for patient,
  ## and all new potential trans times
  other.global.trans.data = convert.to.global(other.person.trans.data)
  person.failure.time = chosen.person.obs.data$time[chosen.person.obs.data$state==3]
  other.global.trans.data.shortened = other.global.trans.data[other.global.trans.data$time < person.failure.time, ]                                                              
  
  ## Add in all the other ppl trans times
  all.trans.times = other.global.trans.data.shortened[,c("time", "num.healthy", "num.ill", "switch1", "switch2", 
                                                         "failure.time")]
  all.trans.times$person.obs.times = FALSE
  all.trans.times$person.obs.states = NA
  all.trans.times$new.trans.time = FALSE
  
  ## Add in prior runs trans times for current person
  which.person.time.include = !is.element(chosen.person.trans.data$time, all.trans.times$time)
  temp.persontrans = chosen.person.trans.data[which.person.time.include,]
  temp.failure.time = as.numeric(temp.persontrans$state == 3)
  temp.persontrans = data.frame(time = temp.persontrans$time)
  if(nrow(temp.persontrans) != 0) {
    temp.persontrans$num.healthy = NA; temp.persontrans$num.ill = NA
    temp.persontrans$switch1 = 0; temp.persontrans$switch2 = 0
    temp.persontrans$person.obs.times = FALSE
    temp.persontrans$person.obs.states = NA
    temp.persontrans$new.trans.time = FALSE
    temp.persontrans$failure.time = temp.failure.time
    all.trans.times = rbind(all.trans.times,temp.persontrans)
  }    
  ## Add in potential new trans times
  if(length(new.potential.trans.times) != 0) {
    temp.newtrans = data.frame(time = new.potential.trans.times)
    temp.newtrans$num.healthy = NA; temp.newtrans$num.ill = NA
    temp.newtrans$switch1 = 0; temp.newtrans$switch2 = 0
    temp.newtrans$person.obs.times = FALSE
    temp.newtrans$person.obs.states = NA
    temp.newtrans$new.trans.time = TRUE
    temp.newtrans$failure.time = FALSE
    all.trans.times = rbind(all.trans.times, temp.newtrans)
  }  
  
  ## Add in obs times for current person
  temp.personobs = data.frame(time = chosen.person.obs.data$time, person.obs.states = chosen.person.obs.data$state)
  temp.personobs$num.healthy = NA; temp.personobs$num.ill = NA
  temp.personobs$switch1 = 0; temp.personobs$switch2 = 0
  temp.personobs$person.obs.times = TRUE
  temp.personobs$new.trans.time = FALSE
  temp.personobs$failure.time = as.numeric(temp.personobs$person.obs.states == 3)
  all.trans.times = rbind(all.trans.times, temp.personobs)
  
  # Correct for temporal ordering
  all.trans.times = all.trans.times[order(all.trans.times$time),]
  
  # NOTE: need to be careful around duplicate times
  # Kept for now.  Will check if necessary later
  return(all.trans.times)
}

filterfwd.backsampl <- function(all.trans.times, current.params, rho, bar.Omega, bar.time) {
  ### Alpha is a matrix 3 rows (for each state), with each column summing to 1
  ### Each column corresponds to a potential trans.time
  unique.trans.times = unique(all.trans.times$time)
  Alpha.matrix = matrix(nrow = 3, ncol = length(unique.trans.times))
  t = 1; current.time = unique.trans.times[t]; current.temp = all.trans.times[all.trans.times$time == current.time,]
  Alpha.matrix[,t] = c(150/250, 100/250, 0)
  list.of.current.Bs = list()
  L.matrix = matrix(nrow = 3, ncol = length(unique.trans.times))
  
  if(any(current.temp$person.obs.times == TRUE)) {
    L.matrix[,t] = as.numeric(1:3 == current.temp$person.obs.states[current.temp$person.obs.times == TRUE])
  } else {
    L.matrix[,t] = rep(1, 3)
  } 
  if ( any(!is.na(current.temp$num.healthy)) ) {
    n1 = current.temp$num.healthy[!is.na(current.temp$num.healthy)]
    n2 = current.temp$num.ill[!is.na(current.temp$num.ill)]
  }
  
  for (t in 2:length(unique.trans.times)) {
    
    current.time = unique.trans.times[t]; 
    current.temp = all.trans.times[all.trans.times$time == current.time,]
    
    
    if ( all(current.temp$person.obs.times == FALSE) | current.time == max(unique.trans.times) ) {
      if( all(current.temp$switch1 == 0 & current.temp$switch2 == 0) &
          (current.time != max(unique.trans.times)) ) {
        current.B = matrix(0, nrow = 3, ncol = 3)
        
        # Need to deal with death times
        temp.which.Omega = min(max(which(current.time >= bar.time)), length(bar.Omega))
        
        current.Omega = bar.Omega[temp.which.Omega]
        
        rate1.baseline = current.params[1] * zeta(n1, n2, beta = current.params[3], rho = rho) * rho
        rate2.baseline = current.params[2] * zeta(n1, n2, beta = current.params[4], rho = rho) * rho
        
        rate1.addhealthy = current.params[1] * zeta(n1+1, n2, beta = current.params[3], rho = rho) * rho
        rate2.addhealthy = current.params[2] * zeta(n1+1, n2, beta = current.params[4], rho = rho) * rho
        total.rate.addhealthy = rate1.addhealthy + rate2.addhealthy - (rate1.baseline + rate2.baseline)
        # total.rate.addhealthy = rate1.addhealthy - rate1.baseline
        
        rate1.addill = current.params[1] * zeta(n1, n2+1, beta = current.params[3], rho = rho) * rho
        rate2.addill = current.params[2] * zeta(n1, n2+1, beta = current.params[4], rho = rho) * rho
        total.rate.addill = rate1.addill + rate2.addill - (rate1.baseline + rate2.baseline)
        # total.rate.addill = rate1.addill - rate1.baseline
        
        diag(current.B) = 1 - c(total.rate.addhealthy, total.rate.addill,0)/current.Omega
        current.B[1,2] = (rate1.addhealthy - rate1.baseline)/(current.Omega)
        current.B[1,3] = (rate2.addhealthy - rate2.baseline)/(current.Omega)
        current.B[2,1] = (rate1.addill - rate1.baseline)/(current.Omega)
        current.B[2,3] = (rate2.addill - rate2.baseline)/(current.Omega)
        list.of.current.Bs[[t-1]] = current.B
      } else if (all(current.temp$failure.time != 1)) {
        ## Transition time is a prior existing one
        ## Need to put up
        temp.switch1 = current.temp$switch1[current.temp$switch1 != 0 | current.temp$switch2 != 0]
        temp.switch2 = current.temp$switch2[current.temp$switch1 != 0 | current.temp$switch2 != 0]
        
        denominator = q.split(r1 = n1 - temp.switch1[1], d1 = temp.switch1[1], 
                              r2 = n2 - temp.switch2[1], d2 = temp.switch2[1], 
                              beta = current.params[3], rho = rho)
        
        current.B = matrix(0, nrow = 3, ncol = 3); current.B[3,3] = 1
        
        numerator11 = q.split(r1 = n1 + 1 - temp.switch1[1], d1 = temp.switch1[1], 
                              r2 = n2 - temp.switch2[1], d2 = temp.switch2[1], 
                              beta = current.params[3], rho = rho)
        current.B[1,1] = numerator11/denominator; current.B[1,2] = 1 - current.B[1,1]
        
        numerator22 = q.split(r1 = n1 - temp.switch1[1], d1 = temp.switch1[1], 
                              r2 = n2 + 1 - temp.switch2[1], d2 = temp.switch2[1], 
                              beta = current.params[3], rho = rho)
        current.B[2,2] = numerator22/denominator; current.B[2,1] = 1 - current.B[2,2]
      
        list.of.current.Bs[[t-1]] = current.B
      } else if ( all(current.temp$person.obs.times == FALSE) &
                  all(current.temp$failure.time == TRUE) ) {
        ## Switch to death times
        temp.switch1 = current.temp$switch1[current.temp$switch1 != 0 | current.temp$switch2 != 0]
        temp.switch2 = current.temp$switch2[current.temp$switch1 != 0 | current.temp$switch2 != 0]
        
        denominator = q.split(r1 = n1 - temp.switch1[1], d1 = temp.switch1[1], 
                              r2 = n2 - temp.switch2[1], d2 = temp.switch2[1], 
                              beta = current.params[4], rho = rho)
        
        current.B = matrix(0, nrow = 3, ncol = 3); current.B[3,3] = 1
        
        numerator11 = q.split(r1 = n1 + 1 - temp.switch1[1], d1 = temp.switch1[1], 
                              r2 = n2 - temp.switch2[1], d2 = temp.switch2[1], 
                              beta = current.params[4], rho = rho)
        current.B[1,1] = numerator11/denominator; current.B[1,3] = 1 - current.B[1,1]
        
        numerator22 = q.split(r1 = n1 - temp.switch1[1], d1 = temp.switch1[1], 
                              r2 = n2 + 1 - temp.switch2[1], d2 = temp.switch2[1], 
                              beta = current.params[4], rho = rho)
        current.B[2,2] = numerator22/denominator; current.B[2,3] = 1 - current.B[2,2]
        list.of.current.Bs[[t-1]] = current.B
      } else {
        current.B = matrix(0, nrow = 3, ncol = 3)
        
        rate1.baseline = current.params[1] * zeta(n1, n2, beta = current.params[3], rho = rho) * rho
        rate2.baseline = current.params[2] * zeta(n1, n2, beta = current.params[4], rho = rho) * rho
        
        rate1.addhealthy = current.params[1] * zeta(n1+1, n2, beta = current.params[3], rho = rho) * rho
        rate2.addhealthy = current.params[2] * zeta(n1+1, n2, beta = current.params[4], rho = rho) * rho
        total.rate.addhealthy = rate1.addhealthy + rate2.addhealthy - (rate1.baseline + rate2.baseline)
        
        rate1.addill = current.params[1] * zeta(n1, n2+1, beta = current.params[3], rho = rho) * rho
        rate2.addill = current.params[2] * zeta(n1, n2+1, beta = current.params[4], rho = rho) * rho
        total.rate.addill = rate1.addill + rate2.addill - (rate1.baseline + rate2.baseline)
        
        current.B[1,2] = (rate1.addhealthy - rate1.baseline) / total.rate.addhealthy
        current.B[1,3] = (rate2.addhealthy - rate2.baseline) / total.rate.addhealthy
        current.B[2,1] = (rate1.addill - rate1.baseline) / total.rate.addill
        current.B[2,3] = (rate2.addill - rate2.baseline) / total.rate.addill
        current.B[3,3] = 1
        list.of.current.Bs[[t-1]] = current.B
      }
    } else {
      current.B = diag(rep(1,3))
      list.of.current.Bs[[t-1]] = current.B
    }
    
    if(any(current.temp$person.obs.times == TRUE)) {
      L.matrix[,t] = as.numeric(1:3 == current.temp$person.obs.states[current.temp$person.obs.times == TRUE])
    } else {
      L.matrix[,t] = rep(1, 3)
    } 
    
    Alpha.matrix[,t] =  (L.matrix[,t-1] * Alpha.matrix[,t-1]) %*% current.B
    
    if ( any(!is.na(current.temp$num.healthy)) ) {
      n1 = current.temp$num.healthy[!is.na(current.temp$num.healthy)]
      n2 = current.temp$num.ill[!is.na(current.temp$num.ill)]
    }
  }
  
  ## And now we backwards sample
  Beta.matrix = Alpha.matrix * 0; sampled.trajectory = vector(length = ncol(Alpha.matrix))
  T = ncol(Alpha.matrix)
  Beta.matrix[,T] = L.matrix[,T]
  # Beta.matrix[,T] = Beta.matrix[,T]/ sum(Beta.matrix[,T])
  sampled.trajectory[T] = sample(1:3, prob = Beta.matrix[,T], size = 1)
  for (t in (T-1):1) {
      B.t = list.of.current.Bs[[t]]
      Beta.matrix[,t] = Alpha.matrix[,t] * L.matrix[,t] * B.t[,sampled.trajectory[t+1]]
      Beta.matrix[,t] = Beta.matrix[,t]/ sum(Beta.matrix[,t])
      sampled.trajectory[t] = sample(1:3, prob = Beta.matrix[,t], size = 1)
      # if(sampled.trajectory[t] != sampled.trajectory[t+1]) {
      #   print( paste("Switch from state",sampled.trajectory[t+1], "to", sampled.trajectory[t], "at time", 
      #                unique.trans.times[t+1], "which is t =", t))
      # }
  }  
  
  is.nonobs.time = which(!is.element(unique.trans.times, all.trans.times$time[all.trans.times$person.obs.times]))
  is.nonobs.time = c(1,is.nonobs.time,length(unique.trans.times))
  
  
  ## NOTE: Removed obs times, but kept all others. Will remove non-trans times in next step.
  return ( 
    list( "trans.times" = unique.trans.times[is.nonobs.time], "sampled.trajectory" = sampled.trajectory[is.nonobs.time])
  )
  
}


build.new.all.person.trans.data <- function(new.trajectory, person.id, other.person.trans.data) {
  ## Build this new person.trans.times
  ## Append to other person.trans.times
  ## New all.persons.trans.data
  len.sampledtrajectory = length(new.trajectory$sampled.trajectory)
  which.trans.times = which(new.trajectory$sampled.trajectory[2:len.sampledtrajectory] != 
                              new.trajectory$sampled.trajectory[1:(len.sampledtrajectory-1)]) + 1
  which.trans.times = c(1,which.trans.times)
  
  temp = data.frame(user = person.id, time = new.trajectory$trans.times[which.trans.times], 
                    state = new.trajectory$sampled.trajectory[which.trans.times])
  
  ## Keep only true transitions and initial value at 0 if available
  return ( 
    rbind(other.person.trans.data, temp)
  )
}

sample.all.new.users <- function(all.person.obs.data, old.all.person.trans.data, 
                                 current.params, current.multiple, rho, 
                                 initial.iter = FALSE, max.attempts = 20) {
  all.users = unique(all.person.obs.data$user)
  if (initial.iter) { 
    keep.users = c(); new.all.person.trans.data = c()
    for(patient.id in all.users) {
      # print(patient.id)
      keep.users = c(patient.id, keep.users)
      new.all.person.trans.data = rbind(new.all.person.trans.data, old.all.person.trans.data[is.element(old.all.person.trans.data$user,patient.id),])
      new.all.person.obs.data = all.person.obs.data[is.element(all.person.obs.data$user,keep.users),]
      all.users = keep.users
      complete = FALSE; how.many.to.completion = 0; attempt = 1
      while(complete == FALSE) {
        how.many.to.completion = how.many.to.completion + 1
        temp = try(sample.one.new.user(patient.id, all.users, new.all.person.obs.data, 
                                       new.all.person.trans.data, current.params, 
                                       current.multiple, rho, initial.iter), TRUE)
        if(!is(temp,"try-error")) {
          new.all.person.trans.data = temp
          complete = TRUE
        } else {
          if (attempt == max.attempts) {
            complete = TRUE  
            print("Hit max attempts, so moving on.")
          } else {
            attempt = attempt + 1
            print("In a fail-safe loop, hold on one moment please. Thanks for your patience.")
          }
        }
      }
      # if(how.many.to.completion > 1 ) { print(paste("Had to run for unit", patient.id, how.many.to.completion, "many times"))}
    } 
  } else {
    new.all.person.trans.data = old.all.person.trans.data
    new.all.person.obs.data = all.person.obs.data
    for(patient.id in all.users) {
      print(patient.id)
      complete = FALSE; how.many.to.completion = 0
      max.attempts = 10; attempt = 1
      while(complete == FALSE) {
        how.many.to.completion = how.many.to.completion + 1
        temp = try(sample.one.new.user(patient.id, all.users, new.all.person.obs.data, 
                                       new.all.person.trans.data, current.params, current.multiple, rho), TRUE)
        if(!is(temp,"try-error")) {
          new.all.person.trans.data = temp
          complete = TRUE
        } else {
          if (attempt == max.attempts) {
            complete = TRUE  
            print("Hit max attempts, so moving on.")
          } else {
            attempt = attempt + 1
            print("In a fail-safe loop, hold on one moment please. Thanks for your patience.")
          }
        }
      }
    } 
  }
  return(new.all.person.trans.data)
}

sample.one.new.user <- function(patient.id, all.users, all.person.obs.data, 
                                all.person.trans.data, current.params, 
                                current.multiple, rho, initial.iter = FALSE) {
  # print(patient.id)
  temp.extract = extract.patient(all.person.obs.data, all.person.trans.data, patient.id)
  chosen.person.obs.data = temp.extract$chosen.person.obs.data
  chosen.person.trans.data = temp.extract$chosen.person.trans.data
  other.person.trans.data = temp.extract$other.person.trans.data
  
  temp.extractOmega.and.A = extract.Omega.and.A(chosen.person.obs.data, chosen.person.trans.data, other.person.trans.data, 
                                    current.params, current.multiple = current.multiple , rho) 
  
  temp.new.potential.trans.times = generate.points(temp.extractOmega.and.A$bar.Omega, temp.extractOmega.and.A$bar.Omega.time, 
                                                   temp.extractOmega.and.A$bar.A, temp.extractOmega.and.A$bar.Omega.time,
                                                   initial.iter)
  
  temp.constructinfo = construct.info(temp.new.potential.trans.times, chosen.person.trans.data, other.person.trans.data, 
                                      chosen.person.obs.data, current.params)
  
  temp.newsample <- filterfwd.backsampl(temp.constructinfo, current.params, rho, temp.extractOmega.and.A$bar.Omega, temp.extractOmega.and.A$bar.Omega.time)
  
  all.person.trans.data <- build.new.all.person.trans.data(temp.newsample, patient.id, other.person.trans.data) 
  
  return(all.person.trans.data)
}

nu.gibbs.updates <- function(data, current.params, rho, prior.params) {
  total.zeta.integral = c(0,0); switch.count = c(0,0)
  for(t in 2:nrow(data)) {
    isdeathtime = data$failure.time[t]
    n1 = data$num.healthy[t-1]; n2 = data$num.ill[t-1]
    d1 = data$switch1[t]; d2 = data$switch2[t]
    window.length = data$time[t] - data$time[t-1]
    rate1 = zeta(n1, n2, current.params[3], rho) * rho
    rate2 = zeta(n1, n2, current.params[4], rho) * rho
    current.rates = c(rate1, rate2)
    switch.count = switch.count + as.numeric(0:1 == isdeathtime)
    total.zeta.integral = total.zeta.integral + current.rates * window.length
  }
  new.alpha = prior.params$alpha + switch.count
  new.beta = prior.params$beta + total.zeta.integral
  return( c(
    rgamma(1, shape  = new.alpha[1], rate = new.beta[1]), 
    rgamma(1, shape  = new.alpha[2], rate = new.beta[2])
  ) )
}

beta.gibbs.updates <- function(data, new.nu.params, current.params, rho, prior.params) {
  temp.params = c(new.nu.params, current.params[3:4])
  proposal.beta.params = log(temp.params[3:4]) + rnorm(2, sd = prior.params$step.size)
  proposal.params = c(temp.params[1:2], exp(proposal.beta.params))
  current.llik = total.llik(data, rho)(temp.params)
  # current.prior.llik = -sum(log(dcauchy(proposal.beta.params, location = prior.params$log.mu, scale = prior.params$log.sd)))
  # proposal.prior.llik = -sum(log(dcauchy(log(temp.params[3:4]), location = prior.params$log.mu, scale = prior.params$log.sd)))
  current.prior.llik = -sum(log(dnorm(proposal.beta.params, mean = prior.params$log.mu, sd = prior.params$log.sd)))
  proposal.prior.llik = -sum(log(dnorm(log(temp.params[3:4]), mean = prior.params$log.mu, sd = prior.params$log.sd)))
  proposal.llik = total.llik(data, rho)(proposal.params)
  acceptance.rate = min(1, exp(current.llik + current.prior.llik - proposal.llik - proposal.prior.llik)  )
  accept.proposal = as.logical(rbinom(n = 1, size = 1, prob = acceptance.rate) ==1)
  new.beta.params = accept.proposal * proposal.params + (1-accept.proposal) * temp.params
  return(
    list( "new.params" = new.beta.params, "accept" = accept.proposal, "accept.rate" = acceptance.rate)
  )
}

Survival.fit <- function(current.global, current.params, grid.length = 0.01) {
  
  unique.trans.times = unique(current.global$time)
  Survival.inithealthy = c(1, 0, 0); Survival.initill = c(0, 1, 0)
  computed.at.times = c(0) 
  Survival.inithealthy.matrix = as.matrix(Survival.inithealthy, ncol = 1)
  Survival.initill.matrix = as.matrix(Survival.initill, ncol = 1)
  
  
  
  for (t in 2:length(unique.trans.times)) {
    
    n1 = current.global$num.healthy[t-1]; n2 = current.global$num.ill[t-1]
    rate1.baseline = current.params[1] * zeta(n1, n2, beta = current.params[3], rho = rho) * rho
    rate2.baseline = current.params[2] * zeta(n1, n2, beta = current.params[4], rho = rho) * rho
    
    rate1.addhealthy = current.params[1] * zeta(n1+1, n2, beta = current.params[3], rho = rho) * rho
    rate2.addhealthy = current.params[2] * zeta(n1+1, n2, beta = current.params[4], rho = rho) * rho
    total.rate.addhealthy = rate1.addhealthy + rate2.addhealthy - (rate1.baseline + rate2.baseline)
    current.Q = matrix(nrow = 3, ncol = 3)
    current.Q[1,1] = - total.rate.addhealthy
    current.Q[1,2] = (rate1.addhealthy - rate1.baseline) 
    current.Q[1,3] = (rate2.addhealthy - rate2.baseline)
    
    rate1.addill = current.params[1] * zeta(n1, n2+1, beta = current.params[3], rho = rho) * rho
    rate2.addill = current.params[2] * zeta(n1, n2+1, beta = current.params[4], rho = rho) * rho
    total.rate.addill = rate1.addill + rate2.addill - (rate1.baseline + rate2.baseline)
    current.Q[2,2] = - total.rate.addill
    current.Q[2,1] = (rate1.addill - rate1.baseline) 
    current.Q[2,3] = (rate2.addill - rate2.baseline)
    current.Q[3,1:3] = 0
    
    # gap = current.global$time[t] - current.global$time[t-1]
    
    ## Expand on a 0.01 grid 
    ## To compute survival at these time
    temp.times = c(seq(current.global$time[t-1], current.global$time[t], grid.length), current.global$time[t])
    diff.gaps = diff(temp.times)
    computed.at.times = c(computed.at.times, temp.times[2:length(temp.times)])
    for (i in 1:length(diff.gaps)) {
      current.B = as.matrix(expm(current.Q*diff.gaps[i]))
      Survival.inithealthy = Survival.inithealthy%*%current.B
      Survival.inithealthy.matrix = cbind(Survival.inithealthy.matrix, t(Survival.inithealthy))
      Survival.initill = Survival.initill%*%current.B
      Survival.initill.matrix = cbind(Survival.initill.matrix, t(Survival.initill))
    }
    
    if (current.global$failure.time[t] != 1) { 
      denominator = q.split(r1 = n1 - current.global$switch1[t], d1 = current.global$switch1[t], 
                            r2 = n2 - current.global$switch2[t], d2 = current.global$switch2[t], 
                            beta = current.params[3], rho = rho)
      
      current.B = matrix(0, nrow = 3, ncol = 3); current.B[3,3] = 1
      
      numerator11 = q.split(r1 = n1 + 1 - current.global$switch1[t], d1 = current.global$switch1[t], 
                            r2 = n2 - current.global$switch2[t], d2 = current.global$switch2[t], 
                            beta = current.params[3], rho = rho)
      current.B[1,1] = numerator11/denominator; current.B[1,2] = 1 - current.B[1,1]
      
      numerator22 = q.split(r1 = n1 - current.global$switch1[t], d1 = current.global$switch1[t], 
                            r2 = n2 + 1 - current.global$switch2[t], d2 = current.global$switch2[t], 
                            beta = current.params[3], rho = rho)
      current.B[2,2] = numerator22/denominator; current.B[2,1] = 1 - current.B[2,2]
      
      Survival.inithealthy = Survival.inithealthy%*%current.B
      Survival.inithealthy.matrix = cbind(Survival.inithealthy.matrix, t(Survival.inithealthy))
      Survival.initill = Survival.initill%*%current.B
      Survival.initill.matrix = cbind(Survival.initill.matrix, t(Survival.initill))
      computed.at.times = c(computed.at.times, current.global$time[t])
      
    } else {
      denominator = q.split(r1 = n1 - current.global$switch1[t], d1 = current.global$switch1[t], 
                            r2 = n2 - current.global$switch2[t], d2 = current.global$switch2[t], 
                            beta = current.params[4], rho = rho)
      
      current.B = matrix(0, nrow = 3, ncol = 3); current.B[3,3] = 1
      
      numerator11 = q.split(r1 = n1 + 1 - current.global$switch1[t], d1 = current.global$switch1[t], 
                            r2 = n2 - current.global$switch2[t], d2 = current.global$switch2[t], 
                            beta = current.params[4], rho = rho)
      current.B[1,1] = numerator11/denominator; current.B[1,3] = 1 - current.B[1,1]
      
      numerator22 = q.split(r1 = n1 - current.global$switch1[t], d1 = current.global$switch1[t], 
                            r2 = n2 + 1 - current.global$switch2[t], d2 = current.global$switch2[t], 
                            beta = current.params[4], rho = rho)
      current.B[2,2] = numerator22/denominator; current.B[2,3] = 1 - current.B[2,2]
      
      Survival.inithealthy = Survival.inithealthy%*%current.B
      Survival.inithealthy.matrix = cbind(Survival.inithealthy.matrix, t(Survival.inithealthy))
      Survival.initill = Survival.initill%*%current.B
      Survival.initill.matrix = cbind(Survival.initill.matrix, t(Survival.initill))
      computed.at.times = c(computed.at.times, current.global$time[t])
      
    }
    
  }
  
  return(cbind(
    computed.at.times,
    Survival.inithealthy.matrix[3,],
    Survival.initill.matrix[3,]
  ))
  
}