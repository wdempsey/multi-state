compute.pdf <- function(trans.matrix,beta.vector, rho) {
  # Trans.matrix (diag is r_11, ...., r_kk; off-diagonal are those who transition))
  f2 <- function(x) {
    log.total = (rho-1)* log(x) - log(1-x)
    for (i in 1:nrow(trans.matrix)) {
      log.total = log.total + sum(trans.matrix[i,-i]) * log(1-x^beta.vector[i]) + trans.matrix[i,i]*beta.vector[i]*log(x)
    }
    return( exp(log.total) )
  }
  return(f2)
}

q.split <- function(trans.matrix, beta.vector, alpha.matrix, rho, tol = 0.000001) {
  temp = compute.pdf(trans.matrix,beta.vector, rho)
  temp.integrate = integrate(temp,0.00,tol)$value + integrate(temp,tol,1)$value
  log.temp.total = 0
  for (i in 1:nrow(trans.matrix)) {
    log.multi.coef = lfactorial(sum(trans.matrix[i,-i])) - sum(lfactorial(trans.matrix[i,-i]))  
    log.temp.total = log.temp.total + log.multi.coef + sum(trans.matrix[i,-i] * log(alpha.matrix[i,-i]), na.rm = TRUE)
  }
  return(exp(log.temp.total)*temp.integrate)
}

compute.normalizing <- function(n.vector, beta.vector, rho) {
  f2 <- function(x) {
    log.total = (rho-1)* log(x) - log(1-x)
    log.total = log.total + log(1-x^(sum(beta.vector*n.vector)))
    return( exp(log.total) )
  }
  return(f2)
}

zeta <- function(n.vector, beta.vector, rho) {
  temp = compute.normalizing(n.vector, beta.vector, rho)
  return(integrate(temp,0,1)$value)
}

zeta.mapply <- function(beta.vector, rho) {
  interior.zeta <- function(num.nocav, num.mild, num.severe) {
    n.vector = c(num.nocav, num.mild, num.severe)
    temp = compute.normalizing(n.vector, beta.vector, rho)
    return(integrate(temp,0,1)$value)
  }
  return(interior.zeta)
}

lik.component <- function(params, trans.matrix, isdeathtime, window.time, rho) {
  n.vector = rowSums(trans.matrix)
  # Alternative just for mapply functionality
  nu = params[1:2]; beta.alive = c(1,params[3:4]); beta.dead = c(1,params[5:6]); 
  alpha.alive = params[7]
  # Setup up alive matrix
  alpha.alive.matrix = matrix(0, nrow = 3, ncol = 3)
  alpha.alive.matrix[1,2] = 1; alpha.alive.matrix[2,1] = alpha.alive
  alpha.alive.matrix[2,3] = (1-alpha.alive); alpha.alive.matrix[3,2] = 1
  diag(alpha.alive.matrix) = 1
  # Set up death time matrix
  alpha.death.matrix = matrix(0,nrow = 3, ncol = 4) 
  alpha.death.matrix[1:3,4] = 1; diag(alpha.death.matrix) = 1
  
  rate.alive = nu[1] * zeta(n.vector, beta.alive, rho) * rho
  rate.dead  = nu[2] * zeta(n.vector, beta.dead, rho) * rho
  total.rate =  rate.alive + rate.dead
  if(isdeathtime) {
    nu = params[2]*rho; beta = beta.dead; alpha = alpha.death.matrix
  } else {
    nu = params[1]*rho; beta = beta.alive; alpha = alpha.alive.matrix
  }
  
  if( all(trans.matrix[row(trans.matrix)!=col(trans.matrix)] == 0) ) {
    return(
      log(nu) - total.rate*window.time
    )
  } else {
    return(
      log(nu) - total.rate*window.time + log(q.split(trans.matrix, beta, alpha, rho))
    )
  }
}

total.llik <- function(data, rho) {
  subfunction <- function(params) {
    total = 0
    for(t in 2:nrow(data)) {
      isdeathtime = data$failure.time[t]
      num.nocav = data$num.nocav[t-1]; num.mild = data$num.mild[t-1]
      num.severe = data$num.severe[t-1]
      d1 = data$switch1[t]; d2a = data$switch2a[t]
      d2b = data$switch2b[t]; d3 = data$switch3[t]
      r.nocav = num.nocav - d1; r.mild = num.mild - d2a - d2b;
      r.severe = num.severe - d3;
      if(r.mild < 0) {r.mild = 0}
      if(!isdeathtime) {
        trans.matrix = matrix(0, nrow = 3, ncol = 3)
        diag(trans.matrix) = c(r.nocav, r.mild, r.severe);
        trans.matrix[1,2] = d1; trans.matrix[2,1] = d2a; 
        trans.matrix[2,3] = d2b; trans.matrix[3,2] = d3
      } else {
        trans.matrix = matrix(0, nrow = 3, ncol = 4)
        diag(trans.matrix) = c(r.nocav, r.mild, r.severe);
        trans.matrix[1,4] = d1; trans.matrix[2,4] = d2a; 
        trans.matrix[3,4] = d3
      }
      window.length = data$time[t] - data$time[t-1]
      total = total + lik.component(params, trans.matrix, isdeathtime, window.length, rho)
    }
    return(-total)
  }
  return(subfunction)
}

convert.to.global <- function(person.data, censoring.data) {
  ## Convert data to global configuration
  ## Assume each participant starts with initial obs
  ## At t = 0
  current.time = 0.0
  states = person.data$state[person.data$time == current.time]
  num.nocav = sum(states == 1); num.mild = sum(states == 2)
  num.severe = sum(states == 3); num.dead = sum(states == 4)
  switch1 = 0; switch2a = switch2b = 0; switch3 = 0; failure.time = 0
  init = matrix(c(current.time, num.nocav, num.mild, num.severe, num.dead, 
                  switch1, switch2a, switch2b, switch3, failure.time),
                nrow = 1)
  simdata = data.frame(init)
  ## I will use switch1, switch3 cause only one type per failure block
  ## I will use switch2a for 2-> 1 or 2-> 4
  ## I will use switch2b for 2-> 3.
  names(simdata) = c("time", "num.nocav", "num.mild", "num.severe", 
                     "num.dead", "switch1", "switch2a", 
                     "switch2b", "switch3", "failure.time")
  if(nrow(person.data) == 0) {max.time = 0} else {max.time = max(person.data$time, censoring.data$time)}
  total.censored = 0
  
  while (current.time < max.time) {
    next.time = min(person.data$time[person.data$time > current.time], 
                    censoring.data$time[censoring.data$time > current.time])
    trans.users.list = person.data$user[person.data$time == next.time]
    next.states = person.data$state[person.data$time == next.time]
    
    censoring.users.list = censoring.data$user[censoring.data$time == next.time & 
                                                 censoring.data$status == 0]
    
    ## Grab all data prior to the next.time for the users who are observed
    ## at that time. This is so we can compute how they transitioned!
    temp = person.data[is.element(person.data$user,trans.users.list) & person.data$time < next.time,]
    switch1 = 0; switch2a = 0; switch2b = 0; switch3 = 0; failure.time = as.numeric(any(next.states==4))
    
    for (user in trans.users.list) {
      ## Iterate over all users in cohort with data at next.time
      ## Either all failures or all transitions at time t
      temp.time = temp$time[temp$user == user]
      temp.state = temp$state[temp$user == user]
      previous.state = temp.state[temp.time == max(temp.time)]
      next.state = next.states[trans.users.list == user]
      user.data = person.data[person.data$user == user,]
      if(previous.state == 1 & is.element(next.state, c(2,4)) ) {
        ## Switch out of no Cav 
        switch1 = switch1 + 1
      } else if (previous.state == 2 & is.element(next.state, c(1,4)) ) {
        ## Switch out of mild to no cav or death
        switch2a = switch2a + 1
      } else if (previous.state == 2 & next.state == 3) {
        ## Switch out of mild to severe
        switch2b = switch2b + 1
      } else if (previous.state == 3 & is.element(next.state, c(2,4)) ) {
        ## Switch out of severe
        switch3 = switch3 + 1
      } else if (previous.state == 1 & next.state == 3) {
        ## Account for a double jump from no cav to severe
        switch1 = switch1 + 1
        switch2b = switch2b + 1
      } else if (previous.state == 3 & next.state == 1) {
        ## Account for a double jump from severe to no cav
        switch3 = switch3 + 1
        switch2a = switch2a + 1
      }
    }
        
    if (failure.time == 1) {
      num.nocav = num.nocav - switch1; num.mild = num.mild - switch2a - switch2b; 
      num.severe = num.severe - switch3; 
      num.dead = num.dead + switch1 + switch2a + switch2b + switch3
    } else {
      num.nocav = num.nocav - switch1 + switch2a; 
      num.mild = num.mild - switch2a - switch2b + switch1 + switch3; 
      num.severe = num.severe - switch3 + switch2b
      num.dead = num.dead
    }
    
    for (user in censoring.users.list) {
      temp.data = person.data[person.data$user == user,]
      final.state = temp.data$state[temp.data$time == max(temp.data$time)]
      if (final.state == 1) {
        num.nocav = num.nocav - 1; total.censored = total.censored + 1
      } else if (final.state == 2) {
        num.mild = num.mild - 1; total.censored = total.censored + 1
      } else if (final.state == 3) {
        num.severe = num.severe - 1; total.censored = total.censored + 1
      }
    }
    # if(length(censoring.users.list) > 0 ) {print("This was a censoring time")}
    # print(paste("total censored =",total.censored))
    # print(paste("total dead =",num.dead))
    simdata = rbind(simdata, c(next.time, num.nocav, num.mild, 
                               num.severe, num.dead, switch1, switch2a,
                               switch2b, switch3, failure.time)
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

extract.Omega.and.A <- function(chosen.person.obs.data, chosen.person.trans.data, other.person.trans.data, 
                                censoring.data, current.params, current.multiple = 5, rho = 10) {
  ## Generate the vector of Omega's for each window of time
  ## This changes at each transition time
  other.global.trans.data = convert.to.global(other.person.trans.data, censoring.data)
  max.time = max(chosen.person.obs.data$time)
  other.global.trans.data.shortened = other.global.trans.data[other.global.trans.data$time < max.time, ]
  
  # Alternative just for mapply functionality
  nu = current.params[1:2]; beta.alive = c(1,current.params[3:4]); beta.dead = c(1,current.params[5:6]); 
  alpha.alive = current.params[7]
  # Setup up alive matrix
  alpha.alive.matrix = matrix(0, nrow = 3, ncol = 3)
  alpha.alive.matrix[1,2] = 1; alpha.alive.matrix[2,1] = alpha.alive
  alpha.alive.matrix[2,3] = (1-alpha.alive); alpha.alive.matrix[3,2] = 1
  diag(alpha.alive.matrix) = 1
  # Set up death time matrix
  alpha.death.matrix = matrix(0,nrow = 3, ncol = 4) 
  alpha.death.matrix[1:3,4] = 1; diag(alpha.death.matrix) = 1
  
  ## Construct rate functions 
  rate.alive.addhealthy = nu[1] * mapply(zeta.mapply(beta.alive, rho), num.nocav = other.global.trans.data.shortened$num.nocav+1, 
                                         num.mild = other.global.trans.data.shortened$num.mild,
                                         num.severe = other.global.trans.data.shortened$num.severe) * rho
  
  rate.dead.addhealthy = nu[2] * mapply(zeta.mapply(beta.dead,rho), num.nocav = other.global.trans.data.shortened$num.nocav+1, 
                                        num.mild = other.global.trans.data.shortened$num.mild,
                                        num.severe = other.global.trans.data.shortened$num.severe) * rho
  
  rate.alive.baseline = nu[1] * mapply(zeta.mapply(beta.alive, rho), num.nocav = other.global.trans.data.shortened$num.nocav, 
                                       num.mild = other.global.trans.data.shortened$num.mild,
                                       num.severe = other.global.trans.data.shortened$num.severe) * rho
  
  rate.dead.baseline = nu[2] * mapply(zeta.mapply(beta.dead, rho), num.nocav = other.global.trans.data.shortened$num.nocav, 
                                      num.mild = other.global.trans.data.shortened$num.mild,
                                      num.severe = other.global.trans.data.shortened$num.severe) * rho
  
  total.rate.addhealthy = rate.alive.addhealthy + rate.dead.addhealthy - rate.alive.baseline - rate.dead.baseline
  
  rate.alive.addmild = nu[1] * mapply(zeta.mapply(beta.alive, rho), num.nocav = other.global.trans.data.shortened$num.nocav, 
                                         num.mild = other.global.trans.data.shortened$num.mild+1,
                                         num.severe = other.global.trans.data.shortened$num.severe) * rho
  
  rate.dead.addmild = nu[2] * mapply(zeta.mapply(beta.dead, rho), num.nocav = other.global.trans.data.shortened$num.nocav, 
                                        num.mild = other.global.trans.data.shortened$num.mild+1,
                                        num.severe = other.global.trans.data.shortened$num.severe) * rho
  
  total.rate.addmild = rate.alive.addmild + rate.dead.addmild - rate.alive.baseline - rate.dead.baseline
  
  rate.alive.addsevere = nu[1] * mapply(zeta.mapply(beta.alive, rho), num.nocav = other.global.trans.data.shortened$num.nocav, 
                                      num.mild = other.global.trans.data.shortened$num.mild,
                                      num.severe = other.global.trans.data.shortened$num.severe+1) * rho
  
  rate.dead.addsevere = nu[2] * mapply(zeta.mapply(beta.dead, rho), num.nocav = other.global.trans.data.shortened$num.nocav, 
                                     num.mild = other.global.trans.data.shortened$num.mild,
                                     num.severe = other.global.trans.data.shortened$num.severe+1) * rho
  
  total.rate.addsevere = rate.alive.addsevere + rate.dead.addsevere - rate.alive.baseline - rate.dead.baseline
  
  all.rates = cbind(total.rate.addhealthy, total.rate.addmild, total.rate.addsevere)
  bar.Omega = current.multiple * apply(all.rates,1,max)
  bar.Omega.time = other.global.trans.data.shortened$time
  bar.Omega.time = c(bar.Omega.time, max.time)
  
  ## NOW DO FOR ALL PEOPLE TO CONSTRUCT A.t
  all.global.trans.data = convert.to.global(rbind(other.person.trans.data, chosen.person.trans.data),
                                            censoring.data)
  all.global.trans.data.shortened = all.global.trans.data[all.global.trans.data$time < max.time, ]
  
  ## Construct rate functions 
  rate1 = nu[1] * mapply(zeta.mapply(beta.alive, rho), num.nocav = all.global.trans.data.shortened$num.nocav, 
                  num.mild = all.global.trans.data.shortened$num.mild,
                  num.severe = all.global.trans.data.shortened$num.severe) * rho
  
  rate2 = nu[2] * mapply(zeta.mapply(beta.dead, rho), num.nocav = all.global.trans.data.shortened$num.nocav, 
                         num.mild = all.global.trans.data.shortened$num.mild,
                         num.severe = all.global.trans.data.shortened$num.severe) * rho
  
  ## Remove chosen person from all.global.trans.data.shortened
  temp.trans.data.shortened = all.global.trans.data.shortened
  for (i in 1:nrow(temp.trans.data.shortened)) {
    temp.which.obs = max(which(chosen.person.trans.data$time <= temp.trans.data.shortened$time[i]))
    temp.state = chosen.person.trans.data$state[temp.which.obs]
    if( temp.state == 1 ) {
      temp.trans.data.shortened$num.nocav[i] = temp.trans.data.shortened$num.nocav[i] - 1
    } else if (temp.state == 2) {
      temp.trans.data.shortened$num.mild[i] = temp.trans.data.shortened$num.mild[i] - 1
    } else if (temp.state == 3) {
      temp.trans.data.shortened$num.severe[i] = temp.trans.data.shortened$num.severe[i] - 1
    }
  }
  
  rate1.At.baseline = nu[1] * mapply(zeta.mapply(beta.alive, rho), num.nocav = temp.trans.data.shortened$num.nocav, 
                                     num.mild = temp.trans.data.shortened$num.mild,
                                     num.severe = temp.trans.data.shortened$num.severe) * rho
  
  rate2.At.baseline = nu[2] * mapply(zeta.mapply(beta.dead, rho), num.nocav = temp.trans.data.shortened$num.nocav, 
                                     num.mild = temp.trans.data.shortened$num.mild,
                                     num.severe = temp.trans.data.shortened$num.severe) * rho
  
  bar.A = rate1 + rate2 - rate1.At.baseline - rate2.At.baseline
  bar.A.time = all.global.trans.data.shortened$time
  bar.A.time = c(bar.A.time, max.time)

  return(list("bar.Omega" = bar.Omega, "bar.Omega.time" = bar.Omega.time, 
              "bar.A" = bar.A, "bar.A.time" = bar.A.time))
  
}

generate.points <- function(bar.Omega, bar.Omega.time, bar.A, bar.A.time, max.time, initial.iter) {
  ## Generate the vector of potential "new" transition times
  all.times = unique(c(bar.Omega.time, bar.A.time, max.time))
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

construct.info <- function(new.potential.trans.times, chosen.person.trans.data, other.person.trans.data, 
                           chosen.person.obs.data, censoring.data, current.params) {
  ## Combine all the information necessary to perform 
  ## the filter forward, backward sampling algorithm.
  ## These are all trans times, obs times for patient,
  ## and all new potential trans times
  other.global.trans.data = convert.to.global(other.person.trans.data,censoring.data)
  max.time = max(chosen.person.obs.data$time)
  other.global.trans.data.shortened = other.global.trans.data[other.global.trans.data$time <= max.time, ]
  ## Add in all the other ppl trans times
  all.trans.times = other.global.trans.data.shortened[,c("time", "num.nocav", "num.mild", "num.severe", "num.dead",
                                                         "switch1", "switch2a", "switch2b", "switch3", "failure.time")]
  all.trans.times$person.obs.times = FALSE
  all.trans.times$person.obs.states = NA
  all.trans.times$new.trans.time = FALSE
  
  ## Add in prior runs trans times for current person
  which.person.time.include = !is.element(chosen.person.trans.data$time, all.trans.times$time)
  temp.persontrans = chosen.person.trans.data[which.person.time.include,]
  temp.failure.time = as.numeric(temp.persontrans$state == 4)
  temp.persontrans = data.frame(time = temp.persontrans$time)
  if(nrow(temp.persontrans) != 0) {
    temp.persontrans$num.nocav = temp.persontrans$num.mild = temp.persontrans$num.severe = temp.persontrans$num.dead = NA
    temp.persontrans$switch1 = 0; temp.persontrans$switch2a = 0
    temp.persontrans$switch2b = 0; temp.persontrans$switch3 = 0; temp.persontrans$failure.time = 0
    temp.persontrans$person.obs.times = FALSE
    temp.persontrans$person.obs.states = NA
    temp.persontrans$new.trans.time = FALSE
    all.trans.times = rbind(all.trans.times,temp.persontrans)
  }    
  ## Add in potential new trans times
  if(length(new.potential.trans.times) != 0) {
    temp.newtrans = data.frame(time = new.potential.trans.times)
    temp.newtrans$num.nocav = temp.newtrans$num.mild = temp.newtrans$num.severe = temp.newtrans$num.dead = NA
    temp.newtrans$switch1 = temp.newtrans$switch2a = temp.newtrans$switch2b = temp.newtrans$switch3 = 0
    temp.newtrans$failure.time = FALSE
    temp.newtrans$person.obs.times = FALSE
    temp.newtrans$person.obs.states = NA
    temp.newtrans$new.trans.time = TRUE
    all.trans.times = rbind(all.trans.times, temp.newtrans)
  }  
  
  ## Add in obs times for current person
  temp.personobs = data.frame(time = chosen.person.obs.data$time, person.obs.states = chosen.person.obs.data$state)
  temp.personobs$num.nocav = temp.personobs$num.mild = temp.personobs$num.severe = temp.personobs$num.dead = NA
  temp.personobs$switch1 = temp.personobs$switch2a = temp.personobs$switch2b = temp.personobs$switch3 = 0
  temp.personobs$person.obs.times = TRUE; temp.personobs$new.trans.time = FALSE
  temp.personobs$failure.time = as.numeric(temp.personobs$person.obs.states == 4)
  all.trans.times = rbind(all.trans.times, temp.personobs)
  
  # Correct for temporal ordering
  all.trans.times = all.trans.times[order(all.trans.times$time),]
  
  # Remove non-obs columns with 0s in all switches
  remove.obs = all.trans.times$switch1 == 0 & all.trans.times$switch2a == 0 & 
    all.trans.times$switch2b == 0 & all.trans.times$switch3 == 0 &
    all.trans.times$person.obs.times == FALSE & all.trans.times$time > 0 & all.trans.times$new.trans.time == FALSE

  # NOTE: need to be careful around duplicate times
  # Kept for now.  Will check if necessary later
  return(all.trans.times[!remove.obs,])
}

filterfwd.backsampl <- function(all.trans.times, current.params, rho, bar.Omega, bar.time) {
  ### Alpha is a matrix 3 rows (for each state), with each column summing to 1
  ### Each column corresponds to a potential trans.time
  unique.trans.times = unique(all.trans.times$time)
  Alpha.matrix = matrix(nrow = 4, ncol = length(unique.trans.times))
  t = 1; current.time = unique.trans.times[t]; current.temp = all.trans.times[all.trans.times$time == current.time,]
  Alpha.matrix[,t] = c(1, 0, 0, 0)
  list.of.current.Bs = list()
  L.matrix = matrix(nrow = 4, ncol = length(unique.trans.times))
  
  ## Construct beta.matrices
  nu = params[1:2]; beta.alive = c(1,params[3:4]); beta.dead = c(1,params[5:6]); 
  alpha.alive = params[7]
  # Setup up alive matrix
  alpha.alive.matrix = matrix(0, nrow = 3, ncol = 3)
  alpha.alive.matrix[1,2] = 1; alpha.alive.matrix[2,1] = alpha.alive
  alpha.alive.matrix[2,3] = 1-alpha.alive; alpha.alive.matrix[3,2] = 1
  diag(alpha.alive.matrix) = 1
  # Set up death time matrix
  alpha.death.matrix = matrix(0,nrow = 3, ncol = 4) 
  alpha.death.matrix[1:3,4] = 1; diag(alpha.death.matrix) = 1
  
  if(any(current.temp$person.obs.times == TRUE)) {
    L.matrix[,t] = as.numeric(1:4 == current.temp$person.obs.states[current.temp$person.obs.times == TRUE])
  } else {
    L.matrix[,t] = c(rep(1, 3),0)
  } 
  if ( any(!is.na(current.temp$num.nocav)) ) {
    num.nocav = current.temp$num.nocav[!is.na(current.temp$num.nocav)]
    num.mild = current.temp$num.mild[!is.na(current.temp$num.mild)]
    num.severe = current.temp$num.severe[!is.na(current.temp$num.severe)]
    num.dead = current.temp$num.dead[!is.na(current.temp$num.dead)]
  }
  
  for (t in 2:length(unique.trans.times)) {
    
    current.time = unique.trans.times[t]; 
    current.temp = all.trans.times[all.trans.times$time == current.time,]
    
    if ( all(current.temp$person.obs.times == FALSE) | 
         (current.time == max(unique.trans.times) & 
          any(current.temp$person.obs.states == 4, na.rm = TRUE) ) 
         ) {
      if( all(current.temp$switch1 == 0 & current.temp$switch2a == 0 & 
              current.temp$switch2b == 0 & current.temp$switch3 == 0 ) &
          (current.time != max(unique.trans.times)) ) {
        current.B = matrix(0, nrow = 4, ncol = 4)
        
        # Need to deal with death times
        temp.which.Omega = min(max(which(current.time >= bar.time)), length(bar.Omega))
        
        current.Omega = bar.Omega[temp.which.Omega]
        n.vector = c(num.nocav, num.mild, num.severe)
        
        rate.alive.baseline = nu[1] * zeta(n.vector, beta.alive, rho = rho) * rho
        rate.death.baseline = nu[2] * zeta(n.vector, beta.dead, rho = rho) * rho
        
        n.vector.addnocav = n.vector; n.vector.addnocav[1] = n.vector.addnocav[1] + 1
        rate.alive.addnocav = nu[1] * zeta(n.vector.addnocav, beta.alive, rho = rho) * rho
        rate.death.addnocav = nu[2] * zeta(n.vector.addnocav, beta.dead, rho = rho) * rho
        total.rate.addnocav = rate.alive.addnocav + rate.death.addnocav - (rate.alive.baseline + rate.death.baseline)
        # total.rate.addhealthy = rate1.addhealthy - rate1.baseline
        
        n.vector.addmild = n.vector; n.vector.addmild[2] = n.vector.addmild[2] + 1
        rate.alive.addmild = nu[1] * zeta(n.vector.addmild, beta.alive, rho = rho) * rho
        rate.death.addmild = nu[2] * zeta(n.vector.addmild, beta.dead, rho = rho) * rho
        total.rate.addmild = rate.alive.addmild + rate.death.addmild - (rate.alive.baseline + rate.death.baseline)
        # total.rate.addhealthy = rate1.addhealthy - rate1.baseline
        trans.matrix.mild = matrix(0,nrow = 3, ncol = 3)
        diag(trans.matrix.mild) = n.vector; trans.matrix.mild[2,1] = 1
        temp.mildtonocav = q.split(trans.matrix.mild, beta.alive, alpha.alive.matrix, rho)
        trans.matrix.mild[2,1] = 0; trans.matrix.mild[2,3] = 1
        temp.mildtosevere = q.split(trans.matrix.mild, beta.alive, alpha.alive.matrix, rho)
        frac.mildtonocav = temp.mildtonocav/ (temp.mildtonocav + temp.mildtosevere)
        
        n.vector.addsevere = n.vector; n.vector.addsevere[3] = n.vector.addsevere[3] + 1
        rate.alive.addsevere = nu[1] * zeta(n.vector.addsevere, beta.alive, rho = rho) * rho
        rate.death.addsevere = nu[2] * zeta(n.vector.addsevere, beta.dead, rho = rho) * rho
        total.rate.addsevere = rate.alive.addsevere + rate.death.addsevere - (rate.alive.baseline + rate.death.baseline)
        # total.rate.addhealthy = rate1.addhealthy - rate1.baseline
        
        diag(current.B) = 1 - c(total.rate.addnocav, total.rate.addmild, total.rate.addsevere,0)/current.Omega
        current.B[1,2] = (rate.alive.addnocav - rate.alive.baseline)/(current.Omega)
        current.B[1,4] = (rate.death.addnocav - rate.death.baseline)/(current.Omega)
        current.B[2,1] =  frac.mildtonocav * (rate.alive.addmild - rate.alive.baseline)/(current.Omega)
        current.B[2,3] = (1-frac.mildtonocav) * (rate.alive.addmild - rate.alive.baseline)/(current.Omega)
        current.B[2,4] = (rate.death.addsevere - rate.death.baseline)/(current.Omega)
        current.B[3,2] = (rate.alive.addsevere - rate.alive.baseline)/(current.Omega)
        current.B[3,4] = (rate.death.addsevere - rate.death.baseline)/(current.Omega)
        list.of.current.Bs[[t-1]] = current.B
      } else if (all(current.temp$failure.time != 1)) {
        ## Transition time is a prior existing one
        ## Need to put up
        temp.obs = which(current.temp$switch1 != 0 | current.temp$switch2a != 0 |
                           current.temp$switch2b != 0 | current.temp$switch3 != 0)
        temp.switch1 = current.temp$switch1[temp.obs]; temp.switch2a = current.temp$switch2a[temp.obs]
        temp.switch2b = current.temp$switch2b[temp.obs]; temp.switch3 = current.temp$switch3[temp.obs]
        
        ## Construct trans.matrix
        r.nocav = num.nocav - temp.switch1; r.mild = num.mild - temp.switch2a - temp.switch2b;
        r.severe = num.severe - temp.switch3;
        if(r.mild < 0) {r.mild = 0}
        trans.matrix = matrix(0, nrow = 3, ncol = 3)
        diag(trans.matrix) = c(r.nocav, r.mild, r.severe);
        trans.matrix[1,2] = temp.switch1; trans.matrix[2,1] = temp.switch2a; 
        trans.matrix[2,3] = temp.switch2b; trans.matrix[3,2] = temp.switch3
        
        denominator = q.split(trans.matrix, beta.alive, alpha.alive.matrix, rho)
        
        current.B = matrix(0, nrow = 4, ncol = 4); current.B[4,4] = 1
        
        trans.matrix.addnocav = trans.matrix; trans.matrix.addnocav[1,2] = trans.matrix.addnocav[1,2] + 1
        numerator12 = q.split(trans.matrix.addnocav, beta.alive, alpha.alive.matrix, rho)
        
        current.B[1,2] = numerator12/denominator; current.B[1,1] = 1 - current.B[1,2]
        
        trans.matrix.addmild.2 = trans.matrix.addmild.2a = trans.matrix.addmild.2b = trans.matrix; 
        trans.matrix.addmild.2[2,2] = trans.matrix.addmild.2[2,2] + 1
        trans.matrix.addmild.2a[2,1] = trans.matrix.addmild.2a[2,1] + 1
        trans.matrix.addmild.2b[2,3] = trans.matrix.addmild.2b[2,3] + 1
        numerator22 = q.split(trans.matrix.addmild.2, beta.alive, alpha.alive.matrix, rho)
        temp.input21 = q.split(trans.matrix.addmild.2a, beta.alive, alpha.alive.matrix, rho)
        temp.input23 = q.split(trans.matrix.addmild.2b, beta.alive, alpha.alive.matrix, rho)
        
        current.B[2,1] = temp.input21/denominator
        current.B[2,3] = temp.input23/denominator
        current.B[2,2] = 1-current.B[2,1] - current.B[2,3]
        
        trans.matrix.addsevere = trans.matrix; trans.matrix.addsevere[3,2] = trans.matrix.addsevere[3,2] + 1
        numerator32 = q.split(trans.matrix.addsevere, beta.alive, alpha.alive.matrix, rho)
        
        current.B[3,2] = numerator32/denominator; current.B[3,3] = 1 - current.B[3,2]
        
        list.of.current.Bs[[t-1]] = current.B
      } else if ( all(current.temp$person.obs.times == FALSE) &
                  all(current.temp$failure.time == TRUE) ) {
        ## Switch to death times
        temp.obs = which(current.temp$switch1 != 0 | current.temp$switch2a != 0 |
                           current.temp$switch2b != 0 | current.temp$switch3 != 0)
        temp.switch1 = current.temp$switch1[temp.obs]; temp.switch2a = current.temp$switch2a[temp.obs]
        temp.switch2b = current.temp$switch2b[temp.obs]; temp.switch3 = current.temp$switch3[temp.obs]
        
        ## Construct trans.matrix
        r.nocav = num.nocav - temp.switch1; r.mild = num.mild - temp.switch2a - temp.switch2b;
        r.severe = num.severe - temp.switch3;
        trans.matrix = matrix(0, nrow = 3, ncol = 4)
        diag(trans.matrix) = c(r.nocav, r.mild, r.severe);
        trans.matrix[1,4] = temp.switch1; trans.matrix[2,4] = temp.switch2a; trans.matrix[3,4] = temp.switch3
        
        denominator = q.split(trans.matrix, beta.dead, alpha.death.matrix, rho)
        
        current.B = matrix(0, nrow = 4, ncol = 4); current.B[4,4] = 1
        
        trans.matrix.addnocav = trans.matrix; trans.matrix.addnocav[1,4] = trans.matrix.addnocav[1,4] + 1
        numerator14 = q.split(trans.matrix.addnocav, beta.dead, alpha.death.matrix, rho)
        current.B[1,4] = numerator14/denominator; current.B[1,1] = 1 - current.B[1,4]
        
        trans.matrix.addmild = trans.matrix; trans.matrix.addmild[2,4] = trans.matrix.addnocav[2,4] + 1
        numerator24 = q.split(trans.matrix.addmild, beta.dead, alpha.death.matrix, rho)
        current.B[2,4] = numerator24/denominator; current.B[2,2] = 1 - current.B[2,4]
        
        trans.matrix.addsevere = trans.matrix; trans.matrix.addsevere[3,4] = trans.matrix.addsevere[3,4] + 1
        numerator34 = q.split(trans.matrix.addsevere, beta.dead, alpha.death.matrix, rho)
        current.B[3,4] = numerator34/denominator; current.B[3,3] = 1 - current.B[3,4]
        
        list.of.current.Bs[[t-1]] = current.B
      } else {
        current.B = matrix(0, nrow = 4, ncol = 4)
        
        # Need to deal with death times
        temp.which.Omega = min(max(which(current.time >= bar.time)), length(bar.Omega))
        
        current.Omega = bar.Omega[temp.which.Omega]
        n.vector = c(num.nocav, num.mild, num.severe)
        
        rate.alive.baseline = nu[1] * zeta(n.vector, beta.alive, rho = rho) * rho
        rate.death.baseline = nu[2] * zeta(n.vector, beta.dead, rho = rho) * rho
        
        n.vector.addnocav = n.vector; n.vector.addnocav[1] = n.vector.addnocav[1] + 1
        rate.alive.addnocav = nu[1] * zeta(n.vector.addnocav, beta.alive, rho = rho) * rho
        rate.death.addnocav = nu[2] * zeta(n.vector.addnocav, beta.dead, rho = rho) * rho
        total.rate.addnocav = rate.alive.addnocav + rate.death.addnocav - (rate.alive.baseline + rate.death.baseline)
        # total.rate.addhealthy = rate1.addhealthy - rate1.baseline
        
        n.vector.addmild = n.vector; n.vector.addmild[2] = n.vector.addmild[2] + 1
        rate.alive.addmild = nu[1] * zeta(n.vector.addmild, beta.alive, rho = rho) * rho
        rate.death.addmild = nu[2] * zeta(n.vector.addmild, beta.dead, rho = rho) * rho
        total.rate.addmild = rate.alive.addmild + rate.death.addmild - (rate.alive.baseline + rate.death.baseline)
        # total.rate.addhealthy = rate1.addhealthy - rate1.baseline
        trans.matrix.mild = matrix(0,nrow = 3, ncol = 3)
        diag(trans.matrix.mild) = n.vector; trans.matrix.mild[2,1] = 1
        temp.mildtonocav = q.split(trans.matrix.mild, beta.alive, alpha.alive.matrix, rho)
        trans.matrix.mild[2,1] = 0; trans.matrix.mild[2,3] = 1
        temp.mildtosevere = q.split(trans.matrix.mild, beta.alive, alpha.alive.matrix, rho)
        frac.mildtonocav = temp.mildtonocav/ (temp.mildtonocav + temp.mildtosevere)
        
        n.vector.addsevere = n.vector; n.vector.addsevere[3] = n.vector.addsevere[3] + 1
        rate.alive.addsevere = nu[1] * zeta(n.vector.addsevere, beta.alive, rho = rho) * rho
        rate.death.addsevere = nu[2] * zeta(n.vector.addsevere, beta.dead, rho = rho) * rho
        total.rate.addsevere = rate.alive.addsevere + rate.death.addsevere - (rate.alive.baseline + rate.death.baseline)
        # total.rate.addhealthy = rate1.addhealthy - rate1.baseline
        
        diag(current.B) = 1 - c(total.rate.addnocav, total.rate.addmild, total.rate.addsevere,0)/current.Omega
        current.B[1,2] = (rate.alive.addnocav - rate.alive.baseline)/(current.Omega)
        current.B[1,4] = (rate.death.addnocav - rate.death.baseline)/(current.Omega)
        current.B[2,1] =  frac.mildtonocav * (rate.alive.addmild - rate.alive.baseline)/(current.Omega)
        current.B[2,3] = (1-frac.mildtonocav) * (rate.alive.addmild - rate.alive.baseline)/(current.Omega)
        current.B[2,4] = (rate.death.addsevere - rate.death.baseline)/(current.Omega)
        current.B[3,2] = (rate.alive.addsevere - rate.alive.baseline)/(current.Omega)
        current.B[3,4] = (rate.death.addsevere - rate.death.baseline)/(current.Omega)
        list.of.current.Bs[[t-1]] = current.B
      }
    } else {
      current.B = diag(rep(1,4))
      list.of.current.Bs[[t-1]] = current.B
    }
    
    if(any(current.temp$person.obs.times == TRUE)) {
      L.matrix[,t] = as.numeric(1:4 == current.temp$person.obs.states[current.temp$person.obs.times == TRUE])
    } else {
      L.matrix[,t] = c(rep(1, 3),0)
    } 
    
    Alpha.matrix[,t] =  (L.matrix[,t-1] * Alpha.matrix[,t-1]) %*% current.B
    
    if ( any(!is.na(current.temp$num.nocav)) ) {
      num.nocav = current.temp$num.nocav[!is.na(current.temp$num.nocav)]
      num.mild = current.temp$num.mild[!is.na(current.temp$num.mild)]
      num.severe = current.temp$num.severe[!is.na(current.temp$num.severe)]
      num.dead = current.temp$num.dead[!is.na(current.temp$num.dead)]
    }
    if(any(current.B <0)) {print(t)}
  }
  
  ## And now we backwards sample
  Beta.matrix = Alpha.matrix * 0; sampled.trajectory = vector(length = ncol(Alpha.matrix))
  T = ncol(Alpha.matrix)
  Beta.matrix[,T] = L.matrix[,T]
  # Beta.matrix[,T] = Beta.matrix[,T]/ sum(Beta.matrix[,T])
  sampled.trajectory[T] = sample(1:4, prob = Beta.matrix[,T], size = 1)
  for (t in (T-1):1) {
      B.t = list.of.current.Bs[[t]]
      Beta.matrix[,t] = Alpha.matrix[,t] * L.matrix[,t] * B.t[,sampled.trajectory[t+1]]
      Beta.matrix[,t] = Beta.matrix[,t]/ sum(Beta.matrix[,t])
      if(t == T-1 & sampled.trajectory[T] == 4) {
        sampled.trajectory[t] = sample(1:3, prob = Beta.matrix[1:3,t], size = 1)
      } else {
        sampled.trajectory[t] = sample(1:4, prob = Beta.matrix[,t], size = 1)
      }
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

sample.all.new.users <- function(all.person.obs.data, old.all.person.trans.data, censoring.data,
                                 current.params, current.multiple, rho, 
                                 initial.iter = FALSE, max.attempts = 20) {
  all.users = unique(all.person.obs.data$user)
  if (initial.iter) { 
    keep.users = c(); new.all.person.trans.data = c(); temp.censoring.data = c();
    for(patient.id in all.users) {
      # print(patient.id)
      keep.users = c(patient.id, keep.users)
      new.all.person.trans.data = rbind(new.all.person.trans.data, old.all.person.trans.data[is.element(old.all.person.trans.data$user,patient.id),])
      temp.censoring.data = rbind(temp.censoring.data, censoring.data[is.element(censoring.data$user,patient.id),])
      new.all.person.obs.data = all.person.obs.data[is.element(all.person.obs.data$user,keep.users),]
      complete = FALSE; how.many.to.completion = 0; attempt = 1; temp.multiple = current.multiple
      while(complete == FALSE) {
        how.many.to.completion = how.many.to.completion + 1
        temp = try(sample.one.new.user(patient.id, keep.users, new.all.person.obs.data, new.all.person.trans.data, 
                                       temp.censoring.data, current.params, temp.multiple, rho, initial.iter), TRUE)
        if(!is(temp,"try-error")) {
          new.all.person.trans.data = temp
          complete = TRUE
        } else {
          if (attempt == max.attempts) {
            complete = TRUE  
            print("Hit max attempts, so moving on.")
            break
          } else {
            attempt = attempt + 1
            temp.multiple = 2*temp.multiple
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
      # print(patient.id)
      complete = FALSE; how.many.to.completion = 0
      max.attempts = 10; attempt = 1; temp.multiple = current.multiple
      while(complete == FALSE) {
        how.many.to.completion = how.many.to.completion + 1
        temp = try(sample.one.new.user(patient.id, all.users, new.all.person.obs.data, new.all.person.trans.data, 
                                       censoring.data, current.params, temp.multiple, rho), TRUE)
        if(!is(temp,"try-error")) {
          new.all.person.trans.data = temp
          complete = TRUE
        } else {
          if (attempt == max.attempts) {
            complete = TRUE  
            print("Hit max attempts, so moving on.")
          } else {
            attempt = attempt + 1
            temp.multiple = 2*temp.multiple
            print("In a fail-safe loop, hold on one moment please. Thanks for your patience.")
          }
        }
      }
    } 
  }
  return(new.all.person.trans.data)
}

sample.one.new.user <- function(patient.id, all.users, all.person.obs.data, 
                                all.person.trans.data, censoring.data, 
                                current.params, current.multiple, rho, initial.iter = FALSE) {
  # print(patient.id)
  temp.extract = extract.patient(all.person.obs.data, all.person.trans.data, patient.id)
  chosen.person.obs.data = temp.extract$chosen.person.obs.data
  chosen.person.trans.data = temp.extract$chosen.person.trans.data
  other.person.trans.data = temp.extract$other.person.trans.data
  max.time = max(chosen.person.obs.data$time)
  
  temp.censoring.data = censoring.data[is.element(censoring.data$user,other.person.trans.data$user),]
  
  temp.extractOmega.and.A = extract.Omega.and.A(chosen.person.obs.data, chosen.person.trans.data, other.person.trans.data, 
                                                temp.censoring.data, current.params, current.multiple = current.multiple , rho) 
  
  temp.new.potential.trans.times = generate.points(temp.extractOmega.and.A$bar.Omega, temp.extractOmega.and.A$bar.Omega.time, 
                                                   temp.extractOmega.and.A$bar.A, temp.extractOmega.and.A$bar.Omega.time,
                                                   max.time, initial.iter)
  
  temp.constructinfo = construct.info(temp.new.potential.trans.times, chosen.person.trans.data, other.person.trans.data, 
                                      chosen.person.obs.data, temp.censoring.data, current.params)
  
  ## REMOVE BAD OBS - Due to censoring
  good.obs = (temp.constructinfo$new.trans.time == TRUE) | (temp.constructinfo$person.obs.times == TRUE) | 
    (any(is.na(temp.constructinfo[,2:4]) == TRUE))
  temp.constructinfo = temp.constructinfo[good.obs,]
  
  # NEED TO DEAL WITH CAN'T TRANS TO 4 PRIOR TO FAILURE TIME
  temp.newsample <- filterfwd.backsampl(temp.constructinfo, current.params, rho, temp.extractOmega.and.A$bar.Omega, temp.extractOmega.and.A$bar.Omega.time)
  
  all.person.trans.data <- build.new.all.person.trans.data(temp.newsample, patient.id, other.person.trans.data) 
  
  return(all.person.trans.data)
}

nu.gibbs.updates <- function(data, current.params, rho, prior.params) {
  total.zeta.integral = c(0,0); switch.count = c(0,0)
  nu = current.params[1:2]; beta.alive = c(1,current.params[3:4]); beta.dead = c(1,current.params[5:6]); 
  for(t in 2:nrow(data)) {
    isdeathtime = data$failure.time[t]
    num.nocav = data$num.nocav[t-1]; num.mild = data$num.mild[t-1]
    num.severe = data$num.severe[t-1]
    n.vector = c(num.nocav, num.mild, num.severe)
    window.length = data$time[t] - data$time[t-1]
    rate.alive = zeta(n.vector, beta.alive, rho) * rho
    rate.dead  = zeta(n.vector, beta.dead, rho) * rho
    current.rates = c(rate.alive, rate.dead)
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

beta.gibbs.updates <- function(data, new.nu.params, new.alpha.alive, current.params, rho, prior.params) {
  temp.params = c(new.nu.params, current.params[3:(length(current.params)-1)], new.alpha.alive)
  proposal.beta.params = log(temp.params[3:6]) + rnorm(4, sd = prior.params$step.size)
  # proposal.alpha.params = rbeta(1, shape1 = temp.params[7] * prior.params$alpha.step.size, 
  #                               shape2 = (1 - temp.params[7]) * prior.params$alpha.step.size)
  proposal.params = c(temp.params[1:2], exp(proposal.beta.params), temp.params[7])
  temp.proposal.params = new.params = temp.params
  all.acceptance.rate = vector(length = 4)
  for (which.beta in 1:4) {
    # Sequentially perform accept/reject to avoid stickiness
    current.llik = total.llik(data, rho)(new.params)
    current.prior.llik = -sum(log(dnorm(log(temp.params[2+which.beta]), mean = prior.params$log.mu, sd = prior.params$log.sd)))
    
    proposal.prior.llik = -sum(log(dnorm(proposal.beta.params[which.beta], mean = prior.params$log.mu, sd = prior.params$log.sd)))
    temp.proposal.params[2 + which.beta]  = proposal.params[2 + which.beta] 
    proposal.llik = total.llik(data, rho)(temp.proposal.params)
    
    acceptance.rate = min(1, exp(current.llik + current.prior.llik - proposal.llik - proposal.prior.llik)  )
    accept.proposal = as.logical(rbinom(n = 1, size = 1, prob = acceptance.rate) ==1)
    new.params = accept.proposal * temp.proposal.params + (1-accept.proposal) * temp.params
    all.acceptance.rate[which.beta] = acceptance.rate
  }
  
  return(
    list( "new.params" = new.params, "accept.rate" = all.acceptance.rate)
  )
}


alpha.alive.gibbs.updates <- function(data, prior.params) {
  new.alpha.alive = sum(data$switch2a[data$failure.time!=1]) + prior.params$alpha.alive
  new.beta.alive = sum(data$switch2b[data$failure.time!=1]) + prior.params$beta.alive
  return( 
    rbeta(1, shape1  = new.alpha.alive, shape2 = new.beta.alive)
  )
}

Survival.fit <- function(current.global, current.params, rho, init.time, grid.length = 0.01) {
  
  unique.trans.times = unique(current.global$time)
  Survival.init.nocav = c(1, 0, 0, 0); computed.at.times = c(init.time) 
  Survival.initnocav.matrix = as.matrix(Survival.init.nocav, ncol = 1)
  
  Survival.init.mild = c(0, 1, 0, 0); computed.at.times = c(init.time) 
  Survival.initmild.matrix = as.matrix(Survival.init.mild, ncol = 1)
  
  Survival.init.severe = c(0, 0, 1, 0); computed.at.times = c(init.time) 
  Survival.initsevere.matrix = as.matrix(Survival.init.severe, ncol = 1)
  
  nu = current.params[1:2]; beta.alive = c(1,current.params[3:4]); beta.dead = c(1,current.params[5:6]); 
  alpha.alive = current.params[7]
  # Setup up alive matrix
  alpha.alive.matrix = matrix(0, nrow = 3, ncol = 3)
  alpha.alive.matrix[1,2] = 1; alpha.alive.matrix[2,1] = alpha.alive
  alpha.alive.matrix[2,3] = 1-alpha.alive; alpha.alive.matrix[3,2] = 1
  diag(alpha.alive.matrix) = 1
  # Set up death time matrix
  alpha.death.matrix = matrix(0,nrow = 3, ncol = 4) 
  alpha.death.matrix[1:3,4] = 1; diag(alpha.death.matrix) = 1
  
  for (t in 2:length(unique.trans.times)) {
    
    num.nocav = current.global$num.nocav[t-1]; num.mild = current.global$num.mild[t-1]; num.severe = current.global$num.severe[t-1]
    n.vector = c(num.nocav, num.mild, num.severe)
    
    rate.alive.baseline = nu[1] * zeta(n.vector, beta.alive, rho) * rho
    rate.death.baseline  = nu[2] * zeta(n.vector, beta.dead, rho) * rho
    
    n.vector.addnocav = n.vector; n.vector.addnocav[1] = n.vector.addnocav[1] + 1
    rate.alive.addnocav = nu[1] * zeta(n.vector.addnocav, beta.alive, rho = rho) * rho
    rate.death.addnocav = nu[2] * zeta(n.vector.addnocav, beta.dead, rho = rho) * rho
    total.rate.addnocav = rate.alive.addnocav + rate.death.addnocav - (rate.alive.baseline + rate.death.baseline)
    # total.rate.addhealthy = rate1.addhealthy - rate1.baseline
    
    n.vector.addmild = n.vector; n.vector.addmild[2] = n.vector.addmild[2] + 1
    rate.alive.addmild = nu[1] * zeta(n.vector.addmild, beta.alive, rho = rho) * rho
    rate.death.addmild = nu[2] * zeta(n.vector.addmild, beta.dead, rho = rho) * rho
    total.rate.addmild = rate.alive.addmild + rate.death.addmild - (rate.alive.baseline + rate.death.baseline)
    # total.rate.addhealthy = rate1.addhealthy - rate1.baseline
    trans.matrix.mild = matrix(0,nrow = 3, ncol = 3)
    diag(trans.matrix.mild) = n.vector; trans.matrix.mild[2,1] = 1
    temp.mildtonocav = q.split(trans.matrix.mild, beta.alive, alpha.alive.matrix, rho)
    trans.matrix.mild[2,1] = 0; trans.matrix.mild[2,3] = 1
    temp.mildtosevere = q.split(trans.matrix.mild, beta.alive, alpha.alive.matrix, rho)
    frac.mildtonocav = temp.mildtonocav/ (temp.mildtonocav + temp.mildtosevere)
    
    n.vector.addsevere = n.vector; n.vector.addsevere[3] = n.vector.addsevere[3] + 1
    rate.alive.addsevere = nu[1] * zeta(n.vector.addsevere, beta.alive, rho = rho) * rho
    rate.death.addsevere = nu[2] * zeta(n.vector.addsevere, beta.dead, rho = rho) * rho
    total.rate.addsevere = rate.alive.addsevere + rate.death.addsevere - (rate.alive.baseline + rate.death.baseline)
    # total.rate.addhealthy = rate1.addhealthy - rate1.baseline
    
    current.Q = matrix(nrow = 4, ncol = 4)
    current.Q[1,1] = - total.rate.addnocav
    current.Q[1,2] = (rate.alive.addnocav - rate.alive.baseline) 
    current.Q[1,3] = 0
    current.Q[1,4] = (rate.death.addnocav - rate.death.baseline)
    
    current.Q[2,2] = - total.rate.addmild
    current.Q[2,1] = (rate.alive.addmild - rate.alive.baseline) * frac.mildtonocav
    current.Q[2,3] = (rate.alive.addmild - rate.alive.baseline) * (1-frac.mildtonocav)
    current.Q[2,4] = (rate.death.addmild - rate.death.baseline)
    
    current.Q[3,3] = - total.rate.addsevere
    current.Q[3,2] = (rate.alive.addsevere - rate.alive.baseline)
    current.Q[3,1] = 0
    current.Q[3,4] = (rate.death.addsevere - rate.death.baseline)
    
    current.Q[4,1:4] = 0
    
    ## Expand on a 0.01 grid 
    ## To compute survival at these time
    temp.times = c(seq(current.global$time[t-1], current.global$time[t], grid.length), current.global$time[t])
    diff.gaps = diff(temp.times)
    computed.at.times = c(computed.at.times, temp.times[2:length(temp.times)])
    for (i in 1:length(diff.gaps)) {
      current.B = as.matrix(expm(current.Q*diff.gaps[i]))
      Survival.init.nocav = Survival.init.nocav%*%current.B
      Survival.initnocav.matrix = cbind(Survival.initnocav.matrix, t(Survival.init.nocav))
      Survival.init.mild = Survival.init.mild%*%current.B
      Survival.initmild.matrix = cbind(Survival.initmild.matrix, t(Survival.init.mild))
      Survival.init.severe = Survival.init.severe%*%current.B
      Survival.initsevere.matrix = cbind(Survival.initsevere.matrix, t(Survival.init.severe))
    }
    
    
    if (
      ( (current.global$switch1[t] != 0) | current.global$switch2a[t] != 0 |
        (current.global$switch2b[t] != 0) | (current.global$switch3[t] != 0) ) &  
      (current.global$failure.time[t] == 0)
    ) { 
      
      temp.switch1 = current.global$switch1[t]; temp.switch2a = current.global$switch2a[t]
      temp.switch2b = current.global$switch2b[t]; temp.switch3 = current.global$switch3[t]
      ## Construct trans.matrix
      r.nocav = num.nocav - temp.switch1; r.mild = num.mild - temp.switch2a - temp.switch2b;
      r.severe = num.severe - temp.switch3;
      if(r.mild < 0) {r.mild = 0}
      trans.matrix = matrix(0, nrow = 3, ncol = 3)
      diag(trans.matrix) = c(r.nocav, r.mild, r.severe);
      trans.matrix[1,2] = temp.switch1; trans.matrix[2,1] = temp.switch2a; 
      trans.matrix[2,3] = temp.switch2b; trans.matrix[3,2] = temp.switch3
      
      denominator = q.split(trans.matrix, beta.alive, alpha.alive.matrix, rho = rho)
      
      current.B = matrix(0, nrow = 4, ncol = 4); current.B[4,4] = 1
      
      trans.matrix.addnocav = trans.matrix; trans.matrix.addnocav[1,2] = trans.matrix.addnocav[1,2] + 1
      numerator12 = q.split(trans.matrix.addnocav, beta.alive, alpha.alive.matrix, rho)
      
      current.B[1,2] = numerator12/denominator; current.B[1,1] = 1 - current.B[1,2]
      
      trans.matrix.addmild.2 = trans.matrix.addmild.2a = trans.matrix.addmild.2b = trans.matrix; 
      trans.matrix.addmild.2[2,2] = trans.matrix.addmild.2[2,2] + 1
      trans.matrix.addmild.2a[2,1] = trans.matrix.addmild.2a[2,1] + 1
      trans.matrix.addmild.2b[2,3] = trans.matrix.addmild.2b[2,3] + 1
      numerator22 = q.split(trans.matrix.addmild.2, beta.alive, alpha.alive.matrix, rho)
      temp.input21 = q.split(trans.matrix.addmild.2a, beta.alive, alpha.alive.matrix, rho)
      temp.input23 = q.split(trans.matrix.addmild.2b, beta.alive, alpha.alive.matrix, rho)
      
      current.B[2,1] = temp.input21/denominator
      current.B[2,3] = temp.input23/denominator
      current.B[2,2] = 1-current.B[2,1] - current.B[2,3]
      
      trans.matrix.addsevere = trans.matrix; trans.matrix.addsevere[3,2] = trans.matrix.addsevere[3,2] + 1
      numerator32 = q.split(trans.matrix.addsevere, beta.alive, alpha.alive.matrix, rho)
      
      current.B[3,2] = numerator32/denominator; current.B[3,3] = 1 - current.B[3,2]
      
      Survival.init.nocav = Survival.init.nocav%*%current.B
      Survival.initnocav.matrix = cbind(Survival.initnocav.matrix, t(Survival.init.nocav))
      
      Survival.init.mild = Survival.init.mild%*%current.B
      Survival.initmild.matrix = cbind(Survival.initmild.matrix, t(Survival.init.mild))
      
      Survival.init.severe = Survival.init.severe%*%current.B
      Survival.initsevere.matrix = cbind(Survival.initsevere.matrix, t(Survival.init.severe))
      
      computed.at.times = c(computed.at.times, current.global$time[t])
      
    } else if(current.global$failure.time[t] == 1) {
      
      temp.switch1 = current.global$switch1[t]; temp.switch2a = current.global$switch2a[t]
      temp.switch2b = current.global$switch2b[t]; temp.switch3 = current.global$switch3[t]
      ## Construct trans.matrix
      r.nocav = num.nocav - temp.switch1; r.mild = num.mild - temp.switch2a - temp.switch2b;
      r.severe = num.severe - temp.switch3;
      if(r.mild < 0) {r.mild = 0}
      trans.matrix = matrix(0, nrow = 3, ncol = 4)
      diag(trans.matrix) = c(r.nocav, r.mild, r.severe);
      trans.matrix[1,4] = temp.switch1; trans.matrix[2,4] = temp.switch2a + temp.switch2b; 
      trans.matrix[3,4] = temp.switch3
      
      denominator = q.split(trans.matrix, beta.dead, alpha.death.matrix, rho = rho)
      
      current.B = matrix(0, nrow = 4, ncol = 4); current.B[4,4] = 1
      
      trans.matrix.addnocav = trans.matrix; trans.matrix.addnocav[1,4] = trans.matrix.addnocav[1,4] + 1
      numerator14 = q.split(trans.matrix.addnocav, beta.dead, alpha.death.matrix, rho)
      
      current.B[1,4] = numerator14/denominator; current.B[1,1] = 1 - current.B[1,4]
      
      trans.matrix.addmild.2 = trans.matrix; 
      trans.matrix.addmild.2[2,4] = trans.matrix.addmild.2[2,4] + 1
      numerator24 = q.split(trans.matrix.addmild.2, beta.dead, alpha.death.matrix, rho)

      current.B[2,4] = numerator24/denominator
      current.B[2,2] = 1-current.B[2,4]
      
      trans.matrix.addsevere = trans.matrix; trans.matrix.addsevere[3,4] = trans.matrix.addsevere[3,4] + 1
      numerator34 = q.split(trans.matrix.addsevere, beta.dead, alpha.death.matrix, rho)
      
      current.B[3,4] = numerator34/denominator; current.B[3,3] = 1 - current.B[3,4]
      
      Survival.init.nocav = Survival.init.nocav%*%current.B
      Survival.initnocav.matrix = cbind(Survival.initnocav.matrix, t(Survival.init.nocav))
      
      Survival.init.mild = Survival.init.mild%*%current.B
      Survival.initmild.matrix = cbind(Survival.initmild.matrix, t(Survival.init.mild))
      
      Survival.init.severe = Survival.init.severe%*%current.B
      Survival.initsevere.matrix = cbind(Survival.initsevere.matrix, t(Survival.init.severe))
      
      computed.at.times = c(computed.at.times, current.global$time[t])
      
    }
    
  }
  
  return(cbind(
    computed.at.times,
    Survival.initnocav.matrix[4,],
    Survival.initmild.matrix[4,],
    Survival.initsevere.matrix[4,]
  ))
  
}



logLik.iid <- function(person.data) {
  internal.fn <- function(parameters) {
    ## Log Likelihood for iid data
    ## Sum over people 
    Q = matrix(0, nrow = 4, ncol = 4)
    Q[1,c(2,4)] = parameters[1:2]; Q[1,1] = - sum(Q[1,2:4])
    Q[2,c(1,3:4)] = parameters[3:5]; Q[2,2] = - sum(Q[2,c(1,3:4)])
    Q[3,c(2,4)] = parameters[6:7]; Q[3,3] = - sum(Q[3, c(1:2,4)])
    Q[4,1:4] = 0
    loglik.total = 0
    for (user in unique(person.data$user)) {
      user.data = person.data[person.data$user == user,]
      for (row in 1:(nrow(user.data)-1)) {
        gap = user.data$time[row+1] - user.data$time[row]
        current.state = user.data$state[row]
        next.state = user.data$state[row+1]
        trans.P = expm(Q*gap)
        if(next.state != 4) {
          loglik.total = loglik.total + log(trans.P[current.state, next.state])
        } else {
          
          loglik.total = loglik.total + log(sum(-trans.P[current.state, 1:3] * Q[1:3,4]/Q[1,1]))
        }
      }
    }
    return(-loglik.total)
  }
  return(internal.fn)
}


nonparametric_estimator <- function(init.distribution, global.data) {
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
  return(result.distribution)
}
