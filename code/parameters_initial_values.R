
migration_situation = function(set, b0, b1){
  if(set == "set1"){
    return(list(beta0 = c(rep(-9,N.v),rep(-8.5,N.r)) ,
                beta1 = c(rep(-.0005,N.v),rep(-.0003,N.r))))
  }
  if(set == "set2"){
    return(list(beta0 = c(rep(-5,N.v),rep(-8,N.r)) ,
                beta1 = c(rep(-.1,N.v),rep(-.15,N.r))))
  }
  if(set == "set3"){
    return(list(beta0 = c(rep(-9.8,N.v),rep(-9.5,N.r)) ,
                beta1 = c(rep(-.003,N.v),rep(-.002,N.r))))
  }
  if(set == "set4"){
    return(list(beta0 = c(rep(-5,N.v),rep(-5,N.r)) ,
                beta1 = c(rep(-.05,N.v),rep(-.05,N.r))))
  }
  else{
    return(list(beta0 = c(rep(b0[1],N.v),rep(b0[2],N.r)) ,
                beta1 = c(rep(b1[1],N.v),rep(b1[2],N.r))))
  }
}


testrate_situation = function(set,testrate){
  if(set == "set1"){
    return(1.01)
  }
  if(set == "set2"){
    return(1.03)
  }
  if(set == "set3"){
    return(1.05)
  }
  else{
    return(testrate)
  }
}


vaccine_situation = function(set,vaccine){
  if(set == "set1"){
    return(c(rep(0,N.v),rep(0,N.r)))
  }
  if(set == "set2"){
    return(c(rep(0.003,N.v),rep(0.003,N.r)))
  }
  if(set == "set3"){
    return(c(rep(0.003,N.v),rep(0,N.r)))
  }
  if(set == "set4"){
    #return(c(rep(0.005,N.v),rep(0.003,N.r)))
    return(c(rep(0.01,N.v),rep(0.005,N.r))) ## smaller vac rate for the robust 0.002/0.003
  }
  if(set == "set5"){
    return(c(rep(0.02,N.v),rep(0.01,N.r))) ## smaller vac rate for the robust 0.002/0.003
  }
  else{
    return(vaccine)
  }
}


N.v = 5 # number of nodes with more vulnerable population
N.r = 5 # number of nodes with less vulnerable (robust) population
N = N.v + N.r # total number of nodes

T = 400 # number of days (time steps)
T.Q = 14 # length of quarantine period
Max = 5*10^5 # maximum population size

psi.true = 0.99 # true positive detection probability of the test
psi.false = 0.02 # false positive detection probability of the test

tau = matrix(0,N,T) # true positive probability for an asymptomatic person to be detected and quarantined
phi = matrix(0,N,T) # false positive probability for an uninfected person to be (falsely) detected and quarantined
iota = matrix(0,N,T)  # basic infection probability

init.states = list()
init.states$omega = matrix(0,N,N)  # omega[j,k]= rate of migration from node j to node k

# States of the system

init.states$S = matrix(0,N,T)  # susceptible (uninfected) population
init.states$Gp= matrix(0,N,T)  # daily added uninfected and quarantined population
init.states$G = matrix(0,N,T)  # uninfected and quarantined population
init.states$J = matrix(0,N,T)  # infected but asymptomatic population
init.states$I = matrix(0,N,T)  # infected and symptomatic population
init.states$Q = matrix(0,N,T)  # asymptomatic and quarantined population
init.states$H = matrix(0,N,T)  # hospitalized population
init.states$R = matrix(0,N,T)  # recovered (cured) population
init.states$D = matrix(0,N,T)  # number of deaths
init.states$V = matrix(0,N,T)  # number of vaccinations
init.states$IV = matrix(0,N,T) # number of immunized due to vaccination
init.states$curr_pop = matrix(0, nrow = N, ncol = T)
init.states$curr_max = matrix(0, nrow = N, ncol = T)

init.states$mu = matrix(0,N,T)  # number of people an average individual meets on a given day
init.states$Test = matrix(0,N,T)  # number of tests done on any day
init.states$cost = matrix(0,N,T) #Cost associated with the lockdown
##################################
## System dynamics
###############

## Initialization of states
init.states$curr_max =  matrix(0, nrow = N, ncol = T)
init.states$J[,1] = c(0,100, rep(0,N-2))
init.states$S[1:3,1] = Max*0.9 - init.states$J[1:3,1]
init.states$S[4:8,1] = Max*0.7 - init.states$J[4:8,1]
init.states$S[9:10,1] = Max*0.5 - init.states$J[9:10,1]

init.states$mu[,1] = c(rep(3,N.v), rep(5,N.r)) # c(rep(3,N.v), rep(6,N.r))//c(rep(5,N.v), rep(10,N.r)) #
init.states$Test[,1] = c(rep(1000,N.v), rep(500,N.r))
# beta0 = c(rep(-10,N.v),rep(-9,N.r)) 
# beta1 = c(rep(-.05,N.v),rep(-.01,N.r))

init.states$omega = matrix(0,N,N)  

##epidemiological parameters

epi_para = list()
epi_para$iota.c = .1
epi_para$gamma = c(rep(0.03,N.v), rep(0.01, N.r)) # probability of transition from asymptomatic to symptomatic state
epi_para$alpha = c(rep(0.005,N.v), rep(0.02, N.r)) # probability of transition from asymptomatic to recovered state
epi_para$zeta = c(rep(0.05,N.v), rep(0.1,N.r))  # probability of transition from symptomatic to recovered state
epi_para$rho  = c(rep(0.07,N.v), rep(0.15,N.r))   # probability of transition from hospitalized to recovered state
epi_para$xi = c(rep(0.2,N.v), rep(0.1,N.r))   # probability of transition from symptomatic to hospitalized state
epi_para$delta = c(rep(0.03,N.v), rep(0.01,N.r))  # probability that a hospitalized person dies on any given day

