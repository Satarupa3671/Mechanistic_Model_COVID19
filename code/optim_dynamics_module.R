update_states = function(t, states, epi_para, migration, test_rate, vac_rate, t_vac,intervention_para, vac_efficacy=0.8, testbound){
  iota.c = epi_para$iota.c
  gamma = epi_para$gamma
  alpha = epi_para$alpha
  zeta = epi_para$zeta
  xi = epi_para$xi
  rho = epi_para$rho
  delta = epi_para$delta
  kappa = intervention_para$kappa
  threshold = intervention_para$threshold
  
  beta0 = migration$beta0
  beta1 = migration$beta1
  
  for(k in 1:N){
    states$omega[,k] = pmax(exp(beta0[k] + beta1[k] * delta * states$H[,t])/(1 + exp(beta0[k] + beta1[k] * delta * states$H[,t])),0.0000001)
  } 
  #states$omega[,1] = pmax(exp(beta0[1] + beta1[1] * delta * states$H[,t])/(1 + exp(beta0[1] + beta1[1] * delta * states$H[,t])),0.0000001)
  diag(states$omega) = 0
  
  mu_lockdown <- {
    ###### 4/19/21 replace || with |
    cond <- ((states$H[,t]/states$S[,1]) > threshold) | ((states$I[,t]/states$S[,1]) > .004)
    
    ######
    val1 <- c(rep(kappa^2, times = c(N.v,N.r)) * states$mu[,1], 
              rep((1 - kappa), times = c(N.v,N.r)) *(states$S[,t]+ states$J[,t] + 
                                                       states$R[,t] + states$IV[,t]))
    
    val2 <- c(states$mu[,1], rep(0,N))
    matrix(mapply(ifelse, cond, val1, val2), nrow = 2, byrow = TRUE)
  }
  
  states$mu[,t] <- mu_lockdown[1,]
  states$cost[,(t+1)] <- mu_lockdown[2,] #+ states$cost[,t] 
  
  # infection rate equation
  iota[,t] = iota.c * states$mu[,t] * states$J[,t]/(states$S[,t] + states$J[,t] + states$R[,t] + states$IV[,t])
  # updating true positive rate
  tau[,t] = psi.true * states$Test[,t]/(states$S[,t] + states$J[,t] + states$IV[,t] +states$R[,t])
  
  # updating false positive rate
  phi[,t] = psi.false * states$Test[,t]/(states$S[,t] + states$J[,t] + states$IV[,t] +states$R[,t])
  
  ## updating active states
  for(k in 1:N){ 
    
    states$S[k,t+1] = states$S[k,t] + sum(states$omega[,k] * states$S[,t]) -
      sum(states$omega[k,]* states$S[k,t]) -
      (iota[k,t] + phi[k,t] ) * states$S[k,t] + states$Gp[k,max(t-T.Q+1,1)] 
    
    states$J[k,t+1] = states$J[k,t] + sum(states$omega[,k] * states$J[,t]) - 
      sum(states$omega[k,]* states$J[k,t]) - 
      (tau[k,t] + gamma[k] + alpha[k]) * states$J[k,t] + iota[k,t] * states$S[k,t] 
    
  } 
  
  states$Gp[,t+1]= phi[,t] * states$S[,t]
  states$G[,t+1] = states$G[,t] + states$Gp[,t+1] -  states$Gp[,max(t-T.Q+1,1)]   
  
  states$I[,t+1] = states$I[,t] + gamma * (states$J[,t] + states$Q[,t]) - 
    (zeta + xi) * states$I[,t]
  
  states$Q[,t+1] = states$Q[,t] - gamma * states$Q[,t] + tau[,t] * states$J[,t] - 
    alpha * states$Q[,t]
  
  states$H[,t+1] = states$H[,t] + xi * states$I[,t] - (rho + delta) * states$H[,t]
  states$D[,t+1] = states$D[,t] + delta * states$H[,t]
  states$R[,t+1] = states$R[,t] + alpha * (states$J[,t] + states$Q[,t]) + 
    zeta * states$I[,t] + rho * states$H[,t] 
  
  if(t > t_vac){
    
    states$V[,t+1] = states$V[,t] + vac_rate * (states$S[,t] + states$J[,t] + states$R[,t])
    states$IV[,t+1] = states$IV[,t] + vac_efficacy * vac_rate * states$S[,t] + vac_rate * states$R[,t] #(states$V[,t+1] -states$V[,t])
    states$R[,t+1] = states$R[,t] + alpha * (states$J[,t] + states$Q[,t]) + 
      zeta * states$I[,t] + rho * states$H[,t] - vac_rate * states$R[,t]
    for(k in 1:N){ 
      
      states$S[k,t+1] = states$S[k,t] + sum(states$omega[,k] * states$S[,t]) -
        sum(states$omega[k,]* states$S[k,t]) -(iota[k,t] + phi[k,t] + vac_efficacy * vac_rate[k]) * states$S[k,t] + 
        states$Gp[k,max(t-T.Q+1,1)] 
    }
    
  }
  
  ## ensuring feasible values of the active states
  for(k in 1:N){
    states$curr_max[k,t] = states$S[k,t] + 
      states$J[k,t] + states$I[k,t] + states$Q[k,t] + states$G[k,t] +
      states$H[k,t] + states$R[k,t] + states$D[k,t] + states$IV[k,t]
  }
  
  for(k in 1:N){
    states$S[k,t+1] = pmin(pmax(states$S[k,t+1],0), states$curr_max[k,t]); #states$S[states$S[k,t+1]<0.05,t+1]=0  
    states$R[k,t+1] = pmin(pmax(states$R[k,t+1],0), states$curr_max[k,t]); #R[R[,t+1]<0.05,t+1]=0
    states$J[k,t+1] = pmin(pmax(states$J[k,t+1],0), states$curr_max[k,t]); #states$J[states$J[k,t+1]<0.05,t+1]=0
    states$I[k,t+1] = pmin(pmax(states$I[k,t+1],0), states$curr_max[k,t]); #states$I[states$I[k,t+1]<0.05,t+1]=0
    states$Q[k,t+1] = pmin(pmax(states$Q[k,t+1],0), states$curr_max[k,t]); #states$Q[states$Q[k,t+1]<0.05,t+1]=0
    states$G[k,t+1] = pmin(pmax(states$G[k,t+1],0), states$curr_max[k,t]); #states$G[states$G[k,t+1]<0.05,t+1]=0
    states$H[k,t+1] = pmin(pmax(states$H[k,t+1],0), states$curr_max[k,t]); #states$H[states$H[k,t+1]<0.05,t+1]=0
    states$Test[k,t+1] = pmin((test_rate * states$Test[k,t]),(states$J[k,t+1] +states$S[k,t+1]), testbound, na.rm = TRUE)
    states$V[k,t+1] = pmin(pmax(states$V[k,t+1],0), states$curr_max[k,t]); #states$V[states$V[,t+1]<0.05,t+1]=0
    states$IV[k,t+1] = pmin(pmax(states$IV[k,t+1],0), states$curr_max[k,t]); #states$IV[states$IV[,t+1]<0.05,t+1]=0
    
  }  
  return(states)
  
  ##update states
  ##return the updated values based on the model for the state P-A-D (observed )
}


obj_fn_cost_30 = function(interv_para){
  kappa =  c(interv_para[1],interv_para[2])
  kappa_trans = exp(kappa)/(1 + exp(kappa))
  # threshold = interv_para[3]
  threshold = interv_para[3]
  threshold_trans = exp(interv_para[3])/(1 + exp(interv_para[3]))
  intervention_para = list(kappa = kappa_trans, threshold = threshold_trans)
  #print(interv_para)
  ans = init.states
  for(t in 1: (T-1)){
    states = ans
    ans = update_states(t, states,epi_para, migration, test_rate, vac_rate, t_vac, intervention_para, vac_efficacy, testbound)
  }
  dth = ans$D
  cst = ans$cost
  
  cond = sapply(1:10, function(k) (max(states$H[k,])/states$S[k,1]) > 0.005 )
  cost_lockdown = sum(sapply(1:N, function(k){
    #dth[k,T]/sum(states$S[,1])+ sum(cst[k,])/(sum(states$S[,1])*365*30)
    30*dth[k,T]/sum(states$S[,1])+ sum(cst[k,])/(sum(states$S[,1])*365)
  }))
  
  cond_mult =  1-prod(1-cond)
  #cost_lockdown = (cost_lockdown + cost_lockdown*cond_mult)
  cost_lockdown = (cost_lockdown*10^3 + cost_lockdown*10^8 *cond_mult)
  
  #print(cost_lockdown)
  return(cost_lockdown)
}


obj_fn_cost_1 = function(interv_para){
  kappa =  c(interv_para[1],interv_para[2])
  kappa_trans = exp(kappa)/(1 + exp(kappa))
  # threshold = interv_para[3]
  threshold = interv_para[3]
  threshold_trans = exp(interv_para[3])/(1 + exp(interv_para[3]))
  intervention_para = list(kappa = kappa_trans, threshold = threshold_trans)
  #print(interv_para)
  ans = init.states
  for(t in 1: (T-1)){
    states = ans
    ans = update_states(t, states,epi_para, migration, test_rate, vac_rate, t_vac, intervention_para, vac_efficacy, testbound)
  }
  dth = ans$D
  cst = ans$cost
  
  cond = sapply(1:10, function(k) (max(states$H[k,])/states$S[k,1]) > 0.005 )
  cost_lockdown = sum(sapply(1:N, function(k){
    #dth[k,T]/sum(states$S[,1])+ sum(cst[k,])/(sum(states$S[,1])*365*30)
    #dth[k,T]/sum(states$S[,1])+ sum(cst[k,])/(sum(states$S[,1])*365)
    dth[k,T]/sum(states$S[,1])+ sum(cst[k,])/(sum(states$S[,1])*365)
  }))
  
  cond_mult =  1-prod(1-cond)
  #cost_lockdown = (cost_lockdown + cost_lockdown*cond_mult)
  cost_lockdown = (cost_lockdown*10^3 + cost_lockdown*10^8 *cond_mult)
  
  #print(cost_lockdown)
  return(cost_lockdown)
}


