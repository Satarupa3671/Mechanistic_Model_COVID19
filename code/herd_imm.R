rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr)
library(ggplot2)
library(parallel)
#library(DMwR)
library(patchwork)
library(nloptr)
library(optimx)
source("./parameters_initial_values.R")
source("./optim_dynamics_module.R")
source("./plots_and_outputs_051821.R")
source("./len_and_cost_lockdown.R")
source("./find_peaks_output.R")



vac_efficacy = 0.95#0.8
t_vac =10
testbound = 20000#rep(c(5000, 20000), each = 6)
migration = migration_situation("set1")
test_rate = testrate_situation("set2")
vac_rate = vaccine_situation("set1")

###herd Imm
#No Interv + high mig (set1) + no vacc (set1) -- setE
#No Interv + high mig (set1) +  vacc (set4) -- setF
 #: no vacc
interv_para_opt = list(kappa = c(1,1) , threshold = 5000)

final_dyn = list()
final_dyn = init.states
for(t in 1: (T-1)){
  states = final_dyn
  final_dyn = update_states(t, states = states, epi_para = epi_para,
                            migration = migration, test_rate = test_rate,
                            vac_rate = vac_rate, t_vac = t_vac,
                            intervention_para = interv_para_opt, 
                            vac_efficacy = vac_efficacy,
                            testbound = testbound)
}
plot_func(final_dyn, vac_rate = vac_rate, herd_imm = T)
ggsave("./herd_imm/no_vac_herd_imm.eps",
       width = 15, height = 10, units='in', device = cairo_ps)

len_cost_lockdown = len_and_cost_lockdown(final_state=final_dyn, herd_imm = T)
save(interv_para_opt, final_dyn, vac_rate, migration, 
     test_rate, vac_efficacy, len_cost_lockdown,
     file = "./herd_imm/no_vac_herd_imm.Rda")

sink("./herd_imm/no_vac_herd_imm.txt")
summary_func(final_dyn, interv_para_opt)
cat("\n")
find_peaks(final_dyn)
cat("\n")
sapply(1:10, function(k) (max(final_dyn$H[k,])/final_dyn$S[k,1]) > 0.005 ) #are the cond met
sapply(1:10, function(k) (max(final_dyn$H[k,])/final_dyn$S[k,1]) ) # what are the max rates
sink()

## high vacc
migration = migration_situation("set1") #:high mig
test_rate = testrate_situation("set2")
vac_rate = vaccine_situation("set5") #: high vacc
t_vac = 10
interv_para_opt = list(kappa = c(1,1) , threshold = 5000)



final_dyn = list()
final_dyn = init.states
for(t in 1: (T-1)){
  states = final_dyn
  final_dyn = update_states(t, states = states, epi_para = epi_para,
                            migration = migration, test_rate = test_rate,
                            vac_rate = vac_rate, t_vac = t_vac,
                            intervention_para = interv_para_opt, 
                            vac_efficacy = vac_efficacy,
                            testbound = testbound)
}
plot_func(final_dyn, vac_rate = vac_rate, herd_imm = T)
ggsave("./herd_imm/vac_herd_imm.eps",
       width = 15, height = 10, units='in', device = cairo_ps)

len_cost_lockdown = len_and_cost_lockdown(final_state=final_dyn, herd_imm = T)
save(interv_para_opt, final_dyn, vac_rate, migration, 
     test_rate, vac_efficacy, len_cost_lockdown,
     file = "./herd_imm/vac_herd_imm.Rda")

sink("./herd_imm/vac_herd_imm.txt")
summary_func(final_dyn, interv_para_opt)
cat("\n")
find_peaks(final_dyn)
cat("\n")
sapply(1:10, function(k) (max(final_dyn$H[k,])/final_dyn$S[k,1]) > 0.005 ) #are the cond met
sapply(1:10, function(k) (max(final_dyn$H[k,])/final_dyn$S[k,1]) ) # what are the max rates
sink()


##Output
{
  vac_rate = vaccine_situation("set5")
  load("./herd_imm/vac_herd_imm.Rda")
  print(round(final_dyn$D[,T])) #Max D
  print(paste0("sum_D= ", round(sum(final_dyn$D[,T])))) #Max D
  cat("\n")
  print(round(
    sapply(1:N, function(k) max(final_dyn$H[k,]) ))) #Max H
  print(paste0("sum_H= ", round(
    sum(
      sapply(1:N, function(k) max(final_dyn$H[k,])) )))) 
  cat("\n")
  print(round(
    sapply(1:N, function(k) max(final_dyn$Q[k,]) ))) #Max Q
  print(paste0("sum_Q= ", round(
    sum(
      sapply(1:N, function(k) max(final_dyn$Q[k,])) )))) 
  cat("\n")
  print(round(
    sapply(1:N, function(k) max(final_dyn$J[k,]) ))) #Max J
  print(paste0("sum_J= ", round(
    sum(
      sapply(1:N, function(k) max(final_dyn$J[k,])) )))) 
  cat("\n")
  print(round(
    sapply(1:N, function(k) max(final_dyn$I[k,]) ))) #Max I
  print(paste0("sum_I= ", round(
    sum(
      sapply(1:N, function(k) max(final_dyn$I[k,])) )))) 
  cat("\n")
  print(round(
    sapply(1:N, function(k) max(final_dyn$IV[k,]) ))) #Max IV
  print(paste0("sum_IV= ", round(
    sum(
      sapply(1:N, function(k) max(final_dyn$IV[k,])) )))) 
  cat("\n")
  print(round(
    sapply(1:N, function(k) max(final_dyn$R[k,]) ))) #Max R
  print(paste0("sum_R= ", round(
    sum(
      sapply(1:N, function(k) max(final_dyn$R[k,])) )))) 
  
  ##
  cat("\n")
  print(paste0("opt= ", interv_para_opt)) 
  #economic cost  (annualized) 
  econ_cost_lockdown = sum(sapply(1:N, function(k){
    (sum(final_dyn$cost[k,])/(sum(final_dyn$S[,1])*365)) *(365/400)#/10^3#/10^3 *365/400
  }))
  cat("\n")
  print(paste0("econ_cost= ", round(econ_cost_lockdown,4))) 
  interv_para_opt = unlist(interv_para_opt)
  print(sapply(1:10, function(k) (max(final_dyn$H[k,])/final_dyn$S[k,1]) > 0.005 )) #are the cond met
  print(sapply(1:10, function(k) (max(final_dyn$H[k,])/final_dyn$S[k,1]) )) # what are the max rates
  cost_lockdown = sum(sapply(1:N, function(k){
    final_dyn$D[k,T]/sum(final_dyn$S[,1])
  }))
 print(paste0("cost= ", cost_lockdown))
}






