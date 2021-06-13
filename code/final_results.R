
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(dplyr)
library(ggplot2)
library(parallel)
#library(DMwR)
library(patchwork)
#library(nloptr)
#library(optimx)
source("./parameters_initial_values.R")
source("./optim_dynamics_module.R")
source("./plots_and_outputs.R")
source("./len_and_cost_lockdown.R")
source("./find_peaks_output.R")

###results included in the main article

vac_efficacy = 0.95#0.8
t_vac = 100
test_bound = rep(c(5000, 20000), each = 6)
vac = rep(rep(sprintf("set%d",c(1,4,5)), each = 2),2)
cost_weight = rep(c(30,1),12)

migration = migration_situation("set1")
test_rate = testrate_situation("set2")


para.init = c(0.1,.1,.0001)
para.init = log(para.init/(1 - para.init))

for(i in seq(1,11,2)){
  vac_rate = vaccine_situation(vac[i])
  testbound = test_bound[i]
  opt = optim(para.init, obj_fn_cost_30, gr = NULL,method = "BFGS")
  interv_para_opt_trans = exp(opt$par)/(1 + exp(opt$par))
  interv_para_opt =   list(kappa = c(interv_para_opt_trans[1],interv_para_opt_trans[2]) , 
                         threshold = interv_para_opt_trans[3])
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
  plot_func(final_dyn, vac_rate = vac_rate)
  ggsave(sprintf("./results/optim_12_scenarios/case_%d_plots.eps",i),
         width = 15, height = 10, units='in', device = cairo_ps)
  
  len_cost_lockdown = len_and_cost_lockdown(final_state=final_dyn, herd_imm = F)
  save(interv_para_opt, final_dyn, vac_rate, migration, 
       test_rate, vac_efficacy, len_cost_lockdown,
       file = sprintf("./results/optim_12_scenarios/case_%d_Results.Rda",i))
  
  sink(sprintf("./results/optim_12_scenarios/case_%d_capacity_peaks.txt",i))
  summary_func(final_dyn, interv_para_opt)
  cat("\n")
  find_peaks(final_dyn)
  cat("\n")
  print(sapply(1:10, function(k) (max(final_dyn$H[k,])/final_dyn$S[k,1]) > 0.005 )) #are the cond met
  print(sapply(1:10, function(k) (max(final_dyn$H[k,])/final_dyn$S[k,1]) )) # what are the max rates
  sink()
}



for(i in seq(2,12,2)){
  vac_rate = vaccine_situation(vac[i])
  testbound = test_bound[i]
  opt = optim(para.init, obj_fn_cost_1, gr = NULL,method = "BFGS")
  interv_para_opt_trans = exp(opt$par)/(1 + exp(opt$par))
  interv_para_opt =   list(kappa = c(interv_para_opt_trans[1],interv_para_opt_trans[2]) , 
                                threshold = interv_para_opt_trans[3])
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
  plot_func(final_dyn, vac_rate = vac_rate)
  ggsave(sprintf("./results/optim_12_scenarios/case_%d_plots.eps",i),
         width = 15, height = 10, units='in', device = cairo_ps)
  
  len_cost_lockdown = len_and_cost_lockdown(final_state=final_dyn, herd_imm = F)
  save(interv_para_opt, final_dyn, vac_rate, migration, 
       test_rate, vac_efficacy, len_cost_lockdown,
       file = sprintf("./results/optim_12_scenarios/case_%d_Results.Rda",i))
  
  sink(sprintf("./results/optim_12_scenarios/case_%d_capacity_peaks.txt",i))
  summary_func(final_dyn, interv_para_opt)
  cat("\n")
  find_peaks(final_dyn)
  cat("\n")
  print(sapply(1:10, function(k) (max(final_dyn$H[k,])/final_dyn$S[k,1]) > 0.005 )) #are the cond met
  print(sapply(1:10, function(k) (max(final_dyn$H[k,])/final_dyn$S[k,1]) )) # what are the max rates
  sink()
}



