

#### 6/10
#### more scenarios in supplmentary
#### scenario S1-S4: case c(1,5,7,11) 
#### corresponding to scenario 3-6 in main body but with more weight on death (obj_fun_1)
#### scenario S5: case 12 but with migration set3 (restrictive); 



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
source("./plots_and_outputs.R") ## 3 nodes: node (2,1,7)
source("./len_and_cost_lockdown.R")
source("./find_peaks_output.R")



### run the scenario S5 for results
### results are saved in ../results/supple_plots/
############

## setting:
vac_efficacy = 0.95#0.8
t_vac = 100
testbound = 20000
vac_rate = vaccine_situation('set5')
cost_weight = 1
migration = migration_situation("set3")
test_rate = testrate_situation("set2")
para.init = c(0.1,.1,.0001)
para.init = log(para.init/(1 - para.init))

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


len_cost_lockdown = len_and_cost_lockdown(final_state=final_dyn, herd_imm = F)
save(interv_para_opt, final_dyn, vac_rate, migration, 
     test_rate, vac_efficacy, len_cost_lockdown,
     file = sprintf("../results/supple_plots/Scenario_S%d_Results.Rda",5))

sink(sprintf("../results/supple_plots/Scenario_S%d_capacity_peaks.txt",5))
summary_func(final_dyn, interv_para_opt)
cat("\n")
find_peaks(final_dyn)
cat("\n")
print(sapply(1:10, function(k) (max(final_dyn$H[k,])/final_dyn$S[k,1]) > 0.005 )) #are the cond met
print(sapply(1:10, function(k) (max(final_dyn$H[k,])/final_dyn$S[k,1]) )) # what are the max rates
sink()



#############



##### figures
###########

## scenario S1-S4 (from /initA: case 1, 5, 7, 11 actually)
case = c(1,5,7,11)
for(s in 1:length(case)){
  
  load(sprintf("../results/optim_12_scenarios/case_%d_Results.Rda",case[s]))
  plot_func(final_dyn, cav_rate,herd_imm = F)
  ggsave(sprintf("../results/supple_plots/case_%d_scen_S%d.eps",case[s],s),
         width = 15, height = 10, units='in', device = cairo_ps)
  
  
}


## scenario S5
load(sprintf("../results/supple_plots/Scenario_S%d_Results.Rda",5))
plot_func(final_dyn, vac_rate = vac_rate)
ggsave(sprintf("../results/supple_plots/scen_S%d.eps",5),
       width = 15, height = 10, units='in', device = cairo_ps)

###########




##### tables
###########

## table 1 -- report maxH(vul (%), rob (%)), maxQ(vul (%), rob (%)), 
## maxD(vul (%), rob (%)) 5 scenarios \times 3 columns
######
sink("../results/supple_plots/table_1.txt")
case = c(1,5,7,11)

# scen S1-S4
for(i in 1:4){
  
  print(sprintf("###### start scenario S%d", i))
  cat("\n")

  load(sprintf("../results/optim_12_scenarios/case_%d_Results.Rda",case[i]))
  
  pop_day_one <- sapply(node, function(i){
    return(final_dyn$D[i,1] + final_dyn$H[i,1] +
             final_dyn$I[i,1] + final_dyn$J[i,1] + 
             final_dyn$Q[i,1] + final_dyn$R[i,1] + 
             final_dyn$S[i,1] + final_dyn$IV[i,1] +
             final_dyn$G[i,1])})
  
  maxH_each = round(
    sapply(1:N, function(k) max(final_dyn$H[k,]) )) #Max H
  vul.ind = which.max(maxH_each[1:5])
  rob.ind = which.max(maxH_each[6:10])+5
  print(paste0("maxH = ", "[ ", maxH_each[vul.ind], '(', 
               round(maxH_each[vul.ind]/pop_day_one[vul.ind]*100,2), '%)',' , ',
               maxH_each[rob.ind], '(', round(maxH_each[rob.ind]/pop_day_one[rob.ind]*100,2), '%)',' ]'))
  cat("\n")
  
  maxD_each = round(
    sapply(1:N, function(k) max(final_dyn$D[k,]) )) #Max H
  vul.ind = which.max(maxD_each[1:5])
  rob.ind = which.max(maxD_each[6:10])+5
  print(paste0("maxD = ", "[ ", maxD_each[vul.ind], '(', 
               round(maxD_each[vul.ind]/pop_day_one[vul.ind]*100,2), '%)',' , ',
               maxD_each[rob.ind], '(', round(maxD_each[rob.ind]/pop_day_one[rob.ind]*100,2), '%)',' ]'))
  cat("\n")
  
  
  maxQ_each = round(
    sapply(1:N, function(k) max(final_dyn$Q[k,]) )) #Max H
  vul.ind = which.max(maxQ_each[1:5])
  rob.ind = which.max(maxQ_each[6:10])+5
  print(paste0("maxQ = ", "[ ", maxQ_each[vul.ind], '(', 
               round(maxQ_each[vul.ind]/pop_day_one[vul.ind]*100,2), '%)',' , ',
               maxQ_each[rob.ind], '(', round(maxQ_each[rob.ind]/pop_day_one[rob.ind]*100,2), '%)',' ]'))
  cat("\n")
  
  print(sprintf("###### end scenario S%d", i))
  cat("\n")
}

## scen S5
load(sprintf("../results/supple_plots/Scenario_S%d_Results.Rda",5))

print(sprintf("###### start scenario S%d", 5))
cat("\n")

pop_day_one <- sapply(node, function(i){
  return(final_dyn$D[i,1] + final_dyn$H[i,1] +
           final_dyn$I[i,1] + final_dyn$J[i,1] + 
           final_dyn$Q[i,1] + final_dyn$R[i,1] + 
           final_dyn$S[i,1] + final_dyn$IV[i,1] +
           final_dyn$G[i,1])})

maxH_each = round(
  sapply(1:N, function(k) max(final_dyn$H[k,]) )) #Max H
vul.ind = which.max(maxH_each[1:5])
rob.ind = which.max(maxH_each[6:10])+5
print(paste0("maxH = ", "[ ", maxH_each[vul.ind], '(', 
             round(maxH_each[vul.ind]/pop_day_one[vul.ind]*100,2), '%)',' , ',
             maxH_each[rob.ind], '(', round(maxH_each[rob.ind]/pop_day_one[rob.ind]*100,2), '%)',' ]'))
cat("\n")

maxD_each = round(
  sapply(1:N, function(k) max(final_dyn$D[k,]) )) #Max H
vul.ind = which.max(maxD_each[1:5])
rob.ind = which.max(maxD_each[6:10])+5
print(paste0("maxD = ", "[ ", maxD_each[vul.ind], '(', 
             round(maxD_each[vul.ind]/pop_day_one[vul.ind]*100,2), '%)',' , ',
             maxD_each[rob.ind], '(', round(maxD_each[rob.ind]/pop_day_one[rob.ind]*100,2), '%)',' ]'))
cat("\n")


maxQ_each = round(
  sapply(1:N, function(k) max(final_dyn$Q[k,]) )) #Max H
vul.ind = which.max(maxQ_each[1:5])
rob.ind = which.max(maxQ_each[6:10])+5
print(paste0("maxQ = ", "[ ", maxQ_each[vul.ind], '(', 
             round(maxQ_each[vul.ind]/pop_day_one[vul.ind]*100,2), '%)',' , ',
             maxQ_each[rob.ind], '(', round(maxQ_each[rob.ind]/pop_day_one[rob.ind]*100,2), '%)',' ]'))
cat("\n")

print(sprintf("###### end scenario S%d", 5))
cat("\n")

sink()
######

## table 2 -- report opt_para (vul, rob), eco_cost(vul (%), rob (%)), 
## total_lockdown_length(vul (%), rob (%)), total_death (vul (%), rob (%))

sink("../results/supple_plots/table_2.txt")

case = c(1,5,7,11)

# scen 3-6
for(i in 1:4){
  
  print(sprintf("###### start scenario S%d", i))
  cat("\n")

  load(sprintf("../results/optim_12_scenarios/case_%d_Results.Rda",case[i]))
  
  print(paste0("opt= ", interv_para_opt))
  cat("\n")
  
  #economic cost  (annualized)
  econ_cost_lockdown = sum(sapply(1:N, function(k){
    (sum(final_dyn$cost[k,])/(sum(final_dyn$S[,1])*365)) *(365/400)#/10^3#/10^3 *365/400
  }))
  print(paste0("econ_cost= ", round(econ_cost_lockdown,4)))
  cat("\n")
  
  # max total length of lockdown
  len_each = round(
    sapply(1:N, function(k) sum(len_cost_lockdown$len_cost_lockdown[[k]][,3]) )) 
  print(paste0("max_total_len = ", "[ ", max(len_each[1:5]),' , ',
               max(len_each[6:10]),' ]'))
  cat("\n")
  
  
  print(paste0("sum_D= ", round(sum(final_dyn$D[,T])))) 
  cat("\n")
  
  print(sprintf("###### end scenario S%d", i))
  cat("\n")
}

## scen S5
load(sprintf("../results/supple_plots/Scenario_S%d_Results.Rda",5))

print(sprintf("###### start scenario S%d", 5))
cat("\n")

print(paste0("opt= ", interv_para_opt))
cat("\n")

#economic cost  (annualized)
econ_cost_lockdown = sum(sapply(1:N, function(k){
  (sum(final_dyn$cost[k,])/(sum(final_dyn$S[,1])*365)) *(365/400)#/10^3#/10^3 *365/400
}))
print(paste0("econ_cost= ", round(econ_cost_lockdown,4)))
cat("\n")

# max total length of lockdown
len_each = round(
  sapply(1:N, function(k) sum(len_cost_lockdown$len_cost_lockdown[[k]][,3]) )) 
print(paste0("max_total_len = ", "[ ", max(len_each[1:5]),' , ',
             max(len_each[6:10]),' ]'))
cat("\n")


print(paste0("sum_D= ", round(sum(final_dyn$D[,T])))) 
cat("\n")

print(sprintf("###### end scenario S%d", 5))
cat("\n")

sink()

###########


