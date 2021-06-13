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



vac_efficacy = 0.95#0.8
t_vac = 100

migration = migration_situation("set1")
test_rate = testrate_situation("set2")




##### figures
###########
## scenario 1 (from /herd_imm)

s=1

load("../results/herd_imm/no_vac_herd_imm.Rda")

plot_func(final_dyn, vac_rate,herd_imm = T)
ggsave(sprintf("../results/main_article_plots/no_vac_herd_imm_scen_%d.eps",s),
       width = 15, height = 10, units='in', device = cairo_ps)


## scenario 2 (from /herd_imm)

s=2

load("../results/herd_imm/vac_herd_imm.Rda")

plot_func(final_dyn, vac_rate,herd_imm = T)
ggsave(sprintf("../results/main_article_plots/vac_herd_imm_scen_%d.eps",s),
       width = 15, height = 10, units='in', device = cairo_ps)

## scenario 3-6 (from /initA case 2, 6, 8, 12 actually)
case = c(2,6,8,12)
for(s in 1:length(case)){
  
  load(sprintf("../results/optim_12_scenarios/case_%d_Results.Rda",case[s]))
  plot_func(final_dyn, cav_rate,herd_imm = F)
  ggsave(sprintf("./../results/main_article_plots/case_%d_scen_%d.eps",case[s],s+2),
         width = 15, height = 10, units='in', device = cairo_ps)
  
  
}

###########




##### tables
###########

## table 1 -- report maxH(vul (%), rob (%)), maxQ(vul (%), rob (%)), 
## maxD(vul (%), rob (%)) 6 scenarios \times 3 columns
######
sink("../results/main_article_plots/table_1_test.txt")
case = c(2,6,8,12)

# scen 1
i=1
print(sprintf("###### start scenario %d", i))
cat("\n")

load("../results/herd_imm/no_vac_herd_imm.Rda")

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

print(sprintf("###### end scenario %d", i))
cat("\n")

# scen 2
i=2
print(sprintf("###### start scenario %d", i))
cat("\n")

load("../results/herd_imm/vac_herd_imm.Rda")

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

print(sprintf("###### end scenario %d", i))
cat("\n")

# scen 3-6
for(i in 3:6){
  
  print(sprintf("###### start scenario %d", i))
  cat("\n")

  load(sprintf("../results/optim_12_scenarios/case_%d_Results.Rda",case[i-2]))
  
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
  
  print(sprintf("###### end scenario %d", i))
  cat("\n")
}
sink()
######

## table 2 -- report opt_para (vul, rob), eco_cost(vul (%), rob (%)), 
## total_lockdown_length(vul (%), rob (%)), total_death (vul (%), rob (%))

sink("../results/main_article_plots/table_2.txt")

case = c(2,6,8,12)

# scen 3-6
for(i in 3:6){
  
  print(sprintf("###### start scenario %d", i))
  cat("\n")

  load(sprintf("../results/optim_12_scenarios/case_%d_Results.Rda",case[i-2]))
  
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
  
  print(sprintf("###### end scenario %d", i))
  cat("\n")
}
sink()

###########


