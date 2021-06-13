library(tidyverse)

## length of each lockdown
## cost and relative cost of each lockdown of each node

node = 1:10
T=400


len_and_cost_lockdown <- function(final_state, herd_imm = FALSE){
  if(herd_imm==T){print('no lockdowns and costs in herd immunity')}else{
    
    lockdown.period <- list()
    final.state.lockdown <- sapply(1:T, function(t){
      cond <- final_state$H[,t]/final_state$S[,1] > interv_para_opt$threshold | final_state$I[,t]/final_state$S[,1] > .004
      return(cond)
    })
    for(i in 1:nrow(final.state.lockdown)){
      node.lock <- which(final.state.lockdown[i,])
      ind.start.lock <- which(diff(node.lock)>1)+1
      lockdown.period[[i]] <- as.data.frame(
        cbind(node.lock[c(1,ind.start.lock)+1],
              # ifelse(is.na(node.lock[c(1,ind.start.lock)+1]),400,
              #            node.lock[c(1,ind.start.lock)+1]),
              c(node.lock[c(ind.start.lock-1)]+1,node.lock[length(node.lock)]))) %>%
        rename(start.day=V1, end.day=V2) %>%
        mutate(length = end.day-start.day+1) %>%
        na.omit()
      
    }
    
    for(j in 1:length(node)){
      
      lockdown.period[[j]]$cost <- sapply(1:nrow(lockdown.period[[j]]), function(s){
        
        daily_cost <- final_state$cost[j,] 
        return(
          sum(
            daily_cost[(lockdown.period[[j]][s,1]):(lockdown.period[[j]][s,2])]))
        
      })
      
      lockdown.period[[j]]$relative.econ.cost <- lockdown.period[[j]]$cost/(sum(final_state$S[,1])*365)*365/400 #sum(final_state$S[,1]*365) # DP: remove 30
      
      
    }
    
    return(list(len_cost_lockdown = lockdown.period))}}


# ss2=len_and_cost_lockdown(final_state=final_dyn, herd_imm = F)
# ss2$len_cost_lockdown[[1]]


#sum(final_dyn$cost)/(sum(states$S[,1])*365*30)
