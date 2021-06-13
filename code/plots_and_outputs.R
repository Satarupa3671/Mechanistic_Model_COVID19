
## dual y-axis
## show daily death -> max is less than 100, hard to see in any scales of other states
## so I keep the previous

## 051821
## make states are displayed in the same panels for both nodes accordingly 
## DHI + S(IV)J + RQT
## DHIQ + S(IV) + RJT

library(ggrepel)
library(RColorBrewer)

color <- c(brewer.pal(8,'Dark2'), 'lightpink1')
color[5] <- '#5AA6BD'




plot_func = function(final_state, vac_rate, herd_imm = FALSE){
  
  node = c(2,1,7)
  ## will create one new object to store the lockdown period for two nodes
  ###
  lockdown.period <- list()
  final.state.lockdown <- sapply(1:T, function(t){
    cond <- final_state$H[,t]/final_state$S[,1] > interv_para_opt$threshold | final_state$I[,t]/final_state$S[,1] > .004
    return(cond[node])
  })
  
  for(i in 1:length(node)){
    node.lock <- which(final.state.lockdown[i,])
    ind.start.lock <- which(diff(node.lock)>1)+1
    
    lockdown.period[[i]] <- as.data.frame(
      cbind(node.lock[c(1,ind.start.lock)+1],
            c(node.lock[c(ind.start.lock-1)]+1,node.lock[length(node.lock)]))) %>%
      rename(start.day=V1, end.day=V2) %>%
      mutate(length = end.day-start.day+1)
    
  }
  
  pop_day_one <- sapply(node, function(i){
    return(final_dyn$D[i,1] + final_dyn$H[i,1] +
             final_dyn$I[i,1] + final_dyn$J[i,1] + 
             final_dyn$Q[i,1] + final_dyn$R[i,1] + 
             final_dyn$S[i,1] + final_dyn$IV[i,1] +
             final_dyn$G[i,1])})
  ###

    
    df1 <- data.frame(Days=rep(1:T,9)) %>%
      mutate(value = c(final_state$I[node[1],],final_state$H[node[1],],
                       final_state$D[node[1],],#epi_para$delta[1]*final_state$H[1,], 
                       final_state$Q[node[1],],
                       final_state$S[node[1],],final_state$IV[node[1],],
                       final_state$J[node[1],],final_state$R[node[1],],final_state$Test[node[1],]),
             states = c(rep('I',T),rep('H',T),
                        rep('D',T),#rep('delD',T), 
                        rep('Q',T),rep('S',T), rep('IV',T), 
                        rep('J',T), rep('R',T), rep('Test',T)),
             type = c(rep("state_type1_node2",4*T),rep("state_type2_node2",2*T),
                      rep("state_type3_node2",3*T)))#,
    #node = rep("Node1_type1", 8*T))
    df1$states = as.factor(df1$states)
    
    df2 <- data.frame(Days=rep(1:T,9)) %>%
      mutate(value = c(final_state$I[node[2],],final_state$H[node[2],],
                       final_state$D[node[2],],#epi_para$delta[7]*final_state$H[7,], 
                       final_state$Q[node[2],],
                       final_state$S[node[2],],final_state$IV[node[2],],
                       final_state$J[node[2],],final_state$R[node[2],],final_state$Test[node[2],]
                       ),
             states = c(rep('I',T),rep('H',T),
                        rep('D',T),#rep('delD',T), 
                        rep('Q',T),rep('S',T), rep('IV',T), 
                        rep('J',T), rep('R',T), rep('Test',T)),
             type = c(rep("state_type1_node1",4*T),rep("state_type2_node1",2*T),
                      rep("state_type3_node1",3*T)))#,
    #  node = rep("Node7",8*T))
    df2$states = as.factor(df2$states)
    
    df3 <- data.frame(Days=rep(1:T,9)) %>%
      mutate(value = c(final_state$I[node[3],],final_state$H[node[3],],
                       final_state$D[node[3],],#epi_para$delta[7]*final_state$H[7,], 
                       final_state$Q[node[3],],
                       final_state$S[node[3],],final_state$IV[node[3],],
                       final_state$J[node[3],],final_state$R[node[3],],final_state$Test[node[3],]
      ),
      states = c(rep('I',T),rep('H',T),
                 rep('D',T),#rep('delD',T), 
                 rep('Q',T),rep('S',T), rep('IV',T), 
                 rep('J',T), rep('R',T), rep('Test',T)),
      type = c(rep("state_type1_node7",4*T),rep("state_type2_node7",2*T),
               rep("state_type3_node7",3*T)))#,
    #  node = rep("Node7",8*T))
    df3$states = as.factor(df3$states)
  
  
  ## 
  if(herd_imm == FALSE){
    
    p1 <- ggplot(df1, aes(x=Days,y=value,col=states,group=states)) + 
      geom_line(size=1.2,aes(linetype=states))+
      labs(#title = paste('node ',k), 
        x = "Days") +theme_bw()+
      theme(axis.text=element_text(size=12),
            axis.title.y = element_text(size=18),
            axis.title=element_text(size=12),
            strip.text = element_text(size=18),
            legend.title = element_text(size=15),
            legend.text=element_text(size=15),
            legend.key.size = unit(3,"line")) + 
      facet_wrap(~ type, ncol = 3, dir="v",
                 labeller = labeller(type = 
                                       c("state_type1_node2" = "Node 2",
                                         "state_type2_node2" = "Node 2",
                                         "state_type3_node2" = "Node 2")),
                 scales = "free") +
      scale_linetype_manual(values = c("solid", "22", "11", "4211", "13", "1343", "73", "2262", "3152"))+
      scale_color_manual(values = color) +
      annotate("rect", xmin = lockdown.period[[1]]$start.day, 
               xmax = lockdown.period[[1]]$end.day, ymin = -Inf, ymax = Inf,
               alpha = .15) +
      
      # Custom the Y scales:
      scale_y_continuous(
        
        # Features of the first axis
        name = "Numbers",
        
        # Add a second axis and specify its features
        sec.axis = sec_axis( trans=~./pop_day_one[1]*100, name="Percent",
                             labels = function(b){paste0(b,'%')})
      ) +
      ### remove x-axis label and reduce bottom margin
      theme(
        #axis.text.x = element_blank(),
        axis.title.x = element_blank()
        #axis.ticks.x = element_blank()
      ) 
    
    
    p2 <- ggplot(df2, aes(x=Days,y=value,col=states,group=states)) + 
      geom_line(size=1.2,aes(linetype=states))+
      labs(#title = paste('node ',k), 
        x = "Days") +theme_bw()+
      theme(axis.text=element_text(size=12),
            axis.title.y = element_text(size=18),
            axis.title=element_text(size=12),
            strip.text = element_text(size=18),
            legend.title = element_text(size=15),
            legend.text=element_text(size=15),
            legend.key.size = unit(3,"line")) + 
      facet_wrap(~ type, ncol = 3, dir="v",
                 labeller = labeller(type = 
                                       c("state_type1_node1" = "Node 1",
                                         "state_type2_node1" = "Node 1",
                                         "state_type3_node1" = "Node 1")),
                 scales = "free") +
      annotate("rect", xmin = lockdown.period[[2]]$start.day, 
               xmax = lockdown.period[[2]]$end.day, ymin = -Inf, ymax = Inf,
               alpha = .15)+
      scale_linetype_manual(values = c("solid", "22", "11", "4211", "13", "1343", "73", "2262", "3152"))+
      scale_color_manual(values = color)+
      
      # Custom the Y scales:
      scale_y_continuous(
        
        # Features of the first axis
        name = "Numbers",
        
        # Add a second axis and specify its features
        sec.axis = sec_axis( trans=~./pop_day_one[2]*100, name="Percent",
                             labels = function(b){paste0(b,'%')})
      ) 
    
    
    p3 <- ggplot(df3, aes(x=Days,y=value,col=states,group=states)) + 
      geom_line(size=1.2,aes(linetype=states))+
      labs(#title = paste('node ',k), 
        x = "Days") +theme_bw()+
      theme(axis.text=element_text(size=12),
            axis.title.y = element_text(size=18),
            axis.title=element_text(size=12),
            strip.text = element_text(size=18),
            legend.title = element_text(size=15),
            legend.text=element_text(size=15),
            legend.key.size = unit(3,"line")) + 
      facet_wrap(~ type, ncol = 3, dir="v",
                 labeller = labeller(type = 
                                       c("state_type1_node7" = "Node 7",
                                         "state_type2_node7" = "Node 7",
                                         "state_type3_node7" = "Node 7")),
                 scales = "free") +
      annotate("rect", xmin = lockdown.period[[3]]$start.day, 
               xmax = lockdown.period[[3]]$end.day, ymin = -Inf, ymax = Inf,
               alpha = .15)+
      scale_linetype_manual(values = c("solid", "22", "11", "4211", "13", "1343", "73", "2262", "3152"))+
      scale_color_manual(values = color)+
      
      # Custom the Y scales:
      scale_y_continuous(
        
        # Features of the first axis
        name = "Numbers",
        
        # Add a second axis and specify its features
        sec.axis = sec_axis( trans=~./pop_day_one[3]*100, name="Percent",
                             labels = function(b){paste0(b,'%')})
      ) 
    
    
    (p1) / (p2) /(p3) +
      plot_layout(guides = "collect")&theme(legend.position = 'right')
  } else {
    
    p1 <-  ggplot(df1, aes(x=Days,y=value,col=states,group=states)) + 
      geom_line(size=1.2,aes(linetype=states))+
      labs(#title = paste('node ',k), 
        x = "Days") +theme_bw()+
      theme(axis.text=element_text(size=12),
            axis.title.y = element_text(size=18),
            axis.title=element_text(size=12),
            strip.text = element_text(size=18),
            legend.title = element_text(size=15),
            legend.text=element_text(size=15),
            legend.key.size = unit(3,"line")) + 
      facet_wrap(~ type, ncol = 3, dir="v",
                 labeller = labeller(type = 
                                       c("state_type1_node2" = "Node 2",
                                         "state_type2_node2" = "Node 2",
                                         "state_type3_node2" = "Node 2")),
                 scales = "free") +
      scale_linetype_manual(values = c("solid", "22", "11", "4211", "13", "1343", "73", "2262", "3152"))+
      scale_color_manual(values = color) +
      # Custom the Y scales:
      scale_y_continuous(
        
        # Features of the first axis
        name = "Numbers",
        
        # Add a second axis and specify its features
        sec.axis = sec_axis( trans=~./pop_day_one[1]*100, name="Percent",
                             labels = function(b){paste0(b,'%')})
      ) +
      ### remove x-axis label and reduce bottom margin
      theme(
        #axis.text.x = element_blank(),
        axis.title.x = element_blank()
        #axis.ticks.x = element_blank()
      ) 
    
    p2 <-  ggplot(df2, aes(x=Days,y=value,col=states,group=states)) + 
      geom_line(size=1.2,aes(linetype=states))+
      labs(#title = paste('node ',k), 
        x = "Days") +theme_bw()+
      theme(axis.text=element_text(size=12),
            axis.title.y = element_text(size=18),
            axis.title=element_text(size=12),
            strip.text = element_text(size=18),
            legend.title = element_text(size=15),
            legend.text=element_text(size=15),
            legend.key.size = unit(3,"line")) + 
      facet_wrap(~ type, ncol = 3, dir="v",
                 labeller = labeller(type = 
                                       c("state_type1_node1" = "Node 1",
                                         "state_type2_node1" = "Node 1",
                                         "state_type3_node1" = "Node 1")),
                 scales = "free") +
      scale_linetype_manual(values = c("solid", "22", "11", "4211", "13", "1343", "73", "2262", "3152"))+
      scale_color_manual(values = color) +
      # Custom the Y scales:
      scale_y_continuous(
        
        # Features of the first axis
        name = "Numbers",
        
        # Add a second axis and specify its features
        sec.axis = sec_axis( trans=~./pop_day_one[2]*100, name="Percent",
                             labels = function(b){paste0(b,'%')})
      ) +
      ### remove x-axis label and reduce bottom margin
      theme(
        #axis.text.x = element_blank(),
        axis.title.x = element_blank()
        #axis.ticks.x = element_blank()
      ) 
    
    
    p3 <- ggplot(df3, aes(x=Days,y=value,col=states,group=states)) + 
      geom_line(size=1.2,aes(linetype=states))+
      labs(#title = paste('node ',k), 
        x = "Days") +theme_bw()+
      theme(axis.text=element_text(size=12),
            axis.title.y = element_text(size=18),
            axis.title=element_text(size=12),
            strip.text = element_text(size=18),
            legend.title = element_text(size=15),
            legend.text=element_text(size=15),
            legend.key.size = unit(3,"line")) + 
      facet_wrap(~ type, ncol = 3, dir="v",
                 labeller = labeller(type = 
                                       c("state_type1_node7" = "Node 7",
                                         "state_type2_node7" = "Node 7",
                                         "state_type3_node7" = "Node 7")),
                 scales = "free") +
      scale_linetype_manual(values = c("solid", "22", "11", "4211", "13", "1343", "73", "2262", "3152"))+
      scale_color_manual(values = color)+
      
      # Custom the Y scales:
      scale_y_continuous(
        
        # Features of the first axis
        name = "Numbers",
        
        # Add a second axis and specify its features
        sec.axis = sec_axis( trans=~./pop_day_one[3]*100, name="Percent",
                             labels = function(b){paste0(b,'%')})
      ) 
    (p1) / (p2)/(p3) +
      plot_layout(guides = "collect")&theme(legend.position = 'right')
  }
  
}


summary_func = function(final_states, interv_para_opt){
  print(paste("optimal kappa = ", interv_para_opt[1],
              "optimal threshold parameter = ", interv_para_opt[2],
              "optimal vaccination rate = ", interv_para_opt[3]))
  print(paste("final mortality rate =", round(sum(final_states$D[,T])/sum(final_states$S[,1]),5)))
  print(paste("final recovery rate = ", round(sum(final_states$R[,T] + final_states$IV[,T])/sum(final_states$S[,1]),5)))
}











