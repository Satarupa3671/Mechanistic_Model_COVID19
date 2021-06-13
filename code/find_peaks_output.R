## return the peaks of H D and daily Death (delD) in .txt files

find_peaks <- function(final_state){
  
  
  h_peak = c(max(final_dyn$H[1,]),
             max(final_dyn$H[7,]))
  d_peak = c(max(final_dyn$D[1,]),
             max(final_dyn$D[7,]))
  d_daily_peak = c(max(diff(final_dyn$D[1,])),
                   max(diff(final_dyn$D[7,])))
  print(c(paste("H_peaks_node1 = ", h_peak[1]),
          paste("H_peaks_node7 = ", h_peak[2]),
          paste("D_peaks_node1 = ", d_peak[1]),
          paste("D_peaks_node7 = ", d_peak[2]),
          paste("D_daily_peaks_node1 = ", d_daily_peak[1]),
          paste("D_daily_peaks_node7 = ", d_daily_peak[2])))
}
