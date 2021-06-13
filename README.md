# Mechanistic_Model_COVID19

----------------Description


We propose a multi-compartment dynamic model 
as a realistic emulator for the COVID-19 pandemic by taking into account the major factors that influence its progression. Specifically, our model incorporates various epidemiological factors, and different policy intervention options such as lockdown, and testing and vaccination rates. Furthermore, the model involves a network of nodes representing different population centers or strata, with possible migrations across nodes. The focus is on optimal policy decisions to 
minimize the impact of the pandemic, by making a meaningful assessment of the costs due to deaths and lockdowns, while also taking into account the capacity of the healthcare system. We compare several scenarios with different levels of intervention through extensive simulation studies. The results suggest that a high rate of testing coupled with rapid vaccination, can be very effective in bringing the pandemic under control quickly and economically.

----------------The relevant files on the code folder are as follows

1. parameters_initial_values.R -- defines the initial values for the epidemiological and control parameters. These can usually be estimated from the observed compartments of an epidemic.

2. optim_dynamics_module.R -- 
  (A) Function dynamics_module() -- a detailed and interpretable mechanistic modeling of a viral epidemic, explicitly considering a multi-compartment model through a set of difference equation is considered. The updating equations roughly capture the mean behaviour of the compartments in discrete time and can be interpreted as the current value of the compartments in the dynamics of the transmission. The evolution of the pandemic can be expressed mechanistically as follows. 
  (B) obj_fn() -- An overall lost function that minimizes the overall death rate, while penalizing for the economic cost of lockdown is described. 
  
3. plots_and_outputs.R -- for a given pandemic trajectories computes the relevant summary statistics and outputs the ggplots.

4. len_and_cost_lockdown.R -- outputs the length and economic costs of a lockdown for any given state/node.

5. find_peaks_output.R -- calculates the maximum values of the various compartments of the pandemic.

6. final_results.R -- making use of the above R functions, computes the optimal intervention strategy using Optim() function, outputs the plots and results for the final trajectories, for various combination of the control parameters (testing rate, vaccination, migration).

7. herd_imm.R -- computes the results for a herd immunity situation, if no intervention measure or restriction is implemented at all and the pandemic is allowed to run its course. 
8. main_article_plots.R -- generates all the plots and tables included in the main article.
9. supple_plots.R -- generates all the plots and tables included in the supplementary material to the main article.

----------------Results

1. The folder main_article_plots contains all the plots and tables included in the main article.
2. supple_plots -- generates all the plots and tables included in the supplementary material to the main article.
