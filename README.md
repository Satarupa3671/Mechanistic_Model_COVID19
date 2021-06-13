# Mechanistic_Model_COVID19

The relevant files on the code folder are as follows

1. parameters_initial_values.R -- defines the initial values for the epidemiological and control parameters. These can usually be estimated from the observed compartments of an epidemic.

2. optim_dynamics_module.R -- 
  (A) Function dynamics_module() -- a detailed and interpretable mechanistic modeling of a viral epidemic, explicitly considering a multi-compartment model through a set of difference equation is considered. The updating equations roughly capture the mean behaviour of the compartments in discrete time and can be interpreted as the current value of the compartments in the dynamics of the transmission. The evolution of the pandemic can be expressed mechanistically as follows. 
  (B) obj_fn() -- An overall lost function that minimizes the overall death rate, while penalizing for the economic cost of lockdown is described. 
  
3. plots_and_outputs.R -- for a given pandemic trajectories computes the relevant summary statistics and outputs the ggplots.

4. len_and_cost_lockdown.R -- outputs the length and economic costs of a lockdown for any given state/node.

5. find_peaks_output.R -- calculates the maximum values of the various compartments of the pandemic.

6. final_results.R -- making use of the above R functions, computes the optimal intervention strategy using Optim() function, outputs the plots and results for the final trajectories, for various combination of the control parameters (testing rate, vaccination, migration).

7. herd_imm.R -- computes the results for a herd immunity situation, if no intervention measure or restriction is implemented at all and the pandemic is allowed to run its course. 
