This is the information file for the code used in the paper "Adaptive use of information during growth can explain long-term effects of early-life experiences". 

Any comments or questions, please contact Sinead English: sineadenglish[at]cantab[dot]net

Last updated: 11 December 2015


### CODE ###
The following files are provided in the 'code' folder:

run_code.txt -- instructions to build the programmes and run them.

growthfun.h -- header file
growthfun.cpp -- file describing the functions used in the main code
growth_info.cpp -- file describing the backwards iterations to produce the optimal decision, reproductive value and information value arrays
run_sim_all.cpp -- file describing the forward simulations

plot_results.r -- R code for plotting results


### DATA ###
Simulation data used in the paper are provided in the 'data' folder

The following files are produced as a result of the dynamic optimization procedure, with the number following the title indicating that assigned to the combination of parameter values, or scenario considered (i.e. 1 = same predation in both environments; 2 = higher predation in low-food environment; 3 = higher predation in high-food environment): 

ConvergeVector -- a vector showing the difference
Decision -- the optimal foraging (alpha) array (across all size and belief states) 
FitnessAlpha -- reproductive value at maturity for all size, belief and alpha combinations (to calculate the fitness value of information)
FitnessValue -- reproductive value at maturity for all size and belief combinations, when alpha is optimized
InfoVal -- value of information for all size and belief combinations
Params -- list of parameters used for a particular scenario

The simulation data files are named using the following rule: 

[Whether results are just at maturity (All), or for each time-step of development (Ind)].[Scenario].[Mortality no (0) or yes (1)].[Starting prior (0.1, 0.5, or 0.9)].[Probability of finding food in time window (0.0, 0.1, 0.9)].[Age at which experimental manipulation started] : 

thus "SimOutputAll.1.0.0.1.0.0.0.txt" indicates the data from the end-point of the simulation (i.e. at maturity or death) for individuals in scenario 1, which did not experience mortality, had a starting prior of 0.1, and did not have experimentally altered foraging, hence the age at starting the foraging manipulation is 0; 

and "SimOutputInd.2.0.0.5.0.1.20.txt" is the age-specific output for individuals in scenario 2, which did not experience predation, had a starting prior of 0.5, a likelihood of finding food in the experimental window of 0.1 and the window started at age 20. 

