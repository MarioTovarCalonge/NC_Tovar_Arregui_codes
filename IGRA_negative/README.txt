System requirements:

-GCC (checked with Apple LLVM version 7.0.2)
-R (checked with R version 3.4.1 (2017-06-30) -- "Single Candle")
-Python3 (checked with python version 3.6.8)
-Python3 scipy, matplotlib and numpy packages

Once python3, R and gcc are available, no explicit installation, or special hardware is needed. When running the master script as described below, C codes will be compiled and required R packages (survival, bbmle) will be installed automatically from CRAN, or locally. If some error raised please install manually 'survival' and 'bbmle' packages for R.

If no python suite is installed or some package is missing, installation can be done via terminal and pip3 
(checked with MacOS high sierra and Ubuntu 18.04LS)
-Python packages can be installed by writing: sudo pip3 install numpy, scipy, matplotlib

PLEASE REMEMBER THAT PYHON DEPENDENCIES WON'T INSTALL AUTOMATICALLY AND SHOULD BE MANUALLY INSTALLED.

##########################################################################################

Executing the following command from this folder:

bash Master_script.txt 3000 4.00 0.00 0.00 0.74

will produce an example of a simulation of a TB vaccine clinical trial formally equivalent to the one published on the main text for a 100-0 vaccine. Also, it will infer the three vaccine descriptors proposed in our manuscript using the methodology explained there. 

also
bash Master_script 3000 4.00 0.00 0.50 0.00

will produce an example of a simulation of a TB vaccine clinical trial formally equivalent to the one published on the main text for a 0-100 vaccine.

-Standard output will be written, by default in 'data' folder.

PLEASE DON'T CHANGE ANY FOLDER NAME, NOR 'data' NEITHER 'IGRA_negative'. Otherwise errors will be raised as scripts need the correct paths for properly working.

-The script's call arguments represent the following:

##########################################################################################

3000 : N = Number of individuals in each cohort 

4.00 : Follow-up period length (years).

0.00 0.00 0.74
or
0.00 0.50 0.00

are the ground-truth characteristics of the vaccine, assumed to be known a priori.

##########################################################################################


In the Master_script file, one block of input values are found:
- Epidemiological parameters: calibrated to reproduce the infection and disease incidence levels reported in the setting of Worcester, South Africa.

Beta 0.068527: Force of infection.
p 0.375000: Fraction of newly infected individuals that undergo fast progression to disease.
q 0.21: Protection of slow progressors against TB reinfections.

A Second block of parameters can be found hard-coded in several scripts that need them. In the example provided, these are:

r_f 0.972160: Inverse of incubation time for fast progressors (years^-1)
r_s 0.000750: Inverse of incubation time for slow progressors (years^-1)

Both of them stands for the rates at which individuals develop TB from latency, where the subscript f means fast latency and s, slow latency.

##########################################################################################

The Master_script works calling sequentially five different scripts, coded in C, R and python:

1st. trial_simulator.c : Simulate one stochastic trial in-silico, given its dimensions, (cohort size and duration), and the ground-truth characteristics of the vaccine, assumed to be known a priori.

2nd. discrete_times.c Discretise the times-to-infection and times-to-disease that the previous step produced. This is to take into account the fact that in an actual trial, time resolution for these times is typically limited (e.g. QFT/diagnostic tests ran once every 3 months)

3rd. Hazard_inf.R Algorithms to estimate, from the series of (discretised) times-to-infection and times-to-disease, the efficacy of the vaccine at preventing infection (epsilon_beta).

4rd. mle.R Algorithms to estimate, from the series of (discretised) times-to-infection and times-to-disease the transition rates "r" and "rv" from which the protection and at delaying fast_progression (epsilon_r) can be obtained. (as explained in the methods section of the the main text) 

5th. Merge.py The efficacy of the vaccine against disease is estimated (VE_dis, or POD ). Then, from the analytic expression that binds VE_dis to the three vaccine descriptors epsilon_beta, epsilon_r and epsilon_p, the last one is isolated and estimated using the values of the others. Finally, all the results of the simulations are merged to give one global estimator with CI for each one of the vaccine parameters.

##########################################################################################

Expected output: as a final result of the simulation+vaccine inference procedure, the algorithm will write one file:

-Final_values_N_T.txt

where N and T are the dimension of the cohort and the follow-up period (in years). Here the inferred values of epsilon_r, epsilon_p and epsilon_beta, along with their confidence intervals associated to the trial simulated are provided. 

Also if matplotlib is installed, one image containing the violin plots of the inferred parameters and the global estimators with its CI will pop.

This violin plot is formally equivalent to the ones included on the main text.

