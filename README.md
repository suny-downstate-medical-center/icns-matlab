# Intrinsic Cardiac Nervous System - Principal Neuron model - MATLAB
Version 0.1

# Description
Principle Neuron model adapted to run in MATLAB, with graphical output of voltage v. time. 

# Contents
LoadInitialConditions1.m; 
PN_model.m; 
plot_results.m;
README.md; 
example_plot.png

# Usage
## Download 3 scripts from repository: 
    LoadInitialConditions1.m, PN_model.m, plot_results.m 

## From MATLAB command line, enter [y0,t] = LoadInitialConditions
    Outputs: y0 = vector of voltages at time t0; t = t0

## From command line, enter dydt = PN_model(t, y0)
    Inputs: y0 = vector of initial voltages; t0 = initial time;  
    Outputs: dydt: vector length y0; 
    Functionality: Runs ODE
    
## From command line, enter plot_results
    Functionality: generates plot of voltage (mV) v time (s) (example_plot.png)
