# Intrinsic Cardiac Nervous System - Principal Neuron model
Version 0.1

# Description
Principle Neuron model adapted to run in MATLAB, with graphical output of voltage v. time. 

# Contents
LoadInitialConditions1.m
PN_model.m
Plot_results.m
README_matlab.md
example_plot.png

# Usage
1. Download 3 scripts from repository: 
    LoadInitialConditions1.m
    PN_model.m 
    plot_results.m
2. From MATLAB command line, type [y0,t] = LoadInitialConditions1
    Outputs:
      y0 = vector of voltages at time t0
      t = t0
3. From command line, type dydt = PN_model(t, y0)
    Inputs:
      y0 = vector of initial voltages
      t0 = initial time 
    Outputs:
      dydt: vector length y0
    Functionality: 
      Runs ODE
4. From command line, type plot_results
    Functionality: generates plot of voltage (mV) v time (s) (example_plot.png)
