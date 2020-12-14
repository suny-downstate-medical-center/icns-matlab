# Intrinsic Cardiac Nervous System-Right Atrial Ganglionic Plexus: Principal Neuron model - MATLAB

# Description
Principle Neuron model adapted to run in MATLAB, with graphical output of voltage v. time. 

# Contents
LoadInitialConditions1.m; 
PN_model.m; 
plot_results.m;
README.md; 
example_plot.png

# Usage
## Clone repository: 
    git clone https://github.com/suny-downstate-medical-center/icns-matlab/
  
## Open MATLAB, cd to icns-matlab folder, and add to path. Then, from MATLAB command line, enter [y0,t] = LoadInitialConditions1
    Outputs: y0 = vector of voltages at time t0; t = t0

## From command line, enter dydt = PN_model(t, y0)
    Inputs: y0 = vector of initial voltages; t0 = initial time;  
    Outputs: dydt: vector length y0; 
    Functionality: Runs ODEs
    
## From command line, enter plot_results
    Functionality: generates and saves plot of voltage (mV) v time (s) (see example_plot.png)
