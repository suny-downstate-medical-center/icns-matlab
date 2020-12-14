  
%PLOTTING
%========
% relevant outputs
% H = current axes
  format long eng
  options = odeset( 'RelTol', 1e-9, 'AbsTol', 1e-9 );


%Time steps- simulation
%========================

ti = 0;
dt = 0.5*10^(-3);       % time step (0.0005 s)
tf = 1;                 % 800 for step response

tspan = ti:dt:tf;

% load initial conditons
%========================

[y0, t] = LoadInitialConditions1;
y1 = abs(y0);


% run simuations
%===============

[t, y] = ode15s(@PN_model, tspan, y0, options);


%Plotting Results
%=====================
H = gcf;
hold on;
plot(tspan, y(:,1),'LineWidth',0.8, 'Color', [0 0 1]);
xlabel('Time (s)')
ylabel('Membrane Voltage (mV)')
hold off;

% save figure not in local folder so as to not overwrite example)
saveas(H, 'test_plot_matlab', 'png')



 
 
 






