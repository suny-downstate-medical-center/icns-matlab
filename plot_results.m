  
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

plot(tspan, y(:,1),'LineWidth',0.5);
H = gca;
hold on;
legend('V')
xlabel('Time (s)')
ylabel('Voltage (mV)')
hold off;





 
 
 






