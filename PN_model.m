
function dydt = PN_model(t, y)

%Parameters
%==========

Cm       = 0.025;   % membrane capacitance  (nano-farad)   */
gNa_bar  = 3     ;   % maximal conductance of Na  (micro-siemen) */
gDR_bar  = 0.9    ;  % maximal conductance of Kdr (delayed rectifier) (micro-siemen) */
gA_bar   = 0.15   ; % maximal conductance of KA (transient)- micro-siemen */
gAHP_bar = 0.15   ;% maximal conductance of KAHP(Ca dep K ) - (micro-siemen) */
gCaL_bar = 0.0015  ; % maximal conductance of Ca_L (micro-siemen) */
gCaN_bar = 0.002  ; % maximal conductance of Ca_N (micro-siemen) */
gL       = 0.05    ; % maximal conductance of leak channel (micro-siemen) */


ENa      = 55  ;     % reversal potential of Na+ (milli-volt)   */
EK       = -94 ;     % reversal potential of K+ (milli-volt)   */
ESynE    = -10  ;     % reversal potential of (excitatory Synapse) (milli-volt)  */
ESynI    =- 94  ;     % reversal potential of (inhibitory Synapse) (milli-volt)  */
EL       = -43 ;     % reversal potential of leak channel (milli-volt - (-41.31 for 1Hz; -52.57 for 2Hz; -52.56 for 3Hz))   */;


v        = 2.5e-4;  % Volume of the shell for Ca channel (nano-liter)
Btot     = 0.03   ;  % Total concentration of Calcium (free + bound) milli-molar
K        = 0.001  ;  % kb/kf (ratio of backward to forward rate constatnts (milli-molar)
Cai0     = 5e-5   ;  % Equilibrium Ca2+ concentration (milli-molar)
tauE     = 30     ;  % Time constatnt for excitatory synapse (milli-second)
tauI     = 30     ;  % Time constatnt for inhibitory synapse (milli-second)
                             
F        = 9.648e4;  % Parameter for inwardly rectofying K+ current (coulomb/mol)
T        = 308    ;  % Temp - parameter for inwardly rectofying K+ current (kelvin)
KR       = 1.5    ;  % Boltzmann binding constatnt for  cAMP modulaiton
KRcAMP   = 1e-3   ; % Boltzmann binding constatnt for  cAMP modulaiton (milli-molar)
DRcAMP   = 0.4e-3 ;  % Boltzmann binding constatnt for  cAMP modulaiton (milli-molar)
gR_bar   = (0.18/3); % maximal unmodulated conductance (micro-siemen)
Z        = 2       ; % Parameter for inwardly rectofying K+ current (coulomb/mol)
RT       = (8314*308); % Joule/Kmol
                             
                             
k_adc    = 0.6e-6;   % Kinetic parameter for cAMP modulation (mM/m-sec)
K_mod    = 1.5    ;  % Kinetic parameter for cAMP modulation
K_5HT    = 6e-3   ;  % Kinetic parameter for cAMP modulation (milli-molar)
v_pde    = 2.4e-6 ;  % Kinetic parameter for cAMP modulation (mM/m-sec)
K_pde    = 3e-3   ;  % Kinetic parameter for cAMP modulation (milli-molar)
                             

%variables
%=========

V     = y(1); % Membrane volatge from interg over Vm
mNa   = y(2); % activation variable of Na+
hNa   = y(3); % inactivation variable of Na+
mDR   = y(4); % activation variable of Kdr
mA1   = y(5); % activation variable of transient potassium K_A1
hA1   = y(6); % inactivation variable of transient potassium K_A1
mA2   = y(7); % activation variable of transient potassium K_A2
hA2   = y(8); % inactivation variable of transient potassium K_A2
mAHP  = y(9); % activation variable of Calcium dependent potassium K_AH
mCaL  = y(10); % activation variable of L-type calcium
mCaN  = y(11); % m1-activation variable of N-type calcium
hCaN1 = y(12); %h1-inactivation variable of N1-type calcium
hCaN2 = y(13); %h2-inactivation variable of N2-type calcium
Cai   = y(14); % Calcium concentration
gSynE = y(15); % pre synaptic excitation
gSynI = y(16); % pre synaptic inhibition
cAMP  = y(17); % cyclic AMP modulation


%Input currents
%==============
                              
u= [
 0.012
 0.00
 0.012
 0.0
];

InE  = u(1); % Excitatory Synapse input current
InI  = u(2); % Inhibitory Synapse input current
Iinj = u(3); % Injected input current
                              
HT5  = u(4); % Serotonin concentration


                              
                              
    function a1= tau_pump(x)
        a1  = (17.7*exp(x/35));
    end

    function b1= PB(x)
        b1       = (Btot/(x+Btot+K));
    end

    function c1= ECa(x)
        c1     = (13.27*log(4/x));
    end



% Current from channels
%=====================
                              
INa   = -gNa_bar*mNa^3*hNa*(ENa-V); %Fast sodium
IDR   = -gDR_bar*mDR^4*(EK-V);  %Potassium delayed rectifier
IA    = -gA_bar*( 0.6*mA1^4*hA1 + 0.4*mA2^4*hA2 )*(EK-V); %Transient potassium
IAHP  = -gAHP_bar*mAHP*mAHP*(EK-V);     %Calcium dependent potassium
FR    = 1+KR/( 1+exp((KRcAMP-cAMP)/DRcAMP) );   %cAMP regulation
IR    = gR_bar*FR*(V-EK+5.66)/( 1+exp((V-EK-15.3)*Z*F/RT) );    %Inwardly rectifying potassium, KAR
ICa   = gCaL_bar*mCaL^2*(ECa(Cai)-V);   %Calcium-L type
                              
ISynE = gSynE*(ESynE-V);   %Excitatory Synapse
ISynI = gSynI*(ESynI-V);   %Inhibitory Synapse
IL    = gL*(EL-V);         %Leak channel current
                              

%Steady state values
%=====================
                              
% Fast Sodium, Na_fast */

miNa = ( (V+38)/(1-exp(-(V+38)/5)) );
tmNa = 1/( 0.091*miNa+0.062*miNa*exp(-(V+38)/5) );
miNa = 0.091*miNa*tmNa;

hiNa = 0.016*exp(-(V+55)/15);
thNa = 1/( hiNa+2.07/(1+exp(-(V-17)/21)) );
hiNa = hiNa*thNa;

gNa  = gNa_bar*mNa^3*hNa;

mNa_dot = (miNa-mNa)/tmNa;
hNa_dot = (hiNa-hNa)/thNa;


% Potassium-delayed rectifier, K_DR */

miDR = ( 0.01*(V+45)/(1-exp(-(V+45)/5)) );
tmDR = 1/( miDR+0.17*exp(-(V+50)/40) );
miDR = miDR*tmDR;

gDR  = gDR_bar*mDR^4;

mDR_dot = (miDR-mDR)/tmDR;


% Transient Potassium-A, K_A */

miA1 = 1/( 1+exp(-(V+60)/8.5) );
tmA1 = 1/( exp((V+35.82)/19.69) + exp(-(V+79.69)/12.7) + 0.37 );
hiA1 = 1/( 1+exp((V+78)/6) );
thA1 = 1/( 1+exp((V+46.05)/5) + exp(-(V+238.4)/37.45) );
miA2 = 1/( 1+exp(-(V+36)/20) );
tmA2 = tmA1;
hiA2 = hiA1;
thA2 = thA1;

if V > -63
    thA1 = 19;
elseif V > -73 
    thA2 = 60;
end    
    
  
gA = gA_bar*( 0.6*mA1^4*hA1 + 0.4*mA2^4*hA2 );

mA1_dot = (miA1-mA1)/tmA1;
hA1_dot = (hiA1-hA1)/thA1;
mA2_dot = (miA2-mA2)/tmA2;
hA2_dot = (hiA2-hA2)/thA2;


% Calcium-dependent potassium, K_AHP */

miAHP = 1.25e8*Cai*Cai;
tmAHP = 1e3/(miAHP+2.5);
miAHP = miAHP*1e-3*tmAHP;

gAHP  = gAHP_bar*mAHP*mAHP;
mAHP_dot = (miAHP-mAHP)/tmAHP;



% Inwardly rectifying (Anomalous Rectifier), KR */

%FR = 1+KR/( 1+exp((KRcAMP-cAMP)/DRcAMP) );
FR = 1;
IR = gR_bar*FR*(V-EK+5.66)/( 1+exp((V-EK-15.3)*Z*F/RT) );


% High-threshold calcium, CaL */

miCaL = 1.6/( 1+exp(-0.072*(V-5)) );
tmCaL = 1/( miCaL + 0.02*(V-1.31)/(exp((V-1.31)/5.36)-1) );
miCaL = miCaL*tmCaL;

gCaL  = gCaL_bar*mCaL*mCaL;

mCaL_dot = (miCaL-mCaL)/tmCaL;

                              
% Low-threshold N-type calcium, CaN */
                              
miCaN = 1.0/( 1+ exp(-(V+20)/4.5) );
tmCaN = 0.364*exp(-(0.042^2)*(V+31)^2) +0.442 ;

hiCaN1 = 1.0/( 1+ exp(V+20)/25);
thCaN1 = 3.752*exp(-(0.0395^2)*(V+30)^2) +0.56;
                              
hiCaN2 = 0.2/( 1+ exp(-(V+40)/10)) +  1.0/( 1+ exp(-(V+20)/40));
thCaN2 = 25.2*exp(-(0.0275^2)*(V+40)^2) + 8.4;
                              
gCaN  = gCaN_bar*mCaN*(0.55*hCaN1+0.45*hCaN2);
                              
mCaN_dot = (miCaN-mCaN)/tmCaN;
hCaN1_dot = (hiCaN1-hCaN1)/thCaN1;
hCaN2_dot = (hiCaN2-hCaN2)/thCaN2;
                              
                              
% cAMP balance in the cell */
                              
cAMP_dot = 0;
%cAMP_dot = k_adc*( 1+K_mod*(HT5/(HT5+K_5HT)) ) - v_pde*(cAMP/(cAMP+K_pde));
  
                              
% Calcium concentration  - with buffer*/
ICaL  = gCaL*(ECa(Cai)-V); %Calcium-L type
ICaN  = gCaN*(ECa(Cai)-V); %Calcium-L type
ICa = ICaL + ICaN;
Cai_dot = ( ICa*(1-PB(Cai))/(2*F*v) ) + ( (Cai0-Cai)/tau_pump(V) );
                              
                              
% Synaptic conductances */

gSynE_dot = (InE-gSynE)/tauE;
gSynI_dot = (InI-gSynI)/tauI;

% Membrane potential */

V_dot = gNa*(ENa-V) + (gDR+gA+gAHP)*(EK-V)- IR + ICa + gL*(EL-V)+ gSynE*(ESynE-V) + gSynI*(ESynI-V) + Iinj;
                             
V_dot = V_dot/Cm;

%differential equations
%======================

dydt=[
 (1e3*V_dot)
 (1e3*mNa_dot)
 (1e3*hNa_dot)
 (1e3*mDR_dot)
 (1e3*mA1_dot)
 (1e3*hA1_dot)
 (1e3*mA2_dot)
 (1e3*hA2_dot)
 (1e3*mAHP_dot)
 (1e3*mCaL_dot)
 (1e3*mCaN_dot)
 (1e3*hCaN1_dot)
 (1e3*hCaN2_dot)
 (1e3*Cai_dot)
 (1e3*gSynE_dot)
 (1e3*gSynI_dot)
 (1e3*cAMP_dot)
   ];

  end
  
  
