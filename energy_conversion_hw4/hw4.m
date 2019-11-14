% Eliminate all extraneous variables and plots
close all
clearvars

% First run with air temperature of 0°C
Tair = 0;
% Set initial conditions
IC = [0, 0, Tair];
% Generate 50 seconds of rest before first cycle
Idlc = 0;
[t0, V0] = ode45(@(t, V)odefcn(t, V, Idlc, Tair), [-50, 0], IC);
I0 = zeros(length(t0),1);
% t0, V0, and I0 refer to time, differential equation results, and 
% current at 0°C. They are constantly added to in each cycle.

% Run 5 cycles
for i = 1:5
    % Start with charging art of cycle
    Idlc = -100;
    % Use event function to stop when Vdlc reaches rated voltage.
    % The anonymous function @(t, V)event(t,V,Idlc) allows ode45
    % to interpret using it as a function of t and V, while still
    % using Idlc to calculate Vdlc.
    options = odeset('Events',@(t, V)event(t,V,Idlc));
    % Gain initial conditions from end of V0
    IC = V0(length(t0),:);
    % Run from last time value to arbitrary large time.
    % The event function should stop it well before this time.
    tRange = [max(t0),1000000];
    % Run ode45 solver on this part of the cycle. The anonymous
    % function @(t,V)odefcn(t,V,Idlc,Tair) allows ode45 to interpret
    % the ODE function as a function of just t and V, while still
    % using the Idlc and Tair parameters.
    [t1,V1] = ode45(@(t,V)odefcn(t,V,Idlc,Tair), tRange, IC, options);
    % Store array of current values
    I1 = zeros(length(t1),1) + Idlc;
    
    % Continue with first resting step
    Idlc = 0;
    % Initial conditions are from the end of charging
    IC = V1(length(t1),:);
    % Go from end of charging to 15 seconds past it
    tRange = [max(t1), max(t1)+15];
    [t2,V2] = ode45(@(t,V)odefcn(t,V,Idlc,Tair), tRange, IC);
    % Store array of current values
    I2 = zeros(length(t2),1);
    
    % Discharge step
    Idlc = 100;
    % Refresh event function with new current
    options = odeset('Events',@(t, V)event(t,V,Idlc));
    % Initial conditions taken from end of rest
    IC = V2(length(t2),:);
    % Run from end of rest to arbitrary large time.
    % The event function should stop it well before this time.
    tRange = [max(t2),1000000];
    [t3,V3] = ode45(@(t,V)odefcn(t,V,Idlc,Tair), tRange, IC, options);
    % Store array of current values
    I3 = zeros(length(t3),1) + Idlc;
    
    % Second resting step
    Idlc = 0;
    % Initial conditions are from the end of discharging
    IC = V3(length(t3),:);
    % Go from end of discharge to 15 seconds past it
    tRange = [max(t3), max(t3)+15];
    [t4, V4] = ode45(@(t, V)odefcn(t, V, Idlc, Tair), tRange, IC);
    % Store array of current values
    I4 = zeros(length(t4),1);
    
    % Concatenate the results to the t0, V0, and I0 arrays
    t0 = vertcat(t0, t1, t2, t3, t4);
    V0 = vertcat(V0, V1, V2, V3, V4);
    I0 = vertcat(I0, I1, I2, I3, I4);
end

% Go through each time value and calculate the Vdlc at that time
Vdlc0 = zeros(length(t1), 1);
for i = 1:length(t0)
    Vdlc0(i) = voltage(V0(i,:), I0(i));
end

% Repeat process with air temperature of 40°C
Tair = 40;
IC = [0, 0, Tair];

Idlc = 0;
[t40, V40] = ode45(@(t, V)odefcn(t, V, Idlc, Tair), [-50, 0], IC);
I40 = zeros(length(t40),1);

for i = 1:5
    Idlc = -100;
    options = odeset('Events',@(t, V)event(t,V,Idlc));
    IC = V40(length(t40),:);
    tRange = [max(t40),1000000];
    [t1,V1] = ode45(@(t,V)odefcn(t,V,Idlc,Tair), tRange, IC, options);
    I1 = zeros(length(t1),1) + Idlc;
    
    Idlc = 0;
    IC = V1(length(t1),:);
    tRange = [max(t1), max(t1)+15];
    [t2, V2] = ode45(@(t, V)odefcn(t, V, Idlc, Tair), tRange, IC);
    I2 = zeros(length(t2),1);
    
    Idlc = 100;
    options = odeset('Events',@(t, V)event(t,V,Idlc));
    IC = V2(length(t2),:);
    tRange = [max(t2),1000000];
    [t3,V3] = ode45(@(t,V)odefcn(t,V,Idlc,Tair), tRange, IC, options);
    I3 = zeros(length(t3),1) + Idlc;
    
    Idlc = 0;
    IC = V3(length(t3),:);
    tRange = [max(t3), max(t3)+15];
    [t4, V4] = ode45(@(t, V)odefcn(t, V, Idlc, Tair), tRange, IC);
    I4 = zeros(length(t4),1);
    
    t40 = vertcat(t40, t1, t2, t3, t4);
    V40 = vertcat(V40, V1, V2, V3, V4);
    I40 = vertcat(I40, I1, I2, I3, I4);
end

Vdlc40 = zeros(length(t1), 1);
for i = 1:length(t40)
    Vdlc40(i) = voltage(V40(i,:), I40(i));
end

% Plot and save results
% Vdlc at 0°C
figure
yyaxis left
plot(t0, I0)
ylim([-160,160])
xlabel("Time (s)")
ylabel("Current (A)")
grid on
yyaxis right
plot(t0,Vdlc0)
ylim([-0.1,3.1])
ylabel("Voltage (V)")
saveas(gcf,'Vdlc_0.png')
% Voltage at 40°C
figure
yyaxis left
plot(t40, I40)
ylim([-160,160])
xlabel("Time (s)")
ylabel("Current (A)")
grid on
yyaxis right
plot(t40,Vdlc40)
ylim([-0.1,3.1])
ylabel("Voltage (V)")
saveas(gcf,'Vdlc_40.png')
% Temperature at 0°C
figure
yyaxis left
plot(t0, I0)
ylim([-160, 160])
xlabel("Time (s)")
ylabel("Current (A)")
grid on
yyaxis right
plot(t0,V0(:,3))
ylim([-0.4, 12.4])
ylabel("Temperature (°C)")
saveas(gcf,'Tdlc_0.png')
% Temperature at 40°C
figure
yyaxis left
plot(t40,I40)
ylim([-160, 160])
xlabel("Time (s)")
ylabel("Current (A)")
grid on
yyaxis right
plot(t40,V40(:,3))
ylim([39.6, 52.4])
ylabel("Temperature (°C)")
saveas(gcf,'Tdlc_40.png')

% ODE function to determine how Vcf, Vcd, and Tdlc change.
% It takes in a time value, an array containing the Vcf, Vcd,
% and Tdlc, Idlc, and Tair. It returns a vertical array showing
% the derivatives of Vcf, Vcd, and Tdlc.
function dVdt = odefcn(t, V, Idlc, Tair)
Vcf = V(1);
Vcd = V(2);
Tdlc = V(3);
% Obtain electrical parameters at the current temperature
[Rd, Rf, Cd, Cf] = Eparams(Tdlc);
% Solve for Id and If
Id = (Idlc * Rf - Vcf + Vcd) / (Rd + Rf);
If = Idlc - Id;
% Solve for derivatives of Vcf and Vcd
dVcf = -If / Cf;
dVcd = -Id / Cd;
% Obtain temperature parameters
[R, eta, Rth, Cth] = Tparams;
% Determine heat generated
if Idlc >= 0 % if discharging
    Qgen = Idlc^2 * R;
else % if charging
    Qgen = Idlc^2 * R - Idlc * voltage(V, Idlc) * (1 - eta);
end
% Solve for derivative of temperature
dTdlc = (Tair / Rth + Qgen - Tdlc / Rth) / Cth;
% Place derivatives in vertical array to return
dVdt = [dVcf; dVcd; dTdlc];
end

% Event function to stop charge and discharge cycles
function [ev,s,dir] = event(t, V, Idlc)
% Determine Vdlc
Vdlc = voltage(V, Idlc);
% Trigger when Vdlc is rated voltage or half rated voltage
ev = [Vdlc - 2.7; Vdlc - 1.35];
% Stop in both cases
s = [1; 1];
% Only trigger when rising to 2.7V or falling to 1.35V
dir = [1; -1];
end

% Function to calculate Vdlc with given Vcf, Vcd, and Idlc
function Vdlc = voltage(V, Idlc)
Vcf = V(1);
Vcd = V(2);
Tdlc = V(3);
% Obtain electrical parameters at current temperature
[Rd, Rf, Cd, Cf] = Eparams(Tdlc);
% Split Idlc into Id and If
Id = (Idlc * Rf - Vcf + Vcd) / (Rd + Rf);
% Calculate Vdlc
Vdlc = -Id * Rd + Vcd;
end

% Function to determine electrical parameters at a given temperature
function [Rd, Rf, Cd, Cf] = Eparams(Tdlc)
% Take parameter values at each temperature and interpolate
Rd = interp1([-18, 25, 50], [0.238, 0.287, 0.226], Tdlc);
Rf = interp1([-18, 25, 50], [0.549, 0.533, 0.496], Tdlc)*10^-3;
Cd = interp1([-18, 25, 50], [64.3, 64.3, 112.1], Tdlc);
Cf = interp1([-18, 25, 50], [1355, 1368, 1481], Tdlc);
end

% Function to store temperature parameters (to keep in one location)
function [R, eta, Rth, Cth] = Tparams
R = 0.00047;
eta = 0.9;
Rth = 4.5;
Cth = 320;
end