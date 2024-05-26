% MAE 119 Introduction to Renewable Energy
% Spring 2024
% Prof. Hidalgo-Gonzalez
% Homework 3

%%
clear all
clc
% interest rate for Net Present Value calculations
rate = 4/100; % 4.0%
% data for generator dispatch problem
D = 365; % days in a year
Y = 20; % years simulated
T = 24; % hours simulated in each day
t= (1:T); % vector of time
 d = 40*[0.4 0.4 0.4 0.65 ...
     0.65 0.6  0.6  0.65  ...
     0.8  0.85  0.8  0.8  ...
     0.75  0.6  0.5  0.3  ...
     0.3  0.4  0.54  0.6  ...
     0.75  0.7  0.6  0.5 ]; % daily demand in MWh, demand per hour for 24 hrs

n = 4; % number of dispatchable generators (2 hard coal and 2 OCGT plants)
% to choose from, in addition to  a solar plant

solar_cap_fact = [zeros(1,6) sin(pi*[0:12]/12) zeros(1,length(d)-6-length(0:12))]; 
% capacity factor of the potential solar plant = power_available/installed_capacity

% (To be used in Part II only: capacity factor of the potential wind plant = 
% power_available/installed_capacity)
wind_cap_fact = [0.3277    0.2865    0.3303    0.3073    0.2994...
    0.3071    0.2980    0.2988    0.3149    0.3141    0.3142 ...
    0.5*0.3067   0.5*0.2879   0.5*0.3072   0.5*0.3163   0.5*0.3049   0.5*0.3103 ...
    0.3914    0.4008    0.3879    0.3889    0.3999 0.4153    0.3923]; 


Pmax = [10 5 10 15];  % Maximum generator capacities in MW 
% [hard coal, hard coal, OCGT, OCGT] due to land constraints
Pmax_solar = 40; % Maximum capacity in MW for solar plant due to land constraints
% in part 2, no solar capacity constraints
Pmin = [0 0 0 0]; % generator minimum capacities in MW
Pmin_solar = 0;
%pt. 2
Pmin_wind = 0; 

ramp_hard_coal = 0.015 % 1.5% [% Pnom/minute]
hourly_ramp_hard_coal = ramp_hard_coal*60; % [% Pnom/h]

ramp_ocgt = 0.08 % 1.5% [% Pnom/min]
hourly_ramp_ocgt = ramp_ocgt*60; % [% Pnom/h]

R = [hourly_ramp_hard_coal*Pmax(1) hourly_ramp_hard_coal*Pmax(2) ...
    hourly_ramp_ocgt*Pmax(3) hourly_ramp_ocgt*Pmax(4)];  % ramp-rate limits [MWh] for non PV

% source:
% https://www.eia.gov/electricity/annual/html/epa_08_04.html
% 
fuel_cost_per_MWh_2018 = [25.4 25.4 27.35 27.35]; % $/MWh, 2018 for non PV or wind
%
%fuel_cost_per_MWh_2019 = [24.28 24.28 23.11 23.11]; % $/MWh, 2019 This
%cost won't be used

% Investment costs (or capital costs per MW of built capacity). 
% From https://atb.nrel.gov/
%wind relevant in pt.2
capital_cost_solar_plant_per_MW = 1200000; % dollars/MW, denoted as M_i in the pdf
capital_cost_wind_plant_per_MW = 1300000;  % dollars/MW, denoted as M_i in the pdf
capital_cost_coal_plant = 4200000; % dollars/MW, denoted as M_i in the pdf
capital_cost_gas_plant = 2600000; % dollars/MW, denoted as M_i in the pdf

%pt.2 relevant
% Emission rates of pounds of CO_2 per MWh for each tech 
% (https://www.eia.gov/tools/faqs/faq.php?id=74&t=11)
co2_coal = 2.30 * 10^3; % pounds per MWh
co2_gas = 0.97 * 10^3; % pounds per MWh
co2_solar = 0;

%co2_emissions = [co2_coal co2_coal co2_gas co2_gas]

%% Solution - Part 1

%% Optimal capacity and hourly dispatch

cvx_begin

variables p(n,T) p_solar(1,T) installed_thermal_capacity(1,n) ...
    installed_solar_capacity(1);

% To-do: add constraints for upper and lower bounds for installed capacity
installed_thermal_capacity <= Pmax;
installed_solar_capacity <= Pmax_solar;
installed_thermal_capacity >= Pmin;
installed_solar_capacity >= Pmin_solar;

% To-do: add constraints for upper and lower bounds for dispatched electricity

sum(p,1) + p_solar >= d;
p >= 0;
for i = 1:n
    p(i,:) <= installed_thermal_capacity(i);
end

% To-do: add constraint that will make the output of the solar plant
% (p_solar) equal to its capacity factor times the installed capacity
% (similar to lecture, but now the installed capacity is the variable
% installed_solar_capacity)
p_solar == solar_cap_fact * installed_solar_capacity;

% Implement ramp rate control for thermal plants
abs(p(:,2:T)-p(:,1:T-1)) <= R'*ones(1,T-1);

% To-do: Calculate annual fuel costs, then calculate its net present value
years = 1:Y;
annual_fuel_cost = D*sum(fuel_cost_per_MWh_2018 * p);
NPV_total_fuel_cost = sum(annual_fuel_cost./((1+rate).^years));


% To-do: Calculate the investment cost
Investment_cost = capital_cost_solar_plant_per_MW * installed_solar_capacity ...
    + capital_cost_coal_plant * installed_thermal_capacity(1) ...
    + capital_cost_coal_plant * installed_thermal_capacity(2) ...
    + capital_cost_gas_plant * installed_thermal_capacity(3) ...
    + capital_cost_gas_plant * installed_thermal_capacity(4);

% State minimization objective
minimize (NPV_total_fuel_cost + Investment_cost)

cvx_end


% white background for plots:
set(gcf,'color','w');
   subplot(3,1,1)
   plot(t,d, 'LineWidth',2.0);
   title('Demand')
   xlabel('Hours')
   ylabel('MWh')
   subplot(3,1,2)
   plot(t,p, t, p_solar, 'LineWidth',2.0);
   title('Generation')
   xlabel('Hours')
   ylabel('MWh')
   legend('coal 1', 'coal 2', 'gas 1', 'gas 2', 'solar')
   subplot(3,1,3)
   plot(t, sum(p)+p_solar, t, d, 'LineWidth',2.0)
   title('Gen and demand')
   xlabel('Hours')
   ylabel('MWh')
   legend('generation', 'demand')
print -depsc gen_dispatch

installed_thermal_capacity;
installed_solar_capacity;

%% LCOE for each generator from using the optimal dispatch obtained above

% Set up useful arrays for LCOE and fuel costs; input values for costs
LCOE = zeros(1,n+1);
% Fuel_Costs = zeros(n,Y);
Fuel_Costs = zeros(1,n);
for i = 1:n
    Fuel_Costs(i) = D * fuel_cost_per_MWh_2018(i) * sum(p(i,:));
end

% Coal LCOE
for i = 1:2
    LCOE(i) = (capital_cost_coal_plant * installed_thermal_capacity(i) ...
        + sum(Fuel_Costs(i) .* (1+rate).^-years)) ...
        / sum(D*sum(p(i,:)) * (1+rate).^-years);
end

% Gas LCOE
for i = 3:n
    LCOE(i) = (capital_cost_gas_plant * installed_thermal_capacity(i) ...
        + sum(Fuel_Costs(i) .* (1+rate).^-years)) ...
        / sum(D*sum(p(i,:)) * (1+rate).^-years);
end

% Solar LCOE
LCOE(n+1) = capital_cost_solar_plant_per_MW * installed_solar_capacity ...
    / sum(D*sum(p_solar) * (1+rate).^-years);

LCOE;

%% Part II:
%{
%}
co2 = [co2_coal co2_coal co2_gas co2_gas];
emission = sum(co2*p)*D; 
emissionBaseline = sum(emission); %baseline parameter

% size curtailment vector
curtailment = zeros(6,24);

for v = 1:6  %6 scenarios starting at 100 co2
    perc = [1 .8 .6 .4 .2 0]; %percentages 100 to 0
    
cvx_begin

variables p(n,T) p_solar(1,T) installed_thermal_capacity(1,n) ...
    installed_solar_capacity(1) p_wind(1,T) installed_wind_capacity(1);

% To-do: add constraints for upper and lower bounds for installed capacity
%and CO2 emissions

installed_thermal_capacity <= Pmax;
%installed_solar_capacity <= Pmax_solar; No max capacity for solar or wind in Pt.2
installed_wind_capacity >= Pmin_wind;
installed_thermal_capacity >= Pmin;
installed_solar_capacity >= Pmin_solar;

% To-do: add constraints for upper and lower bounds for dispatched electricity
% adjust for carbon

sum(p,1) + p_solar + p_wind >= d; %greater than or eq to demand (a double)
p >= 0;  %p is vector of (4,24) each plant with 24 hrs of power gen
emission = zeros(4,24);
for i = 1:n  %for all thermal plants, p is less than/= installed capacity
    p(i,:) <= installed_thermal_capacity(i);
end 
co2 = [co2_coal co2_coal co2_gas co2_gas]; %co2_gas and coal are just doubles
emission = sum(co2*p)*D*20;  
emissionTotal = sum(emission);
emissionTotal >= 0
emissionTotal <= perc(v)*emissionBaseline

% To-do: add constraint that will make the output of the solar plant
% (p_solar) equal to its capacity factor times the installed capacity
% (similar to lecture, but now the installed capacity is the variable
% installed_solar_capacity)

p_solar == solar_cap_fact * installed_solar_capacity; %in pt.2, only minimums
p_wind == wind_cap_fact * installed_wind_capacity; 

% Implement ramp rate control for thermal plants

abs(p(:,2:T)-p(:,1:T-1)) <= R'*ones(1,T-1);

% To-do: Calculate the investment cost (pt.2)
%This varies per scenario
Investment_cost2 = capital_cost_solar_plant_per_MW * installed_solar_capacity ...
    + capital_cost_coal_plant * installed_thermal_capacity(1) ...
    + capital_cost_coal_plant * installed_thermal_capacity(2) ...
    + capital_cost_gas_plant * installed_thermal_capacity(3) ...
    + capital_cost_gas_plant * installed_thermal_capacity(4) ...
    + capital_cost_wind_plant_per_MW * installed_wind_capacity;

% To-do: Calculate annual fuel costs, then calculate its net present value

years = 1:Y;
annual_fuel_cost2 = D*sum(fuel_cost_per_MWh_2018 * p); %constraint is p
NPV_total_fuel_cost2 = sum(annual_fuel_cost2./((1+rate).^years)); 
totalNPV = NPV_total_fuel_cost2 + Investment_cost2;

% co2 emissions total //this code is not actually used
co2_emissions = co2_coal*installed_thermal_capacity(1,1) +...
    co2_coal*installed_thermal_capacity(1,2) + ...
    co2_gas*installed_thermal_capacity(1,3) + ...
    co2_gas*installed_thermal_capacity(1,4) ;

% %Curtailment
% % curtailment(v) <= ( abs( d.*D - sum( p(1,:) + p(2,:) + p(3,:) + p(4,:) + p_wind + p_solar ) ) ); %d=daily demand D=365 days
% abs(d.*D - sum(p(1,:) + p(2,:) + p(3,:) + p(4,:) + p_wind + p_solar)) >= curtailment(v) ; %d=daily demand D=365 days

% State minimization objective
 minimize (totalNPV) 

cvx_end
%should return 100% 80% 60% etc of the baseline (6 values in each vector)
scenarioTotals(v) = emissionTotal; %vector of 6
npvTotals(v) = totalNPV;
fuelTotals(v) = NPV_total_fuel_cost2;
investTotals(v) = Investment_cost2;
installedCoal(v) = installed_thermal_capacity(1) + installed_thermal_capacity(2);
installedGas(v) = installed_thermal_capacity(3) + installed_thermal_capacity(4);
installedSolar(v) = installed_solar_capacity;
installedWind(v) = installed_wind_capacity;

curtailment(v,:)=(abs( d - sum( p(1,:) + p(2,:) + p(3,:) + p(4,:) + p_wind + p_solar ) ) );

end 
%format long
npvTotals;
scenarioTotals; %carbon dioxide emissions per scenario
NPV_total_fuel_cost + Investment_cost; %from Part 1
   for i = 1:6
       scenesTh{i} = [sum(curtailment(i,:))]; %add up each power plant's curtailment
   end
scenesTh;

%% Baseline, no CO2 cap
%% CO2 cap


%% Summary figs
% white background for plots:
figure;
set(gcf,'color','w');
   tiledlayout(1,3)
   nexttile
   scenesO = cell(1,6);
   for i = 1:6
       scenesO{i} = [fuelTotals(i) investTotals(i) npvTotals(i)];
   end
   bar([scenesO{1}; scenesO{2}; scenesO{3}; scenesO{4}; scenesO{5}; scenesO{6}])
   name={'Unconstrained';'<80%';'<60%';'<40%';'<20%';'0%'};
   set(gca,'xticklabel',name);
   grid on;
   title('Total cost in NPV')
   xlabel('CO2 Scenario Compared to Unconstrained (Baseline)')
   ylabel('Cost in NPV')
   legend('Fuel Cost', 'Investment Cost', 'Total Cost');
   
   nexttile %per energy type (coal gas solar wind)
   scenes = cell(1,6);
   for i = 1:6
       scenes{i} = [installedCoal(i) installedGas(i) installedSolar(i) installedWind(i)];
   end
   bar([scenes{1}; scenes{2}; scenes{3}; scenes{4}; scenes{5}; scenes{6}])
   names={'Unconstrained';'<80%';'<60%';'<40%';'<20%';'0%'};
   set(gca,'xticklabel',names);
   title('Installed Capacity by Energy Type')
   xlabel('CO2 Scenario Compared to Unconstrained (Baseline)')
   ylabel('Installed Capacity [MW]')
   legend('Coal', 'Gas', 'Solar', 'Wind')
   
   nexttile %curtailment = abs( demand - (p + p_wind + p_solar) )
   scenesTh = cell(1,6);
   for i = 1:6
       scenesTh{i} = [sum(curtailment(i,:))]; %add up each plant's curtailment
   end
   bar([scenesTh{1}; scenesTh{2}; scenesTh{3}; scenesTh{4}; scenesTh{5}; scenesTh{6}])
   name={'Unconstrained';'<80%';'<60%';'<40%';'<20%';'0%'};
   set(gca,'xticklabel',name);
   title('Curtailment')
   xlabel('CO2 Scenario Compared to Unconstrained (Baseline)')
   ylabel('Curtailment [MWh]')
print -depsc gen_dispatch

%format long
npvTotals(1);
 scenes{6};
scenarioTotals; %carbon dioxide emissions per scenario
NPV_total_fuel_cost + Investment_cost; %from Part 1
curtailment;

%{
%}