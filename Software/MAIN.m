%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Optimization Strategies for Pump-hydro Storage and Wind Farm Coordination Including Wind Power Uncertainty
%  Jorge Filipe, July 2014
%  ee07300@fe.up.pt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRINT message
fprintf('\n');
fprintf('                   PHS Optimization 2014       \n');
fprintf('              Jorge Filipe, ee07300@fe.up.pt   \n');
fprintf('\n');
fprintf('Notas:\n');
fprintf('- P_waste(max)=0\n');
fprintf('- NAIVE model\n');
fprintf('\n ********************************************************************** \n\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tStart_total=tic;
global fval_vector          %initializtion of the global vector used in the output function

%% INITIALIZATION
sd = 1;         % starting day
nd = 92;%23;        % n. days [1;92]
np = 24;        % n. periods [1;24]
t = 1;          % n. hours 1
s = 50;         % n. scenarios 
prob=1/s;       % probability of each scenario
a=3;            % decison maker risk factor prone:a<0; averse:a>0
k_profit=0.66029;  % profit utility function parameter; k_waste=1-k_profit;
%0.66029; %Decision Maker A
%0.71415  %Decision Maker B
%0.85361  %Decision Maker C

%Output configuration
arredondar =1;          % rounds output to zero if value is below 'to_round'(day ahead)
arredondar_opt =1;      % rounds output to zero if value is below 'to_round'(operational strategy)
to_round=0.1;           % sets precision 

%% HYDRO WITH PUMPED STORAGE STATION CHARACTERIZATION 

%High and lower hydro generation capacity (MW)
hydro_max = 190;                %high level -> Max capacity 190MW
hydro_min = 0.25*hydro_max;     %lower level -> 25% of PhM

%Max pumping capacity (MW)
pump_max = 190;

%Efficiencies
pump_eff = 0.84;      %pumping (%)
hydro_eff = 0.92;     %generation (%)

%Storage
e_max = 630;            %Available capacity (MW)
e_begin = 0.0*e_max;    %initial & final level (MW) 

%Operation costs
hydro_cost = 0;      %pumping costs 
pump_cost = 0;       %hydro generation costs


%% WIND DATA

wind_max=246;        %Wind Farm capacity

%For the Day Ahead Strategy 
% wind_data=xlsread('wind_data.xlsx',1,'C2:ALN2209'); %Normalized Wind power scenarios
% wind_date=xlsread('wind_data.xlsx',1,'A2:B2209');   %Date
% save('wind_data.mat','wind_data');
% save('wind_date.mat','wind_date');
load wind_data
load wind_date

pw=wind_data*wind_max/e_max; %Wind power scenarios normalized by the storage capacity

% For the Operational Management Strategy
% real_wind=xlsread('wind_data.xlsx',2,'A2:A2209');       %wind power measured
% forecast_wind=xlsread('wind_data.xlsx',2,'C2:H2209');   %updated forecasts
% point_wind=xlsread('wind_data.xlsx',2,'B2:B2209');      %point forecasts
% save('real_wind.mat','real_wind');
% save('forecast_wind.mat','forecast_wind');
% save('point_wind.mat','point_wind');

load real_wind
load point_wind
load forecast_wind

real_wind=real_wind*wind_max;                           %no need to normalize
point_wind=point_wind*wind_max;                         %no need to normalize
forecast_wind=forecast_wind*wind_max;                   %no need to normalize

%% PRICE DATA

% p=xlsread('price_data.xlsx',1,'C2:C2409');          %market price
% p_plus=xlsread('price_data.xlsx',1,'E2:E2409');     %positive imbalance price (real > schedule)
% p_minus=xlsread('price_data.xlsx',1,'D2:D2409');    %negative imbalance price (real < schedule)
% save('p.mat','p');
% save('p_plus.mat','p_plus');
% save('p_minus.mat','p_minus');

load p
load p_plus
load p_minus

%% UTILITY FUNCTION DATA

C_max=wind_max/e_max*24;                      %MAX wind power wasted - normalized by the storage capacity
C_min=0*24;                                   %MIN wind power wasted
L_max=(hydro_max/e_max+wind_max/e_max)*24;    %MAX daily profit
L_min=-(hydro_max/e_max+wind_max/e_max)*24;   %MIN daily profit

%% Algorithms

% ###################################################
%       Point Forecast
% ###################################################

[data_EV_PF,data_opt_EV_PF,profit_EV_PF]=point_f_EV_function( nd, sd, np, t, prob, a, k_profit, hydro_max, pump_max, ...
    pump_eff, hydro_eff, e_max, hydro_cost, pump_cost,pw, p, p_plus, p_minus, C_max, C_min, L_max, L_min,wind_max, arredondar, ...
    to_round,wind_date,real_wind,forecast_wind,point_wind);

% ###################################################
%       Expected Value
% ###################################################

[data_EV_U,data_opt_EV_U,profit_EV_U]=scenarios_EV_function( nd, sd, np, t, s, prob, a, k_profit, hydro_max, pump_max, ...
    pump_eff, hydro_eff, e_max, hydro_cost, pump_cost,pw, p, p_plus, p_minus, C_max, C_min, L_max, L_min,wind_max, arredondar, ...
    to_round,wind_date,real_wind,forecast_wind,point_wind);

% ###################################################
%       Unidimensinal Utility Function
% ###################################################

%Averse to risk
% beta=a;
% [data_simple_U_Averse,data_opt_simple_U_Averse,profit_simple_U_Averse]=scenarios_simple_U_function( nd, sd, np, t, s, prob,...
%     beta, k_profit, hydro_max, pump_max, pump_eff, hydro_eff, e_max, hydro_cost, pump_cost,pw, p, p_plus, p_minus, C_max,...
%     C_min, L_max, L_min,wind_max, arredondar, to_round,wind_date,real_wind,forecast_wind,point_wind);

%Prone to risk
% beta=-a;
% [data_simple_U_Prone,data_opt_simple_U_Prone,profit_simple_U_Prone]=scenarios_simple_U_function( nd, sd, np, t, s, prob, beta,...
%     k_profit, hydro_max, pump_max, pump_eff, hydro_eff, e_max, hydro_cost, pump_cost,pw, p, p_plus, p_minus, C_max,...
%     C_min, L_max, L_min,wind_max, arredondar, to_round,wind_date,real_wind,forecast_wind,point_wind);

% ###################################################
%       Multi-attribute Utility Function
% ###################################################

%Averse to curtailement of wind power
% beta=a;
% [data_multi_U_Averse,data_opt_multi_U_Averse,profit_multi_U_Averse]=scenarios_multi_U_function( nd, sd, np, t, s, prob, beta, ...
%     k_profit, hydro_max, pump_max, pump_eff, hydro_eff, e_max, hydro_cost, pump_cost,pw, p, p_plus, p_minus, C_max, ...
%     C_min, L_max, L_min,wind_max, arredondar, to_round,wind_date,real_wind,forecast_wind,point_wind);

%Prone to curtailement of wind power
% beta=-a;
% [data_multi_U_Prone,data_opt_multi_U_Prone,profit_multi_U_Prone]=scenarios_multi_U_function( nd, sd, np, t, s, prob, beta, ...
%     k_profit, hydro_max, pump_max, pump_eff, hydro_eff, e_max, hydro_cost, pump_cost,pw, p, p_plus, p_minus, C_max, ...
%     C_min, L_max, L_min,wind_max, arredondar, to_round,wind_date,real_wind,forecast_wind,point_wind);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = msgbox('Operation Completed: MAIN 1');

tEnd = toc(tStart_total);
fprintf('Total: %d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));
diary off