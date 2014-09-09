%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Optimization Strategies for Pump-hydro Storage and Wind Farm Coordination Including Wind Power Uncertainty
%  Jorge Filipe, July 2014
%  ee07300@fe.up.pt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
close all
clear 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRINT message
fprintf('\n');
fprintf('                   PHS Optimization 2014       \n');
fprintf('              Jorge Filipe, ee07300@fe.up.pt   \n');
fprintf('\n');
fprintf('Notas:\n');
fprintf('- P_waste(max)=0\n');
fprintf('- d_23,d_24(max)=0\n');
fprintf('\n ********************************************************************** \n\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tStart_total=tic;
% global fval_vector          %initialization of the global vector used in the output function

%% INITIALIZATION

main_data.sd = 1;               % starting day
main_data.nd = 1;%23;          % n. days [1;92]
main_data.np = 24;              % n. periods [1;24]
main_data.t = 1;                % n. hours 1
main_data.s = 10;               % n. scenarios 
main_data.prob = 1/main_data.s; % probability of each scenario
main_data.a = 3;                % decison maker risk factor prone:a<0; averse:a>0
main_data.DM_A=0.66029;  %Decision Maker A
main_data.DM_B=0.71415;  %Decision Maker B
main_data.DM_C=0.85361;  %Decision Maker C
main_data.k_profit=main_data.DM_A;     % profit utility function parameter; k_waste=1-k_profit;

%Output configuration
main_data.arredondar =1;          % rounds output to zero if value is below 'to_round'(day ahead)
main_data.arredondar_opt =1;      % rounds output to zero if value is below 'to_round'(operational strategy)
main_data.to_round=0.1;           % sets precision 

%% HYDRO WITH PUMPED STORAGE STATION CHARACTERIZATION 

%Max hydro generation capacity (MW)
system_data.hydro_max = 190;                %high level -> Max capacity 190MW

%Max pumping capacity (MW)
system_data.pump_max = 190;

%Efficiencies
system_data.pump_eff = 0.84;      %pumping (%)
system_data.hydro_eff = 0.92;     %generation (%)

%Storage
system_data.e_max = 630;                        %Available capacity (MW)
system_data.e_begin = 0.0*system_data.e_max;    %initial & final level (MW) 

%Operation costs
system_data.hydro_cost = 0;      %pumping costs 
system_data.pump_cost = 0;       %hydro generation costs


%% WIND DATA

wind_data.wind_max=246;        %Wind Farm capacity

%For the Day Ahead Strategy 
% wind_dados=xlsread('wind_data.xlsx',1,'A2:ALN2209'); 								
% wind_data.scenarios=wind_dados(:,3:end)*wind_data.wind_max/system_data.e_max; 		%Wind power scenarios normalized by the storage capacity
% wind_data.date=wind_dados(:,1:2);													%Wind power scenarios date and hour
% clear wind_dados

% For the Operational Management Strategy
% wind_dados=xlsread('wind_data.xlsx',2,'C2:J2209'); 
% wind_data.real_wind=wind_dados(:,1)*wind_data.wind_max;								%wind power measured, without normalization
% wind_data.point_forecast=wind_dados(:,2)*wind_data.wind_max;						%point forecasts, without normalization
% wind_data.updated=wind_dados(:,3:end)*wind_data.wind_max;							%updated forecasts, without normalization
% clear wind_dados

% save('wind_data.mat','wind_data');
load wind_data

%% PRICE DATA

% market_dados=xlsread('price_data.xlsx','Matlab_input','A3:H6770');         
% [r,c] = size(market_dados);
% market_dados_new=zeros(r/2,c);
% k=1;
% for j=1:r/2
%     market_dados_new(j,:)=market_dados(k,:);
%     k=k+2;
% end
% clear market_dados
% market_dados=market_dados_new;
%
% market_data.price_forecast=market_dados(:,1);
% market_data.p_plus_forecast=market_dados(:,2);
% market_data.p_minus_forecast=market_dados(:,3);
% market_data.price_real=market_dados(:,4);
% market_data.p_plus_real=market_dados(:,5);
% market_data.p_minus_real=market_dados(:,6);
% market_data.p_plus_prob=market_dados(:,7);
% market_data.p_minus_prob=market_dados(:,8);
% clear market_dados
% 
% save('market_data.mat','market_data');
load market_data

market_data.positive_imbalance=market_data.p_plus_forecast.*market_data.p_plus_prob;   %generation bigger than schedule
market_data.negative_imbalance=market_data.p_minus_forecast.*market_data.p_minus_prob; %generation smaller than schedule
market_data.price_normalization=max(market_data.price_forecast);

market_data.price_forecast=market_data.price_forecast/market_data.price_normalization;
market_data.positive_imbalance=market_data.positive_imbalance/market_data.price_normalization;
market_data.negative_imbalance=market_data.negative_imbalance/market_data.price_normalization;

%% UTILITY FUNCTION DATA

main_data.C_max=wind_data.wind_max/system_data.e_max*24;              						          %MAX wind power wasted - normalized by the storage capacity
main_data.C_min=0*24;                                   											  %MIN wind power wasted
main_data.L_max=(system_data.hydro_max/system_data.e_max+wind_data.wind_max/system_data.e_max)*24;    %MAX daily profit
main_data.L_min=-(system_data.hydro_max/system_data.e_max+wind_data.wind_max/system_data.e_max)*24;   %MIN daily profit

%% Algorithms

% ###################################################
%       Point Forecast Optimization
% ###################################################

[PointForecast.Day_Ahead,PointForecast.Operational_Strategy,PointForecast.profit]=PointForecast_optimization(system_data, wind_data, market_data, main_data);

% ###################################################
%       Expected Value Optimization
% ###################################################

[ExpectedValue.Day_Ahead,ExpectedValue.Operational_Strategy,ExpectedValue.profit]=ExpectedValue_optimization(system_data, wind_data, market_data, main_data);

gamma=(sum(ExpectedValue.profit(:,2))-sum(PointForecast.profit(:,2)))/sum(PointForecast.profit(:,2))*100

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
clear h
tEnd = toc(tStart_total);
fprintf('Total: %d minutes and %f seconds\n',floor(tEnd/60),rem(tEnd,60));
clear tEnd tStart_total
diary off