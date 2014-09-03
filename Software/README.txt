%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Optimization Strategies for Pump-hydro Storage and Wind Farm Coordination Including Wind Power Uncertainty
%  Jorge Filipe, July 2014
%  ee07300@fe.up.pt
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# MAIN.mat - Main file, where all setup options are defined

scenarios_multi_U_function.m - Multi-attribute Utility Function
scenarios_EV_function.m - Expected Value
scenarios_simple_U_function.m - Unidimensional Utility Function
point_f_EV_function - Point Forecast Optimization

# nl_obj_fun.m - gradient and objective function of the Multi-attribute Utility Function
# nl_obj_fun_EV.m - gradient and objective function of the Expected Value
# nl_obj_fun_simple_utility.m - gradient and objective function of the Unidimensional Utility Function
# nl_obj_fun_PF.m - gradient and objective function of the Point Forecast

# hessianfcn.m - Hessian Matrix of the Multi-attribute Utility Function
# hessianfcn_EV.m - Hessian Matrix of the Expected Value and Point Forecast
# hessianfcn_simply_utility.m - Hessian Matrix of the Unidimensional Utility Function

# outfun.m - Output fun used in the fmincon algorithm aiming to stop the optimization if fval is increasing

# operational_strat.m - Operational Management Strategy Optimization

#price_data.xlsx - market and imbalance prices 

# wind_data.xlsx - wind power scenarios, wind power measured, wind power updated forecasts and point forecasts


