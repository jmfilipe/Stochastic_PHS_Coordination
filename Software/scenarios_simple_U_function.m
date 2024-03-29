function [data, data_opt,profits_opt] = scenarios_simple_U_function( nd, sd, np, t, s, prob, a, k_profit, hydro_max, pump_max, pump_eff, hydro_eff, e_max, hydro_cost, pump_cost,pw, p, p_plus, p_minus, C_max, C_min, L_max, L_min,wind_max, arredondar, to_round, wind_date,real_wind,forecast_wind,point_wind)
data=[];
data_opt=[];
profits_opt=[];
e_storage_initial = 0.0*e_max;     %initial & final level (MW) 
global fval_vector

%modiefiers used for the naive function
p_plus_modiefied=p*0.87;
p_minus_modiefied=p*1.09;
p_modiefied=p;

% p(1:168)=[];
% p_plus(1:168)=[];
% p_minus(1:168)=[];

%Creates log file with date and time as the filename
% DATO=datestr(now);
% filename_log=['matlab-log-',DATO];
% filename_log = strrep(filename_log,':','_');
% filename_log = strrep(filename_log,' ','_');
% filename_log = sprintf('%s%s%s%d%s%d','C:\Users\cpsilva\OneDrive\Tese\Software\Inicio\Results\Logs\',filename_log,'_____',sd,'_',nd);
% eval(['diary ',filename_log]);

%% ALGORITH
for d = sd:nd   %simulates for the time horizon defined - nd: n. of days
    tStart = tic;
    %==================================================================
    % Optimization
    %==================================================================

    % LINEAR PROGRAM PARAMETERS

    n_var = (1+7*s)*np;            % num. variables 

    % equal equations
    neq = (4*s+1)*np;           %n equal equations
    Aeq = sparse(neq,n_var); 
    beq = sparse(1,neq);

    % inequalities
    nde = (4*s+1)*np;             %n inequal equations              
    Ad = sparse(nde,n_var);
    bd = sparse(1,nde);

    % variable limits
    lb = sparse(1,n_var);     %lower limits
    ub = sparse(1,n_var);     %upper limits


    % ALGORITHM

            for j = 1:np  %simulates for every hour of the day

                %==============================================================
                %EQUAL EQ.
                %==============================================================
                
                %Eq. 0: Market Bid for period i
                % 
                BID=0;
                P_hydro=1;
                P_grid=s+1;
                for eq=(0):(s-1)
                    Aeq(eq*np+j,P_hydro*np+j) = -1;     %PH
                    Aeq(eq*np+j,P_grid*np+j) = -1;      %PG
                    Aeq(eq*np+j,BID*np+j) = 1;          %BID
                    P_hydro=P_hydro+1;
                    P_grid=P_grid+1;
                end
                clear P_grid
                clear P_hydro
                clear BID
                
                %Eq. 1 & 2 Reservoir level
                %
                if j == 1 %1st hour
                    
                    %Eq. 1 Initial Reservoir Level
                    E_storage=3*s+1;
                    for eq=(s):(2*s-1)
                        Aeq(eq*np+j,E_storage*np+j) = 1;
                        beq(eq*np+j) = e_storage_initial;                
                        E_storage=E_storage+1;
                    end
                    clear E_storage
                    clear eq
                           
                
                
                else %remaining hours
                    
                    %Eq. 2 Storage level
                    E_storage=3*s+1;
                    P_hydro=1;
                    P_pump=2*s+1;
                    for eq=(2*s):(3*s-1)
                        Aeq(eq*np+j,P_hydro*np+j-1) = t/hydro_eff;
                        Aeq(eq*np+j,P_pump*np+j-1) = -t*pump_eff;                    
                        Aeq(eq*np+j,E_storage*np+j) = 1;
                        Aeq(eq*np+j,E_storage*np+j-1 ) = -1;
                        E_storage=E_storage+1;
                        P_hydro=P_hydro+1;
                        P_pump=P_pump+1;
                    end
                    clear P_pump
                    clear E_storage
                    clear eq
                    clear P_hydro
                end
                             
                %Eq. 3: Imbalance for every scenario 
                P_grid=s+1;
                P_pump=2*s+1;
                d_imbalance=4*s+1;
                P_waste=6*s+1;
                for eq=(3*s):(4*s-1)
                    Aeq(eq*np+j,P_grid*np+j) = 1;
                    Aeq(eq*np+j,P_pump*np+j) = 1;
                    Aeq(eq*np+j,d_imbalance*np+j) = 1;
                    Aeq(eq*np+j,P_waste*np+j) = 1;
                    %(P_grid-s) equals from 1 to s, this way the use of another incremental variable is avoided
                    beq(eq*np+j) = pw((d-1)*np+j,P_grid-s);
                    P_grid=P_grid+1;
                    P_waste=P_waste+1;
                    d_imbalance=d_imbalance+1;
                    P_pump=1+P_pump;
                end
                clear P_pump
                clear P_grid
                clear d_imbalance
                clear P_waste
                clear eq
                                              
                %==============================================================
                %INEQUAL EQ.
                %==============================================================

                %Ineq. 0: Hydro Generation
                P_hydro=1;
                P_pump=2*s+1;
                E_storage=3*s+1;
                for eq=(0):(s-1)
                    Ad(eq*np+j,P_hydro*np+j) = 1;
                    Ad(eq*np+j,E_storage*np+j) = -hydro_eff/t;  
                    Ad(eq*np+j,P_pump*np+j) = -t*hydro_eff*pump_eff;
                    P_pump=1+P_pump;
                    P_hydro=1+P_hydro;
                    E_storage=E_storage+1;
                end
                clear P_pump
                clear E_storage
                clear eq
                clear P_hydro
                              
                %Ineq. 1: Epigraph Form (1)
                d_imbalance=4*s+1;
                T_costs=5*s+1;
                for eq=(s):(2*s-1)
                    Ad(eq*np+j,(T_costs)*np+j) = -1;           %Ts,i
                    Ad(eq*np+j,(d_imbalance)*np+j) = (p_modiefied((d-1)*np+j)-p_plus_modiefied((d-1)*np+j));     %d,si
                    bd((eq)*np+j)=0;
                
                    %Ineq. 2: Epigraph Form (2)
                    Ad((eq+s)*np+j,(T_costs)*np+j) = -1;           %Ts,i
                    Ad((eq+s)*np+j,(d_imbalance)*np+j) = -(p_minus_modiefied((d-1)*np+j)-p_modiefied((d-1)*np+j));     %d,si
                    bd((eq+s)*np+j)=0; 
                    
                    d_imbalance=d_imbalance+1;
                    T_costs=T_costs+1;
                end
                clear T_costs
                clear d_imbalance
                clear eq
                
                % Ineq. 3: Hydro + Pump Constraint
                P_hydro=1;
                P_pump=2*s+1;
                for eq=(3*s):(4*s-1)
                    Ad(eq*np+j,P_hydro*np+j) = 1;
                    Ad(eq*np+j,P_pump*np+j) = 1;
                    bd(eq*np+j) = hydro_max/e_max;
                    P_pump=P_pump+1;
                    P_hydro=P_hydro+1;
                end
                clear P_hydro
                clear P_pump
                clear eq
                              
% %                 Ineq. 6: Objective Function Bounds (1)
%                 
%                 Ad((5*s)*np+j,0*np+j) = p_modiefied((d-1)*np+j);
%                 Ad((5*s)*np+j,(2*s+1)*np+j) = p_modiefied((d-1)*np+j);
%                 Ad((5*s)*np+j,(2*s+2)*np+j) = -1;
%                 bd((5*s)*np+j) = L_max/24;
%                 
% %                 Ineq. 7: Objective Function Bounds (1)
%                 
%                 Ad((5*s+1)*np+j,0*np+j) = -p_modiefied((d-1)*np+j);
%                 Ad((5*s+1)*np+j,(2*s+1)*np+j) = -p_modiefied((d-1)*np+j);
%                 Ad((5*s+1)*np+j,(2*s+2)*np+j) = 1;
%                 bd((5*s+1)*np+j) = -L_min/24;
%                
                %==============================================================
                %BOUNDS
                %==============================================================

                %upper bound
                ub(0*np+j) = Inf;                                   %BID
                for incre=1:s 
                    ub((incre)*np+j) = hydro_max/e_max;             %Hydro
                    ub((2*s+incre)*np+j) = pump_max/e_max;          %Pump
                    ub((s+incre)*np+j) = max(pw((d-1)*np+j,1:s));   %Pgrid
                    ub((3*s+incre)*np+j) = e_max/e_max;             %E_storage
                    ub((4*s+incre)*np+j) = Inf;                     %d_imbalance
                    ub((5*s+incre)*np+j) = Inf;                     %T_costs
                    ub((6*s+incre)*np+j) = Inf;                     %P_waste
                end
                clear incre
                
                %lower bound
                for incre=1:s
                    lb((4*s+incre)*np+j) =-Inf;                     %d_imbalance
                end
                clear incre
            end
            
            
        disp_message=['Utility Function, day: ', num2str(d)];
        disp(disp_message)
        
        % Optimization Funcion
 
        x0= sparse(n_var,1);
        options = optimset('Algorithm','interior-point','Display','iter-detailed','GradObj','on','DerivativeCheck','on','FunValCheck','on','UseParallel','always',...
    'Hessian','on','HessFcn',@(x,lambda)hessianfcn_simple_utility(x,lambda,s,a,k_profit,np,L_min, L_max,n_var,d,p_modiefied,hydro_cost,pump_cost),'OutputFcn',@outfun);%notify iter
            
        options.MaxFunEvals = (s*50000000); 
        options.MaxIter = 250;
        options.TolFun = 1e-9;
        options.TolCon = 1e-7;
        options.TolX = 1e-9;
        
        [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x)nl_obj_fun_simple_utility(x,np,d,s,p_modiefied,hydro_cost,pump_cost,a,k_profit,wind_max,n_var, hydro_max, e_max, C_min, C_max, L_min, L_max),x0,Ad,bd,Aeq,beq,lb,ub,[],options);
        
        % SE DER EXIT FLAG NEGATIVA
        
        if 1;%(exitflag ~= 1) && (exitflag ~= 2 )
          
            x_error=x;
            fval_error=fval;
            clear x
            clear x0
            clear output
            clear lambda
            clear grad
            clear hessian
            clear exitflag
            
            % CALCULA PONTO INICIAL
            
            x0= sparse(n_var,1);
            %f_initial = zeros(size(x0)); % assumes x0 is the initial point        

            options_lin = optimset('Algorithm','interior-point','simplex','on','Display','notify','GradObj','on','DerivativeCheck','off','UseParallel','always',...
                'Hessian','on','HessFcn',@(x,lambda)hessianfcn_EV(x,lambda,n_var));
                    options.MaxFunEvals = (s*1000000); 
                options_lin.MaxIter = 50;
                options_lin.TolFun = 1e-8;
                options_lin.TolCon = 1e-8;
                options_lin.TolX = 1e-8;
            xnew = fmincon(@(x)initial_x0(x0),x0,Ad,bd,Aeq,beq,lb,ub,[],options_lin);
            x0= xnew;

            %SEGUNDA FUNCAO USANDO O PONTO INICIAL
            
            options = optimset('Algorithm','interior-point','Display','iter-detailed','GradObj','on','DerivativeCheck','off','FunValCheck','on','UseParallel','always',...
    'Hessian','on','HessFcn',@(x,lambda)hessianfcn_simple_utility(x,lambda,s,a,k_profit,np,L_min, L_max,n_var,d,p_modiefied,hydro_cost,pump_cost),'OutputFcn',@outfun);%notify iter

            options.MaxFunEvals = (s*1000000); 
            options.MaxIter = 250;
            options.TolFun = 1e-9;
            options.TolCon = 1e-7;
            options.TolX = 1e-9;
            
            [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x)nl_obj_fun_simple_utility(x,np,d,s,p_modiefied,hydro_cost,pump_cost,a,k_profit,wind_max,n_var, hydro_max, e_max, C_min, C_max, L_min, L_max),x0,Ad,bd,Aeq,beq,lb,ub,[],options);
            
            if fval_error < fval
                x=x_error;
                fprintf ('************\n...isto e estranho\n**********\n');
            else
                fprintf ('************\n...tudo na boa\n**********\n');
            end
        
        end
        
        %save x;
        %% output
        x=x;
        xx = x;%round(x*10000)/10000;
        dados=zeros(np,n_var/np);
        for i=1:n_var/np
            dados(1:np,i)=xx(np*(i-1)+1:i*np,1);
        end
        clear i
        dados=dados*e_max;
        arredondar=1;
        if arredondar == 1    
            [r,c] = size (dados);
            for j=1:c
                for i=1:r
                    if(abs(dados(i,j))<to_round)
                        dados(i,j)=0;
                    end
                end
            end
        end
        
        deviations=dados(:,3*s+3+1:4*s+2+1);
        power_to_grid=dados(:,2*s+1+1);
        for j=1:24
            vetor_a_testar=deviations(j,:);
            if all(vetor_a_testar < -0.01)
                maximum=max(vetor_a_testar);
                for i=1:s
                    deviations(j,i)=deviations(j,i)+abs(maximum);  
                end
                power_to_grid(j,1)=power_to_grid(j,1)+maximum;
            else
                %do something else
            end
        end
        dados(:,3*s+3+1:4*s+2+1)=deviations;
        dados(:,2*s+1+1)=power_to_grid(:,1);
        
        %% OPERATIONAL STRATEGY
        day_ahead_data=horzcat(dados(:,1),dados(:,2*s+2));
        [dados_opt,profit_opt,taxas_opt]=operational_strat( np,d, day_ahead_data, e_max, e_storage_initial, hydro_eff, pump_eff, t, real_wind, forecast_wind,point_wind, hydro_max, pump_max, to_round,p,p_plus,p_minus);
        %######
        dados_opt=horzcat(wind_date((d-1)*np+1:(d-1)*np+24,1:2),dados_opt);
        data_opt=vertcat(data_opt,dados_opt);
        %######
        profit_opt=horzcat(wind_date((d-1)*np+1:(d-1)*np+1,1),profit_opt,taxas_opt);
        profits_opt=vertcat(profits_opt,profit_opt);
        %######
        dados=horzcat(wind_date((d-1)*np+1:(d-1)*np+24,1:2),dados);
        data=vertcat(data,dados);   
        
        if a>0.5
            filename_all=sprintf('%s%s%d%s%d','C:\Users\cpsilva\OneDrive\Tese\Software\Inicio\Results\Final_Results\','Single_AV_Scena_from_day_', sd, '_to_day_', d);
        elseif a<-0.5
            filename_all=sprintf('%s%s%d%s%d','C:\Users\cpsilva\OneDrive\Tese\Software\Inicio\Results\Final_Results\','Single_PR_Scena_from_day_', sd, '_to_day_', d);
        else
            filename_all=sprintf('%s%s%d%s%d','C:\Users\cpsilva\OneDrive\Tese\Software\Inicio\Results\Final_Results\','Single_EX_Scena_from_day_', sd, '_to_day_', d);
        end
        save(filename_all,'data','data_opt','profits_opt');
        
        diary off
end
