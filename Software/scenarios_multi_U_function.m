function [data, data_opt,profits_opt] = scenarios_multi_U_function( nd, sd, np, t, s, prob, a, k_profit, hydro_max, pump_max, pump_eff, hydro_eff, e_max, hydro_cost, pump_cost,pw, p, p_plus, p_minus, C_max, C_min, L_max, L_min,wind_max, arredondar, to_round, wind_date,real_wind,forecast_wind,point_wind)
data=[]; %Data obtained from the day-ahead strategy
data_opt=[]; %Data obtained from the operational management strategy
profits_opt=[];
e_begin = 0.0*e_max;     %initial & final level of the reservoir (MW) 
global fval_vector       %global vector for the output function

%% ALGORITH
for d = sd:nd   %simulates for the time horizon defined - nd: n. of days
    tStart = tic;
    %==================================================================
    % Optimization
    %==================================================================

    % LINEAR PROGRAM PARAMETERS

    n_var = (4+5*s+1)*np;            % num. variables 

    % equal equations
    neq = (5*s+3+1)*np;           %n equal equations
    Aeq = sparse(neq,n_var); 
    beq = sparse(1,neq);

    % inequalities
    nde = (4*s-1+1+2)*np;             %n inequal equations              
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
                
                %Eq. 0: Regulation costs of period i
                    Aeq(0*np+j,(2*s+2)*np+j) = 1;         %Ti
                for Ts_Ts=(2*s+3):(3*s+2)
                    Aeq(0*np+j,Ts_Ts*np+j) = -1*prob;    %Ts,i
                end
                clear Ts_Ts
                
                %Eq. 1 Hydro Pumped Storage 
                if j == 1 %1st hour
                    
                    %Eq. 1 Hydro Pumped Storage (1)
                    q_q=1;
                    for Es_Es=(s+1):(2*s)
                        Aeq(q_q*np+j,Es_Es*np+j) = 1;
                        beq(q_q*np+j) = e_begin;                
                        q_q=q_q+1;
                    end
                    clear Es_Es
                    clear q_q
                    
                    %Eq. 2 Hydro Pumped Storage (2)
                    Aeq((s+1)*np+j,0*np+j) = 1;
                    beq((s+1)*np+j) = hydro_eff/t * e_begin;        
                
                %Eq. 3 Hydro Pumped Storage (3)
                else %remaining hours
                    q_q=s+2;
                    Es_Es=s+1;
                    for Pps_Pps=1:s
                        Aeq(q_q*np+j,0*np+j-1) = t/hydro_eff;
                        Aeq(q_q*np+j,Pps_Pps*np+j-1) = -t*pump_eff;                    
                        Aeq(q_q*np+j,Es_Es*np+j) = 1;
                        Aeq(q_q*np+j,Es_Es*np+j-1 ) = -1;
                        q_q=q_q+1;
                        Es_Es=Es_Es+1;
                    end
                end
                clear Pps_Pps
                clear Es_Es
                clear q_q                
                
                %Eq. 4: Imbalance for every scenario 
                q_q=2*s+2;
                Pg_Pg=2*s+1;
                d_d=3*s+3;
                Psps_Psps=4*s+4;
                for Pps_Pps=1:s
                    Aeq(q_q*np+j,Pps_Pps*np+j) = 1;
                    Aeq(q_q*np+j,Pg_Pg*np+j) = 1;
                    Aeq(q_q*np+j,d_d*np+j) = 1;
                    Aeq(q_q*np+j,Psps_Psps*np+j) = 1;
                    beq(q_q*np+j) = pw((d-1)*np+j,Pps_Pps);
                    d_d=d_d+1;
                    Psps_Psps=Psps_Psps+1;
                    q_q=q_q+1;
                end
                clear Pps_Pps
                clear Pg_Pg
                clear d_d
                clear Psps_Psps
                clear q_q
                
               %Eq. 5: Spill for period i
                    Aeq((3*s+2)*np+j,(4*s+3)*np+j) = 1;         %Pspi
                for Psps_Psps=(4*s+4):(5*s+3)
                    Aeq((3*s+2)*np+j,Psps_Psps*np+j) = -1*prob;    %Psps,i
                end
                clear Psps_Psps
                
                %Eq. 6: Pump Power for period i
                    Aeq((3*s+3)*np+j,(5*s+4)*np+j) = 1;         %Ppi
                for Pps_Pps=(1):(s)
                    Aeq((3*s+3)*np+j,Pps_Pps*np+j) = -1*prob;    %Pps,i
                end
                clear Pps_Pps
                
                %Eq. 7: Pump Power for period 24, PP_24 = 0
                if j == np
                    q_q=3*s+4;
                    for Pps_Pps=1:s
                         Aeq((q_q)*np+j,Pps_Pps*np+j) = 1;
                         q_q=q_q+1;
                    end
                end
                clear Pps_Pps
                clear q_q
                
                %Eq. 8      : Energy for period 24, E_24 = 0
                if j == np
                    q_q=4*s+4;
                    for E_E=s+1:2*s
                         Aeq((q_q)*np+j,E_E*np+j) = 1;
                         Aeq((q_q)*np+j,0*np+j) = -1/hydro_eff;
                         q_q=q_q+1;
                    end
                end
                clear E_E
                clear q_q
                               
                %==============================================================
                %INEQUAL EQ.
                %==============================================================

                %Ineq. 0: Hydro Generation
                q_q=0;
                Es_Es=s+1;
                for Pps_Pps=1:s
                    Ad(q_q*np+j,0*np+j) = 1;
                    Ad(q_q*np+j,Es_Es*np+j) = -hydro_eff/t;  
                    Ad(q_q*np+j,Pps_Pps*np+j) = -t*hydro_eff*pump_eff;
                    q_q=q_q+1;
                    Es_Es=Es_Es+1;
                end
                clear Pps_Pps
                clear Es_Es
                clear q_q
                
               
                %Ineq. 1: Epigraph Form (1)
                q_q=s;
                d_d=3*s+3;
                for Ts_Ts=(2*s+3):(3*s+2)
                    Ad(q_q*np+j,(Ts_Ts)*np+j) = -1;           %Ts,i
                    Ad(q_q*np+j,(d_d)*np+j) = (p((d-1)*np+j)-p_plus((d-1)*np+j));     %d,si
                    bd((q_q)*np+j)=0;
                
                    %Ineq. 2: Epigraph Form (2)
                    Ad((q_q+s)*np+j,(Ts_Ts)*np+j) = -1;           %Ts,i
                    Ad((q_q+s)*np+j,(d_d)*np+j) = -(p_minus((d-1)*np+j)-p((d-1)*np+j));     %d,si
                    bd((q_q+s)*np+j)=0; 
                    
                    d_d=d_d+1;
                    q_q=q_q+1;
                end
                clear q_q
                clear d_d
                clear Ts_Ts
                
                % Ineq. 3: Hydro + Pump Constraint
                q_q=3*s;
                for Pps_Pps=1:s
                    Ad(q_q*np+j,0*np+j) = 1;
                    Ad(q_q*np+j,Pps_Pps*np+j) = 1;
                    bd(q_q*np+j) = hydro_max/e_max;
                    q_q=q_q+1;
                end
                clear q_q
                clear Pps_Pps
                
%                 %Ineq. 4: Wind Power Constraint
%                 q_q=4*s;
%                 Pg_Pg=2*s+1;
%                 Psps_Psps=4*s+4;
%                 for Pps_Pps=1:s
%                     Ad(q_q*np+j,Pps_Pps*np+j) = -1;
%                     Ad(q_q*np+j,Pg_Pg*np+j) = -1;
%                     Ad(q_q*np+j,Psps_Psps*np+j) = -1;
%                     bd(q_q*np+j)=-pw((d-1)*np+j,Pps_Pps);
%                     Psps_Psps=Psps_Psps+1;
%                     q_q=q_q+1;
%                 end
%                 clear Psps_Psps
%                 clear Pg_Pg
%                 clear q_q
%                 clear Pps_Pps

%                 
                % Ineq. 5: Pump + Spill Constraint
                q_q=4*s;
                Psps_Psps=4*s+4;
                for Pps_Pps=1:s
                    Ad(q_q*np+j,Psps_Psps*np+j) = 1;
                    Ad(q_q*np+j,Pps_Pps*np+j) = 1;
                    bd(q_q*np+j) = pw((d-1)*np+j,Pps_Pps);
                    q_q=q_q+1;
                    Psps_Psps=Psps_Psps+1;
                end
                clear q_q
                clear Pps_Pps
                clear Psps_Psps
                
%                 Ineq. 6: Objective Function Bounds (1)
                
                Ad((5*s)*np+j,0*np+j) = p((d-1)*np+j);
                Ad((5*s)*np+j,(2*s+1)*np+j) = p((d-1)*np+j);
                Ad((5*s)*np+j,(2*s+2)*np+j) = -1;
                bd((5*s)*np+j) = L_max/24;
                
%                 Ineq. 7: Objective Function Bounds (1)
                
                Ad((5*s+1)*np+j,0*np+j) = -p((d-1)*np+j);
                Ad((5*s+1)*np+j,(2*s+1)*np+j) = -p((d-1)*np+j);
                Ad((5*s+1)*np+j,(2*s+2)*np+j) = 1;
                bd((5*s+1)*np+j) = -L_min/24;
               
                %==============================================================
                %BOUNDS
                %==============================================================

                %upper bound
                ub(0*np+j) = hydro_max/e_max;         %PHi
                ub((2*s+1)*np+j) = max(pw((d-1)*np+j,1:s));    %Pgrid
                ub((2*s+2)*np+j) = Inf;         %Ti
                ub((4*s+3)*np+j) = Inf;         %Pspi
                ub((5*s+4)*np+j) = Inf;         %Ppi
                %ub((5*s+4)*np+j) = Inf;     %Di
                for incre=1:s                         
                    ub((incre)*np+j) = pump_max/e_max;      %Pps                        
                    ub((s+incre)*np+j) = e_max/e_max;       %Es 
                    ub((2*s+2+incre)*np+j) = Inf;     %Ts
                    ub((3*s+2+incre)*np+j) = Inf;     %ds
                    ub((4*s+3+incre)*np+j) = Inf/e_max;     %Psps
                    %ub((5*s+5+incre)*np+j) = Inf;     %Ds
                    
                end
                clear incre
                
                %lower bound
                for incre=1:s
                    lb((3*s+2+incre)*np+j) = -Inf;    %ds
                end
                clear incre
            end
            
            
        disp_message=['Additive Function, day: ', num2str(d)];
        disp(disp_message)
        
         % FUNCAO INICIAL 
         x0= sparse(n_var,1);
    
        options = optimset('Algorithm','interior-point','Display','iter-detailed','GradObj','on','DerivativeCheck','off','FunValCheck','on','UseParallel','always',...
    'Hessian','on','HessFcn',@(x,lambda)hessianfcn(x,lambda,s,a,k_profit,np,C_min, C_max,n_var),'OutputFcn',@outfun);%notify iter
            
        options.MaxFunEvals = (s*5000); 
        options.MaxIter = 250;
        options.TolFun = 1e-9;
        options.TolCon = 1e-7;
        options.TolX = 1e-9;
        
        [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x)nl_obj_fun(x,np,d,s,p,hydro_cost,pump_cost,a,k_profit,wind_max,n_var, hydro_max, e_max, C_min, C_max, L_min, L_max),x0,Ad,bd,Aeq,beq,lb,ub,[],options);
        %exitflag
        
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
                    options.MaxFunEvals = (s*1000); 
                options_lin.MaxIter = 50;
                options_lin.TolFun = 1e-8;
                options_lin.TolCon = 1e-8;
                options_lin.TolX = 1e-8;
            xnew = fmincon(@(x)initial_x0(x0),x0,Ad,bd,Aeq,beq,lb,ub,[],options_lin);
            x0= xnew;

            %SEGUNDA FUNCAO USANDO O PONTO INICIAL
            
            options = optimset('Algorithm','interior-point','Display','iter-detailed','GradObj','on','DerivativeCheck','off','FunValCheck','on','UseParallel','always',...
    'Hessian','on','HessFcn',@(x,lambda)hessianfcn(x,lambda,s,a,k_profit,np,C_min, C_max,n_var),'OutputFcn',@outfun);%notify iter

            options.MaxFunEvals = (s*100); 
            options.MaxIter = 250;
            options.TolFun = 1e-9;
            options.TolCon = 1e-7;
            options.TolX = 1e-9;
            fval_vector=[9999999 0];
            [x,fval,exitflag,output,lambda,grad,hessian] = fmincon(@(x)nl_obj_fun(x,np,d,s,p,hydro_cost,pump_cost,a,k_profit,wind_max,n_var, hydro_max, e_max, C_min, C_max, L_min, L_max),x0,Ad,bd,Aeq,beq,lb,ub,[],options);
            
            if fval_error < fval
                x=x_error;
                fprintf ('************\n...isto e estranho\n**********\n');
            else
                fprintf ('************\n...tudo na boa\n**********\n');
            end
        
        end
        

        
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
        [dados_opt,profit_opt,taxas_opt]=operational_strat( np,d, day_ahead_data, e_max, e_begin, hydro_eff, pump_eff, t, real_wind, forecast_wind,point_wind, hydro_max, pump_max, to_round,p,p_plus,p_minus);
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
            filename_all=sprintf('%s%s%d%s%d','C:\Users\cpsilva\OneDrive\Tese\Software\Inicio\Results\Final_Results\','Multi_AV_Scena_from_day_', sd, '_to_day_', d);
        elseif a<-0.5
            filename_all=sprintf('%s%s%d%s%d','C:\Users\cpsilva\OneDrive\Tese\Software\Inicio\Results\Final_Results\','Multi_PR_Scena_from_day_', sd, '_to_day_', d);
        else
            filename_all=sprintf('%s%s%d%s%d','C:\Users\cpsilva\OneDrive\Tese\Software\Inicio\Results\Final_Results\','Multi_EX_Scena_from_day_', sd, '_to_day_', d);
        end
        save(filename_all,'data','data_opt','profits_opt');
        
end
