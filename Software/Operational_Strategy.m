function [ final_data,profit, taxas, final_storage, x,fval,exitflag,output,lambda ] = Operational_Strategy( d,dados,system_data, wind_data, market_data, main_data)
np=main_data.np; %n. periods 

%% INITIALIZATION
day_ahead=dados;
final_data=zeros(np,7);

    for h=1:np
        %==================================================================
        % Optimization
        %==================================================================

        % LINEAR PROGRAM PARAMETERS

        n_var = (7)*np;            % num. variables 
        
        % objective funcion vector
        f=sparse(1,n_var);
        
        % equal equations
        neq = (8)*np;           %n equal equations
        Aeq = sparse(neq,n_var); 
        beq = sparse(1,neq);

        % inequalities
        nde = (4)*np;             %n inequal equations              
        Ad = sparse(nde,n_var);
        bd = sparse(1,nde);

        % variable limits
        lb = sparse(1,n_var);     %lower limits
        ub = sparse(1,n_var);     %upper limits
        
        for j=h:np
            
            %==============================================================
            % O.F.
            %==============================================================

            if j==h
                f(4*np+j) = 100;           
                f(6*np+j) = 99;
                f(0*np+j) = 99;
            else
                f(4*np+j) = 1;
                f(6*np+j) = 0.99;
                f(0*np+j) = 0.99;
            end
                      
            %==============================================================
            %EQUAL EQ.
            %==============================================================

            %eq. 0: Real - Forecasted
            Aeq(0*np+j,0*np+j) = 1;
            Aeq(0*np+j,3*np+j) = 1;
            Aeq(0*np+j,5*np+j) = -1;
            %Aeq(0*np+j,6*np+j) = -1;
            beq(0*np+j) = day_ahead(j,1);

            %Eq. 1Hydro Pumped Storage
            if j == h %1st hour of operational strategy
                if j==1 %first hour of the day
                    Aeq(1*np+j,2*np+j) = 1;
                    beq(1*np+j) = system_data.final_storage;                %initial level has to be e_begin
%                     Aeq(3*np+j,0*np+j) = 1;
%                     beq(3*np+j) = hydro_eff/t * e_begin;
                   
                else
                    Aeq(5*np+j,2*np+j) = 1;
                    %Aeq(5*np+j,6*np+j) = 1;
                    beq(5*np+j)=final_data(h-1,3)+main_data.t*system_data.pump_eff*final_data(h-1,2)-main_data.t/system_data.hydro_eff*final_data(h-1,1);
                end                      

            else %remaining hours                   
                Aeq(2*np+j,1*np+j-1) = -main_data.t*system_data.pump_eff;
                Aeq(2*np+j,0*np+j-1) = main_data.t/system_data.hydro_eff;
                Aeq(2*np+j,2*np+j) = 1;
                Aeq(2*np+j,2*np+j-1) = -1;
                %Aeq(2*np+j,6*np+j-1) = 1; 
            end

            %eq. 4 Wind Power constraint
            Aeq(4*np+j,3*np+j) = 1;
            Aeq(4*np+j,1*np+j) = 1;
            Aeq(4*np+j,6*np+j) = 1;
            if j==h
                beq(4*np+j)=wind_data.real_wind((d-1)*np+h);
            elseif (j>h && j<(h+5))
                beq(4*np+j)=wind_data.updated((d-1)*np+h,j-h);
            else
                beq(4*np+j)=wind_data.point_forecast((d-1)*np+h);
            end
            
            %Eq. 5: Energy for period 24, E_24 = 0
            if j == np
%                  Aeq(6*np+j,2*np+j) = 1;
%                  Aeq(6*np+j,0*np+j) = -1/hydro_eff;
%                  Aeq(6*np+j,1*np+j) = 1*pump_eff;
            end
            
            %Eq. 5: Pump for period 24, E_24 = 0
            if j == np
                 %Aeq(5*np+j,2*np+j) = 1;
                 %Aeq(5*np+j,0*np+j) = -1/hydro_eff;
%                  Aeq(7*np+j,1*np+j) = 1;
            end
            %==============================================================
            %INEQUAL EQ.
            %==============================================================

            %ineq. 0:Hydro Generation Operation
            Ad(0*np+j,0*np+j) = 1;
            Ad(0*np+j,2*np+j) = -system_data.hydro_eff/main_data.t;  
            %Ad(0*np+j,1*np+j) = -t*hydro_eff*pump_eff;

            %ineq. 1:Imbalances p+*d-T=<0
            Ad(1*np+j,(4)*np+j) = -1;           %Ti
            Ad(1*np+j,(5)*np+j) = 1;%*(p((d-1)*np+j)-p_plus((d-1)*np+j));     %di
            %ineq. 2:Imbalances p-*d-T=<0
            Ad(2*np+j,(4)*np+j) = -1;           %Ti
            Ad(2*np+j,(5)*np+j) = -1;%*(p_minus((d-1)*np+j)-p((d-1)*np+j));     %di

            %max hydro
            Ad(3*np+j,0*np+j) = 1;
            Ad(3*np+j,2*np+j) = -system_data.hydro_eff/main_data.t;
            Ad(3*np+j,1*np+j) = -main_data.t*system_data.hydro_eff*system_data.pump_eff;

            %pump/hydro coordination per hour
            Ad(4*np+j,0*np+j) = 1;
            Ad(4*np+j,1*np+j) = 1;
            bd(4*np+j) = system_data.hydro_max;


            %==============================================================
            %BOUNDS
            %==============================================================

            %upper bound 
            ub(0*np+j) = system_data.hydro_max;                           
            ub(1*np+j) = system_data.pump_max;                        
            ub(2*np+j) = system_data.e_max;                             
            ub(3*np+j) = Inf;                           
            ub(4*np+j) = Inf;
            ub(5*np+j) = Inf;  
            ub(6*np+j) = Inf;
            
            %lower bound
            lb(5*np+j) = -Inf;
        end
        
%         options = optimset('Display','iter','MaxIter',10000,'TolFun',.000001, 'TolCon', .000001, ...
%           'TolX',.000001, 'LargeScale', 'off','simplex','on');
        options = optimset('Display','off','MaxIter',10000,'TolFun',.00001, 'TolCon', .00001, ...
          'TolX',.00001, 'LargeScale', 'on','simplex','off');
        %options = optimset('Display','final','MaxIter',10000,'LargeScale', 'off','simplex','on');
        % Optimization Funcion
        [x,fval,exitflag,output,lambda] = linprog(f,Ad,bd,Aeq,beq,lb,ub,[],options);
        xx = x;%round(x*10000)/10000; 
        
        final_data(h,1)=x(0*np+h);%PH
        final_data(h,2)=x(1*np+h);%PP
        final_data(h,3)=x(2*np+h);%E
        final_data(h,4)=x(3*np+h);%Grid
        final_data(h,5)=x(4*np+h);%T
        final_data(h,6)=x(5*np+h);%d
        final_data(h,7)=x(6*np+h);%spill
        
        
    end
    
    final_data(:,6)=final_data(:,6)+final_data(:,7); %accounts for positive imbalances
%     final_data(:,6)=-final_data(:,6);
    arredondar=1;
    if arredondar == 1    
        [r,c] = size (final_data);
        for j=1:c
            for i=1:r
                if(abs(final_data(i,j))<main_data.to_round)
                    final_data(i,j)=0;
                end
            end
        end
    end
    
    %profit=sum(final_data(:,1))+sum(final_data(:,4))-sum(final_data(:,5));
%     p((d-1)*np+j)
%     p_plus((d-1)*np+j)
%     p_minus((d-1)*np+j)
    
    profit=0;
    lucro=0;
    taxas=0;
    taxas_2=zeros(np,1);
    
    %adds the positive imbalance to the power injected into the grid
    for j=1:24
        if final_data(j,6) > 0 
            final_data(j,4)=final_data(j,4)+final_data(j,6);
        end
    end
        
        
    for k=1:np
       lucro=(final_data(k,1)+final_data(k,4))*market_data.price_real((d-1)*np+k);
       if final_data(k,6)<0
           taxas_2(k,1)=-(final_data(k,6))*market_data.p_minus_real((d-1)*np+k);
       elseif final_data(k,6)>0
           taxas_2(k,1)=(final_data(k,6))*market_data.p_plus_real((d-1)*np+k);        
       else
           taxas_2(k,1)=0;
       end
       profit = profit + lucro - taxas_2(k,1); 
       taxas = taxas + taxas_2(k,1);
    end
    
    final_storage=final_data(24,3)+final_data(24,2)*system_data.pump_eff-final_data(24,1)/system_data.hydro_eff;
    
    final_data(:,5)=final_data(:,6);
    final_data(:,7)=[];
    final_data(:,6)=taxas_2;
end

