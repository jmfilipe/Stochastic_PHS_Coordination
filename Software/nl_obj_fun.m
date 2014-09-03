function [f,g] = nl_obj_fun(x,np,d,s,p,hydro_cost,pump_cost,a,k_profit,wind_max,n_var, hydro_max, e_max,C_min, C_max, L_min, L_max)                

f=0;
f_profit=sparse(s,1); 
utility_profit=zeros(s,1);
utility_cost=zeros(s,1);
f_cost=sparse(s,1);

k_costs=1-k_profit;  

for incre = 0:s-1
    for j=1:np  
        f_profit(incre+1,1) = f_profit(incre+1,1) + (p((d-1)*np+j))*(x(0*np+j)+x((2*s+1)*np+j)+x((3*s+3+incre)*np+j))- x((2*s+3+incre)*np+j);
        f_cost(incre+1,1)=f_cost(incre+1,1)+(x((4*s+4+incre)*np+j));
    end
    utility_profit(incre+1,1)=(L_max-f_profit(incre+1,1))/(L_max-L_min);
    utility_cost(incre+1,1)=1/(1-exp(a))*(1-exp((a*(f_cost(incre+1,1)-C_min)/(C_max-C_min))));
end

for incre = 1:s
    f=f+(1/s*(k_profit*utility_profit(incre,1)+k_costs*utility_cost(incre,1)));
end

clear j
clear aux

% g=[];
% soma_Pspi=0;
% for j=1:np
%     soma_Pspi=soma_Pspi + x(4*s+3)*np+j;
% end

for j=1:np
    g(0*np+j,1)=(s*k_profit*p((d-1)*np+j))/(s*(L_min - L_max));
    g(1*np+j:s*np+j,1)=0;
    g((s+1)*np+j:(2*s)*np+j,1)=0;    
    g((2*s+1)*np+j:(2*s+1)*np+j,1)=(s*k_profit*p((d-1)*np+j))/(s*(L_min - L_max));    
    g((2*s+2)*np+j:(2*s+2)*np+j,1)=0;    
    g((2*s+3)*np+j:(3*s+2)*np+j,1)=-k_profit/(s*(L_min - L_max));
    g((3*s+3)*np+j:(4*s+2)*np+j,1)=(k_profit*p((d-1)*np+j))/(s*(L_min - L_max)); %ds   
    g((4*s+3)*np+j:(4*s+3)*np+j,1)=0;
    for incre =0:s-1
            g((4*s+4+incre)*np+j,1)=-(a*k_costs*exp(-(a*(- C_min + f_cost(incre+1,1)))/(C_min - C_max)))/(s*(exp(a) - 1)*(C_min - C_max));
            
    end
    g((5*s+4)*np+j:(5*s+4)*np+j,1)=0;
end

end





    
    
