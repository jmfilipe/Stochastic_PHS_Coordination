function h = hessianfcn(x,lambda,s,a,k_profit,np,C_min, C_max,n_var)
k_costs=1-k_profit;  
f_cost=sparse(s,1);

for incre = 0:s-1
    for j=1:np  
        %f_profit(incre+1,1) = f_profit(incre+1,1) + (p((d-1)*np+j))*(x(0*np+j)+x((2*s+1)*np+j))- x((2*s+3+incre)*np+j)-(x(0*np+j))*hydro_cost-x((1+incre)*np+j)*pump_cost;
        f_cost(incre+1,1)=f_cost(incre+1,1)+(x((4*s+4+incre)*np+j));
    end
end
h=sparse(n_var,n_var);

% for j=1:np
%     h((4*s+3)*np+j,(4*s+3)*np+j)=(a^2*k_costs*exp(-(a*(- C_min + f_cost))/(C_min - C_max)))/((exp(a) - 1)*(C_min - C_max)^2);
% end

% for var = 0:s-1
%     h((2*s+3+var)*np+1:(2*s+3+var)*np+np,(2*s+3+var)*np+1:(2*s+3+var)*np+np)=(a^2*k_costs*exp(-(a*(- C_min+f_cost(var+1,1)))/(C_min - C_max)))/(s*(exp(a) - 1)*(C_min - C_max)^2);
% end

% for var = 0:s-1
%     h((2*s+3+var)*np+1:(2*s+3+var)*np+np,(2*s+3+var)*np+1:(2*s+3+var)*np+np)=(a^2*k_costs*exp((a*(- C_min+f_cost(var+1,1)))/(C_min - C_max)))/(s*(exp(-a) - 1)*(C_min - C_max)^2);
% end

end

