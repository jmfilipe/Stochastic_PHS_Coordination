function [f,g] = nl_obj_fun_EV(x,np,d,s,p,hydro_cost,pump_cost,a,k_profit,wind_max, n_var, hydro_max, e_max,C_min, C_max, L_min, L_max)                

f=0;
% f_profit=0; 
f_profit=sparse(s,1); 


for incre = 0:s-1
    for j=1:np  
%         f_profit(incre+1,1)=f_profit(incre+1,1)+(p((d-1)*np+j))*(x(0*np+j)+x((2*s+1)*np+j)+x((3*s+3+incre)*np+j))-x((2*s+3+incre)*np+j);
           f_profit(incre+1,1)=f_profit(incre+1,1)+(p((d-1)*np+j))*(x(0*np+j)+x((2*s+1)*np+j)+x((3*s+3+incre)*np+j))-x((2*s+3+incre)*np+j);
    end
end
f_profit;
for incre = 1:s
    f=f+1/s*((f_profit(incre,1)));
    
end


f=-f;

clear j
clear aux

for j=1:np
    g(0*np+j,1)=- p((d-1)*np+j);
    g(1*np+j:s*np+j,1)=0;
    g((s+1)*np+j:(2*s)*np+j,1)=0;    
    g((2*s+1)*np+j:(2*s+1)*np+j,1)=-p((d-1)*np+j);    
    g((2*s+2)*np+j:(2*s+2)*np+j,1)=0;
    for incre=0:(s-1)
    g((2*s+3+incre)*np+j,1)=1/s; %Ts
    g((3*s+3+incre)*np+j,1)=-p((d-1)*np+j)/s;
    end
    
    g((4*s+3)*np+j:(4*s+3)*np+j,1)=0;
    g((4*s+4)*np+j:(5*s+3)*np+j,1)=0;
    g((5*s+4)*np+j:(5*s+4)*np+j,1)=pump_cost;
    
end

end
