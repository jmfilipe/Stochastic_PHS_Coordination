function [f,g] = nl_obj_fun_EV(x,np,d,s,p,hydro_cost,pump_cost,a,k_profit,wind_max, n_var, hydro_max, e_max,C_min, C_max, L_min, L_max)                

f=0;
% f_profit=0; 
f_profit=sparse(s,1); 


for incre = 0:s-1
    for j=1:np  
%         f_profit(incre+1,1)=f_profit(incre+1,1)+(p((d-1)*np+j))*(x(0*np+j)+x((2*s+1)*np+j)+x((3*s+3+incre)*np+j))-x((2*s+3+incre)*np+j);
           f_profit(incre+1,1)=f_profit(incre+1,1)+(p((d-1)*np+j))*(x((0*s+1+incre)*np+j)+x((1*s+1+incre)*np+j)+x((4*s+1+incre)*np+j))-x((5*s+1+incre)*np+j);
    end
end

for incre = 1:s
    f=f+1/s*((f_profit(incre,1)));  
end


f=-f;

clear j
clear aux
g=zeros(n_var,1);
for j=1:np
    g((0)*np+j,1)=0; %BID
    for incre=0:(s-1)
        g((0*s+1+incre)*np+j,1)=-p((d-1)*np+j)/s;%PH
        g((1*s+1+incre)*np+j,1)=-p((d-1)*np+j)/s;%PG
        g((2*s+1+incre)*np+j,1)=0; %PP
        g((3*s+1+incre)*np+j,1)=0; %E
        g((4*s+1+incre)*np+j,1)=-p((d-1)*np+j)/s;%d
        g((5*s+1+incre)*np+j,1)=1/s;%T
        g((6*s+1+incre)*np+j,1)=0; %P_DL
    end
end

end
