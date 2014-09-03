function [f,g] = nl_obj_fun_simple_utility(x,np,d,s,p,hydro_cost,pump_cost,a,k_profit,wind_max, n_var, hydro_max, e_max,C_min, C_max, L_min, L_max)                
L_min=L_min;
L_max=L_max;
f=0;
f_profit=sparse(s,1); 
utility_profit=zeros(s,1);

spill_cost=0;

for incre = 0:s-1
    for j=1:np
        f_profit(incre+1,1)=f_profit(incre+1,1)+(p((d-1)*np+j))*(x(0*np+j)+x((2*s+1)*np+j)+x((3*s+3+incre)*np+j))- x((2*s+3+incre)*np+j);
    end
    utility_profit(incre+1,1)=(1/s)*((1/(1-exp(a))*(1-exp((a*(f_profit(incre+1,1)-L_max)/(L_min-L_max))))));
end

for incre =1:s
    f=f+(utility_profit(incre,1));
end
f_profit;
utility_profit;
clear j
clear aux
clear incre


for j=1:np 
    aux_PH=0;
    aux_Pgrid=0;
    for s_aux=1:s
        aux_PH=aux_PH+(a*p((d-1)*np+j)*exp(-(a*(L_max - f_profit(s_aux,1)))/(L_min - L_max)))/(s*(exp(a) - 1)*(L_min - L_max));
        aux_Pgrid=aux_Pgrid+(a*p((d-1)*np+j)*exp(-(a*(L_max - f_profit(s_aux,1)))/(L_min - L_max)))/(s*(exp(a) - 1)*(L_min - L_max));
    end
    %aux_PH
    %aux_Pgrid
    g(0*np+j,1)=aux_PH; %Ph
    g((2*s+1)*np+j,1)=aux_Pgrid; %Pgrid    
    g((2*s+2)*np+j,1)=0; %Ti
    g((4*s+3)*np+j,1)=0; %pspi
    g((5*s+4)*np+j,1)=0; %PPi
    for incre=0:s-1
        g((1+incre)*np+j,1)=0;%-(a*pump_cost*exp(-(a*(L_max -f_profit(incre+1)))/(L_min - L_max)))/(s*(exp(a) - 1)*(L_min - L_max)); %PP
        g((s+1+incre)*np+j,1)=0; %Es      
        g((2*s+3+incre)*np+j,1)=-(a*exp(-(a*(L_max - f_profit(incre+1)))/(L_min - L_max)))/(s*(exp(a) - 1)*(L_min - L_max)); %Ts
        g((3*s+3+incre)*np+j,1)=(a*p((d-1)*np+j)*exp(-(a*(L_max - f_profit(s_aux,1)))/(L_min - L_max)))/(s*(exp(a) - 1)*(L_min - L_max)); %ds
        g((4*s+4+incre)*np+j,1)=0; %Psps
    end
end

end
