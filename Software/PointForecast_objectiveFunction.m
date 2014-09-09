function [f,g] = PointForecast_objectiveFunction(x,d,n_var,system_data, wind_data, market_data, main_data)                
np=main_data.np; %n. periods 
s=1; %n. scenarios
f=0;
% f_profit=0; 
f_profit=sparse(s,1); 


for incre = 0:s-1
    for j=1:np  
           f_profit(incre+1,1)=f_profit(incre+1,1)+(market_data.price_forecast((d-1)*np+j))*(x((0*s+1+incre)*np+j)+x((1*s+1+incre)*np+j)+x((4*s+1+incre)*np+j))-x((5*s+1+incre)*np+j);
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
        g((0*s+1+incre)*np+j,1)=-market_data.price_forecast((d-1)*np+j)/s;          %PH
        g((1*s+1+incre)*np+j,1)=-market_data.price_forecast((d-1)*np+j)/s;          %PG
        g((2*s+1+incre)*np+j,1)=0;                                                  %PP
        g((3*s+1+incre)*np+j,1)=0;                                                  %E
        g((4*s+1+incre)*np+j,1)=-market_data.price_forecast((d-1)*np+j)/s;          %d
        g((5*s+1+incre)*np+j,1)=1/s;                                                %T
        g((6*s+1+incre)*np+j,1)=0;                                                  %P_DL
    end
end

end
