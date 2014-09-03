function h = hessianfcn_simple_utility(x,lambda,s,a,k_profit,np,L_min, L_max,n_var,d,p,hydro_cost,pump_cost)

f_profit=sparse(s,1);
spill_cost=0;
for incre = 0:s-1
    for j=1:np
        f_profit(incre+1,1)=f_profit(incre+1,1)+(p((d-1)*np+j))*(x(0*np+j)+x((2*s+1)*np+j))- x((2*s+3+incre)*np+j)-(x(0*np+j))*hydro_cost-x((1+incre)*np+j)*pump_cost-x((4*s+4+incre)*np+j)*spill_cost;
    end
    %utility_profit=utility_profit+(1/s)*((1/(1-exp(a))*(1-exp((a*(f_profit(incre+1,1)-L_max)/(L_min-L_max))))));
end

h=sparse(n_var,n_var);
% 
% price_row=1;
% price_column=1;
% 
% for row=1:np %percorre todas as linhas
% %     if price_row == np+1 %reset ao indice do preço, que so varia de 1 a np
% %         price_row = 1;
% %     end
%     for column=1:n_var %percorre todas as colunas
%         if price_column == np+1 %reset ao indice do preço, que so varia de 1 a np
%             price_column = 1;
%         end
%         if column <= (0)*np+np
%             aux_PH=0;
%             aux_Pgrid=0;
%             for s_aux=1:s
%                 aux_PH=aux_PH+(a^2*p((d-1)*np+row)*p((d-1)*np+column)*exp(-(a*(L_max - f_profit(s_aux,1)))/(L_min - L_max)))/(s*(exp(a) - 1)*(L_min - L_max)^2);
%                 aux_Pgrid=aux_Pgrid+(a^2*p((d-1)*np+row)*p((d-1)*np+column)*exp(-(a*(L_max - f_profit(s_aux,1)))/(L_min - L_max)))/(s*(exp(a) - 1)*(L_min - L_max)^2);
%             end
%             h((0)*np+row,(0)*np+column) = aux_PH; %PH
%             h((2*s+1)*np+row,(0)*np+column) = aux_Pgrid; %Pgrid
% 
%         elseif (column >= (2*s+1)*np+1) && (column <= (2*s+1)*np+np)
%             aux_PH=0;
%             aux_Pgrid=0;
%             for s_aux=1:s
%                 aux_PH=aux_PH+(a^2*p((d-1)*np+row)*p((d-1)*np+price_column)*exp(-(a*(L_max - f_profit(s_aux,1)))/(L_min - L_max)))/(s*(exp(a) - 1)*(L_min - L_max)^2);
%                 aux_Pgrid=aux_Pgrid+(a^2*p((d-1)*np+row)*p((d-1)*np+price_column)*exp(-(a*(L_max - f_profit(s_aux,1)))/(L_min - L_max)))/(s*(exp(a) - 1)*(L_min - L_max)^2);
%             end
%             h((0)*np+row,column) = aux_PH; %PH
%             h((2*s+1)*np+row,column) = aux_Pgrid; %Pgrid
%             
% %h((0)*np+row,(0)*np+column) =-(a^2*p((d-1)*np+row)*exp(-(a*(L_max -f_profix(xx)))/(L_min - L_max)))/(s*(exp(a) - 1)*(L_min - L_max)^2);
%         end
%         price_column=price_column+1;
%     end
%     price_row=price_row+1;
% end
% 
% for var = 0:s-1
%     h((2*s+3+var)*np+1:(2*s+3+var)*np+np,(2*s+3+var)*np+1:(2*s+3+var)*np+np)=(a^2*exp(-(a*(L_max -f_profit(var+1,1)))/(L_min - L_max)))/(s*(exp(a) - 1)*(L_min - L_max)^2); %Ts em relacao a Ts
% end
% clear var
% 
% for j=1:np
%     for var = 0:s-1
%         h((0)*np+j,(2*s+3+var)*np+1:(2*s+3+var)*np+np) =-(a^2*p((d-1)*np+j)*exp(-(a*(L_max -f_profit(var+1,1)))/(L_min - L_max)))/(s*(exp(a) - 1)*(L_min - L_max)^2);%PH em relação a Ts
%         h((2*s+1)*np+j,(2*s+3+var)*np+1:(2*s+3+var)*np+np) =-(a^2*p((d-1)*np+j)*exp(-(a*(L_max -f_profit(var+1,1)))/(L_min - L_max)))/(s*(exp(a) - 1)*(L_min - L_max)^2);%Pgrid em relação a TS        
%     end    
% end

end


