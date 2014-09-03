function stop = outfun(x,optimValues,state) 
    global fval_vector
    
    stop = false;
%     new_fval=optimValues.fval;
%     old_fval=fval_vector(1,1);
%     
%     if new_fval<old_fval 
%         fval_vector(1,2)=0;
%         fval_vector(1,1)=new_fval;
% %         disp('menor')
%     else
%         fval_vector(1,2)=fval_vector(1,2)+1;
%         fval_vector(1,1)=new_fval;
% %         disp('maior')
%     end
%     
%     if fval_vector(1,2) > 50
%         disp('ERROR: fval increasing')
%         stop = true;
%     end
%     
%     diary off

end