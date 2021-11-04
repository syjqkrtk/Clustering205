function [out check cg_result1 cg_result_g] = f_all_clustering_process_BIC_181016(CM,n_m,n_m_g) 






%% outlier pre-screening

if length(n_m) == 1
check = 0;
cg_result1 = n_m;
cg_result_g = n_m_g;
out = 0;

else
result1 = 0;
result2 = 0;
check = 1;
out = 0;

o = 0;
mean_CM = mean(mean(CM));

% 
CM_t = CM;

mean_CM = mean(mean(CM));

CM_outlier = sum(CM) - 1;

mean_CM_out = mean(CM_outlier);

sd = sqrt(var(CM_outlier));

CM_out_nd = (CM_outlier-mean_CM_out)/sd;


[a b] = find(CM_out_nd > 3 | CM_out_nd < -3);

outlier_i = b;

if length(b) > 0
    result1 = 1;
    
end
  
 

%% Clsutering 할지 말지?
    CM_test1 = 1-CM_t;
    CM_test2 = CM_t;
    CM_test3 = 1-CM_t;
    
    
    %% BIC k = 1
    
            [m u_i] = min(sum(CM_test1));
            R = length(CM);
            Rn = R;
            K = 1;      
            o2 = sum(CM_test1(:,u_i).^2)/(R-K);
            
            l_D = -Rn/2*log2(2*pi)-(Rn*R)/2*log2(o2)-(Rn-K)/2+Rn*log2(Rn)-Rn*log2(R);
            BIC = l_D-R/2*log2(R);
            
            
            
    %% BIC k = 2
    
      
       

C_n = 2;

li = length(CM_test1);


        


    
    
if li == 2
    
   CM_m{1} = 1;
   CM_m{2} = 1;
   
   
               for i = 1 : C_n
                   
                    [m cen(i)] = min(sum(CM_m{i}));
                    
                    
                    
                    
                     Rn = length(CM_m{i});
                     K = 2;      
                     o2 = sum(CM_m{i}(:,cen(i)).^2)/(R-K);
            
                     l_D = -Rn/2*log2(2*pi)-(Rn*R)/2*log2(o2)-(Rn-K)/2+Rn*log2(Rn)-Rn*log2(R);
                     BIC2_p(i) = l_D-R/2*log2(R);
                    
                    
               end
             
               BIC2 = mean(BIC2_p);

    
else

            [a1 b1] = min(min(CM_test2));
            [a2 b2] = min(CM_test2(b1,:));
           
             cl_c_p_p = [b1 b2];
             
             itable = [CM_test2(b1,:); CM_test2(b2,:)];
             itable(3,:) = 1:li;
             
             itable(:,cl_c_p_p) = [];
             
             for i = 1 : C_n
                    Cluster_result{i}(1) = cl_c_p_p(i);
             end
           
             [si1 si2] = size(itable);
             
             for i = 1 : si2
                if itable(1,i) > itable(2,i)
                    Cluster_result{1}(length(Cluster_result{1})+1) = itable(3,i);
                else
                    Cluster_result{2}(length(Cluster_result{2})+1) = itable(3,i);
                end
                 
             end
             
             
             n_ele_clu = 0;

             for i = 1 : C_n
                    CM_m{i} = 1-CM_test2(Cluster_result{i},Cluster_result{i});
                    [m cen(i)] = min(sum(CM_m{i}));
                    
                    
                    
                    
                     Rn = length(CM_m{i});
                     K = 2;      
                     o2 = sum(CM_m{i}(:,cen(i)).^2)/(R-K);
            
                     l_D = -Rn/2*log2(2*pi)-(Rn*R)/2*log2(o2)-(Rn-K)/2+Rn*log2(Rn)-Rn*log2(R);
                     BIC2_p(i) = l_D-R/2*log2(R);
                    
                    
             end
             
             
             BIC2 = mean(BIC2_p);

                
      

end


if BIC2 > BIC
    result2 = 1;
else
    result2 = 0;
end

% 
% r1 = abs(raw_result(1,2)-raw_result(1,6))/raw_result(1,2);
% r2 = abs(raw_result(1,2)-raw_result(2,6))/raw_result(1,2);
% 
% 
% r1_1 = double(r1 > ri)+double(r2 > ri);
% 
% result2 = double(r1_1 > 0);



%% clustering 

if result1 == 0 && result2 == 0
    
      cg_result1 = n_m;
      cg_result_g = n_m_g;
      check = 0;
      
elseif result1 == 1 && result2 == 0
   
    [cg_result1 cg_result_g] = f_clustering_algorithm1_real_180404_1(CM,n_m,n_m_g);
    out = 1;
elseif result1 == 0 && result2 == 1
    
    [cg_result1 cg_result_g] = f_clustering_algorithm1_real_180404_2(CM,n_m,n_m_g);
    
elseif result1 == 1 && result2 == 1
    out = 1;
    [cg_result1 cg_result_g] = f_clustering_algorithm1_real_180404_3(CM,n_m,n_m_g);

    
end





% end


end

end
