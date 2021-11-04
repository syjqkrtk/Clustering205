function [out check cg_result1 cg_result_g] = f_all_clustering_process_out_181017(CM,n_m,n_m_g,ri) 



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
    CM_test1 = CM_t;
    CM_test2 = CM_t;
    CM_test3 = CM_t;
    
    CM_t2 = CM_t.^2;
    
    
    size_CM = length(CM_t2)*length(CM_t2);
    
    
    v1_p = (sum(sum(CM_t2))-length(CM_t2))/(size_CM-length(CM_t2));
    v2_p = ((sum(sum(CM_t))-length(CM_t))/(size_CM-length(CM_t)))^2;
    
    var_clu_all2 = abs(sqrt(v1_p-v2_p));
    
     mean_CM_all = mean(mean(CM_t));
     var_clu_all = sqrt(var(double(CM_t(:))));                  
       

C_n = 2;

li = length(CM_test1);


raw_result(1,1) = length(CM_test1);
raw_result(1,2) = mean(mean(CM_test1));
raw_result(1,3) = var(double(CM_test1(:)));

        


    
    
if li == 2
    
   CM_m{1} = 1;
   CM_m{2} = 1;
   
   
             for i = 1 : C_n
 
                    
                    raw_result(i,5) = length(CM_m{i});
                    raw_result(i,6) = mean(mean(CM_m{i}));
                    raw_result(i,7) = var(double(CM_m{i}(:)));
             
             
                    
             end
   

    
else

            [a1 b1] = min(min(CM_test1));
            [a2 b2] = min(CM_test1(b1,:));
           
             cl_c_p_p = [b1 b2];
             
             itable = [CM_test1(b1,:); CM_test1(b2,:)];
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
                    CM_m{i} = CM_test1(Cluster_result{i},Cluster_result{i});
             end
             
             
             for i = 1 : C_n
 
                    
                    raw_result(i,5) = length(CM_m{i});
                    raw_result(i,6) = mean(mean(CM_m{i}));
                    raw_result(i,7) = var(double(CM_m{i}(:)));
             
             
                    
             end
                
      

end


r1 = abs(raw_result(1,2)-raw_result(1,6))/raw_result(1,2);
r2 = abs(raw_result(1,2)-raw_result(2,6))/raw_result(1,2);


r1_1 = double(r1 > ri)+double(r2 > ri);

result2 = double(r1_1 > 0);



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





end


end
    
