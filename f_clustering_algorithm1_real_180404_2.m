function [result1 result2] = f_clustering_algorithm1_real_180404_2(CM,n_m,n_m_g)

l = length(n_m);

if l == 2
    result1{1} = n_m(1);
    result1{2} = n_m(2);
    
    result2{1} = n_m_g(1);
    result2{2} = n_m_g(2);
    
else
    

%% 중심 노드 결정
    CM_test1 = CM;
    CM_test2 = CM;
    CM_test3 = CM;
    CM_test4 = CM;
    
    mean_CM_all = mean(mean(CM));
    var_clu_all = var(double(CM(:)));                  
       
       li = length(CM);
%        li2 = round(li*(li-1)*0.02/2);
%      
       li2 = round(li*0.02);
       m_C_n = round(li/5);
       li3 = ceil(li*0.02);
       
       if m_C_n < 2
          
           m_C_n = 2;
           
       end
       
       

%        C_n = 2;
 
       if li2 == 0
           for C_n = 2 : m_C_n
               
               if C_n == 2
                   [a1 b1] = min(min(CM_test2));
                   [a2 b2] = min(CM_test2(b1,:));
                   
                   cl_c_p_p(1) = b1;
                   cl_c_p_p(2) = b2;
                   
                   CM_test2(b1,b2) = 1;
                   CM_test2(b2,b1) = 1;
                   
                   for i = 1 : C_n
                       Cluster_result{i}(1) = cl_c_p_p(i);
                   end
                   
                   complete_n = cl_c_p_p;
                   l_complete = length(complete_n);
                   
         
                   
                   while l_complete ~= l
                       clear test_c
                       for i = 1 : C_n
                           test_c{i}(1,:) = zeros(1,l);
                           
                           for j = 1 : length(Cluster_result{i})
                               test_c{i}(1,:) = test_c{i}(1,:)+CM(Cluster_result{i}(j),:);
                           end
                           
                           test_c{i}(1,:) = test_c{i}(1,:)/length(Cluster_result{i});
                           
                           test_c{i}(2,:) = 1:1:l;
                       end
                       
                       
                       
                       for i = 1 : C_n
                           test_c{i}(:,complete_n) = [];
                           [a b] = max(test_c{i}(1,:));
                           test_c2{i}(1) = a;
                           test_c2{i}(2) = test_c{i}(2,b);
                           test_c3(i) = a;
                       end
                       
                       
                       [a b] = max(test_c3);
                       
                       index = test_c2{b}(2);
                       Cluster_result{b}(length(Cluster_result{b})+1) = index;
                       complete_n(length(complete_n)+1) = index;
                       
                       l_complete = length(complete_n);
                       clear test_ccategorization
                       
                   end
                   
                   
                   
                   Cluster_result_s{C_n} = Cluster_result;
                   
                   
                   
                   clear Cluster_result
               else
                   
                   CM_test2 = CM;
                    for i = 1 : length(cl_c_p_p)
                        for j = 1 : length(cl_c_p_p)
                            if i ~= j
                                CM_test2(cl_c_p_p(i),cl_c_p_p(j)) = 1;
                                CM_test2(cl_c_p_p(j),cl_c_p_p(i)) = 1;
                            end
                        end
                   end
                   
                   
                   t = zeros(1,l);
                   for i = 1 : length(cl_c_p_p)
                       t = t + CM_test2(cl_c_p_p(i),:);
                   end
                   
                   
                   [a3 b3] = min(t);
                   
                   for i = 1 : length(cl_c_p_p)
                       CM_test2(cl_c_p_p(i),b3) = 1;
                       CM_test2(b3,cl_c_p_p(i)) = 1;
                   end
                   
                   cl_c_p_p(C_n) = b3;
                   
                   for i = 1 : C_n
                       Cluster_result{i}(1) = cl_c_p_p(i);
                   end
                   
                   complete_n = cl_c_p_p;
                   
%                    clear cl_c_p_p;
                   
                   l_complete = length(complete_n);
                   
                   while l_complete ~= l
                       clear test_c
                       for i = 1 : C_n
                           test_c{i}(1,:) = zeros(1,l);
                           
                           for j = 1 : length(Cluster_result{i})
                               test_c{i}(1,:) = test_c{i}(1,:)+CM(Cluster_result{i}(j),:);
                           end
                           
                           test_c{i}(1,:) = test_c{i}(1,:)/length(Cluster_result{i});
                           
                           test_c{i}(2,:) = 1:1:l;
                       end
                       
                       
                       
                       for i = 1 : C_n
                           test_c{i}(:,complete_n) = [];
                           [a b] = max(test_c{i}(1,:));
                           test_c2{i}(1) = a;
                           test_c2{i}(2) = test_c{i}(2,b);
                           test_c3(i) = a;
                       end
                       
                       
                       [a b] = max(test_c3);
                       
                       index = test_c2{b}(2);
                       Cluster_result{b}(length(Cluster_result{b})+1) = index;
                       complete_n(length(complete_n)+1) = index;
                       
                       l_complete = length(complete_n);
                       clear test_c
                       
                   end
                   
                   
                   
                   Cluster_result_s{C_n} = Cluster_result;
                   
                   
                   clear Cluster_result
                   
               end
           end
           
           
           
       else
           
           for C_n = 2 : m_C_n
              
               if C_n == 2
                   
                   h = 1;
                   m1 = 0;
                   m2 = 10^-10;
%                    
                       while(m1 < m2)
                           m1 = m2;
                        
                           for r = (h-1)*li2+1 : h*li2
                               
                               [a1 b1] = min(min(CM_test3));
                               [a2 b2] = min(CM_test3(b1,:));

                               min_corr(r) = a2;

                               CM_test3(b1,b2) = 1;
                               CM_test3(b2,b1) = 1;

                               cl_c_p_all{r}(1) = b1;
                               cl_c_p_all{r}(2) = b2;

                               cl_c_p_p = cl_c_p_all{r};
                               %                    cl_c_p_p = [200 100];
                               for i = 1 : C_n
                                   Cluster_result{i}(1) = cl_c_p_p(i);
                               end

                               complete_n = cl_c_p_p;
                               l_complete = length(complete_n);
                               while l_complete ~= l
                                   %         for re = 1 : 1
                                   for i = 1 : C_n
                                       test_c{i}(1,:) = zeros(1,l);

                                       for j = 1 : length(Cluster_result{i})
                                           test_c{i}(1,:) = test_c{i}(1,:)+CM(Cluster_result{i}(j),:);
                                       end

                                       test_c{i}(1,:) = test_c{i}(1,:)/length(Cluster_result{i});

                                       test_c{i}(2,:) = 1:1:l;
                                   end



                                   for i = 1 : C_n
                                       test_c{i}(:,complete_n) = [];
                                       [a b] = max(test_c{i}(1,:));
                                       test_c2{i}(1) = a;
                                       test_c2{i}(2) = test_c{i}(2,b);
                                       test_c3(i) = a;
                                   end


                                   [a b] = max(test_c3);

                                   index = test_c2{b}(2);
                                   Cluster_result{b}(length(Cluster_result{b})+1) = index;
                                   complete_n(length(complete_n)+1) = index;

                                   l_complete = length(complete_n);
                                   clear test_c

                               end

                               for i = 1 : C_n
                                   CM_m{i} = CM_test1(Cluster_result{i},Cluster_result{i});
                               end
                               sum_all_cm_pur = sum(sum(CM_test1));
                               n_ele_all = length(CM_test1)^2;

                               n_ele_clu = 0;

                               for i = 1:length(Cluster_result)

                                   n_ele_clu = n_ele_clu + length(Cluster_result{i})^2;
                                   sum_clu_p(i) = sum(sum(CM_m{i}));
                                   var_clu_p(i) = var(double(CM_m{i}(:)));
                               end

                               sum_clu = sum(sum_clu_p);
                               mean_clu(r) = (sum_clu-l)/(n_ele_clu-l);

                               n_ele_not_clu = n_ele_all - n_ele_clu;

                               sum_not_clu = sum_all_cm_pur - sum_clu;

                               mean_not_clu(r) = sum_not_clu/n_ele_not_clu;
                               var_clu(r) = mean(var_clu_p);


                               Cluster_result_spp{r} = Cluster_result;


                               clear Cluster_result
                           end
                            pur = (mean_clu-mean_CM_all).*(mean_CM_all-mean_not_clu); 
                           
                            [m2 i2_p] = max(pur);
                            
                            Cluster_result_sp_m(h) = Cluster_result_spp(i2_p);
                            cl_c_p_m(h) = cl_c_p_all(i2_p);
                            
                            clear cl_c_p_all mean_clu mean_not_clu var_clu
                            
                            h = h + 1;
                           
                           
                       end


                        cl_c_p_p = cl_c_p_m{h-2};
                        Cluster_result_s{C_n} =  Cluster_result_sp_m{h-2};
                   
                  
               else
                   
   
                    CM_test4 = CM;
                    for i = 1 : length(cl_c_p_p)
                        for j = 1 : length(cl_c_p_p)
                            if i ~= j
                                CM_test4(cl_c_p_p(i),cl_c_p_p(j)) = 1;
                                CM_test4(cl_c_p_p(j),cl_c_p_p(i)) = 1;
                            end
                        end
                   end
                     

                 
                   
                   

                       t = zeros(1,l);
                       for i = 1 : length(cl_c_p_p)
                           t = t + CM_test4(cl_c_p_p(i),:);
                       end
                       
                       
                       
                       [a3 b3] = sort(t);
                       
                       b3_2 = b3((h-1)*li3+1 : h*li3);
                       
                       
                       for i = 1 : length(cl_c_p_p)
                           CM_test4(cl_c_p_p(i),b3_2) = 1;
                           CM_test4(b3_2,cl_c_p_p(i)) = 1;
                       end
                       
                       
                         h = 1;
                         m1 = 0;
                         m2 = 10^-7;
                   
%                        
%                        
%                        
                       while (m1 < m2)
                       
                       m1 = m2;
%                        
                       for r = 1 : li3
                           clear Cluster_result
                           
                            cl_c_p_p2 = cl_c_p_p;
                            cl_c_p_p2(C_n) = b3_2(r);
                            cl_c_p_all{C_n}{r} = cl_c_p_p2;
                            for i = 1 : C_n
                               Cluster_result{i}(1) = cl_c_p_p2(i);
                            end

                             complete_n = cl_c_p_p2;
                             l_complete = length(complete_n);
                            
                             
                         while l_complete ~= l
                             for i = 1 : C_n
                               test_c{i}(1,:) = zeros(1,l);

                               for j = 1 : length(Cluster_result{i})
                                   test_c{i}(1,:) = test_c{i}(1,:)+CM(Cluster_result{i}(j),:);
                               end

                               test_c{i}(1,:) = test_c{i}(1,:)/length(Cluster_result{i});

                               test_c{i}(2,:) = 1:1:l;
                             end
                           
                         for i = 1 :C_n
                           test_c{i}(:,complete_n) = [];
                           [a b] = max(test_c{i}(1,:));
                           test_c2{i}(1) = a;
                           test_c2{i}(2) = test_c{i}(2,b);
                           test_c3(i) = a;
                         end
                             
                         [a b] = max(test_c3);
                       
                       index = test_c2{b}(2);
                       Cluster_result{b}(length(Cluster_result{b})+1) = index;
                       complete_n(length(complete_n)+1) = index;
                       
                       l_complete = length(complete_n);
                       clear test_c test_c2 test_c3

                         
                         end
                      
                         for i = 1 :C_n
                       CM_m{i} = CM_test1(Cluster_result{i},Cluster_result{i});
                   end
                   sum_all_cm_pur = sum(sum(CM_test1));
                   n_ele_all = length(CM_test1)^2;
                   
                   n_ele_clu = 0;
                   
                   for i = 1:length(Cluster_result)
                       
                       n_ele_clu = n_ele_clu + length(Cluster_result{i})^2;
                       sum_clu_p(i) = sum(sum(CM_m{i}));
                       var_clu_p(i) = var(double(CM_m{i}(:)));
                       
                   end
                   sum_clu = sum(sum_clu_p);
                   mean_clu(r) = (sum_clu-l)/(n_ele_clu-l);
                   
                   n_ele_not_clu = n_ele_all - n_ele_clu;
                   
                   sum_not_clu = sum_all_cm_pur - sum_clu;
                   
                   mean_not_clu(r) = sum_not_clu/n_ele_not_clu;
                    var_clu(r) = mean(var_clu_p);
               
                   
                   Cluster_result_spp{r} = Cluster_result;
                   
                   
                   
                   clear Cluster_result CM_m
                   
                       end 
                       
                        pur = (mean_clu-mean_CM_all).*(mean_CM_all-mean_not_clu); 
                        
                        [m2 i2_p] = max(pur);
                            
                       Cluster_result_sp_m(h) = Cluster_result_spp(i2_p);
                       cl_c_p_m(h) = cl_c_p_all{C_n}(i2_p);
                            
                        clear cl_c_p_all mean_clu mean_not_clu var_clu
                            
                      h = h + 1;
                       end
                      
                        cl_c_p_p = cl_c_p_m{h-2};
                        Cluster_result_s{C_n} =  Cluster_result_sp_m{h-2};
                        
                         clear cl_c_p_m cl_c_p_all

                   
                             



               end
               
               
               
               
               
               
           end
           
           
           
           
           
           
           
           
       end
           
           
           

           
           


Cluster_result_c = Cluster_result_s;


cl_c_p = cl_c_p_p;
   
        


%% 중심 노드 보정

clear mean_clu mean_not_clu mean_corr Cluster_result_s

for C_n = 2 : m_C_n
%     C_n
%     Cluster_result_c
%     Cluster_result_c{1}
%     Cluster_result_c{2}
    
     for i = 1 : C_n
        CM_m{i} = CM_test1(Cluster_result_c{C_n}{i},Cluster_result_c{C_n}{i});
        [cen1 cen2] = max(sum(CM_m{i}));
        cen_m(i) = Cluster_result_c{C_n}{i}(cen2);
     end
    
      tf = isequal(cl_c_p(1:C_n),cen_m);
      hh = 0;
      
      if tf == 1 
          Cluster_result = Cluster_result_c{C_n};
          
      else
        while tf < 1
            clear CM_m
 
          hh = hh + 1;
          
          
        clear test_c test_c2 test_c3 Cluster_result test
        cen{1} = cen_m;
        cl_c_p2 = cen_m;
        l_complete = length(cen_m);
        complete_n = cl_c_p2;
        
%         cl_c_p2
        
        for i = 1 : C_n
            Cluster_result{i}(1) = cl_c_p2(i);
        end


        while l_complete ~= l
            %     for re = 1: 3
            for i = 1 : C_n
                test_c{i}(1,:) = zeros(1,l);

                for j = 1 : length(Cluster_result{i})
                    test_c{i}(1,:) = test_c{i}(1,:)+CM(Cluster_result{i}(j),:);
                end

                test_c{i}(1,:) = test_c{i}(1,:)/length(Cluster_result{i});

                test_c{i}(2,:) = 1 : 1 : l;
            end



            for i = 1 : C_n
                test_c{i}(:,complete_n) = [];
                [a b] = max(test_c{i}(1,:));
                test_c2{i}(1) = a;
                test_c2{i}(2) = test_c{i}(2,b);
                test_c3(i) = a;
            end


            [a b] = max(test_c3);

            index = test_c2{b}(2);
            Cluster_result{b}(length(Cluster_result{b})+1) = index;
            complete_n(length(complete_n)+1) = index;

            l_complete = length(complete_n);
            clear test_c

        end
  

        
        Cluster_result2 = Cluster_result;

       
        for i = 1 : C_n
            CM_m{i} = CM_test1(Cluster_result2{i},Cluster_result2{i});
            [cen1 cen2] = max(sum(CM_m{i}));
            cen_m2(i) = Cluster_result2{i}(cen2);
         end
             for i = 1: length(cen)
                 tf_p(i) = isequal(cen{i},cen_m2); 
             end
         cen{length(cen)+1} = cen_m2;
        cen_m = cen_m2;
         tf = sum(tf_p);
         
        end
      end
        
        
        clear tf_p cen 
    

    sum_all_cm_pur = sum(sum(CM_test1));
    
    n_ele_all = length(CM_test1)^2;
    
    n_ele_clu = 0;
    
    
    for i = 1:length(Cluster_result)
       
        n_ele_clu = n_ele_clu + length(Cluster_result{i})^2;
        sum_clu_p(i) = sum(sum(CM_m{i}));
                
    end
    
    Cluster_result_final_p{C_n} = Cluster_result;
    for i = 1: length(Cluster_result)

        for j = 1:length(Cluster_result{i})
            Cluster_result_n{C_n}{i}(j) = n_m(Cluster_result{i}(j));
            Cluster_result_g{C_n}{i}(j) = n_m_g(Cluster_result_n{C_n}{i}(j));
            
        end
        
        

   end
    
    
    clear Cluster_result
    
    sum_clu = sum(sum_clu_p);
    mean_clu(C_n) = (sum_clu-l)/(n_ele_clu-l);
    
    n_ele_not_clu = n_ele_all - n_ele_clu;
    
    sum_not_clu = sum_all_cm_pur - sum_clu;
    
    mean_not_clu(C_n) = sum_not_clu/n_ele_not_clu;
%     Cluster_result_n{C_n}
end



mean_clu(1) = mean_CM_all;
mean_not_clu(1) = 0;

 
pur = (mean_clu-mean_CM_all).*(mean_CM_all-mean_not_clu); 
pur(1) = -100;

[mp1 mp2] = max(pur);


result1 = Cluster_result_n{mp2};
% pur
% mp2
% result1

result2 = Cluster_result_g{mp2};
end
end


