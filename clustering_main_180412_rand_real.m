clc
clear all
close all

temp = load('../../CM1_9_info');
temp_r = load('../../genus_r_result_0523');
temp_c = load('../../color_index_genus');
temp_t = load('../../genus_info');
r_result = temp_r.r_result;

mkdir result0524
folder_name = 'result0524\';

for w = 1: 10
w
clearvars -except w temp r_result temp_c temp_t out_p_pro out_p_conven folder_name

w_c = num2str(w);
% 

% 
CM_pp = temp.CM{w};
matrix_mi_name = temp.matrix_mi_name;



genus_name_real = temp_t.genus_name_real;
genus_info_real = temp_t.genus_info_real;


v_r = r_result(w);

% ii = [18 20 68 75 79 84 184]; 

ii = [18 20];

gg_s = length(ii);

genus_selection = ii;
genus_name_selection = genus_name_real(ii);
% 
     selection_index_id = [];
     n_m_g = [];
     mi_name_selection = [];
     
    for i = 1 : length(genus_selection)

        a = find(genus_info_real == genus_selection(i));
        selection_index{i} = a;
        selection_index_id = [selection_index_id a];
        n_m_g(length(n_m_g)+1:length(n_m_g)+length(a)) = i;
        mi_name_selection = [mi_name_selection matrix_mi_name(a)];

    end
    CM = CM_pp(selection_index_id,selection_index_id);

    
    
    
    
    l = length(n_m_g);
    l_all = l;
    n_m = 1 : 1 : l;

    l = length(n_m_g);
    l_all = l;
    n_m = 1 : 1 : l;

    [out check cg_result1 cg_result_g] = f_all_clustering_process_180404(CM,n_m,n_m_g,v_r);
    
    for i = 1: length(cg_result1)
        pre_index(i) = check;
        out_index(i) = 0;
        
    end
    
    if out == 1
        out_index(length(out_index)) = 1;
    end
        
    
    r = 1;
    
    while sum(pre_index) > 0
%      for j = 1 : 1
%         clear pre_index
        
        r
         for i = 1 : length(cg_result1)
             
             if pre_index(i) > 0
             
                CM_pp1{i} = CM(cg_result1{i},cg_result1{i});
               
                [out_index_p(i) pre_index_p(i) cg_result1_p{i} cg_result_g_p{i}] = f_all_clustering_process_180404(CM_pp1{i},cg_result1{i},n_m_g,v_r);
            
                
             else
                         out_index_p(i) = out_index(i);
                         pre_index_p(i) = 0;
                         cg_result1_p{i} = cg_result1{i};
                        cg_result_g_p{i} = n_m_g;
   
             end
         end
         
         
         clear pre_index cg_result1 cg_result_g CM_pp1
         
         h = 1;
         
        
         
         for i = 1 : length(cg_result1_p) 
            
             if pre_index_p(i) > 0
            
                 for k = 1 : length(cg_result1_p{i})

                     cg_result1{h} = cg_result1_p{i}{k};
                     pre_index(h) = pre_index_p(i);
                     cg_result_g{h} = cg_result_g_p{i}{k};
                     h = h + 1;
                 end
             
             else
                 
                
                 
                     cg_result1{h} =  cg_result1_p{i};
                     pre_index(h) = pre_index_p(i);
                     cg_result_g{h} = cg_result_g_p{i};
                     h = h + 1;
                
             
                     
             end
         end
         
         for i = 1 : length(out_index)
             if pre_index_p(i) == 0
                 
                 out_index_a{i} = out_index(i);
             else
             
             if out_index(i) == 0
                 
                 if out_index_p(i) == 0
                 
                         out_index_a{i} = zeros(1,length(cg_result1_p{i}));
                    
                 else
                     
                     out_index_a{i} = zeros(1,length(cg_result1_p{i}));
                     out_index_a{i}(length(out_index_a{i})) = 1;
                     
                 end
                 
                 
             else
                 
                  out_index_a{i} = zeros(1,length(cg_result1_p{i}))+1;
                 
             end         
             
             end
             
             
         end
         
         clear out_index
         out_index = [];
         for i = 1 : length(out_index_a)
            out_index  = [out_index out_index_a{i}];
         end
        
         r = r + 1;
%       
         clear pre_index_p cg_result1_p cg_result_g_p out_index_a


    end

    

  







%%  


cg_result_p = cg_result1;

clear cg_result1

aa = find(out_index == 1);

cg_result_out_p = cg_result_p(aa);

cg_result_out = [];
for i = 1 : length(cg_result_out_p)
   
    cg_result_out = [cg_result_out cg_result_out_p{i}];
    
end

cg_result_p(aa) = [];

for i = 1:length(cg_result_p)
    l_cg(i) = length(cg_result_p{i});
end
    
[lcs1 lcs2] = sort(l_cg,'descend');
%     
for i = 1:length(cg_result_p)
	cg_result1{i} = cg_result_p{lcs2(i)};
end

cg_result1{length(cg_result1)+1} = cg_result_out;

  
k = 1;

for i = 1: length(cg_result1)
            
            c_size(i) = length(cg_result1{i});
            for j = 1:length(cg_result1{i})
                cc_ss_result(k) = cg_result1{i}(j);
                group_re_pro(k) = n_m_g(cg_result1{i}(j));
                k = k + 1;
            end



end
 
S1 = 1 : 1 : length(cc_ss_result);
S2 = cc_ss_result;

CM_re = mat_ind_change(CM, S2, S1);
Chart_name = ' ';

figure(1)

visual_dmat(CM_re,Chart_name)

    al = 0.5; 

    for i = 1 : length(cg_result1)

        l3 = length(cg_result1{i});
        al_t(i) = al;
        l_t(i) = l3;

        rectangle('Position',[al al l3 l3],'EdgeColor','r','LineWidth',2);
        hold on;

        al = al + l3;
    end
%     gg_s_c = int2str(gg_s);
    v_r_c = int2str(v_r*1000);
    tt0 = ['240_result_','w=',w_c,', r=',v_r_c];
%         title([tt3,'%']);
%     
   tt1 = [folder_name,tt0,'_corr.png'];

   saveas(gcf,tt1)

%%

    group_m = zeros(length(group_re_pro),length( group_re_pro));


    
    color_sp = temp_c.color_s;
    
    
     color_s = color_sp(ii,:);

    for i = 1: length(group_re_pro)
        for j = 1:length(group_re_pro)

            if group_re_pro(i) == group_re_pro(j)
                group_m(i,j) = group_re_pro(i);
                group_c1(i,j) = color_s(group_re_pro(i),1);
                group_c2(i,j) = color_s(group_re_pro(i),2);
                group_c3(i,j) = color_s(group_re_pro(i),3);
                
%                 group_m2(i,j) = group_re_pro(i)*0.1;
            end
        end
    end
    

    
    group_mi(:,:,1) = 1- group_c1;
    group_mi(:,:,2) = 1- group_c2;
    group_mi(:,:,3) = 1- group_c3;
          
figure(2)
          
    image(group_mi);
    axis image;
    
            % % % % % 클러스터 그리기
      al = 0.5; 

    for i = 1 : length(cg_result1)
        l3 = length(cg_result1{i});
        al_t(i) = al;
        l_t(i) = l3;

        rectangle('Position',[al al l3 l3],'EdgeColor','r','LineWidth',2);
        hold on;

        al = al + l3;
    end
    
   tt2 = [folder_name,tt0,'_mi_result.png'];

   saveas(gcf,tt2)

%     
    
%%


     labeling_p1 = zeros(5*gg_s,25);
     labeling_p2 = zeros(5*gg_s,25);
     labeling_p3 = zeros(5*gg_s,25);
     
    
     
     
    for i = 1 : gg_s
        labeling_p1(5*(i-1)+1:5*i,1:5) = color_s(i,1);
        labeling_p2(5*(i-1)+1:5*i,1:5) = color_s(i,2);
        labeling_p3(5*(i-1)+1:5*i,1:5) = color_s(i,3);
    end
    
     labeling_mi(:,:,1) = 1 -labeling_p1;
    
     labeling_mi(:,:,2) = 1 -labeling_p2;
    
     labeling_mi(:,:,3) = 1 -labeling_p3;
     
   figure(3)
    
    image(labeling_mi)
    
    axis image;
   
    po = [0 0];
    tt = 'red';
    
    xt = zeros(gg_s,1)+7;
    for i = 1 : gg_s
        yt(i) = 5*(i-1)+2.5;
    end
    
text(xt,yt,genus_name_selection,'fontsize',5);

set(gca,'ytick',[])
set(gca,'xtick',[])

    
   tt3 = [folder_name,tt0,'_mi_labeling.png'];

   saveas(gcf,tt3)
        
        
        
        %%

    
    [c_s1 c_s2] = sort(c_size,'descend');
    
    for i = 1 : length(c_s2)
       
        cluster_name_result(1:c_s1(i),c_s2(i)) = mi_name_selection(cg_result1{c_s2(i)});
        
    end
    
    xlswrite([folder_name,tt0,'_cluster_name'],cluster_name_result);
    
clear cg_result_g 
sum_l = 0;

 for i = 1: length(cg_result1)
   
    cg_result_g{i} = n_m_g(cg_result1{i});
    sum_l = sum_l + length(cg_result_g{i});
 end

 for g = 1: gg_s
    for cg = 1: length(cg_result_g)
       
        g_num_temp(cg) = length(find(cg_result_g{cg}==g));
        
        
    end    
    [group_num_info(g) group_mode(g)] = max(g_num_temp);
     group_num_real(g) = length(selection_index{g});
    clear g_num_temp
 end
 
 nc = nchoosek(1:gg_s,2);
 
 [s_n1 s_n2] = size(nc);
 
 for i = 1 : s_n1
    a1 = nc(i,1);
    a2 = nc(i,2);
    
    g1 = group_mode(a1);
    g2 = group_mode(a2);
    
    if g1 == g2
        
        gn = max([group_num_info(a1) group_num_info(a2)]);
        
        rand_p(i) = gn / (group_num_real(a1)+group_num_real(a2));
    else
        
         
        rand_p(i) = (group_num_info(a1)+group_num_info(a2))/(group_num_real(a1)+group_num_real(a2));
        
    end

   
     
     
     
     
 end
 
%     
    out_p_pro(w) = mean(rand_p);
    

 
 
% 
% gg_s2 = gg_s;
% 
% for i = 1 : length(cg_result_g)
%       
%     h1(i) = mode(cg_result_g{i});
%     c1(i) = length(find(cg_result_g{i} == h1(i)));
%     
% end
% 
% h_hist = hist(h1,1:gg_s);
% h_f1 = find(h_hist > 1);
% r = 1;
% 
% while length(h_f1) > 0
%    
%     for i = 1 : length(h_f1)
%        
%         h_ff = find(h1 == h_f1(i));
%         [h_fm1 h_fm2] = max(c1(h_ff));
%         
%         h_ff(h_fm2) = [];
%         
%         for j = 1 : length(h_ff)
%              cg_result_g{h_ff(j)}(cg_result_g{h_ff(j)} == h_f1(i)) = [];
%              if length(cg_result_g{h_ff(j)}) > 0
%              h1(h_ff(j)) = mode(cg_result_g{h_ff(j)});
%              c1(h_ff(j)) = length(find(cg_result_g{h_ff(j)} == h1(h_ff(j))));
%              else
%                   gg_s2 = gg_s2+1;
%                   h1(h_ff(j)) = gg_s2;
%                   c1(h_ff(j)) = 0;
%              end
%              
%              
%         end
%        
%         
%         
%     end
%    
%     h_hist = hist(h1,1: gg_s2);
%     h_f1 = find(h_hist > 1);
%     
%     close all
%     
% 
% end
% 
%     out_p_p = sum(c1)/sum_l;
%     
%     out_p_conven(w) = mean(out_p_p);
end

% 
%         
% 
% save([folder_name,'out_p_pro.mat'],'out_p_pro','out_p_conven')
% 
%         
% 


