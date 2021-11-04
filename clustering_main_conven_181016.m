clc
clear all
close all

temp = load('CM_205_2018');
% temp_r = load('../../genus_r_result_0523');
temp_c = load('color_index_genus');
% temp_t = load('../../genus_info');
% r_result = temp_r.r_result;

mkdir result1018_c
folder_name = 'result1018_c\';

for w = 1: 1
    tic    
w
clearvars -except w temp r_result temp_c temp_t out_p_pro out_p_conven folder_name

w_c = num2str(w);
% 
CM = 1-temp.DM{w};

        n_m_g = temp.group205;
        mi_name_selection = temp.name205;
% 
temp5 = load('pdf205_w110');

pdf_p = temp5.pdf{w};

clear pdf_pp

for i = 1 : length(pdf_p)
    pdf_pp(i,:) = pdf_p{i};
end

    tt0 = ['240_result_','w=',w_c,];

close all


pdf = pdf_pp;

       m_C_n = round(205/5);
       if m_C_n <= 2
           m_C_n = 2;
       end

for i = 2 : m_C_n
     clear idx_p
     idx_p= kmeans(pdf,i);
     
     for j = 1 : i
        fp{j} = find(idx_p == j);
        clustering_result_pp{i}{j} = fp{j};
        CM_intra{j} = CM(fp{j},fp{j});
     end
     
     h = 1 ;
     
     for j = 1 : i
         for k = 1 : length(fp{j})
        
%              fp{j}(k) 
             w_p{j}{k} = CM_intra{j}(k,:);
             
             if length(w_p{j}{k}) == 1
                 s(h) = 1;
                 
             else
             w_p{j}{k}(k) = [];
             w_value_p{j}(k) = mean(w_p{j}{k});
             w_value(h) = w_value_p{j}(k);
             
             hh = 1;
             
             for jj = 1 : i
                if jj ~= j
                B_p{j}{k}{hh} = CM(fp{j}(k),fp{jj});
                B_value_p{j}{k}(hh) = mean(B_p{j}{k}{hh});
                hh = hh + 1;
                end
             end
             
             b_value_p{j}(k) = min(B_value_p{j}{k});
             b_value(h) = b_value_p{j}(k);
             
             s_pp{h} = [b_value(h) w_value(h)];
             s_p(h) = max(s_pp{h});
             
             s(h) = (b_value(h)-w_value(h))/s_p(h);
             

             
             
                
             end
    
              h = h + 1;
             
         end
     end
     
     silhouettes(i) = mean(s);
     clear s
         
     
    
end

[ma mb] = max(abs(silhouettes));

cg_result1 = clustering_result_pp{mb};

cg_result_conven = cg_result1;

k = 1;

for i = 1: length(cg_result1)

            for j = 1:length(cg_result1{i})
                cc_ss_result(k) = cg_result1{i}(j);
                group_re2(k) = n_m_g(cg_result1{i}(j));
                k = k + 1;
            end



end

    S1 = 1 : 1 : length(cc_ss_result);
      S2 = cc_ss_result;
      
      clear cc_ss_result

      CM_re = mat_ind_change(CM, S2, S1);
      Chart_name = ' ';
      
        figure(1)
      

      visual_dmat(1-CM_re,Chart_name)

      al = 0.5; 

    for i = 1 : length(cg_result1)

        l3 = length(cg_result1{i});
        al_t(i) = al;
        l_t(i) = l3;

        rectangle('Position',[al al l3 l3],'EdgeColor','r','LineWidth',2);
        hold on;

        al = al + l3;
    end

    
    tt1 = [folder_name,tt0,'_corr.png'];

   saveas(gcf,tt1)

 
    group_m = zeros(length(group_re2),length(group_re2));
    group_o = zeros(length(group_re2),length(group_re2));
    
    temp_c = load('color_index_genus');
    
    color_sp = temp_c.color_s;
    
    
     color_s = color_sp(1:10,:);
   gg_s = 10;

    for i = 1: length(group_re2)
        for j = 1:length(group_re2)

            if group_re2(i) == group_re2(j)
                group_m(i,j) = group_re2(i);
                group_o(i,j) = 1;
                group_c1(i,j) = color_s(group_re2(i),1);
                group_c2(i,j) = color_s(group_re2(i),2);
                group_c3(i,j) = color_s(group_re2(i),3);
                
%                 group_m2(i,j) = group_re2(i)*0.1;
            end
        end
    end
    
    group_o_result{w} = group_o;
    
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
        yt(i) = 5*(i-1)+3;
    end
    
       temp_g = load('group_name');
   genus_name_selection_p = temp_g.group;
    
   temp_g = load('group_name');
   genus_name_selection_p = temp_g.group;
   for i  = 1 : gg_s
       
       genus_name_selection{i} = genus_name_selection_p{i};
       
   end
       
text(xt,yt,genus_name_selection,'fontsize',8);

set(gca,'ytick',[])
set(gca,'xtick',[])



tt3 = [folder_name,tt0,'_mi_labeling.png'];


saveas(gcf,tt3)


close all
clear cg_result1

if gg_s ~= 1

for x = 2 : 2
    if x == 1
      cg_result1  = cg_result_pro;
    else
       cg_result1  = cg_result_conven;
    end
    
    for i = 1: length(cg_result1)
   
    cluster_info(cg_result1{i}) = i;
    
end


cc = nchoosek(1:gg_s,2); 
[llc1 llc2] = size(cc);

for i = 1 : llc1
    cc1 = cc(i,:);
    a1 = find(group_re2 == cc1(1));
    a2 = find(group_re2 == cc1(2));
    
    a1_c = cluster_info(a1);
    a2_c = cluster_info(a2);
    
    h1 = mode(a1_c);
    h1_f = length(find(a1_c == h1));
    
    
    h2 = mode(a2_c);
    h2_f = length(find(a2_c == h2));
    
    if h1 ~= h2
       
        out_p_p(i) = (h1_f+h2_f)/(length(a1)+length(a2));
    elseif h1 == h2
        
        oo = max([h1_f h2_f]);
        
        out_p_p(i) = oo/(length(a1)+length(a2));
    end
    
    
end

    out_p(x) = mean(out_p_p);
    

    
    if x == 1
      out_p_conven{w}(oi) = mean(out_p_p);
    else
      out_p_conven(w)  =  mean(out_p_p);
    end
    
    clear out_p_p
end    
    

 
 
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

toc
end

% 
%         
% 
% save([folder_name,'out_p_pro.mat'],'out_p_pro','out_p_conven')
% 
%         
% 


