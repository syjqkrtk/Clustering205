clc
clear all
close all

temp = load('pdf_w1');

pdf = temp.pdf;


temp = load('CM_205_2018');

n_m_g = temp.group205;
mi_name_selection = temp.name205;

unique_g = unique(n_m_g);

for i = 1 : 10
   
    
    a = find(n_m_g == i);
    
    
    group_name_p = mi_name_selection{a(1)};
    
    b = strfind(group_name_p,' ');
    
    group_name{i} = group_name_p(1:b(1)-1);
    
    
    
    
end


for i = 1 : 10
   
    a = find(n_m_g == i);
    
    for j = 1 : length(a)
       
        cg_ratio_p{i}(j) = pdf{a(j)}(2) + pdf{a(j)}(3);
        
        
    end
    
    cg_ratio_r(i) = mean(cg_ratio_p{i});
    
    
    
end


[cg_s1 cg_s2] = sort(cg_ratio_r);
group_s = group_name(cg_s2);

bar(1:10,cg_s1)
set(gca,'xtick',[1:10],'xticklabel',group_s )
xtickangle(90)
xlim([0.5 10.5])
ylabel('CG ratio of genus')