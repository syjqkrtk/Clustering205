clc
clear all
close all

temp = load('CM_205_2018');

CM_all = temp.DM;
group = temp.group205;

g1 = find(group == 2);
d = [9 12 14 15 35 91];
% for i = 1 : length(d)
% g1(g1 == d(i)) = [];
% end
g2 = find(group == 7);
g2(11) = [];
g2(1) = [];

g2 = g2(1:10);

g = [g1 g2];


for w = 3 : 6
   
    CM = CM_all{w};
    
    m_CM = mean(mean(CM));
    v_CM = var(reshape(CM,1,length(CM).^2));
    
    CM_g = CM(g,g);
    
    c1 = mean(mean(CM(g1,g1)));
    c2 = mean(mean(CM(g2,g2)));
    
    c3 = mean(mean(CM(g1,g2)));
    
    a(w) = abs(m_CM-c3)
    
end
% temp = load('pdf205_w110');
% 
% 
% for w = 1 : 10
%    w
%    ffp_p = temp.pdf{w};
%    for i = 6 : 6
%    ffp = ffp_p{i};
%    ffp_all_pp = ffp(0<ffp);
%    information_result_p = -sum(log2(ffp_all_pp));
%    end
%    
%    information_result(w) = information_result_p;
%    clear information_result_p
%     
% end
% 
% plot(information_result)