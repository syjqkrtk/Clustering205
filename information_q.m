clc
clear all
close all

temp = load('pdf205_w110');


for w = 1 : 10
   w
   ffp_p = temp.pdf{w};
   for i = 6 : 6
   ffp = ffp_p{i};
   ffp_all_pp = ffp(0<ffp);
   information_result_p = -sum(log2(ffp_all_pp));
   end
   
   information_result(w) = information_result_p;
   clear information_result_p
    
end

plot(information_result)