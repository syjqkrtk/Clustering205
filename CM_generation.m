clc
clear all
close all

temp = load('pdf205_w110');
pdf = temp.pdf;


for i = 1 : 10
   tic 
    
   
   for j = 1 : 205
       for k = 1 : 205
           
      
           CM{i}(j,k) = 1-abs(corr(pdf{i}{j}',pdf{i}{k}'));
           
       end
   end
   
    
   
   toc
end