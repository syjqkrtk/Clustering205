clc
clear all
close all

temp = load('pdf1443');

pdf1443 = temp.pdf18_1443;
mi_name1443 = temp.matrix_mi_name';

temp = load('name205_del');

name205 = temp.name;


for i = 1 : length(name205)
   
    
    a = find(double(ismember(mi_name1443,name205{i}))==1);
    
    if length(a) > 0
    index(i,1) = a(1);
    
    end
    
    
    
    
    
    
    
    
end