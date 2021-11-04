clc
clear all
close all

temp = load('proposed_out_p');

out_p_pro = temp.out_p_pro;

temp = load('conventional_out_p');

out_p_conven = temp.out_p_pro;

temp = load('kmeans_out_p');

out_p_kmeans = temp.out_p_pro;

w = 1 : 10;

subplot(1,2,1)

plot(w,out_p_pro,'ro-');
hold on
grid on
% 
% plot(w,out_p_conven,'bs-');
% plot(w,out_p_kmeans,'gx-');

ylim([0 1])

xlabel('Word length')
ylabel('Mean Rand measure for genus')
legend('Proposed LAC system')
% legend('Proposed LAC system','Silhouettes clustering','K-means + elbow method')


c1 = [9.31 8.32 7.99 7.91 8.24 9.5 14.297 41.504 202.11920 897.037988];

for i = 1: 10
c2(i) = 30*(1.1)^(i+1);
end

c3 = [3.25 3.38 3.40 3.45 3.65 6.15 10.8555 21.8886 197.56150 748.91648];

p1 = [85.099664 92.066267 86.767036 76.642242 77.61923 89.734130 87.575806 86.745136 89.753182 86.089803];

c = c1 + c2 + c3;

p = p1 + c2 + c3;

w = 1 : 10;
subplot(1,2,2)
plot(w,p,'ro-');
hold on
grid on

plot(w,c,'bs-');

% ylim([0 1])

xlabel('Word length')
ylabel('Execution time (s)')
legend('Proposed LAC system','Silhouettes clustering')