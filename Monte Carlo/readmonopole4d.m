clear all;
%close all;
clc;
format long;

hold on

Q=dlmread('monopoledensity4d.txt');

for i=1:20
    M(i)=sum(Q(((i-1)*30)+1:i*30))/30;
    x(i)=1.02-0.0015*i;
end
% 
% for int i=1:20
%     density(i)=sum(M(i,:));
% end
figure
scatter(x,M)
% set(gca,'yscale','log')
% ylim([10^(-3) 0.3])
% xlim([0 2.1])