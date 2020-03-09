clear all;
%close all;
clc;
format long;

hold on

Q=dlmread('monopoledensity.txt');

for i=1:20
    M(i)=sum(Q(((i-1)*100)+1:i*100))/100;
    x(i)=0.01+0.1*i;
end
% 
% for int i=1:20
%     density(i)=sum(M(i,:));
% end
figure
scatter(x,M)
set(gca,'yscale','log')
ylim([10^(-3) 0.3])
xlim([0 2.1])