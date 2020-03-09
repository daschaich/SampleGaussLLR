clear all;
%close all;
clc;
format long;

hold on

Q=dlmread('Wilsonloop1x1.txt');
P=dlmread('Wilsonloop2x2.txt');
for i=1:41
    M(i)=sum(Q(((i-1)*30)+1:i*30))/30;
    N(i)=sum(P(((i-1)*30)+1:i*30))/30;
    x(i)=0.5+0.025*i;
end
hold on
% x = [1.006, 1.008, 1.010, 1.012, 1.014];
scatter(x,M)
hold on
scatter(x,N)