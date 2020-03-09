clear all;
%close all;
clc;
format long;

hold on

Q=dlmread('Actionrandom.txt')/6;
P=dlmread('Actionorder.txt')/6;
figure
%plot(Q)
plot([1:1000],Q,[1:1000],P)
%semilogx([1:500],Q(1:500),[1:500],P(1:500))
% 
% x = 1:500;
% l(1:9)=Q(1:9);
% u(1:9)=P(1:9);
% l(10:19)=Q(10:10:100)
% u(10:19)=P(10:10:100)
% l(20:23)=Q(110:100:500)
% u(20:23)=P(110:100:500)
% y(1:9)=x(1:9)
% y(10:19)=x(10:10:100)
% y(20:23)=x(110:100:500)
% figure
% scatter(y,l)
% hold on
% scatter(y,u)
% set(gca,'xscale','log')