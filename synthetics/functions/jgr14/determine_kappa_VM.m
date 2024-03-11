% how large is kappa in a von mises distribution, if x% of samples should
% be beyond +/- delta
clear all

%sigma = 30;                  % Nodal plane uncertainty from catalog;
%kappa = 1/(deg2rad(sigma)^2)   % equivalent to 1/sigma^2 in normal distribution
kappa = 1.31;                 % equivalent to 1/sigma^2 in normal distribution

delta = 30;                 % error threshold from FM generation in [deg]

addpath(genpath('../../fct_downloads/CircStat2010e/'));
n = 1E5;

%x=linspace(-pi,pi,n);

% Generate $n random samples that follow a VM distribution
p_stoch=circ_vmrnd(0,kappa,n);
p_stoch = rad2deg(p_stoch);


% % Analytic VM-distribution for comparison, scaled approximately to
% % stochastic VM-distribution
% p_anal = 800*circ_vmpdf(x,0,kappa);

figure(123)
clf
x = (-180:1:180);
[h,xbin] = hist(p_stoch,x);
bar(xbin,h)
hold on
line([delta,delta],[0,n/50],'color','r')
line([-delta,-delta],[0,n/50],'color','r')
title(['Kappa = ',num2str(kappa)])
tail_up  = p_stoch >  delta;
tail_low = p_stoch < -delta;

tail_p100 = (sum(tail_up) + sum(tail_low))/n*100