%% ENEL 671 Project3
% Pouyan Keshavarzian
% FALL 2016
%% Effect of Eignevalue Spread
clear all 
close all
clc

M = 10;
delta = 0.01;
N = 600;
K = 100;
h = [0.2194 1.0 0.2194;0.2798 1.0 0.2798;0.3365 1.0 0.3365;0.3887 1.0 0.3887];
lamda = 1;
for k=1:K
      a = BPSK(N);
      u = filterinput(a,h); 
      [e11, gamma_s] = RLSL_algorithm(M,N,lamda,delta,a,u(:,1));
      %e11 = e(11,:);
      %e11 = e11(:);
      ed11(:,k) = e11.^2;
      MSEE11 = sum(ed11,2)/K;
      MSEE11n = MSEE11./gamma_s(11,:)';
end
    figure(1)
    plot(1:N,gamma_s(11,:),'LineWidth',2)
    legend('Channel 1','Channel 2','Channel 3','Channel 4')
    grid on
    xlabel('Time (s)');
    ylabel('Gamma'); 
    title('Likilihood');

    figure(2)
    semilogy(1:N,MSEE11n,'LineWidth',2)
    title('Normalized MSEE');